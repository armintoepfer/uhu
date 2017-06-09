// Copyright (c) 2017, Pacific Biosciences of California, Inc.
//
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted (subject to the limitations in the
// disclaimer below) provided that the following conditions are met:
//
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//
//  * Redistributions in binary form must reproduce the above
//    copyright notice, this list of conditions and the following
//    disclaimer in the documentation and/or other materials provided
//    with the distribution.
//
//  * Neither the name of Pacific Biosciences nor the names of its
//    contributors may be used to endorse or promote products derived
//    from this software without specific prior written permission.
//
// NO EXPRESS OR IMPLIED LICENSES TO ANY PARTY'S PATENT RIGHTS ARE
// GRANTED BY THIS LICENSE. THIS SOFTWARE IS PROVIDED BY PACIFIC
// BIOSCIENCES AND ITS CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED
// WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES
// OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
// DISCLAIMED. IN NO EVENT SHALL PACIFIC BIOSCIENCES OR ITS
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
// SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT
// LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF
// USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
// ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
// OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT
// OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF
// SUCH DAMAGE.

// Author: Armin Töpfer

#include <algorithm>
#include <atomic>
#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <numeric>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include <ssw_cpp.h>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamReader.h>
#include <pbbam/BamWriter.h>
#include <pbbam/Cigar.h>
#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>

#include <uhu/threadpool/ThreadPool.h>

namespace PacBio {
namespace Demux {
struct AlignerConfig
{
    AlignerConfig(uint8_t matchScore, uint8_t mismatchPenalty, uint8_t gapOpenPenalty,
                  uint8_t gapExtPenalty)
        : MatchScore(matchScore)
        , MismatchPenalty(mismatchPenalty)
        , GapOpenPenalty(gapOpenPenalty)
        , GapExtPenalty(gapExtPenalty)
    {
    }

    uint8_t MatchScore;
    uint8_t MismatchPenalty;
    uint8_t GapOpenPenalty;
    uint8_t GapExtPenalty;
};
struct Barcode
{
    Barcode(const std::string& name, const std::string& bases) : Name(name), Bases(bases) {}
    std::string Name;
    std::string Bases;
};

struct BarcodeHit
{
    BarcodeHit(int idx, int bq, int clipLeft, int clipRight)
        : IdxL(idx), IdxR(idx), Bq(bq), ClipStart(clipLeft), ClipRight(clipRight)
    {
    }
    BarcodeHit(int idxL, int idxR, int bq, int clipLeft, int clipRight)
        : IdxL(idxL), IdxR(idxR), Bq(bq), ClipStart(clipLeft), ClipRight(clipRight)
    {
    }
    uint16_t IdxL;
    uint16_t IdxR;
    uint8_t Bq;
    int ClipStart;
    int ClipRight;

    operator std::string() const
    {
        std::stringstream out;
        out << static_cast<int>(IdxL) << "\t" << static_cast<int>(IdxR) << "\t"
            << static_cast<int>(Bq) << "\t" << ClipStart << "\t" << ClipRight;
        return out.str();
    }

    friend std::ostream& operator<<(std::ostream& stream, const BarcodeHit& bh)
    {
        stream << std::string(bh);
        return stream;
    }
};

enum class Mode : int
{
    SYMMETRIC = 0,
    ASYMMETRIC,
    UNKNOWN
};

static Mode StringToMode(const std::string& mode)
{
    if (mode == "symmetric")
        return Mode::SYMMETRIC;
    else if (mode == "asymmetric")
        return Mode::ASYMMETRIC;
    return Mode::UNKNOWN;
}

static PacBio::CLI::Interface CreateCLI()
{
    using Option = PacBio::CLI::Option;

    PacBio::CLI::Interface i{"lima", "Demultiplex Barcoded CCS Data and Clip Barcodes", "0.3.0"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    // clang-format off
    i.AddOptions({});
    i.AddGroup("Barcode Configuration", {
        {"mode",      {"m","mode"},       "Barcoding mode. Available: symmetric",      Option::StringType("symmetric"), {"symmetric", "asymmetric"}},
        {"tryRC",     {"t","try-rc"},     "Try barcodes also as reverse complements.", Option::BoolType()},
    });
    i.AddGroup("Tuning", {
        {"windowSizeMult",  {"w","window-size-mult"},  "The candidate region size multiplier: barcode_length * multiplier.", Option::FloatType(1.2)},
        {"minScore",        {"s","min-score"},         "Minimum barcode score.",                                             Option::IntType(51)},
        {"minLength",       {"l","min-length"},        "Minimum sequence length after clipping.",                            Option::IntType(50)}
    });
    i.AddGroup("Aligner Configuration", {
        {"matchScore",      {"A", "match-score"},      "Score for a sequence match.",                           Option::IntType(2)},
        {"mismatchPenalty", {"B", "mismatch-penalty"}, "Penalty for a mismatch.",                               Option::IntType(2)},
        {"gapOpenPenalty",  {"O", "gap-open-penalty"}, "Gap open penalties for deletions and insertions.",      Option::IntType(3)},
        {"gapExtPenalty",   {"e", "gap-ext-penalty"},  "Gap extension penalties for deletions and insertions.", Option::IntType(1)}
    });
    // clang-format on

    i.AddPositionalArguments(
        {{"bam", "Source BAM", "BAM_FILE"}, {"fasta", "Barcode file", "FASTA_FILE"}});

    return i;
}

static void ParsePositionalArgs(const std::vector<std::string>& args,
                                std::vector<std::string>* datasetPaths,
                                std::vector<Barcode>* barcodes)
{
    std::vector<std::string> fastaPaths;
    for (const auto& i : args) {
        const bool fileExist = PacBio::Utility::FileExists(i);
        if (!fileExist) {
            throw std::runtime_error("File does not exist: " + i);
        }
        BAM::DataSet ds(i);

        switch (ds.Type()) {
            case BAM::DataSet::TypeEnum::SUBREAD:
            case BAM::DataSet::TypeEnum::ALIGNMENT:
            case BAM::DataSet::TypeEnum::CONSENSUS_ALIGNMENT:
            case BAM::DataSet::TypeEnum::CONSENSUS_READ:
                datasetPaths->push_back(i);
                break;
            case BAM::DataSet::TypeEnum::BARCODE:
            case BAM::DataSet::TypeEnum::REFERENCE:
                fastaPaths.push_back(i);
                break;
            default:
                throw std::runtime_error("Unsupported input file: " + i + " of type " +
                                         BAM::DataSet::TypeToName(ds.Type()));
        }
    }

    for (const auto& fasta : fastaPaths) {
        BAM::DataSet ds(fasta);
        for (const auto& fastaFile : ds.FastaFiles()) {
            BAM::FastaReader msaReader(fastaFile);
            BAM::FastaSequence f;
            while (msaReader.GetNext(f)) {
                barcodes->emplace_back(f.Name(), f.Bases());
            }
        }
    }
}

char Complement(char base)
{
    switch (base) {
        case 'A':
            return 'T';
        case 'a':
            return 't';
        case 'C':
            return 'G';
        case 'c':
            return 'g';
        case 'G':
            return 'C';
        case 'g':
            return 'c';
        case 'T':
            return 'A';
        case 't':
            return 'a';
        case '-':
            return '-';
        default:
            throw std::invalid_argument("invalid base!");
    }
}

static std::string RCerseComplement(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (auto it = input.crbegin(); it != input.crend(); ++it)
        output.push_back(Complement(*it));
    return output;
}

static StripedSmithWaterman::Alignment Align(StripedSmithWaterman::Aligner& aligner,
                                             const char* bases)
{
    StripedSmithWaterman::Filter filter;
    StripedSmithWaterman::Alignment alignment;
    aligner.Align(bases, filter, &alignment);
    return alignment;
};

static StripedSmithWaterman::Alignment AlignForward(StripedSmithWaterman::Aligner& aligner,
                                                    const Barcode& query)
{
    return Align(aligner, query.Bases.c_str());
};

static StripedSmithWaterman::Alignment AlignRC(StripedSmithWaterman::Aligner& aligner,
                                               const Barcode& query)
{
    auto revComp = RCerseComplement(query.Bases);
    return Align(aligner, revComp.c_str());
};

BarcodeHit SimdNeedleWunschAlignment(const AlignerConfig& ac, const std::string& target,
                                     const std::vector<Barcode>& queries, const Mode& mode,
                                     const bool& tryRC, const double& windowSizeMult)
{
    if (mode == Mode::UNKNOWN) throw std::runtime_error("Unsupported barcoding mode");

    int barcodeLength = 0;
    for (const auto& q : queries)
        barcodeLength = std::max(barcodeLength, static_cast<int>(q.Bases.size()));
    int barcodeLengthWSpacing = barcodeLength * windowSizeMult;
    int targetLength = target.size();

    StripedSmithWaterman::Aligner alignerLeft(ac.MatchScore, ac.MismatchPenalty, ac.GapOpenPenalty,
                                              ac.GapExtPenalty);
    alignerLeft.SetReferenceSequence(target.c_str(), std::min(targetLength, barcodeLengthWSpacing));

    StripedSmithWaterman::Aligner alignerRight(ac.MatchScore, ac.MismatchPenalty, ac.GapOpenPenalty,
                                               ac.GapExtPenalty);
    auto alignerRightBegin = std::max(targetLength - barcodeLengthWSpacing, 0);
    alignerRight.SetReferenceSequence(target.c_str() + alignerRightBegin,
                                      targetLength - alignerRightBegin);

    auto AlignTo = [&](StripedSmithWaterman::Aligner& aligner) {
        std::vector<int> scores(queries.size(), 0);
        std::vector<int> scoresRC(queries.size(), 0);

        for (size_t i = 0; i < queries.size(); ++i) {
            scores[i] = AlignForward(aligner, queries[i]).sw_score;
            scoresRC[i] = AlignRC(aligner, queries[i]).sw_score;
        }

        return std::make_pair(scores, scoresRC);
    };

    auto GetBestIndex = [&](std::vector<int>& v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });

        int bq = std::round(100.0 * v.at(idx.front()) / (barcodeLength * ac.MatchScore));
        return std::make_pair(idx.front(), bq);
    };

    if (mode == Mode::ASYMMETRIC) {

        auto Compute = [&](StripedSmithWaterman::Aligner& aligner, int* idx, int* score, int* clip,
                           bool left) {
            std::vector<int> scores;
            std::vector<int> scoresRC;
            std::tie(scores, scoresRC) = AlignTo(aligner);

            int idxFwd;
            int scoreFwd;
            std::tie(idxFwd, scoreFwd) = GetBestIndex(scores);

            int scoreRev;
            int idxRev;
            std::tie(idxRev, scoreRev) = GetBestIndex(scoresRC);

            if (scoreFwd > scoreRev) {
                *score = scoreFwd;
                *idx = idxFwd;
                if (left)
                    *clip = AlignForward(aligner, queries[*idx]).ref_end;
                else
                    *clip = alignerRightBegin + AlignForward(alignerRight, queries[*idx]).ref_begin;
            } else {
                *score = scoreRev;
                *idx = idxRev;
                if (left)
                    *clip = AlignRC(aligner, queries[*idx]).ref_end;
                else
                    *clip = alignerRightBegin + AlignRC(alignerRight, queries[*idx]).ref_begin;
            }
        };

        int leftIdx;
        int rightIdx;
        int leftScore;
        int rightScore;
        int clipLeft;
        int clipRight;
        Compute(alignerLeft, &leftIdx, &leftScore, &clipLeft, true);
        Compute(alignerRight, &rightIdx, &rightScore, &clipRight, false);

        return BarcodeHit(leftIdx, rightIdx, (leftScore + rightScore) / 2.0, clipLeft, clipRight);
    } else if (mode == Mode::SYMMETRIC && tryRC) {
        auto ComputeCombinedScore = [&](std::vector<int>* scores, std::vector<int>* scoresRC) {
            std::vector<int> scoresLeft;
            std::vector<int> scoresRCLeft;
            std::tie(scoresLeft, scoresRCLeft) = AlignTo(alignerLeft);

            std::vector<int> scoresRight;
            std::vector<int> scoresRCRight;
            std::tie(scoresRight, scoresRCRight) = AlignTo(alignerRight);

            assert(scoresLeft.size() == scoresRCRight.size());
            for (size_t i = 0; i < scoresLeft.size(); ++i)
                scores->emplace_back((scoresLeft.at(i) + scoresRCRight.at(i)) / 2.0);

            assert(scoresRCLeft.size() == scoresRCRight.size());
            for (size_t i = 0; i < scoresLeft.size(); ++i)
                scoresRC->emplace_back((scoresRCLeft.at(i) + scoresRight.at(i)) / 2.0);
        };

        std::vector<int> scores;
        std::vector<int> scoresRC;
        ComputeCombinedScore(&scores, &scoresRC);

        int scoreFwd;
        int idxFwd;
        std::tie(idxFwd, scoreFwd) = GetBestIndex(scores);
        int scoreRev;
        int idxRev;
        std::tie(idxRev, scoreRev) = GetBestIndex(scoresRC);

        int idx;
        int score;
        int clipLeft;
        int clipRight;
        if (scoreFwd > scoreRev) {
            score = scoreFwd;
            idx = idxFwd;
            clipLeft = AlignForward(alignerLeft, queries[idx]).ref_end;
            clipRight = alignerRightBegin + AlignRC(alignerRight, queries[idx]).ref_begin;
        } else {
            score = scoreRev;
            idx = idxRev;
            clipLeft = AlignRC(alignerLeft, queries[idx]).ref_end;
            clipRight = alignerRightBegin + AlignForward(alignerRight, queries[idx]).ref_begin;
        }

        return BarcodeHit(idx, score, clipLeft, clipRight);
    } else if (mode == Mode::SYMMETRIC && !tryRC) {
        auto ComputeCombinedForwardScore = [&](std::vector<int>* scores) {
            for (size_t i = 0; i < queries.size(); ++i) {
                (*scores)[i] = (AlignForward(alignerLeft, queries[i]).sw_score +
                                AlignRC(alignerRight, queries[i]).sw_score) /
                               ac.MatchScore;
            }
        };

        std::vector<int> scores(queries.size(), 0);
        ComputeCombinedForwardScore(&scores);

        int score;
        int idx;
        std::tie(idx, score) = GetBestIndex(scores);
        int clipLeft = std::max(0, AlignForward(alignerLeft, queries[idx]).ref_end);
        int clipRight = std::max(targetLength,
                                 alignerRightBegin + AlignRC(alignerRight, queries[idx]).ref_begin);

        return BarcodeHit(idx, score, clipLeft, clipRight);
    }

    throw std::runtime_error("Unknown barcoding mode");
}

static int Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM input, see --help" << std::endl;
        return EXIT_FAILURE;
    }

    const double windowSizeMult = options["windowSizeMult"];
    const bool tryRC = options["tryRC"];
    const std::string mode = options["mode"];
    const int minScore = options["minScore"];
    const int minLength = options["minLength"];

    AlignerConfig ac(options["matchScore"], options["mismatchPenalty"], options["gapOpenPenalty"],
                     options["gapExtPenalty"]);

    std::vector<std::string> datasetPaths;
    std::vector<Barcode> barcodes;
    ParsePositionalArgs(options.PositionalArguments(), &datasetPaths, &barcodes);

    auto BamQuery = [](const std::string& filePath) {
        BAM::DataSet ds(filePath);
        const auto filter = BAM::PbiFilter::FromDataSet(ds);
        std::unique_ptr<BAM::internal::IQuery> query(nullptr);
        if (filter.IsEmpty())
            query.reset(new BAM::EntireFileQuery(ds));
        else
            query.reset(new BAM::PbiFilterQuery(filter, ds));
        return query;
    };

    auto FilePrefixInfix = [](const std::string& path) -> std::string {
        size_t fileStart = path.find_last_of("/");

        if (fileStart == std::string::npos) fileStart = -1;

        // increment beyond the '/'
        ++fileStart;

        size_t extStart = path.substr(fileStart, path.length() - fileStart).find_last_of(".");

        if (extStart == std::string::npos) return "";

        auto suffix = path.substr(fileStart, extStart);
        return suffix;
    };

    std::unique_ptr<BAM::BamWriter> writer;
    std::map<int, BAM::BamRecord> map;
    for (const auto& datasetPath : datasetPaths) {
        auto query = BamQuery(datasetPath);
        std::vector<
            Uhu::Threadpool::ThreadPool::TaskFuture<std::tuple<BAM::BamRecord, std::string, bool>>>
            v;
        std::string prefix = FilePrefixInfix(datasetPath);
        std::atomic_int belowMinLength(0);
        std::atomic_int belowMinScore(0);
        std::atomic_int belowBoth(0);
        std::atomic_int aboveThresholds(0);
        for (auto& r : *query) {
            if (!writer) {
                writer.reset(new BAM::BamWriter(prefix + ".demux.bam", r.Header().DeepCopy()));
            }
            v.push_back(Uhu::Threadpool::DefaultThreadPool::submitJob(
                [&](BAM::BamRecord r) {
                    BAM::BamRecord recordOut;
                    std::string report;
                    BarcodeHit bh = SimdNeedleWunschAlignment(
                        ac, r.Sequence(), barcodes, StringToMode(mode), tryRC, windowSizeMult);
                    bool aboveMinLength = (bh.ClipRight - bh.ClipStart) >= minLength;
                    bool aboveMinScore = bh.Bq >= minScore;
                    report = r.FullName() + "\t" + std::string(bh);
                    if (aboveMinLength && aboveMinScore) {
                        r.Clip(BAM::ClipType::CLIP_TO_QUERY, bh.ClipStart, bh.ClipRight);
                        r.Barcodes(std::make_pair(bh.IdxL, bh.IdxR));
                        r.BarcodeQuality(bh.Bq);
                        recordOut = std::move(r);
                        ++aboveThresholds;
                    } else if (!aboveMinLength && !aboveMinScore) {
                        ++belowBoth;
                    } else if (!aboveMinLength) {
                        ++belowMinLength;
                    } else if (!aboveMinScore) {
                        ++belowMinScore;
                    }
                    return std::make_tuple(std::move(recordOut), report,
                                           aboveMinLength && aboveMinScore);
                },
                r));
        }

        std::ofstream report(prefix + ".demux.report");
        report << "ZMW\tBcLeft\tBcRight\tScore\tClipLeft\tClipRight" << std::endl;

        for (auto& item : v) {
            auto p = item.get();
            if (std::get<2>(p)) writer->Write(std::get<0>(p));
            report << std::get<1>(p) << std::endl;
        }

        std::ofstream summary(prefix + ".demux.summary");
        summary << "Above length and score threshold : " << aboveThresholds << std::endl;
        summary << "Below length and score threshold : " << belowBoth << std::endl;
        summary << "Below length threshold           : " << belowMinLength << std::endl;
        summary << "Below score threshold            : " << belowMinScore << std::endl;
        writer.reset(nullptr);
    }

    return EXIT_SUCCESS;
}
}
};

// Entry point
int main(int argc, char* argv[])
{
    return PacBio::CLI::Run(argc, argv, PacBio::Demux::CreateCLI(), &PacBio::Demux::Runner);
}