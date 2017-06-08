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
struct Barcode
{
    Barcode(const std::string& name, const std::string& bases) : Name(name), Bases(bases) {}
    std::string Name;
    std::string Bases;
};

struct BarcodeHit
{
    BarcodeHit(int idx, int bq, int clipStart, int clipEnd)
        : Idx(idx), Bq(bq), ClipStart(clipStart), ClipEnd(clipEnd)
    {
    }
    uint16_t Idx;
    uint8_t Bq;
    int ClipStart;
    int ClipEnd;

    operator std::string() const
    {
        std::stringstream out;
        out << static_cast<int>(Idx) << "\t" << static_cast<int>(Bq) << "\t" << ClipStart << "\t"
            << ClipEnd;
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
    SYMMETRIC_BOTH
};

static Mode StringToMode(const std::string& mode)
{
    if (mode == "symmetric")
        return Mode::SYMMETRIC;
    else if (mode == "symmetric_both")
        return Mode::SYMMETRIC_BOTH;
    return Mode::SYMMETRIC;
}

static PacBio::CLI::Interface CreateCLI()
{
    using Option = PacBio::CLI::Option;

    PacBio::CLI::Interface i{"demux_ccs", "Demultiplex Barcoded CCS Data and Clip Barcodes",
                             "0.1.0"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    // clang-format off
    i.AddOptions({
        {"mode", {"m","mode"},
         "Barcoding mode. Suffix \"_both\" indicates that barcodes are also tested as reverse complements."
         " Available: symmetric, symmetric_both",
         Option::StringType("symmetric"), {"symmetric", "symmetric_both"}},
        {"minScore",  {"s","min-score"},  "Minimum barcode score.", Option::IntType(51)},
        {"minLength", {"l","min-length"}, "Minimum sequence length after clipping.", Option::IntType(50)}
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

static std::string ReverseComplement(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (auto it = input.crbegin(); it != input.crend(); ++it)
        output.push_back(Complement(*it));
    return output;
}

BarcodeHit SimdNeedleWunschAlignment(const std::string& target, const std::vector<Barcode>& queries,
                                     Mode mode)
{
    int barcodeLength = queries.front().Bases.size();
    int barcodeLengthWSpacing = barcodeLength * 1.2;
    int targetLength = target.size();

    StripedSmithWaterman::Aligner alignerBegin;
    alignerBegin.SetReferenceSequence(target.c_str(),
                                      std::min(targetLength, barcodeLengthWSpacing));

    StripedSmithWaterman::Aligner alignerEnd;
    auto alignerEndBegin = std::max(targetLength - barcodeLengthWSpacing, 0);
    alignerEnd.SetReferenceSequence(target.c_str() + alignerEndBegin,
                                    targetLength - alignerEndBegin);
    StripedSmithWaterman::Filter filter;

    auto AlignForward = [&filter](StripedSmithWaterman::Aligner& aligner, const Barcode& query) {
        StripedSmithWaterman::Alignment alignment;
        aligner.Align(query.Bases.c_str(), filter, &alignment);
        return alignment;
    };

    auto AlignRC = [&filter](StripedSmithWaterman::Aligner& aligner, const Barcode& query) {
        StripedSmithWaterman::Alignment alignment;
        auto revComp = ReverseComplement(query.Bases);
        aligner.Align(revComp.c_str(), filter, &alignment);
        return alignment;
    };

    auto AlignTo = [&AlignForward, &AlignRC, &queries, &filter,
                    &barcodeLength](StripedSmithWaterman::Aligner& aligner) {
        std::vector<int> scores(queries.size(), 0);
        std::vector<int> scoresRev(queries.size(), 0);

        for (size_t i = 0; i < queries.size(); ++i) {
            scores[i] = AlignForward(aligner, queries[i]).sw_score;
            scoresRev[i] = AlignRC(aligner, queries[i]).sw_score;
        }

        return std::make_pair(scores, scoresRev);
    };

    auto GetBestIndex = [&barcodeLength](std::vector<int>& v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });

        int bq = std::round(100.0 * v.at(idx.front()) / (barcodeLength * 2.0));
        return std::make_pair(idx.front(), bq);
    };

    if (mode == Mode::SYMMETRIC_BOTH) {

        auto ComputeCombinedScore = [&AlignForward, &AlignRC, &AlignTo, &alignerBegin, &alignerEnd](
            std::vector<int>* scores, std::vector<int>* scoresRev) {
            std::vector<int> scoresBegin;
            std::vector<int> scoresRevBegin;
            std::tie(scoresBegin, scoresRevBegin) = AlignTo(alignerBegin);

            std::vector<int> scoresEnd;
            std::vector<int> scoresRevEnd;
            std::tie(scoresEnd, scoresRevEnd) = AlignTo(alignerEnd);

            assert(scoresBegin.size() == scoresRevEnd.size());
            for (size_t i = 0; i < scoresBegin.size(); ++i)
                scores->emplace_back((scoresBegin.at(i) + scoresRevEnd.at(i)) / 2);

            assert(scoresRevBegin.size() == scoresRevEnd.size());
            for (size_t i = 0; i < scoresBegin.size(); ++i)
                scoresRev->emplace_back((scoresRevBegin.at(i) + scoresEnd.at(i)) / 2);
        };

        std::vector<int> scores;
        std::vector<int> scoresRev;
        ComputeCombinedScore(&scores, &scoresRev);

        int forwardScore;
        int forwardIdx;
        std::tie(forwardIdx, forwardScore) = GetBestIndex(scores);
        int revScore;
        int revIdx;
        std::tie(revIdx, revScore) = GetBestIndex(scoresRev);

        int idx;
        int score;
        int clipStart;
        int clipEnd;
        StripedSmithWaterman::Alignment alignmentBegin;
        StripedSmithWaterman::Alignment alignmentEnd;
        if (forwardScore > revScore) {
            score = forwardScore;
            idx = forwardIdx;
            alignmentBegin = AlignForward(alignerBegin, queries[idx]);
            alignmentEnd = AlignRC(alignerEnd, queries[idx]);
        } else {
            score = revScore;
            idx = revIdx;
            alignmentBegin = AlignRC(alignerBegin, queries[idx]);
            alignmentEnd = AlignForward(alignerEnd, queries[idx]);
        }
        clipStart = alignmentBegin.ref_end;
        clipEnd = alignerEndBegin + alignmentEnd.ref_begin;

        return BarcodeHit(idx, score, clipStart, clipEnd);
    } else if (mode == Mode::SYMMETRIC) {

        auto ComputeCombinedForwardScore = [&AlignForward, &AlignRC, &queries, &alignerBegin,
                                            &alignerEnd](std::vector<int>* scores) {
            for (size_t i = 0; i < queries.size(); ++i) {
                (*scores)[i] = (AlignForward(alignerBegin, queries[i]).sw_score +
                                AlignRC(alignerEnd, queries[i]).sw_score) /
                               2;
            }
        };

        std::vector<int> scores(queries.size(), 0);
        ComputeCombinedForwardScore(&scores);

        int score;
        int idx;
        std::tie(idx, score) = GetBestIndex(scores);
        int clipStart = std::max(0, AlignForward(alignerBegin, queries[idx]).ref_end);
        int clipEnd =
            std::max(targetLength, alignerEndBegin + AlignRC(alignerEnd, queries[idx]).ref_begin);

        return BarcodeHit(idx, score, clipStart, clipEnd);
    }
}

static int Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM input, see --help" << std::endl;
        return EXIT_FAILURE;
    }

    const bool tryRC = options["tryRC"];
    const std::string mode = options["mode"];
    const int minScore = options["minScore"];
    const int minLength = options["minLength"];

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
        std::vector<Uhu::Threadpool::ThreadPool::TaskFuture<std::pair<BAM::BamRecord, std::string>>>
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
                    BarcodeHit bh =
                        SimdNeedleWunschAlignment(r.Sequence(), barcodes, StringToMode(mode));
                    bool aboveMinLength = (bh.ClipEnd - bh.ClipStart) >= minLength;
                    bool aboveMinScore = bh.Bq >= minScore;
                    if (aboveMinLength && aboveMinScore) {
                        r.Clip(BAM::ClipType::CLIP_TO_QUERY, bh.ClipStart, bh.ClipEnd);
                        r.Barcodes(std::make_pair(bh.Idx, bh.Idx));
                        r.BarcodeQuality(bh.Bq);
                        report = r.FullName() + "\t" + std::string(bh);
                        recordOut = std::move(r);
                        ++aboveThresholds;
                    } else if (!aboveMinLength && !aboveMinScore) {
                        ++belowBoth;
                    } else if (!aboveMinLength) {
                        ++belowMinLength;
                    } else if (!aboveMinScore) {
                        ++belowMinScore;
                    }
                    return std::make_pair(std::move(recordOut), report);
                },
                r));
        }

        std::ofstream report(prefix + ".demux.report");
        report << "ZMW\tIndex\tScore\tClipStar\tClipEnd" << std::endl;

        for (auto& item : v) {
            auto p = item.get();
            if (!p.second.empty()) {
                writer->Write(p.first);
                report << p.second << std::endl;
            }
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