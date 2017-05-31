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
    Barcode(std::string name, std::string bases) : Name(name), Bases(bases) {}
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

    operator std::string()
    {
        std::stringstream out;
        out << static_cast<int>(Idx);
        out << "\t";
        out << static_cast<int>(Bq);
        out << "\t";
        out << ClipStart;
        out << "\t";
        out << ClipEnd;
        return out.str();
    }

    friend std::ostream& operator<<(std::ostream& stream, const BarcodeHit& bh)
    {
        stream << bh.Idx << "\t" << static_cast<int>(bh.Bq) << "\t" << bh.ClipStart << "\t"
               << bh.ClipEnd;
        return stream;
    }
};

static PacBio::CLI::Interface CreateCLI()
{
    using Option = PacBio::CLI::Option;

    PacBio::CLI::Interface i{"demux_ccs", "Demultiplex Barcoded CCS Data", "0.0.1"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    i.AddPositionalArguments(
        {{"bam", "Source BAM", "BAM_FILE"}, {"fasta", "Barcode file", "FASTA_FILE"}});

    return i;
}

static void ParsePositionalArgs(const std::vector<std::string>& args,
                                std::vector<std::string>* datasetPaths,
                                std::vector<Barcode>* barcodes, std::string* outputFile)
{
    std::vector<std::string> fastaPaths;
    for (const auto& i : args) {
        const bool fileExist = PacBio::Utility::FileExists(i);
        if (!fileExist) {
            if (!outputFile->empty())
                throw std::runtime_error(
                    "Only one output file allowed. Following files do not exist: " + *outputFile +
                    " and " + i);
            *outputFile = i;
            continue;
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

BarcodeHit SimdNeedleWunschAlignment(const std::string& target, const std::vector<Barcode>& queries)
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

    auto GetBestIndex = [&barcodeLength](std::vector<int>& v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });

        int bq = std::round(100.0 * v.at(idx.front()) / (barcodeLength * 2.0));
        return std::make_pair(idx.front(), bq);
    };

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
    if (forwardScore > revScore) {
        score = forwardScore;
        idx = forwardIdx;
        auto alignmentTmp = AlignForward(alignerBegin, queries[idx]);
        clipStart = alignmentTmp.ref_end;

        alignmentTmp = AlignRC(alignerEnd, queries[idx]);
        clipEnd = alignerEndBegin + alignmentTmp.ref_begin;

#if 0
        std::cerr << "Best Smith-Waterman score:\t" << alignmentTmp.sw_score << std::endl
                  << "Next-best Smith-Waterman score:\t" << alignmentTmp.sw_score_next_best
                  << std::endl
                  << "Reference start:\t" << alignmentTmp.ref_begin << std::endl
                  << "Reference end:\t" << alignmentTmp.ref_end << std::endl
                  << "Query start:\t" << alignmentTmp.query_begin << std::endl
                  << "Query end:\t" << alignmentTmp.query_end << std::endl
                  << "Next-best reference end:\t" << alignmentTmp.ref_end_next_best << std::endl
                  << "Number of mismatches:\t" << alignmentTmp.mismatches << std::endl
                  << "Cigar: " << alignmentTmp.cigar_string << std::endl;
#endif
    } else {
        score = revScore;
        idx = revIdx;
        auto alignmentTmp = AlignRC(alignerBegin, queries[idx]);
        clipStart = alignmentTmp.ref_end;

        alignmentTmp = AlignForward(alignerEnd, queries[idx]);
        clipEnd = alignerEndBegin + alignmentTmp.ref_begin;
    }

    return BarcodeHit(idx, score, clipStart, clipEnd);
}

static int Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM input, see --help" << std::endl;
        return EXIT_FAILURE;
    }

    std::vector<std::string> datasetPaths;
    std::vector<Barcode> barcodes;
    std::string outputFile;
    ParsePositionalArgs(options.PositionalArguments(), &datasetPaths, &barcodes, &outputFile);

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

    std::unique_ptr<BAM::BamWriter> writer;
    int counter = 0;
    std::cout << "ZMW\tIndex\tScore\tClipStar\tClipEnd" << std::endl;
    std::map<int, BAM::BamRecord> map;
    for (const auto& datasetPath : datasetPaths) {
        auto query = BamQuery(datasetPath);
        std::vector<Uhu::Threadpool::ThreadPool::TaskFuture<std::pair<BAM::BamRecord, std::string>>>
            v;
        for (auto& r : *query) {
            if (!writer)
                writer.reset(new BAM::BamWriter("out-" + std::to_string(counter++) + ".bam",
                                                r.Header().DeepCopy()));
            v.push_back(Uhu::Threadpool::DefaultThreadPool::submitJob(
                [&barcodes](BAM::BamRecord r) {
                    BAM::BamRecord recordOut;
                    std::string report;
                    BarcodeHit bh = SimdNeedleWunschAlignment(r.Sequence(), barcodes);
                    if (bh.ClipEnd - bh.ClipStart >= 50) {
                        r.Clip(BAM::ClipType::CLIP_TO_QUERY, bh.ClipStart, bh.ClipEnd);
                        r.Barcodes(std::make_pair(bh.Idx, bh.Idx));
                        r.BarcodeQuality(bh.Bq);
                        // writer->Write(r);
                        report = r.FullName() + "\t" + std::string(bh);
                        recordOut = std::move(r);
                    }
                    return std::make_pair(std::move(recordOut), report);
                },
                r));
        }

        for (auto& item : v) {
            auto p = item.get();
            if (!p.second.empty()) {
                writer->Write(p.first);
                std::cout << p.second << std::endl;
            }
        }
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