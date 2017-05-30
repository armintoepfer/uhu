// Copyright (c) 2016-2017, Pacific Biosciences of California, Inc.
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
#include <stdexcept>
#include <string>
#include <vector>

#include <ssw_cpp.h>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamReader.h>
#include <pbbam/Cigar.h>
#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>

namespace PacBio {
namespace Cleric {

struct Barcode
{
    Barcode(std::string name, std::string bases) : Name(name), Bases(bases) {}
    std::string Name;
    std::string Bases;
};

static PacBio::CLI::Interface CreateCLI()
{
    using Option = PacBio::CLI::Option;

    PacBio::CLI::Interface i{"demux_ccs", "Demultiplex Barcoded CCS Data", ""};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    i.AddPositionalArguments({{"bam", "Source BAM", "FILE"}});

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

void SimdNeedleWunschAlignment(const std::string& target, const std::vector<Barcode>& queries)
{
    int barcodeLength = queries.front().Bases.size() * 1.5;
    int targetLength = target.size();

    StripedSmithWaterman::Aligner alignerBegin;
    alignerBegin.SetReferenceSequence(target.c_str(), std::min(targetLength, barcodeLength));

    StripedSmithWaterman::Aligner alignerEnd;
    auto alignerEndBegin = std::max(targetLength - barcodeLength, 0);
    alignerEnd.SetReferenceSequence(target.c_str() + alignerEndBegin,
                                    targetLength - alignerEndBegin);
    StripedSmithWaterman::Filter filter;

    auto AlignTo = [&queries, &filter, &barcodeLength](StripedSmithWaterman::Aligner& aligner) {
        std::vector<int> scores;
        scores.resize(queries.size());
        std::vector<int> scoresRev;
        scoresRev.resize(queries.size());

        for (size_t i = 0; i < queries.size(); ++i) {
            StripedSmithWaterman::Alignment alignment;
            aligner.Align(queries[i].Bases.c_str(), filter, &alignment);
            scores[i] = alignment.sw_score;

            StripedSmithWaterman::Alignment alignmentRev;
            auto revComp = ReverseComplement(queries[i].Bases);
            aligner.Align(revComp.c_str(), filter, &alignmentRev);
            scoresRev[i] = alignmentRev.sw_score;
        }

        auto GetBestIndex = [](std::vector<int>& v) {
            // initialize original index locations
            std::vector<size_t> idx(v.size());
            std::iota(idx.begin(), idx.end(), 0);

            // sort indexes based on comparing values in v
            std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });

            std::cerr << idx.front() << " => " << v.at(idx.front()) << std::endl;
        };
        GetBestIndex(scores);
        GetBestIndex(scoresRev);
    };

    AlignTo(alignerBegin);
    AlignTo(alignerEnd);

    std::cerr << barcodeLength << std::endl;
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

    for (const auto& datasetPath : datasetPaths) {
        auto query = BamQuery(datasetPath);
        for (const auto& r : *query) {
            SimdNeedleWunschAlignment(r.Sequence(), barcodes);
        }
    }

    return EXIT_SUCCESS;
}
}
};

// Entry point
int main(int argc, char* argv[])
{
    return PacBio::CLI::Run(argc, argv, PacBio::Cleric::CreateCLI(), &PacBio::Cleric::Runner);
}