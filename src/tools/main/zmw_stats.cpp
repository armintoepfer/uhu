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

#include <exception>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <pbcopper/cli/CLI.h>
#include <pbcopper/utility/FileUtils.h>

#include <pbbam/BamReader.h>
#include <pbbam/DataSet.h>
#include <pbbam/EntireFileQuery.h>
#include <pbbam/FastaReader.h>
#include <pbbam/PbiFilter.h>
#include <pbbam/PbiFilterQuery.h>
#include <pbbam/PbiIndex.h>

namespace PacBio {
namespace Cleric {

static PacBio::CLI::Interface CreateCLI()
{
    using Option = PacBio::CLI::Option;

    PacBio::CLI::Interface i{"zmw_stats", "Extracts per ZMW stats", "0.0.1"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    i.AddPositionalArguments({{"bam", "Source BAM", "FILE"}});

    return i;
}

static int Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM input, see --help" << std::endl;
        return EXIT_FAILURE;
    }

    std::map<int, std::vector<uint8_t>> cxPerZmw;

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

    auto query = BamQuery(options.PositionalArguments().front());

    std::cerr << "zmw";
    for (int i = 0; i < 16; ++i)
        std::cerr << "\t" << i;
    std::cerr << "\tbcf\tbq\tseq_count\tseq_mean\tseq_median\tseq_sd" << std::endl;

    int curZmw = -1;
    std::array<int, 16> cxUniq;
    cxUniq.fill(0);
    std::vector<int> seqLengths;

    using namespace boost::accumulators;
    using StatsAcc =
        accumulator_set<int, features<tag::count, tag::mean, tag::median, tag::variance>>;
    StatsAcc acc;

    std::vector<int> subreadLengths;
    int bq = -1;
    int bcf = -1;

    for (const auto& r : *query) {
        if (curZmw == -1) {
            curZmw = r.HoleNumber();
        } else if (curZmw != r.HoleNumber()) {
            std::cerr << curZmw;
            for (const auto& cx : cxUniq)
                std::cerr << "\t" << cx;
            std::cerr << "\t" << bcf << "\t" << bq;
            std::cerr << "\t" << count(acc) << "\t" << mean(acc) << "\t" << median(acc) << "\t"
                      << std::sqrt(variance(acc));
            std::cerr << std::endl;

            std::ofstream out(std::to_string(curZmw) + ".subreads");
            for (const auto& s : subreadLengths)
                out << s << std::endl;

            cxUniq.fill(0);
            curZmw = r.HoleNumber();
            acc = StatsAcc();
            subreadLengths.clear();
        }
        ++cxUniq[static_cast<uint8_t>(r.LocalContextFlags())];
        acc(r.Sequence().size());
        subreadLengths.emplace_back(r.Sequence().size());
        if (r.HasBarcodeQuality()) {
            bcf = r.BarcodeForward();
            bq = static_cast<int>(r.BarcodeQuality());
        } else {
            bq = -1;
            bcf = -1;
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