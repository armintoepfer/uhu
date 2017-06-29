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
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

#include <pacbio/data/PlainOption.h>

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
namespace OptionNames {
using PlainOption = Data::PlainOption;
// clang-format off
const PlainOption Mapping{
    "mapping",
    {"m","mapping"},
    "Mapping",
    "Mappings of bcs to ref names, example: 10=02.A,11=23.C,143=21.A",
    CLI::Option::StringType()
};

const PlainOption MinLength{
    "min_length",
    {"l","min-length"},
    "MinLength",
    "Minimum reference span to score a read.",
    CLI::Option::IntType(0)
};

const PlainOption Percentiles{
    "percentiles",
    {"p", "percentiles"},
    "Percentiles",
    "Number of percentiles between [0, 100] to compute.",
    CLI::Option::IntType(0)
};

const PlainOption Tailed{
    "tailed",
    {"t", "tailed"},
    "Tailed",
    "Flag to analyze in tailed mode.",
    CLI::Option::BoolType()
};

const PlainOption ZMW{
    "zmw",
    {"z", "zmw"},
    "Zmw",
    "Flag to analyze only one subread per ZMW.",
    CLI::Option::BoolType()
};

const PlainOption NumBC{
    "num_bc",
    {"b","num-barcodes"},
    "NumBC",
    "Number of barcodes used; 0 means, don't compute FN rate.",
    CLI::Option::IntType(0)
};

const PlainOption MinPPV{
    "min_ppv",
    {"v","min-ppv"},
    "MinPPV",
    "Compute the minimal Barcode Score for a given PPV.",
    CLI::Option::FloatType(0)
};

// clang-format on
}  // namespace OptionNames

static PacBio::CLI::Interface CreateCLI()
{
    using Option = PacBio::CLI::Option;

    PacBio::CLI::Interface i{"zmw_to_ref", "Maps Record to Reference", "0.1.0"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    i.AddPositionalArguments({{"bam", "Source BAM", "FILE"}});

    // clang-format off
    i.AddOptions({
        OptionNames::Mapping,
        OptionNames::MinLength,
        OptionNames::Percentiles,
        OptionNames::Tailed,
        OptionNames::ZMW,
        OptionNames::NumBC,
        OptionNames::MinPPV
    });
    // clang-format on

    return i;
}

static int Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM input, see --help" << std::endl;
        return EXIT_FAILURE;
    }

    std::string mapping = options[OptionNames::Mapping];
    const int minLength = options[OptionNames::MinLength];
    const int nPercentiles = options[OptionNames::Percentiles];
    const int numBC = options[OptionNames::NumBC];
    const double minPPV = options[OptionNames::MinPPV];
    const bool tailed = options[OptionNames::Tailed];
    const bool zmwMode = options[OptionNames::ZMW];
    const bool computeMinBQ = minPPV != 0;

    if (mapping.empty()) {
        mapping =
            "379=068.B,119=068.B,350=068.B,349=068.B,51=068.B,205=068.B,97=068.B,159=068.B,275=055."
            "A,245=055.A,103=055.A,347=055.A,121=055.A,95=055.A,207=055.A,5=055.A,297=001.C,2=001."
            "C,71=001.C,72=001.C,378=001.C,284=001.C,41=001.C,1=001.C,177=002.A,122=002.A,9=002.A,"
            "281=002.A,20=002.A,253=002.A,155=002.A,360=002.A,63=002.B,70=002.B,271=002.B,326=002."
            "B,310=002.B,161=002.B,293=002.B,212=002.B,373=002.C,52=002.C,346=002.C,214=002.C,201="
            "002.C,12=002.C,8=002.C,90=002.C,302=003.A,225=003.A,303=003.A,185=003.A,140=003.A,294="
            "003.A,114=003.A,239=003.A,199=003.B,328=003.B,311=003.B,380=003.B,15=003.B,105=003.B,"
            "46=003.B,92=003.B,62=003.C,333=003.C,22=003.C,255=003.C,123=003.C,163=003.C,258=003.C,"
            "96=003.C,260=004.A,115=004.A,299=004.A,25=004.A,314=004.A,203=004.A,158=004.A,40=004."
            "A,132=004.B,168=004.B,175=004.B,300=004.B,179=004.B,60=004.B,45=004.B,220=004.B,54="
            "004.C,370=004.C,127=004.C,305=004.C,137=004.C,355=004.C,234=004.C,43=004.C,150=005.A,"
            "341=005.A,312=005.A,200=005.A,50=005.A,206=005.A,156=005.A,21=005.A,130=005.B,182=005."
            "B,216=005.B,268=005.B,222=005.B,354=005.B,10=005.B,30=005.B,236=005.C,316=005.C,218="
            "005.C,6=005.C,330=005.C,88=005.C,336=005.C,94=005.C,93=006.A,377=006.A,215=006.A,106="
            "006.A,323=006.A,375=006.A,231=006.A,35=006.A,331=006.B,143=006.B,169=006.B,285=006.B,"
            "198=006.B,69=006.B,28=006.B,102=006.B,280=006.C,348=006.C,306=006.C,295=006.C,53=006."
            "C,99=006.C,221=006.C,345=006.C,286=007.A,219=007.A,320=007.A,190=007.A,75=007.A,186="
            "007.A,224=007.A,153=007.A,149=007.B,229=007.B,116=007.B,194=007.B,309=007.B,76=007.B,"
            "107=007.B,26=007.B,363=007.C,256=007.C,352=007.C,87=007.C,329=007.C,217=007.C,213=007."
            "C,73=007.C,172=008.A,342=008.A,84=008.A,261=008.A,151=008.A,176=008.A,364=008.A,246="
            "008.A,283=008.B,269=008.B,237=008.B,301=008.B,282=008.B,324=008.B,125=008.B,59=008.B,"
            "232=008.C,353=008.C,47=008.C,335=008.C,33=008.C,29=008.C,56=008.C,37=008.C,4=009.B,"
            "242=009.B,164=009.B,251=009.B,57=009.B,238=009.B,191=009.B,113=009.B,384=009.C,291="
            "009.C,304=009.C,154=009.C,79=009.C,296=009.C,288=009.C,49=009.C,371=011.A,279=011.A,"
            "313=011.A,368=011.A,78=011.A,148=011.A,170=011.A,298=011.A,139=011.B,289=011.B,367="
            "011.B,204=011.B,357=011.B,129=011.B,274=011.B,209=011.B,133=012.A,187=012.A,66=012.A,"
            "152=012.A,146=012.A,356=012.A,273=012.A,189=012.A,165=013.A,16=013.A,358=013.A,262="
            "013.A,267=013.A,166=013.A,257=013.A,7=013.A,68=013.B,265=013.B,240=013.B,91=013.B,83="
            "013.B,383=013.B,89=013.B,58=013.B,80=013.C,366=013.C,202=013.C,351=013.C,42=013.C,111="
            "013.C,77=013.C,292=013.C,337=014.A,372=014.A,17=014.A,18=014.A,197=014.A,278=014.A,81="
            "014.A,39=014.A,228=015.A,319=015.A,277=015.A,85=015.A,74=015.A,131=015.A,248=015.A,"
            "241=015.A,100=015.C,365=015.C,361=015.C,108=015.C,264=015.C,321=015.C,259=015.C,86="
            "015.C,376=016.A,250=016.A,263=016.A,162=016.A,266=016.A,24=016.A,244=016.A,183=016.A,"
            "233=018.B,136=018.B,362=018.B,226=018.B,112=018.B,38=018.B,272=018.B,193=018.B,223="
            "019.C,227=019.C,374=019.C,171=019.C,322=019.C,276=019.C,120=019.C,178=019.C,109=021.A,"
            "101=021.A,196=021.A,192=021.A,167=021.A,359=021.A,135=021.A,180=021.A,290=022.A,340="
            "022.A,252=022.A,332=022.A,243=022.A,325=022.A,307=022.A,36=022.A,173=067.B,124=067.B,"
            "235=067.B,381=067.B,208=067.B,188=067.B,65=067.B,134=067.B,32=023.B,145=023.B,160=023."
            "B,317=023.B,343=023.B,14=023.B,3=023.B,339=023.B,327=024.A,55=024.A,144=024.A,211=024."
            "A,249=024.A,318=024.A,210=024.A,174=024.A,254=063.A,118=063.A,181=063.A,138=063.A,308="
            "063.A,64=063.A,195=063.A,34=063.A,31=051.C,344=051.C,338=051.C,13=051.C,147=051.C,157="
            "051.C,11=051.C,27=051.C,184=052.A,67=052.A,334=052.A,110=052.A,369=052.A,82=052.A,126="
            "052.A,287=052.A,270=053.B,141=053.B,19=053.B,230=053.B,128=053.B,247=053.B,44=053.B,"
            "117=053.B,104=059.A,48=059.A,61=059.A,142=059.A,315=059.A,382=059.A,98=059.A,23=059.A";
        // std::cerr << "ERROR: --mapping must be set" << std::endl;
        // return EXIT_FAILURE;
    }

    if (minLength < 0) {
        std::cerr << "ERROR: --min-length must be >= 0" << std::endl;
        return EXIT_FAILURE;
    }

    if (nPercentiles < 0) {
        std::cerr << "ERROR: --percentiles must be >= 0" << std::endl;
        return EXIT_FAILURE;
    }

    auto SplitMappingString = [&mapping]() {
        std::stringstream ss(mapping);
        std::string item;
        std::map<int, std::string> barcodeMapping;
        while (std::getline(ss, item, ',')) {
            std::stringstream ss2(item);
            std::string item2;
            if (!std::getline(ss2, item2, '='))
                throw std::runtime_error("Could not parse mapping string.");
            int providedBC = std::stoi(item2);
            if (!std::getline(ss2, item2, '='))
                throw std::runtime_error("Could not parse mapping string.");
            barcodeMapping[providedBC] = item2;
        }
        return barcodeMapping;
    };
    std::map<int, std::string> barcodeMapping = SplitMappingString();

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

    std::ofstream report("report");
    report << "ReadName,HoleNumber,RefName,RefStart,RefEnd,RefLength,MapQuality,MappedID,"
              "BarcodedID,BarcodeFwd,BarcodeRev,BarcodeQuality"
           << std::endl;

    std::vector<std::pair<int, bool>> lengthsMatch;
    std::vector<std::pair<int, bool>> bqsMatch;
    std::map<int, std::vector<bool>> barcodeHits;
    for (int i = 1; i <= numBC; ++i)
        barcodeHits[i] = std::vector<bool>();
    std::map<int, int> zmwSubreads;
    std::map<int, int> zmsSubreadsMeasured;
    for (const auto& r : *query) {
        if (r.Impl().IsSupplementaryAlignment()) continue;
        if (!r.Impl().IsPrimaryAlignment()) continue;
        if (!r.HasBarcodes() || !r.HasBarcodeQuality()) continue;

        const auto zmwNum = r.HoleNumber();
        const int length = r.ReferenceEnd() - r.ReferenceStart();
        ++zmwSubreads[zmwNum];
        if (length < minLength) {
            continue;
        }

        if ((zmwMode && zmsSubreadsMeasured.find(zmwNum) == zmsSubreadsMeasured.cend()) ||
            !zmwMode) {

            std::stringstream ss(r.ReferenceName());
            std::string item;
            std::vector<std::string> elems;
            int i = 0;
            std::string refName;
            while (std::getline(ss, item, '.')) {
                if (i == 1 || i == 2) refName += item;
                if (i == 1) refName += ".";
                ++i;
            }

            const int idx =
                tailed ? std::floor(r.BarcodeForward() / 2.0) + 1 : r.BarcodeForward() + 1;
            const std::string bcRef = barcodeMapping.at(idx);

            const bool positive = refName == bcRef;
            barcodeHits[idx].emplace_back(positive);

            ++zmsSubreadsMeasured[zmwNum];

            report << r.FullName() << ',' << zmwNum << ',' << r.ReferenceName() << ','
                   << r.ReferenceStart() << ',' << r.ReferenceEnd() << ',' << length << ','
                   << (int)r.MapQuality() << ',' << refName << ',' << bcRef << ','
                   << r.BarcodeForward() << ',' << r.BarcodeReverse() << ','
                   << (int)r.BarcodeQuality() << std::endl;

            if (computeMinBQ) bqsMatch.emplace_back(std::make_pair(r.BarcodeQuality(), positive));
            if (nPercentiles > 1) lengthsMatch.emplace_back(std::make_pair(length, positive));
        }
    }

    double ppvSum = 0.0;
    double ppvCounter = 0;
    int missingBC = 0;
    for (const auto& bc_hits : barcodeHits) {
        if (std::any_of(bc_hits.second.cbegin(), bc_hits.second.cend(), [](bool x) { return x; })) {
            ppvSum += 1.0 * std::accumulate(bc_hits.second.cbegin(), bc_hits.second.cend(), 0) /
                      bc_hits.second.size();
            ++ppvCounter;
        } else {
            ++missingBC;
        }
    }

    // clang-format off
    std::cerr << "#Subreads input        : "
              << std::accumulate(zmwSubreads.cbegin(), zmwSubreads.cend(), 0,
                                 [](const auto& sum, const auto& next) { return sum + next.second; })
              << std::endl
              << "#Subreads BC & >" << std::setw(7) << std::left << std::to_string(minLength) + "bp" << ": "
              << std::accumulate(zmsSubreadsMeasured.cbegin(), zmsSubreadsMeasured.cend(), 0,
                                 [](const auto& sum, const auto& next) { return sum + next.second; })
              << std::endl
              << std::endl
              << "#ZMWs input            : " << zmwSubreads.size() << std::endl
              << "#ZMWs BC & >" << std::setw(11) << std::left << std::to_string(minLength) + "bp"
              << ": " << zmsSubreadsMeasured.size() << std::endl
              << std::endl;
    if (numBC > 0) std::cerr << "Barcode FN rate        : " << 1.0 * missingBC / numBC << std::endl;
    std::cerr << "PPV                    : " << ppvSum / ppvCounter << std::endl;
    // clang-format on

    // Determine minimal BQ needed for given minimal PPV
    const int nBQs = bqsMatch.size();
    if (nBQs > 1 && computeMinBQ) {
        std::sort(bqsMatch.begin(), bqsMatch.end(),
                  [](const auto& a, const auto& b) { return a > b; });
        double ppv = 1;
        int i = 0;
        int pos = 0;
        int neg = 0;
        do {
            if (bqsMatch[i].second)
                ++pos;
            else
                ++neg;
            ppv = 1.0 * pos / (pos + neg);
            ++i;
        } while (ppv > minPPV && i < nBQs);
        std::cerr << "Min BQ for PPV " << std::setw(8) << std::left << minPPV << ": "
                  << bqsMatch[i].first << std::endl;
    }

    // Percentiles for length
    const int nObs = lengthsMatch.size();
    if (nObs > 1 && nPercentiles > 1) {
        std::cerr << std::endl;
        std::sort(lengthsMatch.begin(), lengthsMatch.end());
        // Second Variant, C = 0 from
        //  https://en.wikipedia.org/wiki/Percentile#The_Linear_Interpolation_Between_Closest_Ranks_method
        std::vector<int> breaks;
        for (int p = 1; p <= nPercentiles; ++p) {
            // no + 1 because we're 0-indexing here
            const double x = static_cast<double>(p) / nPercentiles * (nObs - 1);
            breaks.push_back(static_cast<int>(std::floor(x) + 1));
        }
        int b = 0;
        int i = 0;
        int pos = 0;
        int neg = 0;
        for (const auto obs : lengthsMatch) {
            if (breaks[b] == i) {
                const double perc = static_cast<double>(b + 1) / nPercentiles;
                const double x = perc * (nObs - 1);
                const double v = x - i;
                const double l =
                    lengthsMatch[i - 1].first + v * (obs.first - lengthsMatch[i - 1].first);
                std::cerr << "PPV(" << perc * 100 << ", " << l << ") : " << pos << " ("
                          << (1.0 * pos / (pos + neg)) << ")" << std::endl;
                ++b;
                pos = 0;
                neg = 0;
            }
            if (obs.second)
                ++pos;
            else
                ++neg;
            ++i;
        }
        std::cerr << "PPV(100, " << lengthsMatch.back().first << ") : " << pos << " ("
                  << (1.0 * pos / (pos + neg)) << ")" << std::endl;
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
