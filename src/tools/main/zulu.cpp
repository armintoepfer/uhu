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
#include <iomanip>
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
namespace Zulu {
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

    PacBio::CLI::Interface i{"Zulu", "Barcode Validation Tool", "0.2.0"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    i.AddPositionalArguments({{"bam", "Source BAM", "FILE"}});

    // clang-format off
    i.AddOptions({
        OptionNames::Mapping,
        OptionNames::MinLength,
        OptionNames::Percentiles,
        OptionNames::Tailed,
        OptionNames::NumBC,
        OptionNames::MinPPV
    });
    // clang-format on

    return i;
}

// clang-format off
static std::string defaultMapping_ = "379=2kb_Gold.068.B.1501_3500.For,119=2kb_Gold.068.B.1501_3500.For,350=2kb_Gold.068.B.1501_3500.For,349=2kb_Gold.068.B.1501_3500.For,51=2kb_Gold.068.B.1501_3500.For,205=2kb_Gold.068.B.1501_3500.For,97=2kb_Gold.068.B.1501_3500.For,159=2kb_Gold.068.B.1501_3500.For,245=2kb_Gold.055.A.0001_2000.For,103=2kb_Gold.055.A.0001_2000.For,347=2kb_Gold.055.A.0001_2000.For,121=2kb_Gold.055.A.0001_2000.For,95=2kb_Gold.055.A.0001_2000.For,207=2kb_Gold.055.A.0001_2000.For,5=2kb_Gold.055.A.0001_2000.For,297=2kb_Gold.001.C.3001_5000.For,2=2kb_Gold.001.C.3001_5000.For,71=2kb_Gold.001.C.3001_5000.For,72=2kb_Gold.001.C.3001_5000.For,378=2kb_Gold.001.C.3001_5000.For,284=2kb_Gold.001.C.3001_5000.For,41=2kb_Gold.001.C.3001_5000.For,1=2kb_Gold.001.C.3001_5000.For,177=2kb_Gold.002.A.0001_2000.For,122=2kb_Gold.002.A.0001_2000.For,9=2kb_Gold.002.A.0001_2000.For,281=2kb_Gold.002.A.0001_2000.For,20=2kb_Gold.002.A.0001_2000.For,253=2kb_Gold.002.A.0001_2000.For,155=2kb_Gold.002.A.0001_2000.For,360=2kb_Gold.002.A.0001_2000.For,63=2kb_Gold.002.B.1501_3500.For,70=2kb_Gold.002.B.1501_3500.For,271=2kb_Gold.002.B.1501_3500.For,326=2kb_Gold.002.B.1501_3500.For,310=2kb_Gold.002.B.1501_3500.For,161=2kb_Gold.002.B.1501_3500.For,293=2kb_Gold.002.B.1501_3500.For,212=2kb_Gold.002.B.1501_3500.For,373=2kb_Gold.002.C.3001_5000.For,52=2kb_Gold.002.C.3001_5000.For,346=2kb_Gold.002.C.3001_5000.For,214=2kb_Gold.002.C.3001_5000.For,201=2kb_Gold.002.C.3001_5000.For,12=2kb_Gold.002.C.3001_5000.For,8=2kb_Gold.002.C.3001_5000.For,90=2kb_Gold.002.C.3001_5000.For,302=2kb_Gold.003.A.0001_2000.For,225=2kb_Gold.003.A.0001_2000.For,303=2kb_Gold.003.A.0001_2000.For,185=2kb_Gold.003.A.0001_2000.For,140=2kb_Gold.003.A.0001_2000.For,294=2kb_Gold.003.A.0001_2000.For,114=2kb_Gold.003.A.0001_2000.For,239=2kb_Gold.003.A.0001_2000.For,199=2kb_Gold.003.B.1501_3500.For,328=2kb_Gold.003.B.1501_3500.For,311=2kb_Gold.003.B.1501_3500.For,380=2kb_Gold.003.B.1501_3500.For,15=2kb_Gold.003.B.1501_3500.For,105=2kb_Gold.003.B.1501_3500.For,46=2kb_Gold.003.B.1501_3500.For,92=2kb_Gold.003.B.1501_3500.For,62=2kb_Gold.003.C.3001_5000.For,333=2kb_Gold.003.C.3001_5000.For,22=2kb_Gold.003.C.3001_5000.For,255=2kb_Gold.003.C.3001_5000.For,123=2kb_Gold.003.C.3001_5000.For,163=2kb_Gold.003.C.3001_5000.For,258=2kb_Gold.003.C.3001_5000.For,96=2kb_Gold.003.C.3001_5000.For,260=2kb_Gold.004.A.0001_2000.For,115=2kb_Gold.004.A.0001_2000.For,299=2kb_Gold.004.A.0001_2000.For,25=2kb_Gold.004.A.0001_2000.For,314=2kb_Gold.004.A.0001_2000.For,203=2kb_Gold.004.A.0001_2000.For,158=2kb_Gold.004.A.0001_2000.For,40=2kb_Gold.004.A.0001_2000.For,132=2kb_Gold.004.B.1501_3500.For,168=2kb_Gold.004.B.1501_3500.For,175=2kb_Gold.004.B.1501_3500.For,300=2kb_Gold.004.B.1501_3500.For,179=2kb_Gold.004.B.1501_3500.For,60=2kb_Gold.004.B.1501_3500.For,45=2kb_Gold.004.B.1501_3500.For,220=2kb_Gold.004.B.1501_3500.For,54=2kb_Gold.004.C.3001_5000.For,370=2kb_Gold.004.C.3001_5000.For,127=2kb_Gold.004.C.3001_5000.For,305=2kb_Gold.004.C.3001_5000.For,137=2kb_Gold.004.C.3001_5000.For,355=2kb_Gold.004.C.3001_5000.For,234=2kb_Gold.004.C.3001_5000.For,43=2kb_Gold.004.C.3001_5000.For,150=2kb_Gold.005.A.0001_2000.For,341=2kb_Gold.005.A.0001_2000.For,312=2kb_Gold.005.A.0001_2000.For,200=2kb_Gold.005.A.0001_2000.For,50=2kb_Gold.005.A.0001_2000.For,206=2kb_Gold.005.A.0001_2000.For,156=2kb_Gold.005.A.0001_2000.For,21=2kb_Gold.005.A.0001_2000.For,130=2kb_Gold.005.B.1501_3500.For,182=2kb_Gold.005.B.1501_3500.For,216=2kb_Gold.005.B.1501_3500.For,268=2kb_Gold.005.B.1501_3500.For,222=2kb_Gold.005.B.1501_3500.For,354=2kb_Gold.005.B.1501_3500.For,10=2kb_Gold.005.B.1501_3500.For,30=2kb_Gold.005.B.1501_3500.For,236=2kb_Gold.005.C.3001_5000.For,316=2kb_Gold.005.C.3001_5000.For,218=2kb_Gold.005.C.3001_5000.For,6=2kb_Gold.005.C.3001_5000.For,330=2kb_Gold.005.C.3001_5000.For,88=2kb_Gold.005.C.3001_5000.For,336=2kb_Gold.005.C.3001_5000.For,94=2kb_Gold.005.C.3001_5000.For,93=2kb_Gold.006.A.0001_2000.For,377=2kb_Gold.006.A.0001_2000.For,215=2kb_Gold.006.A.0001_2000.For,106=2kb_Gold.006.A.0001_2000.For,323=2kb_Gold.006.A.0001_2000.For,375=2kb_Gold.006.A.0001_2000.For,231=2kb_Gold.006.A.0001_2000.For,35=2kb_Gold.006.A.0001_2000.For,331=2kb_Gold.006.B.1501_3500.For,143=2kb_Gold.006.B.1501_3500.For,169=2kb_Gold.006.B.1501_3500.For,285=2kb_Gold.006.B.1501_3500.For,198=2kb_Gold.006.B.1501_3500.For,69=2kb_Gold.006.B.1501_3500.For,28=2kb_Gold.006.B.1501_3500.For,102=2kb_Gold.006.B.1501_3500.For,280=2kb_Gold.006.C.3001_5000.For,348=2kb_Gold.006.C.3001_5000.For,306=2kb_Gold.006.C.3001_5000.For,295=2kb_Gold.006.C.3001_5000.For,53=2kb_Gold.006.C.3001_5000.For,99=2kb_Gold.006.C.3001_5000.For,221=2kb_Gold.006.C.3001_5000.For,345=2kb_Gold.006.C.3001_5000.For,286=2kb_Gold.007.A.0001_2000.For,219=2kb_Gold.007.A.0001_2000.For,320=2kb_Gold.007.A.0001_2000.For,190=2kb_Gold.007.A.0001_2000.For,75=2kb_Gold.007.A.0001_2000.For,186=2kb_Gold.007.A.0001_2000.For,224=2kb_Gold.007.A.0001_2000.For,153=2kb_Gold.007.A.0001_2000.For,149=2kb_Gold.007.B.1501_3500.For,229=2kb_Gold.007.B.1501_3500.For,116=2kb_Gold.007.B.1501_3500.For,194=2kb_Gold.007.B.1501_3500.For,309=2kb_Gold.007.B.1501_3500.For,76=2kb_Gold.007.B.1501_3500.For,107=2kb_Gold.007.B.1501_3500.For,26=2kb_Gold.007.B.1501_3500.For,363=2kb_Gold.007.C.3001_5000.For,256=2kb_Gold.007.C.3001_5000.For,352=2kb_Gold.007.C.3001_5000.For,87=2kb_Gold.007.C.3001_5000.For,329=2kb_Gold.007.C.3001_5000.For,217=2kb_Gold.007.C.3001_5000.For,213=2kb_Gold.007.C.3001_5000.For,73=2kb_Gold.007.C.3001_5000.For,172=2kb_Gold.008.A.0001_2000.For,342=2kb_Gold.008.A.0001_2000.For,84=2kb_Gold.008.A.0001_2000.For,261=2kb_Gold.008.A.0001_2000.For,151=2kb_Gold.008.A.0001_2000.For,176=2kb_Gold.008.A.0001_2000.For,364=2kb_Gold.008.A.0001_2000.For,246=2kb_Gold.008.A.0001_2000.For,283=2kb_Gold.008.B.1501_3500.For,269=2kb_Gold.008.B.1501_3500.For,237=2kb_Gold.008.B.1501_3500.For,301=2kb_Gold.008.B.1501_3500.For,282=2kb_Gold.008.B.1501_3500.For,324=2kb_Gold.008.B.1501_3500.For,125=2kb_Gold.008.B.1501_3500.For,59=2kb_Gold.008.B.1501_3500.For,232=2kb_Gold.008.C.3001_5000.For,353=2kb_Gold.008.C.3001_5000.For,47=2kb_Gold.008.C.3001_5000.For,335=2kb_Gold.008.C.3001_5000.For,33=2kb_Gold.008.C.3001_5000.For,29=2kb_Gold.008.C.3001_5000.For,56=2kb_Gold.008.C.3001_5000.For,37=2kb_Gold.008.C.3001_5000.For,4=2kb_Gold.009.B.1501_3500.For,242=2kb_Gold.009.B.1501_3500.For,164=2kb_Gold.009.B.1501_3500.For,251=2kb_Gold.009.B.1501_3500.For,57=2kb_Gold.009.B.1501_3500.For,238=2kb_Gold.009.B.1501_3500.For,191=2kb_Gold.009.B.1501_3500.For,113=2kb_Gold.009.B.1501_3500.For,384=2kb_Gold.009.C.3001_5000.For,291=2kb_Gold.009.C.3001_5000.For,304=2kb_Gold.009.C.3001_5000.For,154=2kb_Gold.009.C.3001_5000.For,79=2kb_Gold.009.C.3001_5000.For,296=2kb_Gold.009.C.3001_5000.For,288=2kb_Gold.009.C.3001_5000.For,49=2kb_Gold.009.C.3001_5000.For,371=2kb_Gold.011.A.0001_2000.For,279=2kb_Gold.011.A.0001_2000.For,313=2kb_Gold.011.A.0001_2000.For,368=2kb_Gold.011.A.0001_2000.For,78=2kb_Gold.011.A.0001_2000.For,148=2kb_Gold.011.A.0001_2000.For,170=2kb_Gold.011.A.0001_2000.For,298=2kb_Gold.011.A.0001_2000.For,139=2kb_Gold.011.B.1501_3500.For,289=2kb_Gold.011.B.1501_3500.For,367=2kb_Gold.011.B.1501_3500.For,204=2kb_Gold.011.B.1501_3500.For,357=2kb_Gold.011.B.1501_3500.For,129=2kb_Gold.011.B.1501_3500.For,274=2kb_Gold.011.B.1501_3500.For,209=2kb_Gold.011.B.1501_3500.For,133=2kb_Gold.012.A.0001_2000.For,187=2kb_Gold.012.A.0001_2000.For,66=2kb_Gold.012.A.0001_2000.For,152=2kb_Gold.012.A.0001_2000.For,146=2kb_Gold.012.A.0001_2000.For,356=2kb_Gold.012.A.0001_2000.For,273=2kb_Gold.012.A.0001_2000.For,189=2kb_Gold.012.A.0001_2000.For,165=2kb_Gold.013.A.0001_2000.For,16=2kb_Gold.013.A.0001_2000.For,358=2kb_Gold.013.A.0001_2000.For,262=2kb_Gold.013.A.0001_2000.For,267=2kb_Gold.013.A.0001_2000.For,166=2kb_Gold.013.A.0001_2000.For,257=2kb_Gold.013.A.0001_2000.For,7=2kb_Gold.013.A.0001_2000.For,68=2kb_Gold.013.B.1501_3500.For,265=2kb_Gold.013.B.1501_3500.For,240=2kb_Gold.013.B.1501_3500.For,91=2kb_Gold.013.B.1501_3500.For,83=2kb_Gold.013.B.1501_3500.For,383=2kb_Gold.013.B.1501_3500.For,89=2kb_Gold.013.B.1501_3500.For,58=2kb_Gold.013.B.1501_3500.For,80=2kb_Gold.013.C.3001_5000.For,366=2kb_Gold.013.C.3001_5000.For,202=2kb_Gold.013.C.3001_5000.For,351=2kb_Gold.013.C.3001_5000.For,42=2kb_Gold.013.C.3001_5000.For,111=2kb_Gold.013.C.3001_5000.For,77=2kb_Gold.013.C.3001_5000.For,292=2kb_Gold.013.C.3001_5000.For,337=2kb_Gold.014.A.0001_2000.For,372=2kb_Gold.014.A.0001_2000.For,17=2kb_Gold.014.A.0001_2000.For,18=2kb_Gold.014.A.0001_2000.For,197=2kb_Gold.014.A.0001_2000.For,278=2kb_Gold.014.A.0001_2000.For,81=2kb_Gold.014.A.0001_2000.For,39=2kb_Gold.014.A.0001_2000.For,228=2kb_Gold.015.A.0001_2000.For,319=2kb_Gold.015.A.0001_2000.For,277=2kb_Gold.015.A.0001_2000.For,85=2kb_Gold.015.A.0001_2000.For,74=2kb_Gold.015.A.0001_2000.For,131=2kb_Gold.015.A.0001_2000.For,248=2kb_Gold.015.A.0001_2000.For,241=2kb_Gold.015.A.0001_2000.For,100=2kb_Gold.015.C.3001_5000.For,365=2kb_Gold.015.C.3001_5000.For,361=2kb_Gold.015.C.3001_5000.For,108=2kb_Gold.015.C.3001_5000.For,264=2kb_Gold.015.C.3001_5000.For,321=2kb_Gold.015.C.3001_5000.For,259=2kb_Gold.015.C.3001_5000.For,86=2kb_Gold.015.C.3001_5000.For,376=2kb_Gold.016.A.0001_2000.For,250=2kb_Gold.016.A.0001_2000.For,263=2kb_Gold.016.A.0001_2000.For,162=2kb_Gold.016.A.0001_2000.For,266=2kb_Gold.016.A.0001_2000.For,24=2kb_Gold.016.A.0001_2000.For,244=2kb_Gold.016.A.0001_2000.For,183=2kb_Gold.016.A.0001_2000.For,233=2kb_Gold.018.B.1501_3500.For,136=2kb_Gold.018.B.1501_3500.For,362=2kb_Gold.018.B.1501_3500.For,226=2kb_Gold.018.B.1501_3500.For,112=2kb_Gold.018.B.1501_3500.For,38=2kb_Gold.018.B.1501_3500.For,272=2kb_Gold.018.B.1501_3500.For,193=2kb_Gold.018.B.1501_3500.For,223=2kb_Gold.019.C.3001_5000.For,227=2kb_Gold.019.C.3001_5000.For,374=2kb_Gold.019.C.3001_5000.For,171=2kb_Gold.019.C.3001_5000.For,322=2kb_Gold.019.C.3001_5000.For,276=2kb_Gold.019.C.3001_5000.For,120=2kb_Gold.019.C.3001_5000.For,178=2kb_Gold.019.C.3001_5000.For,109=2kb_Gold.021.A.0001_2000.For,101=2kb_Gold.021.A.0001_2000.For,196=2kb_Gold.021.A.0001_2000.For,192=2kb_Gold.021.A.0001_2000.For,167=2kb_Gold.021.A.0001_2000.For,359=2kb_Gold.021.A.0001_2000.For,135=2kb_Gold.021.A.0001_2000.For,180=2kb_Gold.021.A.0001_2000.For,290=2kb_Gold.022.A.0001_2000.For,340=2kb_Gold.022.A.0001_2000.For,252=2kb_Gold.022.A.0001_2000.For,332=2kb_Gold.022.A.0001_2000.For,243=2kb_Gold.022.A.0001_2000.For,325=2kb_Gold.022.A.0001_2000.For,307=2kb_Gold.022.A.0001_2000.For,36=2kb_Gold.022.A.0001_2000.For,173=2kb_Gold.067.B.1501_3500.For,124=2kb_Gold.067.B.1501_3500.For,235=2kb_Gold.067.B.1501_3500.For,381=2kb_Gold.067.B.1501_3500.For,208=2kb_Gold.067.B.1501_3500.For,188=2kb_Gold.067.B.1501_3500.For,65=2kb_Gold.067.B.1501_3500.For,134=2kb_Gold.067.B.1501_3500.For,32=2kb_Gold.023.B.1501_3500.For,145=2kb_Gold.023.B.1501_3500.For,160=2kb_Gold.023.B.1501_3500.For,317=2kb_Gold.023.B.1501_3500.For,343=2kb_Gold.023.B.1501_3500.For,14=2kb_Gold.023.B.1501_3500.For,3=2kb_Gold.023.B.1501_3500.For,339=2kb_Gold.023.B.1501_3500.For,327=2kb_Gold.024.A.0001_2000.For,55=2kb_Gold.024.A.0001_2000.For,144=2kb_Gold.024.A.0001_2000.For,211=2kb_Gold.024.A.0001_2000.For,249=2kb_Gold.024.A.0001_2000.For,318=2kb_Gold.024.A.0001_2000.For,210=2kb_Gold.024.A.0001_2000.For,174=2kb_Gold.024.A.0001_2000.For,254=2kb_Gold.063.A.0001_2000.For,118=2kb_Gold.063.A.0001_2000.For,181=2kb_Gold.063.A.0001_2000.For,138=2kb_Gold.063.A.0001_2000.For,308=2kb_Gold.063.A.0001_2000.For,64=2kb_Gold.063.A.0001_2000.For,195=2kb_Gold.063.A.0001_2000.For,34=2kb_Gold.063.A.0001_2000.For,31=2kb_Gold.051.C.3001_5000.For,344=2kb_Gold.051.C.3001_5000.For,338=2kb_Gold.051.C.3001_5000.For,13=2kb_Gold.051.C.3001_5000.For,147=2kb_Gold.051.C.3001_5000.For,157=2kb_Gold.051.C.3001_5000.For,11=2kb_Gold.051.C.3001_5000.For,27=2kb_Gold.051.C.3001_5000.For,184=2kb_Gold.052.A.0001_2000.For,67=2kb_Gold.052.A.0001_2000.For,334=2kb_Gold.052.A.0001_2000.For,110=2kb_Gold.052.A.0001_2000.For,369=2kb_Gold.052.A.0001_2000.For,82=2kb_Gold.052.A.0001_2000.For,126=2kb_Gold.052.A.0001_2000.For,287=2kb_Gold.052.A.0001_2000.For,270=2kb_Gold.053.B.1501_3500.For,141=2kb_Gold.053.B.1501_3500.For,19=2kb_Gold.053.B.1501_3500.For,230=2kb_Gold.053.B.1501_3500.For,128=2kb_Gold.053.B.1501_3500.For,247=2kb_Gold.053.B.1501_3500.For,44=2kb_Gold.053.B.1501_3500.For,117=2kb_Gold.053.B.1501_3500.For,104=2kb_Gold.059.A.0001_2000.For,48=2kb_Gold.059.A.0001_2000.For,61=2kb_Gold.059.A.0001_2000.For,142=2kb_Gold.059.A.0001_2000.For,315=2kb_Gold.059.A.0001_2000.For,382=2kb_Gold.059.A.0001_2000.For,98=2kb_Gold.059.A.0001_2000.For,23=2kb_Gold.059.A.0001_2000.For";
// clang-format on

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
    const bool computeMinBQ = minPPV != 0;

    if (mapping.empty()) {
        mapping = defaultMapping_;
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

    std::ofstream report("report.uhu");
    report << "ReadName,HoleNumber,RefName,RefStart,RefEnd,RefLength,MapQuality,MappedID,"
              "BarcodedID,BarcodeFwd,BarcodeRev,BarcodeQuality"
           << std::endl;

    std::vector<std::pair<int, bool>> lengthsMatch;
    std::vector<std::pair<int, bool>> bqsMatch;
    std::map<int, std::vector<bool>> barcodeHits;
    std::map<std::string, int> zmwSubreads;
    std::map<std::string, int> zmwSubreadsMeasured;
    // a datastructure for each subread as a (barcodeCall, barcodeQuality) pair:
    //   (holeNumber) -> (mappedRefId) -> [subread<bc, bq>]
    std::map<std::string, std::map<std::string, std::vector<std::pair<int, int>>>> readsByZmw;
    for (const auto& r : *query) {
        if (r.Impl().IsSupplementaryAlignment()) continue;
        if (!r.Impl().IsPrimaryAlignment()) continue;
        if (!r.HasBarcodes() || !r.HasBarcodeQuality()) continue;

        const auto zmwNum = r.MovieName() + "/" + std::to_string(r.HoleNumber());
        const int length = r.ReferenceEnd() - r.ReferenceStart();
        ++zmwSubreads[zmwNum];
        if (length < minLength) {
            continue;
        }

        const std::string refName = r.ReferenceName();

        const int idx = tailed ? std::floor(r.BarcodeForward() / 2.0) + 1 : r.BarcodeForward() + 1;
        if (barcodeMapping.find(idx) == barcodeMapping.cend()) continue;
        const std::string bcRef = barcodeMapping.at(idx);

        const bool positive = refName == bcRef;
        barcodeHits[idx].emplace_back(positive);

        ++zmwSubreadsMeasured[zmwNum];

        report << r.FullName() << ',' << r.HoleNumber() << ',' << r.ReferenceName() << ','
               << r.ReferenceStart() << ',' << r.ReferenceEnd() << ',' << length << ','
               << (int)r.MapQuality() << ',' << refName << ',' << bcRef << ',' << r.BarcodeForward()
               << ',' << r.BarcodeReverse() << ',' << (int)r.BarcodeQuality() << std::endl;

        if (computeMinBQ) bqsMatch.emplace_back(std::make_pair(r.BarcodeQuality(), positive));
        if (nPercentiles > 1) lengthsMatch.emplace_back(std::make_pair(length, positive));

        readsByZmw[zmwNum][refName].emplace_back(std::make_pair(idx, (int)r.BarcodeQuality()));
    }

    double ppvSum = 0.0;
    double ppvCounter = 0;
    int missingBC = 0;
    std::ofstream barcodePPvStream("barcode_ppv.uhu");
    barcodePPvStream << "BC COUNTS PPV" << std::endl;
    for (const auto& bc_hits : barcodeHits) {
        if (std::any_of(bc_hits.second.cbegin(), bc_hits.second.cend(), [](bool x) { return x; })) {
            double bcPpv = 1.0 *
                           std::accumulate(bc_hits.second.cbegin(), bc_hits.second.cend(), 0) /
                           bc_hits.second.size();
            barcodePPvStream << bc_hits.first << " " << bc_hits.second.size() << " " << bcPpv
                             << std::endl;
            ppvSum += bcPpv;
            ++ppvCounter;
        } else {
            ++missingBC;
        }
    }

    // a datastructure for each bc as the mode of the mappedRefId to subread:
    //   (mappedRefId) -> [subread<bc, bq>]
    std::map<std::string, std::vector<std::pair<int, int>>> countsByRefName;
    // (holeNumber) -> (mappedRefId) -> [subread<bc, bq>]
    double byZmwAgreement = 0;
    for (const auto& zmwNum_refName_reads : readsByZmw) {
        std::string maxRefName;
        int maxReads = 0;

        // (mappedRefId) -> [subread<bc, bq>]
        for (const auto& refName_reads : zmwNum_refName_reads.second) {
            const int nReads = refName_reads.second.size();
            if (nReads > maxReads) {
                maxRefName = refName_reads.first;
                maxReads = nReads;
            }
        }

        // everybody from this holeNumber goes into the mode barcode
        int p = 0, n = 0;
        for (const auto& refName_reads : zmwNum_refName_reads.second) {
            for (const auto& r : refName_reads.second) {
                countsByRefName[maxRefName].emplace_back(r);
                if (barcodeMapping.at(r.first) == maxRefName) ++p;
                ++n;
            }
        }
        byZmwAgreement += (1.0 * p / n);
    }
    byZmwAgreement /= readsByZmw.size();

    double byBcPPV = 0;
    int truePositive2 = 0;
    for (const auto& refName_subreads : countsByRefName) {
        int p = 0, n = 0;

        for (const auto& r : refName_subreads.second) {
            if (barcodeMapping.at(r.first) == refName_subreads.first) {
                ++truePositive2;
                ++p;
            }
            ++n;
        }
        byBcPPV += (1.0 * p / n);
    }
    byBcPPV /= countsByRefName.size();

    const int nSubreads =
        std::accumulate(zmwSubreads.cbegin(), zmwSubreads.cend(), 0,
                        [](const int sum, const auto& next) { return sum + next.second; });
    const int nMeasured =
        std::accumulate(zmwSubreadsMeasured.cbegin(), zmwSubreadsMeasured.cend(), 0,
                        [](const int sum, const auto& next) { return sum + next.second; });
    // clang-format off
    std::cerr << "#Subreads input        : " << nSubreads << std::endl
              << "#Subreads BC & >" << std::setw(7) << std::left << std::to_string(minLength) + "bp" << ": "
              << nMeasured << std::endl
              << std::endl
              << "#ZMWs input            : " << zmwSubreads.size() << std::endl
              << "#ZMWs BC & >" << std::setw(11) << std::left << std::to_string(minLength) + "bp"
              << ": " << zmwSubreadsMeasured.size() << std::endl
              << std::endl;
    if (numBC > 0) std::cerr << "Barcode FN rate        : " << 1.0 * missingBC / numBC << std::endl;
    std::cerr << "PPV                    : " << std::round(1000.0 * ppvSum / ppvCounter) / 1000.0 << std::endl
              << std::endl
              << "%Mode/zmw              : " << std::round(1000.0 * byZmwAgreement) / 1000.0 << std::endl
              << "PPV/bc                 : " << std::round(1000.0 * byBcPPV) / 1000.0 << std::endl
              << "PPV/sr                 : " << std::round(1000.0 * (1.0 * truePositive2 / nMeasured)) / 1000.0 << std::endl
              << "#Refs                  : " << std::round(1000.0 * countsByRefName.size()) / 1000.0 << std::endl;
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
    return PacBio::CLI::Run(argc, argv, PacBio::Zulu::CreateCLI(), &PacBio::Zulu::Runner);
}