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

#pragma once

#include <string>
#include <vector>

#include <pacbio/lima/Lima.h>

namespace PacBio {
namespace Lima {
struct LimaWorkflow
{
    static BarcodeHitPair Tag(const std::vector<BAM::BamRecord> records,
                              const std::vector<Barcode>& queries, const LimaSettings& settings,
                              const AlignParameters& alignParameters);

    static void Process(
        const LimaSettings& settings,
        const std::vector<std::pair<std::string, BAM::DataSet::TypeEnum>>& datasetPaths,
        const std::vector<Barcode>& barcodes,
        const std::pair<std::string, BAM::DataSet::TypeEnum>& barcodePath);

    static int Runner(const PacBio::CLI::Results& options);

    static void ParsePositionalArgs(
        const std::vector<std::string>& args,
        std::vector<std::pair<std::string, BAM::DataSet::TypeEnum>>* datasetPaths,
        std::vector<Barcode>* barcodes,
        std::pair<std::string, BAM::DataSet::TypeEnum>* barcodePath);
};
}
}  // ::PacBio::Lima