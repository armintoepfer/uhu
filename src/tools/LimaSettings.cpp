// Copyright (c) 2014-2017, Pacific Biosciences of California, Inc.
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

#include <pacbio/data/PlainOption.h>

#include <pacbio/lima/LimaSettings.h>

namespace PacBio {
namespace Lima {
namespace OptionNames {
using PlainOption = Data::PlainOption;
// clang-format off
const PlainOption KeepSymmetric{
    "KeepSymmetric",
    {"s","symmetric"},
    "KeepSymmetric",
    "Only keep symmetric barcodes in BAM output.",
    CLI::Option::BoolType()
};

const PlainOption TryRC{
    "tryRC",
    {"t","try-rc"},
    "TryRC",
    "Try barcodes also as reverse complements.",
    CLI::Option::BoolType()
};

const PlainOption WindowSizeMult{
    "windowSizeMult",
    {"w","window-size-mult"},
    "WindowSizeMult",
    "The candidate region size multiplier: barcode_length * multiplier.",
    CLI::Option::FloatType(1.2)
};

const PlainOption MinScore{
    "minScore",
    {"m","min-score"},
    "MinScore",
    "Minimum barcode score.",
    CLI::Option::IntType(51)
};

const PlainOption MinLength{
    "minLength",
    {"l","min-length"},
    "MinLength",
    "Minimum sequence length after clipping.",
    CLI::Option::IntType(50)
};

const PlainOption MatchScore{
    "matchScore",
    {"A", "match-score"},
    "MatchScore",
    "Score for a sequence match.",
    CLI::Option::IntType(4)
};

const PlainOption MismatchPenalty{
    "mismatchPenalty",
    {"B", "mismatch-penalty"},
    "MismatchPenalty",
    "Penalty for a mismatch.",
    CLI::Option::IntType(13)
};

const PlainOption GapOpenPenalty{
    "gapOpenPenalty",
    {"O", "gap-open-penalty"},
    "GapOpenPenalty",
    "Gap open penalties for deletions and insertions.",
    CLI::Option::IntType(7)
};

const PlainOption GapExtPenalty{
    "gapExtPenalty",
    {"E", "gap-ext-penalty"},
    "GapExtPenalty",
    "Gap extension penalties for deletions and insertions.",
    CLI::Option::IntType(7)
};

const PlainOption NoBam{
    "NoBam",
    {"no-bam"},
    "NoBam",
    "Do not generate BAM output.",
    CLI::Option::BoolType()
};

const PlainOption NoReports{
    "NoReports",
    {"no-reports"},
    "NoReports",
    "Do not generate reports.",
    CLI::Option::BoolType()
};

const PlainOption SplitBam{
    "SplitBam",
    {"split-bam"},
    "SplitBam",
    "Split BAM output by barcode pair.",
    CLI::Option::BoolType()
};

const PlainOption CCS{
    "CCS",
    {"ccs"},
    "CCS",
    "CCS mode, use optimal alignment options -A 4 -B 1 -O 3 -E 1.",
    CLI::Option::BoolType()
};
// clang-format on
}  // namespace OptionNames

LimaSettings::LimaSettings(const PacBio::CLI::Results& options)
    : CLI(options.InputCommandLine())
    , InputFiles(options.PositionalArguments())
    , WindowSizeMult(options[OptionNames::WindowSizeMult])
    , KeepSymmetric(options[OptionNames::KeepSymmetric])
    , MinScore(options[OptionNames::MinScore])
    , MinLength(options[OptionNames::MinLength])
    , NoBam(options[OptionNames::NoBam])
    , NoReports(options[OptionNames::NoReports])
    , SplitBam(options[OptionNames::SplitBam])
{
    if (SplitBam && NoBam)
        throw std::runtime_error("Options --split-bam and --no-bam are mutually exclusive!");

    if (options[OptionNames::SplitBam]) {
        MatchScore = 4;
        MismatchPenalty = 1;
        GapOpenPenalty = 3;
        GapExtPenalty = 1;
    }

    if (static_cast<int>(options[OptionNames::MatchScore]) !=
        static_cast<int>(OptionNames::MatchScore.defaultValue))
        MatchScore = options[OptionNames::MatchScore];
    if (static_cast<int>(options[OptionNames::MismatchPenalty]) !=
        static_cast<int>(OptionNames::MismatchPenalty.defaultValue))
        MismatchPenalty = options[OptionNames::MismatchPenalty];
    if (static_cast<int>(options[OptionNames::GapOpenPenalty]) !=
        static_cast<int>(OptionNames::GapOpenPenalty.defaultValue))
        GapOpenPenalty = options[OptionNames::GapOpenPenalty];
    if (static_cast<int>(options[OptionNames::GapExtPenalty]) !=
        static_cast<int>(OptionNames::GapExtPenalty.defaultValue))
        GapExtPenalty = options[OptionNames::GapExtPenalty];
}

PacBio::CLI::Interface LimaSettings::CreateCLI()
{
    using Option = PacBio::CLI::Option;
    using Task = PacBio::CLI::ToolContract::Task;

    PacBio::CLI::Interface i{"lima", "Lima, Demultiplex Barcoded PacBio Data and Clip Barcodes ",
                             "0.8.0"};

    i.AddHelpOption();     // use built-in help output
    i.AddVersionOption();  // use built-in version output

    // clang-format off
    i.AddPositionalArguments({
        {"bam", "Source BAM", "BAM_FILE"},
        {"fasta", "Barcode file", "FASTA_FILE"}
    });

    i.AddGroup("Tuning",
    {
        OptionNames::KeepSymmetric,
        OptionNames::WindowSizeMult,
        OptionNames::MinLength,
        OptionNames::MinScore
    });

    i.AddGroup("Aligner Configuration",
    {
        OptionNames::CCS,
        OptionNames::MatchScore,
        OptionNames::MismatchPenalty,
        OptionNames::GapOpenPenalty,
        OptionNames::GapExtPenalty
    });

    i.AddGroup("Output Restrictions",
    {
        OptionNames::NoBam,
        OptionNames::SplitBam,
        OptionNames::NoReports
    });
    // clang-format on

    return i;
}
}
}  // ::PacBio::Lima
