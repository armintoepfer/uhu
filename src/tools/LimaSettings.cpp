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

const PlainOption MaxScoredReads{
    "maxScoredReads",
    {"n","max-scored-reads"},
    "MaxScoredReads",
    "Only use up to N reads to find the barcode, 0 means use all.",
    CLI::Option::IntType(0)
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

const PlainOption DeletionPenalty{
    "deletionPenalty",
    {"D", "deletion-penalty"},
    "DeletionPenalty",
    "Deletions penalty.",
    CLI::Option::IntType(7)
};

const PlainOption InsertionPenalty{
    "insertionPenalty",
    {"I", "insertion-penalty"},
    "InsertionPenalty",
    "Insertion penalty.",
    CLI::Option::IntType(7)
};

const PlainOption BranchPenalty{
    "branchPenalty",
    {"X", "branch-penalty"},
    "BranchPenalty",
    "Branch penalty.",
    CLI::Option::IntType(4)
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
    "CCS mode, use optimal alignment options -A 4 -B 1 -D 3 -I 3 -X 4.",
    CLI::Option::BoolType()
};
const PlainOption NumThreads{
    "NumThreads",
    { "j", "numThreads" },
    "Number of Threads",
    "Number of threads to use, 0 means autodetection.",
    CLI::Option::IntType(0)
};
const PlainOption Chunks{
    "Chunks",
    { "c", "chunk-size" },
    "Size of Chunks",
    "Size of Chunks.",
    CLI::Option::IntType(10)
};
const PlainOption PerSubread{
    "PerSubread",
    { "p", "per-subread" },
    "Tag per subread",
    "Do not tag per ZMW, but per subread.",
    CLI::Option::BoolType()
};
const PlainOption MinPasses{
    "minPasses",
    { "u", "min-passes" },
    "Minimal Number Passes",
    "Minimal number of full passes.",
    CLI::Option::IntType(1)
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
    , MatchScore(options[OptionNames::MatchScore])
    , MismatchPenalty(options[OptionNames::MismatchPenalty])
    , DeletionPenalty(options[OptionNames::DeletionPenalty])
    , InsertionPenalty(options[OptionNames::InsertionPenalty])
    , BranchPenalty(options[OptionNames::BranchPenalty])
    , NoBam(options[OptionNames::NoBam])
    , NoReports(options[OptionNames::NoReports])
    , SplitBam(options[OptionNames::SplitBam])
    , MaxScoredReads(options[OptionNames::MaxScoredReads])
    , Chunks(options[OptionNames::Chunks])
    , PerSubread(options[OptionNames::PerSubread])
    , MinPasses(options[OptionNames::MinPasses])
{
    if (SplitBam && NoBam)
        throw std::runtime_error("Options --split-bam and --no-bam are mutually exclusive!");

    if (options[OptionNames::CCS]) {
        MatchScore = 4;
        MismatchPenalty = 1;
        DeletionPenalty = 3;
        InsertionPenalty = 3;
        BranchPenalty = 3;
    }

    if (static_cast<int>(options[OptionNames::MatchScore]) !=
        static_cast<int>(OptionNames::MatchScore.defaultValue))
        MatchScore = options[OptionNames::MatchScore];
    if (static_cast<int>(options[OptionNames::MismatchPenalty]) !=
        static_cast<int>(OptionNames::MismatchPenalty.defaultValue))
        MismatchPenalty = options[OptionNames::MismatchPenalty];
    if (static_cast<int>(options[OptionNames::DeletionPenalty]) !=
        static_cast<int>(OptionNames::DeletionPenalty.defaultValue))
        DeletionPenalty = options[OptionNames::DeletionPenalty];
    if (static_cast<int>(options[OptionNames::InsertionPenalty]) !=
        static_cast<int>(OptionNames::InsertionPenalty.defaultValue))
        InsertionPenalty = options[OptionNames::InsertionPenalty];
    if (static_cast<int>(options[OptionNames::BranchPenalty]) !=
        static_cast<int>(OptionNames::BranchPenalty.defaultValue))
        BranchPenalty = options[OptionNames::BranchPenalty];

    int requestedNThreads;
    if (options.IsFromRTC()) {
        requestedNThreads = options.NumProcessors();
    } else {
        requestedNThreads = options[OptionNames::NumThreads];
    }
    NumThreads = ThreadCount(requestedNThreads);
}

size_t LimaSettings::ThreadCount(int n)
{
    const int m = std::thread::hardware_concurrency();

    if (n < 1) return std::max(1, m + n);

    return std::min(m, n);
}

PacBio::CLI::Interface LimaSettings::CreateCLI()
{
    using Option = PacBio::CLI::Option;
    using Task = PacBio::CLI::ToolContract::Task;

    PacBio::CLI::Interface i{"lima", "Lima, Demultiplex Barcoded PacBio Data and Clip Barcodes ",
                             "0.13.0"};

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
        OptionNames::MinScore,
        OptionNames::MaxScoredReads,
        OptionNames::MinPasses,
        OptionNames::Chunks,
        OptionNames::PerSubread
    });

    i.AddGroup("Aligner Configuration",
    {
        OptionNames::CCS,
        OptionNames::MatchScore,
        OptionNames::MismatchPenalty,
        OptionNames::DeletionPenalty,
        OptionNames::InsertionPenalty,
        OptionNames::BranchPenalty
    });

    i.AddGroup("Output Restrictions",
    {
        OptionNames::NoBam,
        OptionNames::SplitBam,
        OptionNames::NoReports
    });
    i.AddOptions({
        OptionNames::NumThreads
    });
    // clang-format on

    return i;
}
}
}  // ::PacBio::Lima
