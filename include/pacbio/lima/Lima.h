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

#include <atomic>
#include <memory>
#include <string>
#include <vector>

namespace BAM {
namespace internal {
class IQuery;
}
class BamRecord;
}

namespace StripedSmithWaterman {
struct Alignment;
class Aligner;
}

namespace PacBio {
namespace Lima {
struct LimaSettings;

struct SequenceUtils
{
    static char Complement(char base);
    static std::string ReverseComplement(const std::string& input);
};

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
    Barcode(const std::string& name, const std::string& bases)
        : Name(name), Bases(bases), BasesRC(SequenceUtils::ReverseComplement(bases))
    {
    }
    std::string Name;
    std::string Bases;
    std::string BasesRC;
};

struct BarcodeHit
{
    BarcodeHit() {}
    BarcodeHit(size_t reserveSize)
    {
        Scores.reserve(reserveSize);
        Clips.reserve(reserveSize);
    }

    uint16_t Idx = 0;
    uint8_t Score = 0;
    double ScoreSum = 0;
    std::vector<int> Scores;
    std::vector<int> Clips;

    void Add(int score, int clip);
    void AddWithSumScore(int score, int clip);
    void Normalize(int denominator);
};

struct BarcodeHitPair
{
    BarcodeHitPair() = delete;
    BarcodeHitPair(const BarcodeHit& left, const BarcodeHit& right)
        : Left(left), Right(right), MeanScore((Left.Score + Right.Score) / 2)
    {
    }
    BarcodeHitPair(BarcodeHit&& left, BarcodeHit&& right)
        : Left(std::forward<BarcodeHit>(left))
        , Right(std::forward<BarcodeHit>(right))
        , MeanScore((Left.Score + Right.Score) / 2)
    {
    }

    const BarcodeHit Left;
    const BarcodeHit Right;
    const uint8_t MeanScore;

    operator std::string() const;
};

struct Summary
{
    std::atomic_int BelowMinLength{0};
    std::atomic_int BelowMinScore{0};
    std::atomic_int BelowBoth{0};
    std::atomic_int AboveThresholds{0};
    int SymmetricCounts = 0;
    int AsymmetricCounts = 0;
    std::atomic_int SubreadBelowMinLength{0};
    std::atomic_int SubreadAboveMinLength{0};

    operator std::string() const;
};

struct AdvancedFileUtils
{
    static std::string FilePrefixInfix(const std::string& path);
    static std::unique_ptr<BAM::internal::IQuery> BamQuery(const std::string& filePath);
};

struct AlignParameters
{
    AlignParameters(int32_t matchScore, int32_t mismatchPenalty, int32_t deletionPenalty,
                    int32_t insertionPenalty, int32_t branchPenalty)
        : MatchScore(matchScore)
        , MismatchPenalty(mismatchPenalty)
        , DeletionPenalty(deletionPenalty)
        , InsertionPenalty(insertionPenalty)
        , BranchPenalty(branchPenalty)
    {
    }

    int32_t MatchScore;
    int32_t MismatchPenalty;
    int32_t DeletionPenalty;
    int32_t InsertionPenalty;
    int32_t BranchPenalty;
};

struct AlignUtils
{
    // Fills out a supplied SW matrix.
    static void SWComputeMatrix(const char* const query, const int32_t M, const char* const read,
                                const int32_t N, const bool globalInQuery,
                                std::vector<int32_t>& matrix,
                                const AlignParameters& parameters) noexcept;

    // Traverse the last row of an SW matrix (i.e. representing
    // alignments terminating with the last base of the query
    // sequence) and return the max score and it's position
    static std::pair<int32_t, int32_t> SWLastRowMax(const std::vector<int32_t>& matrix,
                                                    const int32_t queryLength,
                                                    const int32_t readLength) noexcept;

    static std::pair<int32_t, int32_t> Align(const std::string& bcBases, const char* target,
                                             const int targetSize, std::vector<int32_t>& matrix,
                                             const AlignParameters& parameters) noexcept;
};
}
}

#include "pacbio/lima/internal/Lima.inl"