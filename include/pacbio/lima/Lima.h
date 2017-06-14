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
    Barcode(const std::string& name, const std::string& bases) : Name(name), Bases(bases) {}
    std::string Name;
    std::string Bases;
};

struct ScoreClip
{
    ScoreClip(size_t reserveSize)
    {
        Scores.reserve(reserveSize);
        Clips.reserve(reserveSize);
    }

    double ScoreSum = 0;
    std::vector<int> Scores;
    std::vector<int> Clips;

    void Add(int score, int clip);
};

struct BarcodeHit
{
    BarcodeHit() = default;
    BarcodeHit(int idx, int score, int clip) : Idx(idx), Score(score), Clip(clip) {}
    BarcodeHit(int idx, int score, std::vector<int>& clips) : Idx(idx), Score(score), Clips(clips)
    {
    }
    BarcodeHit(int idx, int score, std::vector<int>& scores, std::vector<int>& clips)
        : Idx(idx), Score(score), Scores(scores), Clips(clips)
    {
    }

    uint16_t Idx = 0;
    uint8_t Score = 0;
    int Clip = 0;
    std::vector<int> Scores;
    std::vector<int> Clips;
};

struct BarcodeHitPair
{
    BarcodeHitPair() = default;
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

    operator std::string() const;
};

struct SequenceUtils
{
    static char Complement(char base);
    static std::string ReverseComplement(const std::string& input);
};

struct AdvancedFileUtils
{
    static std::string FilePrefixInfix(const std::string& path);
    static std::unique_ptr<BAM::internal::IQuery> BamQuery(const std::string& filePath);
};

struct AlignUtils
{
    static StripedSmithWaterman::Alignment Align(StripedSmithWaterman::Aligner& aligner,
                                                 const char* bases);

    static StripedSmithWaterman::Alignment AlignForward(StripedSmithWaterman::Aligner& aligner,
                                                        const Barcode& query);

    static StripedSmithWaterman::Alignment AlignRC(StripedSmithWaterman::Aligner& aligner,
                                                   const Barcode& query);
};
}
}

#include "pacbio/lima/internal/Lima.inl"