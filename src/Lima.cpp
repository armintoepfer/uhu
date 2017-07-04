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

#include <threadpool/ThreadPool.h>

#include <pacbio/lima/Lima.h>

namespace PacBio {
namespace Lima {
namespace {
const int leftAdapterFlag = static_cast<int>(BAM::LocalContextFlags::ADAPTER_BEFORE);
const int rightAdapterFlag = static_cast<int>(BAM::LocalContextFlags::ADAPTER_AFTER);
}

// ## BarcodeHit ##
BarcodeHitPair::operator std::string() const
{
    std::stringstream out;
    out << static_cast<int>(Left.Idx) << "\t" << static_cast<int>(Right.Idx) << "\t"
        << static_cast<int>(Left.Score) << "\t" << static_cast<int>(Right.Score) << "\t"
        << static_cast<int>(MeanScore) << "\t";
    if (Left.Clips.empty()) {
        out << "-";
    } else {
        out << Left.Clips.at(0);
        for (size_t i = 1; i < Left.Clips.size(); ++i)
            out << "," << Left.Clips.at(i);
    }
    out << "\t";
    if (Right.Clips.empty()) {
        out << "-";
    } else {
        out << Right.Clips.at(0);
        for (size_t i = 1; i < Right.Clips.size(); ++i)
            out << "," << Right.Clips.at(i);
        out << "\t";
    }
    if (Left.Scores.empty()) {
        out << "-";
    } else {
        out << static_cast<int>(Left.Scores.at(0));
        for (size_t i = 1; i < Left.Scores.size(); ++i)
            out << "," << static_cast<int>(Left.Scores.at(i));
    }
    out << "\t";
    if (Right.Scores.empty()) {
        out << "-";
    } else {
        out << static_cast<int>(Right.Scores.at(0));
        for (size_t i = 1; i < Right.Scores.size(); ++i)
            out << "," << static_cast<int>(Right.Scores.at(i));
    }
    return out.str();
}

// ## SequenceUtils ##
char SequenceUtils::Complement(char base)
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

std::string SequenceUtils::ReverseComplement(const std::string& input)
{
    std::string output;
    output.reserve(input.length());
    for (auto it = input.crbegin(); it != input.crend(); ++it)
        output.push_back(Complement(*it));
    return output;
}

// ## AdvancedFileUtils ##
std::string AdvancedFileUtils::FilePrefixInfix(const std::string& path)
{
    size_t fileStart = path.find_last_of("/");

    if (fileStart == std::string::npos) fileStart = -1;

    // increment beyond the '/'
    ++fileStart;

    size_t extStart = path.substr(fileStart, path.length() - fileStart).find_last_of(".");

    if (extStart == std::string::npos) return "";

    auto suffix = path.substr(fileStart, extStart);
    return suffix;
};

std::unique_ptr<BAM::internal::IQuery> AdvancedFileUtils::BamQuery(const std::string& filePath)
{
    BAM::DataSet ds(filePath);
    const auto filter = BAM::PbiFilter::FromDataSet(ds);
    std::unique_ptr<BAM::internal::IQuery> query(nullptr);
    if (filter.IsEmpty())
        query.reset(new BAM::EntireFileQuery(ds));
    else
        query.reset(new BAM::PbiFilterQuery(filter, ds));
    return query;
};

// ## AlignUtils ##
std::pair<int32_t, int32_t> AlignUtils::Align(const std::string& bcBases, const char* target,
                                              const int targetSize, std::vector<int32_t>& matrix,
                                              const AlignParameters& parameters) noexcept
{
    SWComputeMatrix(bcBases.c_str(), bcBases.size() + 1, target, targetSize + 1, false, matrix,
                    parameters);
    return SWLastRowMax(matrix, bcBases.size(), targetSize);
}

std::pair<int32_t, int32_t> AlignUtils::SWLastRowMax(const std::vector<int32_t>& matrix,
                                                     const int32_t queryLength,
                                                     const int32_t readLength) noexcept
{
    // Calculate the starting position of the last row
    const int32_t M = queryLength + 1;
    const int32_t N = readLength + 1;
    const int32_t beginLastRow = (M - 1) * N;

    // Find maximal score in last row and it's position
    int32_t maxScore = -1;
    int32_t endPos = 0;
    for (int32_t j = 0; j < N; ++j) {
        if (matrix[beginLastRow + j] > maxScore) {
            maxScore = matrix[beginLastRow + j];
            endPos = j;
        }
    }

    // Return the maximum score and position as a pair
    return std::make_pair(maxScore, endPos);
}

void AlignUtils::SWComputeMatrix(const char* const query, const int32_t M, const char* const read,
                                 const int32_t N, const bool globalInQuery,
                                 std::vector<int32_t>& matrix,
                                 const AlignParameters& parameters) noexcept
{
    matrix[0] = 0;

    if (globalInQuery)
        for (int32_t i = 1; i < M; ++i)
            matrix[i * N] = i * parameters.DeletionPenalty;
    else
        for (int32_t i = 1; i < M; ++i)
            matrix[i * N] = 0;

    for (int32_t j = 1; j < N; ++j)
        matrix[j] = 0;

    char iQuery;
    char iBeforeQuery;
    int32_t mismatchDelta = parameters.MatchScore - parameters.MismatchPenalty;
    int32_t insertionDelta = parameters.BranchPenalty - parameters.InsertionPenalty;
    for (int32_t i = 1; __builtin_expect(i < M, 1); ++i) {
        iQuery = query[i];
        iBeforeQuery = query[i - 1];
        if (__builtin_expect(i < M - 1, 1)) {
            for (int32_t j = 1; __builtin_expect(j < N, 1); ++j) {
                int32_t a = matrix[(i - 1) * N + j - 1] + parameters.MatchScore;
                int32_t b = matrix[i * N + j - 1] + parameters.BranchPenalty;
                const int32_t c = matrix[(i - 1) * N + j] + parameters.DeletionPenalty;
                if (read[j - 1] != iBeforeQuery) a -= mismatchDelta;
                if (read[j - 1] != iQuery) b -= insertionDelta;
                matrix[i * N + j] = std::max(a, std::max(b, c));
            }
        } else {
            for (int32_t j = 1; __builtin_expect(j < N, 1); ++j) {
                int32_t a = matrix[(i - 1) * N + j - 1] + parameters.MatchScore;
                int32_t b = matrix[i * N + j - 1] + parameters.InsertionPenalty;
                const int32_t c = matrix[(i - 1) * N + j] + parameters.DeletionPenalty;
                if (read[j - 1] != iBeforeQuery) a -= mismatchDelta;
                matrix[i * N + j] = std::max(a, std::max(b, c));
            }
        }
    }
}

// ## Summary ##
Summary::operator std::string() const
{
    std::stringstream summaryStream;
    summaryStream << "ZMWs input                    : " << NumZMWs << std::endl;
    summaryStream << "ZMWs above all thresholds (A) : " << AboveThresholds << std::endl;
    summaryStream << "ZMWs below any threshold  (B) : " << NumZMWs - AboveThresholds << std::endl;
    summaryStream << std::endl;
    summaryStream << "Marginals for (B)" << std::endl;
    summaryStream << "ZMWs below length threshold   : " << BelowMinLength << std::endl;
    summaryStream << "ZMWs below score threshold    : " << BelowMinScore << std::endl;
    summaryStream << "ZMWs below passes threshold   : " << BelowNumPasses << std::endl;
    summaryStream << std::endl;
    summaryStream << "For (A)" << std::endl;
    summaryStream << "ZMWs symmetric                : " << SymmetricCounts << std::endl;
    summaryStream << "ZMWs asymmetric               : " << AsymmetricCounts << std::endl;
    summaryStream << std::endl;
    summaryStream << "For (A)" << std::endl;
    summaryStream << "Reads above length            : " << SubreadAboveMinLength << std::endl;
    summaryStream << "Reads below length            : " << SubreadBelowMinLength << std::endl;
    return summaryStream.str();
}
}
}