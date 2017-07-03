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
#include <chrono>
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

#include <parasail.h>
#include <ssw_cpp.h>

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

#include <pacbio/parallel/WorkQueue.h>
#include <threadpool/ThreadPool.h>

#include <pacbio/lima/LimaSettings.h>

#include <pacbio/lima/LimaWorkflow.h>

namespace PacBio {
namespace Lima {
namespace {
const int leftAdapterFlag = static_cast<int>(BAM::LocalContextFlags::ADAPTER_BEFORE);
const int rightAdapterFlag = static_cast<int>(BAM::LocalContextFlags::ADAPTER_AFTER);
}

/// Fills out a supplied SW matrix.
/// \param  query       char* to the query
/// \param  queryLength Length of the query array
/// \param  read        char* to the read
/// \param  readLength  Length of the read array
/// \param  scoring     ScoringScheme for DP algorithm
/// \param  matrix      int32_t* to the SW matrix
static void SWComputeMatrix(const char* const query, const int32_t M, const char* const read,
                            const int32_t N, const bool globalInQuery, int32_t*& matrix,
                            const int32_t matchScore = 4, const int32_t mismatchPenalty = -13,
                            const int32_t deletionPenalty = -7, const int32_t insertionPenalty = -7,
                            const int32_t branchPenalty = -4) noexcept
{
    matrix[0] = 0;

    if (globalInQuery)
        for (int32_t i = 1; i < M; ++i)
            matrix[i * N] = i * deletionPenalty;
    else
        for (int32_t i = 1; i < M; ++i)
            matrix[i * N] = 0;

    for (int32_t j = 1; j < N; ++j)
        matrix[j] = 0;

    char iQuery;
    char iBeforeQuery;
    int32_t mismatchDelta = matchScore - mismatchPenalty;
    int32_t insertionDelta = branchPenalty - insertionPenalty;
    for (int32_t i = 1; __builtin_expect(i < M, 1); ++i) {
        iQuery = query[i];
        iBeforeQuery = query[i - 1];
        if (__builtin_expect(i < M - 1, 1)) {
            for (int32_t j = 1; __builtin_expect(j < N, 1); ++j) {
                // branch = match && read[j - 2] == read[j - 1];
                int32_t a = matrix[(i - 1) * N + j - 1] + matchScore;
                int32_t b = matrix[i * N + j - 1] + branchPenalty;
                int32_t c = matrix[(i - 1) * N + j] + deletionPenalty;
                if (read[j - 1] != iBeforeQuery) a -= mismatchDelta;
                if (read[j - 1] != iQuery) b -= insertionDelta;
                matrix[i * N + j] =
                    std::max(a, std::max(b, c));  //(a > b) ? ((a > c) ? a : c) : ((c > b) ? c : b);
            }
        } else {
            for (int32_t j = 1; __builtin_expect(j < N, 1); ++j) {
                // branch = match && read[j - 2] == read[j - 1];
                int32_t a = matrix[(i - 1) * N + j - 1] + matchScore;
                int32_t b = matrix[i * N + j - 1] + insertionPenalty;
                int32_t c = matrix[(i - 1) * N + j] + deletionPenalty;
                if (read[j - 1] != iBeforeQuery) a -= mismatchDelta;
                matrix[i * N + j] = std::max(a, std::max(b, c));
            }
        }
    }
}

/// Traverse the last row of an SW matrix (i.e. representing
///     alignments terminating with the last base of the query
///     sequence) and return the max score and it's position
///
/// \param  matrix       pointer to SW matrix
/// \param  queryLength  length of the query sequence
/// \param  readLength   length of the read sequence
///
/// \return  A std::pair of the max score and it's position
static std::pair<int32_t, int32_t> SWLastRowMax(const int32_t* matrix, const int32_t queryLength,
                                                const int32_t readLength)
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

static std::pair<int32_t, int32_t> AlignBarcode(const std::string& bcBases, const char* target,
                                                const int targetSize, int32_t*& matrix)
{
    SWComputeMatrix(bcBases.c_str(), bcBases.size() + 1, target, targetSize + 1, false, matrix);
    return SWLastRowMax(matrix, bcBases.size(), targetSize);
}

BarcodeHitPair LimaWorkflow::Tag(const std::vector<BAM::BamRecord> records,
                                 const std::vector<Barcode>& queries, const LimaSettings& settings)
{
    int barcodeLength = 0;
    for (const auto& q : queries)
        barcodeLength = std::max(barcodeLength, static_cast<int>(q.Bases.size()));
    int barcodeLengthWSpacing = barcodeLength * settings.WindowSizeMult;

    size_t numBarcodes = queries.size();
    size_t numRecords = records.size();

    int counterLeft = 0;
    int counterFullLeft = 0;
    std::vector<BarcodeHit> left(numBarcodes, numRecords);

    int counterRight = 0;
    int counterFullRight = 0;
    std::vector<BarcodeHit> right(numBarcodes, numRecords);

    // TODO: Inc those to fit fasta names
    for (size_t i = 0; i < numBarcodes; ++i) {
        left[i].Idx = i;
        right[i].Idx = i;
    }

    auto NormalizeScore = [&](const double& score) {
        return std::round(100.0 * score) / (barcodeLength * settings.MatchScore);
    };

    const int maxScoredReads = settings.MaxScoredReads;
    const bool maxScoring = settings.MaxScoredReads > 0;

    for (const auto& r : records) {
        // Activate if there is no context flag or if the left/right adapter is present
        const bool hasCX = r.HasLocalContextFlags();
        const bool hasAdapterLeft =
            ((hasCX && (r.LocalContextFlags() & leftAdapterFlag)) || !hasCX);
        const bool hasAdapterRight =
            ((hasCX && (r.LocalContextFlags() & rightAdapterFlag)) || !hasCX);
        const bool isFull = hasAdapterLeft && hasAdapterRight;

        const auto target = r.Sequence().c_str();
        const int targetLength = r.Sequence().size();

        if (hasAdapterLeft) {
            const auto targetSize = std::min(targetLength, barcodeLengthWSpacing);
            int32_t* matrix = new int32_t[(targetSize + 1) * (barcodeLength + 1)];

            for (size_t i = 0; i < queries.size(); ++i) {
                auto pair = AlignBarcode(queries[i].Bases, target, targetSize, matrix);
                const auto score = NormalizeScore(pair.first);
                const auto refEnd = pair.second;

                pair = AlignBarcode(queries[i].BasesRC, target, targetSize, matrix);
                const auto scoreRC = NormalizeScore(pair.first);
                const auto refEndRC = pair.second;

                if ((maxScoring && isFull && counterFullLeft < maxScoredReads) || !maxScoring) {
                    if (score > scoreRC)
                        left[i].AddWithSumScore(score, refEnd);
                    else
                        left[i].AddWithSumScore(scoreRC, refEndRC);
                } else {
                    if (score > scoreRC)
                        left[i].Add(score, refEnd);
                    else
                        left[i].Add(scoreRC, refEndRC);
                }
            }
            delete[] matrix;
            if ((maxScoring && isFull && counterFullLeft < maxScoredReads) || !maxScoring)
                ++counterFullLeft;
            ++counterLeft;
        } else {
            for (size_t i = 0; i < queries.size(); ++i) {
                left[i].Add(-1, 0);
            }
        }

        if (hasAdapterRight) {
            // Set reference as the last few bases
            int alignerRightBegin = std::max(targetLength - barcodeLengthWSpacing, 0);
            const auto targetSize = targetLength - alignerRightBegin;
            int32_t* matrix = new int32_t[(targetSize + 1) * (barcodeLength + 1)];
            for (size_t i = 0; i < queries.size(); ++i) {
                auto pair =
                    AlignBarcode(queries[i].Bases, target + alignerRightBegin, targetSize, matrix);
                const auto score = NormalizeScore(pair.first);
                const auto refEnd = pair.second;

                pair = AlignBarcode(queries[i].BasesRC, target + alignerRightBegin, targetSize,
                                    matrix);
                const auto scoreRC = NormalizeScore(pair.first);
                const auto refEndRC = pair.second;

                if ((maxScoring && isFull && counterFullRight < maxScoredReads) || !maxScoring) {
                    if (score > scoreRC)
                        right[i].AddWithSumScore(score, alignerRightBegin + refEnd);
                    else
                        right[i].AddWithSumScore(scoreRC, alignerRightBegin + refEndRC);
                } else {
                    if (score > scoreRC)
                        right[i].Add(score, alignerRightBegin + refEnd);
                    else
                        right[i].Add(scoreRC, alignerRightBegin + refEndRC);
                }
            }
            delete[] matrix;
            if ((maxScoring && isFull && counterFullRight < maxScoredReads) || !maxScoring)
                ++counterFullRight;
            ++counterRight;
        } else {
            for (size_t i = 0; i < queries.size(); ++i) {
                right[i].Add(-1, targetLength);
            }
        }
    }

    auto Compute = [&](std::vector<BarcodeHit>& fwd, int denominator) {
        if (denominator == 0) denominator = 1;
        for (size_t i = 0; i < numBarcodes; ++i) {
            fwd[i].Normalize(denominator);
        }
        auto cmp = [](const BarcodeHit& l, const BarcodeHit& r) { return l.Score < r.Score; };
        auto bestFwd = std::max_element(fwd.begin(), fwd.end(), cmp);

        return *bestFwd;
    };

    BarcodeHit leftBH;
    BarcodeHit rightBH;

    if (counterLeft > 0)
        leftBH = Compute(left, maxScoring ? counterFullLeft : counterLeft);
    else
        leftBH.Clips.resize(numRecords);

    if (counterRight > 0)
        rightBH = Compute(right, maxScoring ? counterFullRight : counterRight);
    else
        rightBH.Clips.resize(numRecords);

    return BarcodeHitPair(std::move(leftBH), std::move(rightBH));
}

struct TaskResult
{
    TaskResult(BarcodeHitPair&& bhp) : BHP(std::forward<BarcodeHitPair>(bhp)) {}
    std::vector<BAM::BamRecord> Records;
    std::string Report;
    BarcodeHitPair BHP;
    bool PassingFilters;
};

void WorkerThread(PacBio::Parallel::WorkQueue<std::vector<TaskResult>>& queue,
                  std::unique_ptr<BAM::BamWriter>& writer, const LimaSettings& settings,
                  const std::string& prefix, Summary& summary, BAM::BamHeader& header)
{
    std::map<uint16_t, std::map<uint16_t, int>> barcodePairCounts;

    std::ofstream report;
    if (!settings.NoReports) {
        report.open(prefix + ".demux.report");
        report << "ZMW\tIndexLeft\tIndexRight\tMeanScoreLeft\tMeanScoreRight\tMeanScore\tClipsL"
                  "eft\tClipsRight\tScoresLeft\tScoresRight"
               << std::endl;
    }

    std::map<std::pair<uint16_t, uint16_t>, std::vector<BAM::BamRecord>> barcodeToRecords;

    auto LambdaWorker = [&](std::vector<TaskResult>&& ps) {
        for (auto&& p : ps) {
            if (p.PassingFilters) {
                const auto leftIdx = p.BHP.Left.Idx;
                const auto rightIdx = p.BHP.Right.Idx;
                if (leftIdx == rightIdx)
                    ++summary.SymmetricCounts;
                else
                    ++summary.AsymmetricCounts;

                if ((settings.KeepSymmetric && leftIdx == rightIdx) || !settings.KeepSymmetric) {
                    if (settings.SplitBam)
                        for (auto&& r : p.Records)
                            barcodeToRecords[std::make_pair(leftIdx, rightIdx)].emplace_back(
                                std::move(r));
                    else if (!settings.NoBam)
                        for (auto& r : p.Records) {
                            writer->Write(r);
                        }
                    if (!settings.NoReports) ++barcodePairCounts[leftIdx][rightIdx];
                }
            }
            if (!settings.NoReports) report << p.Report << std::endl;
        }
    };

    while (queue.ConsumeWith(LambdaWorker)) {
    }

    if (settings.SplitBam) {
        for (const auto& bc_records : barcodeToRecords) {
            std::stringstream fileName;
            fileName << prefix;
            fileName << ".";
            fileName << static_cast<int>(bc_records.first.first);
            fileName << "-";
            fileName << static_cast<int>(bc_records.first.second);
            fileName << ".demux.bam";
            BAM::BamWriter writer(fileName.str(), header);
            for (const auto& r : bc_records.second)
                writer.Write(r);
        }
    }

    if (!settings.NoReports) {
        std::ofstream summaryStream(prefix + ".demux.summary");
        summaryStream << summary;

        std::ofstream counts(prefix + ".demux.counts");
        counts << "IndexLeft\tIndexRight\tCounts" << std::endl;
        for (const auto& left_right_counts : barcodePairCounts)
            for (const auto& right_counts : left_right_counts.second)
                counts << static_cast<int>(left_right_counts.first) << "\t"
                       << static_cast<int>(right_counts.first) << "\t" << right_counts.second
                       << std::endl;
    }
}

void LimaWorkflow::Process(const LimaSettings& settings,
                           const std::vector<std::string>& datasetPaths,
                           const std::vector<Barcode>& barcodes)
{

    // Single writer for non-split mode
    std::unique_ptr<BAM::BamWriter> writer;
    // Header can be used for split mode
    BAM::BamHeader header;
    // Treat every dataset as an individual entity
    for (const auto& datasetPath : datasetPaths) {
        writer.reset(nullptr);
        std::atomic_int threadCount(0);

        std::string prefix = AdvancedFileUtils::FilePrefixInfix(datasetPath);
        Summary summary;
        // Individual queue per dataset
        PacBio::Parallel::WorkQueue<std::vector<TaskResult>> workQueue(settings.NumThreads);
        std::future<void> workerThread =
            std::async(std::launch::async, WorkerThread, std::ref(workQueue), std::ref(writer),
                       std::ref(settings), std::ref(prefix), std::ref(summary), std::ref(header));
        // Get a query to the underlying BAM files, respecting filters
        auto query = AdvancedFileUtils::BamQuery(datasetPath);
        auto Submit = [&barcodes, &settings, &summary,
                       &threadCount](std::vector<std::vector<BAM::BamRecord>> chunk) {
            ++threadCount;
            std::vector<TaskResult> results;
            for (const auto& records : chunk) {
                TaskResult result{LimaWorkflow::Tag(records, barcodes, settings)};
                const auto& bhp = result.BHP;

                bool aboveMinScore = bhp.MeanScore >= settings.MinScore;
                bool aboveMinLength = false;

                if (bhp.Right.Clips.size() != bhp.Left.Clips.size() ||
                    bhp.Right.Clips.size() != records.size())
                    throw std::runtime_error("Internal error, clips sizes not equal! " +
                                             records.at(0).FullName() + " " +
                                             std::to_string(bhp.Left.Clips.size()) + " " +
                                             std::to_string(bhp.Right.Clips.size()));
                for (size_t i = 0; i < bhp.Right.Clips.size(); ++i) {
                    if (bhp.Right.Clips.at(i) - bhp.Left.Clips.at(i) > settings.MinLength) {
                        aboveMinLength = true;
                        break;
                    }
                }
                result.PassingFilters = aboveMinScore && aboveMinLength;

                if (!settings.NoReports)
                    result.Report =
                        std::to_string(records.at(0).HoleNumber()) + "\t" + std::string(bhp);

                if (aboveMinScore && aboveMinLength) {
                    if (!settings.NoBam) {
                        for (size_t i = 0; i < records.size(); ++i) {
                            int clipLeft = bhp.Left.Clips.at(i);
                            int clipRight = bhp.Right.Clips.at(i);
                            if (clipRight - clipLeft > settings.MinLength) {
                                auto r = records[i];
                                if (r.HasQueryStart()) {
                                    clipLeft += r.QueryStart();
                                    clipRight += r.QueryStart();
                                }
                                r.Clip(BAM::ClipType::CLIP_TO_QUERY, clipLeft, clipRight);
                                r.Barcodes(std::make_pair(bhp.Left.Idx, bhp.Right.Idx));
                                r.BarcodeQuality(bhp.MeanScore);
                                result.Records.emplace_back(std::move(r));
                                ++summary.SubreadAboveMinLength;
                            } else {
                                ++summary.SubreadBelowMinLength;
                            }
                        }
                    }
                    ++summary.AboveThresholds;
                } else if (!aboveMinLength && !aboveMinScore) {
                    ++summary.BelowBoth;
                } else if (!aboveMinLength) {
                    ++summary.BelowMinLength;
                } else if (!aboveMinScore) {
                    ++summary.BelowMinScore;
                }
                results.emplace_back(std::move(result));
            }
            --threadCount;
            return results;
        };

        int zmwNum = -1;

        std::vector<std::vector<BAM::BamRecord>> chunk;
        std::vector<BAM::BamRecord> records;
        for (auto& r : *query) {
            if (!writer)
                if (!settings.NoBam && !settings.SplitBam)
                    writer.reset(new BAM::BamWriter(
                        prefix + ".demux.bam", r.Header().DeepCopy(),
                        BAM::BamWriter::CompressionLevel::CompressionLevel_0, settings.NumThreads));
            if (settings.SplitBam) header = r.Header().DeepCopy();

            if (zmwNum == -1) {
                zmwNum = r.HoleNumber();
            } else if (zmwNum != r.HoleNumber()) {
                if (!records.empty()) chunk.emplace_back(std::move(records));
                if (chunk.size() == 1) {
                    workQueue.ProduceWith(Submit, chunk);
                    chunk.clear();
                }
                zmwNum = r.HoleNumber();
                records = std::vector<BAM::BamRecord>();
            }
            records.push_back(r);
        }
        if (!chunk.empty()) workQueue.ProduceWith(Submit, chunk);
        workQueue.Finalize();
        while (threadCount > 0) {
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        }
    }
}

int LimaWorkflow::Runner(const PacBio::CLI::Results& options)
{
    // Check args size, as pbcopper does not enforce the correct number
    if (options.PositionalArguments().empty()) {
        std::cerr << "ERROR: Please provide BAM and Barcode input, see --help" << std::endl;
        return EXIT_FAILURE;
    }

    const LimaSettings settings(options);
    std::vector<std::string> datasetPaths;
    std::vector<Barcode> barcodes;
    ParsePositionalArgs(options.PositionalArguments(), &datasetPaths, &barcodes);

    Process(settings, datasetPaths, barcodes);

    return EXIT_SUCCESS;
}

void LimaWorkflow::ParsePositionalArgs(const std::vector<std::string>& args,
                                       std::vector<std::string>* datasetPaths,
                                       std::vector<Barcode>* barcodes)
{
    std::vector<std::string> fastaPaths;
    for (const auto& i : args) {
        const bool fileExist = PacBio::Utility::FileExists(i);
        if (!fileExist) {
            throw std::runtime_error("File does not exist: " + i);
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
}
}  // ::PacBio::Lima