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

BarcodeHitPair LimaWorkflow::Tag(const std::vector<BAM::BamRecord> records,
                                 const std::vector<Barcode>& queries, const LimaSettings& settings,
                                 const AlignParameters& alignParameters)
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

    std::vector<int32_t> matrix;
    for (const auto& r : records) {
        // Activate if there is no context flag or if the left/right adapter is present
        const bool hasCX = r.HasLocalContextFlags();
        const bool hasAdapterLeft =
            ((hasCX && (r.LocalContextFlags() & leftAdapterFlag)) || !hasCX);
        const bool hasAdapterRight =
            ((hasCX && (r.LocalContextFlags() & rightAdapterFlag)) || !hasCX);
        const bool isFull = hasAdapterLeft && hasAdapterRight;

        const auto seq = r.Sequence();
        const auto target = seq.c_str();
        const int targetLength = seq.size();
        const auto targetSizeLeft = std::min(targetLength, barcodeLengthWSpacing);
        if (hasAdapterLeft && targetSizeLeft > 0) {
            matrix.resize((targetSizeLeft + 1) * (barcodeLength + 1));
            for (size_t i = 0; i < queries.size(); ++i) {
                auto pair = AlignUtils::Align(queries[i].Bases, target, targetSizeLeft, matrix,
                                              alignParameters);
                const auto score = NormalizeScore(pair.first);
                const auto refEnd = pair.second;

                pair = AlignUtils::Align(queries[i].BasesRC, target, targetSizeLeft, matrix,
                                         alignParameters);
                const auto scoreRC = NormalizeScore(pair.first);
                const auto refEndRC = pair.second;

                if ((!maxScoring || (maxScoring && isFull && counterFullLeft < maxScoredReads)) &&
                    (!settings.PerSubread || (settings.PerSubread && isFull))) {
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
            if (!maxScoring || (maxScoring && isFull && counterFullLeft < maxScoredReads))
                ++counterFullLeft;
            ++counterLeft;
        } else {
            for (size_t i = 0; i < queries.size(); ++i) {
                left[i].Add(-1, 0);
            }
        }

        int alignerRightBegin = std::max(targetLength - barcodeLengthWSpacing, 0);
        const auto targetSizeRight = targetLength - alignerRightBegin;
        if (hasAdapterRight && targetSizeRight) {
            // Set reference as the last few bases
            matrix.resize((targetSizeRight + 1) * (barcodeLength + 1));
            for (size_t i = 0; i < queries.size(); ++i) {
                auto pair = AlignUtils::Align(queries[i].Bases, target + alignerRightBegin,
                                              targetSizeRight, matrix, alignParameters);
                const auto score = NormalizeScore(pair.first);
                const auto refEnd = pair.second;

                pair = AlignUtils::Align(queries[i].BasesRC, target + alignerRightBegin,
                                         targetSizeRight, matrix, alignParameters);
                const auto scoreRC = NormalizeScore(pair.first);
                const auto refEndRC = pair.second;
                if ((!maxScoring || (maxScoring && isFull && counterFullRight < maxScoredReads)) &&
                    (!settings.PerSubread || (settings.PerSubread && isFull))) {
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
            if (!maxScoring || (maxScoring && isFull && counterFullRight < maxScoredReads))
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
    else {
        leftBH.Clips.resize(numRecords);
        leftBH.Scores.resize(numRecords, -1);
    }

    if (counterRight > 0)
        rightBH = Compute(right, maxScoring ? counterFullRight : counterRight);
    else {
        rightBH.Clips.resize(numRecords);
        rightBH.Scores.resize(numRecords, -1);
    }

    return BarcodeHitPair(std::move(leftBH), std::move(rightBH));
}

struct TaskResult
{
    TaskResult(BarcodeHitPair&& bhp) : BHP(std::forward<BarcodeHitPair>(bhp)) {}
    std::vector<BAM::BamRecord> Records;
    std::string Report;
    BarcodeHitPair BHP;
    bool PassingFilters;
    int NumPasses;
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
                  "eft\tClipsRight\tScoresLeft\tScoresRight\tNumPasses\tPassing"
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

                if (((settings.KeepSymmetric && p.BHP.Left.Idx == p.BHP.Right.Idx) ||
                     !settings.KeepSymmetric) &&
                    (!settings.PerSubread ||
                     (settings.PerSubread && p.BHP.Left.Score > 0 && p.BHP.Right.Score > 0))) {
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
            if (!settings.NoReports)
                report << p.Report << "\t" << p.NumPasses << "\t"
                       << (p.PassingFilters &&
                           ((settings.KeepSymmetric && p.BHP.Left.Idx == p.BHP.Right.Idx) ||
                            !settings.KeepSymmetric) &&
                           (!settings.PerSubread ||
                            (settings.PerSubread && p.BHP.Left.Score > 0 && p.BHP.Right.Score > 0)))
                       << std::endl;
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
    AlignParameters alignParameters(settings.MatchScore, -settings.MismatchPenalty,
                                    -settings.DeletionPenalty, -settings.InsertionPenalty,
                                    -settings.BranchPenalty);

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
        auto Submit = [&barcodes, &settings, &summary, &threadCount,
                       &alignParameters](std::vector<std::vector<BAM::BamRecord>> chunk) {
            ++threadCount;
            std::vector<TaskResult> results;
            for (const auto& records : chunk) {
                TaskResult result{LimaWorkflow::Tag(records, barcodes, settings, alignParameters)};
                ++summary.NumZMWs;
                const auto& bhp = result.BHP;

                bool aboveMinScore = bhp.MeanScore >= settings.MinScore;
                bool aboveMinLength = false;

                if (bhp.Right.Clips.size() != bhp.Left.Clips.size() ||
                    bhp.Right.Clips.size() != records.size()) {
                    std::cerr << "Internal error, clips sizes not equal!" << std::endl;
                    exit(1);
                }
                for (size_t i = 0; i < bhp.Right.Clips.size(); ++i) {
                    if (bhp.Right.Clips.at(i) - bhp.Left.Clips.at(i) > settings.MinLength) {
                        aboveMinLength = true;
                        break;
                    }
                }
                int numPasses = 0;
                if (bhp.Right.Scores.size() != bhp.Left.Scores.size()) {
                    std::cerr << "Internal error, scores sizes not equal!" << std::endl;
                    exit(1);
                }
                for (size_t i = 0; i < bhp.Right.Scores.size(); ++i) {
                    numPasses += bhp.Right.Scores.at(i) != -1 && bhp.Left.Scores.at(i) != -1;
                }
                bool aboveNumPasses = numPasses >= settings.MinPasses;
                result.PassingFilters = aboveMinScore && aboveMinLength && aboveNumPasses;
                result.NumPasses = numPasses;

                if (!settings.NoReports)
                    result.Report =
                        std::to_string(records.at(0).HoleNumber()) + "\t" + std::string(bhp);

                if (aboveMinScore && aboveMinLength && aboveNumPasses) {
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
                } else {
                    if (!aboveMinLength) ++summary.BelowMinLength;
                    if (!aboveMinScore) ++summary.BelowMinScore;
                    if (!aboveNumPasses) ++summary.BelowNumPasses;
                }
                results.emplace_back(std::move(result));
            }
            --threadCount;
            return results;
        };

        int zmwNum = -1;

        std::vector<std::vector<BAM::BamRecord>> chunk;
        std::vector<BAM::BamRecord> records;
        int chunkNum = 0;
        for (auto& r : *query) {
            if (!writer)
                if (!settings.NoBam && !settings.SplitBam)
                    writer.reset(new BAM::BamWriter(
                        prefix + ".demux.bam", r.Header().DeepCopy(),
                        BAM::BamWriter::CompressionLevel::CompressionLevel_0, settings.NumThreads));
            if (settings.SplitBam) header = r.Header().DeepCopy();

            if (zmwNum == -1) {
                zmwNum = r.HoleNumber();
            } else if (settings.PerSubread || (!settings.PerSubread && zmwNum != r.HoleNumber())) {
                if (!records.empty()) chunk.emplace_back(std::move(records));
                if (static_cast<int>(chunk.size()) >= settings.Chunks) {
                    workQueue.ProduceWith(Submit, chunk);
                    chunk.clear();
                }
                zmwNum = r.HoleNumber();
                records = std::vector<BAM::BamRecord>();
                ++chunkNum;
            }
            records.push_back(r);
        }
        if (!records.empty()) chunk.emplace_back(std::move(records));
        if (!chunk.empty()) workQueue.ProduceWith(Submit, chunk);
        workQueue.Finalize();
        while (threadCount > 0)
            std::this_thread::sleep_for(std::chrono::milliseconds(100));
        workerThread.wait();
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