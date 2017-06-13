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

#include <threadpool/ThreadPool.h>

#include <pacbio/lima/LimaRawSettings.h>

#include <pacbio/lima/LimaRawWorkflow.h>

namespace PacBio {
namespace Lima {
namespace RAW {
namespace {
const int leftAdapterFlag = static_cast<int>(BAM::LocalContextFlags::ADAPTER_BEFORE);
const int rightAdapterFlag = static_cast<int>(BAM::LocalContextFlags::ADAPTER_AFTER);
}
// ## Lima #
BarcodeHitPair LimaWorkflow::Tag(const std::vector<BAM::BamRecord> records,
                                 const std::vector<Barcode>& queries, const LimaSettings& settings)
{
    int barcodeLength = 0;
    for (const auto& q : queries)
        barcodeLength = std::max(barcodeLength, static_cast<int>(q.Bases.size()));
    int barcodeLengthWSpacing = barcodeLength * settings.WindowSizeMult;

    size_t numBarcodes = queries.size();

    int counterLeft = 0;
    std::vector<double> scoresLeft(numBarcodes, 0);
    std::vector<double> scoresLeftRC(numBarcodes, 0);
    std::vector<std::vector<int>> scoresLeftV(numBarcodes);
    std::vector<std::vector<int>> scoresLeftRCV(numBarcodes);
    std::vector<std::vector<int>> clipsLeftV(numBarcodes);
    std::vector<std::vector<int>> clipsLeftRCV(numBarcodes);

    int counterRight = 0;
    std::vector<double> scoresRight(numBarcodes, 0);
    std::vector<double> scoresRightRC(numBarcodes, 0);
    std::vector<std::vector<int>> scoresRightV(numBarcodes);
    std::vector<std::vector<int>> scoresRightRCV(numBarcodes);
    std::vector<std::vector<int>> clipsRightV(numBarcodes);
    std::vector<std::vector<int>> clipsRightRCV(numBarcodes);

    auto NormalizeScore = [&](const double& score) {
        return std::round(100.0 * score) / (barcodeLength * settings.MatchScore);
    };

    int recordCounter = 0;
    for (const auto& r : records) {
        const int cx = r.LocalContextFlags();
        bool hasAdapterLeft = cx & leftAdapterFlag;
        bool hasAdapterRight = cx & rightAdapterFlag;

        const auto target = r.Sequence();
        const int targetLength = target.size();

        StripedSmithWaterman::Aligner alignerLeft(settings.MatchScore, settings.MismatchPenalty,
                                                  settings.GapOpenPenalty, settings.GapExtPenalty);
        StripedSmithWaterman::Aligner alignerRight(settings.MatchScore, settings.MismatchPenalty,
                                                   settings.GapOpenPenalty, settings.GapExtPenalty);

        if (hasAdapterLeft) {
            alignerLeft.SetReferenceSequence(target.c_str(),
                                             std::min(targetLength, barcodeLengthWSpacing));
            for (size_t i = 0; i < queries.size(); ++i) {
                const auto align = AlignUtils::AlignForward(alignerLeft, queries[i]);
                scoresLeft[i] += align.sw_score;
                scoresLeftV[i].push_back(NormalizeScore(align.sw_score));
                clipsLeftV[i].push_back(align.ref_end);

                const auto alignRC = AlignUtils::AlignRC(alignerLeft, queries[i]);
                scoresLeftRC[i] += alignRC.sw_score;
                scoresLeftRCV[i].push_back(NormalizeScore(alignRC.sw_score));
                clipsLeftRCV[i].push_back(alignRC.ref_end);
            }
            ++counterLeft;
        } else {
            for (size_t i = 0; i < queries.size(); ++i) {
                scoresLeftV[i].push_back(-1);
                scoresLeftRCV[i].push_back(-1);
                clipsLeftV[i].push_back(0);
                clipsLeftRCV[i].push_back(0);
            }
        }
        if (hasAdapterRight) {
            int alignerRightBegin = std::max(targetLength - barcodeLengthWSpacing, 0);
            alignerRight.SetReferenceSequence(target.c_str() + alignerRightBegin,
                                              targetLength - alignerRightBegin);
            for (size_t i = 0; i < queries.size(); ++i) {
                const auto align = AlignUtils::AlignForward(alignerRight, queries[i]);
                scoresRight[i] += align.sw_score;
                scoresRightV[i].push_back(NormalizeScore(align.sw_score));
                clipsRightV[i].push_back(alignerRightBegin + align.ref_begin);

                const auto alignRC = AlignUtils::AlignRC(alignerRight, queries[i]);
                scoresRightRC[i] += alignRC.sw_score;
                scoresRightRCV[i].push_back(NormalizeScore(alignRC.sw_score));
                clipsRightRCV[i].push_back(alignerRightBegin + alignRC.ref_begin);
            }
            ++counterRight;
        } else {
            for (size_t i = 0; i < queries.size(); ++i) {
                scoresRightV[i].push_back(-1);
                scoresRightRCV[i].push_back(-1);
                clipsRightV[i].push_back(targetLength);
                clipsRightRCV[i].push_back(targetLength);
            }
        }
        ++recordCounter;
    }

    auto GetBestIndex = [&](std::vector<double>& v) {
        std::vector<size_t> idx(v.size());
        std::iota(idx.begin(), idx.end(), 0);
        std::sort(idx.begin(), idx.end(), [&v](size_t i1, size_t i2) { return v[i1] > v[i2]; });

        return std::make_pair(idx.front(), NormalizeScore(v.at(idx.front())));
    };

    auto Compute = [&](std::vector<double>& scores, std::vector<double>& scoresRC, bool left) {
        BarcodeHit bc;

        int idxFwd;
        int scoreFwd;
        std::tie(idxFwd, scoreFwd) = GetBestIndex(scores);

        int scoreRev;
        int idxRev;
        std::tie(idxRev, scoreRev) = GetBestIndex(scoresRC);

        int counter = 0;
        if (scoreFwd > scoreRev) {
            bc.Score = scoreFwd;
            bc.Idx = idxFwd;

            if (left) {
                bc.Clips = clipsLeftV.at(bc.Idx);
                bc.Scores = scoresLeftV.at(bc.Idx);
            } else {
                bc.Clips = clipsRightV.at(bc.Idx);
                bc.Scores = scoresRightV.at(bc.Idx);
            }
        } else {
            bc.Score = scoreRev;
            bc.Idx = idxRev;

            if (left) {
                bc.Clips = clipsLeftRCV.at(bc.Idx);
                bc.Scores = scoresLeftRCV.at(bc.Idx);
            } else {
                bc.Clips = clipsRightRCV.at(bc.Idx);
                bc.Scores = scoresRightRCV.at(bc.Idx);
            }
        }

        return bc;
    };

    BarcodeHit left;
    BarcodeHit right;

    if (counterLeft > 0) {
        for (size_t i = 0; i < numBarcodes; ++i) {
            scoresLeft[i] /= counterLeft;
            scoresLeftRC[i] /= counterLeft;
        }
        left = Compute(scoresLeft, scoresLeftRC, true);
    }

    if (counterRight > 0) {
        for (size_t i = 0; i < numBarcodes; ++i) {
            scoresRight[i] /= counterRight;
            scoresRightRC[i] /= counterRight;
        }
        right = Compute(scoresRight, scoresRightRC, false);
    }

    return BarcodeHitPair(std::move(left), std::move(right));
}

void LimaWorkflow::Process(const LimaSettings& settings,
                           const std::vector<std::string>& datasetPaths,
                           const std::vector<Barcode>& barcodes)
{
    std::unique_ptr<BAM::BamWriter> writer;
    std::map<int, std::vector<BAM::BamRecord>> map;
    BAM::BamHeader header;
    for (const auto& datasetPath : datasetPaths) {
        auto query = AdvancedFileUtils::BamQuery(datasetPath);
        std::vector<Uhu::Threadpool::ThreadPool::TaskFuture<
            std::tuple<std::vector<BAM::BamRecord>, std::string, BarcodeHitPair, bool>>>
            v;
        Summary summary;
        auto Submit = [&](const std::vector<BAM::BamRecord>& records) {
            std::vector<BAM::BamRecord> recordOut;
            int cumRecordLength = 0;
            std::string report;
            BarcodeHitPair bh = LimaWorkflow::Tag(records, barcodes, settings);

            bool aboveMinScore = bh.MeanScore >= settings.MinScore;
            if (!settings.NoReports)
                report = std::to_string(records.at(0).HoleNumber()) + "\t" + std::string(bh);
            if (aboveMinScore) {
                if (!settings.NoBam) {
                    for (size_t i = 0; i < records.size(); ++i) {
                        auto r = records[i];
                        int clipLeft = bh.Left.Clips.at(i);
                        int clipRight = bh.Right.Clips.at(i);

                        r.Clip(BAM::ClipType::CLIP_TO_QUERY, r.QueryStart() + clipLeft,
                               r.QueryStart() + clipRight);
                        r.Barcodes(std::make_pair(bh.Left.Idx, bh.Right.Idx));
                        r.BarcodeQuality(bh.MeanScore);
                        recordOut.emplace_back(std::move(r));
                    }
                }
                ++summary.AboveThresholds;
            } else if (!aboveMinScore) {
                ++summary.BelowMinScore;
            }
            return std::make_tuple(std::move(records), report, bh, aboveMinScore);
        };

        std::string prefix = AdvancedFileUtils::FilePrefixInfix(datasetPath);

        int zmwNum = -1;
        std::vector<BAM::BamRecord> records;
        for (auto& r : *query) {
            if (!writer && !settings.NoBam && !settings.SplitBam)
                writer.reset(new BAM::BamWriter(prefix + ".demux.bam", r.Header().DeepCopy()));
            if (settings.SplitBam) header = r.Header().DeepCopy();

            if (zmwNum == -1) {
                zmwNum = r.HoleNumber();
            } else if (zmwNum != r.HoleNumber()) {
                if (!records.empty())
                    v.push_back(Uhu::Threadpool::DefaultThreadPool::submitJob(Submit, records));
                zmwNum = r.HoleNumber();
                records.clear();
            }
            records.push_back(r);
        }
        if (!records.empty())
            v.push_back(Uhu::Threadpool::DefaultThreadPool::submitJob(Submit, records));

        std::map<uint8_t, std::map<uint8_t, int>> barcodePairCounts;

        std::ofstream report;
        if (!settings.NoReports) {
            report.open(prefix + ".demux.report");
            report << "ZMW\tIndexLeft\tIndexRight\tMeanScoreLeft\tMeanScoreRight\tMeanScore\tClipsL"
                      "eft\tClipsRight\tScoresLeft\tScoresRight"
                   << std::endl;
        }

        std::map<std::pair<uint8_t, uint8_t>, std::vector<BAM::BamRecord>> barcodeToRecords;

        for (auto& item : v) {
            auto p = item.get();
            if (std::get<3>(p)) {
                const auto leftIdx = std::get<2>(p).Left.Idx;
                const auto rightIdx = std::get<2>(p).Right.Idx;
                if (leftIdx == rightIdx)
                    ++summary.SymmetricCounts;
                else
                    ++summary.AsymmetricCounts;

                if ((settings.KeepSymmetric && leftIdx == rightIdx) || !settings.KeepSymmetric) {
                    if (settings.SplitBam)
                        for (auto&& r : std::get<0>(p))
                            barcodeToRecords[std::make_pair(leftIdx, rightIdx)].emplace_back(
                                std::move(r));
                    else if (!settings.NoBam)
                        for (auto& r : std::get<0>(p))
                            writer->Write(r);
                    if (!settings.NoReports) ++barcodePairCounts[leftIdx][rightIdx];
                }
            }
            if (!settings.NoReports) report << std::get<1>(p) << std::endl;
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
        writer.reset(nullptr);
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
}
}