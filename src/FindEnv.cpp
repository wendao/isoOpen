#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <numeric>
#include "Envelope.h"

/**
 * Candidate envelope with score information for non-greedy selection.
 */
struct CandidateEnvelope {
    Envelope* envelope;           // The envelope object
    std::vector<size_t> indices;  // Peak indices used by this envelope
    double score;                 // Combined score for ranking
};

/**
 * Calculate spacing error (average ppm deviation from expected deltaM).
 */
double CalculateSpacingError(const std::vector<double>& mz, int charge) {
    if (mz.size() < 2) return 0.0;

    const double C13C12 = 1.003354835;
    const double expectedDelta = C13C12 / charge;
    double totalError = 0.0;

    for (size_t i = 1; i < mz.size(); ++i) {
        double actualDelta = mz[i] - mz[i-1];
        double errorPPM = std::abs(actualDelta - expectedDelta) / expectedDelta * 1e6;
        totalError += errorPPM;
    }

    return totalError / (mz.size() - 1);
}

/**
 * Find all candidate envelopes without greedy assignment.
 *
 * Instead of immediately marking peaks as used, we find ALL possible envelopes
 * and score them by quality metrics.
 */
void FindAllCandidates(std::vector<CandidateEnvelope>& candidates,
                       const std::vector<double>& mz,
                       const std::vector<double>& intensity,
                       int lenMin, int lenMax)
{
    const double C13C12 = 1.003354835;
    const double tol = 0.01;
    const std::vector<int> charges = {7, 6, 5, 4, 3, 2, 1};
    size_t n = mz.size();

    for (int charge : charges) {
        double deltaM = C13C12 / charge;

        for (size_t i = 0; i < n; ++i) {
            // Early termination if we can't reach minimum length
            if (mz[i] + (lenMin - 1) * deltaM > mz.back() + tol)
                break;

            std::vector<size_t> idxs;
            std::vector<double> envelopeMzs;
            idxs.reserve(lenMax);
            envelopeMzs.reserve(lenMax);
            idxs.push_back(i);
            envelopeMzs.push_back(mz[i]);

            double target = mz[i] + deltaM;
            size_t searchStart = idxs.back() + 1;  // Start searching from last peak + 1

            while ((int)idxs.size() < lenMax && searchStart < n) {
                // Linear search with early termination - cache-friendly for small tolerance
                // Envelopes grow monotonically, so we only search forward
                const double tolSquared = tol * tol;  // Precompute for squared comparison

                bool found = false;
                while (searchStart < n) {
                    double mzVal = mz[searchStart];
                    double diff = mzVal - target;

                    // If we've passed the target beyond tolerance, stop
                    if (diff > tol) {
                        break;
                    }

                    // Check if within tolerance using squared comparison (avoids sqrt)
                    if (diff >= -tol && diff * diff <= tolSquared) {
                        idxs.push_back(searchStart);
                        envelopeMzs.push_back(mzVal);
                        searchStart++;  // Move past this peak for next iteration
                        found = true;
                        break;
                    }

                    searchStart++;
                }

                if (!found) {
                    break;  // No matching peak found within tolerance
                }

                target += deltaM;
            }

            if ((int)idxs.size() >= lenMin) {
                CandidateEnvelope cand;
                cand.indices = std::move(idxs);

                // Create envelope
                Envelope* e = new Envelope(charge);
                for (size_t idx : cand.indices) {
                    e->AddPeak(mz[idx], intensity[idx]);
                }
                cand.envelope = e;

                // Calculate score based on multiple metrics
                double totalIntensity = 0.0;
                for (size_t idx : cand.indices) {
                    totalIntensity += intensity[idx];
                }

                // Spacing quality: lower error is better
                double spacingErrorPPM = CalculateSpacingError(envelopeMzs, charge);

                // Length score: longer envelopes are generally better (diminishing returns)
                double lengthScore = std::log2(static_cast<double>(cand.indices.size()));

                // Combined score: prioritize intensity, penalize spacing error, reward length
                // Higher is better
                cand.score = std::log(totalIntensity + 1.0) * 2.0
                           - spacingErrorPPM * 0.001
                           + lengthScore * 0.5;

                candidates.push_back(std::move(cand));
            }
        }
    }
}

/**
 * Check if two candidate envelopes overlap significantly.
 * Returns true if they share more than threshold fraction of peaks.
 */
bool OverlapSignificant(const CandidateEnvelope& a,
                        const CandidateEnvelope& b,
                        double overlapThreshold = 0.8)
{
    // Use two-pointer technique for sorted indices
    size_t i = 0, j = 0;
    int common = 0;

    // Both indices are sorted, use merge-like approach
    while (i < a.indices.size() && j < b.indices.size()) {
        if (a.indices[i] == b.indices[j]) {
            ++common;
            ++i;
            ++j;
        } else if (a.indices[i] < b.indices[j]) {
            ++i;
        } else {
            ++j;
        }
    }

    double overlapRatio = static_cast<double>(common) / std::min(a.indices.size(), b.indices.size());
    return overlapRatio >= overlapThreshold;
}

/**
 * Check overlap with lower threshold (0.5) for building conflict graph.
 */
bool OverlapSignificantLight(const CandidateEnvelope& a,
                            const CandidateEnvelope& b)
{
    size_t i = 0, j = 0;
    int common = 0;

    while (i < a.indices.size() && j < b.indices.size()) {
        if (a.indices[i] == b.indices[j]) {
            ++common;
            ++i;
            ++j;
        } else if (a.indices[i] < b.indices[j]) {
            ++i;
        } else {
            ++j;
        }
    }

    double overlapRatio = static_cast<double>(common) / std::min(a.indices.size(), b.indices.size());
    return overlapRatio >= 0.5;
}

/**
 * Check if candidate can be selected given current selection and conflict graph.
 */
bool CanSelect(size_t candIdx,
               const std::vector<bool>& selected,
               const std::vector<std::vector<int>>& adj) {
    for (int neighbor : adj[candIdx]) {
        if (selected[neighbor]) {
            return false;
        }
    }
    return true;
}

/**
 * Find a candidate that conflicts with given candidate and has lower score.
 */
int FindConflictingWithLowerScore(size_t candIdx,
                                   const std::vector<bool>& selected,
                                   const std::vector<CandidateEnvelope>& candidates,
                                   const std::vector<std::vector<int>>& adj) {
    for (int neighbor : adj[candIdx]) {
        if (selected[neighbor] && candidates[neighbor].score < candidates[candIdx].score) {
            return neighbor;
        }
    }
    return -1;
}

/**
 * Improved non-greedy selection using local search optimization.
 *
 * Algorithm:
 * 1. Sort candidates by score (descending)
 * 2. Build conflict graph (overlap >= 50%)
 * 3. Initial greedy selection
 * 4. Local search: try to replace lower-scoring selected with higher-scoring unselected
 */
void ImprovedSelection(std::vector<Envelope*>& envelopes,
                      std::vector<CandidateEnvelope>& candidates) {
    if (candidates.empty()) return;

    // Sort candidates by score (descending)
    std::sort(candidates.begin(), candidates.end(),
              [](const CandidateEnvelope& a, const CandidateEnvelope& b) {
                  return a.score > b.score;
              });

    size_t n = candidates.size();
    size_t peakCount = 0;
    for (const auto& c : candidates) {
        peakCount = std::max(peakCount, *max_element(c.indices.begin(), c.indices.end()));
    }
    ++peakCount;

    // Build conflict graph (O(n^2) but with early termination)
    std::vector<std::vector<int>> adj(n);
    for (size_t i = 0; i < n; ++i) {
        // Skip if this candidate has very low score relative to best
        if (candidates[i].score < candidates[0].score * 0.1) {
            continue;  // Likely too poor to be selected
        }

        for (size_t j = i + 1; j < n; ++j) {
            if (candidates[j].score < candidates[0].score * 0.1) {
                continue;
            }

            // Quick overlap check: if m/z ranges don't overlap, skip
            if (candidates[i].envelope->GetLastMz() < candidates[j].envelope->GetMz(0) ||
                candidates[j].envelope->GetLastMz() < candidates[i].envelope->GetMz(0)) {
                continue;
            }

            if (OverlapSignificantLight(candidates[i], candidates[j])) {
                adj[i].push_back(j);
                adj[j].push_back(i);
            }
        }
    }

    // Initial greedy selection
    std::vector<bool> selected(n, false);
    double currentScore = 0;

    for (size_t i = 0; i < n; ++i) {
        if (CanSelect(i, selected, adj)) {
            selected[i] = true;
            currentScore += candidates[i].score;
        }
    }

    // Local search optimization
    int improvements = 0;
    for (int iter = 0; iter < 500; ++iter) {
        bool madeImprovement = false;

        // Try to find a replacement opportunity
        for (size_t i = 0; i < n; ++i) {
            if (selected[i]) continue;  // Skip already selected

            // Find a conflicting candidate that has lower score
            for (int neighbor : adj[i]) {
                if (selected[neighbor] && candidates[i].score > candidates[neighbor].score) {
                    // Check if removing neighbor allows selecting i
                    selected[neighbor] = false;
                    currentScore -= candidates[neighbor].score;

                    if (CanSelect(i, selected, adj)) {
                        selected[i] = true;
                        currentScore += candidates[i].score;
                        ++improvements;
                        madeImprovement = true;
                        break;
                    } else {
                        // Revert
                        selected[neighbor] = true;
                        currentScore += candidates[neighbor].score;
                    }
                }
            }

            if (madeImprovement) break;
        }

        if (!madeImprovement) break;
    }

    // Collect selected envelopes
    for (size_t i = 0; i < n; ++i) {
        if (selected[i]) {
            envelopes.push_back(candidates[i].envelope);
        } else {
            delete candidates[i].envelope;  // Free memory
        }
    }
}

/**
 * Non-greedy envelope selection using improved local search algorithm.
 *
 * This approach:
 * 1. Finds ALL candidate envelopes without marking peaks as used
 * 2. Scores each envelope by quality metrics (intensity, spacing, length)
 * 3. Selects non-overlapping envelopes using greedy + local search optimization
 *
 * @param envelopes  Output vector of selected Envelope pointers
 * @param mz         m/z values (ascending)
 * @param intensity  Corresponding intensity values
 * @param lenMin     Minimum envelope length (inclusive)
 * @param lenMax     Maximum envelope length (inclusive)
 * @param useImprovedSelection If true, use improved selection (default: true)
 */
void FindEnvelope(std::vector<Envelope*>& envelopes,
                  const std::vector<double>& mz,
                  const std::vector<double>& intensity,
                  int lenMin = 6,
                  int lenMax = 14,
                  bool useImprovedSelection = true)
{
    envelopes.clear();
    size_t n = mz.size();
    if (n == 0 || intensity.size() != n) return;

    // Step 1: Find all candidate envelopes
    std::vector<CandidateEnvelope> candidates;
    FindAllCandidates(candidates, mz, intensity, lenMin, lenMax);

    if (candidates.empty()) return;

    // Sort indices by score (descending) for both algorithms
    std::sort(candidates.begin(), candidates.end(),
              [](const CandidateEnvelope& a, const CandidateEnvelope& b) {
                  return a.score > b.score;
              });

    if (useImprovedSelection) {
        // Use improved selection with local search
        ImprovedSelection(envelopes, candidates);
    } else {
        // Original greedy selection (for comparison)
        std::vector<bool> selected(candidates.size(), false);

        for (size_t i = 0; i < candidates.size(); ++i) {
            if (selected[i]) continue;

            // Check if this envelope overlaps with any already selected one
            bool canSelect = true;
            for (size_t j = 0; j < i; ++j) {
                if (selected[j] && OverlapSignificant(candidates[i], candidates[j])) {
                    canSelect = false;
                    break;
                }
            }

            if (canSelect) {
                selected[i] = true;
                envelopes.push_back(candidates[i].envelope);
            } else {
                delete candidates[i].envelope;
            }
        }
    }
}

int main(int argc, char* argv[]) {
    using std::cerr;
    using std::cout;
    using std::endl;
    using std::ifstream;
    using std::istringstream;
    using std::vector;
    using std::string;

    if (argc < 2) {
        cerr << "Usage: " << argv[0] << " <ms1_data.txt>\n";
        return EXIT_FAILURE;
    }

    ifstream infile(argv[1]);
    if (!infile) {
        cerr << "Error: cannot open file " << argv[1] << "\n";
        return EXIT_FAILURE;
    }

    vector<double> mzs;
    vector<double> ints;
    string line;
    int spectrum_idx = 0;
    vector<Envelope*> envelopes;

    auto process_spectrum = [&]() {
        if (mzs.empty()) return;
        ++spectrum_idx;

        // 释放上一谱的 Envelope 对象
        for (auto *e : envelopes) delete e;
        envelopes.clear();

        cout << "Spectrum " << spectrum_idx 
             << ": points=" << mzs.size() << endl;

        FindEnvelope(envelopes, mzs, ints, 6, 15);

        // 这里可以遍历 envelopes，打印或者处理每个 envelope
        for (auto *e : envelopes) {
            cout << "  Charge=" << e->GetCharge()
                 << "  Length=" << e->GetLength()
                 << "  Last m/z=" << e->GetLastMz()
                 << endl;
        }

        mzs.clear();
        ints.clear();
    };

    while (getline(infile, line)) {
        if (line.empty()) {
            continue;
        }
        // 如果是新谱开始标志，先处理上一谱
        if (line[0] == 'S') {
            process_spectrum();
            continue;
        }
        // 跳过所有以 'I' 开头的注释行
        if (line[0] == 'I') {
            continue;
        }
        // 否则应为数据行：四列数值，取前两列 mz 和 intensity
        istringstream iss(line);
        double mz, intensity;
        // 读两列，忽略后面两列
        if (iss >> mz >> intensity) {
            mzs.push_back(mz);
            ints.push_back(intensity);
        }
    }
    // 文件末尾，处理最后一谱
    process_spectrum();

    return EXIT_SUCCESS;
}

