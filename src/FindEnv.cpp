#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <limits>
#include <map>
#include <memory>
#include "Envelope.h"

// ---------------------------------------------------------------------------
// Constants
// ---------------------------------------------------------------------------
#ifndef PPMTOL
#define PPMTOL 12.0   // dynamic m/z match tolerance (ppm); override with -DPPMTOL=
#endif
#ifndef MAXCHARGE
#define MAXCHARGE 6   // highest charge state tried; override with -DMAXCHARGE=
#endif
#ifndef LENMIN
#define LENMIN 6      // minimum envelope length (# isotope peaks); override with -DLENMIN=
#endif
#ifndef LENMAX
#define LENMAX 15     // maximum envelope length before splitting; override with -DLENMAX=
#endif
namespace {
constexpr double kC13C12 = 1.003354835;  // 13C - 12C mass difference (Da)
double g_matchPpm = PPMTOL;
const std::vector<int> kCharges = [] {
    std::vector<int> c;
    for (int z = MAXCHARGE; z >= 1; --z) c.push_back(z);
    return c;
}();
}  // namespace

/**
 * Candidate envelope plus the data needed to rank and de-conflict it.
 */
struct CandidateEnvelope {
    std::unique_ptr<Envelope> envelope;  // the envelope object
    std::vector<size_t> indices;         // peak indices it occupies (ascending)
    double score;                        // ranking score
};

// ===========================================================================
// Peak matching
// ===========================================================================

/**
 * Find the best peak near `target` within a ppm-scaled window.
 * "Best" = smallest |m/z - target|, ties broken by higher intensity.
 * Assumes `mz` is sorted ascending. Returns -1 if no peak matches.
 */
static int FindBestPeak(const std::vector<double>& mz,
                        const std::vector<double>& intensity,
                        double target) {
    const double tol = target * g_matchPpm / 1e6;
    auto lo = std::lower_bound(mz.begin(), mz.end(), target - tol);
    auto hi = std::upper_bound(mz.begin(), mz.end(), target + tol);

    int best = -1;
    double bestDiff = tol + 1.0;
    double bestInt = -1.0;
    for (auto it = lo; it != hi; ++it) {
        size_t k = static_cast<size_t>(it - mz.begin());
        double diff = std::abs(mz[k] - target);
        if (diff < bestDiff - 1e-12 ||
            (diff <= bestDiff + 1e-12 && intensity[k] > bestInt)) {
            bestDiff = diff;
            bestInt = intensity[k];
            best = static_cast<int>(k);
        }
    }
    return best;
}

/**
 * Grow an isotope envelope from a seed peak in both directions.
 * Targets march along the ideal isotope grid anchored on the seed.
 * Returns peak indices in ascending m/z order (includes the seed).
 */
static std::vector<size_t> GrowEnvelope(const std::vector<double>& mz,
                                        const std::vector<double>& intensity,
                                        size_t seed, double deltaM) {
    std::vector<size_t> back;   // lower-m/z peaks, collected high -> low
    std::vector<size_t> fwd;    // higher-m/z peaks, collected low -> high

    for (double t = mz[seed] - deltaM; ; t -= deltaM) {
        int idx = FindBestPeak(mz, intensity, t);
        if (idx < 0) break;
        back.push_back(static_cast<size_t>(idx));
    }
    for (double t = mz[seed] + deltaM; ; t += deltaM) {
        int idx = FindBestPeak(mz, intensity, t);
        if (idx < 0) break;
        fwd.push_back(static_cast<size_t>(idx));
    }

    std::vector<size_t> idxs;
    idxs.reserve(back.size() + 1 + fwd.size());
    idxs.insert(idxs.end(), back.rbegin(), back.rend());  // reverse -> ascending
    idxs.push_back(seed);
    idxs.insert(idxs.end(), fwd.begin(), fwd.end());
    return idxs;
}

// ===========================================================================
// Scoring helpers
// ===========================================================================

/**
 * Median ppm deviation of the peak spacing from the expected deltaM.
 * Median is more robust to outliers than the mean.
 */
static double CalculateSpacingError(const std::vector<double>& mz, int charge) {
    if (mz.size() < 2) return 0.0;
    const double expectedDelta = kC13C12 / charge;

    std::vector<double> errors;
    errors.reserve(mz.size() - 1);
    for (size_t i = 1; i < mz.size(); ++i) {
        double err = std::abs((mz[i] - mz[i - 1]) - expectedDelta) / expectedDelta * 1e6;
        errors.push_back(err);
    }

    size_t mid = errors.size() / 2;
    std::nth_element(errors.begin(), errors.begin() + mid, errors.end());
    double hiVal = errors[mid];
    if (errors.size() % 2 != 0) return hiVal;
    // even count: average the two central order statistics
    std::nth_element(errors.begin(), errors.begin() + mid - 1, errors.end());
    return (hiVal + errors[mid - 1]) / 2.0;
}

/**
 * Combined ranking score: rewards intensity and length, penalizes spacing error.
 */
static double ScoreEnvelope(const std::vector<double>& mz,
                            const std::vector<double>& intensity,
                            const std::vector<size_t>& idxs, int charge) {
    double totalIntensity = 0.0;
    std::vector<double> mzs;
    mzs.reserve(idxs.size());
    for (size_t idx : idxs) {
        totalIntensity += intensity[idx];
        mzs.push_back(mz[idx]);
    }
    double spacingErrorPPM = CalculateSpacingError(mzs, charge);
    return std::log(totalIntensity + 1.0) * 2.0
         - spacingErrorPPM * 0.001
         + std::log2(static_cast<double>(idxs.size())) * 0.5;
}

/**
 * Recursively split an oversized envelope at its lowest-intensity interior peak.
 *  - length < lenMin            -> discard
 *  - lenMin <= length <= lenMax -> keep
 *  - length > lenMax            -> split at min-intensity interior peak, recurse
 */
static void SplitAndValidate(const std::vector<size_t>& idxs,
                             const std::vector<double>& intensity,
                             int lenMin, int lenMax,
                             std::vector<std::vector<size_t>>& result) {
    int len = static_cast<int>(idxs.size());
    if (len < lenMin) return;
    if (len <= lenMax) {
        result.push_back(idxs);
        return;
    }

    int splitIdx = 1;
    double minInt = intensity[idxs[1]];
    for (int k = 2; k < len - 1; ++k) {
        if (intensity[idxs[k]] < minInt) {
            minInt = intensity[idxs[k]];
            splitIdx = k;
        }
    }

    std::vector<size_t> left(idxs.begin(), idxs.begin() + splitIdx);
    std::vector<size_t> right(idxs.begin() + splitIdx + 1, idxs.end());
    SplitAndValidate(left, intensity, lenMin, lenMax, result);
    SplitAndValidate(right, intensity, lenMin, lenMax, result);
}

// ===========================================================================
// Candidate enumeration
// ===========================================================================

/**
 * Enumerate every plausible envelope (no greedy peak assignment) and score it.
 * Optionally pre-filters out peaks below `minIntensity`.
 */
static void FindAllCandidates(std::vector<CandidateEnvelope>& candidates,
                              const std::vector<double>& mz,
                              const std::vector<double>& intensity,
                              int lenMin, int lenMax,
                              double minIntensity,
                              double minTotalIntensity,
                              double maxSpacingPpm) {
    // Optional intensity pre-filter. `work*` are the arrays we search in;
    // `origIdx` maps a working index back to the original peak index.
    std::vector<double> fMz, fInt;
    std::vector<size_t> origIdx;
    const bool useFilteredPeaks = minIntensity > 0.0;
    if (useFilteredPeaks) {
        for (size_t k = 0; k < mz.size(); ++k) {
            if (intensity[k] >= minIntensity) {
                fMz.push_back(mz[k]);
                fInt.push_back(intensity[k]);
                origIdx.push_back(k);
            }
        }
    }
    const std::vector<double>& workMz = useFilteredPeaks ? fMz : mz;
    const std::vector<double>& workInt = useFilteredPeaks ? fInt : intensity;

    size_t n = workMz.size();
    if (n == 0) return;

    for (int charge : kCharges) {
        double deltaM = kC13C12 / charge;

        for (size_t i = 0; i < n; ++i) {
            std::vector<size_t> idxs = GrowEnvelope(workMz, workInt, i, deltaM);
            if (static_cast<int>(idxs.size()) < lenMin) continue;

            // Map working indices back to original peak indices.
            if (useFilteredPeaks) {
                for (auto& idx : idxs) idx = origIdx[idx];
            }

            // Oversized envelopes are split into valid sub-envelopes.
            std::vector<std::vector<size_t>> subEnvelopes;
            SplitAndValidate(idxs, intensity, lenMin, lenMax, subEnvelopes);

            for (auto& subIdx : subEnvelopes) {
                double totalIntensity = 0.0;
                std::vector<double> subMzs;
                subMzs.reserve(subIdx.size());
                for (size_t idx : subIdx) {
                    totalIntensity += intensity[idx];
                    subMzs.push_back(mz[idx]);
                }
                if (totalIntensity < minTotalIntensity ||
                    CalculateSpacingError(subMzs, charge) > maxSpacingPpm) {
                    continue;
                }

                CandidateEnvelope cand;
                cand.score = ScoreEnvelope(mz, intensity, subIdx, charge);
                cand.envelope = std::make_unique<Envelope>(charge);
                cand.envelope->Reserve(subIdx.size());
                for (size_t idx : subIdx) {
                    cand.envelope->AddPeak(mz[idx], intensity[idx]);
                }
                cand.indices = std::move(subIdx);
                candidates.push_back(std::move(cand));
            }
        }
    }
}

// ===========================================================================
// Conflict-graph selection
// ===========================================================================

/**
 * Drop exact-duplicate candidates (same charge + same peak set), keeping the
 * highest-scoring one. Every peak of a real envelope acts as a seed, so the
 * same envelope is enumerated many times; removing duplicates shrinks the set
 * fed into the conflict graph without changing the selection outcome
 * (duplicates fully overlap and only the best would ever be chosen).
 */
static void DedupCandidates(std::vector<CandidateEnvelope>& cands) {
    std::map<std::pair<int, std::vector<size_t>>, size_t> best;  // key -> kept index
    std::vector<CandidateEnvelope> kept;
    kept.reserve(cands.size());
    for (auto& c : cands) {
        auto key = std::make_pair(c.envelope->GetCharge(), c.indices);
        auto it = best.find(key);
        if (it == best.end()) {
            best.emplace(std::move(key), kept.size());
            kept.push_back(std::move(c));
        } else if (c.score > kept[it->second].score) {
            kept[it->second] = std::move(c);
        }
    }
    cands = std::move(kept);
}

static bool CanSelect(size_t cand, const std::vector<bool>& selected,
                      const std::vector<std::vector<int>>& adj) {
    for (int neighbor : adj[cand]) {
        if (selected[neighbor]) return false;
    }
    return true;
}

/**
 * Select a non-overlapping subset of candidates that maximizes total score:
 *  1. sort by score descending
 *  2. build a conflict graph (peak overlap >= 50%)
 *  3. greedy initial pick
 *  4. local search: swap in a higher-scoring candidate for a conflicting lower one
 */
static void ImprovedSelection(std::vector<std::unique_ptr<Envelope>>& envelopes,
                              std::vector<CandidateEnvelope>& candidates,
                              size_t nPeaks) {
    if (candidates.empty()) return;

    std::sort(candidates.begin(), candidates.end(),
              [](const CandidateEnvelope& a, const CandidateEnvelope& b) {
                  return a.score > b.score;
              });

    size_t n = candidates.size();

    // Conflict graph via an inverted index. Two candidates can only conflict if
    // they share a peak, so we list candidates per peak and only examine pairs
    // that actually co-occur, instead of comparing all O(n^2) pairs.
    std::vector<std::vector<int>> peakToCands(nPeaks);
    for (size_t i = 0; i < n; ++i) {
        for (size_t idx : candidates[i].indices) {
            peakToCands[idx].push_back(static_cast<int>(i));
        }
    }

    std::vector<std::vector<int>> adj(n);
    std::vector<int> shared(n, 0);   // # peaks shared with current candidate i
    std::vector<int> touched;        // which j were touched (to reset cheaply)
    for (size_t i = 0; i < n; ++i) {
        touched.clear();
        for (size_t idx : candidates[i].indices) {
            for (int j : peakToCands[idx]) {
                if (static_cast<size_t>(j) <= i) continue;  // add each edge once
                if (shared[j]++ == 0) touched.push_back(j);
            }
        }
        size_t lenI = candidates[i].indices.size();
        for (int j : touched) {
            double ratio = static_cast<double>(shared[j]) /
                           std::min(lenI, candidates[j].indices.size());
            if (ratio >= 0.5) {
                adj[i].push_back(j);
                adj[j].push_back(static_cast<int>(i));
            }
            shared[j] = 0;  // reset for next i
        }
    }

    // Greedy initial selection (candidates already sorted by score).
    std::vector<bool> selected(n, false);
    for (size_t i = 0; i < n; ++i) {
        if (CanSelect(i, selected, adj)) selected[i] = true;
    }

    // Local search: replace a selected candidate with a higher-scoring conflict.
    for (int iter = 0; iter < 500; ++iter) {
        bool improved = false;
        for (size_t i = 0; i < n && !improved; ++i) {
            if (selected[i]) continue;
            for (int neighbor : adj[i]) {
                if (!selected[neighbor] ||
                    candidates[i].score <= candidates[neighbor].score) {
                    continue;
                }
                selected[neighbor] = false;
                if (CanSelect(i, selected, adj)) {
                    selected[i] = true;
                    improved = true;
                    break;
                }
                selected[neighbor] = true;  // revert
            }
        }
        if (!improved) break;
    }

    for (size_t i = 0; i < n; ++i) {
        if (selected[i]) envelopes.push_back(std::move(candidates[i].envelope));
    }
}

// ===========================================================================
// Public entry point
// ===========================================================================

/**
 * Detect isotope envelopes in one spectrum.
 * @param mz/intensity  parallel peak arrays, m/z ascending
 * @param lenMin/lenMax min / max envelope length
 * @param minIntensity  optional peak pre-filter (0 = keep all)
 */
void FindEnvelope(std::vector<std::unique_ptr<Envelope>>& envelopes,
                  const std::vector<double>& mz,
                  const std::vector<double>& intensity,
                  int lenMin = 6, int lenMax = 15, double minIntensity = 0.0,
                  double minTotalIntensity = 0.0,
                  double maxSpacingPpm = std::numeric_limits<double>::infinity()) {
    envelopes.clear();
    if (mz.empty() || intensity.size() != mz.size()) return;

    std::vector<CandidateEnvelope> candidates;
    FindAllCandidates(candidates, mz, intensity, lenMin, lenMax, minIntensity,
                      minTotalIntensity, maxSpacingPpm);
    DedupCandidates(candidates);
    ImprovedSelection(envelopes, candidates, mz.size());
}

// ===========================================================================
// CLI: read an .ms1 file and print detected envelopes per spectrum
// ===========================================================================

// One MS1 spectrum: scan number, retention time (minutes, as stored in .ms1),
// and the peak lists.
struct Spectrum {
    int scan;
    double rtMin;
    std::vector<double> mz;
    std::vector<double> intensity;
};

static double g_basePeakMin = 0.0;  // drop envelope if its base (max) peak < this
static double g_peakFloor   = 0.0;  // ignore peaks below this when building envelopes
static double g_totalIntensityMin = 0.0;  // reject weak candidates before selection
static double g_maxSpacingPpm = std::numeric_limits<double>::infinity();
static int g_lenMin = LENMIN;
static int g_lenMax = LENMAX;

// Detect envelopes in a batch of spectra in parallel, then emit one CSV row per
// detected envelope in file order. Spectra are independent -> scales on cores.
// Columns: scan,rt_sec,charge,mono_mz,length,base_peak,intensity,spacing_ppm
static void ProcessBatch(std::vector<Spectrum>& batch) {
    if (batch.empty()) return;
    std::vector<std::string> out(batch.size());

    #pragma omp parallel for schedule(dynamic)
    for (size_t k = 0; k < batch.size(); ++k) {
        std::vector<std::unique_ptr<Envelope>> envelopes;
        FindEnvelope(envelopes, batch[k].mz, batch[k].intensity, g_lenMin, g_lenMax,
                     g_peakFloor, g_totalIntensityMin, g_maxSpacingPpm);

        double rtSec = batch[k].rtMin * 60.0;  // GT files use seconds
        std::ostringstream os;
        os << std::fixed;
        for (const auto& e : envelopes) {
            double totalIntensity = 0.0, basePeak = 0.0;
            for (std::size_t p = 0; p < e->GetLength(); ++p) {
                double v = e->GetIntensity(p);
                totalIntensity += v;
                if (v > basePeak) basePeak = v;   // base (most intense) peak
            }
            if (basePeak < g_basePeakMin) continue;
            double spacingPpm = CalculateSpacingError(e->GetMzs(), e->GetCharge());
            os << batch[k].scan << ','
               << std::setprecision(4) << rtSec << ','
               << e->GetCharge() << ','
               << std::setprecision(6) << e->GetMz(0) << ','   // monoisotopic (lowest m/z)
               << e->GetLength() << ','
               << std::setprecision(2) << basePeak << ','
               << std::setprecision(2) << totalIntensity << ','
               << std::setprecision(3) << spacingPpm << '\n';
        }
        out[k] = os.str();
    }

    for (const auto& s : out) std::cout << s;
    batch.clear();
}

static void PrintUsage(const char* program) {
    std::cerr
        << "Usage: " << program << " <ms1_data.txt> [options]\n"
        << "\n"
        << "Options:\n"
        << "  --profile ups2-gt       validated balanced GT profile\n"
        << "  --min-base-peak N       output filter on strongest peak (default 0)\n"
        << "  --peak-floor N          ignore input peaks below N (default 0)\n"
        << "  --min-total-intensity N reject weak candidates before selection (default 0)\n"
        << "  --max-spacing-ppm N     reject irregular candidates before selection (default inf)\n"
        << "  --ppm N                 isotope peak-match tolerance (default " << PPMTOL << ")\n"
        << "  --len-min N             minimum isotope peaks (default " << LENMIN << ")\n"
        << "  --len-max N             split candidates longer than N (default " << LENMAX << ")\n"
        << "\n"
        << "Legacy positional form is still accepted:\n"
        << "  " << program << " file [min_base_peak] [peak_floor]"
           " [min_total_intensity] [max_spacing_ppm]\n";
}

static bool ParseOptions(int argc, char* argv[]) {
    int legacyPosition = 0;
    for (int i = 2; i < argc; ++i) {
        const std::string arg = argv[i];
        if (arg == "--help" || arg == "-h") {
            PrintUsage(argv[0]);
            return false;
        }

        auto value = [&](const char* option) -> const char* {
            if (++i >= argc) {
                std::cerr << "Error: missing value for " << option << "\n";
                return nullptr;
            }
            return argv[i];
        };

        if (arg == "--profile") {
            const char* name = value("--profile");
            if (!name) return false;
            if (std::string(name) != "ups2-gt") {
                std::cerr << "Error: unknown profile '" << name << "'\n";
                return false;
            }
            // Full UPS2 validation: +4.5 points GT recall while removing 59%
            // of unmatched predicted features relative to the generic baseline.
            g_matchPpm = 6.0;
            g_lenMin = 3;
            g_totalIntensityMin = 3500.0;
            g_maxSpacingPpm = 15000.0;
        } else if (arg == "--min-base-peak") {
            const char* v = value("--min-base-peak"); if (!v) return false;
            g_basePeakMin = std::atof(v);
        } else if (arg == "--peak-floor") {
            const char* v = value("--peak-floor"); if (!v) return false;
            g_peakFloor = std::atof(v);
        } else if (arg == "--min-total-intensity") {
            const char* v = value("--min-total-intensity"); if (!v) return false;
            g_totalIntensityMin = std::atof(v);
        } else if (arg == "--max-spacing-ppm") {
            const char* v = value("--max-spacing-ppm"); if (!v) return false;
            g_maxSpacingPpm = std::atof(v);
        } else if (arg == "--ppm") {
            const char* v = value("--ppm"); if (!v) return false;
            g_matchPpm = std::atof(v);
        } else if (arg == "--len-min") {
            const char* v = value("--len-min"); if (!v) return false;
            g_lenMin = std::atoi(v);
        } else if (arg == "--len-max") {
            const char* v = value("--len-max"); if (!v) return false;
            g_lenMax = std::atoi(v);
        } else if (!arg.empty() && arg[0] == '-') {
            std::cerr << "Error: unknown option '" << arg << "'\n";
            return false;
        } else {
            const double v = std::atof(arg.c_str());
            switch (legacyPosition++) {
                case 0: g_basePeakMin = v; break;
                case 1: g_peakFloor = v; break;
                case 2: g_totalIntensityMin = v; break;
                case 3: g_maxSpacingPpm = v; break;
                default:
                    std::cerr << "Error: too many positional arguments\n";
                    return false;
            }
        }
    }

    if (g_matchPpm <= 0.0 || g_lenMin < 2 || g_lenMax < g_lenMin ||
        g_basePeakMin < 0.0 || g_peakFloor < 0.0 || g_totalIntensityMin < 0.0 ||
        g_maxSpacingPpm <= 0.0) {
        std::cerr << "Error: invalid detector parameters\n";
        return false;
    }
    return true;
}

int main(int argc, char* argv[]) {
    if (argc < 2 || (argc == 2 && (std::string(argv[1]) == "--help" ||
                                   std::string(argv[1]) == "-h"))) {
        PrintUsage(argv[0]);
        return argc < 2 ? EXIT_FAILURE : EXIT_SUCCESS;
    }
    if (!ParseOptions(argc, argv)) return EXIT_FAILURE;

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Error: cannot open file " << argv[1] << "\n";
        return EXIT_FAILURE;
    }

    constexpr size_t kBatchSize = 256;  // spectra processed per parallel batch
    std::vector<Spectrum> batch;
    std::vector<double> mzs, ints;
    int curScan = 0;
    double curRtMin = 0.0;

    auto finalizeSpectrum = [&]() {
        if (mzs.empty()) return;
        batch.push_back({curScan, curRtMin, std::move(mzs), std::move(ints)});
        mzs.clear();  // moved-from vectors: restore to known empty state
        ints.clear();
        if (batch.size() >= kBatchSize) ProcessBatch(batch);
    };

    std::cout << "scan,rt_sec,charge,mono_mz,length,base_peak,intensity,spacing_ppm\n";

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == 'S') {                 // new spectrum
            finalizeSpectrum();
            std::istringstream iss(line);
            std::string tag;
            iss >> tag >> curScan;            // "S <scan> <scan_id>"
            curRtMin = 0.0;
            continue;
        }
        if (line[0] == 'I') {                 // info line; capture RTime if present
            if (line.rfind("I\tRTime", 0) == 0 || line.rfind("I RTime", 0) == 0) {
                std::istringstream iss(line);
                std::string a, b;
                double v;
                if (iss >> a >> b >> v) curRtMin = v;
            }
            continue;
        }
        if (line[0] == 'H') continue;         // header line

        std::istringstream iss(line);
        double mz, intensity;
        if (iss >> mz >> intensity) {  // first two columns only
            mzs.push_back(mz);
            ints.push_back(intensity);
        }
    }
    finalizeSpectrum();  // final spectrum
    ProcessBatch(batch);  // flush remainder

    return EXIT_SUCCESS;
}
