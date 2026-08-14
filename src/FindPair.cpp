#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <string>
#include <cstdlib>
#include <cmath>
#include <algorithm>
#include <iomanip>
#include <map>
#include <memory>
#include <limits>
#include "Envelope.h"

// ===========================================================================
// FindPair — detect LIGHT/HEAVY isotope-envelope PAIRS in one MS1 spectrum.
//
// Written in the same spirit as FindEnv.cpp, but the unit of detection is a
// *pair* of envelopes (a light-labeled peptide and its heavy-labeled partner)
// rather than a single envelope. Requiring a co-eluting partner at an exact
// label-defined m/z offset is a strong physical constraint that buys precision
// the single-envelope detector cannot get.
//
// The labelling chemistry is DELIBERATELY NOT part of this algorithm. The only
// chemistry-dependent quantity, the light->heavy mass difference per labelled
// residue (delta), is supplied at run time via --tag-delta, so one binary works
// for any reagent and a decoy run is simply a non-physical delta -- no rebuild.
// A peptide carrying nLabel labels has its heavy partner at
//   m/z shift = nLabel * delta / charge.
// Only the 13C-12C spacing (1.003354835) is hard-coded: that is physics, and it
// defines the isotope grid itself rather than the labelling scheme.
// Reference deltas: IA-TEV 6.0138091 | 6x13C 6.020129 | SILAC K8 8.014199.
//
// In the recommended forward mode, each lower-m/z light anchor is probed above
// (+n*delta) for its heavy partner. SEARCHBOTH=1 additionally probes below when
// reverse anchoring is required. Output is always oriented light (lower m/z) to
// heavy (higher m/z), with an MS1 light/heavy intensity ratio.
// ===========================================================================

#ifndef PPMTOL
#define PPMTOL 12.0     // dynamic m/z match tolerance (ppm); override -DPPMTOL=
#endif
#ifndef MAXCHARGE
#define MAXCHARGE 6     // highest charge state tried; override -DMAXCHARGE=
#endif
#ifndef LENMIN
#define LENMIN 3        // min envelope length (# isotope peaks); override -DLENMIN=
#endif
#ifndef LENMAX
#define LENMAX 15       // max envelope length before splitting; override -DLENMAX=
#endif
#ifndef HEAVYLENMIN
#define HEAVYLENMIN 3   // min length required of the HEAVY partner; override -DHEAVYLENMIN=
#endif
#ifndef MAXLABEL
#define MAXLABEL 4      // DEFAULT only; set at runtime with --max-label
#endif
#ifndef PAIRPPM
#define PAIRPPM 8.0     // tolerance on the heavy-mono position vs its expected m/z
#endif
#ifndef MAXGAP
#define MAXGAP 0        // # consecutive missing isotope peaks tolerated in growth.
                        // Kept at 0: bridging gaps added far more random pairs
                        // than real recall on R1 (frontier got strictly worse).
#endif
#ifndef MINCORR
#define MINCORR 0.6     // reject a pair whose light/heavy patterns disagree below this
#endif
#ifndef CONFLICT_OVERLAP
#define CONFLICT_OVERLAP 0.5  // two pairs conflict if peak overlap >= this
#endif
#ifndef SEARCHBOTH
#define SEARCHBOTH 0    // 0 (default/recommended): forward-only — anchor on the
                        //    light (lower-m/z) member; fewer random partners ->
                        //    higher precision. Still reports light/heavy + ratio.
                        // 1: probe partner both above (+n*delta) and below
                        //    (-n*delta) and keep the best-shape one, so any
                        //    isolated envelope self-labels light/heavy at MS1
                        //    (use when reverse ratios / heavy>light can occur).
#endif

namespace {
constexpr double kC13C12   = 1.003354835;   // 13C - 12C (Da) -- physics, not chemistry
constexpr double kPpmTol   = PPMTOL;
constexpr double kPairPpm  = PAIRPPM;
const std::vector<int> kCharges = [] {
    std::vector<int> c;
    for (int z = MAXCHARGE; z >= 1; --z) c.push_back(z);
    return c;
}();

// The labelling chemistry is NOT part of the algorithm. The light/heavy mass
// difference per labelled residue is supplied at run time (--tag-delta), so one
// binary serves any reagent (IA-TEV 6.0138, 6x13C 6.0201, SILAC, dimethyl, ...)
// and decoy runs are just a different value -- never a recompile.
double gTagDelta = 0.0;            // REQUIRED at run time; no default on purpose

// Runtime candidate gates.  Keeping these at candidate-generation time is
// important: a weak pair that will be discarded later must not occupy peaks
// and block a stronger overlapping pair during conflict selection.
double gMinCandidateCorr = MINCORR;
double gMaxCandidateMassErr = kPairPpm;
double gMinLightTotal = 0.0;
int gMaxLabel = MAXLABEL;
}  // namespace

// ---------------------------------------------------------------------------
// A candidate light/heavy pair plus what we need to rank and de-conflict it.
// ---------------------------------------------------------------------------
struct CandidatePair {
    int charge;
    int nLabel;                          // number of tag copies (Cys count)
    std::vector<size_t> lightIdx;        // light-envelope peak indices (ascending)
    std::vector<size_t> heavyIdx;        // heavy-envelope peak indices (ascending)
    std::vector<size_t> allIdx;          // union, for conflict graph
    double lightInt, heavyInt;           // summed intensities
    double lightMono, heavyMono;         // monoisotopic m/z of each member
    double massErrPpm;                   // heavy-mono deviation from expected
    double spacingPpm;                   // worse of the two spacing errors
    double patternCorr;                  // cosine similarity of the two patterns
    double score;
};

// ===========================================================================
// Peak matching / envelope growth  (same primitives as FindEnv)
// ===========================================================================

static int FindBestPeak(const std::vector<double>& mz,
                        const std::vector<double>& intensity, double target) {
    const double tol = target * kPpmTol / 1e6;
    auto lo = std::lower_bound(mz.begin(), mz.end(), target - tol);
    auto hi = std::upper_bound(mz.begin(), mz.end(), target + tol);
    int best = -1;
    double bestDiff = tol + 1.0, bestInt = -1.0;
    for (auto it = lo; it != hi; ++it) {
        size_t k = static_cast<size_t>(it - mz.begin());
        double diff = std::abs(mz[k] - target);
        if (diff < bestDiff - 1e-12 ||
            (diff <= bestDiff + 1e-12 && intensity[k] > bestInt)) {
            bestDiff = diff; bestInt = intensity[k]; best = static_cast<int>(k);
        }
    }
    return best;
}

// Grow an isotope envelope from a seed in both directions along the ideal grid.
// Up to MAXGAP consecutive missing peaks are bridged: low-abundance envelopes
// routinely drop an interior isotope, and a hard stop at the first gap
// truncates them below lenMin and loses real signal (observed on R1).
static std::vector<size_t> GrowEnvelope(const std::vector<double>& mz,
                                        const std::vector<double>& intensity,
                                        size_t seed, double deltaM) {
    std::vector<size_t> back, fwd;
    int miss = 0;
    for (double t = mz[seed] - deltaM; ; t -= deltaM) {
        int idx = FindBestPeak(mz, intensity, t);
        if (idx < 0) { if (++miss > MAXGAP) break; else continue; }
        miss = 0;
        back.push_back(static_cast<size_t>(idx));
    }
    miss = 0;
    for (double t = mz[seed] + deltaM; ; t += deltaM) {
        int idx = FindBestPeak(mz, intensity, t);
        if (idx < 0) { if (++miss > MAXGAP) break; else continue; }
        miss = 0;
        fwd.push_back(static_cast<size_t>(idx));
    }
    std::vector<size_t> idxs;
    idxs.reserve(back.size() + 1 + fwd.size());
    idxs.insert(idxs.end(), back.rbegin(), back.rend());
    idxs.push_back(seed);
    idxs.insert(idxs.end(), fwd.begin(), fwd.end());
    return idxs;
}

// Grow forward-only from an anchor (used for the heavy partner: the anchor is
// its monoisotope, so we only extend to higher m/z). Also gap-tolerant.
static std::vector<size_t> GrowForward(const std::vector<double>& mz,
                                       const std::vector<double>& intensity,
                                       size_t anchor, double deltaM) {
    std::vector<size_t> idxs;
    idxs.push_back(anchor);
    int miss = 0;
    for (double t = mz[anchor] + deltaM; ; t += deltaM) {
        int idx = FindBestPeak(mz, intensity, t);
        if (idx < 0) { if (++miss > MAXGAP) break; else continue; }
        miss = 0;
        idxs.push_back(static_cast<size_t>(idx));
    }
    return idxs;
}

static double MedianSpacingErr(const std::vector<double>& mz, int charge) {
    if (mz.size() < 2) return 0.0;
    const double expected = kC13C12 / charge;
    std::vector<double> e;
    e.reserve(mz.size() - 1);
    for (size_t i = 1; i < mz.size(); ++i) {
        double diff = mz[i] - mz[i - 1];
        // account for bridged gaps: neighbors may span >1 isotope step
        double steps = std::floor(diff / expected + 0.5);
        if (steps < 1.0) steps = 1.0;
        e.push_back(std::abs(diff - steps * expected) / (steps * expected) * 1e6);
    }
    size_t mid = e.size() / 2;
    std::nth_element(e.begin(), e.begin() + mid, e.end());
    double hi = e[mid];
    if (e.size() % 2) return hi;
    std::nth_element(e.begin(), e.begin() + mid - 1, e.end());
    return (hi + e[mid - 1]) / 2.0;
}

// Cosine similarity of the two isotope patterns over their leading L peaks.
// Real light/heavy pairs share (nearly) the same theoretical isotope envelope,
// so a high correlation is a strong true-pair signal — and it needs no
// averagine model, only that the two observed patterns look alike.
static double PatternCorr(const std::vector<double>& mz,
                          const std::vector<double>& intensity,
                          const std::vector<size_t>& a,
                          const std::vector<size_t>& b) {
    size_t L = std::min(a.size(), b.size());
    if (L < 2) return 0.0;
    double dot = 0, na = 0, nb = 0;
    for (size_t k = 0; k < L; ++k) {
        double x = intensity[a[k]], y = intensity[b[k]];
        dot += x * y; na += x * x; nb += y * y;
    }
    if (na <= 0 || nb <= 0) return 0.0;
    return dot / (std::sqrt(na) * std::sqrt(nb));
}

// ===========================================================================
// Candidate pair enumeration
// ===========================================================================

// Assemble an oriented candidate from a light (lower-m/z) and heavy (higher-m/z)
// envelope. massErr/corr are supplied by the caller.
static CandidatePair MakePair(const std::vector<double>& mz,
                              const std::vector<double>& intensity,
                              int charge, int nLabel,
                              std::vector<size_t> lightIdx,
                              std::vector<size_t> heavyIdx,
                              double massErr, double corr) {
    CandidatePair c;
    c.charge = charge;
    c.nLabel = nLabel;
    c.lightMono = mz[lightIdx.front()];
    c.heavyMono = mz[heavyIdx.front()];
    c.massErrPpm = massErr;
    c.patternCorr = corr;

    double lInt = 0, hInt = 0;
    std::vector<double> lmz, hmz;
    lmz.reserve(lightIdx.size()); hmz.reserve(heavyIdx.size());
    for (size_t idx : lightIdx) { lInt += intensity[idx]; lmz.push_back(mz[idx]); }
    for (size_t idx : heavyIdx) { hInt += intensity[idx]; hmz.push_back(mz[idx]); }
    c.lightInt = lInt; c.heavyInt = hInt;
    c.spacingPpm = std::max(MedianSpacingErr(lmz, charge),
                            MedianSpacingErr(hmz, charge));

    c.allIdx = lightIdx;
    c.allIdx.insert(c.allIdx.end(), heavyIdx.begin(), heavyIdx.end());
    std::sort(c.allIdx.begin(), c.allIdx.end());
    c.allIdx.erase(std::unique(c.allIdx.begin(), c.allIdx.end()), c.allIdx.end());

    int minLen = static_cast<int>(std::min(lightIdx.size(), heavyIdx.size()));
    c.lightIdx = std::move(lightIdx);
    c.heavyIdx = std::move(heavyIdx);
    // mass_err carries a lot of signal here: because the tag shift (6.0138) is
    // ~6x the isotope spacing, a spurious "self-pair" (light paired with its own
    // isotope tail) sits ~n*3-7 ppm off the true tag position, whereas a real
    // heavy is ~0 ppm. Penalising mass_err hard lets the true pair win selection.
    c.score = std::log(lInt + hInt + 1.0) * 2.0
            - c.spacingPpm * 0.001
            - c.massErrPpm * 0.15
            + std::log2(static_cast<double>(minLen)) * 0.5
            + corr * 2.0;
    return c;
}

static void FindAllPairs(std::vector<CandidatePair>& out,
                         const std::vector<double>& mz,
                         const std::vector<double>& intensity,
                         int lenMin, int lenMax, int maxLabel,
                         double minAnchorBase, int partnerLenMin) {
    const size_t n = mz.size();
    if (n == 0) return;

    for (int charge : kCharges) {
        const double deltaIso = kC13C12 / charge;

        for (size_t i = 0; i < n; ++i) {
            // Grow the ANCHOR envelope from the seed and re-anchor on its true
            // monoisotope. We do NOT presume the anchor is light or heavy — that
            // is decided below by which partner (above or below) fits its shape.
            std::vector<size_t> anchorAll = GrowEnvelope(mz, intensity, i, deltaIso);
            if (static_cast<int>(anchorAll.size()) < lenMin) continue;
            size_t anchorMonoIdx = anchorAll.front();
            std::vector<size_t> anchorIdx =
                GrowForward(mz, intensity, anchorMonoIdx, deltaIso);
            if (static_cast<int>(anchorIdx.size()) < lenMin) continue;
            // NB: do NOT cap the anchor to lenMax here. Its usable length is
            // capped per-nLabel to the tag gap below, so a light and its heavy
            // partner never bleed into one another.

            double anchorMono = mz[anchorMonoIdx];
            double anchorBase = 0.0;
            for (size_t idx : anchorIdx)
                if (intensity[idx] > anchorBase) anchorBase = intensity[idx];
            if (anchorBase < minAnchorBase) continue;

            for (int nLabel = 1; nLabel <= maxLabel; ++nLabel) {
                double sep = nLabel * gTagDelta / charge;
                // The light and heavy monos are ~sep/deltaIso isotope steps apart
                // (~6 per label). Cap each member's length to that gap so the two
                // envelopes stay disjoint even when they physically merge into one
                // contiguous peak run (the classic "two envelopes look like one").
                int stepGap = static_cast<int>(std::lround(sep / deltaIso));
                if (stepGap < lenMin) continue;         // partners too close to split
                int truncLen = std::min(stepGap, lenMax);
                // Two hypotheses per nLabel:
                //   partner ABOVE (+sep) -> anchor is the LIGHT member;
                //   partner BELOW (-sep) -> anchor is the HEAVY member.
                // If both exist, keep the one whose isotope pattern best matches
                // the anchor (highest pattern_corr) -> the anchor's light/heavy
                // identity is resolved here, at the MS1 level.
                double bestCorr = -1.0;
                CandidatePair best;
                bool haveBest = false;
                for (int dir = +1; dir >= -1; dir -= 2) {
                    if (!SEARCHBOTH && dir < 0) break;  // forward-only mode
                    double partnerTarget = anchorMono + dir * sep;
                    if (partnerTarget <= 0) continue;
                    int pAnchor = FindBestPeak(mz, intensity, partnerTarget);
                    if (pAnchor < 0) continue;
                    double massErr =
                        std::abs(mz[pAnchor] - partnerTarget) / partnerTarget * 1e6;
                    if (massErr > gMaxCandidateMassErr) continue;
                    std::vector<size_t> partnerIdx = GrowForward(
                        mz, intensity, static_cast<size_t>(pAnchor), deltaIso);
                    if (static_cast<int>(partnerIdx.size()) < partnerLenMin) continue;

                    // Cap both members to the tag gap (disjoint envelopes).
                    std::vector<size_t> aEnv = anchorIdx, pEnv = partnerIdx;
                    if (static_cast<int>(aEnv.size()) > truncLen) aEnv.resize(truncLen);
                    if (static_cast<int>(pEnv.size()) > truncLen) pEnv.resize(truncLen);
                    if (static_cast<int>(aEnv.size()) < lenMin) continue;
                    if (static_cast<int>(pEnv.size()) < partnerLenMin) continue;

                    double corr = PatternCorr(mz, intensity, aEnv, pEnv);
                    if (corr < gMinCandidateCorr ||
                        (haveBest && corr <= bestCorr)) continue;
                    bestCorr = corr;
                    if (dir == +1)  // anchor light, partner heavy
                        best = MakePair(mz, intensity, charge, nLabel,
                                        aEnv, pEnv, massErr, corr);
                    else            // anchor heavy, partner light
                        best = MakePair(mz, intensity, charge, nLabel,
                                        pEnv, aEnv, massErr, corr);
                    haveBest = true;
                }
                if (haveBest && best.lightInt >= gMinLightTotal)
                    out.push_back(std::move(best));
            }
        }
    }
}

// ===========================================================================
// Dedup + conflict-graph selection  (mirrors FindEnv::ImprovedSelection)
// ===========================================================================

static void DedupPairs(std::vector<CandidatePair>& cands) {
    std::map<std::tuple<int, int, std::vector<size_t>>, size_t> best;
    std::vector<CandidatePair> kept;
    kept.reserve(cands.size());
    for (auto& c : cands) {
        auto key = std::make_tuple(c.charge, c.nLabel, c.allIdx);
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
    for (int nb : adj[cand]) if (selected[nb]) return false;
    return true;
}

static void SelectPairs(std::vector<CandidatePair>& picked,
                        std::vector<CandidatePair>& cands, size_t nPeaks) {
    if (cands.empty()) return;
    std::sort(cands.begin(), cands.end(),
              [](const CandidatePair& a, const CandidatePair& b) {
                  // Resolve shared peaks from high charge to low charge. Score
                  // ranks candidates only within a charge state.
                  if (a.charge != b.charge) return a.charge > b.charge;
                  return a.score > b.score;
              });
    size_t n = cands.size();

    std::vector<std::vector<int>> peakToCands(nPeaks);
    for (size_t i = 0; i < n; ++i)
        for (size_t idx : cands[i].allIdx)
            peakToCands[idx].push_back(static_cast<int>(i));

    std::vector<std::vector<int>> adj(n);
    std::vector<int> shared(n, 0), touched;
    for (size_t i = 0; i < n; ++i) {
        touched.clear();
        for (size_t idx : cands[i].allIdx)
            for (int j : peakToCands[idx]) {
                if (static_cast<size_t>(j) <= i) continue;
                if (shared[j]++ == 0) touched.push_back(j);
            }
        size_t lenI = cands[i].allIdx.size();
        for (int j : touched) {
            double ratio = static_cast<double>(shared[j]) /
                           std::min(lenI, cands[j].allIdx.size());
            if (ratio >= CONFLICT_OVERLAP) {
                adj[i].push_back(j);
                adj[j].push_back(static_cast<int>(i));
            }
            shared[j] = 0;
        }
    }

    std::vector<bool> selected(n, false);
    for (size_t i = 0; i < n; ++i)
        if (CanSelect(i, selected, adj)) selected[i] = true;

    for (int iter = 0; iter < 500; ++iter) {
        bool improved = false;
        for (size_t i = 0; i < n && !improved; ++i) {
            if (selected[i]) continue;
            for (int nb : adj[i]) {
                // Never let a lower-charge candidate displace a selected
                // higher-charge pair.  Initial greedy selection already gives
                // higher charges first ownership of ambiguous peaks.
                if (!selected[nb] || cands[i].charge != cands[nb].charge ||
                    cands[i].score <= cands[nb].score) continue;
                selected[nb] = false;
                if (CanSelect(i, selected, adj)) {
                    selected[i] = true; improved = true; break;
                }
                selected[nb] = true;
            }
        }
        if (!improved) break;
    }

    for (size_t i = 0; i < n; ++i)
        if (selected[i]) picked.push_back(std::move(cands[i]));
}

// Public entry: detect light/heavy pairs in one spectrum.
static void FindPairs(std::vector<CandidatePair>& picked,
                      const std::vector<double>& mz,
                      const std::vector<double>& intensity,
                      int lenMin, int lenMax, int maxLabel,
                      double minLightBase) {
    picked.clear();
    if (mz.empty() || intensity.size() != mz.size()) return;
    std::vector<CandidatePair> cands;
    FindAllPairs(cands, mz, intensity, lenMin, lenMax, maxLabel, minLightBase,
                 HEAVYLENMIN);
    DedupPairs(cands);
#if defined(NOSELECT) && NOSELECT
    picked = std::move(cands);   // dump all candidates (generation ceiling)
#else
    SelectPairs(picked, cands, mz.size());
#endif
}

// ===========================================================================
// CLI: read an .ms1 file, print one CSV row per detected pair
// ===========================================================================

#ifndef MERGE
#define MERGE 1   // # consecutive MS1 scans merged before pairing (1 = per-scan).
                  // Merging neighbouring scans lets low-abundance members that
                  // flicker across scans co-occur, lifting the pairing recall
                  // ceiling (per-scan detection otherwise misses ~half the
                  // "both-visible-over-RT-but-never-in-one-scan" pairs).
#endif

struct Spectrum {
    int scan;
    double rtMin;
    std::vector<double> mz;
    std::vector<double> intensity;
};

// Merge a group of consecutive scans into one peak list: pool all peaks, sort by
// m/z, chain-cluster within kPpmTol, and represent each cluster by its
// intensity-weighted mean m/z and its MAX intensity (the per-scan apex).
static Spectrum MergeScans(std::vector<Spectrum>& grp) {
    Spectrum out;
    out.scan = grp[grp.size() / 2].scan;
    double rtSum = 0;
    for (auto& s : grp) rtSum += s.rtMin;
    out.rtMin = rtSum / grp.size();
    if (grp.size() == 1) { out.mz = std::move(grp[0].mz);
                           out.intensity = std::move(grp[0].intensity); return out; }

    std::vector<std::pair<double, double>> all;  // (mz, intensity)
    for (auto& s : grp)
        for (size_t i = 0; i < s.mz.size(); ++i)
            all.emplace_back(s.mz[i], s.intensity[i]);
    std::sort(all.begin(), all.end());

    for (size_t i = 0; i < all.size();) {
        double wSum = all[i].first * all[i].second;
        double iSum = all[i].second, maxI = all[i].second;
        size_t j = i + 1;
        while (j < all.size() &&
               all[j].first - all[j - 1].first <= all[j].first * kPpmTol / 1e6) {
            wSum += all[j].first * all[j].second;
            iSum += all[j].second;
            if (all[j].second > maxI) maxI = all[j].second;
            ++j;
        }
        out.mz.push_back(iSum > 0 ? wSum / iSum : all[i].first);
        out.intensity.push_back(maxI);
        i = j;
    }
    return out;
}

static double g_minAnchorBase = 0.0;  // drop pair if the anchor (stronger,
                                      // detectable) member's base peak < this

static bool ParseDouble(const char* text, double& value) {
    char* end = nullptr;
    value = std::strtod(text, &end);
    return end != text && end && *end == '\0' && std::isfinite(value);
}

static bool ParseInt(const char* text, int& value) {
    char* end = nullptr;
    long parsed = std::strtol(text, &end, 10);
    if (end == text || !end || *end != '\0' || parsed < 1 ||
        parsed > std::numeric_limits<int>::max()) return false;
    value = static_cast<int>(parsed);
    return true;
}

static void PrintUsage(const char* program) {
    std::cerr
        << "Usage: " << program
        << " <ms1_data.txt> --tag-delta <Da> [min_anchor_base_peak] [options]\n"
        << "Required:\n"
        << "  --tag-delta N               light->heavy mass difference per labelled\n"
        << "                              residue (Da). Property of the REAGENT, not of\n"
        << "                              the algorithm, so there is no default.\n"
        << "                              IA-TEV 6.0138091 | 6x13C 6.020129 |\n"
        << "                              SILAC K8 8.014199 | dimethyl 4x 4.025107\n"
        << "                              (a decoy run is just a non-physical value)\n"
        << "Options:\n"
        << "  --profile high-confidence  preselect corr>=0.9, massErr<=3, lightTotal>=1e6\n"
        << "  --min-anchor-base N         minimum anchor envelope base peak\n"
        << "  --min-light-total N         minimum integrated light-envelope intensity\n"
        << "  --min-corr N                candidate pattern correlation gate\n"
        << "  --max-mass-error N          candidate pair mass error in ppm\n"
        << "  --max-label N               largest tag-copy count to enumerate\n";
}

static bool ParseOptions(int argc, char* argv[]) {
    int i = 2;
    if (i < argc && std::string(argv[i]).rfind("--", 0) != 0) {
        if (!ParseDouble(argv[i], g_minAnchorBase)) return false;
        ++i;
    }
    for (; i < argc; ++i) {
        const std::string option = argv[i];
        if (option == "--help" || option == "-h") {
            PrintUsage(argv[0]);
            return false;
        }
        if (option == "--profile") {
            if (++i >= argc || std::string(argv[i]) != "high-confidence") {
                std::cerr << "Unknown or missing profile\n";
                return false;
            }
            gMinCandidateCorr = 0.90;
            gMaxCandidateMassErr = 3.0;
            gMinLightTotal = 1e6;
            continue;
        }
        if (++i >= argc) {
            std::cerr << "Missing value for " << option << '\n';
            return false;
        }
        bool ok = false;
        if (option == "--min-anchor-base")
            ok = ParseDouble(argv[i], g_minAnchorBase);
        else if (option == "--min-light-total")
            ok = ParseDouble(argv[i], gMinLightTotal);
        else if (option == "--min-corr")
            ok = ParseDouble(argv[i], gMinCandidateCorr);
        else if (option == "--max-mass-error")
            ok = ParseDouble(argv[i], gMaxCandidateMassErr);
        else if (option == "--max-label")
            ok = ParseInt(argv[i], gMaxLabel);
        else if (option == "--tag-delta")
            ok = ParseDouble(argv[i], gTagDelta);
        else {
            std::cerr << "Unknown option: " << option << '\n';
            return false;
        }
        if (!ok) {
            std::cerr << "Invalid value for " << option << ": " << argv[i] << '\n';
            return false;
        }
    }
    if (gTagDelta <= 0.0) {
        std::cerr << "Error: --tag-delta <Da> is required (light->heavy mass "
                     "difference per labelled residue). It is a property of the\n"
                     "       labelling reagent, not of the algorithm, so it has "
                     "no built-in default.\n"
                     "       e.g. IA-TEV 6.0138091 | 6x13C 6.020129 | "
                     "SILAC K8 8.014199 | dimethyl 4x 4.025107\n";
        return false;
    }
    return g_minAnchorBase >= 0.0 && gMinLightTotal >= 0.0 &&
           gMinCandidateCorr >= 0.0 && gMinCandidateCorr <= 1.0 &&
           gMaxCandidateMassErr > 0.0 && gMaxLabel >= 1;
}

static void ProcessBatch(std::vector<Spectrum>& batch) {
    if (batch.empty()) return;
    std::vector<std::string> out(batch.size());

    #pragma omp parallel for schedule(dynamic)
    for (size_t k = 0; k < batch.size(); ++k) {
        std::vector<CandidatePair> pairs;
        FindPairs(pairs, batch[k].mz, batch[k].intensity,
                  LENMIN, LENMAX, gMaxLabel, g_minAnchorBase);
        double rtSec = batch[k].rtMin * 60.0;
        std::ostringstream os;
        os << std::fixed;
        for (const auto& p : pairs) {
            double ratio = (p.heavyInt > 0) ? p.lightInt / p.heavyInt : 0.0;
            // log2(L/H); smoothed by +1 to stay finite when a member is ~0.
            double log2Ratio = std::log2((p.lightInt + 1.0) / (p.heavyInt + 1.0));
            os << batch[k].scan << ','
               << std::setprecision(4) << rtSec << ','
               << p.charge << ','
               << p.nLabel << ','
               << std::setprecision(6) << p.lightMono << ','
               << std::setprecision(6) << p.heavyMono << ','
               << p.lightIdx.size() << ','
               << p.heavyIdx.size() << ','
               << std::setprecision(2) << p.lightInt << ','
               << std::setprecision(2) << p.heavyInt << ','
               << std::setprecision(4) << ratio << ','
               << std::setprecision(4) << log2Ratio << ','
               << std::setprecision(3) << p.spacingPpm << ','
               << std::setprecision(3) << p.massErrPpm << ','
               << std::setprecision(4) << p.patternCorr << '\n';
        }
        out[k] = os.str();
    }
    for (const auto& s : out) std::cout << s;
    batch.clear();
}

int main(int argc, char* argv[]) {
    if (argc < 2) {
        PrintUsage(argv[0]);
        return EXIT_FAILURE;
    }
    if (argc == 2 && (std::string(argv[1]) == "--help" ||
                      std::string(argv[1]) == "-h")) {
        PrintUsage(argv[0]);
        return EXIT_SUCCESS;
    }
    if (!ParseOptions(argc, argv)) return EXIT_FAILURE;

    std::ifstream infile(argv[1]);
    if (!infile) {
        std::cerr << "Error: cannot open file " << argv[1] << "\n";
        return EXIT_FAILURE;
    }

    constexpr size_t kBatchSize = 256;
    std::vector<Spectrum> batch;
    std::vector<Spectrum> mergeBuf;     // holds up to MERGE raw scans
    std::vector<double> mzs, ints;
    int curScan = 0;
    double curRtMin = 0.0;

    auto flushMerge = [&]() {
        if (mergeBuf.empty()) return;
        batch.push_back(MergeScans(mergeBuf));
        mergeBuf.clear();
        if (batch.size() >= kBatchSize) ProcessBatch(batch);
    };
    auto finalize = [&]() {
        if (mzs.empty()) return;
        mergeBuf.push_back({curScan, curRtMin, std::move(mzs), std::move(ints)});
        mzs.clear(); ints.clear();
        if (mergeBuf.size() >= static_cast<size_t>(MERGE)) flushMerge();
    };

    std::cout << "scan,rt_sec,charge,n_label,light_mono_mz,heavy_mono_mz,"
                 "light_len,heavy_len,light_int,heavy_int,ratio,log2_ratio,"
                 "spacing_ppm,mass_err_ppm,pattern_corr\n";

    std::string line;
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        if (line[0] == 'S') {
            finalize();
            std::istringstream iss(line);
            std::string tag; iss >> tag >> curScan;
            curRtMin = 0.0;
            continue;
        }
        if (line[0] == 'I') {
            if (line.rfind("I\tRTime", 0) == 0 || line.rfind("I RTime", 0) == 0) {
                std::istringstream iss(line);
                std::string a, b; double v;
                if (iss >> a >> b >> v) curRtMin = v;
            }
            continue;
        }
        if (line[0] == 'H') continue;
        std::istringstream iss(line);
        double m, in;
        if (iss >> m >> in) { mzs.push_back(m); ints.push_back(in); }
    }
    finalize();
    flushMerge();       // merge any trailing scans (< MERGE)
    ProcessBatch(batch);
    return EXIT_SUCCESS;
}
