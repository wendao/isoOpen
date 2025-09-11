import sys
import numpy as np
from dataclasses import dataclass, field
from typing import List, Tuple, Optional, Sequence
from ms_utils import ms1_spectra

@dataclass(order=True)
class IsotopeEnvelope:
    """
    A dataclass to store information about a detected isotopic envelope.
    The `order=True` makes instances sortable based on the fields in order.
    We prioritize sorting by total_intensity (desc) and spacing_error_ppm (asc).
    """
    total_intensity: float = field(init=False)
    spacing_error_ppm: float
    # --- Sort-insensitive fields ---
    charge: int = field(compare=False)
    indices: List[int] = field(compare=False)
    mz: List[float] = field(compare=False)
    intensity: List[float] = field(compare=False)
    start_mz: float = field(compare=False)
    end_mz: float = field(compare=False)
    seed_index: int = field(compare=False)

    def __post_init__(self):
        # We negate total_intensity to achieve descending sort order
        # because default sorting is ascending.
        self.total_intensity = -sum(self.intensity)

    @property
    def positive_total_intensity(self) -> float:
        """Returns the actual, non-negated total intensity."""
        return -self.total_intensity

# -----------------
# Core Logic
# -----------------

def extract_isotope_envelopes(
    mz: np.ndarray,
    intensity: np.ndarray,
    *,
    delta_c13: float = 1.0033548378,
    charges: Sequence[int] = (7, 6, 5, 4, 3, 2),
    ppm_tol: float = 10.0,
    abs_tol: float = 0.001,
    min_len: int = 5,
    max_len: int = 15,
    min_intensity: float = 0.0,
    dedup_ppm: float = 5.0,
) -> List[IsotopeEnvelope]:
    """
    Extracts evenly spaced isotopic envelopes from a single MS1 spectrum.

    The algorithm seeds envelopes from high-intensity peaks and grows them outwards,
    validates their length, splits oversized envelopes, and deduplicates the results.

    Args:
        mz: NumPy array of m/z values, assumed to be sorted ascending.
        intensity: NumPy array of corresponding intensity values.
        delta_c13: Mass difference between ¹³C and ¹²C isotopes in Daltons.
        charges: A sequence of integer charges to search for.
        ppm_tol: Relative tolerance in parts-per-million for matching peaks.
        abs_tol: Absolute m/z tolerance for matching peaks. The effective
                 tolerance is max(abs_tol, expected_mz * ppm_tol / 1e6).
        min_len: The minimum number of peaks required for a valid envelope.
        max_len: The maximum number of peaks allowed. Envelopes longer than this
                 will be split.
        min_intensity: Peaks below this absolute intensity are ignored.
        dedup_ppm: PPM tolerance used to identify duplicate envelopes. Two
                   envelopes are duplicates if they have the same charge and
                   >80% of their peaks match within this tolerance.

    Returns:
        A list of IsotopeEnvelope objects, sorted by descending total intensity.
    """
    if mz.size == 0:
        return []

    # 1. Pre-filtering and setup
    if min_intensity > 0:
        mask = intensity >= min_intensity
        mz = mz[mask]
        intensity = intensity[mask]
        original_indices = np.where(mask)[0]
    else:
        original_indices = np.arange(mz.size)
    
    if mz.size < min_len:
        return []

    seed_indices = np.argsort(intensity)[::-1]

    # 2. Grow all possible envelopes from seeds
    all_found_envelopes = []
    for seed_idx in seed_indices:
        for z in charges:
            # Grow a raw envelope from the current seed and charge
            raw_envelope_peaks = _grow_envelope(
                seed_idx, z, mz, intensity, original_indices,
                delta_c13, ppm_tol, abs_tol
            )

            if len(raw_envelope_peaks) < min_len:
                continue

            # Validate length, splitting if necessary
            validated_peak_lists = _validate_and_split(
                raw_envelope_peaks, min_len, max_len
            )

            # Score and format each valid envelope
            for peak_list in validated_peak_lists:
                env = _score_envelope(peak_list, z, delta_c13, original_indices[seed_idx])
                all_found_envelopes.append(env)
    
    if not all_found_envelopes:
        return []
        
    # 3. Deduplicate and select the best envelopes
    # Sort by quality: high intensity is best, low spacing error is second best.
    all_found_envelopes.sort() # Uses the dataclass's defined order

    final_envelopes = []
    for new_env in all_found_envelopes:
        is_duplicate = any(
            _are_envelopes_duplicate(new_env, final_env, dedup_ppm)
            for final_env in final_envelopes
        )
        if not is_duplicate:
            final_envelopes.append(new_env)

    # 4. Final sort by descending total intensity as per requirement
    final_envelopes.sort(key=lambda env: env.positive_total_intensity, reverse=True)

    return final_envelopes

# -----------------
# Helper Functions
# -----------------

def _find_closest_peak(
    target_mz: float,
    mz_array: np.ndarray,
    intensity_array: np.ndarray,
    original_indices: np.ndarray,
    ppm_tol: float,
    abs_tol: float,
) -> Optional[Tuple[int, float, float]]:
    """Finds the best matching peak within a tolerance window."""
    tolerance = max(abs_tol, target_mz * ppm_tol * 1e-6)
    lower_bound = target_mz - tolerance
    upper_bound = target_mz + tolerance

    start_idx = np.searchsorted(mz_array, lower_bound, side='left')
    end_idx = np.searchsorted(mz_array, upper_bound, side='right')

    if start_idx == end_idx:
        return None

    # Slice candidate peaks
    candidate_indices = np.arange(start_idx, end_idx)
    candidate_mz = mz_array[candidate_indices]
    
    # Find the best match: minimum m/z deviation, tie-break with max intensity
    mz_diffs = np.abs(candidate_mz - target_mz)
    min_diff = np.min(mz_diffs)
    
    # Indices of candidates with the minimum m/z difference
    tie_indices = candidate_indices[mz_diffs == min_diff]
    
    # Select the one with the highest intensity among those tied
    best_candidate_local_idx = np.argmax(intensity_array[tie_indices])
    best_candidate_idx = tie_indices[best_candidate_local_idx]

    return (
        original_indices[best_candidate_idx],
        mz_array[best_candidate_idx],
        intensity_array[best_candidate_idx],
    )

def _grow_envelope(
    seed_idx: int,
    charge: int,
    mz: np.ndarray,
    intensity: np.ndarray,
    original_indices: np.ndarray,
    delta_c13: float,
    ppm_tol: float,
    abs_tol: float,
) -> List[Tuple[int, float, float]]:
    """Grows an envelope outwards from a seed peak."""
    step = delta_c13 / charge
    
    # Initial peak
    seed_peak = (original_indices[seed_idx], mz[seed_idx], intensity[seed_idx])
    
    # Grow forward
    forward_peaks = []
    current_peak_mz = seed_peak[1]
    while True:
        next_peak = _find_closest_peak(
            current_peak_mz + step, mz, intensity, original_indices, ppm_tol, abs_tol
        )
        if next_peak:
            forward_peaks.append(next_peak)
            current_peak_mz = next_peak[1]
        else:
            break

    # Grow backward
    backward_peaks = []
    current_peak_mz = seed_peak[1]
    while True:
        prev_peak = _find_closest_peak(
            current_peak_mz - step, mz, intensity, original_indices, ppm_tol, abs_tol
        )
        if prev_peak:
            backward_peaks.append(prev_peak)
            current_peak_mz = prev_peak[1]
        else:
            break
            
    return backward_peaks[::-1] + [seed_peak] + forward_peaks

def _validate_and_split(
    envelope_peaks: List[Tuple[int, float, float]],
    min_len: int,
    max_len: int,
) -> List[List[Tuple[int, float, float]]]:
    """Recursively validates envelope length, splitting if it exceeds max_len."""
    L = len(envelope_peaks)
    if L < min_len:
        return []
    if L <= max_len:
        return [envelope_peaks]

    # Envelope is too long (L > max_len), must split at the lowest intensity point.
    # We search for the minimum in the "inner" part of the envelope to avoid trivial splits.
    intensities = np.array([p[2] for p in envelope_peaks[1:-1]])
    if intensities.size == 0: # Should not happen if L > max_len and max_len >= min_len >= 2
        return []
    
    # np.argmin finds the first occurrence in case of a tie.
    split_idx = np.argmin(intensities) + 1  # Add 1 to adjust for the [1:-1] slice

    left_sub_envelope = envelope_peaks[:split_idx]
    right_sub_envelope = envelope_peaks[split_idx + 1:]

    # Recursively validate the new sub-envelopes
    return _validate_and_split(left_sub_envelope, min_len, max_len) + \
           _validate_and_split(right_sub_envelope, min_len, max_len)

def _score_envelope(
    peak_list: List[Tuple[int, float, float]],
    charge: int,
    delta_c13: float,
    seed_original_index: int,
) -> IsotopeEnvelope:
    """Calculates metrics for a single valid envelope."""
    indices, mzs, intensities = zip(*peak_list)
    
    # Calculate spacing error
    observed_steps = np.diff(mzs)
    expected_step = delta_c13 / charge
    errors = np.abs((observed_steps - expected_step) / expected_step) * 1e6
    spacing_error_ppm = np.median(errors) if len(errors) > 0 else 0.0

    return IsotopeEnvelope(
        charge=charge,
        indices=list(indices),
        mz=list(mzs),
        intensity=list(intensities),
        spacing_error_ppm=spacing_error_ppm,
        start_mz=mzs[0],
        end_mz=mzs[-1],
        seed_index=seed_original_index
    )

def _are_envelopes_duplicate(
    env1: IsotopeEnvelope,
    env2: IsotopeEnvelope,
    dedup_ppm: float
) -> bool:
    """Checks if two envelopes are duplicates based on charge and peak similarity."""
    if env1.charge != env2.charge:
        return False

    # Identify which list of m/z values is shorter
    if len(env1.mz) <= len(env2.mz):
        shorter_mz, longer_mz = env1.mz, np.array(env2.mz)
    else:
        shorter_mz, longer_mz = env2.mz, np.array(env1.mz)
    
    if len(shorter_mz) == 0:
        return True # Two empty envelopes of same charge are duplicates

    # Count how many peaks in the shorter envelope have a match in the longer one
    common_peaks = 0
    for mz_val in shorter_mz:
        tolerance = mz_val * dedup_ppm * 1e-6
        if np.any(np.abs(longer_mz - mz_val) <= tolerance):
            common_peaks += 1
            
    # Duplicates if > 80% of peaks in the smaller envelope are shared
    overlap_fraction = common_peaks / len(shorter_mz)
    return overlap_fraction > 0.8

if __name__ == '__main__':
    for mz, inten in ms1_spectra(sys.argv[1]):
        # 在这里处理每一张 MS1 谱图
        envelopes = extract_isotope_envelopes(mz, inten)
        for envelope in envelopes:
            print("charge=", envelope.charge, "size=", len(envelope.mz))
