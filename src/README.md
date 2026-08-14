# MS1 Isotope Envelope and Light/Heavy Pair Detection

This repository provides two C++17 command-line detectors for centroided MS1
spectra in `.ms1` text format:

- **FindEnv** (`FindEnv.cpp`) detects isotope envelopes in each MS1 scan.
- **FindPair** (`FindPair.cpp`) detects a light/heavy pair of isotope envelopes
  with a known mass shift per labelled residue. It reports the MS1 light/heavy
  intensity ratio for every selected pair.

Both programs use mass-accuracy-aware peak matching, non-greedy candidate
selection, and OpenMP parallelism across independent spectra.

## Build

Requirements: a C++17 compiler and OpenMP. On macOS, install Homebrew `libomp`.
The build script detects macOS versus Linux automatically.

```bash
git clone <repository-url>
cd src
bash comp.sh
```

This creates `find_envelope` and `find_pair`. Set the worker count when running:

```bash
OMP_NUM_THREADS=8 ./find_envelope --help
OMP_NUM_THREADS=8 ./find_pair --help
```

## Input format

Both tools read an MS1 text file containing spectra such as:

```text
S 1 1
I	RTime	0.00415
653.1970 462398.4
654.2004 104532.7
...
S 2 2
```

- `S` begins a spectrum; the second field is emitted as `scan`.
- `I RTime` is optional and is interpreted as minutes.
- `H` and other `I` lines are ignored.
- Peak rows need at least `m/z intensity`; subsequent columns are ignored.
- Input peaks must be sorted by increasing m/z within each spectrum.

## FindEnv: isotope-envelope detection

### Run

```bash
# Generic defaults: 12 ppm, charges 1--6, at least 6 isotope peaks.
OMP_NUM_THREADS=8 ./find_envelope data.ms1 > envelopes.csv

# UPS2-specific validated operating point (not a universal instrument profile).
OMP_NUM_THREADS=8 ./find_envelope data.ms1 --profile ups2-gt > envelopes_ups2.csv

# Inspect or tune parameters.
./find_envelope --help
./find_envelope data.ms1 --ppm 8 --len-min 3 \
  --min-total-intensity 5000 --max-spacing-ppm 15000 > envelopes.csv
```

The output is one envelope per selected spectrum-level detection:

```text
scan,rt_sec,charge,mono_mz,length,base_peak,intensity,spacing_ppm
```

`mono_mz` is the lowest m/z peak in the detected isotope series; `intensity`
is the sum of envelope intensities, and `base_peak` is its largest peak.

The generic defaults aim to be conservative and portable. The `ups2-gt` profile
uses 6 ppm, `len-min=3`, total intensity >= 3500, and median spacing error <=
15000 ppm. Absolute intensity thresholds must be re-tuned for another
instrument or acquisition scale.

## FindPair: light/heavy envelope-pair detection

`FindPair` jointly detects a lower-m/z light envelope and its heavier partner.
The labelling chemistry is supplied at runtime, rather than being hard-coded.
The required `--tag-delta` value is the light-to-heavy mass difference in Da
**per labelled residue**.

```bash
# Generic example: IA-TEV, delta = 6.0138091 Da per labelled Cys.
OMP_NUM_THREADS=8 ./find_pair data.ms1 \
  --tag-delta 6.0138091 > pairs.csv

# Examples of common deltas: 6x13C = 6.020129; SILAC K8 = 8.014199.
# Limit the number of labelled residues searched if it is known.
OMP_NUM_THREADS=8 ./find_pair data.ms1 \
  --tag-delta 6.020129 --max-label 1 > pairs.csv

# Apply validated high-confidence candidate gates before conflict selection.
OMP_NUM_THREADS=8 ./find_pair data.ms1 \
  --tag-delta 6.0138091 --profile high-confidence > pairs_hc.csv
```

`--profile high-confidence` is equivalent to correlation >= 0.90, pair mass
error <= 3 ppm, and light-envelope total intensity >= 1e6. These values are
data-dependent; inspect the output and re-tune them for a new experiment.

The output is one selected pair per processed spectrum:

```text
scan,rt_sec,charge,n_label,light_mono_mz,heavy_mono_mz,light_len,heavy_len,
light_int,heavy_int,ratio,log2_ratio,spacing_ppm,mass_err_ppm,pattern_corr
```

`n_label` is the inferred label-copy count. `ratio` is `light_int / heavy_int`;
the reported orientation is always light (lower m/z) to heavy (higher m/z).

## Quick validation results

These are reference evaluations, not guarantees for another instrument,
labelling reagent, or concentration range.

### FindEnv: UPS2 MS1 ground truth

Dataset: `18185_REP2_4pmol_UPS2_IDA_1` (6,616 MS1 spectra), evaluated against
7,234 MS1 clusters and 3,412 PSMs. Cluster and PSM matching use 10 ppm.

| Configuration | Detections | GT cluster recall | Detection precision | Correct-charge PSM recall |
|---|---:|---:|---:|---:|
| Generic default | 1,359,840 | 88.6% | 4.75% | 95.8% |
| `--profile ups2-gt` | 748,186 | **93.1%** | **12.75%** | 95.1% |

Reproduce the score when the UPS2 input and ground-truth files are available:

```bash
OMP_NUM_THREADS=8 ./find_envelope 18185_REP2_4pmol_UPS2_IDA_1.ms1 \
  --profile ups2-gt > detections_ups2.csv
python3 GT/eval_feature_finding.py detections_ups2.csv \
  GT/charged_pos_std.csv GT/psm.tsv
```

### FindPair: YJ7 light/heavy ratio series

For the Q Exactive Plus YJ7 single-label, 6x13C dataset, `FindPair` was run
directly per spectrum with `--tag-delta 6.020129 --max-label 1`. The table
uses paired MS2 identifications as ground truth; precision is estimated by a
non-physical mass-shift decoy run.

| Result across L10H1--L1H10 | Value |
|---|---:|
| Pair recall | 98.0%--100.0% |
| Target-decoy precision | 95.2%--99.1% |
| Absolute median log2-ratio error | 0.10--0.28 |

This is a single-instrument, single-label validation. Re-tune `--tag-delta`,
`--max-label`, and any confidence gates for another labelling chemistry or
instrument; it is not a chemistry- or instrument-agnostic FDR guarantee.

## Implementation notes

- **FindEnv** searches charges 6 to 1 by default, grows candidates on the
  isotope grid using ppm-scaled matching, and resolves shared peaks with a
  conflict graph.
- **FindPair** adds an exact label-mass-offset and isotope-pattern-correlation
  constraint before the same type of conflict selection.
- Both outputs are per-spectrum detections. If LC-level entities are required,
  perform an explicit RT-aware grouping/deduplication step.

For algorithmic detail, see [`Envelope_Algorithm.md`](Envelope_Algorithm.md).
