# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

isoOpen (eyes open) is a mass spectrometry data analysis toolkit for proteomics research. It extracts features from MS1 spectra, transfers identification results with target modifications (e.g., outputs from MSFragger or MSGF+), and generates inclusion lists for targeted quantification.

## Build Commands

### C++ Components

From the `src/` directory:

```bash
# Build envelope detection tools
g++ -std=c++17 -O2 -Wall -c Envelope.cpp -o Envelope.o
g++ -std=c++17 -O2 -Wall FindEnv.cpp Envelope.o -o find_envelope

# Build SESTAR (main MS1 processing executable)
# Uncomment in comp.sh:
# g++ -O2 -o SESTAR SESTAR.cpp
```

### Python Dependencies

```bash
pip install numpy scipy matplotlib pyteomics
```

## Architecture

### C++ and Python Toolkit

The toolkit consists of C++ and Python utilities for proteomics workflows.

**C++ Layer (`src/`):**
- `SESTAR.cpp` - Main MS1 spectrum processor for Se (GGJ)
- `Envelope.cpp/.h` - Core isotope envelope representation class
- `FindEnv.cpp` - Envelope detection algorithm
- Compiled binaries: `SESTAR`, `find_envelope`

**Python Layer (`scripts/`):**
- `ms1_isotope_extractor.py` - Core isotope envelope extraction with `IsotopeEnvelope` dataclass
- `detect_mod_shift.py` / `detect_mod_shift_Ascore.py` - Modification detection
- `find_best_mods.py` - Optimal modification identification
- `test_XL.py` - Cross-linking analysis
- `ms_utils.py` - MS1 spectra utilities (imported by extractors)

### Data Flow

1. Raw MS1 spectra (`.ms1` format) processed by C++ tools
2. Python scripts perform isotope envelope extraction using `extract_isotope_envelopes()`
3. Envelopes validated for charge state, spacing (deltaC12C13 = 1.0033548), length
4. Modifications detected via mass shift analysis
5. Results mapped to search engine identifications

### Key Parameters

From `scripts/params/peptide.py`:
- `AA_MonoMassMap` - Monoisotopic amino acid masses
- `deltaC12C13` = 1.0033548 - ¹³C/¹²C mass difference
- `MonoHplus` = 1.0072765 - Proton mass
- `AA_NumCarbonMap` - Carbon count per amino acid (for isotope calculations)

### Isotope Envelope Algorithm

The Python envelope extractor (`ms1_isotope_extractor.py`) works by:
1. Seeding envelopes from high-intensity peaks (sorted descending)
2. Growing envelopes outwards for each charge state (default: charges=[7,6,5,4,3,2])
3. Validating length (min_len=5, max_len=15) and spacing (ppm_tol=10.0)
4. Splitting oversized envelopes
5. Deduplicating (dedup_ppm=5.0, >80% peak overlap)

## Data Formats

- `.ms1` - MS1 spectrum format (H, S, I lines for headers, scans, intensities)
- `scripts/` expects mzML/mzXML files (via pyteomics)
- Output: Envelopes with m/z, intensity, charge, indices

## Directory Structure

- `src/` - C++ source, compiled binaries, MS1 data files
- `scripts/` - Python analysis scripts and `params/` subdirectory
- `scripts/params/` - Atomic and peptide mass constants
- `reference_data/` - Reference datasets and processed results
- `benchmarks/` - Performance benchmarks (feature_detection/, paired_features/)
- `docs/` - Logo and demo images
