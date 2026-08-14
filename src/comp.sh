#!/usr/bin/env bash
set -euo pipefail

CXX="${CXX:-g++}"
compile_omp=(-fopenmp)
link_omp=(-fopenmp)

# Apple Clang ships without OpenMP enabled, while Homebrew libomp provides the
# headers and runtime. Linux GCC continues to use the ordinary -fopenmp path.
if [[ "$(uname -s)" == "Darwin" ]]; then
    libomp="${LIBOMP_PREFIX:-$(brew --prefix libomp)}"
    compile_omp=(-Xpreprocessor -fopenmp -I"$libomp/include")
    link_omp=(-Xpreprocessor -fopenmp -I"$libomp/include" \
              -L"$libomp/lib" -Wl,-rpath,"$libomp/lib" -lomp)
fi

"$CXX" -std=c++17 -O2 -Wall "${compile_omp[@]}" -c Envelope.cpp -o Envelope.o
"$CXX" -std=c++17 -O2 -Wall "${compile_omp[@]}" FindEnv.cpp Envelope.o \
    "${link_omp[@]}" -o find_envelope
"$CXX" -std=c++17 -O2 -Wall "${compile_omp[@]}" FindPair.cpp Envelope.o \
    "${link_omp[@]}" -o find_pair
