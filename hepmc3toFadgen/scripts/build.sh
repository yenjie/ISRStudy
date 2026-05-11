#!/usr/bin/env bash
set -euo pipefail

PKG_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$PKG_ROOT"

CXX="${CXX:-g++}"
CXXFLAGS="${CXXFLAGS:--std=c++17 -O2 -Wall -Wextra}"
mkdir -p bin

"$CXX" $CXXFLAGS src/hepmc3_to_fadgen.cc \
  -I/usr/include \
  -L/usr/lib -L/usr/lib/x86_64-linux-gnu \
  -Wl,-rpath,/usr/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu \
  -lHepMC3 \
  -o bin/hepmc3_to_fadgen

echo "[done] built $PKG_ROOT/bin/hepmc3_to_fadgen"
