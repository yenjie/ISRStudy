#!/usr/bin/env bash
# Build the KKMCee 5 driver against the KKMCee 5 tree, Pythia 8, HepMC3,
# Photos++ and ROOT.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
KK="${KKMCEE5:-/raid5/data/yjlee/ISR/KKMCee5}"
PHOTOS="${PHOTOS_INSTALL:-/raid5/data/yjlee/ISR/external/photos-install}"
ROOTSYS="${ROOTSYS:-/raid5/root/root-v6.34.04/root}"
export ROOTSYS PATH="$ROOTSYS/bin:$PATH" LD_LIBRARY_PATH="$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"

LIBDIRS="$KK/SRCee/.libs $KK/MCdev/.libs $KK/Foam2/.libs $KK/tauola-photos/.libs $PHOTOS/lib $ROOTSYS/lib"
RPATH=""; for d in $LIBDIRS; do RPATH="$RPATH -Wl,-rpath,$d"; done

g++ -std=c++17 -O2 -g "$HERE/kkmcee5_dump.cxx" -o "$HERE/kkmcee5_dump" \
  -I"$KK/SRCee" -I"$KK/MCdev" -I"$KK/Foam2" -I"$PHOTOS/include" -I/usr/include \
  $(root-config --cflags) \
  -L"$KK/SRCee/.libs" -L"$KK/MCdev/.libs" -L"$KK/Foam2/.libs" -L"$KK/tauola-photos/.libs" -L"$PHOTOS/lib" \
  -lKKee -lKKfm -lFOAM -lMCdev -lTauolaPhotos \
  -lPhotospp -lPhotosppHepMC3 -lPhotosppHEPEVT -lHepMC3 -lpythia8 \
  $(root-config --libs) -lgfortran -ldl $RPATH
echo "[done] built $HERE/kkmcee5_dump"
