#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO_ROOT"

ROOT_CONFIG="${ROOT_CONFIG:-/raid5/root/root-v6.34.04/root/bin/root-config}"

g++ -std=c++17 -O2 src/real_isr_ntuple_producer.cc \
  -I/usr/include -I"$("$ROOT_CONFIG" --incdir)" \
  -L/usr/lib -L/usr/lib/x86_64-linux-gnu \
  -Wl,-rpath,/usr/lib -Wl,-rpath,/usr/lib/x86_64-linux-gnu \
  -lpythia8 -lHepMC3 $("$ROOT_CONFIG" --libs) -ldl \
  -o real_isr_ntuple_producer

echo "[done] built $REPO_ROOT/real_isr_ntuple_producer"

g++ -std=c++17 -O2 src/make_endpoint_diagnostics.cc \
  -I"$("$ROOT_CONFIG" --incdir)" \
  $("$ROOT_CONFIG" --libs) \
  -o make_endpoint_diagnostics

echo "[done] built $REPO_ROOT/make_endpoint_diagnostics"
