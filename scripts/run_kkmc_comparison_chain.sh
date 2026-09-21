#!/usr/bin/env bash
#
# Waits for the KKMC production to finish, then runs the downstream chain:
# endpoint diagnostics, the cross-generator comparison macro, and the copy of
# the figures into the Overleaf tree.
#
set -uo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
KKMC_OUT="${KKMC_OUT:-/data2/yjlee/ISRsample/kkmc_1M_20260921}"
OVERLEAF="${OVERLEAF:-/raid5/data/yjlee/ISR/overleaf}"
ROOTSYS_DIR="${ROOTSYS_DIR:-/raid5/root/root-v6.34.04/root}"

export ROOTSYS="$ROOTSYS_DIR"
export PATH="$ROOTSYS/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"

echo "[wait] waiting for the three KKMC drivers to finish"
while [ "$(grep -l 'KKMC production in' "$KKMC_OUT"/drv_*.log 2>/dev/null | wc -l)" -lt 3 ]; do
  sleep 60
done
echo "[wait] all three KKMC drivers reported done"

for f in mc_KKMC424_ISR_OFF.root mc_KKMC424_ISR_ON.root mc_KKMC424_ISR_ON_IFI.root; do
  if [[ ! -s "$KKMC_OUT/$f" ]]; then
    echo "[error] missing $KKMC_OUT/$f" >&2
    exit 1
  fi
done

echo "[step] endpoint diagnostics"
REAL_DIR="$KKMC_OUT" OUTDIR="$KKMC_OUT/endpoint_diagnostics" \
  SAMPLE_SET=KKMC MAX_WORKERS=3 \
  "$REPO_ROOT/scripts/run_endpoint_diagnostics.sh" || exit 1

echo "[step] comparison macro"
root -l -b -q "$REPO_ROOT/macros/plot_isr_model_comparison.C(\"$KKMC_OUT/results\")" || exit 1

echo "[step] copy figures into Overleaf"
mkdir -p "$OVERLEAF/results/kkmc"
cp -f "$KKMC_OUT"/results/isr_model_*.png "$OVERLEAF/results/kkmc/"
cp -f "$KKMC_OUT"/results/isr_model_*.pdf "$OVERLEAF/results/kkmc/" 2>/dev/null
cp -f "$KKMC_OUT"/results/isr_model_*.csv "$OVERLEAF/results/kkmc/"

echo "[done] KKMC comparison chain complete; outputs in $KKMC_OUT/results"
