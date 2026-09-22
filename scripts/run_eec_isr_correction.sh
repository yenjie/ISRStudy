#!/usr/bin/env bash
# Rebuild the charged-EEC ISR correction tables and the double-log figure.
#
#   1. plot_eec_isr_correction.C   event loop -> per-bin and per-region CSV
#   2. plot_eec_isr_doublelog.C    CSV -> figure
#
# The event loop is O(N_charged^2) per event over ~30M events, which is several
# hours in one process, so the six sample pairs are run concurrently and their
# CSV parts concatenated.  Step 2 is cheap: restyling the figure only needs it.
set -euo pipefail

ROOTSYS="${ROOTSYS:-/raid5/root/root-v6.34.04/root}"
export ROOTSYS
export PATH="$ROOTSYS/bin:$PATH"
export LD_LIBRARY_PATH="$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"

OUT="${1:-/data2/yjlee/ISRsample/kkmc_1M_20260921/results}"
NSAMPLES="${NSAMPLES:-6}"
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
mkdir -p "$OUT"

echo "[1/4] event loop, $NSAMPLES concurrent jobs -> $OUT"
pids=()
for i in $(seq 0 $((NSAMPLES - 1))); do
  root -l -b -q "$HERE/macros/plot_eec_isr_correction.C(\"$OUT\",-1,$i)" \
      > "$OUT/eec_job_$i.log" 2>&1 &
  pids+=($!)
done
fail=0
for p in "${pids[@]}"; do wait "$p" || fail=1; done
[[ $fail -eq 0 ]] || { echo "at least one job failed; see $OUT/eec_job_*.log" >&2; exit 1; }

echo "[2/4] concatenating parts"
for base in eec_isr_correction eec_isr_correction_regions; do
  head -1 "$OUT/${base}_s0.csv" > "$OUT/$base.csv"
  for i in $(seq 0 $((NSAMPLES - 1))); do
    tail -n +2 "$OUT/${base}_s$i.csv" >> "$OUT/$base.csv"
  done
  rm -f "$OUT/${base}"_s*.csv
done

echo "[3/4] figure -> $OUT"
root -l -b -q "$HERE/macros/plot_eec_isr_doublelog.C(\"$OUT\")"

# ROOT writes a portrait PDF page whatever the canvas aspect, so the vector
# version needs cropping before it goes into a deck or the combined PDF.
pdfcrop --margins 2 "$OUT/eec_isr_correction_doublelog.pdf" \
        "$OUT/eec_isr_correction_doublelog_cropped.pdf" > /dev/null

# Verification: the closed-form closure and the charged-energy radiation
# measurement.  A linear pass, so all six run concurrently in a couple of
# minutes even at full statistics.
echo "[4/4] verification -> $OUT"
pids=()
for i in $(seq 0 $((NSAMPLES - 1))); do
  root -l -b -q "$HERE/macros/verify_eec_isr_correction.C(\"$OUT\",-1,$i)" \
      > "$OUT/verify_job_$i.log" 2>&1 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p" || { echo "verification job failed" >&2; exit 1; }; done
head -1 "$OUT/eec_verification_s0.csv" > "$OUT/eec_verification.csv"
for i in $(seq 0 $((NSAMPLES - 1))); do
  tail -n +2 "$OUT/eec_verification_s$i.csv" >> "$OUT/eec_verification.csv"
done
rm -f "$OUT"/eec_verification_s*.csv
root -l -b -q "$HERE/macros/plot_eec_verification.C(\"$OUT\")"
pdfcrop --margins 2 "$OUT/eec_isr_verification.pdf" \
        "$OUT/eec_isr_verification_cropped.pdf" > /dev/null

echo "[done] $OUT"
