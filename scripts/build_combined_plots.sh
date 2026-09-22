#!/usr/bin/env bash
# Concatenate the ISR study figures into one browsable PDF.
#
# The page order is fixed here rather than by a glob so that the combined file
# is reproducible and so that a superseded figure cannot creep back in.
set -euo pipefail

RES="${1:-/data2/yjlee/ISRsample/kkmc_1M_20260921/results}"
OUT="$RES/isr_all_plots.pdf"

PAGES=(
  "$RES/isr_model_cisr_thrust.pdf"              # thrust correction, main result
  "$RES/isr_model_photon_energy.pdf"            # radiated photon energy
  "$RES/isr_model_kkmc_ifi.pdf"                 # initial-final interference
  "$RES/eec_isr_correction_doublelog_cropped.pdf"  # charged EEC, double-log
)

for p in "${PAGES[@]}"; do
  [[ -f "$p" ]] || { echo "missing: $p" >&2; exit 1; }
done

pdfunite "${PAGES[@]}" "$OUT"
echo "[done] $OUT ($(pdfinfo "$OUT" | awk '/^Pages/{print $2}') pages)"
