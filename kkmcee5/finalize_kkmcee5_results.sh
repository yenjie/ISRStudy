#!/usr/bin/env bash
# Rebuild the comparison figures and slide tables from completed KKMCee 5 CSVs.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
S="${KKMCEE5_OUT:-/raid5/data/yjlee/ISR/samples/kkmcee5_20261002}"
export KKMCEE5_OUT="$S"
DEST="${OVERLEAF_KKMC_RESULTS:-/raid5/data/yjlee/ISR/overleaf/results/kkmc}"
ALEPH_CORRECTION_ROOT="${ALEPH_CORRECTION_ROOT:-/raid5/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/Isr/isr_corr_precomputed.root}"
ROOTSYS="${ROOTSYS:-/raid5/root/root-v6.34.04/root}"
export ROOTSYS PATH="$ROOTSYS/bin:$PATH" LD_LIBRARY_PATH="$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"

python3 "$HERE/scripts/merge_kkmcee5_eec_results.py" \
  --baseline "$DEST" --new "$S/results" --output "$S/results"

root -l -b -q "$HERE/macros/plot_cisr_tau_kkmcee5.C(\"$DEST/isr_model_cisr_thrust_bins.csv\",\"$S/results/isr_model_cisr_thrust_bins_kkmcee5.csv\",\"$DEST/anthony_cisr_tau_digitized.csv\",\"$ALEPH_CORRECTION_ROOT\",\"$S/results\")"
root -l -b -q "$HERE/macros/plot_eec_isr_doublelog.C(\"$S/results\")"
root -l -b -q "$HERE/macros/plot_eec_verification.C(\"$S/results\")"
pdfcrop --margins 2 "$S/results/isr_model_cisr_tau_with_kkmcee5.pdf" \
  "$S/results/isr_model_cisr_tau_with_kkmcee5_cropped.pdf" > /dev/null
pdfcrop --margins 2 "$S/results/eec_isr_correction_doublelog.pdf" \
  "$S/results/eec_isr_correction_doublelog_cropped.pdf" > /dev/null
pdfcrop --margins 2 "$S/results/eec_isr_verification.pdf" \
  "$S/results/eec_isr_verification_cropped.pdf" > /dev/null

python3 "$HERE/scripts/make_eec_region_table.py" \
  --regions "$S/results/eec_isr_correction_regions.csv" \
  --out "$S/results/eec_isr_region_table.tex"
python3 "$HERE/scripts/make_eec_verification_tables.py" \
  --results "$S/results" --out "$S/results"
python3 "$HERE/scripts/make_kkmcee5_table.py" \
  --main "$DEST/isr_model_cisr_thrust_bins.csv" \
  --kkmcee5 "$S/results/isr_model_cisr_thrust_bins_kkmcee5.csv" \
  --anthony "$DEST/anthony_cisr_tau_digitized.csv" \
  --out "$S/results/kkmcee5_table.tex"

cp "$S/results"/eec_isr_correction.csv "$S/results"/eec_isr_correction_regions.csv \
   "$S/results"/eec_verification.csv "$S/results"/isr_model_summary.csv \
   "$S/results"/eec_isr_correction_doublelog.* \
   "$S/results"/eec_isr_correction_doublelog_cropped.pdf \
   "$S/results"/eec_isr_verification.* \
   "$S/results"/eec_isr_verification_cropped.pdf \
   "$S/results"/eec_isr_region_table.tex "$S/results"/eec_closure_table.tex \
   "$S/results"/eec_energy_table.tex "$S/results"/kkmcee5_table.tex \
   "$S/results"/isr_model_cisr_tau_with_kkmcee5.* \
   "$S/results"/isr_model_cisr_tau_with_kkmcee5_cropped.pdf \
   "$S/results"/isr_model_cisr_thrust_with_kkmcee5.* \
   "$S/results"/isr_model_cisr_thrust_bins_kkmcee5.csv \
   "$S/results"/kkmcee5_xsec.csv "$DEST/"

"$HERE/scripts/build_combined_plots.sh" "$DEST"
echo "[done] KKMCee 5 comparison products in $DEST"
