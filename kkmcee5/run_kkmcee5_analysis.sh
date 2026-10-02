#!/usr/bin/env bash
# Analyse the KKMCee 5.00.02 pair once the ntuples exist: endpoint diagnostics,
# thrust C_ISR, charged EEC, photon summary, and comparison figures.
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
S="${KKMCEE5_OUT:-/raid5/data/yjlee/ISR/samples/kkmcee5_20261002}"
export KKMCEE5_OUT="$S"
ROOTSYS="${ROOTSYS:-/raid5/root/root-v6.34.04/root}"
export ROOTSYS PATH="$ROOTSYS/bin:$PATH" LD_LIBRARY_PATH="$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"
mkdir -p "$S/endpoint_diagnostics" "$S/results"

echo "[1/4] endpoint diagnostics"
pids=()
for st in ON OFF; do
  "$HERE/make_endpoint_diagnostics" --input "$S/mc_KKMCee50002_ISR_$st.root" \
      --output "$S/endpoint_diagnostics/endpoint_diagnostics_KKMCee50002_ISR_$st.root" \
      > "$S/endpoint_diagnostics/endpoint_diagnostics_KKMCee50002_ISR_$st.log" 2>&1 &
  pids+=($!)
done
for p in "${pids[@]}"; do wait "$p"; done

echo "[2/4] thrust C_ISR per bin and definition"
root -l -b -q "$HERE/macros/kkmcee5_thrust_cisr.C(\"$S/endpoint_diagnostics\",\"$S/results\")"

echo "[3/4] charged EEC correction, verification and photon summary (sample index 6)"
root -l -b -q "$HERE/macros/plot_eec_isr_correction.C(\"$S/results\",-1,6)" > "$S/results/eec_job_6.log" 2>&1 &
eec_pid=$!
root -l -b -q "$HERE/macros/verify_eec_isr_correction.C(\"$S/results\",-1,6)" > "$S/results/verify_job_6.log" 2>&1 &
verify_pid=$!
root -l -b -q "$HERE/macros/kkmcee5_photon_summary.C(\"$S\",\"$S/results/isr_model_summary_kkmcee5.csv\")" > "$S/results/photon_summary.log" 2>&1 &
photon_pid=$!
for p in "$eec_pid" "$verify_pid" "$photon_pid"; do wait "$p"; done

echo "[4/4] figures and slide tables"
bash "$HERE/kkmcee5/finalize_kkmcee5_results.sh"
echo "[done] outputs in $S/results"
ls -la "$S/results"
