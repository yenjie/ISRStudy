#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO_ROOT"

REAL_DIR="${REAL_DIR:-/data2/yjlee/ISRsample/real_3M_20260511}"
OUTDIR="${OUTDIR:-$REAL_DIR/endpoint_diagnostics}"
MAX_WORKERS="${MAX_WORKERS:-5}"
MAX_EVENTS="${MAX_EVENTS:-}"
MAX_EXACT_PARTICLES="${MAX_EXACT_PARTICLES:-0}"
FORCE="${FORCE:-0}"
SAMPLE_SET="${SAMPLE_SET:-ON}"

mkdir -p "$OUTDIR"

if [[ ! -x "$REPO_ROOT/make_endpoint_diagnostics" ]]; then
  "$REPO_ROOT/scripts/build_real_isr_tools.sh"
fi

on_samples=(
  "Herwig730_QEDshower:mc_Herwig730_QEDshower.root"
  "Pythia8315:mc_Pythia8315_ISR_ON.root"
  "Sherpa303:mc_Sherpa303_ISR_ON.root"
  "Sherpa303_YFS:mc_Sherpa303_ISR_YFS.root"
  "Pythia8315_Vincia:mc_Pythia8315_Vincia_ISR_ON.root"
)

off_samples=(
  "Herwig730_OFF:mc_Herwig730_ISR_OFF.root"
  "Pythia8315_OFF:mc_Pythia8315_ISR_OFF.root"
  "Sherpa303_OFF:mc_Sherpa303_ISR_OFF.root"
  "Pythia8315_Vincia_OFF:mc_Pythia8315_Vincia_ISR_OFF.root"
)

samples=()
case "$SAMPLE_SET" in
  ON)
    samples=("${on_samples[@]}")
    ;;
  OFF)
    samples=("${off_samples[@]}")
    ;;
  ALL)
    samples=("${on_samples[@]}" "${off_samples[@]}")
    ;;
  *)
    echo "Unknown SAMPLE_SET=$SAMPLE_SET; use ON, OFF, or ALL" >&2
    exit 2
    ;;
esac

declare -a pids=()
declare -A labels=()

wait_for_slot() {
  while (( ${#pids[@]} >= MAX_WORKERS )); do
    local next=("${pids[@]}")
    pids=()
    for pid in "${next[@]}"; do
      if kill -0 "$pid" 2>/dev/null; then
        pids+=("$pid")
      else
        wait "$pid"
      fi
    done
    sleep 2
  done
}

for sample in "${samples[@]}"; do
  tag="${sample%%:*}"
  file="${sample#*:}"
  input="$REAL_DIR/$file"
  output="$OUTDIR/endpoint_diagnostics_${tag}.root"
  log="$OUTDIR/endpoint_diagnostics_${tag}.log"
  if [[ -f "$output" && "$FORCE" != "1" ]]; then
    echo "[reuse] $output"
    continue
  fi
  wait_for_slot
  echo "[launch] $tag -> $output"
  args=(--input "$input" --output "$output" --maxExactParticles "$MAX_EXACT_PARTICLES")
  if [[ -n "$MAX_EVENTS" ]]; then
    args+=(--maxEvents "$MAX_EVENTS")
  fi
  "$REPO_ROOT/make_endpoint_diagnostics" "${args[@]}" > "$log" 2>&1 &
  pid=$!
  pids+=("$pid")
  labels["$pid"]="$tag"
done

status=0
for pid in "${pids[@]}"; do
  if wait "$pid"; then
    echo "[done] ${labels[$pid]}"
  else
    echo "[fail] ${labels[$pid]}" >&2
    status=1
  fi
done

if (( status != 0 )); then
  echo "[error] one or more endpoint diagnostic tasks failed" >&2
  exit "$status"
fi

echo "[done] endpoint diagnostics in $OUTDIR"
