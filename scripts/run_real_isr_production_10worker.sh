#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO_ROOT"

EVENTS="${EVENTS:-2000000}"
ECM="${ECM:-91.1876}"
OUTDIR="${OUTDIR:-/data2/yjlee/ISRsample/real_2M_20260511}"
WORKDIR="${WORKDIR:-$OUTDIR/work}"
MAX_WORKERS="${MAX_WORKERS:-10}"
FORCE="${FORCE:-0}"

mkdir -p "$OUTDIR" "$WORKDIR"

targets=(
  pythia-off
  pythia-on
  pythia-vincia-off
  pythia-vincia-on
  herwig-off
  herwig-qedshower
  sherpa-off
  sherpa-on
  sherpa-yfs
)

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

for target in "${targets[@]}"; do
  wait_for_slot
  log="$WORKDIR/driver_${target}.log"
  echo "[launch] $target -> $log"
  (
    EVENTS="$EVENTS" ECM="$ECM" OUTDIR="$OUTDIR" WORKDIR="$WORKDIR" FORCE="$FORCE" \
      TARGETS="$target" "$REPO_ROOT/scripts/run_real_isr_production.sh"
  ) > "$log" 2>&1 &
  pid=$!
  pids+=("$pid")
  labels["$pid"]="$target"
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
  echo "[error] one or more 2M production tasks failed" >&2
  exit "$status"
fi

echo "[done] 10-worker capped real ISR production in $OUTDIR"
