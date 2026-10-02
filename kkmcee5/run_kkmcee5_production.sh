#!/usr/bin/env bash
# Produce a KKMCee 5 sample: N concurrent jobs with distinct seeds, each in its
# own subdirectory sharing the run card, defaults and DIZET tables, then the
# dumps are concatenated and converted with real_isr_ntuple_producer --mode kkmc.
#
# Usage: run_kkmcee5_production.sh <ISR_ON|ISR_OFF> <events per job> <jobs>
set -euo pipefail
HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
STATE="${1:?ISR_ON or ISR_OFF}"; PERJOB="${2:?events per job}"; NJOBS="${3:?jobs}"
BASE="${KKMCEE5_OUT:-/raid5/data/yjlee/ISR/samples/kkmcee5_20261002}"
RUN="$BASE/$STATE"
ROOTSYS="${ROOTSYS:-/raid5/root/root-v6.34.04/root}"
export ROOTSYS PATH="$ROOTSYS/bin:$PATH" LD_LIBRARY_PATH="$ROOTSYS/lib:${LD_LIBRARY_PATH:-}"
[[ -s "$RUN/DIZET-table1" ]] || { echo "no DIZET tables in $RUN; run DZface/TabMainC there first" >&2; exit 1; }
[[ -x "$HERE/kkmcee5_dump" ]] || { echo "build the driver first: $HERE/build_kkmcee5_dump.sh" >&2; exit 1; }

echo "[1/3] $STATE: $NJOBS jobs x $PERJOB events"
pids=()
for j in $(seq 1 "$NJOBS"); do
  d="$RUN/job$j"; mkdir -p "$d"
  ln -sfn "$RUN/KKMCee_defaults" "$d/KKMCee_defaults"
  ln -sfn "$RUN/pro.input" "$d/pro.input"
  for t in "$RUN"/DIZET-table*; do ln -sfn "$t" "$d/$(basename "$t")"; done
  off=0; [[ "$STATE" == ISR_OFF ]] && off=500
  seed=$((20261002 + 1000 * j + off))
  ( cd "$d" && "$HERE/kkmcee5_dump" "$PERJOB" events.dat "$seed" > job.log 2>&1 ) &
  pids+=($!)
done
fail=0; for p in "${pids[@]}"; do wait "$p" || fail=1; done
[[ $fail -eq 0 ]] || { echo "a job failed; see $RUN/job*/job.log" >&2; exit 1; }

echo "[2/3] concatenating"
head -2 "$RUN/job1/events.dat" > "$RUN/kkmcee5_events.dat"
for j in $(seq 1 "$NJOBS"); do grep -v '^#' "$RUN/job$j/events.dat" >> "$RUN/kkmcee5_events.dat"; done
TOTAL=$((PERJOB * NJOBS))

echo "[3/3] converting to the ISR-study ntuple"
isr=1; [[ "$STATE" == ISR_OFF ]] && isr=0
"$HERE/../real_isr_ntuple_producer" --mode kkmc --kkmcInput "$RUN/kkmcee5_events.dat" \
  --nEvents "$TOTAL" --sqrtS 91.1876 --isrOn "$isr" \
  --generatorName "KKMCee50002_${STATE}" --generatorId 6 \
  --output "$BASE/mc_KKMCee50002_${STATE}.root"
gzip -f "$RUN/kkmcee5_events.dat"
rm -f "$RUN"/job*/events.dat
echo "[done] $BASE/mc_KKMCee50002_${STATE}.root"
