#!/usr/bin/env bash
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
cd "$REPO_ROOT"

EVENTS="${EVENTS:-100000}"
ECM="${ECM:-91.2}"
OUTDIR="${OUTDIR:-/data2/yjlee/ISRsample/real_generators_20260510}"
WORKDIR="${WORKDIR:-$OUTDIR/work}"
FORCE="${FORCE:-0}"
TARGETS="${TARGETS:-pythia,pythia-vincia,herwig,sherpa}"

HERWIG_DIR="${HERWIG_DIR:-/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/herwig-7.3.0}"
HERWIG_BIN="${HERWIG_BIN:-$HERWIG_DIR/bin/Herwig}"
HERWIG_EECOLLIDER="${HERWIG_EECOLLIDER:-$HERWIG_DIR/share/Herwig/snippets/EECollider.in}"

SHERPA_BIN="${SHERPA_BIN:-/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3-build/outputs/bin/Sherpa}"
SHERPA_BUILD_ROOT="${SHERPA_BUILD_ROOT:-/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3-build}"
SHERPA_SOURCE_ROOT="${SHERPA_SOURCE_ROOT:-/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3}"

mkdir -p "$OUTDIR" "$WORKDIR"

if [[ ! -x "$REPO_ROOT/real_isr_ntuple_producer" ]]; then
  "$REPO_ROOT/scripts/build_real_isr_tools.sh"
fi

contains_target() {
  local needle="$1"
  [[ ",$TARGETS," == *",$needle,"* || "$TARGETS" == "all" ]]
}

run_pythia() {
  local tag="$1"
  local name="$2"
  local gid="$3"
  local isr="$4"
  local vincia="$5"
  local seed="$6"
  local mode="OFF"
  [[ "$isr" == "1" ]] && mode="ON"
  local out="$OUTDIR/mc_${tag}_ISR_${mode}.root"
  if [[ -f "$out" && "$FORCE" != "1" ]]; then
    echo "[reuse] $out"
    return
  fi
  "$REPO_ROOT/real_isr_ntuple_producer" \
    --mode pythia \
    --nEvents "$EVENTS" \
    --sqrtS "$ECM" \
    --isrOn "$isr" \
    --vincia "$vincia" \
    --seed "$seed" \
    --generatorName "$name" \
    --generatorId "$gid" \
    --output "$out"
}

run_herwig() {
  local isr="$1"
  local seed="$2"
  local mode="OFF"
  [[ "$isr" == "1" ]] && mode="ON"
  local out="$OUTDIR/mc_Herwig730_ISR_${mode}.root"
  if [[ -f "$out" && "$FORCE" != "1" ]]; then
    echo "[reuse] $out"
    return
  fi
  if [[ ! -x "$HERWIG_BIN" ]]; then
    echo "Missing Herwig binary: $HERWIG_BIN" >&2
    exit 1
  fi
  local work="$WORKDIR/herwig_ISR_${mode}"
  mkdir -p "$work"
  local run_name="LEP_Herwig730_ISR_${mode}"
  local card="$work/${run_name}.in"
  local setup="$work/${run_name}.setup"
  local fifo="$work/${run_name}.hepmc2"
  local template="$REPO_ROOT/cards/herwig_zpole_hepmc_ISR_${mode}.in.template"
  sed \
    -e "s|@EECOLLIDER@|$HERWIG_EECOLLIDER|g" \
    -e "s|@ENERGY@|$ECM|g" \
    -e "s|@RUN_NAME@|$run_name|g" \
    "$template" > "$card"
  cat > "$setup" <<EOF_SETUP
set /Herwig/Analysis/HepMC:Filename $fifo
set /Herwig/Generators/EventGenerator:DebugLevel 0
set /Herwig/Generators/EventGenerator:PrintEvent 0
EOF_SETUP
  rm -f "$fifo"
  mkfifo "$fifo"
  "$REPO_ROOT/real_isr_ntuple_producer" \
    --mode hepmc \
    --hepmcInput "$fifo" \
    --hepmcFormat hepmc2 \
    --nEvents "$EVENTS" \
    --sqrtS "$ECM" \
    --isrOn "$isr" \
    --generatorName "Herwig730" \
    --generatorId 1 \
    --output "$out" > "$work/converter.log" 2>&1 &
  local converter_pid=$!
  (
    cd "$work"
    LD_LIBRARY_PATH="$HERWIG_DIR/lib:${LD_LIBRARY_PATH:-}" "$HERWIG_BIN" read "$card"
  ) > "$work/herwig_read.log" 2>&1
  (
    cd "$work"
    LD_LIBRARY_PATH="$HERWIG_DIR/lib:${LD_LIBRARY_PATH:-}" "$HERWIG_BIN" run "${run_name}.run" -q -x "$setup" -N "$EVENTS" -s "$seed"
  ) > "$work/herwig_run.log" 2>&1
  wait "$converter_pid"
  rm -f "$fifo"
  echo "[done] $out"
}

setup_sherpa_runtime() {
  local runtime_share="$SHERPA_BUILD_ROOT/outputs/share/SHERPA-MC-runtime"
  local runtime_lib="$SHERPA_BUILD_ROOT/outputs/lib/SHERPA-MC"
  local runtime_include="$SHERPA_BUILD_ROOT/outputs/include/SHERPA-MC"
  local makelibs_source="$SHERPA_BUILD_ROOT/AMEGIC++/Main/makelibs"

  mkdir -p "$runtime_share"
  ln -sfn "$SHERPA_SOURCE_ROOT/HADRONS++/Decaydata.yaml" "$runtime_share/Decaydata.yaml"
  ln -sfn "$SHERPA_SOURCE_ROOT/HADRONS++/Current_Library/PhaseSpaceFunctions" "$runtime_share/PhaseSpaceFunctions"
  ln -sfn "$SHERPA_SOURCE_ROOT/ATOOLS/Math/Sobol" "$runtime_share/Sobol"
  if [[ -d "$SHERPA_SOURCE_ROOT/PDF/CT14/CT14Grid" ]]; then
    ln -sfn "$SHERPA_SOURCE_ROOT/PDF/CT14/CT14Grid" "$runtime_share/CT14Grid"
  fi
  if [[ -f "$makelibs_source" ]]; then
    ln -sfn "$makelibs_source" "$runtime_share/makelibs"
  fi

  export SHERPA_SHARE_PATH="$runtime_share"
  export SHERPA_LIBRARY_PATH="$runtime_lib"
  export SHERPA_INCLUDE_PATH="$runtime_include"
  export LD_LIBRARY_PATH="$runtime_lib:${LD_LIBRARY_PATH:-}"
}

run_sherpa() {
  local isr="$1"
  local seed="$2"
  local mode="OFF"
  [[ "$isr" == "1" ]] && mode="ON"
  local out="$OUTDIR/mc_Sherpa303_ISR_${mode}.root"
  if [[ -f "$out" && "$FORCE" != "1" ]]; then
    echo "[reuse] $out"
    return
  fi
  if [[ ! -x "$SHERPA_BIN" ]]; then
    echo "Missing Sherpa binary: $SHERPA_BIN" >&2
    exit 1
  fi
  setup_sherpa_runtime
  local work="$WORKDIR/sherpa_ISR_${mode}"
  mkdir -p "$work"
  local output_basename="events"
  local card="$work/Sherpa.yaml"
  local fifo="$work/$output_basename"
  sed -e "s|@OUTPUT_BASENAME@|$output_basename|g" \
    "$REPO_ROOT/cards/sherpa_zpole_hepmc_ISR_${mode}.yaml.template" > "$card"
  rm -f "$fifo"
  mkfifo "$fifo"
  "$REPO_ROOT/real_isr_ntuple_producer" \
    --mode hepmc \
    --hepmcInput "$fifo" \
    --hepmcFormat hepmc2 \
    --nEvents "$EVENTS" \
    --sqrtS "$ECM" \
    --isrOn "$isr" \
    --generatorName "Sherpa303" \
    --generatorId 3 \
    --output "$out" > "$work/converter.log" 2>&1 &
  local converter_pid=$!
  (
    cd "$work"
    "$SHERPA_BIN" -f "$card" -e "$EVENTS" -R "$seed" -g
  ) > "$work/sherpa.log" 2>&1
  wait "$converter_pid"
  rm -f "$fifo"
  echo "[done] $out"
}

if contains_target pythia; then
  run_pythia Pythia8315 PYTHIA8315 2 0 0 1200510
  run_pythia Pythia8315 PYTHIA8315 2 1 0 1200511
fi

if contains_target pythia-vincia; then
  run_pythia Pythia8315_Vincia Pythia8315_Vincia 4 0 1 1200520
  run_pythia Pythia8315_Vincia Pythia8315_Vincia 4 1 1 1200521
fi

if contains_target herwig; then
  run_herwig 0 1300510
  run_herwig 1 1300511
fi

if contains_target sherpa; then
  run_sherpa 0 1400510
  run_sherpa 1 1400511
fi

echo "[done] real ISR production in $OUTDIR"
