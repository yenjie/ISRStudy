#!/usr/bin/env bash
set -euo pipefail

PKG_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
REPO_ROOT="$(cd "$PKG_ROOT/.." && pwd -P)"

EVENTS="${EVENTS:-1000}"
ECM="${ECM:-91.1876}"
OUTDIR="${OUTDIR:-}"
FORCE="${FORCE:-0}"
TARGETS="${TARGETS:-OFF,PDFESherpa,YFS}"
WORKDIR="${WORKDIR:-$OUTDIR/work}"
HEPMC_FORMAT="${HEPMC_FORMAT:-hepmc2}"
MAX_WORKERS="${MAX_WORKERS:-1}"
FADGEN_OUTPUT_KIND="${FADGEN_OUTPUT_KIND:-delsim-lujets}"

SHERPA_BIN="${SHERPA_BIN:-/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3-build/outputs/bin/Sherpa}"
SHERPA_BUILD_ROOT="${SHERPA_BUILD_ROOT:-/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3-build}"
SHERPA_SOURCE_ROOT="${SHERPA_SOURCE_ROOT:-/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3}"

if [[ -z "$OUTDIR" ]]; then
  echo "Set OUTDIR, e.g. OUTDIR=/data2/yjlee/ISRsample/sherpa_fadgen_smoke" >&2
  exit 2
fi

mkdir -p "$OUTDIR"/{cards,hepmc,fad,root_validation,logs,metadata} "$WORKDIR"

if [[ ! -x "$PKG_ROOT/bin/hepmc3_to_fadgen" ]]; then
  "$PKG_ROOT/scripts/build.sh"
fi
if [[ ! -x "$REPO_ROOT/real_isr_ntuple_producer" ]]; then
  "$REPO_ROOT/scripts/build_real_isr_tools.sh"
fi
if [[ ! -x "$SHERPA_BIN" ]]; then
  echo "Missing Sherpa binary: $SHERPA_BIN" >&2
  exit 1
fi

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

contains_target() {
  local needle="$1"
  [[ ",$TARGETS," == *",$needle,"* || "$TARGETS" == "all" ]]
}

mode_config() {
  local mode="$1"
  case "$mode" in
    OFF)
      echo "sherpa_zpole_fad_ISR_OFF.yaml.template sherpa_OFF 0 1400510"
      ;;
    PDFESherpa)
      echo "sherpa_zpole_fad_ISR_ON.yaml.template sherpa_PDFESherpa 1 1400511"
      ;;
    YFS)
      echo "sherpa_zpole_fad_ISR_YFS.yaml.template sherpa_YFS 1 1400512"
      ;;
    *)
      echo "Unknown mode: $mode" >&2
      exit 2
      ;;
  esac
}

run_mode() {
  local mode="$1"
  read -r template_name sample isr seed <<< "$(mode_config "$mode")"
  local fad="$OUTDIR/fad/${sample}.fadgen"
  local hepmc="$OUTDIR/hepmc/${sample}.hepmc"
  local root="$OUTDIR/root_validation/${sample}.root"
  local metadata="$OUTDIR/metadata/${sample}.json"
  local card_out="$OUTDIR/cards/${sample}.yaml"
  local log_prefix="$OUTDIR/logs/${sample}"
  local work="$WORKDIR/${sample}"
  mkdir -p "$work"

  if [[ -f "$fad" && -f "$root" && "$FORCE" != "1" ]]; then
    echo "[reuse] $sample"
    return
  fi

  local output_basename="${sample}.hepmc"
  sed -e "s|@OUTPUT_BASENAME@|$output_basename|g" \
    "$PKG_ROOT/cards/$template_name" > "$card_out"
  cp "$card_out" "$work/Sherpa.yaml"

  echo "[sherpa] $mode -> $hepmc"
  rm -f "$hepmc" "$fad" "$root" "$metadata"
  (
    cd "$work"
    "$SHERPA_BIN" -f "$work/Sherpa.yaml" -e "$EVENTS" -R "$seed" -g
  ) > "${log_prefix}.sherpa.log" 2>&1
  if [[ -f "$work/$output_basename" ]]; then
    mv "$work/$output_basename" "$hepmc"
  fi
  if [[ ! -s "$hepmc" ]]; then
    echo "Sherpa did not produce expected HepMC output: $hepmc" >&2
    echo "See ${log_prefix}.sherpa.log" >&2
    exit 1
  fi

  echo "[convert] $hepmc -> $fad"
  "$PKG_ROOT/bin/hepmc3_to_fadgen" \
    --input "$hepmc" \
    --output "$fad" \
    --metadata "$metadata" \
    --format "$HEPMC_FORMAT" \
    --outputKind "$FADGEN_OUTPUT_KIND" \
    --generatorName "Sherpa 3.0.3" \
    --mode "$mode" \
    --sqrtS "$ECM" \
    --maxEvents "$EVENTS" > "${log_prefix}.fadgen.log" 2>&1

  echo "[root-validation] $hepmc -> $root"
  "$REPO_ROOT/real_isr_ntuple_producer" \
    --mode hepmc \
    --hepmcInput "$hepmc" \
    --hepmcFormat "$HEPMC_FORMAT" \
    --nEvents "$EVENTS" \
    --sqrtS "$ECM" \
    --isrOn "$isr" \
    --generatorName "Sherpa303_${mode}" \
    --generatorId 30 \
    --output "$root" > "${log_prefix}.root_validation.log" 2>&1

  sha256sum "$card_out" "$hepmc" "$fad" "$root" > "$OUTDIR/metadata/${sample}.sha256"
  echo "[done] $sample"
}

setup_sherpa_runtime

pids=()
for mode in OFF PDFESherpa YFS; do
  if contains_target "$mode"; then
    run_mode "$mode" &
    pids+=("$!")
    while [[ "${#pids[@]}" -ge "$MAX_WORKERS" ]]; do
      wait "${pids[0]}"
      pids=("${pids[@]:1}")
    done
  fi
done
for pid in "${pids[@]}"; do
  wait "$pid"
done

cat > "$OUTDIR/metadata/production_summary.md" <<EOF_SUMMARY
# Sherpa to FADGEN Production

- date: $(date -Iseconds)
- events per mode: $EVENTS
- sqrtS: $ECM
- targets: $TARGETS
- HepMC format: $HEPMC_FORMAT
- FADGEN output kind: $FADGEN_OUTPUT_KIND
- Sherpa binary: $SHERPA_BIN
- output: $OUTDIR

## Validation next step

Run:

\`\`\`bash
SKELEANA_CMD=/path/to/skeleana \\
FADGEN_INPUT=$OUTDIR/fad/sherpa_PDFESherpa.fadgen \\
OUTDIR=$OUTDIR/skeleana \\
$PKG_ROOT/scripts/run_skeleana_validation.sh
\`\`\`
EOF_SUMMARY

echo "[done] production summary: $OUTDIR/metadata/production_summary.md"
