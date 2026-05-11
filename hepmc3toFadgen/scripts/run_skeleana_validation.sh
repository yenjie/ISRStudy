#!/usr/bin/env bash
set -euo pipefail

PKG_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
OUTDIR="${OUTDIR:-$PKG_ROOT/skeleana_validation}"
SKELEANA_CMD="${SKELEANA_CMD:-}"
DELPHI_REFERENCE_INPUT="${DELPHI_REFERENCE_INPUT:-}"
FADGEN_INPUT="${FADGEN_INPUT:-}"
SKELEANA_REFERENCE_ARGS="${SKELEANA_REFERENCE_ARGS:-}"
SKELEANA_FADGEN_ARGS="${SKELEANA_FADGEN_ARGS:-}"

mkdir -p "$OUTDIR"
manifest="$OUTDIR/validation_manifest.md"

{
  echo "# DELPHI skelana Validation Manifest"
  echo
  echo "- package: $PKG_ROOT"
  echo "- date: $(date -Iseconds)"
  echo "- SKELEANA_CMD: ${SKELEANA_CMD:-unset}"
  echo "- DELPHI_REFERENCE_INPUT: ${DELPHI_REFERENCE_INPUT:-unset}"
  echo "- FADGEN_INPUT: ${FADGEN_INPUT:-unset}"
  echo
} > "$manifest"

if [[ -z "$SKELEANA_CMD" ]]; then
  {
    echo "## Status"
    echo
    echo "Skipped: set SKELEANA_CMD to the local DELPHI skelana executable or wrapper."
  } >> "$manifest"
  echo "[skip] SKELEANA_CMD is not set; wrote $manifest" >&2
  exit 2
fi

if [[ ! -x "$SKELEANA_CMD" ]]; then
  {
    echo "## Status"
    echo
    echo "Failed: SKELEANA_CMD is not executable."
  } >> "$manifest"
  echo "[error] SKELEANA_CMD is not executable: $SKELEANA_CMD" >&2
  exit 2
fi

run_skeleana() {
  local label="$1"
  local input="$2"
  local arg_string="$3"
  local log="$OUTDIR/${label}.log"
  echo "[run] $label -> $log"
  if [[ -n "$arg_string" ]]; then
    read -r -a args <<< "$arg_string"
    "$SKELEANA_CMD" "${args[@]}" > "$log" 2>&1
  else
    "$SKELEANA_CMD" "$input" > "$log" 2>&1
  fi
  {
    echo "## $label"
    echo
    echo "- input: $input"
    echo "- log: $log"
    echo "- status: success"
    echo
  } >> "$manifest"
}

if [[ -n "$DELPHI_REFERENCE_INPUT" || -n "$SKELEANA_REFERENCE_ARGS" ]]; then
  run_skeleana "official_delphi_reference" "$DELPHI_REFERENCE_INPUT" "$SKELEANA_REFERENCE_ARGS"
else
  {
    echo "## official_delphi_reference"
    echo
    echo "Skipped: set DELPHI_REFERENCE_INPUT or SKELEANA_REFERENCE_ARGS."
    echo
  } >> "$manifest"
fi

if [[ -n "$FADGEN_INPUT" || -n "$SKELEANA_FADGEN_ARGS" ]]; then
  run_skeleana "sherpa_fadgen" "$FADGEN_INPUT" "$SKELEANA_FADGEN_ARGS"
else
  {
    echo "## sherpa_fadgen"
    echo
    echo "Skipped: set FADGEN_INPUT or SKELEANA_FADGEN_ARGS."
    echo
  } >> "$manifest"
fi

echo "[done] validation manifest: $manifest"
