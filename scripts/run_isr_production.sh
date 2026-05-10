#!/usr/bin/env bash
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"

if [[ -z "${ROOT_BIN:-}" ]]; then
  if [[ -x /snap/root-framework/954/usr/local/bin/root ]]; then
    ROOT_BIN=/snap/root-framework/954/usr/local/bin/root
  else
    ROOT_BIN=root
  fi
fi
if ! command -v "${ROOT_BIN}" >/dev/null 2>&1; then
  echo "ERROR: ROOT executable not found. Set ROOT_BIN=/path/to/root." >&2
  exit 1
fi

EVENTS="${EVENTS:-20000000}"
OUTPUT_DIR="${OUTPUT_DIR:-/data2/yjlee/ISRsample}"
REGENERATE="${REGENERATE:-true}"

mkdir -p "${OUTPUT_DIR}"

echo "ROOT_BIN=${ROOT_BIN}"
echo "EVENTS=${EVENTS}"
echo "OUTPUT_DIR=${OUTPUT_DIR}"
echo "REGENERATE=${REGENERATE}"

cd "${REPO_DIR}"
"${ROOT_BIN}" -l -b -q \
  -e ".L macros/reproduce_isr_plot.C" \
  -e "nEvents=${EVENTS}; regenerateSamples=${REGENERATE}; gOutputDir=\"${OUTPUT_DIR}\"; reproduce_isr_plot();"
