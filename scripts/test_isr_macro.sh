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

EVENTS="${EVENTS:-2000}"
OUTDIR="${OUTDIR:-/tmp/isrstudy_smoke_${USER}_$$}"
export EVENTS OUTDIR
mkdir -p "${OUTDIR}"

echo "Running ISRStudy smoke test"
echo "ROOT_BIN=${ROOT_BIN}"
echo "EVENTS=${EVENTS}"
echo "OUTDIR=${OUTDIR}"

cd "${REPO_DIR}"
"${ROOT_BIN}" -l -b -q \
  -e ".L macros/reproduce_isr_plot.C" \
  -e "nEvents=${EVENTS}; regenerateSamples=true; gOutputDir=\"${OUTDIR}\"; reproduce_isr_plot();"

"${ROOT_BIN}" -l -b -q -e '
const char* out = gSystem->Getenv("OUTDIR");
Long64_t expected = atoll(gSystem->Getenv("EVENTS"));
const char* files[] = {
  "mc_Herwig715_ISR_ON.root",
  "mc_Herwig715_ISR_OFF.root",
  "mc_Pythia8317_ISR_ON.root",
  "mc_Pythia8317_ISR_OFF.root",
  "mc_Sherpa226_ISR_ON.root",
  "mc_Sherpa226_ISR_OFF.root",
  "mc_Pythia8317_Vincia_ISR_ON.root",
  "mc_Pythia8317_Vincia_ISR_OFF.root"
};
bool ok = true;
for (auto name : files) {
  TString path = TString::Format("%s/%s", out, name);
  TFile file(path, "READ");
  TTree* tree = (TTree*)file.Get("Events");
  Long64_t entries = tree ? tree->GetEntries() : -1;
  std::cout << path << " entries=" << entries << std::endl;
  if (entries != expected) ok = false;
}
TString png = TString::Format("%s/isr_correction_studies_reproduction.png", out);
TString pdf = TString::Format("%s/isr_correction_studies_reproduction.pdf", out);
if (gSystem->AccessPathName(png) || gSystem->AccessPathName(pdf)) ok = false;
if (!ok) gSystem->Exit(2);
'

echo "Smoke test passed"
if [[ "${KEEP_ISRTEST:-0}" != "1" ]]; then
  rm -rf "${OUTDIR}"
else
  echo "Kept test output in ${OUTDIR}"
fi
