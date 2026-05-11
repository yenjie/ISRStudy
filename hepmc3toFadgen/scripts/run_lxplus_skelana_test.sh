#!/usr/bin/env bash
set -euo pipefail

PKG_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
REPO_ROOT="$(cd "$PKG_ROOT/.." && pwd -P)"

LOCAL_FADGEN="${LOCAL_FADGEN:-/raid5/data/yjlee/ISR/hepmc3toFadgen_smoke/fad/sherpa_OFF.fadgen}"
LOCAL_EXAMPLES="${LOCAL_EXAMPLES:-/raid5/data/yjlee/ISR/delphi_examples}"
REMOTE="${CERN_REMOTE:-lxplus}"
REMOTE_BASE="${REMOTE_BASE:-/tmp/yjlee_sherpa_fadgen_skelana_test}"
REMOTE_WORK="${REMOTE_BASE}/work_$(date +%Y%m%d_%H%M%S)"
NEVMAX="${NEVMAX:-3}"
VERSION="${VERSION:-v94c}"
LABO="${LABO:-CERN}"
NRUN="${NRUN:-1000}"
EBEAM="${EBEAM:-45.5938}"

if [[ ! -s "$LOCAL_FADGEN" ]]; then
  echo "Missing LOCAL_FADGEN: $LOCAL_FADGEN" >&2
  exit 2
fi

if ! ssh -O check "$REMOTE" >/dev/null 2>&1; then
  cat >&2 <<EOF
No active SSH ControlMaster for $REMOTE.

Open one in your terminal and complete CERN 2FA:

  ssh -MNf $REMOTE

Then rerun this script.
EOF
  exit 2
fi

if [[ ! -d "$LOCAL_EXAMPLES/dump" ]]; then
  git clone https://gitlab.cern.ch/delphi/examples.git "$LOCAL_EXAMPLES"
fi

tmp_tar="$(mktemp /tmp/delphi_dump_example.XXXXXX.tar.gz)"
tar -C "$LOCAL_EXAMPLES" -czf "$tmp_tar" dump
trap 'rm -f "$tmp_tar"' EXIT

ssh "$REMOTE" "mkdir -p '$REMOTE_WORK'"
scp "$LOCAL_FADGEN" "$REMOTE:$REMOTE_WORK/input.fadgen"
scp "$tmp_tar" "$REMOTE:$REMOTE_WORK/dump.tar.gz"

ssh "$REMOTE" bash -s -- "$REMOTE_WORK" "$NEVMAX" "$VERSION" "$LABO" "$NRUN" "$EBEAM" <<'REMOTE_TEST'
set -euo pipefail

work="$1"
nevmax="$2"
version="$3"
labo="$4"
nrun="$5"
ebeam="$6"

cd "$work"
tar -xzf dump.tar.gz

{
  echo "# Sherpa FADGEN to DELPHI skelana validation"
  echo
  echo "- host: $(hostname)"
  echo "- date: $(date -Iseconds)"
  echo "- workdir: $work"
  echo "- input: $work/input.fadgen"
  echo "- version: $version"
  echo "- beam energy: $ebeam"
  echo "- nevmax: $nevmax"
  echo
} > validation_manifest.md

if [[ ! -r /cvmfs/delphi.cern.ch/setup.sh ]]; then
  echo "Missing /cvmfs/delphi.cern.ch/setup.sh" | tee -a validation_manifest.md
  exit 3
fi

. /cvmfs/delphi.cern.ch/setup.sh

{
  echo "## Environment"
  echo
  for cmd in runsim dellib nypatchy cernlib gfortran; do
    printf -- "- %-10s " "$cmd"
    command -v "$cmd" || true
  done
  echo "- DELPHI_DATA_ROOT: ${DELPHI_DATA_ROOT:-unset}"
  echo
} >> validation_manifest.md

echo "[runsim] testing external FADGEN input"
set +e
runsim -VERSION "$version" -LABO "$labo" -NRUN "$nrun" -EBEAM "$ebeam" -gext input.fadgen -NEVMAX "$nevmax" > runsim.log 2> runsim.err
runsim_status=$?
set -e

{
  echo "## runsim"
  echo
  echo "- exit status: $runsim_status"
  for f in simana.fadsim simana.fadana simana.sdst simana.xsdst; do
    if [[ -e "$f" ]]; then
      echo "- produced: $f ($(stat -c%s "$f") bytes)"
    fi
  done
  echo
} >> validation_manifest.md

if [[ "$runsim_status" -ne 0 ]]; then
  echo "runsim failed; skelana was not run" | tee -a validation_manifest.md
  exit "$runsim_status"
fi

dst=""
if [[ -s simana.sdst ]]; then
  dst="simana.sdst"
elif [[ -s simana.xsdst ]]; then
  dst="simana.xsdst"
else
  echo "runsim completed but no simana.sdst/simana.xsdst was produced" | tee -a validation_manifest.md
  exit 4
fi

cd "$work/dump"
printf '*\nFILE = %s/%s\n*\n' "$work" "$dst" > PDLINPUT

echo "[skelana] running dump example on $dst"
set +e
./dump.sh > dump_driver.log 2> dump_driver.err
dump_status=$?
set -e

{
  echo "## skelana dump"
  echo
  echo "- input DST: $work/$dst"
  echo "- exit status: $dump_status"
  if [[ -s dump.log ]]; then
    echo "- dump.log: $work/dump/dump.log"
    echo "- CHECK lines: $(grep -c 'CHECK:' dump.log || true)"
  fi
  echo
} >> "$work/validation_manifest.md"

exit "$dump_status"
REMOTE_TEST

echo "[done] remote workdir: $REMOTE_WORK"
echo "Fetch manifest with:"
echo "  scp $REMOTE:$REMOTE_WORK/validation_manifest.md ./"
