#!/usr/bin/env bash
#
# KKMC production for the ISR study.
#
# Generates e+e- -> q qbar + n gamma at the Z pole with KKMC (CEEX, DIZET 6.45,
# PYTHIA 6.202 hadronization) in three configurations and converts the ASCII
# event dump into the standard ISRStudy `Events` ROOT schema.
#
#   ISR_OFF     KeyISR=0  KeyINT=0   beam ISR off
#   ISR_ON      KeyISR=1  KeyINT=0   beam ISR on, ISR-FSR interference off
#   ISR_ON_IFI  KeyISR=1  KeyINT=2   beam ISR on, interference on (KKMC default)
#
# The ISR_OFF / ISR_ON pair is the ISR toggle that matches the PYTHIA, Herwig and
# Sherpa samples.  ISR_ON_IFI versus ISR_ON isolates the interference term, which
# no other generator in this study implements.
#
set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
KKMC_DIR="${KKMC_DIR:-/raid5/data/yjlee/ISR/KKMCee}"
EVENTS="${EVENTS:-1000000}"
ECM="${ECM:-91.1876}"
OUTDIR="${OUTDIR:-}"
MODES="${MODES:-ISR_OFF,ISR_ON,ISR_ON_IFI}"
FORCE="${FORCE:-0}"

if [[ -z "$OUTDIR" ]]; then
  echo "Set OUTDIR explicitly, for example OUTDIR=/data2/yjlee/ISRsample/kkmc_20260921" >&2
  exit 2
fi
mkdir -p "$OUTDIR"

KKDUMP="$KKMC_DIR/ffbench/KKdump.exe"
if [[ ! -x "$KKDUMP" ]]; then
  echo "Missing KKdump.exe.  Build it with:" >&2
  echo "  cd $KKMC_DIR && ln -sfn dizet-6.45 dizet" >&2
  echo "  cd ffbench && make -f KKMakefile makflag && make -f KKMakefile makprod" >&2
  echo "  make -f KKMakefile EWtables && make -f KKMakefile KKdump.exe" >&2
  exit 1
fi

if [[ ! -x "$REPO_ROOT/real_isr_ntuple_producer" ]]; then
  "$REPO_ROOT/scripts/build_real_isr_tools.sh"
fi

write_input() {
  # write_input <file> <KeyISR> <KeyINT> <nevt>
  local f="$1" keyisr="$2" keyint="$3" nevt="$4"
  cat > "$f" <<EOF
********************** KKMC INPUT FOR THE ISR STUDY ****************************
  $(printf '%8d' "$nevt")      NEVT number of events   <- Input for the main program
----------cccccccccccccccccccccccccccccommenttttttttttttttttttttttttttttttttttt

BeginX
*indx_____data______ccccccccc0cccc__________General_____ccc0ccccccccc0ccccccccc0
    1     $(printf '%10s' "${ECM}e0")      CMSene =xpar( 1) Center of mass energy [GeV]
    2        0.000e0      DelEne =xpar( 2) Beam energy spread [GeV]
*     Hadronization/showering ON so that the final state is hadronic
   50              1      KeyHad=xpar(50)
*indx_____data______ccccccccc0ccccc________Process______ccc0ccccccccc0ccccccccc0
*     Quark final states only, to match e+e- -> q qbar in the other generators
  401              1      KFfin, d
  402              1      KFfin, u
  403              1      KFfin, s
  404              1      KFfin, c
  405              1      KFfin, b
*indx_____data______ccccccccc0ccccc________Radiation____ccc0ccccccccc0ccccccccc0
   20     $(printf '%10d' "$keyisr")      KeyISR=xpar(20)  beam ISR
   21              1      KeyFSR=xpar(21)  final-state radiation, same in all modes
   27     $(printf '%10d' "$keyint")      KeyINT=xpar(27)  ISR-FSR interference
   28              1      KeyGPS=xpar(28)  CEEX exponentiation
   29              1      KeyQSR=xpar(29)  photon emission from final quarks
*indx_____data______ccccccccc0ccccc________Miscel_______ccc0ccccccccc0ccccccccc0
    5              0      LevPri =xpar( 5)  PrintOut Level
    6              1      Ie1Pri =xpar( 6)
    7              0      Ie2Pri =xpar( 7)
********************************************************************************
EndX
EOF
}

run_mode() {
  local mode="$1" keyisr="$2" keyint="$3" seed="$4"
  local out="$OUTDIR/mc_KKMC424_${mode}.root"
  if [[ -f "$out" && "$FORCE" != "1" ]]; then
    echo "[reuse] $out"
    return
  fi
  local work="$KKMC_DIR/ffbench/ISRrun_${mode}"
  rm -rf "$work"
  mkdir -p "$work"
  ln -sfn ../../.KK2f_defaults "$work/.KK2f_defaults"
  ln -sfn ../../dizet "$work/dizet"
  printf '%10d\n%10d\n%10d\n' "$seed" 0 0 > "$work/iniseed"
  write_input "$work/pro.input" "$keyisr" "$keyint" "$EVENTS"

  echo "[run] KKMC $mode  KeyISR=$keyisr KeyINT=$keyint  $EVENTS events"
  ( cd "$work" && "$KKDUMP" ) > "$work/kkmc_run.log" 2>&1

  if [[ ! -s "$work/kkmc_events.dat" ]]; then
    echo "KKMC produced no events for $mode; see $work/kkmc_run.log" >&2
    exit 1
  fi

  "$REPO_ROOT/real_isr_ntuple_producer" \
    --mode kkmc \
    --kkmcInput "$work/kkmc_events.dat" \
    --nEvents "$EVENTS" \
    --sqrtS "$ECM" \
    --isrOn "$keyisr" \
    --generatorName "KKMC424_${mode}" \
    --generatorId 5 \
    --output "$out" > "$work/converter.log" 2>&1

  echo "[done] $out"
  gzip -f "$work/kkmc_events.dat"
}

contains() { [[ ",$MODES," == *",$1,"* ]]; }

contains ISR_OFF     && run_mode ISR_OFF     0 0 1500510
contains ISR_ON      && run_mode ISR_ON      1 0 1500511
contains ISR_ON_IFI  && run_mode ISR_ON_IFI  1 2 1500512

echo "[done] KKMC production in $OUTDIR"
