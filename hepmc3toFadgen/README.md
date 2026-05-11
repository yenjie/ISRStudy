# hepmc3toFadgen

Sherpa-to-FADGEN production and validation package for the ISR study.

The production chain is intentionally staged:

```text
Sherpa 3.0.3 -> HepMC ASCII -> FADGEN exchange file -> validation ROOT ntuple
```

The converter in this directory can write two formats:

- `DELSIM_UNFORMATTED_LUJETS`: native DELPHI external-generator input for
  `runsim -gext`, using the same Fortran unformatted LUJETS record layout as
  the DELPHI Pythia example `WINWRITE` routine.
- `FADGEN_EXCHANGE_ASCII_V1`: a loss-preserving generator-level exchange record
  used for debugging and format audits.

The Sherpa production wrapper writes native `DELSIM_UNFORMATTED_LUJETS` by
default.  Use `FADGEN_OUTPUT_KIND=exchange` only for human-readable converter
debugging.

## Build

```bash
cd /raid5/data/yjlee/ISR/ISRStudy/hepmc3toFadgen
./scripts/build.sh
```

The build uses the system HepMC3 headers and library.

## External References

- Sherpa ISR structure functions: the Sherpa 3.0 manual documents
  `PDF_LIBRARY: PDFESherpa` as the built-in electron structure function,
  `ISR_E_ORDER`, `ISR_E_SCHEME`, and `PDF_LIBRARY: None` for fixed beam
  energy: <https://sherpa-team.gitlab.io/sherpa/v3.0.1/manual/parameters/isr.html>
- Sherpa YFS lepton-collision ISR: the Sherpa QED-corrections manual documents
  the `YFS` settings block and `MODE: ISR` for lepton-lepton scattering:
  <https://sherpa-team.gitlab.io/sherpa/master/manual/parameters/qed-corrections.html>
- HepMC3 event record: the converter reads HepMC2-style and HepMC3 ASCII
  through HepMC3 readers; the library and release tags are documented at
  <https://gitlab.cern.ch/hepmc/HepMC3>.
- DELPHI Open Data external-generator path: CERN Open Data documents running
  external generator output through DELPHI simulation and then analysing the
  result with the `skelana` skeleton analysis framework:
  <https://opendata.cern/docs/delphi-getting-started> and
  <https://opendata.cern/docs/delphi-guide-analysis>.

## Quick Smoke Test

Run 10 events for all Sherpa ISR modes and produce FADGEN exchange files plus
ROOT validation ntuples:

```bash
EVENTS=10 \
OUTDIR=/data2/yjlee/ISRsample/sherpa_fadgen_smoke \
FORCE=1 \
./scripts/run_sherpa_to_fadgen.sh
```

Outputs are written as:

```text
$OUTDIR/
  cards/
  hepmc/
  fad/                 # native DELSIM unformatted LUJETS by default
  root_validation/
  logs/
  metadata/
```

Validate that the native FADGEN generator-level payload matches the original
HepMC final-state content:

```bash
hepmc3toFadgen/bin/validate_fadgen_vs_hepmc \
  --hepmc /data2/yjlee/ISRsample/sherpa_fadgen_1k/hepmc/sherpa_OFF.hepmc \
  --fadgen /data2/yjlee/ISRsample/sherpa_fadgen_1k/fad/sherpa_OFF.fadgen \
  --format hepmc2 \
  --report /data2/yjlee/ISRsample/sherpa_fadgen_1k/metadata/fadgen_vs_hepmc_validation.md
```

The validator decodes each native DELSIM LUJETS Fortran record and compares
the event-by-event final-particle multiplicity, PDG id, and `px, py, pz, E,
mass` against the source HepMC record, allowing only float-rounding tolerance.

## Full Production

```bash
EVENTS=3000000 \
OUTDIR=/data2/yjlee/ISRsample/sherpa_fadgen_3M_YYYYMMDD \
TARGETS=OFF,PDFESherpa,YFS \
MAX_WORKERS=1 \
FORCE=0 \
./scripts/run_sherpa_to_fadgen.sh
```

Run with `MAX_WORKERS=1` first.  Sherpa integration libraries and output files
are simpler to debug one mode at a time.

## DELPHI skelana Validation

The final integration test is not just that this converter can write a file.
The goal is to show that the downstream DELPHI Open Data analysis path can
consume the generated output.

1. Run `skelana` on an official DELPHI Open Data sample to validate the local
   DELPHI runtime.
2. Run the Sherpa-derived FADGEN/FAD output through the same `skelana` path.
3. Compare event counts and basic observables with the ROOT validation ntuple.

Configure the harness with environment variables:

```bash
SKELEANA_CMD=/path/to/skelana \
DELPHI_REFERENCE_INPUT=/path/to/official/delphi/sample \
FADGEN_INPUT=/data2/yjlee/ISRsample/sherpa_fadgen_smoke/fad/sherpa_OFF.fadgen \
OUTDIR=/data2/yjlee/ISRsample/sherpa_fadgen_smoke/skelana \
./scripts/run_skeleana_validation.sh
```

If `skelana` needs more arguments, pass them as strings:

```bash
SKELEANA_REFERENCE_ARGS="-i /path/to/reference -o reference.root" \
SKELEANA_FADGEN_ARGS="-i /path/to/fadgen -o sherpa_fadgen.root" \
./scripts/run_skeleana_validation.sh
```

For the real DELPHI external-generator chain on LXPLUS, use:

```bash
LOCAL_FADGEN=/raid5/data/yjlee/ISR/hepmc3toFadgen_smoke/fad/sherpa_OFF.fadgen \
./scripts/run_lxplus_skelana_test.sh
```

This requires an already authenticated CERN SSH ControlMaster (`ssh -MNf
lxplus`).  The script copies the candidate FADGEN file to LXPLUS, sources
`/cvmfs/delphi.cern.ch/setup.sh`, runs `runsim -gext input.fadgen`, and then
runs the DELPHI `dump` skelana example on the resulting `simana.sdst` or
`simana.xsdst`.

Validation status:

- The old text exchange output starts `runsim`, but `SXRDLU` reports a read
  error and no events are processed.  It is not native FADGEN.
- The native `DELSIM_UNFORMATTED_LUJETS` output passed the LXPLUS chain on
  2026-05-11: `runsim -gext` processed 3 events and wrote 3 DST records, then
  the DELPHI `dump` skelana example read `simana.sdst` and produced 758
  `CHECK:` lines.

The local copy of that validation is stored in:

```text
/raid5/data/yjlee/ISR/hepmc3toFadgen_native_skelana_validation
```

## Native DELSIM LUJETS Format

The native output follows the DELPHI example in
`gitlab.cern.ch/delphi/examples/pythia/delpyt.for`, routine `WINWRITE`:

- one Fortran unformatted sequential record per event
- 32-bit record markers, as written/read by gfortran on x86_64
- payload: `N`, then for each LUJETS particle `K(1:5), P(1:5), V(1:5)`
- `K` entries are 32-bit integers
- `P` and `V` entries are 32-bit floating point values
- final HepMC particles are written as stable LUJETS particles with `K(1)=1`
- a final generator-comment line is appended with `K(1)=21`

## FADGEN Exchange Format

Each event block contains:

- event number
- event weight
- HepMC momentum and length units
- every particle in the HepMC event
- PDG id and generator status
- first/last mother and daughter indices
- four-vector and production vertex
- a conservative beam-collinear ISR-photon candidate flag used only for
  diagnostics

The exact DELPHI FAD/FADGEN binary/native layout still needs to be supplied by
the local FADGEN installation or DELPHI documentation.  Until then, this
package provides a deterministic, testable bridge from Sherpa to a stable
intermediate record plus the `skelana` validation hooks needed to close the
loop.
