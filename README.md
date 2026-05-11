# ISRStudy

Reusable ROOT-based ISR correction study for archived ALEPH-style
`e+e-` events at the `Z` pole, `sqrt(s) = 91.2 GeV`.

The main observable is

```text
C_ISR(T) = N_ISR_OFF(T) / N_ISR_ON(T)
```

computed from generated ROOT ntuples in bins of thrust.  The same ntuples store
enough event-record information to recompute visible energy, charged-particle
multiplicity, ISR photon kinematics, event-shape variables, generator categories,
particle-level observables, and parton-level observables.

## Repository Contents

- `src/real_isr_ntuple_producer.cc`: standalone real-generator ROOT ntuple
  producer.  It directly runs PYTHIA/PYTHIA+Vincia and converts Herwig/Sherpa
  HepMC streams into the same `Events` TTree schema.
- `scripts/run_real_isr_production.sh`: real standalone generator production
  wrapper for PYTHIA 8.315, PYTHIA 8.315 (Vincia), Herwig 7.3.0, and Sherpa
  3.0.3.
- `scripts/run_real_isr_production_10worker.sh`: worker-capped launcher used
  for the 3M parallel production.
- `macros/plot_real_isr_results.C`: real-generator plotting macro for ISR
  corrections, ALEPH thrust comparisons, ratio panels, ISR photons, and visible
  energy.
- `macros/write_real_isr_stats.C`: standalone ROOT statistics writer for the
  scalar sample-summary CSV.
- `src/make_endpoint_diagnostics.cc`: derived per-event diagnostic ntuple
  producer for checking whether high-thrust endpoint behavior is correlated
  with low visible mass, longitudinal boost, or explicit ISR photon activity.
- `macros/plot_endpoint_diagnostics.C`: 2D endpoint-correlation plotter and
  summary-table writer for the derived diagnostic ntuples.
- `hepmc3toFadgen/`: Sherpa 3.0.3 to HepMC to FADGEN-exchange bridge, with a
  documented converter and a DELPHI `skeleana` validation harness.
- `cards/`: Herwig and Sherpa standalone cards for ISR ON/OFF.
- `docs/REAL_GENERATOR_PRODUCTION.md`: current real-generator production status,
  version references, ISR ON/OFF definitions, sample statistics, and output paths.
- `docs/TREE_SCHEMA.md`: event and particle branch definitions.

The multi-GB ROOT ntuples and derived plot files are intentionally not committed
to git.  The current outputs live under `/data2/yjlee/ISRsample/real_3M_20260511`
and are copied to Overleaf when needed.

## Current Production

The current 3M real-generator production is stored on the analysis machine at:

```text
/data2/yjlee/ISRsample/real_3M_20260511
```

It contains nine ROOT files produced from real standalone generators:

- PYTHIA 8.315
- PYTHIA 8.315 with Vincia
- Herwig 7.3.0 baseline and QED-shower comparison
- Sherpa 3.0.3 baseline, PDFESherpa ISR, and YFS ISR alternative
- `3,000,000` events per generated file

The generator-setting definitions are:

- PYTHIA/PYTHIA-Vincia: `PDF:lepton`, `PartonLevel:ISR`, and
  `SpaceShower:QEDshowerByL` are all toggled together; `gm/Z -> q qbar` is
  selected with `WeakSingleBoson:ffbar2gmZ`, `WeakZ0:gmZmode = 0`, and
  `23:onIfAny = 1 2 3 4 5`.
- Herwig 7.3.0: the baseline is `snippets/EECollider.in` plus
  `MEee2gZ2qq`, `91.1876*GeV`, and `Interactions QCD`; the comparison sample
  is explicitly labeled `QEDshower` and uses the installed Herwig spelling
  `Interactions QEDQCD` for QCD+QED showering.  The literal `QCDandQED` token
  is rejected by the installed Herwig binary.  The HepMC `PrintEvent` setting
  follows the requested `EVENTS` value, avoiding the earlier 1M cap.
- Sherpa 3.0.3: nominal ISR uses `PDF_LIBRARY: PDFESherpa`; YFS ISR is a
  separately labeled alternative and is not mixed with `PDFESherpa`.

The old fallback samples and the old 5M archive were removed from
`/data2/yjlee/ISRsample` to avoid confusion.

Audit caveat for the existing 3M ROOT files: the current plots and CSV recompute
visible energy from particle vectors with the `|eta| < 1.74` definition.  The
stored scalar `visibleEnergy` branch in those already-produced files predates
that definition, and the Sherpa `weight` branch contains HepMC bookkeeping
weights that are not used by the current unweighted plots.  The producer has
been patched for future HepMC conversions to write `weight = 1.0`.

## Quick Start

Build the real-generator producer:

```bash
./scripts/build_real_isr_tools.sh
```

Load the ROOT 6.34.04 environment used for production:

```bash
cd /raid5/root/root-v6.34.04/root/bin
. ./thisroot.sh
cd /raid5/data/yjlee/ISR/ISRStudy
```

Run or reuse the 3M real-generator production:

```bash
EVENTS=3000000 MAX_WORKERS=15 OUTDIR=/data2/yjlee/ISRsample/real_3M_20260511 FORCE=0 ./scripts/run_real_isr_production_10worker.sh
```

Regenerate the real-generator plots:

```bash
ISR_REAL_DIR=/data2/yjlee/ISRsample/real_3M_20260511 root -l -b -q macros/plot_real_isr_results.C
```

Build and run the endpoint diagnostic pass for the ISR ON samples:

```bash
./scripts/build_real_isr_tools.sh
REAL_DIR=/data2/yjlee/ISRsample/real_3M_20260511 MAX_WORKERS=5 ./scripts/run_endpoint_diagnostics.sh
ISR_REAL_DIR=/data2/yjlee/ISRsample/real_3M_20260511 root -l -b -q macros/plot_endpoint_diagnostics.C
```

The same runner can produce the ISR OFF derived trees needed for correction
ratios:

```bash
REAL_DIR=/data2/yjlee/ISRsample/real_3M_20260511 OUTDIR=/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics_allfinal_3M SAMPLE_SET=OFF MAX_WORKERS=4 ./scripts/run_endpoint_diagnostics.sh
ISR_REAL_DIR=/data2/yjlee/ISRsample/real_3M_20260511 ISR_ENDPOINT_DIAG_DIR=/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics_allfinal_3M root -l -b -q macros/plot_thrust_definition_corrections.C
```

The diagnostic trees are written to:

```text
/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics_allfinal_3M
```

The endpoint correlation figures and summary tables are written to:

```text
/data2/yjlee/ISRsample/real_3M_20260511/results_allfinal_3M/endpoint_diagnostics
```

The thrust-definition `C_ISR(T)` comparison figures are written to:

```text
/data2/yjlee/ISRsample/real_3M_20260511/results_allfinal_3M/thrust_definition_corrections
```

For the current nominal 3M diagnostic, `T_lab_allFinal_including_ISR_photons`
is computed from all stable final-state particles, including neutrinos and
tagged ISR photons.  The diagnostic tree also stores the visible-particle
cross-checks: `T_lab_excluding_ISR_photons`, `T_lab_including_ISR_photons`,
and `T_visibleCM_excluding_ISR_photons`.

Regenerate only the ALEPH thrust comparison figures with `MC / ALEPH` ratio
panels:

```bash
ISR_REAL_DIR=/data2/yjlee/ISRsample/real_3M_20260511 root -l -b -q -e 'gROOT->LoadMacro("macros/plot_real_isr_results.C"); plot_real_aleph_ratio_figures();'
```

## Generator Status

The current files are not toy fallback samples.  They were produced from real
standalone generator paths:

- Direct PYTHIA 8.315 C++ generation, linked against `/usr/lib/libpythia8.so`.
- Direct PYTHIA 8.315 plus `Pythia8::Vincia`.
- Standalone Herwig 7.3.0 writing HepMC through a FIFO into the ROOT converter.
- Standalone Sherpa 3.0.3 writing HepMC through a FIFO into the ROOT converter.

See `docs/REAL_GENERATOR_PRODUCTION.md` for the exact ISR ON/OFF switches and
version references.

## Output Files

The production writes:

```text
mc_Herwig730_ISR_OFF.root
mc_Herwig730_QEDshower.root
mc_Pythia8315_ISR_ON.root
mc_Pythia8315_ISR_OFF.root
mc_Sherpa303_ISR_ON.root
mc_Sherpa303_ISR_OFF.root
mc_Sherpa303_ISR_YFS.root
mc_Pythia8315_Vincia_ISR_ON.root
mc_Pythia8315_Vincia_ISR_OFF.root
results/real_isr_correction_double_ratio.png
results/real_isr_on_thrust_vs_aleph.png
results/real_isr_on_thrust_vs_aleph_ratio.png
results/real_isr_off_thrust_vs_aleph.png
results/real_isr_off_thrust_vs_aleph_ratio.png
results/real_isr_photon_energy_theta_spectra.png
results/real_generator_sample_statistics.csv
results/isr_correction_studies_reproduction.png
```

Each ROOT file contains a `TTree` named `Events`.
