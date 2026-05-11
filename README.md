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

- `macros/reproduce_isr_plot.C`: legacy C++ ROOT macro for the original
  slide-style ISR correction figure.  It is not used for the current
  real-generator production.
- `src/real_isr_ntuple_producer.cc`: standalone real-generator ROOT ntuple
  producer.  It directly runs PYTHIA/PYTHIA+Vincia and converts Herwig/Sherpa
  HepMC streams into the same `Events` TTree schema.
- `scripts/run_real_isr_production.sh`: real standalone generator production
  wrapper for PYTHIA 8.315, PYTHIA 8.315 (Vincia), Herwig 7.3.0, and Sherpa
  3.0.3.
- `scripts/run_real_isr_production_10worker.sh`: worker-capped launcher used
  for the 100k validation and larger parallel productions.
- `scripts/run_isr_production.sh`: legacy wrapper for plot-only reruns from
  the older macro workflow.
- `scripts/test_isr_macro.sh`: smoke test that generates tiny temporary samples,
  verifies ROOT tree entries, and checks plot creation.
- `cards/`: Herwig and Sherpa standalone cards for ISR ON/OFF.
- `docs/MC_GENERATION.md`: detailed production note, model status, links, sample
  sizes, and 20M statistics.
- `docs/REAL_GENERATOR_PRODUCTION.md`: current real-generator production status,
  version references, ISR ON/OFF definitions, sample statistics, and output paths.
- `docs/TREE_SCHEMA.md`: event and particle branch definitions.
- `docs/VALIDATION.md`: validation commands and checked production numbers.
- `results/`: small reference PNG/PDF outputs from the 20M production.

The multi-GB ROOT ntuples are intentionally not committed to git.

## Current Production

The latest generated 100k validation production is stored on the analysis
machine at:

```text
/data2/yjlee/ISRsample/real_100k_20260511
```

It contains eight ROOT files produced from real standalone generators:

- PYTHIA 8.315
- PYTHIA 8.315 with Vincia
- Herwig 7.3.0
- Sherpa 3.0.3
- ISR ON and ISR OFF for each configuration
- `100,000` events per file in the validation refresh

Those 100k files were produced before the generator-setting requirements were
updated.  The code and cards now define the next production as:

- PYTHIA/PYTHIA-Vincia: `PDF:lepton`, `PartonLevel:ISR`, and
  `SpaceShower:QEDshowerByL` are all toggled together; `gm/Z -> q qbar` is
  selected with `WeakSingleBoson:ffbar2gmZ`, `WeakZ0:gmZmode = 0`, and
  `23:onIfAny = 1 2 3 4 5`.
- Herwig 7.3.0: the baseline is `snippets/EECollider.in` plus
  `MEee2gZ2qq`, `91.1876*GeV`, and `Interactions QCD`; the comparison sample
  is explicitly labeled `QEDshower` and uses the installed Herwig spelling
  `Interactions QEDQCD` for QCD+QED showering.  The literal `QCDandQED` token
  is rejected by the installed Herwig binary.
- Sherpa 3.0.3: nominal ISR uses `PDF_LIBRARY: PDFESherpa`; YFS ISR is a
  separately labeled alternative and is not mixed with `PDFESherpa`.

The old fallback samples and the old 5M archive were removed from
`/data2/yjlee/ISRsample` to avoid confusion.

## Quick Start

Run the checked-in smoke test:

```bash
./scripts/test_isr_macro.sh
```

Build the real-generator producer:

```bash
./scripts/build_real_isr_tools.sh
```

Run the real-generator validation production:

```bash
EVENTS=100000 MAX_WORKERS=15 OUTDIR=/data2/yjlee/ISRsample/real_100k_20260511 FORCE=0 ./scripts/run_real_isr_production_10worker.sh
```

Regenerate the real-generator plots:

```bash
ISR_REAL_DIR=/data2/yjlee/ISRsample/real_100k_20260511 root -l -b -q macros/plot_real_isr_results.C
```

The old slide macro still runs with:

```bash
root -l -b -q macros/reproduce_isr_plot.C
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
mc_Herwig730_ISR_ON.root
mc_Herwig730_ISR_OFF.root
mc_Pythia8315_ISR_ON.root
mc_Pythia8315_ISR_OFF.root
mc_Sherpa303_ISR_ON.root
mc_Sherpa303_ISR_OFF.root
mc_Pythia8315_Vincia_ISR_ON.root
mc_Pythia8315_Vincia_ISR_OFF.root
results/real_isr_correction_double_ratio.png
results/real_isr_on_thrust_vs_aleph.png
results/real_isr_photon_energy_theta_spectra.png
results/isr_correction_studies_reproduction.png
```

Each ROOT file contains a `TTree` named `Events`.
