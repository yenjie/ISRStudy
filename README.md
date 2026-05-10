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

- `macros/reproduce_isr_plot.C`: C++ ROOT macro that can generate MC-like
  toy fallback samples, analyze ISR ON/OFF pairs, and draw the original
  slide-style ISR correction figure.
- `src/real_isr_ntuple_producer.cc`: standalone real-generator ROOT ntuple
  producer.  It directly runs PYTHIA/PYTHIA+Vincia and converts Herwig/Sherpa
  HepMC streams into the same `Events` TTree schema.
- `scripts/run_real_isr_production.sh`: real standalone generator production
  wrapper for PYTHIA 8.315, PYTHIA 8.315 (Vincia), Herwig 7.3.0, and Sherpa
  3.0.3.
- `scripts/run_isr_production.sh`: production wrapper for full 20M-per-sample
  generation or plot-only reruns from the older macro workflow.
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

The current real-generator validation production is stored on the analysis
machine at:

```text
/data2/yjlee/ISRsample/real_20260510
```

It contains eight ROOT files produced from real standalone generators:

- PYTHIA 8.315
- PYTHIA 8.315 with Vincia
- Herwig 7.3.0
- Sherpa 3.0.3
- ISR ON and ISR OFF for each configuration
- `20,000` events per file in the validation refresh

The old toy samples and the old 5M toy archive were removed from
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
EVENTS=20000 OUTDIR=/data2/yjlee/ISRsample/real_20260510 TARGETS=all FORCE=1 ./scripts/run_real_isr_production.sh
```

Regenerate the real-generator plots:

```bash
root -l -b -q macros/plot_real_isr_results.C
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
```

Each ROOT file contains a `TTree` named `Events`.
