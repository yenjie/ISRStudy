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
  samples, analyze ISR ON/OFF pairs, and draw the ISR correction figure.
- `scripts/run_isr_production.sh`: production wrapper for full 20M-per-sample
  generation or plot-only reruns from existing samples.
- `scripts/test_isr_macro.sh`: smoke test that generates tiny temporary samples,
  verifies ROOT tree entries, and checks plot creation.
- `docs/MC_GENERATION.md`: detailed production note, model status, links, sample
  sizes, and 20M statistics.
- `docs/TREE_SCHEMA.md`: event and particle branch definitions.
- `docs/VALIDATION.md`: validation commands and checked production numbers.
- `results/`: small reference PNG/PDF outputs from the 20M production.

The multi-GB ROOT ntuples are intentionally not committed to git.

## Current Production

The completed 20M production is stored on the analysis machine at:

```text
/data2/yjlee/ISRsample
```

It contains eight ROOT files:

- four generator-equivalent configurations
- ISR ON and ISR OFF for each configuration
- `20,000,000` events per file
- `160,000,000` total generated events

The previous 5M production is archived at:

```text
/data2/yjlee/ISRsample/5M_20260509
```

## Quick Start

Run the checked-in smoke test:

```bash
./scripts/test_isr_macro.sh
```

Run a plot-only update from existing 20M samples:

```bash
REGENERATE=false ./scripts/run_isr_production.sh
```

Run the full 20M production:

```bash
EVENTS=20000000 OUTPUT_DIR=/data2/yjlee/ISRsample REGENERATE=true ./scripts/run_isr_production.sh
```

The bare ROOT command also works from the repository root:

```bash
root -l -b -q macros/reproduce_isr_plot.C
```

By default, the macro writes samples and plots to `/data2/yjlee/ISRsample`.

## Generator Status

The macro can use Pythia8 only when a complete, linkable Pythia8 installation is
available.  On the machine used for this production, `pythia8-config --version`
reported `8.315`, but `pythia8-config --cxxflags --libs` failed.  The 20M files
therefore use the macro's labeled fallback toy `e+e-` `Z`-pole ISR generator.

The labels `Herwig 7.1.5`, `PYTHIA 8.317`, `Sherpa 2.2.6`, and
`PYTHIA 8.317 (Vincia)` are generator-equivalent model configurations for
correction-shape comparison.  They are not claims that external Herwig, PYTHIA,
Sherpa, or Vincia executables produced these files.

## Output Files

The production writes:

```text
mc_Herwig715_ISR_ON.root
mc_Herwig715_ISR_OFF.root
mc_Pythia8317_ISR_ON.root
mc_Pythia8317_ISR_OFF.root
mc_Sherpa226_ISR_ON.root
mc_Sherpa226_ISR_OFF.root
mc_Pythia8317_Vincia_ISR_ON.root
mc_Pythia8317_Vincia_ISR_OFF.root
isr_correction_studies_reproduction.png
isr_correction_studies_reproduction.pdf
```

Each ROOT file contains a `TTree` named `Events`.

