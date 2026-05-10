# ISR Z-Pole MC Generation Note

Generated/updated: 2026-05-10 11:34 EDT

## Goal

Produce reusable MC-like ROOT ntuples for studying ISR effects in archived ALEPH-style
e+e- events at the Z pole, sqrt(s) = 91.2 GeV. The immediate analysis target is the
ISR correction factor

```text
C_ISR(T) = ISR OFF / ISR ON
```

as a function of thrust T. The same files are intended to be reusable for later
studies of visible energy, charged-particle multiplicity, ISR photon kinematics,
event-shape observables, generator/process categories, particle-level quantities,
and parton-level quantities.

## Generation Command

The 20M-event-per-sample MC production and plot were generated from:

```bash
/snap/root-framework/954/usr/local/bin/root -l -b -q -e '.L macros/reproduce_isr_plot.C' -e 'nEvents=20000000; regenerateSamples=true; gOutputDir="/data2/yjlee/ISRsample"; reproduce_isr_plot();'
```

The requested command remains:

```bash
root -l -b -q macros/reproduce_isr_plot.C
```

On this machine, the direct ROOT binary was used because the Snap wrapper has path
visibility issues in some mounted analysis directories. The macro default is now
`nEvents = 20000000` and `gOutputDir = "/data2/yjlee/ISRsample"`, so the requested
command will write future samples and plot outputs under `/data2/yjlee/ISRsample`
unless edited.

## Software Versions

- ROOT: 6.36.04, built for `linuxx8664gcc` on 2025-08-27.
- Macro: `macros/reproduce_isr_plot.C`
- Macro SHA256:

```text
2e09ecb2d91d36c61a0e4da6c5cf7a89ccf2225e11581ce4b56581963d8c1e15  macros/reproduce_isr_plot.C
```

## Generator Model Used For These Files

The generated files in this directory were produced with the macro's labeled
fallback toy e+e- Z-pole ISR generator.

Reason: a Pythia8 installation is partially visible, but not valid for the macro's
configured Pythia path:

```bash
pythia8-config --version
# 8.315

pythia8-config --cxxflags --libs
# Error: cannot find valid configuration for Pythia 8
```

Therefore this production did not use external Pythia events. It used the internal
toy generator in `macros/reproduce_isr_plot.C`, with deterministic `TRandom3` seeds and
the same output branch structure as the Pythia-backed path. The toy generator
simulates e+e- -> hadrons at sqrt(s) = 91.2 GeV, creates parton-level and
particle-level records, injects explicit ISR photons for ISR ON samples, and
stores event-shape observables computed from final-state visible particles.

The toy model has no external web model page. Its effective version is the macro
checksum above. The generated samples are four configured toy-generator variants
named after Herwig 7.1.5, PYTHIA 8.317, Sherpa 2.2.6, and PYTHIA 8.317 Vincia
for correction-shape comparison; those labels are not claims that external Herwig,
Pythia, Sherpa, or Vincia executables generated these files.

## External Model Path Documented For Reproducibility

The macro can use Pythia8 when Pythia is fully configured. The official Pythia
website describes Pythia as a general-purpose high-energy physics event generator
including hard/soft interactions, showers, fragmentation, and decays:

- Official Pythia site: https://pythia.org/
- Pythia release/download page: https://www.pythia.org/releases/
- Pythia documentation page: https://pythia.org/documentation/
- Local version detected here: Pythia 8.315
- Online manual for Pythia 8.315: https://pythia.org/manuals/pythia8315/Welcome.html
- Main Pythia 8.3 reference: https://arxiv.org/abs/2203.11601
- SciPost Pythia 8.3 manual page: https://scipost.org/10.21468/SciPostPhysCodeb.8
- DOI for the Pythia 8.3 manual: https://doi.org/10.21468/SciPostPhysCodeb.8

Web sources checked on 2026-05-09: the official Pythia release page lists
PYTHIA 8.317 as the latest version, and the documentation page includes online
manuals for PYTHIA 8.317, 8.316, 8.315, and older releases. The locally detected
8.315 is therefore not the latest upstream version, and it was not usable for this
production because `pythia8-config --cxxflags --libs` failed.

## Output Samples

Each file contains a `TTree` named `Events`. The completed 20M production is stored
under `/data2/yjlee/ISRsample`; all eight top-level files were verified to contain
exactly 20,000,000 entries. The previous 5M production was moved to
`/data2/yjlee/ISRsample/5M_20260509` and occupies about 32G. The full sample area,
including the 5M archive, occupies about 156G.

| File | Entries | Size |
|---|---:|---:|
| `mc_Herwig715_ISR_ON.root` | 20,000,000 | 16G |
| `mc_Herwig715_ISR_OFF.root` | 20,000,000 | 15G |
| `mc_Pythia8317_ISR_ON.root` | 20,000,000 | 16G |
| `mc_Pythia8317_ISR_OFF.root` | 20,000,000 | 16G |
| `mc_Sherpa226_ISR_ON.root` | 20,000,000 | 17G |
| `mc_Sherpa226_ISR_OFF.root` | 20,000,000 | 16G |
| `mc_Pythia8317_Vincia_ISR_ON.root` | 20,000,000 | 16G |
| `mc_Pythia8317_Vincia_ISR_OFF.root` | 20,000,000 | 16G |

ISR ON summary statistics recomputed from the 20M files:

| Configuration | Emit fraction | Mean ISR photons/event | Mean ISR photon energy/event | Mean ISR photon energy/emitting event |
|---|---:|---:|---:|---:|
| Herwig 7.1.5 equivalent | 43.23% | 0.473 | 1.857 GeV | 4.294 GeV |
| PYTHIA 8.317 equivalent | 47.01% | 0.518 | 2.167 GeV | 4.609 GeV |
| Sherpa 2.2.6 equivalent | 50.77% | 0.564 | 2.505 GeV | 4.934 GeV |
| PYTHIA 8.317 Vincia equivalent | 46.05% | 0.507 | 2.087 GeV | 4.532 GeV |

The final slide-style outputs are:

```text
isr_correction_studies_reproduction.png
isr_correction_studies_reproduction.pdf
```

## Event Record

The `Events` tree stores event-level observables:

```text
eventId, generatorId, generatorName, isrOn, processId, processName,
sqrtS, weight, visibleEnergy, missingEnergy, thrust, thrustMajor,
thrustMinor, sphericity, cParameter, nFinalParticles,
nChargedFinalParticles, nPartons, nISRPhotons, totalISRPhotonEnergy
```

It also stores particle-level vectors:

```text
pdgId, status, mother1, mother2, daughter1, daughter2,
px, py, pz, energy, mass, charge,
isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton
```

Default observable selection excludes neutrinos and excludes ISR photons from
visible energy and event shapes. Charged-only selections and ISR-photon-inclusive
selections are configurable through `SelectionConfig`.

## Analysis Chain

The final ISR correction plot is made by reading the ROOT ntuples, filling matched
ISR OFF and ISR ON histograms in thrust bins, dividing them bin-by-bin, propagating
statistical uncertainties, and converting the ratios to `TGraphErrors`. The plot
does not use hard-coded final graph points.
