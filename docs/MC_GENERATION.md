# ISR Z-Pole MC Generation Note

Generated/updated: 2026-05-10 14:35 EDT

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
/snap/root-framework/954/usr/local/bin/root -l -b -q -e '.L reproduce_isr_plot.C' -e 'nEvents=20000000; regenerateSamples=true; gOutputDir="/data2/yjlee/ISRsample"; reproduce_isr_plot();'
```

The requested command remains:

```bash
root -l -b -q reproduce_isr_plot.C
```

On this machine, the direct ROOT binary was used because the Snap wrapper has path
visibility issues in `/data/yjlee/ISR`. The macro default is now `nEvents = 20000000`
and `gOutputDir = "/data2/yjlee/ISRsample"`, so the requested command will write
future samples and plot outputs under `/data2/yjlee/ISRsample` unless edited.

## Software Versions

- ROOT: 6.36.04, built for `linuxx8664gcc` on 2025-08-27.
- Macro used for the completed 20M production: `reproduce_isr_plot.C`
- Production macro SHA256:

```text
2e09ecb2d91d36c61a0e4da6c5cf7a89ccf2225e11581ce4b56581963d8c1e15  reproduce_isr_plot.C
```

- Current checked-in macro SHA256 after the ALEPH thrust-alignment patch:

```text
c981516005c3dedb3ce9e38d4f82733f405129f513fac62064452c3169a76412  reproduce_isr_plot.C
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
toy generator in `reproduce_isr_plot.C`, with deterministic `TRandom3` seeds and
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

The stale stored-branch thrust plots were removed. The current thrust-dependent
outputs use all-visible thrust recomputed from particle vectors, with charged
plus neutral final-state particles kept and neutrinos plus flagged ISR photons
excluded. The corrected diagnostic uses 50,000 entries per ISR ON/OFF file.

```text
isr_correction_all_visible_recomputed.png
isr_correction_all_visible_recomputed.pdf
thrust_stored_vs_all_visible_recomputed_diagnostic.png
thrust_stored_vs_all_visible_recomputed_diagnostic.pdf
thrust_all_visible_recomputed_metrics.csv
thrust_all_visible_recomputed_histograms.csv
update_thrust_all_visible_results.C
```

Stored-branch versus recomputed thrust summary:

| Configuration | Category | Mean recomputed T | Mean stored T | Mean stored-recomputed | Mean abs delta | Max abs delta |
|---|---|---:|---:|---:|---:|---:|
| Herwig 7.1.5 equivalent | ISR ON | 0.8805 | 0.8783 | -0.00224 | 0.00224 | 0.0498 |
| PYTHIA 8.317 equivalent | ISR ON | 0.8829 | 0.8807 | -0.00222 | 0.00222 | 0.0579 |
| Sherpa 2.2.6 equivalent | ISR ON | 0.8778 | 0.8755 | -0.00229 | 0.00229 | 0.0437 |
| PYTHIA 8.317 Vincia equivalent | ISR ON | 0.8827 | 0.8805 | -0.00225 | 0.00225 | 0.0520 |

ISR photon energy and angular spectra were produced from the ISR ON particle
arrays by selecting entries with `isISRPhoton=1`. The spectra are normalized per
generated event, so the integral keeps the ISR photon rate difference across
generator-equivalent samples. The angular figure includes both the full polar
angle `theta` and the folded beam angle `min(theta, pi - theta)`.

```text
isr_photon_energy_theta_spectra.png
isr_photon_energy_theta_spectra.pdf
isr_photon_energy_theta_spectra.csv
plot_isr_photon_spectra.C
```

Summary from the 20M ISR ON files:

| Configuration | Emit fraction | Photons/event | Mean ISR photon energy | Mean folded angle |
|---|---:|---:|---:|---:|
| Herwig 7.1.5 equivalent | 0.432 | 0.473 | 3.92 GeV | 0.0246 rad |
| PYTHIA 8.317 equivalent | 0.470 | 0.518 | 4.18 GeV | 0.0247 rad |
| Sherpa 2.2.6 equivalent | 0.508 | 0.564 | 4.44 GeV | 0.0247 rad |
| PYTHIA 8.317 Vincia equivalent | 0.461 | 0.507 | 4.12 GeV | 0.0247 rad |

Visible-energy ISR ON/OFF comparison plots were produced for each
generator-equivalent sample from the stored event-level `visibleEnergy` branch.
The branch corresponds to the default selection: neutrinos are excluded and ISR
photons are not counted as visible energy. Each plotted spectrum is normalized
to unit area and includes a bottom `ISR ON / ISR OFF` shape ratio.

```text
visible_energy_Herwig715_ISR_ON_OFF.png
visible_energy_Herwig715_ISR_ON_OFF.pdf
visible_energy_Pythia8317_ISR_ON_OFF.png
visible_energy_Pythia8317_ISR_ON_OFF.pdf
visible_energy_Sherpa226_ISR_ON_OFF.png
visible_energy_Sherpa226_ISR_ON_OFF.pdf
visible_energy_Pythia8317_Vincia_ISR_ON_OFF.png
visible_energy_Pythia8317_Vincia_ISR_ON_OFF.pdf
visible_energy_isr_on_off_summary.csv
visible_energy_isr_on_off_bins.csv
plot_visible_energy_isr_on_off.C
```

Weighted visible-energy summary:

| Configuration | Mean ON | RMS ON | Mean OFF | RMS OFF |
|---|---:|---:|---:|---:|
| Herwig 7.1.5 equivalent | 88.48 GeV | 4.27 GeV | 90.32 GeV | 2.24 GeV |
| PYTHIA 8.317 equivalent | 88.18 GeV | 4.57 GeV | 90.32 GeV | 2.22 GeV |
| Sherpa 2.2.6 equivalent | 87.87 GeV | 4.86 GeV | 90.35 GeV | 2.12 GeV |
| PYTHIA 8.317 Vincia equivalent | 88.27 GeV | 4.48 GeV | 90.33 GeV | 2.19 GeV |

The recomputed all-visible ISR ON spectra were compared with the published
ALEPH archived-data thrust distribution at sqrt(s) = 91.20 GeV:

- Source: HEPData Table 54, DOI `10.17182/hepdata.12794.v1/t54`
- Local CSV: `/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/ALEPH/HEPData-ins636645-v1-Table_54.csv`
- MC input: four ISR ON ROOT ntuples in `/data2/yjlee/ISRsample`, 50k entries
  per file recomputed from particle vectors
- Normalization: each MC histogram normalized to unit area over `0.58 < T < 1.00`
  using the published ALEPH binning and stored event weights
- ALEPH uncertainty: statistical, `sys_1`, and `sys_2` components summed in quadrature

Output comparison artifacts:

```text
aleph_isr_on_thrust_comparison_all_visible_recomputed.png
aleph_isr_on_thrust_comparison_all_visible_recomputed.pdf
aleph_isr_on_thrust_comparison_all_visible_recomputed_metrics.csv
aleph_isr_on_thrust_comparison_all_visible_recomputed_bins.csv
```

The comparison shows that the stored branch was stale but not the dominant
cause of the bad absolute thrust spectrum. The MC/ALEPH discrepancy remains
large after all-visible recomputation:

| Configuration | chi2/ndf | Mean abs frac. | Central mean abs frac. |
|---|---:|---:|---:|
| Herwig 7.1.5 equivalent | 3567.1 | 7.32 | 2.41 |
| PYTHIA 8.317 equivalent | 3089.0 | 6.53 | 2.30 |
| Sherpa 2.2.6 equivalent | 3369.2 | 6.90 | 2.47 |
| PYTHIA 8.317 Vincia equivalent | 2994.5 | 6.25 | 2.27 |

The corrected spectra still overpopulate low `T`, are closest near `T ~= 0.90`,
and strongly underpopulate `0.93 < T < 0.98`. This points to the fallback toy
generator's final-state topology and thrust-shape model as the main issue, not
only the thrust algorithm. A physics-level agreement with ALEPH requires a
retuned toy model or fully configured external Pythia/Herwig/Sherpa production.

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

## ALEPH Thrust/Selection Cross-Check

The ISR selection and thrust definitions were compared against
`/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/DataProcessing/include`.
For generator-level visible event shapes, the selections agree in intent:
both use visible final-state particles, include charged and neutral visible
particles, and exclude ISR photons from thrust. The ISR macro does this through
PDG/status flags and `isISRPhoton`; the ALEPH code uses `pwflag` categories and
drops `pwflag==4` photons in the PYTHIA8 path. The ISR toy study is particle
level and does not attempt to reproduce the ALEPH detector-level track, neutral,
event-selection, or optional MET variations.

The original ISR macro used a bounded candidate-axis thrust approximation. A
100k-event comparison against ALEPH `getThrust(..., THRUST::OPTIMAL)` gave mean
absolute differences of about `2.1e-3` and maxima of `4.6e-2` to `5.8e-2`,
depending on ISR category. The macro was therefore updated to use an
ALEPH-compatible Herwig-style hemisphere enumeration. After the update, the same
100k-event PYTHIA 8.317-equivalent check gave:

| Sample | Entries | Mean ISR-ALEPH | Mean absolute difference | Max absolute difference |
|---|---:|---:|---:|---:|
| ISR ON | 100,000 | -1.6e-9 | 6.1e-8 | 3.5e-7 |
| ISR OFF | 100,000 | -1.3e-9 | 6.1e-8 | 3.9e-7 |

No checked event had `|delta T| > 1e-6` after the update. The completed 20M ROOT
files retain the pre-update stored `thrust` branch, but their particle-level
arrays are sufficient to recompute ALEPH-aligned thrust. Future regenerated
samples will store the aligned thrust definition directly. Existing thrust
results in this note now use recomputed all-visible thrust and the old
stored-branch thrust plots have been deleted to avoid confusion.

## Analysis Chain

The corrected ISR correction plot is made by reading the ROOT ntuples,
recomputing all-visible thrust from the stored final-state vectors, filling
matched ISR OFF and ISR ON histograms, dividing them bin-by-bin, propagating
statistical uncertainties, and converting the ratios to `TGraphErrors`. The plot
does not use hard-coded final graph points.
