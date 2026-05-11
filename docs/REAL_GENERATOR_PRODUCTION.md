# Real-generator ISR production status

Date: 2026-05-11  
Latest production directory: `/data2/yjlee/ISRsample/real_3M_20260511`  

## Goal

Replace the previous toy fallback samples with reusable ROOT ntuples produced from real standalone event generators for `e+e- -> hadrons` at `sqrt(s) = 91.2 GeV`.  The ntuples are used to recompute visible energy, thrust, event-shape observables, ISR photon summaries, and ISR correction factors.

## Important cleanup

All old toy `mc_*.root` samples and toy-derived plots were removed from `/data2/yjlee/ISRsample`.  The old `5M_20260509` toy archive was also removed.  The current real-generator production samples live under:

```text
/data2/yjlee/ISRsample/real_3M_20260511/
```

The checked-in `ISRStudy` repository no longer tracks toy-derived plots, old
100k/20k result snapshots, or legacy macros that used the Herwig715,
Pythia8317, and Sherpa226 toy file names.  The plotting macro defaults to the
3M directory above, and the low-level production wrapper now requires an
explicit `OUTDIR` instead of silently writing to an old default directory.

## Real generators and references

- PYTHIA 8.315, using the local library `/usr/lib/libpythia8.so` and XML docs at `/usr/share/Pythia8/xmldoc`.  Reference: PYTHIA 8.3 manual/release, SciPost Phys. Codebases 8 and 8-r8.3, <https://scipost.org/SciPostPhysCodeb.8> and <https://scipost.org/SciPostPhysCodeb.8-r8.3>.
- PYTHIA 8.315 with Vincia, using `Pythia8::Vincia` through `Pythia::setShowerModelPtr` plus `/usr/share/Pythia8/tunes/VinciaDefault.cmnd`.  Reference: PYTHIA Vincia manual, <https://pythia.org/latest-manual/Vincia.html>.
- Herwig 7.3.0, using the local standalone binary at `/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/herwig-7.3.0/bin/Herwig`.  Reference: Herwig 7.3 release note, EPJC 84, 1053 (2024), <https://doi.org/10.1140/epjc/s10052-024-13211-9>.
- Sherpa 3.0.3, using the local standalone binary at `/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3-build/outputs/bin/Sherpa`.  Reference: Sherpa 3 release paper, JHEP 2024, 156, <https://doi.org/10.1007/JHEP12(2024)156>; Sherpa release/download documentation, <https://sherpa-team.gitlab.io/changelog.html>.

## ISR on/off definitions

The generator/process, beam energy, hadronization, decays, output schema, and
event selection must be held fixed within each generator pair.  The required
ON/OFF changes are the generator-specific controls below.

- PYTHIA ISR OFF:
  - `PDF:lepton = off`
  - `PartonLevel:ISR = off`
  - `SpaceShower:QEDshowerByL = off`
- PYTHIA ISR ON:
  - `PDF:lepton = on`
  - `PartonLevel:ISR = on`
  - `SpaceShower:QEDshowerByL = on`
- PYTHIA hard process and decay restriction:
  - `WeakSingleBoson:ffbar2gmZ = on`
  - `WeakZ0:gmZmode = 0`
  - `23:onMode = off`
  - `23:onIfAny = 1 2 3 4 5`
- Vincia uses the same PYTHIA hard-process, decay, beam, and ISR settings, with the Vincia shower model attached.
- Herwig 7.3.0 baseline:
  - read `snippets/EECollider.in`
  - insert `MEee2gZ2qq`
  - set the luminosity-function energy to `91.1876*GeV`
  - set `/Herwig/Shower/ShowerHandler:Interactions QCD`
- Herwig 7.3.0 QED-shower comparison:
  - same hard setup as the baseline
  - requested physics mode: QCD plus QED showering
  - the installed Herwig 7.3.0 switch spelling is
    `/Herwig/Shower/ShowerHandler:Interactions QEDQCD`; the literal
    `QCDandQED` token is rejected by the installed binary
  - label the output `QEDshower`, not structure-function ISR
  - do not claim this is a direct equivalent of PYTHIA `PDF:lepton` or Sherpa
    `PDFESherpa` without verifying that from installed Herwig documentation or source.
- Sherpa 3.0.3 ISR OFF:
  - `BEAMS: [11, -11]`
  - `BEAM_ENERGIES: 45.5938`
  - `PDF_LIBRARY: None`
  - `YFS: MODE: None`
- Sherpa 3.0.3 nominal ISR ON, structure-function mode:
  - `PDF_LIBRARY: PDFESherpa`
  - `PDF_SET: Default`
  - `ISR_E_ORDER: 1`
  - `ISR_E_SCHEME: 2`
  - `YFS: MODE: None`
- Sherpa 3.0.3 YFS alternative:
  - `PDF_LIBRARY: None`
  - `YFS: MODE: ISR`
  - `YFS: BETA: 1`
- Never enable `PDFESherpa` and YFS ISR in the same nominal ISR sample.

The local ALEPH agentic ISR correction uses precomputed PYTHIA 8 files,
`Isr/isr0_ALL.root` and `Isr/isr1_ALL.root`, with `tgenBefore/thrust`.  It is
not a Sherpa ISR correction.  Those files contain 2.5M entries each and show
`sqrt_sHat` fixed at `91.1876 GeV` for ISR OFF, while ISR ON has reduced
`sqrt_sHat` from photon radiation.

## Settings smoke test

The updated settings were checked with tiny local samples after the requirements
change:

- PYTHIA ISR OFF and ISR ON each initialized and wrote 10 events with the
  requested `PDF:lepton`, `PartonLevel:ISR`, and `SpaceShower:QEDshowerByL`
  settings.
- Herwig OFF, Herwig `QEDshower`, Sherpa OFF, Sherpa `PDFESherpa` ISR, and
  Sherpa YFS alternative each wrote two-event ROOT files under
  `/tmp/isr_generator_requirement_smoke`.
- The Herwig smoke test verified that the installed `Interactions` switch
  accepts `QEDQCD`, not `QCDandQED`.

## Production command

```bash
cd /raid5/root/root-v6.34.04/root/bin
. ./thisroot.sh
cd /raid5/data/yjlee/ISR
EVENTS=3000000 MAX_WORKERS=15 OUTDIR=/data2/yjlee/ISRsample/real_3M_20260511 FORCE=0 ./scripts/run_real_isr_production_10worker.sh
```

The event count is intentionally configurable.  The wrapper name still contains
`10worker`, but the current invocation honors `MAX_WORKERS=15`.  With the
updated requirements, `all` includes the Sherpa YFS alternative in addition to
the nominal Sherpa `PDFESherpa` ISR sample.

## Produced ROOT ntuples

Each file contains a TTree named `Events` with event-level branches and particle-level vectors.

```text
mc_Herwig730_ISR_OFF.root
mc_Herwig730_QEDshower.root
mc_Pythia8315_ISR_OFF.root
mc_Pythia8315_ISR_ON.root
mc_Sherpa303_ISR_OFF.root
mc_Sherpa303_ISR_ON.root
mc_Sherpa303_ISR_YFS.root
mc_Pythia8315_Vincia_ISR_OFF.root
mc_Pythia8315_Vincia_ISR_ON.root
```

The 2026-05-11 3M production initially exposed a Herwig-only issue: the Herwig
HepMC card had `set /Herwig/Analysis/HepMC:PrintEvent 1000000`, so the first
large Herwig files stopped at 1M converted `Events` entries.  The card now uses
`@N_EVENTS@`, and `scripts/run_real_isr_production.sh` substitutes the requested
`EVENTS` value before running Herwig.

## 3M sample validation

All current ROOT ntuples in `/data2/yjlee/ISRsample/real_3M_20260511` contain
exactly 3,000,000 entries in the `Events` TTree:

```text
mc_Herwig730_ISR_OFF.root             3000000
mc_Herwig730_QEDshower.root          3000000
mc_Pythia8315_ISR_OFF.root           3000000
mc_Pythia8315_ISR_ON.root            3000000
mc_Pythia8315_Vincia_ISR_OFF.root    3000000
mc_Pythia8315_Vincia_ISR_ON.root     3000000
mc_Sherpa303_ISR_OFF.root            3000000
mc_Sherpa303_ISR_ON.root             3000000
mc_Sherpa303_ISR_YFS.root            3000000
```

The refreshed plots are in:

```text
/data2/yjlee/ISRsample/real_3M_20260511/results
```

The thrust-to-ALEPH comparison is written both as spectrum-only figures and as
separate ratio-panel figures:

```text
real_isr_on_thrust_vs_aleph.png
real_isr_on_thrust_vs_aleph_ratio.png
real_isr_off_thrust_vs_aleph.png
real_isr_off_thrust_vs_aleph_ratio.png
```

The full plotting macro wrote all figures and then ROOT crashed during
end-of-process cleanup before completing the CSV.  The CSV was regenerated with
the standalone ROOT macro `write_real_isr_stats.C`.

Summary file:

```text
/data2/yjlee/ISRsample/real_3M_20260511/results/real_generator_sample_statistics.csv
```

Key means from the 3M production after recomputing visible energy from
final-state non-neutrino particles within `|eta| < 1.74`, including photons.
The `mean thrust` column here is the original stored `Events.thrust` branch,
which was the pre-diagnostic visible-particle thrust with tagged ISR photons
removed.  The current nominal all-final thrust is stored in the derived
diagnostic tree described below.

| sample | mode | entries | mean thrust | mean visible energy [GeV] | mean ISR photons | mean ISR photon energy [GeV] |
|---|---:|---:|---:|---:|---:|---:|
| Herwig 7.3.0 QED shower | OFF | 3000000 | 0.94003 | 82.04 | 0.000 | 0.000 |
| Herwig 7.3.0 QED shower | QEDshower | 3000000 | 0.94040 | 82.05 | 0.545 | 0.518 |
| PYTHIA 8.315 | OFF | 3000000 | 0.93197 | 82.11 | 0.000 | 0.000 |
| PYTHIA 8.315 | ISR ON | 3000000 | 0.93163 | 81.99 | 2.730 | 0.202 |
| Sherpa 3.0.3 PDFESherpa | OFF | 3000000 | 0.93395 | 82.22 | 0.000 | 0.000 |
| Sherpa 3.0.3 PDFESherpa | ISR ON | 3000000 | 0.93598 | 82.04 | 17.826 | 22.785 |
| Sherpa 3.0.3 YFS | OFF | 3000000 | 0.93395 | 82.22 | 0.000 | 0.000 |
| Sherpa 3.0.3 YFS | ISR ON | 3000000 | 0.93593 | 82.07 | 17.932 | 22.836 |
| PYTHIA 8.315 (Vincia) | OFF | 3000000 | 0.93277 | 82.16 | 0.000 | 0.000 |
| PYTHIA 8.315 (Vincia) | ISR ON | 3000000 | 0.93265 | 81.96 | 2.003 | 0.202 |

## Branch-level audit notes

- The current figures and the CSV above recompute visible energy from the
  particle vectors using final-state non-neutrinos within `|eta| < 1.74`,
  including photons.  The existing 3M ROOT files were produced before that
  acceptance-definition update, so their stored event-level `visibleEnergy`
  branch still contains the production-time definition.  Use the recomputed
  values in the current analysis code, or regenerate the samples with the
  current producer before relying on the scalar `visibleEnergy` branch.
- The existing Sherpa 3M files carry large HepMC generator bookkeeping weights
  in the `weight` branch.  The current ISR correction, ALEPH comparison,
  visible-energy, and photon plots use unweighted event counts.  The producer
  has been patched so future HepMC conversions write `weight = 1.0` for these
  event-shape samples.
- A stored-vs-recomputed thrust spot check over the first 100 events in each
  3M file exactly reproduced the `thrust` branch from the stored final-state
  non-neutrino particles with tagged ISR photons excluded.  The thrust branch
  itself is therefore consistent with the current thrust algorithm and
  selection used at production time.  The producer default has since been
  updated so future `Events.thrust` branches use the requested all-stable-final
  convention, including neutrinos and tagged ISR photons.

## Endpoint mass and boost diagnostic

An additional derived-ntuple pass was added to test whether Sherpa endpoint
behavior is driven by low visible mass, a longitudinal boost, or unusually
large explicit ISR photon activity.  The pass does not modify the original
generator ROOT files.  It reads the stored particle vectors and writes one
compact ROOT file per ISR ON/OFF sample.  The current all-final 3M diagnostic
is stored under:

```text
/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics_allfinal_3M
```

Each derived `EndpointDiagnostics` tree stores:

```text
T_lab_excluding_ISR_photons
T_lab_including_ISR_photons
T_visibleCM_excluding_ISR_photons
T_lab_allFinal_including_ISR_photons
T_lab_allFinal_excluding_ISR_photons
T_allFinalCM_including_ISR_photons
cosTheta_thrust_lab_allFinal_including_ISR_photons
absCosTheta_thrust_lab_allFinal_including_ISR_photons
Evis_excluding_ISR_photons
Evis_including_ISR_photons
Mvis_excluding_ISR_photons
Mvis_including_ISR_photons
sPrime_vis_over_s
betaZ_vis
Eall_including_ISR_photons
Eall_excluding_ISR_photons
Mall_including_ISR_photons
Mall_excluding_ISR_photons
sAll_including_ISR_over_s
betaZ_all_including_ISR
leading_ISR_photon_energy
leading_ISR_photon_theta
N_ISR_photons
sum_ISR_photon_energy
event_weight
```

For the current nominal 3M diagnostic, `T_lab_allFinal_including_ISR_photons`
uses all stable final-state particles, including neutrinos and tagged ISR
photons.  The visible-particle variants keep the earlier ALEPH-style
cross-checks that exclude neutrinos and either remove or include tagged ISR
photons.  The CM variants boost the selected four-vector sum to its own rest
frame before recomputing thrust.

The 2D correlation plots are written to:

```text
/data2/yjlee/ISRsample/real_3M_20260511/results_allfinal_3M/endpoint_diagnostics
```

The same diagnostic producer was also run on the ISR OFF samples to compare
the nominal thrust definition and the visible-particle cross-check definitions
in the ISR correction ratio:

```text
A: T_lab_excluding_ISR_photons
B: T_lab_including_ISR_photons
C: T_visibleCM_excluding_ISR_photons
N: T_lab_allFinal_including_ISR_photons
```

Those `C_ISR(T)` plots and bin-summary CSV files are written to:

```text
/data2/yjlee/ISRsample/real_3M_20260511/results_allfinal_3M/thrust_definition_corrections
```

The 3M all-final endpoint summary is definition-dependent.  In
the last ALEPH thrust bin, `0.99 < T < 1.00`, Sherpa PDFESherpa has
`C_ISR = 1.150 +/- 0.007` for definition N and
`C_ISR = 1.149 +/- 0.007` for definition B, where explicit tagged ISR photons
are kept in the thrust input.  The same sample gives
`C_ISR = 0.492 +/- 0.002` for definition A and
`C_ISR = 0.494 +/- 0.002` for definition C, where tagged ISR photons are
removed from the visible-particle thrust input.  The endpoint sign therefore
changes with the photon/neutrino convention:

```text
N: T_lab_allFinal_including_ISR_photons
A: T_lab_excluding_ISR_photons
B: T_lab_including_ISR_photons
C: T_visibleCM_excluding_ISR_photons
```

## Repository and Overleaf sync status

- The 3M ROOT ntuples and refreshed plots were produced with ROOT 6.34.04 from
  `/raid5/root/root-v6.34.04/root/bin/thisroot.sh`.
- The Overleaf analysis note references only the current 3M result files.
- A byte-for-byte check confirmed that the Overleaf result files match
  `/data2/yjlee/ISRsample/real_3M_20260511/results`.
- Pushing `ISRStudy` to GitHub is blocked because the configured remote is `git@github.com:yenjie/ISRStudy.git` and GitHub returns `Repository not found`.
- Once the correct GitHub remote or access permission is available, push the local commit with:

```bash
git -C /raid5/data/yjlee/ISR/ISRStudy push --set-upstream origin main
```

## Code entry points

- `src/real_isr_ntuple_producer.cc`: direct PYTHIA/Vincia production and HepMC-to-ROOT conversion for Herwig/Sherpa.
- `scripts/build_real_isr_tools.sh`: build command for the producer.
- `scripts/run_real_isr_production.sh`: standalone real-generator production wrapper.
- `cards/`: real-generator Herwig and Sherpa cards for ISR ON/OFF.
- `macros/plot_real_isr_results.C`: plotting and sample-statistics macro.
