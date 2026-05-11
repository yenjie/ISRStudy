# Real-generator ISR production status

Date: 2026-05-11  
Latest production directory: `/data2/yjlee/ISRsample/real_3M_20260511`  
Previous 20k directory: `/data2/yjlee/ISRsample/real_20260510`

## Goal

Replace the previous toy fallback samples with reusable ROOT ntuples produced from real standalone event generators for `e+e- -> hadrons` at `sqrt(s) = 91.2 GeV`.  The ntuples are used to recompute visible energy, thrust, event-shape observables, ISR photon summaries, and ISR correction factors.

## Important cleanup

All old toy `mc_*.root` samples and toy-derived plots were removed from `/data2/yjlee/ISRsample`.  The old `5M_20260509` toy archive was also removed.  The current real-generator production samples live under:

```text
/data2/yjlee/ISRsample/real_3M_20260511/
```

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

The full plotting macro wrote all figures and then ROOT crashed during
end-of-process cleanup before completing the CSV.  The CSV was regenerated with
the standalone ROOT macro `write_real_isr_stats.C`.

Summary file:

```text
/data2/yjlee/ISRsample/real_3M_20260511/results/real_generator_sample_statistics.csv
```

Key means from the 3M production after recomputing visible energy from
final-state non-neutrino particles within `|eta| < 1.74`, including photons:

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

## Superseded 100k sample statistics

The existing 100k files in `/data2/yjlee/ISRsample/real_100k_20260511` were
generated before the requirements update in this note.  In particular, Herwig
used a lepton-PDF restoration attempt and Sherpa nominal ISR used YFS.  Those
outputs remain useful only as a historical validation snapshot; rerun production
before using the plots or statistics as the current generator comparison.

Summary file:

```text
/data2/yjlee/ISRsample/real_100k_20260511/results/real_generator_sample_statistics.csv
```

Key means from that superseded 100k/event-mode validation set:

| sample | ISR | entries | mean thrust | mean visible energy [GeV] | mean ISR photons | mean ISR photon energy [GeV] |
|---|---:|---:|---:|---:|---:|---:|
| Herwig 7.3.0 | OFF | 100000 | 0.93958 | 89.67 | 0.000 | 0.000 |
| Herwig 7.3.0 | ON | 100000 | 0.94030 | 89.18 | 0.543 | 0.509 |
| PYTHIA 8.315 | OFF | 100000 | 0.93162 | 89.49 | 0.000 | 0.000 |
| PYTHIA 8.315 | ON | 100000 | 0.93173 | 89.36 | 2.733 | 0.206 |
| Sherpa 3.0.3 | OFF | 100000 | 0.93394 | 89.63 | 0.000 | 0.000 |
| Sherpa 3.0.3 | ON | 100000 | 0.93606 | 66.82 | 17.929 | 22.822 |
| PYTHIA 8.315 (Vincia) | OFF | 100000 | 0.93267 | 89.52 | 0.000 | 0.000 |
| PYTHIA 8.315 (Vincia) | ON | 100000 | 0.93245 | 89.33 | 2.003 | 0.206 |

Sherpa YFS ISR was visibly much stronger than the PYTHIA and Herwig configurations in this superseded validation sample.  The updated nominal Sherpa ISR sample uses `PDFESherpa` instead, with YFS kept as a separately labeled alternative.

## Superseded refreshed plots

The existing plots were made by reading the superseded 100k ROOT ntuples:

```text
/data2/yjlee/ISRsample/real_100k_20260511/results/isr_correction_studies_reproduction.png
/data2/yjlee/ISRsample/real_100k_20260511/results/real_isr_correction_double_ratio.png
/data2/yjlee/ISRsample/real_100k_20260511/results/real_isr_on_thrust_vs_aleph.png
/data2/yjlee/ISRsample/real_100k_20260511/results/real_isr_photon_energy_theta_spectra.png
/data2/yjlee/ISRsample/real_100k_20260511/results/real_visible_energy_Herwig730_ISR_ON_OFF.png
/data2/yjlee/ISRsample/real_100k_20260511/results/real_visible_energy_Pythia8315_ISR_ON_OFF.png
/data2/yjlee/ISRsample/real_100k_20260511/results/real_visible_energy_Sherpa303_ISR_ON_OFF.png
/data2/yjlee/ISRsample/real_100k_20260511/results/real_visible_energy_Pythia8315_Vincia_ISR_ON_OFF.png
```

## Repository and Overleaf sync status

Status from the 100k validation update:

- The 100k ROOT ntuples and refreshed plots were produced with ROOT 6.34.04 from `/raid5/root/root-v6.34.04/root/bin/thisroot.sh`.
- The Overleaf analysis note was updated and pushed successfully to `https://git.overleaf.com/69ff6df5f5dd70ea72cedce4`; the generator-requirements markdown update starts at commit `b8feff0`.
- The local `ISRStudy` generator-setting implementation commit is `c9a2868` with message `Apply generator-specific ISR settings`; later status-only commits may also be present.
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
- `plot_real_isr_results.C`: plotting and sample-statistics macro.
