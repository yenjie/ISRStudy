# Real-generator ISR production status

Date: 2026-05-10  
Latest validation directory: `/data2/yjlee/ISRsample/real_100k_20260511`  
Previous 20k directory: `/data2/yjlee/ISRsample/real_20260510`

## Goal

Replace the previous toy fallback samples with reusable ROOT ntuples produced from real standalone event generators for `e+e- -> hadrons` at `sqrt(s) = 91.2 GeV`.  The ntuples are used to recompute visible energy, thrust, event-shape observables, ISR photon summaries, and ISR correction factors.

## Important cleanup

All old toy `mc_*.root` samples and toy-derived plots were removed from `/data2/yjlee/ISRsample`.  The old `5M_20260509` toy archive was also removed.  The current real-generator validation samples live under:

```text
/data2/yjlee/ISRsample/real_100k_20260511/
```

## Real generators and references

- PYTHIA 8.315, using the local library `/usr/lib/libpythia8.so` and XML docs at `/usr/share/Pythia8/xmldoc`.  Reference: PYTHIA 8.3 manual/release, SciPost Phys. Codebases 8 and 8-r8.3, <https://scipost.org/SciPostPhysCodeb.8> and <https://scipost.org/SciPostPhysCodeb.8-r8.3>.
- PYTHIA 8.315 with Vincia, using `Pythia8::Vincia` through `Pythia::setShowerModelPtr` plus `/usr/share/Pythia8/tunes/VinciaDefault.cmnd`.  Reference: PYTHIA Vincia manual, <https://pythia.org/latest-manual/Vincia.html>.
- Herwig 7.3.0, using the local standalone binary at `/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/herwig-7.3.0/bin/Herwig`.  Reference: Herwig 7.3 release note, EPJC 84, 1053 (2024), <https://doi.org/10.1140/epjc/s10052-024-13211-9>.
- Sherpa 3.0.3, using the local standalone binary at `/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3-build/outputs/bin/Sherpa`.  Reference: Sherpa 3 release paper, JHEP 2024, 156, <https://doi.org/10.1007/JHEP12(2024)156>; Sherpa release/download documentation, <https://sherpa-team.gitlab.io/changelog.html>.

## ISR on/off definitions

The production was audited so the generator/process, beam energy, hadronization,
decays, output schema, and event selection are held fixed within each generator
pair.  The only intended ON/OFF changes are the generator-specific ISR controls
listed below.

- PYTHIA ISR OFF:
  - `PDF:lepton = off`
  - `SpaceShower:QEDshowerByL = off`
- PYTHIA ISR ON:
  - `PDF:lepton = on`
  - `SpaceShower:QEDshowerByL = on`
- Vincia uses the same PYTHIA beam/ISR settings, with the Vincia shower model attached.
- Herwig ISR OFF uses the standard `snippets/EECollider.in`, which sets the lepton shower PDFs to `NULL`.
- Herwig ISR ON restores the built-in lepton PDFs for the shower:
  - `/Herwig/Shower/ShowerHandler:PDFA = /Herwig/Partons/LeptonPDF`
  - `/Herwig/Shower/ShowerHandler:PDFB = /Herwig/Partons/LeptonPDF`
- Sherpa ISR OFF uses `PDF_LIBRARY: None`, `11 -11 -> 93 93`, and the same
  `PARTICLE_DATA: 11: Massive: true` setting as ISR ON.
- Sherpa ISR ON uses the installed LEP YFS pattern:
  - `YFS: MODE: ISR`

After the audit, the actual Sherpa OFF/ON card diff is only:

```diff
+YFS:
+  MODE: ISR
```

The local ALEPH agentic ISR correction uses precomputed PYTHIA 8 files,
`Isr/isr0_ALL.root` and `Isr/isr1_ALL.root`, with `tgenBefore/thrust`.  It is
not a Sherpa ISR correction.  Those files contain 2.5M entries each and show
`sqrt_sHat` fixed at `91.1876 GeV` for ISR OFF, while ISR ON has reduced
`sqrt_sHat` from photon radiation.

## Production command

```bash
EVENTS=100000 MAX_WORKERS=15 OUTDIR=/data2/yjlee/ISRsample/real_100k_20260511 FORCE=0 ./scripts/run_real_isr_production_10worker.sh
```

The event count is intentionally configurable.  The 100k/event-mode sample is a validation refresh; larger 2M, 5M, or 20M production should use the same ISR definitions after accepting the generator-specific ISR settings.  The wrapper name still contains `10worker`, but the current invocation honors `MAX_WORKERS=15`.

## Produced ROOT ntuples

Each file contains a TTree named `Events` with event-level branches and particle-level vectors.

```text
mc_Herwig730_ISR_OFF.root
mc_Herwig730_ISR_ON.root
mc_Pythia8315_ISR_OFF.root
mc_Pythia8315_ISR_ON.root
mc_Sherpa303_ISR_OFF.root
mc_Sherpa303_ISR_ON.root
mc_Pythia8315_Vincia_ISR_OFF.root
mc_Pythia8315_Vincia_ISR_ON.root
```

## Current sample statistics

Summary file:

```text
/data2/yjlee/ISRsample/real_100k_20260511/results/real_generator_sample_statistics.csv
```

Key means from the 100k/event-mode validation set:

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

Sherpa YFS ISR is visibly much stronger than the PYTHIA and Herwig configurations in this validation sample.  Treat that as a configuration item to review before committing CPU/storage to a 20M-event run.

## Refreshed plots

All plots are made by reading the new ROOT ntuples:

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

## Code entry points

- `src/real_isr_ntuple_producer.cc`: direct PYTHIA/Vincia production and HepMC-to-ROOT conversion for Herwig/Sherpa.
- `scripts/build_real_isr_tools.sh`: build command for the producer.
- `scripts/run_real_isr_production.sh`: standalone real-generator production wrapper.
- `cards/`: real-generator Herwig and Sherpa cards for ISR ON/OFF.
- `plot_real_isr_results.C`: plotting and sample-statistics macro.
