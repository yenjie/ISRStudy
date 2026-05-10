# Real-generator ISR production status

Date: 2026-05-10  
Output directory: `/data2/yjlee/ISRsample/real_20260510`

## Goal

Replace the previous toy fallback samples with reusable ROOT ntuples produced from real standalone event generators for `e+e- -> hadrons` at `sqrt(s) = 91.2 GeV`.  The ntuples are used to recompute visible energy, thrust, event-shape observables, ISR photon summaries, and ISR correction factors.

## Important cleanup

All old toy `mc_*.root` samples and toy-derived plots were removed from `/data2/yjlee/ISRsample`.  The old `5M_20260509` toy archive was also removed.  The current real-generator validation samples live under:

```text
/data2/yjlee/ISRsample/real_20260510/
```

## Real generators and references

- PYTHIA 8.315, using the local library `/usr/lib/libpythia8.so` and XML docs at `/usr/share/Pythia8/xmldoc`.  Reference: PYTHIA 8.3 manual/release, SciPost Phys. Codebases 8 and 8-r8.3, <https://scipost.org/SciPostPhysCodeb.8> and <https://scipost.org/SciPostPhysCodeb.8-r8.3>.
- PYTHIA 8.315 with Vincia, using `Pythia8::Vincia` through `Pythia::setShowerModelPtr` plus `/usr/share/Pythia8/tunes/VinciaDefault.cmnd`.  Reference: PYTHIA Vincia manual, <https://pythia.org/latest-manual/Vincia.html>.
- Herwig 7.3.0, using the local standalone binary at `/raid5/data/yjlee/strangeness/Strangeness/MainAnalysis/20260218_KtoPiInversion/external/herwig-7.3.0/bin/Herwig`.  Reference: Herwig 7.3 release note, EPJC 84, 1053 (2024), <https://doi.org/10.1140/epjc/s10052-024-13211-9>.
- Sherpa 3.0.3, using the local standalone binary at `/data/yjlee/OnePointChargeCorrelator/external/sherpa-3.0.3-build/outputs/bin/Sherpa`.  Reference: Sherpa 3 release paper, JHEP 2024, 156, <https://doi.org/10.1007/JHEP12(2024)156>; Sherpa release/download documentation, <https://sherpa-team.gitlab.io/changelog.html>.

## ISR on/off definitions

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
- Sherpa ISR OFF uses `PDF_LIBRARY: None` and `11 -11 -> 93 93`.
- Sherpa ISR ON uses the installed LEP YFS pattern:
  - `YFS: MODE: ISR`
  - electron mass enabled in `PARTICLE_DATA`

## Production command

```bash
EVENTS=20000 OUTDIR=/data2/yjlee/ISRsample/real_20260510 TARGETS=all FORCE=1 ./scripts/run_real_isr_production.sh
```

The event count is intentionally configurable.  The 20k/event-mode sample is a validation refresh; larger 5M or 20M production should use the same script after accepting the generator-specific ISR settings.

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
/data2/yjlee/ISRsample/real_20260510/results/real_generator_sample_statistics.csv
```

Key means from the 20k/event-mode validation set:

| sample | ISR | entries | mean thrust | mean visible energy [GeV] | mean ISR photons | mean ISR photon energy [GeV] |
|---|---:|---:|---:|---:|---:|---:|
| Herwig 7.3.0 | OFF | 20000 | 0.93946 | 89.65 | 0.000 | 0.000 |
| Herwig 7.3.0 | ON | 20000 | 0.93983 | 89.17 | 0.543 | 0.513 |
| PYTHIA 8.315 | OFF | 20000 | 0.93159 | 89.48 | 0.000 | 0.000 |
| PYTHIA 8.315 | ON | 20000 | 0.93205 | 89.32 | 2.731 | 0.212 |
| Sherpa 3.0.3 | OFF | 20000 | 0.93408 | 89.59 | 0.000 | 0.000 |
| Sherpa 3.0.3 | ON | 20000 | 0.93565 | 66.69 | 17.972 | 22.98 |
| PYTHIA 8.315 (Vincia) | OFF | 20000 | 0.93297 | 89.51 | 0.000 | 0.000 |
| PYTHIA 8.315 (Vincia) | ON | 20000 | 0.93258 | 89.33 | 2.004 | 0.198 |

Sherpa YFS ISR is visibly much stronger than the PYTHIA and Herwig configurations in this validation sample.  Treat that as a configuration item to review before committing CPU/storage to a 20M-event run.

## Refreshed plots

All plots are made by reading the new ROOT ntuples:

```text
/data2/yjlee/ISRsample/real_20260510/results/real_isr_correction_double_ratio.png
/data2/yjlee/ISRsample/real_20260510/results/real_isr_on_thrust_vs_aleph.png
/data2/yjlee/ISRsample/real_20260510/results/real_isr_photon_energy_theta_spectra.png
/data2/yjlee/ISRsample/real_20260510/results/real_visible_energy_Herwig730_ISR_ON_OFF.png
/data2/yjlee/ISRsample/real_20260510/results/real_visible_energy_Pythia8315_ISR_ON_OFF.png
/data2/yjlee/ISRsample/real_20260510/results/real_visible_energy_Sherpa303_ISR_ON_OFF.png
/data2/yjlee/ISRsample/real_20260510/results/real_visible_energy_Pythia8315_Vincia_ISR_ON_OFF.png
```

## Code entry points

- `src/real_isr_ntuple_producer.cc`: direct PYTHIA/Vincia production and HepMC-to-ROOT conversion for Herwig/Sherpa.
- `scripts/build_real_isr_tools.sh`: build command for the producer.
- `scripts/run_real_isr_production.sh`: standalone real-generator production wrapper.
- `cards/`: real-generator Herwig and Sherpa cards for ISR ON/OFF.
- `plot_real_isr_results.C`: plotting and sample-statistics macro.
