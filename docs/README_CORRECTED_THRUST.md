# Corrected thrust status

Updated: 2026-05-10 14:35 EDT

The existing 20M ROOT ntuples in `/data2/yjlee/ISRsample` remain the particle-level
source samples. Their stored event-level `thrust` branch was produced with the
older bounded candidate-axis approximation and should not be used for current
thrust-spectrum plots.

Current thrust-dependent results use all-visible thrust recomputed from the
stored final-state particle vectors:

- keep charged and neutral visible final-state particles
- exclude neutrinos
- exclude particles flagged with `isISRPhoton=1`
- use the ALEPH-compatible Herwig-style hemisphere enumeration implemented in
  `update_thrust_all_visible_results.C`

The stale stored-branch thrust plot artifacts were removed. Current corrected
outputs are:

```text
isr_correction_all_visible_recomputed.png
isr_correction_all_visible_recomputed.pdf
aleph_isr_on_thrust_comparison_all_visible_recomputed.png
aleph_isr_on_thrust_comparison_all_visible_recomputed.pdf
thrust_stored_vs_all_visible_recomputed_diagnostic.png
thrust_stored_vs_all_visible_recomputed_diagnostic.pdf
thrust_all_visible_recomputed_metrics.csv
aleph_isr_on_thrust_comparison_all_visible_recomputed_metrics.csv
```

The 50k-entry-per-file diagnostic shows the stored branch is biased low by about
0.0021--0.0023 in mean thrust, with rare event-level shifts up to about 0.04--0.058.
This stale branch is a real issue, but it does not explain the large disagreement
with ALEPH. The remaining absolute thrust-spectrum disagreement is caused mainly
by the fallback toy generator topology and requires a retuned toy model or real
configured Pythia/Herwig/Sherpa production.

## Hadron-only / no-photon check

Updated: 2026-05-10 15:55 EDT

The thrust algorithm has been checked against the local ALEPH agentic
`THRUST::OPTIMAL` implementation and agrees at the `1e-7` level. A separate
selection diagnostic was added in `macros/diagnose_hadron_no_photon_thrust.C` to
test whether ordinary final-state photons explain the bad ALEPH thrust-spectrum
agreement.

Selections tested on 50k PYTHIA 8.317-equivalent ISR ON events:

- all visible: neutrinos and flagged ISR photons excluded
- no photons: all visible plus every PDG 22 photon excluded
- hadron-only, no photons: neutrinos, flagged ISR photons, ordinary photons, and
  leptons excluded; require at least 5 selected hadrons and `E_had > 20 GeV`

Key result: removing photons changes the spectrum only modestly. For the
PYTHIA-equivalent sample, the full-range mean absolute MC/ALEPH fractional
difference improves from 6.53 to 6.13 with photons removed, and to 5.72 with the
hadron-only/no-photon selection. The central `0.65<T<0.95` metric remains about
2.2--2.3, so the dominant problem is still the fallback toy generator topology,
not the thrust algorithm or photon inclusion alone.

New outputs:

```text
thrust_selection_mode_diagnostic_pythia.png
aleph_isr_on_thrust_comparison_hadron_no_photon.png
isr_correction_hadron_no_photon.png
thrust_hadron_no_photon_metrics.csv
thrust_hadron_no_photon_bins.csv
```
