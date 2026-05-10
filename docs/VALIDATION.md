# Validation

## Smoke Test

Run from the repository root:

```bash
./scripts/test_isr_macro.sh
```

The smoke test:

1. Generates eight small ROOT files under `/tmp`.
2. Builds the ISR correction plot from those files.
3. Verifies that every `Events` tree has the requested number of entries.
4. Verifies that both PNG and PDF plot outputs exist.

Use a larger temporary sample if needed:

```bash
EVENTS=10000 ./scripts/test_isr_macro.sh
```

Keep the temporary output for inspection:

```bash
KEEP_ISRTEST=1 ./scripts/test_isr_macro.sh
```

## Full 20M Production Checks

The completed production was verified with ROOT:

```text
/data2/yjlee/ISRsample/mc_Herwig715_ISR_ON.root entries=20000000
/data2/yjlee/ISRsample/mc_Herwig715_ISR_OFF.root entries=20000000
/data2/yjlee/ISRsample/mc_Pythia8317_ISR_ON.root entries=20000000
/data2/yjlee/ISRsample/mc_Pythia8317_ISR_OFF.root entries=20000000
/data2/yjlee/ISRsample/mc_Sherpa226_ISR_ON.root entries=20000000
/data2/yjlee/ISRsample/mc_Sherpa226_ISR_OFF.root entries=20000000
/data2/yjlee/ISRsample/mc_Pythia8317_Vincia_ISR_ON.root entries=20000000
/data2/yjlee/ISRsample/mc_Pythia8317_Vincia_ISR_OFF.root entries=20000000
```

The regenerated reference image was checked as:

```text
PNG: 1596x872, 8-bit sRGB
PDF: one page, ROOT 6.36.04
```

## ISR ON Summary From 20M Files

| Configuration | Emit fraction | Mean ISR photons/event | Mean ISR photon energy/event | Mean ISR photon energy/emitting event |
|---|---:|---:|---:|---:|
| Herwig 7.1.5 equivalent | 43.23% | 0.473 | 1.857 GeV | 4.294 GeV |
| PYTHIA 8.317 equivalent | 47.01% | 0.518 | 2.167 GeV | 4.609 GeV |
| Sherpa 2.2.6 equivalent | 50.77% | 0.564 | 2.505 GeV | 4.934 GeV |
| PYTHIA 8.317 Vincia equivalent | 46.05% | 0.507 | 2.087 GeV | 4.532 GeV |

