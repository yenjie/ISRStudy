# ISR Modelling: KKMC Compared with PYTHIA, Herwig and Sherpa

This note has two parts.

1. What KKMC/KKMCee actually does for initial-state radiation, and how that
   compares with the ISR models used by the PYTHIA 8.315, Herwig 7.3.0 and
   Sherpa 3.0.3 samples in this repository.
2. An audit of the existing 3M study, with the incorrect statements and results
   identified and corrected, and the source/card fixes that were applied.

All numbers labelled "measured" come from the existing
`/data2/yjlee/ISRsample/real_3M_20260511` ntuples (3,000,000 events per file)
unless a smaller sample is stated explicitly.

---

## 1. KKMC / KKMCee

### 1.1 What it is

KKMC (also written `KK MC`, originally `KK2f`) is the precision multi-photon
Monte Carlo for `e+e- -> f fbar + n gamma` by Jadach, Ward and Was, with Yost
and Siodmok joining for the modern versions.  It covers
`f = mu, tau, d, u, s, c, b` from the tau threshold up to about 1 TeV, and it
was the reference ISR/FSR generator for the LEP and SLC fermion-pair and
Z-lineshape measurements.

- Original Fortran program: Comput. Phys. Commun. 130 (2000) 260,
  arXiv:hep-ph/9912214.
- KKMC 4.22, with quark-initiated collisions and updated electroweak input:
  Comput. Phys. Commun. 260 (2021) 107734, arXiv:2007.07964.
- KKMCee, the complete C++ rewrite, version 5.00 reproducing the Fortran
  benchmarks: Comput. Phys. Commun. 283 (2023) 108556, arXiv:2204.11949.
  Current release 5.02 at <https://github.com/KrakowHEPSoft/KKMCee>.

### 1.2 The ISR model: CEEX

The distinguishing feature is **CEEX**, coherent exclusive exponentiation.  YFS
(Yennie-Frautschi-Suura) soft-photon resummation is applied at the level of the
**amplitude** rather than at the level of the cross section:

- All real photons are generated **exclusively**, with exact multi-photon phase
  space.  Photons carry their true transverse momentum; nothing is collapsed
  into a collinear `x` fraction.
- Hard-photon residuals are included with exact `O(alpha)` and `O(alpha^2)`
  next-to-leading-log matrix elements.
- Because the sum over photon multiplicities is done coherently inside
  `|sum_i M_i|^2`, ISR, FSR and their **interference (IFI)** are treated
  consistently to all orders, including the automatic suppression of IFI near a
  narrow resonance such as the Z.
- The older **EEX** mode is the cross-section-level YFS exponentiation.  EEX
  reaches `O(alpha^3 L^3)` leading logarithms, which CEEX does not, so the two
  modes are complementary cross-checks of the same program.

Electroweak form factors come from the DIZET library (6.45 in KKMCee), so KKMC
carries genuine one-loop electroweak corrections in addition to QED, which none
of the three general-purpose generators used here do in this configuration.
Beamstrahlung (CIRCE1) and beam energy spread are available.

Quoted precision is a 0.2% tag at LEP2 energies for the total cross section; at
the Z pole KKMC/ZFITTER-class calculations are what the LEP experiments used to
unfold ISR from the lineshape.

Run control is through a `KKMCee_defaults` file overridden by a short user
`pro.input`.  The switches relevant here are `KeyISR`, `KeyFSR`, `KeyINT`
(IFI; default on) and `KeyGPS` (CEEX level; default on), plus `KeyFix` for the
beamstrahlung/beam-spread mode.  Exact `xpar` indices should be read from
`SRCee/KKMCee_defaults` in the release rather than assumed.

### 1.3 What KKMC is not

KKMC is a fermion-pair generator, not a general-purpose event generator.  For a
hadronic event-shape study the quark pair plus photons still has to be showered
and hadronized externally (JETSET/PYTHIA 6 in the Fortran versions; KKMCee
writes HepMC3 for interfacing to a modern shower).  So a KKMC thrust
distribution inherits the QCD shower and hadronization model of whatever it is
interfaced to, and only the QED/electroweak part is KKMC's own.

---

## 2. Side-by-side ISR comparison

| | PYTHIA 8.315 | PYTHIA 8.315 Vincia | Herwig 7.3.0 | Sherpa 3.0.3 `PDFESherpa` | Sherpa 3.0.3 YFS | KKMC / KKMCee |
|---|---|---|---|---|---|---|
| ISR formalism | collinear QED electron structure function (`PDF:lepton`) plus backward-evolution QED spacelike shower | same structure function, coherent QED antenna shower | collinear QED electron structure function `ThePEG::LeptonLeptonPDF` | collinear QED electron structure function, built-in `PDFESherpa` | YFS soft-photon resummation for the initial state | YFS resummation at amplitude level (CEEX); EEX available |
| Resummation | LL collinear, exponentiated | LL collinear, exponentiated | `O(alpha^2)` soft-photon-exponentiated structure function | exponentiated, `ISR_E_ORDER`/`ISR_E_SCHEME` control the non-leading terms (defaults 1 and 2, the "beta choice") | all-orders soft, with fixed-order improvement up to `O(alpha^3 L^3)` | exact `O(alpha)` and `O(alpha^2)` NLL hard residuals on top of all-orders soft |
| Photon kinematics | shower gives the ISR photons real `pT` | antenna shower gives real `pT` | strictly collinear remnant photons | strictly collinear remnant photons | explicit exclusive photons with exact phase space | explicit exclusive photons with exact phase space |
| Photons in the record | `~2` tagged final photons/event | `~2` | `~1.3` exactly beam-collinear photons/event | exactly 2 beam-collinear remnant photons/event | explicit multi-photon final state | explicit multi-photon final state |
| ISR-FSR interference | no | no | no | no | no (noted as a future target for Sherpa) | **yes**, coherently to all orders |
| Electroweak loops | effective couplings only | effective couplings only | effective couplings only | effective couplings only | effective couplings only | DIZET 6.45 form factors |
| Hadronization | built in | built in | built in | built in | built in | external |

The practical consequence for a Z-pole event-shape correction:

- PYTHIA, Herwig and Sherpa `PDFESherpa` all implement the *same class* of
  model, a collinear QED electron structure function.  They should, and in the
  measurements below do, agree on the mean energy removed from the beams.
- Sherpa's YFS mode and KKMC's CEEX differ from that class in that the photon
  transverse momentum is real, so the thrust axis and the thrust-axis polar
  angle see the photons differently.  For a thrust spectrum integrated over
  photon angle this is a small effect; for an oriented thrust analysis or any
  tagged-photon selection it is not.
- Only KKMC gives IFI and genuine electroweak loop corrections.  IFI is
  suppressed at the Z peak but is the standard systematic for forward-backward
  asymmetries, and it is the main reason KKMC is still the reference.

### 2.1 Measured ISR energy loss in the existing samples

Analytic expectation first.  With
`beta = (2 alpha / pi) (ln(s/m_e^2) - 1) ~ 0.108` at `sqrt(s) = 91.1876 GeV`,
folding the leading-log radiator with the Z Breit-Wigner at the peak gives
`<1 - s'/s> ~ 0.006`, hence a mean ISR photon energy of roughly
**0.2 to 0.3 GeV per event**.  Radiative return is not available at the peak,
so the large inclusive logarithm `ln(m_Z^2/m_e^2) = 24.18` does *not* translate
into a large mean energy loss for this observable.

Measured mean energy in exactly beam-collinear final photons
(`|cos theta| > 0.9999`), which is the model-independent proxy for
structure-function ISR:

| Sample | ISR OFF [GeV] | ISR ON [GeV] |
|---|---|---|
| PYTHIA 8.315 | 0.0036 | 0.077 |
| PYTHIA 8.315 Vincia | 0.0034 | 0.201 |
| Herwig 7.3.0 | **0.194** | 0.197 |
| Sherpa 3.0.3 `PDFESherpa` | 0.0030 | 0.216 |
| Sherpa 3.0.3 YFS | 0.0030 | 0.119 |

and the same conclusion from the ISR-independent side, the mean invariant mass
of the final state after removing neutrinos and beam-collinear photons:

| Sample | `<M_vis>` OFF [GeV] | `<M_vis>` ON [GeV] | shift |
|---|---|---|---|
| PYTHIA 8.315 | 89.355 | 89.281 | -0.074 |
| PYTHIA 8.315 Vincia | 89.410 | 89.216 | -0.194 |
| Herwig 7.3.0 | 89.359 | 89.339 | -0.020 |
| Sherpa 3.0.3 `PDFESherpa` | 89.524 | 89.304 | -0.220 |
| Sherpa 3.0.3 YFS | 89.524 | 89.392 | -0.132 |

The two columns agree per generator, as they must, because
`sqrt(s') ~ sqrt(s) - E_gamma^ISR`.  Sherpa `PDFESherpa` (0.22 GeV) and Vincia
(0.20 GeV) sit on the analytic expectation; default PYTHIA (0.08 GeV) is a
factor of about three lower, which is a genuine and so far unexplained
shower-model difference within PYTHIA, since Vincia uses the same
`PDF:lepton` structure function with a different QED shower.  The Herwig row is
discussed in section 3.3.

---

## 3. Audit of the existing 3M study

### 3.1 The ISR-photon tag was catching every photon in the Sherpa events

**Statement being corrected.** The note and slides report
`<N_gamma^ISR> = 17.8` and `<E_gamma^ISR> = 22.8 GeV` per event for both Sherpa
samples, and describe this as "a large explicit ISR photon system" written by
Sherpa into HepMC.

**What is actually happening.** The Sherpa cards set `HEPMC3_SHORT: true`.  That
option collapses the whole event onto a **single vertex** whose incoming
particles are the two beam leptons.  The converter's HepMC ISR rule in
`src/real_isr_ntuple_producer.cc` tagged a photon when its production vertex had
an incoming `e+-`, so with a short record that rule fires for **every** final
photon in the event, including all pi0 decay photons.

Measured, 3M events, final-state particles only:

| ISR ON sample | `<N_gamma>` all | `<E_gamma>` all [GeV] | `<N>` tagged | `<E>` tagged [GeV] |
|---|---|---|---|---|
| PYTHIA 8.315 | 23.00 | 23.84 | 2.57 | 0.14 |
| PYTHIA 8.315 Vincia | 22.75 | 24.09 | 2.00 | 0.20 |
| Herwig 7.3.0 QED shower | 20.49 | 23.49 | 0.54 | 0.52 |
| Sherpa 3.0.3 `PDFESherpa` | 20.38 | 22.90 | **17.79** | **22.86** |
| Sherpa 3.0.3 YFS | 19.71 | 22.92 | **17.89** | **22.89** |

For Sherpa the tagged photon energy is 99.8% of the *total* photon energy in the
event.  About 23 GeV of photons per Z -> hadrons event is simply the normal pi0
decay content, and roughly 18 photons above 50 MeV is the normal pi0 photon
multiplicity.  The identical result for `PDFESherpa` and YFS, which are very
different ISR models, was itself the clue: the number has nothing to do with
ISR.

**Corrected value.** The true Sherpa ISR photon energy is about
**0.22 GeV/event**, not 22.8 GeV.  Verified directly: a 200-event Sherpa
`PDFESherpa` run with `HEPMC3_SHORT` removed produces a proper HepMC tree
(44 vertices per event instead of 1), and the *unchanged* tagging rule then
returns `<N_ISR> = 2.00` and `<E_ISR> = 0.218 GeV` -- exactly the two collinear
structure-function remnant photons, one per beam.

### 3.2 The tag was applied to ISR ON samples only

`real_isr_ntuple_producer.cc` computed `isrPhoton = b.isrOn && ...`.  Any
observable that *removes* tagged photons therefore removed them from the
numerator and denominator of `C_ISR` asymmetrically: from the ON sample only.

This invalidates thrust definitions **A** and **C** as ISR corrections, for
every generator, independently of the tagging bug in 3.1:

| Last ALEPH bin, `0.99 < T < 1.00` | N (valid) | A (invalid) | B (valid) | C (invalid) |
|---|---|---|---|---|
| PYTHIA 8.315 | 1.166 +- 0.009 | 1.055 | 1.167 +- 0.009 | 1.051 |
| PYTHIA 8.315 Vincia | 1.170 +- 0.010 | 1.018 | 1.170 +- 0.010 | 1.014 |
| Sherpa 3.0.3 `PDFESherpa` | 1.150 +- 0.007 | **0.492** | 1.149 +- 0.007 | **0.494** |
| Sherpa 3.0.3 YFS | 1.141 +- 0.007 | **0.493** | 1.140 +- 0.007 | **0.497** |
| Herwig 7.3.0 QED shower | 1.021 +- 0.005 | 0.903 | 1.021 +- 0.005 | 0.900 |

Definitions N and B are unaffected because they keep all photons, so the tag
cancels out of them.  **The `C_ISR = 0.49` Sherpa endpoint deficit, which the
note and slides treat as the central physics puzzle, is entirely an artefact of
stripping every photon out of the Sherpa ISR ON events and none out of the ISR
OFF events.**

Once the invalid definitions are dropped, there is no Sherpa endpoint anomaly.
The four genuine ISR ON/OFF pairs agree:

```
C_ISR(N), last ALEPH bin:   PYTHIA 1.166   Vincia 1.170   Sherpa 1.150   Sherpa YFS 1.141
```

spread `0.03`, comparable with the `0.007`-`0.010` statistical errors times the
number of comparisons.  Herwig's 1.021 is the outlier, for the reason in 3.3.

### 3.3 Herwig ISR was never switched off

**Statement being corrected.** The note describes the Herwig pair as a
`QCD` versus `QEDQCD` shower comparison that is "the least ISR-only comparison
here".  That is directionally right but understates the problem: the Herwig pair
contains **no ISR toggle at all**.

`snippets/EECollider.in` leaves the ThePEG default e+/e- PDF in place, which is
`ThePEG::LeptonLeptonPDF`.  Reading the implementation
(`ThePEG-2.3.0/PDF/LeptonLeptonPDF.cc`), that object is the
soft-photon-exponentiated collinear QED electron structure function with
`beta/2 = (alpha/pi)(ln(Q^2/m^2) - 1)`, i.e. real beam ISR.  The cards only
change `ShowerHandler:Interactions`, which is the shower, not the beam.

Measured confirmation: the two 3M Herwig samples have the **same** beam-collinear
photon content, 1.29 photons/event carrying 0.194 GeV (OFF) and 0.197 GeV
("ON").  Compare PYTHIA ISR OFF, 0.0036 GeV, where ISR really is off.

**Consequence.** The Herwig curve in every `C_ISR` figure is a QED-shower ratio
at fixed ISR, and should not be read as a weak-ISR generator.  Its near-unity
value is exactly what an ISR-unchanged ratio should give.

**Fix, validated.** Setting the beam PDF explicitly is what toggles ISR.  A
2000-event test at `sqrt(s) = 91.1876 GeV` with the QED shower held fixed at
`QEDQCD`:

| Herwig configuration | beam-collinear photons/event | energy [GeV] |
|---|---|---|
| `e+-:PDF = /Herwig/Partons/LeptonPDF` (ISR ON) | 1.292 | 0.205 |
| `e+-:PDF = /Herwig/Partons/NoPDF` (ISR OFF) | 0.003 | 0.004 |

### 3.4 The "Weizsaecker-Williams" description of older Herwig ISR is wrong

The slide on older reference setups says Herwig 7.1.5 used "an effective
Weizsaecker-Williams style energy-loss treatment".  Herwig/ThePEG does ship a
`ThePEG::WeizsackerWilliamsPDF`, but it is a separate object describing the
**photon** content of a lepton, used for photon-initiated (gamma-gamma)
processes, and it is not attached to the e+/e- beams for `MEee2gZ2qq`.  The
object actually used is `LeptonLeptonPDF`, the QED electron-in-electron
structure function -- the same class of model as PYTHIA's `PDF:lepton` and
Sherpa's `PDFESherpa`, not an equivalent-photon approximation.

### 3.5 The visible-energy explanation is wrong

The note explains the shrinking of the Sherpa/PYTHIA visible-energy gap from
23 GeV to below 0.2 GeV by saying that "most of that explicit Sherpa ISR photon
energy is forward enough to fail the `|eta| < 1.74` acceptance".

That cannot be true, and the numbers say so: with `<E_vis> = 82.04 GeV` inside
`|eta| < 1.74` for Sherpa and 91.19 GeV of total energy, there is only about
9 GeV outside the acceptance plus neutrinos, so 22.8 GeV of forward photons does
not fit.  Only about 0.5 GeV of the tagged photons is at `|cos theta| > 0.99`.
The mis-tagged photons are ordinary central hadronic photons, comfortably inside
the acceptance.  The 23 GeV gap under the old definition came from *removing*
them via `isISRPhoton`; the fix was to stop removing them, not an acceptance
effect.

### 3.6 The quoted mean-thrust shifts are the wrong definition and the wrong sign

The note quotes mean-thrust shifts "for the 3M ntuples" in a paragraph about the
nominal all-final thrust.  Those numbers come from the `thrust` branch of the
`Events` tree, which in these already-produced files is definition **A**, the
broken one, not definition N.  For Sherpa the sign flips.

| Sample | quoted (definition A) | corrected `<T_N>` OFF | corrected `<T_N>` ON | corrected shift |
|---|---|---|---|---|
| PYTHIA 8.315 | 0.93197 -> 0.93163 | 0.932656 | 0.931820 | -0.00084 |
| PYTHIA 8.315 Vincia | 0.93277 -> 0.93265 | 0.933416 | 0.932582 | -0.00083 |
| Herwig 7.3.0 | 0.94003 -> 0.94040 | 0.940606 | 0.940286 | -0.00032 |
| Sherpa `PDFESherpa` | 0.93395 -> **0.93598** | 0.934525 | 0.933686 | **-0.00084** |
| Sherpa YFS | 0.93395 -> **0.93593** | 0.934525 | 0.933702 | **-0.00082** |

Corrected, all four genuine ISR pairs shift the mean nominal thrust by the same
`-0.00083`.  That is a much stronger and much simpler result than the one
currently written up, and it is consistent with the common structure-function
ISR model underlying PYTHIA, Vincia and Sherpa.  Herwig's smaller shift is the
residual QED-shower effect at fixed ISR.

### 3.7 Everything downstream of `M_vis` for Sherpa is contaminated

`Mvis_excluding_ISR_photons` is built by removing the tagged photons, so for
Sherpa it is the hadronic system with all photons deleted:

| Sample, all ISR ON events | `<M_vis>` as stored [GeV] | corrected `<M_vis>` [GeV] |
|---|---|---|
| PYTHIA 8.315 | 89.23 | 89.28 |
| Herwig 7.3.0 | 89.02 | 89.34 |
| Sherpa `PDFESherpa` | **65.76** | 89.30 |
| Sherpa YFS | **65.71** | 89.39 |

Consequently these statements in the note are all artefacts and must be
withdrawn:

- "the last ALEPH thrust bin is still a low-visible-mass ISR ON population",
  with `<M_vis>` around 66 GeV and `<s'_vis/s>` around 0.54;
- "`s'_vis/s > 0.95` gives a last-bin `PDFESherpa` ratio of about 36.7 with
  OFF/ON survival fractions of 0.839 and 0.0052" -- the ON sample fails the cut
  because its photons were deleted, not because its `s'` is low;
- "Sherpa ISR ON ... has a much lower visible mass spectrum in the current
  ntuples, so the OFF/ON ratio probes a different selected phase space";
- the slide table giving `<E_vis>` as 89.3 GeV for PYTHIA and 66.8 GeV for
  Sherpa.  Under the current definition both are 82.0 GeV.

The whole cut-scan section, and the `endpoint_diagnostics` and `cut_scans`
figures, are unreliable for Sherpa and mildly biased for Herwig and PYTHIA until
the diagnostics are regenerated.

### 3.8 Smaller items

- `nISRPhotons` and `totalISRPhotonEnergy` in the PYTHIA path counted
  non-final-state photons as well, inflating them from 2.57 / 0.140 GeV
  (final only) to 2.73 / 0.202 GeV.  Fixed by requiring `isFinal`.
- The `|cos theta| > 0.990` leg of the HepMC rule has a large hadronic
  background: in ISR OFF samples it selects about 0.33 GeV/event of ordinary
  forward photons, larger than the ISR signal itself.  Tightened to
  `|cos theta| > 0.9999`.
- `ISR_E_ORDER: 1` and `ISR_E_SCHEME: 2` are described as chosen settings; they
  are the Sherpa defaults (order in alpha = 1, the "beta choice" scheme).
- Reference tidy-ups: Sherpa 3 is JHEP 12 (2024) 156, arXiv:2410.22148; Herwig
  7.3 is Eur. Phys. J. C 84 (2024) 1053, arXiv:2312.05175; arXiv:2203.12557 is
  correctly cited (Frixione, Laenen et al., initial-state QED radiation for
  future e+e- colliders).

---

## 4. Fixes applied in this repository

| File | Change |
|---|---|
| `src/real_isr_ntuple_producer.cc` | `hepmcLooksLikeISRPhoton` now ignores the incoming-lepton test on vertices with more than 4 outgoing particles, so a flat/short HepMC record cannot tag the whole event; the loose `\|cos theta\| > 0.990` leg is tightened to `0.9999`; the tag is evaluated for ISR OFF samples too; the PYTHIA tag requires `isFinal`. |
| `cards/sherpa_zpole_hepmc_ISR_{OFF,ON,YFS}.yaml.template` | `HEPMC3_SHORT: true` removed, with a comment explaining why. |
| `cards/herwig_zpole_hepmc_ISR_OFF.in.template` | Now genuinely ISR OFF: `e+-:PDF = /Herwig/Partons/NoPDF`, QED shower held at `QEDQCD`. |
| `cards/herwig_zpole_hepmc_ISR_ON.in.template` | New: genuine ISR ON partner, `e+-:PDF = /Herwig/Partons/LeptonPDF`, same shower settings. |
| `cards/herwig_zpole_hepmc_QEDshower.in.template` | Kept, with a header warning that it is not an ISR toggle. |
| `scripts/run_real_isr_production.sh`, `scripts/run_real_isr_production_10worker.sh` | `run_herwig` takes `OFF`/`ON`/`QEDshower`; a `herwig-on` target now exists. |

Validation of the fixes (small samples, this session):

```
Sherpa PDFESherpa, HEPMC3_SHORT removed, 200 events:
    <N_ISR> = 2.000   <E_ISR> = 0.218 GeV      (was 17.79 / 22.86)
Herwig, QED shower fixed at QEDQCD, 2000 events:
    ISR ON  (LeptonPDF)  <N_ISR> = 1.292   <E_ISR> = 0.2050 GeV
    ISR OFF (NoPDF)      <N_ISR> = 0.003   <E_ISR> = 0.0042 GeV
```

## 5. What still has to be re-run

The HepMC intermediates for the 3M production were not kept, so the Sherpa and
Herwig ROOT ntuples cannot be re-tagged offline; they have to be regenerated.

1. Regenerate the Sherpa OFF/ON/YFS and Herwig OFF/ON samples with the corrected
   cards (PYTHIA and Vincia are unaffected except for the `isFinal` fix, which
   changes only the `nISRPhotons`/`totalISRPhotonEnergy` bookkeeping branches).
2. Re-run `make_endpoint_diagnostics` and the plotting macros.
3. Results that survive unchanged: the nominal `C_ISR(T_N)` and `C_ISR(T_B)`
   curves and the ALEPH comparisons, because they do not use the tag.  Results
   that must be regenerated before being quoted: everything using
   `Mvis_excluding_ISR_photons`, `sPrime_vis_over_s`, `betaZ_vis`, the tagged
   photon spectra, the cut scans, and definitions A and C.

## 6. KKMC production

KKMC is now built and run in this study.  The build is the Fortran line of the
KKMCee repository, which carries the same CEEX physics as the C++ KKMCee 5
release and has JETSET/PYTHIA hadronization built in, so no external
hadronization step is needed.

### 6.1 Build

```bash
git clone https://github.com/KrakowHEPSoft/KKMCee.git
cd KKMCee && ln -sfn dizet-6.45 dizet
cd ffbench
# gfortran 13 needs the legacy dialect; the master KKMakefile FFLAGS line was
# extended with -std=legacy -fallow-argument-mismatch
make -f KKMakefile makflag
make -f KKMakefile makprod
make -f KKMakefile EWtables
make -f KKMakefile KKdump.exe
```

`ffbench/KKdump.f` is an addition of this study.  It runs the standard KKMC
event loop and writes, per event, KKMC's own initial-state photon list taken
from `KarLud_GetPhotons` followed by the stable final state from the PYTHIA
6.202 `PYJETS` record.  That makes KKMC the only sample in this study whose ISR
photon set is the generator's own definition rather than an analysis-side tag.

### 6.2 Configurations

`scripts/run_kkmc_production.sh` drives three 1M-event configurations at
`sqrt(s) = 91.1876 GeV`, `e+e- -> q qbar` with `d,u,s,c,b` only:

| Mode | `KeyISR` | `KeyINT` | Purpose |
|---|---|---|---|
| `ISR_OFF` | 0 | 0 | ISR toggle partner |
| `ISR_ON` | 1 | 0 | nominal, comparable with the other generators |
| `ISR_ON_IFI` | 1 | 2 | ISR-FSR interference, KKMC default |

`KeyGPS = 1` (CEEX), `KeyFSR = 1`, `KeyQSR = 1` and `KeyHad = 1` are identical in
all three, so `ISR_OFF` versus `ISR_ON` isolates beam radiation and
`ISR_ON_IFI` versus `ISR_ON` isolates interference.  Interference is switched
off in the nominal pair because no other generator in this study implements it,
so leaving it on would confound the comparison.

`--mode kkmc` was added to `real_isr_ntuple_producer` to read the dump into the
standard `Events` schema; `SAMPLE_SET=KKMC` was added to
`scripts/run_endpoint_diagnostics.sh`.

### 6.3 Validation

- Cross section: KKMC returns about 30.4 nb for the hadronic cross section at
  the peak.  The LEP pole cross section is `sigma^0_had = 41.48 nb` and
  radiation reduces the observed peak by roughly a quarter, so this is the
  expected value.
- Energy: the summed final-state energy is 91.1876 GeV event by event.
- ISR energy: about 0.18 GeV per event in a 2000-event test, matching the
  analytic Z-peak expectation of 0.2 to 0.3 GeV and the other generators.

### 6.4 Comparison outputs

`macros/plot_isr_model_comparison.C` builds the cross-generator comparison.  It
deliberately avoids the `isISRPhoton` tag everywhere: radiation is measured from
beam-collinear photons and from the shift in the reduced hadronic mass, and
`C_ISR(T)` uses the all-particle thrust.  For KKMC it additionally reads the
generator's own ISR photon list, which calibrates how much of the real radiated
energy the geometric estimator captures.

Slides: `overleaf/slides/20260921-isr_model_comparison.tex`, with the KKMC
numbers substituted from the macro outputs by
`scripts/fill_kkmc_slide_numbers.py` so that the deck cannot drift from the
files it was made from.

### 6.5 Results

Validation, 1M events per configuration:

| Quantity | KKMC | Expected |
|---|---|---|
| cross section, ISR off | 41.33 +- 0.01 nb | LEP pole 41.48 nb |
| cross section, ISR on | 30.38 +- 0.01 nb | 26.5% radiative reduction |
| summed final-state energy | 91.1876 GeV per event | exact |

Radiation:

| Quantity | KKMC ISR off | KKMC ISR on |
|---|---|---|
| generator ISR photon energy | 0 | **0.220 GeV** |
| beam-collinear estimator | 0.0026 GeV | 0.136 GeV |
| `<M_vis>` | 89.545 GeV | 89.423 GeV |
| `<T_N>` | 0.93282 | 0.93212 |

KKMC agrees with Sherpa `PDFESherpa` (0.216 GeV) on the amount of energy
radiated.  The beam-collinear estimator recovers only 62% of it, and that 38%
deficit is the radiation carrying transverse momentum, which a collinear
structure function cannot produce.  The estimator is therefore unbiased for
PYTHIA, Vincia, Herwig and Sherpa `PDFESherpa` and a lower limit for Sherpa YFS
and KKMC, which explains why Sherpa YFS appeared to radiate only 0.119 GeV.

Endpoint correction, last ALEPH thrust bin, definition N:

| Sample | `C_ISR` | `<T>` shift |
|---|---|---|
| PYTHIA 8.315 | 1.166 +- 0.009 | -0.00084 |
| PYTHIA 8.315 Vincia | 1.170 +- 0.010 | -0.00083 |
| Sherpa 3.0.3 `PDFESherpa` | 1.150 +- 0.007 | -0.00084 |
| Sherpa 3.0.3 YFS | 1.141 +- 0.007 | -0.00082 |
| **KKMC 4.30 CEEX** | **1.097 +- 0.015** | **-0.00070** |
| Herwig 7.3.0, ISR not toggled | 1.021 +- 0.005 | -0.00032 |

**This is the main physics result of the comparison.**  The four collinear
models agree with each other; KKMC sits 3 to 4 standard deviations below them.
The KKMC to PYTHIA difference is `0.069 +- 0.017`, about 6% of the correction,
and it is the first defensible radiation-model uncertainty in this study.

It is not a hadronization effect.  The KKMC and PYTHIA ISR OFF samples agree on
the mean thrust to `2e-4` and on the endpoint bin population to 3%, despite
KKMC using PYTHIA 6.202 and the other samples PYTHIA 8.315.  The models agree on
how much energy is radiated and disagree on where it goes: with an all-particle
thrust definition, photons carrying transverse momentum migrate events out of
the two-jet endpoint differently from photons pinned to the beam axis.

Below `T = 0.97` every model including KKMC agrees within 1%.

ISR-FSR interference, measured for the first time in this study:

| Quantity | `N(KeyINT=2) / N(KeyINT=0)` |
|---|---|
| integrated over the thrust spectrum | 1.0000 +- 0.0014 |
| last ALEPH bin | 0.987 +- 0.014 |

Interference is consistent with no effect on the thrust spectrum at the Z peak,
which is the expected suppression on a narrow resonance.  It can be dropped from
the uncertainty budget for this observable.

## 7. Consistency audit against ALEPH_Agentic_Event_Shape_Analysis

Performed after the first version of the comparison deck, at the request to
check the code, the ratios, the KKMC version, and the thrust definition.

### 7.1 Thrust definition: was inconsistent, now fixed

The ALEPH event-shape analysis corrects its unfolded thrust with
`C_ISR(j) = N_noISR(j) / N_ISR(j)` built from the `tgenBefore/thrust` branch of
`Isr/isr0_ALL.root` and `Isr/isr1_ALL.root` (2.5M events each).  That branch is
not documented in the note, so the particle selection was determined by
recomputing thrust from the stored particles under several definitions and
comparing:

| Particle set | events matching the stored branch | mean abs. difference |
|---|---|---|
| all final-state particles | 2358 / 3000 | 9.6e-4 |
| **final state excluding neutrinos** | **3000 / 3000** | **1.5e-8** |

So the ALEPH analysis thrust is all stable final-state particles **excluding
neutrinos**, photons included.  Final state means Pythia status > 0; the tree
also stores decayed entries, which must be skipped.

The first version of the comparison used definition N, all stable final-state
particles **including** neutrinos, and said so on the slides.  That was the
wrong observable.  The comparison now uses definition B, which is the ALEPH
definition.  Numerically the change is small, at most 0.002 in `C_ISR`, so no
conclusion moves, but the stated definition was wrong and is corrected.

| Sample | `C_ISR` last bin, def N (was) | def B (ALEPH, now) |
|---|---|---|
| PYTHIA 8.315 | 1.166 +- 0.009 | 1.167 +- 0.009 |
| PYTHIA Vincia | 1.170 +- 0.010 | 1.170 +- 0.010 |
| Sherpa PDFESherpa | 1.150 +- 0.007 | 1.149 +- 0.007 |
| Sherpa YFS | 1.141 +- 0.007 | 1.140 +- 0.007 |
| KKMC CEEX | 1.097 +- 0.015 | 1.094 +- 0.015 |
| Herwig, ISR not toggled | 1.021 +- 0.005 | 1.021 +- 0.005 |

Binning was also checked: the published ALEPH thrust table
(HEPData ins636645 Table 54) is uniform from 0.58 to 1.00 in steps of 0.01, which
is what the macro assumes.

### 7.2 Cross-check against the correction the analysis applies

Recomputing `C_ISR` from their own two files reproduces their cached
`Isr/isr_corr_precomputed.root` exactly, including 1.1996 in the last bin, so
their cache is built by binning the branch as T directly.

Comparing their Pythia8 correction with ours, in the same definition and bins:

| bin | ALEPH analysis | this study, PYTHIA 8.315 | difference |
|---|---|---|---|
| 0.94-0.95 | 0.9651 +- 0.0031 | 0.9754 +- 0.0029 | -2.5 sigma |
| 0.95-0.96 | 0.9756 +- 0.0027 | 0.9766 +- 0.0025 | -0.3 sigma |
| 0.96-0.97 | 0.9859 +- 0.0024 | 0.9872 +- 0.0022 | -0.4 sigma |
| 0.97-0.98 | 1.0204 +- 0.0022 | 1.0189 +- 0.0020 | +0.5 sigma |
| 0.98-0.99 | 1.0905 +- 0.0030 | 1.0840 +- 0.0027 | +1.6 sigma |
| 0.99-1.00 | 1.1996 +- 0.0108 | 1.1673 +- 0.0092 | +2.3 sigma |

Agreement is good through the bulk but drifts at the endpoint, in the same
direction in the last two bins.  That is larger than statistics comfortably
explains and is unresolved: the generation settings for their samples are not in
that repository, so the Pythia version, tune, ISR-off switches and flavour
selection cannot be compared.  This matters for the headline number, because
measured against the correction actually in use the KKMC difference grows from
0.073 +- 0.017 to 0.105 +- 0.019.

### 7.3 Two documentation problems in that repository, for its owner

- `sections/Datasets.tex` states the ISR correction "is not yet implemented in
  the current analysis", while `sections/AnalysisMethod.tex` and
  `sections/SystematicUncertainties.tex` describe it as applied.  One of these
  is stale.
- `sections/AnalysisMethod.tex` says the `tgenBefore/thrust` branch is "stored
  as `1-T`, transformed back to T".  The branch holds T directly (mean 0.931,
  range 0.56 to 0.999), and the cached correction is consistent with binning T
  without any transformation.  The computation looks right and the sentence
  looks wrong, but it should be reconciled so nobody applies the transformation.

### 7.4 KKMC version: mislabelled, and not the newest release

The banner of the built program reads `Version 4.30, October 2020`, not 4.24 as
`configure.ac` suggests and as the first version of the deck claimed.  All
labels are corrected to 4.30.

The current release is **KKMCee v5.00.02**, a C++ rewrite.  Its release notes
describe the rewrite, FOAM, HepMC3 output and speed, with no change of physics
content, and the CPC paper states that 5.00 reproduces the Fortran benchmarks.
The repository master branch carries only the Fortran tree, which is why 4.30 is
what got built.  Repeating the measurement on v5.00.02 remains an open check.

CEEX was verified to be active for quarks rather than silently falling back to
EEX: `KK2f/KK2f.f` gates it on `m_KeyGPS != 0 && SvarQ > MminCEEX^2` with
`MminCEEX = 20 GeV` from the quark entries of `.KK2f_defaults`, and at the Z pole
`SvarQ` is about 8300 GeV^2.

### 7.5 Checks that passed

- The `C_ISR` ratio and its independent-Poisson error were reproduced by a
  separate macro from the same inputs.
- KKMC cross sections bracket the radiative correction correctly:
  41.33 +- 0.01 nb with radiation off against the LEP pole value 41.48 nb, and
  30.38 +- 0.01 nb with radiation on, a 26.5% reduction.
- Each KKMC run directory's `pro.input` was checked against its own output
  banner and dump header, confirming that `KeyISR` and `KeyINT` are what the
  sample names claim.

## 8. References

- S. Jadach, B.F.L. Ward, Z. Was, *The precision Monte Carlo event generator KK
  for two-fermion final states in e+e- collisions*, Comput. Phys. Commun. 130
  (2000) 260, arXiv:hep-ph/9912214.
- S. Jadach, B.F.L. Ward, Z. Was, S.A. Yost, *KKMC 4.22 electroweak updates*,
  Comput. Phys. Commun. 260 (2021) 107734, arXiv:2007.07964.
- S. Jadach, B.F.L. Ward, Z. Was, S.A. Yost, A. Siodmok, *Multi-photon Monte
  Carlo event generator KKMCee*, Comput. Phys. Commun. 283 (2023) 108556,
  arXiv:2204.11949; release notes arXiv:2410.07294.
- C. Bierlich et al., *A comprehensive guide to the physics and usage of PYTHIA
  8.3*, SciPost Phys. Codebases 8 (2022).
- G. Bewick et al., *Herwig 7.3 release note*, Eur. Phys. J. C 84 (2024) 1053,
  arXiv:2312.05175.
- E. Bothmann et al., *Event generation with Sherpa 3*, JHEP 12 (2024) 156,
  arXiv:2410.22148.
- B.F.L. Ward et al. / F. Krauss, M. Schoenherr, P. Price, *YFS resummation for
  future lepton-lepton colliders*, arXiv:2203.10948.
- S. Frixione, E. Laenen et al., *Initial state QED radiation aspects for future
  e+e- colliders*, arXiv:2203.12557.
- Sherpa ISR parameters:
  <https://sherpa-team.gitlab.io/sherpa/v3.0.0/manual/parameters/isr.html>
