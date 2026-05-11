# ROOT Tree Schema

Each generated ROOT file contains one `TTree` named `Events`.

## Event-Level Branches

| Branch | Type | Meaning |
|---|---|---|
| `eventId` | `int` | Event index within the sample. |
| `generatorId` | `int` | Integer generator-equivalent configuration ID. |
| `generatorName` | `std::string` | Human-readable generator-equivalent label. |
| `isrOn` | `bool` | ISR category flag. |
| `processId` | `int` | Toy hard-process category. |
| `processName` | `std::string` | Human-readable process category. |
| `sqrtS` | `double` | Nominal collision energy in GeV. |
| `weight` | `double` | Event weight used in histograms. |
| `visibleEnergy` | `double` | Scalar visible energy after the default selection. |
| `missingEnergy` | `double` | Nominal missing energy from invisible particles and ISR. |
| `thrust` | `double` | Thrust from selected final-state visible particles. |
| `thrustMajor` | `double` | Approximate thrust-major observable. |
| `thrustMinor` | `double` | Approximate thrust-minor observable. |
| `sphericity` | `double` | Sphericity from the visible final state. |
| `cParameter` | `double` | C-parameter from the visible final state. |
| `nFinalParticles` | `int` | Stable final-state particle count. |
| `nChargedFinalParticles` | `int` | Charged stable final-state particle count. |
| `nPartons` | `int` | Stored hard/shower parton count. |
| `nISRPhotons` | `int` | Number of flagged ISR photons. |
| `totalISRPhotonEnergy` | `double` | Scalar sum of flagged ISR photon energies. |

## Particle-Level Branches

| Branch | Type | Meaning |
|---|---|---|
| `pdgId` | `std::vector<int>` | PDG particle ID. |
| `status` | `std::vector<int>` | Event-record status code. |
| `mother1`, `mother2` | `std::vector<int>` | Mother index range. |
| `daughter1`, `daughter2` | `std::vector<int>` | Daughter index range. |
| `px`, `py`, `pz` | `std::vector<float>` | Momentum components in GeV. |
| `energy` | `std::vector<float>` | Energy in GeV. |
| `mass` | `std::vector<float>` | Particle mass in GeV. |
| `charge` | `std::vector<float>` | Electric charge in units of `e`. |
| `isFinal` | `std::vector<char>` | Stable final-state particle flag. |
| `isCharged` | `std::vector<char>` | Charged particle flag. |
| `isParton` | `std::vector<char>` | Hard/shower parton flag. |
| `isHadron` | `std::vector<char>` | Hadron flag. |
| `isLepton` | `std::vector<char>` | Lepton flag. |
| `isPhoton` | `std::vector<char>` | Photon flag. |
| `isISRPhoton` | `std::vector<char>` | Explicit ISR photon flag. |

## Default Analysis Selection

The default `SelectionConfig` is:

```cpp
includeISRPhotonsInVisibleEnergy = true;
includeISRPhotonsInEventShapes = true;
includeNeutrinosInVisibleEnergy = false;
includeNeutrinosInEventShapes = true;
useChargedParticlesOnly = false;
visibleEtaMax = 1.74;
```

Visible energy is the scalar energy sum of final-state particles within
`|eta| < 1.74`, excluding neutrinos and including photons.  The current
nominal thrust quick check uses all stable final-state particles with no
`eta` restriction, including neutrinos and tagged ISR photons.  ISR photons
are stored in the event record so the visible-particle cross-checks can remove
or include them explicitly.

## Derived Endpoint Diagnostics Tree

The endpoint diagnostic pass writes separate ROOT files with a tree named
`EndpointDiagnostics`.  These files are derived from the original `Events`
trees and leave the generator ntuples unchanged.  The nominal diagnostic is
all stable final-state particles, including neutrinos and tagged ISR photons;
the visible-system branches are retained for ALEPH-style cross-checks.

| Branch | Type | Meaning |
|---|---|---|
| `eventId` | `int` | Event index copied from the source tree. |
| `generatorName` | `std::string` | Source generator label. |
| `isrOn` | `bool` | Source ISR category flag. |
| `T_lab_allFinal_including_ISR_photons` | `double` | Nominal lab-frame thrust from all stable final-state particles, including neutrinos and tagged ISR photons. |
| `T_lab_allFinal_excluding_ISR_photons` | `double` | Lab-frame thrust from all stable final-state particles after removing tagged ISR photons. |
| `T_allFinalCM_including_ISR_photons` | `double` | Thrust after boosting the all-final-particle system to its own rest frame. |
| `cosTheta_thrust_lab_allFinal_including_ISR_photons` | `double` | Cosine of the nominal all-final lab-frame thrust-axis polar angle. |
| `absCosTheta_thrust_lab_allFinal_including_ISR_photons` | `double` | Absolute value of the nominal all-final lab-frame thrust-axis cosine. |
| `Eall_including_ISR_photons` | `double` | Scalar energy of all stable final-state particles, including tagged ISR photons. |
| `Eall_excluding_ISR_photons` | `double` | Scalar energy of all stable final-state particles after removing tagged ISR photons. |
| `Mall_including_ISR_photons` | `double` | Invariant mass of all stable final-state particles, including tagged ISR photons. |
| `Mall_excluding_ISR_photons` | `double` | Invariant mass of all stable final-state particles after removing tagged ISR photons. |
| `sAll_including_ISR_over_s` | `double` | `Mall_including_ISR_photons^2 / s`. |
| `betaZ_all_including_ISR` | `double` | `abs(Pz_all) / Eall` for the nominal all-final system. |
| `T_lab_excluding_ISR_photons` | `double` | Visible-particle lab-frame thrust: final-state non-neutrinos excluding tagged ISR photons. |
| `T_lab_including_ISR_photons` | `double` | Visible-particle lab-frame thrust after adding tagged ISR photons. |
| `T_visibleCM_excluding_ISR_photons` | `double` | Thrust after boosting the visible non-neutrino, ISR-photon-removed system to its rest frame. |
| `T_lab_hadron_excluding_ISR_photons` | `double` | Lab-frame thrust from final-state hadrons only, excluding tagged ISR photons and leptons. |
| `T_lab_hadron_including_ISR_photons` | `double` | Lab-frame thrust from final-state hadrons plus tagged ISR photons. |
| `T_hadronicCM_excluding_ISR_photons` | `double` | Thrust after boosting the final-state-hadron, ISR-photon-removed system to its rest frame. |
| `Evis_excluding_ISR_photons` | `double` | Scalar energy of the visible non-neutrino system after removing tagged ISR photons. |
| `Evis_including_ISR_photons` | `double` | Scalar energy of visible non-neutrino particles plus tagged ISR photons. |
| `Mvis_excluding_ISR_photons` | `double` | Invariant mass of the visible non-neutrino system after removing tagged ISR photons. |
| `Mvis_including_ISR_photons` | `double` | Invariant mass of visible non-neutrino particles plus tagged ISR photons. |
| `sPrime_vis_over_s` | `double` | `Mvis_excluding_ISR_photons^2 / s`. |
| `betaZ_vis` | `double` | `abs(Pz_vis) / Evis` for the final-state-hadron system. |
| `leading_ISR_photon_energy` | `double` | Leading tagged ISR photon energy in GeV, or `-1` if absent. |
| `leading_ISR_photon_theta` | `double` | Leading tagged ISR photon polar angle in radians, or `-1` if absent. |
| `N_ISR_photons` | `int` | Number of tagged ISR photons. |
| `sum_ISR_photon_energy` | `double` | Scalar sum of tagged ISR photon energies in GeV. |
| `event_weight` | `double` | Source event weight copied for reference. |
