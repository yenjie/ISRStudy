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
includeISRPhotonsInVisibleEnergy = false;
includeISRPhotonsInEventShapes = false;
useChargedParticlesOnly = false;
```

Visible energy and event-shape calculations exclude neutrinos.  ISR photons are
stored in the event record but are excluded from visible energy and event shapes
unless the selection flags are changed.

