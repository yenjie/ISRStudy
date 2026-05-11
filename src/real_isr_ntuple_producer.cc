#include <Pythia8/Pythia.h>
#include <Pythia8/Vincia.h>

#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/ReaderAsciiHepMC2.h>

#include <TFile.h>
#include <TLorentzVector.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string mode = "pythia";
    std::string hepmcFormat = "hepmc2";
    std::string hepmcInput;
    std::string output = "events.root";
    std::string generatorName = "PYTHIA 8.315";
    int generatorId = 2;
    int nEvents = 100000;
    int seed = 1200510;
    bool isrOn = false;
    bool vincia = false;
    double sqrtS = 91.1876;
};

struct SelectionConfig {
    bool includeISRPhotonsInVisibleEnergy = true;
    bool includeISRPhotonsInEventShapes = false;
    bool useChargedParticlesOnly = false;
    double visibleEtaMax = 1.74;
};

struct EventBranches {
    int eventId = 0;
    int generatorId = 0;
    std::string generatorName;
    bool isrOn = false;
    int processId = 0;
    std::string processName;
    double sqrtS = 91.2;
    double weight = 1.0;
    double visibleEnergy = 0.0;
    double missingEnergy = 0.0;
    double thrust = 0.0;
    double thrustMajor = 0.0;
    double thrustMinor = 0.0;
    double sphericity = 0.0;
    double cParameter = 0.0;
    int nFinalParticles = 0;
    int nChargedFinalParticles = 0;
    int nPartons = 0;
    int nISRPhotons = 0;
    double totalISRPhotonEnergy = 0.0;

    std::vector<int> pdgId;
    std::vector<int> status;
    std::vector<int> mother1;
    std::vector<int> mother2;
    std::vector<int> daughter1;
    std::vector<int> daughter2;

    std::vector<float> px;
    std::vector<float> py;
    std::vector<float> pz;
    std::vector<float> energy;
    std::vector<float> mass;
    std::vector<float> charge;

    std::vector<char> isFinal;
    std::vector<char> isCharged;
    std::vector<char> isParton;
    std::vector<char> isHadron;
    std::vector<char> isLepton;
    std::vector<char> isPhoton;
    std::vector<char> isISRPhoton;
};

bool isNeutrino(int pdgId)
{
    const int apdg = std::abs(pdgId);
    return apdg == 12 || apdg == 14 || apdg == 16;
}

bool isLeptonPdg(int pdgId)
{
    const int apdg = std::abs(pdgId);
    return apdg == 11 || apdg == 12 || apdg == 13 || apdg == 14 || apdg == 15 || apdg == 16;
}

bool isHadronPdg(int pdgId)
{
    const int apdg = std::abs(pdgId);
    return apdg > 100 || apdg == 2212 || apdg == 2112;
}

double particleCharge(int pdgId)
{
    const int apdg = std::abs(pdgId);
    if (apdg == 2 || apdg == 4 || apdg == 6) return pdgId > 0 ? 2.0 / 3.0 : -2.0 / 3.0;
    if (apdg == 1 || apdg == 3 || apdg == 5) return pdgId > 0 ? -1.0 / 3.0 : 1.0 / 3.0;
    if (apdg == 11 || apdg == 13 || apdg == 15) return pdgId > 0 ? -1.0 : 1.0;
    if (apdg == 211 || apdg == 321 || apdg == 2212) return pdgId > 0 ? 1.0 : -1.0;
    if (apdg == 24) return pdgId > 0 ? 1.0 : -1.0;
    return 0.0;
}

void reset(EventBranches& b)
{
    const int keepEventId = b.eventId;
    const int keepGeneratorId = b.generatorId;
    const std::string keepGeneratorName = b.generatorName;
    const bool keepIsrOn = b.isrOn;
    const double keepSqrtS = b.sqrtS;

    b = EventBranches();
    b.eventId = keepEventId;
    b.generatorId = keepGeneratorId;
    b.generatorName = keepGeneratorName;
    b.isrOn = keepIsrOn;
    b.sqrtS = keepSqrtS;
}

void bookTree(TTree* tree, EventBranches& b)
{
    tree->Branch("eventId", &b.eventId);
    tree->Branch("generatorId", &b.generatorId);
    tree->Branch("generatorName", &b.generatorName);
    tree->Branch("isrOn", &b.isrOn);
    tree->Branch("processId", &b.processId);
    tree->Branch("processName", &b.processName);
    tree->Branch("sqrtS", &b.sqrtS);
    tree->Branch("weight", &b.weight);
    tree->Branch("visibleEnergy", &b.visibleEnergy);
    tree->Branch("missingEnergy", &b.missingEnergy);
    tree->Branch("thrust", &b.thrust);
    tree->Branch("thrustMajor", &b.thrustMajor);
    tree->Branch("thrustMinor", &b.thrustMinor);
    tree->Branch("sphericity", &b.sphericity);
    tree->Branch("cParameter", &b.cParameter);
    tree->Branch("nFinalParticles", &b.nFinalParticles);
    tree->Branch("nChargedFinalParticles", &b.nChargedFinalParticles);
    tree->Branch("nPartons", &b.nPartons);
    tree->Branch("nISRPhotons", &b.nISRPhotons);
    tree->Branch("totalISRPhotonEnergy", &b.totalISRPhotonEnergy);
    tree->Branch("pdgId", &b.pdgId);
    tree->Branch("status", &b.status);
    tree->Branch("mother1", &b.mother1);
    tree->Branch("mother2", &b.mother2);
    tree->Branch("daughter1", &b.daughter1);
    tree->Branch("daughter2", &b.daughter2);
    tree->Branch("px", &b.px);
    tree->Branch("py", &b.py);
    tree->Branch("pz", &b.pz);
    tree->Branch("energy", &b.energy);
    tree->Branch("mass", &b.mass);
    tree->Branch("charge", &b.charge);
    tree->Branch("isFinal", &b.isFinal);
    tree->Branch("isCharged", &b.isCharged);
    tree->Branch("isParton", &b.isParton);
    tree->Branch("isHadron", &b.isHadron);
    tree->Branch("isLepton", &b.isLepton);
    tree->Branch("isPhoton", &b.isPhoton);
    tree->Branch("isISRPhoton", &b.isISRPhoton);
}

std::vector<TLorentzVector> selectedVisibleParticles(const EventBranches& b, const SelectionConfig& config, bool forEventShapes)
{
    std::vector<TLorentzVector> selected;
    selected.reserve(b.pdgId.size());
    for (size_t i = 0; i < b.pdgId.size(); ++i) {
        if (!b.isFinal[i]) continue;
        if (isNeutrino(b.pdgId[i])) continue;
        if (forEventShapes && !config.includeISRPhotonsInEventShapes && b.isISRPhoton[i]) continue;
        if (!forEventShapes && !config.includeISRPhotonsInVisibleEnergy && b.isISRPhoton[i]) continue;
        if (config.useChargedParticlesOnly && !b.isCharged[i]) continue;
        TLorentzVector p4(b.px[i], b.py[i], b.pz[i], b.energy[i]);
        if (p4.P() <= 1e-9) continue;
        selected.push_back(p4);
    }
    return selected;
}

double computeVisibleEnergy(const EventBranches& b, const SelectionConfig& config)
{
    double evis = 0.0;
    for (size_t i = 0; i < b.pdgId.size(); ++i) {
        if (!b.isFinal[i]) continue;
        if (isNeutrino(b.pdgId[i])) continue;
        if (!config.includeISRPhotonsInVisibleEnergy && b.isISRPhoton[i]) continue;
        if (config.useChargedParticlesOnly && !b.isCharged[i]) continue;
        TLorentzVector p4(b.px[i], b.py[i], b.pz[i], b.energy[i]);
        if (p4.P() <= 1e-9) continue;
        if (config.visibleEtaMax > 0.0 && std::abs(p4.Eta()) >= config.visibleEtaMax) continue;
        evis += p4.E();
    }
    return evis;
}

TVector3 computeThrustAxis(const std::vector<TLorentzVector>& particles, double* thrustValue = nullptr)
{
    double pSum = 0.0;
    std::vector<TVector3> momenta;
    momenta.reserve(particles.size());
    for (const TLorentzVector& p4 : particles) {
        const TVector3 p = p4.Vect();
        const double mag = p.Mag();
        if (mag <= 1e-9) continue;
        pSum += mag;
        momenta.push_back(p);
    }

    if (pSum <= 1e-9 || momenta.empty()) {
        if (thrustValue) *thrustValue = 0.0;
        return TVector3(0.0, 0.0, 1.0);
    }

    TVector3 bestVector(0.0, 0.0, 0.0);
    double bestMag2 = -1.0;
    const size_t n = momenta.size();

    if (n <= 3) {
        const unsigned int nMasks = 1u << n;
        for (unsigned int mask = 0; mask < nMasks; ++mask) {
            TVector3 sum(0.0, 0.0, 0.0);
            for (size_t i = 0; i < n; ++i) sum += ((mask >> i) & 1u) ? momenta[i] : -momenta[i];
            const double mag2 = sum.Mag2();
            if (mag2 > bestMag2) {
                bestMag2 = mag2;
                bestVector = sum;
            }
        }
    } else {
        for (size_t i = 1; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                const TVector3 cross = momenta[i].Cross(momenta[j]);
                TVector3 ptot(0.0, 0.0, 0.0);
                for (size_t k = 0; k < n; ++k) {
                    if (k == i || k == j) continue;
                    ptot += (momenta[k].Dot(cross) > 0.0) ? momenta[k] : -momenta[k];
                }
                const TVector3 candidates[4] = {
                    ptot - momenta[j] - momenta[i],
                    ptot - momenta[j] + momenta[i],
                    ptot + momenta[j] - momenta[i],
                    ptot + momenta[j] + momenta[i]
                };
                for (const TVector3& candidate : candidates) {
                    const double mag2 = candidate.Mag2();
                    if (mag2 > bestMag2) {
                        bestMag2 = mag2;
                        bestVector = candidate;
                    }
                }
            }
        }
    }

    const double best = std::min(1.0, std::max(0.0, bestVector.Mag() / pSum));
    if (thrustValue) *thrustValue = best;
    return bestVector.Mag() > 1e-12 ? bestVector.Unit() : momenta.front().Unit();
}

double computeThrust(const std::vector<TLorentzVector>& particles)
{
    double thrust = 0.0;
    computeThrustAxis(particles, &thrust);
    return thrust;
}

double computeThrustMajor(const std::vector<TLorentzVector>& particles)
{
    double thrust = 0.0;
    const TVector3 axis = computeThrustAxis(particles, &thrust);
    double denom = 0.0;
    std::vector<TVector3> transverse;
    for (const TLorentzVector& p4 : particles) {
        TVector3 p = p4.Vect();
        const double mag = p.Mag();
        if (mag <= 1e-9) continue;
        denom += mag;
        p -= p.Dot(axis) * axis;
        if (p.Mag() > 1e-9) transverse.push_back(p);
    }
    if (denom <= 1e-9 || transverse.empty()) return 0.0;
    double best = 0.0;
    for (const TVector3& candidate : transverse) {
        const TVector3 n = candidate.Unit();
        double sum = 0.0;
        for (const TVector3& p : transverse) sum += std::abs(p.Dot(n));
        best = std::max(best, sum / denom);
    }
    return std::min(1.0, std::max(0.0, best));
}

double computeThrustMinor(const std::vector<TLorentzVector>& particles)
{
    double thrust = 0.0;
    const TVector3 axis = computeThrustAxis(particles, &thrust);
    double denom = 0.0;
    std::vector<TVector3> transverse;
    for (const TLorentzVector& p4 : particles) {
        TVector3 p = p4.Vect();
        const double mag = p.Mag();
        if (mag <= 1e-9) continue;
        denom += mag;
        p -= p.Dot(axis) * axis;
        if (p.Mag() > 1e-9) transverse.push_back(p);
    }
    if (denom <= 1e-9 || transverse.empty()) return 0.0;
    double bestMajor = -1.0;
    TVector3 majorAxis(1.0, 0.0, 0.0);
    for (const TVector3& candidate : transverse) {
        const TVector3 n = candidate.Unit();
        double sum = 0.0;
        for (const TVector3& p : transverse) sum += std::abs(p.Dot(n));
        if (sum > bestMajor) {
            bestMajor = sum;
            majorAxis = n;
        }
    }
    TVector3 minorAxis = axis.Cross(majorAxis);
    if (minorAxis.Mag() <= 1e-9) return 0.0;
    minorAxis = minorAxis.Unit();
    double minor = 0.0;
    for (const TVector3& p : transverse) minor += std::abs(p.Dot(minorAxis));
    return std::min(1.0, std::max(0.0, minor / denom));
}

std::vector<double> sortedEigenvalues(const TMatrixDSym& matrix)
{
    TMatrixDSymEigen eigen(matrix);
    TVectorD values = eigen.GetEigenValues();
    std::vector<double> out(3, 0.0);
    for (int i = 0; i < 3; ++i) out[i] = values[i];
    std::sort(out.begin(), out.end(), std::greater<double>());
    return out;
}

double computeSphericity(const std::vector<TLorentzVector>& particles)
{
    TMatrixDSym tensor(3);
    tensor.Zero();
    double norm = 0.0;
    for (const TLorentzVector& p4 : particles) {
        const TVector3 p = p4.Vect();
        const double p2 = p.Mag2();
        if (p2 <= 1e-12) continue;
        norm += p2;
        tensor(0, 0) += p.X() * p.X();
        tensor(0, 1) += p.X() * p.Y();
        tensor(0, 2) += p.X() * p.Z();
        tensor(1, 1) += p.Y() * p.Y();
        tensor(1, 2) += p.Y() * p.Z();
        tensor(2, 2) += p.Z() * p.Z();
    }
    if (norm <= 1e-12) return 0.0;
    tensor(1, 0) = tensor(0, 1);
    tensor(2, 0) = tensor(0, 2);
    tensor(2, 1) = tensor(1, 2);
    tensor *= 1.0 / norm;
    const std::vector<double> ev = sortedEigenvalues(tensor);
    return std::min(1.0, std::max(0.0, 1.5 * (ev[1] + ev[2])));
}

double computeCParameter(const std::vector<TLorentzVector>& particles)
{
    TMatrixDSym tensor(3);
    tensor.Zero();
    double norm = 0.0;
    for (const TLorentzVector& p4 : particles) {
        const TVector3 p = p4.Vect();
        const double pmag = p.Mag();
        if (pmag <= 1e-12) continue;
        norm += pmag;
        tensor(0, 0) += p.X() * p.X() / pmag;
        tensor(0, 1) += p.X() * p.Y() / pmag;
        tensor(0, 2) += p.X() * p.Z() / pmag;
        tensor(1, 1) += p.Y() * p.Y() / pmag;
        tensor(1, 2) += p.Y() * p.Z() / pmag;
        tensor(2, 2) += p.Z() * p.Z() / pmag;
    }
    if (norm <= 1e-12) return 0.0;
    tensor(1, 0) = tensor(0, 1);
    tensor(2, 0) = tensor(0, 2);
    tensor(2, 1) = tensor(1, 2);
    tensor *= 1.0 / norm;
    const std::vector<double> ev = sortedEigenvalues(tensor);
    const double c = 3.0 * (ev[0] * ev[1] + ev[0] * ev[2] + ev[1] * ev[2]);
    return std::min(1.0, std::max(0.0, c));
}

void finalizeObservables(EventBranches& b)
{
    SelectionConfig config;
    b.visibleEnergy = computeVisibleEnergy(b, config);
    b.missingEnergy = std::max(0.0, b.sqrtS - b.visibleEnergy);
    const std::vector<TLorentzVector> selected = selectedVisibleParticles(b, config, true);
    b.thrust = computeThrust(selected);
    b.thrustMajor = computeThrustMajor(selected);
    b.thrustMinor = computeThrustMinor(selected);
    b.sphericity = computeSphericity(selected);
    b.cParameter = computeCParameter(selected);
}

void addParticle(EventBranches& b, int pdg, int stat, int m1, int m2, int d1, int d2,
                 const TLorentzVector& p4, bool finalState, bool parton, bool hadron,
                 bool lepton, bool photon, bool isrPhoton, double q)
{
    b.pdgId.push_back(pdg);
    b.status.push_back(stat);
    b.mother1.push_back(m1);
    b.mother2.push_back(m2);
    b.daughter1.push_back(d1);
    b.daughter2.push_back(d2);
    b.px.push_back(static_cast<float>(p4.Px()));
    b.py.push_back(static_cast<float>(p4.Py()));
    b.pz.push_back(static_cast<float>(p4.Pz()));
    b.energy.push_back(static_cast<float>(p4.E()));
    b.mass.push_back(static_cast<float>(p4.M()));
    b.charge.push_back(static_cast<float>(q));
    b.isFinal.push_back(finalState ? 1 : 0);
    b.isCharged.push_back(std::abs(q) > 0.1 ? 1 : 0);
    b.isParton.push_back(parton ? 1 : 0);
    b.isHadron.push_back(hadron ? 1 : 0);
    b.isLepton.push_back(lepton ? 1 : 0);
    b.isPhoton.push_back(photon ? 1 : 0);
    b.isISRPhoton.push_back(isrPhoton ? 1 : 0);

    if (finalState) {
        ++b.nFinalParticles;
        if (std::abs(q) > 0.1) ++b.nChargedFinalParticles;
    }
    if (parton) ++b.nPartons;
    if (isrPhoton) {
        ++b.nISRPhotons;
        b.totalISRPhotonEnergy += p4.E();
    }
}

bool hasImmediateLeptonMother(const Pythia8::Event& event, int index)
{
    const int m1 = event[index].mother1();
    const int m2 = event[index].mother2();
    if (m1 > 0 && m1 < event.size() && std::abs(event[m1].id()) == 11) return true;
    if (m2 > 0 && m2 < event.size() && std::abs(event[m2].id()) == 11) return true;
    return false;
}

void fillFromPythia(const Pythia8::Pythia& pythia, EventBranches& b)
{
    b.processId = pythia.info.code();
    b.processName = "e+e- -> gamma*/Z -> hadrons";
    b.weight = pythia.info.weight();
    for (int i = 0; i < pythia.event.size(); ++i) {
        const Pythia8::Particle& p = pythia.event[i];
        const int pdg = p.id();
        const int apdg = std::abs(pdg);
        const bool finalState = p.isFinal();
        const bool parton = apdg <= 6 || apdg == 21;
        const bool photon = pdg == 22;
        const bool isrPhoton = photon && b.isrOn && hasImmediateLeptonMother(pythia.event, i);
        TLorentzVector p4(p.px(), p.py(), p.pz(), p.e());
        addParticle(b, pdg, p.status(), p.mother1(), p.mother2(), p.daughter1(), p.daughter2(),
                    p4, finalState, parton, p.isHadron(), p.isLepton(), photon, isrPhoton, p.charge());
    }
}

int runPythia(const Options& opt)
{
    TFile output(opt.output.c_str(), "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("Could not open output file " + opt.output);
    TTree tree("Events", "Real generator event record for ISR studies");
    EventBranches b;
    b.generatorId = opt.generatorId;
    b.generatorName = opt.generatorName;
    b.isrOn = opt.isrOn;
    b.sqrtS = opt.sqrtS;
    bookTree(&tree, b);

    Pythia8::Pythia pythia("/usr/share/Pythia8/xmldoc", false);
    if (opt.vincia) {
        Pythia8::ShowerModelPtr showerModelPtr = std::make_shared<Pythia8::Vincia>();
        if (!pythia.setShowerModelPtr(showerModelPtr)) {
            throw std::runtime_error("Failed to attach Vincia shower model");
        }
        pythia.readFile("/usr/share/Pythia8/tunes/VinciaDefault.cmnd");
    }
    pythia.readString("Print:quiet = on");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Stat:showProcessLevel = off");
    pythia.readString("Stat:showErrors = off");
    pythia.readString("Random:setSeed = on");
    pythia.readString("Random:seed = " + std::to_string(opt.seed));
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("Beams:eCM = " + std::to_string(opt.sqrtS));
    pythia.readString(opt.isrOn ? "PDF:lepton = on" : "PDF:lepton = off");
    pythia.readString(opt.isrOn ? "PartonLevel:ISR = on" : "PartonLevel:ISR = off");
    pythia.readString(opt.isrOn ? "SpaceShower:QEDshowerByL = on" : "SpaceShower:QEDshowerByL = off");
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("WeakZ0:gmZmode = 0");
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");

    if (!pythia.init()) throw std::runtime_error("PYTHIA initialization failed");

    int nAccepted = 0;
    int nFailed = 0;
    for (int iev = 0; iev < opt.nEvents; ++iev) {
        if (!pythia.next()) {
            ++nFailed;
            continue;
        }
        b.eventId = nAccepted;
        reset(b);
        fillFromPythia(pythia, b);
        finalizeObservables(b);
        tree.Fill();
        ++nAccepted;
    }

    tree.Write();
    output.Close();
    std::cout << "Wrote " << nAccepted << " events to " << opt.output
              << " (" << nFailed << " generator failures)" << std::endl;
    return 0;
}

int firstOrZero(const std::vector<HepMC3::ConstGenParticlePtr>& particles)
{
    return particles.empty() ? 0 : particles.front()->id();
}

int lastOrZero(const std::vector<HepMC3::ConstGenParticlePtr>& particles)
{
    return particles.empty() ? 0 : particles.back()->id();
}

bool hepmcLooksLikeISRPhoton(const HepMC3::ConstGenParticlePtr& p)
{
    if (!p || p->pid() != 22 || p->status() != 1) return false;
    const HepMC3::FourVector& v = p->momentum();
    const double pabs = std::sqrt(v.px() * v.px() + v.py() * v.py() + v.pz() * v.pz());
    if (pabs <= 1e-9 || v.e() < 0.05) return false;
    const double absCosTheta = std::abs(v.pz() / pabs);
    if (absCosTheta > 0.990) return true;
    auto prod = p->production_vertex();
    if (!prod) return false;
    for (const auto& mother : prod->particles_in()) {
        if (mother && std::abs(mother->pid()) == 11) return true;
    }
    return false;
}

void fillFromHepMC(const HepMC3::GenEvent& event, EventBranches& b)
{
    b.processId = 1;
    b.processName = "e+e- -> hadrons";
    b.weight = event.weights().empty() ? 1.0 : event.weights().front();

    for (const auto& p : event.particles()) {
        const HepMC3::FourVector& v = p->momentum();
        TLorentzVector p4(v.px(), v.py(), v.pz(), v.e());
        const int pdg = p->pid();
        const int apdg = std::abs(pdg);
        const bool finalState = p->status() == 1 || !p->end_vertex();
        const bool parton = apdg <= 6 || apdg == 21;
        const bool photon = pdg == 22;
        const bool isrPhoton = b.isrOn && hepmcLooksLikeISRPhoton(p);

        std::vector<HepMC3::ConstGenParticlePtr> mothers;
        std::vector<HepMC3::ConstGenParticlePtr> daughters;
        if (p->production_vertex()) mothers = p->production_vertex()->particles_in();
        if (p->end_vertex()) daughters = p->end_vertex()->particles_out();
        const int m1 = firstOrZero(mothers);
        const int m2 = lastOrZero(mothers);
        const int d1 = firstOrZero(daughters);
        const int d2 = lastOrZero(daughters);

        addParticle(b, pdg, p->status(), m1, m2, d1, d2, p4, finalState, parton,
                    isHadronPdg(pdg), isLeptonPdg(pdg), photon, isrPhoton, particleCharge(pdg));
    }
}

int runHepMC(const Options& opt)
{
    if (opt.hepmcInput.empty()) throw std::runtime_error("--hepmcInput is required for --mode hepmc");
    TFile output(opt.output.c_str(), "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("Could not open output file " + opt.output);
    TTree tree("Events", "Real generator event record for ISR studies");
    EventBranches b;
    b.generatorId = opt.generatorId;
    b.generatorName = opt.generatorName;
    b.isrOn = opt.isrOn;
    b.sqrtS = opt.sqrtS;
    bookTree(&tree, b);

    std::unique_ptr<HepMC3::Reader> reader;
    if (opt.hepmcFormat == "hepmc3") {
        reader.reset(new HepMC3::ReaderAscii(opt.hepmcInput));
    } else {
        reader.reset(new HepMC3::ReaderAsciiHepMC2(opt.hepmcInput));
    }

    int nRead = 0;
    while (!reader->failed() && nRead < opt.nEvents) {
        HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
        reader->read_event(event);
        if (reader->failed()) break;
        b.eventId = nRead;
        reset(b);
        fillFromHepMC(event, b);
        finalizeObservables(b);
        tree.Fill();
        ++nRead;
    }

    tree.Write();
    output.Close();
    std::cout << "Converted " << nRead << " HepMC events to " << opt.output << std::endl;
    return 0;
}

Options parseOptions(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto requireValue = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + name);
            return argv[++i];
        };
        if (key == "--mode") opt.mode = requireValue(key);
        else if (key == "--hepmcFormat") opt.hepmcFormat = requireValue(key);
        else if (key == "--hepmcInput") opt.hepmcInput = requireValue(key);
        else if (key == "--output") opt.output = requireValue(key);
        else if (key == "--generatorName") opt.generatorName = requireValue(key);
        else if (key == "--generatorId") opt.generatorId = std::atoi(requireValue(key).c_str());
        else if (key == "--nEvents") opt.nEvents = std::atoi(requireValue(key).c_str());
        else if (key == "--seed") opt.seed = std::atoi(requireValue(key).c_str());
        else if (key == "--sqrtS") opt.sqrtS = std::atof(requireValue(key).c_str());
        else if (key == "--isrOn") opt.isrOn = std::atoi(requireValue(key).c_str()) != 0;
        else if (key == "--vincia") opt.vincia = std::atoi(requireValue(key).c_str()) != 0;
        else if (key == "--help") {
            std::cout
                << "Usage:\n"
                << "  real_isr_ntuple_producer --mode pythia --nEvents N --isrOn 0|1 --output file.root [--vincia 0|1]\n"
                << "  real_isr_ntuple_producer --mode hepmc --hepmcInput file_or_fifo --hepmcFormat hepmc2|hepmc3 --output file.root\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown option " + key);
        }
    }
    return opt;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        const Options opt = parseOptions(argc, argv);
        if (opt.mode == "pythia") return runPythia(opt);
        if (opt.mode == "hepmc") return runHepMC(opt);
        throw std::runtime_error("Unsupported mode " + opt.mode);
    } catch (const std::exception& e) {
        std::cerr << "real_isr_ntuple_producer: " << e.what() << std::endl;
        return 1;
    }
}
