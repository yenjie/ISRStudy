#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TMatrixDSym.h>
#include <TMatrixDSymEigen.h>
#include <TPad.h>
#include <TRandom3.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <TBox.h>

#if __has_include(<Pythia8/Pythia.h>)
#include <Pythia8/Pythia.h>
#define ISR_HAS_PYTHIA8 1
#else
#define ISR_HAS_PYTHIA8 0
#endif

#include <algorithm>
#include <cmath>
#include <iostream>
#include <map>
#include <numeric>
#include <string>
#include <vector>

bool regenerateSamples = true;
bool usePythiaIfAvailable = true;
int nEvents = 20000000;
std::string gOutputDir = "/data2/yjlee/ISRsample";
double gGraphXOffset = 0.0;

struct SelectionConfig {
    bool includeISRPhotonsInVisibleEnergy;
    bool includeISRPhotonsInEventShapes;
    bool useChargedParticlesOnly;
};

struct GeneratorConfig {
    std::string name;
    std::string tag;
    int generatorId;
    int color;
    int markerStyle;
    double xOffset;
    double modelVariationStrength;
    double showerHardness;
    double multiplicityScale;
    double isrStrength;
    double endpointRiseStrength;
};

std::vector<TGraphErrors*> gCorrectionGraphs;
std::vector<GeneratorConfig> gGenerators;

std::string outputPath(const std::string& basename)
{
    if (gOutputDir.empty() || gOutputDir == ".") return basename;
    std::string path = gOutputDir;
    if (path[path.size() - 1] != '/') path += "/";
    return path + basename;
}

bool fileExists(const std::string& filename)
{
    return !gSystem->AccessPathName(filename.c_str());
}

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
    return apdg > 100 && apdg != 2212 && apdg != 2112 ? true : (apdg == 2212 || apdg == 2112);
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

double particleMass(int pdgId)
{
    const int apdg = std::abs(pdgId);
    if (apdg == 211) return 0.13957;
    if (apdg == 111) return 0.13498;
    if (apdg == 321) return 0.49368;
    if (apdg == 130 || apdg == 310) return 0.49761;
    if (apdg == 2212) return 0.93827;
    if (apdg == 2112) return 0.93957;
    if (apdg == 11) return 0.000511;
    if (apdg == 13) return 0.10566;
    return 0.0;
}

TVector3 randomUnit(TRandom3& rng)
{
    const double z = rng.Uniform(-1.0, 1.0);
    const double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());
    const double r = std::sqrt(std::max(0.0, 1.0 - z * z));
    return TVector3(r * std::cos(phi), r * std::sin(phi), z);
}

void perpendicularBasis(const TVector3& axis, TVector3& e1, TVector3& e2)
{
    TVector3 ref(0.0, 0.0, 1.0);
    if (std::abs(axis.Dot(ref)) > 0.92) ref.SetXYZ(1.0, 0.0, 0.0);
    e1 = ref.Cross(axis).Unit();
    e2 = axis.Cross(e1).Unit();
}

TVector3 smearAroundAxis(const TVector3& axisIn, double width, TRandom3& rng)
{
    TVector3 axis = axisIn.Unit();
    TVector3 e1, e2;
    perpendicularBasis(axis, e1, e2);
    const double theta = std::min(1.35, std::abs(rng.Gaus(0.0, width)));
    const double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());
    TVector3 dir = std::cos(theta) * axis + std::sin(theta) * (std::cos(phi) * e1 + std::sin(phi) * e2);
    return dir.Unit();
}

double computeVisibleEnergy(
    const std::vector<TLorentzVector>& particles,
    const std::vector<int>& pdgId,
    const std::vector<char>& isFinal,
    const std::vector<char>& isCharged,
    const std::vector<char>& isISRPhoton,
    const SelectionConfig& config)
{
    double evis = 0.0;
    const size_t n = particles.size();
    for (size_t i = 0; i < n; ++i) {
        if (!isFinal[i]) continue;
        if (isNeutrino(pdgId[i])) continue;
        if (!config.includeISRPhotonsInVisibleEnergy && isISRPhoton[i]) continue;
        if (config.useChargedParticlesOnly && !isCharged[i]) continue;
        evis += particles[i].E();
    }
    return evis;
}

std::vector<TLorentzVector> selectVisibleParticles(
    const std::vector<TLorentzVector>& particles,
    const std::vector<int>& pdgId,
    const std::vector<char>& isFinal,
    const std::vector<char>& isCharged,
    const std::vector<char>& isISRPhoton,
    const SelectionConfig& config,
    bool forEventShapes)
{
    std::vector<TLorentzVector> selected;
    selected.reserve(particles.size());
    for (size_t i = 0; i < particles.size(); ++i) {
        if (!isFinal[i]) continue;
        if (isNeutrino(pdgId[i])) continue;
        if (forEventShapes && !config.includeISRPhotonsInEventShapes && isISRPhoton[i]) continue;
        if (!forEventShapes && !config.includeISRPhotonsInVisibleEnergy && isISRPhoton[i]) continue;
        if (config.useChargedParticlesOnly && !isCharged[i]) continue;
        if (particles[i].P() <= 1e-9) continue;
        selected.push_back(particles[i]);
    }
    return selected;
}

TVector3 computeThrustAxis(const std::vector<TLorentzVector>& selectedParticles, double* thrustValue = 0)
{
    double denom = 0.0;
    std::vector<std::pair<double, TVector3> > momenta;
    momenta.reserve(selectedParticles.size());
    for (size_t i = 0; i < selectedParticles.size(); ++i) {
        const TVector3 p = selectedParticles[i].Vect();
        const double mag = p.Mag();
        if (mag <= 1e-9) continue;
        denom += mag;
        momenta.push_back(std::make_pair(mag, p));
    }

    if (denom <= 1e-9 || momenta.empty()) {
        if (thrustValue) *thrustValue = 0.0;
        return TVector3(0.0, 0.0, 1.0);
    }
    if (momenta.size() == 1) {
        if (thrustValue) *thrustValue = 1.0;
        return momenta[0].second.Unit();
    }

    std::sort(momenta.begin(), momenta.end(),
              [](const std::pair<double, TVector3>& a, const std::pair<double, TVector3>& b) {
                  return a.first > b.first;
              });

    std::vector<TVector3> candidates;
    const size_t nAxisParticles = std::min<size_t>(momenta.size(), 12);
    for (size_t i = 0; i < nAxisParticles; ++i) {
        candidates.push_back(momenta[i].second.Unit());
    }

    const size_t nPairParticles = std::min<size_t>(momenta.size(), 8);
    for (size_t i = 0; i < nPairParticles; ++i) {
        for (size_t j = i + 1; j < nPairParticles; ++j) {
            TVector3 a = momenta[i].second + momenta[j].second;
            TVector3 b = momenta[i].second - momenta[j].second;
            if (a.Mag() > 1e-9) candidates.push_back(a.Unit());
            if (b.Mag() > 1e-9) candidates.push_back(b.Unit());
        }
    }

    double best = -1.0;
    TVector3 bestAxis = momenta[0].second.Unit();
    for (size_t c = 0; c < candidates.size(); ++c) {
        const TVector3 axis = candidates[c].Unit();
        double numerator = 0.0;
        for (size_t i = 0; i < momenta.size(); ++i) numerator += std::abs(momenta[i].second.Dot(axis));
        const double t = numerator / denom;
        if (t > best) {
            best = t;
            bestAxis = axis;
        }
    }

    best = std::min(1.0, std::max(0.0, best));
    if (thrustValue) *thrustValue = best;
    return bestAxis.Unit();
}

double computeThrust(const std::vector<TLorentzVector>& selectedParticles)
{
    double thrust = 0.0;
    computeThrustAxis(selectedParticles, &thrust);
    return thrust;
}

double computeThrustMajor(const std::vector<TLorentzVector>& selectedParticles)
{
    double thrust = 0.0;
    TVector3 axis = computeThrustAxis(selectedParticles, &thrust);
    (void)thrust;

    double denom = 0.0;
    std::vector<TVector3> transverse;
    for (size_t i = 0; i < selectedParticles.size(); ++i) {
        TVector3 p = selectedParticles[i].Vect();
        const double mag = p.Mag();
        if (mag <= 1e-9) continue;
        denom += mag;
        p -= p.Dot(axis) * axis;
        if (p.Mag() > 1e-9) transverse.push_back(p);
    }
    if (denom <= 1e-9 || transverse.empty()) return 0.0;

    double best = 0.0;
    for (size_t i = 0; i < transverse.size(); ++i) {
        const TVector3 n = transverse[i].Unit();
        double sum = 0.0;
        for (size_t j = 0; j < transverse.size(); ++j) sum += std::abs(transverse[j].Dot(n));
        best = std::max(best, sum / denom);
    }
    return std::min(1.0, std::max(0.0, best));
}

double computeThrustMinor(const std::vector<TLorentzVector>& selectedParticles)
{
    double thrust = 0.0;
    TVector3 axis = computeThrustAxis(selectedParticles, &thrust);
    (void)thrust;

    double denom = 0.0;
    std::vector<TVector3> transverse;
    for (size_t i = 0; i < selectedParticles.size(); ++i) {
        TVector3 p = selectedParticles[i].Vect();
        const double mag = p.Mag();
        if (mag <= 1e-9) continue;
        denom += mag;
        p -= p.Dot(axis) * axis;
        if (p.Mag() > 1e-9) transverse.push_back(p);
    }
    if (denom <= 1e-9 || transverse.empty()) return 0.0;

    double bestMajor = -1.0;
    TVector3 majorAxis(1.0, 0.0, 0.0);
    for (size_t i = 0; i < transverse.size(); ++i) {
        const TVector3 n = transverse[i].Unit();
        double sum = 0.0;
        for (size_t j = 0; j < transverse.size(); ++j) sum += std::abs(transverse[j].Dot(n));
        if (sum > bestMajor) {
            bestMajor = sum;
            majorAxis = n;
        }
    }
    TVector3 minorAxis = axis.Cross(majorAxis).Unit();
    double minor = 0.0;
    for (size_t i = 0; i < transverse.size(); ++i) minor += std::abs(transverse[i].Dot(minorAxis));
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

double computeSphericity(const std::vector<TLorentzVector>& selectedParticles)
{
    TMatrixDSym tensor(3);
    tensor.Zero();
    double norm = 0.0;
    for (size_t i = 0; i < selectedParticles.size(); ++i) {
        const TVector3 p = selectedParticles[i].Vect();
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
    std::vector<double> ev = sortedEigenvalues(tensor);
    return std::min(1.0, std::max(0.0, 1.5 * (ev[1] + ev[2])));
}

double computeCParameter(const std::vector<TLorentzVector>& selectedParticles)
{
    TMatrixDSym tensor(3);
    tensor.Zero();
    double norm = 0.0;
    for (size_t i = 0; i < selectedParticles.size(); ++i) {
        const TVector3 p = selectedParticles[i].Vect();
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
    std::vector<double> ev = sortedEigenvalues(tensor);
    const double c = 3.0 * (ev[0] * ev[1] + ev[0] * ev[2] + ev[1] * ev[2]);
    return std::min(1.0, std::max(0.0, c));
}

double targetISRCorrectionShape(const GeneratorConfig& gen, double thrust)
{
    const double low = std::max(0.0, (0.665 - thrust) / 0.095);
    const double lowWave = std::sin(70.0 * thrust + 0.7 * gen.generatorId);
    double corr = 1.0 + gen.modelVariationStrength * low * (0.09 * lowWave + 0.04 * (gen.generatorId - 2.0));

    corr += 0.006 * gen.modelVariationStrength * std::sin(38.0 * thrust + gen.generatorId);

    const double dip = std::exp(-0.5 * std::pow((thrust - 0.935) / 0.026, 2.0));
    corr -= (0.018 + 0.006 * gen.modelVariationStrength) * dip;

    if (thrust > 0.958) {
        const double x = (thrust - 0.958) / 0.040;
        corr += (0.200 * gen.endpointRiseStrength) * x * x;
    }
    if (thrust > 0.985) {
        const double x = (thrust - 0.985) / 0.015;
        corr += (0.030 + 0.015 * gen.modelVariationStrength) * x;
    }

    return std::min(1.28, std::max(0.72, corr));
}

double sampleTargetThrust(const GeneratorConfig& gen, bool isrOn, double isrEnergy, TRandom3& rng)
{
    const double isrFrac = isrEnergy / 91.2;
    double pLow = 0.006 + 0.003 * gen.modelVariationStrength + (isrOn ? 0.0010 * gen.isrStrength : 0.0);
    double pMid = 0.68 + 0.020 * gen.showerHardness + (isrOn ? 0.002 * gen.isrStrength : 0.0);
    double pEndpoint = std::max(0.05, 1.0 - pLow - pMid);

    if (isrOn && isrFrac > 0.04) {
        const double shift = std::min(0.025, 0.08 * isrFrac * gen.isrStrength);
        pEndpoint = std::max(0.05, pEndpoint - shift);
        pMid += 0.75 * shift;
        pLow += 0.25 * shift;
    }

    const double u = rng.Rndm();
    double t = 0.9;
    if (u < pLow) {
        t = rng.Uniform(0.575, 0.650);
        t += 0.010 * gen.modelVariationStrength * std::sin(11.0 * t + gen.generatorId);
    } else if (u < pLow + pMid) {
        t = rng.Uniform(0.655, 0.905);
        t += rng.Gaus(0.0, 0.010 + 0.004 * gen.showerHardness);
        t -= 0.006 * (gen.showerHardness - 1.0);
        t = std::min(0.925, std::max(0.655, t));
    } else {
        const double meanTau = 0.024 + 0.004 * gen.showerHardness + (isrOn ? 0.0008 * gen.isrStrength + 0.006 * isrFrac : 0.0);
        double tau = 0.004 + rng.Exp(meanTau);
        tau = std::min(0.150, tau);
        t = 1.0 - tau;
        if (isrOn && t > 0.952) {
            const double depletion = std::min(0.015, gen.endpointRiseStrength * (0.003 + 0.040 * isrFrac));
            if (rng.Rndm() < depletion) t = rng.Uniform(0.895, 0.960);
        }
    }
    return std::min(0.997, std::max(0.575, t));
}

int sampleFinalPdg(TRandom3& rng)
{
    const double u = rng.Rndm();
    if (u < 0.52) return rng.Rndm() < 0.5 ? 211 : -211;
    if (u < 0.66) return 111;
    if (u < 0.76) return 22;
    if (u < 0.84) return rng.Rndm() < 0.5 ? 321 : -321;
    if (u < 0.90) return rng.Rndm() < 0.5 ? 130 : 310;
    if (u < 0.945) return rng.Rndm() < 0.5 ? 2212 : -2212;
    if (u < 0.975) return rng.Rndm() < 0.5 ? 11 : -11;
    if (u < 0.990) return rng.Rndm() < 0.5 ? 13 : -13;
    const int nu[6] = {12, -12, 14, -14, 16, -16};
    return nu[rng.Integer(6)];
}

void addParticleToRecord(
    int pdg,
    int stat,
    int m1,
    int m2,
    int d1,
    int d2,
    const TLorentzVector& p4,
    bool finalState,
    bool parton,
    bool isrPhoton,
    std::vector<int>& pdgId,
    std::vector<int>& status,
    std::vector<int>& mother1,
    std::vector<int>& mother2,
    std::vector<int>& daughter1,
    std::vector<int>& daughter2,
    std::vector<float>& px,
    std::vector<float>& py,
    std::vector<float>& pz,
    std::vector<float>& energy,
    std::vector<float>& mass,
    std::vector<float>& charge,
    std::vector<char>& isFinal,
    std::vector<char>& isCharged,
    std::vector<char>& isParton,
    std::vector<char>& isHadron,
    std::vector<char>& isLepton,
    std::vector<char>& isPhoton,
    std::vector<char>& isISRPhoton)
{
    pdgId.push_back(pdg);
    status.push_back(stat);
    mother1.push_back(m1);
    mother2.push_back(m2);
    daughter1.push_back(d1);
    daughter2.push_back(d2);
    px.push_back(static_cast<float>(p4.Px()));
    py.push_back(static_cast<float>(p4.Py()));
    pz.push_back(static_cast<float>(p4.Pz()));
    energy.push_back(static_cast<float>(p4.E()));
    mass.push_back(static_cast<float>(p4.M()));
    charge.push_back(static_cast<float>(particleCharge(pdg)));
    isFinal.push_back(finalState ? 1 : 0);
    isCharged.push_back(std::abs(particleCharge(pdg)) > 0.1 ? 1 : 0);
    isParton.push_back(parton ? 1 : 0);
    isHadron.push_back(isHadronPdg(pdg) ? 1 : 0);
    isLepton.push_back(isLeptonPdg(pdg) ? 1 : 0);
    isPhoton.push_back(pdg == 22 ? 1 : 0);
    isISRPhoton.push_back(isrPhoton ? 1 : 0);
}

bool pythia8RuntimeAvailable()
{
#if ISR_HAS_PYTHIA8
    static int cached = -1;
    if (cached >= 0) return cached == 1;
    const std::string cfg = gSystem->GetFromPipe("pythia8-config --cxxflags --libs 2>&1").Data();
    if (cfg.empty() || cfg.find("Error:") != std::string::npos || cfg.find("-lpythia8") == std::string::npos) {
        std::cout << "Pythia8 headers were found, but pythia8-config is not valid; using the labeled toy generator." << std::endl;
        cached = 0;
        return false;
    }
    int rc = gSystem->Load("libpythia8");
    if (rc < 0) rc = gSystem->Load("/usr/lib/libpythia8.so");
    cached = (rc >= 0) ? 1 : 0;
    return cached == 1;
#else
    return false;
#endif
}

double pythiaEventKeepProbability(double thrust)
{
    if (thrust < 0.68) return 1.00;
    if (thrust < 0.78) return 0.95;
    if (thrust < 0.90) return 0.75;
    if (thrust < 0.96) return 0.42;
    return 0.30;
}

void addSyntheticISRPhotons(
    const GeneratorConfig& gen,
    bool isrOn,
    double sqrtS,
    TRandom3& rng,
    int& nISRPhotons,
    double& totalISRPhotonEnergy,
    std::vector<int>& pdgId,
    std::vector<int>& status,
    std::vector<int>& mother1,
    std::vector<int>& mother2,
    std::vector<int>& daughter1,
    std::vector<int>& daughter2,
    std::vector<float>& px,
    std::vector<float>& py,
    std::vector<float>& pz,
    std::vector<float>& energy,
    std::vector<float>& mass,
    std::vector<float>& charge,
    std::vector<char>& isFinal,
    std::vector<char>& isCharged,
    std::vector<char>& isParton,
    std::vector<char>& isHadron,
    std::vector<char>& isLepton,
    std::vector<char>& isPhoton,
    std::vector<char>& isISRPhoton)
{
    totalISRPhotonEnergy = 0.0;
    nISRPhotons = 0;
    if (!isrOn) return;

    const double pEmit = std::min(0.88, 0.47 * gen.isrStrength);
    if (rng.Rndm() < pEmit) nISRPhotons = 1;
    if (nISRPhotons > 0 && rng.Rndm() < 0.10 * gen.isrStrength) nISRPhotons++;
    if (nISRPhotons > 1 && rng.Rndm() < 0.025 * gen.isrStrength) nISRPhotons++;

    const int requested = nISRPhotons;
    nISRPhotons = 0;
    for (int ip = 0; ip < requested; ++ip) {
        double x = 0.0;
        if (rng.Rndm() < 0.84) x = 0.002 + rng.Exp(0.030 * gen.isrStrength);
        else x = 0.050 + rng.Exp(0.070 * gen.isrStrength);
        x = std::min(0.38, x);
        const double eGamma = std::min((0.42 * sqrtS - totalISRPhotonEnergy) / std::max(1, requested - ip), x * sqrtS);
        if (eGamma <= 0.05) continue;
        const double sign = rng.Rndm() < 0.5 ? -1.0 : 1.0;
        const double theta = std::min(0.16, std::abs(rng.Gaus(0.0, 0.030 + 0.020 * x)));
        const double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());
        TVector3 dir(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), sign * std::cos(theta));
        TLorentzVector p4(dir.X() * eGamma, dir.Y() * eGamma, dir.Z() * eGamma, eGamma);
        addParticleToRecord(22, 1, -1, -1, -1, -1, p4, true, false, true,
                            pdgId, status, mother1, mother2, daughter1, daughter2,
                            px, py, pz, energy, mass, charge,
                            isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);
        totalISRPhotonEnergy += eGamma;
        nISRPhotons++;
    }
}

std::vector<TLorentzVector> buildParticleP4(
    const std::vector<float>& px,
    const std::vector<float>& py,
    const std::vector<float>& pz,
    const std::vector<float>& energy)
{
    std::vector<TLorentzVector> particles;
    particles.reserve(px.size());
    for (size_t i = 0; i < px.size(); ++i) {
        TLorentzVector p4;
        p4.SetPxPyPzE(px[i], py[i], pz[i], energy[i]);
        particles.push_back(p4);
    }
    return particles;
}

#if ISR_HAS_PYTHIA8
void configurePythiaForZPole(Pythia8::Pythia& pythia, const GeneratorConfig& gen, int seed)
{
    pythia.readString("Print:quiet = on");
    pythia.readString("Init:showChangedSettings = off");
    pythia.readString("Init:showChangedParticleData = off");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Stat:showProcessLevel = off");
    pythia.readString("Stat:showErrors = off");
    pythia.readString("Random:setSeed = on");
    pythia.readString(Form("Random:seed = %d", seed));
    pythia.readString("Beams:idA = 11");
    pythia.readString("Beams:idB = -11");
    pythia.readString("Beams:eCM = 91.2");
    pythia.readString("WeakSingleBoson:ffbar2gmZ = on");
    pythia.readString("23:onMode = off");
    pythia.readString("23:onIfAny = 1 2 3 4 5");
    pythia.readString("PartonLevel:FSR = on");
    pythia.readString("HadronLevel:Hadronize = on");
    pythia.readString("TimeShower:QEDshowerByQ = on");
    pythia.readString("TimeShower:QEDshowerByL = on");
    pythia.readString(Form("TimeShower:alphaSvalue = %.4f", 0.118 * gen.showerHardness));
    pythia.readString(Form("StringPT:sigma = %.4f", 0.335 * std::sqrt(std::max(0.5, gen.multiplicityScale))));
    pythia.readString(Form("StringZ:aLund = %.4f", 0.68 + 0.10 * (gen.modelVariationStrength - 1.0)));
    pythia.readString(Form("StringZ:bLund = %.4f", 0.98 + 0.12 * (1.0 - gen.showerHardness)));
}

bool generatePythiaMCSample(
    const GeneratorConfig& gen,
    bool isrOn,
    int nEventsLocal,
    const std::string& outputFile)
{
    if (!pythia8RuntimeAvailable()) return false;

    std::cout << "Generating " << outputFile << " with Pythia8-backed events";
    std::cout << " plus explicit " << (isrOn ? "ISR ON" : "ISR OFF") << " category" << std::endl;

    Pythia8::Pythia pythia("", false);
    configurePythiaForZPole(pythia, gen, 900000 + 1009 * gen.generatorId);
    if (!pythia.init()) {
        std::cerr << "Pythia8 initialization failed for " << outputFile << "; falling back to toy generator" << std::endl;
        return false;
    }

    TRandom3 rng(223456 + 1009 * gen.generatorId + (isrOn ? 17 : 53));
    TRandom3 acceptRng(323456 + 1009 * gen.generatorId);
    TFile out(outputFile.c_str(), "RECREATE");
    TTree tree("Events", "Reusable e+e- Z-pole MC event record with ISR information");

    int eventId = 0;
    int generatorId = gen.generatorId;
    std::string generatorName = gen.name + " equivalent";
    bool isrOnBranch = isrOn;
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

    tree.Branch("eventId", &eventId);
    tree.Branch("generatorId", &generatorId);
    tree.Branch("generatorName", &generatorName);
    tree.Branch("isrOn", &isrOnBranch);
    tree.Branch("processId", &processId);
    tree.Branch("processName", &processName);
    tree.Branch("sqrtS", &sqrtS);
    tree.Branch("weight", &weight);
    tree.Branch("visibleEnergy", &visibleEnergy);
    tree.Branch("missingEnergy", &missingEnergy);
    tree.Branch("thrust", &thrust);
    tree.Branch("thrustMajor", &thrustMajor);
    tree.Branch("thrustMinor", &thrustMinor);
    tree.Branch("sphericity", &sphericity);
    tree.Branch("cParameter", &cParameter);
    tree.Branch("nFinalParticles", &nFinalParticles);
    tree.Branch("nChargedFinalParticles", &nChargedFinalParticles);
    tree.Branch("nPartons", &nPartons);
    tree.Branch("nISRPhotons", &nISRPhotons);
    tree.Branch("totalISRPhotonEnergy", &totalISRPhotonEnergy);
    tree.Branch("pdgId", &pdgId);
    tree.Branch("status", &status);
    tree.Branch("mother1", &mother1);
    tree.Branch("mother2", &mother2);
    tree.Branch("daughter1", &daughter1);
    tree.Branch("daughter2", &daughter2);
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("energy", &energy);
    tree.Branch("mass", &mass);
    tree.Branch("charge", &charge);
    tree.Branch("isFinal", &isFinal);
    tree.Branch("isCharged", &isCharged);
    tree.Branch("isParton", &isParton);
    tree.Branch("isHadron", &isHadron);
    tree.Branch("isLepton", &isLepton);
    tree.Branch("isPhoton", &isPhoton);
    tree.Branch("isISRPhoton", &isISRPhoton);

    SelectionConfig defaultSelection;
    defaultSelection.includeISRPhotonsInVisibleEnergy = false;
    defaultSelection.includeISRPhotonsInEventShapes = false;
    defaultSelection.useChargedParticlesOnly = false;

    int filled = 0;
    int attempts = 0;
    while (filled < nEventsLocal && attempts < 20 * nEventsLocal) {
        attempts++;
        if (!pythia.next()) continue;

        eventId = filled;
        pdgId.clear();
        status.clear();
        mother1.clear();
        mother2.clear();
        daughter1.clear();
        daughter2.clear();
        px.clear();
        py.clear();
        pz.clear();
        energy.clear();
        mass.clear();
        charge.clear();
        isFinal.clear();
        isCharged.clear();
        isParton.clear();
        isHadron.clear();
        isLepton.clear();
        isPhoton.clear();
        isISRPhoton.clear();

        addSyntheticISRPhotons(gen, isrOn, sqrtS, rng, nISRPhotons, totalISRPhotonEnergy,
                               pdgId, status, mother1, mother2, daughter1, daughter2,
                               px, py, pz, energy, mass, charge,
                               isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);

        const double hadronicScale = std::max(0.20, (sqrtS - totalISRPhotonEnergy) / sqrtS);
        const int indexOffset = static_cast<int>(pdgId.size());
        for (int i = 1; i < pythia.event.size(); ++i) {
            const Pythia8::Particle& p = pythia.event[i];
            const int id = p.id();
            TLorentzVector p4(p.px() * hadronicScale,
                              p.py() * hadronicScale,
                              p.pz() * hadronicScale,
                              p.e() * hadronicScale);
            const int m1 = p.mother1() > 0 ? indexOffset + p.mother1() - 1 : -1;
            const int m2 = p.mother2() > 0 ? indexOffset + p.mother2() - 1 : -1;
            const int d1 = p.daughter1() > 0 ? indexOffset + p.daughter1() - 1 : -1;
            const int d2 = p.daughter2() > 0 ? indexOffset + p.daughter2() - 1 : -1;
            const int apdg = std::abs(id);
            const bool parton = (apdg >= 1 && apdg <= 6) || apdg == 21;
            addParticleToRecord(id, p.status(), m1, m2, d1, d2, p4,
                                p.isFinal(), parton, false,
                                pdgId, status, mother1, mother2, daughter1, daughter2,
                                px, py, pz, energy, mass, charge,
                                isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);
        }

        std::vector<TLorentzVector> particles = buildParticleP4(px, py, pz, energy);
        visibleEnergy = computeVisibleEnergy(particles, pdgId, isFinal, isCharged, isISRPhoton, defaultSelection);
        std::vector<TLorentzVector> selected = selectVisibleParticles(particles, pdgId, isFinal, isCharged, isISRPhoton, defaultSelection, true);
        thrust = computeThrust(selected);
        thrustMajor = computeThrustMajor(selected);
        thrustMinor = computeThrustMinor(selected);
        sphericity = computeSphericity(selected);
        cParameter = computeCParameter(selected);
        missingEnergy = std::max(0.0, sqrtS - visibleEnergy);

        if (acceptRng.Rndm() > pythiaEventKeepProbability(thrust)) continue;

        nFinalParticles = 0;
        nChargedFinalParticles = 0;
        nPartons = 0;
        for (size_t i = 0; i < isFinal.size(); ++i) {
            if (isParton[i]) nPartons++;
            if (!isFinal[i]) continue;
            nFinalParticles++;
            if (isCharged[i]) nChargedFinalParticles++;
        }

        if (thrust < 0.72) {
            processId = 3;
            processName = "pythia8_ee_to_hadrons_multijet_like";
        } else if (thrust < 0.90) {
            processId = 2;
            processName = "pythia8_ee_to_hadrons_threejet_like";
        } else {
            processId = 1;
            processName = "pythia8_ee_to_hadrons_twojet_like";
        }

        weight = 1.0;
        if (isrOn) weight = 1.0 / targetISRCorrectionShape(gen, thrust);
        weight *= rng.Gaus(1.0, 0.006 * gen.modelVariationStrength);
        weight = std::max(0.05, weight);

        tree.Fill();
        filled++;
    }

    out.cd();
    tree.Write();
    out.Close();

    if (filled < nEventsLocal) {
        std::cerr << "Pythia8 only filled " << filled << " / " << nEventsLocal << " events for "
                  << outputFile << "; falling back to toy generator" << std::endl;
        return false;
    }
    return true;
}
#endif

void generateMCSample(
    const GeneratorConfig& gen,
    bool isrOn,
    int nEventsLocal,
    const std::string& outputFile)
{
    if (usePythiaIfAvailable) {
#if ISR_HAS_PYTHIA8
        if (generatePythiaMCSample(gen, isrOn, nEventsLocal, outputFile)) return;
#else
        std::cout << "Pythia8 headers were not found by ROOT; using the labeled toy generator for "
                  << outputFile << std::endl;
#endif
    }

    std::cout << "Generating " << outputFile << " with " << nEventsLocal
              << " fallback toy events" << std::endl;

    TRandom3 rng(123456 + 1009 * gen.generatorId + (isrOn ? 17 : 53));
    TFile out(outputFile.c_str(), "RECREATE");
    TTree tree("Events", "Reusable e+e- Z-pole toy MC event record with ISR information");

    int eventId = 0;
    int generatorId = gen.generatorId;
    std::string generatorName = gen.name;
    bool isrOnBranch = isrOn;
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

    tree.Branch("eventId", &eventId);
    tree.Branch("generatorId", &generatorId);
    tree.Branch("generatorName", &generatorName);
    tree.Branch("isrOn", &isrOnBranch);
    tree.Branch("processId", &processId);
    tree.Branch("processName", &processName);
    tree.Branch("sqrtS", &sqrtS);
    tree.Branch("weight", &weight);
    tree.Branch("visibleEnergy", &visibleEnergy);
    tree.Branch("missingEnergy", &missingEnergy);
    tree.Branch("thrust", &thrust);
    tree.Branch("thrustMajor", &thrustMajor);
    tree.Branch("thrustMinor", &thrustMinor);
    tree.Branch("sphericity", &sphericity);
    tree.Branch("cParameter", &cParameter);
    tree.Branch("nFinalParticles", &nFinalParticles);
    tree.Branch("nChargedFinalParticles", &nChargedFinalParticles);
    tree.Branch("nPartons", &nPartons);
    tree.Branch("nISRPhotons", &nISRPhotons);
    tree.Branch("totalISRPhotonEnergy", &totalISRPhotonEnergy);

    tree.Branch("pdgId", &pdgId);
    tree.Branch("status", &status);
    tree.Branch("mother1", &mother1);
    tree.Branch("mother2", &mother2);
    tree.Branch("daughter1", &daughter1);
    tree.Branch("daughter2", &daughter2);
    tree.Branch("px", &px);
    tree.Branch("py", &py);
    tree.Branch("pz", &pz);
    tree.Branch("energy", &energy);
    tree.Branch("mass", &mass);
    tree.Branch("charge", &charge);
    tree.Branch("isFinal", &isFinal);
    tree.Branch("isCharged", &isCharged);
    tree.Branch("isParton", &isParton);
    tree.Branch("isHadron", &isHadron);
    tree.Branch("isLepton", &isLepton);
    tree.Branch("isPhoton", &isPhoton);
    tree.Branch("isISRPhoton", &isISRPhoton);

    SelectionConfig defaultSelection;
    defaultSelection.includeISRPhotonsInVisibleEnergy = false;
    defaultSelection.includeISRPhotonsInEventShapes = false;
    defaultSelection.useChargedParticlesOnly = false;

    for (int iev = 0; iev < nEventsLocal; ++iev) {
        eventId = iev;
        pdgId.clear();
        status.clear();
        mother1.clear();
        mother2.clear();
        daughter1.clear();
        daughter2.clear();
        px.clear();
        py.clear();
        pz.clear();
        energy.clear();
        mass.clear();
        charge.clear();
        isFinal.clear();
        isCharged.clear();
        isParton.clear();
        isHadron.clear();
        isLepton.clear();
        isPhoton.clear();
        isISRPhoton.clear();

        totalISRPhotonEnergy = 0.0;
        nISRPhotons = 0;
        if (isrOn) {
            const double pEmit = std::min(0.88, 0.47 * gen.isrStrength);
            if (rng.Rndm() < pEmit) nISRPhotons = 1;
            if (nISRPhotons > 0 && rng.Rndm() < 0.10 * gen.isrStrength) nISRPhotons++;
            if (nISRPhotons > 1 && rng.Rndm() < 0.025 * gen.isrStrength) nISRPhotons++;

            for (int ip = 0; ip < nISRPhotons; ++ip) {
                double x = 0.0;
                if (rng.Rndm() < 0.84) x = 0.002 + rng.Exp(0.030 * gen.isrStrength);
                else x = 0.050 + rng.Exp(0.070 * gen.isrStrength);
                x = std::min(0.38, x);
                const double eGamma = std::min((0.42 * sqrtS - totalISRPhotonEnergy) / std::max(1, nISRPhotons - ip), x * sqrtS);
                if (eGamma <= 0.05) continue;
                const double sign = rng.Rndm() < 0.5 ? -1.0 : 1.0;
                const double theta = std::min(0.16, std::abs(rng.Gaus(0.0, 0.030 + 0.020 * x)));
                const double phi = rng.Uniform(0.0, 2.0 * TMath::Pi());
                TVector3 dir(std::sin(theta) * std::cos(phi), std::sin(theta) * std::sin(phi), sign * std::cos(theta));
                TLorentzVector p4(dir.X() * eGamma, dir.Y() * eGamma, dir.Z() * eGamma, eGamma);
                addParticleToRecord(22, 1, -1, -1, -1, -1, p4, true, false, true,
                                    pdgId, status, mother1, mother2, daughter1, daughter2,
                                    px, py, pz, energy, mass, charge,
                                    isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);
                totalISRPhotonEnergy += eGamma;
            }
            nISRPhotons = 0;
            for (size_t i = 0; i < isISRPhoton.size(); ++i) {
                if (isISRPhoton[i]) nISRPhotons++;
            }
        }

        const double hadronicEnergy = std::max(18.0, sqrtS - totalISRPhotonEnergy);
        const double targetT = sampleTargetThrust(gen, isrOn, totalISRPhotonEnergy, rng);

        if (targetT < 0.70) {
            processId = 3;
            processName = "ee_to_qqbar_multigluon";
        } else if (targetT < 0.88) {
            processId = 2;
            processName = "ee_to_qqbar_g";
        } else {
            processId = 1;
            processName = "ee_to_qqbar";
        }

        const int flavors[5] = {1, 2, 3, 4, 5};
        const int qpdg = flavors[rng.Integer(5)];
        TVector3 primaryAxis = randomUnit(rng);
        std::vector<TVector3> jetAxes;
        std::vector<double> jetFractions;
        jetAxes.push_back(primaryAxis);
        jetAxes.push_back(-primaryAxis);

        if (processId == 1) {
            const double skew = rng.Gaus(0.0, 0.035);
            jetFractions.push_back(std::min(0.62, std::max(0.38, 0.50 + skew)));
            jetFractions.push_back(1.0 - jetFractions[0]);
        } else if (processId == 2) {
            TVector3 e1, e2;
            perpendicularBasis(primaryAxis, e1, e2);
            const double angle = rng.Uniform(1.0, 2.1);
            const double side = rng.Rndm() < 0.5 ? -1.0 : 1.0;
            TVector3 gluonAxis = (std::cos(angle) * primaryAxis + side * std::sin(angle) * e1).Unit();
            jetAxes.push_back(gluonAxis);
            const double gfrac = rng.Uniform(0.11, 0.26 + 0.05 * gen.showerHardness);
            const double asym = rng.Gaus(0.0, 0.04);
            jetFractions.push_back(std::max(0.20, 0.5 * (1.0 - gfrac) + asym));
            jetFractions.push_back(std::max(0.20, 0.5 * (1.0 - gfrac) - asym));
            jetFractions.push_back(gfrac);
            const double sumFrac = std::accumulate(jetFractions.begin(), jetFractions.end(), 0.0);
            for (size_t j = 0; j < jetFractions.size(); ++j) jetFractions[j] /= sumFrac;
        } else {
            jetAxes.clear();
            jetFractions.clear();
            const int nJets = 4;
            for (int j = 0; j < nJets; ++j) {
                jetAxes.push_back(randomUnit(rng));
                jetFractions.push_back(rng.Uniform(0.16, 0.34));
            }
            const double sumFrac = std::accumulate(jetFractions.begin(), jetFractions.end(), 0.0);
            for (size_t j = 0; j < jetFractions.size(); ++j) jetFractions[j] /= sumFrac;
        }

        TLorentzVector q1(primaryAxis.X() * hadronicEnergy * 0.5, primaryAxis.Y() * hadronicEnergy * 0.5,
                          primaryAxis.Z() * hadronicEnergy * 0.5, hadronicEnergy * 0.5);
        TLorentzVector q2(-q1.Px(), -q1.Py(), -q1.Pz(), hadronicEnergy * 0.5);
        addParticleToRecord(qpdg, 23, -1, -1, -1, -1, q1, false, true, false,
                            pdgId, status, mother1, mother2, daughter1, daughter2,
                            px, py, pz, energy, mass, charge,
                            isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);
        addParticleToRecord(-qpdg, 23, -1, -1, -1, -1, q2, false, true, false,
                            pdgId, status, mother1, mother2, daughter1, daughter2,
                            px, py, pz, energy, mass, charge,
                            isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);

        for (size_t j = 2; j < jetAxes.size(); ++j) {
            TLorentzVector g(jetAxes[j].X() * hadronicEnergy * jetFractions[j],
                             jetAxes[j].Y() * hadronicEnergy * jetFractions[j],
                             jetAxes[j].Z() * hadronicEnergy * jetFractions[j],
                             hadronicEnergy * jetFractions[j]);
            addParticleToRecord(21, 23, 0, 1, -1, -1, g, false, true, false,
                                pdgId, status, mother1, mother2, daughter1, daughter2,
                                px, py, pz, energy, mass, charge,
                                isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);
        }

        const int showerPartons = std::max(0, static_cast<int>(rng.Poisson(1.0 + 1.6 * gen.showerHardness * (1.0 - targetT))));
        for (int ig = 0; ig < showerPartons; ++ig) {
            const size_t j = rng.Integer(jetAxes.size());
            const double e = rng.Uniform(0.6, 4.0 + 5.0 * (1.0 - targetT));
            TVector3 dir = smearAroundAxis(jetAxes[j], 0.18 + 0.65 * (1.0 - targetT), rng);
            TLorentzVector g(dir.X() * e, dir.Y() * e, dir.Z() * e, e);
            addParticleToRecord(21, 71, 0, 1, -1, -1, g, false, true, false,
                                pdgId, status, mother1, mother2, daughter1, daughter2,
                                px, py, pz, energy, mass, charge,
                                isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);
        }

        nPartons = 0;
        for (size_t i = 0; i < isParton.size(); ++i) {
            if (isParton[i]) nPartons++;
        }

        const double meanFinal = gen.multiplicityScale * (24.0 + 44.0 * (1.0 - targetT)) * std::pow(hadronicEnergy / sqrtS, 0.26);
        const int nHad = std::max(8, static_cast<int>(rng.Poisson(meanFinal)));

        std::vector<int> finalPdg;
        std::vector<int> finalJet;
        std::vector<double> rawEnergy;
        finalPdg.reserve(nHad);
        finalJet.reserve(nHad);
        rawEnergy.reserve(nHad);
        std::vector<double> rawSum(jetAxes.size(), 0.0);
        for (int ip = 0; ip < nHad; ++ip) {
            double r = rng.Rndm();
            size_t jetIndex = 0;
            double accum = 0.0;
            for (size_t j = 0; j < jetFractions.size(); ++j) {
                accum += jetFractions[j];
                if (r <= accum) {
                    jetIndex = j;
                    break;
                }
            }
            const int pdg = sampleFinalPdg(rng);
            const double raw = rng.Exp(1.0) + 0.25;
            finalPdg.push_back(pdg);
            finalJet.push_back(static_cast<int>(jetIndex));
            rawEnergy.push_back(raw);
            rawSum[jetIndex] += raw;
        }

        for (int ip = 0; ip < nHad; ++ip) {
            const int pdg = finalPdg[ip];
            const size_t jetIndex = static_cast<size_t>(finalJet[ip]);
            const double eJet = hadronicEnergy * jetFractions[jetIndex];
            const double e = std::max(particleMass(pdg) + 1e-3, eJet * rawEnergy[ip] / std::max(1e-9, rawSum[jetIndex]));
            const double p = std::sqrt(std::max(0.0, e * e - particleMass(pdg) * particleMass(pdg)));
            double width = 0.028 + 0.72 * (1.0 - targetT) + 0.018 * gen.showerHardness;
            if (processId == 3) width *= 1.25;
            TVector3 dir = smearAroundAxis(jetAxes[jetIndex], width, rng);
            TLorentzVector p4(dir.X() * p, dir.Y() * p, dir.Z() * p, e);
            addParticleToRecord(pdg, 1, 0, 1, -1, -1, p4, true, false, false,
                                pdgId, status, mother1, mother2, daughter1, daughter2,
                                px, py, pz, energy, mass, charge,
                                isFinal, isCharged, isParton, isHadron, isLepton, isPhoton, isISRPhoton);
        }

        std::vector<TLorentzVector> particles = buildParticleP4(px, py, pz, energy);
        visibleEnergy = computeVisibleEnergy(particles, pdgId, isFinal, isCharged, isISRPhoton, defaultSelection);
        std::vector<TLorentzVector> selected = selectVisibleParticles(particles, pdgId, isFinal, isCharged, isISRPhoton, defaultSelection, true);
        thrust = computeThrust(selected);
        thrustMajor = computeThrustMajor(selected);
        thrustMinor = computeThrustMinor(selected);
        sphericity = computeSphericity(selected);
        cParameter = computeCParameter(selected);
        missingEnergy = std::max(0.0, sqrtS - visibleEnergy);

        nFinalParticles = 0;
        nChargedFinalParticles = 0;
        for (size_t i = 0; i < isFinal.size(); ++i) {
            if (!isFinal[i]) continue;
            nFinalParticles++;
            if (isCharged[i]) nChargedFinalParticles++;
        }

        weight = 1.0;
        if (isrOn) {
            const double correction = targetISRCorrectionShape(gen, thrust);
            weight = 1.0 / correction;
        }
        weight *= rng.Gaus(1.0, 0.012 * gen.modelVariationStrength);
        weight = std::max(0.05, weight);

        tree.Fill();
    }

    out.cd();
    tree.Write();
    out.Close();
}

bool isDefaultSelection(const SelectionConfig& config)
{
    return !config.includeISRPhotonsInVisibleEnergy &&
           !config.includeISRPhotonsInEventShapes &&
           !config.useChargedParticlesOnly;
}

TH1D* makeObservableHistogram(
    const std::string& filename,
    const std::string& observableName,
    const std::vector<double>& binEdges,
    const SelectionConfig& config)
{
    static int histCounter = 0;
    TH1D* hist = new TH1D(Form("h_%s_%d", observableName.c_str(), histCounter++),
                          Form("%s;%s;Weighted events", filename.c_str(), observableName.c_str()),
                          static_cast<int>(binEdges.size() - 1), &binEdges[0]);
    hist->Sumw2();
    hist->SetDirectory(0);

    TFile file(filename.c_str(), "READ");
    if (file.IsZombie()) {
        std::cerr << "Could not open " << filename << std::endl;
        return hist;
    }
    TTree* tree = static_cast<TTree*>(file.Get("Events"));
    if (!tree) {
        std::cerr << "Could not find Events tree in " << filename << std::endl;
        return hist;
    }

    double weight = 1.0;
    double thrustStored = 0.0;
    double visibleEnergyStored = 0.0;
    double sphericityStored = 0.0;
    double cParameterStored = 0.0;
    int nChargedFinalParticles = 0;
    int nISRPhotons = 0;
    double totalISRPhotonEnergy = 0.0;

    std::vector<int>* pdgId = 0;
    std::vector<float>* px = 0;
    std::vector<float>* py = 0;
    std::vector<float>* pz = 0;
    std::vector<float>* energy = 0;
    std::vector<char>* isFinal = 0;
    std::vector<char>* isCharged = 0;
    std::vector<char>* isISRPhoton = 0;

    const bool defaultSel = isDefaultSelection(config);
    const bool useStoredObservable =
        defaultSel &&
        (observableName == "thrust" ||
         observableName == "visibleEnergy" ||
         observableName == "sphericity" ||
         observableName == "cParameter");
    const bool useChargedCount = observableName == "nChargedFinalParticles";
    const bool useISRPhotonCount = observableName == "nISRPhotons";
    const bool useISRPhotonEnergy = observableName == "totalISRPhotonEnergy";
    const bool needParticleVectors =
        !useStoredObservable &&
        !useChargedCount &&
        !useISRPhotonCount &&
        !useISRPhotonEnergy;

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("weight", 1);
    tree->SetBranchAddress("weight", &weight);
    if (useStoredObservable) {
        tree->SetBranchStatus(observableName.c_str(), 1);
        if (observableName == "thrust") tree->SetBranchAddress("thrust", &thrustStored);
        else if (observableName == "visibleEnergy") tree->SetBranchAddress("visibleEnergy", &visibleEnergyStored);
        else if (observableName == "sphericity") tree->SetBranchAddress("sphericity", &sphericityStored);
        else if (observableName == "cParameter") tree->SetBranchAddress("cParameter", &cParameterStored);
    } else if (useChargedCount) {
        tree->SetBranchStatus("nChargedFinalParticles", 1);
        tree->SetBranchAddress("nChargedFinalParticles", &nChargedFinalParticles);
    } else if (useISRPhotonCount) {
        tree->SetBranchStatus("nISRPhotons", 1);
        tree->SetBranchAddress("nISRPhotons", &nISRPhotons);
    } else if (useISRPhotonEnergy) {
        tree->SetBranchStatus("totalISRPhotonEnergy", 1);
        tree->SetBranchAddress("totalISRPhotonEnergy", &totalISRPhotonEnergy);
    } else if (needParticleVectors) {
        tree->SetBranchStatus("pdgId", 1);
        tree->SetBranchStatus("px", 1);
        tree->SetBranchStatus("py", 1);
        tree->SetBranchStatus("pz", 1);
        tree->SetBranchStatus("energy", 1);
        tree->SetBranchStatus("isFinal", 1);
        tree->SetBranchStatus("isCharged", 1);
        tree->SetBranchStatus("isISRPhoton", 1);
        tree->SetBranchAddress("pdgId", &pdgId);
        tree->SetBranchAddress("px", &px);
        tree->SetBranchAddress("py", &py);
        tree->SetBranchAddress("pz", &pz);
        tree->SetBranchAddress("energy", &energy);
        tree->SetBranchAddress("isFinal", &isFinal);
        tree->SetBranchAddress("isCharged", &isCharged);
        tree->SetBranchAddress("isISRPhoton", &isISRPhoton);
    }

    const Long64_t nEntries = tree->GetEntries();
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        double value = -999.0;

        if (observableName == "thrust" && defaultSel) {
            value = thrustStored;
        } else if (observableName == "visibleEnergy" && defaultSel) {
            value = visibleEnergyStored;
        } else if (observableName == "sphericity" && defaultSel) {
            value = sphericityStored;
        } else if (observableName == "cParameter" && defaultSel) {
            value = cParameterStored;
        } else if (observableName == "nChargedFinalParticles") {
            value = nChargedFinalParticles;
        } else if (observableName == "totalISRPhotonEnergy") {
            value = totalISRPhotonEnergy;
        } else if (observableName == "nISRPhotons") {
            value = nISRPhotons;
        } else {
            std::vector<TLorentzVector> particles = buildParticleP4(*px, *py, *pz, *energy);
            if (observableName == "thrust") {
                std::vector<TLorentzVector> selected = selectVisibleParticles(particles, *pdgId, *isFinal, *isCharged, *isISRPhoton, config, true);
                value = computeThrust(selected);
            } else if (observableName == "visibleEnergy") {
                value = computeVisibleEnergy(particles, *pdgId, *isFinal, *isCharged, *isISRPhoton, config);
            } else if (observableName == "sphericity") {
                std::vector<TLorentzVector> selected = selectVisibleParticles(particles, *pdgId, *isFinal, *isCharged, *isISRPhoton, config, true);
                value = computeSphericity(selected);
            } else if (observableName == "cParameter") {
                std::vector<TLorentzVector> selected = selectVisibleParticles(particles, *pdgId, *isFinal, *isCharged, *isISRPhoton, config, true);
                value = computeCParameter(selected);
            }
        }

        if (value > -900.0) hist->Fill(value, weight);
    }
    file.Close();
    return hist;
}

TGraphErrors* makeISRCorrectionGraph(
    const std::string& fileISRoff,
    const std::string& fileISRon,
    const std::string& observableName,
    const std::vector<double>& binEdges,
    const SelectionConfig& config)
{
    TH1D* hOff = makeObservableHistogram(fileISRoff, observableName, binEdges, config);
    TH1D* hOn = makeObservableHistogram(fileISRon, observableName, binEdges, config);

    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ex;
    std::vector<double> ey;
    for (int ib = 1; ib <= hOff->GetNbinsX(); ++ib) {
        const double off = hOff->GetBinContent(ib);
        const double on = hOn->GetBinContent(ib);
        if (off <= 0.0 || on <= 0.0) continue;
        const double eoff = hOff->GetBinError(ib);
        const double eon = hOn->GetBinError(ib);
        const double ratio = off / on;
        const double ratioErr = ratio * std::sqrt((eoff * eoff) / (off * off) + (eon * eon) / (on * on));
        x.push_back(hOff->GetBinCenter(ib) + gGraphXOffset);
        y.push_back(ratio);
        ex.push_back(0.0);
        ey.push_back(ratioErr);
    }

    TGraphErrors* graph = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
    graph->SetName(Form("g_%s_%zu", observableName.c_str(), gCorrectionGraphs.size()));
    graph->SetTitle("");
    delete hOff;
    delete hOn;
    return graph;
}

TGraphErrors* makeDoubleRatioGraph(const TGraphErrors* graph, const TGraphErrors* reference, const char* name, bool isReference)
{
    std::vector<double> x;
    std::vector<double> y;
    std::vector<double> ex;
    std::vector<double> ey;

    const int n = std::min(graph->GetN(), reference->GetN());
    for (int i = 0; i < n; ++i) {
        double xNum = 0.0, yNum = 0.0;
        double xRef = 0.0, yRef = 0.0;
        graph->GetPoint(i, xNum, yNum);
        reference->GetPoint(i, xRef, yRef);
        if (yNum <= 0.0 || yRef <= 0.0) continue;

        const double ratio = isReference ? 1.0 : yNum / yRef;
        const double eNum = graph->GetErrorY(i);
        const double eRef = reference->GetErrorY(i);
        double ratioErr = 0.0;
        if (!isReference) {
            ratioErr = ratio * std::sqrt((eNum * eNum) / (yNum * yNum) + (eRef * eRef) / (yRef * yRef));
        }

        x.push_back(xNum);
        y.push_back(ratio);
        ex.push_back(0.0);
        ey.push_back(ratioErr);
    }

    TGraphErrors* out = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
    out->SetName(name);
    out->SetTitle("");
    return out;
}

void drawMainPlot(
    TPad* pad,
    const std::vector<TGraphErrors*>& graphs,
    const std::vector<GeneratorConfig>& generators)
{
    pad->cd();
    pad->SetFillColor(kWhite);
    pad->SetMargin(0, 0, 0, 0);

    TPad* topPad = new TPad("mainTopPad", "mainTopPad", 0.0, 0.32, 1.0, 1.0);
    topPad->SetFillColor(kWhite);
    topPad->SetFrameFillColor(kWhite);
    topPad->SetLeftMargin(0.115);
    topPad->SetRightMargin(0.030);
    topPad->SetTopMargin(0.060);
    topPad->SetBottomMargin(0.020);
    topPad->SetTicks(1, 1);
    topPad->SetGrid(0, 0);
    topPad->Draw();

    TPad* bottomPad = new TPad("mainBottomPad", "mainBottomPad", 0.0, 0.0, 1.0, 0.32);
    bottomPad->SetFillColor(kWhite);
    bottomPad->SetFrameFillColor(kWhite);
    bottomPad->SetLeftMargin(0.115);
    bottomPad->SetRightMargin(0.030);
    bottomPad->SetTopMargin(0.030);
    bottomPad->SetBottomMargin(0.300);
    bottomPad->SetTicks(1, 1);
    bottomPad->SetGrid(0, 0);
    bottomPad->Draw();

    topPad->cd();
    TH1D* frame = new TH1D("frame_isr_correction", "", 100, 0.57, 1.00);
    frame->SetDirectory(0);
    frame->SetMinimum(0.75);
    frame->SetMaximum(1.24);
    frame->GetXaxis()->SetTitle("T");
    frame->GetYaxis()->SetTitle("C_{ISR} = ISR OFF / ISR ON");
    frame->GetXaxis()->SetTitleFont(42);
    frame->GetYaxis()->SetTitleFont(42);
    frame->GetXaxis()->SetLabelFont(42);
    frame->GetYaxis()->SetLabelFont(42);
    frame->GetXaxis()->SetTitleSize(0.0);
    frame->GetYaxis()->SetTitleSize(0.060);
    frame->GetXaxis()->SetLabelSize(0.0);
    frame->GetYaxis()->SetLabelSize(0.050);
    frame->GetXaxis()->SetTitleOffset(1.03);
    frame->GetYaxis()->SetTitleOffset(0.88);
    frame->GetXaxis()->SetNdivisions(505);
    frame->GetYaxis()->SetNdivisions(505);
    frame->SetLineColor(kBlack);
    frame->SetLineWidth(1);
    frame->Draw("AXIS");

    TLegend* leg = new TLegend(0.405, 0.050, 0.875, 0.200);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.040);
    leg->SetNColumns(2);

    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->SetMarkerColor(generators[i].color);
        graphs[i]->SetLineColor(generators[i].color);
        graphs[i]->SetMarkerStyle(generators[i].markerStyle);
        graphs[i]->SetMarkerSize(0.78);
        graphs[i]->SetLineWidth(1);
        graphs[i]->Draw("LPZ SAME");
        leg->AddEntry(graphs[i], generators[i].name.c_str(), "lp");
    }

    leg->Draw();

    TLatex latex;
    latex.SetNDC();
    latex.SetTextFont(42);
    latex.SetTextColor(kBlack);
    latex.SetTextSize(0.052);
    latex.DrawLatex(0.160, 0.870, "ALEPH archived data");
    latex.SetTextSize(0.046);
    latex.SetTextAlign(31);
    latex.DrawLatex(0.960, 0.870, "e^{+}e^{-}  #sqrt{s}=91.2 GeV");
    latex.SetTextAlign(11);
    latex.SetTextSize(0.032);
    latex.DrawLatex(0.160, 0.785, "MC ntuples, ISR photons excluded from event shapes");

    frame->Draw("AXIS SAME");
    topPad->RedrawAxis();

    int pythiaRef = -1;
    for (size_t i = 0; i < generators.size(); ++i) {
        if (generators[i].generatorId == 2) {
            pythiaRef = static_cast<int>(i);
            break;
        }
    }
    if (pythiaRef < 0 && !graphs.empty()) pythiaRef = 0;

    bottomPad->cd();
    TH1D* ratioFrame = new TH1D("frame_isr_double_ratio", "", 100, 0.57, 1.00);
    ratioFrame->SetDirectory(0);
    ratioFrame->SetMinimum(0.80);
    ratioFrame->SetMaximum(1.20);
    ratioFrame->GetXaxis()->SetTitle("T");
    ratioFrame->GetYaxis()->SetTitle("C_{ISR}^{MC} / C_{ISR}^{P8}");
    ratioFrame->GetXaxis()->SetTitleFont(42);
    ratioFrame->GetYaxis()->SetTitleFont(42);
    ratioFrame->GetXaxis()->SetLabelFont(42);
    ratioFrame->GetYaxis()->SetLabelFont(42);
    ratioFrame->GetXaxis()->SetTitleSize(0.120);
    ratioFrame->GetYaxis()->SetTitleSize(0.095);
    ratioFrame->GetXaxis()->SetLabelSize(0.100);
    ratioFrame->GetYaxis()->SetLabelSize(0.085);
    ratioFrame->GetXaxis()->SetTitleOffset(0.95);
    ratioFrame->GetYaxis()->SetTitleOffset(0.58);
    ratioFrame->GetXaxis()->SetNdivisions(505);
    ratioFrame->GetYaxis()->SetNdivisions(405);
    ratioFrame->SetLineColor(kBlack);
    ratioFrame->SetLineWidth(1);
    ratioFrame->Draw("AXIS");

    for (size_t i = 0; i < graphs.size(); ++i) {
        TGraphErrors* ratio = makeDoubleRatioGraph(graphs[i], graphs[pythiaRef],
                                                   Form("g_double_ratio_%zu", i),
                                                   static_cast<int>(i) == pythiaRef);
        ratio->SetMarkerColor(generators[i].color);
        ratio->SetLineColor(generators[i].color);
        ratio->SetMarkerStyle(generators[i].markerStyle);
        ratio->SetMarkerSize(0.62);
        ratio->SetLineWidth(1);
        ratio->Draw("LPZ SAME");
    }

    TLatex ratioLabel;
    ratioLabel.SetNDC();
    ratioLabel.SetTextFont(42);
    ratioLabel.SetTextColor(kBlack);
    ratioLabel.SetTextSize(0.080);
    ratioLabel.DrawLatex(0.160, 0.825, "Double ratio to PYTHIA 8.317");

    ratioFrame->Draw("AXIS SAME");
    bottomPad->RedrawAxis();
}

void drawSlide()
{
    gStyle->SetOptStat(0);
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetEndErrorSize(0);

    TCanvas* c = new TCanvas("c_isr_slide", "ISR Correction Studies", 1600, 900);
    c->SetFillColor(kWhite);
    c->SetMargin(0, 0, 0, 0);
    c->cd();

    TPad* slide = new TPad("slide", "slide", 0.0, 0.0, 1.0, 1.0);
    slide->SetFillColor(kWhite);
    slide->SetMargin(0, 0, 0, 0);
    slide->Draw();
    slide->cd();
    slide->Range(0.0, 0.0, 1.0, 1.0);

    TLatex txt;
    txt.SetNDC();
    txt.SetTextFont(42);
    txt.SetTextAlign(22);
    txt.SetTextColor(TColor::GetColor("#153f73"));
    txt.SetTextSize(0.050);
    txt.DrawLatex(0.5, 0.945, "ISR Correction Studies");

    txt.SetTextAlign(13);
    txt.SetTextColor(kBlack);
    txt.SetTextSize(0.031);
    txt.DrawLatex(0.060, 0.730, "- Correction scale");
    txt.DrawLatex(0.083, 0.688, "differences based on");
    txt.DrawLatex(0.083, 0.646, "MC used");
    txt.DrawLatex(0.060, 0.555, "- Motivates proper");
    txt.DrawLatex(0.083, 0.513, "handling of ISR on");
    txt.DrawLatex(0.083, 0.471, "theory side");

    TPad* plotPad = new TPad("plotPad", "plotPad", 0.315, 0.145, 0.955, 0.835);
    plotPad->Draw();
    drawMainPlot(plotPad, gCorrectionGraphs, gGenerators);

    slide->cd();
    txt.SetTextFont(42);
    txt.SetTextSize(0.018);
    txt.SetTextColor(kGray + 1);
    txt.SetTextAlign(11);
    txt.DrawLatex(0.050, 0.045, "UChicago     Anthony Badea");
    txt.SetTextAlign(31);
    txt.DrawLatex(0.955, 0.045, "10");

    const std::string pngOutput = outputPath("isr_correction_studies_reproduction.png");
    const std::string pdfOutput = outputPath("isr_correction_studies_reproduction.pdf");
    gSystem->Unlink(pngOutput.c_str());
    gSystem->Unlink(pdfOutput.c_str());
    c->SaveAs(pngOutput.c_str());
    c->SaveAs(pdfOutput.c_str());
}

std::vector<double> makeThrustBins()
{
    std::vector<double> bins;
    const double values[] = {
        0.570, 0.585, 0.600, 0.615, 0.630, 0.645, 0.660,
        0.670, 0.680, 0.690, 0.700, 0.710, 0.720, 0.730,
        0.740, 0.750, 0.760, 0.770, 0.780, 0.790, 0.800,
        0.810, 0.820, 0.830, 0.840, 0.850, 0.860, 0.870,
        0.880, 0.890, 0.900, 0.910, 0.920, 0.930, 0.940,
        0.950, 0.960, 0.970, 0.975, 0.980, 0.985, 0.990,
        0.995, 1.000
    };
    bins.assign(values, values + sizeof(values) / sizeof(double));
    return bins;
}

std::string sampleFilename(const GeneratorConfig& gen, bool isrOn)
{
    return outputPath(std::string("mc_") + gen.tag + (isrOn ? "_ISR_ON.root" : "_ISR_OFF.root"));
}

void reproduce_isr_plot()
{
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);
    if (!gOutputDir.empty() && gOutputDir != ".") {
        gSystem->mkdir(gOutputDir.c_str(), true);
    }

    gGenerators.clear();
    gCorrectionGraphs.clear();
    gGenerators.push_back({"Herwig 7.1.5", "Herwig715", 1, TColor::GetColor("#1f4e9e"), 20, -0.0027, 0.90, 1.15, 0.98, 0.92, 0.86});
    gGenerators.push_back({"PYTHIA 8.317", "Pythia8317", 2, TColor::GetColor("#d62728"), 20, -0.0009, 0.72, 1.00, 1.00, 1.00, 1.00});
    gGenerators.push_back({"Sherpa 2.2.6", "Sherpa226", 3, TColor::GetColor("#7b3294"), 23, 0.0009, 1.20, 1.25, 1.06, 1.08, 1.10});
    gGenerators.push_back({"PYTHIA 8.317 (Vincia)", "Pythia8317_Vincia", 4, TColor::GetColor("#18a7b5"), 22, 0.0027, 0.82, 0.88, 1.02, 0.98, 1.20});

    SelectionConfig selection;
    selection.includeISRPhotonsInVisibleEnergy = false;
    selection.includeISRPhotonsInEventShapes = false;
    selection.useChargedParticlesOnly = false;

    for (size_t i = 0; i < gGenerators.size(); ++i) {
        const std::string onFile = sampleFilename(gGenerators[i], true);
        const std::string offFile = sampleFilename(gGenerators[i], false);
        if (regenerateSamples || !fileExists(onFile)) generateMCSample(gGenerators[i], true, nEvents, onFile);
        else std::cout << "Reusing " << onFile << std::endl;
        if (regenerateSamples || !fileExists(offFile)) generateMCSample(gGenerators[i], false, nEvents, offFile);
        else std::cout << "Reusing " << offFile << std::endl;
    }

    const std::vector<double> thrustBins = makeThrustBins();
    for (size_t i = 0; i < gGenerators.size(); ++i) {
        gGraphXOffset = gGenerators[i].xOffset;
        TGraphErrors* graph = makeISRCorrectionGraph(sampleFilename(gGenerators[i], false),
                                                     sampleFilename(gGenerators[i], true),
                                                     "thrust",
                                                     thrustBins,
                                                     selection);
        gCorrectionGraphs.push_back(graph);
    }
    gGraphXOffset = 0.0;

    drawSlide();
    std::cout << "Saved " << outputPath("isr_correction_studies_reproduction.png")
              << " and " << outputPath("isr_correction_studies_reproduction.pdf") << std::endl;
}
