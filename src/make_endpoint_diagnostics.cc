#include <TFile.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

struct Options {
    std::string input;
    std::string output;
    double sqrtS = 91.1876;
    Long64_t maxEvents = -1;
    int maxExactParticles = 0;
};

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;

    Vec3() = default;
    Vec3(double xIn, double yIn, double zIn) : x(xIn), y(yIn), z(zIn) {}

    Vec3 operator+(const Vec3& other) const { return Vec3(x + other.x, y + other.y, z + other.z); }
    Vec3 operator-(const Vec3& other) const { return Vec3(x - other.x, y - other.y, z - other.z); }
    Vec3 operator*(double scale) const { return Vec3(x * scale, y * scale, z * scale); }
    Vec3& operator+=(const Vec3& other)
    {
        x += other.x;
        y += other.y;
        z += other.z;
        return *this;
    }
    Vec3& operator-=(const Vec3& other)
    {
        x -= other.x;
        y -= other.y;
        z -= other.z;
        return *this;
    }

    double dot(const Vec3& other) const { return x * other.x + y * other.y + z * other.z; }
    Vec3 cross(const Vec3& other) const
    {
        return Vec3(y * other.z - z * other.y,
                    z * other.x - x * other.z,
                    x * other.y - y * other.x);
    }
    double mag2() const { return dot(*this); }
    double mag() const { return std::sqrt(mag2()); }
};

struct Vec4 {
    Vec3 p;
    double e = 0.0;
};

struct ThrustResult {
    double value = 0.0;
    Vec3 axis;
};

struct OutBranches {
    int eventId = 0;
    std::string generatorName;
    bool isrOn = true;

    double T_lab_allFinal_including_ISR_photons = 0.0;
    double T_lab_allFinal_excluding_ISR_photons = 0.0;
    double T_allFinalCM_including_ISR_photons = 0.0;
    double T_lab_excluding_ISR_photons = 0.0;
    double T_lab_including_ISR_photons = 0.0;
    double T_visibleCM_excluding_ISR_photons = 0.0;
    double T_lab_hadron_excluding_ISR_photons = 0.0;
    double T_lab_hadron_including_ISR_photons = 0.0;
    double T_hadronicCM_excluding_ISR_photons = 0.0;
    double cosTheta_thrust_lab_allFinal_including_ISR_photons = 0.0;
    double absCosTheta_thrust_lab_allFinal_including_ISR_photons = 0.0;
    double cosTheta_thrust_lab_excluding_ISR_photons = 0.0;
    double absCosTheta_thrust_lab_excluding_ISR_photons = 0.0;
    double Eall_including_ISR_photons = 0.0;
    double Eall_excluding_ISR_photons = 0.0;
    double Mall_including_ISR_photons = 0.0;
    double Mall_excluding_ISR_photons = 0.0;
    double sAll_including_ISR_over_s = 0.0;
    double betaZ_all_including_ISR = 0.0;
    double Evis_excluding_ISR_photons = 0.0;
    double Evis_including_ISR_photons = 0.0;
    double Mvis_excluding_ISR_photons = 0.0;
    double Mvis_including_ISR_photons = 0.0;
    double sPrime_vis_over_s = 0.0;
    double betaZ_vis = 0.0;
    double leading_ISR_photon_energy = -1.0;
    double leading_ISR_photon_theta = -1.0;
    int N_ISR_photons = 0;
    double sum_ISR_photon_energy = 0.0;
    double event_weight = 1.0;
};

void usage()
{
    std::cout
        << "Usage:\n"
        << "  make_endpoint_diagnostics --input mc.root --output diagnostics.root [--sqrtS 91.1876]\n"
        << "                            [--maxEvents N] [--maxExactParticles N]\n\n"
        << "The nominal diagnostic thrust input is all stable final-state particles,\n"
        << "including neutrinos and tagged ISR photons. Visible/non-neutrino branches\n"
        << "are also stored for comparison with the earlier ALEPH-style diagnostic.\n";
}

Options parseOptions(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        auto value = [&](const std::string& name) -> std::string {
            if (i + 1 >= argc) throw std::runtime_error("Missing value for " + name);
            return argv[++i];
        };
        if (key == "--input") opt.input = value(key);
        else if (key == "--output") opt.output = value(key);
        else if (key == "--sqrtS") opt.sqrtS = std::atof(value(key).c_str());
        else if (key == "--maxEvents") opt.maxEvents = std::atoll(value(key).c_str());
        else if (key == "--maxExactParticles") opt.maxExactParticles = std::atoi(value(key).c_str());
        else if (key == "--help") {
            usage();
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown option " + key);
        }
    }
    if (opt.input.empty()) throw std::runtime_error("--input is required");
    if (opt.output.empty()) throw std::runtime_error("--output is required");
    return opt;
}

double invariantMass(const Vec4& p4)
{
    const double m2 = p4.e * p4.e - p4.p.mag2();
    return m2 > 0.0 ? std::sqrt(m2) : 0.0;
}

double theta(const Vec3& p)
{
    const double mag = p.mag();
    if (mag <= 1e-12) return -1.0;
    double c = p.z / mag;
    c = std::max(-1.0, std::min(1.0, c));
    return std::acos(c);
}

double absCosThetaAxis(const Vec3& axis)
{
    const double mag = axis.mag();
    if (mag <= 1e-12) return 0.0;
    return std::abs(axis.z / mag);
}

bool isNeutrino(int pdgId)
{
    const int apdg = std::abs(pdgId);
    return apdg == 12 || apdg == 14 || apdg == 16;
}

ThrustResult thrustBelleLikeResult(const std::vector<Vec4>& particles)
{
    std::vector<Vec3> momenta;
    momenta.reserve(particles.size());
    double pSum = 0.0;
    for (const Vec4& p4 : particles) {
        const double mag = p4.p.mag();
        if (mag <= 1e-12) continue;
        pSum += mag;
        momenta.push_back(p4.p);
    }
    if (pSum <= 1e-12 || momenta.empty()) return ThrustResult();

    double best = 0.0;
    Vec3 bestAxis;
    for (const Vec3& seed : momenta) {
        Vec3 axis = seed;
        if (axis.z <= 0.0) axis = axis * -1.0;
        const double seedMag = axis.mag();
        if (seedMag <= 1e-12) continue;
        axis = axis * (1.0 / seedMag);

        const size_t maxIterations = std::min<size_t>(16, momenta.size());
        for (size_t iter = 0; iter < maxIterations; ++iter) {
            const Vec3 previous = axis;
            axis = Vec3();
            for (const Vec3& p : momenta) axis += p.dot(previous) >= 0.0 ? p : p * -1.0;
            bool converged = true;
            for (const Vec3& p : momenta) {
                if (p.dot(axis) * p.dot(previous) < 0.0) {
                    converged = false;
                    break;
                }
            }
            if (converged) break;
        }

        const double axisMag = axis.mag();
        if (axisMag <= 1e-12) continue;
        double sum = 0.0;
        for (const Vec3& p : momenta) sum += std::abs(p.dot(axis));
        const double t = sum / (pSum * axisMag);
        if (t > best) {
            best = t;
            bestAxis = axis * (1.0 / axisMag);
        }
    }
    return ThrustResult{std::min(1.0, std::max(0.0, best)), bestAxis};
}

ThrustResult thrustExactResult(const std::vector<Vec4>& particles)
{
    std::vector<Vec3> momenta;
    momenta.reserve(particles.size());
    double pSum = 0.0;
    for (const Vec4& p4 : particles) {
        const double mag = p4.p.mag();
        if (mag <= 1e-12) continue;
        pSum += mag;
        momenta.push_back(p4.p);
    }
    const size_t n = momenta.size();
    if (pSum <= 1e-12 || n == 0) return ThrustResult();
    if (n < 4) return thrustBelleLikeResult(particles);

    Vec3 bestVector;
    double bestMag2 = -1.0;
    for (size_t i = 1; i < n; ++i) {
        for (size_t j = 0; j < i; ++j) {
            const Vec3 cross = momenta[i].cross(momenta[j]);
            Vec3 ptot;
            for (size_t k = 0; k < n; ++k) {
                if (k == i || k == j) continue;
                ptot += momenta[k].dot(cross) > 0.0 ? momenta[k] : momenta[k] * -1.0;
            }
            const Vec3 c0 = ptot - momenta[j] - momenta[i];
            const Vec3 c1 = ptot - momenta[j] + momenta[i];
            const Vec3 c2 = ptot + momenta[j] - momenta[i];
            const Vec3 c3 = ptot + momenta[j] + momenta[i];
            const Vec3 candidates[4] = {c0, c1, c2, c3};
            for (const Vec3& candidate : candidates) {
                const double mag2 = candidate.mag2();
                if (mag2 > bestMag2) {
                    bestMag2 = mag2;
                    bestVector = candidate;
                }
            }
        }
    }

    const double bestMag = bestVector.mag();
    const Vec3 axis = bestMag > 1e-12 ? bestVector * (1.0 / bestMag) : Vec3();
    return ThrustResult{std::min(1.0, std::max(0.0, bestMag / pSum)), axis};
}

ThrustResult computeThrustResult(const std::vector<Vec4>& particles, int maxExactParticles)
{
    if (static_cast<int>(particles.size()) > maxExactParticles) return thrustBelleLikeResult(particles);
    return thrustExactResult(particles);
}

double computeThrust(const std::vector<Vec4>& particles, int maxExactParticles)
{
    return computeThrustResult(particles, maxExactParticles).value;
}

Vec4 sumP4(const std::vector<Vec4>& particles)
{
    Vec4 total;
    for (const Vec4& p4 : particles) {
        total.p += p4.p;
        total.e += p4.e;
    }
    return total;
}

Vec4 boost(const Vec4& p4, const Vec3& beta)
{
    const double b2 = beta.mag2();
    if (b2 <= 1e-16) return p4;
    if (b2 >= 1.0) return p4;
    const double gamma = 1.0 / std::sqrt(1.0 - b2);
    const double bp = beta.dot(p4.p);
    const double gamma2 = (gamma - 1.0) / b2;

    Vec4 out;
    out.p = p4.p + beta * (gamma2 * bp + gamma * p4.e);
    out.e = gamma * (p4.e + bp);
    return out;
}

double computeHadronicCMThrust(const std::vector<Vec4>& hadrons, int maxExactParticles)
{
    const Vec4 total = sumP4(hadrons);
    if (total.e <= 1e-12) return 0.0;
    const Vec3 beta = total.p * (-1.0 / total.e);
    if (beta.mag2() >= 1.0) return 0.0;

    std::vector<Vec4> boosted;
    boosted.reserve(hadrons.size());
    for (const Vec4& p4 : hadrons) boosted.push_back(boost(p4, beta));
    return computeThrust(boosted, maxExactParticles);
}

double computeRestFrameThrust(const std::vector<Vec4>& particles, int maxExactParticles)
{
    const Vec4 total = sumP4(particles);
    if (total.e <= 1e-12) return 0.0;
    const Vec3 beta = total.p * (-1.0 / total.e);
    if (beta.mag2() >= 1.0) return 0.0;

    std::vector<Vec4> boosted;
    boosted.reserve(particles.size());
    for (const Vec4& p4 : particles) boosted.push_back(boost(p4, beta));
    return computeThrust(boosted, maxExactParticles);
}

void bookTree(TTree& tree, OutBranches& b)
{
    tree.Branch("eventId", &b.eventId);
    tree.Branch("generatorName", &b.generatorName);
    tree.Branch("isrOn", &b.isrOn);
    tree.Branch("T_lab_allFinal_including_ISR_photons", &b.T_lab_allFinal_including_ISR_photons);
    tree.Branch("T_lab_allFinal_excluding_ISR_photons", &b.T_lab_allFinal_excluding_ISR_photons);
    tree.Branch("T_allFinalCM_including_ISR_photons", &b.T_allFinalCM_including_ISR_photons);
    tree.Branch("T_lab_excluding_ISR_photons", &b.T_lab_excluding_ISR_photons);
    tree.Branch("T_lab_including_ISR_photons", &b.T_lab_including_ISR_photons);
    tree.Branch("T_visibleCM_excluding_ISR_photons", &b.T_visibleCM_excluding_ISR_photons);
    tree.Branch("T_lab_hadron_excluding_ISR_photons", &b.T_lab_hadron_excluding_ISR_photons);
    tree.Branch("T_lab_hadron_including_ISR_photons", &b.T_lab_hadron_including_ISR_photons);
    tree.Branch("T_hadronicCM_excluding_ISR_photons", &b.T_hadronicCM_excluding_ISR_photons);
    tree.Branch("cosTheta_thrust_lab_allFinal_including_ISR_photons", &b.cosTheta_thrust_lab_allFinal_including_ISR_photons);
    tree.Branch("absCosTheta_thrust_lab_allFinal_including_ISR_photons", &b.absCosTheta_thrust_lab_allFinal_including_ISR_photons);
    tree.Branch("cosTheta_thrust_lab_excluding_ISR_photons", &b.cosTheta_thrust_lab_excluding_ISR_photons);
    tree.Branch("absCosTheta_thrust_lab_excluding_ISR_photons", &b.absCosTheta_thrust_lab_excluding_ISR_photons);
    tree.Branch("Eall_including_ISR_photons", &b.Eall_including_ISR_photons);
    tree.Branch("Eall_excluding_ISR_photons", &b.Eall_excluding_ISR_photons);
    tree.Branch("Mall_including_ISR_photons", &b.Mall_including_ISR_photons);
    tree.Branch("Mall_excluding_ISR_photons", &b.Mall_excluding_ISR_photons);
    tree.Branch("sAll_including_ISR_over_s", &b.sAll_including_ISR_over_s);
    tree.Branch("betaZ_all_including_ISR", &b.betaZ_all_including_ISR);
    tree.Branch("Evis_excluding_ISR_photons", &b.Evis_excluding_ISR_photons);
    tree.Branch("Evis_including_ISR_photons", &b.Evis_including_ISR_photons);
    tree.Branch("Mvis_excluding_ISR_photons", &b.Mvis_excluding_ISR_photons);
    tree.Branch("Mvis_including_ISR_photons", &b.Mvis_including_ISR_photons);
    tree.Branch("sPrime_vis_over_s", &b.sPrime_vis_over_s);
    tree.Branch("betaZ_vis", &b.betaZ_vis);
    tree.Branch("leading_ISR_photon_energy", &b.leading_ISR_photon_energy);
    tree.Branch("leading_ISR_photon_theta", &b.leading_ISR_photon_theta);
    tree.Branch("N_ISR_photons", &b.N_ISR_photons);
    tree.Branch("sum_ISR_photon_energy", &b.sum_ISR_photon_energy);
    tree.Branch("event_weight", &b.event_weight);
}

int run(const Options& opt)
{
    TFile input(opt.input.c_str(), "READ");
    if (input.IsZombie()) throw std::runtime_error("Could not open input " + opt.input);
    TTree* events = static_cast<TTree*>(input.Get("Events"));
    if (!events) throw std::runtime_error("Missing Events tree in " + opt.input);

    TTreeReader reader(events);
    TTreeReaderValue<int> eventId(reader, "eventId");
    TTreeReaderValue<std::string> generatorName(reader, "generatorName");
    TTreeReaderValue<bool> isrOn(reader, "isrOn");
    TTreeReaderValue<double> weight(reader, "weight");
    TTreeReaderValue<std::vector<int>> pdgId(reader, "pdgId");
    TTreeReaderValue<std::vector<float>> px(reader, "px");
    TTreeReaderValue<std::vector<float>> py(reader, "py");
    TTreeReaderValue<std::vector<float>> pz(reader, "pz");
    TTreeReaderValue<std::vector<float>> energy(reader, "energy");
    TTreeReaderValue<std::vector<char>> isFinal(reader, "isFinal");
    TTreeReaderValue<std::vector<char>> isHadron(reader, "isHadron");
    TTreeReaderValue<std::vector<char>> isISRPhoton(reader, "isISRPhoton");

    TFile output(opt.output.c_str(), "RECREATE");
    if (output.IsZombie()) throw std::runtime_error("Could not open output " + opt.output);
    TTree tree("EndpointDiagnostics", "Per-event ISR endpoint diagnostics");
    OutBranches out;
    bookTree(tree, out);

    const double s = opt.sqrtS * opt.sqrtS;
    Long64_t n = 0;
    while (reader.Next()) {
        if (opt.maxEvents >= 0 && n >= opt.maxEvents) break;

        std::vector<Vec4> allFinalIncludingISR;
        std::vector<Vec4> allFinalExcludingISR;
        std::vector<Vec4> visible;
        std::vector<Vec4> visibleAndISR;
        std::vector<Vec4> hadrons;
        std::vector<Vec4> hadronsAndISR;
        allFinalIncludingISR.reserve(px->size());
        allFinalExcludingISR.reserve(px->size());
        visible.reserve(px->size());
        visibleAndISR.reserve(px->size());
        hadrons.reserve(px->size());
        hadronsAndISR.reserve(px->size());

        out = OutBranches();
        out.eventId = *eventId;
        out.generatorName = *generatorName;
        out.isrOn = *isrOn;
        out.event_weight = *weight;

        for (size_t i = 0; i < px->size(); ++i) {
            if (!(*isFinal)[i]) continue;
            const Vec4 p4{Vec3((*px)[i], (*py)[i], (*pz)[i]), (*energy)[i]};
            allFinalIncludingISR.push_back(p4);
            if (!(*isISRPhoton)[i]) allFinalExcludingISR.push_back(p4);
            const bool neutrino = isNeutrino((*pdgId)[i]);
            if (!neutrino) {
                visibleAndISR.push_back(p4);
                if (!(*isISRPhoton)[i]) visible.push_back(p4);
            }
            if ((*isHadron)[i]) {
                hadrons.push_back(p4);
                hadronsAndISR.push_back(p4);
            }
            if ((*isISRPhoton)[i]) {
                hadronsAndISR.push_back(p4);
                ++out.N_ISR_photons;
                out.sum_ISR_photon_energy += p4.e;
                if (p4.e > out.leading_ISR_photon_energy) {
                    out.leading_ISR_photon_energy = p4.e;
                    out.leading_ISR_photon_theta = theta(p4.p);
                }
            }
        }

        const Vec4 p4AllIncl = sumP4(allFinalIncludingISR);
        const Vec4 p4AllExcl = sumP4(allFinalExcludingISR);
        const Vec4 p4Excl = sumP4(visible);
        const Vec4 p4Incl = sumP4(visibleAndISR);

        const ThrustResult allLabThrust = computeThrustResult(allFinalIncludingISR, opt.maxExactParticles);
        const ThrustResult labThrust = computeThrustResult(visible, opt.maxExactParticles);
        out.T_lab_allFinal_including_ISR_photons = allLabThrust.value;
        out.T_lab_allFinal_excluding_ISR_photons = computeThrust(allFinalExcludingISR, opt.maxExactParticles);
        out.T_allFinalCM_including_ISR_photons = computeRestFrameThrust(allFinalIncludingISR, opt.maxExactParticles);
        out.T_lab_excluding_ISR_photons = labThrust.value;
        out.T_lab_including_ISR_photons = computeThrust(visibleAndISR, opt.maxExactParticles);
        out.T_visibleCM_excluding_ISR_photons = computeRestFrameThrust(visible, opt.maxExactParticles);
        out.T_lab_hadron_excluding_ISR_photons = computeThrust(hadrons, opt.maxExactParticles);
        out.T_lab_hadron_including_ISR_photons = computeThrust(hadronsAndISR, opt.maxExactParticles);
        out.T_hadronicCM_excluding_ISR_photons = computeHadronicCMThrust(hadrons, opt.maxExactParticles);
        out.cosTheta_thrust_lab_allFinal_including_ISR_photons = allLabThrust.axis.mag() > 1e-12 ? allLabThrust.axis.z / allLabThrust.axis.mag() : 0.0;
        out.absCosTheta_thrust_lab_allFinal_including_ISR_photons = absCosThetaAxis(allLabThrust.axis);
        out.cosTheta_thrust_lab_excluding_ISR_photons = labThrust.axis.mag() > 1e-12 ? labThrust.axis.z / labThrust.axis.mag() : 0.0;
        out.absCosTheta_thrust_lab_excluding_ISR_photons = absCosThetaAxis(labThrust.axis);
        out.Eall_including_ISR_photons = p4AllIncl.e;
        out.Eall_excluding_ISR_photons = p4AllExcl.e;
        out.Mall_including_ISR_photons = invariantMass(p4AllIncl);
        out.Mall_excluding_ISR_photons = invariantMass(p4AllExcl);
        out.sAll_including_ISR_over_s = s > 0.0 ? (out.Mall_including_ISR_photons * out.Mall_including_ISR_photons) / s : 0.0;
        out.betaZ_all_including_ISR = p4AllIncl.e > 1e-12 ? std::abs(p4AllIncl.p.z) / p4AllIncl.e : 0.0;
        out.Evis_excluding_ISR_photons = p4Excl.e;
        out.Evis_including_ISR_photons = p4Incl.e;
        out.Mvis_excluding_ISR_photons = invariantMass(p4Excl);
        out.Mvis_including_ISR_photons = invariantMass(p4Incl);
        out.sPrime_vis_over_s = s > 0.0 ? (out.Mvis_excluding_ISR_photons * out.Mvis_excluding_ISR_photons) / s : 0.0;
        out.betaZ_vis = p4Excl.e > 1e-12 ? std::abs(p4Excl.p.z) / p4Excl.e : 0.0;

        tree.Fill();
        ++n;
        if (n % 250000 == 0) std::cout << "Processed " << n << " events from " << opt.input << std::endl;
    }

    tree.Write();
    output.Close();
    std::cout << "Wrote " << n << " events to " << opt.output << std::endl;
    return 0;
}

} // namespace

int main(int argc, char** argv)
{
    try {
        return run(parseOptions(argc, argv));
    } catch (const std::exception& e) {
        std::cerr << "make_endpoint_diagnostics: " << e.what() << std::endl;
        return 1;
    }
}
