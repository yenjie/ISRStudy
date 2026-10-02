// Is the KKMC endpoint difference the DIRECTION of its ISR photons?
//
// Recompute the ALEPH-definition thrust for the KKMC pair twice: once with the
// generator's ISR photons as they are, and once with each ISR photon replaced
// by a photon of the same energy exactly along the beam axis (the sign of its
// pz kept).  Everything else in the event is untouched.  If the collinearised
// C_ISR lands where the collinear-structure-function generators sit, the
// difference is the photon direction and nothing else.
//
// The first pass also reproduces the stored diagnostics value, which checks
// this thrust implementation against the producer's.
//
// Usage:
//   root -l -b -q 'macros/kkmc_collinearize_isr_test.C("<out dir>")'

#include <cmath>
#include <cstdio>
#include <fstream>
#include <string>
#include <vector>

namespace {

const char* kKkmcNtup = "/data2/yjlee/ISRsample/kkmc_1M_20260921";

struct V3 { double x, y, z; };
inline double dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }
inline double mag(const V3& a) { return std::sqrt(dot(a, a)); }

// Iterative thrust from every particle direction as a seed; converges to the
// global maximum for the multiplicities here and matches the producer to the
// precision that matters for bin counts.
double thrustOf(const std::vector<V3>& p)
{
    if (p.size() < 2) return -1;
    double sumMag = 0;
    for (const auto& q : p) sumMag += mag(q);
    double best = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        const double m = mag(p[i]);
        if (!(m > 0)) continue;
        V3 ax{p[i].x / m, p[i].y / m, p[i].z / m};
        for (int it = 0; it < 20; ++it) {
            V3 nw{0, 0, 0};
            for (const auto& q : p) {
                const double s = dot(q, ax) > 0 ? 1.0 : -1.0;
                nw.x += s * q.x; nw.y += s * q.y; nw.z += s * q.z;
            }
            const double nm = mag(nw);
            if (!(nm > 0)) break;
            nw.x /= nm; nw.y /= nm; nw.z /= nm;
            const bool same = std::fabs(nw.x - ax.x) + std::fabs(nw.y - ax.y) + std::fabs(nw.z - ax.z) < 1e-12;
            ax = nw;
            if (same) break;
        }
        double num = 0;
        for (const auto& q : p) num += std::fabs(dot(q, ax));
        if (sumMag > 0 && num / sumMag > best) best = num / sumMag;
    }
    return best;
}

bool isNeutrino(int pdg) { const int a = std::abs(pdg); return a == 12 || a == 14 || a == 16; }

struct Counts { Long64_t events = 0; double sumW = 0, lastActual = 0, lastColl = 0, lastNoIsr = 0;
                double lastActual2 = 0, lastColl2 = 0, lastNoIsr2 = 0;  // sum w^2 for errors
                std::vector<double> hA, hC, hN; };

Counts scan(const std::string& path, Long64_t maxEvents)
{
    Counts r;
    r.hA.assign(42, 0); r.hC.assign(42, 0); r.hN.assign(42, 0);
    TFile* f = TFile::Open(path.c_str());
    TTree* t = static_cast<TTree*>(f->Get("Events"));
    std::vector<char>* isFinal = nullptr; std::vector<char>* isIsr = nullptr;
    std::vector<int>* pdg = nullptr;
    std::vector<float>*px = nullptr, *py = nullptr, *pz = nullptr, *en = nullptr;
    double weight = 1.0;
    t->SetBranchStatus("*", 0);
    for (const char* b : {"isFinal", "isISRPhoton", "pdgId", "px", "py", "pz", "energy", "weight"}) t->SetBranchStatus(b, 1);
    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("isFinal", &isFinal); t->SetBranchAddress("isISRPhoton", &isIsr);
    t->SetBranchAddress("pdgId", &pdg);
    t->SetBranchAddress("px", &px); t->SetBranchAddress("py", &py); t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("energy", &en);
    Long64_t n = t->GetEntries();
    if (maxEvents > 0 && n > maxEvents) n = maxEvents;
    std::vector<V3> actual, coll, noIsr;
    auto binOf = [](double T) { const int b = static_cast<int>(std::floor((T - 0.58) / 0.01)); return (b >= 0 && b < 42) ? b : -1; };
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        actual.clear(); coll.clear(); noIsr.clear();
        for (size_t j = 0; j < isFinal->size(); ++j) {
            if (!(*isFinal)[j] || isNeutrino((*pdg)[j])) continue;
            const V3 p{(*px)[j], (*py)[j], (*pz)[j]};
            actual.push_back(p);
            if ((*isIsr)[j]) {
                coll.push_back(V3{0, 0, (p.z >= 0 ? 1.0 : -1.0) * static_cast<double>((*en)[j])});
            } else {
                coll.push_back(p);
                noIsr.push_back(p);
            }
        }
        const double tA = thrustOf(actual), tC = thrustOf(coll), tN = thrustOf(noIsr);
        const double w = weight;
        int b;
        if ((b = binOf(tA)) >= 0) r.hA[b] += w;
        if ((b = binOf(tC)) >= 0) r.hC[b] += w;
        if ((b = binOf(tN)) >= 0) r.hN[b] += w;
        if (tA >= 0.99) { r.lastActual += w; r.lastActual2 += w * w; }
        if (tC >= 0.99) { r.lastColl += w; r.lastColl2 += w * w; }
        if (tN >= 0.99) { r.lastNoIsr += w; r.lastNoIsr2 += w * w; }
        r.sumW += w;
        ++r.events;
        if ((i + 1) % 200000 == 0) printf("    %lld/%lld\n", i + 1, n);
    }
    f->Close();
    return r;
}

}  // namespace

void kkmc_collinearize_isr_test(const char* outDir = "/data2/yjlee/ISRsample/kkmc_1M_20260921/results",
                                Long64_t maxEvents = -1)
{
    printf("[off]\n");
    Counts off = scan(std::string(kKkmcNtup) + "/mc_KKMC424_ISR_OFF.root", maxEvents);
    printf("[on]\n");
    Counts on = scan(std::string(kKkmcNtup) + "/mc_KKMC424_ISR_ON.root", maxEvents);
    // Weighted shape ratio: each state normalised to its own sum of weights.
    auto ratio = [&](double a, double a2, double c, double c2, double& r, double& e) {
        const double fa = a / off.sumW, fc = c / on.sumW;
        r = fc > 0 ? fa / fc : 0;
        e = (a > 0 && c > 0) ? r * std::sqrt(a2 / (a * a) + c2 / (c * c)) : 0; };
    double rA, eA, rC, eC, rN, eN;
    ratio(off.lastActual, off.lastActual2, on.lastActual, on.lastActual2, rA, eA);
    ratio(off.lastColl, off.lastColl2, on.lastColl, on.lastColl2, rC, eC);
    ratio(off.lastNoIsr, off.lastNoIsr2, on.lastNoIsr, on.lastNoIsr2, rN, eN);
    printf("  sum of weights: OFF %.1f over %lld events, ON %.1f over %lld\n", off.sumW, off.events, on.sumW, on.events);
    printf("\n=== KKMC 4.30, 0.99 < T < 1.00, ALEPH thrust definition ===\n");
    printf("  ISR photons as generated         C_ISR = %.4f +- %.4f   (unweighted diagnostics gave 1.094)\n", rA, eA);
    printf("  ISR photons forced beam-collinear C_ISR = %.4f +- %.4f\n", rC, eC);
    printf("  ISR photons removed               C_ISR = %.4f +- %.4f   (unweighted diagnostics A gave 0.958)\n", rN, eN);
    std::ofstream csv(std::string(outDir) + "/kkmc_collinearize_isr_test.csv");
    csv << "variant,sumw_off_lastbin,sumw_on_lastbin,c_isr,c_isr_err\n";
    csv << "as generated," << off.lastActual << "," << on.lastActual << "," << rA << "," << eA << "\n";
    csv << "forced collinear," << off.lastColl << "," << on.lastColl << "," << rC << "," << eC << "\n";
    csv << "removed," << off.lastNoIsr << "," << on.lastNoIsr << "," << rN << "," << eN << "\n";
    csv.close();
    std::ofstream bins(std::string(outDir) + "/kkmc_collinearize_isr_test_bins.csv");
    bins << "T_low,T_high,off_actual,on_actual,off_coll,on_coll,off_noisr,on_noisr\n";
    for (int b = 0; b < 42; ++b)
        bins << 0.58 + 0.01 * b << "," << 0.59 + 0.01 * b << "," << off.hA[b] << "," << on.hA[b] << ","
             << off.hC[b] << "," << on.hC[b] << "," << off.hN[b] << "," << on.hN[b] << "\n";
    bins.close();
    printf("[done] %s/kkmc_collinearize_isr_test.csv\n", outDir);
}
