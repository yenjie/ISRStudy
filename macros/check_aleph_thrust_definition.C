// Does the ISR-study thrust definition match what the ALEPH thrust analysis needs?
//
// The ALEPH event-shape analysis applies
//   C_ISR(j) = N_noISR(j) / N_ISR(j)
// built from the `tgenBefore/thrust` branch of its own standalone Pythia8
// samples.  That branch's particle selection is not written down anywhere, so
// this macro determines it by brute force: recompute thrust from the stored
// particles under every plausible particle-level definition and see which one
// reproduces the stored branch.
//
// It then quantifies how much the choice actually matters, by computing
// C_ISR in the last ALEPH thrust bin under each definition.
//
// Usage:
//   root -l -b -q macros/check_aleph_thrust_definition.C

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

namespace {

const char* kIsrOff = "/raid5/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/Isr/isr0_ALL.root";
const char* kIsrOn = "/raid5/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/Isr/isr1_ALL.root";

// Candidate particle-level definitions.
enum Def {
    kAllFinal = 0,      // every status>0 particle, neutrinos included
    kNoNeutrino,        // status>0, neutrinos removed        <-- expected match
    kChargedOnly,       // status>0, charged only
    kNoNuAcc094,        // status>0, no neutrinos, |cos(theta)| < 0.94
    kNoNuP200,          // status>0, no neutrinos, p > 0.2 GeV
    kNDef
};
const char* kDefName[kNDef] = {
    "all final state (neutrinos in)",
    "final state, no neutrinos",
    "charged only",
    "no neutrinos, |cos t| < 0.94",
    "no neutrinos, p > 0.2 GeV"};

bool isNeutrino(int pdg)
{
    const int a = std::abs(pdg);
    return a == 12 || a == 14 || a == 16;
}

double thrustOf(const std::vector<TVector3>& p)
{
    if (p.size() < 2) return -1;
    double best = 0;
    for (size_t i = 0; i < p.size(); ++i) {
        TVector3 ax = p[i].Unit();
        for (int it = 0; it < 10; ++it) {
            TVector3 nw(0, 0, 0);
            for (size_t j = 0; j < p.size(); ++j) nw += (p[j].Dot(ax) > 0 ? p[j] : -p[j]);
            if (nw.Mag() < 1e-12) break;
            nw = nw.Unit();
            if ((nw - ax).Mag() < 1e-12) { ax = nw; break; }
            ax = nw;
        }
        double num = 0, den = 0;
        for (size_t j = 0; j < p.size(); ++j) { num += std::fabs(p[j].Dot(ax)); den += p[j].Mag(); }
        if (den > 0 && num / den > best) best = num / den;
    }
    return best;
}

}  // namespace

void check_aleph_thrust_definition(Long64_t nMatch = 3000)
{
    const int MAXN = 10000;
    int nPart;
    float thr;
    static float px[MAXN], py[MAXN], pz[MAXN];
    static int pdg[MAXN], st[MAXN];
    static bool isC[MAXN];

    // ---------------------------------------------------------------- part 1
    // Which definition reproduces the stored branch?
    TFile* f = TFile::Open(kIsrOn);
    TTree* t = static_cast<TTree*>(f->Get("tgenBefore"));
    t->SetBranchAddress("nParticle", &nPart);
    t->SetBranchAddress("thrust", &thr);
    t->SetBranchAddress("px", px);
    t->SetBranchAddress("py", py);
    t->SetBranchAddress("pz", pz);
    t->SetBranchAddress("pdgid", pdg);
    t->SetBranchAddress("status", st);
    t->SetBranchAddress("isCharged", isC);

    int matched[kNDef] = {0};
    double sumAbsDiff[kNDef] = {0};
    for (Long64_t i = 0; i < nMatch; ++i) {
        t->GetEntry(i);
        std::vector<TVector3> v[kNDef];
        for (int j = 0; j < nPart; ++j) {
            if (st[j] <= 0) continue;  // Pythia: positive status = final state
            const TVector3 p(px[j], py[j], pz[j]);
            const bool nu = isNeutrino(pdg[j]);
            const double mag = p.Mag();
            const double absCos = mag > 0 ? std::fabs(pz[j] / mag) : 0.0;
            v[kAllFinal].push_back(p);
            if (!nu) v[kNoNeutrino].push_back(p);
            if (isC[j]) v[kChargedOnly].push_back(p);
            if (!nu && absCos < 0.94) v[kNoNuAcc094].push_back(p);
            if (!nu && mag > 0.2) v[kNoNuP200].push_back(p);
        }
        for (int d = 0; d < kNDef; ++d) {
            const double x = thrustOf(v[d]);
            sumAbsDiff[d] += std::fabs(x - thr);
            if (std::fabs(x - thr) < 1e-5) matched[d]++;
        }
    }
    f->Close();

    printf("\n=== Which particle-level definition is the ALEPH `tgenBefore/thrust` branch? ===\n");
    printf("%-32s %12s %14s\n", "definition", "matched", "mean |diff|");
    for (int d = 0; d < kNDef; ++d)
        printf("%-32s %7d/%-4lld %14.2e\n", kDefName[d], matched[d], nMatch, sumAbsDiff[d] / nMatch);

    // ---------------------------------------------------------------- part 2
    // How much does the choice change the correction they apply?
    std::vector<double> ed;
    for (int i = 0; i <= 42; ++i) ed.push_back(0.58 + 0.01 * i);
    const int nb = static_cast<int>(ed.size()) - 1;
    TH1D* h[2][kNDef];
    for (int k = 0; k < 2; ++k)
        for (int d = 0; d < kNDef; ++d) {
            h[k][d] = new TH1D(Form("h_%d_%d", k, d), "", nb, &ed[0]);
            h[k][d]->SetDirectory(nullptr);
        }

    // Thrust is recomputed O(N^2) per definition, so this is the cost driver.
    // 150k events per state gives about 2% on the endpoint bin, which is ample
    // for a definition-sensitivity statement.
    const Long64_t nScan = 150000;
    for (int k = 0; k < 2; ++k) {
        TFile* fk = TFile::Open(k == 0 ? kIsrOff : kIsrOn);
        TTree* tk = static_cast<TTree*>(fk->Get("tgenBefore"));
        tk->SetBranchAddress("nParticle", &nPart);
        tk->SetBranchAddress("thrust", &thr);
        tk->SetBranchAddress("px", px);
        tk->SetBranchAddress("py", py);
        tk->SetBranchAddress("pz", pz);
        tk->SetBranchAddress("pdgid", pdg);
        tk->SetBranchAddress("status", st);
        tk->SetBranchAddress("isCharged", isC);
        const Long64_t n = std::min(nScan, tk->GetEntries());
        for (Long64_t i = 0; i < n; ++i) {
            tk->GetEntry(i);
            std::vector<TVector3> v[kNDef];
            for (int j = 0; j < nPart; ++j) {
                if (st[j] <= 0) continue;
                const TVector3 p(px[j], py[j], pz[j]);
                const bool nu = isNeutrino(pdg[j]);
                const double mag = p.Mag();
                const double absCos = mag > 0 ? std::fabs(pz[j] / mag) : 0.0;
                v[kAllFinal].push_back(p);
                if (!nu) v[kNoNeutrino].push_back(p);
                if (isC[j]) v[kChargedOnly].push_back(p);
                if (!nu && absCos < 0.94) v[kNoNuAcc094].push_back(p);
                if (!nu && mag > 0.2) v[kNoNuP200].push_back(p);
            }
            for (int d = 0; d < kNDef; ++d) h[k][d]->Fill(thrustOf(v[d]));
        }
        fk->Close();
        printf("[scan] %s done (%lld events)\n", k == 0 ? "ISR OFF" : "ISR ON ", n);
    }

    printf("\n=== C_ISR in the last ALEPH bin, 0.99 < T < 1.00, by definition ===\n");
    printf("(ALEPH samples, %lld events per state)\n", nScan);
    printf("%-32s %18s\n", "definition", "C_ISR last bin");
    for (int d = 0; d < kNDef; ++d) {
        const double a = h[0][d]->GetBinContent(nb), c = h[1][d]->GetBinContent(nb);
        const double r = c > 0 ? a / c : 0;
        const double e = (a > 0 && c > 0) ? r * std::sqrt(1. / a + 1. / c) : 0;
        printf("%-32s %9.4f +- %.4f\n", kDefName[d], r, e);
    }
}
