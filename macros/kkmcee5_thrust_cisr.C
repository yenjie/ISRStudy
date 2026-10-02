// Thrust ISR correction for the KKMCee 5.00.02 sample, per ALEPH bin and per
// stored thrust definition, from its EndpointDiagnostics trees.  Writes rows
// in the same layout as isr_model_cisr_thrust_bins.csv so the two can be
// compared directly.  The last bin of each definition is printed as a check.
//
// Usage:
//   root -l -b -q 'macros/kkmcee5_thrust_cisr.C("<diag dir>", "<out dir>")'

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {
struct Def { std::string branch, tag; };
const std::vector<Def> kDefs = {
    {"T_lab_including_ISR_photons", "B"}, {"T_lab_excluding_ISR_photons", "A"},
    {"T_lab_allFinal_including_ISR_photons", "N"}, {"T_lab_hadron_including_ISR_photons", "HB"},
    {"T_lab_hadron_excluding_ISR_photons", "HA"}, {"T_visibleCM_excluding_ISR_photons", "C"},
    {"T_hadronicCM_excluding_ISR_photons", "HC"}};

std::map<std::string, TH1D*> fillAll(const std::string& path, const std::string& tag)
{
    std::map<std::string, TH1D*> out;
    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) { std::cerr << "cannot open " << path << std::endl; return out; }
    TTree* t = static_cast<TTree*>(f->Get("EndpointDiagnostics"));
    t->SetBranchStatus("*", 0);
    double w = 1.0;
    t->SetBranchStatus("event_weight", 1);
    t->SetBranchAddress("event_weight", &w);
    std::vector<double> val(kDefs.size(), 0.0);
    for (size_t d = 0; d < kDefs.size(); ++d) {
        t->SetBranchStatus(kDefs[d].branch.c_str(), 1);
        t->SetBranchAddress(kDefs[d].branch.c_str(), &val[d]);
        TH1D* h = new TH1D((tag + "_" + kDefs[d].tag).c_str(), "", 42, 0.58, 1.00);
        h->SetDirectory(nullptr); h->Sumw2();
        out[kDefs[d].tag] = h;
    }
    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        for (size_t d = 0; d < kDefs.size(); ++d) out[kDefs[d].tag]->Fill(val[d], w);
    }
    for (auto& kv : out)
        if (kv.second->GetSumOfWeights() > 0) kv.second->Scale(static_cast<double>(n) / kv.second->GetSumOfWeights());
    std::cout << "  [fill] " << tag << " " << n << " events" << std::endl;
    f->Close();
    return out;
}
}  // namespace

void kkmcee5_thrust_cisr(const char* diagDir = "/raid5/data/yjlee/ISR/samples/kkmcee5_20261002/endpoint_diagnostics",
                         const char* outDir = "/raid5/data/yjlee/ISR/samples/kkmcee5_20261002/results")
{
    gSystem->mkdir(outDir, kTRUE);
    auto off = fillAll(std::string(diagDir) + "/endpoint_diagnostics_KKMCee50002_ISR_OFF.root", "off");
    auto on = fillAll(std::string(diagDir) + "/endpoint_diagnostics_KKMCee50002_ISR_ON.root", "on");
    if (off.empty() || on.empty()) return;

    // Seven ON events failed hadronization.  Put both spectra on the same
    // accepted-event exposure before forming the shape ratio.
    const double nOff = off.at("B")->GetEntries();
    const double nOn = on.at("B")->GetEntries();
    if (nOff <= 0 || nOn <= 0) return;
    for (auto& kv : on) kv.second->Scale(nOff / nOn);

    std::ofstream csv(std::string(outDir) + "/isr_model_cisr_thrust_bins_kkmcee5.csv");
    csv << std::setprecision(8);
    csv << "sample,definition,T_low,T_high,tau_low,tau_high,n_off,n_on,c_isr,c_isr_err\n";
    printf("\n=== KKMCee 5.00.02 + Pythia 8.316, C_ISR in 0.99 < T < 1.00 by definition ===\n");
    for (const auto& d : kDefs) {
        TH1D* ho = off[d.tag]; TH1D* hn = on[d.tag];
        for (int b = 1; b <= 42; ++b) {
            const double a = ho->GetBinContent(b), c = hn->GetBinContent(b);
            const double ea = ho->GetBinError(b), ec = hn->GetBinError(b);
            const double tl = 1.0 - ho->GetBinLowEdge(b + 1), th = 1.0 - ho->GetBinLowEdge(b);
            double r = 0, er = 0;
            if (a > 0 && c > 0) { r = a / c; er = r * std::sqrt((ea / a) * (ea / a) + (ec / c) * (ec / c)); }
            csv << "\"KKMCee 5.00.02 + Pythia 8.316\"," << d.tag << "," << ho->GetBinLowEdge(b) << ","
                << ho->GetBinLowEdge(b + 1) << "," << tl << "," << th << "," << a << "," << c << ","
                << r << "," << er << "\n";
            if (b == 42) printf("  %-3s  %.4f +- %.4f   (OFF %.0f, ON %.0f)\n", d.tag.c_str(), r, er, a, c);
        }
    }
    csv.close();
    std::cout << "[done] " << outDir << "/isr_model_cisr_thrust_bins_kkmcee5.csv" << std::endl;
}
