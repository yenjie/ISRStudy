// Pass-1 photon summary for the KKMCee 5 pair, in the layout of
// isr_model_summary.csv, so the EEC verification table and figure can show the
// beam-collinear estimator for the new sample too.
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>
#include <vector>
void kkmcee5_photon_summary(const char* base = "/raid5/data/yjlee/ISR/samples/kkmcee5_20261002",
                            const char* outCsv = "/raid5/data/yjlee/ISR/samples/kkmcee5_20261002/results/isr_model_summary_kkmcee5.csv")
{
    std::ofstream csv(outCsv);
    csv << std::setprecision(6);
    csv << "sample,isr_state,entries_scanned,mean_n_gamma,mean_e_gamma,mean_n_beam_collinear_gamma,"
           "mean_e_beam_collinear_gamma,mean_e_isr_generator_truth,mean_m_vis\n";
    for (const char* st : {"OFF", "ON"}) {
        TFile* f = TFile::Open(Form("%s/mc_KKMCee50002_ISR_%s.root", base, st));
        TTree* t = (TTree*)f->Get("Events");
        std::vector<int>* pdg = nullptr; std::vector<char>* isFinal = nullptr; std::vector<char>* isIsr = nullptr;
        std::vector<float>*px = nullptr, *py = nullptr, *pz = nullptr, *en = nullptr;
        double w = 1;
        t->SetBranchStatus("*", 0);
        for (const char* b : {"pdgId", "isFinal", "isISRPhoton", "px", "py", "pz", "energy", "weight"}) t->SetBranchStatus(b, 1);
        t->SetBranchAddress("pdgId", &pdg); t->SetBranchAddress("isFinal", &isFinal); t->SetBranchAddress("isISRPhoton", &isIsr);
        t->SetBranchAddress("px", &px); t->SetBranchAddress("py", &py); t->SetBranchAddress("pz", &pz); t->SetBranchAddress("energy", &en);
        t->SetBranchAddress("weight", &w);
        double sumW = 0, nG = 0, eG = 0, nBC = 0, eBC = 0, eIsr = 0, mVis = 0;
        const Long64_t n = t->GetEntries();
        for (Long64_t i = 0; i < n; ++i) {
            t->GetEntry(i);
            double E = 0, PX = 0, PY = 0, PZ = 0;
            for (size_t j = 0; j < pdg->size(); ++j) {
                if (!(*isFinal)[j]) continue;
                const int a = std::abs((*pdg)[j]);
                const double p = std::sqrt((*px)[j]*(*px)[j] + (*py)[j]*(*py)[j] + (*pz)[j]*(*pz)[j]);
                const bool bc = ((*pdg)[j] == 22 && p > 0 && std::fabs((*pz)[j] / p) > 0.9999);
                if ((*pdg)[j] == 22) { nG += w; eG += w * (*en)[j]; }
                if (bc) { nBC += w; eBC += w * (*en)[j]; }
                if ((*isIsr)[j]) eIsr += w * (*en)[j];
                if (!(a == 12 || a == 14 || a == 16) && !bc) { E += (*en)[j]; PX += (*px)[j]; PY += (*py)[j]; PZ += (*pz)[j]; }
            }
            const double m2 = E*E - PX*PX - PY*PY - PZ*PZ;
            mVis += w * (m2 > 0 ? std::sqrt(m2) : 0);
            sumW += w;
        }
        csv << "\"KKMCee 5.00.02 + Pythia 8.316\"," << st << "," << n << "," << nG/sumW << "," << eG/sumW << ","
            << nBC/sumW << "," << eBC/sumW << "," << eIsr/sumW << "," << mVis/sumW << "\n";
        printf("%s: n_gamma %.3f  E_gamma %.3f  n_bc %.4f  E_bc %.4f  E_isr_truth %.4f  M_vis %.3f\n",
               st, nG/sumW, eG/sumW, nBC/sumW, eBC/sumW, eIsr/sumW, mVis/sumW);
        f->Close();
    }
    csv.close();
    printf("[done] %s\n", outCsv);
}
