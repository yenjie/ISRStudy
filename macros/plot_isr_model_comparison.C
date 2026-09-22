// Cross-generator ISR model comparison, including KKMC.
//
// Pass 1 reads the Events trees and measures, generator-neutrally, the energy
// carried by beam-collinear final-state photons (|cos theta| > 0.9999) and the
// invariant mass of the final state after removing neutrinos and those photons.
// Neither quantity uses the isISRPhoton tag, so the comparison is immune to the
// tagging defects documented in docs/ISR_MODELLING_AND_KKMC.md.  For KKMC the
// generator's own ISR photon list is read as well, which gives the only
// available ground truth for what an "ISR photon" is.
//
// Pass 2 reads the endpoint-diagnostics trees and builds
//   C_ISR(T) = N_ISR_OFF(T) / N_ISR_ON(T)
// in the ALEPH thrust binning, using the tag-independent definition N
// (lab-frame thrust from all stable final-state particles).
//
// Usage:
//   root -l -b -q macros/plot_isr_model_comparison.C

#include <fstream>
#include <iostream>
#include <map>
#include <string>
#include <vector>

namespace {

const char* kRealDiag =
    "/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics_allfinal_3M";
const char* kRealNtup = "/data2/yjlee/ISRsample/real_3M_20260511";
const char* kKkmcDiag = "/data2/yjlee/ISRsample/kkmc_1M_20260921/endpoint_diagnostics";
const char* kKkmcNtup = "/data2/yjlee/ISRsample/kkmc_1M_20260921";

// Beam-collinear photon definition used as the generator-neutral ISR proxy.
const double kBeamCollinearCos = 0.9999;

struct Sample {
    std::string label;       // legend label
    std::string ntupleOff;   // Events tree, ISR OFF
    std::string ntupleOn;    // Events tree, ISR ON
    std::string diagOff;     // diagnostics tree, ISR OFF
    std::string diagOn;      // diagnostics tree, ISR ON
    int color;
    int marker;
    bool genuineIsrToggle;   // false for the Herwig QED-shower pair
};

struct PassOne {
    double nGamma = 0;       // all final photons per event
    double eGamma = 0;       // their total energy
    double eBeamCol = 0;     // energy in |cos theta| > kBeamCollinearCos photons
    double nBeamCol = 0;
    double eIsrTruth = -1;   // KKMC only: generator's own ISR photon energy
    double mVis = 0;         // mass excluding neutrinos and beam-collinear photons
    Long64_t entries = 0;
    TH1D* hBeamCol = nullptr;   // per-event beam-collinear photon energy
    TH1D* hIsrTruth = nullptr;  // KKMC only: per-event generator ISR photon energy
};

std::vector<double> alephThrustEdges()
{
    std::vector<double> edges;
    for (int i = 0; i <= 42; ++i) edges.push_back(0.58 + 0.01 * i);
    return edges;
}

bool isNeutrino(int pdg)
{
    const int a = std::abs(pdg);
    return a == 12 || a == 14 || a == 16;
}

PassOne scanEvents(const std::string& path, Long64_t maxEvents, const std::string& tag)
{
    PassOne r;
    r.hBeamCol = new TH1D(("hBeamCol_" + tag).c_str(), "", 120, 0.0, 12.0);
    r.hBeamCol->SetDirectory(nullptr);
    r.hIsrTruth = new TH1D(("hIsrTruth_" + tag).c_str(), "", 120, 0.0, 12.0);
    r.hIsrTruth->SetDirectory(nullptr);
    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "  [skip] cannot open " << path << std::endl;
        return r;
    }
    TTree* t = static_cast<TTree*>(f->Get("Events"));
    if (!t) {
        std::cerr << "  [skip] no Events tree in " << path << std::endl;
        f->Close();
        return r;
    }
    std::vector<int>* pdg = nullptr;
    std::vector<char>* isFinal = nullptr;
    std::vector<char>* isIsr = nullptr;
    std::vector<float>*px = nullptr, *py = nullptr, *pz = nullptr, *en = nullptr;
    t->SetBranchAddress("pdgId", &pdg);
    t->SetBranchAddress("isFinal", &isFinal);
    t->SetBranchAddress("isISRPhoton", &isIsr);
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("energy", &en);

    const bool isKkmc = path.find("KKMC") != std::string::npos;
    Long64_t n = t->GetEntries();
    if (maxEvents > 0 && n > maxEvents) n = maxEvents;
    double eIsrTruth = 0;

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        double E = 0, PX = 0, PY = 0, PZ = 0;
        double eBeamColEvt = 0, eIsrTruthEvt = 0;
        for (size_t j = 0; j < pdg->size(); ++j) {
            if (!(*isFinal)[j]) continue;
            const double p = std::sqrt((*px)[j] * (*px)[j] + (*py)[j] * (*py)[j] +
                                       (*pz)[j] * (*pz)[j]);
            const double absCos = p > 0 ? std::fabs((*pz)[j] / p) : 0.0;
            const bool beamCol = ((*pdg)[j] == 22 && absCos > kBeamCollinearCos);
            if ((*pdg)[j] == 22) {
                r.nGamma += 1;
                r.eGamma += (*en)[j];
            }
            if (beamCol) {
                r.nBeamCol += 1;
                r.eBeamCol += (*en)[j];
                eBeamColEvt += (*en)[j];
            }
            if (isKkmc && (*isIsr)[j]) {
                eIsrTruth += (*en)[j];
                eIsrTruthEvt += (*en)[j];
            }
            if (!isNeutrino((*pdg)[j]) && !beamCol) {
                E += (*en)[j];
                PX += (*px)[j];
                PY += (*py)[j];
                PZ += (*pz)[j];
            }
        }
        const double m2 = E * E - PX * PX - PY * PY - PZ * PZ;
        r.mVis += m2 > 0 ? std::sqrt(m2) : 0.0;
        r.hBeamCol->Fill(eBeamColEvt);
        if (isKkmc) r.hIsrTruth->Fill(eIsrTruthEvt);
    }
    r.entries = n;
    if (n > 0) {
        r.nGamma /= n;
        r.eGamma /= n;
        r.nBeamCol /= n;
        r.eBeamCol /= n;
        r.mVis /= n;
        if (isKkmc) r.eIsrTruth = eIsrTruth / n;
    }
    f->Close();
    return r;
}

TH1D* thrustHist(const std::string& path, const std::string& name,
                 const std::vector<double>& edges)
{
    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) {
        std::cerr << "  [skip] cannot open " << path << std::endl;
        return nullptr;
    }
    TTree* t = static_cast<TTree*>(f->Get("EndpointDiagnostics"));
    if (!t) {
        std::cerr << "  [skip] no EndpointDiagnostics tree in " << path << std::endl;
        f->Close();
        return nullptr;
    }
    // Fill by reading the branch directly.  TTree::Draw with a ">>name" target
    // resolves the histogram through gDirectory, which is fragile once the
    // histogram has been detached from a directory.
    TH1D* h = new TH1D(name.c_str(), "", static_cast<int>(edges.size()) - 1, &edges[0]);
    h->SetDirectory(nullptr);
    // Definition B: all stable final-state particles EXCLUDING neutrinos.  This
    // was verified to reproduce the `tgenBefore/thrust` branch of
    // ALEPH_Agentic_Event_Shape_Analysis to 1.5e-8 over 3000 events, so the
    // correction here is for the same observable that analysis corrects.
    double thrust = 0.0;
    t->SetBranchStatus("*", 0);
    t->SetBranchStatus("T_lab_including_ISR_photons", 1);
    t->SetBranchAddress("T_lab_including_ISR_photons", &thrust);
    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        h->Fill(thrust);
    }
    for (int b = 1; b <= h->GetNbinsX(); ++b) {
        h->SetBinError(b, std::sqrt(std::max(0.0, h->GetBinContent(b))));
    }
    std::cout << "  [thrust] " << name << " entries " << n << std::endl;
    f->Close();
    return h;
}

// C_ISR = OFF / ON with independent-Poisson errors.
TGraphErrors* ratioGraph(TH1D* off, TH1D* on, double xShift)
{
    if (!off || !on) return nullptr;
    std::vector<double> x, y, ex, ey;
    for (int b = 1; b <= off->GetNbinsX(); ++b) {
        const double nOff = off->GetBinContent(b);
        const double nOn = on->GetBinContent(b);
        if (nOff < 1 || nOn < 1) continue;
        const double r = nOff / nOn;
        x.push_back(off->GetBinCenter(b) + xShift);
        y.push_back(r);
        ex.push_back(0.0);
        ey.push_back(r * std::sqrt(1.0 / nOff + 1.0 / nOn));
    }
    if (x.empty()) return nullptr;
    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
    return g;
}

void styleGraph(TGraphErrors* g, int color, int marker)
{
    if (!g) return;
    g->SetLineColor(color);
    g->SetMarkerColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(1.0);
    g->SetLineWidth(2);
}

}  // namespace

void plot_isr_model_comparison(const char* outDir =
                                   "/data2/yjlee/ISRsample/kkmc_1M_20260921/results")
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gSystem->mkdir(outDir, kTRUE);

    std::vector<Sample> samples = {
        {"Pythia 8.315",
         std::string(kRealNtup) + "/mc_Pythia8315_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Pythia8315_ISR_ON.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Pythia8315_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Pythia8315.root",
         kBlack, 20, true},
        {"Pythia 8.315 Vincia",
         std::string(kRealNtup) + "/mc_Pythia8315_Vincia_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Pythia8315_Vincia_ISR_ON.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Pythia8315_Vincia_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Pythia8315_Vincia.root",
         kGray + 2, 24, true},
        {"Sherpa 3.0.3 PDFESherpa",
         std::string(kRealNtup) + "/mc_Sherpa303_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Sherpa303_ISR_ON.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Sherpa303_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Sherpa303.root",
         kAzure + 2, 21, true},
        {"Sherpa 3.0.3 YFS",
         std::string(kRealNtup) + "/mc_Sherpa303_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Sherpa303_ISR_YFS.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Sherpa303_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Sherpa303_YFS.root",
         kAzure + 7, 25, true},
        {"Herwig 7.3.0 QED shower, ISR unchanged",
         std::string(kRealNtup) + "/mc_Herwig730_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Herwig730_QEDshower.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Herwig730_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Herwig730_QEDshower.root",
         kOrange + 7, 22, false},
        {"KKMC 4.30 CEEX",
         std::string(kKkmcNtup) + "/mc_KKMC424_ISR_OFF.root",
         std::string(kKkmcNtup) + "/mc_KKMC424_ISR_ON.root",
         std::string(kKkmcDiag) + "/endpoint_diagnostics_KKMC424_ISR_OFF.root",
         std::string(kKkmcDiag) + "/endpoint_diagnostics_KKMC424_ISR_ON.root",
         kRed + 1, 29, true},
    };

    // ---------------------------------------------------------------------
    // Pass 1: generator-neutral radiation measurements from the Events trees.
    // ---------------------------------------------------------------------
    const Long64_t kScanEvents = 1000000;
    std::ofstream csv(std::string(outDir) + "/isr_model_summary.csv");
    csv << "sample,isr_state,entries_scanned,mean_n_gamma,mean_e_gamma,"
           "mean_n_beam_collinear_gamma,mean_e_beam_collinear_gamma,"
           "mean_e_isr_generator_truth,mean_m_vis\n";

    std::map<std::string, PassOne> offScan, onScan;
    for (size_t si = 0; si < samples.size(); ++si) {
        const Sample& s = samples[si];
        std::cout << "[pass1] " << s.label << std::endl;
        PassOne o = scanEvents(s.ntupleOff, kScanEvents, Form("off%d", static_cast<int>(si)));
        PassOne n = scanEvents(s.ntupleOn, kScanEvents, Form("on%d", static_cast<int>(si)));
        offScan[s.label] = o;
        onScan[s.label] = n;
        for (int k = 0; k < 2; ++k) {
            const PassOne& r = k == 0 ? o : n;
            csv << "\"" << s.label << "\"," << (k == 0 ? "OFF" : "ON") << ","
                << r.entries << "," << r.nGamma << "," << r.eGamma << ","
                << r.nBeamCol << "," << r.eBeamCol << "," << r.eIsrTruth << ","
                << r.mVis << "\n";
        }
    }
    csv.close();

    // Per-event ISR energy spectrum, shape normalized.
    TCanvas* c0 = new TCanvas("c0", "", 1000, 750);
    c0->SetTopMargin(0.09);
    c0->SetLeftMargin(0.13);
    c0->SetLogy();
    TH1D* f0 = new TH1D("frame0", "", 1, 0.0, 10.0);
    f0->SetMinimum(2e-5);
    f0->SetMaximum(5.0);
    f0->GetXaxis()->SetTitle("Per-event beam-collinear photon energy [GeV]");
    f0->GetYaxis()->SetTitle("Fraction of events / 0.1 GeV");
    f0->GetXaxis()->SetTitleSize(0.045);
    f0->GetYaxis()->SetTitleSize(0.045);
    f0->GetXaxis()->SetLabelSize(0.040);
    f0->GetYaxis()->SetLabelSize(0.040);
    f0->Draw();
    TLegend* leg0 = new TLegend(0.40, 0.58, 0.92, 0.90);
    leg0->SetBorderSize(0);
    leg0->SetFillStyle(0);
    leg0->SetTextSize(0.032);
    for (size_t i = 0; i < samples.size(); ++i) {
        TH1D* h = onScan[samples[i].label].hBeamCol;
        if (!h || h->Integral() <= 0) continue;
        h->Scale(1.0 / h->Integral());
        h->SetLineColor(samples[i].color);
        h->SetLineWidth(2);
        h->Draw("HIST SAME");
        leg0->AddEntry(h, samples[i].label.c_str(), "l");
    }
    TH1D* hTruth = onScan["KKMC 4.30 CEEX"].hIsrTruth;
    if (hTruth && hTruth->Integral() > 0) {
        hTruth->Scale(1.0 / hTruth->Integral());
        hTruth->SetLineColor(kRed + 1);
        hTruth->SetLineStyle(2);
        hTruth->SetLineWidth(2);
        hTruth->Draw("HIST SAME");
        leg0->AddEntry(hTruth, "KKMC generator ISR photon list", "l");
    }
    leg0->Draw();
    TLatex wip0;
    wip0.SetNDC();
    wip0.SetTextSize(0.032);
    wip0.SetTextColor(kGray + 2);
    wip0.DrawLatex(0.13, 0.935, "ISR study, work in progress");
    c0->SaveAs(Form("%s/isr_model_photon_energy.png", outDir));
    c0->SaveAs(Form("%s/isr_model_photon_energy.pdf", outDir));

    // ---------------------------------------------------------------------
    // Pass 2: C_ISR(T) in the ALEPH thrust binning, definition N.
    // ---------------------------------------------------------------------
    const std::vector<double> edges = alephThrustEdges();
    std::vector<TGraphErrors*> graphs;
    std::vector<TH1D*> offHists, onHists;
    for (size_t i = 0; i < samples.size(); ++i) {
        const Sample& s = samples[i];
        std::cout << "[pass2] " << s.label << std::endl;
        TH1D* off = thrustHist(s.diagOff, Form("h_off_%zu", i), edges);
        TH1D* on = thrustHist(s.diagOn, Form("h_on_%zu", i), edges);
        offHists.push_back(off);
        onHists.push_back(on);
        TGraphErrors* g = ratioGraph(off, on, 0.0015 * (static_cast<double>(i) - 2.5));
        styleGraph(g, s.color, s.marker);
        graphs.push_back(g);
    }

    TCanvas* c1 = new TCanvas("c1", "", 1000, 800);
    c1->SetTopMargin(0.09);
    c1->SetLeftMargin(0.13);
    c1->cd();
    // Bins below 0.70 are omitted: they carry very large statistical
    // uncertainties and compress the region the correction actually matters in.
    TH1D* frame = new TH1D("frame", "", 1, 0.70, 1.0);
    frame->SetMinimum(0.92);
    frame->SetMaximum(1.28);
    frame->GetXaxis()->SetTitle("Thrust T");
    frame->GetYaxis()->SetTitle("C_{ISR} = N_{ISR OFF} / N_{ISR ON}");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->Draw();
    TLine* unity = new TLine(0.70, 1.0, 1.0, 1.0);
    unity->SetLineStyle(2);
    unity->SetLineColor(kGray + 1);
    unity->Draw();
    TLegend* leg = new TLegend(0.17, 0.60, 0.66, 0.87);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    for (size_t i = 0; i < graphs.size(); ++i) {
        if (!graphs[i]) continue;
        graphs[i]->Draw("P SAME");
        leg->AddEntry(graphs[i], samples[i].label.c_str(), "lp");
    }
    // Reference: the ISR correction the ALEPH analysis itself uses, from its own
    // standalone Pythia8 samples and its own thrust branch.
    {
        TH1D* aOff = new TH1D("a_off", "", static_cast<int>(edges.size()) - 1, &edges[0]);
        TH1D* aOn = new TH1D("a_on", "", static_cast<int>(edges.size()) - 1, &edges[0]);
        aOff->SetDirectory(nullptr);
        aOn->SetDirectory(nullptr);
        const char* ap[2] = {
            "/raid5/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/Isr/isr0_ALL.root",
            "/raid5/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/Isr/isr1_ALL.root"};
        bool ok = true;
        for (int k = 0; k < 2; ++k) {
            TFile* af = TFile::Open(ap[k]);
            if (!af || af->IsZombie()) { ok = false; break; }
            TTree* at = static_cast<TTree*>(af->Get("tgenBefore"));
            if (!at) { ok = false; af->Close(); break; }
            float thr = 0;
            at->SetBranchStatus("*", 0);
            at->SetBranchStatus("thrust", 1);
            at->SetBranchAddress("thrust", &thr);
            const Long64_t an = at->GetEntries();
            for (Long64_t i = 0; i < an; ++i) { at->GetEntry(i); (k == 0 ? aOff : aOn)->Fill(thr); }
            af->Close();
        }
        if (ok) {
            TGraphErrors* ga = ratioGraph(aOff, aOn, 0.0);
            if (ga) {
                ga->SetLineColor(kGreen + 3);
                ga->SetMarkerColor(kGreen + 3);
                ga->SetMarkerStyle(34);
                ga->SetMarkerSize(1.2);
                ga->SetLineWidth(2);
                ga->Draw("P SAME");
                leg->AddEntry(ga, "ALEPH analysis Pythia8 (in use)", "lp");
            }
        }
    }
    leg->Draw();
    TLatex wip1;
    wip1.SetNDC();
    wip1.SetTextSize(0.032);
    wip1.SetTextColor(kGray + 2);
    wip1.DrawLatex(0.13, 0.935, "ISR study, work in progress");
    wip1.DrawLatex(0.62, 0.935, "stat. uncertainties only");
    c1->SaveAs(Form("%s/isr_model_cisr_thrust.png", outDir));
    c1->SaveAs(Form("%s/isr_model_cisr_thrust.pdf", outDir));

    // ---------------------------------------------------------------------
    // KKMC-only: effect of ISR-FSR interference on the same observable.
    // ---------------------------------------------------------------------
    TH1D* kkOn = thrustHist(std::string(kKkmcDiag) + "/endpoint_diagnostics_KKMC424_ISR_ON.root",
                            "h_kk_on", edges);
    TH1D* kkIfi = thrustHist(std::string(kKkmcDiag) + "/endpoint_diagnostics_KKMC424_ISR_ON_IFI.root",
                             "h_kk_ifi", edges);
    if (kkOn && kkIfi) {
        TGraphErrors* gIfi = ratioGraph(kkIfi, kkOn, 0.0);
        styleGraph(gIfi, kRed + 1, 29);
        TCanvas* c2 = new TCanvas("c2", "", 1000, 700);
        c2->SetTopMargin(0.10);
        c2->SetLeftMargin(0.13);
        TH1D* f2 = new TH1D("frame2", "", 1, 0.70, 1.0);
        f2->SetMinimum(0.95);
        f2->SetMaximum(1.05);
        f2->GetXaxis()->SetTitle("Thrust T");
        f2->GetYaxis()->SetTitle("N(interference on) / N(interference off)");
        f2->GetXaxis()->SetTitleSize(0.048);
        f2->GetYaxis()->SetTitleSize(0.044);
        f2->GetXaxis()->SetLabelSize(0.042);
        f2->GetYaxis()->SetLabelSize(0.042);
        f2->Draw();
        TLine* u2 = new TLine(0.70, 1.0, 1.0, 1.0);
        u2->SetLineStyle(2);
        u2->SetLineColor(kGray + 1);
        u2->Draw();
        if (gIfi) gIfi->Draw("P SAME");
        TLatex wip2;
        wip2.SetNDC();
        wip2.SetTextSize(0.036);
        wip2.SetTextColor(kGray + 2);
        wip2.DrawLatex(0.13, 0.935, "KKMC 4.30, ISR study, work in progress");
        wip2.DrawLatex(0.70, 0.935, "stat. only");
        c2->SaveAs(Form("%s/isr_model_kkmc_ifi.png", outDir));
        c2->SaveAs(Form("%s/isr_model_kkmc_ifi.pdf", outDir));

        // Integrated and endpoint interference effect.
        std::ofstream ifi(std::string(outDir) + "/isr_model_kkmc_ifi.csv");
        ifi << "quantity,value,error\n";
        const double totOn = kkOn->Integral();
        const double totIfi = kkIfi->Integral();
        const double rTot = totOn > 0 ? totIfi / totOn : 0.0;
        const double eTot = (totOn > 0 && totIfi > 0)
                                ? rTot * std::sqrt(1.0 / totOn + 1.0 / totIfi)
                                : 0.0;
        ifi << "integrated_ratio," << rTot << "," << eTot << "\n";
        const int nb = kkOn->GetNbinsX();
        const double lOn = kkOn->GetBinContent(nb);
        const double lIfi = kkIfi->GetBinContent(nb);
        const double rLast = lOn > 0 ? lIfi / lOn : 0.0;
        const double eLast = (lOn > 0 && lIfi > 0)
                                 ? rLast * std::sqrt(1.0 / lOn + 1.0 / lIfi)
                                 : 0.0;
        ifi << "lastbin_ratio," << rLast << "," << eLast << "\n";
        ifi << "mean_thrust_ifi_off," << kkOn->GetMean() << ",0\n";
        ifi << "mean_thrust_ifi_on," << kkIfi->GetMean() << ",0\n";
        ifi.close();
    }

    // ---------------------------------------------------------------------
    // Endpoint and mean-thrust table.
    // ---------------------------------------------------------------------
    std::ofstream tab(std::string(outDir) + "/isr_model_endpoint.csv");
    tab << "sample,genuine_isr_toggle,mean_T_off,mean_T_on,delta_mean_T,"
           "lastbin_off,lastbin_on,lastbin_cisr,lastbin_err\n";
    for (size_t i = 0; i < samples.size(); ++i) {
        TH1D* off = offHists[i];
        TH1D* on = onHists[i];
        if (!off || !on) continue;
        const int nb = off->GetNbinsX();
        const double nOff = off->GetBinContent(nb);
        const double nOn = on->GetBinContent(nb);
        const double r = (nOn > 0) ? nOff / nOn : 0.0;
        const double er = (nOff > 0 && nOn > 0)
                              ? r * std::sqrt(1.0 / nOff + 1.0 / nOn)
                              : 0.0;
        tab << "\"" << samples[i].label << "\","
            << (samples[i].genuineIsrToggle ? "yes" : "no") << ","
            << off->GetMean() << "," << on->GetMean() << ","
            << (on->GetMean() - off->GetMean()) << ","
            << nOff << "," << nOn << "," << r << "," << er << "\n";
    }
    tab.close();

    std::cout << "[done] wrote comparison outputs to " << outDir << std::endl;
}
