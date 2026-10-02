// Thrust C_ISR against tau = 1 - T for every sample, from the per-bin CSVs:
// the six samples of the main study, the new KKMCee 5.00.02 + Pythia 8 pair,
// the correction in use, and Anthony's values read off his plot.  Reads only
// CSVs, so it does not need the sample files.
//
// Usage:
//   root -l -b -q 'macros/plot_cisr_tau_kkmcee5.C("<main bins csv>","<kkmcee5 bins csv>","<anthony csv>","<aleph precomputed root>","<out dir>")'

#include <array>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {
std::vector<std::string> splitCsv(const std::string& line)
{
    std::vector<std::string> out; std::string cur; bool q = false;
    for (char c : line) {
        if (c == '"') { q = !q; continue; }
        if (c == ',' && !q) { out.push_back(cur); cur.clear(); continue; }
        cur.push_back(c);
    }
    out.push_back(cur); return out;
}
struct Pt { double tau, c, e; };
// per-bin rows of definition B -> sample -> points
std::map<std::string, std::vector<Pt>> readBins(const char* path)
{
    std::map<std::string, std::vector<Pt>> out;
    std::ifstream in(path); std::string line; std::getline(in, line);
    while (std::getline(in, line)) {
        auto f = splitCsv(line);
        if (f.size() < 10 || f[1] != "B") continue;
        const double tl = std::atof(f[4].c_str()), th = std::atof(f[5].c_str());
        const double c = std::atof(f[8].c_str()), e = std::atof(f[9].c_str());
        if (c > 0) out[f[0]].push_back({0.5 * (tl + th), c, e});
    }
    return out;
}
}  // namespace

void plot_cisr_tau_kkmcee5(
    const char* mainCsv = "/raid5/data/yjlee/ISR/overleaf/results/kkmc/isr_model_cisr_thrust_bins.csv",
    const char* kkmcee5Csv = "/raid5/data/yjlee/ISR/samples/kkmcee5_20261002/results/isr_model_cisr_thrust_bins_kkmcee5.csv",
    const char* anthonyCsv = "/raid5/data/yjlee/ISR/overleaf/results/kkmc/anthony_cisr_tau_digitized.csv",
    const char* alephRoot = "/raid5/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/Isr/isr_corr_precomputed.root",
    const char* outDir = "/raid5/data/yjlee/ISR/samples/kkmcee5_20261002/results")
{
    gStyle->SetOptStat(0); gStyle->SetPadTickX(1); gStyle->SetPadTickY(1);
    gStyle->SetTextFont(42); gStyle->SetLabelFont(42, "XYZ"); gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetEndErrorSize(0); gStyle->SetErrorX(0);

    auto bins = readBins(mainCsv);
    auto kk5 = readBins(kkmcee5Csv);
    for (auto& kv : kk5) bins[kv.first] = kv.second;

    struct Style { std::string label, legend; int color, marker; double shift; };
    const std::vector<Style> styles = {
        {"Pythia 8.315", "Pythia 8.315", kBlack, 20, -0.0030},
        {"Pythia 8.315 Vincia", "Pythia 8.315 Vincia", kGray + 2, 24, -0.0018},
        {"Sherpa 3.0.3 PDFESherpa", "Sherpa 3.0.3 PDFESherpa", kAzure + 2, 21, -0.0006},
        {"Sherpa 3.0.3 YFS", "Sherpa 3.0.3 YFS", kAzure + 7, 25, 0.0006},
        {"Herwig 7.3.0 QED shower, ISR unchanged", "Herwig 7.3.0 (pair does not switch ISR)", kOrange + 7, 22, 0.0018},
        {"KKMC 4.30 CEEX", "KKMC 4.30 CEEX + Pythia 6.202", kRed + 1, 29, 0.0030},
        {"KKMCee 5.00.02 + Pythia 8.316", "KKMCee 5.00.02 CEEX + Pythia 8.316", kMagenta + 2, 33, 0.0042},
    };
    std::map<std::string, TGraphErrors*> g;
    for (const auto& s : styles) {
        if (!bins.count(s.label)) continue;
        std::vector<double> x, y, ex, ey;
        for (const auto& p : bins[s.label]) { x.push_back(p.tau + s.shift); y.push_back(p.c); ex.push_back(0); ey.push_back(p.e); }
        TGraphErrors* gr = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
        gr->SetMarkerColor(s.color); gr->SetLineColor(s.color); gr->SetMarkerStyle(s.marker);
        gr->SetMarkerSize(s.marker == 33 ? 1.5 : 0.9); gr->SetLineWidth(1);
        g[s.label] = gr;
    }
    // Anthony
    std::map<std::string, TGraphErrors*> anthony;
    {
        std::ifstream in(anthonyCsv); std::string line; std::getline(in, line);
        std::map<std::string, std::vector<std::array<double, 3>>> pts;
        while (std::getline(in, line)) {
            auto f = splitCsv(line);
            if (f.size() < 4) continue;
            if (f[0] != "KKMCee 5.00.02" && f[0] != "collinear models") continue;
            pts[f[0]].push_back({std::atof(f[1].c_str()), std::atof(f[2].c_str()), std::atof(f[3].c_str())});
        }
        for (auto& kv : pts) {
            std::vector<double> x, y, ex, ey;
            for (auto& p : kv.second) { x.push_back(p[0]); y.push_back(p[1]); ex.push_back(0); ey.push_back(p[2]); }
            TGraphErrors* gr = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
            const bool kk = kv.first == "KKMCee 5.00.02";
            gr->SetMarkerColor(kk ? kRed + 1 : kGray + 1); gr->SetLineColor(kk ? kRed + 1 : kGray + 1);
            gr->SetMarkerStyle(kk ? 30 : 27); gr->SetMarkerSize(kk ? 2.6 : 2.4);
            anthony[kv.first] = gr;
        }
    }
    // in-use correction
    TGraph* ga = nullptr;
    if (TFile* fa = TFile::Open(alephRoot)) {
        TIter it(fa->GetListOfKeys()); TKey* k;
        while ((k = static_cast<TKey*>(it()))) {
            TObject* o = k->ReadObj(); if (!o->InheritsFrom("TH1")) continue;
            TH1* h = static_cast<TH1*>(o); std::vector<double> x, y;
            for (int b = 1; b <= h->GetNbinsX(); ++b) { if (h->GetBinLowEdge(b) < 0.58 - 1e-9) continue; x.push_back(1.0 - h->GetBinCenter(b)); y.push_back(h->GetBinContent(b)); }
            ga = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
            ga->SetMarkerStyle(34); ga->SetMarkerSize(1.2); ga->SetMarkerColor(kGreen + 2); ga->SetLineColor(kGreen + 2);
            break;
        }
    }

    TCanvas* c1 = new TCanvas("c_tau5", "", 1600, 720);
    TPad* pL = new TPad("pL", "", 0.00, 0.0, 0.50, 1.0);
    TPad* pR = new TPad("pR", "", 0.50, 0.0, 1.00, 1.0);
    for (TPad* p : {pL, pR}) { p->SetLeftMargin(0.14); p->SetRightMargin(0.03); p->SetTopMargin(0.08); p->SetBottomMargin(0.13); p->Draw(); }
    auto panel = [&](TPad* pad, double xmax, double ymin, double ymax, const char* title, bool legend) {
        pad->cd();
        TH1D* fr = new TH1D(Form("fr_%s", pad->GetName()), "", 1, 0.0, xmax);
        fr->SetMinimum(ymin); fr->SetMaximum(ymax);
        fr->GetXaxis()->SetTitle("#tau = 1 - T"); fr->GetYaxis()->SetTitle("C_{ISR} = N_{ISR OFF} / N_{ISR ON}");
        fr->GetXaxis()->SetTitleSize(0.052); fr->GetYaxis()->SetTitleSize(0.052);
        fr->GetXaxis()->SetLabelSize(0.044); fr->GetYaxis()->SetLabelSize(0.044); fr->GetYaxis()->SetTitleOffset(1.25);
        fr->Draw();
        TLine* u = new TLine(0.0, 1.0, xmax, 1.0); u->SetLineStyle(2); u->SetLineColor(kGray + 1); u->Draw();
        for (const auto& s : styles) if (g.count(s.label)) g[s.label]->Draw("P SAME");
        if (ga) ga->Draw("P SAME");
        for (auto& kv : anthony) kv.second->Draw("P SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.040); tx.DrawLatex(0.14, 0.94, title);
        if (!legend) return;
        TLegend* leg = new TLegend(0.28, 0.42, 0.97, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.032);
        for (const auto& s : styles) if (g.count(s.label)) leg->AddEntry(g[s.label], s.legend.c_str(), "p");
        if (ga) leg->AddEntry(ga, "ALEPH analysis Pythia8, in use", "p");
        if (anthony.count("KKMCee 5.00.02")) leg->AddEntry(anthony["KKMCee 5.00.02"], "Anthony: KKMCee 5.00.02", "p");
        if (anthony.count("collinear models")) leg->AddEntry(anthony["collinear models"], "Anthony: mean of his 4 non-KKMC models", "p");
        leg->Draw();
    };
    panel(pL, 0.42, 0.75, 1.25, "Full range, Anthony's axes", false);
    panel(pR, 0.12, 0.95, 1.20, "Zoom on the bins where the correction varies", true);
    c1->cd();
    TLatex tf; tf.SetNDC(); tf.SetTextFont(42); tf.SetTextSize(0.022); tf.SetTextColor(kGray + 2); tf.SetTextAlign(31);
    tf.DrawLatex(0.985, 0.972, "ISR study, work in progress.  ALEPH thrust definition, ALEPH bins, stat. bars.  Anthony's values read off his plot, #pm0.004");
    c1->SaveAs(Form("%s/isr_model_cisr_tau_with_kkmcee5.png", outDir));
    c1->SaveAs(Form("%s/isr_model_cisr_tau_with_kkmcee5.pdf", outDir));

    TCanvas* cT = new TCanvas("c_thrust5", "", 1000, 760);
    cT->SetLeftMargin(0.13); cT->SetRightMargin(0.04);
    cT->SetTopMargin(0.09); cT->SetBottomMargin(0.13); cT->SetTicks(1, 1);
    TH1D* ft = new TH1D("fr_thrust5", "", 1, 0.70, 1.00);
    ft->SetMinimum(0.92); ft->SetMaximum(1.27);
    ft->GetXaxis()->SetTitle("Thrust T");
    ft->GetYaxis()->SetTitle("C_{ISR} = N_{ISR OFF} / N_{ISR ON}");
    ft->GetXaxis()->SetTitleSize(0.050); ft->GetYaxis()->SetTitleSize(0.050);
    ft->GetXaxis()->SetLabelSize(0.043); ft->GetYaxis()->SetLabelSize(0.043);
    ft->GetYaxis()->SetTitleOffset(1.25);
    ft->Draw();
    TLine* oneT = new TLine(0.70, 1.0, 1.0, 1.0);
    oneT->SetLineStyle(2); oneT->SetLineColor(kGray + 1); oneT->Draw();
    TLegend* lt = new TLegend(0.17, 0.50, 0.78, 0.88);
    lt->SetBorderSize(0); lt->SetFillStyle(0); lt->SetTextSize(0.029);
    for (const auto& s : styles) {
        if (!g.count(s.label)) continue;
        std::vector<double> x, y, ex, ey;
        TGraphErrors* src = g[s.label];
        for (int i = 0; i < src->GetN(); ++i) {
            double tau, c; src->GetPoint(i, tau, c);
            const double thrust = 1.0 - tau;
            if (thrust < 0.70 || thrust > 1.00) continue;
            x.push_back(thrust); y.push_back(c);
            ex.push_back(0.0); ey.push_back(src->GetErrorY(i));
        }
        if (x.empty()) continue;
        TGraphErrors* gr = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
        gr->SetLineColor(s.color); gr->SetMarkerColor(s.color);
        gr->SetMarkerStyle(s.marker); gr->SetMarkerSize(s.marker == 33 ? 1.2 : 0.75);
        gr->SetLineWidth(1); gr->Draw("PL SAME");
        lt->AddEntry(gr, s.legend.c_str(), "lp");
    }
    if (ga) {
        std::vector<double> x, y;
        for (int i = 0; i < ga->GetN(); ++i) {
            double tau, c; ga->GetPoint(i, tau, c);
            if (1.0 - tau < 0.70) continue;
            x.push_back(1.0 - tau); y.push_back(c);
        }
        if (!x.empty()) {
            TGraph* gAleph = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
            gAleph->SetLineColor(kGreen + 2); gAleph->SetMarkerColor(kGreen + 2);
            gAleph->SetMarkerStyle(34); gAleph->SetMarkerSize(1.1);
            gAleph->SetLineWidth(1); gAleph->Draw("PL SAME");
            lt->AddEntry(gAleph, "ALEPH analysis Pythia8 (in use)", "lp");
        }
    }
    lt->Draw();
    TLatex tT; tT.SetNDC(); tT.SetTextFont(42); tT.SetTextSize(0.030);
    tT.DrawLatex(0.13, 0.955, "ALEPH visible-particle thrust, #sqrt{s} = 91.1876 GeV");
    cT->SaveAs(Form("%s/isr_model_cisr_thrust_with_kkmcee5.png", outDir));
    cT->SaveAs(Form("%s/isr_model_cisr_thrust_with_kkmcee5.pdf", outDir));
    std::cout << "[done] thrust and tau figures in " << outDir << std::endl;
}
