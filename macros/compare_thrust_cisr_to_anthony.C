// Our thrust C_ISR per ALEPH bin, drawn against tau = 1 - T in the conventions
// of Anthony's plot (x from 0 to 0.42, y from 0.75 to 1.25), with his values
// overlaid when a digitised CSV is supplied.
//
// Also writes every bin for every available thrust definition, so that the
// sensitivity of the correction to how ISR photons enter the thrust can be
// read off directly.  That matters for the comparison: the one model on which
// the two studies disagree is KKMC, and the two studies need not treat the
// generator's ISR photons the same way.
//
// Inputs are the EndpointDiagnostics trees, one entry per event with the
// thrust under each definition already computed.
//
// Usage:
//   root -l -b -q 'macros/compare_thrust_cisr_to_anthony.C("<out>", "<anthony csv or empty>")'

#include <array>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

const char* kRealDiag =
    "/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics_allfinal_3M";
const char* kKkmcDiag = "/data2/yjlee/ISRsample/kkmc_1M_20260921/endpoint_diagnostics";
const char* kAlephPrecomputed =
    "/raid5/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/Isr/isr_corr_precomputed.root";

struct Sample {
    std::string label, off, on;
    int color, marker;
    double xShift;   // in tau, for visibility only
};

// The thrust definitions stored in the diagnostics tree.  B is the ALEPH
// analysis definition and the nominal one here.
struct Def { std::string branch, tag, text; };
const std::vector<Def> kDefs = {
    {"T_lab_including_ISR_photons",          "B", "lab, final state no neutrinos, ISR photons in  (ALEPH)"},
    {"T_lab_excluding_ISR_photons",          "A", "lab, final state no neutrinos, ISR photons out"},
    {"T_lab_allFinal_including_ISR_photons", "N", "lab, all final state incl. neutrinos, ISR photons in"},
    {"T_lab_hadron_including_ISR_photons",   "HB", "lab, hadrons only + ISR photons"},
    {"T_lab_hadron_excluding_ISR_photons",   "HA", "lab, hadrons only"},
    {"T_visibleCM_excluding_ISR_photons",    "C", "visible rest frame, ISR photons out"},
    {"T_hadronicCM_excluding_ISR_photons",   "HC", "hadronic rest frame, ISR photons out"},
};

std::vector<double> alephEdges()
{
    std::vector<double> e;
    for (int i = 0; i <= 42; ++i) e.push_back(0.58 + 0.01 * i);
    return e;
}

// One pass over a file fills one histogram per definition.
std::map<std::string, TH1D*> fillAll(const std::string& path, const std::string& tag,
                                     const std::vector<double>& edges)
{
    std::map<std::string, TH1D*> out;
    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) { std::cerr << "  [skip] " << path << std::endl; return out; }
    TTree* t = static_cast<TTree*>(f->Get("EndpointDiagnostics"));
    if (!t) { std::cerr << "  [skip] no tree " << path << std::endl; f->Close(); return out; }
    t->SetBranchStatus("*", 0);
    std::vector<double> val(kDefs.size(), 0.0);
    std::vector<bool> have(kDefs.size(), false);
    for (size_t d = 0; d < kDefs.size(); ++d) {
        if (!t->GetBranch(kDefs[d].branch.c_str())) continue;
        have[d] = true;
        t->SetBranchStatus(kDefs[d].branch.c_str(), 1);
        t->SetBranchAddress(kDefs[d].branch.c_str(), &val[d]);
        TH1D* h = new TH1D((tag + "_" + kDefs[d].tag).c_str(), "",
                           static_cast<int>(edges.size()) - 1, &edges[0]);
        h->SetDirectory(nullptr);
        out[kDefs[d].tag] = h;
    }
    const Long64_t n = t->GetEntries();
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        for (size_t d = 0; d < kDefs.size(); ++d)
            if (have[d]) out[kDefs[d].tag]->Fill(val[d]);
    }
    for (auto& kv : out)
        for (int b = 1; b <= kv.second->GetNbinsX(); ++b)
            kv.second->SetBinError(b, std::sqrt(std::max(0.0, kv.second->GetBinContent(b))));
    std::cout << "  [fill] " << tag << "  " << n << " events" << std::endl;
    f->Close();
    return out;
}

}  // namespace

void compare_thrust_cisr_to_anthony(
    const char* outDir = "/data2/yjlee/ISRsample/kkmc_1M_20260921/results",
    const char* anthonyCsv = "")
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetEndErrorSize(0);
    gStyle->SetErrorX(0);

    const std::vector<Sample> samples = {
        {"Pythia 8.315", std::string(kRealDiag) + "/endpoint_diagnostics_Pythia8315_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Pythia8315.root", kBlack, 20, -0.0030},
        {"Pythia 8.315 Vincia", std::string(kRealDiag) + "/endpoint_diagnostics_Pythia8315_Vincia_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Pythia8315_Vincia.root", kGray + 2, 24, -0.0018},
        {"Sherpa 3.0.3 PDFESherpa", std::string(kRealDiag) + "/endpoint_diagnostics_Sherpa303_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Sherpa303.root", kAzure + 2, 21, -0.0006},
        {"Sherpa 3.0.3 YFS", std::string(kRealDiag) + "/endpoint_diagnostics_Sherpa303_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Sherpa303_YFS.root", kAzure + 7, 25, 0.0006},
        {"Herwig 7.3.0 QED shower, ISR unchanged", std::string(kRealDiag) + "/endpoint_diagnostics_Herwig730_OFF.root",
         std::string(kRealDiag) + "/endpoint_diagnostics_Herwig730_QEDshower.root", kOrange + 7, 22, 0.0018},
        {"KKMC 4.30 CEEX", std::string(kKkmcDiag) + "/endpoint_diagnostics_KKMC424_ISR_OFF.root",
         std::string(kKkmcDiag) + "/endpoint_diagnostics_KKMC424_ISR_ON.root", kRed + 1, 29, 0.0030},
    };
    const std::vector<double> edges = alephEdges();
    const int nb = static_cast<int>(edges.size()) - 1;

    // Fill each distinct file once.
    std::map<std::string, std::map<std::string, TH1D*>> byFile;
    for (const auto& s : samples)
        for (const std::string& p : {s.off, s.on})
            if (!byFile.count(p)) byFile[p] = fillAll(p, Form("f%zu", byFile.size()), edges);

    // ------------------------------------------------------------ per-bin CSV
    std::ofstream csv(std::string(outDir) + "/isr_model_cisr_thrust_bins.csv");
    csv << std::setprecision(8);
    csv << "sample,definition,T_low,T_high,tau_low,tau_high,n_off,n_on,c_isr,c_isr_err\n";
    std::map<std::string, std::map<std::string, TGraphErrors*>> graphs;  // [def][sample]
    for (const auto& s : samples) {
        for (const auto& d : kDefs) {
            auto& fo = byFile[s.off];
            auto& fn = byFile[s.on];
            if (!fo.count(d.tag) || !fn.count(d.tag)) continue;
            TH1D* off = fo[d.tag];
            TH1D* on = fn[d.tag];
            std::vector<double> x, y, ex, ey;
            for (int b = 1; b <= nb; ++b) {
                const double a = off->GetBinContent(b), c = on->GetBinContent(b);
                const double tl = 1.0 - edges[b], th = 1.0 - edges[b - 1];
                double r = 0, er = 0;
                if (a > 0 && c > 0) {
                    r = a / c;
                    er = r * std::sqrt(1.0 / a + 1.0 / c);
                    x.push_back(0.5 * (tl + th) + s.xShift);
                    y.push_back(r); ex.push_back(0.0); ey.push_back(er);
                }
                csv << "\"" << s.label << "\"," << d.tag << "," << edges[b - 1] << "," << edges[b]
                    << "," << tl << "," << th << "," << a << "," << c << "," << r << "," << er << "\n";
            }
            if (!x.empty()) {
                TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
                g->SetMarkerColor(s.color); g->SetLineColor(s.color);
                g->SetMarkerStyle(s.marker); g->SetMarkerSize(0.9); g->SetLineWidth(1);
                graphs[d.tag][s.label] = g;
            }
        }
    }
    csv.close();

    // ------------------------------------------------ last-bin definition table
    std::ofstream tab(std::string(outDir) + "/isr_model_cisr_lastbin_by_definition.csv");
    tab << std::setprecision(6);
    tab << "sample,definition,description,n_off,n_on,c_isr,c_isr_err\n";
    printf("\n=== C_ISR in 0.99 < T < 1.00 by thrust definition ===\n");
    printf("%-40s", "sample");
    for (const auto& d : kDefs) printf(" %10s", d.tag.c_str());
    printf("\n");
    for (const auto& s : samples) {
        printf("%-40s", s.label.substr(0, 40).c_str());
        for (const auto& d : kDefs) {
            auto& fo = byFile[s.off];
            auto& fn = byFile[s.on];
            if (!fo.count(d.tag) || !fn.count(d.tag)) { printf(" %10s", "-"); continue; }
            const double a = fo[d.tag]->GetBinContent(nb), c = fn[d.tag]->GetBinContent(nb);
            const double r = c > 0 ? a / c : 0, er = (a > 0 && c > 0) ? r * std::sqrt(1 / a + 1 / c) : 0;
            printf(" %10.3f", r);
            tab << "\"" << s.label << "\"," << d.tag << ",\"" << d.text << "\"," << a << "," << c
                << "," << r << "," << er << "\n";
        }
        printf("\n");
    }
    tab.close();

    // ------------------------------------------------------------------ plot
    // Two panels: the full tau range in Anthony's conventions, and a zoom on
    // tau < 0.12 where the correction actually varies.  The legend lives in the
    // zoom panel's empty upper right; in the full-range panel the large error
    // bars beyond tau = 0.3 reach the top of the frame and would run through it.
    TCanvas* c1 = new TCanvas("c_tau", "", 1600, 720);
    TPad* pL = new TPad("pL", "", 0.00, 0.0, 0.50, 1.0);
    TPad* pR = new TPad("pR", "", 0.50, 0.0, 1.00, 1.0);
    for (TPad* p : {pL, pR}) {
        p->SetLeftMargin(0.14); p->SetRightMargin(0.03);
        p->SetTopMargin(0.08);  p->SetBottomMargin(0.13);
        p->Draw();
    }

    // Anthony's values, read off his plot.  Only the two series that matter
    // for the comparison are drawn: KKMCee, and the mean of his four collinear
    // models per bin.  Columns: series,tau,c_isr,reading_precision.
    std::map<std::string, TGraphErrors*> anthony;
    if (anthonyCsv && *anthonyCsv) {
        std::ifstream in(anthonyCsv);
        std::map<std::string, std::vector<std::array<double, 3>>> pts;
        std::string line;
        std::getline(in, line);
        while (std::getline(in, line)) {
            std::stringstream ss(line);
            std::string series, a, b, c;
            std::getline(ss, series, ','); std::getline(ss, a, ','); std::getline(ss, b, ','); std::getline(ss, c, ',');
            if (series != "KKMCee 5.00.02" && series != "collinear models") continue;
            pts[series].push_back({std::atof(a.c_str()), std::atof(b.c_str()), std::atof(c.c_str())});
        }
        for (auto& kv : pts) {
            std::vector<double> x, y, ex, ey;
            for (auto& p : kv.second) { x.push_back(p[0]); y.push_back(p[1]); ex.push_back(0); ey.push_back(p[2]); }
            TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
            const bool kk = kv.first == "KKMCee 5.00.02";
            g->SetMarkerColor(kk ? kRed + 1 : kGray + 1); g->SetLineColor(kk ? kRed + 1 : kGray + 1);
            g->SetMarkerStyle(kk ? 30 : 27); g->SetMarkerSize(kk ? 2.6 : 2.4); g->SetLineWidth(1);
            anthony[kv.first] = g;
        }
    }

    // The correction the ALEPH analysis applies today.  The cache file stores
    // no uncertainties; our own recount from its source files gives +-0.011 in
    // the last bin.
    TGraph* ga = nullptr;
    TFile* fa = TFile::Open(kAlephPrecomputed);
    if (fa && !fa->IsZombie()) {
        TIter it(fa->GetListOfKeys());
        TKey* k;
        while ((k = static_cast<TKey*>(it()))) {
            TObject* o = k->ReadObj();
            if (!o->InheritsFrom("TH1")) continue;
            TH1* h = static_cast<TH1*>(o);
            std::vector<double> x, y;
            for (int b = 1; b <= h->GetNbinsX(); ++b) {
                if (h->GetBinLowEdge(b) < 0.58 - 1e-9) continue;
                x.push_back(1.0 - h->GetBinCenter(b));
                y.push_back(h->GetBinContent(b));
            }
            ga = new TGraph(static_cast<int>(x.size()), &x[0], &y[0]);
            ga->SetMarkerStyle(34); ga->SetMarkerSize(1.2);
            ga->SetMarkerColor(kGreen + 2); ga->SetLineColor(kGreen + 2);
            break;
        }
    }

    auto drawPanel = [&](TPad* pad, double xmax, double ymin, double ymax, const char* title, bool withLegend) {
        pad->cd();
        TH1D* frame = new TH1D(Form("frame_%s", pad->GetName()), "", 1, 0.0, xmax);
        frame->SetMinimum(ymin); frame->SetMaximum(ymax);
        frame->GetXaxis()->SetTitle("#tau = 1 - T");
        frame->GetYaxis()->SetTitle("C_{ISR} = N_{ISR OFF} / N_{ISR ON}");
        frame->GetXaxis()->SetTitleSize(0.052); frame->GetYaxis()->SetTitleSize(0.052);
        frame->GetXaxis()->SetLabelSize(0.044); frame->GetYaxis()->SetLabelSize(0.044);
        frame->GetYaxis()->SetTitleOffset(1.25);
        frame->Draw();
        TLine* unity = new TLine(0.0, 1.0, xmax, 1.0);
        unity->SetLineStyle(2); unity->SetLineColor(kGray + 1); unity->Draw();
        for (const auto& s : samples)
            if (graphs["B"].count(s.label)) graphs["B"][s.label]->Draw("P SAME");
        if (ga) ga->Draw("P SAME");
        for (auto& kv : anthony) kv.second->Draw("P SAME");
        TLatex tx; tx.SetNDC(); tx.SetTextFont(42); tx.SetTextSize(0.040);
        tx.DrawLatex(0.14, 0.94, title);
        if (!withLegend) return;
        TLegend* leg = new TLegend(0.30, 0.44, 0.97, 0.88);
        leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextSize(0.033);
        for (const auto& s : samples) {
            if (!graphs["B"].count(s.label)) continue;
            leg->AddEntry(graphs["B"][s.label],
                          s.label == "Herwig 7.3.0 QED shower, ISR unchanged"
                              ? "Herwig 7.3.0 (pair does not switch ISR)" : s.label.c_str(), "p");
        }
        if (ga) leg->AddEntry(ga, "ALEPH analysis Pythia8, in use", "p");
        if (anthony.count("KKMCee 5.00.02"))
            leg->AddEntry(anthony["KKMCee 5.00.02"], "Anthony: KKMCee 5.00.02", "p");
        if (anthony.count("collinear models"))
            leg->AddEntry(anthony["collinear models"], "Anthony: mean of his 4 collinear models", "p");
        leg->Draw();
    };
    drawPanel(pL, 0.42, 0.75, 1.25, "Full range, Anthony's axes", false);
    drawPanel(pR, 0.12, 0.95, 1.20, "Zoom on the bins where the correction varies", true);
    c1->cd();
    // Canvas-level note above both frames, clear of the pad headers.
    TLatex tf; tf.SetNDC(); tf.SetTextFont(42); tf.SetTextSize(0.022); tf.SetTextColor(kGray + 2);
    tf.SetTextAlign(31);
    tf.DrawLatex(0.985, 0.972, "ISR study, work in progress.  Ours: ALEPH thrust definition, ALEPH bins, stat. bars.  "
                               "Anthony's values read off his plot, #pm0.004");

    c1->SaveAs(Form("%s/isr_model_cisr_tau_vs_anthony.png", outDir));
    c1->SaveAs(Form("%s/isr_model_cisr_tau_vs_anthony.pdf", outDir));
    std::cout << "[done] wrote per-bin table, definition table and tau plot to " << outDir << std::endl;
}
