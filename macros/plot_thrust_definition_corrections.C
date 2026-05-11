#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct ThrustDef {
    std::string tag;
    std::string branch;
    std::string label;
    std::string outputStem;
};

struct CorrSample {
    std::string label;
    std::string tag;
    std::string offFile;
    std::string onFile;
    int color;
    int marker;
    double xOffset;
};

std::string gCorrRealDir = "/data2/yjlee/ISRsample/real_3M_20260511";
std::string gCorrDiagDir = "/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics";
std::string gCorrOutDir = "/data2/yjlee/ISRsample/real_3M_20260511/results/thrust_definition_corrections";

std::string corrJoin(const std::string& dir, const std::string& name)
{
    return dir + "/" + name;
}

std::string corrAlephThrustCsv()
{
    std::string alephCsv = "/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/ALEPH/HEPData-ins636645-v1-Table_54.csv";
    if (const char* envAleph = gSystem->Getenv("ALEPH_THRUST_CSV")) alephCsv = envAleph;
    return alephCsv;
}

std::vector<double> corrAlephThrustEdges()
{
    std::vector<double> edges;
    std::ifstream in(corrAlephThrustCsv().c_str());
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.find("THRUST,THRUST LOW") != std::string::npos) continue;
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> cols;
        while (std::getline(ss, item, ',')) cols.push_back(item);
        if (cols.size() < 3) continue;
        const double low = std::stod(cols[1]);
        const double high = std::stod(cols[2]);
        if (edges.empty()) edges.push_back(low);
        edges.push_back(high);
    }
    if (edges.size() > 1) return edges;
    for (int i = 0; i <= 42; ++i) edges.push_back(0.58 + 0.01 * i);
    return edges;
}

TH1D* makeDiagHist(const std::string& filename,
                   const std::string& name,
                   const std::string& branch,
                   const std::vector<double>& edges)
{
    TFile f(filename.c_str(), "READ");
    if (f.IsZombie()) {
        std::cerr << "Cannot open " << filename << std::endl;
        return nullptr;
    }
    TTree* t = static_cast<TTree*>(f.Get("EndpointDiagnostics"));
    if (!t) {
        std::cerr << "Missing EndpointDiagnostics tree in " << filename << std::endl;
        return nullptr;
    }
    gROOT->cd();
    TH1D* h = new TH1D(name.c_str(), "", static_cast<int>(edges.size()) - 1, edges.data());
    h->Sumw2();
    h->SetDirectory(gROOT);
    t->Draw((branch + ">>" + name).c_str(), "", "goff");
    h->SetDirectory(nullptr);
    return h;
}

TGraphErrors* makeRatioGraph(const CorrSample& sample,
                             const ThrustDef& def,
                             const std::vector<double>& edges,
                             std::ofstream* lastBinSummary,
                             std::ofstream* allBinSummary)
{
    TH1D* hOff = makeDiagHist(sample.offFile, "h_off_" + sample.tag + "_" + def.tag, def.branch, edges);
    TH1D* hOn = makeDiagHist(sample.onFile, "h_on_" + sample.tag + "_" + def.tag, def.branch, edges);
    if (!hOff || !hOn) return nullptr;

    std::vector<double> x, y, ex, ey;
    int lastValidBin = -1;
    double lastOff = 0.0;
    double lastOn = 0.0;
    double lastRatio = 0.0;
    double lastError = 0.0;
    for (int ib = 1; ib <= hOff->GetNbinsX(); ++ib) {
        const double off = hOff->GetBinContent(ib);
        const double on = hOn->GetBinContent(ib);
        if (off <= 0.0 || on <= 0.0) continue;
        if (off < 10.0 || on < 10.0) continue;
        const double r = off / on;
        const double e = r * std::sqrt(1.0 / off + 1.0 / on);
        x.push_back(hOff->GetBinCenter(ib) + sample.xOffset);
        y.push_back(r);
        ex.push_back(0.0);
        ey.push_back(e);
        if (allBinSummary) {
            (*allBinSummary) << sample.label << "," << def.tag << "," << hOff->GetBinLowEdge(ib) << ","
                             << hOff->GetBinLowEdge(ib + 1) << "," << off << "," << on << ","
                             << r << "," << e << "\n";
        }
        lastValidBin = ib;
        lastOff = off;
        lastOn = on;
        lastRatio = r;
        lastError = e;
    }
    if (lastBinSummary && lastValidBin > 0) {
        (*lastBinSummary) << sample.label << "," << def.tag << "," << hOff->GetBinLowEdge(lastValidBin) << ","
                          << hOff->GetBinLowEdge(lastValidBin + 1) << "," << lastOff << "," << lastOn << ","
                          << lastRatio << "," << lastError << "\n";
    }
    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetName(("g_" + sample.tag + "_" + def.tag).c_str());
    g->SetLineColor(sample.color);
    g->SetMarkerColor(sample.color);
    g->SetMarkerStyle(sample.marker);
    g->SetMarkerSize(0.74);
    g->SetLineWidth(2);
    delete hOff;
    delete hOn;
    return g;
}

void drawOneDefinition(const ThrustDef& def,
                       const std::vector<CorrSample>& samples,
                       const std::vector<double>& edges,
                       std::ofstream& lastBinSummary,
                       std::ofstream& allBinSummary)
{
    std::vector<TGraphErrors*> graphs;
    for (const CorrSample& sample : samples) graphs.push_back(makeRatioGraph(sample, def, edges, &lastBinSummary, &allBinSummary));
    const double xmin = edges.front();
    const double xmax = edges.back();

    TCanvas* c = new TCanvas(("c_cisr_" + def.tag).c_str(), "CISR thrust definition", 1050, 780);
    c->SetLeftMargin(0.125);
    c->SetRightMargin(0.035);
    c->SetTopMargin(0.070);
    c->SetBottomMargin(0.105);
    c->SetTicks(1, 1);

    TH1D* frame = new TH1D(("frame_" + def.tag).c_str(), "", 100, xmin, xmax);
    frame->SetMinimum(0.55);
    frame->SetMaximum(1.45);
    frame->GetXaxis()->SetTitle("Thrust T");
    frame->GetYaxis()->SetTitle("C_{ISR} = ISR OFF / ISR ON");
    frame->GetYaxis()->SetTitleOffset(1.05);
    frame->Draw("AXIS");

    TLine* one = new TLine(xmin, 1.0, xmax, 1.0);
    one->SetLineColor(kGray + 2);
    one->SetLineStyle(2);
    one->SetLineWidth(2);
    one->Draw();

    for (TGraphErrors* g : graphs) {
        if (g) g->Draw("LP SAME");
    }

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.036);
    lat.DrawLatex(0.135, 0.890, def.label.c_str());
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.135, 0.845, "#sqrt{s}=91.1876 GeV, derived from EndpointDiagnostics trees");

    TLegend* leg = new TLegend(0.45, 0.18, 0.90, 0.39);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    for (size_t i = 0; i < samples.size(); ++i) {
        if (graphs[i]) leg->AddEntry(graphs[i], samples[i].label.c_str(), "lp");
    }
    leg->AddEntry(one, "C_{ISR}=1", "l");
    leg->Draw();
    frame->Draw("AXIS SAME");

    c->SaveAs(corrJoin(gCorrOutDir, "cisr_" + def.outputStem + ".png").c_str());
    c->SaveAs(corrJoin(gCorrOutDir, "cisr_" + def.outputStem + ".pdf").c_str());
}

void drawThreeDefinitionOverlay(const std::vector<CorrSample>& samples,
                                const std::vector<ThrustDef>& defs,
                                const std::vector<double>& edges)
{
    const bool fourPanel = defs.size() > 3;
    TCanvas* c = new TCanvas("c_cisr_definitions_panel", "CISR definition panel",
                             fourPanel ? 1350 : 1650,
                             fourPanel ? 1000 : 620);
    c->Divide(fourPanel ? 2 : 3, fourPanel ? 2 : 1, 0.002, 0.002);
    const double xmin = edges.front();
    const double xmax = edges.back();

    for (size_t id = 0; id < defs.size(); ++id) {
        c->cd(static_cast<int>(id + 1));
        gPad->SetLeftMargin((!fourPanel && id == 0) || (fourPanel && id % 2 == 0) ? 0.16 : 0.11);
        gPad->SetRightMargin(0.030);
        gPad->SetTopMargin(0.085);
        gPad->SetBottomMargin(0.115);
        gPad->SetTicks(1, 1);

        TH1D* frame = new TH1D(("frame_panel_" + defs[id].tag).c_str(), "", 100, xmin, xmax);
        frame->SetMinimum(0.55);
        frame->SetMaximum(1.45);
        frame->GetXaxis()->SetTitle("Thrust T");
        frame->GetYaxis()->SetTitle("C_{ISR}");
        frame->Draw("AXIS");
        TLine* one = new TLine(xmin, 1.0, xmax, 1.0);
        one->SetLineColor(kGray + 2);
        one->SetLineStyle(2);
        one->Draw();

        std::vector<TGraphErrors*> graphs;
        for (const CorrSample& sample : samples) {
            TGraphErrors* g = makeRatioGraph(sample, defs[id], edges, nullptr, nullptr);
            graphs.push_back(g);
            if (g) g->Draw("LP SAME");
        }

        TLatex lat;
        lat.SetNDC();
        lat.SetTextFont(42);
        lat.SetTextSize(0.042);
        lat.DrawLatex(0.16, 0.905, defs[id].tag.c_str());
        lat.SetTextSize(0.030);
        lat.DrawLatex(0.16, 0.855, defs[id].label.c_str());

        if (id == 0) {
            TLegend* leg = new TLegend(0.18, 0.18, 0.92, 0.39);
            leg->SetBorderSize(0);
            leg->SetFillStyle(0);
            leg->SetTextFont(42);
            leg->SetTextSize(0.030);
            leg->SetNColumns(1);
            for (size_t i = 0; i < samples.size(); ++i) {
                if (graphs[i]) leg->AddEntry(graphs[i], samples[i].label.c_str(), "lp");
            }
            leg->Draw();
        }
        frame->Draw("AXIS SAME");
    }

    c->SaveAs(corrJoin(gCorrOutDir, "cisr_thrust_definition_panel.png").c_str());
    c->SaveAs(corrJoin(gCorrOutDir, "cisr_thrust_definition_panel.pdf").c_str());
}

void plot_thrust_definition_corrections()
{
    if (const char* envReal = gSystem->Getenv("ISR_REAL_DIR")) {
        gCorrRealDir = envReal;
        gCorrDiagDir = corrJoin(gCorrRealDir, "endpoint_diagnostics");
        gCorrOutDir = corrJoin(corrJoin(gCorrRealDir, "results"), "thrust_definition_corrections");
    }
    if (const char* envDiag = gSystem->Getenv("ISR_ENDPOINT_DIAG_DIR")) {
        gCorrDiagDir = envDiag;
    }
    if (const char* envOut = gSystem->Getenv("ISR_THRUST_DEF_OUT_DIR")) {
        gCorrOutDir = envOut;
    }
    gSystem->mkdir(gCorrOutDir.c_str(), true);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    const int colorBlue = TColor::GetColor("#0072B2");
    const int colorVermillion = TColor::GetColor("#D55E00");
    const int colorGreen = TColor::GetColor("#009E73");
    const int colorPurple = TColor::GetColor("#CC79A7");
    const int colorOrange = TColor::GetColor("#E69F00");

    const std::vector<CorrSample> samples = {
        {"Herwig 7.3.0 QED shower", "Herwig730_QEDshower",
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Herwig730_OFF.root"),
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Herwig730_QEDshower.root"),
         colorBlue, 20, -0.0018},
        {"PYTHIA 8.315", "Pythia8315",
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Pythia8315_OFF.root"),
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Pythia8315.root"),
         colorVermillion, 21, -0.0009},
        {"Sherpa 3.0.3 PDFESherpa", "Sherpa303",
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Sherpa303_OFF.root"),
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Sherpa303.root"),
         colorGreen, 23, 0.0000},
        {"Sherpa 3.0.3 YFS", "Sherpa303_YFS",
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Sherpa303_OFF.root"),
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Sherpa303_YFS.root"),
         colorPurple, 34, 0.0009},
        {"PYTHIA 8.315 (Vincia)", "Pythia8315_Vincia",
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Pythia8315_Vincia_OFF.root"),
         corrJoin(gCorrDiagDir, "endpoint_diagnostics_Pythia8315_Vincia.root"),
         colorOrange, 22, 0.0018}
    };

    const std::vector<ThrustDef> defs = {
        {"N", "T_lab_allFinal_including_ISR_photons", "N: nominal all final particles, including #nu and ISR #gamma", "N_lab_allFinal_including_ISR"},
        {"A", "T_lab_excluding_ISR_photons", "A: lab visible thrust, excluding tagged ISR photons", "A_lab_visible_excluding_ISR"},
        {"B", "T_lab_including_ISR_photons", "B: lab visible thrust, including tagged ISR photons", "B_lab_visible_including_ISR"},
        {"C", "T_visibleCM_excluding_ISR_photons", "C: visible-CM thrust, excluding tagged ISR photons", "C_visibleCM_excluding_ISR"}
    };

    const std::vector<double> edges = corrAlephThrustEdges();

    std::ofstream lastBinSummary(corrJoin(gCorrOutDir, "cisr_thrust_definition_last_bin_summary.csv").c_str());
    std::ofstream allBinSummary(corrJoin(gCorrOutDir, "cisr_thrust_definition_all_bins.csv").c_str());
    lastBinSummary << std::setprecision(10);
    allBinSummary << std::setprecision(10);
    lastBinSummary << "sample,definition,bin_low,bin_high,off_count,on_count,ratio,error\n";
    allBinSummary << "sample,definition,bin_low,bin_high,off_count,on_count,ratio,error\n";
    for (const ThrustDef& def : defs) drawOneDefinition(def, samples, edges, lastBinSummary, allBinSummary);
    drawThreeDefinitionOverlay(samples, defs, edges);
}
