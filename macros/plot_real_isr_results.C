#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TMath.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct SampleDef {
    std::string label;
    std::string tag;
    std::string offFile;
    std::string onFile;
    std::string onLabel;
    int color;
    int marker;
    double xOffset;
};

struct AlephPoint {
    double center;
    double low;
    double high;
    double value;
    double err;
};

std::string gRealDir = "/data2/yjlee/ISRsample/real_3M_20260511";
std::string gOutDir = "/data2/yjlee/ISRsample/real_3M_20260511/results";
std::string gDiagDir = "/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics";
std::string gNominalThrustBranch = "T_lab_allFinal_including_ISR_photons";

std::string joinPath(const std::string& dir, const std::string& name)
{
    return dir + "/" + name;
}

constexpr double kVisibleEtaMax = 1.74;

bool isNeutrinoPdg(int pdg)
{
    const int apdg = std::abs(pdg);
    return apdg == 12 || apdg == 14 || apdg == 16;
}

double etaFromMomentum(double px, double py, double pz)
{
    const double p = std::sqrt(px * px + py * py + pz * pz);
    if (p <= 1e-12) return 1e9;
    const double plus = p + pz;
    const double minus = p - pz;
    if (plus <= 0.0) return -1e9;
    if (minus <= 0.0) return 1e9;
    return 0.5 * std::log(plus / minus);
}

double acceptedVisibleEnergy(const std::vector<int>& pdgId,
                             const std::vector<float>& px,
                             const std::vector<float>& py,
                             const std::vector<float>& pz,
                             const std::vector<float>& energy,
                             const std::vector<char>& isFinal)
{
    double evis = 0.0;
    const size_t n = pdgId.size();
    for (size_t i = 0; i < n; ++i) {
        if (!isFinal[i]) continue;
        if (isNeutrinoPdg(pdgId[i])) continue;
        const double eta = etaFromMomentum(px[i], py[i], pz[i]);
        if (std::abs(eta) >= kVisibleEtaMax) continue;
        evis += energy[i];
    }
    return evis;
}

std::vector<AlephPoint> readAleph(const std::string& filename)
{
    std::vector<AlephPoint> points;
    std::ifstream in(filename.c_str());
    std::string line;
    while (std::getline(in, line)) {
        if (line.empty() || line[0] == '#') continue;
        if (line.find("THRUST,THRUST LOW") != std::string::npos) continue;
        std::stringstream ss(line);
        std::string item;
        std::vector<std::string> cols;
        while (std::getline(ss, item, ',')) cols.push_back(item);
        if (cols.size() < 10) continue;
        AlephPoint p;
        p.center = std::stod(cols[0]);
        p.low = std::stod(cols[1]);
        p.high = std::stod(cols[2]);
        p.value = std::stod(cols[3]);
        const double stat = std::abs(std::stod(cols[4]));
        const double sys1 = std::abs(std::stod(cols[6]));
        const double sys2 = std::abs(std::stod(cols[8]));
        p.err = std::sqrt(stat * stat + sys1 * sys1 + sys2 * sys2);
        points.push_back(p);
    }
    return points;
}

std::vector<double> alephEdges(const std::vector<AlephPoint>& points)
{
    std::vector<double> edges;
    if (points.empty()) return edges;
    edges.push_back(points.front().low);
    for (const AlephPoint& p : points) edges.push_back(p.high);
    return edges;
}

std::string alephThrustCsv()
{
    std::string alephCsv = "/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/ALEPH/HEPData-ins636645-v1-Table_54.csv";
    if (const char* envAleph = gSystem->Getenv("ALEPH_THRUST_CSV")) alephCsv = envAleph;
    return alephCsv;
}

std::vector<double> alephThrustEdges()
{
    std::vector<double> edges = alephEdges(readAleph(alephThrustCsv()));
    if (edges.size() > 1) return edges;
    for (int i = 0; i <= 42; ++i) edges.push_back(0.58 + 0.01 * i);
    return edges;
}

TGraphErrors* makeAlephGraph(const std::vector<AlephPoint>& aleph, const std::string& name)
{
    std::vector<double> ax, ay, aex, aey;
    for (const AlephPoint& p : aleph) {
        ax.push_back(p.center);
        ay.push_back(p.value);
        aex.push_back(0.5 * (p.high - p.low));
        aey.push_back(p.err);
    }
    TGraphErrors* g = new TGraphErrors(ax.size(), ax.data(), ay.data(), aex.data(), aey.data());
    g->SetName(name.c_str());
    g->SetMarkerStyle(20);
    g->SetMarkerSize(0.65);
    g->SetMarkerColor(kBlack);
    g->SetLineColor(kBlack);
    g->SetLineWidth(1);
    return g;
}

TGraphErrors* makeMCAlephRatioGraph(const TH1D* h,
                                    const std::vector<AlephPoint>& aleph,
                                    const SampleDef& s,
                                    const std::string& name)
{
    std::vector<double> x, y, ex, ey;
    const int n = std::min(h ? h->GetNbinsX() : 0, static_cast<int>(aleph.size()));
    for (int i = 0; i < n; ++i) {
        const AlephPoint& p = aleph[i];
        const double mc = h->GetBinContent(i + 1);
        const double emc = h->GetBinError(i + 1);
        if (mc <= 0.0 || p.value <= 0.0) continue;
        const double r = mc / p.value;
        const double relMC = emc / mc;
        const double relData = p.err / p.value;
        x.push_back(p.center + s.xOffset);
        y.push_back(r);
        ex.push_back(0.0);
        ey.push_back(r * std::sqrt(relMC * relMC + relData * relData));
    }
    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetName(name.c_str());
    g->SetLineColor(s.color);
    g->SetMarkerColor(s.color);
    g->SetMarkerStyle(s.marker);
    g->SetMarkerSize(0.62);
    g->SetLineWidth(2);
    return g;
}

void normalizeDensity(TH1D* h)
{
    const double integral = h->Integral();
    if (integral <= 0.0) return;
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
        const double width = h->GetBinWidth(ib);
        h->SetBinContent(ib, h->GetBinContent(ib) / (integral * width));
        h->SetBinError(ib, h->GetBinError(ib) / (integral * width));
    }
}

TH1D* makeHist(const std::string& filename, const std::string& name, const std::string& expr,
               int nbins, double xmin, double xmax, const std::string& selection = "")
{
    TFile f(filename.c_str());
    TTree* t = static_cast<TTree*>(f.Get("Events"));
    if (!t) {
        std::cerr << "Missing Events tree in " << filename << std::endl;
        return nullptr;
    }
    gROOT->cd();
    TH1D* h = new TH1D(name.c_str(), "", nbins, xmin, xmax);
    h->Sumw2();
    h->SetDirectory(gROOT);
    t->Draw((expr + ">>" + name).c_str(), selection.c_str(), "goff");
    h->SetDirectory(0);
    return h;
}

TH1D* makeVisibleEnergyHist(const std::string& filename, const std::string& name,
                            int nbins, double xmin, double xmax)
{
    TFile f(filename.c_str());
    TTree* t = static_cast<TTree*>(f.Get("Events"));
    if (!t) {
        std::cerr << "Missing Events tree in " << filename << std::endl;
        return nullptr;
    }
    gROOT->cd();
    TH1D* h = new TH1D(name.c_str(), "", nbins, xmin, xmax);
    h->Sumw2();
    h->SetDirectory(gROOT);

    TTreeReader reader(t);
    TTreeReaderValue<std::vector<int>> pdgId(reader, "pdgId");
    TTreeReaderValue<std::vector<float>> px(reader, "px");
    TTreeReaderValue<std::vector<float>> py(reader, "py");
    TTreeReaderValue<std::vector<float>> pz(reader, "pz");
    TTreeReaderValue<std::vector<float>> energy(reader, "energy");
    TTreeReaderValue<std::vector<char>> isFinal(reader, "isFinal");
    while (reader.Next()) {
        h->Fill(acceptedVisibleEnergy(*pdgId, *px, *py, *pz, *energy, *isFinal));
    }
    h->SetDirectory(0);
    return h;
}

TH1D* makeHistEdges(const std::string& filename, const std::string& name, const std::string& expr,
                    const std::vector<double>& edges, const std::string& selection = "")
{
    TFile f(filename.c_str());
    TTree* t = static_cast<TTree*>(f.Get("Events"));
    if (!t) {
        std::cerr << "Missing Events tree in " << filename << std::endl;
        return nullptr;
    }
    gROOT->cd();
    TH1D* h = new TH1D(name.c_str(), "", static_cast<int>(edges.size()) - 1, edges.data());
    h->Sumw2();
    h->SetDirectory(gROOT);
    t->Draw((expr + ">>" + name).c_str(), selection.c_str(), "goff");
    h->SetDirectory(0);
    return h;
}

std::string endpointTag(const SampleDef& s, bool useOn)
{
    if (useOn) return s.tag;
    if (s.tag == "Herwig730_QEDshower") return "Herwig730_OFF";
    if (s.tag == "Pythia8315") return "Pythia8315_OFF";
    if (s.tag == "Sherpa303" || s.tag == "Sherpa303_YFS") return "Sherpa303_OFF";
    if (s.tag == "Pythia8315_Vincia") return "Pythia8315_Vincia_OFF";
    return s.tag + "_OFF";
}

std::string endpointFile(const SampleDef& s, bool useOn)
{
    return joinPath(gDiagDir, "endpoint_diagnostics_" + endpointTag(s, useOn) + ".root");
}

TH1D* makeEndpointThrustHistEdges(const SampleDef& s,
                                  bool useOn,
                                  const std::string& name,
                                  const std::vector<double>& edges,
                                  const std::string& branch = "")
{
    const std::string filename = endpointFile(s, useOn);
    TFile f(filename.c_str());
    TTree* t = static_cast<TTree*>(f.Get("EndpointDiagnostics"));
    if (!t) {
        std::cerr << "Missing EndpointDiagnostics tree in " << filename << std::endl;
        return nullptr;
    }
    gROOT->cd();
    TH1D* h = new TH1D(name.c_str(), "", static_cast<int>(edges.size()) - 1, edges.data());
    h->Sumw2();
    h->SetDirectory(gROOT);
    const std::string drawBranch = branch.empty() ? gNominalThrustBranch : branch;
    t->Draw((drawBranch + ">>" + name).c_str(), "", "goff");
    h->SetDirectory(0);
    return h;
}

TGraphErrors* makeCorrectionGraph(const SampleDef& s, const std::vector<double>& edges)
{
    TH1D* hOff = makeEndpointThrustHistEdges(s, false, "h_off_" + s.tag, edges);
    TH1D* hOn = makeEndpointThrustHistEdges(s, true, "h_on_" + s.tag, edges);
    std::vector<double> x, y, ex, ey;
    for (int ib = 1; ib <= hOff->GetNbinsX(); ++ib) {
        const double a = hOff->GetBinContent(ib);
        const double b = hOn->GetBinContent(ib);
        if (b <= 0.0 || a <= 0.0) continue;
        if (a < 10.0 || b < 10.0) continue;
        const double r = a / b;
        const double err = r * std::sqrt(1.0 / a + 1.0 / b);
        x.push_back(hOff->GetBinCenter(ib) + s.xOffset);
        y.push_back(r);
        ex.push_back(0.0);
        ey.push_back(err);
    }
    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetName(("g_corr_" + s.tag).c_str());
    g->SetLineColor(s.color);
    g->SetMarkerColor(s.color);
    g->SetMarkerStyle(s.marker);
    g->SetMarkerSize(0.82);
    g->SetLineWidth(3);
    return g;
}

TGraphErrors* makeDoubleRatio(TGraphErrors* g, TGraphErrors* ref, const SampleDef& s)
{
    std::vector<double> x, y, ex, ey;
    const int n = std::min(g->GetN(), ref->GetN());
    for (int i = 0; i < n; ++i) {
        double gx, gy, rx, ry;
        g->GetPoint(i, gx, gy);
        ref->GetPoint(i, rx, ry);
        if (ry <= 0.0) continue;
        x.push_back(gx);
        y.push_back(gy / ry);
        ex.push_back(0.0);
        const double eg = g->GetErrorY(i);
        const double er = ref->GetErrorY(i);
        ey.push_back((gy / ry) * std::sqrt((eg / gy) * (eg / gy) + (er / ry) * (er / ry)));
    }
    TGraphErrors* out = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    out->SetName(("g_dr_" + s.tag).c_str());
    out->SetLineColor(s.color);
    out->SetMarkerColor(s.color);
    out->SetMarkerStyle(s.marker);
    out->SetMarkerSize(0.74);
    out->SetLineWidth(3);
    return out;
}

void drawISRCorrection(const std::vector<SampleDef>& samples)
{
    const std::vector<double> edges = alephThrustEdges();
    const double xmin = edges.front();
    const double xmax = edges.back();

    std::vector<TGraphErrors*> graphs;
    for (const SampleDef& s : samples) graphs.push_back(makeCorrectionGraph(s, edges));
    std::vector<TGraphErrors*> dr;
    for (size_t i = 0; i < samples.size(); ++i) dr.push_back(makeDoubleRatio(graphs[i], graphs[1], samples[i]));

    TCanvas* c = new TCanvas("c_real_isr_correction", "real ISR correction", 1000, 900);
    TPad* top = new TPad("top", "top", 0.0, 0.30, 1.0, 1.0);
    TPad* bot = new TPad("bot", "bot", 0.0, 0.0, 1.0, 0.30);
    top->SetLeftMargin(0.145); top->SetRightMargin(0.035); top->SetBottomMargin(0.020); top->SetTopMargin(0.055); top->SetTicks(1, 1);
    bot->SetLeftMargin(0.145); bot->SetRightMargin(0.035); bot->SetBottomMargin(0.310); bot->SetTopMargin(0.020); bot->SetTicks(1, 1);
    top->Draw(); bot->Draw();

    top->cd();
    TH1D* frame = new TH1D("frame_corr", "", 100, xmin, xmax);
    frame->SetMinimum(0.55);
    frame->SetMaximum(1.42);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("C_{ISR} = ISR OFF / ISR ON");
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetLabelSize(0.044);
    frame->Draw("AXIS");
    TLine* one = new TLine(xmin, 1.0, xmax, 1.0);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 2);
    one->Draw();
    for (TGraphErrors* g : graphs) g->Draw("LP SAME");
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.040);
    lat.DrawLatex(0.145, 0.875, "ALEPH archived data");
    lat.DrawLatex(0.630, 0.875, "e^{+}e^{-}  #sqrt{s}=91.1876 GeV");
    lat.SetTextSize(0.028);
    lat.DrawLatex(0.145, 0.815, "Nominal truth thrust: all final particles incl. neutrinos and tagged ISR photons");
    TLegend* leg = new TLegend(0.43, 0.18, 0.91, 0.38);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetNColumns(2);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    for (size_t i = 0; i < samples.size(); ++i) leg->AddEntry(graphs[i], samples[i].label.c_str(), "lp");
    leg->AddEntry(one, "C_{ISR}=1", "l");
    leg->Draw();
    frame->Draw("AXIS SAME");

    bot->cd();
    TH1D* rframe = new TH1D("frame_dr", "", 100, xmin, xmax);
    rframe->SetMinimum(0.50);
    rframe->SetMaximum(1.50);
    rframe->GetXaxis()->SetTitle("Thrust T");
    rframe->GetYaxis()->SetTitle("C_{ISR}/PYTHIA");
    rframe->GetXaxis()->SetTitleSize(0.105);
    rframe->GetXaxis()->SetLabelSize(0.085);
    rframe->GetYaxis()->SetTitleSize(0.080);
    rframe->GetYaxis()->SetTitleOffset(0.82);
    rframe->GetYaxis()->SetLabelSize(0.075);
    rframe->Draw("AXIS");
    TLine* one2 = new TLine(xmin, 1.0, xmax, 1.0);
    one2->SetLineStyle(2);
    one2->SetLineColor(kGray + 2);
    one2->Draw();
    for (TGraphErrors* g : dr) g->Draw("LP SAME");
    rframe->Draw("AXIS SAME");

    c->SaveAs(joinPath(gOutDir, "real_isr_correction_double_ratio.png").c_str());
    c->SaveAs(joinPath(gOutDir, "real_isr_correction_double_ratio.pdf").c_str());
    c->SaveAs(joinPath(gOutDir, "isr_correction_studies_reproduction.png").c_str());
    c->SaveAs(joinPath(gOutDir, "isr_correction_studies_reproduction.pdf").c_str());
}

void drawThrustVsAleph(const std::vector<SampleDef>& samples)
{
    std::string alephCsv = alephThrustCsv();
    const std::vector<AlephPoint> aleph = readAleph(alephCsv);
    const std::vector<double> edges = alephEdges(aleph);
    std::vector<double> ax, ay, aex, aey;
    for (const AlephPoint& p : aleph) {
        ax.push_back(p.center);
        ay.push_back(p.value);
        aex.push_back(0.5 * (p.high - p.low));
        aey.push_back(p.err);
    }
    TGraphErrors* gAleph = new TGraphErrors(ax.size(), ax.data(), ay.data(), aex.data(), aey.data());
    gAleph->SetMarkerStyle(20);
    gAleph->SetMarkerSize(0.65);
    gAleph->SetMarkerColor(kBlack);
    gAleph->SetLineColor(kBlack);

    TCanvas* c = new TCanvas("c_real_thrust_aleph", "real thrust vs ALEPH", 1050, 900);
    c->SetLogy();
    c->SetLeftMargin(0.145);
    c->SetRightMargin(0.035);
    c->SetTopMargin(0.055);
    c->SetBottomMargin(0.105);
    c->SetTicks(1, 1);
    TH1D* frame = new TH1D("frame_aleph", "", 100, 0.58, 1.0);
    frame->SetMinimum(5e-4);
    frame->SetMaximum(45.0);
    frame->GetXaxis()->SetTitle("Thrust T");
    frame->GetYaxis()->SetTitle("(1/N) dN/dT");
    frame->GetYaxis()->SetTitleOffset(1.10);
    frame->Draw("AXIS");
    gAleph->Draw("P SAME");
    TLegend* leg = new TLegend(0.14, 0.68, 0.58, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry(gAleph, "ALEPH published data", "pe");
    for (const SampleDef& s : samples) {
        TH1D* h = makeEndpointThrustHistEdges(s, true, "h_aleph_" + s.tag, edges);
        normalizeDensity(h);
        h->SetLineColor(s.color);
        h->SetMarkerColor(s.color);
        h->SetLineWidth(3);
        h->SetMarkerStyle(s.marker);
        h->SetMarkerSize(0.48);
        h->Draw("HIST SAME");
        leg->AddEntry(h, s.label.c_str(), "l");
    }
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.58, 0.86, "Truth-level all final particles, incl. #nu and ISR #gamma");
    leg->Draw();
    frame->Draw("AXIS SAME");
    c->SaveAs(joinPath(gOutDir, "real_isr_on_thrust_vs_aleph.png").c_str());
    c->SaveAs(joinPath(gOutDir, "real_isr_on_thrust_vs_aleph.pdf").c_str());
}

std::string offLegendLabel(const SampleDef& s)
{
    if (s.tag == "Herwig730_QEDshower") return "Herwig 7.3.0 ISR OFF";
    if (s.tag == "Sherpa303") return "Sherpa 3.0.3 ISR OFF";
    return s.label + " ISR OFF";
}

void drawThrustOffVsAleph(const std::vector<SampleDef>& samples)
{
    std::string alephCsv = alephThrustCsv();
    const std::vector<AlephPoint> aleph = readAleph(alephCsv);
    const std::vector<double> edges = alephEdges(aleph);
    std::vector<double> ax, ay, aex, aey;
    for (const AlephPoint& p : aleph) {
        ax.push_back(p.center);
        ay.push_back(p.value);
        aex.push_back(0.5 * (p.high - p.low));
        aey.push_back(p.err);
    }
    TGraphErrors* gAleph = new TGraphErrors(ax.size(), ax.data(), ay.data(), aex.data(), aey.data());
    gAleph->SetMarkerStyle(20);
    gAleph->SetMarkerSize(0.65);
    gAleph->SetMarkerColor(kBlack);
    gAleph->SetLineColor(kBlack);

    TCanvas* c = new TCanvas("c_real_thrust_aleph_off", "real ISR OFF thrust vs ALEPH", 1050, 900);
    c->SetLogy();
    c->SetLeftMargin(0.145);
    c->SetRightMargin(0.035);
    c->SetTopMargin(0.055);
    c->SetBottomMargin(0.105);
    c->SetTicks(1, 1);
    TH1D* frame = new TH1D("frame_aleph_off", "", 100, 0.58, 1.0);
    frame->SetMinimum(5e-4);
    frame->SetMaximum(45.0);
    frame->GetXaxis()->SetTitle("Thrust T");
    frame->GetYaxis()->SetTitle("(1/N) dN/dT");
    frame->GetYaxis()->SetTitleOffset(1.10);
    frame->Draw("AXIS");
    gAleph->Draw("P SAME");
    TLegend* leg = new TLegend(0.14, 0.68, 0.58, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.030);
    leg->AddEntry(gAleph, "ALEPH published data", "pe");

    std::vector<std::string> drawnOffFiles;
    for (const SampleDef& s : samples) {
        if (std::find(drawnOffFiles.begin(), drawnOffFiles.end(), s.offFile) != drawnOffFiles.end()) continue;
        drawnOffFiles.push_back(s.offFile);
        TH1D* h = makeEndpointThrustHistEdges(s, false, "h_aleph_off_" + s.tag, edges);
        normalizeDensity(h);
        h->SetLineColor(s.color);
        h->SetMarkerColor(s.color);
        h->SetLineWidth(3);
        h->SetMarkerStyle(s.marker);
        h->SetMarkerSize(0.48);
        h->Draw("HIST SAME");
        leg->AddEntry(h, offLegendLabel(s).c_str(), "l");
    }
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.58, 0.86, "Truth-level all final particles, incl. #nu and ISR #gamma");
    lat.SetTextSize(0.026);
    lat.DrawLatex(0.58, 0.815, "Sherpa PDFESherpa and YFS share the same OFF file");
    leg->Draw();
    frame->Draw("AXIS SAME");
    c->SaveAs(joinPath(gOutDir, "real_isr_off_thrust_vs_aleph.png").c_str());
    c->SaveAs(joinPath(gOutDir, "real_isr_off_thrust_vs_aleph.pdf").c_str());
}

void drawThrustVsAlephRatioPanels(const std::vector<SampleDef>& samples, bool useISRLikeSample)
{
    std::string alephCsv = alephThrustCsv();
    const std::vector<AlephPoint> aleph = readAleph(alephCsv);
    const std::vector<double> edges = alephEdges(aleph);
    TGraphErrors* gAleph = makeAlephGraph(aleph, useISRLikeSample ? "g_aleph_ratio_on" : "g_aleph_ratio_off");

    const std::string suffix = useISRLikeSample ? "on" : "off";
    TCanvas* c = new TCanvas(("c_real_thrust_aleph_ratio_" + suffix).c_str(),
                             "real thrust vs ALEPH with ratio", 1050, 900);
    TPad* top = new TPad(("top_aleph_ratio_" + suffix).c_str(), "top", 0.0, 0.30, 1.0, 1.0);
    TPad* bot = new TPad(("bot_aleph_ratio_" + suffix).c_str(), "bottom", 0.0, 0.0, 1.0, 0.30);
    top->SetLeftMargin(0.145);
    top->SetRightMargin(0.035);
    top->SetTopMargin(0.055);
    top->SetBottomMargin(0.020);
    top->SetTicks(1, 1);
    top->SetLogy();
    bot->SetLeftMargin(0.145);
    bot->SetRightMargin(0.035);
    bot->SetTopMargin(0.025);
    bot->SetBottomMargin(0.310);
    bot->SetTicks(1, 1);
    top->Draw();
    bot->Draw();

    top->cd();
    TH1D* frame = new TH1D(("frame_aleph_ratio_" + suffix).c_str(), "", 100, 0.58, 1.0);
    frame->SetMinimum(5e-4);
    frame->SetMaximum(45.0);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("(1/N) dN/dT");
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetLabelSize(0.044);
    frame->Draw("AXIS");

    TLegend* leg = new TLegend(0.14, 0.67, 0.58, 0.89);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.029);
    leg->AddEntry(gAleph, "ALEPH published data", "pe");

    std::vector<TGraphErrors*> ratios;
    std::vector<std::string> drawnOffFiles;
    for (const SampleDef& s : samples) {
        if (!useISRLikeSample) {
            if (std::find(drawnOffFiles.begin(), drawnOffFiles.end(), s.offFile) != drawnOffFiles.end()) continue;
            drawnOffFiles.push_back(s.offFile);
        }
        const std::string histName = useISRLikeSample ? ("h_aleph_ratio_on_" + s.tag)
                                                      : ("h_aleph_ratio_off_" + s.tag);
        TH1D* h = makeEndpointThrustHistEdges(s, useISRLikeSample, histName, edges);
        normalizeDensity(h);
        h->SetLineColor(s.color);
        h->SetMarkerColor(s.color);
        h->SetLineWidth(3);
        h->SetMarkerStyle(s.marker);
        h->SetMarkerSize(0.48);
        h->Draw("HIST SAME");
        leg->AddEntry(h, useISRLikeSample ? s.label.c_str() : offLegendLabel(s).c_str(), "l");

        const std::string ratioName = useISRLikeSample ? ("g_mc_aleph_ratio_on_" + s.tag)
                                                       : ("g_mc_aleph_ratio_off_" + s.tag);
        ratios.push_back(makeMCAlephRatioGraph(h, aleph, s, ratioName));
    }

    gAleph->Draw("P SAME");
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.034);
    if (useISRLikeSample) {
        lat.DrawLatex(0.58, 0.86, "Truth-level all final particles, incl. #nu and ISR #gamma");
    } else {
        lat.DrawLatex(0.58, 0.86, "Truth-level all final particles, incl. #nu and ISR #gamma");
        lat.SetTextSize(0.026);
        lat.DrawLatex(0.58, 0.815, "Sherpa PDFESherpa and YFS share the same OFF file");
    }
    leg->Draw();
    frame->Draw("AXIS SAME");

    bot->cd();
    TH1D* rframe = new TH1D(("frame_mc_aleph_ratio_" + suffix).c_str(), "", 100, 0.58, 1.0);
    rframe->SetMinimum(0.35);
    rframe->SetMaximum(1.65);
    rframe->GetXaxis()->SetTitle("Thrust T");
    rframe->GetYaxis()->SetTitle("MC / ALEPH");
    rframe->GetXaxis()->SetTitleSize(0.105);
    rframe->GetXaxis()->SetLabelSize(0.085);
    rframe->GetYaxis()->SetTitleSize(0.080);
    rframe->GetYaxis()->SetTitleOffset(0.72);
    rframe->GetYaxis()->SetLabelSize(0.075);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->Draw("AXIS");
    TLine* one = new TLine(0.58, 1.0, 1.0, 1.0);
    one->SetLineStyle(2);
    one->SetLineColor(kGray + 2);
    one->SetLineWidth(2);
    one->Draw();
    for (TGraphErrors* g : ratios) g->Draw("LP SAME");
    rframe->Draw("AXIS SAME");

    c->SaveAs(joinPath(gOutDir, "real_isr_" + suffix + "_thrust_vs_aleph_ratio.png").c_str());
    c->SaveAs(joinPath(gOutDir, "real_isr_" + suffix + "_thrust_vs_aleph_ratio.pdf").c_str());
}

void drawISRPhotonSpectra(const std::vector<SampleDef>& samples)
{
    TCanvas* c = new TCanvas("c_isr_photons", "ISR photon spectra", 1200, 520);
    c->Divide(2, 1);
    c->cd(1);
    gPad->SetTicks(1, 1);
    gPad->SetLogy();
    TLegend* leg = new TLegend(0.44, 0.68, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.030);
    bool first = true;
    for (const SampleDef& s : samples) {
        TH1D* h = makeHist(s.onFile, "h_isr_e_" + s.tag, "energy", 80, 0.0, 50.0, "isISRPhoton");
        normalizeDensity(h);
        h->SetLineColor(s.color);
        h->SetLineWidth(3);
        h->SetMinimum(1e-5);
        h->SetMaximum(4.0);
        h->GetXaxis()->SetTitle("ISR photon E [GeV]");
        h->GetYaxis()->SetTitle("normalized entries");
        h->Draw(first ? "HIST" : "HIST SAME");
        leg->AddEntry(h, s.label.c_str(), "l");
        first = false;
    }
    leg->Draw();
    c->cd(2);
    gPad->SetTicks(1, 1);
    gPad->SetLogy();
    first = true;
    for (const SampleDef& s : samples) {
        TH1D* h = makeHist(s.onFile, "h_isr_theta_" + s.tag,
                           "TMath::ACos(pz/TMath::Sqrt(px*px+py*py+pz*pz))*180.0/TMath::Pi()",
                           90, 0.0, 180.0, "isISRPhoton");
        normalizeDensity(h);
        h->SetLineColor(s.color);
        h->SetLineWidth(3);
        h->SetMinimum(1e-5);
        h->SetMaximum(4.0);
        h->GetXaxis()->SetTitle("ISR photon #theta [deg]");
        h->GetYaxis()->SetTitle("normalized entries");
        h->Draw(first ? "HIST" : "HIST SAME");
        first = false;
    }
    c->SaveAs(joinPath(gOutDir, "real_isr_photon_energy_theta_spectra.png").c_str());
    c->SaveAs(joinPath(gOutDir, "real_isr_photon_energy_theta_spectra.pdf").c_str());
}

void drawVisibleEnergy(const std::vector<SampleDef>& samples)
{
    for (const SampleDef& s : samples) {
        TCanvas* c = new TCanvas(("c_evis_" + s.tag).c_str(), "visible energy", 850, 650);
        c->SetTicks(1, 1);
        TH1D* hOff = makeVisibleEnergyHist(s.offFile, "h_evis_off_" + s.tag, 80, 0.0, 100.0);
        TH1D* hOn = makeVisibleEnergyHist(s.onFile, "h_evis_on_" + s.tag, 80, 0.0, 100.0);
        normalizeDensity(hOff);
        normalizeDensity(hOn);
        hOff->SetLineColor(kBlack);
        hOff->SetLineWidth(3);
        hOn->SetLineColor(s.color);
        hOn->SetLineWidth(3);
        hOff->SetMaximum(std::max(hOff->GetMaximum(), hOn->GetMaximum()) * 1.25);
        hOff->GetXaxis()->SetTitle("Visible energy E_{vis} (|#eta|<1.74) [GeV]");
        hOff->GetYaxis()->SetTitle("normalized entries");
        hOff->Draw("HIST");
        hOn->Draw("HIST SAME");
        TLegend* leg = new TLegend(0.15, 0.74, 0.58, 0.88);
        leg->SetBorderSize(0);
        leg->SetFillStyle(0);
        leg->SetTextFont(42);
        leg->SetTextSize(0.034);
        leg->AddEntry(hOff, (s.label + " ISR OFF").c_str(), "l");
        leg->AddEntry(hOn, (s.label + " " + s.onLabel).c_str(), "l");
        leg->Draw();
        c->SaveAs(joinPath(gOutDir, "real_visible_energy_" + s.tag + "_ISR_ON_OFF.png").c_str());
        c->SaveAs(joinPath(gOutDir, "real_visible_energy_" + s.tag + "_ISR_ON_OFF.pdf").c_str());
    }
}

double treeMean(const std::string& filename, const std::string& branch)
{
    TFile f(filename.c_str());
    TTree* t = static_cast<TTree*>(f.Get("Events"));
    t->Draw(branch.c_str(), "", "goff");
    return TMath::Mean(t->GetSelectedRows(), t->GetV1());
}

void writeStats(const std::vector<SampleDef>& samples)
{
    std::ofstream out(joinPath(gOutDir, "real_generator_sample_statistics.csv").c_str());
    out << std::setprecision(10);
    out << "sample,isr,entries,mean_thrust,mean_visibleEnergy,mean_nISRPhotons,mean_totalISRPhotonEnergy,file\n";
    for (const SampleDef& s : samples) {
        for (const auto& modeFile : {std::make_pair(std::string("OFF"), s.offFile), std::make_pair(s.onLabel, s.onFile)}) {
            TFile f(modeFile.second.c_str());
            TTree* t = static_cast<TTree*>(f.Get("Events"));
            TTreeReader reader(t);
            TTreeReaderValue<double> thrust(reader, "thrust");
            TTreeReaderValue<int> nISRPhotons(reader, "nISRPhotons");
            TTreeReaderValue<double> totalISRPhotonEnergy(reader, "totalISRPhotonEnergy");
            TTreeReaderValue<std::vector<int>> pdgId(reader, "pdgId");
            TTreeReaderValue<std::vector<float>> px(reader, "px");
            TTreeReaderValue<std::vector<float>> py(reader, "py");
            TTreeReaderValue<std::vector<float>> pz(reader, "pz");
            TTreeReaderValue<std::vector<float>> energy(reader, "energy");
            TTreeReaderValue<std::vector<char>> isFinal(reader, "isFinal");
            Long64_t entries = 0;
            long double sumThrust = 0.0;
            long double sumVisibleEnergy = 0.0;
            long double sumNISRPhotons = 0.0;
            long double sumISRPhotonEnergy = 0.0;
            while (reader.Next()) {
                ++entries;
                sumThrust += *thrust;
                sumVisibleEnergy += acceptedVisibleEnergy(*pdgId, *px, *py, *pz, *energy, *isFinal);
                sumNISRPhotons += *nISRPhotons;
                sumISRPhotonEnergy += *totalISRPhotonEnergy;
            }
            const long double denom = entries > 0 ? entries : 1;
            out << s.label << "," << modeFile.first << "," << t->GetEntries() << ","
                << static_cast<double>(sumThrust / denom) << ","
                << static_cast<double>(sumVisibleEnergy / denom) << ","
                << static_cast<double>(sumNISRPhotons / denom) << ","
                << static_cast<double>(sumISRPhotonEnergy / denom) << ","
                << modeFile.second << "\n";
        }
    }
}

void plot_real_isr_results()
{
    if (const char* envDir = gSystem->Getenv("ISR_REAL_DIR")) {
        gRealDir = envDir;
        gOutDir = joinPath(gRealDir, "results");
        gDiagDir = joinPath(gRealDir, "endpoint_diagnostics");
    }
    if (const char* envDiag = gSystem->Getenv("ISR_ENDPOINT_DIAG_DIR")) {
        gDiagDir = envDiag;
    }
    if (const char* envOut = gSystem->Getenv("ISR_OUT_DIR")) {
        gOutDir = envOut;
    }
    gSystem->mkdir(gOutDir.c_str(), true);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    const int colorBlue = TColor::GetColor("#0072B2");
    const int colorVermillion = TColor::GetColor("#D55E00");
    const int colorGreen = TColor::GetColor("#009E73");
    const int colorPurple = TColor::GetColor("#CC79A7");
    const int colorOrange = TColor::GetColor("#E69F00");

    std::vector<SampleDef> samples = {
        {"Herwig 7.3.0 QED shower", "Herwig730_QEDshower", joinPath(gRealDir, "mc_Herwig730_ISR_OFF.root"), joinPath(gRealDir, "mc_Herwig730_QEDshower.root"), "QEDshower", colorBlue, 20, -0.0018},
        {"PYTHIA 8.315", "Pythia8315", joinPath(gRealDir, "mc_Pythia8315_ISR_OFF.root"), joinPath(gRealDir, "mc_Pythia8315_ISR_ON.root"), "ISR ON", colorVermillion, 21, -0.0009},
        {"Sherpa 3.0.3 PDFESherpa", "Sherpa303", joinPath(gRealDir, "mc_Sherpa303_ISR_OFF.root"), joinPath(gRealDir, "mc_Sherpa303_ISR_ON.root"), "ISR ON (PDFESherpa)", colorGreen, 23, 0.0000},
        {"Sherpa 3.0.3 YFS", "Sherpa303_YFS", joinPath(gRealDir, "mc_Sherpa303_ISR_OFF.root"), joinPath(gRealDir, "mc_Sherpa303_ISR_YFS.root"), "ISR ON (YFS)", colorPurple, 34, 0.0009},
        {"PYTHIA 8.315 (Vincia)", "Pythia8315_Vincia", joinPath(gRealDir, "mc_Pythia8315_Vincia_ISR_OFF.root"), joinPath(gRealDir, "mc_Pythia8315_Vincia_ISR_ON.root"), "ISR ON", colorOrange, 22, 0.0018}
    };

    drawISRCorrection(samples);
    drawThrustVsAleph(samples);
    drawThrustOffVsAleph(samples);
    drawThrustVsAlephRatioPanels(samples, true);
    drawThrustVsAlephRatioPanels(samples, false);
    if (gSystem->Getenv("ISR_THRUST_ONLY")) return;
    drawISRPhotonSpectra(samples);
    drawVisibleEnergy(samples);
    if (!gSystem->Getenv("ISR_SKIP_STATS")) writeStats(samples);
}

void plot_real_thrust_figures()
{
    gSystem->Setenv("ISR_THRUST_ONLY", "1");
    plot_real_isr_results();
}

void plot_real_aleph_ratio_figures()
{
    if (const char* envDir = gSystem->Getenv("ISR_REAL_DIR")) {
        gRealDir = envDir;
        gOutDir = joinPath(gRealDir, "results");
        gDiagDir = joinPath(gRealDir, "endpoint_diagnostics");
    }
    if (const char* envDiag = gSystem->Getenv("ISR_ENDPOINT_DIAG_DIR")) {
        gDiagDir = envDiag;
    }
    if (const char* envOut = gSystem->Getenv("ISR_OUT_DIR")) {
        gOutDir = envOut;
    }
    gSystem->mkdir(gOutDir.c_str(), true);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    const int colorBlue = TColor::GetColor("#0072B2");
    const int colorVermillion = TColor::GetColor("#D55E00");
    const int colorGreen = TColor::GetColor("#009E73");
    const int colorPurple = TColor::GetColor("#CC79A7");
    const int colorOrange = TColor::GetColor("#E69F00");

    std::vector<SampleDef> samples = {
        {"Herwig 7.3.0 QED shower", "Herwig730_QEDshower", joinPath(gRealDir, "mc_Herwig730_ISR_OFF.root"), joinPath(gRealDir, "mc_Herwig730_QEDshower.root"), "QEDshower", colorBlue, 20, -0.0018},
        {"PYTHIA 8.315", "Pythia8315", joinPath(gRealDir, "mc_Pythia8315_ISR_OFF.root"), joinPath(gRealDir, "mc_Pythia8315_ISR_ON.root"), "ISR ON", colorVermillion, 21, -0.0009},
        {"Sherpa 3.0.3 PDFESherpa", "Sherpa303", joinPath(gRealDir, "mc_Sherpa303_ISR_OFF.root"), joinPath(gRealDir, "mc_Sherpa303_ISR_ON.root"), "ISR ON (PDFESherpa)", colorGreen, 23, 0.0000},
        {"Sherpa 3.0.3 YFS", "Sherpa303_YFS", joinPath(gRealDir, "mc_Sherpa303_ISR_OFF.root"), joinPath(gRealDir, "mc_Sherpa303_ISR_YFS.root"), "ISR ON (YFS)", colorPurple, 34, 0.0009},
        {"PYTHIA 8.315 (Vincia)", "Pythia8315_Vincia", joinPath(gRealDir, "mc_Pythia8315_Vincia_ISR_OFF.root"), joinPath(gRealDir, "mc_Pythia8315_Vincia_ISR_ON.root"), "ISR ON", colorOrange, 22, 0.0018}
    };

    drawThrustVsAlephRatioPanels(samples, true);
    drawThrustVsAlephRatioPanels(samples, false);
}

void plot_real_isr_off_thrust_vs_aleph()
{
    if (const char* envDir = gSystem->Getenv("ISR_REAL_DIR")) {
        gRealDir = envDir;
        gOutDir = joinPath(gRealDir, "results");
        gDiagDir = joinPath(gRealDir, "endpoint_diagnostics");
    }
    if (const char* envDiag = gSystem->Getenv("ISR_ENDPOINT_DIAG_DIR")) {
        gDiagDir = envDiag;
    }
    if (const char* envOut = gSystem->Getenv("ISR_OUT_DIR")) {
        gOutDir = envOut;
    }
    gSystem->mkdir(gOutDir.c_str(), true);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    const int colorBlue = TColor::GetColor("#0072B2");
    const int colorVermillion = TColor::GetColor("#D55E00");
    const int colorGreen = TColor::GetColor("#009E73");
    const int colorPurple = TColor::GetColor("#CC79A7");
    const int colorOrange = TColor::GetColor("#E69F00");

    std::vector<SampleDef> samples = {
        {"Herwig 7.3.0 QED shower", "Herwig730_QEDshower", joinPath(gRealDir, "mc_Herwig730_ISR_OFF.root"), joinPath(gRealDir, "mc_Herwig730_QEDshower.root"), "QEDshower", colorBlue, 20, -0.0018},
        {"PYTHIA 8.315", "Pythia8315", joinPath(gRealDir, "mc_Pythia8315_ISR_OFF.root"), joinPath(gRealDir, "mc_Pythia8315_ISR_ON.root"), "ISR ON", colorVermillion, 21, -0.0009},
        {"Sherpa 3.0.3 PDFESherpa", "Sherpa303", joinPath(gRealDir, "mc_Sherpa303_ISR_OFF.root"), joinPath(gRealDir, "mc_Sherpa303_ISR_ON.root"), "ISR ON (PDFESherpa)", colorGreen, 23, 0.0000},
        {"Sherpa 3.0.3 YFS", "Sherpa303_YFS", joinPath(gRealDir, "mc_Sherpa303_ISR_OFF.root"), joinPath(gRealDir, "mc_Sherpa303_ISR_YFS.root"), "ISR ON (YFS)", colorPurple, 34, 0.0009},
        {"PYTHIA 8.315 (Vincia)", "Pythia8315_Vincia", joinPath(gRealDir, "mc_Pythia8315_Vincia_ISR_OFF.root"), joinPath(gRealDir, "mc_Pythia8315_Vincia_ISR_ON.root"), "ISR ON", colorOrange, 22, 0.0018}
    };

    drawThrustOffVsAleph(samples);
}
