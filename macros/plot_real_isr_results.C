#include <TCanvas.h>
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

std::string gRealDir = "/data2/yjlee/ISRsample/real_20260510";
std::string gOutDir = "/data2/yjlee/ISRsample/real_20260510/results";

std::string joinPath(const std::string& dir, const std::string& name)
{
    return dir + "/" + name;
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

TGraphErrors* makeCorrectionGraph(const SampleDef& s, const std::vector<double>& edges)
{
    TH1D* hOff = makeHistEdges(s.offFile, "h_off_" + s.tag, "thrust", edges);
    TH1D* hOn = makeHistEdges(s.onFile, "h_on_" + s.tag, "thrust", edges);
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
    g->SetMarkerSize(0.62);
    g->SetLineWidth(2);
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
    out->SetMarkerSize(0.55);
    out->SetLineWidth(2);
    return out;
}

void drawISRCorrection(const std::vector<SampleDef>& samples)
{
    std::vector<double> edges;
    for (int i = 0; i <= 32; ++i) edges.push_back(0.57 + i * (0.43 / 32.0));

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
    TH1D* frame = new TH1D("frame_corr", "", 100, 0.57, 1.0);
    frame->SetMinimum(0.55);
    frame->SetMaximum(1.42);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("C_{ISR} = ISR OFF / ISR ON");
    frame->GetYaxis()->SetTitleSize(0.055);
    frame->GetYaxis()->SetTitleOffset(1.05);
    frame->GetYaxis()->SetLabelSize(0.044);
    frame->Draw("AXIS");
    TLine* one = new TLine(0.57, 1.0, 1.0, 1.0);
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
    lat.DrawLatex(0.145, 0.815, "Real standalone generators; Herwig uses QCD vs QCD+QED shower");
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
    TH1D* rframe = new TH1D("frame_dr", "", 100, 0.57, 1.0);
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
    TLine* one2 = new TLine(0.57, 1.0, 1.0, 1.0);
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
    std::string alephCsv = "/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/ALEPH/HEPData-ins636645-v1-Table_54.csv";
    if (const char* envAleph = gSystem->Getenv("ALEPH_THRUST_CSV")) {
        alephCsv = envAleph;
    }
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
        TH1D* h = makeHistEdges(s.onFile, "h_aleph_" + s.tag, "thrust", edges);
        normalizeDensity(h);
        h->SetLineColor(s.color);
        h->SetMarkerColor(s.color);
        h->SetLineWidth(2);
        h->SetMarkerStyle(s.marker);
        h->SetMarkerSize(0.35);
        h->Draw("HIST SAME");
        leg->AddEntry(h, s.label.c_str(), "l");
    }
    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.58, 0.86, "#sqrt{s}=91.1876 GeV, ISR/QED-shower comparison");
    leg->Draw();
    frame->Draw("AXIS SAME");
    c->SaveAs(joinPath(gOutDir, "real_isr_on_thrust_vs_aleph.png").c_str());
    c->SaveAs(joinPath(gOutDir, "real_isr_on_thrust_vs_aleph.pdf").c_str());
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
        h->SetLineWidth(2);
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
        h->SetLineWidth(2);
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
        TH1D* hOff = makeHist(s.offFile, "h_evis_off_" + s.tag, "visibleEnergy", 80, 0.0, 100.0);
        TH1D* hOn = makeHist(s.onFile, "h_evis_on_" + s.tag, "visibleEnergy", 80, 0.0, 100.0);
        normalizeDensity(hOff);
        normalizeDensity(hOn);
        hOff->SetLineColor(kBlack);
        hOff->SetLineWidth(2);
        hOn->SetLineColor(s.color);
        hOn->SetLineWidth(2);
        hOff->SetMaximum(std::max(hOff->GetMaximum(), hOn->GetMaximum()) * 1.25);
        hOff->GetXaxis()->SetTitle("Visible energy E_{vis} [GeV]");
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
            auto mean = [&](const char* x) {
                t->Draw(x, "", "goff");
                return TMath::Mean(t->GetSelectedRows(), t->GetV1());
            };
            out << s.label << "," << modeFile.first << "," << t->GetEntries() << ","
                << mean("thrust") << "," << mean("visibleEnergy") << ","
                << mean("nISRPhotons") << "," << mean("totalISRPhotonEnergy") << ","
                << modeFile.second << "\n";
        }
    }
}

void plot_real_isr_results()
{
    if (const char* envDir = gSystem->Getenv("ISR_REAL_DIR")) {
        gRealDir = envDir;
        gOutDir = joinPath(gRealDir, "results");
    }
    if (const char* envOut = gSystem->Getenv("ISR_OUT_DIR")) {
        gOutDir = envOut;
    }
    gSystem->mkdir(gOutDir.c_str(), true);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    std::vector<SampleDef> samples = {
        {"Herwig 7.3.0 QED shower", "Herwig730_QEDshower", joinPath(gRealDir, "mc_Herwig730_ISR_OFF.root"), joinPath(gRealDir, "mc_Herwig730_QEDshower.root"), "QEDshower", kBlue + 2, 20, -0.0015},
        {"PYTHIA 8.315", "Pythia8315", joinPath(gRealDir, "mc_Pythia8315_ISR_OFF.root"), joinPath(gRealDir, "mc_Pythia8315_ISR_ON.root"), "ISR ON", kRed + 1, 20, -0.0005},
        {"Sherpa 3.0.3 PDFESherpa", "Sherpa303", joinPath(gRealDir, "mc_Sherpa303_ISR_OFF.root"), joinPath(gRealDir, "mc_Sherpa303_ISR_ON.root"), "ISR ON (PDFESherpa)", kMagenta + 2, 23, 0.0005},
        {"PYTHIA 8.315 (Vincia)", "Pythia8315_Vincia", joinPath(gRealDir, "mc_Pythia8315_Vincia_ISR_OFF.root"), joinPath(gRealDir, "mc_Pythia8315_Vincia_ISR_ON.root"), "ISR ON", kCyan + 2, 22, 0.0015}
    };

    drawISRCorrection(samples);
    drawThrustVsAleph(samples);
    drawISRPhotonSpectra(samples);
    drawVisibleEnergy(samples);
    writeStats(samples);
}
