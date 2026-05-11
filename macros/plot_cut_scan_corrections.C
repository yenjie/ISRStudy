#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH1D.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLine.h>
#include <TPad.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <algorithm>
#include <cmath>
#include <cctype>
#include <cstdio>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

namespace {

struct ScanSample {
    std::string label;
    std::string tag;
    std::string offFile;
    std::string onFile;
    int color;
    int marker;
    double xOffset;
};

enum class Var {
    SPrime,
    Mvis,
    Evis,
    BetaZ,
    LeadingISR,
    AbsCosThetaThrust
};

struct EventValues {
    double thrust = 0.0;
    double sPrime = 0.0;
    double mvis = 0.0;
    double evis = 0.0;
    double betaZ = 0.0;
    double leadingISR = -1.0;
    double absCosThetaThrust = 0.0;
};

struct Condition {
    Var var;
    bool greaterThan;
    double threshold;
    std::string tag;
    std::string label;
};

struct CutScenario {
    std::string tag;
    std::string label;
    bool cumulative = false;
    std::vector<Condition> conditions;
};

struct HistStats {
    TH1D* hist = nullptr;
    long long total = 0;
    long long selected = 0;
};

struct RatioPoint {
    double low = 0.0;
    double high = 0.0;
    double center = 0.0;
    double ratio = 0.0;
    double error = 0.0;
    double off = 0.0;
    double on = 0.0;
    bool valid = false;
};

std::string gCutRealDir = "/data2/yjlee/ISRsample/real_3M_20260511";
std::string gCutDiagDir = "/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics";
std::string gCutOutDir = "/data2/yjlee/ISRsample/real_3M_20260511/results/cut_scans";
constexpr double kSqrtS = 91.1876;

std::string joinPath(const std::string& dir, const std::string& name)
{
    return dir + "/" + name;
}

std::string safeName(std::string s)
{
    for (char& c : s) {
        if (!std::isalnum(static_cast<unsigned char>(c)) && c != '_') c = '_';
    }
    return s;
}

std::string csvQuote(const std::string& s)
{
    std::string out = "\"";
    for (char c : s) {
        if (c == '"') out += "\"\"";
        else out += c;
    }
    out += "\"";
    return out;
}

double valueFor(const EventValues& e, Var var)
{
    switch (var) {
    case Var::SPrime: return e.sPrime;
    case Var::Mvis: return e.mvis;
    case Var::Evis: return e.evis;
    case Var::BetaZ: return e.betaZ;
    case Var::LeadingISR: return e.leadingISR;
    case Var::AbsCosThetaThrust: return e.absCosThetaThrust;
    }
    return 0.0;
}

bool passCondition(const EventValues& e, const Condition& c)
{
    const double v = valueFor(e, c.var);
    return c.greaterThan ? (v > c.threshold) : (v < c.threshold);
}

bool passScenario(const EventValues& e, const CutScenario& scenario)
{
    for (const Condition& c : scenario.conditions) {
        if (!passCondition(e, c)) return false;
    }
    return true;
}

std::vector<double> thrustEdges()
{
    std::vector<double> edges;
    std::string alephCsv = "/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/ALEPH/HEPData-ins636645-v1-Table_54.csv";
    if (const char* envAleph = gSystem->Getenv("ALEPH_THRUST_CSV")) alephCsv = envAleph;
    std::ifstream in(alephCsv.c_str());
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

std::vector<Condition> baseConditions()
{
    return {
        {Var::SPrime, true, 0.90, "sPrime_gt_0p90", "s'_{vis}/s > 0.90"},
        {Var::SPrime, true, 0.95, "sPrime_gt_0p95", "s'_{vis}/s > 0.95"},
        {Var::SPrime, true, 0.98, "sPrime_gt_0p98", "s'_{vis}/s > 0.98"},
        {Var::Mvis, true, 80.0, "Mvis_gt_80", "M_{vis} > 80 GeV"},
        {Var::Mvis, true, 85.0, "Mvis_gt_85", "M_{vis} > 85 GeV"},
        {Var::Mvis, true, 88.0, "Mvis_gt_88", "M_{vis} > 88 GeV"},
        {Var::Evis, true, 0.80 * kSqrtS, "Evis_gt_0p80sqrtS", "E_{vis} > 0.80 #sqrt{s}"},
        {Var::Evis, true, 0.90 * kSqrtS, "Evis_gt_0p90sqrtS", "E_{vis} > 0.90 #sqrt{s}"},
        {Var::BetaZ, false, 0.30, "betaZ_lt_0p30", "#beta_{z,vis} < 0.30"},
        {Var::BetaZ, false, 0.20, "betaZ_lt_0p20", "#beta_{z,vis} < 0.20"},
        {Var::BetaZ, false, 0.10, "betaZ_lt_0p10", "#beta_{z,vis} < 0.10"},
        {Var::LeadingISR, false, 5.0, "leadingISR_lt_5", "leading ISR photon E < 5 GeV"},
        {Var::LeadingISR, false, 2.0, "leadingISR_lt_2", "leading ISR photon E < 2 GeV"},
        {Var::AbsCosThetaThrust, false, 0.90, "absCosThetaThrust_lt_0p9", "|cos#theta_{thrust}| < 0.9"},
        {Var::AbsCosThetaThrust, false, 0.70, "absCosThetaThrust_lt_0p7", "|cos#theta_{thrust}| < 0.7"}
    };
}

std::vector<CutScenario> buildScenarios()
{
    const std::vector<Condition> cuts = baseConditions();
    std::vector<CutScenario> scenarios;
    scenarios.push_back({"baseline_no_cut", "baseline: no additional cut", false, {}});
    for (const Condition& c : cuts) {
        scenarios.push_back({"ind_" + c.tag, "individual: " + c.label, false, {c}});
    }
    std::vector<Condition> running;
    for (size_t i = 0; i < cuts.size(); ++i) {
        running.push_back(cuts[i]);
        char buf[16];
        std::snprintf(buf, sizeof(buf), "%02zu", i + 1);
        scenarios.push_back({"cum_" + std::string(buf) + "_" + cuts[i].tag,
                             "cumulative through: " + cuts[i].label,
                             true,
                             running});
    }
    return scenarios;
}

std::vector<HistStats> fillScenarios(const std::string& filename,
                                     const std::string& namePrefix,
                                     const std::vector<CutScenario>& scenarios,
                                     const std::vector<double>& edges)
{
    std::vector<HistStats> out(scenarios.size());
    for (size_t i = 0; i < scenarios.size(); ++i) {
        const std::string hname = "h_" + safeName(namePrefix + "_" + scenarios[i].tag);
        out[i].hist = new TH1D(hname.c_str(), "", static_cast<int>(edges.size()) - 1, edges.data());
        out[i].hist->Sumw2();
        out[i].hist->SetDirectory(nullptr);
    }

    TFile f(filename.c_str(), "READ");
    if (f.IsZombie()) {
        std::cerr << "Cannot open " << filename << std::endl;
        return out;
    }
    TTree* t = static_cast<TTree*>(f.Get("EndpointDiagnostics"));
    if (!t) {
        std::cerr << "Missing EndpointDiagnostics tree in " << filename << std::endl;
        return out;
    }

    TTreeReader reader(t);
    TTreeReaderValue<double> thrust(reader, "T_lab_allFinal_including_ISR_photons");
    TTreeReaderValue<double> sPrime(reader, "sPrime_vis_over_s");
    TTreeReaderValue<double> mvis(reader, "Mvis_excluding_ISR_photons");
    TTreeReaderValue<double> evis(reader, "Evis_excluding_ISR_photons");
    TTreeReaderValue<double> betaZ(reader, "betaZ_vis");
    TTreeReaderValue<double> leadingISR(reader, "leading_ISR_photon_energy");
    TTreeReaderValue<double> absCosTheta(reader, "absCosTheta_thrust_lab_allFinal_including_ISR_photons");

    while (reader.Next()) {
        EventValues e;
        e.thrust = *thrust;
        e.sPrime = *sPrime;
        e.mvis = *mvis;
        e.evis = *evis;
        e.betaZ = *betaZ;
        e.leadingISR = *leadingISR;
        e.absCosThetaThrust = *absCosTheta;

        for (size_t i = 0; i < scenarios.size(); ++i) {
            ++out[i].total;
            if (!passScenario(e, scenarios[i])) continue;
            ++out[i].selected;
            out[i].hist->Fill(e.thrust);
        }
    }
    return out;
}

RatioPoint ratioForBin(const TH1D* hOff, const TH1D* hOn, int ib)
{
    RatioPoint p;
    p.low = hOff->GetBinLowEdge(ib);
    p.high = hOff->GetBinLowEdge(ib + 1);
    p.center = hOff->GetBinCenter(ib);
    p.off = hOff->GetBinContent(ib);
    p.on = hOn->GetBinContent(ib);
    if (p.off > 0.0 && p.on > 0.0 && p.off >= 10.0 && p.on >= 10.0) {
        p.valid = true;
        p.ratio = p.off / p.on;
        p.error = p.ratio * std::sqrt(1.0 / p.off + 1.0 / p.on);
    }
    return p;
}

std::vector<RatioPoint> ratioPoints(const HistStats& off, const HistStats& on)
{
    std::vector<RatioPoint> points;
    for (int ib = 1; ib <= off.hist->GetNbinsX(); ++ib) points.push_back(ratioForBin(off.hist, on.hist, ib));
    return points;
}

TGraphErrors* makeGraph(const std::vector<RatioPoint>& points, const ScanSample& sample)
{
    std::vector<double> x, y, ex, ey;
    for (const RatioPoint& p : points) {
        if (!p.valid) continue;
        x.push_back(p.center + sample.xOffset);
        y.push_back(p.ratio);
        ex.push_back(0.0);
        ey.push_back(p.error);
    }
    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetLineColor(sample.color);
    g->SetMarkerColor(sample.color);
    g->SetMarkerStyle(sample.marker);
    g->SetMarkerSize(0.62);
    g->SetLineWidth(2);
    return g;
}

TGraphErrors* makeDoubleRatioGraph(const std::vector<RatioPoint>& points,
                                   const std::vector<RatioPoint>& ref,
                                   const ScanSample& sample)
{
    std::vector<double> x, y, ex, ey;
    const size_t n = std::min(points.size(), ref.size());
    for (size_t i = 0; i < n; ++i) {
        if (!points[i].valid || !ref[i].valid || ref[i].ratio <= 0.0) continue;
        const double r = points[i].ratio / ref[i].ratio;
        const double rel = std::sqrt(std::pow(points[i].error / points[i].ratio, 2) +
                                     std::pow(ref[i].error / ref[i].ratio, 2));
        x.push_back(points[i].center + sample.xOffset);
        y.push_back(r);
        ex.push_back(0.0);
        ey.push_back(r * rel);
    }
    TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), x.data(), y.data(), ex.data(), ey.data());
    g->SetLineColor(sample.color);
    g->SetMarkerColor(sample.color);
    g->SetMarkerStyle(sample.marker);
    g->SetMarkerSize(0.58);
    g->SetLineWidth(2);
    return g;
}

void drawScenario(const CutScenario& scenario,
                  const std::vector<ScanSample>& samples,
                  const std::vector<std::vector<RatioPoint>>& ratios,
                  const std::vector<double>& edges)
{
    const double xmin = edges.front();
    const double xmax = edges.back();
    TCanvas* c = new TCanvas(("c_cutscan_" + scenario.tag).c_str(), "cut scan", 1000, 900);
    TPad* top = new TPad(("top_" + scenario.tag).c_str(), "top", 0.0, 0.30, 1.0, 1.0);
    TPad* bot = new TPad(("bot_" + scenario.tag).c_str(), "bottom", 0.0, 0.0, 1.0, 0.30);
    top->SetLeftMargin(0.135); top->SetRightMargin(0.035); top->SetBottomMargin(0.020); top->SetTopMargin(0.075); top->SetTicks(1, 1);
    bot->SetLeftMargin(0.135); bot->SetRightMargin(0.035); bot->SetBottomMargin(0.315); bot->SetTopMargin(0.025); bot->SetTicks(1, 1);
    top->Draw();
    bot->Draw();

    top->cd();
    TH1D* frame = new TH1D(("frame_" + scenario.tag).c_str(), "", 100, xmin, xmax);
    frame->SetMinimum(0.25);
    frame->SetMaximum(2.45);
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitle("C_{ISR} = OFF / ON");
    frame->GetYaxis()->SetTitleSize(0.050);
    frame->GetYaxis()->SetTitleOffset(1.10);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->Draw("AXIS");
    TLine* one = new TLine(xmin, 1.0, xmax, 1.0);
    one->SetLineColor(kGray + 2);
    one->SetLineStyle(2);
    one->Draw();

    std::vector<TGraphErrors*> graphs;
    for (size_t i = 0; i < samples.size(); ++i) {
        TGraphErrors* g = makeGraph(ratios[i], samples[i]);
        graphs.push_back(g);
        g->Draw("LP SAME");
    }

    TLatex lat;
    lat.SetNDC();
    lat.SetTextFont(42);
    lat.SetTextSize(0.034);
    lat.DrawLatex(0.145, 0.900, scenario.label.c_str());
    lat.SetTextSize(0.026);
    lat.DrawLatex(0.145, 0.855, "T = all stable final particles, including #nu and tagged ISR photons");

    TLegend* leg = new TLegend(0.48, 0.55, 0.91, 0.86);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.027);
    leg->SetNColumns(1);
    for (size_t i = 0; i < samples.size(); ++i) leg->AddEntry(graphs[i], samples[i].label.c_str(), "lp");
    leg->Draw();
    frame->Draw("AXIS SAME");

    bot->cd();
    TH1D* rframe = new TH1D(("rframe_" + scenario.tag).c_str(), "", 100, xmin, xmax);
    rframe->SetMinimum(0.25);
    rframe->SetMaximum(2.45);
    rframe->GetXaxis()->SetTitle("Thrust T");
    rframe->GetYaxis()->SetTitle("C_{ISR}/PYTHIA");
    rframe->GetXaxis()->SetTitleSize(0.100);
    rframe->GetXaxis()->SetLabelSize(0.080);
    rframe->GetYaxis()->SetTitleSize(0.078);
    rframe->GetYaxis()->SetTitleOffset(0.72);
    rframe->GetYaxis()->SetLabelSize(0.070);
    rframe->GetYaxis()->SetNdivisions(505);
    rframe->Draw("AXIS");
    TLine* oneBot = new TLine(xmin, 1.0, xmax, 1.0);
    oneBot->SetLineColor(kGray + 2);
    oneBot->SetLineStyle(2);
    oneBot->Draw();
    const std::vector<RatioPoint>& ref = ratios[1];
    for (size_t i = 0; i < samples.size(); ++i) {
        TGraphErrors* g = makeDoubleRatioGraph(ratios[i], ref, samples[i]);
        g->Draw("LP SAME");
    }
    rframe->Draw("AXIS SAME");

    c->SaveAs(joinPath(gCutOutDir, "cutscan_" + safeName(scenario.tag) + ".png").c_str());
    c->SaveAs(joinPath(gCutOutDir, "cutscan_" + safeName(scenario.tag) + ".pdf").c_str());
    delete c;
}

} // namespace

void plot_cut_scan_corrections()
{
    if (const char* envReal = gSystem->Getenv("ISR_REAL_DIR")) {
        gCutRealDir = envReal;
        gCutDiagDir = joinPath(gCutRealDir, "endpoint_diagnostics");
        gCutOutDir = joinPath(joinPath(gCutRealDir, "results"), "cut_scans");
    }
    if (const char* envDiag = gSystem->Getenv("ISR_ENDPOINT_DIAG_DIR")) gCutDiagDir = envDiag;
    if (const char* envOut = gSystem->Getenv("ISR_CUTSCAN_OUT_DIR")) gCutOutDir = envOut;
    gSystem->mkdir(gCutOutDir.c_str(), true);
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);

    const int colorBlue = TColor::GetColor("#0072B2");
    const int colorVermillion = TColor::GetColor("#D55E00");
    const int colorGreen = TColor::GetColor("#009E73");
    const int colorPurple = TColor::GetColor("#CC79A7");
    const int colorOrange = TColor::GetColor("#E69F00");

    const std::vector<ScanSample> samples = {
        {"Herwig 7.3.0 QED shower", "Herwig730_QEDshower",
         joinPath(gCutDiagDir, "endpoint_diagnostics_Herwig730_OFF.root"),
         joinPath(gCutDiagDir, "endpoint_diagnostics_Herwig730_QEDshower.root"), colorBlue, 20, -0.0018},
        {"PYTHIA 8.315", "Pythia8315",
         joinPath(gCutDiagDir, "endpoint_diagnostics_Pythia8315_OFF.root"),
         joinPath(gCutDiagDir, "endpoint_diagnostics_Pythia8315.root"), colorVermillion, 21, -0.0009},
        {"Sherpa 3.0.3 PDFESherpa", "Sherpa303",
         joinPath(gCutDiagDir, "endpoint_diagnostics_Sherpa303_OFF.root"),
         joinPath(gCutDiagDir, "endpoint_diagnostics_Sherpa303.root"), colorGreen, 23, 0.0000},
        {"Sherpa 3.0.3 YFS", "Sherpa303_YFS",
         joinPath(gCutDiagDir, "endpoint_diagnostics_Sherpa303_OFF.root"),
         joinPath(gCutDiagDir, "endpoint_diagnostics_Sherpa303_YFS.root"), colorPurple, 34, 0.0009},
        {"PYTHIA 8.315 (Vincia)", "Pythia8315_Vincia",
         joinPath(gCutDiagDir, "endpoint_diagnostics_Pythia8315_Vincia_OFF.root"),
         joinPath(gCutDiagDir, "endpoint_diagnostics_Pythia8315_Vincia.root"), colorOrange, 22, 0.0018}
    };

    const std::vector<CutScenario> scenarios = buildScenarios();
    const std::vector<double> edges = thrustEdges();

    std::vector<std::vector<HistStats>> offStats(samples.size());
    std::vector<std::vector<HistStats>> onStats(samples.size());
    for (size_t is = 0; is < samples.size(); ++is) {
        std::cout << "Filling OFF " << samples[is].label << std::endl;
        offStats[is] = fillScenarios(samples[is].offFile, samples[is].tag + "_off", scenarios, edges);
        std::cout << "Filling ON " << samples[is].label << std::endl;
        onStats[is] = fillScenarios(samples[is].onFile, samples[is].tag + "_on", scenarios, edges);
    }

    std::ofstream summary(joinPath(gCutOutDir, "cutscan_summary.csv").c_str());
    std::ofstream allBins(joinPath(gCutOutDir, "cutscan_all_bins.csv").c_str());
    std::ofstream sherpa(joinPath(gCutOutDir, "cutscan_sherpa_endpoint_summary.csv").c_str());
    summary << std::setprecision(10);
    allBins << std::setprecision(10);
    sherpa << std::setprecision(10);
    summary << "scenario,scenario_type,cut_label,sample,off_selected,off_total,off_survival,on_selected,on_total,on_survival,"
            << "prev_bin_low,prev_bin_high,prev_ratio,prev_error,prev_off_count,prev_on_count,"
            << "last_bin_low,last_bin_high,last_ratio,last_error,last_off_count,last_on_count\n";
    allBins << "scenario,scenario_type,cut_label,sample,bin_low,bin_high,off_count,on_count,ratio,error,double_ratio_to_pythia,double_ratio_error\n";
    sherpa << "scenario,scenario_type,cut_label,sample,prev_bin_ratio,prev_bin_error,last_bin_ratio,last_bin_error,"
           << "off_survival,on_survival\n";

    for (size_t ic = 0; ic < scenarios.size(); ++ic) {
        std::vector<std::vector<RatioPoint>> ratios(samples.size());
        for (size_t is = 0; is < samples.size(); ++is) ratios[is] = ratioPoints(offStats[is][ic], onStats[is][ic]);
        drawScenario(scenarios[ic], samples, ratios, edges);

        const std::string scenarioType = scenarios[ic].tag == "baseline_no_cut" ? "baseline" : (scenarios[ic].cumulative ? "cumulative" : "individual");
        const std::vector<RatioPoint>& ref = ratios[1];
        for (size_t is = 0; is < samples.size(); ++is) {
            const RatioPoint& prev = ratios[is][ratios[is].size() - 2];
            const RatioPoint& last = ratios[is][ratios[is].size() - 1];
            const double offSurv = offStats[is][ic].total > 0 ? static_cast<double>(offStats[is][ic].selected) / offStats[is][ic].total : 0.0;
            const double onSurv = onStats[is][ic].total > 0 ? static_cast<double>(onStats[is][ic].selected) / onStats[is][ic].total : 0.0;
            summary << csvQuote(scenarios[ic].tag) << "," << csvQuote(scenarioType) << "," << csvQuote(scenarios[ic].label) << ","
                    << csvQuote(samples[is].label) << ","
                    << offStats[is][ic].selected << "," << offStats[is][ic].total << "," << offSurv << ","
                    << onStats[is][ic].selected << "," << onStats[is][ic].total << "," << onSurv << ","
                    << prev.low << "," << prev.high << "," << (prev.valid ? prev.ratio : 0.0) << "," << (prev.valid ? prev.error : 0.0) << ","
                    << prev.off << "," << prev.on << ","
                    << last.low << "," << last.high << "," << (last.valid ? last.ratio : 0.0) << "," << (last.valid ? last.error : 0.0) << ","
                    << last.off << "," << last.on << "\n";

            if (samples[is].tag.find("Sherpa") != std::string::npos) {
                sherpa << csvQuote(scenarios[ic].tag) << "," << csvQuote(scenarioType) << "," << csvQuote(scenarios[ic].label) << ","
                       << csvQuote(samples[is].label) << ","
                       << (prev.valid ? prev.ratio : 0.0) << "," << (prev.valid ? prev.error : 0.0) << ","
                       << (last.valid ? last.ratio : 0.0) << "," << (last.valid ? last.error : 0.0) << ","
                       << offSurv << "," << onSurv << "\n";
            }

            for (size_t ib = 0; ib < ratios[is].size(); ++ib) {
                const RatioPoint& p = ratios[is][ib];
                if (!p.valid) continue;
                double dr = 0.0;
                double der = 0.0;
                if (ib < ref.size() && ref[ib].valid && ref[ib].ratio > 0.0) {
                    dr = p.ratio / ref[ib].ratio;
                    der = dr * std::sqrt(std::pow(p.error / p.ratio, 2) + std::pow(ref[ib].error / ref[ib].ratio, 2));
                }
                allBins << csvQuote(scenarios[ic].tag) << "," << csvQuote(scenarioType) << "," << csvQuote(scenarios[ic].label) << ","
                        << csvQuote(samples[is].label) << "," << p.low << "," << p.high << ","
                        << p.off << "," << p.on << "," << p.ratio << "," << p.error << ","
                        << dr << "," << der << "\n";
            }
        }
    }
}
