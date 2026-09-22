// Charged-EEC ISR correction in the project's double-logarithmic style.
//
// The EEC binning is already double-logarithmic: 100 logarithmic theta_L bins
// from 0.002 to pi/2, mirror-reflected about pi/2, so both the collinear limit
// (z -> 0) and the back-to-back limit (z -> 1) are resolved.  The spectrum is
// therefore drawn against the BIN INDEX on a linear frame -- not with SetLogx,
// which would squash the back-to-back side -- and physical z values are written
// on as bin labels at the nearest index.  Only the spectrum panel uses log y.
//
// This mirrors ConfigureIndexAxis in
// ALEPH_EEC/reproduction/plotting/MakeReportPlots.C (lines ~209-240).
//
// Input is the CSV written by macros/plot_eec_isr_correction.C, so restyling
// does not require rerunning the event loop.
//
// Usage:
//   root -l -b -q 'macros/plot_eec_isr_doublelog.C("<results dir>")'

#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

constexpr int kAngleBins = 200;

struct Series {
    std::vector<double> eecOn, eecOnErr, ratio, ratioErr, zCenter;
    int color = 1;
    int marker = 20;
};

// Split a CSV line, honouring the quoted sample field.
std::vector<std::string> splitCsv(const std::string& line)
{
    std::vector<std::string> out;
    std::string cur;
    bool inQuote = false;
    for (char c : line) {
        if (c == '"') { inQuote = !inQuote; continue; }
        if (c == ',' && !inQuote) { out.push_back(cur); cur.clear(); continue; }
        cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

double toD(const std::string& s)
{
    if (s.empty()) return 0.0;
    return std::atof(s.c_str());
}

// Index of the bin whose z centre is closest to a target z.
int closestIndex(const std::vector<double>& centre, double target)
{
    int best = 0;
    double bestD = 1e30;
    for (int i = 0; i < static_cast<int>(centre.size()); ++i) {
        const double d = std::fabs(centre[i] - target);
        if (d < bestD) { bestD = d; best = i; }
    }
    return best;
}

void configureIndexAxis(TH1D* frame, const std::vector<double>& centre, bool bottom)
{
    TAxis* a = frame->GetXaxis();
    a->SetTitle("z = (1 - cos#theta_{L})/2");
    a->CenterTitle();
    a->SetTickLength(bottom ? 0.060 : 0.025);
    a->SetLabelSize(bottom ? 0.075 : 0.0);
    a->SetTitleSize(bottom ? 0.090 : 0.0);
    a->SetTitleOffset(1.02);
    if (!bottom) return;
    const double pi = std::acos(-1.0);
    const auto z = [](double t) { return (1.0 - std::cos(t)) / 2.0; };
    const std::vector<double> target = {z(0.01), z(0.1), 0.5, z(pi - 0.1), z(pi - 0.01)};
    const std::vector<std::string> label = {
        "2.5#times10^{-5}", "2.5#times10^{-3}", "1/2",
        "1-2.5#times10^{-3}", "1-2.5#times10^{-5}"};
    for (size_t i = 0; i < target.size(); ++i)
        a->SetBinLabel(closestIndex(centre, target[i]) + 1, label[i].c_str());
    a->LabelsOption("h");
}

}  // namespace

void plot_eec_isr_doublelog(const char* dir =
                                "/data2/yjlee/ISRsample/kkmc_1M_20260921/results")
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetErrorX(0);

    const std::string csvPath = std::string(dir) + "/eec_isr_correction.csv";
    std::ifstream in(csvPath);
    if (!in) { std::cerr << "cannot open " << csvPath << std::endl; return; }

    // Preserve the order the samples appear in, and give each a style.
    const std::vector<std::string> order = {
        "Pythia 8.315", "Pythia 8.315 Vincia", "Sherpa 3.0.3 PDFESherpa",
        "Sherpa 3.0.3 YFS", "Herwig 7.3.0 QED shower, ISR unchanged",
        "KKMC 4.30 CEEX"};
    const std::vector<int> colors = {kBlack, kGray + 2, kAzure + 2, kAzure + 7,
                                     kOrange + 7, kRed + 1};
    const std::vector<int> markers = {20, 24, 21, 25, 22, 29};

    std::map<std::string, Series> series;
    std::vector<double> centre;

    std::string line;
    std::getline(in, line);  // header
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        const std::vector<std::string> f = splitCsv(line);
        if (f.size() < 10) continue;
        if (f[1] == "integral") continue;
        const std::string sample = f[0];
        const double zLow = toD(f[2]), zHigh = toD(f[3]);
        Series& s = series[sample];
        s.eecOn.push_back(toD(f[6]));
        s.eecOnErr.push_back(toD(f[7]));
        s.ratio.push_back(toD(f[8]));
        s.ratioErr.push_back(toD(f[9]));
        s.zCenter.push_back(0.5 * (zLow + zHigh));
    }
    in.close();
    if (series.empty()) { std::cerr << "no rows parsed" << std::endl; return; }
    centre = series.begin()->second.zCenter;
    std::cout << "parsed " << series.size() << " samples, "
              << centre.size() << " bins" << std::endl;

    for (size_t i = 0; i < order.size(); ++i) {
        auto it = series.find(order[i]);
        if (it == series.end()) continue;
        it->second.color = colors[i];
        it->second.marker = markers[i];
    }

    TCanvas* c = new TCanvas("c_eec_dl", "", 1200, 1200);
    TPad* top = new TPad("top", "", 0, 0.34, 1, 1);
    TPad* bot = new TPad("bot", "", 0, 0, 1, 0.34);
    top->SetLeftMargin(0.13); top->SetRightMargin(0.03);
    top->SetTopMargin(0.06);  top->SetBottomMargin(0.015);
    bot->SetLeftMargin(0.13); bot->SetRightMargin(0.03);
    bot->SetTopMargin(0.015); bot->SetBottomMargin(0.30);
    top->SetLogy();
    top->Draw();
    bot->Draw();

    // ---------------------------------------------------------------- top
    top->cd();
    TH1D* ftop = new TH1D("ftop", "", kAngleBins, 0, kAngleBins);
    ftop->SetMinimum(2e-5);
    ftop->SetMaximum(0.2);
    ftop->GetYaxis()->SetTitle("EEC per event, ISR on");
    ftop->GetYaxis()->SetTitleSize(0.050);
    ftop->GetYaxis()->SetLabelSize(0.042);
    ftop->GetYaxis()->SetTitleOffset(1.20);
    configureIndexAxis(ftop, centre, false);
    ftop->Draw();

    // The EEC peaks near 2.5e-3 per event, so the whole upper band of the pad
    // is empty; the legend goes there rather than over the curves.
    TLegend* leg = new TLegend(0.17, 0.52, 0.72, 0.82);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.033);

    std::vector<TGraphErrors*> keep;
    for (const std::string& name : order) {
        auto it = series.find(name);
        if (it == series.end()) continue;
        const Series& s = it->second;
        std::vector<double> x(s.eecOn.size()), ex(s.eecOn.size(), 0.0);
        for (size_t i = 0; i < x.size(); ++i) x[i] = i + 0.5;
        TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), &x[0],
            const_cast<double*>(&s.eecOn[0]), &ex[0],
            const_cast<double*>(&s.eecOnErr[0]));
        g->SetLineColor(s.color);
        g->SetMarkerColor(s.color);
        g->SetLineWidth(2);
        g->SetMarkerStyle(s.marker);
        g->SetMarkerSize(0.7);
        g->Draw("L SAME");
        leg->AddEntry(g, name.c_str(), "l");
        keep.push_back(g);
    }
    leg->Draw();

    TLatex tx;
    tx.SetNDC();
    tx.SetTextFont(42);
    tx.SetTextSize(0.040);
    tx.DrawLatex(0.17, 0.91, "Charged EEC, particle level");
    tx.SetTextSize(0.031);
    tx.SetTextColor(kGray + 2);
    tx.DrawLatex(0.17, 0.865, "ISR study, work in progress");

    // ---------------------------------------------------------------- bottom
    bot->cd();
    TH1D* fbot = new TH1D("fbot", "", kAngleBins, 0, kAngleBins);
    fbot->SetMinimum(0.94);
    fbot->SetMaximum(1.06);
    fbot->GetYaxis()->SetTitle("C_{ISR} = OFF / ON");
    fbot->GetYaxis()->SetTitleSize(0.095);
    fbot->GetYaxis()->SetLabelSize(0.080);
    fbot->GetYaxis()->SetTitleOffset(0.62);
    fbot->GetYaxis()->SetNdivisions(505);
    configureIndexAxis(fbot, centre, true);
    fbot->Draw();

    TLine* unity = new TLine(0, 1.0, kAngleBins, 1.0);
    unity->SetLineStyle(2);
    unity->SetLineColor(kGray + 1);
    unity->Draw();

    for (const std::string& name : order) {
        auto it = series.find(name);
        if (it == series.end()) continue;
        const Series& s = it->second;
        std::vector<double> x, y, ex, ey;
        for (size_t i = 0; i < s.ratio.size(); ++i) {
            if (!(s.ratio[i] > 0)) continue;
            x.push_back(i + 0.5);
            y.push_back(s.ratio[i]);
            ex.push_back(0.0);
            ey.push_back(s.ratioErr[i]);
        }
        if (x.empty()) continue;
        TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0],
                                           &ex[0], &ey[0]);
        g->SetLineColor(s.color);
        g->SetMarkerColor(s.color);
        g->SetMarkerStyle(s.marker);
        g->SetMarkerSize(0.6);
        g->SetLineWidth(1);
        g->Draw("P SAME");
        keep.push_back(g);
    }

    c->cd();
    // ROOT writes a portrait PDF page regardless of canvas aspect, leaving white
    // space below the drawing.  The PNG is already tight; the PDF is cropped by
    // the caller.
    c->SaveAs(Form("%s/eec_isr_correction_doublelog.png", dir));
    c->SaveAs(Form("%s/eec_isr_correction_doublelog.pdf", dir));
    std::cout << "[done] wrote double-log EEC plots to " << dir << std::endl;
}
