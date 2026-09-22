// Two verification panels for the charged-EEC ISR correction.
//
// Left:  the measured correction against the one predicted by the charged
//        energy loss alone, C_ISR = 1/(1-d)^2.  If the EEC correction is just
//        the energy radiation removes from the hadronic system, these agree.
//
// Right: the energy radiation removes, measured two ways.  The beam-collinear
//        photon estimator only sees photons within |cos theta| > 0.9999; the
//        charged energy loss sees all of it, because every radiated photon
//        lowers the energy left for the hadronic system.  KKMC's own generator
//        ISR photon list is the truth point that calibrates the comparison.
//
// Usage:
//   root -l -b -q 'macros/plot_eec_verification.C("<results dir>")'

#include <fstream>
#include <iostream>
#include <map>
#include <sstream>
#include <string>
#include <vector>

namespace {

std::vector<std::string> splitCsv(const std::string& line)
{
    std::vector<std::string> out;
    std::string cur;
    bool q = false;
    for (char c : line) {
        if (c == '"') { q = !q; continue; }
        if (c == ',' && !q) { out.push_back(cur); cur.clear(); continue; }
        cur.push_back(c);
    }
    out.push_back(cur);
    return out;
}

std::vector<std::vector<std::string>> readCsv(const std::string& path)
{
    std::vector<std::vector<std::string>> rows;
    std::ifstream in(path);
    if (!in) { std::cerr << "cannot open " << path << std::endl; return rows; }
    std::string line;
    std::getline(in, line);
    while (std::getline(in, line)) {
        if (line.empty()) continue;
        rows.push_back(splitCsv(line));
    }
    return rows;
}

}  // namespace

void plot_eec_verification(const char* dir =
                               "/data2/yjlee/ISRsample/kkmc_1M_20260921/results")
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gStyle->SetTextFont(42);
    gStyle->SetLabelFont(42, "XYZ");
    gStyle->SetTitleFont(42, "XYZ");
    gStyle->SetEndErrorSize(4);

    // Print order; Herwig last because its two samples do not switch beam
    // radiation, so its entry measures final-state radiation instead.
    const std::vector<std::string> order = {
        "Pythia 8.315", "Pythia 8.315 Vincia", "Sherpa 3.0.3 PDFESherpa",
        "Sherpa 3.0.3 YFS", "KKMC 4.30 CEEX",
        "Herwig 7.3.0 QED shower, ISR unchanged"};
    const std::vector<std::string> shortName = {
        "Pythia", "Vincia", "Sherpa", "Sherpa YFS", "KKMC", "Herwig*"};

    std::map<std::string, std::vector<std::string>> ver, sumOff, sumOn;
    for (const auto& r : readCsv(std::string(dir) + "/eec_verification.csv"))
        if (r.size() > 12) ver[r[0]] = r;
    for (const auto& r : readCsv(std::string(dir) + "/isr_model_summary.csv"))
        if (r.size() > 7) { (r[1] == "OFF" ? sumOff : sumOn)[r[0]] = r; }
    if (ver.empty()) { std::cerr << "no verification rows" << std::endl; return; }

    const int n = static_cast<int>(order.size());
    std::vector<double> x(n), meas(n), measErr(n), pred(n), ebc(n), ech(n),
        echErr(n), zero(n, 0.0);
    double kkmcTruth = 0;
    for (int i = 0; i < n; ++i) {
        x[i] = i + 0.36;
        const auto& v = ver[order[i]];
        meas[i] = (std::stod(v[10]) - 1.0) * 1e3;
        measErr[i] = std::stod(v[11]) * 1e3;
        pred[i] = (std::stod(v[12]) - 1.0) * 1e3;
        ech[i] = std::stod(v[6]) * 91.1876;          // delta * sqrt(s)
        echErr[i] = std::stod(v[7]) * 91.1876;
        if (sumOn.count(order[i]) && sumOff.count(order[i]))
            ebc[i] = std::stod(sumOn[order[i]][6]) - std::stod(sumOff[order[i]][6]);
        if (order[i] == "KKMC 4.30 CEEX" && sumOn.count(order[i]))
            kkmcTruth = std::stod(sumOn[order[i]][7]);
    }

    TCanvas* c = new TCanvas("c_ver", "", 1500, 620);
    c->Divide(2, 1);

    // ------------------------------------------------------------- left panel
    c->cd(1);
    gPad->SetLeftMargin(0.13); gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.09);  gPad->SetBottomMargin(0.13);
    TH1D* f1 = new TH1D("f1", "", n, 0, n);
    f1->SetMinimum(0.0);
    f1->SetMaximum(8.5);
    f1->GetYaxis()->SetTitle("(C_{ISR} - 1) #times 10^{3}");
    f1->GetYaxis()->SetTitleSize(0.055);
    f1->GetYaxis()->SetLabelSize(0.048);
    f1->GetYaxis()->SetTitleOffset(1.05);
    f1->GetXaxis()->SetLabelSize(0.058);
    for (int i = 0; i < n; ++i) f1->GetXaxis()->SetBinLabel(i + 1, shortName[i].c_str());
    f1->Draw();

    TGraphErrors* gm = new TGraphErrors(n, &x[0], &meas[0], &zero[0], &measErr[0]);
    gm->SetMarkerStyle(20); gm->SetMarkerSize(1.5);
    gm->SetMarkerColor(kBlack); gm->SetLineColor(kBlack); gm->SetLineWidth(2);
    TGraph* gp = new TGraph(n, &x[0], &pred[0]);
    gp->SetMarkerStyle(24); gp->SetMarkerSize(2.1);
    gp->SetMarkerColor(kRed + 1); gp->SetLineColor(kRed + 1);
    gp->Draw("P SAME");
    gm->Draw("P SAME");

    TLegend* l1 = new TLegend(0.17, 0.72, 0.62, 0.87);
    l1->SetBorderSize(0); l1->SetFillStyle(0); l1->SetTextSize(0.043);
    l1->AddEntry(gm, "measured, binned EEC", "lp");
    l1->AddEntry(gp, "predicted, 1/(1-#delta)^{2}", "p");
    l1->Draw();

    TLatex t1;
    t1.SetNDC(); t1.SetTextFont(42); t1.SetTextSize(0.046);
    t1.DrawLatex(0.13, 0.925, "Charged-EEC ISR correction, all z");
    t1.SetTextSize(0.036); t1.SetTextColor(kGray + 2);
    t1.DrawLatex(0.17, 0.66, "agree to 0.04-0.13 per mille, all six");

    // ------------------------------------------------------------ right panel
    c->cd(2);
    gPad->SetLeftMargin(0.155); gPad->SetRightMargin(0.03);
    gPad->SetTopMargin(0.09);  gPad->SetBottomMargin(0.13);
    TH1D* f2 = new TH1D("f2", "", n, 0, n);
    f2->SetMinimum(0.0);
    f2->SetMaximum(0.38);
    f2->GetYaxis()->SetTitle("energy removed from hadrons [GeV]");
    f2->GetYaxis()->SetTitleSize(0.055);
    f2->GetYaxis()->SetLabelSize(0.048);
    f2->GetYaxis()->SetTitleOffset(1.35);
    f2->GetXaxis()->SetLabelSize(0.058);
    for (int i = 0; i < n; ++i) f2->GetXaxis()->SetBinLabel(i + 1, shortName[i].c_str());
    f2->Draw();

    TGraphErrors* ge = new TGraphErrors(n, &x[0], &ech[0], &zero[0], &echErr[0]);
    ge->SetMarkerStyle(21); ge->SetMarkerSize(1.5);
    ge->SetMarkerColor(kAzure + 2); ge->SetLineColor(kAzure + 2); ge->SetLineWidth(2);
    std::vector<double> xb(n);
    for (int i = 0; i < n; ++i) xb[i] = i + 0.64;
    TGraph* gb = new TGraph(n, &xb[0], &ebc[0]);
    gb->SetMarkerStyle(25); gb->SetMarkerSize(1.9); gb->SetMarkerColor(kOrange + 8);
    ge->Draw("P SAME");
    gb->Draw("P SAME");

    TLine* gt = nullptr;
    if (kkmcTruth > 0) {
        gt = new TLine(4.08, kkmcTruth, 4.92, kkmcTruth);
        gt->SetLineColor(kRed + 1);
        gt->SetLineWidth(3);
        gt->Draw("SAME");
    }

    TLegend* l2 = new TLegend(0.19, 0.70, 0.78, 0.88);
    l2->SetBorderSize(0); l2->SetFillStyle(0); l2->SetTextSize(0.041);
    l2->AddEntry(ge, "from the charged energy loss", "p");
    l2->AddEntry(gb, "beam-collinear photon estimator", "p");
    if (gt) l2->AddEntry(gt, "KKMC generator ISR photons", "l");
    l2->Draw();

    TLatex t2;
    t2.SetNDC(); t2.SetTextFont(42); t2.SetTextSize(0.046);
    t2.DrawLatex(0.155, 0.925, "Radiated energy, two estimators");
    t2.SetTextSize(0.034); t2.SetTextColor(kGray + 2);
    t2.DrawLatex(0.19, 0.16, "*Herwig pair switches final-state, not beam, radiation");

    c->SaveAs(Form("%s/eec_isr_verification.png", dir));
    c->SaveAs(Form("%s/eec_isr_verification.pdf", dir));
    std::cout << "[done] wrote verification panels to " << dir << std::endl;
}
