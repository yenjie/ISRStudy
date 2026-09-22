// ISR correction for the charged-particle EEC.
//
//   C_ISR(z) = EEC_ISR_OFF(z) / EEC_ISR_ON(z)
//
// This is a different observable from the thrust correction: the ALEPH
// inclusive EEC analysis is built on CHARGED particles only, energy-weighted
// pairs, in z = (1 - cos theta_ij)/2.  The thrust correction cannot be
// rescaled into this one, so it is recomputed here from the same ISR ON/OFF
// generator samples.
//
// Definitions are taken from the EEC reproduction code so that the correction
// is directly usable there:
//   particles  reproduction/src/BuildThrustSlicedGeneratorEEC.cpp:112
//              (generator level: |charge| > 0, |p| > 0, E > 0; no pT and no
//              |cos theta| cut, and no highPurity, which is a detector concept)
//   pairs      unordered, j = i+1, counted once
//   coordinate z = (1 - cos theta_ij)/2
//   weight     E_i E_j / s with sqrt(s) = 91.1876 GeV
//   binning    reproduction/include/EECBinning.h, 200 bins: 100 logarithmic in
//              theta from 0.002 to pi/2, reflected about pi/2, transformed to z
//
// Statistical uncertainties are computed from the spread of the PER-EVENT EEC
// across events, not from sqrt(sum w^2) over pairs.  Pairs inside one event are
// correlated, so the pair-level estimate would be too small.
//
// Usage:
//   root -l -b -q 'macros/plot_eec_isr_correction.C("/output/dir")'

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

const char* kRealNtup = "/data2/yjlee/ISRsample/real_3M_20260511";
const char* kKkmcNtup = "/data2/yjlee/ISRsample/kkmc_1M_20260921";

constexpr int kHalfAngleBins = 100;
constexpr int kAngleBins = 2 * kHalfAngleBins;
constexpr double kThetaMin = 0.002;
constexpr double kSqrtS = 91.1876;

// Exact copy of EEC::ThetaEdges() / EEC::ZEdges().
const std::vector<double>& thetaEdges()
{
    static const std::vector<double> edges = []() {
        std::vector<double> r(kAngleBins + 1);
        const double middle = std::acos(-1.0) / 2.0;
        for (int i = 0; i <= kHalfAngleBins; ++i) {
            const double low = std::exp(std::log(kThetaMin) +
                (std::log(middle) - std::log(kThetaMin)) * i / kHalfAngleBins);
            r[i] = low;
            r[kAngleBins - i] = 2.0 * middle - low;
        }
        return r;
    }();
    return edges;
}

const std::vector<double>& zEdges()
{
    static const std::vector<double> edges = []() {
        std::vector<double> r(kAngleBins + 1);
        const auto& th = thetaEdges();
        for (int i = 0; i <= kAngleBins; ++i) r[i] = (1.0 - std::cos(th[i])) / 2.0;
        return r;
    }();
    return edges;
}

int findBin(double v, const std::vector<double>& e)
{
    if (!std::isfinite(v) || v < e.front() || v >= e.back()) return -1;
    const auto up = std::upper_bound(e.begin(), e.end(), v);
    const int i = static_cast<int>(up - e.begin()) - 1;
    return (i >= 0 && i + 1 < static_cast<int>(e.size())) ? i : -1;
}

struct EecResult {
    std::vector<double> sumW;    // sum over events of the per-event EEC
    std::vector<double> sumW2;   // sum of squares of the per-event EEC
    Long64_t events = 0;
    bool ok = false;
};

// Accumulate the charged-particle EEC, per event, from an Events tree.
EecResult buildEec(const std::string& path, Long64_t maxEvents)
{
    EecResult r;
    r.sumW.assign(kAngleBins, 0.0);
    r.sumW2.assign(kAngleBins, 0.0);

    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) { std::cerr << "  [skip] cannot open " << path << std::endl; return r; }
    TTree* t = static_cast<TTree*>(f->Get("Events"));
    if (!t) { std::cerr << "  [skip] no Events tree in " << path << std::endl; f->Close(); return r; }

    std::vector<char>* isFinal = nullptr;
    std::vector<float>*px = nullptr, *py = nullptr, *pz = nullptr, *en = nullptr, *ch = nullptr;
    t->SetBranchStatus("*", 0);
    for (const char* b : {"isFinal", "px", "py", "pz", "energy", "charge"})
        t->SetBranchStatus(b, 1);
    t->SetBranchAddress("isFinal", &isFinal);
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("energy", &en);
    t->SetBranchAddress("charge", &ch);

    Long64_t n = t->GetEntries();
    if (maxEvents > 0 && n > maxEvents) n = maxEvents;

    const std::vector<double>& ze = zEdges();
    std::vector<double> evt(kAngleBins, 0.0);
    std::vector<double> qx, qy, qz, qe;

    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        qx.clear(); qy.clear(); qz.clear(); qe.clear();
        for (size_t j = 0; j < isFinal->size(); ++j) {
            if (!(*isFinal)[j]) continue;
            if (std::fabs((*ch)[j]) <= 0) continue;
            const double m = std::sqrt((*px)[j] * (*px)[j] + (*py)[j] * (*py)[j] +
                                       (*pz)[j] * (*pz)[j]);
            if (!(m > 0) || !((*en)[j] > 0)) continue;
            qx.push_back((*px)[j]); qy.push_back((*py)[j]);
            qz.push_back((*pz)[j]); qe.push_back((*en)[j]);
        }
        std::fill(evt.begin(), evt.end(), 0.0);
        const size_t np = qx.size();
        for (size_t a = 0; a < np; ++a) {
            const double ma = std::sqrt(qx[a]*qx[a] + qy[a]*qy[a] + qz[a]*qz[a]);
            for (size_t b = a + 1; b < np; ++b) {
                const double mb = std::sqrt(qx[b]*qx[b] + qy[b]*qy[b] + qz[b]*qz[b]);
                double c = (qx[a]*qx[b] + qy[a]*qy[b] + qz[a]*qz[b]) / (ma * mb);
                c = std::clamp(c, -1.0, 1.0);
                const double z = (1.0 - c) / 2.0;
                const int bin = findBin(z, ze);
                if (bin < 0) continue;
                evt[bin] += qe[a] * qe[b] / (kSqrtS * kSqrtS);
            }
        }
        for (int b = 0; b < kAngleBins; ++b) {
            r.sumW[b] += evt[b];
            r.sumW2[b] += evt[b] * evt[b];
        }
        ++r.events;
        if ((i + 1) % 500000 == 0) std::cout << "    " << i + 1 << "/" << n << std::endl;
    }
    f->Close();
    r.ok = r.events > 0;
    std::cout << "  [eec] " << path << "  events " << r.events << std::endl;
    return r;
}

// Mean per-event EEC in a bin and the uncertainty on that mean.
void meanAndError(const EecResult& r, int b, double& mean, double& err)
{
    mean = err = 0;
    if (!r.ok) return;
    const double n = static_cast<double>(r.events);
    mean = r.sumW[b] / n;
    const double var = r.sumW2[b] / n - mean * mean;
    err = var > 0 ? std::sqrt(var / n) : 0.0;
}

}  // namespace

void plot_eec_isr_correction(const char* outDir =
                                 "/data2/yjlee/ISRsample/kkmc_1M_20260921/results",
                             Long64_t maxEvents = -1)
{
    gStyle->SetOptStat(0);
    gStyle->SetPadTickX(1);
    gStyle->SetPadTickY(1);
    gSystem->mkdir(outDir, kTRUE);

    struct Sample { std::string label; std::string off; std::string on; int color; int marker; };
    std::vector<Sample> samples = {
        {"Pythia 8.315", std::string(kRealNtup) + "/mc_Pythia8315_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Pythia8315_ISR_ON.root", kBlack, 20},
        {"Pythia 8.315 Vincia", std::string(kRealNtup) + "/mc_Pythia8315_Vincia_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Pythia8315_Vincia_ISR_ON.root", kGray + 2, 24},
        {"Sherpa 3.0.3 PDFESherpa", std::string(kRealNtup) + "/mc_Sherpa303_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Sherpa303_ISR_ON.root", kAzure + 2, 21},
        {"Sherpa 3.0.3 YFS", std::string(kRealNtup) + "/mc_Sherpa303_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Sherpa303_ISR_YFS.root", kAzure + 7, 25},
        {"Herwig 7.3.0 QED shower, ISR unchanged",
         std::string(kRealNtup) + "/mc_Herwig730_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Herwig730_QEDshower.root", kOrange + 7, 22},
        {"KKMC 4.30 CEEX", std::string(kKkmcNtup) + "/mc_KKMC424_ISR_OFF.root",
         std::string(kKkmcNtup) + "/mc_KKMC424_ISR_ON.root", kRed + 1, 29},
    };

    const std::vector<double>& ze = zEdges();
    std::ofstream csv(std::string(outDir) + "/eec_isr_correction.csv");
    // The z bin edges crowd towards 1 on the back-to-back side; at the default
    // 6 significant digits several of the highest bins print with z_low equal to
    // z_high, which makes the file useless for forming a density there.
    csv << std::setprecision(12);
    csv << "sample,bin,z_low,z_high,eec_off,eec_off_err,eec_on,eec_on_err,c_isr,c_isr_err\n";

    std::vector<TGraphErrors*> graphs;
    TCanvas* c1 = new TCanvas("c_eec", "", 1000, 800);
    c1->SetTopMargin(0.09);
    c1->SetLeftMargin(0.13);
    c1->SetLogx();
    TH1D* frame = new TH1D("frame_eec", "", 1, ze[1], 1.0);
    frame->SetMinimum(0.90);
    frame->SetMaximum(1.20);
    frame->GetXaxis()->SetTitle("z = (1 - cos#theta)/2");
    frame->GetYaxis()->SetTitle("C_{ISR} = EEC_{ISR OFF} / EEC_{ISR ON}");
    frame->GetXaxis()->SetTitleSize(0.045);
    frame->GetYaxis()->SetTitleSize(0.045);
    frame->GetXaxis()->SetLabelSize(0.040);
    frame->GetYaxis()->SetLabelSize(0.040);
    frame->Draw();
    TLine* unity = new TLine(ze[1], 1.0, 1.0, 1.0);
    unity->SetLineStyle(2);
    unity->SetLineColor(kGray + 1);
    unity->Draw();
    TLegend* leg = new TLegend(0.17, 0.62, 0.66, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextSize(0.029);

    for (const auto& s : samples) {
        std::cout << "[eec] " << s.label << std::endl;
        EecResult off = buildEec(s.off, maxEvents);
        EecResult on = buildEec(s.on, maxEvents);
        if (!off.ok || !on.ok) continue;

        std::vector<double> x, y, ex, ey;
        double intOff = 0, intOn = 0;
        for (int b = 0; b < kAngleBins; ++b) {
            double mo, eo, mn, eN;
            meanAndError(off, b, mo, eo);
            meanAndError(on, b, mn, eN);
            intOff += mo;
            intOn += mn;
            double r = 0, er = 0;
            if (mo > 0 && mn > 0) {
                r = mo / mn;
                er = r * std::sqrt((eo / mo) * (eo / mo) + (eN / mn) * (eN / mn));
                const double zc = 0.5 * (ze[b] + ze[b + 1]);
                x.push_back(zc); y.push_back(r); ex.push_back(0.0); ey.push_back(er);
            }
            csv << "\"" << s.label << "\"," << b << "," << ze[b] << "," << ze[b + 1] << ","
                << mo << "," << eo << "," << mn << "," << eN << "," << r << "," << er << "\n";
        }
        std::cout << "  integral EEC: OFF " << intOff << "  ON " << intOn
                  << "  ratio " << (intOn > 0 ? intOff / intOn : 0) << std::endl;
        csv << "\"" << s.label << "\",integral,,," << intOff << ",," << intOn << ","
            << "," << (intOn > 0 ? intOff / intOn : 0) << ",\n";

        if (x.empty()) continue;
        TGraphErrors* g = new TGraphErrors(static_cast<int>(x.size()), &x[0], &y[0], &ex[0], &ey[0]);
        g->SetLineColor(s.color);
        g->SetMarkerColor(s.color);
        g->SetMarkerStyle(s.marker);
        g->SetMarkerSize(0.9);
        g->SetLineWidth(2);
        g->Draw("P SAME");
        leg->AddEntry(g, s.label.c_str(), "lp");
        graphs.push_back(g);
    }
    csv.close();

    leg->Draw();
    TLatex tx;
    tx.SetNDC();
    tx.SetTextSize(0.032);
    tx.SetTextColor(kGray + 2);
    tx.DrawLatex(0.13, 0.935, "Charged EEC, ISR study, work in progress");
    tx.DrawLatex(0.62, 0.935, "stat. uncertainties only");
    c1->SaveAs(Form("%s/eec_isr_correction.png", outDir));
    c1->SaveAs(Form("%s/eec_isr_correction.pdf", outDir));

    std::cout << "[done] wrote EEC ISR correction to " << outDir << std::endl;
}
