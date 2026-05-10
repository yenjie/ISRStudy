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
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TVector3.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <vector>

Long64_t gMaxEntriesPerFile = 100000;

struct Sample {
    std::string label;
    std::string tag;
    std::string onFile;
    std::string offFile;
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

struct HistResult {
    TH1D* corrected = nullptr;
    TH1D* stored = nullptr;
    Long64_t entries = 0;
    double sumW = 0.0;
    double meanCorrected = 0.0;
    double meanStored = 0.0;
    double meanDelta = 0.0;
    double meanAbsDelta = 0.0;
    double maxAbsDelta = 0.0;
    double meanNSelected = 0.0;
};

bool isNeutrino(int pdg)
{
    const int a = std::abs(pdg);
    return a == 12 || a == 14 || a == 16;
}

TVector3 computeThrustAxis(const std::vector<TVector3>& momenta, double* thrustValue)
{
    double pSum = 0.0;
    std::vector<TVector3> pvec;
    pvec.reserve(momenta.size());
    for (const TVector3& p : momenta) {
        const double mag = p.Mag();
        if (mag <= 1e-9) continue;
        pSum += mag;
        pvec.push_back(p);
    }
    if (pSum <= 1e-9 || pvec.empty()) {
        if (thrustValue) *thrustValue = 0.0;
        return TVector3(0, 0, 1);
    }

    TVector3 best(0, 0, 0);
    double bestMag2 = -1.0;
    const size_t n = pvec.size();
    if (n <= 3) {
        const unsigned int nMasks = 1u << n;
        for (unsigned int mask = 0; mask < nMasks; ++mask) {
            TVector3 sum(0, 0, 0);
            for (size_t i = 0; i < n; ++i) sum += ((mask >> i) & 1u) ? pvec[i] : -pvec[i];
            const double mag2 = sum.Mag2();
            if (mag2 > bestMag2) {
                bestMag2 = mag2;
                best = sum;
            }
        }
    } else {
        for (size_t i = 1; i < n; ++i) {
            for (size_t j = 0; j < i; ++j) {
                const TVector3 cross = pvec[i].Cross(pvec[j]);
                TVector3 ptot(0, 0, 0);
                for (size_t k = 0; k < n; ++k) {
                    if (k == i || k == j) continue;
                    ptot += (pvec[k].Dot(cross) > 0.0) ? pvec[k] : -pvec[k];
                }
                const TVector3 candidates[4] = {
                    ptot - pvec[j] - pvec[i],
                    ptot - pvec[j] + pvec[i],
                    ptot + pvec[j] - pvec[i],
                    ptot + pvec[j] + pvec[i]
                };
                for (int c = 0; c < 4; ++c) {
                    const double mag2 = candidates[c].Mag2();
                    if (mag2 > bestMag2) {
                        bestMag2 = mag2;
                        best = candidates[c];
                    }
                }
            }
        }
    }
    const double t = std::min(1.0, std::max(0.0, best.Mag() / pSum));
    if (thrustValue) *thrustValue = t;
    return best.Mag() > 1e-12 ? best.Unit() : pvec[0].Unit();
}

std::vector<double> makeThrustBins()
{
    return {
        0.570, 0.585, 0.600, 0.615, 0.630, 0.645, 0.660,
        0.670, 0.680, 0.690, 0.700, 0.710, 0.720, 0.730,
        0.740, 0.750, 0.760, 0.770, 0.780, 0.790, 0.800,
        0.810, 0.820, 0.830, 0.840, 0.850, 0.860, 0.870,
        0.880, 0.890, 0.900, 0.910, 0.920, 0.930, 0.940,
        0.950, 0.960, 0.970, 0.975, 0.980, 0.985, 0.990,
        0.995, 1.000
    };
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

HistResult fillThrustHist(const std::string& filename, const std::vector<double>& edges,
                          const std::string& name)
{
    HistResult r;
    r.corrected = new TH1D((name + "_corrected").c_str(), "", edges.size() - 1, edges.data());
    r.stored = new TH1D((name + "_stored").c_str(), "", edges.size() - 1, edges.data());
    r.corrected->Sumw2();
    r.stored->Sumw2();
    r.corrected->SetDirectory(0);
    r.stored->SetDirectory(0);

    TFile file(filename.c_str(), "READ");
    TTree* tree = static_cast<TTree*>(file.Get("Events"));
    if (!tree) {
        std::cerr << "Missing Events tree in " << filename << std::endl;
        return r;
    }

    double weight = 1.0;
    double storedThrust = 0.0;
    std::vector<int>* pdgId = nullptr;
    std::vector<float>* px = nullptr;
    std::vector<float>* py = nullptr;
    std::vector<float>* pz = nullptr;
    std::vector<char>* isFinal = nullptr;
    std::vector<char>* isISRPhoton = nullptr;

    tree->SetBranchStatus("*", 0);
    tree->SetBranchStatus("weight", 1);
    tree->SetBranchStatus("thrust", 1);
    tree->SetBranchStatus("pdgId", 1);
    tree->SetBranchStatus("px", 1);
    tree->SetBranchStatus("py", 1);
    tree->SetBranchStatus("pz", 1);
    tree->SetBranchStatus("isFinal", 1);
    tree->SetBranchStatus("isISRPhoton", 1);
    tree->SetBranchAddress("weight", &weight);
    tree->SetBranchAddress("thrust", &storedThrust);
    tree->SetBranchAddress("pdgId", &pdgId);
    tree->SetBranchAddress("px", &px);
    tree->SetBranchAddress("py", &py);
    tree->SetBranchAddress("pz", &pz);
    tree->SetBranchAddress("isFinal", &isFinal);
    tree->SetBranchAddress("isISRPhoton", &isISRPhoton);

    const Long64_t nEntries = gMaxEntriesPerFile > 0 ? std::min(gMaxEntriesPerFile, tree->GetEntries()) : tree->GetEntries();
    double sumWStored = 0.0, sumWCorr = 0.0, sumDelta = 0.0, sumAbsDelta = 0.0, sumN = 0.0;
    for (Long64_t i = 0; i < nEntries; ++i) {
        tree->GetEntry(i);
        std::vector<TVector3> selected;
        selected.reserve(px->size());
        for (size_t j = 0; j < pdgId->size(); ++j) {
            if (!(*isFinal)[j]) continue;
            if ((*isISRPhoton)[j]) continue;
            if (isNeutrino((*pdgId)[j])) continue;
            TVector3 p((*px)[j], (*py)[j], (*pz)[j]);
            if (p.Mag() <= 1e-9) continue;
            selected.push_back(p);
        }
        double corrected = 0.0;
        computeThrustAxis(selected, &corrected);
        r.corrected->Fill(corrected, weight);
        r.stored->Fill(storedThrust, weight);
        r.entries++;
        r.sumW += weight;
        sumWCorr += weight * corrected;
        sumWStored += weight * storedThrust;
        const double delta = storedThrust - corrected;
        sumDelta += weight * delta;
        sumAbsDelta += weight * std::abs(delta);
        r.maxAbsDelta = std::max(r.maxAbsDelta, std::abs(delta));
        sumN += selected.size();
    }
    file.Close();
    if (r.sumW > 0.0) {
        r.meanCorrected = sumWCorr / r.sumW;
        r.meanStored = sumWStored / r.sumW;
        r.meanDelta = sumDelta / r.sumW;
        r.meanAbsDelta = sumAbsDelta / r.sumW;
    }
    if (r.entries > 0) r.meanNSelected = sumN / r.entries;
    return r;
}

void normalizeDensity(TH1D* h)
{
    double integral = 0.0;
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) integral += h->GetBinContent(ib);
    if (integral <= 0.0) return;
    for (int ib = 1; ib <= h->GetNbinsX(); ++ib) {
        const double width = h->GetBinWidth(ib);
        h->SetBinContent(ib, h->GetBinContent(ib) / (integral * width));
        h->SetBinError(ib, h->GetBinError(ib) / (integral * width));
    }
}

TGraphErrors* makeRatioGraph(const TH1D* hOff, const TH1D* hOn, int color, int marker, double xOffset)
{
    std::vector<double> x, y, ex, ey;
    for (int ib = 1; ib <= hOff->GetNbinsX(); ++ib) {
        const double off = hOff->GetBinContent(ib);
        const double on = hOn->GetBinContent(ib);
        if (off <= 0.0 || on <= 0.0) continue;
        const double ratio = off / on;
        const double eoff = hOff->GetBinError(ib);
        const double eon = hOn->GetBinError(ib);
        x.push_back(hOff->GetBinCenter(ib) + xOffset);
        y.push_back(ratio);
        ex.push_back(0.0);
        ey.push_back(ratio * std::sqrt((eoff/off)*(eoff/off) + (eon/on)*(eon/on)));
    }
    TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
    g->SetMarkerColor(color);
    g->SetLineColor(color);
    g->SetMarkerStyle(marker);
    g->SetMarkerSize(0.72);
    g->SetLineWidth(2);
    return g;
}

void drawCorrectedISRPlot(const std::string& outDir, const std::vector<Sample>& samples,
                          const std::vector<TGraphErrors*>& graphs)
{
    TCanvas* c = new TCanvas("c_corrected_isr", "Corrected all-visible ISR correction", 1400, 850);
    c->SetFillColor(kWhite);
    TPad* top = new TPad("top", "top", 0.0, 0.30, 1.0, 1.0);
    TPad* bot = new TPad("bot", "bot", 0.0, 0.0, 1.0, 0.30);
    top->SetLeftMargin(0.105); top->SetRightMargin(0.030); top->SetBottomMargin(0.020); top->SetTopMargin(0.055); top->SetTicks(1, 1);
    bot->SetLeftMargin(0.105); bot->SetRightMargin(0.030); bot->SetBottomMargin(0.290); bot->SetTopMargin(0.020); bot->SetTicks(1, 1);
    top->Draw(); bot->Draw();
    top->cd();
    TH1D* frame = new TH1D("frame_corr", "", 100, 0.57, 1.0);
    frame->SetMinimum(0.70); frame->SetMaximum(1.35);
    frame->GetYaxis()->SetTitle("C_{ISR} = ISR OFF / ISR ON");
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitleSize(0.050);
    frame->GetYaxis()->SetLabelSize(0.042);
    frame->Draw("AXIS");
    TLine* one = new TLine(0.57, 1.0, 1.0, 1.0); one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();
    TLegend* leg = new TLegend(0.42, 0.08, 0.88, 0.25);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.034); leg->SetNColumns(2);
    for (size_t i = 0; i < graphs.size(); ++i) {
        graphs[i]->Draw("LP SAME");
        leg->AddEntry(graphs[i], samples[i].label.c_str(), "lp");
    }
    leg->Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.040);
    lat.DrawLatex(0.145, 0.875, "All visible final-state thrust, recomputed");
    lat.SetTextSize(0.030);
    lat.DrawLatex(0.145, 0.825, "charged + neutral visible particles; neutrinos and ISR photons excluded");
    lat.DrawLatex(0.145, 0.775, Form("diagnostic sample: %lld entries/file", gMaxEntriesPerFile));
    frame->Draw("AXIS SAME");

    bot->cd();
    TH1D* rframe = new TH1D("rframe_corr", "", 100, 0.57, 1.0);
    rframe->SetMinimum(0.75); rframe->SetMaximum(1.25);
    rframe->GetXaxis()->SetTitle("Thrust T");
    rframe->GetYaxis()->SetTitle("MC / P8");
    rframe->GetXaxis()->SetTitleSize(0.105); rframe->GetXaxis()->SetLabelSize(0.085);
    rframe->GetYaxis()->SetTitleSize(0.085); rframe->GetYaxis()->SetLabelSize(0.075); rframe->GetYaxis()->SetTitleOffset(0.55);
    rframe->Draw("AXIS");
    TLine* rone = new TLine(0.57, 1.0, 1.0, 1.0); rone->SetLineStyle(2); rone->SetLineColor(kGray+2); rone->Draw();
    const int ref = 1;
    for (size_t i = 0; i < graphs.size(); ++i) {
        std::vector<double> x, y, ex, ey;
        const int n = std::min(graphs[i]->GetN(), graphs[ref]->GetN());
        for (int j = 0; j < n; ++j) {
            double xi, yi, xr, yr;
            graphs[i]->GetPoint(j, xi, yi); graphs[ref]->GetPoint(j, xr, yr);
            if (yi <= 0.0 || yr <= 0.0) continue;
            x.push_back(xi); y.push_back(i == ref ? 1.0 : yi / yr); ex.push_back(0.0); ey.push_back(0.0);
        }
        TGraphErrors* gr = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
        gr->SetMarkerColor(samples[i].color); gr->SetLineColor(samples[i].color); gr->SetMarkerStyle(samples[i].marker); gr->SetMarkerSize(0.60); gr->SetLineWidth(2);
        gr->Draw("LP SAME");
    }
    rframe->Draw("AXIS SAME");
    c->SaveAs((outDir + "/isr_correction_all_visible_recomputed.png").c_str());
    c->SaveAs((outDir + "/isr_correction_all_visible_recomputed.pdf").c_str());
}

void drawAlephComparison(const std::string& outDir, const std::vector<Sample>& samples,
                         const std::vector<HistResult>& onResults,
                         const std::vector<AlephPoint>& aleph)
{
    std::vector<TH1D*> hists;
    for (size_t i = 0; i < onResults.size(); ++i) {
        TH1D* h = static_cast<TH1D*>(onResults[i].corrected->Clone(Form("h_aleph_%zu", i)));
        h->SetDirectory(0);
        normalizeDensity(h);
        h->SetLineColor(samples[i].color); h->SetMarkerColor(samples[i].color); h->SetMarkerStyle(samples[i].marker); h->SetLineWidth(2); h->SetMarkerSize(0.45);
        hists.push_back(h);
    }
    std::vector<double> ax, ay, aex, aey;
    for (const AlephPoint& p : aleph) {
        ax.push_back(p.center); ay.push_back(p.value); aex.push_back(0.5*(p.high-p.low)); aey.push_back(p.err);
    }
    TGraphErrors* gAleph = new TGraphErrors(ax.size(), ax.data(), ay.data(), aex.data(), aey.data());
    gAleph->SetMarkerStyle(20); gAleph->SetMarkerSize(0.75); gAleph->SetMarkerColor(kBlack); gAleph->SetLineColor(kBlack);

    TCanvas* c = new TCanvas("c_aleph_corr", "ALEPH comparison corrected", 1100, 950);
    c->SetFillColor(kWhite);
    TPad* top = new TPad("topA", "topA", 0.0, 0.30, 1.0, 1.0);
    TPad* bot = new TPad("botA", "botA", 0.0, 0.0, 1.0, 0.30);
    top->SetLeftMargin(0.105); top->SetRightMargin(0.030); top->SetBottomMargin(0.020); top->SetTopMargin(0.055); top->SetTicks(1, 1); top->SetLogy();
    bot->SetLeftMargin(0.105); bot->SetRightMargin(0.030); bot->SetBottomMargin(0.290); bot->SetTopMargin(0.020); bot->SetTicks(1, 1); bot->SetLogy();
    top->Draw(); bot->Draw();
    top->cd();
    TH1D* frame = new TH1D("frameA", "", 100, 0.58, 1.0);
    frame->SetMinimum(5e-4); frame->SetMaximum(40);
    frame->GetYaxis()->SetTitle("(1/N) dN/dT");
    frame->GetXaxis()->SetLabelSize(0);
    frame->GetYaxis()->SetTitleSize(0.050); frame->GetYaxis()->SetLabelSize(0.042);
    frame->Draw("AXIS");
    gAleph->Draw("P SAME");
    for (TH1D* h : hists) h->Draw("HIST SAME");
    for (TH1D* h : hists) h->Draw("E0 SAME");
    gAleph->Draw("P SAME");
    TLegend* leg = new TLegend(0.14, 0.62, 0.58, 0.90);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.031);
    leg->AddEntry(gAleph, "ALEPH published data", "pe");
    for (size_t i = 0; i < hists.size(); ++i) leg->AddEntry(hists[i], (samples[i].label + " ISR ON").c_str(), "lp");
    leg->Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.034);
    lat.DrawLatex(0.60, 0.86, "All-visible thrust recomputed");
    lat.DrawLatex(0.60, 0.81, Form("%lld entries/file", gMaxEntriesPerFile));
    frame->Draw("AXIS SAME");

    bot->cd();
    TH1D* rframe = new TH1D("rframeA", "", 100, 0.58, 1.0);
    rframe->SetMinimum(0.05); rframe->SetMaximum(100);
    rframe->GetXaxis()->SetTitle("Thrust T");
    rframe->GetYaxis()->SetTitle("MC / ALEPH");
    rframe->GetXaxis()->SetTitleSize(0.105); rframe->GetXaxis()->SetLabelSize(0.085);
    rframe->GetYaxis()->SetTitleSize(0.085); rframe->GetYaxis()->SetLabelSize(0.075); rframe->GetYaxis()->SetTitleOffset(0.58);
    rframe->Draw("AXIS");
    TLine* one = new TLine(0.58, 1.0, 1.0, 1.0); one->SetLineStyle(2); one->SetLineColor(kGray+2); one->Draw();
    for (size_t ih = 0; ih < hists.size(); ++ih) {
        std::vector<double> x, y, ex, ey;
        for (int ib = 1; ib <= hists[ih]->GetNbinsX() && ib <= static_cast<int>(aleph.size()); ++ib) {
            if (aleph[ib-1].value <= 0.0) continue;
            x.push_back(aleph[ib-1].center);
            y.push_back(hists[ih]->GetBinContent(ib) / aleph[ib-1].value);
            ex.push_back(0.0); ey.push_back(0.0);
        }
        TGraphErrors* g = new TGraphErrors(x.size(), x.data(), y.data(), ex.data(), ey.data());
        g->SetMarkerColor(samples[ih].color); g->SetLineColor(samples[ih].color); g->SetMarkerStyle(samples[ih].marker); g->SetMarkerSize(0.55); g->SetLineWidth(2);
        g->Draw("LP SAME");
    }
    rframe->Draw("AXIS SAME");
    c->SaveAs((outDir + "/aleph_isr_on_thrust_comparison_all_visible_recomputed.png").c_str());
    c->SaveAs((outDir + "/aleph_isr_on_thrust_comparison_all_visible_recomputed.pdf").c_str());
}

void writeAlephMetrics(const std::string& outDir, const std::vector<Sample>& samples,
                       const std::vector<HistResult>& onResults,
                       const std::vector<AlephPoint>& aleph)
{
    std::ofstream metrics((outDir + "/aleph_isr_on_thrust_comparison_all_visible_recomputed_metrics.csv").c_str());
    metrics << std::setprecision(8);
    metrics << "sample,entries,chi2_full,ndf_full,chi2ndf_full,mean_abs_frac_full,"
            << "chi2_central,ndf_central,chi2ndf_central,mean_abs_frac_central\n";

    std::ofstream bins((outDir + "/aleph_isr_on_thrust_comparison_all_visible_recomputed_bins.csv").c_str());
    bins << std::setprecision(8);
    bins << "sample,T_low,T_center,T_high,ALEPH,ALEPH_total_unc,MC_density,MC_over_ALEPH\n";

    for (size_t is = 0; is < samples.size(); ++is) {
        TH1D* h = static_cast<TH1D*>(onResults[is].corrected->Clone(Form("h_metric_%zu", is)));
        h->SetDirectory(0);
        normalizeDensity(h);

        double chi2Full = 0.0, chi2Central = 0.0;
        double meanAbsFull = 0.0, meanAbsCentral = 0.0;
        int nFull = 0, nCentral = 0;
        for (int ib = 1; ib <= h->GetNbinsX() && ib <= static_cast<int>(aleph.size()); ++ib) {
            const AlephPoint& p = aleph[ib - 1];
            if (p.value <= 0.0) continue;
            const double mc = h->GetBinContent(ib);
            const double frac = (mc - p.value) / p.value;
            if (p.err > 0.0) chi2Full += (mc - p.value) * (mc - p.value) / (p.err * p.err);
            meanAbsFull += std::abs(frac);
            nFull++;
            if (p.center >= 0.65 && p.center < 0.95) {
                if (p.err > 0.0) chi2Central += (mc - p.value) * (mc - p.value) / (p.err * p.err);
                meanAbsCentral += std::abs(frac);
                nCentral++;
            }
            bins << samples[is].label << "," << p.low << "," << p.center << "," << p.high << ","
                 << p.value << "," << p.err << "," << mc << "," << mc / p.value << "\n";
        }
        if (nFull > 0) meanAbsFull /= nFull;
        if (nCentral > 0) meanAbsCentral /= nCentral;
        metrics << samples[is].label << "," << onResults[is].entries << ","
                << chi2Full << "," << nFull << "," << (nFull > 0 ? chi2Full / nFull : 0.0) << ","
                << meanAbsFull << "," << chi2Central << "," << nCentral << ","
                << (nCentral > 0 ? chi2Central / nCentral : 0.0) << "," << meanAbsCentral << "\n";
        delete h;
    }
}

void drawStoredVsCorrectedDiagnostic(const std::string& outDir, const Sample& sample,
                                     const HistResult& result)
{
    TH1D* stored = static_cast<TH1D*>(result.stored->Clone("stored_diag"));
    TH1D* corr = static_cast<TH1D*>(result.corrected->Clone("corr_diag"));
    stored->SetDirectory(0); corr->SetDirectory(0);
    normalizeDensity(stored);
    normalizeDensity(corr);
    stored->SetLineColor(kGray+2); stored->SetLineStyle(2); stored->SetLineWidth(2);
    corr->SetLineColor(sample.color); corr->SetMarkerColor(sample.color); corr->SetMarkerStyle(20); corr->SetLineWidth(2);
    TCanvas* c = new TCanvas("c_stored_vs_corrected", "Stored vs corrected thrust", 1000, 800);
    c->SetFillColor(kWhite);
    c->SetLogy();
    c->SetTicks(1, 1);
    TH1D* frame = new TH1D("frameS", "", 100, 0.57, 1.0);
    frame->SetMinimum(5e-4); frame->SetMaximum(60);
    frame->GetXaxis()->SetTitle("Thrust T");
    frame->GetYaxis()->SetTitle("(1/N) dN/dT");
    frame->Draw("AXIS");
    stored->Draw("HIST SAME");
    corr->Draw("HIST SAME");
    TLegend* leg = new TLegend(0.18, 0.73, 0.62, 0.88);
    leg->SetBorderSize(0); leg->SetFillStyle(0); leg->SetTextFont(42); leg->SetTextSize(0.034);
    leg->AddEntry(stored, "stored 20M thrust branch", "l");
    leg->AddEntry(corr, "recomputed all-visible thrust", "l");
    leg->Draw();
    TLatex lat; lat.SetNDC(); lat.SetTextFont(42); lat.SetTextSize(0.034);
    lat.DrawLatex(0.18, 0.68, (sample.label + " ISR ON diagnostic").c_str());
    lat.DrawLatex(0.18, 0.63, Form("mean |stored - recomputed| = %.4f", result.meanAbsDelta));
    c->SaveAs((outDir + "/thrust_stored_vs_all_visible_recomputed_diagnostic.png").c_str());
    c->SaveAs((outDir + "/thrust_stored_vs_all_visible_recomputed_diagnostic.pdf").c_str());
}

void update_thrust_all_visible_results()
{
    gStyle->SetOptStat(0);
    TH1::AddDirectory(kFALSE);
    const std::string outDir = "/data2/yjlee/ISRsample";
    const std::vector<double> correctionBins = makeThrustBins();
    const std::vector<AlephPoint> aleph = readAleph("/data/yjlee/ALEPH_Agentic_Event_Shape_Analysis/ALEPH/HEPData-ins636645-v1-Table_54.csv");
    const std::vector<double> alephBinEdges = alephEdges(aleph);

    std::vector<Sample> samples = {
        {"Herwig 7.1.5", "Herwig715", outDir + "/mc_Herwig715_ISR_ON.root", outDir + "/mc_Herwig715_ISR_OFF.root", TColor::GetColor("#1f4e9e"), 20, -0.0027},
        {"PYTHIA 8.317", "Pythia8317", outDir + "/mc_Pythia8317_ISR_ON.root", outDir + "/mc_Pythia8317_ISR_OFF.root", TColor::GetColor("#d62728"), 20, -0.0009},
        {"Sherpa 2.2.6", "Sherpa226", outDir + "/mc_Sherpa226_ISR_ON.root", outDir + "/mc_Sherpa226_ISR_OFF.root", TColor::GetColor("#7b3294"), 23, 0.0009},
        {"PYTHIA 8.317 (Vincia)", "Pythia8317_Vincia", outDir + "/mc_Pythia8317_Vincia_ISR_ON.root", outDir + "/mc_Pythia8317_Vincia_ISR_OFF.root", TColor::GetColor("#18a7b5"), 22, 0.0027}
    };

    std::vector<HistResult> corrOn, corrOff, alephOn;
    std::vector<TGraphErrors*> graphs;
    std::ofstream metrics((outDir + "/thrust_all_visible_recomputed_metrics.csv").c_str());
    metrics << std::setprecision(8);
    metrics << "sample,category,entries,sumw,mean_corrected,mean_stored,mean_stored_minus_corrected,mean_abs_delta,max_abs_delta,mean_n_selected\n";

    for (const Sample& s : samples) {
        std::cout << "Recomputing all-visible thrust for " << s.label << " ISR ON" << std::endl;
        HistResult on = fillThrustHist(s.onFile, correctionBins, s.tag + "_on_corrbins");
        std::cout << "Recomputing all-visible thrust for " << s.label << " ISR OFF" << std::endl;
        HistResult off = fillThrustHist(s.offFile, correctionBins, s.tag + "_off_corrbins");
        HistResult onAleph = fillThrustHist(s.onFile, alephBinEdges, s.tag + "_on_alephbins");

        metrics << s.label << ",ISR_ON," << on.entries << "," << on.sumW << "," << on.meanCorrected << "," << on.meanStored << ","
                << on.meanDelta << "," << on.meanAbsDelta << "," << on.maxAbsDelta << "," << on.meanNSelected << "\n";
        metrics << s.label << ",ISR_OFF," << off.entries << "," << off.sumW << "," << off.meanCorrected << "," << off.meanStored << ","
                << off.meanDelta << "," << off.meanAbsDelta << "," << off.maxAbsDelta << "," << off.meanNSelected << "\n";

        graphs.push_back(makeRatioGraph(off.corrected, on.corrected, s.color, s.marker, s.xOffset));
        corrOn.push_back(on);
        corrOff.push_back(off);
        alephOn.push_back(onAleph);
    }
    metrics << "\nNotes\n";
    metrics << "nominal_selection,all visible final-state particles: charged plus neutral; neutrinos and flagged ISR photons excluded\n";
    metrics << "max_entries_per_file," << gMaxEntriesPerFile << "\n";
    metrics << "reason,existing 20M files are retained as particle-level sources; stale thrust-result artifacts are replaced by recomputed all-visible outputs\n";
    metrics.close();

    drawCorrectedISRPlot(outDir, samples, graphs);
    drawAlephComparison(outDir, samples, alephOn, aleph);
    writeAlephMetrics(outDir, samples, alephOn, aleph);
    drawStoredVsCorrectedDiagnostic(outDir, samples[1], corrOn[1]);

    std::ofstream binout((outDir + "/thrust_all_visible_recomputed_histograms.csv").c_str());
    binout << std::setprecision(8);
    binout << "sample,category,T_low,T_center,T_high,corrected_count,stored_count\n";
    for (size_t is = 0; is < samples.size(); ++is) {
        for (int cat = 0; cat < 2; ++cat) {
            const HistResult& r = cat == 0 ? corrOn[is] : corrOff[is];
            for (int ib = 1; ib <= r.corrected->GetNbinsX(); ++ib) {
                binout << samples[is].label << "," << (cat == 0 ? "ISR_ON" : "ISR_OFF") << ","
                       << r.corrected->GetBinLowEdge(ib) << "," << r.corrected->GetBinCenter(ib) << ","
                       << r.corrected->GetBinLowEdge(ib + 1) << ","
                       << r.corrected->GetBinContent(ib) << "," << r.stored->GetBinContent(ib) << "\n";
            }
        }
    }
    binout.close();

    std::cout << "Saved corrected all-visible thrust outputs in " << outDir << std::endl;
}
