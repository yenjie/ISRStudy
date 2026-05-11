#include <TCanvas.h>
#include <TColor.h>
#include <TFile.h>
#include <TH1.h>
#include <TH2D.h>
#include <TLatex.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>
#include <TString.h>

#include <array>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

struct EndpointSample {
    std::string label;
    std::string tag;
    std::string file;
};

struct EndpointSelection {
    std::string tag;
    std::string label;
    double minT;
};

struct VarDef {
    std::string tag;
    std::string title;
    int bins;
    double xmin;
    double xmax;
};

struct EndpointSummary {
    long long entries = 0;
    long double sumT = 0.0;
    long double sumTIncl = 0.0;
    long double sumTCM = 0.0;
    long double sumMvis = 0.0;
    long double sumSprime = 0.0;
    long double sumBetaZ = 0.0;
    long double sumLeadE = 0.0;
    long double sumISR = 0.0;
    long double sumNISR = 0.0;
};

std::string gEndpointRealDir = "/data2/yjlee/ISRsample/real_3M_20260511";
std::string gEndpointDiagDir = "/data2/yjlee/ISRsample/real_3M_20260511/endpoint_diagnostics";
std::string gEndpointOutDir = "/data2/yjlee/ISRsample/real_3M_20260511/results/endpoint_diagnostics";

std::string epJoin(const std::string& dir, const std::string& name)
{
    return dir + "/" + name;
}

std::string epSanitize(const std::string& name)
{
    std::string out = name;
    for (char& c : out) {
        if (!(std::isalnum(c) || c == '_')) c = '_';
    }
    return out;
}

void addSummary(EndpointSummary& s, double t, double tIncl, double tCM, double mvis, double sprime,
                double betaZ, double leadE, double sumISR, int nISR)
{
    ++s.entries;
    s.sumT += t;
    s.sumTIncl += tIncl;
    s.sumTCM += tCM;
    s.sumMvis += mvis;
    s.sumSprime += sprime;
    s.sumBetaZ += betaZ;
    s.sumLeadE += std::max(0.0, leadE);
    s.sumISR += sumISR;
    s.sumNISR += nISR;
}

double meanOrZero(long double sum, long long n)
{
    return n > 0 ? static_cast<double>(sum / n) : 0.0;
}

void drawEndpointCanvas(const EndpointSample& sample,
                        const EndpointSelection& selection,
                        const std::vector<VarDef>& vars,
                        const std::vector<TH2D*>& hists)
{
    TCanvas* c = new TCanvas(("c_endpoint_" + sample.tag + "_" + selection.tag).c_str(),
                             "endpoint diagnostics", 1650, 1000);
    c->Divide(3, 2, 0.002, 0.002);
    for (size_t i = 0; i < vars.size(); ++i) {
        c->cd(static_cast<int>(i + 1));
        gPad->SetTicks(1, 1);
        gPad->SetRightMargin(0.145);
        gPad->SetLeftMargin(0.115);
        gPad->SetBottomMargin(0.115);
        gPad->SetTopMargin(0.105);
        gPad->SetLogz();
        hists[i]->SetTitle("");
        hists[i]->GetXaxis()->SetTitle(vars[i].title.c_str());
        hists[i]->GetYaxis()->SetTitle("T_{lab}^{all final}, incl. #nu and ISR #gamma");
        hists[i]->GetXaxis()->SetTitleSize(0.046);
        hists[i]->GetYaxis()->SetTitleSize(0.046);
        hists[i]->GetYaxis()->SetTitleOffset(1.05);
        hists[i]->Draw("COLZ");

        TLatex lat;
        lat.SetNDC();
        lat.SetTextFont(42);
        lat.SetTextSize(0.038);
        lat.DrawLatex(0.13, 0.915, sample.label.c_str());
        lat.SetTextSize(0.030);
        lat.DrawLatex(0.13, 0.865, selection.label.c_str());
    }
    const std::string base = epJoin(gEndpointOutDir, "endpoint_diagnostics_" + sample.tag + "_" + selection.tag);
    c->SaveAs((base + ".png").c_str());
    c->SaveAs((base + ".pdf").c_str());
    delete c;
}

void writeSummaryFiles(const std::string& csvName,
                       const std::string& mdName,
                       const std::vector<EndpointSample>& samples,
                       const std::vector<EndpointSelection>& selections,
                       const std::vector<std::vector<EndpointSummary>>& summaries)
{
    std::ofstream csv(csvName.c_str());
    csv << std::setprecision(10);
    csv << "sample,selection,entries,fraction,mean_T_lab_allFinal_including_ISR_photons,mean_T_lab_visible_including_ISR_photons,"
        << "mean_T_allFinalCM_including_ISR_photons,"
        << "mean_Mvis_excluding_ISR_photons,mean_sPrime_vis_over_s,mean_betaZ_vis,"
        << "mean_leading_ISR_photon_energy,mean_sum_ISR_photon_energy,mean_N_ISR_photons\n";

    for (size_t is = 0; is < samples.size(); ++is) {
        const long long total = summaries[is][0].entries > 0 ? summaries[is][0].entries : 1;
        for (size_t ic = 0; ic < selections.size(); ++ic) {
            const EndpointSummary& s = summaries[is][ic];
            csv << samples[is].label << "," << selections[ic].tag << "," << s.entries << ","
                << static_cast<double>(s.entries) / total << ","
                << meanOrZero(s.sumT, s.entries) << ","
                << meanOrZero(s.sumTIncl, s.entries) << ","
                << meanOrZero(s.sumTCM, s.entries) << ","
                << meanOrZero(s.sumMvis, s.entries) << ","
                << meanOrZero(s.sumSprime, s.entries) << ","
                << meanOrZero(s.sumBetaZ, s.entries) << ","
                << meanOrZero(s.sumLeadE, s.entries) << ","
                << meanOrZero(s.sumISR, s.entries) << ","
                << meanOrZero(s.sumNISR, s.entries) << "\n";
        }
    }

    std::ofstream md(mdName.c_str());
    md << "# Endpoint ISR diagnostic summary\n\n";
    md << "Diagnostic thrust input: all stable final-state particles, including neutrinos and tagged ISR photons. "
       << "The plotted T_lab variable is `T_lab_allFinal_including_ISR_photons`.\n\n";
    md << "| sample | selection | entries | fraction | <T_lab all final incl ISR> | <T_lab visible incl ISR> | <T_allFinalCM incl ISR> | <Mvis excl ISR> [GeV] | <sPrime_vis/s> | <betaZ_vis> | <sum E_ISR> [GeV] | <N_ISR> |\n";
    md << "|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|\n";
    md << std::setprecision(4);
    for (size_t is = 0; is < samples.size(); ++is) {
        const long long total = summaries[is][0].entries > 0 ? summaries[is][0].entries : 1;
        for (size_t ic = 0; ic < selections.size(); ++ic) {
            const EndpointSummary& s = summaries[is][ic];
            md << "| " << samples[is].label
               << " | " << selections[ic].tag
               << " | " << s.entries
               << " | " << static_cast<double>(s.entries) / total
               << " | " << meanOrZero(s.sumT, s.entries)
               << " | " << meanOrZero(s.sumTIncl, s.entries)
               << " | " << meanOrZero(s.sumTCM, s.entries)
               << " | " << meanOrZero(s.sumMvis, s.entries)
               << " | " << meanOrZero(s.sumSprime, s.entries)
               << " | " << meanOrZero(s.sumBetaZ, s.entries)
               << " | " << meanOrZero(s.sumISR, s.entries)
               << " | " << meanOrZero(s.sumNISR, s.entries)
               << " |\n";
        }
    }
}

void plot_endpoint_diagnostics()
{
    if (const char* envReal = gSystem->Getenv("ISR_REAL_DIR")) {
        gEndpointRealDir = envReal;
        gEndpointDiagDir = epJoin(gEndpointRealDir, "endpoint_diagnostics");
        gEndpointOutDir = epJoin(epJoin(gEndpointRealDir, "results"), "endpoint_diagnostics");
    }
    if (const char* envDiag = gSystem->Getenv("ISR_ENDPOINT_DIAG_DIR")) {
        gEndpointDiagDir = envDiag;
    }
    if (const char* envOut = gSystem->Getenv("ISR_ENDPOINT_OUT_DIR")) {
        gEndpointOutDir = envOut;
    }
    gSystem->mkdir(gEndpointOutDir.c_str(), true);
    gStyle->SetOptStat(0);
    gStyle->SetPalette(kBird);
    TH1::AddDirectory(kFALSE);

    const std::vector<EndpointSample> samples = {
        {"Herwig 7.3.0 QED shower", "Herwig730_QEDshower", epJoin(gEndpointDiagDir, "endpoint_diagnostics_Herwig730_QEDshower.root")},
        {"PYTHIA 8.315", "Pythia8315", epJoin(gEndpointDiagDir, "endpoint_diagnostics_Pythia8315.root")},
        {"Sherpa 3.0.3 PDFESherpa", "Sherpa303", epJoin(gEndpointDiagDir, "endpoint_diagnostics_Sherpa303.root")},
        {"Sherpa 3.0.3 YFS", "Sherpa303_YFS", epJoin(gEndpointDiagDir, "endpoint_diagnostics_Sherpa303_YFS.root")},
        {"PYTHIA 8.315 (Vincia)", "Pythia8315_Vincia", epJoin(gEndpointDiagDir, "endpoint_diagnostics_Pythia8315_Vincia.root")}
    };

    const double lastBinLow = 0.99;
    const std::vector<EndpointSelection> selections = {
        {"all", "all ISR ON events", -1.0},
        {"Tgt0p97", "T_{lab} > 0.97", 0.97},
        {"Tgt0p98", "T_{lab} > 0.98", 0.98},
        {"lastBin", Form("last plotted thrust bin: T_{lab} > %.4f", lastBinLow), lastBinLow}
    };

    const std::vector<VarDef> vars = {
        {"Mvis", "M_{vis} excluding ISR photons [GeV]", 80, 0.0, 100.0},
        {"sPrime", "s'_{vis}/s = M_{vis}^{2}/s", 80, 0.0, 1.05},
        {"betaZ", "#beta_{z,vis}=|P_{z,vis}|/E_{vis}", 80, 0.0, 1.0},
        {"leadISR", "leading ISR photon E [GeV]", 80, 0.0, 80.0},
        {"sumISR", "#Sigma E_{ISR #gamma} [GeV]", 80, 0.0, 100.0},
        {"nISR", "N_{ISR #gamma}", 60, -0.5, 59.5}
    };

    std::vector<std::vector<EndpointSummary>> summaries(samples.size(), std::vector<EndpointSummary>(selections.size()));

    for (size_t is = 0; is < samples.size(); ++is) {
        TFile f(samples[is].file.c_str(), "READ");
        if (f.IsZombie()) {
            std::cerr << "Cannot open " << samples[is].file << std::endl;
            continue;
        }
        TTree* tree = static_cast<TTree*>(f.Get("EndpointDiagnostics"));
        if (!tree) {
            std::cerr << "Missing EndpointDiagnostics tree in " << samples[is].file << std::endl;
            continue;
        }

        std::vector<std::vector<TH2D*>> hists(selections.size());
        for (size_t ic = 0; ic < selections.size(); ++ic) {
            const double ymin = selections[ic].minT > 0.0 ? selections[ic].minT : 0.55;
            hists[ic].reserve(vars.size());
            for (const VarDef& v : vars) {
                const std::string hname = "h_" + samples[is].tag + "_" + selections[ic].tag + "_" + v.tag;
                TH2D* h = new TH2D(hname.c_str(), "", v.bins, v.xmin, v.xmax, 72, ymin, 1.005);
                hists[ic].push_back(h);
            }
        }

        TTreeReader reader(tree);
        TTreeReaderValue<double> tLab(reader, "T_lab_allFinal_including_ISR_photons");
        TTreeReaderValue<double> tLabIncl(reader, "T_lab_including_ISR_photons");
        TTreeReaderValue<double> tCM(reader, "T_allFinalCM_including_ISR_photons");
        TTreeReaderValue<double> mvis(reader, "Mvis_excluding_ISR_photons");
        TTreeReaderValue<double> sprime(reader, "sPrime_vis_over_s");
        TTreeReaderValue<double> betaZ(reader, "betaZ_vis");
        TTreeReaderValue<double> leadE(reader, "leading_ISR_photon_energy");
        TTreeReaderValue<double> sumISR(reader, "sum_ISR_photon_energy");
        TTreeReaderValue<int> nISR(reader, "N_ISR_photons");

        while (reader.Next()) {
            const double varsValue[6] = {
                *mvis,
                *sprime,
                *betaZ,
                std::max(0.0, *leadE),
                *sumISR,
                static_cast<double>(*nISR)
            };
            for (size_t ic = 0; ic < selections.size(); ++ic) {
                if (selections[ic].minT > 0.0 && *tLab <= selections[ic].minT) continue;
                addSummary(summaries[is][ic], *tLab, *tLabIncl, *tCM, *mvis, *sprime, *betaZ, *leadE, *sumISR, *nISR);
                for (size_t iv = 0; iv < vars.size(); ++iv) hists[ic][iv]->Fill(varsValue[iv], *tLab);
            }
        }

        for (size_t ic = 0; ic < selections.size(); ++ic) {
            drawEndpointCanvas(samples[is], selections[ic], vars, hists[ic]);
            for (TH2D* h : hists[ic]) delete h;
        }
    }

    writeSummaryFiles(epJoin(gEndpointOutDir, "endpoint_diagnostics_summary.csv"),
                      epJoin(gEndpointOutDir, "endpoint_diagnostics_summary.md"),
                      samples, selections, summaries);
}
