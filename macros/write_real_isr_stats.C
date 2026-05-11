#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <cmath>
#include <string>
#include <vector>

struct StatSample {
    std::string label;
    std::string mode;
    std::string file;
};

std::string statJoinPath(const std::string& dir, const std::string& name)
{
    return dir + "/" + name;
}

constexpr double kVisibleEtaMax = 1.74;

bool statIsNeutrinoPdg(int pdg)
{
    const int apdg = std::abs(pdg);
    return apdg == 12 || apdg == 14 || apdg == 16;
}

double statEtaFromMomentum(double px, double py, double pz)
{
    const double p = std::sqrt(px * px + py * py + pz * pz);
    if (p <= 1e-12) return 1e9;
    const double plus = p + pz;
    const double minus = p - pz;
    if (plus <= 0.0) return -1e9;
    if (minus <= 0.0) return 1e9;
    return 0.5 * std::log(plus / minus);
}

double statAcceptedVisibleEnergy(const std::vector<int>& pdgId,
                                 const std::vector<float>& px,
                                 const std::vector<float>& py,
                                 const std::vector<float>& pz,
                                 const std::vector<float>& energy,
                                 const std::vector<char>& isFinal)
{
    double evis = 0.0;
    for (size_t i = 0; i < pdgId.size(); ++i) {
        if (!isFinal[i]) continue;
        if (statIsNeutrinoPdg(pdgId[i])) continue;
        if (std::abs(statEtaFromMomentum(px[i], py[i], pz[i])) >= kVisibleEtaMax) continue;
        evis += energy[i];
    }
    return evis;
}

void write_real_isr_stats()
{
    std::string realDir = "/data2/yjlee/ISRsample/real_3M_20260511";
    if (const char* envDir = gSystem->Getenv("ISR_REAL_DIR")) realDir = envDir;
    std::string outDir = statJoinPath(realDir, "results");
    if (const char* envOut = gSystem->Getenv("ISR_OUT_DIR")) outDir = envOut;
    gSystem->mkdir(outDir.c_str(), true);

    const std::vector<StatSample> samples = {
        {"Herwig 7.3.0 QED shower", "OFF", statJoinPath(realDir, "mc_Herwig730_ISR_OFF.root")},
        {"Herwig 7.3.0 QED shower", "QEDshower", statJoinPath(realDir, "mc_Herwig730_QEDshower.root")},
        {"PYTHIA 8.315", "OFF", statJoinPath(realDir, "mc_Pythia8315_ISR_OFF.root")},
        {"PYTHIA 8.315", "ISR ON", statJoinPath(realDir, "mc_Pythia8315_ISR_ON.root")},
        {"Sherpa 3.0.3 PDFESherpa", "OFF", statJoinPath(realDir, "mc_Sherpa303_ISR_OFF.root")},
        {"Sherpa 3.0.3 PDFESherpa", "ISR ON (PDFESherpa)", statJoinPath(realDir, "mc_Sherpa303_ISR_ON.root")},
        {"Sherpa 3.0.3 YFS", "OFF", statJoinPath(realDir, "mc_Sherpa303_ISR_OFF.root")},
        {"Sherpa 3.0.3 YFS", "ISR ON (YFS)", statJoinPath(realDir, "mc_Sherpa303_ISR_YFS.root")},
        {"PYTHIA 8.315 (Vincia)", "OFF", statJoinPath(realDir, "mc_Pythia8315_Vincia_ISR_OFF.root")},
        {"PYTHIA 8.315 (Vincia)", "ISR ON", statJoinPath(realDir, "mc_Pythia8315_Vincia_ISR_ON.root")}
    };

    const std::string outName = statJoinPath(outDir, "real_generator_sample_statistics.csv");
    std::ofstream out(outName.c_str());
    out << std::setprecision(10);
    out << "sample,isr,entries,mean_thrust,mean_visibleEnergy,mean_nISRPhotons,mean_totalISRPhotonEnergy,file\n";

    for (const StatSample& sample : samples) {
        TFile* file = TFile::Open(sample.file.c_str(), "READ");
        if (!file || file->IsZombie()) {
            std::cerr << "Cannot open " << sample.file << std::endl;
            continue;
        }
        TTree* tree = static_cast<TTree*>(file->Get("Events"));
        if (!tree) {
            std::cerr << "Missing Events tree in " << sample.file << std::endl;
            file->Close();
            delete file;
            continue;
        }

        TTreeReader reader(tree);
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
            sumVisibleEnergy += statAcceptedVisibleEnergy(*pdgId, *px, *py, *pz, *energy, *isFinal);
            sumNISRPhotons += *nISRPhotons;
            sumISRPhotonEnergy += *totalISRPhotonEnergy;
        }

        const long double denom = entries > 0 ? entries : 1;
        out << sample.label << "," << sample.mode << "," << entries << ","
            << static_cast<double>(sumThrust / denom) << ","
            << static_cast<double>(sumVisibleEnergy / denom) << ","
            << static_cast<double>(sumNISRPhotons / denom) << ","
            << static_cast<double>(sumISRPhotonEnergy / denom) << ","
            << sample.file << "\n";
        file->Close();
        delete file;
    }

    std::cout << "Wrote " << outName << std::endl;
}
