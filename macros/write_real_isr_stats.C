#include <TFile.h>
#include <TSystem.h>
#include <TTree.h>
#include <TTreeReader.h>
#include <TTreeReaderValue.h>

#include <fstream>
#include <iomanip>
#include <iostream>
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
        TTreeReaderValue<double> visibleEnergy(reader, "visibleEnergy");
        TTreeReaderValue<int> nISRPhotons(reader, "nISRPhotons");
        TTreeReaderValue<double> totalISRPhotonEnergy(reader, "totalISRPhotonEnergy");

        Long64_t entries = 0;
        long double sumThrust = 0.0;
        long double sumVisibleEnergy = 0.0;
        long double sumNISRPhotons = 0.0;
        long double sumISRPhotonEnergy = 0.0;
        while (reader.Next()) {
            ++entries;
            sumThrust += *thrust;
            sumVisibleEnergy += *visibleEnergy;
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
