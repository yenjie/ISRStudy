// Independent verification of the charged-EEC ISR correction.
//
// The per-event EEC summed over ALL pairs has a closed form,
//
//   sum_{i<j} E_i E_j / s  =  ( (sum_i E_i)^2 - sum_i E_i^2 ) / (2 s),
//
// which needs only a linear pass over the particles.  Comparing it with the
// number the binned O(N^2) pair loop produces tests the pair loop, the z
// binning and the bin lookup at once, using no shared code beyond reading the
// branches.  It also gives the charged energy directly, so the size of the
// correction can be predicted rather than only described:
//
//   EEC ~ (sum E_ch)^2, so losing a fraction d of the charged energy gives
//   C_ISR = EEC_OFF / EEC_ON ~ 1 / (1-d)^2 ~ 1 + 2d.
//
// Also splits each sample in half and compares the two halves, to check the
// quoted statistical uncertainty against an independent estimate.
//
// Usage:
//   root -l -b -q 'macros/verify_eec_isr_correction.C("/out/dir",-1,<sample>)'

#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

namespace {

const char* kRealNtup = "/data2/yjlee/ISRsample/real_3M_20260511";
const char* kKkmcNtup = "/data2/yjlee/ISRsample/kkmc_1M_20260921";
constexpr double kSqrtS = 91.1876;

// Weighted means: m = A/W, var(m) = (B - 2mC + m^2 D) / W^2 with
// A = sum w e, B = sum w^2 e^2, C = sum w^2 e, D = sum w^2.  The KKMC ISR off
// sample is weighted (mean 1.18); the others have unit weights.
struct WMean {
    double A = 0, B = 0, C = 0, D = 0;
    void add(double e, double w) { A += w * e; B += w * w * e * e; C += w * w * e; D += w * w; }
    double mean(double W) const { return W > 0 ? A / W : 0; }
    double err(double W) const { const double m = mean(W), v = B - 2 * m * C + m * m * D;
                                 return (W > 0 && v > 0) ? std::sqrt(v) / W : 0; }
};

struct Scan {
    Long64_t events = 0;
    double sumW = 0;
    WMean sumE, eec;              // per-event charged energy and exact EEC
    double nch = 0;
    double eecHalf[2] = {0, 0}, wHalf[2] = {0, 0};
    bool ok = false;
};

Scan scan(const std::string& path, Long64_t maxEvents)
{
    Scan r;
    TFile* f = TFile::Open(path.c_str());
    if (!f || f->IsZombie()) { std::cerr << "  [skip] " << path << std::endl; return r; }
    TTree* t = static_cast<TTree*>(f->Get("Events"));
    if (!t) { std::cerr << "  [skip] no tree " << path << std::endl; f->Close(); return r; }

    std::vector<char>* isFinal = nullptr;
    std::vector<float>*px = nullptr, *py = nullptr, *pz = nullptr, *en = nullptr, *ch = nullptr;
    double weight = 1.0;
    t->SetBranchStatus("*", 0);
    for (const char* b : {"isFinal", "px", "py", "pz", "energy", "charge", "weight"})
        t->SetBranchStatus(b, 1);
    t->SetBranchAddress("weight", &weight);
    t->SetBranchAddress("isFinal", &isFinal);
    t->SetBranchAddress("px", &px);
    t->SetBranchAddress("py", &py);
    t->SetBranchAddress("pz", &pz);
    t->SetBranchAddress("energy", &en);
    t->SetBranchAddress("charge", &ch);

    Long64_t n = t->GetEntries();
    if (maxEvents > 0 && n > maxEvents) n = maxEvents;
    for (Long64_t i = 0; i < n; ++i) {
        t->GetEntry(i);
        double sE = 0, sE2 = 0;
        int nc = 0;
        for (size_t j = 0; j < isFinal->size(); ++j) {
            if (!(*isFinal)[j]) continue;
            if (std::fabs((*ch)[j]) <= 0) continue;
            const double m = std::sqrt((*px)[j] * (*px)[j] + (*py)[j] * (*py)[j] +
                                       (*pz)[j] * (*pz)[j]);
            if (!(m > 0) || !((*en)[j] > 0)) continue;
            sE += (*en)[j];
            sE2 += (*en)[j] * (*en)[j];
            ++nc;
        }
        const double eec = (sE * sE - sE2) / (2.0 * kSqrtS * kSqrtS);
        const double w = weight;
        r.sumE.add(sE, w); r.eec.add(eec, w); r.nch += w * nc; r.sumW += w;
        const int h = (i < n / 2) ? 0 : 1;
        r.eecHalf[h] += w * eec; r.wHalf[h] += w;
        ++r.events;
    }
    f->Close();
    r.ok = r.events > 0;
    return r;
}

void meanErr(const Scan& s, double& m, double& e)
{
    m = s.eec.mean(s.sumW);
    e = s.eec.err(s.sumW);
}

}  // namespace

void verify_eec_isr_correction(const char* outDir =
                                   "/data2/yjlee/ISRsample/kkmc_1M_20260921/results",
                               Long64_t maxEvents = -1,
                               int sampleIndex = -1)
{
    struct Sample { std::string label; std::string off; std::string on; };
    std::vector<Sample> samples = {
        {"Pythia 8.315", std::string(kRealNtup) + "/mc_Pythia8315_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Pythia8315_ISR_ON.root"},
        {"Pythia 8.315 Vincia", std::string(kRealNtup) + "/mc_Pythia8315_Vincia_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Pythia8315_Vincia_ISR_ON.root"},
        {"Sherpa 3.0.3 PDFESherpa", std::string(kRealNtup) + "/mc_Sherpa303_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Sherpa303_ISR_ON.root"},
        {"Sherpa 3.0.3 YFS", std::string(kRealNtup) + "/mc_Sherpa303_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Sherpa303_ISR_YFS.root"},
        {"Herwig 7.3.0 QED shower, ISR unchanged",
         std::string(kRealNtup) + "/mc_Herwig730_ISR_OFF.root",
         std::string(kRealNtup) + "/mc_Herwig730_QEDshower.root"},
        {"KKMC 4.30 CEEX", std::string(kKkmcNtup) + "/mc_KKMC424_ISR_OFF.root",
         std::string(kKkmcNtup) + "/mc_KKMC424_ISR_ON.root"},
        {"KKMCee 5.00.02 + Pythia 8.316", "/raid5/data/yjlee/ISR/samples/kkmcee5_20261002/mc_KKMCee50002_ISR_OFF.root",
         "/raid5/data/yjlee/ISR/samples/kkmcee5_20261002/mc_KKMCee50002_ISR_ON.root"},
    };
    if (sampleIndex >= 0) {
        if (sampleIndex >= static_cast<int>(samples.size())) return;
        samples = {samples[sampleIndex]};
    }
    const std::string tag = sampleIndex >= 0 ? Form("_s%d", sampleIndex) : "";

    std::ofstream csv(std::string(outDir) + "/eec_verification" + tag + ".csv");
    csv << std::setprecision(10);
    csv << "sample,events,nch_off,nch_on,sumE_off,sumE_on,delta,delta_err,"
           "eec_exact_off,eec_exact_on,c_isr_exact,c_isr_exact_err,"
           "c_isr_predicted,half_diff_off,half_diff_on\n";

    for (const auto& s : samples) {
        std::cout << "[verify] " << s.label << std::endl;
        Scan off = scan(s.off, maxEvents);
        Scan on = scan(s.on, maxEvents);
        if (!off.ok || !on.ok) continue;

        const double nOff = off.sumW, nOn = on.sumW;   // sums of weights
        const double eOff = off.sumE.mean(nOff), eOn = on.sumE.mean(nOn);
        const double delta = (eOff - eOn) / eOff;
        // Statistical error on delta from the spread of the per-event charged
        // energy; the two states are independent samples.
        const double sOffE = off.sumE.err(nOff), sOnE = on.sumE.err(nOn);
        const double deltaErr = (eOn / eOff) *
            std::sqrt((sOffE / eOff) * (sOffE / eOff) + (sOnE / eOn) * (sOnE / eOn));
        double mo, so, mn, sn;
        meanErr(off, mo, so);
        meanErr(on, mn, sn);
        const double r = mo / mn;
        const double er = r * std::sqrt((so / mo) * (so / mo) + (sn / mn) * (sn / mn));
        const double pred = 1.0 / ((1.0 - delta) * (1.0 - delta));

        const double ho = off.eecHalf[0] / off.wHalf[0] - off.eecHalf[1] / off.wHalf[1];
        const double hn = on.eecHalf[0] / on.wHalf[0] - on.eecHalf[1] / on.wHalf[1];

        printf("  events        %lld / %lld\n", off.events, on.events);
        printf("  <N_ch>        %.3f (off)  %.3f (on)\n", off.nch / nOff, on.nch / nOn);
        printf("  <sum E_ch>    %.4f (off)  %.4f (on)  GeV,  delta = %.5f +- %.5f\n",
               eOff, eOn, delta, deltaErr);
        printf("  E_rad         %.4f +- %.4f GeV\n", delta * kSqrtS, deltaErr * kSqrtS);
        printf("  EEC exact     %.8f (off)  %.8f (on)\n", mo, mn);
        printf("  C_ISR exact   %.6f +- %.6f\n", r, er);
        printf("  C_ISR pred    %.6f   (1/(1-delta)^2)\n", pred);
        printf("  half split    off %+.2e   on %+.2e  (quoted sigma %.2e / %.2e)\n",
               ho, hn, so, sn);

        csv << "\"" << s.label << "\"," << off.events << "," << off.nch / nOff << ","
            << on.nch / nOn << "," << eOff << "," << eOn << "," << delta << ","
            << deltaErr << ","
            << mo << "," << mn << "," << r << "," << er << "," << pred << ","
            << ho << "," << hn << "\n";
    }
    csv.close();
    std::cout << "[done] wrote " << outDir << "/eec_verification" << tag << ".csv" << std::endl;
}
