// Cross-check: does our generator charged-EEC density have the same shape as
// the one the ALEPH inclusive EEC analysis reports?
//
// Both are charged-particle EECs in z = (1 - cos theta_L)/2 on the same 200-bin
// double-logarithmic grid, so after normalising each to unit integral the
// shapes should agree up to genuine generator-vs-data differences.  A gross
// mismatch would mean our EEC is wrong, not merely differently normalised.

#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

void compare_eec_shape_to_aleph(
    const char* alephFile =
        "/raid5/data/yjlee/ALEPH_EEC/reproduction/output/report_20260910/"
        "joint_3d_direct_final/final_nominal.root",
    const char* alephHist = "reported_iter1_inclusive_z",
    const char* ourCsv =
        "/data2/yjlee/ISRsample/kkmc_1M_20260921/results/eec_isr_correction.csv",
    const char* ourSample = "Pythia 8.315")
{
    TFile* f = TFile::Open(alephFile);
    if (!f || f->IsZombie()) { printf("cannot open %s\n", alephFile); return; }
    TH1* h = dynamic_cast<TH1*>(f->Get(alephHist));
    if (!h) {
        printf("histogram %s not found; keys are:\n", alephHist);
        f->ls();
        return;
    }
    const int nb = h->GetNbinsX();
    printf("ALEPH histogram %s: %d bins, x from %g to %g\n",
           alephHist, nb, h->GetXaxis()->GetXmin(), h->GetXaxis()->GetXmax());

    // Our per-bin sums and exact edges.
    // Our file always has the full 200-bin grid, whatever nb the reported
    // histogram was merged down to.
    const int kOurBins = 200;
    std::vector<double> ours(kOurBins, 0.0), zlo(kOurBins, 0.0), zhi(kOurBins, 0.0);
    FILE* fp = fopen(ourCsv, "r");
    if (!fp) { printf("cannot open %s\n", ourCsv); return; }
    char line[4096];
    fgets(line, sizeof(line), fp);
    while (fgets(line, sizeof(line), fp)) {
        std::string l(line);
        const std::string key = std::string(ourSample) + ",";
        const size_t p = l.find(key);
        if (p == std::string::npos) continue;
        std::string rest = l.substr(p + key.size());
        int bin = 0;
        double a = 0, b = 0, eoff = 0, eofferr = 0, eon = 0;
        if (sscanf(rest.c_str(), "%d,%lf,%lf,%lf,%lf,%lf", &bin, &a, &b, &eoff,
                   &eofferr, &eon) != 6) continue;
        if (bin < 0 || bin >= kOurBins) continue;
        ours[bin] = eon;
        zlo[bin] = a;
        zhi[bin] = b;
    }
    fclose(fp);

    // Their reported spectrum is merged, so group our 200 bins to match.
    const int merge = kOurBins / nb;
    printf("merging our 200 bins in groups of %d to match\n\n", merge);
    std::vector<double> mo(nb, 0.0), mlo(nb, 0.0), mhi(nb, 0.0);
    for (int i = 0; i < nb; ++i) {
        mlo[i] = zlo[i * merge];
        mhi[i] = zhi[i * merge + merge - 1];
        for (int k = 0; k < merge; ++k) mo[i] += ours[i * merge + k];
    }

    // Two readings of their histogram: bin content already a density, or content
    // a per-bin sum that still needs dividing by the width.  Normalise each to
    // unit integral over z and see which one tracks ours.
    double normAdens = 0, normAsum = 0, sumO = 0;
    for (int i = 0; i < nb; ++i) {
        const double w = mhi[i] - mlo[i];
        normAdens += h->GetBinContent(i + 1) * w;   // if content is a density
        normAsum += h->GetBinContent(i + 1);        // if content is a sum
        sumO += mo[i];
    }
    printf("%4s %12s %12s %12s %12s %9s %9s\n", "bin", "z_centre", "width",
           "ours dens", "A dens(raw)", "r_raw", "r_div");
    for (int i = 0; i < nb; ++i) {
        const double w = mhi[i] - mlo[i];
        if (!(w > 0) || sumO <= 0) continue;
        const double dour = mo[i] / sumO / w;
        const double aRaw = h->GetBinContent(i + 1) / normAdens;       // content = density
        const double aDiv = h->GetBinContent(i + 1) / normAsum / w;    // content = sum
        printf("%4d %12.4e %12.4e %12.4e %12.4e %9.3f %9.3f\n", i,
               0.5 * (mlo[i] + mhi[i]), w, dour, aRaw,
               aRaw > 0 ? dour / aRaw : 0.0, aDiv > 0 ? dour / aDiv : 0.0);
    }
    f->Close();
}
