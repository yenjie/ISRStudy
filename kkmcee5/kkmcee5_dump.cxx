// Drive KKMCee 5 for e+e- -> q qbar (+ n gamma) at the Z pole, hadronise the
// quark pair with Pythia 8, and write the KKMC ASCII event dump that
// real_isr_ntuple_producer --mode kkmc already reads:
//
//   E  <iev> <WtMain> <nISR> <nFinal>
//   I  px py pz e                      ISR photons, generator's own list
//   P  kf px py pz e m                 stable final-state particles
//
// KKMCee 5 has no hadronisation of its own: it delivers the bare quark pair
// after ISR and FSR, with the photons listed separately.  Only the quark pair
// goes into Pythia, as a colour-singlet hard process with the shower and
// hadronisation run on top (ProcessLevel:all = off).  The KKMCee photons are
// appended to the final state unchanged, so the ISR photons remain exactly the
// generator's, and the converter matches them back by four-momentum as it does
// for the Fortran 4.30 dump.
//
// Pythia's QED emission off the quarks is switched off because KKMCee already
// radiates from the final quarks (KeyFSR = 1, KeyQSR = 1), as in our 4.30 runs.
//
// Reads ./KKMCee_defaults, ./pro.input and the DIZET tables from the current
// directory, exactly like the stock KKMCee main programs.
//
// Usage:  kkmcee5_dump <nEvents> <outFile> <seed>

#include <cstdio>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include "TH1D.h"
#include "TRandom3.h"
#include "TLorentzVector.h"

#include "KKee2f.h"
#include "KKevent.h"

#include "Pythia8/Pythia.h"

int main(int argc, char** argv)
{
    if (argc < 4) {
        std::cerr << "usage: kkmcee5_dump <nEvents> <outFile> <seed>" << std::endl;
        return 2;
    }
    const long nEvents = std::atol(argv[1]);
    const std::string outName = argv[2];
    const long seed = std::atol(argv[3]);

    std::ofstream log("./pro.output", std::ios::out);
    TRandom3* rn = new TRandom3(seed);
    TH1D* hNorma = new TH1D("HST_KKMC_NORMA", "Normalization histo", 10000, 0, 10000);

    KKee2f* gen = new KKee2f("MCgen");
    gen->Initialize(rn, &log, hNorma);

    Pythia8::Pythia pythia("/usr/share/Pythia8/xmldoc", false);
    pythia.readString("ProcessLevel:all = off");
    pythia.readString("PartonLevel:ISR = off");
    pythia.readString("TimeShower:QEDshowerByQ = off");      // KKMCee does the quark FSR
    pythia.readString("TimeShower:QEDshowerByGamma = off");
    pythia.readString("Check:event = on");
    pythia.readString("Next:numberShowEvent = 0");
    pythia.readString("Next:numberShowInfo = 0");
    pythia.readString("Next:numberShowProcess = 0");
    pythia.readString("Next:numberCount = 0");
    pythia.readString("Random:setSeed = on");
    pythia.readString(("Random:seed = " + std::to_string((seed % 900000000) + 1)).c_str());
    if (!pythia.init()) { std::cerr << "Pythia init failed" << std::endl; return 1; }

    std::ofstream out(outName);
    out << "# KKMCee 5.00.02 ISR-study event dump v1, quark pair hadronised by Pythia "
        << PYTHIA_VERSION << "\n";
    out << "# KeyISR " << gen->m_xpar[20] << " KeyFSR " << gen->m_xpar[21]
        << " KeyINT " << gen->m_xpar[27] << " KeyGPS " << gen->m_xpar[28]
        << " KeyQSR " << gen->m_xpar[29] << " KeyWgt " << gen->m_xpar[10]
        << " CMSene " << gen->m_xpar[1] << "\n";

    long nWritten = 0, nPythiaFail = 0;
    double sumWt = 0;
    char buf[256];
    for (long iev = 1; iev <= nEvents; ++iev) {
        gen->Generate();
        KKevent* ev = gen->m_Event;
        const double wt = gen->m_WtMain;
        const int kf = ev->m_KFfin;             // negative if the first is the antifermion
        const TLorentzVector& q1 = ev->m_Qf1;
        const TLorentzVector& q2 = ev->m_Qf2;

        // Hadronise the quark pair.  The pair is entered as the decay of a
        // colour-singlet Z (status -22) with the quarks as its daughters
        // (status 23), and the event scale set to the pair mass.  Without the
        // scale Pythia starts no shower at all, and without the resonance
        // mother it applies no matrix-element correction to the first gluon
        // emission; either way the multiplicity comes out far too low.  With
        // both, this is the same shower configuration Pythia uses for its own
        // e+e- -> Z -> q qbar, which is what the other samples in this study
        // were made with.  Retry on a (rare) shower/hadronisation failure.
        const TLorentzVector Z = q1 + q2;
        const double scale = Z.M();
        bool ok = false;
        for (int tries = 0; tries < 5 && !ok; ++tries) {
            pythia.event.reset();
            pythia.event.append(23, -22, 0, 0, 2, 3, 0, 0,
                                Z.Px(), Z.Py(), Z.Pz(), Z.E(), Z.M(), scale);
            pythia.event.append( kf, 23, 1, 0, 0, 0, kf > 0 ? 101 : 0, kf > 0 ? 0 : 101,
                                q1.Px(), q1.Py(), q1.Pz(), q1.E(), q1.M(), scale);
            pythia.event.append(-kf, 23, 1, 0, 0, 0, kf > 0 ? 0 : 101, kf > 0 ? 101 : 0,
                                q2.Px(), q2.Py(), q2.Pz(), q2.E(), q2.M(), scale);
            pythia.event.scale(scale);
            ok = pythia.next();
        }
        if (!ok) { ++nPythiaFail; continue; }

        std::vector<int> fin;
        for (int j = 0; j < pythia.event.size(); ++j)
            if (pythia.event[j].isFinal()) fin.push_back(j);
        const int nIsr = ev->m_nPhotISR, nFsr = ev->m_nPhotFSR;
        const int nFinal = static_cast<int>(fin.size()) + nIsr + nFsr;

        std::snprintf(buf, sizeof(buf), "E %10ld %16.9E %4d %6d\n", iev, wt, nIsr, nFinal);
        out << buf;
        for (int i = 1; i <= nIsr; ++i) {
            const TLorentzVector& g = ev->m_PhotISR[i];
            std::snprintf(buf, sizeof(buf), "I  %16.9E %16.9E %16.9E %16.9E\n", g.Px(), g.Py(), g.Pz(), g.E());
            out << buf;
        }
        for (int j : fin) {
            const Pythia8::Particle& p = pythia.event[j];
            std::snprintf(buf, sizeof(buf), "P %8d %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                          p.id(), p.px(), p.py(), p.pz(), p.e(), p.m());
            out << buf;
        }
        for (int i = 1; i <= nIsr; ++i) {
            const TLorentzVector& g = ev->m_PhotISR[i];
            std::snprintf(buf, sizeof(buf), "P %8d %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                          22, g.Px(), g.Py(), g.Pz(), g.E(), 0.0);
            out << buf;
        }
        for (int i = 1; i <= nFsr; ++i) {
            const TLorentzVector& g = ev->m_PhotFSR[i];
            std::snprintf(buf, sizeof(buf), "P %8d %16.9E %16.9E %16.9E %16.9E %16.9E\n",
                          22, g.Px(), g.Py(), g.Pz(), g.E(), 0.0);
            out << buf;
        }
        ++nWritten;
        sumWt += wt;
        if (iev % 50000 == 0)
            std::cout << "  " << iev << "/" << nEvents << "  <wt> so far " << sumWt / nWritten << std::endl;
    }
    out.close();
    gen->Finalize();
    std::cout << "[done] " << nWritten << " events written to " << outName
              << ", Pythia failures " << nPythiaFail
              << ", mean WtMain " << (nWritten ? sumWt / nWritten : 0) << std::endl;
    return 0;
}
