#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/ReaderAsciiHepMC2.h>
#include <HepMC3/Units.h>

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstdlib>
#include <cstring>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

struct Options {
    std::string hepmc;
    std::string fadgen;
    std::string report;
    std::string format = "auto";
    long long maxEvents = -1;
    double absTolerance = 2.0e-5;
    double relTolerance = 2.0e-5;
};

struct LujetsParticle {
    int32_t k[5] = {0, 0, 0, 0, 0};
    float p[5] = {0, 0, 0, 0, 0};
    float v[5] = {0, 0, 0, 0, 0};
};

struct LujetsEvent {
    int32_t n = 0;
    std::vector<LujetsParticle> particles;
};

struct Summary {
    long long eventsCompared = 0;
    long long hepmcFinalParticles = 0;
    long long fadgenParticles = 0;
    long long pdgMismatches = 0;
    long long multiplicityMismatches = 0;
    long long fourVectorMismatches = 0;
    double maxAbsDelta = 0.0;
    double maxRelDelta = 0.0;
};

std::string requireValue(int& i, int argc, char** argv)
{
    if (i + 1 >= argc) throw std::runtime_error(std::string("Missing value for ") + argv[i]);
    return argv[++i];
}

Options parseArgs(int argc, char** argv)
{
    Options opt;
    for (int i = 1; i < argc; ++i) {
        const std::string key = argv[i];
        if (key == "--hepmc") opt.hepmc = requireValue(i, argc, argv);
        else if (key == "--fadgen") opt.fadgen = requireValue(i, argc, argv);
        else if (key == "--report") opt.report = requireValue(i, argc, argv);
        else if (key == "--format") opt.format = requireValue(i, argc, argv);
        else if (key == "--maxEvents") opt.maxEvents = std::atoll(requireValue(i, argc, argv).c_str());
        else if (key == "--absTolerance") opt.absTolerance = std::atof(requireValue(i, argc, argv).c_str());
        else if (key == "--relTolerance") opt.relTolerance = std::atof(requireValue(i, argc, argv).c_str());
        else if (key == "--help" || key == "-h") {
            std::cout
                << "Usage:\n"
                << "  validate_fadgen_vs_hepmc --hepmc events.hepmc --fadgen sample.fadgen [options]\n\n"
                << "Options:\n"
                << "  --format auto|hepmc2|hepmc3\n"
                << "  --report validation.md\n"
                << "  --maxEvents N\n"
                << "  --absTolerance VALUE\n"
                << "  --relTolerance VALUE\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + key);
        }
    }
    if (opt.hepmc.empty()) throw std::runtime_error("--hepmc is required");
    if (opt.fadgen.empty()) throw std::runtime_error("--fadgen is required");
    return opt;
}

std::string detectFormat(const std::string& filename, const std::string& requested)
{
    if (requested != "auto") return requested;
    std::ifstream in(filename);
    std::string line;
    for (int n = 0; n < 20 && std::getline(in, line); ++n) {
        if (line.find("HepMC::Version 2") != std::string::npos) return "hepmc2";
        if (line.find("HepMC::IO_GenEvent") != std::string::npos) return "hepmc2";
        if (line.find("HepMC::Version 3") != std::string::npos) return "hepmc3";
        if (line.find("HepMC::Asciiv3") != std::string::npos) return "hepmc3";
    }
    return "hepmc3";
}

std::unique_ptr<HepMC3::Reader> makeReader(const std::string& filename, const std::string& format)
{
    if (format == "hepmc2") return std::unique_ptr<HepMC3::Reader>(new HepMC3::ReaderAsciiHepMC2(filename));
    if (format == "hepmc3") return std::unique_ptr<HepMC3::Reader>(new HepMC3::ReaderAscii(filename));
    throw std::runtime_error("Unsupported HepMC format: " + format);
}

double invariantMass(const HepMC3::FourVector& p)
{
    const double m2 = p.e() * p.e() - p.px() * p.px() - p.py() * p.py() - p.pz() * p.pz();
    return m2 >= 0.0 ? std::sqrt(m2) : -std::sqrt(-m2);
}

bool isFinalForDelsim(const HepMC3::ConstGenParticlePtr& p)
{
    return p && p->status() == 1 && !p->end_vertex();
}

std::vector<HepMC3::ConstGenParticlePtr> selectedHepMCParticles(const HepMC3::GenEvent& event)
{
    std::vector<HepMC3::ConstGenParticlePtr> selected;
    for (const auto& p : event.particles()) {
        if (isFinalForDelsim(p)) selected.push_back(p);
    }
    return selected;
}

template <typename T>
T readScalar(const std::vector<char>& data, size_t& offset)
{
    if (offset + sizeof(T) > data.size()) throw std::runtime_error("Malformed FADGEN record");
    T value;
    std::memcpy(&value, data.data() + offset, sizeof(T));
    offset += sizeof(T);
    return value;
}

bool readFortranRecord(std::istream& in, std::vector<char>& payload)
{
    int32_t nbytes = 0;
    in.read(reinterpret_cast<char*>(&nbytes), sizeof(nbytes));
    if (!in) return false;
    if (nbytes <= 0 || nbytes > 10000000) {
        throw std::runtime_error("Invalid Fortran record size in FADGEN file");
    }
    payload.resize(static_cast<size_t>(nbytes));
    in.read(payload.data(), nbytes);
    if (!in) throw std::runtime_error("Truncated FADGEN payload");
    int32_t trailer = 0;
    in.read(reinterpret_cast<char*>(&trailer), sizeof(trailer));
    if (!in) throw std::runtime_error("Missing FADGEN record trailer");
    if (trailer != nbytes) throw std::runtime_error("Mismatched Fortran record markers in FADGEN file");
    return true;
}

LujetsEvent decodeLujetsRecord(const std::vector<char>& payload)
{
    LujetsEvent event;
    size_t offset = 0;
    event.n = readScalar<int32_t>(payload, offset);
    if (event.n <= 0 || event.n > 4000) throw std::runtime_error("Invalid LUJETS particle count");
    const size_t expected = sizeof(int32_t) + static_cast<size_t>(event.n) * 15U * sizeof(int32_t);
    if (payload.size() != expected) throw std::runtime_error("Unexpected LUJETS record byte count");
    event.particles.resize(event.n);
    for (int i = 0; i < event.n; ++i) {
        for (int j = 0; j < 5; ++j) event.particles[i].k[j] = readScalar<int32_t>(payload, offset);
        for (int j = 0; j < 5; ++j) event.particles[i].p[j] = readScalar<float>(payload, offset);
        for (int j = 0; j < 5; ++j) event.particles[i].v[j] = readScalar<float>(payload, offset);
    }
    return event;
}

double relativeDelta(double a, double b)
{
    const double scale = std::max({1.0, std::abs(a), std::abs(b)});
    return std::abs(a - b) / scale;
}

bool outsideTolerance(double a, double b, const Options& opt, Summary& summary)
{
    const double ad = std::abs(a - b);
    const double rd = relativeDelta(a, b);
    summary.maxAbsDelta = std::max(summary.maxAbsDelta, ad);
    summary.maxRelDelta = std::max(summary.maxRelDelta, rd);
    return ad > opt.absTolerance && rd > opt.relTolerance;
}

bool compareEvent(const HepMC3::GenEvent& hepmcEvent,
                  const LujetsEvent& fadgenEvent,
                  const Options& opt,
                  Summary& summary,
                  std::ostream& detail)
{
    const auto hepmcParticles = selectedHepMCParticles(hepmcEvent);
    const int expectedN = static_cast<int>(hepmcParticles.size()) + 1;
    bool ok = true;

    summary.eventsCompared += 1;
    summary.hepmcFinalParticles += static_cast<long long>(hepmcParticles.size());
    summary.fadgenParticles += fadgenEvent.n - 1;

    if (fadgenEvent.n != expectedN) {
        summary.multiplicityMismatches += 1;
        detail << "- Event " << summary.eventsCompared << ": multiplicity mismatch, HepMC final "
               << hepmcParticles.size() << ", FADGEN payload " << (fadgenEvent.n - 1) << "\n";
        ok = false;
    }

    const size_t nCompare = std::min(hepmcParticles.size(), fadgenEvent.particles.size() > 0 ? fadgenEvent.particles.size() - 1 : 0);
    for (size_t i = 0; i < nCompare; ++i) {
        const auto& hp = hepmcParticles[i];
        const auto& fp = fadgenEvent.particles[i];
        if (fp.k[0] != 1 || fp.k[1] != hp->pid()) {
            summary.pdgMismatches += 1;
            if (summary.pdgMismatches <= 20) {
                detail << "- Event " << summary.eventsCompared << ", particle " << (i + 1)
                       << ": K mismatch, FADGEN (" << fp.k[0] << ", " << fp.k[1]
                       << "), HepMC pid " << hp->pid() << "\n";
            }
            ok = false;
        }
        const HepMC3::FourVector mom = hp->momentum();
        const double hepmcVals[5] = {mom.px(), mom.py(), mom.pz(), mom.e(), invariantMass(mom)};
        for (int j = 0; j < 5; ++j) {
            if (outsideTolerance(static_cast<double>(fp.p[j]), hepmcVals[j], opt, summary)) {
                summary.fourVectorMismatches += 1;
                if (summary.fourVectorMismatches <= 20) {
                    detail << "- Event " << summary.eventsCompared << ", particle " << (i + 1)
                           << ", component " << j << ": FADGEN " << std::setprecision(10) << fp.p[j]
                           << ", HepMC " << hepmcVals[j] << "\n";
                }
                ok = false;
            }
        }
    }

    if (!fadgenEvent.particles.empty()) {
        const auto& comment = fadgenEvent.particles.back();
        if (comment.k[0] != 21 || comment.k[1] != 0 || comment.k[2] != 303) {
            detail << "- Event " << summary.eventsCompared << ": missing/invalid generator comment line\n";
            ok = false;
        }
    }
    return ok;
}

void writeReport(const Options& opt,
                 const std::string& format,
                 const Summary& summary,
                 bool pass,
                 const std::string& details)
{
    if (opt.report.empty()) return;
    std::ofstream out(opt.report);
    if (!out) throw std::runtime_error("Failed to write report: " + opt.report);
    out << "# FADGEN vs HepMC Generator-Level Validation\n\n"
        << "- status: " << (pass ? "PASS" : "FAIL") << "\n"
        << "- HepMC: `" << opt.hepmc << "`\n"
        << "- FADGEN: `" << opt.fadgen << "`\n"
        << "- HepMC format: " << format << "\n"
        << "- events compared: " << summary.eventsCompared << "\n"
        << "- HepMC final particles: " << summary.hepmcFinalParticles << "\n"
        << "- FADGEN particles excluding comment lines: " << summary.fadgenParticles << "\n"
        << "- multiplicity mismatches: " << summary.multiplicityMismatches << "\n"
        << "- PDG/status mismatches: " << summary.pdgMismatches << "\n"
        << "- four-vector/mass mismatches: " << summary.fourVectorMismatches << "\n"
        << "- max absolute delta: " << std::setprecision(10) << summary.maxAbsDelta << "\n"
        << "- max relative delta: " << std::setprecision(10) << summary.maxRelDelta << "\n"
        << "- absolute tolerance: " << opt.absTolerance << "\n"
        << "- relative tolerance: " << opt.relTolerance << "\n\n";
    if (!details.empty()) out << "## First Differences\n\n" << details;
}

int main(int argc, char** argv)
{
    try {
        const Options opt = parseArgs(argc, argv);
        const std::string format = detectFormat(opt.hepmc, opt.format);
        std::unique_ptr<HepMC3::Reader> reader = makeReader(opt.hepmc, format);
        std::ifstream fadgen(opt.fadgen, std::ios::binary);
        if (!fadgen) throw std::runtime_error("Failed to open FADGEN file: " + opt.fadgen);

        Summary summary;
        bool pass = true;
        std::ostringstream details;
        long long eventLimit = opt.maxEvents;

        while (!reader->failed()) {
            if (eventLimit > 0 && summary.eventsCompared >= eventLimit) break;
            HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
            reader->read_event(event);
            if (reader->failed()) break;
            std::vector<char> payload;
            if (!readFortranRecord(fadgen, payload)) {
                details << "- Missing FADGEN record for HepMC event " << event.event_number() << "\n";
                pass = false;
                break;
            }
            LujetsEvent lujets = decodeLujetsRecord(payload);
            if (!compareEvent(event, lujets, opt, summary, details)) pass = false;
        }

        std::vector<char> extra;
        if ((eventLimit < 0 || summary.eventsCompared < eventLimit) && readFortranRecord(fadgen, extra)) {
            details << "- FADGEN contains extra records after HepMC stream ended\n";
            pass = false;
        }
        reader->close();

        writeReport(opt, format, summary, pass, details.str());
        std::cout << (pass ? "PASS" : "FAIL")
                  << " events=" << summary.eventsCompared
                  << " hepmcFinal=" << summary.hepmcFinalParticles
                  << " fadgenParticles=" << summary.fadgenParticles
                  << " maxAbsDelta=" << std::setprecision(10) << summary.maxAbsDelta
                  << " maxRelDelta=" << std::setprecision(10) << summary.maxRelDelta
                  << std::endl;
        return pass ? 0 : 1;
    } catch (const std::exception& e) {
        std::cerr << "validate_fadgen_vs_hepmc: " << e.what() << std::endl;
        return 1;
    }
}
