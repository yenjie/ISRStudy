#include <HepMC3/FourVector.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenParticle.h>
#include <HepMC3/GenVertex.h>
#include <HepMC3/Reader.h>
#include <HepMC3/ReaderAscii.h>
#include <HepMC3/ReaderAsciiHepMC2.h>
#include <HepMC3/Units.h>

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

struct Options {
    std::string input;
    std::string output;
    std::string metadata;
    std::string format = "auto";
    std::string generatorName = "Sherpa 3.0.3";
    std::string mode = "UNKNOWN";
    double sqrtS = 91.1876;
    long long maxEvents = -1;
};

struct Summary {
    long long events = 0;
    long long particles = 0;
    long long finalParticles = 0;
    long long isrPhotonCandidates = 0;
    double sumWeights = 0.0;
    double sumWeight2 = 0.0;
    double sumEnergyFinal = 0.0;
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
        if (key == "--input") opt.input = requireValue(i, argc, argv);
        else if (key == "--output") opt.output = requireValue(i, argc, argv);
        else if (key == "--metadata") opt.metadata = requireValue(i, argc, argv);
        else if (key == "--format") opt.format = requireValue(i, argc, argv);
        else if (key == "--generatorName") opt.generatorName = requireValue(i, argc, argv);
        else if (key == "--mode") opt.mode = requireValue(i, argc, argv);
        else if (key == "--sqrtS") opt.sqrtS = std::atof(requireValue(i, argc, argv).c_str());
        else if (key == "--maxEvents") opt.maxEvents = std::atoll(requireValue(i, argc, argv).c_str());
        else if (key == "--help" || key == "-h") {
            std::cout
                << "Usage:\n"
                << "  hepmc3_to_fadgen --input events.hepmc --output sample.fadgen [options]\n\n"
                << "Options:\n"
                << "  --format auto|hepmc2|hepmc3\n"
                << "  --metadata sample.json\n"
                << "  --generatorName NAME\n"
                << "  --mode OFF|PDFESherpa|YFS\n"
                << "  --sqrtS GEV\n"
                << "  --maxEvents N\n\n"
                << "Output is FADGEN_EXCHANGE_ASCII_V1, a documented intermediate\n"
                << "record pending the exact local DELPHI/FADGEN FAD writer spec.\n";
            std::exit(0);
        } else {
            throw std::runtime_error("Unknown argument: " + key);
        }
    }
    if (opt.input.empty()) throw std::runtime_error("--input is required");
    if (opt.output.empty()) throw std::runtime_error("--output is required");
    return opt;
}

std::string jsonEscape(const std::string& s)
{
    std::ostringstream out;
    for (char c : s) {
        switch (c) {
            case '\\': out << "\\\\"; break;
            case '"': out << "\\\""; break;
            case '\n': out << "\\n"; break;
            case '\r': out << "\\r"; break;
            case '\t': out << "\\t"; break;
            default: out << c; break;
        }
    }
    return out.str();
}

std::string unitName(HepMC3::Units::MomentumUnit u)
{
    return u == HepMC3::Units::MEV ? "MEV" : "GEV";
}

std::string unitName(HepMC3::Units::LengthUnit u)
{
    switch (u) {
        case HepMC3::Units::MM: return "MM";
        case HepMC3::Units::CM: return "CM";
        default: return "UNKNOWN";
    }
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

double mass(const HepMC3::FourVector& p)
{
    const double m2 = p.e() * p.e() - p.px() * p.px() - p.py() * p.py() - p.pz() * p.pz();
    return m2 >= 0.0 ? std::sqrt(m2) : -std::sqrt(-m2);
}

bool isISRPhotonCandidate(const HepMC3::ConstGenParticlePtr& p)
{
    if (!p || p->pid() != 22 || p->status() != 1) return false;
    const HepMC3::FourVector& v = p->momentum();
    const double pabs = std::sqrt(v.px() * v.px() + v.py() * v.py() + v.pz() * v.pz());
    if (pabs <= 1e-12 || v.e() < 0.05) return false;
    return std::abs(v.pz() / pabs) > 0.990;
}

std::pair<int, int> firstLastIndices(const std::vector<HepMC3::ConstGenParticlePtr>& particles,
                                     const std::unordered_map<const HepMC3::GenParticle*, int>& index)
{
    int first = 0;
    int last = 0;
    for (const auto& p : particles) {
        const auto it = index.find(p.get());
        if (it == index.end()) continue;
        if (first == 0) first = it->second;
        last = it->second;
    }
    return {first, last};
}

void writeHeader(std::ostream& out, const Options& opt, const std::string& format)
{
    out << "# FADGEN_EXCHANGE_ASCII_V1\n"
        << "# This is a generator-level exchange file for the FADGEN/FAD adapter.\n"
        << "# Replace this writer with the exact DELPHI FAD writer once the local spec is available.\n"
        << "BEGIN_RUN\n"
        << "generator_name " << std::quoted(opt.generatorName) << "\n"
        << "mode " << std::quoted(opt.mode) << "\n"
        << "sqrtS_GeV " << std::setprecision(12) << opt.sqrtS << "\n"
        << "source_hepmc " << std::quoted(opt.input) << "\n"
        << "source_format " << format << "\n"
        << "columns index barcode pdg status mother1 mother2 daughter1 daughter2"
        << " px py pz energy mass vx vy vz vt is_final is_isr_photon_candidate\n"
        << "END_RUN\n";
}

void writeEvent(std::ostream& out, const HepMC3::GenEvent& event, Summary& summary)
{
    const std::vector<HepMC3::ConstGenParticlePtr> particles = event.particles();
    std::unordered_map<const HepMC3::GenParticle*, int> index;
    index.reserve(particles.size());
    for (size_t i = 0; i < particles.size(); ++i) index[particles[i].get()] = static_cast<int>(i + 1);

    const double weight = event.weights().empty() ? 1.0 : event.weights().front();
    summary.events += 1;
    summary.particles += static_cast<long long>(particles.size());
    summary.sumWeights += weight;
    summary.sumWeight2 += weight * weight;

    out << "BEGIN_EVENT "
        << event.event_number() << " "
        << std::setprecision(17) << weight << " "
        << particles.size() << " "
        << unitName(event.momentum_unit()) << " "
        << unitName(event.length_unit()) << "\n";

    for (size_t i = 0; i < particles.size(); ++i) {
        const auto& p = particles[i];
        const HepMC3::FourVector& mom = p->momentum();
        auto prod = p->production_vertex();
        auto end = p->end_vertex();
        const auto mothers = prod ? firstLastIndices(prod->particles_in(), index) : std::make_pair(0, 0);
        const auto daughters = end ? firstLastIndices(end->particles_out(), index) : std::make_pair(0, 0);
        HepMC3::FourVector vtx(0.0, 0.0, 0.0, 0.0);
        if (prod && prod->has_set_position()) vtx = prod->position();
        const bool isFinal = p->status() == 1 || !end;
        const bool isISR = isISRPhotonCandidate(p);
        if (isFinal) {
            summary.finalParticles += 1;
            summary.sumEnergyFinal += mom.e();
        }
        if (isISR) summary.isrPhotonCandidates += 1;

        out << "P "
            << (i + 1) << " "
            << p->id() << " "
            << p->pid() << " "
            << p->status() << " "
            << mothers.first << " "
            << mothers.second << " "
            << daughters.first << " "
            << daughters.second << " "
            << std::setprecision(17)
            << mom.px() << " "
            << mom.py() << " "
            << mom.pz() << " "
            << mom.e() << " "
            << mass(mom) << " "
            << vtx.x() << " "
            << vtx.y() << " "
            << vtx.z() << " "
            << vtx.t() << " "
            << (isFinal ? 1 : 0) << " "
            << (isISR ? 1 : 0) << "\n";
    }
    out << "END_EVENT\n";
}

void writeMetadata(const Options& opt, const std::string& format, const Summary& summary)
{
    if (opt.metadata.empty()) return;
    std::ofstream out(opt.metadata);
    if (!out) throw std::runtime_error("Failed to open metadata output: " + opt.metadata);
    out << "{\n"
        << "  \"format\": \"FADGEN_EXCHANGE_ASCII_V1\",\n"
        << "  \"generatorName\": \"" << jsonEscape(opt.generatorName) << "\",\n"
        << "  \"mode\": \"" << jsonEscape(opt.mode) << "\",\n"
        << "  \"sqrtS_GeV\": " << std::setprecision(12) << opt.sqrtS << ",\n"
        << "  \"input\": \"" << jsonEscape(opt.input) << "\",\n"
        << "  \"output\": \"" << jsonEscape(opt.output) << "\",\n"
        << "  \"sourceHepMCFormat\": \"" << jsonEscape(format) << "\",\n"
        << "  \"events\": " << summary.events << ",\n"
        << "  \"particles\": " << summary.particles << ",\n"
        << "  \"finalParticles\": " << summary.finalParticles << ",\n"
        << "  \"isrPhotonCandidates\": " << summary.isrPhotonCandidates << ",\n"
        << "  \"sumWeights\": " << std::setprecision(17) << summary.sumWeights << ",\n"
        << "  \"sumWeight2\": " << std::setprecision(17) << summary.sumWeight2 << ",\n"
        << "  \"sumEnergyFinal\": " << std::setprecision(17) << summary.sumEnergyFinal << "\n"
        << "}\n";
}

int main(int argc, char** argv)
{
    try {
        const Options opt = parseArgs(argc, argv);
        const std::string format = detectFormat(opt.input, opt.format);
        std::unique_ptr<HepMC3::Reader> reader = makeReader(opt.input, format);
        std::ofstream out(opt.output);
        if (!out) throw std::runtime_error("Failed to open output: " + opt.output);

        writeHeader(out, opt, format);
        Summary summary;
        while (!reader->failed()) {
            HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
            reader->read_event(event);
            if (reader->failed()) break;
            writeEvent(out, event, summary);
            if (opt.maxEvents > 0 && summary.events >= opt.maxEvents) break;
        }
        out << "END_FADGEN_EXCHANGE\n";
        reader->close();
        out.close();
        writeMetadata(opt, format, summary);

        std::cout << "Converted " << summary.events << " events, "
                  << summary.particles << " particles, "
                  << summary.isrPhotonCandidates << " ISR photon candidates -> "
                  << opt.output << std::endl;
        return 0;
    } catch (const std::exception& e) {
        std::cerr << "hepmc3_to_fadgen: " << e.what() << std::endl;
        return 1;
    }
}
