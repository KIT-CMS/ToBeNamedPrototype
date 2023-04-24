#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "RooTrace.h"
#include "TStopwatch.h"
#include "include/genparticles.hxx"
#include "include/htxs.hxx"
#include "include/jets.hxx"
#include "include/lorentzvectors.hxx"
#include "include/met.hxx"
#include "include/metfilter.hxx"
#include "include/pairselection.hxx"
#include "include/tripleselection.hxx"
#include "include/physicsobjects.hxx"
#include "include/quantities.hxx"
#include "include/reweighting.hxx"
#include "include/scalefactors.hxx"
#include "include/triggers.hxx"
#include "include/utility/Logger.hxx"
#include <ROOT/RLogger.hxx>
#include <TFile.h>
#include <TTree.h>
#include <regex>
#include <string>

// {INCLUDES}

int main(int argc, char *argv[]) {
    // {DEBUGLEVEL}
    // ROOT logging
    if (debug) {
        auto verbosity = ROOT::Experimental::RLogScopedVerbosity(
            ROOT::Detail::RDF::RDFLogChannel(),
            ROOT::Experimental::ELogLevel::kInfo);
        RooTrace::verbose(kTRUE);
        Logger::setLevel(Logger::LogLevel::DEBUG);
    } else {
        RooTrace::verbose(kFALSE);
        Logger::setLevel(Logger::LogLevel::INFO);
        gErrorIgnoreLevel = 6001; // ignore all ROOT errors
    }
    if (argc < 3) {
        Logger::get("main")->critical(
            "Require at least two arguments: a single output file and N input "
            "files \n"
            "Example:\n"
            "./analysis output.root /path/to/inputfiles/*.root");
        return 1;
    }
    std::vector<std::string> input_files;
    int nevents = 0;
    Logger::get("main")->info("Checking input files");
    for (int i = 2; i < argc; i++) {
        input_files.push_back(std::string(argv[i]));
        // Check if the input file exists and is readable, also get the number
        // of events, using a timeout of 30 seconds
        TFile *f1 = TFile::Open(argv[i],"TIMEOUT=30");
        if (!f1 || f1->IsZombie()) {
            Logger::get("main")->critical("File {} does not exist or is not "
                                          "readable",
                                          argv[i]);
            return 1;
        }
        TTree *t1 = (TTree *)f1->Get("Events");
        nevents += t1->GetEntries();
        Logger::get("main")->info("input_file {}: {} - {} Events", i - 1,
                                  argv[i], t1->GetEntries());
    }
    const auto output_path = argv[1];
    Logger::get("main")->info("Output directory: {}", output_path);
    TStopwatch timer;
    timer.Start();
    int quantile = 10000;

    // file logging
    Logger::enableFileLogging("logs/main.txt");

    // {MULTITHREADING}

    // initialize df
    ROOT::RDataFrame df0("Events", input_files);
    Logger::get("main")->info("Starting Setup of Dataframe with {} events",
                              nevents);

    // {CODE_GENERATION}

    ROOT::RDF::RSnapshotOptions dfconfig;
    dfconfig.fLazy = true;

    // {RUN_COMMANDS}

    // Add meta-data
    // clang-format off
    const std::map<std::string, std::vector<std::string>> output_quanties = {OUTPUT_QUANTITIES};
    const std::map<std::string, std::vector<std::string>> variations = {SYSTEMATIC_VARIATIONS};
    // clang-format on
    const std::string analysis = {ANALYSISTAG};
    const std::string era = {ERATAG};
    const std::string sample = {SAMPLETAG};
    const std::string commit_hash = {COMMITHASH};
    bool setup_clean = {SETUP_IS_CLEAN};
    for (auto const &x : output_quanties) {
        TFile outputfile(x.first.c_str(), "UPDATE");
        TTree quantities_meta = TTree("quantities", "quantities");
        for (auto const &quantity : x.second) {
            quantities_meta.Branch(quantity.c_str(), &setup_clean);
        }
        quantities_meta.Write();
        TTree variations_meta = TTree("variations", "variations");
        for (auto const &variation : variations.at(x.first)) {
            variations_meta.Branch(variation.c_str(), &setup_clean);
        }
        variations_meta.Write();
        TTree conditions_meta = TTree("conditions", "conditions");
        conditions_meta.Branch(analysis.c_str(), &setup_clean);
        conditions_meta.Branch(era.c_str(), &setup_clean);
        conditions_meta.Branch(sample.c_str(), &setup_clean);
        conditions_meta.Write();
        TTree commit_meta = TTree("commit", "commit");
        commit_meta.Branch(commit_hash.c_str(), &setup_clean);
        commit_meta.Fill();
        commit_meta.Write();
        outputfile.Close();
    }

    Logger::get("main")->info("Finished Evaluation");
    Logger::get("main")->info(
        "Overall runtime (real time: {0:.2f}, CPU time: {1:.2f})",
        timer.RealTime(), timer.CpuTime());
}
