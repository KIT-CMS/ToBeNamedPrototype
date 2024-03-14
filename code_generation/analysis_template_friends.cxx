#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "RooTrace.h"
#include "TStopwatch.h"
#include "include/fakefactors.hxx"
#include "include/genparticles.hxx"
#include "include/htxs.hxx"
#include "include/jets.hxx"
#include "include/lorentzvectors.hxx"
#include "include/met.hxx"
#include "include/metfilter.hxx"
#include "include/ml.hxx"
#include "include/pairselection.hxx"
#include "include/physicsobjects.hxx"
#include "include/quantities.hxx"
#include "include/reweighting.hxx"
#include "include/scalefactors.hxx"
#include "include/topreco.hxx"
#include "include/triggers.hxx"
#include "include/tripleselection.hxx"
#include "include/utility/Logger.hxx"
#include <ROOT/RLogger.hxx>
#include <TFile.h>
#include <TMap.h>
#include <TObjString.h>
#include <TTree.h>
#include <TVector.h>
#include <filesystem>
#include <regex>
#include <string>
// {INCLUDES}

int validate_rootfile(std::string file, std::string &basetree) {
    int nevents = 0;
    TFile *f1 = TFile::Open(file.c_str(), "TIMEOUT=30");
    if (!f1 || f1->IsZombie()) {
        Logger::get("main")->critical("File {} does not exist or is not "
                                      "readable",
                                      file);
        return -1;
    }
    // Get a list of all keys in the file
    TList *list = f1->GetListOfKeys();
    // Check if the Events tree exists
    if (list->FindObject("Events")) {
        TTree *t1 = (TTree *)f1->Get("Events");
        nevents += t1->GetEntries();
        Logger::get("main")->info("NanoAOD input_file: {} - {} Events", file,
                                  t1->GetEntries());
        return nevents;
    } else if (list->FindObject("ntuple")) {
        TTree *t1 = (TTree *)f1->Get("ntuple");
        nevents += t1->GetEntries();
        basetree = "ntuple";
        Logger::get("main")->info("CROWN input_file: {} - {} Events", file,
                                  t1->GetEntries());
        return nevents;
    } else {
        Logger::get("main")->critical("File {} does not contain a tree "
                                      "named 'Events' or 'ntuple'",
                                      file);
        return -1;
    }
}

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
            "Require at least two arguments: a single output file a single "
            "input file and N additioal friend files"
            "files \n"
            "Example:\n"
            "./analysis output.root input.root input_friend_1.root "
            "input_friend_2.root");
        return 1;
    }
    // check if CROWN is run from the correct directory, if the folder "data"
    // does not exist, exit
    if (!std::filesystem::exists("data")) {
        Logger::get("main")->critical(
            "CROWN is not run from the correct directory, "
            "data folder does not exist. Did you run CROWN from the correct "
            "directory?");
        return 1;
    }
    std::vector<std::string> input_friends;
    int nevents = 0;
    Logger::get("main")->info("Checking if all input files exist, are readable "
                              "and have the same number of events");
    std::string basetree = "Events";
    std::string input_file = argv[2];
    nevents = validate_rootfile(input_file, basetree);
    int nfriends = argc - 3;
    if (nevents == -1) {
        return 1;
    }

    for (int i = 3; i < argc; i++) {
        input_friends.push_back(std::string(argv[i]));
        int nevents_friend = validate_rootfile(input_friends.back(), basetree);
        if (nevents_friend == -1) {
            return 1;
        }
        if (nevents != nevents_friend) {
            Logger::get("main")->critical(
                "Input file {} has {} events, but input friend {} has {} "
                "events. All input files must have the same number of events.",
                input_file, nevents, input_friends.back(), nevents_friend);
            return 1;
        }
    }
    const auto output_path = argv[1];
    Logger::get("main")->info("Output directory: {}", output_path);
    TStopwatch timer;
    timer.Start();
    int quantile = 10000;

    // file logging
    Logger::enableFileLogging("logs/main.txt");

    // {MULTITHREADING}

    // build a tchain from input file with all friends
    TChain dataset(basetree.c_str());
    std::vector<TChain *> friends;
    dataset.Add(input_file.c_str());
    for (int index = 0; index < nfriends; index++) {
        std::string friend_file = input_friends[index];
        std::string chain_name = "friend_" + std::to_string(index);
        // each friend needs to be a separate TChain
        TChain *friend_chain = new TChain(basetree.c_str());
        Logger::get("main")->debug("Adding friend file {}", friend_file);
        friend_chain->Add(friend_file.c_str());
        // then add the friend to the main dataset
        dataset.AddFriend(friend_chain, chain_name.c_str());
        friends.push_back(friend_chain);
    }
    // initialize df
    ROOT::RDataFrame df0(dataset);
    // ROOT::RDF::Experimental::AddProgressBar(df0); ROOT 6.30 not available for CS8 on lcg
    // print all available branches to the log
    Logger::get("main")->debug("Available branches:");
    for (auto const &branch : df0.GetColumnNames()) {
        Logger::get("main")->debug("{}", branch);
    }
    Logger::get("main")->info(
        "Starting Setup of Dataframe with {} events and {} friends", nevents,
        nfriends);
    std::vector<ROOT::RDF::RResultPtr<ROOT::RDF::RCutFlowReport>> cutReports;
    // setup output files
    // {OUTPUT_PATHS}
    if (nevents == 0) {
        // {ZERO_EVENTS_FALLBACK}
    } else {
        // {CODE_GENERATION}
        ROOT::RDF::RSnapshotOptions dfconfig;
        dfconfig.fLazy = true;
        // {RUN_COMMANDS}
    }
    // Add meta-data
    // clang-format off
    const std::map<std::string, std::vector<std::string>> output_quanties = {OUTPUT_QUANTITIES};
    const std::map<std::string, std::vector<std::string>> variations = {SYSTEMATIC_VARIATIONS};
    std::map<std::string, std::map<std::string, std::vector<std::string>>> shift_quantities_map = {SHIFT_QUANTITIES_MAP};
    std::map<std::string, std::map<std::string, std::vector<std::string>>> quantities_shift_map = {QUANTITIES_SHIFT_MAP};
    // clang-format on
    const std::string analysis = {ANALYSISTAG};
    const std::string config = {CONFIGTAG};
    const std::string era = {ERATAG};
    const std::string sample = {SAMPLETAG};
    const std::string commit_hash = {COMMITHASH};
    bool setup_clean = {CROWN_IS_CLEAN};
    const std::string analysis_commit_hash = {ANALYSIS_COMMITHASH};
    bool analysis_setup_clean = {ANALYSIS_IS_CLEAN};
    int scope_counter = 0;
    for (auto const &output : output_quanties) {
        // output.first is the output file name
        // output.second is the list of quantities
        TFile outputfile(output.first.c_str(), "UPDATE");
        TTree quantities_meta = TTree("quantities", "quantities");
        for (auto const &quantity : output.second) {
            quantities_meta.Branch(quantity.c_str(), &setup_clean);
        }
        quantities_meta.Write();
        TTree variations_meta = TTree("variations", "variations");
        for (auto const &variation : variations.at(output.first)) {
            variations_meta.Branch(variation.c_str(), &setup_clean);
        }
        variations_meta.Write();
        TTree conditions_meta = TTree("conditions", "conditions");
        conditions_meta.Branch(analysis.c_str(), &setup_clean);
        conditions_meta.Branch(config.c_str(), &setup_clean);
        conditions_meta.Branch(era.c_str(), &setup_clean);
        conditions_meta.Branch(sample.c_str(), &setup_clean);
        conditions_meta.Write();
        TTree commit_meta = TTree("commit", "commit");
        commit_meta.Branch(commit_hash.c_str(), &setup_clean);
        commit_meta.Fill();
        commit_meta.Write();
        TTree analysis_commit_meta =
            TTree("analysis_commit", "analysis_commit");
        analysis_commit_meta.Branch(analysis_commit_hash.c_str(),
                                    &analysis_setup_clean);
        analysis_commit_meta.Fill();
        analysis_commit_meta.Write();
        if (nevents != 0) {
            TH1D cutflow;
            cutflow.SetName("cutflow");
            cutflow.SetTitle("cutflow");
            // iterate through the cutflow vector and fill the histogram with
            // the .GetPass() values
            if (cutReports.size() < scope_counter || cutReports.empty()) {
                Logger::get("main")->critical(
                    "cutReports vector is too small, this should not happen");
                return 1;
            }
            for (auto cut = cutReports[scope_counter].begin();
                 cut != cutReports[scope_counter].end(); cut++) {
                cutflow.SetBinContent(
                    std::distance(cutReports[scope_counter].begin(), cut) + 1,
                    cut->GetPass());
                cutflow.GetXaxis()->SetBinLabel(
                    std::distance(cutReports[scope_counter].begin(), cut) + 1,
                    cut->GetName().c_str());
            }
            // store it in the output file
            cutflow.Write();
        }
        outputfile.Close();
        TFile *fout = TFile::Open(output.first.c_str(), "UPDATE");
        Logger::get("main")->info("Writing quantities map to {}", output.first);
        fout->WriteObject(&shift_quantities_map.at(output.first),
                          "shift_quantities_map");
        fout->WriteObject(&quantities_shift_map.at(output.first),
                          "quantities_shift_map");
        fout->Close();
        scope_counter++;
    }

    Logger::get("main")->info("Finished Evaluation");
    Logger::get("main")->info(
        "Overall runtime (real time: {0:.2f}, CPU time: {1:.2f})",
        timer.RealTime(), timer.CpuTime());
}
