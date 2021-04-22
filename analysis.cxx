#include "ROOT/RDataFrame.hxx"
#include "ROOT/RDFHelpers.hxx"

#include "physicsobjects.hxx"
#include "metfilter.hxx"
#include "pairselection.hxx"

static std::vector<std::string> varSet = {};

int main(){

    ROOT::RDataFrame df("Events", "/ceph/htautau/nanoaod/testing/DYJetsToLL_SUmmer19_UL18.root");
    // for testing, we limit to 1000 events only
    // 1st stage: Good object selection
    std::cout << "Starting Setup of Dataframe \n";
    auto df2 = df.Range(100);
    // MET Filters
    auto df3 = metfilter::ApplyMetFilter(df2, "Flag_goodVertices");

    // Electron

    // Taus
    // Steps for good taus in the mt channel:
    // tau_pt > 30 GeV
    // tau_eta < 2.3
    // Tau_decayMode != 5,6,7,8,9
    // Tau_dz < 0.2
    // Tau_idDeepTau2017v2p1VSjet > 4
    // Tau_idDeepTau2017v2p1VSe > 4
    // Tau_idDeepTau2017v2p1VSmu > 1
    auto df4 = physicsobject::CutPt(df3, "Tau_pt", "tau_maskpt", 30);
    auto df5 = physicsobject::CutEta(df4, "Tau_eta", "tau_masketa", 2.3);
    auto df6 = physicsobject::tau::FilterDecayModes(df5, "tau_maskdecaymodes", {1,2,3,4,10});
    auto df7 = physicsobject::tau::FilterTauID(df6, "tau_maskIDvsJet", "Tau_idDeepTau2017v2p1VSjet", 4);
    auto df8 = physicsobject::tau::FilterTauID(df7, "tau_maskIDvsEle", "Tau_idDeepTau2017v2p1VSe", 4);
    auto df9 = physicsobject::tau::FilterTauID(df8, "tau_maskIDvsMu", "Tau_idDeepTau2017v2p1VSmu", 1);
    auto df10 = physicsobject::CutDz(df9, "Tau_dz", "tau_maskdz", 0.2);

    auto df12 = physicsobject::CombineMasks(df10, "good_taus_mask",{"tau_maskpt", "tau_masketa", "tau_maskdecaymodes", "tau_maskIDvsJet", "tau_maskIDvsEle", "tau_maskIDvsMu","tau_maskdz"});
    // muons
    auto df20 = physicsobject::CutPt(df12, "Muon_pt", "muon_maskpt", 23);
    auto df21 = physicsobject::CutEta(df20, "Muon_eta", "muon_masketa", 2.5);
    auto df22 = physicsobject::muon::FilterID(df21, "muon_maskid", "Muon_mediumId");
    auto df23 = physicsobject::muon::FilterIsolation(df22, "muon_maskiso", "Muon_pfRelIso04_all", 0.15);
    auto df24 = physicsobject::CombineMasks(df23, "good_muons_mask",{"muon_maskpt", "muon_masketa", "muon_maskid", "muon_maskiso"});


    // Build the pair !
    auto df24_1 = physicsobject::FilterObjects(df24, "nTau", 1);
    auto df24_2 = physicsobject::FilterObjects(df24_1, "nMuon", 1);
    auto df25 = pairselection::mutau::PairSelection(df24_2, "good_taus_mask", "good_muons_mask", "mtpair", {""});
    auto df_final = df25;
    auto cutReport = df_final.Report();
    varSet.push_back("good_muons_mask");
    varSet.push_back("good_taus_mask");
    varSet.push_back("mtpair");
    std::cout << "Finished Setup \n";
    std::cout << "Starting Evaluation \n";

    cutReport->Print();
    df_final.Snapshot("ntuple", "test.root", varSet);
    std::cout << "Finished Evaluation \n";
    // as a first testcase, we work on selecting good muons
    return 0;
}