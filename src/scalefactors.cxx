#ifndef GUARD_SCALEFACTORS_H
#define GUARD_SCALEFACTORS_H

#include "../include/basefunctions.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/RooFunctorThreadsafe.hxx"
#include "ROOT/RDataFrame.hxx"
#include "RooFunctor.h"
#include "RooWorkspace.h"
#include "TFile.h"
#include "correction.h"
/// namespace used for scale factor related functions
namespace scalefactor {
namespace muon {
/**
 * @brief Function used to evaluate id scale factors from muons
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param id_output name of the id scale factor column
 * @param workspace_name path to the Rooworkspace
 * @param id_functor_name name of the function from the workspace
 * @param id_arguments arguments of the function
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id_rooworkspace(ROOT::RDF::RNode df, const std::string &pt,
                                 const std::string &eta,
                                 const std::string &id_output,
                                 const std::string &workspace_name,
                                 const std::string &id_functor_name,
                                 const std::string &id_arguments) {

    Logger::get("muonsf")->debug("Setting up functions for muon sf");
    Logger::get("muonsf")->debug("ID - Function {} // argset {}",
                                 id_functor_name, id_arguments);

    const std::shared_ptr<RooFunctorThreadsafe> id_function =
        loadFunctor(workspace_name, id_functor_name, id_arguments);
    auto df1 = basefunctions::evaluateWorkspaceFunction(df, id_output,
                                                        id_function, pt, eta);
    return df1;
}
/**
 * @brief Function used to evaluate iso scale factors from muons
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param iso muon iso
 * @param iso_output name of the iso scale factor column
 * @param workspace_name path to the Rooworkspace
 * @param iso_functor_name name of the function from the workspace
 * @param iso_arguments arguments of the function
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode iso_rooworkspace(ROOT::RDF::RNode df, const std::string &pt,
                                  const std::string &eta,
                                  const std::string &iso,
                                  const std::string &iso_output,
                                  const std::string &workspace_name,
                                  const std::string &iso_functor_name,
                                  const std::string &iso_arguments) {

    Logger::get("muonsf")->debug("Setting up functions for muon sf");
    Logger::get("muonsf")->debug("Iso - Function {} // argset {}",
                                 iso_functor_name, iso_arguments);

    const std::shared_ptr<RooFunctorThreadsafe> iso_function =
        loadFunctor(workspace_name, iso_functor_name, iso_arguments);
    auto df1 = basefunctions::evaluateWorkspaceFunction(
        df, iso_output, iso_function, pt, eta, iso);
    return df1;
}
/**
 * @brief Function used to evaluate id scale factors from muons with
 * correctionlib. Configuration:
 * - [UL2018 Muon
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/MUO_muon_Z_Run2_UL/MUO_muon_Z_2018_UL.html)
 * - [UL2017 Muon
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/MUO_muon_Z_Run2_UL/MUO_muon_Z_2017_UL.html)
 * - [UL2016preVFP Muon
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/MUO_muon_Z_Run2_UL/MUO_muon_Z_2016preVFP_UL.html)
 * - [UL2016postVFP Muon
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/MUO_muon_Z_Run2_UL/MUO_muon_Z_2016postVFP_UL.html)
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param year_id id for the year of data taking and mc compaign
 * @param variation id for the variation of the scale factor "sf" for nominal
 * and "systup"/"systdown" for up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the muon scale factors
 * @param idAlgorithm name of the muon id scale factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &pt,
                    const std::string &eta, const std::string &year_id,
                    const std::string &variation, const std::string &id_output,
                    const std::string &sf_file,
                    const std::string &idAlgorithm) {

    Logger::get("muonIdSF")->debug("Setting up functions for muon id sf");
    Logger::get("muonIdSF")->debug("ID - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        id_output,
        [evaluator, year_id, variation](const float &pt, const float &eta) {
            Logger::get("muonIdSF")->debug("ID - pt {}, eta {}", pt, eta);
            double sf = 1.;
            // preventing muons with default values due to tau energy correction
            // shifts below good tau pt selection
            if (pt >= 0.0 && std::abs(eta) >= 0.0) {
                sf = evaluator->evaluate(
                    {year_id, std::abs(eta), pt, variation});
            }
            return sf;
        },
        {pt, eta});
    return df1;
}
/**
 * @brief Function used to evaluate iso scale factors from muons with
 * correctionlib. Configurations:
 * - [UL2018 Muon
 * Iso](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/MUO_muon_Z_Run2_UL/MUO_muon_Z_2018_UL.html)
 * - [UL2017 Muon
 * Iso](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/MUO_muon_Z_Run2_UL/MUO_muon_Z_2017_UL.html)
 * - [UL2016preVFP Muon
 * Iso](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/MUO_muon_Z_Run2_UL/MUO_muon_Z_2016preVFP_UL.html)
 * - [UL2016postVFP Muon
 * Iso](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/MUO_muon_Z_Run2_UL/MUO_muon_Z_2016postVFP_UL.html)
 *
 * @param df The input dataframe
 * @param pt muon pt
 * @param eta muon eta
 * @param year_id id for the year of data taking and mc compaign
 * @param variation id for the variation of the scale factor "sf" for nominal
 * and "systup"/"systdown" the up/down variation
 * @param iso_output name of the iso scale factor column
 * @param sf_file path to the file with the muon scale factors
 * @param idAlgorithm name of the muon iso scale factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode iso(ROOT::RDF::RNode df, const std::string &pt,
                     const std::string &eta, const std::string &year_id,
                     const std::string &variation,
                     const std::string &iso_output, const std::string &sf_file,
                     const std::string &idAlgorithm) {

    Logger::get("muonIsoSF")->debug("Setting up functions for muon iso sf");
    Logger::get("muonIsoSF")->debug("ISO - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        iso_output,
        [evaluator, year_id, variation](const float &pt, const float &eta) {
            Logger::get("muonIsoSF")->debug("ISO - pt {}, eta {}", pt, eta);
            double sf = 1.;
            // preventing muons with default values due to tau energy correction
            // shifts below good tau pt selection
            if (pt >= 0.0 && std::abs(eta) >= 0.0) {
                sf = evaluator->evaluate(
                    {year_id, std::abs(eta), pt, variation});
            }
            return sf;
        },
        {pt, eta});
    return df1;
}
} // namespace muon
namespace tau {
/**
 * @brief Function used to evaluate vsJets tau id scale factors in the lt
channel with
 * correctionlib

Description of the bit map used to define the tau id working points of the
DeepTau2017v2p1 tagger.
vsJets                              | Value | Bit (value used in the config)
------------------------------------|-------|-------
no ID selection (takes every tau)   |  0    | -
VVVLoose                            |  1    | 1
VVLoose                             |  2    | 2
VLoose                              |  4    | 3
Loose                               |  8    | 4
Medium                              |  16   | 5
Tight                               |  32   | 6
VTight                              |  64   | 7
VVTight                             |  128  | 8
 * @param df The input dataframe
 * @param pt tau pt
 * @param decayMode decay mode of the tau
 * @param genMatch column with genmatch values (from prompt e, prompt mu,
 * tau->e, tau->mu, had. tau)
 * @param selectedDMs list of allowed decay modes for which a scale factor
 * should be calculated
 * @param wp working point of the ID cut
 * @param sf_vsjet_tau30to35 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau35to40 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau40to500 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau500to1000 id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau1000toinf id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_dependence "pt", "dm" or "eta" based scale factors
 * @param vsele_wp working point of the vsEle cut
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param idAlgorithm name of the tau id scale factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id_vsJet_lt(
    ROOT::RDF::RNode df, const std::string &pt, const std::string &decayMode,
    const std::string &genMatch, const std::vector<int> &selectedDMs,
    const std::string &wp, const std::string &sf_vsjet_tau30to35,
    const std::string &sf_vsjet_tau35to40,
    const std::string &sf_vsjet_tau40to500,
    const std::string &sf_vsjet_tau500to1000,
    const std::string &sf_vsjet_tau1000toinf, const std::string &sf_dependence,
    const std::string &vsele_wp, const std::string &id_output,
    const std::string &sf_file, const std::string &idAlgorithm) {

    Logger::get("TauIDvsJet_lt_SF")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("TauIDvsJet_lt_SF")->debug("ID - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto idSF_calculator = [evaluator, wp, vsele_wp, sf_vsjet_tau30to35,
                            sf_vsjet_tau35to40, sf_vsjet_tau40to500,
                            sf_vsjet_tau500to1000, sf_vsjet_tau1000toinf,
                            sf_dependence, selectedDMs,
                            idAlgorithm](const float &pt, const int &decayMode,
                                         const UChar_t &genMatch) {
        Logger::get("TauIDvsJet_lt_SF")->debug("ID - decayMode {}", decayMode);
        // only calculate SFs for allowed tau decay modes (also excludes default
        // values due to tau energy correction shifts below good tau pt
        // selection)
        double sf = 1.;
        if (std::find(selectedDMs.begin(), selectedDMs.end(), decayMode) !=
            selectedDMs.end()) {
            Logger::get("TauIDvsJet_lt_SF")
                ->debug("ID {} - pt {}, decayMode {}, genMatch {}, wp {}, "
                        "vsele_wp {}"
                        "sf_vsjet_tau30to35 {}, sf_vsjet_tau35to40 {}, "
                        "sf_vsjet_tau40to500{}, sf_vsjet_tau500to1000 {}, "
                        "sf_vsjet_tau1000toinf {}, sf_dependence {}",
                        idAlgorithm, pt, decayMode, genMatch, wp, vsele_wp,
                        sf_vsjet_tau30to35, sf_vsjet_tau35to40,
                        sf_vsjet_tau40to500, sf_vsjet_tau500to1000,
                        sf_vsjet_tau1000toinf, sf_dependence);
            if (pt >= 30.0 && pt < 35.0) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tau30to35, sf_dependence});
            } else if (pt >= 35.0 && pt < 40.0) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tau35to40, sf_dependence});
            } else if (pt >= 40.0 && pt < 500.0) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tau40to500, sf_dependence});
            } else if (pt >= 500.0 && pt < 1000.0) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tau500to1000, sf_dependence});
            } else if (pt >= 1000.0 && pt < 2000.0) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tau1000toinf, sf_dependence});
            } else {
                sf = 1.;
            }
        }
        Logger::get("TauIDvsJet_lt_SF")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {pt, decayMode, genMatch});
    return df1;
}
/**
 * @brief Function used to evaluate vsJets tau id scale factors in the lt
channel with the correctionlib for tauembedded samples


 * @param df The input dataframe
 * @param pt tau pt
 * @param wp working point of the ID cut
 * @param sf_vsjet_tau20to25 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau25to30 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau30to35 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau35to40 id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tau40toInf id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param correctionset name of the correction set containing the tau id scale
factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id_vsJet_lt_embedding(
    ROOT::RDF::RNode df, const std::string &pt, const std::string &wp,
    const std::string &sf_vsjet_tau20to25,
    const std::string &sf_vsjet_tau25to30,
    const std::string &sf_vsjet_tau30to35,
    const std::string &sf_vsjet_tau35to40,
    const std::string &sf_vsjet_tau40toInf, const std::string &id_output,
    const std::string &sf_file, const std::string &correctionset) {

    Logger::get("TauIDvsJet_lt_SF_embedding")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("TauIDvsJet_lt_SF_embedding")
        ->debug("ID - Name {}", correctionset);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(correctionset);
    auto idSF_calculator = [evaluator, wp, sf_vsjet_tau20to25,
                            sf_vsjet_tau25to30, sf_vsjet_tau30to35,
                            sf_vsjet_tau35to40, sf_vsjet_tau40toInf,
                            correctionset](const float &pt) {
        double sf = 1.;
        Logger::get("TauIDvsJet_lt_SF_embedding")
            ->debug("ID {} - pt {}, wp {} "
                    "sf_vsjet_tau20to25 {}, sf_vsjet_tau25to30 {}, "
                    "sf_vsjet_tau30to35{}, sf_vsjet_tau35to40 {}, "
                    "sf_vsjet_tau40toInf {},",
                    correctionset, pt, wp, sf_vsjet_tau20to25,
                    sf_vsjet_tau25to30, sf_vsjet_tau30to35, sf_vsjet_tau35to40,
                    sf_vsjet_tau40toInf);
        if (pt >= 20.0 && pt < 25.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau20to25, wp});
        } else if (pt >= 25.0 && pt < 30.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau25to30, wp});
        } else if (pt >= 30.0 && pt < 35.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau30to35, wp});
        } else if (pt >= 35.0 && pt < 40.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau35to40, wp});
        } else if (pt >= 40.0 && pt < 10000.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_tau40toInf, wp});
        } else {
            sf = 1.;
        }
        Logger::get("TauIDvsJet_lt_SF_embedding")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {pt});
    return df1;
}
/**
 * @brief Function used to evaluate vsJets tau id scale factors in the lt
channel with the correctionlib for tauembedded samples


 * @param df The input dataframe
 * @param pt tau pt
 * @param wp working point of the ID cut
 * @param vsele_wp working point of the vsEle cut
 * @param sf_vsjet_vsele_tau20to25 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_vsele_tau25to30 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_vsele_tau30to35 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_vsele_tau35to40 id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_vsele_tau40toInf id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param correctionset name of the correction set containing the tau id scale
factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id_vsJet_vsEle_lt_embedding(
    ROOT::RDF::RNode df, const std::string &pt, const std::string &wp,
    const std::string &vsele_wp,
    const std::string &sf_vsjet_vsele_tau20to25,
    const std::string &sf_vsjet_vsele_tau25to30,
    const std::string &sf_vsjet_vsele_tau30to35,
    const std::string &sf_vsjet_vsele_tau35to40,
    const std::string &sf_vsjet_vsele_tau40toInf, const std::string &id_output,
    const std::string &sf_file, const std::string &correctionset) {

    Logger::get("Tau_2_VsJet_VsEle_TauID_lt_SF_embedding")
        ->debug("Setting up function for tau id vsJet vsEle sf");
    Logger::get("Tau_2_VsJet_VsEle_TauID_lt_SF_embedding")
        ->debug("ID - Name {}", correctionset);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(correctionset);
    auto idSF_calculator = [evaluator, wp, vsele_wp, sf_vsjet_vsele_tau20to25,
                            sf_vsjet_vsele_tau25to30, sf_vsjet_vsele_tau30to35,
                            sf_vsjet_vsele_tau35to40, sf_vsjet_vsele_tau40toInf,
                            correctionset](const float &pt) {
        double sf = 1.;
        Logger::get("Tau_2_VsJet_VsEle_TauID_lt_SF_embedding")
            ->debug("ID {} - pt {}, wp {} "
                    "vsele_wp {}"
                    "sf_vsjet_vsele_tau20to25 {}, sf_vsjet_vsele_tau25to30 {}, "
                    "sf_vsjet_vsele_tau30to35{}, sf_vsjet_vsele_tau35to40 {}, "
                    "sf_vsjet_vsele_tau40toInf {},",
                    correctionset, pt, wp, vsele_wp, sf_vsjet_vsele_tau20to25,
                    sf_vsjet_vsele_tau25to30, sf_vsjet_vsele_tau30to35, sf_vsjet_vsele_tau35to40,
                    sf_vsjet_vsele_tau40toInf);
        if (pt >= 20.0 && pt < 25.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_vsele_tau20to25, wp, vsele_wp});
        } else if (pt >= 25.0 && pt < 30.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_vsele_tau25to30, wp, vsele_wp});
        } else if (pt >= 30.0 && pt < 35.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_vsele_tau30to35, wp, vsele_wp});
        } else if (pt >= 35.0 && pt < 40.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_vsele_tau35to40, wp, vsele_wp});
        } else if (pt >= 40.0 && pt < 10000.0) {
            sf = evaluator->evaluate({pt, sf_vsjet_vsele_tau40toInf, wp, vsele_wp});
        } else {
            sf = 1.;
        }
        Logger::get("Tau_2_VsJet_VsEle_TauID_lt_SF_embedding")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {pt});
    return df1;
}

/**
 * @brief Function used to evaluate vsJets tau id scale factors in the tt
channel with the correctionlib for tauembedded samples

 * @param df The input dataframe
 * @param decaymode decay mode of the tau
 * @param wp working point of the ID cut
 * @param vsele_wp working point of the vsEle cut
 * @param sf_vsjet_vsele_tauDM0 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_vsele_tauDM1 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_vsele_tauDM10 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_vsele_tauDM11 id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param correctionset name of the correction set containing the tau id scale
factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id_vsJet_vsEle_lt_dm_embedding(
    ROOT::RDF::RNode df, const std::string &decaymode, const std::string &wp,
    const std::string &vsele_wp,
    const std::string &sf_vsjet_vsele_tauDM0, const std::string &sf_vsjet_vsele_tauDM1,
    const std::string &sf_vsjet_vsele_tauDM10, const std::string &sf_vsjet_vsele_tauDM11,
    const std::string &id_output, const std::string &sf_file,
    const std::string &correctionset) {

    Logger::get("Tau_VsJet_VsEle_TauID_lt_SF_DM_embedding")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("Tau_VsJet_VsEle_TauID_lt_SF_DM_embedding")
        ->debug("ID - Name {}", correctionset);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(correctionset);
    auto idSF_calculator = [evaluator, wp, vsele_wp,
                            sf_vsjet_vsele_tauDM0, sf_vsjet_vsele_tauDM1,
                            sf_vsjet_vsele_tauDM10, sf_vsjet_vsele_tauDM11,
                            correctionset](const int &decaymode) {
        double sf = 1.;
        Logger::get("Tau_VsJet_VsEle_TauID_lt_SF_DM_embedding")
            ->debug("ID {} - decaymode {}, wp {}, "
                    "vsele_wp {},"
                    "sf_vsjet_vsele_tauDM0 {}, sf_vsjet_vsele_tauDM1 {}, "
                    "sf_vsjet_vsele_tauDM10{}, sf_vsjet_vsele_tauDM11 {}, ",
                    correctionset, decaymode, wp, vsele_wp, sf_vsjet_vsele_tauDM0,
                    sf_vsjet_vsele_tauDM1, sf_vsjet_vsele_tauDM10, sf_vsjet_vsele_tauDM11);
        if (decaymode == 0) {
            sf = evaluator->evaluate({decaymode, sf_vsjet_vsele_tauDM0, wp, vsele_wp});
        } else if (decaymode == 1) {
            sf = evaluator->evaluate({decaymode, sf_vsjet_vsele_tauDM1, wp, vsele_wp});
        } else if (decaymode == 10) {
            sf = evaluator->evaluate({decaymode, sf_vsjet_vsele_tauDM10, wp, vsele_wp});
        } else if (decaymode == 11) {
            sf = evaluator->evaluate({decaymode, sf_vsjet_vsele_tauDM11, wp, vsele_wp});
        } else {
            sf = 1.;
        }
        Logger::get("Tau_VsJet_VsEle_TauID_lt_SF_DM_embedding")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {decaymode});
    return df1;
}

/**
 * @brief Function used to evaluate vsJets tau id scale factors in the tt
channel with the correctionlib for tauembedded samples

 * @param df The input dataframe
 * @param decaymode decay mode of the tau
 * @param wp working point of the ID cut
 * @param sf_vsjet_tauDM0 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM1 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM10 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM11 id for the variation of the scale factor "sf"
for nominal
 * and "systup"/"systdown" the up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param correctionset name of the correction set containing the tau id scale
factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id_vsJet_tt_embedding(
    ROOT::RDF::RNode df, const std::string &decaymode, const std::string &wp,
    const std::string &sf_vsjet_tauDM0, const std::string &sf_vsjet_tauDM1,
    const std::string &sf_vsjet_tauDM10, const std::string &sf_vsjet_tauDM11,
    const std::string &id_output, const std::string &sf_file,
    const std::string &correctionset) {

    Logger::get("TauIDvsJet_tt_SF_embedding")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("TauIDvsJet_tt_SF_embedding")
        ->debug("ID - Name {}", correctionset);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(correctionset);
    auto idSF_calculator = [evaluator, wp, sf_vsjet_tauDM0, sf_vsjet_tauDM1,
                            sf_vsjet_tauDM10, sf_vsjet_tauDM11,
                            correctionset](const int &decaymode) {
        double sf = 1.;
        Logger::get("TauIDvsJet_tt_SF_embedding")
            ->debug("ID {} - decaymode {}, wp {} "
                    "sf_vsjet_tauDM0 {}, sf_vsjet_tauDM1 {}, "
                    "sf_vsjet_tauDM10{}, sf_vsjet_tauDM11 {}, ",
                    correctionset, decaymode, wp, sf_vsjet_tauDM0,
                    sf_vsjet_tauDM1, sf_vsjet_tauDM10, sf_vsjet_tauDM11);
        if (decaymode == 0) {
            sf = evaluator->evaluate({decaymode, sf_vsjet_tauDM0, wp});
        } else if (decaymode == 1) {
            sf = evaluator->evaluate({decaymode, sf_vsjet_tauDM1, wp});
        } else if (decaymode == 10) {
            sf = evaluator->evaluate({decaymode, sf_vsjet_tauDM10, wp});
        } else if (decaymode == 11) {
            sf = evaluator->evaluate({decaymode, sf_vsjet_tauDM11, wp});
        } else {
            sf = 1.;
        }
        Logger::get("TauIDvsJet_tt_SF_embedding")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {decaymode});
    return df1;
}
/**
 * @brief Function used to evaluate vsJets tau id scale factors in the tt
channel with
 * correctionlib

Description of the bit map used to define the tau id working points of the
DeepTau2017v2p1 tagger.
vsJets                              | Value | Bit (value used in the config)
------------------------------------|-------|-------
no ID selection (takes every tau)   |  0    | -
VVVLoose                            |  1    | 1
VVLoose                             |  2    | 2
VLoose                              |  4    | 3
Loose                               |  8    | 4
Medium                              |  16   | 5
Tight                               |  32   | 6
VTight                              |  64   | 7
VVTight                             |  128  | 8
 * @param df The input dataframe
 * @param pt tau pt
 * @param decayMode decay mode of the tau
 * @param genMatch column with genmatch values (from prompt e, prompt mu,
 * tau->e, tau->mu, had. tau)
 * @param selectedDMs list of allowed decay modes for which a scale factor
 * should be calculated
 * @param wp working point of the ID cut
 * @param sf_vsjet_tauDM0 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM1 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM10 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsjet_tauDM11 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_dependence "pt", "dm" or "eta" based scale factors
 * @param vsele_wp working point of the vsEle cut
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param idAlgorithm name of the tau id scale factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id_vsJet_tt(
    ROOT::RDF::RNode df, const std::string &pt, const std::string &decayMode,
    const std::string &genMatch, const std::vector<int> &selectedDMs,
    const std::string &wp, const std::string &sf_vsjet_tauDM0,
    const std::string &sf_vsjet_tauDM1, const std::string &sf_vsjet_tauDM10,
    const std::string &sf_vsjet_tauDM11, const std::string &sf_dependence,
    const std::string &vsele_wp, const std::string &id_output,
    const std::string &sf_file, const std::string &idAlgorithm) {

    Logger::get("TauIDvsJet_tt_SF")
        ->debug("Setting up function for tau id vsJet sf");
    Logger::get("TauIDvsJet_tt_SF")->debug("ID - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto idSF_calculator = [evaluator, wp, vsele_wp, sf_vsjet_tauDM0,
                            sf_vsjet_tauDM1, sf_vsjet_tauDM10, sf_vsjet_tauDM11,
                            sf_dependence, selectedDMs,
                            idAlgorithm](const float &pt, const int &decayMode,
                                         const UChar_t &genMatch) {
        Logger::get("TauIDvsJet_tt_SF")->debug("ID - decayMode {}", decayMode);
        // only calculate SFs for allowed tau decay modes (also excludes default
        // values due to tau energy correction shifts below good tau pt
        // selection)
        double sf = 1.;
        if (std::find(selectedDMs.begin(), selectedDMs.end(), decayMode) !=
            selectedDMs.end()) {
            Logger::get("TauIDvsJet_tt_SF")
                ->debug("ID {} - pt {}, decayMode {}, genMatch {}, wp {}, "
                        "vsele_wp {}, "
                        "sf_vsjet_tauDM0 {}, sf_vsjet_tauDM1 {}, "
                        "sf_vsjet_tauDM1 {}, sf_vsjet_tauDM10{}, "
                        "sf_vsjet_tauDM11 {}, sf_dependence {}",
                        idAlgorithm, pt, decayMode, genMatch, wp, vsele_wp,
                        sf_vsjet_tauDM0, sf_vsjet_tauDM1, sf_vsjet_tauDM10,
                        sf_vsjet_tauDM11, sf_dependence);
            if (decayMode == 0) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tauDM0, sf_dependence});
            } else if (decayMode == 1) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tauDM1, sf_dependence});
            } else if (decayMode == 10) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tauDM10, sf_dependence});
            } else if (decayMode == 11) {
                sf = evaluator->evaluate(
                    {pt, decayMode, static_cast<int>(genMatch), wp, vsele_wp,
                     sf_vsjet_tauDM11, sf_dependence});
            } else {
                sf = 1.;
            }
        }
        Logger::get("TauIDvsJet_tt_SF")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, idSF_calculator, {pt, decayMode, genMatch});
    return df1;
}
/**
 * @brief Function used to evaluate vsEle tau id scale factors with
 * correctionlib

Description of the bit map used to define the tau id working points of the
DeepTau2017v2p1 tagger.
vsElectrons                         | Value | Bit (value used in the config)
------------------------------------|-------|-------
no ID selection (takes every tau)   |  0    | -
VVVLoose                            |  1    | 1
VVLoose                             |  2    | 2
VLoose                              |  4    | 3
Loose                               |  8    | 4
Medium                              |  16   | 5
Tight                               |  32   | 6
VTight                              |  64   | 7
VVTight                             |  128  | 8

vsMuons                             | Value | Bit (value used in the config)
------------------------------------|-------|-------
no ID selection (takes every tau)   |  0    | -
VLoose                              |  1    | 1
Loose                               |  2    | 2
Medium                              |  4    | 3
Tight                               |  8    | 4
 * @param df The input dataframe
 * @param eta tau eta
 * @param decayMode decay mode of the tau
 * @param genMatch column with genmatch values (from prompt e, prompt mu,
 * tau->e, tau->mu, had. tau)
 * @param selectedDMs list of allowed decay modes for which a scale factor
 * should be calculated
 * @param wp working point of the ID cut
 * @param sf_vsele_barrel id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsele_endcap id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param idAlgorithm name of the tau id scale factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
id_vsEle(ROOT::RDF::RNode df, const std::string &eta,
         const std::string &decayMode, const std::string &genMatch,
         const std::vector<int> &selectedDMs, const std::string &wp,
         const std::string &sf_vsele_barrel, const std::string &sf_vsele_endcap,
         const std::string &id_output, const std::string &sf_file,
         const std::string &idAlgorithm) {

    Logger::get("TauIDvsEleSF")
        ->debug("Setting up function for tau id vsEle sf");
    Logger::get("TauIDvsEleSF")->debug("ID - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto idSF_calculator = [evaluator, wp, sf_vsele_barrel, sf_vsele_endcap,
                            selectedDMs,
                            idAlgorithm](const float &eta, const int &decayMode,
                                         const UChar_t &genMatch) {
        double sf = 1.;
        Logger::get("TauIDvsEleSF")->debug("ID - decayMode {}", decayMode);
        // only calculate SFs for allowed tau decay modes (also excludes
        // default values due to tau energy correction shifts below good tau
        // pt selection)
        if (std::find(selectedDMs.begin(), selectedDMs.end(), decayMode) !=
            selectedDMs.end()) {
            Logger::get("TauIDvsEleSF")
                ->debug("ID {} - eta {}, genMatch {}, wp {}, sf_vsele_barrel "
                        "{}, sf_vsele_endcap {}",
                        idAlgorithm, eta, genMatch, wp, sf_vsele_barrel,
                        sf_vsele_endcap);
            if (std::abs(eta) < 1.46) {
                sf = evaluator->evaluate(
                    {eta, static_cast<int>(genMatch), wp, sf_vsele_barrel});
            } else if (std::abs(eta) >= 1.558 && std::abs(eta) < 2.3) {
                sf = evaluator->evaluate(
                    {eta, static_cast<int>(genMatch), wp, sf_vsele_endcap});
            } else {
                sf = 1.;
            }
        }
        Logger::get("TauIDvsEleSF")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 =
        df.Define(id_output, idSF_calculator, {eta, decayMode, genMatch});
    return df1;
}
/**
 * @brief Function used to evaluate vsMu tau id scale factors with
 * correctionlib

Description of the bit map used to define the tau id working points of the
DeepTau2017v2p1 tagger.
vsElectrons                         | Value | Bit (value used in the config)
------------------------------------|-------|-------
no ID selection (takes every tau)   |  0    | -
VVVLoose                            |  1    | 1
VVLoose                             |  2    | 2
VLoose                              |  4    | 3
Loose                               |  8    | 4
Medium                              |  16   | 5
Tight                               |  32   | 6
VTight                              |  64   | 7
VVTight                             |  128  | 8

vsMuons                             | Value | Bit (value used in the config)
------------------------------------|-------|-------
no ID selection (takes every tau)   |  0    | -
VLoose                              |  1    | 1
Loose                               |  2    | 2
Medium                              |  4    | 3
Tight                               |  8    | 4
 * @param df The input dataframe
 * @param eta tau eta
 * @param decayMode decay mode of the tau
 * @param genMatch column with genmatch values (from prompt e, prompt mu,
 * tau->e, tau->mu, had. tau)
 * @param selectedDMs list of allowed decay modes for which a scale factor
 * should be calculated
 * @param wp working point of the ID cut
 * @param sf_vsmu_wheel1 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsmu_wheel2 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsmu_wheel3 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsmu_wheel4 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param sf_vsmu_wheel5 id for the variation of the scale factor "sf" for
nominal
 * and "systup"/"systdown" the up/down variation
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param idAlgorithm name of the tau id scale factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
id_vsMu(ROOT::RDF::RNode df, const std::string &eta,
        const std::string &decayMode, const std::string &genMatch,
        const std::vector<int> &selectedDMs, const std::string &wp,
        const std::string &sf_vsmu_wheel1, const std::string &sf_vsmu_wheel2,
        const std::string &sf_vsmu_wheel3, const std::string &sf_vsmu_wheel4,
        const std::string &sf_vsmu_wheel5, const std::string &id_output,
        const std::string &sf_file, const std::string &idAlgorithm) {

    Logger::get("TauIDvsMuSF")->debug("Setting up function for tau id vsMu sf");
    Logger::get("TauIDvsMuSF")->debug("ID - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto idSF_calculator = [evaluator, wp, sf_vsmu_wheel1, sf_vsmu_wheel2,
                            sf_vsmu_wheel3, sf_vsmu_wheel4, sf_vsmu_wheel5,
                            selectedDMs,
                            idAlgorithm](const float &eta, const int &decayMode,
                                         const UChar_t &genMatch) {
        double sf = 1.;
        Logger::get("TauIDvsMuSF")->debug("ID - decayMode {}", decayMode);
        // only calculate SFs for allowed tau decay modes (also excludes
        // default values due to tau energy correction shifts below good tau
        // pt selection)
        if (std::find(selectedDMs.begin(), selectedDMs.end(), decayMode) !=
            selectedDMs.end()) {
            Logger::get("TauIDvsMuSF")
                ->debug("ID {} - eta {}, genMatch {}, wp {}, sf_vsmu_wheel1 "
                        "{}, sf_vsmu_wheel2 {}, sf_vsmu_wheel3 {}, "
                        "sf_vsmu_wheel4 {}, sf_vsmu_wheel5 {}",
                        idAlgorithm, eta, genMatch, wp, sf_vsmu_wheel1,
                        sf_vsmu_wheel2, sf_vsmu_wheel3, sf_vsmu_wheel4,
                        sf_vsmu_wheel5);
            if (std::abs(eta) < 0.4) {
                sf = evaluator->evaluate(
                    {eta, static_cast<int>(genMatch), wp, sf_vsmu_wheel1});
            } else if (std::abs(eta) >= 0.4 && std::abs(eta) < 0.8) {
                sf = evaluator->evaluate(
                    {eta, static_cast<int>(genMatch), wp, sf_vsmu_wheel2});
            } else if (std::abs(eta) >= 0.8 && std::abs(eta) < 1.2) {
                sf = evaluator->evaluate(
                    {eta, static_cast<int>(genMatch), wp, sf_vsmu_wheel3});
            } else if (std::abs(eta) >= 1.2 && std::abs(eta) < 1.7) {
                sf = evaluator->evaluate(
                    {eta, static_cast<int>(genMatch), wp, sf_vsmu_wheel4});
            } else if (std::abs(eta) >= 1.7 && std::abs(eta) < 2.3) {
                sf = evaluator->evaluate(
                    {eta, static_cast<int>(genMatch), wp, sf_vsmu_wheel5});
            } else {
                sf = 1.0;
            }
        }
        Logger::get("TauIDvsMuSF")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 =
        df.Define(id_output, idSF_calculator, {eta, decayMode, genMatch});
    return df1;
}

/**
 * @brief Function to evaluate the tau trigger scale factor from a xpog file
 *
 * @param df the input dataframe
 * @param decaymode the name of the column containing the tau decay mode
 * variable
 * @param wp the name of the the tau id working point
 * @param type the type of the tau trigger scale factor, available are emb and
 * mc
 * @param pt the name of the column containing the tau pt variable
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the tau scale factors
 * @param correctionset name of the tau trigger scale factor
 * @return ROOT::RDF::RNode a new dataframe containing the new sf column
 */

ROOT::RDF::RNode
tau_trigger_sf(ROOT::RDF::RNode df, const std::string &decaymode,
               const std::string &pt, const std::string &wp,
               const std::string &type, const std::string &id_output,
               const std::string &sf_file, const std::string &correctionset) {

    Logger::get("tau_trigger_sf")
        ->debug("Setting up function for tau trigger sf");
    Logger::get("tau_trigger_sf")
        ->debug("ID - Name {}, file {}", correctionset, sf_file);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(correctionset);
    Logger::get("tau_trigger_sf")->debug("WP {} - type {}", wp, type);
    auto trigger_sf_calculator = [evaluator, wp, type, correctionset](
                                     const int &decaymode, const float &pt) {
        float sf = 1.;
        Logger::get("tau_trigger_sf")
            ->debug("ID {} - decaymode {}, wp {} "
                   "pt {}, type {}, ",
                   correctionset, decaymode, wp, pt, type);
        sf = evaluator->evaluate({decaymode, wp, type, pt});
        Logger::get("tau_trigger_sf")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(id_output, trigger_sf_calculator, {decaymode, pt});
    return df1;
}
} // namespace tau

namespace electron {
/**
 * @brief Function used to evaluate id scale factors of electrons with
 * correctionlib, configurations:
 * - [UL2018 Electron
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2018_UL.html)
 * - [UL2017 Electron
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2017_UL.html)
 * - [UL2016preVFP Electron
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2016preVFP_UL.html)
 * - [UL2016postVFP Electron
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/EGM_electron_Run2_UL/EGM_electron_2016postVFP_UL.html)
 * @param df The input dataframe
 * @param pt electron pt
 * @param eta electron eta
 * @param year_id id for the year of data taking and mc compaign
 * @param wp wp of the electron id
 * @param variation id for the variation of the scale factor. Available Values:
 * sf, sfdown, sfup
 * @param id_output name of the id scale factor column
 * @param sf_file path to the file with the electron scale factors
 * @param idAlgorithm name of the electron id scale factor
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &pt,
                    const std::string &eta, const std::string &year_id,
                    const std::string &wp, const std::string &variation,
                    const std::string &id_output, const std::string &sf_file,
                    const std::string &idAlgorithm) {

    Logger::get("electronIDSF")
        ->debug("Setting up functions for electron id sf with correctionlib");
    Logger::get("electronIDSF")->debug("ID - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        id_output,
        [evaluator, year_id, idAlgorithm, wp, variation](const float &pt,
                                                         const float &eta) {
            Logger::get("electronIDSF")
                ->debug("Year {}, Name {}, WP {}", year_id, idAlgorithm, wp);
            Logger::get("electronIDSF")->debug("ID - pt {}, eta {}", pt, eta);
            double sf = 1.;
            if (pt >= 0.0) {
                sf = evaluator->evaluate({year_id, variation, wp, eta, pt});
            }
            Logger::get("electronIDSF")->debug("Scale Factor {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}

} // namespace electron
namespace jet {
/**
 * @brief Function used to evaluate b-tagging scale factors of jets with
 * correctionlib, configurations:
 * - [UL2018 b-tagging
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/BTV_btagging_Run2_UL/BTV_btagging_2018_UL.html)
 * - [UL2017 b-tagging
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/BTV_btagging_Run2_UL/BTV_btagging_2017_UL.html)
 * - [UL2016preVFP b-tagging
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/BTV_btagging_Run2_UL/BTV_btagging_2016preVFP_UL.html)
 * - [UL2016postVFP b-tagging
 * ID](https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/BTV_btagging_Run2_UL/BTV_btagging_2016postVFP_UL.html)
 * @param df The input dataframe
 * @param pt jet pt
 * @param eta jet eta
 * @param btag_discr btag value of a jet based on a b-jet tagger (e.g. DeepJet)
 * @param flavor flavor of the jet
 * @param jet_mask mask for good/selected jets
 * @param bjet_mask mask for good/selected b jets
 * @param jet_veto_mask veto mask for overlapping jets
 * @param variation id for the variation of the scale factor. Available Values:
 * central, down_*, up_* (* name of variation)
 * @param sf_output name of the scale factor column
 * @param sf_file path to the file with the btagging scale factors
 * @param corr_algorithm name of the btagging correction algorithm
 * @return a new dataframe containing the new column
 */
ROOT::RDF::RNode
btagSF(ROOT::RDF::RNode df, const std::string &pt, const std::string &eta,
       const std::string &btag_discr, const std::string &flavor,
       const std::string &jet_mask, const std::string &bjet_mask,
       const std::string &jet_veto_mask, const std::string &variation,
       const std::string &sf_output, const std::string &sf_file,
       const std::string &corr_algorithm) {
    Logger::get("btagSF")->debug(
        "Setting up functions for b-tag sf with correctionlib");
    Logger::get("btagSF")->debug("Correction algorithm - Name {}",
                                 corr_algorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(corr_algorithm);

    auto btagSF_lambda = [evaluator,
                          variation](const ROOT::RVec<float> &pt_values,
                                     const ROOT::RVec<float> &eta_values,
                                     const ROOT::RVec<float> &btag_values,
                                     const ROOT::RVec<int> &flavors,
                                     const ROOT::RVec<int> &jet_mask,
                                     const ROOT::RVec<int> &bjet_mask,
                                     const ROOT::RVec<int> &jet_veto_mask) {
        Logger::get("btagSF")->debug("Vatiation - Name {}", variation);
        float sf = 1.;
        for (int i = 0; i < pt_values.size(); i++) {
            Logger::get("btagSF")->debug(
                "jet masks - jet {}, bjet {}, jet veto {}", jet_mask.at(i),
                bjet_mask.at(i), jet_veto_mask.at(i));
            // considering only good jets/bjets, this is needed since jets and
            // bjets might have different cuts depending on the analysis
            if ((jet_mask.at(i) || bjet_mask.at(i)) && jet_veto_mask.at(i)) {
                Logger::get("btagSF")->debug(
                    "SF - pt {}, eta {}, btag value {}, flavor {}",
                    pt_values.at(i), eta_values.at(i), btag_values.at(i),
                    flavors.at(i));
                float jet_sf = 1.;
                // considering only phase space where the scale factors are
                // defined
                if (pt_values.at(i) >= 20.0 && pt_values.at(i) < 10000.0 &&
                    std::abs(eta_values.at(i)) < 2.5) {
                    // for c jet related uncertainties only scale factors of
                    // c-jets are varied, the rest is nominal/central
                    if (variation.find("cferr") != std::string::npos) {
                        // flavor=4 means c-flavor
                        if (flavors.at(i) == 4) {
                            jet_sf = evaluator->evaluate(
                                {variation, flavors.at(i),
                                 std::abs(eta_values.at(i)), pt_values.at(i),
                                 btag_values.at(i)});
                        } else {
                            jet_sf = evaluator->evaluate(
                                {"central", flavors.at(i),
                                 std::abs(eta_values.at(i)), pt_values.at(i),
                                 btag_values.at(i)});
                        }
                    }
                    // for nominal/central and all other uncertainties c-jets
                    // have a scale factor of 1 (only for central defined in
                    // json file from BTV)
                    else {
                        if (flavors.at(i) != 4) {
                            jet_sf = evaluator->evaluate(
                                {variation, flavors.at(i),
                                 std::abs(eta_values.at(i)), pt_values.at(i),
                                 btag_values.at(i)});
                        } else {
                            jet_sf = evaluator->evaluate(
                                {"central", flavors.at(i),
                                 std::abs(eta_values.at(i)), pt_values.at(i),
                                 btag_values.at(i)});
                        }
                    }
                }
                Logger::get("btagSF")->debug("Jet Scale Factor {}", jet_sf);
                sf *= jet_sf;
            }
        };
        Logger::get("btagSF")->debug("Event Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(
        sf_output, btagSF_lambda,
        {pt, eta, btag_discr, flavor, jet_mask, bjet_mask, jet_veto_mask});
    return df1;
}
} // namespace jet

namespace embedding {
/**
 * @brief Function used to readout the embedding selection trigger scalefactors
 *
 * @param df the input dataframe
 * @param pt_1 the pt of the leading generator particle in the event. This
 * corresponds to the leading muon selected by the embedding selection
 * @param eta_1 the eta of the leading generator particle in the event. This
 * corresponds to the leading muon selected by the embedding selection
 * @param pt_2 the pt of the subleading generator particle in the event. This
 * corresponds to the subleading muon selected by the embedding selection
 * @param eta_2 the eta of the subleading generator particle in the event. This
 * corresponds to the subleading muon selected by the embedding selection
 * @param output name of the output column
 * @param sf_file path to the correctionlib file containing the scale factor
 * @param idAlgorithm name of the scale factor in the correctionlib file
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode
selection_trigger(ROOT::RDF::RNode df, const std::string &pt_1,
                  const std::string &eta_1, const std::string &pt_2,
                  const std::string &eta_2, const std::string &output,
                  const std::string &sf_file, const std::string &idAlgorithm) {

    Logger::get("EmbeddingSelectionTriggerSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator](const float &pt_1, const float &eta_1, const float &pt_2,
                    const float &eta_2) {
            Logger::get("EmbeddingSelectionTriggerSF")
                ->debug(" pt_1 {}, eta_1 {}, pt_2 {}, eta_2 {}", pt_1, eta_1,
                        pt_2, eta_2);
            double sf = 1.;
            sf = evaluator->evaluate(
                {pt_1, std::abs(eta_1), pt_2, std::abs(eta_2)});
            Logger::get("EmbeddingSelectionTriggerSF")->debug("sf {}", sf);
            return sf;
        },
        {pt_1, eta_1, pt_2, eta_2});
    return df1;
}
/**
 * @brief Function used to readout the embedding selection trigger scalefactors.
 *
 * @param df the input dataframe
 * @param pt the pt of the generator particle in the event. This corresponds to
 * one of the muons selected by the embedding selection
 * @param eta the eta of the generator particle in the event. This corresponds
 * to one of the muons selected by the embedding selection
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode selection_id(ROOT::RDF::RNode df, const std::string &pt,
                              const std::string &eta, const std::string &output,
                              const std::string &sf_file,
                              const std::string &idAlgorithm) {

    Logger::get("EmbeddingSelectionIDSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 =
        df.Define(output,
                  [evaluator](const float &pt, const float &eta) {
                      Logger::get("EmbeddingSelectionIDSF")
                          ->debug(" pt {}, eta {},", pt, eta);
                      double sf = 1.;
                      sf = evaluator->evaluate({pt, std::abs(eta)});
                      Logger::get("EmbeddingSelectionIDSF")->debug("sf {}", sf);
                      return sf;
                  },
                  {pt, eta});
    return df1;
}
/**
 * @brief Function used to readout SFs from the muon scale factor measurements
 *
 * @param df the input dataframe
 * @param pt the pt of the muon
 * @param eta the eta of the muon
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param correctiontype the type of the correction. Use `emb` for embedding and
 * `mc` for monte carlo
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @param extrapolation_factor The extrapolation factor to be used for the scale
 * factor, defaults to 1.
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode muon_sf(ROOT::RDF::RNode df, const std::string &pt,
                         const std::string &eta, const std::string &output,
                         const std::string &sf_file,
                         const std::string correctiontype,
                         const std::string &idAlgorithm,
                         const float &extrapolation_factor = 1.0) {

    Logger::get("EmbeddingMuonSF")->debug("Correction - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator, correctiontype, extrapolation_factor](const float &pt,
                                                          const float &eta) {
            Logger::get("EmbeddingMuonSF")
                ->debug(" pt {}, eta {}, correctiontype {}, extrapolation "
                        "factor {}",
                        pt, eta, correctiontype, extrapolation_factor);
            double sf = 1.;
            sf = extrapolation_factor *
                 evaluator->evaluate({pt, std::abs(eta), correctiontype});
            Logger::get("EmbeddingMuonSF")->debug("sf {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
/**
 * @brief Function used to readout SFs from the electron scale factor
 * measurements
 *
 * @param df the input dataframe
 * @param pt the pt of the electron
 * @param eta the eta of the electron
 * @param output the name of the output column
 * @param sf_file the path to the correctionlib file containing the scale factor
 * @param correctiontype the type of the correction. Use `emb` for embedding and
 * `mc` for monte carlo
 * @param idAlgorithm the name of the scale factor in the correctionlib
 * file
 * @param extrapolation_factor The extrapolation factor to be used for the scale
 * factor, defaults to 1.
 * @return ROOT::RDF::RNode
 */
ROOT::RDF::RNode electron_sf(ROOT::RDF::RNode df, const std::string &pt,
                             const std::string &eta, const std::string &output,
                             const std::string &sf_file,
                             const std::string correctiontype,
                             const std::string &idAlgorithm,
                             const float &extrapolation_factor = 1.0) {

    Logger::get("EmbeddingElectronSF")
        ->debug("Correction - Name {}", idAlgorithm);
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(idAlgorithm);
    auto df1 = df.Define(
        output,
        [evaluator, correctiontype, extrapolation_factor](const float &pt,
                                                          const float &eta) {
            Logger::get("EmbeddingElectronSF")
                ->debug(" pt {}, eta {}, correctiontype {}, extrapolation "
                        "factor {}",
                        pt, eta, correctiontype, extrapolation_factor);
            double sf = 1.;
            sf = extrapolation_factor *
                 evaluator->evaluate({pt, eta, correctiontype});
            Logger::get("EmbeddingElectronSF")->debug("sf {}", sf);
            return sf;
        },
        {pt, eta});
    return df1;
}
/**
 * @brief Function to evaluate the di-tau trigger or etau/mutau cross trigger
 * scale factor for embedded events from a xpog file
 *
 * @param df the input dataframe
 * @param pt the name of the column containing the tau pt variable
 * @param decaymode the name of the column containing the tau decay mode
 * variable
 * @param output name of the scale factor column
 * @param wp the name of the the tau id working point VVVLoose-VVTight
 * @param sf_file path to the file with the tau trigger scale factors
 * @param type the type of the tau trigger, available are "ditau", "etau",
 * "mutau", "ditauvbf"
 * @param corrtype name of the tau trigger correction type, available are
 * "eff_data", "eff_mc", "sf"
 * @param syst name of the systematic variation, options are "nom", "up", "down"
 * @return ROOT::RDF::RNode a new dataframe containing the new sf column
 */

ROOT::RDF::RNode
ditau_trigger_sf(ROOT::RDF::RNode df, const std::string &pt,
                 const std::string &decaymode, const std::string &output,
                 const std::string &wp, const std::string &sf_file,
                 const std::string &type, const std::string &corrtype,
                 const std::string &syst) {

    Logger::get("ditau_trigger")
        ->debug("Setting up function for di-tau trigger sf");
    Logger::get("ditau_trigger")
        ->debug("correction type {}, file {}", corrtype, sf_file);
    // tauTriggerSF is the only correction set in the file for now, might change
    // with official sf release -> change into additional input parameter
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at("tauTriggerSF");
    Logger::get("ditau_trigger")->debug("WP {} - trigger type {}", wp, type);
    auto trigger_sf_calculator = [evaluator, wp, type, corrtype,
                                  syst](const float &pt, const int &decaymode) {
        float sf = 1.;
        Logger::get("ditau_trigger")
            ->debug("decaymode {}, pt {}", decaymode, pt);
        if (pt > 0) {
            if (decaymode == 0 || decaymode == 1 || decaymode == 10 ||
                decaymode == 11) {
                sf = evaluator->evaluate(
                    {pt, decaymode, type, wp, corrtype, syst});
            } else {
                sf = evaluator->evaluate({pt, -1, type, wp, corrtype, syst});
            }
        }
        Logger::get("ditau_trigger")->debug("Scale Factor {}", sf);
        return sf;
    };
    auto df1 = df.Define(output, trigger_sf_calculator, {pt, decaymode});
    return df1;
}
} // namespace embedding
} // namespace scalefactor
#endif /* GUARD_SCALEFACTORS_H */