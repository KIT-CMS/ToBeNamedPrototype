#ifndef GUARD_TRIGGERS_H
#define GUARD_TRIGGERS_H

#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "bitset"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
#include <fstream>
#include <nlohmann/json.hpp>
#include <regex>

typedef std::bitset<30> IntBits;

namespace trigger {

/**
 * @brief Function used to try and match a object with a trigger object. An
object is successfully matched, if they overlap within the given deltaR cone
(matchDeltaR), the pt and eta thresholds are met, and the triggerobject id and
the triggerobject bit are set to the correct value.

For the trigger objects, two additional quantities are stored, the id of the
particle and some further information on the trigger path encoded in a bitmap.
The triggerobject ID is encoded as follows:
    @code
    1 = Jet
    2 = MET
    3 = HT
    4 = MHT
    6 = FatJet
    11 = Electron (PixelMatched e/gamma)
    22 = Photon (PixelMatch-vetoed e/gamma)
    13 = Muon
    15 = Tau
    @endcode

The additional trigger information is encoded as follows. A more detailed
description can be found in the corresponding NanoAOD producer
[link](https://github.com/cms-sw/cmssw/blob/master/PhysicsTools/NanoAOD/python/triggerObjects_cff.py#L17)


Electrons                           | Value | Bit (value used in the config)
------------------------------------|-------|-------
ignore matching                     |  -    | -1
no match                            |  0    | -
CaloIdL_TrackIdL_IsoVL              |  1    | 0
1e (WPTight)                        |  2    | 1
1e (WPLoose)                        |  4    | 2
OverlapFilter PFTau                 |  8    | 3
2e                                  |  16   | 4
1e-1mu                              |  32   | 5
1e-1tau                             |  64   | 6
3e                                  |  128  | 7
2e-1mu                              |  256  | 8
1e-2mu                              |  512  | 9
1e (32_L1DoubleEG_AND_L1SingleEGOr) |  1024 | 10

Muons               | Value | Bit (value used in the config)
--------------------|-------|-------
ignore matching     |  -    | -1
no match            |  0    | -
TrkIsoVVL           |  1    | 0
Iso                 |  2    | 1
OverlapFilter PFTau |  4    | 2
1mu                 |  8    | 3
2mu                 |  16   | 4
1mu-1e              |  32   | 5
1mu-1tau            |  64   | 6
3mu                 |  128  | 7
2mu-1e              |  256  | 8
1mu-2e              |  512  | 9

Taus                 | Value | Bit (value used in the config)
---------------------|-------|-------
ignore matching      |  -    | -1
no match             |  0    | -
LooseChargedIso      |  1    | 0
MediumChargedIso     |  2    | 1
TightChargedIso      |  4    | 2
TightID OOSC photons |  8    | 3
HPS                  |  16   | 4
single-tau + tau+MET |  32   | 5
di-tau               |  64   | 6
e-tau                |  128  | 7
mu-tau               |  256  | 8
VBF+di-tau for Tau   |  512  | 9

jets                                    | Value | Bit (value used in the config)
----------------------------------------|-------|-------
no match                                |  0    | -
VBF cross-cleaned from loose iso PFTau  |  1    | 0
 * @param particle the `ROOT::Math::PtEtaPhiMVector` vector of the object to
match
 * @param triggerobject_pts `ROOT::RVec<float>` of trigger object pts
 * @param triggerobject_etas `ROOT::RVec<float>` of trigger object etas
 * @param triggerobject_phis `ROOT::RVec<float>` of trigger object phis
 * @param triggerobject_bits `ROOT::RVec<float>` of trigger object bits.
Depending on the trigger object id, this bitmap has a different meaning as
listed in th table above.
 * @param triggerobject_ids `ROOT::RVec<float>` of trigger object ids. The ID
encoding is explained above.
 * @param matchDeltaR The maximum deltaR value used for the match
 * @param pt_cut pt cut value on the trigger object pt
 * @param eta_cut eta cut value on the trigger object eta
 * @param trigger_particle_id_cut the triggerobject id required
 * @param triggerbit_cut the triggerobject filter bit position required
 * @return true, if all criteria are met, false otherwise
 */

bool matchParticle(const ROOT::Math::PtEtaPhiMVector &particle,
                   ROOT::RVec<float> &triggerobject_pts,
                   ROOT::RVec<float> &triggerobject_etas,
                   ROOT::RVec<float> &triggerobject_phis,
                   ROOT::RVec<int> &triggerobject_bits,
                   ROOT::RVec<int> &triggerobject_ids, const float &matchDeltaR,
                   const float &pt_cut, const float &eta_cut,
                   const int &trigger_particle_id_cut,
                   const int &triggerbit_cut) {
    Logger::get("CheckTriggerMatch")->debug("Checking Triggerobjects");
    Logger::get("CheckTriggerMatch")
        ->debug("Total number of triggerobjects: {}", triggerobject_pts.size());
    for (std::size_t idx = 0; idx < triggerobject_pts.size(); ++idx) {
        Logger::get("CheckTriggerMatch")->debug("Triggerobject Nr. {}", idx);
        Logger::get("CheckTriggerMatch")
            ->debug("bit Value: {}", IntBits(triggerobject_bits[idx]));
        Logger::get("CheckTriggerMatch")
            ->debug("bit Value: {}", triggerobject_bits[idx]);
        auto triggerobject = ROOT::Math::RhoEtaPhiVectorF(
            0, triggerobject_etas[idx], triggerobject_phis[idx]);
        // We check the deltaR match as well as that the pt and eta of the
        // triggerobject are above the given thresholds
        bool deltaR = ROOT::Math::VectorUtil::DeltaR(triggerobject, particle) <
                      matchDeltaR;
        // if we don't want to do any matching here, the triggerbut_cut value is
        // -1
        Logger::get("CheckTriggerMatch")
            ->debug("bit Value: {}", triggerobject_bits[idx]);
        bool bit = (triggerbit_cut == -1) ||
                   (IntBits(triggerobject_bits[idx]).test(triggerbit_cut));
        bool id = triggerobject_ids[idx] == trigger_particle_id_cut;
        bool pt = particle.pt() > pt_cut;
        bool eta = abs(particle.eta()) < eta_cut;
        Logger::get("CheckTriggerMatch")
            ->debug("Partice Lorentz Vector: {}, {}, {}, {}", particle.pt(),
                    particle.eta(), particle.phi(), particle.mass());
        Logger::get("CheckTriggerMatch")
            ->debug("-------------------------------------------------------");
        Logger::get("CheckTriggerMatch")
            ->debug("deltaR/matchDeltaR Check: {}/{}", deltaR, matchDeltaR);
        Logger::get("CheckTriggerMatch")
            ->debug("deltaR Value: {}",
                    ROOT::Math::VectorUtil::DeltaR(triggerobject, particle));
        Logger::get("CheckTriggerMatch")
            ->debug("id/trigger_particle_id_cut Check: {}/{}", id,
                    trigger_particle_id_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("id Value: {}", triggerobject_ids[idx]);
        Logger::get("CheckTriggerMatch")
            ->debug("bit/triggerbit_cut Check: {}/{}", bit, triggerbit_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("bit Value: {}", IntBits(triggerobject_bits[idx]));
        Logger::get("CheckTriggerMatch")
            ->debug("pt/pt_cut Check: {}/{}", pt, pt_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("pt Value (trg): {}, pt Value (reco): {}",
                    triggerobject_pts[idx], particle.pt());
        Logger::get("CheckTriggerMatch")
            ->debug("eta/eta_cut Check: {}/{}", eta, eta_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("eta (trg) Value: {}, eta (reco) Value: {}",
                    triggerobject_etas[idx], abs(particle.eta()));
        Logger::get("CheckTriggerMatch")
            ->debug("-------------------------------------------------------");
        if (deltaR && bit && id && pt && eta) {
            // remove the matching object from the object vectors so it cant be
            // matched by the next particle as well (if there is one)
            triggerobject_ids.erase(triggerobject_ids.begin() + idx);
            triggerobject_bits.erase(triggerobject_bits.begin() + idx);
            triggerobject_pts.erase(triggerobject_pts.begin() + idx);
            triggerobject_etas.erase(triggerobject_etas.begin() + idx);
            triggerobject_phis.erase(triggerobject_phis.begin() + idx);
            return true;
        }
    }
    return false;
};
/**
 * @brief Function to generate a trigger flag based on an hlt path and trigger
 * object matching for the given object. This relies on the
 * trigger::matchParticle function which does the matching test.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle_p4 `ROOT::Math::PtEtaPhiMVector` of the object to be checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param hltpath name of the hlt path to be checked, this can be a valid regex.
 * If more than one matching HLT path is found, the function will throw an
 * exception, if no matching HLT path is found, the function will return a
 * dataframe with a flag of false for all entries
 * @param pt_cut minimal pt value for the triggerobject
 * @param eta_cut maximal pt value for the triggerobject
 * @param trigger_particle_id_cut trigger id value the triggerobject has to
 * match (details can be found in the documentation of trigger::matchParticle)
 * @param triggerbit_cut trigger bit value the triggerobject has to match
 * (details can be found in the documentation of trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @return a new dataframe containing the trigger flag column
 */

ROOT::RDF::RNode GenerateSingleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &hltpath, const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold) {

    auto triggermatch = [DeltaR_threshold, pt_cut, eta_cut,
                         trigger_particle_id_cut, triggerbit_cut, hltpath](
                            bool hltpath_match,
                            const ROOT::Math::PtEtaPhiMVector &particle_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
        Logger::get("GenerateSingleTriggerFlag")->debug("Checking Trigger");
        Logger::get("CheckTriggerMatch")
            ->debug("Selected trigger: {}", hltpath);
        bool result = false;
        bool match_result = false;
        if (hltpath_match) {
            Logger::get("CheckTriggerMatch")
                ->debug("Checking Triggerobject match with particles ....");
            match_result = matchParticle(
                particle_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, pt_cut, eta_cut, trigger_particle_id_cut,
                triggerbit_cut);
        }
        result = hltpath_match & match_result;
        Logger::get("GenerateSingleTriggerFlag")
            ->debug("---> HLT Match: {}", hltpath_match);
        Logger::get("GenerateSingleTriggerFlag")
            ->debug("---> Total Match: {}", match_result);
        Logger::get("GenerateSingleTriggerFlag")
            ->debug("--->>>> result: {}", result);
        return result;
    };
    auto available_trigger = df.GetColumnNames();
    std::vector<std::string> matched_trigger_names;
    std::regex hltpath_regex = std::regex(hltpath);
    // loop over all available trigger names and check if the hltpath is
    // matching any of them
    for (auto &trigger : available_trigger) {
        if (std::regex_match(trigger, hltpath_regex)) {
            Logger::get("GenerateSingleTriggerFlag")
                ->debug("Found matching trigger: {}", trigger);
            matched_trigger_names.push_back(trigger);
        }
    }
    // if no matching trigger was found return the initial dataframe
    if (matched_trigger_names.size() == 0) {
        Logger::get("GenerateSingleTriggerFlag")
            ->info("No matching trigger for {} found, returning false for "
                   "trigger flag {}",
                   hltpath, triggerflag_name);
        auto df1 = df.Define(triggerflag_name, []() { return false; });
        return df1;
    } else if (matched_trigger_names.size() > 1) {
        Logger::get("GenerateSingleTriggerFlag")
            ->debug(
                "More than one matching trigger found, not implemented yet");
        throw std::invalid_argument(
            "received too many matching trigger paths, not implemented yet");
    } else {
        Logger::get("GenerateSingleTriggerFlag")
            ->debug("Found matching trigger: {}", matched_trigger_names[0]);
        auto df1 =
            df.Define(triggerflag_name, triggermatch,
                      {matched_trigger_names[0], particle_p4,
                       triggerobject_bits, triggerobject_id, triggerobject_pt,
                       triggerobject_eta, triggerobject_phi});
        return df1;
    }
}

/**
 * @brief Function to generate a trigger flag based on an hlt path and trigger
 * object matching for the given object. This relies on the
 * trigger::matchParticle function which does the matching test. The
 * implementation is similar to the trigger::GenerateSingleTriggerFlag function,
 * but here two objects have to be matched for a passing flag.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle1_p4 `ROOT::Math::PtEtaPhiMVector` of the first object to be
 * checked
 * @param particle2_p4 `ROOT::Math::PtEtaPhiMVector` of the second object to be
 * checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param hltpath name of the hlt path to be checked, this can be a valid regex.
 * If more than one matching HLT path is found, the function will throw an
 * exception, if no matching HLT path is found, the function will return a
 * dataframe with a flag of false for all entries
 * @param p1_pt_cut minimal pt value for the triggerobject matching the first
 * object
 * @param p2_pt_cut minimal pt value for the triggerobject matching the second
 * object
 * @param p1_eta_cut maximal pt value for the triggerobject matching the first
 * object
 * @param p2_eta_cut maximal pt value for the triggerobject matching the second
 * object
 * @param p1_trigger_particle_id_cut trigger id value the triggerobject matching
 * the first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_trigger_particle_id_cut trigger id value the triggerobject matching
 * the second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p1_triggerbit_cut trigger bit value the triggerobject matching the
 * first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_triggerbit_cut trigger bit value the triggerobject matching the
 * second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @return a new dataframe containing the trigger flag column
 *
 */

ROOT::RDF::RNode GenerateDoubleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle1_p4, const std::string &particle2_p4,
    const std::string &triggerobject_bits, const std::string &triggerobject_id,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &hltpath,
    const float &p1_pt_cut, const float &p2_pt_cut, const float &p1_eta_cut,
    const float &p2_eta_cut, const int &p1_trigger_particle_id_cut,
    const int &p2_trigger_particle_id_cut, const int &p1_triggerbit_cut,
    const int &p2_triggerbit_cut, const float &DeltaR_threshold) {

    auto triggermatch = [DeltaR_threshold, p1_pt_cut, p2_pt_cut, p1_eta_cut,
                         p2_eta_cut, p1_trigger_particle_id_cut,
                         p2_trigger_particle_id_cut, p1_triggerbit_cut,
                         p2_triggerbit_cut, hltpath](
                            bool hltpath_match,
                            const ROOT::Math::PtEtaPhiMVector &particle1_p4,
                            const ROOT::Math::PtEtaPhiMVector &particle2_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
        Logger::get("GenerateDoubleTriggerFlag")->debug("Checking Trigger");
        Logger::get("CheckTriggerMatch")
            ->debug("Selected trigger: {}", hltpath);
        bool result = false;
        bool match_result_p1 = false;
        bool match_result_p2 = false;
        if (hltpath_match) {
            Logger::get("GenerateDoubleTriggerFlag")
                ->debug("Checking Triggerobject match with particles ....");
            Logger::get("GenerateDoubleTriggerFlag")->debug("First particle");
            match_result_p1 = matchParticle(
                particle1_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, p1_pt_cut, p1_eta_cut,
                p1_trigger_particle_id_cut, p1_triggerbit_cut);
            Logger::get("GenerateDoubleTriggerFlag")->debug("Second particle");
            match_result_p2 = matchParticle(
                particle2_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, p2_pt_cut, p2_eta_cut,
                p2_trigger_particle_id_cut, p2_triggerbit_cut);
        }
        result = hltpath_match & match_result_p1 & match_result_p2;
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("---> HLT Match: {}", hltpath_match);
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("---> Total Match P1: {}", match_result_p1);
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("---> Total Match P2: {}", match_result_p2);
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("--->>>> result: {}", result);
        return result;
    };
    auto available_trigger = df.GetColumnNames();
    std::vector<std::string> matched_trigger_names;
    std::regex hltpath_regex = std::regex(hltpath);
    // loop over all available trigger names and check if the hltpath is
    // matching any of them
    for (auto &trigger : available_trigger) {
        if (std::regex_match(trigger, hltpath_regex)) {
            Logger::get("GenerateDoubleTriggerFlag")
                ->debug("Found matching trigger: {}", trigger);
            matched_trigger_names.push_back(trigger);
        }
    }
    // if no matching trigger was found return the initial dataframe
    if (matched_trigger_names.size() == 0) {
        Logger::get("GenerateDoubleTriggerFlag")
            ->info("No matching trigger for {} found, returning false for "
                   "trigger flag {}",
                   hltpath, triggerflag_name);
        auto df1 = df.Define(triggerflag_name, []() { return false; });
        return df1;
    } else if (matched_trigger_names.size() > 1) {
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug(
                "More than one matching trigger found, not implemented yet");
        throw std::invalid_argument(
            "received too many matching trigger paths, not implemented yet");
    } else {
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("Found matching trigger: {}", matched_trigger_names[0]);
        auto df1 =
            df.Define(triggerflag_name, triggermatch,
                      {matched_trigger_names[0], particle1_p4, particle2_p4,
                       triggerobject_bits, triggerobject_id, triggerobject_pt,
                       triggerobject_eta, triggerobject_phi});
        return df1;
    }
}

/**
 * @brief Function to generate a trigger flag based on a trigger
 * object matching for two given objects. This relies on the
 * trigger::matchParticle function which does the matching test. The
 * implementation is similar to the trigger::GenerateDoubleTriggerFlag function,
 * but here, no hlt path is required, only the matching of trigger objects is
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle1_p4 `ROOT::Math::PtEtaPhiMVector` of the first object to be
 * checked
 * @param particle2_p4 `ROOT::Math::PtEtaPhiMVector` of the second object to be
 * checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param p1_pt_cut minimal pt value for the triggerobject matching the first
 * object
 * @param p2_pt_cut minimal pt value for the triggerobject matching the second
 * object
 * @param p1_eta_cut maximal pt value for the triggerobject matching the first
 * object
 * @param p2_eta_cut maximal pt value for the triggerobject matching the second
 * object
 * @param p1_trigger_particle_id_cut trigger id value the triggerobject matching
 * the first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_trigger_particle_id_cut trigger id value the triggerobject matching
 * the second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p1_triggerbit_cut trigger bit value the triggerobject matching the
 * first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_triggerbit_cut trigger bit value the triggerobject matching the
 * second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @return a new dataframe containing the trigger flag column
 *
 */

ROOT::RDF::RNode MatchDoubleTriggerObject(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle1_p4, const std::string &particle2_p4,
    const std::string &triggerobject_bits, const std::string &triggerobject_id,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const float &p1_pt_cut,
    const float &p2_pt_cut, const float &p1_eta_cut, const float &p2_eta_cut,
    const int &p1_trigger_particle_id_cut,
    const int &p2_trigger_particle_id_cut, const int &p1_triggerbit_cut,
    const int &p2_triggerbit_cut, const float &DeltaR_threshold) {

    auto triggermatch =
        [DeltaR_threshold, p1_pt_cut, p2_pt_cut, p1_eta_cut, p2_eta_cut,
         p1_trigger_particle_id_cut, p2_trigger_particle_id_cut,
         p1_triggerbit_cut,
         p2_triggerbit_cut](const ROOT::Math::PtEtaPhiMVector &particle1_p4,
                            const ROOT::Math::PtEtaPhiMVector &particle2_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
            bool match_result_p1 = false;
            bool match_result_p2 = false;
            Logger::get("MatchDoubleTriggerObject")
                ->debug("Checking Triggerobject match with particles ....");
            Logger::get("MatchDoubleTriggerObject")->debug("First particle");
            match_result_p1 = matchParticle(
                particle1_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, p1_pt_cut, p1_eta_cut,
                p1_trigger_particle_id_cut, p1_triggerbit_cut);
            Logger::get("MatchDoubleTriggerObject")->debug("Second particle");
            match_result_p2 = matchParticle(
                particle2_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, p2_pt_cut, p2_eta_cut,
                p2_trigger_particle_id_cut, p2_triggerbit_cut);
            bool result = match_result_p1 & match_result_p2;
            Logger::get("MatchDoubleTriggerObject")
                ->debug("---> Total Match P1: {}", match_result_p1);
            Logger::get("MatchDoubleTriggerObject")
                ->debug("---> Total Match P2: {}", match_result_p2);
            Logger::get("MatchDoubleTriggerObject")
                ->debug("--->>>> result: {}", result);
            return result;
        };
    auto df1 = df.Define(triggerflag_name, triggermatch,
                         {particle1_p4, particle2_p4, triggerobject_bits,
                          triggerobject_id, triggerobject_pt, triggerobject_eta,
                          triggerobject_phi});
    return df1;
}

/**
 * @brief Function to generate a trigger flag based on a trigger
 * object match only. This relies on the
 * trigger::matchParticle function which does the matching test.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle_p4 `ROOT::Math::PtEtaPhiMVector` of the object to be checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param pt_cut minimal pt value for the triggerobject
 * @param eta_cut maximal pt value for the triggerobject
 * @param trigger_particle_id_cut trigger id value the triggerobject has to
 * match (details can be found in the documentation of trigger::matchParticle)
 * @param triggerbit_cut trigger bit value the triggerobject has to match
 * (details can be found in the documentation of trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @return a new dataframe containing the trigger object match flag column
 */

ROOT::RDF::RNode MatchSingleTriggerObject(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold) {

    auto triggermatch = [DeltaR_threshold, pt_cut, eta_cut,
                         trigger_particle_id_cut, triggerbit_cut](
                            const ROOT::Math::PtEtaPhiMVector &particle_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
        Logger::get("MatchSingleTriggerObject")->debug("Checking Trigger");
        Logger::get("MatchSingleTriggerObject")
            ->debug("Checking Triggerobject match with particles ....");
        bool match_result =
            matchParticle(particle_p4, triggerobject_pts, triggerobject_etas,
                          triggerobject_phis, triggerobject_bits,
                          triggerobject_ids, DeltaR_threshold, pt_cut, eta_cut,
                          trigger_particle_id_cut, triggerbit_cut);
        Logger::get("MatchSingleTriggerObject")
            ->debug("--->>>> match_result: {}", match_result);
        return match_result;
    };
    auto df1 =
        df.Define(triggerflag_name, triggermatch,
                  {particle_p4, triggerobject_bits, triggerobject_id,
                   triggerobject_pt, triggerobject_eta, triggerobject_phi});
    return df1;
}

namespace tagandprobe {
/**
 * @brief Function used to try and match a object with a trigger object. An
object is successfully matched, if they overlap within the given deltaR cone
(matchDeltaR), the pt and eta thresholds are met, and the triggerobject id and
the triggerobject bit are set to the correct value. Additionally, a threshold
on the pT of the trigger object is enforced.

For the trigger objects, two additional quantities are stored, the id of the
particle and some further information on the trigger path encoded in a bitmap.
More information on the triggerobject IDs and the information encoded in the
bitmap is given in the description of the same function in the trigger
namespace.
 * @param particle the `ROOT::Math::PtEtaPhiMVector` vector of the object to
match
 * @param triggerobject_pts `ROOT::RVec<float>` of trigger object pts
 * @param triggerobject_etas `ROOT::RVec<float>` of trigger object etas
 * @param triggerobject_phis `ROOT::RVec<float>` of trigger object phis
 * @param triggerobject_bits `ROOT::RVec<float>` of trigger object bits.
Depending on the trigger object id, this bitmap has a different meaning as
listed in th table above.
 * @param triggerobject_ids `ROOT::RVec<float>` of trigger object ids. The ID
encoding is explained above.
 * @param matchDeltaR The maximum deltaR value used for the match
 * @param pt_cut pt cut value on the selected particle pt
 * @param eta_cut eta cut value on the selected particle eta
 * @param trigger_particle_id_cut the triggerobject id required
 * @param triggerbit_cut the triggerobject filter bit position required
 * @param trigger_particle_pt_cut the pT threshold enforced on the trigger
object
 * @return true, if all criteria are met, false otherwise
 */

bool matchParticle(const ROOT::Math::PtEtaPhiMVector &particle,
                   ROOT::RVec<float> &triggerobject_pts,
                   ROOT::RVec<float> &triggerobject_etas,
                   ROOT::RVec<float> &triggerobject_phis,
                   ROOT::RVec<int> &triggerobject_bits,
                   ROOT::RVec<int> &triggerobject_ids, const float &matchDeltaR,
                   const float &pt_cut, const float &eta_cut,
                   const int &trigger_particle_id_cut,
                   const int &triggerbit_cut,
                   const float &trigger_particle_pt_cut) {
    Logger::get("CheckTriggerMatch")->debug("Checking Triggerobjects");
    Logger::get("CheckTriggerMatch")
        ->debug("Total number of triggerobjects: {}", triggerobject_pts.size());
    for (std::size_t idx = 0; idx < triggerobject_pts.size(); ++idx) {
        Logger::get("CheckTriggerMatch")->debug("Triggerobject Nr. {}", idx);
        Logger::get("CheckTriggerMatch")
            ->debug("bit Value: {}", IntBits(triggerobject_bits[idx]));
        Logger::get("CheckTriggerMatch")
            ->debug("bit Value: {}", triggerobject_bits[idx]);
        auto triggerobject = ROOT::Math::RhoEtaPhiVectorF(
            0, triggerobject_etas[idx], triggerobject_phis[idx]);
        // We check the deltaR match as well as that the pt and eta of the
        // triggerobject are above the given thresholds
        bool deltaR = ROOT::Math::VectorUtil::DeltaR(triggerobject, particle) <
                      matchDeltaR;
        // if we don't want to do any matching here, the triggerbut_cut value is
        // -1
        Logger::get("CheckTriggerMatch")
            ->debug("bit Value: {}", triggerobject_bits[idx]);
        bool bit = (triggerbit_cut == -1) ||
                   (IntBits(triggerobject_bits[idx]).test(triggerbit_cut));
        bool id = triggerobject_ids[idx] == trigger_particle_id_cut;
        bool pt = particle.pt() > pt_cut;
        bool eta = abs(particle.eta()) < eta_cut;
        bool trigger_particle_pt =
            (trigger_particle_pt_cut < 0.) ||
            (triggerobject_pts[idx] > trigger_particle_pt_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("Partice Lorentz Vector: {}, {}, {}, {}", particle.pt(),
                    particle.eta(), particle.phi(), particle.mass());
        Logger::get("CheckTriggerMatch")
            ->debug("-------------------------------------------------------");
        Logger::get("CheckTriggerMatch")
            ->debug("deltaR/matchDeltaR Check: {}/{}", deltaR, matchDeltaR);
        Logger::get("CheckTriggerMatch")
            ->debug("deltaR Value: {}",
                    ROOT::Math::VectorUtil::DeltaR(triggerobject, particle));
        Logger::get("CheckTriggerMatch")
            ->debug("id/trigger_particle_id_cut Check: {}/{}", id,
                    trigger_particle_id_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("id Value: {}", triggerobject_ids[idx]);
        Logger::get("CheckTriggerMatch")
            ->debug("bit/triggerbit_cut Check: {}/{}", bit, triggerbit_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("bit Value: {}", IntBits(triggerobject_bits[idx]));
        Logger::get("CheckTriggerMatch")
            ->debug("pt/pt_cut Check: {}/{}", pt, pt_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("pt Value (trg): {}, pt Value (reco): {}",
                    triggerobject_pts[idx], particle.pt());
        Logger::get("CheckTriggerMatch")
            ->debug("eta/eta_cut Check: {}/{}", eta, eta_cut);
        Logger::get("CheckTriggerMatch")
            ->debug("eta (trg) Value: {}, eta (reco) Value: {}",
                    triggerobject_etas[idx], abs(particle.eta()));
        Logger::get("CheckTriggerMatch")
            ->debug("-------------------------------------------------------");
        if (deltaR && bit && id && pt && eta && trigger_particle_pt) {
            // remove the matching object from the object vectors so it cant be
            // matched by the next particle as well (if there is one)
            triggerobject_ids.erase(triggerobject_ids.begin() + idx);
            triggerobject_bits.erase(triggerobject_bits.begin() + idx);
            triggerobject_pts.erase(triggerobject_pts.begin() + idx);
            triggerobject_etas.erase(triggerobject_etas.begin() + idx);
            triggerobject_phis.erase(triggerobject_phis.begin() + idx);
            return true;
        }
    }
    return false;
};

/**
 * @brief Function to generate a trigger flag based on a trigger
 * object match only. This relies on the
 * trigger::matchParticle function which does the matching test.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle_p4 `ROOT::Math::PtEtaPhiMVector` of the object to be checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param pt_cut minimal pt value for the triggerobject
 * @param eta_cut maximal pt value for the triggerobject
 * @param trigger_particle_id_cut trigger id value the triggerobject has to
 * match (details can be found in the documentation of trigger::matchParticle)
 * @param triggerbit_cut trigger bit value the triggerobject has to match
 * (details can be found in the documentation of trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @param trigger_particle_pt_cut the pT threshold enforced on the trigger
 * object
 * @return a new dataframe containing the trigger object match flag column
 */

ROOT::RDF::RNode MatchSingleTriggerObject(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold, const float &trigger_particle_pt_cut) {

    auto triggermatch = [DeltaR_threshold, pt_cut, eta_cut,
                         trigger_particle_id_cut, triggerbit_cut,
                         trigger_particle_pt_cut](
                            const ROOT::Math::PtEtaPhiMVector &particle_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
        Logger::get("MatchSingleTriggerObject")->debug("Checking Trigger");
        Logger::get("MatchSingleTriggerObject")
            ->debug("Checking Triggerobject match with particles ....");
        bool match_result = matchParticle(
            particle_p4, triggerobject_pts, triggerobject_etas,
            triggerobject_phis, triggerobject_bits, triggerobject_ids,
            DeltaR_threshold, pt_cut, eta_cut, trigger_particle_id_cut,
            triggerbit_cut, trigger_particle_pt_cut);
        Logger::get("MatchSingleTriggerObject")
            ->debug("--->>>> match_result: {}", match_result);
        return match_result;
    };
    auto df1 =
        df.Define(triggerflag_name, triggermatch,
                  {particle_p4, triggerobject_bits, triggerobject_id,
                   triggerobject_pt, triggerobject_eta, triggerobject_phi});
    return df1;
}

/**
 * @brief Function to generate a trigger flag based on an hlt path and trigger
 * object matching for the given object. This relies on the
 * trigger::tagandprobe::matchParticle function which does the matching test.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle_p4 `ROOT::Math::PtEtaPhiMVector` of the object to be checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param hltpath name of the hlt path to be checked, this can be a valid regex.
 * If more than one matching HLT path is found, the function will throw an
 * exception, if no matching HLT path is found, the function will return a
 * dataframe with a flag of false for all entries
 * @param pt_cut minimal pt value for the selected particle
 * @param eta_cut maximal pt value for the selected particle
 * @param trigger_particle_id_cut trigger id value the triggerobject has to
 * match (details can be found in the documentation of trigger::matchParticle)
 * @param triggerbit_cut trigger bit value the triggerobject has to match
 * (details can be found in the documentation of trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @param trigger_particle_pt_cut minimal pt value for the triggerobject
 * @return a new dataframe containing the trigger flag column
 */

ROOT::RDF::RNode GenerateSingleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &hltpath, const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold, const float &trigger_particle_pt_cut) {

    auto triggermatch =
        [DeltaR_threshold, pt_cut, eta_cut, trigger_particle_id_cut,
         triggerbit_cut, trigger_particle_pt_cut, hltpath](
            bool hltpath_match, const ROOT::Math::PtEtaPhiMVector &particle_p4,
            ROOT::RVec<int> triggerobject_bits,
            ROOT::RVec<int> triggerobject_ids,
            ROOT::RVec<float> triggerobject_pts,
            ROOT::RVec<float> triggerobject_etas,
            ROOT::RVec<float> triggerobject_phis) {
            Logger::get("GenerateSingleTriggerFlag")->debug("Checking Trigger");
            Logger::get("CheckTriggerMatch")
                ->debug("Selected trigger: {}", hltpath);
            bool result = false;
            bool match_result = false;
            if (hltpath_match) {
                Logger::get("CheckTriggerMatch")
                    ->debug("Checking Triggerobject match with particles ....");
                match_result = matchParticle(
                    particle_p4, triggerobject_pts, triggerobject_etas,
                    triggerobject_phis, triggerobject_bits, triggerobject_ids,
                    DeltaR_threshold, pt_cut, eta_cut, trigger_particle_id_cut,
                    triggerbit_cut, trigger_particle_pt_cut);
            }
            result = hltpath_match & match_result;
            Logger::get("GenerateSingleTriggerFlag")
                ->debug("---> HLT Match: {}", hltpath_match);
            Logger::get("GenerateSingleTriggerFlag")
                ->debug("---> Total Match: {}", match_result);
            Logger::get("GenerateSingleTriggerFlag")
                ->debug("--->>>> result: {}", result);
            return result;
        };
    auto available_trigger = df.GetColumnNames();
    std::vector<std::string> matched_trigger_names;
    std::regex hltpath_regex = std::regex(hltpath);
    // loop over all available trigger names and check if the hltpath is
    // matching any of them
    for (auto &trigger : available_trigger) {
        if (std::regex_match(trigger, hltpath_regex)) {
            Logger::get("GenerateSingleTriggerFlag")
                ->debug("Found matching trigger: {}", trigger);
            matched_trigger_names.push_back(trigger);
        }
    }
    // if no matching trigger was found return the initial dataframe
    if (matched_trigger_names.size() == 0) {
        Logger::get("GenerateSingleTriggerFlag")
            ->info("No matching trigger for {} found, returning false for "
                   "trigger flag {}",
                   hltpath, triggerflag_name);
        auto df1 = df.Define(triggerflag_name, []() { return false; });
        return df1;
    } else if (matched_trigger_names.size() > 1) {
        Logger::get("GenerateSingleTriggerFlag")
            ->debug(
                "More than one matching trigger found, not implemented yet");
        throw std::invalid_argument(
            "received too many matching trigger paths, not implemented yet");
    } else {
        Logger::get("GenerateSingleTriggerFlag")
            ->debug("Found matching trigger: {}", matched_trigger_names[0]);
        auto df1 =
            df.Define(triggerflag_name, triggermatch,
                      {matched_trigger_names[0], particle_p4,
                       triggerobject_bits, triggerobject_id, triggerobject_pt,
                       triggerobject_eta, triggerobject_phi});
        return df1;
    }
}
/**
 * @brief Function to generate a trigger flag based on an hlt path and trigger
 * object matching for the given object. This relies on the
 * trigger::matchParticle function which does the matching test. The
 * implementation is similar to the trigger::GenerateSingleTriggerFlag function,
 * but here two objects have to be matched for a passing flag.
 *
 * @param df The input dataframe
 * @param triggerflag_name name of the output flag
 * @param particle1_p4 `ROOT::Math::PtEtaPhiMVector` of the first object to be
 * checked
 * @param particle2_p4 `ROOT::Math::PtEtaPhiMVector` of the second object to be
 * checked
 * @param triggerobject_bits name of the trigger object bits column in the
 * inputfile
 * @param triggerobject_id name of the trigger object id column in the inputfile
 * @param triggerobject_pt name of the trigger object pt column in the inputfile
 * @param triggerobject_eta name of the trigger object eta column in the
 * inputfile
 * @param triggerobject_phi name of the trigger object phi column in the
 * inputfile
 * @param hltpath name of the hlt path to be checked, this can be a valid regex.
 * If more than one matching HLT path is found, the function will throw an
 * exception, if no matching HLT path is found, the function will return a
 * dataframe with a flag of false for all entries
 * @param p1_pt_cut minimal pt value for the triggerobject matching the first
 * object
 * @param p2_pt_cut minimal pt value for the triggerobject matching the second
 * object
 * @param p1_eta_cut maximal pt value for the triggerobject matching the first
 * object
 * @param p2_eta_cut maximal pt value for the triggerobject matching the second
 * object
 * @param p1_trigger_particle_id_cut trigger id value the triggerobject matching
 * the first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_trigger_particle_id_cut trigger id value the triggerobject matching
 * the second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p1_triggerbit_cut trigger bit value the triggerobject matching the
 * first object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param p2_triggerbit_cut trigger bit value the triggerobject matching the
 * second object has to match (details can be found in the documentation of
 * trigger::matchParticle)
 * @param DeltaR_threshold maximal value for the deltaR between the
 * triggerobject and the input object to consider a match
 * @param p1_trigger_particle_pt_cut minimal value of the pT of the trigger
 * object matching the first object
 * @param p2_trigger_particle_pt_cut minimal value of the pT of the trigger
 * object matching the second object
 * @return a new dataframe containing the trigger flag column
 *
 */

ROOT::RDF::RNode GenerateDoubleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle1_p4, const std::string &particle2_p4,
    const std::string &triggerobject_bits, const std::string &triggerobject_id,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &hltpath,
    const float &p1_pt_cut, const float &p2_pt_cut, const float &p1_eta_cut,
    const float &p2_eta_cut, const int &p1_trigger_particle_id_cut,
    const int &p2_trigger_particle_id_cut, const int &p1_triggerbit_cut,
    const int &p2_triggerbit_cut, const float &DeltaR_threshold,
    const float &p1_trigger_particle_pt_cut,
    const float &p2_trigger_particle_pt_cut) {

    auto triggermatch = [DeltaR_threshold, p1_pt_cut, p2_pt_cut, p1_eta_cut,
                         p2_eta_cut, p1_trigger_particle_id_cut,
                         p2_trigger_particle_id_cut, p1_triggerbit_cut,
                         p2_triggerbit_cut, p1_trigger_particle_pt_cut,
                         p2_trigger_particle_pt_cut, hltpath](
                            bool hltpath_match,
                            const ROOT::Math::PtEtaPhiMVector &particle1_p4,
                            const ROOT::Math::PtEtaPhiMVector &particle2_p4,
                            ROOT::RVec<int> triggerobject_bits,
                            ROOT::RVec<int> triggerobject_ids,
                            ROOT::RVec<float> triggerobject_pts,
                            ROOT::RVec<float> triggerobject_etas,
                            ROOT::RVec<float> triggerobject_phis) {
        Logger::get("GenerateDoubleTriggerFlag")->debug("Checking Trigger");
        Logger::get("CheckTriggerMatch")
            ->debug("Selected trigger: {}", hltpath);
        bool result = false;
        bool match_result_p1 = false;
        bool match_result_p2 = false;
        if (hltpath_match) {
            Logger::get("GenerateDoubleTriggerFlag")
                ->debug("Checking Triggerobject match with particles ....");
            Logger::get("GenerateDoubleTriggerFlag")->debug("First particle");
            match_result_p1 = matchParticle(
                particle1_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, p1_pt_cut, p1_eta_cut,
                p1_trigger_particle_id_cut, p1_triggerbit_cut,
                p1_trigger_particle_pt_cut);
            Logger::get("GenerateDoubleTriggerFlag")->debug("Second particle");
            match_result_p2 = matchParticle(
                particle2_p4, triggerobject_pts, triggerobject_etas,
                triggerobject_phis, triggerobject_bits, triggerobject_ids,
                DeltaR_threshold, p2_pt_cut, p2_eta_cut,
                p2_trigger_particle_id_cut, p2_triggerbit_cut,
                p2_trigger_particle_pt_cut);
        }
        result = hltpath_match & match_result_p1 & match_result_p2;
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("---> HLT Match: {}", hltpath_match);
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("---> Total Match P1: {}", match_result_p1);
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("---> Total Match P2: {}", match_result_p2);
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("--->>>> result: {}", result);
        return result;
    };
    auto available_trigger = df.GetColumnNames();
    std::vector<std::string> matched_trigger_names;
    std::regex hltpath_regex = std::regex(hltpath);
    // loop over all available trigger names and check if the hltpath is
    // matching any of them
    for (auto &trigger : available_trigger) {
        if (std::regex_match(trigger, hltpath_regex)) {
            Logger::get("GenerateDoubleTriggerFlag")
                ->debug("Found matching trigger: {}", trigger);
            matched_trigger_names.push_back(trigger);
        }
    }
    // if no matching trigger was found return the initial dataframe
    if (matched_trigger_names.size() == 0) {
        Logger::get("GenerateDoubleTriggerFlag")
            ->info("No matching trigger for {} found, returning false for "
                   "trigger flag {}",
                   hltpath, triggerflag_name);
        auto df1 = df.Define(triggerflag_name, []() { return false; });
        return df1;
    } else if (matched_trigger_names.size() > 1) {
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug(
                "More than one matching trigger found, not implemented yet");
        throw std::invalid_argument(
            "received too many matching trigger paths, not implemented yet");
    } else {
        Logger::get("GenerateDoubleTriggerFlag")
            ->debug("Found matching trigger: {}", matched_trigger_names[0]);
        auto df1 =
            df.Define(triggerflag_name, triggermatch,
                      {matched_trigger_names[0], particle1_p4, particle2_p4,
                       triggerobject_bits, triggerobject_id, triggerobject_pt,
                       triggerobject_eta, triggerobject_phi});
        return df1;
    }
}

} // end namespace tagandprobe
/**
 * @brief Function to generate a new column containing the prescale value for a
 * trigger given run and lumiblock, read from an external JSON file
 *
 * @param df The input dataframe
 * @param prescale_columnname name of the output column for the prescale value
 * @param hlt_columnname name of the column for the HLT path
 * @param run_columnname name of the column for the run number
 * @param lumiblock_columnname name of the column for the lumiblock
 * @param prescale_json_file relative path to the JSON containing the values
 * @return a new dataframe containing the prescale column
 */

ROOT::RDF::RNode GetPrescaleValues(ROOT::RDF::RNode df,
                                   const std::string &prescale_columnname,
                                   const std::string &hlt_columnname,
                                   const std::string &run_columnname,
                                   const std::string &lumiblock_columnname,
                                   const std::string &prescale_json_file) {

    Logger::get("prescale")->debug("reading json from {}", prescale_json_file);
    std::ifstream i(prescale_json_file);
    nlohmann::json prescale_json = nlohmann::json::parse(i);

    auto get_prescale = [prescale_json](const Bool_t hlt, const UInt_t run,
                                        const UInt_t lumiblock) {
        int prescale = -1;

        // Logger::get("prescale")->debug("run, lumi: {},{}", run, lumiblock);

        if (hlt == false) {
            prescale = -2;
            // Logger::get("prescale")->debug("no HLT hit,  prescale value:
            // {}",prescale);
            return prescale;
        }

        const std::string s_run = std::to_string(run);

        if (prescale_json.find(s_run) != prescale_json.end()) {

            Logger::get("prescale")->debug("found run in JSON ...");
            unsigned highest_lumi = 1;

            for (auto &[i_key, i_value] : prescale_json[s_run].items()) {
                Logger::get("prescale")
                    ->debug("... checking lumi {}, prescale {} ...",
                            std::stoi(i_key), int(i_value));
                if (lumiblock >= std::stoi(i_key)) {
                    if (std::stoi(i_key) >= highest_lumi) {
                        highest_lumi = std::stoi(i_key);
                        prescale = i_value;
                        Logger::get("prescale")
                            ->debug("... assigning prescale value: {}",
                                    prescale);
                    }
                }
            }

        } else {
            prescale = -3;
            Logger::get("prescale")
                ->debug(
                    "could not find run and lumi in JSON, prescale value: {}",
                    prescale);
        }

        return prescale;
    };

    auto df1 =
        df.Define(prescale_columnname, get_prescale,
                  {hlt_columnname, run_columnname, lumiblock_columnname});

    return df1;
}

} // end namespace trigger
#endif /* GUARD_TRIGGERS_H */
