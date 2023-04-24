#ifndef GUARD_TRIPLESELECTION_H
#define GUARD_TRIPLESELECTION_H

#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "../include/pairselection.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include "bitset"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

typedef std::bitset<15> StatusBits;
namespace whtautau_tripleselection {
/// Function used to sort two particles based on the isolation and the
/// pt of the two particles. The function is used as the ordering
/// function for the
/// [ROOT::VecOps::Sort()](https://root.cern.ch/doc/master/group__vecops.html#ga882439c2ff958157d2990b52dd76f599)
/// algorithm. If two quantities are the same within a given epsilon of
/// 1e-5, the next criterion is applied. The sorting is done using the
/// following criterion odering:
/// -# Isolation of the first particle
/// -# pt of the first particle
/// -# Isolation of the second particle
/// -# pt of the second particle
///
/// \param lep1pt `ROOT::RVec<float>` containing pts of the first
/// particle
/// \param lep1iso `ROOT::RVec<float>` containing isolations of
/// the first particle
/// \param lep2pt `ROOT::RVec<float>` containing pts
/// of the second particle
/// \param lep2iso `ROOT::RVec<float>` containing
/// isolations of the second particle
///
/// \returns true or false based on the particle ordering.
auto compareForTriples(const ROOT::RVec<float> &lep1pt,
                     const ROOT::RVec<float> &lep1iso,
                     const ROOT::RVec<float> &lep2pt,
                     const ROOT::RVec<float> &lep2iso) {
    return [lep1pt, lep1iso, lep2pt, lep2iso](auto value_next,
                                              auto value_previous) {
        Logger::get("PairSelectionCompare")->debug("lep1 Pt: {}", lep1pt);
        Logger::get("PairSelectionCompare")->debug("lep1 Iso: {}", lep1iso);
        Logger::get("PairSelectionCompare")->debug("lep2 Pt: {}", lep2pt);
        Logger::get("PairSelectionCompare")->debug("lep2 Iso: {}", lep2iso);
        bool result = false;
        Logger::get("PairSelectionCompare")
            ->debug("Next pair: {}, {}", std::to_string(value_next.first),
                    std::to_string(value_next.second));
        Logger::get("PairSelectionCompare")
            ->debug("Previous pair: {}, {}",
                    std::to_string(value_previous.first),
                    std::to_string(value_previous.second));
        const auto i1_next = value_next.first;
        const auto i1_previous = value_previous.first;
        Logger::get("PairSelectionCompare")
            ->debug("i1_next: {}, i1_previous : {}", i1_next, i1_previous);
        // start with lep1 isolation
        const auto iso1_next = lep1iso.at(i1_next);
        const auto iso1_previous = lep1iso.at(i1_previous);
        Logger::get("PairSelectionCompare")
            ->debug("Isolations: {}, {}", iso1_next, iso1_previous);
        if (not utility::ApproxEqual(iso1_next, iso1_previous)) {
            result = iso1_next > iso1_previous;
        } else {
            // if too similar, compare lep1 pt
            Logger::get("PairSelectionCompare")
                ->debug("Isolation lep 1 too similar, taking pt 1");
            const auto pt1_next = lep1pt.at(i1_next);
            const auto pt1_previous = lep1pt.at(i1_previous);
            if (not utility::ApproxEqual(pt1_next, pt1_previous)) {
                result = pt1_next > pt1_previous;
            } else {
                // if too similar, compare lep2 iso
                const auto i2_next = value_next.second;
                const auto i2_previous = value_previous.second;
                Logger::get("PairSelectionCompare")
                    ->debug("Pt lep 1 too similar, taking lep2 iso");
                const auto iso2_next = lep2iso.at(i2_next);
                const auto iso2_previous = lep2iso.at(i2_previous);
                if (not utility::ApproxEqual(iso2_next, iso2_previous)) {
                    result = iso2_next > iso2_previous;
                } else {
                    // if too similar, compare lep2 pt
                    Logger::get("PairSelectionCompare")
                        ->debug("Isolation lep 2 too similar, taking pt 2");
                    const auto pt2_next = lep2pt.at(i2_next);
                    const auto pt2_previous = lep2pt.at(i2_previous);
                    result = pt2_next > pt2_previous;
                }
            }
        }
        Logger::get("PairSelectionCompare")
            ->debug("Returning result {}", result);
        return result;
    };
}
/**
 * @brief Function used to build a triple of GenParticles from the selected
 * triple. This uses the references of the reco particles to the gen
 * particles.
 *
 * @param df the Dataframe
 * @param recopair the column containing the DiTauPair vector
 * @param genindex_particle1 the column containing the index of the GenParticle
 * reference for the first pair particle
 * @param genindex_particle2 the column containing the index of the GenParticle
 * reference for the second pair particle
 * @param genindex_particle3 the column containing the index of the GenParticle
 * reference for the second pair particle
 * @param gentriple name of the new column containing the Gentriple
 * @return a new Dataframe with the Gentriple column
 */
ROOT::RDF::RNode buildgentriple(ROOT::RDF::RNode df, const std::string &recotriple,
                              const std::string &genindex_particle1,
                              const std::string &genindex_particle2,
                              const std::string &genindex_particle3,
                              const std::string &gentriple) {
    auto getGenTriple = [](const ROOT::RVec<int> &recotriple,
                         const ROOT::RVec<int> &genindex_particle1,
                         const ROOT::RVec<int> &genindex_particle2, const ROOT::RVec<int> &genindex_particle3) {
        ROOT::RVec<int> gentriple = {-1, -1, -1};
        Logger::get("buildgentriple")->debug("existing Triple: {}", recotriple);
        gentriple[0] = genindex_particle1.at(recotriple.at(0), -1);
        gentriple[1] = genindex_particle2.at(recotriple.at(1), -1);
        gentriple[2] = genindex_particle3.at(recotriple.at(2), -1);
        Logger::get("buildgentriple")
            ->debug("matching GenTriple: {}", gentriple);
        return gentriple;
    };
    return df.Define(gentriple, getGenTriple,
                     {recotriple, genindex_particle1, genindex_particle2, genindex_particle3});
}
/**
 * @brief Function to get the true gen-level triple from the event. The
 * triple is build by searching for the gen mother particles and the three requested
 * daughter particles. For each stable daughter particle found in the collection
 of gen particles, it is checked if the particle is a daughter of the requested
 mother particle.
 *
 * @param df the Dataframe
 * @param statusflags the column containing the status flags of the gen
 particles
 * @param status the column containing the status of the genparticles (status=1
 means stable)
 * @param pdgids the column containing the PDGID of the gen particles
 * @param motherids the column containing the index of the mother particle of
 the gen particles
 * @param pts the column containing the pt of the gen particles (used for
 sorting the particles by pt)
 * @param gentriple the output column containing the index of the two selected gen
 particles
 * @param mother_pdgid_1 the PDGID of the mother particle of the first particle
 * @param mother_pdgid_23 the PDGID of the mother particle of the second and third particle
 * @param daughter_1_pdgid the PDGID of the first daughter particle
 * @param daughter_2_pdgid the PDGID of the second daughter particle
 * @param daughter_3_pdgid the PDGID of the third daughter particle
 * @return auto the new Dataframe with the genpair column
 */
ROOT::RDF::RNode
buildtruegentriple(ROOT::RDF::RNode df, const std::string &statusflags,
                 const std::string &status, const std::string &pdgids,
                 const std::string &motherids, const std::string &pts,
                 const std::string &gentriple, const int mother_pdgid_1, const int mother_pdgid_23,
                 const int daughter_1_pdgid, const int daughter_2_pdgid, const int daughter_3_pdgid) {

    auto getTrueGenTriple = [mother_pdgid_1, mother_pdgid_23, daughter_1_pdgid,
                           daughter_2_pdgid, daughter_3_pdgid](const ROOT::RVec<int> &statusflags,
                                             const ROOT::RVec<int> &status,
                                             const ROOT::RVec<int> &pdgids,
                                             const ROOT::RVec<int> &motherids,
                                             const ROOT::RVec<float> &pts) {
        ROOT::RVec<int> gentriple = {-1, -1, -1};

        // first we build structs, one for each genparticle
        Logger::get("buildtruegentriple")
            ->debug("Starting to build True Gentriple for event");
        ROOT::RVec<GenParticle> genparticles;
        for (int i = 0; i < statusflags.size(); ++i) {
            GenParticle genparticle;
            genparticle.index = i;
            genparticle.status = status.at(i);
            genparticle.statusflag = StatusBits(statusflags.at(i));
            genparticle.pdgid = abs(pdgids.at(i));
            genparticle.motherid = motherids.at(i);
            genparticles.push_back(genparticle);
        }
        Logger::get("buildtruegentriple")->debug("genparticles: ");
        for (auto &genparticle : genparticles) {
            Logger::get("buildtruegentriple")
                ->debug("|--------------------------------------------");
            Logger::get("buildtruegentriple")
                ->debug("|    Index: {}", genparticle.index);
            Logger::get("buildtruegentriple")
                ->debug("|    Status: {}", genparticle.status);
            Logger::get("buildtruegentriple")
                ->debug("|     Statusflag: {}", genparticle.statusflag);
            Logger::get("buildtruegentriple")
                ->debug("|    Pdgid: {}", genparticle.pdgid);
            Logger::get("buildtruegentriple")
                ->debug("|    motherid: {}", genparticle.motherid);
            Logger::get("buildtruegentriple")
                ->debug("|--------------------------------------------");
        }
        auto gen_candidates_1 = ROOT::VecOps::Filter(
            genparticles, [daughter_1_pdgid](const GenParticle &genparticle) {
                return (daughter_1_pdgid == genparticle.pdgid) &&
                       (genparticle.status == 1);
            });
        auto gen_candidates_2 = ROOT::VecOps::Filter(
            genparticles, [daughter_2_pdgid](const GenParticle &genparticle) {
                return (daughter_2_pdgid == genparticle.pdgid) &&
                       (genparticle.status == 1);
            });
        auto gen_candidates_3 = ROOT::VecOps::Filter(
            genparticles, [daughter_3_pdgid](const GenParticle &genparticle) {
                return (daughter_3_pdgid == genparticle.pdgid) &&
                       (genparticle.status == 1);
            });
        // in the case of three different daughter particles, we need find
        // each one seperately
        for (const auto &gen_candidate_1 : gen_candidates_1) {
            if (check_mother(genparticles, gen_candidate_1.index,
                                mother_pdgid_1)) {
                gentriple[0] = gen_candidate_1.index;
                break;
            }
        }
        for (const auto &gen_candidate_2 : gen_candidates_2) {
            if (check_mother(genparticles, gen_candidate_2.index,
                                mother_pdgid_23)) {
                gentriple[1] = gen_candidate_2.index;
                break;
            }
        }
        for (const auto &gen_candidate_3 : gen_candidates_3) {
            if (check_mother(genparticles, gen_candidate_3.index,
                                mother_pdgid_23)) {
                gentriple[2] = gen_candidate_3.index;
                break;
            }
        }
        Logger::get("buildtruegentriple")
            ->debug("Selected Particles: {} {} {}", gentriple[0], gentriple[1], gentriple[2]);
        if (gentriple[0] != -1 && gentriple[1] != -1 && gentriple[2] != -1) {
            Logger::get("buildtruegentriple")
                ->debug("viable daughter particles found");
//            return gentriple;
        }
        return gentriple;
    };
    return df.Define(gentriple, getTrueGenTriple,
                     {statusflags, status, pdgids, motherids, pts});
}
/// This function flags events, where a suitable particle triple is found.
/// A triple is considered suitable, if a TripleSelectionAlgo (like
/// whtautau_tripleselection::three_flavor::TripleSelectionAlgo) returns indices, that are
/// not -1. Events, where any of the particle indices is -1 are vetoed
/// by this filter.
///
/// \param df The input dataframe
/// \param flagname The name of the generated flag column
/// \param pairname The name of the column, containing the indices of the
/// particles in the particle quantity vectors.
/// \returns a dataframe with the
/// new flag
ROOT::RDF::RNode flagGoodTriples(ROOT::RDF::RNode df, const std::string &flagname,
                               const std::string &triplename) {
    using namespace ROOT::VecOps;
    return df.Define(
        flagname,
        [](const ROOT::RVec<int> &triple) { return bool(Min(triple) >= 0); },
        {triplename});
}
namespace three_flavor {
/// Implementation of the triple selection algorithm. First, only events
/// that contain at least two goodleptons and one goodTau are considered.
/// Events contain at least two good leptons and one good tau, if the
/// the two leptonmasks have nonzero elements. These masks are
/// constructed using the functions from the physicsobject namespace
/// (e.g. physicsobject::CutPt). The argument triple gives information wheather the emt or the met channel is considered.
/// Events with two good electrons or two good muons are vetos immediately.
/// For the fake rate estimation, base leptons are also considered.
/// 
/// \returns an `ROOT::RVec<int>` with three values, the first one beeing
/// the lepton from the W index, the second one beeing the lepton from the tau index and the third one the hadronic tau index.
auto TripleSelectionAlgo(const float &mindeltaR_leptau, const float &mindeltaR_leplep, const std::string &triple) {
    Logger::get("three_flavor::TripleSelectionAlgo")
        ->debug("Setting up algorithm");
    return [mindeltaR_leptau, mindeltaR_leplep, triple](const ROOT::RVec<float> &tau_pt,
                       const ROOT::RVec<float> &tau_eta,
                       const ROOT::RVec<float> &tau_phi,
                       const ROOT::RVec<float> &tau_mass,
                       const ROOT::RVec<float> &tau_iso,
                       const ROOT::RVec<float> &electron_pt,
                       const ROOT::RVec<float> &electron_eta,
                       const ROOT::RVec<float> &electron_phi,
                       const ROOT::RVec<float> &electron_mass,
                       const ROOT::RVec<float> &electron_iso,
                       const ROOT::RVec<float> &muon_pt,
                       const ROOT::RVec<float> &muon_eta,
                       const ROOT::RVec<float> &muon_phi,
                       const ROOT::RVec<float> &muon_mass,
                       const ROOT::RVec<float> &muon_iso,
                       const ROOT::RVec<int> &good_electron_mask,
                       const ROOT::RVec<int> &base_electron_mask,
                       const ROOT::RVec<int> &good_muon_mask,
                       const ROOT::RVec<int> &base_muon_mask,
                       const ROOT::RVec<int> &tau_mask) {

        Logger::get("three_flavor::TripleSelectionAlgo")
            ->debug("masks: org tau {}, good ele {}, good mu {}", tau_mask, good_electron_mask, good_muon_mask);    
        Logger::get("three_flavor::TripleSelectionAlgo")
            ->debug("base masks: base ele {}, base mu {}", base_electron_mask, base_muon_mask);  
        ROOT::RVec<int> selected_triple = {-1, -1, -1};
        //get indices of candidates that pass the given mask
        const auto original_tau_indices = ROOT::VecOps::Nonzero(tau_mask);
        const auto ind_good_electrons =
            ROOT::VecOps::Nonzero(good_electron_mask);
        const auto ind_good_muons = ROOT::VecOps::Nonzero(good_muon_mask);
        const auto ind_base_electrons =
            ROOT::VecOps::Nonzero(base_electron_mask);
        const auto ind_base_muons = ROOT::VecOps::Nonzero(base_muon_mask);
        // if the event has more than 2 good leptons or no good tau or no base leptons decline the event
        if (original_tau_indices.size() < 1 or
            ind_good_electrons.size() > 1 or
            ind_good_muons.size() > 1 or ind_base_electrons.size() < 1 or ind_base_muons.size() < 1) {
            Logger::get("three_flavor::TripleSelectionAlgo")
            ->debug("no valid triple found");
            return selected_triple;
        }

        ROOT::RVec<int> selected_ele_indices;
        ROOT::RVec<int> selected_mu_indices;
        const auto base_ele_iso =
            ROOT::VecOps::Take(electron_iso, ind_base_electrons);
        const auto base_mu_iso =
            ROOT::VecOps::Take(muon_iso, ind_base_muons);
        // if the event contains exactly one good lepton from each flavor, these are chosen (signal region). If at least one lepton fails the good lepton requirements, the event is chosen for the fake rate measurements
        if (ind_good_electrons.size() == 1) {
            selected_ele_indices.push_back(ind_good_electrons.at(0));
        }
        else {
            // sort base leptons in their isolation in ascending order
            const auto sorted_base_ele_idx = ROOT::VecOps::Argsort(base_ele_iso);
            for (const auto &ele_cand : sorted_base_ele_idx){
                Logger::get("three_flavor::TripleSelectionAlgo")
                    ->debug("index ele cand {}, sorted ele {}, base ele iso {}, ind_base_electrons {}, base ele iso {}", ele_cand,sorted_base_ele_idx, base_ele_iso, ind_base_electrons, electron_iso);
                selected_ele_indices.push_back(ind_base_electrons.at(ele_cand));
            }
        }
        if (ind_good_muons.size() == 1) {
            selected_mu_indices.push_back(ind_good_muons.at(0));
        }
        else {
            // sort base leptons by their isolation in ascending order
            const auto sorted_base_mu_idx = ROOT::VecOps::Argsort(base_mu_iso);
            for (auto &mu_cand : sorted_base_mu_idx){
                Logger::get("three_flavor::TripleSmuctionAlgo")
                    ->debug("index mu cand {}, sorted mu {}, base_mu_iso {}", mu_cand,sorted_base_mu_idx, base_mu_iso);
                selected_mu_indices.push_back(ind_base_muons.at(mu_cand));
            }
        }       

        // const auto selected_tau_pt =
        //     ROOT::VecOps::Take(tau_pt, original_tau_indices);
        const auto selected_tau_iso =
            ROOT::VecOps::Take(tau_iso, original_tau_indices);

        // sort hadronic taus by their deeptau vs. jets score in descending order
        const auto sorted_tau_idx = ROOT::VecOps::Argsort(
            -1. * selected_tau_iso);
        Logger::get("three_flavor::TripleSelectionAlgo")
            ->debug("sorted taus {}", sorted_tau_idx);

        //Construct four vectors of the three candidates and check the DeltaR requirements. Triples are rejected if the candidates are too close. Since the hadronic tau has to be good also in the jet to lepton fake rate estimation, the for loop starts with the taus, that are sorted by isolation. As muons are required to be either global or tracker, the muon ID is tighter compared to electrons. therefore the muon for loop follows. For the signal region the algorithm has no effect, since we select exaclty one good lepton. 
        for (auto &candidate_tau : sorted_tau_idx) {
            Logger::get("three_flavor::TripleSelectionAlgo")
                ->debug("{} debuging1 {}", original_tau_indices.at(candidate_tau), sorted_tau_idx.at(candidate_tau));
            ROOT::Math::PtEtaPhiMVector tau = ROOT::Math::PtEtaPhiMVector(
            tau_pt.at(original_tau_indices.at(candidate_tau)), tau_eta.at(original_tau_indices.at(candidate_tau)), tau_phi.at(original_tau_indices.at(candidate_tau)),
            tau_mass.at(original_tau_indices.at(candidate_tau)));
            for (auto &candidate_mu : selected_mu_indices) {
                ROOT::Math::PtEtaPhiMVector muon = ROOT::Math::PtEtaPhiMVector(
                muon_pt.at(candidate_mu), muon_eta.at(candidate_mu),
                muon_phi.at(candidate_mu), muon_mass.at(candidate_mu));
                for (auto &candidate_ele : selected_ele_indices) {
                    ROOT::Math::PtEtaPhiMVector electron = ROOT::Math::PtEtaPhiMVector(
                    electron_pt.at(candidate_ele), electron_eta.at(candidate_ele),
                    electron_phi.at(candidate_ele), electron_mass.at(candidate_ele));
                    Logger::get("three_flavor::TripleSelectionAlgo")
                            ->debug("deltaR electron muon {}", ROOT::Math::VectorUtil::DeltaR(electron, muon)); 
                    Logger::get("three_flavor::TripleSelectionAlgo")
                        ->debug("DeltaR(electron, tau): {}", ROOT::Math::VectorUtil::DeltaR(electron, tau));
                    Logger::get("three_flavor::TripleSelectionAlgo")
                        ->debug("DeltaR(muon, tau): {}", ROOT::Math::VectorUtil::DeltaR(muon, tau));
                    if (ROOT::Math::VectorUtil::DeltaR(electron, tau) > mindeltaR_leptau && ROOT::Math::VectorUtil::DeltaR(muon, tau) > mindeltaR_leptau && ROOT::Math::VectorUtil::DeltaR(electron, muon) > mindeltaR_leplep) {
                        Logger::get("three_flavor::TripleSmuctionAlgo")
                            ->debug("indices of selected candidates: muon {}, electron {}, tau {}", candidate_mu, candidate_ele, original_tau_indices.at(candidate_tau));
                        if ((triple == "emt") && (muon.Pt() > electron.Pt())) {
                            // only select triples with the pt_electron > pt_muon
                                return selected_triple;
                        }
                        else if ((triple == "met") && (electron.Pt() > muon.Pt())) {
                            // only select triples with the pt_muon > pt_electron
                            return selected_triple;
                        }
                        else if ((triple == "emt") && (muon.Pt() < electron.Pt())) {
                            selected_triple = {static_cast<int>(candidate_ele), static_cast<int>(candidate_mu), 
                                                static_cast<int>(original_tau_indices.at(candidate_tau))};
                            Logger::get("three_flavor::TripleSelectionAlgo")
                                ->debug("valid triple found {}", selected_triple);
                            return selected_triple;
                        }
                        else if ((triple == "met") && (electron.Pt() < muon.Pt())) {
                            selected_triple = {static_cast<int>(candidate_mu), static_cast<int>(candidate_ele), 
                                                static_cast<int>(original_tau_indices.at(candidate_tau))};
                            Logger::get("three_flavor::TripleSelectionAlgo")
                                ->debug("valid triple found {}", selected_triple);
                            return selected_triple;
                        }
                    }
                    else {
                        Logger::get("three_flavor::TripleSelectionAlgo")
                            ->debug("building failed because of deltaR requirements");
                    }
                }
            }
        }
        Logger::get("three_flavor::TripleSelectionAlgo")
            ->debug("no valid triple found {}", selected_triple);
        return selected_triple;
    };
} // end namespace TripleSelectionAlgo

/// Implementation of the triple selection algorithm without considering electrons to compare the results of the emt channel with the higgs to tau tau algorithm (mt channel). First, only events
/// that contain at least two goodleptons and one goodTau are considered.
/// Events contain at least two good leptons and one good tau, if the
/// tau_mask the two leptonmasks have nonzero elements. These masks are
/// constructed using the functions from the physicsobject namespace
/// (e.g. physicsobject::CutPt).
/// To estimate the jet to lepton fake rate, we select also events with base leptons
/// \returns an `ROOT::RVec<int>` with three values, the first one beeing
/// the lepton from the W index, the second one beeing the lepton from the tau index and the third one the hadronic tau index.
auto TripleSelectionWithoutEleAlgo(const float &mindeltaR_leptau) {
    Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
        ->debug("Setting up algorithm");
    return [mindeltaR_leptau](const ROOT::RVec<float> &tau_pt,
                       const ROOT::RVec<float> &tau_eta,
                       const ROOT::RVec<float> &tau_phi,
                       const ROOT::RVec<float> &tau_mass,
                       const ROOT::RVec<float> &tau_iso,
                       const ROOT::RVec<float> &electron_pt,
                       const ROOT::RVec<float> &electron_eta,
                       const ROOT::RVec<float> &electron_phi,
                       const ROOT::RVec<float> &electron_mass,
                       const ROOT::RVec<float> &electron_iso,
                       const ROOT::RVec<float> &muon_pt,
                       const ROOT::RVec<float> &muon_eta,
                       const ROOT::RVec<float> &muon_phi,
                       const ROOT::RVec<float> &muon_mass,
                       const ROOT::RVec<int> &electron_mask,
                       const ROOT::RVec<int> &muon_mask,
                       const ROOT::RVec<int> &tau_mask) {
        // first entry is the lepton index,
        // second entry is the tau index
        ROOT::RVec<int> selected_triple = {-1, -2, -2};
        const auto original_tau_indices = ROOT::VecOps::Nonzero(tau_mask);
        const auto original_electron_indices =
            ROOT::VecOps::Nonzero(electron_mask);
        const auto original_muon_indices = ROOT::VecOps::Nonzero(muon_mask);
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("masks: org tau {}, org ele {}, org mu {}", tau_mask, electron_mask, muon_mask);

        // only select events with one good muon and one good electron and at least one tau
        if (original_tau_indices.size() < 1 or
            original_muon_indices.size() != 1) {
            return selected_triple;
        }
        const auto selected_tau_pt =
            ROOT::VecOps::Take(tau_pt, original_tau_indices);
        const auto selected_tau_iso =
            ROOT::VecOps::Take(tau_iso, original_tau_indices);
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("Running algorithm on good taus and leptons");
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("original electron indices {}", original_electron_indices);
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("original muon indices {}", original_muon_indices);
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("original tau indices {}", original_tau_indices);

        // sort hadronic taus by their deeptau vs. jets score in descending order
        const auto sorted_tau_idx = ROOT::VecOps::Argsort(
            -1. * selected_tau_iso);
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("sorted taus {}", sorted_tau_idx);

        // construct the four vectors of the selected leptons         
        if (original_electron_indices.size() > 0) {
            // sort electrons by their isolation in descending order
            const auto selected_ele_iso =
            ROOT::VecOps::Take(electron_iso, original_electron_indices);
            const auto sorted_ele_iso = ROOT::VecOps::Argsort(
            selected_ele_iso);
            Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("sorted eles {}", sorted_ele_iso);
            ROOT::Math::PtEtaPhiMVector electron = ROOT::Math::PtEtaPhiMVector(
                electron_pt.at(original_electron_indices.at(sorted_ele_iso.at(0))), electron_eta.at(original_electron_indices.at(sorted_ele_iso.at(0))),
                electron_phi.at(original_electron_indices.at(sorted_ele_iso.at(0))), electron_mass.at(original_electron_indices.at(sorted_ele_iso.at(0))));
            Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
                ->debug("{} electron vector: {}", original_electron_indices, electron);
        }
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("muon pts: {}", muon_pt); 
        ROOT::Math::PtEtaPhiMVector muon = ROOT::Math::PtEtaPhiMVector(
            muon_pt.at(original_muon_indices.at(0)), muon_eta.at(original_muon_indices.at(0)),
            muon_phi.at(original_muon_indices.at(0)), muon_mass.at(original_muon_indices.at(0)));
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("{} muon vector: {}", original_muon_indices, muon);    

        //Construct four vectors of hadronic tau candidates and check the DeltaR requirements to the leptons. Triples are rejected if the candidates are too close
        for (auto &candidate : sorted_tau_idx) {
            Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
                ->debug("{} debuging1 {}", original_tau_indices.at(candidate), sorted_tau_idx.at(candidate));
            ROOT::Math::PtEtaPhiMVector tau = ROOT::Math::PtEtaPhiMVector(
                tau_pt.at(original_tau_indices.at(candidate)), tau_eta.at(original_tau_indices.at(candidate)), tau_phi.at(original_tau_indices.at(candidate)),
                tau_mass.at(original_tau_indices.at(candidate)));
            if (ROOT::Math::VectorUtil::DeltaR(muon, tau) > mindeltaR_leptau) {
                if (original_electron_indices.size() > 0) {                           
                selected_triple = {static_cast<int>(original_electron_indices.at(0)), static_cast<int>(original_muon_indices.at(0)), 
                                static_cast<int>(original_tau_indices.at(candidate))};
                Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
                ->debug("selected triple {}", selected_triple);
                }
                else {                           
                selected_triple = {-1, static_cast<int>(original_muon_indices.at(0)), 
                                static_cast<int>(original_tau_indices.at(candidate))};
                Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
                ->debug("selected triple {}", selected_triple);
                }
                break;
            }
        }
        Logger::get("three_flavor::TripleSelectionWithoutEleAlgo")
            ->debug("Final triple {} {} {}", selected_triple[0], selected_triple[1], selected_triple[2]);
        return selected_triple;
    };  
} // end namespace TripleSelectionWithoutElectronAlgo
} // end namespace three_flavor
namespace two_flavor {

/// Implementation of the triple selection algorithm. First, only events
/// that contain at least two goodleptons with same flavour and one goodTau are considered.
/// Events contain at least two good leptons and one good tau, if the
/// the leptonmask and the taumask have nonzero elements. These masks are
/// constructed using the functions from the physicsobject namespace
/// (e.g. physicsobject::CutPt). The argument triple gives information wheather the emt or the met channel is considered.
/// To estimate the jet to lepton fake rate, we select also events with base leptons
/// \returns an `ROOT::RVec<int>` with three values, the first one beeing
/// the lepton from the W index, the second one beeing the lepton from the tau index and the third one the hadronic tau index.
auto TripleSelectionAlgo(const float &mindeltaR_leptau, const float &mindeltaR_leplep) {
    Logger::get("two_flavor::TripleSelectionAlgo")
        ->debug("Setting up algorithm");
    return [mindeltaR_leptau, mindeltaR_leplep](const ROOT::RVec<float> &tau_pt,
                       const ROOT::RVec<float> &tau_eta,
                       const ROOT::RVec<float> &tau_phi,
                       const ROOT::RVec<float> &tau_mass,
                       const ROOT::RVec<float> &tau_iso,
                       const ROOT::RVec<float> &muon_pt,
                       const ROOT::RVec<float> &muon_eta,
                       const ROOT::RVec<float> &muon_phi,
                       const ROOT::RVec<float> &muon_mass,
                       const ROOT::RVec<float> &muon_iso,
                       const ROOT::RVec<int> &good_muon_mask,
                       const ROOT::RVec<int> &base_muon_mask,
                       const ROOT::RVec<int> &tau_mask) {

        ROOT::RVec<int> selected_triple = {-1, -1, -1};
        const auto original_tau_indices = ROOT::VecOps::Nonzero(tau_mask);
        const auto ind_good_muons = ROOT::VecOps::Nonzero(good_muon_mask);
        const auto ind_base_muons = ROOT::VecOps::Nonzero(base_muon_mask);
        Logger::get("two_flavor::TripleSelectionAlgo")
            ->debug("good masks: tau {}, mu {}", tau_mask, good_muon_mask);
        Logger::get("two_flavor::TripleSelectionAlgo")
            ->debug("base mask: mu {}", base_muon_mask);  
        // reject events with more than two good muons or less than one good tau or less than two base muons
        if (original_tau_indices.size() < 1 or
            ind_good_muons.size() > 2 or ind_base_muons.size() < 2) {
            Logger::get("two_flavor::TripleSelectionAlgo")
            ->debug("no valid triple found");
            return selected_triple;
        }
        ROOT::RVec<int> selected_mu_indices;
        const auto base_mu_iso =
            ROOT::VecOps::Take(muon_iso, ind_base_muons);
        // if the event contains exactly two good muons, these are chosen (for signal region). In the case, that only bad muons or only one good muon are contained in the event, the most isolated muons from the base muon selection are chosen first for the triple (for fake rate estimation)
        if (ind_good_muons.size() == 2) {
            selected_mu_indices.push_back(ind_good_muons.at(0));
            selected_mu_indices.push_back(ind_good_muons.at(1));
            Logger::get("three_flavor::TripleSmuctionAlgo")
                ->debug("two good muons found");
        }
        else {
            // sort base leptons in their isolation in ascending order
            const auto sorted_base_mu_idx = ROOT::VecOps::Argsort(base_mu_iso);
            for (auto &mu_cand : sorted_base_mu_idx){
                Logger::get("three_flavor::TripleSmuctionAlgo")
                    ->debug("mu cand {}, sorted mu {}, base_mu_iso {}", mu_cand,sorted_base_mu_idx, base_mu_iso);
                selected_mu_indices.push_back(ind_base_muons.at(mu_cand));
            }            
        }       
        Logger::get("two_flavour::TripleSelectionAlgo")
            ->debug("Running algorithm on good taus and leptons");
        Logger::get("two_flavour::TripleSelectionAlgo")
            ->debug("selected muon indices {}", selected_mu_indices);
        Logger::get("two_flavour::TripleSelectionAlgo")
            ->debug("original tau indices {}", original_tau_indices);

        const auto selected_tau_iso =
            ROOT::VecOps::Take(tau_iso, original_tau_indices);
        Logger::get("two_flavour::TripleSelectionAlgo")
            ->debug("tau mask {}", tau_mask);
        Logger::get("two_flavour::TripleSelectionAlgo")
            ->debug("selected tau iso {}", selected_tau_iso);

        // sort hadronic taus by their deeptau vs. jets score in descending order
        const auto sorted_tau_idx = ROOT::VecOps::Argsort(
            -1. * selected_tau_iso);
        Logger::get("two_flavour::TripleSelectionAlgo")
            ->debug("sorted taus {}", sorted_tau_idx);

        //Construct four vectors of hadronic tau and muon candidates and check the DeltaR requirements. Triples are rejected if the candidates are too close. For the signal region exactly two good muons are required, therefore the loops over the muons have no influence on the pairing. For the fake rate region, more muons are allowed. The algorithm prefers triples with the most isolated muon. 
        for (auto &candidate_tau : sorted_tau_idx) {
            Logger::get("two_flavour::TripleSelectionAlgo")
                ->debug("original_tau_indices.at(candidate_tau) {},  sorted_tau_idx.at(candidate_tau) {}", original_tau_indices.at(candidate_tau), sorted_tau_idx.at(candidate_tau));
            ROOT::Math::PtEtaPhiMVector tau = ROOT::Math::PtEtaPhiMVector(
                tau_pt.at(original_tau_indices.at(candidate_tau)), tau_eta.at(original_tau_indices.at(candidate_tau)), tau_phi.at(original_tau_indices.at(candidate_tau)),
                tau_mass.at(original_tau_indices.at(candidate_tau)));
            for (int i=0; i<selected_mu_indices.size()-1; i=i+1) {
                ROOT::Math::PtEtaPhiMVector muon_1 = ROOT::Math::PtEtaPhiMVector(
                    muon_pt.at(selected_mu_indices.at(i)), muon_eta.at(selected_mu_indices.at(i)),
                    muon_phi.at(selected_mu_indices.at(i)), muon_mass.at(selected_mu_indices.at(i)));
                for (int j=i+1;j<selected_mu_indices.size(); j=j+1) {
                    ROOT::Math::PtEtaPhiMVector muon_2 = ROOT::Math::PtEtaPhiMVector(
                        muon_pt.at(selected_mu_indices.at(j)), muon_eta.at(selected_mu_indices.at(j)),
                        muon_phi.at(selected_mu_indices.at(j)), muon_mass.at(selected_mu_indices.at(j)));
                        Logger::get("two_flavour::TripleSelectionAlgo")
                            ->debug("DeltaR(muon_1, tau): {}",ROOT::Math::VectorUtil::DeltaR(muon_1, tau));
                        Logger::get("two_flavour::TripleSelectionAlgo")
                            ->debug("DeltaR(muon_2, tau): {}", ROOT::Math::VectorUtil::DeltaR(muon_2, tau));
                        Logger::get("two_flavour::TripleSelectionAlgo")
                            ->debug("DeltaR(muon_1, muon_2): {}", ROOT::Math::VectorUtil::DeltaR(muon_1, muon_2));
                    if (ROOT::Math::VectorUtil::DeltaR(muon_1, tau) > mindeltaR_leptau && ROOT::Math::VectorUtil::DeltaR(muon_2, tau) > mindeltaR_leptau && ROOT::Math::VectorUtil::DeltaR(muon_1, muon_2) > mindeltaR_leplep) {
                        Logger::get("two_flavour::TripleSelectionAlgo")
                            ->debug("pt muon_1 {}, pt muon_2 {}", muon_1.Pt(), muon_2.Pt());
                        if (muon_1.Pt() > muon_2.Pt()) {
                            Logger::get("two_flavour::TripleSelectionAlgo")
                                ->debug("Selected original triple indices: mu_1 = {}, mu_2 = {} , tau = {}", selected_mu_indices.at(i), selected_mu_indices.at(j), original_tau_indices.at(candidate_tau));           
                            selected_triple = {static_cast<int>(selected_mu_indices.at(i)), static_cast<int>(selected_mu_indices.at(j)), static_cast<int>(original_tau_indices.at(candidate_tau))};
                            Logger::get("two_flavour::TripleSelectionAlgo")
                                ->debug("selected triple {}", selected_triple);
                        }
                        else {
                            Logger::get("two_flavour::TripleSelectionAlgo")
                                ->debug("Selected original triple indices: mu_1 = {}, mu_2 = {} , tau = {}", selected_mu_indices.at(j), selected_mu_indices.at(i), original_tau_indices.at(candidate_tau));           
                            selected_triple = {static_cast<int>(selected_mu_indices.at(j)), static_cast<int>(selected_mu_indices.at(i)), static_cast<int>(original_tau_indices.at(candidate_tau))};
                            Logger::get("two_flavour::TripleSelectionAlgo")
                                ->debug("selected triple {}", selected_triple);
                        }
                        return selected_triple;
                    }
                }
            }
        }
        Logger::get("three_flavor::TripleSelectionAlgo")
        ->debug("no valid triple found {}", selected_triple);
        return selected_triple;
    };
} // end namespace TripleSelectionAlgo
} // end namespace two_flavor
namespace lep_tautau {

/// Implementation of the triple selection algorithm. First, only events
/// that contain one good lepton and at least two good taus are considered.
/// Events contain one good lepton and at least two good taus, if the
/// tau_mask and the leptonmask have nonzero elements. These masks are
/// constructed using the functions from the physicsobject namespace
/// (e.g. physicsobject::CutPt). The argument triple gives information wheather the emt or the met channel is considered.
///
/// \returns an `ROOT::RVec<int>` with three values, the first one beeing
/// the lepton from the W index, the second one beeing the lepton from the tau index and the third one the hadronic tau index.
auto TripleSelectionAlgo(const float &mindeltaR_leptau, const float &mindeltaR_tautau) {
    Logger::get("lep_tautau::TripleSelectionAlgo")
        ->debug("Setting up algorithm");
    return [mindeltaR_leptau, mindeltaR_tautau](const ROOT::RVec<float> &tau_pt,
                       const ROOT::RVec<float> &tau_eta,
                       const ROOT::RVec<float> &tau_phi,
                       const ROOT::RVec<float> &tau_mass,
                       const ROOT::RVec<float> &tau_iso,
                       const ROOT::RVec<float> &lepton_pt,
                       const ROOT::RVec<float> &lepton_eta,
                       const ROOT::RVec<float> &lepton_phi,
                       const ROOT::RVec<float> &lepton_mass,
                       const ROOT::RVec<int> &lepton_mask,
                       const ROOT::RVec<int> &tau_mask) {
                        
        ROOT::RVec<int> selected_triple = {-1, -1, -1};
        const auto original_tau_indices = ROOT::VecOps::Nonzero(tau_mask);
        const auto original_lepton_indices = ROOT::VecOps::Nonzero(lepton_mask);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("masks: org tau {}, org lepton {}", tau_mask, lepton_mask);

        // only select events with at least two good taus and excatly one lepton (muon or electron)
        if (original_tau_indices.size() < 2 or
            original_lepton_indices.size() != 1) {
            return selected_triple;
        }
        
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("Running algorithm on good taus and leptons");
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("original lepton indices {}", original_lepton_indices);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("original tau indices {}", original_tau_indices);

        const auto selected_tau_pt =
            ROOT::VecOps::Take(tau_pt, original_tau_indices);
        const auto selected_tau_iso =
            ROOT::VecOps::Take(tau_iso, original_tau_indices);
        const auto selected_lepton_pt =
            ROOT::VecOps::Take(lepton_pt, original_lepton_indices);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("tau mask {}", tau_mask);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("tau pt {}", tau_pt);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("lepton pt {}", lepton_pt);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("selected tau pt {}", selected_tau_pt);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("selected tau iso {}", selected_tau_iso);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("selected lepton pt {}", selected_lepton_pt);

        // sort hadronic taus by their iso in descending order
        const auto sorted_tau_idx = ROOT::VecOps::Argsort(
            -1. * selected_tau_iso);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("sorted taus {}", sorted_tau_idx);

        // construct the four vectors of the selected lepton
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("lepton pts: {}", lepton_pt); 
        ROOT::Math::PtEtaPhiMVector lepton = ROOT::Math::PtEtaPhiMVector(
        lepton_pt.at(original_lepton_indices.at(0)), lepton_eta.at(original_lepton_indices.at(0)),
        lepton_phi.at(original_lepton_indices.at(0)), lepton_mass.at(original_lepton_indices.at(0)));

        const auto pair_indices = ROOT::VecOps::Combinations(
            selected_tau_pt, 2); // Gives indices of tau-tau pairs
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("Hadronic tau pairs: {} {}", pair_indices[0], pair_indices[1]);

        const auto pairs = ROOT::VecOps::Construct<std::pair<UInt_t, UInt_t>>(
            pair_indices[0], pair_indices[1]);
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("Hadronic tau pairs size: {}", pairs.size());
        int counter = 0;
        for (auto &pair : pairs) {
            counter++;
            Logger::get("lep_tautau::TripleSelectionAlgo")
                ->debug("Constituents hadronic tau pair {}. : {} {}", counter, pair.first,
                        pair.second);
        }
        const auto sorted_pairs = ROOT::VecOps::Sort(
            pairs, compareForTriples(selected_tau_pt, selected_tau_iso,
                                   selected_tau_pt, selected_tau_iso));

        // construct the four vectors of the selected taus to check
        // deltaR and reject a pair if the candidates are too close
        for (auto &candidate : sorted_pairs) {
            auto tau_index_1 = original_tau_indices[candidate.first];
            ROOT::Math::PtEtaPhiMVector tau_1 = ROOT::Math::PtEtaPhiMVector(
                tau_pt.at(tau_index_1), tau_eta.at(tau_index_1),
                tau_phi.at(tau_index_1), tau_mass.at(tau_index_1));
            Logger::get("lep_tautau::TripleSelectionAlgo")
                ->debug("{} leadint tau vector: {}", tau_index_1, tau_1);
            auto tau_index_2 = original_tau_indices[candidate.second];
            ROOT::Math::PtEtaPhiMVector tau_2 = ROOT::Math::PtEtaPhiMVector(
                tau_pt.at(tau_index_2), tau_eta.at(tau_index_2),
                tau_phi.at(tau_index_2), tau_mass.at(tau_index_2));
            Logger::get("lep_tautau::TripleSelectionAlgo")
                ->debug("{} tau vector: {}", tau_index_2, tau_2);
            Logger::get("lep_tautau::TripleSelectionAlgo")
                ->debug("DeltaR: {}",
                        ROOT::Math::VectorUtil::DeltaR(tau_1, tau_2));
            if (ROOT::Math::VectorUtil::DeltaR(tau_1, tau_2) > mindeltaR_tautau && ROOT::Math::VectorUtil::DeltaR(tau_1, lepton) > mindeltaR_leptau && ROOT::Math::VectorUtil::DeltaR(tau_2, lepton) > mindeltaR_leptau) {
                Logger::get("lep_tautau::TripleSelectionAlgo")
                    ->debug("Selected original pair indices: tau_1 = {} , "
                            "tau_2 = {}",
                            tau_index_1, tau_index_2);
                Logger::get("lep_tautau::TripleSelectionAlgo")
                    ->debug("Tau_1 Pt = {} , Tau_2 Pt = {} ", tau_1.Pt(),
                            tau_2.Pt());
                Logger::get("lep_tautau::TripleSelectionAlgo")
                    ->debug("Tau_1 Iso = {} , Tau_2 Iso = {} ",
                            tau_iso[tau_index_1], tau_iso[tau_index_2]);
                if (tau_1.Pt()>tau_2.Pt()) {
                   selected_triple = {static_cast<int>(original_lepton_indices.at(0)), static_cast<int>(tau_index_1), static_cast<int>(tau_index_2)}; 
                }
                else {
                    selected_triple = {static_cast<int>(original_lepton_indices.at(0)), static_cast<int>(tau_index_2), static_cast<int>(tau_index_1)};
                }
                break;
            }
        }
        Logger::get("lep_tautau::TripleSelectionAlgo")
            ->debug("Final triple {} {} {}", selected_triple[0], selected_triple[1], selected_triple[2]);

        return selected_triple;
    };
} // end namespace TripleSelectionAlgo
} // end namespace lep_tautau
namespace lep1lep1_lep2 {
/// Implementation of the triple selection algorithm. First, only events
/// that contain two good leptons lep1 and one loose lepton lep2 are considered.
/// Events contain two good lep1 and one loose lep2, if the corresponding
/// leptonmask have nonzero elements. These masks are
/// constructed using the functions from the physicsobject namespace
/// (e.g. physicsobject::CutPt). The argument triple gives information wheather the emt or the met channel is considered.
///
/// \returns an `ROOT::RVec<int>` with three values, the first one beeing
/// the lepton from the W index, the second one beeing the lepton from the tau index and the third one the hadronic tau index.
auto TripleSelectionAlgo(const float &mindeltaR_lep1lep1, const float &mindeltaR_lep1lep2) {
    Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
        ->debug("Setting up algorithm");
    return [mindeltaR_lep1lep1, mindeltaR_lep1lep2](const ROOT::RVec<float> &lep1_pt,
                       const ROOT::RVec<float> &lep1_eta,
                       const ROOT::RVec<float> &lep1_phi,
                       const ROOT::RVec<float> &lep1_mass,
                       const ROOT::RVec<float> &lep1_iso,
                       const ROOT::RVec<float> &lep2_pt,
                       const ROOT::RVec<float> &lep2_eta,
                       const ROOT::RVec<float> &lep2_phi,
                       const ROOT::RVec<float> &lep2_mass,
                       const ROOT::RVec<float> &lep2_iso,
                       const ROOT::RVec<int> &lep1_mask, const ROOT::RVec<int> &lep2_mask) {
        // first 2 entries are lep1 indices,
        // third entry is the lep2 index
        ROOT::RVec<int> selected_triple = {-1, -1, -1};
        const auto ind_good_lep1 = ROOT::VecOps::Nonzero(lep1_mask);
        const auto ind_base_lep2 = ROOT::VecOps::Nonzero(lep2_mask);
        Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
            ->debug("masks: org lep1 {}, org lep2 {}", lep1_mask, lep2_mask);

        // only select events with at least two good lep1 and at least one base lep2
        if (ind_good_lep1.size() < 2 or
            ind_base_lep2.size() < 1) {
            return selected_triple;
        }
        
        Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
            ->debug("Running algorithm on good taus and leptons");
        Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
            ->debug("good lep1 indices {}", ind_good_lep1);
        Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
            ->debug("base lep2 indices {}", ind_base_lep2);
        
        ROOT::RVec<int> selected_lep1_indices;
        ROOT::RVec<int> selected_lep2_indices;
        // sort the leptons by iso
        const auto good_lep1_iso =
            ROOT::VecOps::Take(lep1_iso, ind_good_lep1);
        const auto base_lep2_iso =
            ROOT::VecOps::Take(lep2_iso, ind_base_lep2);    
        // sort base leptons in their isolation in ascending order
        const auto sorted_good_lep1_idx = ROOT::VecOps::Argsort(good_lep1_iso);
        const auto sorted_base_lep2_idx = ROOT::VecOps::Argsort(base_lep2_iso);
        for (auto &lep1_cand : sorted_good_lep1_idx){
            Logger::get("lep1lep1_lep2::TripleSlep1ctionAlgo")
                ->debug("lep1 cand {}, sorted lep1 {}, good_lep1_iso {}", lep1_cand,sorted_good_lep1_idx, good_lep1_iso);
            selected_lep1_indices.push_back(ind_good_lep1.at(lep1_cand));
        }            
        for (auto &lep2_cand : sorted_base_lep2_idx){
            Logger::get("lep1lep1_lep2::TripleSlep2ctionAlgo")
                ->debug("lep2 cand {}, sorted lep2 {}, base_lep2_iso {}", lep2_cand,sorted_base_lep2_idx, base_lep2_iso);
            selected_lep2_indices.push_back(ind_base_lep2.at(lep2_cand));
        }

        for (int i=0; i<selected_lep1_indices.size()-1; i=i+1) {
            ROOT::Math::PtEtaPhiMVector lep1_1 = ROOT::Math::PtEtaPhiMVector(
                lep1_pt.at(selected_lep1_indices.at(i)), lep1_eta.at(selected_lep1_indices.at(i)),
                lep1_phi.at(selected_lep1_indices.at(i)), lep1_mass.at(selected_lep1_indices.at(i)));
            for (int j=i+1;j<selected_lep1_indices.size(); j=j+1) {
                ROOT::Math::PtEtaPhiMVector lep1_2 = ROOT::Math::PtEtaPhiMVector(
                    lep1_pt.at(selected_lep1_indices.at(j)), lep1_eta.at(selected_lep1_indices.at(j)),
                    lep1_phi.at(selected_lep1_indices.at(j)), lep1_mass.at(selected_lep1_indices.at(j)));
                for (int k=0;k<selected_lep2_indices.size(); k=k+1) {
                    ROOT::Math::PtEtaPhiMVector lep2 = ROOT::Math::PtEtaPhiMVector(
                    lep2_pt.at(selected_lep2_indices.at(k)), lep2_eta.at(selected_lep2_indices.at(k)),
                    lep2_phi.at(selected_lep2_indices.at(k)), lep2_mass.at(selected_lep2_indices.at(k)));
                    Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                        ->debug("DeltaR(lep1_1, lep1_2): {}",ROOT::Math::VectorUtil::DeltaR(lep1_1, lep1_2));
                    Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                        ->debug("DeltaR(lep1_1, lep2): {}", ROOT::Math::VectorUtil::DeltaR(lep1_1, lep2));
                    Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                        ->debug("DeltaR(lep1_2, lep2): {}", ROOT::Math::VectorUtil::DeltaR(lep1_2, lep2));       
                    if (ROOT::Math::VectorUtil::DeltaR(lep1_1, lep1_2) > mindeltaR_lep1lep1 && ROOT::Math::VectorUtil::DeltaR(lep1_1, lep2) > mindeltaR_lep1lep2 && ROOT::Math::VectorUtil::DeltaR(lep1_2, lep2) > mindeltaR_lep1lep2) {
                        Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                            ->debug("lep1_1 Pt = {} , lep1_2 Pt = {} ", lep1_1.Pt(),
                                    lep1_2.Pt());
                        if (lep1_1.Pt() > lep1_2.Pt()) {
                            Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                                ->debug("Selected original triple indices: lep1_1 = {}, lep1_2 = {} , lep2 = {}", selected_lep1_indices.at(i), selected_lep1_indices.at(j), selected_lep2_indices.at(k));           
                            selected_triple = {static_cast<int>(selected_lep1_indices.at(i)), static_cast<int>(selected_lep1_indices.at(j)), static_cast<int>(selected_lep2_indices.at(k))};
                            Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                                ->debug("selected triple {}", selected_triple);
                        }
                        else {
                            Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                                ->debug("Selected original triple indices: lep1_1 = {}, lep1_2 = {} , lep2 = {}", selected_lep1_indices.at(j), selected_lep1_indices.at(i), selected_lep2_indices.at(k));           
                            selected_triple = {static_cast<int>(selected_lep1_indices.at(j)), static_cast<int>(selected_lep1_indices.at(i)), static_cast<int>(selected_lep2_indices.at(k))};
                            Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                                ->debug("selected triple {}", selected_triple);
                        }
                        Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                            ->debug("lep1_1 Pt = {} , lep1_2 Pt = {} , lep2 Pt = {} ", lep1_1.Pt(),
                                    lep1_2.Pt(), lep2.Pt());
                        return selected_triple;
                    }
                }
            }
        }
        Logger::get("three_flavor::TripleSelectionAlgo")
        ->debug("no valid triple found {}", selected_triple);
        return selected_triple;
    };
} // end namespace TripleSelectionAlgo
} // end namespace lep1lep1_lep2
namespace elemutau {

/**
 * @brief Function used to select the lepton triple from the W boson and the Higgs boson decay. 
 Only triples are selected if the electron pt is larger than the pt of the muon. 
 The electron is than assigned to the W boson.
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the elemuTau triple selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - electron_pt
    - electron_eta
    - electron_phi
    - electron_mass
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - electron_masks containing the flags whether the electron is a good electron or a base electron
    - tau_mask containing the flags whether the tau is a good tau or not
    - muon_masks containing the flags whether the muon is a good muon or a base muon
 * @param triplename name of the new column containing the triple index
 * @param mindeltaR_leptau the seperation between each lepton and the tau has to be larger
 than
 * this value
 * @param mindeltaR_leplep the seperation between the leptons has to be larger
 than
 * this value
 * @return a new dataframe with the triple index column added
 */
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_leplep) {
    Logger::get("elemutau::TripleSelection")
        ->debug("Setting up EleMuTau Triple building");
    const std::string triple = "emt" ;
    auto df1 = df.Define(
        triplename,
        whtautau_tripleselection::three_flavor::TripleSelectionAlgo(mindeltaR_leptau, mindeltaR_leplep, triple),
        input_vector);
    return df1;
}

/**
 * @brief Function used to select the lepton triple from the W boson and the Higgs boson decay. 
 Only triples are selected if the electron pt is larger than the pt of the muon. 
 The electron is than assigned to the W boson.
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the elemuTau triple selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - electron_pt
    - electron_eta
    - electron_phi
    - electron_mass
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - electron_masks containing the flags whether the electron is a good electron or a base electron
    - tau_mask containing the flags whether the tau is a good tau or not
    - muon_masks containing the flags whether the muon is a good muon or a base muon
 * @param triplename name of the new column containing the triple index
 * @param mindeltaR_leptau the seperation between each lepton and the tau has to be larger
 than
 * this value
 * @param mindeltaR_leplep the seperation between the leptons has to be larger
 than
 * this value
 * @return a new dataframe with the triple index column added
 */
ROOT::RDF::RNode TripleSelectionWOEle(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau) {
    Logger::get("elemutau::TripleSelectionWOEle")
        ->debug("Setting up EleMuTau Triple building");
    auto df1 = df.Define(
        triplename,
        whtautau_tripleselection::three_flavor::TripleSelectionWithoutEleAlgo(mindeltaR_leptau),
        input_vector);
    return df1;
}
} // end namespace elemutau
namespace mueletau {

/**
 * @brief Function used to select the lepton triple from the W boson and the Higgs boson decay. 
 Only triples are selected if the electron pt is larger than the pt of the muon. 
 The electron is than assigned to the W boson.
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the mueleTau triple selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - electron_pt
    - electron_eta
    - electron_phi
    - electron_mass
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - electron_masks containing the flags whether the electron is a good electron or a base electron
    - tau_mask containing the flags whether the tau is a good tau or not
    - muon_masks containing the flags whether the muon is a good muon or a base muon
 * @param triplename name of the new column containing the triple index
 * @param mindeltaR_leptau the seperation between each lepton and the tau has to be larger
 than
 * this value
 * @param mindeltaR_leplep the seperation between the leptons has to be larger
 than
 * this value
 * @return a new dataframe with the triple index column added
 */
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_leplep) {
    Logger::get("mueletau::TripleSelection")
        ->debug("Setting up MuEleTau Triple building");
    const std::string triple = "met" ;
    auto df1 = df.Define(
        triplename,
        whtautau_tripleselection::three_flavor::TripleSelectionAlgo(mindeltaR_leptau, mindeltaR_leplep, triple),
        input_vector);
    return df1;
}
} // end namespace mueletau
namespace mumutau {

/**
 * @brief Function used to select the lepton triple from the W boson and the Higgs boson decay. 
 Only triples are selected if the electron pt is larger than the pt of the muon. 
 The electron is than assigned to the W boson.
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the mumuTau triple selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - tau_mask containing the flags whether the tau is a good tau or not
    - muon_mask containing the flags whether the muons are a good muons or not
 * @param triplename name of the new column containing the triple index
 * @param mindeltaR_leptau the seperation between each lepton and the tau has to be larger
 than
 * this value
 * @param mindeltaR_leplep the seperation between the leptons has to be larger
 than
 * this value
 * @return a new dataframe with the triple index column added
 */
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_leplep) {
    Logger::get("mumutau::TripleSelection")
        ->debug("Setting up MuMuTau Triple building");
    auto df1 = df.Define(
        triplename,
        whtautau_tripleselection::two_flavor::TripleSelectionAlgo(mindeltaR_leptau, mindeltaR_leplep),
        input_vector);
    return df1;
}
} // end namespace mumutau
namespace mumuele {

/**
 * @brief Function used to select a triple from mumu+jets events to estimate the jet to ele fakerate. 
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the mumuele triple selection these values are:
    - muon_pts
    - muon_etas
    - muon_phis
    - muon_masses
    - ele_mask containing the flags whether the ele is a good tau or not
    - muon_mask containing the flags whether the muons are a good muons or not
 * @param triplename name of the new column containing the triple index
 * @param mindeltaR_leptau the seperation between each lepton and the tau has to be larger
 than
 * this value
 * @param mindeltaR_leplep the seperation between the leptons has to be larger
 than
 * this value
 * @return a new dataframe with the triple index column added
 */
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_lep1lep1, const float &mindeltaR_lep1lep2) {
    Logger::get("mumuele::TripleSelection")
        ->debug("Setting up MuMuEle Triple building");
    auto df1 = df.Define(
        triplename,
        whtautau_tripleselection::lep1lep1_lep2::TripleSelectionAlgo(mindeltaR_lep1lep1, mindeltaR_lep1lep2),
        input_vector);
    return df1;
}
} // end namespace mumuele
namespace eleelemu {

/**
 * @brief Function used to select a triple from eleele+jets events to estimate the jet to mu fakerate. 
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the eleelemu triple selection these values are:
    - muon_pts
    - muon_etas
    - muon_phis
    - muon_masses
    - ele_mask containing the flags whether the ele is a good tau or not
    - muon_mask containing the flags whether the muons are a good muons or not
 * @param triplename name of the new column containing the triple index
 * @param mindeltaR_lep1lep1 the seperation between each lepton1  has to be larger
 than
 * this value
 * @param mindeltaR_lep1lep2 the seperation between the leptons has to be larger
 than
 * this value
 * @return a new dataframe with the triple index column added
 */
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_lep1lep1, const float &mindeltaR_lep1lep2) {
    Logger::get("eleelemu::TripleSelection")
        ->debug("Setting up EleEleMu Triple building");
    auto df1 = df.Define(
        triplename,
        whtautau_tripleselection::lep1lep1_lep2::TripleSelectionAlgo(mindeltaR_lep1lep1, mindeltaR_lep1lep2),
        input_vector);
    return df1;
}
} // end namespace eleelemu
namespace mu_tautau {

/**
 * @brief Function used to select the lepton triple from the W boson and the Higgs boson decay. 
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the mutauTau triple selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - tau_mask containing the flags whether the tau is a good tau or not
    - muon_mask containing the flags whether the muons are a good muons or not
 * @param triplename name of the new column containing the triple index
 * @param mindeltaR_leptau the seperation between each lepton and the tau has to be larger
 than
 * this value
 * @param mindeltaR_leplep the seperation between the leptons has to be larger
 than
 * this value
 * @return a new dataframe with the triple index column added
 */
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_tautau) {
    Logger::get("mu_tautau::TripleSelection")
        ->debug("Setting up mu_TauTau Triple building");
    auto df1 = df.Define(
        triplename,
        whtautau_tripleselection::lep_tautau::TripleSelectionAlgo(mindeltaR_leptau, mindeltaR_tautau),
        input_vector);
    return df1;
}
} // end namespace mutautau
namespace ele_tautau {

/**
 * @brief Function used to select the lepton triple from the W boson and the Higgs boson decay. 
 *
 * @param df the input dataframe
 * @param input_vector vector of strings containing the columns needed for the
 * alogrithm. For the eletautau triple selection these values are:
    - tau_pt
    - tau_eta
    - tau_phi
    - tau_mass
    - muon_pt
    - muon_eta
    - muon_phi
    - muon_mass
    - tau_mask containing the flags whether the tau is a good tau or not
    - muon_mask containing the flags whether the muons are a good muons or not
 * @param triplename name of the new column containing the triple index
 * @param mindeltaR_leptau the seperation between each lepton and the tau has to be larger
 than
 * this value
 * @param mindeltaR_leplep the seperation between the leptons has to be larger
 than
 * this value
 * @return a new dataframe with the triple index column added
 */
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_tautau) {
    Logger::get("ele_tautau::TripleSelection")
        ->debug("Setting up ele_TauTau Triple building");
    auto df1 = df.Define(
        triplename,
        whtautau_tripleselection::lep_tautau::TripleSelectionAlgo(mindeltaR_leptau, mindeltaR_tautau),
        input_vector);
    return df1;
}
} // end namespace eletautau
}// end namespace whtautau_tripleselection
#endif /* GUARD_TRIPLESELECTION_H */