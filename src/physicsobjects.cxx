#ifndef GUARD_PHYSICSOBJECTS_H
#define GUARD_PHYSICSOBJECTS_H

#include "../include/RoccoR.hxx"
#include "../include/basefunctions.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "TRandom3.h"
#include "correction.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <iostream>
#include <string>
#include <type_traits>
#include <vector>
/// Namespace containing function to apply cuts on physics objects. The
/// cut results are typically stored within a mask, which is represented by
/// an `ROOT::RVec<int>`.
///    \code
///    In the mask
///    1 --> cut is passed by the object
///    0 --> cut is not passed by the object
///    \endcode
/// multiple cuts can be combined by multiplying masks using
/// physicsobject::CombineMasks.
namespace physicsobject {
/// Function to select objects above a pt threshold, using
/// basefunctions::FilterMin
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the pt column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] ptThreshold minimal pt value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutPt(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const float &ptThreshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMin(ptThreshold), {quantity});
    return df1;
}
/// Function to select objects blow an eta threshold, using
/// basefunctions::FilterAbsMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the eta column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] EtaThreshold maximal eta value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutEta(ROOT::RDF::RNode df, const std::string &quantity,
                        const std::string &maskname,
                        const float &EtaThreshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterAbsMax(EtaThreshold),
                         {quantity});
    return df1;
}
/// Function to select objects below an Dz threshold, using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the Dz column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] Threshold maximal Dz value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutDz(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const float &Threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterAbsMax(Threshold), {quantity});
    return df1;
}
/// Function to select objects below an Dxy threshold, using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] quantity name of the Dxy column in the NanoAOD
/// \param[out] maskname the name of the mask to be added as column to the
/// dataframe \param[in] Threshold maximal Dxy value
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutDxy(ROOT::RDF::RNode df, const std::string &quantity,
                        const std::string &maskname, const float &Threshold) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterAbsMax(Threshold), {quantity});
    return df1;
}
/// Function to select objects with eta dependent upper and lower thesholds
/// for a given variable
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the mask to be added as column to the
/// \param[in] etaColumnName name of the eta column in the NanoAOD dataframe
/// \param[in] cutVarColumnName name of the variable column to apply the
/// selection in the NanoAOD dataframe \param[in] etaBoundary boundary of
/// absolute eta for the barrel and endcap regions of the detector \param[in]
/// lowerThresholdBarrel lower threshold for the barrel \param[in]
/// upperThresholdBarrel upper threshold for the barrel \param[in]
/// lowerThresholdEndcap lower threshold for the endcap \param[in]
/// upperThresholdEndcap upper threshold for the barrel
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutVariableBarrelEndcap(
    ROOT::RDF::RNode df, const std::string &maskname,
    const std::string &etaColumnName, const std::string &cutVarColumnName,
    const float &etaBoundary, const float &lowerThresholdBarrel,
    const float &upperThresholdBarrel, const float &lowerThresholdEndcap,
    const float &upperThresholdEndcap) {
    auto lambda = [etaBoundary, lowerThresholdBarrel, upperThresholdBarrel,
                   lowerThresholdEndcap,
                   upperThresholdEndcap](const ROOT::RVec<float> &eta,
                                         const ROOT::RVec<float> &variable) {
        ROOT::RVec<int> mask =
            (((abs(eta) < etaBoundary) && (variable >= lowerThresholdBarrel) &&
              (variable < upperThresholdBarrel)) ||
             ((abs(eta) >= etaBoundary) && (variable >= lowerThresholdEndcap) &&
              (variable < upperThresholdEndcap)));
        return mask;
    };

    auto df1 = df.Define(maskname, lambda, {etaColumnName, cutVarColumnName});
    return df1;
}

/// Function to take a mask and create a new one where a particle candidate is
/// set to false
///
/// \param[in] df the input dataframe
/// \param[out] outputmaskname the name of the new mask to be added as column to
/// the dataframe \param[in] inputmaskname the name of the input mask \param[in]
/// dileptonpair name of the column of the dileptonpair \param[in] index index
/// of the particle candidate to be ignored by mask
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode VetoCandInMask(ROOT::RDF::RNode df,
                                const std::string &outputmaskname,
                                const std::string &inputmaskname,
                                const std::string &dileptonpair,
                                const int index) {
    return df.Define(outputmaskname,
                     [index, inputmaskname](const ROOT::RVec<int> &mask,
                                            const ROOT::RVec<int> &pair) {
                        Logger::get("VetoCandInMask")
                            ->debug("Vetoing the selected candidate (index "
                                    "{}) from the mask {}",
                                    index, inputmaskname);
                        auto newmask = mask;
                        if (pair.at(index) >= 0)
                            newmask.at(pair.at(index)) = 0;
                        Logger::get("VetoCandInMask")
                             ->debug("Vetoing the selected candidate (index "
                                     "{}) from the mask {}, new mask {}",
                                     index, inputmaskname, newmask);    
                        return newmask;
                     },
                     {inputmaskname, dileptonpair});
}

/// Function to filter events based on a mask. If the mask contains at least
/// one object, the event is filtered
///
///   \param df the input dataframe
///   \param maskname the name of the column containing the vetomap
///    to be used
///
///   \return a new df with the events filtered
ROOT::RDF::RNode FilterMasks(ROOT::RDF::RNode df, const std::string &maskname) {
    auto df1 = df.Filter(
        [](const ROOT::RVec<Int_t> &mask) {
            auto result = Any(mask);
            return result;
        },
        {maskname});
    return df1;
}

/// Function to generate a veto based on a mask. If the mask contains at least
/// one object, the veto flag is set, meaning that this event contains at least
/// one object matching the requirements of the veto map.
///
/// \code
///  In the veto column
///  1 --> the event contains at least one object matching the veto map
///  0 --> the event does not contain any object matching the veto map
/// \endcode
///
/// \param df the input dataframe
/// \param outputname name of the new quantity containing the veto flags
/// \param vetomap the name of the column containing the vetomap to be used
///
/// \return a new df containing the veto flag column
ROOT::RDF::RNode LeptonVetoFlag(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &vetomap) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &mask) {
                         return ROOT::VecOps::Nonzero(mask).size() != 0;
                     },
                     {vetomap});
}

/// Function to create a boolian flag based on the number of non-zero masks.
/// The flag is set to true if the number of non-zero masks is zero.
///
/// \param df the input dataframe
/// \param outputname the name of the output column that is created
/// \param vetomap the name of the column containing the input mask to be used
///
/// \return a new df containing the output flag column
ROOT::RDF::RNode IsEmptyFlag(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &vetomap) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &mask) {
                         return ROOT::VecOps::Nonzero(mask).size() == 0;
                     },
                     {vetomap});
}

/// Function to create a boolian flag based on the number of non-zero masks.
/// The flag is set to true if the number of non-zero masks matches the
/// allowed value provided as an input.
///
/// \param df the input dataframe
/// \param outputname the name of the output column that is created
/// \param map the name of the column containing the input mask to be used
/// \param n the allowed number of non-zero masks
///
/// \return a new df containing the output flag column
ROOT::RDF::RNode CutNFlag(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &map, const int &n) {
    return df.Define(outputname,
                     [n](const ROOT::RVec<int> &mask) {
                         return ROOT::VecOps::Nonzero(mask).size() == n;
                     },
                     {map});
}

/// Function to create a column for a vector of indices of objects based on
/// input masks. Indices of objects with non-zero mask are stored in the output
/// column.
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] inputmaskname the name of input masks
///
/// \return a dataframe containing a vector of indices for the selected objects
ROOT::RDF::RNode SelectedObjects(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &inputmaskname) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &mask) {
                         Logger::get("SelectedObjects")
                             ->debug("size = {}",
                                     ROOT::VecOps::Nonzero(mask).size());
                         return static_cast<ROOT::VecOps::RVec<int>>(
                             ROOT::VecOps::Nonzero(mask));
                     },
                     {inputmaskname});
}

/**
 * @brief Function used to veto a particle, if it is overlapping within a given
 * DeltaR value with another particle within a mask. For the particle to test,
 * the lorentz vector is used, for the veto, the mask and all input quantities
 * needed to calculate the lorentz vector are used. Returns true if the particle
 * is vetoed, false if it is not.
 *
 * @param df The input dataframe
 * @param output_flag The name of the veto flag to be added as column to the
 * dataframe
 * @param p4 The name of the Lorentz vector column to be used for the particle
 * to test
 * @param particle_mask The name of the mask to be used for the veto
 * @param particle_pt The name of the pt column to be used for the particles to
 * test against
 * @param particle_eta The name of the eta column to be used for the particles
 * to test against
 * @param particle_phi The name of the phi column to be used for the particles
 * to test against
 * @param particle_mass The name of the mass column to be used for the particles
 * to test against
 * @param dR_cut The maximum dR to be used for the veto
 * @return ROOT::RDF::RNode a dataframe containing the new veto flag
 */

ROOT::RDF::RNode DeltaRParticleVeto(
    ROOT::RDF::RNode df, const std::string &output_flag, const std::string &p4,
    const std::string &particle_mask, const std::string &particle_pt,
    const std::string &particle_eta, const std::string &particle_phi,
    const std::string &particle_mass, const float dR_cut) {
    auto veto_overlapping_particle =
        [dR_cut](const ROOT::Math::PtEtaPhiMVector &p4,
                 const ROOT::RVec<float> &particle_pt,
                 const ROOT::RVec<float> &particle_eta,
                 const ROOT::RVec<float> &particle_phi,
                 const ROOT::RVec<float> &particle_mass,
                 const ROOT::RVec<int> &particle_mask) {
            // for all particles in the mask, check if they overlap with the
            // particle, if so, return true
            ROOT::RVec<int> valid_particle_indices =
                ROOT::VecOps::Nonzero(particle_mask);
            const auto selected_pt =
                ROOT::VecOps::Take(particle_pt, valid_particle_indices);
            const auto selected_eta =
                ROOT::VecOps::Take(particle_eta, valid_particle_indices);
            const auto selected_phi =
                ROOT::VecOps::Take(particle_phi, valid_particle_indices);
            const auto selected_mass =
                ROOT::VecOps::Take(particle_mass, valid_particle_indices);
            auto selected_p4s =
                ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(
                    selected_pt, selected_eta, selected_phi, selected_mass);
            for (const auto &p4_test : selected_p4s) {
                if (ROOT::Math::VectorUtil::DeltaR(p4_test, p4) < dR_cut) {
                    return true;
                }
            }
            // if no particle is close enough to the p4, return false
            return false;
        };
    auto df1 = df.Define(output_flag, veto_overlapping_particle,
                         {p4, particle_pt, particle_eta, particle_phi,
                          particle_mass, particle_mask});
    return df1;
}

/// Function to correct object mass in alignment with object pt correction
///
/// \param[in] df the input dataframe
/// \param[out] corrected_mass the name of the corrected masses to be determined
/// \param[in] raw_mass name of the input mass \param[in] raw_pt name of the
/// uncorrected object pts \param[in] corrected_pt name of the corrected object
/// pts
///
/// \return a dataframe containing the modified object masses
ROOT::RDF::RNode ObjectMassCorrectionWithPt(ROOT::RDF::RNode df,
                                            const std::string &corrected_mass,
                                            const std::string &raw_mass,
                                            const std::string &raw_pt,
                                            const std::string &corrected_pt) {
    auto mass_correction_lambda =
        [](const ROOT::RVec<float> &mass_values,
           const ROOT::RVec<float> &pt_values,
           const ROOT::RVec<float> &corrected_pt_values) {
            ROOT::RVec<float> corrected_mass_values(mass_values.size());
            for (int i = 0; i < mass_values.size(); i++) {
                corrected_mass_values[i] = mass_values.at(i) *
                                           corrected_pt_values.at(i) /
                                           pt_values.at(i);
            }
            return corrected_mass_values;
        };
    auto df1 = df.Define(corrected_mass, mass_correction_lambda,
                         {raw_mass, raw_pt, corrected_pt});
    return df1;
}

/// Function to check whether at least one lepton pair is present
///
/// \param[in] df the input dataframe
/// \param[out] output_flag the name of the bool column that is created
/// \param[in] leptons_pt name of the input pt column of the lepton collection
/// \param[in] leptons_eta name of the input eta column of the lepton collection
/// \param[in] leptons_phi name of the input phi column of the lepton collection
/// \param[in] leptons_mass name of the input mass column of the lepton
/// collection \param[in] leptons_charge name of the input charge column of the
/// lepton collection \param[in] leptons_mask name of the input mask column of
/// the lepton collection that marks lepton to be taken into account \param[in]
/// dR_cut minimum required angular distance between the leptons
///
/// \return a dataframe containing the new bool column
ROOT::RDF::RNode CheckForDiLeptonPairs(
    ROOT::RDF::RNode df, const std::string &output_flag,
    const std::string &leptons_pt, const std::string &leptons_eta,
    const std::string &leptons_phi, const std::string &leptons_mass,
    const std::string &leptons_charge, const std::string &leptons_mask,
    const float dR_cut) {
    auto pair_finder_lambda = [dR_cut](const ROOT::RVec<float> &pt_values,
                                       const ROOT::RVec<float> &eta_values,
                                       const ROOT::RVec<float> &phi_values,
                                       const ROOT::RVec<float> &mass_values,
                                       const ROOT::RVec<int> &charge_values,
                                       const ROOT::RVec<int> &mask) {
        const auto valid_lepton_indices = ROOT::VecOps::Nonzero(mask);
        for (auto it1 = valid_lepton_indices.begin();
             it1 != valid_lepton_indices.end(); it1++) {
            for (auto it2 = it1 + 1; it2 != valid_lepton_indices.end(); it2++) {
                if (charge_values.at(*it1) != charge_values.at(*it2)) {
                    auto p4_1 = ROOT::Math::PtEtaPhiMVector(
                        pt_values.at(*it1), eta_values.at(*it1),
                        phi_values.at(*it1), mass_values.at(*it1));
                    auto p4_2 = ROOT::Math::PtEtaPhiMVector(
                        pt_values.at(*it2), eta_values.at(*it2),
                        phi_values.at(*it2), mass_values.at(*it2));
                    if (ROOT::Math::VectorUtil::DeltaR(p4_1, p4_2) >= dR_cut)
                        return true;
                }
            }
        }
        return false;
    };
    auto df1 = df.Define(output_flag, pair_finder_lambda,
                         {leptons_pt, leptons_eta, leptons_phi, leptons_mass,
                          leptons_charge, leptons_mask});
    return df1;
}

/// Muon specific functions
namespace muon {
/// Function to cut on muons based on the muon ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
}
/// Function to cut on muons based on the muon isolation using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold),
                         {isolationName});
    return df1;
}
/// Function to cut on muons based on the muon signature: isTracker or isGlobal
///
/// \param[in] df the input dataframe
/// \param[in] isTracker name of the signature column in the NanoAOD
/// \param[in] isGlobal name of the signature column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe 
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsTrackerOrIsGlobal(ROOT::RDF::RNode df, const std::string &isTracker, const std::string &isGlobal, const std::string &maskname) {
        auto lambda = [](const ROOT::RVec<Bool_t>  &tracker, const ROOT::RVec<Bool_t>  &global) {
            ROOT::RVec<int> mask = (tracker == 1 || global == 1);
            Logger::get("lep1lep1_lep2::TripleSelectionAlgo")
                ->debug("istracker {}, isglobal {}, mask {}", tracker, global, mask);
            return mask;
        };
        auto df1 = df.Define(maskname, lambda, {isTracker, isGlobal});
        return df1;
}
/// Function to create a column of vector of random numbers between 0 and 1
/// with size of the input object collection
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] objCollection the name of the input object collection
/// \param[in] seed the seed of the random number generator
///
/// \return a dataframe with the new column
ROOT::RDF::RNode GenerateRndmRVec(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &objCollection, int seed) {
    gRandom->SetSeed(seed);
    auto lambda = [](const ROOT::RVec<int> &objects) {
        const int len = objects.size();
        float rndm[len];
        gRandom->RndmArray(len, rndm);
        ROOT::RVec<float> out = {};
        for (auto &x : rndm) {
            out.push_back(x);
        }
        return out;
    };
    return df.Define(outputname, lambda, {objCollection});
}

/// Function to create a column of Rochester correction applied transverse
/// momentum for data, see https://gitlab.cern.ch/akhukhun/roccor
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] filename the name of Rochester correction file
/// \param[in] position index of the position in the input vector
/// \param[in] objCollection the name of the column containing the input vector
/// \param[in] chargColumn the name of the column containing muon charges
/// \param[in] ptColumn the name of the column containing muon charge values
/// \param[in] etaColumn the name of the column containing muon eta values
/// \param[in] phiColumn the name of the column containing muon phi values
/// \param[in] error_set the error set number
/// \param[in] error_member the error member number
///
/// \return a dataframe with the new column
ROOT::RDF::RNode
applyRoccoRData(ROOT::RDF::RNode df, const std::string &outputname,
                const std::string &filename, const int &position,
                const std::string &objCollection,
                const std::string &chargColumn, const std::string &ptColumn,
                const std::string &etaColumn, const std::string &phiColumn,
                int error_set, int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set,
                   error_member](const ROOT::RVec<int> &objects,
                                 const ROOT::RVec<int> &chargCol,
                                 const ROOT::RVec<float> &ptCol,
                                 const ROOT::RVec<float> &etaCol,
                                 const ROOT::RVec<float> &phiCol) {
        const int index = objects.at(position);
        double pt_rc =
            ptCol.at(index) * rc.kScaleDT(chargCol.at(index), ptCol.at(index),
                                          etaCol.at(index), phiCol.at(index),
                                          error_set, error_member);
        return pt_rc;
    };

    return df.Define(
        outputname, lambda,
        {objCollection, chargColumn, ptColumn, etaColumn, phiColumn});
}

/// Function to create a column of Rochester correction applied transverse
/// momentum for MC, see https://gitlab.cern.ch/akhukhun/roccor
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the output column that is created
/// \param[in] filename the name of Rochester correction file
/// \param[in] position index of the position in the input vector
/// \param[in] objCollection the name of the column containing the input vector
/// \param[in] chargColumn the name of the column containing muon charges
/// \param[in] ptColumn the name of the column containing muon charge values
/// \param[in] etaColumn the name of the column containing muon eta values
/// \param[in] phiColumn the name of the column containing muon phi values
/// \param[in] genPtColumn the name of the column containing gen-level
/// transverse momentum value of the target muon \param[in] nTrackerLayersColumn
/// the name of the column containing number of tracker layers values \param[in]
/// rndmColumn the name of the column containing random number generated for
/// each muon \param[in] error_set the error set number \param[in] error_member
/// the error member number
///
/// \return a dataframe with the new column
ROOT::RDF::RNode
applyRoccoRMC(ROOT::RDF::RNode df, const std::string &outputname,
              const std::string &filename, const int &position,
              const std::string &objCollection, const std::string &chargColumn,
              const std::string &ptColumn, const std::string &etaColumn,
              const std::string &phiColumn, const std::string &genPtColumn,
              const std::string &nTrackerLayersColumn,
              const std::string &rndmColumn, int error_set, int error_member) {
    RoccoR rc(filename);
    auto lambda = [rc, position, error_set, error_member](
                      const ROOT::RVec<int> &objects,
                      const ROOT::RVec<int> &chargCol,
                      const ROOT::RVec<float> &ptCol,
                      const ROOT::RVec<float> &etaCol,
                      const ROOT::RVec<float> &phiCol, const float &genPt,
                      const ROOT::RVec<int> &nTrackerLayersCol,
                      const ROOT::RVec<float> &rndmCol) {
        double pt_rc = default_float;
        const int index = objects.at(position);
        if (genPt > 0.) {
            pt_rc = ptCol.at(index) *
                    rc.kSpreadMC(chargCol.at(index), ptCol.at(index),
                                 etaCol.at(index), phiCol.at(index), genPt,
                                 error_set, error_member);
        } else {
            pt_rc = ptCol.at(index) *
                    rc.kSmearMC(chargCol.at(index), ptCol.at(index),
                                etaCol.at(index), phiCol.at(index),
                                nTrackerLayersCol.at(index),
                                rndmCol.at(position), error_set, error_member);
        }

        return pt_rc;
    };

    return df.Define(outputname, lambda,
                     {objCollection, chargColumn, ptColumn, etaColumn,
                      phiColumn, genPtColumn, nTrackerLayersColumn,
                      rndmColumn});
}
} // end namespace muon
/// Tau specific functions
namespace tau {
/// Function to cut on taus based on the tau decay mode
///
/// \param[in] df the input dataframe
/// \param[in] tau_dms name of the column with tau decay modes
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] SelectedDecayModes a `std::vector<int>` containing the
/// decay modes, that should pass the cut
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutDecayModes(ROOT::RDF::RNode df, const std::string &maskname,
                               const std::string &tau_dms,
                               const std::vector<int> &SelectedDecayModes) {
    auto df1 = df.Define(
        maskname,
        [SelectedDecayModes](const ROOT::RVec<Int_t> &decaymodes) {
            ROOT::RVec<int> mask;
            for (auto n : decaymodes) {
                mask.push_back(int(std::find(SelectedDecayModes.begin(),
                                             SelectedDecayModes.end(),
                                             n) != SelectedDecayModes.end()));
            }
            return mask;
        },
        {tau_dms});
    return df1;
}
/// Function to cut taus based on the tau ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// idxID bitvalue of the WP the has to be passed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutTauID(ROOT::RDF::RNode df, const std::string &maskname,
                          const std::string &nameID, const int &idxID) {
    auto df1 = df.Define(maskname, basefunctions::FilterID(idxID), {nameID});
    return df1;
}
/// Function to correct e to tau fake pt
///
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] df the input dataframe
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] sf_dm0_b scale factor to be applied to taus with decay mode 0 and
/// eta region barrel \param[in] sf_dm1_b scale factor to be applied to taus
/// with decay mode 1 and eta region barrel \param[in] sf_dm0_e scale factor to
/// be applied to taus with decay mode 0 and eta region endcap \param[in]
/// sf_dm1_e scale factor to be applied to taus with decay mode 1 and eta region
/// endcap name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_eleFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const std::string &decayMode, const std::string &genMatch,
                     const std::string &sf_file, const std::string &jsonESname,
                     const std::string &idAlgorithm,
                     const std::string &sf_dm0_b, const std::string &sf_dm1_b,
                     const std::string &sf_dm0_e, const std::string &sf_dm1_e) {
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda = [evaluator, idAlgorithm, sf_dm0_b, sf_dm1_b,
                                     sf_dm0_e, sf_dm1_e](
                                        const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (genmatch.at(i) == 1 || genmatch.at(i) == 3) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0 &&
                    std::abs(eta_values.at(i)) <= 1.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm0_b});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 0 &&
                           std::abs(eta_values.at(i)) > 1.5 &&
                           std::abs(eta_values.at(i)) <= 2.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm0_e});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(eta_values.at(i)) <= 1.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm1_b});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1 &&
                           std::abs(eta_values.at(i)) > 1.5 &&
                           std::abs(eta_values.at(i)) <= 2.5) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_dm1_e});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                }
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
            Logger::get("ptcorrection ele fake")
                ->debug("tau pt before {}, tau pt after {}", pt_values.at(i),
                        corrected_pt_values.at(i));
        }
        return corrected_pt_values;
    };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct mu to tau fake pt
///
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] df the input dataframe
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] sf_es scale factor to be applied to taus faked by muons
/// name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_muFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &sf_es) {
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda =
        [evaluator, idAlgorithm, sf_es](const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (genmatch.at(i) == 2 || genmatch.at(i) == 4) {
                    // only considering wanted tau decay modes
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, sf_es});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                }
                if (genmatch.at(i) == 2 || genmatch.at(i) == 4) {
                    Logger::get("mu fake")->debug(
                        "tau pt before {}, tau pt after {}", pt_values.at(i),
                        corrected_pt_values.at(i));
                }
            }
            return corrected_pt_values;
        };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
/// Function to correct tau pt
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt \param[in] decayMode
/// \param[in] sf_dm0 scale factor to be applied to taus with decay mode 0
/// \param[in] sf_dm1 scale factor to be applied to other 1 prong taus
/// \param[in] sf_dm10 scale factor to be applied to taus with decay mode 10
/// \param[in] sf_dm11 scale factor to be applied to other 3 prong taus
/// name of the tau decay mode quantity
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &decayMode,
                     const float &sf_dm0, const float &sf_dm1,
                     const float &sf_dm10, const float &sf_dm11) {
    auto tau_pt_correction_lambda =
        [sf_dm0, sf_dm1, sf_dm10, sf_dm11](const ROOT::RVec<float> &pt_values,
                                           const ROOT::RVec<int> &decay_modes) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (decay_modes.at(i) == 0)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm0;
                else if (decay_modes.at(i) > 0 && decay_modes.at(i) < 5)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm1;
                else if (decay_modes.at(i) == 10)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm10;
                else if (decay_modes.at(i) > 10 && decay_modes.at(i) < 15)
                    corrected_pt_values[i] = pt_values.at(i) * sf_dm11;
                else
                    corrected_pt_values[i] = pt_values.at(i);
            }
            return corrected_pt_values;
        };
    auto df1 =
        df.Define(corrected_pt, tau_pt_correction_lambda, {pt, decayMode});
    return df1;
}
/// Function to correct tau pt with correctionlib
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected tau pt to be calculated
/// \param[in] pt name of the raw tau pt
/// \param[in] eta name of raw tau eta
/// \param[in] decayMode decay mode of the tau
/// \param[in] genMatch column with genmatch values (from prompt e, prompt mu,
/// tau->e, tau->mu, had. tau) \param[in] sf_file:
///     2018:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2018_UL.html
///     2017:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2017_UL.html
///     2016:
///     https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016preVFP_UL.html
///           https://cms-nanoaod-integration.web.cern.ch/commonJSONSFs/TAU_tau_Run2_UL/TAU_tau_2016postVFP_UL.html
/// \param[in] jsonESname name of the tau energy correction in the json file
/// \param[in] idAlgorithm name of the used tau id algorithm
/// \param[in] DM0 variation decay mode 0
/// \param[in] DM1 variation decay mode 1
/// \param[in] DM10 variation decay mode 10
/// \param[in] DM11 variation decay mode 11
/// DM values: "nom","up","down"
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_genTau(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &DM0,
                    const std::string &DM1, const std::string &DM10,
                    const std::string &DM11) {
    auto evaluator =
        correction::CorrectionSet::from_file(sf_file)->at(jsonESname);
    auto tau_pt_correction_lambda = [evaluator, idAlgorithm, DM0, DM1, DM10,
                                     DM11](
                                        const ROOT::RVec<float> &pt_values,
                                        const ROOT::RVec<float> &eta_values,
                                        const ROOT::RVec<int> &decay_modes,
                                        const ROOT::RVec<UChar_t> &genmatch) {
        ROOT::RVec<float> corrected_pt_values(pt_values.size());
        for (int i = 0; i < pt_values.size(); i++) {
            if (genmatch.at(i) == 5) {
                // only considering wanted tau decay modes
                if (decay_modes.at(i) == 0) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM0});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 1) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM1});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 10) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM10});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                } else if (decay_modes.at(i) == 11) {
                    auto sf = evaluator->evaluate(
                        {pt_values.at(i), std::abs(eta_values.at(i)),
                         decay_modes.at(i), static_cast<int>(genmatch.at(i)),
                         idAlgorithm, DM11});
                    corrected_pt_values[i] = pt_values.at(i) * sf;
                }
            } else {
                corrected_pt_values[i] = pt_values.at(i);
            }
            Logger::get("tauEnergyCorrection")
                ->warn("tau pt before {}, tau pt after {}, decaymode {}",
                        pt_values.at(i), corrected_pt_values.at(i),
                        decay_modes.at(i));
        }
        return corrected_pt_values;
    };
    auto df1 = df.Define(corrected_pt, tau_pt_correction_lambda,
                         {pt, eta, decayMode, genMatch});
    return df1;
}
} // end namespace tau

namespace electron {
/// Function to correct electron pt
///
/// \param[in] df the input dataframe
/// \param[out] corrected_pt name of the corrected electron pt to be calculated
/// \param[in] pt name of the raw electron pt
/// \param[in] eta the name of the raw electron eta
/// \param[in] sf_barrel scale factor to be applied to electrons in the barrel
/// \param[in] sf_endcap scale factor to be applied to electrons in the endcap
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const float &sf_barrel, const float &sf_endcap) {
    auto electron_pt_correction_lambda =
        [sf_barrel, sf_endcap](const ROOT::RVec<float> &pt_values,
                               const ROOT::RVec<float> &eta) {
            ROOT::RVec<float> corrected_pt_values(pt_values.size());
            for (int i = 0; i < pt_values.size(); i++) {
                if (abs(eta.at(i)) <= 1.479) {
                    corrected_pt_values[i] = pt_values.at(i) * sf_barrel;
                } else if (abs(eta.at(i)) > 1.479) {
                    corrected_pt_values[i] = pt_values.at(i) * sf_endcap;
                } else {
                    corrected_pt_values[i] = pt_values.at(i);
                }
            }
            return corrected_pt_values;
        };
    auto df1 =
        df.Define(corrected_pt, electron_pt_correction_lambda, {pt, eta});
    return df1;
}
/// Function to cut electrons based on the electron MVA ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID) {
    auto df1 = df.Define(
        maskname,
        [](const ROOT::RVec<Bool_t> &id) { return (ROOT::RVec<int>)id; },
        {nameID});
    return df1;
} /// Function to cut jets based on the cut based electron ID
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] nameID name of the ID column in the NanoAOD \param[in]
/// IDvalue value of the WP the has to be passed
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                         const std::string &nameID, const int &IDvalue) {
    auto df1 =
        df.Define(maskname, basefunctions::FilterMinInt(IDvalue), {nameID});
    return df1;
}
/// Function to cut electrons based on the electron isolation using
/// basefunctions::FilterMax
///
/// \param[in] df the input dataframe
/// \param[in] isolationName name of the isolation column in the NanoAOD
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe \param[in] Threshold maximal isolation threshold
///
/// \return a dataframe containing the new mask
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold) {
    auto df1 = df.Define(maskname, basefunctions::FilterMax(Threshold),
                         {isolationName});
    return df1;
}

} // end namespace electron

} // namespace physicsobject
#endif /* GUARD_PHYSICSOBJECTS_H */