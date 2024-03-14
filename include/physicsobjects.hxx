#ifndef GUARD_PHYSICSOBJECTS_H
#define GUARD_PHYSICSOBJECTS_H

namespace physicsobject {
ROOT::RDF::RNode CutPt(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const float &ptThreshold);
ROOT::RDF::RNode CutEta(ROOT::RDF::RNode df, const std::string &quantity,
                        const std::string &maskname, const float &EtaThreshold);
ROOT::RDF::RNode CutDz(ROOT::RDF::RNode df, const std::string &quantity,
                       const std::string &maskname, const float &Threshold);
ROOT::RDF::RNode CutDxy(ROOT::RDF::RNode df, const std::string &quantity,
                        const std::string &maskname, const float &Threshold);
ROOT::RDF::RNode CutVariableBarrelEndcap(
    ROOT::RDF::RNode df, const std::string &maskname,
    const std::string &etaColumnName, const std::string &cutVarColumnName,
    const float &etaBoundary, const float &lowerThresholdBarrel,
    const float &upperThresholdBarrel, const float &lowerThresholdEndcap,
    const float &upperThresholdEndcap);

/// Function to combine a list of masks into a single mask. This is done be
/// multiplying all input masks
///
/// \param[in] df the input dataframe
/// \param[out] maskname the name of the new mask to be added as column to the
/// dataframe
/// \param[in] masks a parameter pack containing an arbitrary number of
/// `std::vector<std::string>` objects. Each string is the name of a mask to be
/// combined
///
/// \return a dataframe containing the new mask
template <class... Masks>
inline ROOT::RDF::RNode CombineMasks(ROOT::RDF::RNode df,
                                     const std::string &maskname,
                                     const Masks &...masks) {
    auto multiplyMasks = [](const ROOT::RVec<ROOT::RVec<int>> &x) {
        ROOT::RVec<int> result(x[0].size(), 1);
        for (auto &xx : x) {
            result *= xx;
        }
        return result;
    };
    // std::vector<std::string> MaskList{{masks...}}; does weird things in case
    // of two arguments in masks
    std::vector<std::string> MaskList;
    utility::appendParameterPackToVector(MaskList, masks...);
    const auto nMasks = sizeof...(Masks);
    return df.Define(
        maskname, ROOT::RDF::PassAsVec<nMasks, ROOT::RVec<int>>(multiplyMasks),
        MaskList);
}
ROOT::RDF::RNode VetoCandInMask(ROOT::RDF::RNode df,
                                const std::string &outputmaskname,
                                const std::string &inputmaskname,
                                const std::string &ditaupair, const int index);
ROOT::RDF::RNode FilterMasks(ROOT::RDF::RNode df, const std::string &maskname);
ROOT::RDF::RNode LeptonVetoFlag(ROOT::RDF::RNode df,
                                const std::string &outputname,
                                const std::string &vetomap);
ROOT::RDF::RNode IsEmptyFlag(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &vetomap);
ROOT::RDF::RNode CutNFlag(ROOT::RDF::RNode df, const std::string &outputname,
                          const std::string &map, const int &n);
ROOT::RDF::RNode SelectedObjects(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &inputmaskname);
ROOT::RDF::RNode DeltaRParticleVeto(
    ROOT::RDF::RNode df, const std::string &output_flag, const std::string &p4,
    const std::string &particle_mask, const std::string &particle_pt,
    const std::string &particle_eta, const std::string &particle_phi,
    const std::string &particle_mass, const float dR_cut);
ROOT::RDF::RNode ObjectMassCorrectionWithPt(ROOT::RDF::RNode df,
                                            const std::string &corrected_mass,
                                            const std::string &raw_mass,
                                            const std::string &raw_pt,
                                            const std::string &corrected_pt);
ROOT::RDF::RNode CheckForDiLeptonPairs(
    ROOT::RDF::RNode df, const std::string &output_flag,
    const std::string &leptons_pt, const std::string &leptons_eta,
    const std::string &leptons_phi, const std::string &leptons_mass,
    const std::string &leptons_charge, const std::string &leptons_mask,
    const float dR_cut);
ROOT::RDF::RNode SelectInt(ROOT::RDF::RNode df, const std::string &maskname,
                           const std::string &nameID, const int &IDvalue);
namespace muon {
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID);
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold);
ROOT::RDF::RNode AntiCutIsolation(ROOT::RDF::RNode df,
                                  const std::string &maskname,
                                  const std::string &isolationName,
                                  const float &Threshold);
ROOT::RDF::RNode CutIsTrackerOrIsGlobal(ROOT::RDF::RNode df, const std::string &isTracker, const std::string &isGlobal, const std::string &maskname);
ROOT::RDF::RNode GenerateRndmRVec(ROOT::RDF::RNode df,
                                  const std::string &outputname,
                                  const std::string &objCollection, int seed);
ROOT::RDF::RNode
applyRoccoRData(ROOT::RDF::RNode df, const std::string &outputname,
                const std::string &filename, const int &position,
                const std::string &objCollection,
                const std::string &chargColumn, const std::string &ptColumn,
                const std::string &etaColumn, const std::string &phiColumn,
                int error_set, int error_member);
ROOT::RDF::RNode
applyRoccoRMC(ROOT::RDF::RNode df, const std::string &outputname,
              const std::string &filename, const int &position,
              const std::string &objCollection, const std::string &chargColumn,
              const std::string &ptColumn, const std::string &etaColumn,
              const std::string &phiColumn, const std::string &genPtColumn,
              const std::string &nTrackerLayersColumn,
              const std::string &rndmColumn, int error_set, int error_member);
} // namespace muon
namespace tau {
ROOT::RDF::RNode CutDecayModes(ROOT::RDF::RNode df, const std::string &maskname,
                               const std::string &tau_dms,
                               const std::vector<int> &SelectedDecayModes);
ROOT::RDF::RNode CutTauID(ROOT::RDF::RNode df, const std::string &maskname,
                          const std::string &nameID, const int &idxID);
ROOT::RDF::RNode
PtCorrection_eleFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const std::string &decayMode, const std::string &genMatch,
                     const std::string &sf_file, const std::string &jsonESname,
                     const std::string &idAlgorithm,
                     const std::string &sf_dm0_b, const std::string &sf_dm1_b,
                     const std::string &sf_dm0_e, const std::string &sf_dm1_e);
ROOT::RDF::RNode
PtCorrection_muFake(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &sf_es);
ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &decayMode,
                     const float &sf_dm0, const float &sf_dm1,
                     const float &sf_dm10, const float &sf_dm11);
ROOT::RDF::RNode
PtCorrection_genTau(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &decayMode, const std::string &genMatch,
                    const std::string &sf_file, const std::string &jsonESname,
                    const std::string &idAlgorithm, const std::string &DM0,
                    const std::string &DM1, const std::string &DM10,
                    const std::string &DM11);
} // namespace tau

namespace electron {

ROOT::RDF::RNode
PtCorrection_byValue(ROOT::RDF::RNode df, const std::string &corrected_pt,
                     const std::string &pt, const std::string &eta,
                     const float &sf_barrel, const float &sf_endcap);
ROOT::RDF::RNode PtCorrection(ROOT::RDF::RNode df, const std::string &corrected_pt,
                    const std::string &pt, const std::string &eta,
                    const std::string &sf_barrel, const std::string &sf_endcap,
                    const std::string &sf_file, const std::string &jsonESname);
ROOT::RDF::RNode CutID(ROOT::RDF::RNode df, const std::string &maskname,
                       const std::string &nameID);
ROOT::RDF::RNode CutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                         const std::string &nameID, const int &IDvalue);
ROOT::RDF::RNode AntiCutCBID(ROOT::RDF::RNode df, const std::string &maskname,
                             const std::string &nameID, const int &IDvalue);
ROOT::RDF::RNode CutIsolation(ROOT::RDF::RNode df, const std::string &maskname,
                              const std::string &isolationName,
                              const float &Threshold);
ROOT::RDF::RNode AntiCutIsolation(ROOT::RDF::RNode df,
                                  const std::string &maskname,
                                  const std::string &isolationName,
                                  const float &Threshold);
ROOT::RDF::RNode CutIP(ROOT::RDF::RNode df, const std::string &eta,
                       const std::string &detasc, const std::string &dxy,
                       const std::string &dz, const std::string &maskname,
                       const float &abseta_eb_ee, const float &max_dxy_eb,
                       const float &max_dz_eb, const float &max_dxy_ee,
                       const float &max_dz_ee);

ROOT::RDF::RNode CutGap(ROOT::RDF::RNode df, const std::string &eta,
                        const std::string &detasc, const std::string &maskname,
                        const float &end_eb, const float &start_ee);

ROOT::RDF::RNode CutCBIDBitmapNoIso(ROOT::RDF::RNode df, const std::string &maskname, const int &IDvalue, const std::string &bitmap);

} // end namespace electron
} // namespace physicsobject
#endif /* GUARD_PHYSICSOBJECTS_H */
