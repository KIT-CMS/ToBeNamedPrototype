#ifndef GUARD_SCALEFACTORS_H
#define GUARD_SCALEFACTORS_H

namespace scalefactor {
namespace muon {
ROOT::RDF::RNode id_rooworkspace(ROOT::RDF::RNode df, const std::string &pt,
                                 const std::string &eta,
                                 const std::string &id_output,
                                 const std::string &workspace_name,
                                 const std::string &id_functor_name,
                                 const std::string &id_arguments);
ROOT::RDF::RNode iso_rooworkspace(ROOT::RDF::RNode df, const std::string &pt,
                                  const std::string &eta,
                                  const std::string &iso,
                                  const std::string &iso_output,
                                  const std::string &workspace_name,
                                  const std::string &iso_functor_name,
                                  const std::string &iso_arguments);
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &pt,
                    const std::string &eta, const std::string &year_id,
                    const std::string &variation, const std::string &id_output,
                    const std::string &sf_file, const std::string &idAlgorithm);
ROOT::RDF::RNode iso(ROOT::RDF::RNode df, const std::string &pt,
                     const std::string &eta, const std::string &year_id,
                     const std::string &variation,
                     const std::string &iso_output, const std::string &sf_file,
                     const std::string &idAlgorithm);
} // namespace muon
namespace tau {

ROOT::RDF::RNode
id_vsJet_lt(ROOT::RDF::RNode df, const std::string &pt,
            const std::string &decayMode, const std::string &genMatch,
            const std::vector<int> &selectedDMs, const std::string &wp,
            const std::string &sf_vsjet_tau30to35,
            const std::string &sf_vsjet_tau35to40,
            const std::string &sf_vsjet_tau40to500,
            const std::string &sf_vsjet_tau500to1000,
            const std::string &sf_vsjet_tau1000toinf,
            const std::string &sf_dependence, const std::string &id_output,
            const std::string &sf_file, const std::string &idAlgorithm);
ROOT::RDF::RNode id_vsJet_lt_embedding(
    ROOT::RDF::RNode df, const std::string &pt, const std::string &wp,
    const std::string &sf_vsjet_tau20to25,
    const std::string &sf_vsjet_tau25to30,
    const std::string &sf_vsjet_tau30to35,
    const std::string &sf_vsjet_tau35to40,
    const std::string &sf_vsjet_tau40toInf, const std::string &id_output,
    const std::string &sf_file, const std::string &correctionset);
ROOT::RDF::RNode id_vsJet_tt_embedding(
    ROOT::RDF::RNode df, const std::string &decaymode, const std::string &wp,
    const std::string &sf_vsjet_tauDM0, const std::string &sf_vsjet_tauDM1,
    const std::string &sf_vsjet_tauDM10, const std::string &sf_vsjet_tauDM11,
    const std::string &id_output, const std::string &sf_file,
    const std::string &correctionset);
ROOT::RDF::RNode id_vsJet_tt(
    ROOT::RDF::RNode df, const std::string &pt, const std::string &decayMode,
    const std::string &genMatch, const std::vector<int> &selectedDMs,
    const std::string &wp, const std::string &sf_vsjet_tauDM0,
    const std::string &sf_vsjet_tauDM1, const std::string &sf_vsjet_tauDM10,
    const std::string &sf_vsjet_tauDM11, const std::string &sf_dependence,
    const std::string &id_output, const std::string &sf_file,
    const std::string &idAlgorithm);
ROOT::RDF::RNode
id_vsEle(ROOT::RDF::RNode df, const std::string &eta,
         const std::string &decayMode, const std::string &genMatch,
         const std::vector<int> &selectedDMs, const std::string &wp,
         const std::string &sf_vsele_barrel, const std::string &sf_vsele_endcap,
         const std::string &id_output, const std::string &sf_file,
         const std::string &idAlgorithm);
ROOT::RDF::RNode
id_vsMu(ROOT::RDF::RNode df, const std::string &eta,
        const std::string &decayMode, const std::string &genMatch,
        const std::vector<int> &selectedDMs, const std::string &wp,
        const std::string &sf_vsmu_wheel1, const std::string &sf_vsmu_wheel2,
        const std::string &sf_vsmu_wheel3, const std::string &sf_vsmu_wheel4,
        const std::string &sf_vsmu_wheel5, const std::string &id_output,
        const std::string &sf_file, const std::string &idAlgorithm);
} // namespace tau

namespace electron {

ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &pt,
                    const std::string &eta, const std::string &year_id,
                    const std::string &wp, const std::string &variation,
                    const std::string &id_output, const std::string &sf_file,
                    const std::string &idAlgorithm);
} // namespace electron

namespace jet {

ROOT::RDF::RNode
btagSF(ROOT::RDF::RNode df, const std::string &pt, const std::string &eta,
       const std::string &btag_discr, const std::string &flavor,
       const std::string &jet_mask, const std::string &bjet_mask,
       const std::string &jet_veto_mask, const std::string &variation,
       const std::string &sf_output, const std::string &sf_file,
       const std::string &corr_algorithm);
} // namespace jet
namespace embedding {
ROOT::RDF::RNode
selection_trigger(ROOT::RDF::RNode df, const std::string &pt_1,
                  const std::string &eta_1, const std::string &pt_2,
                  const std::string &eta_2, const std::string &output,
                  const std::string &sf_file, const std::string &idAlgorithm);
ROOT::RDF::RNode
embedding_cleaning_correction(ROOT::RDF::RNode df, const std::string &puppi_met,
                              const std::string &output, const std::string &sf_file,
                              const std::string &idAlgorithm);
ROOT::RDF::RNode selection_id(ROOT::RDF::RNode df, const std::string &pt,
                              const std::string &eta, const std::string &output,
                              const std::string &sf_file,
                              const std::string &idAlgorithm);
ROOT::RDF::RNode muon_sf(ROOT::RDF::RNode df, const std::string &pt,
                         const std::string &eta, const std::string &output,
                         const std::string &sf_file,
                         const std::string correctiontype,
                         const std::string &idAlgorithm,
                         const float &extrapolation_factor = 1.0);
ROOT::RDF::RNode electron_sf(ROOT::RDF::RNode df, const std::string &pt,
                             const std::string &eta, const std::string &output,
                             const std::string &sf_file,
                             const std::string correctiontype,
                             const std::string &idAlgorithm,
                             const float &extrapolation_factor = 1.0);
} // namespace embedding
} // namespace scalefactor
#endif /* GUARD_SCALEFACTORS_H */