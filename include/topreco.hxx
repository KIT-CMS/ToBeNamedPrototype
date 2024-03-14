#ifndef GUARD_TOPRECO_H
#define GUARD_TOPRECO_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

namespace topreco {

ROOT::RDF::RNode LeptonSelection(
    ROOT::RDF::RNode df, const std::string &str_n_loose_mu,
    const std::string &str_n_loose_el, const std::string &str_n_tight_mu,
    const std::string &str_n_tight_el, const std::string &str_tight_muons_mask,
    const std::string &str_tight_electrons_mask,
    const std::string &str_n_antitight_mu,
    const std::string &str_n_antitight_el,
    const std::string &str_antitight_muons_mask,
    const std::string &str_antitight_electrons_mask,
    const std::string &str_mu_pt, const std::string &str_mu_eta,
    const std::string &str_mu_phi, const std::string &str_mu_mass,
    const std::string &str_mu_charge, const std::string &str_el_pt,
    const std::string &str_el_eta, const std::string &str_el_detasc,
    const std::string &str_el_phi, const std::string &str_el_mass,
    const std::string &str_el_charge, const std::string &str_n_loose_lep,
    const std::string &str_n_tight_lep, const std::string &str_n_antitight_lep,
    const std::string &str_is_mu, const std::string &str_is_el,
    const std::string &str_is_iso, const std::string &str_lep_p4,
    const std::string &str_lep_sceta, const std::string &str_lep_charge,
    const std::string &str_mu_index, const std::string &str_el_index);

double rad_py(double x, double lep_px);
double min_fplus(double *par);
double min_fminus(double *par);
void fcn_plus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
              Int_t iflag);
void fcn_minus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
               Int_t iflag);

ROOT::RDF::RNode ReconstructLeptonicW(ROOT::RDF::RNode df,
                                      const std::string &str_lep_p4,
                                      const std::string &str_met_p4,
                                      const std::string &str_lepw_p4);

ROOT::RDF::RNode JetSelection(ROOT::RDF::RNode df, const int &njets,
                              const int &nbjets,
                              const std::string &str_good_jets_mask,
                              const std::string &str_good_bjets_mask);

ROOT::RDF::RNode
TopReco(ROOT::RDF::RNode df, const std::string &str_wlep_p4,
        const std::string &str_n_nonbjets, const std::string &str_nonbjet_p4_1,
        const std::string &str_nonbjet_btag_1,
        const std::string &str_nonbjet_p4_2,
        const std::string &str_nonbjet_btag_2, const std::string &str_n_bjets,
        const std::string &str_bjet_p4_1, const std::string &str_bjet_btag_1,
        const std::string &str_bjet_p4_2, const std::string &str_bjet_btag_2,
        const std::string &str_is_reco, const std::string &str_is_jj,
        const std::string &str_is_jjb, const std::string &str_is_jjbb,
        const std::string &str_is_jjjb, const std::string &str_is_jjjbb,
        const std::string &str_reco_p4s, const std::string &str_top_p4,
        const std::string &str_tb_p4, const std::string &str_sb_p4);

ROOT::RDF::RNode DNNQuantities(
    ROOT::RDF::RNode df, const std::string &str_is_reco,
    const std::string &str_lep_p4, const std::string &str_met_p4,
    const std::string &str_wlep_p4, const std::string &str_nonbjet_p4_1,
    const std::string &str_nonbjet_btag_1, const std::string &str_nonbjet_p4_2,
    const std::string &str_nonbjet_btag_2, const std::string &str_bjet_p4_1,
    const std::string &str_bjet_btag_1, const std::string &str_bjet_p4_2,
    const std::string &str_bjet_btag_2, const std::string &str_top_p4,
    const std::string &str_tb_p4, const std::string &str_sb_p4,
    const std::string &str_good_jetslowpt_mask, const std::string &str_jet_pt,
    const std::string &str_jet_eta, const std::string &str_jet_phi,
    const std::string &str_jet_mass, const std::string &str_dphi_top_b1,
    const std::string &str_deta_top_sb, const std::string &str_dphi_b1_b2,
    const std::string &str_deta_lep_b1, const std::string &str_m_lep_b2,
    const std::string &str_pt_b1_b2, const std::string &str_costhetastar,
    const std::string &str_sumht, const std::string &str_wolfram,
    const std::string &str_deta_topb2_b1);

ROOT::RDF::RNode LeptonScaleFactors(
    ROOT::RDF::RNode df, const std::string &str_lep_pt,
    const std::string &str_lep_eta, const std::string &str_lep_sceta,
    const std::string &str_lep_is_mu, const std::string &str_lep_is_el,
    const std::string &str_lep_is_iso,
    const std::string &str_lep_sf_mu_trigger_nom,
    const std::string &str_lep_sf_mu_trigger_up,
    const std::string &str_lep_sf_mu_trigger_down,
    const std::string &str_lep_sf_mu_iso_nom,
    const std::string &str_lep_sf_mu_iso_up,
    const std::string &str_lep_sf_mu_iso_down,
    const std::string &str_lep_sf_mu_id_nom,
    const std::string &str_lep_sf_mu_id_up,
    const std::string &str_lep_sf_mu_id_down,
    const std::string &str_lep_sf_el_trigger_nom,
    const std::string &str_lep_sf_el_trigger_up,
    const std::string &str_lep_sf_el_trigger_down,
    const std::string &str_lep_sf_el_id_nom,
    const std::string &str_lep_sf_el_id_up,
    const std::string &str_lep_sf_el_id_down,
    const std::string &str_lep_sf_el_reco_nom,
    const std::string &str_lep_sf_el_reco_up,
    const std::string &str_lep_sf_el_reco_down, const std::string &mu_sf_era,
    const std::string &mu_trigger_sf_file,
    const std::string &mu_trigger_sf_file_syst,
    const std::string &mu_trigger_sf_name,
    const std::string &mu_trigger_sf_name_syst,
    const std::string &mu_iso_sf_file, const std::string &mu_iso_sf_file_syst,
    const std::string &mu_iso_sf_name, const std::string &mu_iso_sf_name_syst,
    const std::string &mu_sf_file, const std::string &mu_id_sf_name,
    const std::string &el_sf_era, const std::string &el_trigger_sf_file,
    const std::string &el_trigger_sf_file_syst,
    const std::string &el_trigger_sf_name,
    const std::string &el_trigger_sf_name_syst, const std::string &el_sf_file,
    const std::string &el_id_sf_name);

ROOT::RDF::RNode BTagScaleFactors(
    ROOT::RDF::RNode df, const std::string &str_is_iso,
    const std::string &str_is_reco, const std::string &str_is_jj,
    const std::string &str_is_jjb, const std::string &str_is_jjbb,
    const std::string &str_is_jjjb, const std::string &str_is_jjjbb,
    const std::string &str_nonbjet_pt_1, const std::string &str_nonbjet_eta_1,
    const std::string &str_nonbjet_btag_1,
    const std::string &str_nonbjet_flavor_1,
    const std::string &str_nonbjet_pt_2, const std::string &str_nonbjet_eta_2,
    const std::string &str_nonbjet_btag_2,
    const std::string &str_nonbjet_flavor_2, const std::string &str_bjet_pt_1,
    const std::string &str_bjet_eta_1, const std::string &str_bjet_btag_1,
    const std::string &str_bjet_flavor_1, const std::string &str_bjet_pt_2,
    const std::string &str_bjet_eta_2, const std::string &str_bjet_btag_2,
    const std::string &str_bjet_flavor_2, const std::string &str_btag_sf_vec,
    const std::string &str_btagw_nom, const std::string &str_btagw_HFup_corr,
    const std::string &str_btagw_HFup_uncorr,
    const std::string &str_btagw_HFdown_corr,
    const std::string &str_btagw_HFdown_uncorr,
    const std::string &str_btagw_LFup_corr,
    const std::string &str_btagw_LFup_uncorr,
    const std::string &str_btagw_LFdown_corr,
    const std::string &str_btagw_LFdown_uncorr, const std::string &btag_sf_file,
    const std::string &btag_corr_algo_HF, const std::string &btag_corr_algo_LF,
    const std::string &btag_eff_file, const std::string &btag_eff_type,
    const std::string &btag_wp, const float &max_bjet_eta_sf);

ROOT::RDF::RNode BTagScaleFactorsGeneric(
    ROOT::RDF::RNode df, const std::string &str_jet_pt,
    const std::string &str_jet_eta, const std::string &str_jet_btag,
    const std::string &str_jet_flavor, const std::string &str_jet_collection,
    const std::string &str_btag_sf_vec, const std::string &str_btagw_nom,
    const std::string &str_btagw_HFup_corr,
    const std::string &str_btagw_HFup_uncorr,
    const std::string &str_btagw_HFdown_corr,
    const std::string &str_btagw_HFdown_uncorr,
    const std::string &str_btagw_LFup_corr,
    const std::string &str_btagw_LFup_uncorr,
    const std::string &str_btagw_LFdown_corr,
    const std::string &str_btagw_LFdown_uncorr, const std::string &btag_sf_file,
    const std::string &btag_corr_algo_HF, const std::string &btag_corr_algo_LF,
    const std::string &btag_eff_file, const std::string &btag_eff_type,
    const std::string &btag_wp, const float &btag_cut,
    const float &max_bjet_eta_sf);

ROOT::RDF::RNode
CombineDNNOutputs(ROOT::RDF::RNode df, const std::string &str_lep_is_mu,
                  const std::string &str_lep_is_el,
                  const std::string &str_lep_is_iso,
                  const std::string &str_is_reco, const std::string &str_dnn_mu,
                  const std::string &str_dnn_el, const std::string &str_dnn);

} // end namespace topreco

#endif /* GUARD_TOPRECO_H */
