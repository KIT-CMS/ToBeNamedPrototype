#ifndef GUARDFAKEFACTORS_H
#define GUARDFAKEFACTORS_H
/// The namespace that contains the fake factor function.
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "correction.h"

namespace fakefactors {
/**
 * @brief Function to calculate raw fake factors without corrections with
 * correctionlib for the semileptonic channels
 *
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param nbtags number of good b-tagged jets in the event
 * @param variation name of the uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
raw_fakefactor_nmssm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &tau_pt, const std::string &njets,
                        const std::string &lep_mt, const std::string &nbtags,
                        const std::string &variation,
                        const std::string &ff_file) {
    Logger::get("RawFakeFactor")
        ->debug("Setting up functions for raw fake factor (without "
                "corrections) evaluation with correctionlib");
    Logger::get("RawFakeFactor")->debug("Variation - Name {}", variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");
    auto calc_fake_factor = [variation, qcd, wjets, ttbar,
                             fractions](const float &pt_2, const int &njets,
                                        const float &mt_1, const int &nbtag) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("RawFakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff = wjets->evaluate({pt_2, (float)njets, variation});
            Logger::get("RawFakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("RawFakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("RawFakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("RawFakeFactor")->debug("N b-jets - value {}", nbtag);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)nbtag, variation});
            Logger::get("RawFakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)nbtag, variation});
            Logger::get("RawFakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)nbtag, variation});
            Logger::get("RawFakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);

            ff = qcd_frac * qcd_ff + wjets_frac * wjets_ff +
                 ttbar_frac * ttbar_ff;
        }

        Logger::get("RawFakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor,
                         {tau_pt, njets, lep_mt, nbtags});
    return df1;
}
/**
 * @brief Function to calculate raw fake factors without corrections with
 * correctionlib for the full hadronic channel
 *
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_idx index of the tau, leading/subleading
 * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
 * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param variation name of the uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
raw_fakefactor_nmssm_tt(ROOT::RDF::RNode df, const std::string &outputname,
                        const int &tau_idx, const std::string &tau_pt_1,
                        const std::string &tau_pt_2, const std::string &njets,
                        const std::string &variation,
                        const std::string &ff_file) {

    Logger::get("RawFakeFactor")
        ->debug("Setting up functions for raw fake factor (without "
                "corrections) evaluation with correctionlib");
    Logger::get("RawFakeFactor")->debug("Variation - Name {}", variation);

    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");

    auto qcd_subleading = correction::CorrectionSet::from_file(ff_file)->at(
        "QCD_subleading_fake_factors");

    auto calc_fake_factor = [tau_idx, variation, qcd, qcd_subleading](
                                const float &pt_1, const float &pt_2,
                                const int &njets) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("RawFakeFactor")
                ->debug("Leading Tau pt - value {}", pt_1);
            Logger::get("RawFakeFactor")
                ->debug("Subleading Tau pt - value {}", pt_2);
            Logger::get("RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = -1.;
            if (tau_idx == 0) {
                float qcd_ff = qcd->evaluate({pt_1, (float)njets, variation});
                Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
                ff = qcd_ff;
            } else if (tau_idx == 1) {
                float qcd_ff =
                    qcd_subleading->evaluate({pt_2, (float)njets, variation});
                Logger::get("RawFakeFactor")->debug("QCD - value {}", qcd_ff);
                ff = qcd_ff;
            }
        }

        Logger::get("RawFakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 =
        df.Define(outputname, calc_fake_factor, {tau_pt_1, tau_pt_2, njets});
    return df1;
}
/**
 * @brief Function to calculate fake factors with correctionlib for the
 * semileptonic channels
 *
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param nbtags number of good b-tagged jets in the event
 * @param lep_pt pt of the leptonic tau in the tau pair
 * @param lep_iso isolation of the leptonic tau in the tau pair
 * @param m_vis visible mass of the tau pair
 * @param variation name of the uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
fakefactor_nmssm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tau_pt, const std::string &njets,
                    const std::string &lep_mt, const std::string &nbtags,
                    const std::string &lep_pt, const std::string &lep_iso,
                    const std::string &m_vis, const std::string &variation,
                    const std::string &ff_file,
                    const std::string &ff_corr_file) {

    Logger::get("FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("FakeFactor")->debug("Variation - Name {}", variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto qcd_lep_pt_closure = correction::CorrectionSet::from_file(ff_corr_file)
                                  ->at("QCD_non_closure_lep_pt_correction");
    auto qcd_lep_iso_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_lep_iso_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");
    auto wjets_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("Wjets_non_closure_lep_pt_correction");
    auto wjets_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                           ->at("Wjets_DR_SR_correction");
    auto ttbar_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_lep_pt_correction");
    auto calc_fake_factor = [variation, qcd, wjets, ttbar, fractions,
                             qcd_lep_pt_closure, qcd_lep_iso_closure, qcd_DR_SR,
                             wjets_lep_pt_closure, wjets_DR_SR,
                             ttbar_lep_pt_closure](
                                const float &pt_2, const int &njets,
                                const float &mt_1, const int &nbtag,
                                const float &pt_1, const float &iso_1,
                                const float &m_vis) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff = wjets->evaluate({pt_2, (float)njets, variation});
            Logger::get("FakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("FakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("FakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("FakeFactor")->debug("N b-jets - value {}", nbtag);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)nbtag, variation});
            Logger::get("FakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)nbtag, variation});
            Logger::get("FakeFactor")->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)nbtag, variation});
            Logger::get("FakeFactor")->debug("ttbar - fraction {}", ttbar_frac);

            Logger::get("FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("FakeFactor")->debug("m_vis - value {}", m_vis);

            float qcd_lep_pt_corr =
                qcd_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("FakeFactor")
                ->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            float qcd_lep_iso_corr =
                qcd_lep_iso_closure->evaluate({iso_1, variation});
            Logger::get("FakeFactor")
                ->debug("QCD - lep iso correction {}", qcd_lep_iso_corr);
            float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, variation});
            Logger::get("FakeFactor")
                ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
            float wjets_lep_pt_corr =
                wjets_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            float wjets_DR_SR_corr = wjets_DR_SR->evaluate({m_vis, variation});
            Logger::get("FakeFactor")
                ->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);
            float ttbar_lep_pt_corr =
                ttbar_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("FakeFactor")
                ->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);

            ff = qcd_frac * qcd_ff * qcd_lep_pt_corr * qcd_lep_iso_corr *
                     qcd_DR_SR_corr +
                 wjets_frac * wjets_ff * wjets_lep_pt_corr * wjets_DR_SR_corr +
                 ttbar_frac * ttbar_ff * ttbar_lep_pt_corr;
        }

        Logger::get("FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 =
        df.Define(outputname, calc_fake_factor,
                  {tau_pt, njets, lep_mt, nbtags, lep_pt, lep_iso, m_vis});
    return df1;
}
/**
 * @brief Function to calculate fake factors with correctionlib for the full
 * hadronic channel
 *
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_idx index of the tau, leading/subleading
 * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
 * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param m_vis visible mass of the tau pair
 * @param variation name of the uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
fakefactor_nmssm_tt(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &tau_idx, const std::string &tau_pt_1,
                    const std::string &tau_pt_2, const std::string &njets,
                    const std::string &m_vis, const std::string &variation,
                    const std::string &ff_file,
                    const std::string &ff_corr_file) {

    Logger::get("FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("FakeFactor")->debug("Variation - Name {}", variation);

    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");

    auto qcd_tau_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_subleading_lep_pt_correction");
    auto qcd_m_vis_closure = correction::CorrectionSet::from_file(ff_corr_file)
                                 ->at("QCD_non_closure_m_vis_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");

    auto qcd_subleading = correction::CorrectionSet::from_file(ff_file)->at(
        "QCD_subleading_fake_factors");

    auto qcd_tau_pt_closure_subleading =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_subleading_non_closure_leading_lep_pt_correction");
    auto qcd_m_vis_closure_subleading =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_subleading_non_closure_m_vis_correction");
    auto qcd_DR_SR_subleading =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_subleading_DR_SR_correction");

    auto calc_fake_factor = [tau_idx, variation, qcd, qcd_tau_pt_closure,
                             qcd_m_vis_closure, qcd_DR_SR, qcd_subleading,
                             qcd_tau_pt_closure_subleading,
                             qcd_m_vis_closure_subleading,
                             qcd_DR_SR_subleading](
                                const float &pt_1, const float &pt_2,
                                const int &njets, const float &m_vis) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("FakeFactor")->debug("Leading Tau pt - value {}", pt_1);
            Logger::get("FakeFactor")
                ->debug("Subleading Tau pt - value {}", pt_2);
            Logger::get("FakeFactor")->debug("m_vis - value {}", m_vis);
            Logger::get("FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = -1.;
            float qcd_tau_pt_corr = -1.;
            float qcd_m_vis_corr = -1.;
            float qcd_DR_SR_corr = -1.;
            if (tau_idx == 0) {
                float qcd_ff = qcd->evaluate({pt_1, (float)njets, variation});
                Logger::get("FakeFactor")->debug("QCD - value {}", qcd_ff);
                float qcd_tau_pt_corr =
                    qcd_tau_pt_closure->evaluate({pt_2, variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - lep pt correction {}", qcd_tau_pt_corr);
                float qcd_m_vis_corr =
                    qcd_m_vis_closure->evaluate({m_vis, variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - visible mass correction {}", qcd_m_vis_corr);
                float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, variation});
                Logger::get("FakeFactor")
                    ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
                ff = qcd_ff * qcd_tau_pt_corr * qcd_m_vis_corr * qcd_DR_SR_corr;
            } else if (tau_idx == 1) {
                float qcd_ff =
                    qcd_subleading->evaluate({pt_2, (float)njets, variation});
                Logger::get("FakeFactor")
                    ->debug("QCD(subleading) - value {}", qcd_ff);
                float qcd_tau_pt_corr =
                    qcd_tau_pt_closure_subleading->evaluate({pt_1, variation});
                Logger::get("FakeFactor")
                    ->debug("QCD(subleading) - lep pt correction {}",
                            qcd_tau_pt_corr);
                float qcd_m_vis_corr =
                    qcd_m_vis_closure_subleading->evaluate({m_vis, variation});
                Logger::get("FakeFactor")
                    ->debug("QCD(subleading) - visible mass correction {}",
                            qcd_m_vis_corr);
                float qcd_DR_SR_corr =
                    qcd_DR_SR_subleading->evaluate({m_vis, variation});
                Logger::get("FakeFactor")
                    ->debug("QCD(subleading) - DR to SR correction {}",
                            qcd_DR_SR_corr);
                ff = qcd_ff * qcd_tau_pt_corr * qcd_m_vis_corr * qcd_DR_SR_corr;
            }
        }

        Logger::get("FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor,
                         {tau_pt_1, tau_pt_2, njets, m_vis});
    return df1;
}
/**
 * @brief Function to calculate raw fake factors without corrections with
 * correctionlib. In difference to the NMSSM version, njets is used for the
 * fraction binning, and an additional split in deltaR is applied for wjets.
 *
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param delta_r delta R between the two taus
 * @param variation name of the uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
raw_fakefactor_sm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &tau_pt, const std::string &njets,
                     const std::string &lep_mt, const std::string &delta_r,
                     const std::string &variation, const std::string &ff_file) {
    Logger::get("SM RawFakeFactor")
        ->debug("Setting up functions for raw fake factor (without "
                "corrections) evaluation with correctionlib");
    Logger::get("SM RawFakeFactor")->debug("Variation - Name {}", variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");
    auto calc_fake_factor = [variation, qcd, wjets, ttbar, fractions](
                                const float &pt_2, const int &njets,
                                const float &mt_1, const float &delta_r) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM RawFakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM RawFakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff =
                wjets->evaluate({pt_2, (float)njets, delta_r, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("ttbar - value {}", ttbar_ff);

            Logger::get("SM RawFakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("SM RawFakeFactor")->debug("N jets - value {}", njets);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)njets, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)njets, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)njets, variation});
            Logger::get("SM RawFakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);

            ff = qcd_frac * qcd_ff + wjets_frac * wjets_ff +
                 ttbar_frac * ttbar_ff;
        }

        Logger::get("SM RawFakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor,
                         {tau_pt, njets, lep_mt, delta_r});
    return df1;
}
/**
 * @brief Function to calculate fake factors with correctionlib. In difference
 * to the NMSSM version, njets is used for the fraction binning, and an
 * additional split in deltaR is applied for wjets.
 * @param df the input dataframe
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param lep_pt pt of the leptonic tau in the tau pair
 * @param lep_iso isolation of the leptonic tau in the tau pair
 * @param m_vis visible mass of the tau pair
 * @param delta_r distance in eta-phi between the two taus
 * @param variation name of the uncertainty variation or nominal
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
fakefactor_sm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                 const std::string &tau_pt, const std::string &njets,
                 const std::string &lep_mt, const std::string &lep_pt,
                 const std::string &lep_iso, const std::string &m_vis,
                 const std::string &delta_r, const std::string &variation,
                 const std::string &ff_file, const std::string &ff_corr_file) {
    Logger::get("SM FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("SM FakeFactor")->debug("Variation - Name {}", variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto qcd_lep_pt_closure = correction::CorrectionSet::from_file(ff_corr_file)
                                  ->at("QCD_non_closure_lep_pt_correction");
    auto qcd_lep_iso_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_lep_iso_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");
    auto wjets_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("Wjets_non_closure_lep_pt_correction");
    auto wjets_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                           ->at("Wjets_DR_SR_correction");
    auto ttbar_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_lep_pt_correction");
    auto calc_fake_factor = [variation, qcd, wjets, ttbar, fractions,
                             qcd_lep_pt_closure, qcd_lep_iso_closure, qcd_DR_SR,
                             wjets_lep_pt_closure, wjets_DR_SR,
                             ttbar_lep_pt_closure](
                                const float &pt_2, const int &njets,
                                const float &mt_1, const float &pt_1,
                                const float &iso_1, const float &m_vis,
                                const float &delta_r) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("QCD - value {}", qcd_ff);
            float wjets_ff =
                wjets->evaluate({pt_2, (float)njets, delta_r, variation});
            Logger::get("SM FakeFactor")->debug("Wjets - value {}", wjets_ff);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("SM FakeFactor")->debug("Lep mt - value {}", mt_1);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);

            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("QCD - fraction {}", qcd_frac);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);

            Logger::get("SM FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FakeFactor")->debug("m_vis - value {}", m_vis);

            float qcd_lep_pt_corr =
                qcd_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            float qcd_lep_iso_corr =
                qcd_lep_iso_closure->evaluate({iso_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - lep iso correction {}", qcd_lep_iso_corr);
            float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
            float wjets_lep_pt_corr =
                wjets_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            float wjets_DR_SR_corr = wjets_DR_SR->evaluate({m_vis, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);
            float ttbar_lep_pt_corr =
                ttbar_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);

            ff = qcd_frac * qcd_ff * qcd_lep_pt_corr * qcd_lep_iso_corr *
                     qcd_DR_SR_corr +
                 wjets_frac * wjets_ff * wjets_lep_pt_corr * wjets_DR_SR_corr +
                 ttbar_frac * ttbar_ff * ttbar_lep_pt_corr;
        }

        Logger::get("SM FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 =
        df.Define(outputname, calc_fake_factor,
                  {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis, delta_r});
    return df1;
}

/// Function to calculate raw fake factors without corrections with
/// correctionlib. In difference to the NMSSM version, njets is used for the
/// fraction binning, and an additional split in deltaR is applied for wjets.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the output column for the fake factor
/// \param qcd_outputname name of the output column for the qcd fake factor
/// \param wjets_outputname name of the output column for the wjets fake factor
/// \param ttbar_outputname name of the output column for the ttbar fake factor
/// \param tau_pt pt of the hadronic tau in the tau pair
/// \param njets number of good jets in the event
/// \param lep_mt transverse mass of the leptonic tau in the tau pair
/// \param delta_r delta R between the two taus
/// \param variation name of the uncertainty variation or nominal
/// \param ff_file correctionlib json file with the fake factors
///
/// \returns a dataframe with the fake factors
ROOT::RDF::RNode
raw_fakefactor_sm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &qcd_outputname,
                     const std::string &wjets_outputname,
                     const std::string &ttbar_outputname,
                     const std::string &tau_pt, const std::string &njets,
                     const std::string &lep_mt, const std::string &delta_r,
                     const std::string &variation, const std::string &ff_file) {
    Logger::get("SM RawFakeFactor")
        ->debug("Setting up functions for raw fake factor (without "
                "corrections) evaluation with correctionlib");
    Logger::get("SM RawFakeFactor")->debug("Variation - Name {}", variation);
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto calc_raw_fake_factor_qcd =
        [variation, qcd, fractions](const float &pt_2, const int &njets,
                                    const float &mt_1) {
            float ff_qcd = -1.;
            if (pt_2 >= 0.) {
                float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
                float qcd_frac =
                    fractions->evaluate({"QCD", mt_1, (float)njets, variation});

                ff_qcd = qcd_frac * qcd_ff;
            }
            return ff_qcd;
        };
    auto calc_raw_fake_factor_wjets =
        [variation, wjets, fractions](const float &pt_2, const int &njets,
                                      const float &mt_1, const float &delta_r) {
            float ff_wjets = -1.;
            if (pt_2 >= 0.) {
                float wjets_ff =
                    wjets->evaluate({pt_2, (float)njets, delta_r, variation});
                float wjets_frac = fractions->evaluate(
                    {"Wjets", mt_1, (float)njets, variation});
                ff_wjets = wjets_frac * wjets_ff;
            }
            return ff_wjets;
        };
    auto calc_raw_fake_factor_ttbar = [variation, ttbar, fractions](
                                          const float &pt_2, const int &njets,
                                          const float &mt_1) {
        float ff_ttbar = -1.;
        if (pt_2 >= 0.) {
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)njets, variation});

            ff_ttbar = ttbar_frac * ttbar_ff;
        }
        return ff_ttbar;
    };
    auto calc_full_raw_fake_factor = [](const float &qcd, const float &wjets,
                                        const float &ttbar) {
        float ff = qcd + wjets + ttbar;
        Logger::get("SM FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(qcd_outputname, calc_raw_fake_factor_qcd,
                         {tau_pt, njets, lep_mt});
    auto df2 = df1.Define(wjets_outputname, calc_raw_fake_factor_wjets,
                          {tau_pt, njets, lep_mt, delta_r});
    auto df3 = df2.Define(ttbar_outputname, calc_raw_fake_factor_ttbar,
                          {
                              tau_pt,
                              njets,
                              lep_mt,
                          });
    auto df4 = df3.Define(outputname, calc_full_raw_fake_factor,
                          {qcd_outputname, wjets_outputname, ttbar_outputname});
    return df4;
}

/// Function to calculate fake factors with correctionlib.
/// In difference to the NMSSM version, njets is used for the
/// fraction binning, and an additional split in deltaR is
/// applied for wjets.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the output column for the fake factor
/// \param qcd_outputname name of the output column for the qcd fake factor
/// \param wjets_outputname name of the output column for the wjets fake factor
/// \param ttbar_outputname name of the output column for the ttbar fake factor
/// \param tau_pt pt of the hadronic tau in the tau pair
/// \param njets number of good jets in the event
/// \param lep_mt transverse mass of the leptonic tau in the tau pair
/// \param lep_pt pt of the leptonic tau in the tau pair
/// \param lep_iso isolation of the leptonic tau in the tau pair
/// \param m_vis visible mass of the tau pair
/// \param delta_r distance in eta-phi between the two taus
/// \param variation name of the uncertainty variation or nominal
/// \param ff_file correctionlib json file with the fake factors
/// \param ff_corr_file correctionlib json file with corrections for the fake
/// factors
///
/// \returns a dataframe with the fake factors
ROOT::RDF::RNode
fakefactor_sm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                 const std::string &qcd_outputname,
                 const std::string &wjets_outputname,
                 const std::string &ttbar_outputname, const std::string &tau_pt,
                 const std::string &njets, const std::string &lep_mt,
                 const std::string &lep_pt, const std::string &lep_iso,
                 const std::string &m_vis, const std::string &delta_r,
                 const std::string &variation, const std::string &ff_file,
                 const std::string &ff_corr_file) {
    Logger::get("SM FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("SM FakeFactor")->debug("Variation - Name {}", variation);
    // load the corrections from the correctionlib json file
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto qcd_lep_pt_closure = correction::CorrectionSet::from_file(ff_corr_file)
                                  ->at("QCD_non_closure_lep_pt_correction");
    auto qcd_lep_iso_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_lep_iso_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");
    auto wjets_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("Wjets_non_closure_lep_pt_correction");
    auto wjets_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                           ->at("Wjets_DR_SR_correction");
    auto ttbar_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_lep_pt_correction");
    // setup lambda functions to calculate the fake factors split by origin
    auto calc_fake_factor_qcd = [variation, qcd, fractions, qcd_lep_pt_closure,
                                 qcd_lep_iso_closure, qcd_DR_SR](
                                    const float &pt_2, const int &njets,
                                    const float &mt_1, const float &pt_1,
                                    const float &iso_1, const float &m_vis,
                                    const float &delta_r) {
        float ff_qcd = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("QCD - value {}", qcd_ff);
            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("QCD - fraction {}", qcd_frac);
            Logger::get("SM FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FakeFactor")->debug("m_vis - value {}", m_vis);

            float qcd_lep_pt_corr =
                qcd_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            float qcd_lep_iso_corr =
                qcd_lep_iso_closure->evaluate({iso_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - lep iso correction {}", qcd_lep_iso_corr);
            float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
            ff_qcd = qcd_frac * qcd_ff * qcd_lep_pt_corr * qcd_lep_iso_corr * qcd_DR_SR_corr;
        }

        Logger::get("SM FakeFactor")
            ->debug("Event Fake Factor for QCD {}", ff_qcd);
        return ff_qcd;
    };
    auto calc_fake_factor_wjets = [variation, wjets, fractions,
                                   wjets_lep_pt_closure, wjets_DR_SR,
                                   ttbar_lep_pt_closure](
                                      const float &pt_2, const int &njets,
                                      const float &mt_1, const float &pt_1,
                                      const float &iso_1, const float &m_vis,
                                      const float &delta_r) {
        float ff_wjets = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);
            float wjets_ff =
                wjets->evaluate({pt_2, (float)njets, delta_r, variation});
            Logger::get("SM FakeFactor")->debug("Wjets - value {}", wjets_ff);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            Logger::get("SM FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FakeFactor")->debug("m_vis - value {}", m_vis);

            float wjets_lep_pt_corr =
                wjets_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            float wjets_DR_SR_corr = wjets_DR_SR->evaluate({m_vis, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);

            ff_wjets =
                wjets_frac * wjets_ff * wjets_lep_pt_corr * wjets_DR_SR_corr;
        }

        Logger::get("SM FakeFactor")
            ->debug("Event Fake Factor for Wjets {}", ff_wjets);
        return ff_wjets;
    };
    auto calc_fake_factor_ttbar = [variation, ttbar, fractions,
                                   ttbar_lep_pt_closure](
                                      const float &pt_2, const int &njets,
                                      const float &mt_1, const float &pt_1,
                                      const float &iso_1, const float &m_vis,
                                      const float &delta_r) {
        float ff_ttbar = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("ttbar - value {}", ttbar_ff);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);
            Logger::get("SM FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FakeFactor")->debug("m_vis - value {}", m_vis);
            float ttbar_lep_pt_corr =
                ttbar_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);

            ff_ttbar = ttbar_frac * ttbar_ff * ttbar_lep_pt_corr;
        }

        Logger::get("SM FakeFactor")
            ->debug("Event Fake Factor fot ttbar {}", ff_ttbar);
        return ff_ttbar;
    };
    auto calc_full_fake_factor = [](const float &qcd, const float &wjets,
                                    const float &ttbar) {
        float ff = qcd + wjets + ttbar;
        Logger::get("SM FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 =
        df.Define(qcd_outputname, calc_fake_factor_qcd,
                  {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis, delta_r});
    auto df2 =
        df1.Define(wjets_outputname, calc_fake_factor_wjets,
                   {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis, delta_r});
    auto df3 =
        df2.Define(ttbar_outputname, calc_fake_factor_ttbar,
                   {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis, delta_r});
    auto df4 = df3.Define(outputname, calc_full_fake_factor,
                          {qcd_outputname, wjets_outputname, ttbar_outputname});
    return df4;
}

/// Function to calculate fake factors with correctionlib.
/// In difference to the NMSSM version, njets is used for the
/// fraction binning, and an additional split in deltaR is
/// applied for wjets.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the output column for the fake factor
/// \param tau_pt pt of the hadronic tau in the tau pair
/// \param njets number of good jets in the event
/// \param lep_mt transverse mass of the leptonic tau in the tau pair
/// \param lep_pt pt of the leptonic tau in the tau pair
/// \param lep_iso isolation of the leptonic tau in the tau pair
/// \param m_vis visible mass of the tau pair
/// \param delta_r distance in eta-phi between the two taus
/// \param variation name of the uncertainty variation or nominal
/// \param ff_file correctionlib json file with the fake factors
/// \param ff_corr_file correctionlib json file with corrections for the fake
/// factors
///
/// \returns a dataframe with the fake factors
ROOT::RDF::RNode fakefactor_sm_lt_no_deltaR(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &qcd_outputname, const std::string &wjets_outputname,
    const std::string &ttbar_outputname, const std::string &tau_pt,
    const std::string &njets, const std::string &lep_mt,
    const std::string &lep_pt, const std::string &lep_iso,
    const std::string &m_vis, const std::string &variation,
    const std::string &ff_file, const std::string &ff_corr_file) {
    Logger::get("SM FakeFactor")
        ->debug("Setting up functions for fake factor evaluation with "
                "correctionlib");
    Logger::get("SM FakeFactor")->debug("Variation - Name {}", variation);
    // load the corrections from the correctionlib json file
    auto qcd =
        correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    auto wjets =
        correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    auto ttbar =
        correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    auto fractions =
        correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto qcd_lep_pt_closure = correction::CorrectionSet::from_file(ff_corr_file)
                                  ->at("QCD_non_closure_lep_pt_correction");
    auto qcd_lep_iso_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("QCD_non_closure_lep_iso_correction");
    auto qcd_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                         ->at("QCD_DR_SR_correction");
    auto wjets_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("Wjets_non_closure_lep_pt_correction");
    auto wjets_DR_SR = correction::CorrectionSet::from_file(ff_corr_file)
                           ->at("Wjets_DR_SR_correction");
    auto ttbar_lep_pt_closure =
        correction::CorrectionSet::from_file(ff_corr_file)
            ->at("ttbar_non_closure_lep_pt_correction");
    // setup lambda functions to calculate the fake factors split by origin
    auto calc_fake_factor_qcd = [variation, qcd, fractions, qcd_lep_pt_closure,
                                 qcd_lep_iso_closure, qcd_DR_SR](
                                    const float &pt_2, const int &njets,
                                    const float &mt_1, const float &pt_1,
                                    const float &iso_1, const float &m_vis) {
        float ff_qcd = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("QCD - value {}", qcd_ff);
            float qcd_frac =
                fractions->evaluate({"QCD", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("QCD - fraction {}", qcd_frac);
            Logger::get("SM FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FakeFactor")->debug("m_vis - value {}", m_vis);

            float qcd_lep_pt_corr =
                qcd_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            float qcd_lep_iso_corr =
                qcd_lep_iso_closure->evaluate({iso_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - lep iso correction {}", qcd_lep_iso_corr);
            float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, variation});
            Logger::get("SM FakeFactor")
                ->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
            ff_qcd = qcd_frac * qcd_ff * qcd_lep_pt_corr * qcd_lep_iso_corr *
                     qcd_DR_SR_corr;
        }

        Logger::get("SM FakeFactor")
            ->debug("Event Fake Factor for QCD {}", ff_qcd);
        return ff_qcd;
    };
    auto calc_fake_factor_wjets = [variation, wjets, fractions,
                                   wjets_lep_pt_closure, wjets_DR_SR,
                                   ttbar_lep_pt_closure](
                                      const float &pt_2, const int &njets,
                                      const float &mt_1, const float &pt_1,
                                      const float &iso_1, const float &m_vis) {
        float ff_wjets = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);
            float wjets_ff = wjets->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("Wjets - value {}", wjets_ff);
            float wjets_frac =
                fractions->evaluate({"Wjets", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - fraction {}", wjets_frac);
            Logger::get("SM FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FakeFactor")->debug("m_vis - value {}", m_vis);

            float wjets_lep_pt_corr =
                wjets_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            float wjets_DR_SR_corr = wjets_DR_SR->evaluate({m_vis, variation});
            Logger::get("SM FakeFactor")
                ->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);

            ff_wjets =
                wjets_frac * wjets_ff * wjets_lep_pt_corr * wjets_DR_SR_corr;
        }

        Logger::get("SM FakeFactor")
            ->debug("Event Fake Factor for Wjets {}", ff_wjets);
        return ff_wjets;
    };
    auto calc_fake_factor_ttbar = [variation, ttbar, fractions,
                                   ttbar_lep_pt_closure](
                                      const float &pt_2, const int &njets,
                                      const float &mt_1, const float &pt_1,
                                      const float &iso_1, const float &m_vis) {
        float ff_ttbar = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor")->debug("N jets - value {}", njets);
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, variation});
            Logger::get("SM FakeFactor")->debug("ttbar - value {}", ttbar_ff);
            float ttbar_frac =
                fractions->evaluate({"ttbar", mt_1, (float)njets, variation});
            Logger::get("SM FakeFactor")
                ->debug("ttbar - fraction {}", ttbar_frac);
            Logger::get("SM FakeFactor")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FakeFactor")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FakeFactor")->debug("m_vis - value {}", m_vis);
            float ttbar_lep_pt_corr =
                ttbar_lep_pt_closure->evaluate({pt_1, variation});
            Logger::get("SM FakeFactor")
                ->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);

            ff_ttbar = ttbar_frac * ttbar_ff * ttbar_lep_pt_corr;
        }

        Logger::get("SM FakeFactor")
            ->debug("Event Fake Factor fot ttbar {}", ff_ttbar);
        return ff_ttbar;
    };
    auto calc_full_fake_factor = [](const float &qcd, const float &wjets,
                                    const float &ttbar) {
        float ff = qcd + wjets + ttbar;
        Logger::get("SM FakeFactor")->debug("Event Fake Factor {}", ff);
        return ff;
    };
    auto df1 = df.Define(qcd_outputname, calc_fake_factor_qcd,
                         {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis});
    auto df2 = df1.Define(wjets_outputname, calc_fake_factor_wjets,
                          {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis});
    auto df3 = df2.Define(ttbar_outputname, calc_fake_factor_ttbar,
                          {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis});
    auto df4 = df3.Define(outputname, calc_full_fake_factor,
                          {qcd_outputname, wjets_outputname, ttbar_outputname});
    return df4;
}

} // namespace fakefactors
#endif /* GUARDFAKEFACTORS_H */