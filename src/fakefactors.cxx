#ifndef GUARDFAKEFACTORS_H
#define GUARDFAKEFACTORS_H
/// The namespace that contains the fake factor function.
#include "../include/utility/CorrectionManager.hxx"
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
 * @param correctionManager The CorrectionManager object
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param delta_r delta R between the two taus
 * @param fraction_variation name of the fraction variation
 * @param qcd_variation name of the QCD variation
 * @param wjets_variation name of the Wjets variation
 * @param ttbar_variation name of the ttbar variation
 * @param ff_file correctionlib json file with the fake factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
raw_fakefactor_sm_lt(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correctionManager,
    const std::string &outputname,
    const std::string &tau_pt,
    const std::string &njets,
    const std::string &lep_mt,
    const std::string &delta_r,
    //
    const std::string &fraction_variation,
    const std::string &qcd_variation,
    const std::string &wjets_variation,
    const std::string &ttbar_variation,
    //
    const std::string &ff_file
) {

    Logger::get("SM RawFakeFactor (lt)")->debug("Setting up functions for raw fake factor (without corrections) evaluation with correctionlib");

    Logger::get("SM RawFakeFactor (lt)")->debug("Fraction variation - Name {}", fraction_variation);
    Logger::get("SM RawFakeFactor (lt)")->debug("QCD variation - Name {}", qcd_variation);
    Logger::get("SM RawFakeFactor (lt)")->debug("Wjets variation - Name {}", wjets_variation);
    Logger::get("SM RawFakeFactor (lt)")->debug("ttbar variation - Name {}", ttbar_variation);

    auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
    auto wjets = correctionManager.loadCorrection(ff_file, "Wjets_fake_factors");
    auto ttbar = correctionManager.loadCorrection(ff_file, "ttbar_fake_factors");
    auto fractions = correctionManager.loadCorrection(ff_file, "process_fractions");


    // auto qcd = correction::CorrectionSet::from_file(ff_file)->at("QCD_fake_factors");
    // auto wjets = correction::CorrectionSet::from_file(ff_file)->at("Wjets_fake_factors");
    // auto ttbar = correction::CorrectionSet::from_file(ff_file)->at("ttbar_fake_factors");
    // auto fractions = correction::CorrectionSet::from_file(ff_file)->at("process_fractions");

    auto calc_fake_factor = [
        qcd, wjets, ttbar, fractions,
        qcd_variation, wjets_variation, ttbar_variation, fraction_variation](
        const float &pt_2, const int &njets, const float &mt_1, const float &delta_r) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM RawFakeFactor (lt)")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM RawFakeFactor (lt)")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, qcd_variation});
            float wjets_ff = wjets->evaluate({pt_2,  delta_r, (float)njets, wjets_variation});
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, ttbar_variation});

            Logger::get("SM RawFakeFactor (lt)")->debug("QCD - value {}", qcd_ff);
            Logger::get("SM RawFakeFactor (lt)")->debug("Wjets - value {}", wjets_ff);
            Logger::get("SM RawFakeFactor (lt)")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("SM RawFakeFactor (lt)")->debug("Lep mt - value {}", mt_1);
            Logger::get("SM RawFakeFactor (lt)")->debug("N jets - value {}", njets);

            float qcd_frac = fractions->evaluate({"QCD", mt_1, (float)njets, fraction_variation});
            float wjets_frac = fractions->evaluate({"Wjets", mt_1, (float)njets, fraction_variation});
            float ttbar_frac = fractions->evaluate({"ttbar", mt_1, (float)njets, fraction_variation});

            Logger::get("SM RawFakeFactor (lt)")->debug("QCD - fraction {}", qcd_frac);
            Logger::get("SM RawFakeFactor (lt)")->debug("Wjets - fraction {}", wjets_frac);
            Logger::get("SM RawFakeFactor (lt)")->debug("ttbar - fraction {}", ttbar_frac);


            ff = qcd_frac * std::max(qcd_ff, (float)0.) + 
                 wjets_frac * std::max(wjets_ff, (float)0.) +
                 ttbar_frac * std::max(ttbar_ff, (float)0.);
        }

        Logger::get("SM RawFakeFactor (lt)")->debug("Event Fake Factor {}", ff);

        return ff;
    };

    auto df1 = df.Define(outputname, calc_fake_factor, {tau_pt, njets, lep_mt, delta_r});

    return df1;
}
/**
 * @brief Function to calculate fake factors with correctionlib. In difference
 * to the NMSSM version, njets is used for the fraction binning, and an
 * additional split in deltaR is applied for wjets.
 * @param df the input dataframe
 * @param correctionManager The CorrectionManager object
 * @param outputname name of the output column for the fake factor
 * @param tau_pt pt of the hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param lep_mt transverse mass of the leptonic tau in the tau pair
 * @param lep_pt pt of the leptonic tau in the tau pair
 * @param lep_iso isolation of the leptonic tau in the tau pair
 * @param m_vis visible mass of the tau pair
 * @param delta_r distance in eta-phi between the two taus
 * @param fraction_variation name of the fraction variation
 * @param qcd_variation name of the QCD variation
 * @param wjets_variation name of the Wjets variation
 * @param ttbar_variation name of the ttbar variation
 * @param qcd_corr_leppt_variation name of the QCD lep pt correction variation
 * @param wjets_corr_leppt_variation name of the Wjets lep pt correction variation
 * @param ttbar_corr_leppt_variation name of the ttbar lep pt correction variation
 * @param qcd_corr_drsr_variation name of the QCD DR to SR correction variation
 * @param wjets_corr_drsr_variation name of the Wjets DR to SR correction variation
 * @param qcd_corr_lep_iso_variation name of the QCD lep iso correction variation
 * @param wjets_corr_lep_iso_variation name of the Wjets lep iso correction variation
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * factors
 * @returns a dataframe with the fake factors
 */
ROOT::RDF::RNode
fakefactor_sm_lt(
    ROOT::RDF::RNode df, 
    correctionManager::CorrectionManager &correctionManager,
    const std::string &outputname,
    const std::string &tau_pt,
    const std::string &njets,
    const std::string &lep_mt,
    const std::string &delta_r,
    const std::string &m_vis,
    const std::string &lep_pt,
    const std::string &lep_iso,
    //
    const std::string &fraction_variation,
    const std::string &qcd_variation,
    const std::string &wjets_variation,
    const std::string &ttbar_variation,
    //
    const std::string &qcd_corr_leppt_variation,
    const std::string &wjets_corr_leppt_variation,
    const std::string &ttbar_corr_leppt_variation,
    //
    const std::string &qcd_corr_drsr_variation,
    const std::string &wjets_corr_drsr_variation,
    //
    const std::string &qcd_corr_lep_iso_variation,
    const std::string &wjets_corr_lep_iso_variation,
    //
    const std::string &ff_file,
    const std::string &ff_corr_file
    ) {

    Logger::get("SM FaceFactor (lt)")->debug("Setting up functions for fake factor evaluation with correctionlib");
    
    Logger::get("SM FaceFactor (lt)")->debug("Fraction variation - Name {}", fraction_variation);
    Logger::get("SM FaceFactor (lt)")->debug("QCD variation - Name {}", qcd_variation);
    Logger::get("SM FaceFactor (lt)")->debug("Wjets variation - Name {}", wjets_variation);
    Logger::get("SM FaceFactor (lt)")->debug("ttbar variation - Name {}", ttbar_variation);

    Logger::get("SM FaceFactor (lt)")->debug("QCD lep pt correction variation - Name {}", qcd_corr_leppt_variation);
    Logger::get("SM FaceFactor (lt)")->debug("Wjets lep pt correction variation - Name {}", wjets_corr_leppt_variation);
    Logger::get("SM FaceFactor (lt)")->debug("ttbar lep pt correction variation - Name {}", ttbar_corr_leppt_variation);

    Logger::get("SM FaceFactor (lt)")->debug("QCD DR to SR correction variation - Name {}", qcd_corr_drsr_variation);
    Logger::get("SM FaceFactor (lt)")->debug("Wjets DR to SR correction variation - Name {}", wjets_corr_drsr_variation);

    Logger::get("SM FaceFactor (lt)")->debug("QCD lep iso correction variation - Name {}", qcd_corr_lep_iso_variation);
    Logger::get("SM FaceFactor (lt)")->debug("Wjets lep iso correction variation - Name {}", wjets_corr_lep_iso_variation);

    auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
    auto wjets = correctionManager.loadCorrection(ff_file, "Wjets_fake_factors");
    auto ttbar = correctionManager.loadCorrection(ff_file, "ttbar_fake_factors");
    
    auto fractions = correctionManager.loadCorrection(ff_file, "process_fractions");

    auto qcd_lep_pt_closure = correctionManager.loadCorrection(ff_corr_file, "QCD_non_closure_leading_lep_pt_correction");
    auto wjets_lep_pt_closure = correctionManager.loadCorrection(ff_corr_file, "Wjets_non_closure_leading_lep_pt_correction");
    auto ttbar_lep_pt_closure = correctionManager.loadCorrection(ff_corr_file, "ttbar_non_closure_leading_lep_pt_correction");
    
    auto qcd_DR_SR = correctionManager.loadCorrection(ff_corr_file, "QCD_DR_SR_correction");
    auto wjets_DR_SR = correctionManager.loadCorrection(ff_corr_file, "Wjets_DR_SR_correction");

    auto qcd_lep_iso_closure = correctionManager.loadCorrection(ff_corr_file, "QCD_non_closure_lep_iso_correction");
    auto wjets_lep_iso_closure = correctionManager.loadCorrection(ff_corr_file, "Wjets_non_closure_lep_iso_correction");
    
    
    auto calc_fake_factor = [
        qcd, wjets, ttbar, fractions,
        qcd_lep_pt_closure, wjets_lep_pt_closure, ttbar_lep_pt_closure,
        qcd_DR_SR, wjets_DR_SR,
        qcd_lep_iso_closure, wjets_lep_iso_closure,
        qcd_variation, wjets_variation, ttbar_variation, fraction_variation,        
        qcd_corr_leppt_variation, wjets_corr_leppt_variation, ttbar_corr_leppt_variation,
        qcd_corr_drsr_variation, wjets_corr_drsr_variation,
        qcd_corr_lep_iso_variation, wjets_corr_lep_iso_variation](
        const float &pt_2, const int &njets,
        const float &mt_1, const float &pt_1,
        const float &iso_1, const float &m_vis,
        const float &delta_r) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM FaceFactor (lt)")->debug("Tau pt - value {}", pt_2);
            Logger::get("SM FaceFactor (lt)")->debug("N jets - value {}", njets);

            float qcd_ff = qcd->evaluate({pt_2, (float)njets, qcd_variation});
            float wjets_ff = wjets->evaluate({pt_2, delta_r, (float)njets, wjets_variation});
            float ttbar_ff = ttbar->evaluate({pt_2, (float)njets, ttbar_variation});

            Logger::get("SM FaceFactor (lt)")->debug("QCD - value {}", qcd_ff);
            Logger::get("SM FaceFactor (lt)")->debug("Wjets - value {}", wjets_ff);
            Logger::get("SM FaceFactor (lt)")->debug("ttbar - value {}", ttbar_ff);

            Logger::get("SM FaceFactor (lt)")->debug("Lep mt - value {}", mt_1);
            Logger::get("SM FaceFactor (lt)")->debug("N jets - value {}", njets);

            float qcd_frac = fractions->evaluate({"QCD", mt_1, (float)njets, fraction_variation});
            float wjets_frac = fractions->evaluate({"Wjets", mt_1, (float)njets, fraction_variation});
            float ttbar_frac = fractions->evaluate({"ttbar", mt_1, (float)njets, fraction_variation});
            
            Logger::get("SM FaceFactor (lt)")->debug("QCD - fraction {}", qcd_frac);
            Logger::get("SM FaceFactor (lt)")->debug("Wjets - fraction {}", wjets_frac);
            Logger::get("SM FaceFactor (lt)")->debug("ttbar - fraction {}", ttbar_frac);

            Logger::get("SM FaceFactor (lt)")->debug("Lep pt - value {}", pt_1);
            Logger::get("SM FaceFactor (lt)")->debug("Lep iso - value {}", iso_1);
            Logger::get("SM FaceFactor (lt)")->debug("m_vis - value {}", m_vis);

            float qcd_lep_pt_corr = qcd_lep_pt_closure->evaluate({pt_1, qcd_corr_leppt_variation});
            float wjets_lep_pt_corr = wjets_lep_pt_closure->evaluate({pt_1, wjets_corr_leppt_variation});
            float ttbar_lep_pt_corr = ttbar_lep_pt_closure->evaluate({pt_1, ttbar_corr_leppt_variation});
            
            Logger::get("SM FaceFactor (lt)")->debug("QCD - lep pt correction {}", qcd_lep_pt_corr);
            Logger::get("SM FaceFactor (lt)")->debug("Wjets - lep pt correction {}", wjets_lep_pt_corr);
            Logger::get("SM FaceFactor (lt)")->debug("ttbar - lep pt correction {}", ttbar_lep_pt_corr);

            float qcd_lep_iso_corr = qcd_lep_iso_closure->evaluate({iso_1, qcd_corr_lep_iso_variation});
            float wjets_lep_iso_corr = wjets_lep_iso_closure->evaluate({iso_1, wjets_corr_lep_iso_variation});
            
            Logger::get("SM FaceFactor (lt)")->debug("QCD - lep iso correction {}", qcd_lep_iso_corr);
            Logger::get("SM FaceFactor (lt)")->debug("Wjets - lep iso correction {}", wjets_lep_iso_corr);

            float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, qcd_corr_drsr_variation});
            float wjets_DR_SR_corr = wjets_DR_SR->evaluate({m_vis, wjets_corr_drsr_variation});
            
            Logger::get("SM FaceFactor (lt)")->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);            
            Logger::get("SM FaceFactor (lt)")->debug("Wjets - DR to SR correction {}", wjets_DR_SR_corr);

            ff = qcd_frac * std::max(qcd_ff, (float)0.) * qcd_lep_pt_corr * qcd_lep_iso_corr * qcd_DR_SR_corr +
                 wjets_frac * std::max(wjets_ff, (float)0.) * wjets_lep_pt_corr * wjets_lep_iso_corr * wjets_DR_SR_corr +
                 ttbar_frac * std::max(ttbar_ff, (float)0.) * ttbar_lep_pt_corr;
        }

        Logger::get("SM FaceFactor (lt)")->debug("Event Fake Factor {}", ff);
        
        return ff;
    };
    
    auto df1 = df.Define(outputname,  calc_fake_factor,  {tau_pt, njets, lep_mt, lep_pt, lep_iso, m_vis, delta_r});
    
    return df1;
}
/** 
 * @brief Function to calculate raw fake factors with corrections. In difference
 * to the NMSSM version, njets is used for the fraction binning, and only QCD
 * fake factors are considered (for now).
 * @param df the input dataframe
 * @param correctionManager The CorrectionManager object
 * @param outputname name of the output column for the fake factor
 * @param tau_idx index of the tau, leading/subleading
 * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
 * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param m_vis visible mass of the tau pair
 * @param fraction_variation name of the fraction variation
 * @param qcd_variation name of the QCD variation
 * @param ff_file correctionlib json file with the fake factors
 * @param ff_corr_file correctionlib json file with corrections for the fake
 * @returns a dataframe with the raw fake factor
 */
ROOT::RDF::RNode
raw_fakefactor_sm_tt(
    ROOT::RDF::RNode df,
    correctionManager::CorrectionManager &correctionManager,
    const std::string &outputname,
    const int &tau_idx,
    const std::string &tau_pt_1,
    const std::string &tau_pt_2,
    const std::string &njets,
    const std::string &m_vis,
    //
    const std::string &fraction_variation,
    const std::string &qcd_variation,
    //
    const std::string &ff_file
) {

    Logger::get("SM RawFakeFactor (tt)")->debug("Setting up functions for raw fake factor (without corrections) evaluation with correctionlib");

    Logger::get("SM RawFakeFactor (tt)")->debug("QCD variation - Name {}", qcd_variation);
    Logger::get("SM RawFakeFactor (tt)")->debug("Fraction variation - Name {}", fraction_variation);

    auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
    auto qcd_subleading = correctionManager.loadCorrection(ff_file, "QCD_subleading_fake_factors");

    auto fractions = correctionManager.loadCorrection(ff_file, "process_fractions");
    auto fractions_subleading = correctionManager.loadCorrection(ff_file, "process_fractions_subleading");

    auto calc_fake_factor = [
        tau_idx, 
        qcd, qcd_subleading, fractions, fractions_subleading, 
        qcd_variation, fraction_variation](
        const float &pt_1, const float &pt_2,
        const int &njets, const float &m_vis) {
        float ff = -1.;
        if (pt_2 >= 0.) {
            Logger::get("SM RawFakeFactor (tt)")->debug("Leading Tau pt - value {}", pt_1);
            Logger::get("SM RawFakeFactor (tt)")->debug("Subleading Tau pt - value {}", pt_2);
            Logger::get("SM RawFakeFactor (tt)")->debug("N jets - value {}", njets);

            float qcd_ff = -1.;
            float qcd_frac = -1.;
            if (tau_idx == 0) {
                float qcd_ff = qcd->evaluate({pt_1, (float)njets, qcd_variation});
                float qcd_frac = fractions->evaluate({"QCD", m_vis, (float)njets, fraction_variation});
                
                Logger::get("SM RawFakeFactor (tt)")->debug("qcd_variation {}", qcd_variation);
                Logger::get("SM RawFakeFactor (tt)")->debug("fraction_variation {}", fraction_variation);
                Logger::get("SM RawFakeFactor (tt)")->debug("tau_idx {}", tau_idx);
                Logger::get("SM RawFakeFactor (tt)")->debug("QCD - value {}", qcd_ff);
                Logger::get("SM RawFakeFactor (tt)")->debug("QCD - fraction {}", qcd_frac);

                ff = qcd_frac * std::max(qcd_ff, (float)0.);
            } else if (tau_idx == 1) {
                float qcd_ff = qcd_subleading->evaluate({pt_2, (float)njets, qcd_variation});
                float qcd_frac = fractions_subleading->evaluate({"QCD", m_vis, (float)njets, fraction_variation});
                
                Logger::get("SM RawFakeFactor (tt)")->debug("qcd_variation {}", qcd_variation);
                Logger::get("SM RawFakeFactor (tt)")->debug("fraction_variation {}", fraction_variation);
                Logger::get("SM RawFakeFactor (tt)")->debug("tau_idx {}", tau_idx);
                Logger::get("SM RawFakeFactor (tt)")->debug("QCD - value {}", qcd_ff);
                Logger::get("SM RawFakeFactor (tt)")->debug("QCD - fraction {}", qcd_frac);

                ff = qcd_frac * std::max(qcd_ff, (float)0.);
            }
        }

        Logger::get("SM RawFakeFactor (tt)")->debug("Event Fake Factor {}", ff);
        
        return ff;
    };
    auto df1 = df.Define(outputname, calc_fake_factor, {tau_pt_1, tau_pt_2, njets, m_vis});
    
    return df1;
}
/** 
 * @brief Function to calculate raw fake factors with corrections. In difference
 * to the NMSSM version, njets is used for the fraction binning, and only QCD
 * fake factors are considered (for now).
 * @param df the input dataframe
 * @param correctionManager The CorrectionManager object
 * @param outputname name of the output column for the fake factor
 * @param tau_idx index of the tau, leading/subleading
 * @param tau_pt_1 pt of the leading hadronic tau in the tau pair
 * @param tau_pt_2 pt of the subleading hadronic tau in the tau pair
 * @param njets number of good jets in the event
 * @param m_vis visible mass of the tau pair
 * @param fraction_variation name of the fraction variation
 * @param qcd_variation name of the QCD variation
 * @param qcd_corr_leppt_variation name of the QCD lep pt correction variation
 * @param qcd_corr_drsr_variation name of the QCD DR to SR correction variation
 * @param qcd_corr_m_vis_variation name of the QCD m_vis correction variation
 * @param ff_file correctionlib json file with the fake factors
 * @returns a dataframe with the raw fake factor
 */
ROOT::RDF::RNode
fakefactor_sm_tt(
    ROOT::RDF::RNode df, 
    correctionManager::CorrectionManager &correctionManager,
    const std::string &outputname,
    const int &tau_idx,
    const std::string &tau_pt_1,
    const std::string &tau_pt_2,
    const std::string &njets,
    const std::string &m_vis,
    //
    const std::string &fraction_variation,
    const std::string &qcd_variation,
    //
    const std::string &qcd_corr_leppt_variation,
    //
    const std::string &qcd_corr_drsr_variation,
    //
    const std::string &qcd_corr_m_vis_variation,
    //
    const std::string &ff_file,
    const std::string &ff_corr_file
) {

    Logger::get("SM FakeFactor (tt)")->debug("Setting up functions for fake factor evaluation with correctionlib");
    
    Logger::get("SM FakeFactor (tt)")->debug("Fraction Variation - Name {}", fraction_variation);
    Logger::get("SM FakeFactor (tt)")->debug("QCD Variation - Name {}", qcd_variation);
    Logger::get("SM FakeFactor (tt)")->debug("QCD lep pt correction Variation - Name {}", qcd_corr_leppt_variation);
    Logger::get("SM FakeFactor (tt)")->debug("QCD DR to SR correction Variation - Name {}", qcd_corr_drsr_variation);
    Logger::get("SM FakeFactor (tt)")->debug("QCD m_vis correction Variation - Name {}", qcd_corr_m_vis_variation);

    auto qcd = correctionManager.loadCorrection(ff_file, "QCD_fake_factors");
    auto qcd_subleading = correctionManager.loadCorrection(ff_file, "QCD_subleading_fake_factors");

    auto qcd_tau_pt_closure = correctionManager.loadCorrection(ff_corr_file, "QCD_non_closure_subleading_lep_pt_correction");
    auto qcd_tau_pt_closure_subleading = correctionManager.loadCorrection(ff_corr_file, "QCD_subleading_non_closure_leading_lep_pt_correction");
    
    auto qcd_m_vis_closure = correctionManager.loadCorrection(ff_corr_file, "QCD_non_closure_m_vis_correction");
    auto qcd_m_vis_closure_subleading = correctionManager.loadCorrection(ff_corr_file, "QCD_subleading_non_closure_m_vis_correction");
    
    auto qcd_DR_SR = correctionManager.loadCorrection(ff_corr_file, "QCD_DR_SR_correction");
    auto qcd_DR_SR_subleading = correctionManager.loadCorrection(ff_corr_file, "QCD_subleading_DR_SR_correction");

    auto fractions = correctionManager.loadCorrection(ff_file, "process_fractions");
    auto fractions_subleading = correctionManager.loadCorrection(ff_file, "process_fractions_subleading");

    auto calc_fake_factor = [
        tau_idx, 
        qcd, qcd_subleading, fractions, fractions_subleading,
        qcd_tau_pt_closure, qcd_tau_pt_closure_subleading,
        qcd_m_vis_closure, qcd_m_vis_closure_subleading,
        qcd_DR_SR, qcd_DR_SR_subleading,
        qcd_variation, fraction_variation,                     
        qcd_corr_leppt_variation, 
        qcd_corr_drsr_variation, 
        qcd_corr_m_vis_variation](
        const float &pt_1, const float &pt_2,
        const int &njets, const float &m_vis) {
        float ff = 0.;
        if (pt_2 >= 0.) {
            Logger::get("SM FakeFactor (tt)")->debug("Leading Tau pt - value {}", pt_1);
            Logger::get("SM FakeFactor (tt)")->debug("Subleading Tau pt - value {}", pt_2);
            Logger::get("SM FakeFactor (tt)")->debug("m_vis - value {}", m_vis);
            Logger::get("SM FakeFactor (tt)")->debug("N jets - value {}", njets);

            float qcd_ff = -1.;
            float qcd_tau_pt_corr = -1.;
            float qcd_m_vis_corr = -1.;
            float qcd_DR_SR_corr = -1.;
            float qcd_frac = -1.;
            if (tau_idx == 0) {
                float qcd_ff = qcd->evaluate({pt_1, (float)njets, qcd_variation});
                float qcd_frac = fractions->evaluate({"QCD", m_vis, (float)njets, fraction_variation});
                
                float qcd_tau_pt_corr = qcd_tau_pt_closure->evaluate({pt_2, qcd_corr_leppt_variation});
                float qcd_m_vis_corr = qcd_m_vis_closure->evaluate({m_vis, qcd_corr_m_vis_variation});
                float qcd_DR_SR_corr = qcd_DR_SR->evaluate({m_vis, qcd_corr_drsr_variation});
                
                Logger::get("SM FakeFactor (tt)")->debug("QCD - leading value {}", qcd_ff);
                Logger::get("SM FakeFactor (tt)")->debug("QCD - leading fraction {}", qcd_frac);

                Logger::get("SM FakeFactor (tt)")->debug("QCD - leading lep pt correction {}", qcd_tau_pt_corr);
                Logger::get("SM FakeFactor (tt)")->debug("QCD - visible mass correction {}", qcd_m_vis_corr);
                Logger::get("SM FakeFactor (tt)")->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);
                
                ff = qcd_frac * std::max(qcd_ff, (float)0.) * qcd_tau_pt_corr * qcd_m_vis_corr * qcd_DR_SR_corr;

            } else if (tau_idx == 1) {
                float qcd_ff = qcd_subleading->evaluate({pt_2, (float)njets, qcd_variation});
                float qcd_frac = fractions_subleading->evaluate({"QCD", m_vis, (float)njets, fraction_variation});
                
                float qcd_tau_pt_corr = qcd_tau_pt_closure_subleading->evaluate({pt_1, qcd_corr_leppt_variation});
                float qcd_m_vis_corr = qcd_m_vis_closure_subleading->evaluate({m_vis, qcd_corr_m_vis_variation});
                float qcd_DR_SR_corr = qcd_DR_SR_subleading->evaluate({m_vis, qcd_corr_drsr_variation});

                Logger::get("SM FakeFactor (tt)")->debug("QCD - subleading value {}", qcd_ff);
                Logger::get("SM FakeFactor (tt)")->debug("QCD - subleading fraction {}", qcd_frac);

                Logger::get("SM FakeFactor (tt)")->debug("QCD - subleading lep pt correction {}", qcd_tau_pt_corr);
                Logger::get("SM FakeFactor (tt)")->debug("QCD - visible mass correction {}", qcd_m_vis_corr);
                Logger::get("SM FakeFactor (tt)")->debug("QCD - DR to SR correction {}", qcd_DR_SR_corr);

                ff = qcd_frac * std::max(qcd_ff, (float)0.) * qcd_tau_pt_corr * qcd_m_vis_corr * qcd_DR_SR_corr;
            }
        }

        Logger::get("SM FakeFactor (tt)")->debug("Event Fake Factor {}", ff);
        
        return ff;
    };

    auto df1 = df.Define(outputname, calc_fake_factor, {tau_pt_1, tau_pt_2, njets, m_vis});
    
    return df1;
}
} // namespace fakefactors
#endif /* GUARDFAKEFACTORS_H */
