#ifndef GUARD_TOPRECO_H
#define GUARD_TOPRECO_H

#include "../include/topreco.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TMinuit.h"
#include "TVector2.h"
#include "correction.h"
#include <Math/Boost.h>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

const float W_MASS = 80.377;  // PDG value as of 10/22
const float TOP_MASS = 172.5; // gen mass

namespace topreco {

/**
 * Function to create columns for isolated and antiisolated leptons quantities,
 * as used for a semileptonic top selection.
 *
 * @param[in] df The input dataframe.
 * @param[in] str_n_loose_mu Name of the column containing the number of loose
 * muons.
 * @param[in] str_n_loose_el Name of the column containing the number of loose
 * electrons.
 * @param[in] str_n_tight_mu Name of the column containing the number of tight
 * muons.
 * @param[in] str_n_tight_el Name of the column containing the number of tight
 * electrons.
 * @param[in] str_tight_muons_mask Name of the column containing the mask for
 * tight muons.
 * @param[in] str_tight_electrons_mask Name of the column containing the mask
 * for tight electrons.
 * @param[in] str_n_antitight_mu Name of the column containing the number of
 * antitight muons.
 * @param[in] str_n_antitight_el Name of the column containing the number of
 * antitight electrons.
 * @param[in] str_antitight_muons_mask Name of the column containing the mask
 * for antitight muons.
 * @param[in] str_antitight_electrons_mask Name of the column containing the
 * mask for antitight electrons.
 * @param[in] str_mu_pt Name of the column for the muon pT in nanoAOD.
 * @param[in] str_mu_eta Name of the column for the muon eta in nanoAOD.
 * @param[in] str_mu_phi Name of the column for the muon phi in nanoAOD.
 * @param[in] str_mu_mass Name of the column for the muon mass in nanoAOD.
 * @param[in] str_mu_charge Name of the column for the muon charge in nanoAOD.
 * @param[in] str_el_pt Name of the column for the electron pT in nanoAOD.
 * @param[in] str_el_eta Name of the column for the electron eta in nanoAOD.
 * @param[in] str_el_detasc Name of the column for the electron delta eta eta_sc
 * in nanoAOD.
 * @param[in] str_el_phi Name of the column for the electron phi in nanoAOD.
 * @param[in] str_el_mass Name of the column for the electron mass in nanoAOD.
 * @param[in] str_el_charge Name of the column for the electron charge in
 * nanoAOD.
 * @param[out] str_n_loose_lep Name of the output column for the number of loose
 * leptons.
 * @param[out] str_n_tight_lep Name of the output column for the number of tight
 * leptons.
 * @param[out] str_n_antitight_lep Name of the output column for the number of
 * antitight leptons.
 * @param[out] str_is_mu Name of the output column for a muon flag (tight and
 * antitight).
 * @param[out] str_is_el Name of the output column for an electron flag (tight
 * and antitight).
 * @param[out] str_is_iso Name of the output column for an isolation flag (muon
 * and electron).
 * @param[out] str_lep_p4 Name of the output column for the lepton four-momenta.
 * @param[out] str_lep_sceta Name of the output column for the lepton
 * supercluster eta.
 * @param[out] str_lep_charge Name of the output column for the lepton charge.
 * @param[out] str_mu_index Name of the output column for the muon index with
 * respect to the nanoAOD input branches.
 * @param[out] str_el_index Name of the output column for the electron index
 * with respect to the nanoAOD input branches.
 *
 * @return A dataframe containing the new columns.
 */
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
    const std::string &str_mu_index, const std::string &str_el_index) {

    auto df1 = df.Define(str_n_loose_lep,
                         [](const int &n_loose_mu, const int &n_loose_el) {
                             Logger::get("lepsel")->debug(
                                 "size of n_loose_mu and n_loose_el: {} {}",
                                 n_loose_mu, n_loose_el);
                             return n_loose_mu + n_loose_el;
                         },
                         {str_n_loose_mu, str_n_loose_el});

    auto df2 = df1.Define(str_n_tight_lep,
                          [](const int &n_tight_mu, const int &n_tight_el) {
                              Logger::get("lepsel")->debug(
                                  "size of n_tight_mu and n_tight_el: {} {}",
                                  n_tight_mu, n_tight_el);
                              return n_tight_mu + n_tight_el;
                          },
                          {str_n_tight_mu, str_n_tight_el});

    auto df3 =
        df2.Define(str_n_antitight_lep,
                   [](const int &n_antitight_mu, const int &n_antitight_el) {
                       Logger::get("lepsel")->debug(
                           "size of n_antitight_mu and n_antitight_el: {} {}",
                           n_antitight_mu, n_antitight_el);
                       return n_antitight_mu + n_antitight_el;
                   },
                   {str_n_antitight_mu, str_n_antitight_el});

    auto df4 = df3.Define(
        str_is_mu,
        [](const int &n_tight_mu, const int &n_antitight_mu,
           const int &n_loose_lep, const int &n_tight_lep,
           const int &n_antitight_lep) {
            return int(((n_tight_mu == 1) && (n_loose_lep == 1) &&
                        (n_tight_lep == 1)) ||
                       ((n_antitight_mu == 1) && (n_loose_lep == 0) &&
                        (n_antitight_lep == 1)));
        },
        {str_n_tight_mu, str_n_antitight_mu, str_n_loose_lep, str_n_tight_lep,
         str_n_antitight_lep});

    auto df5 = df4.Define(
        str_is_el,
        [](const int &n_tight_el, const int &n_antitight_el,
           const int &n_loose_lep, const int &n_tight_lep,
           const int &n_antitight_lep) {
            return int(((n_tight_el == 1) && (n_loose_lep == 1) &&
                        (n_tight_lep == 1)) ||
                       ((n_antitight_el == 1) && (n_loose_lep == 0) &&
                        (n_antitight_lep == 1)));
        },
        {str_n_tight_el, str_n_antitight_el, str_n_loose_lep, str_n_tight_lep,
         str_n_antitight_lep});

    auto is_iso = [](const int &is_mu, const int &is_el, const int &n_tight_lep,
                     const int &n_antitight_lep) {
        if (is_mu || is_el) {
            if (n_tight_lep == 1)
                return +1;
            if (n_antitight_lep == 1)
                return -1;
        }

        return 0;
    };

    auto df6 = df5.Define(
        str_is_iso, is_iso,
        {str_is_mu, str_is_el, str_n_tight_lep, str_n_antitight_lep});

    auto lep_mu_index = [](const int is_mu, const int is_iso,
                           const ROOT::RVec<int> &tight_muons_mask,
                           const ROOT::RVec<int> &antitight_muons_mask) {
        Logger::get("lep_mu_index")
            ->debug("mask mu {}, mask size mu {}, max in mu {}, index of max "
                    "in mu {}",
                    tight_muons_mask, tight_muons_mask.size(),
                    ROOT::VecOps::Max(tight_muons_mask),
                    ROOT::VecOps::ArgMax(tight_muons_mask));

        int lep_mu_index = -1;

        if (is_iso == 0 || is_mu == 0) {
            return lep_mu_index;
        } else if (is_iso == +1 && is_mu == 1) {
            Logger::get("lep_mu_index")->debug("---> should reco iso mu...");
            lep_mu_index = ROOT::VecOps::ArgMax(tight_muons_mask);
        } else if (is_iso == -1 && is_mu == 1) {
            Logger::get("lep_mu_index")
                ->debug("---> should reco antiiso mu...");
            lep_mu_index = ROOT::VecOps::ArgMax(antitight_muons_mask);
        }

        return lep_mu_index;
    };

    auto df7 = df6.Define(str_mu_index, lep_mu_index,
                          {str_is_mu, str_is_iso, str_tight_muons_mask,
                           str_antitight_muons_mask});

    auto lep_el_index = [](const int is_el, const int is_iso,
                           const ROOT::RVec<int> &tight_electrons_mask,
                           const ROOT::RVec<int> &antitight_electrons_mask) {
        Logger::get("lep_el_index")
            ->debug("mask el {}, mask size el {}, max in el {}, index of max "
                    "in el {}",
                    tight_electrons_mask, tight_electrons_mask.size(),
                    ROOT::VecOps::Max(tight_electrons_mask),
                    ROOT::VecOps::ArgMax(tight_electrons_mask));

        int lep_el_index = -1;

        if (is_iso == 0 || is_el == 0) {
            return lep_el_index;
        } else if (is_iso == +1 && is_el == 1) {
            Logger::get("lep_el_index")->debug("---> should reco iso el...");
            lep_el_index = ROOT::VecOps::ArgMax(tight_electrons_mask);
        } else if (is_iso == -1 && is_el == 1) {
            Logger::get("lep_el_index")
                ->debug("---> should reco antiiso el...");
            lep_el_index = ROOT::VecOps::ArgMax(antitight_electrons_mask);
        }

        return lep_el_index;
    };

    auto df7a = df7.Define(str_el_index, lep_el_index,
                           {str_is_el, str_is_iso, str_tight_electrons_mask,
                            str_antitight_electrons_mask});

    auto lep_p4 =
        [](const int is_mu, const int is_el, const int is_iso,
           const int mu_index, const int el_index,
           const ROOT::RVec<float> &mu_pt, const ROOT::RVec<float> &mu_eta,
           const ROOT::RVec<float> &mu_phi, const ROOT::RVec<float> &mu_mass,
           const ROOT::RVec<float> &el_pt, const ROOT::RVec<float> &el_eta,
           const ROOT::RVec<float> &el_phi, const ROOT::RVec<float> &el_mass) {
            auto lep = ROOT::Math::PtEtaPhiMVector(-10, -10, -10, -10);

            if (is_iso == 0) {
                return lep;
            } else if (is_iso == +1) {
                if (is_mu) {
                    lep = ROOT::Math::PtEtaPhiMVector(
                        mu_pt.at(mu_index, -2), mu_eta.at(mu_index, -2),
                        mu_phi.at(mu_index, -2), mu_mass.at(mu_index, -2));
                } else if (is_el) {
                    lep = ROOT::Math::PtEtaPhiMVector(
                        el_pt.at(el_index, -3), el_eta.at(el_index, -3),
                        el_phi.at(el_index, -3), el_mass.at(el_index, -3));
                }
            } else if (is_iso == -1) {
                if (is_mu) {
                    lep = ROOT::Math::PtEtaPhiMVector(
                        mu_pt.at(mu_index, -4), mu_eta.at(mu_index, -4),
                        mu_phi.at(mu_index, -4), mu_mass.at(mu_index, -4));
                } else if (is_el) {
                    lep = ROOT::Math::PtEtaPhiMVector(
                        el_pt.at(el_index, -5), el_eta.at(el_index, -5),
                        el_phi.at(el_index, -5), el_mass.at(el_index, -5));
                }
            } else {
                lep = ROOT::Math::PtEtaPhiMVector(-6, -6, -6, -6);
            }

            Logger::get("final_lep")
                ->debug("building p4 from lepton with {} {} {} {}", lep.Pt(),
                        lep.Eta(), lep.Phi(), lep.M());

            return lep;
        };

    auto df7b = df7a.Define(str_lep_p4, lep_p4,
                            {str_is_mu, str_is_el, str_is_iso, str_mu_index,
                             str_el_index, str_mu_pt, str_mu_eta, str_mu_phi,
                             str_mu_mass, str_el_pt, str_el_eta, str_el_phi,
                             str_el_mass});

    auto lep_sceta = [](const int is_el, const int is_iso, const int el_index,
                        const ROOT::RVec<float> &el_eta,
                        const ROOT::RVec<float> &el_detasc) {
        float lep_sceta = -10;

        if (is_iso == 0) {
            return lep_sceta;
        } else if (is_iso == +1) {
            if (is_el) {
                lep_sceta =
                    el_eta.at(el_index, -5) + el_detasc.at(el_index, -5);
            }
        } else if (is_iso == -1) {
            if (is_el) {
                lep_sceta =
                    el_eta.at(el_index, -5) + el_detasc.at(el_index, -5);
            }
        }
        return lep_sceta;
    };

    auto df7c = df7b.Define(
        str_lep_sceta, lep_sceta,
        {str_is_el, str_is_iso, str_el_index, str_el_eta, str_el_detasc});

    auto lep_charge = [](const int is_mu, const int is_el, const int is_iso,
                         const int mu_index, const int el_index,
                         const ROOT::RVec<int> &mu_charge,
                         const ROOT::RVec<int> &el_charge) {
        int charge = -10;

        if (is_iso == 0) {
            return charge;
        } else if (is_iso == +1) {
            if (is_mu) {
                charge = mu_charge.at(mu_index, -2);
            } else if (is_el) {
                charge = el_charge.at(el_index, -3);
            }
        } else if (is_iso == -1) {
            if (is_mu) {
                charge = mu_charge.at(mu_index, -4);
            } else if (is_el) {
                charge = el_charge.at(el_index, -5);
            }
        } else {
            charge = -6;
        }

        return charge;
    };

    auto df8 = df7c.Define(str_lep_charge, lep_charge,
                           {str_is_mu, str_is_el, str_is_iso, str_mu_index,
                            str_el_index, str_mu_charge, str_el_charge});

    return df8;
}

// helper function for minimizer constraint
double rad_py(double x, double lep_px) {
    return W_MASS * W_MASS + 4 * lep_px * x;
}

// the delta plus function with the py nu plus solution
double min_fplus(double *par) {
    // par[0] = x, par[1] = lep_px, par[2] = lep_py, par[3] = lep_pt, par[4] =
    // px_miss, par[5] = py_miss
    double r = rad_py(par[0], par[1]);
    double y = 0;
    // double res = 99999;
    // if (r>=0) {
    y = (W_MASS * W_MASS * par[2] + 2 * par[1] * par[2] * par[0] +
         W_MASS * par[3] * sqrt(r)) /
        (2 * par[1] * par[1]);
    double res = sqrt((par[0] - par[4]) * (par[0] - par[4]) +
                      (y - par[5]) * (y - par[5]));
    // }
    // else // FIXME: proper constraint in TMinuit?
    // res = 99999;
    return res;
}

// the delta minus function with the py nu minus solution
double min_fminus(double *par) {
    // par[0] = x, par[1] = lep_px, par[2] = lep_py, par[3] = lep_pt, par[4] =
    // px_miss, par[5] = py_miss
    double r = rad_py(par[0], par[1]);
    double y = 0;
    double res = 99999;
    if (r >= 0) {
        y = (W_MASS * W_MASS * par[2] + 2 * par[1] * par[2] * par[0] -
             W_MASS * par[3] * sqrt(r)) /
            (2 * par[1] * par[1]);
        res = sqrt((par[0] - par[4]) * (par[0] - par[4]) +
                   (y - par[5]) * (y - par[5]));
    } else
        res = 99999;
    return res;
}

// TMinuit fit function for the py nu plus solution
void fcn_plus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
              Int_t iflag) {
    f = min_fplus(par);
}

// TMinuit fit function for the py nu minus solution
void fcn_minus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par,
               Int_t iflag) {
    f = min_fminus(par);
}

/**
 * Function to reconstruct a leptonic W boson based on mass constraints and a
 * numerical quadratic equation solver.
 *
 * @param[in] df The input dataframe.
 * @param[in] str_lep_p4 Name of the column containing the lepton four-momenta.
 * @param[in] str_met_p4 Name of the column containing the missing transverse
 * momentum four-momenta.
 * @param[out] str_wlep_p4 Name of the output column for reconstructed W boson
 * four-momenta.
 *
 * @return A dataframe containing the new columns.
 */
ROOT::RDF::RNode ReconstructLeptonicW(ROOT::RDF::RNode df,
                                      const std::string &str_lep_p4,
                                      const std::string &str_met_p4,
                                      const std::string &str_wlep_p4) {

    auto leptonicW = [](const ROOT::Math::PtEtaPhiMVector lep_p4,
                        const ROOT::Math::PtEtaPhiMVector met_p4) {
        auto wlep_p4 = ROOT::Math::PtEtaPhiMVector(-10, -10, -10, -10);

        if (lep_p4.Pt() < 0)
            return wlep_p4;

        double lep_e = lep_p4.E();
        double lep_pt = lep_p4.Pt();
        double lep_px = lep_p4.Px();
        double lep_py = lep_p4.Py();
        double lep_pz = lep_p4.Pz();

        ROOT::Math::PtEtaPhiMVector nu_p4 = met_p4;
        double nu_e = met_p4.Pt();
        double nu_px = met_p4.Px();
        double nu_py = met_p4.Py();

        Logger::get("wlep")->debug(
            "building wlep p4 from lepton with E: {} px: {} py: {} pz: {}",
            lep_e, lep_px, lep_py, lep_pz);
        Logger::get("wlep")->debug(
            "building wlep p4 from pTmiss with E: {} px: {} py: {} pz: ???",
            nu_e, nu_px, nu_py);

        // definition of the constant mu in Eq. 4.5 (here called alpha to not
        // confuse mu and nu) also converting p_T and cos dphi into p_x and p_y
        double alpha =
            (W_MASS * W_MASS) / 2 + (lep_px * nu_px) + (lep_py * nu_py);

        // for p_z,nu there is a quadratic equation with up to two solutions as
        // shown in Eq. 4.6 and A.7 (NOTE: there is a 'power of two' missing in
        // the first denominator of Eq. 4.6) first, check if we have complex
        // solution, i.e. if the radicand is negative
        double rad =
            ((alpha * alpha * lep_pz * lep_pz) /
             (lep_pt * lep_pt * lep_pt * lep_pt)) -
            ((lep_e * lep_e * nu_e * nu_e - alpha * alpha) / (lep_pt * lep_pt));

        if (rad < 0) {
            // complex solutions, in around 30% of all cases
            //    cout << "Complex neutrino p_z" << endl;

            // assumption: p_T^miss does not directly correspond to p_T,nu
            // set m_T^W to m^W, result is a quadratic equation for p_(x,y) that
            // depends on p_(y,x)

            // save p_x^miss and p_y^miss as we need them later to determine the
            // better solution

            Logger::get("wlep")->debug("complex solution");

            double px_miss = nu_px;
            double py_miss = nu_py;

            Logger::get("wlep")->debug("complex debug point 1");

            // initialize TMinuit with a maximum of 6 params for py nu plus and
            // minus solution
            TMinuit *gMinuit_plus = new TMinuit(5);
            TMinuit *gMinuit_minus = new TMinuit(5);

            Logger::get("wlep")->debug("complex debug point 2");

            //    TMinuit *gMinuit_plus = new TMinuit(5);
            gMinuit_plus->SetFCN(fcn_plus);

            //    TMinuit *gMinuit_minus = new TMinuit(5);
            gMinuit_minus->SetFCN(fcn_minus);

            Logger::get("wlep")->debug("complex debug point 3");

            double arglist[10];
            int ierflg = 0;

            // no print
            arglist[0] = -1;
            gMinuit_plus->mnexcm("SET PRI", arglist, 1, ierflg);
            gMinuit_minus->mnexcm("SET PRI", arglist, 1, ierflg);

            // no warnings
            arglist[0] = -1;
            gMinuit_plus->mnexcm("SET NOW", arglist, 1, ierflg);
            gMinuit_minus->mnexcm("SET NOW", arglist, 1, ierflg);
            //    gMinuit_plus->mnexcm("SET WAR", arglist, 1, ierflg);
            //    gMinuit_minus->mnexcm("SET WAR", arglist, 1, ierflg);

            Logger::get("wlep")->debug("complex debug point 4");

            /*
            // set accuracy
            arglist[0] = 1.E-5L;
            gMinuit_plus->mnexcm("SET EPS", arglist, 1, ierflg);
            gMinuit_minus->mnexcm("SET EPS", arglist, 1, ierflg);
            */
            // set error
            arglist[0] = 1; // 0.5 for lnL, 1 for chi2
            // set the 1 sigma tolerance for the change in FCN
            // that determines when a function has been minimized
            gMinuit_plus->mnexcm("SET ERR", arglist, 1, ierflg);
            gMinuit_minus->mnexcm("SET ERR", arglist, 1, ierflg);

            // set strategy (0 is less accurate, 1 is default, 2 is more
            // accurate)
            arglist[0] = 0;
            gMinuit_plus->mnexcm("SET STR", arglist, 1, ierflg);
            gMinuit_minus->mnexcm("SET STR", arglist, 1, ierflg);

            Logger::get("wlep")->debug("complex debug point 5");

            // set start value and range of x and fix all other params
            static double start_val = 0;
            static double step = 0.00001;
            double lower = -9999;
            double upper = 9999;
            if (lep_px > 0) {
                lower = -W_MASS * W_MASS / (4 * lep_px) + 1e-5;
                start_val = lower + 1;
                //      cout << "lower: " << lower << endl;
            }
            if (lep_px < 0) {
                upper = -W_MASS * W_MASS / (4 * lep_px) - 1e-5;
                start_val = upper - 1;
                //      cout << "upper: " << upper << endl;
            }
            gMinuit_plus->mnparm(0, "x", start_val, step, lower, upper, ierflg);
            gMinuit_plus->mnparm(1, "lep_px", lep_px, step, -999, 999, ierflg);
            gMinuit_plus->mnparm(2, "lep_py", lep_py, step, -999, 999, ierflg);
            gMinuit_plus->mnparm(3, "lep_pt", lep_pt, step, 0, 999, ierflg);
            gMinuit_plus->mnparm(4, "px_miss", px_miss, step, -999, 999,
                                 ierflg);
            gMinuit_plus->mnparm(5, "py_miss", py_miss, step, -999, 999,
                                 ierflg);
            gMinuit_plus->FixParameter(1);
            gMinuit_plus->FixParameter(2);
            gMinuit_plus->FixParameter(3);
            gMinuit_plus->FixParameter(4);
            gMinuit_plus->FixParameter(5);

            gMinuit_minus->mnparm(0, "x", start_val, step, lower, upper,
                                  ierflg);
            gMinuit_minus->mnparm(1, "lep_px", lep_px, step, -999, 999, ierflg);
            gMinuit_minus->mnparm(2, "lep_py", lep_py, step, -999, 999, ierflg);
            gMinuit_minus->mnparm(3, "lep_pt", lep_pt, step, 0, 999, ierflg);
            gMinuit_minus->mnparm(4, "px_miss", px_miss, step, -999, 999,
                                  ierflg);
            gMinuit_minus->mnparm(5, "py_miss", py_miss, step, -999, 999,
                                  ierflg);
            gMinuit_minus->FixParameter(1);
            gMinuit_minus->FixParameter(2);
            gMinuit_minus->FixParameter(3);
            gMinuit_minus->FixParameter(4);
            gMinuit_minus->FixParameter(5);

            Logger::get("wlep")->debug("complex debug point 6");

            // now ready for minimization step
            arglist[0] = 500000; // maximum number of iterations
            arglist[1] = 1.;     // related to errors
            gMinuit_plus->mnexcm("MIGARD", arglist, 2, ierflg);
            gMinuit_minus->mnexcm("MIGRAD", arglist, 2, ierflg);

            // obtain fit results and calculate values of delta minus and delta
            // plus functions choose solution that leads to a smaller delta
            // value
            double x_plus, x_pluserr;
            double d_plus;
            gMinuit_plus->GetParameter(0, x_plus, x_pluserr);
            double par_plus[6] = {x_plus, lep_px,  lep_py,
                                  lep_pt, px_miss, py_miss};
            d_plus = min_fplus(par_plus);
            //    cout << "Fit result plus: x=" << x_plus << " " << "d(x)=" <<
            //    d_plus << endl;

            double x_minus, x_minuserr;
            double d_minus;
            gMinuit_minus->GetParameter(0, x_minus, x_minuserr);
            double par_minus[6] = {x_minus, lep_px,  lep_py,
                                   lep_pt,  px_miss, py_miss};
            d_minus = min_fminus(par_minus);
            //    cout << "Fit result minus: x=" << x_minus << " d(x)=" <<
            //    d_minus << endl;

            Logger::get("wlep")->debug("complex debug point 7");

            double nu_pxnew, nu_pynew, r_new;
            if (d_plus < d_minus) {
                nu_pxnew = x_plus;
                r_new = rad_py(nu_pxnew, lep_px);
                nu_pynew =
                    (W_MASS * W_MASS * lep_py + 2 * lep_px * lep_py * nu_pxnew +
                     W_MASS * lep_pt * sqrt(r_new)) /
                    (2 * lep_px * lep_px);
            } else {
                nu_pxnew = x_minus;
                r_new = rad_py(nu_pxnew, lep_px);
                nu_pynew =
                    (W_MASS * W_MASS * lep_py + 2 * lep_px * lep_py * nu_pxnew -
                     W_MASS * lep_pt * sqrt(r_new)) /
                    (2 * lep_px * lep_px);
            }
            // calculate new nu pz (only one solution with fixed px and py)
            double nu_pznew = lep_pz / (lep_pt * lep_pt) *
                              ((W_MASS * W_MASS / 2) + (lep_px * nu_pxnew) +
                               (lep_py * nu_pynew));
            //    cout << "new nu px: " << nu_pxnew << ", new nu py: " <<
            //    nu_pynew << ", new nu pz: " << nu_pznew << endl;

            // set 4 momenta of neutrino and W boson
            nu_p4.SetPxPyPzE(nu_pxnew, nu_pynew, nu_pznew,
                             sqrt(nu_pxnew * nu_pxnew + nu_pynew * nu_pynew +
                                  nu_pznew * nu_pznew));

            Logger::get("wlep")->debug("complex debug point 8");

            delete gMinuit_plus;
            delete gMinuit_minus;

        }

        else { // two real solutions for pz nu
            //    cout << "Two neutrino pz solutions" << endl;
            Logger::get("wlep")->debug("real solution");
            double sol1, sol2, nu_pz;
            sol1 = lep_pz * alpha / (lep_pt * lep_pt) + sqrt(rad);
            sol2 = lep_pz * alpha / (lep_pt * lep_pt) - sqrt(rad);

            // choose the smaller pz solution
            if (abs(sol1) < abs(sol2)) {
                nu_pz = sol1;
            } else {
                nu_pz = sol2;
            }

            // set 4 momenta of neutrino and W boson
            nu_p4.SetPxPyPzE(
                nu_px, nu_py, nu_pz,
                sqrt(nu_px * nu_px + nu_py * nu_py + nu_pz * nu_pz));
        }

        wlep_p4 = lep_p4 + nu_p4;

        Logger::get("wlep")->debug(
            "final leptonic W boson pT: {} eta: {} phi: {} mass: {}",
            wlep_p4.Pt(), wlep_p4.Eta(), wlep_p4.Phi(), wlep_p4.M());

        return wlep_p4;
    };

    return df.Define(str_wlep_p4, leptonicW, {str_lep_p4, str_met_p4});
}

/**
 * Selects events with a given number of jets and b-jets.
 *
 * @param[in] df The input dataframe.
 * @param[in] njets The number of jets to select.
 * @param[in] nbjets The number of b-jets to select.
 * @param[in] str_good_jets_mask The name of the column containing the mask for
 * good jets.
 * @param[in] str_good_bjets_mask The name of the column containing the mask for
 * good b-jets.
 *
 * @return A new dataframe containing the selected events.
 */
ROOT::RDF::RNode JetSelection(ROOT::RDF::RNode df, const int &njets,
                              const int &nbjets,
                              const std::string &str_good_jets_mask,
                              const std::string &str_good_bjets_mask) {

    auto df2 = df.Filter(
        [njets, nbjets](const ROOT::RVec<int> &good_jets_mask,
                        const ROOT::RVec<int> &good_bjets_mask) {
            int nj = ROOT::VecOps::Sum(good_jets_mask, 0);
            int nbj = ROOT::VecOps::Sum(good_bjets_mask, 0);

            return (nj == njets) && (nbj == nbjets);
        },
        {str_good_jets_mask, str_good_bjets_mask},
        "jet selection and b jet selection");

    return df2;
}

/**
 * Function to reconstruct a leptonic top quark from a b jet and a W boson based
 * on the number of jets and b jets in the event.
 *
 * @param[in] df The input dataframe.
 * @param[in] str_wlep_p4 Name of the column containing the W boson
 * four-momenta.
 * @param[in] str_n_nonbjets Name of the column containing the number of
 * untagged jets in the event.
 * @param[in] str_nonbjet_p4_1 Name of the column containing the hardest
 * untagged jet four-momenta.
 * @param[in] str_nonbjet_btag_1 Name of the column containing the hardest
 * untagged jet b tagging score.
 * @param[in] str_nonbjet_p4_2 Name of the column containing the second hardest
 * untagged jet four-momenta.
 * @param[in] str_nonbjet_btag_2 Name of the column containing the second
 * hardest untagged jet b tagging score.
 * @param[in] str_n_bjets Name of the column containing the number of b-tagged
 * jets in the event.
 * @param[in] str_bjet_p4_1 Name of the column containing the hardest b-tagged
 * jet four-momenta.
 * @param[in] str_bjet_btag_1 Name of the column containing the hardest b-tagged
 * jet b tagging score.
 * @param[in] str_bjet_p4_2 Name of the column containing the second hardest
 * b-tagged jet four-momenta.
 * @param[in] str_bjet_btag_2 Name of the column containing the second hardest
 * b-tagged jet b tagging score.
 * @param[out] str_is_reco Name of the output column for the flag if the top
 * quark is reconstructable.
 * @param[out] str_is_jj Name of the output column for the flag if the event
 * falls into the 2j0b category.
 * @param[out] str_is_jjb Name of the output column for the flag if the event
 * falls into the 2j1b category.
 * @param[out] str_is_jjbb Name of the output column for the flag if the event
 * falls into the 2j2b category.
 * @param[out] str_is_jjjb Name of the output column for the flag if the event
 * falls into the 3j1b category.
 * @param[out] str_is_jjjbb Name of the output column for the flag if the event
 * falls into the 3j2b category.
 * @param[out] str_reco_p4s Name of the output column for the vector of
 * reconstructed particles four-momenta.
 * @param[out] str_top_p4 Name of the output column for the reconstructed top
 * quark four-momenta.
 * @param[out] str_tb_p4 Name of the output column for the reconstructed b quark
 * (top quark decay) four-momenta.
 * @param[out] str_sb_p4 Name of the output column for the reconstructed b quark
 * (top quark production) four-momenta.
 *
 * @return A dataframe containing the new columns.
 */
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
        const std::string &str_tb_p4, const std::string &str_sb_p4) {

    auto df2a =
        df.Define(str_is_jj,
                  [](const int n_nonbjets, const int n_bjets) {
                      return int((n_nonbjets + n_bjets) == 2 && n_bjets == 0);
                  },
                  {str_n_nonbjets, str_n_bjets});

    auto df2 =
        df2a.Define(str_is_jjb,
                    [](const int n_nonbjets, const int n_bjets) {
                        return int((n_nonbjets + n_bjets) == 2 && n_bjets == 1);
                    },
                    {str_n_nonbjets, str_n_bjets});

    auto df3 =
        df2.Define(str_is_jjbb,
                   [](const int n_nonbjets, const int n_bjets) {
                       return int((n_nonbjets + n_bjets) == 2 && n_bjets == 2);
                   },
                   {str_n_nonbjets, str_n_bjets});

    auto df4 =
        df3.Define(str_is_jjjb,
                   [](const int n_nonbjets, const int n_bjets) {
                       return int((n_nonbjets + n_bjets) == 3 && n_bjets == 1);
                   },
                   {str_n_nonbjets, str_n_bjets});

    auto df5 =
        df4.Define(str_is_jjjbb,
                   [](const int n_nonbjets, const int n_bjets) {
                       return int((n_nonbjets + n_bjets) == 3 && n_bjets == 2);
                   },
                   {str_n_nonbjets, str_n_bjets});

    auto df6 = df5.Define(str_is_reco,
                          [](const int is_jjb, const int is_jjbb,
                             const int is_jjjb, const int is_jjjbb) {
                              return is_jjb + is_jjbb + is_jjjb + is_jjjbb;
                          },
                          {str_is_jjb, str_is_jjbb, str_is_jjjb, str_is_jjjbb});

    auto top_reco = [](const int is_reco, const int is_jj, const int is_jjb,
                       const int is_jjbb, const int is_jjjb, const int is_jjjbb,
                       const ROOT::Math::PtEtaPhiMVector wlep_p4,
                       const ROOT::Math::PtEtaPhiMVector nonbjet_p4_1,
                       const float nonbjet_btag_1,
                       const ROOT::Math::PtEtaPhiMVector nonbjet_p4_2,
                       const float nonbjet_btag_2,
                       const ROOT::Math::PtEtaPhiMVector bjet_p4_1,
                       const float bjet_btag_1,
                       const ROOT::Math::PtEtaPhiMVector bjet_p4_2,
                       const float bjet_btag_2) {
        ROOT::RVec<ROOT::Math::PtEtaPhiMVector> reco_vec{
            ROOT::Math::PtEtaPhiMVector(-10, -10, -10, -10), // top_p4
            ROOT::Math::PtEtaPhiMVector(-10, -10, -10, -10), // tb_p4
            ROOT::Math::PtEtaPhiMVector(-10, -10, -10, -10)  // sb_p4
        };

        if (wlep_p4.Pt() < 0)
            return reco_vec;

        if (!is_reco)
            return reco_vec;

        if (is_jjb) { // 2j1b
            reco_vec[0] = wlep_p4 + bjet_p4_1;
            reco_vec[1] = bjet_p4_1;
            reco_vec[2] = nonbjet_p4_1;
        } else if (is_jjbb) { // 2j2b
            auto cand1_p4 = wlep_p4 + bjet_p4_1;
            auto cand2_p4 = wlep_p4 + bjet_p4_2;
            if (abs(cand1_p4.M() - TOP_MASS) < abs(cand2_p4.M() - TOP_MASS)) {
                reco_vec[0] = cand1_p4;
                reco_vec[1] = bjet_p4_1;
                reco_vec[2] = bjet_p4_2;
            } else {
                reco_vec[0] = cand2_p4;
                reco_vec[1] = bjet_p4_2;
                reco_vec[2] = bjet_p4_1;
            }
        } else if (is_jjjb) { // 3j1b
            reco_vec[0] = wlep_p4 + bjet_p4_1;
            reco_vec[1] = bjet_p4_1;
            if (nonbjet_btag_1 > nonbjet_btag_2)
                reco_vec[2] = nonbjet_p4_1;
            else
                reco_vec[2] = nonbjet_p4_2;
        } else if (is_jjjbb) { // 3j2b
            auto cand1_p4 = wlep_p4 + bjet_p4_1;
            auto cand2_p4 = wlep_p4 + bjet_p4_2;
            if (abs(cand1_p4.M() - TOP_MASS) < abs(cand2_p4.M() - TOP_MASS)) {
                reco_vec[0] = cand1_p4;
                reco_vec[1] = bjet_p4_1;
                reco_vec[2] = bjet_p4_2;
            } else {
                reco_vec[0] = cand2_p4;
                reco_vec[1] = bjet_p4_2;
                reco_vec[2] = bjet_p4_1;
            }
        }

        return reco_vec;
    };

    auto df7 = df6.Define(str_reco_p4s, top_reco,
                          {str_is_reco, str_is_jj, str_is_jjb, str_is_jjbb,
                           str_is_jjjb, str_is_jjjbb, str_wlep_p4,
                           str_nonbjet_p4_1, str_nonbjet_btag_1,
                           str_nonbjet_p4_2, str_nonbjet_btag_2, str_bjet_p4_1,
                           str_bjet_btag_1, str_bjet_p4_2, str_bjet_btag_2});

    auto df8 =
        df7.Define(str_top_p4,
                   [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &vec) {
                       return vec[0];
                   },
                   {str_reco_p4s});

    auto df9 =
        df8.Define(str_tb_p4,
                   [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &vec) {
                       return vec[1];
                   },
                   {str_reco_p4s});

    auto df10 =
        df9.Define(str_sb_p4,
                   [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &vec) {
                       return vec[2];
                   },
                   {str_reco_p4s});

    return df10;
}

/**
 * Function to create additional variables useful in a DNN training
 *
 * @param[in] df The input dataframe.
 * @param[in] str_is_reco The name of the column containing the flag if a top
 *                        quark was reconstructed.
 * @param[in] str_lep_p4 The name of the column containing the lepton
 * four-momenta.
 * @param[in] str_met_p4 The name of the column containing the missing
 * transverse momentum four-momenta.
 * @param[in] str_wlep_p4 The name of the column containing the W boson
 * four-momenta.
 * @param[in] str_nonbjet_p4_1 The name of the column containing the hardest
 * untagged jet four-momenta.
 * @param[in] str_nonbjet_btag_1 The name of the column containing the hardest
 * untagged jet b tagging score.
 * @param[in] str_nonbjet_p4_2 The name of the column containing the second
 * hardest untagged jet four-momenta.
 * @param[in] str_nonbjet_btag_2 The name of the column containing the second
 * hardest untagged jet b tagging score.
 * @param[in] str_bjet_p4_1 The name of the column containing the hardest
 * b-tagged jet four-momenta.
 * @param[in] str_bjet_btag_1 The name of the column containing the hardest
 * b-tagged jet b tagging score.
 * @param[in] str_bjet_p4_2 The name of the column containing the second hardest
 * b-tagged jet four-momenta.
 * @param[in] str_bjet_btag_2 The name of the column containing the second
 * hardest b-tagged jet b tagging score.
 * @param[in] str_top_p4 The name of the column containing the reconstructed top
 * quark four-momenta.
 * @param[in] str_tb_p4 The name of the column containing the reconstructed b
 * quark (top quark decay) four-momenta.
 * @param[in] str_sb_p4 The name of the column containing the reconstructed b
 * quark (top quark production) four-momenta.
 * @param[in] str_good_jetslowpt_mask The name of the column containing the mask
 * for low pT jets.
 * @param[in] str_jet_pt The name of the column for the jet pT.
 * @param[in] str_jet_eta The name of the column for the jet eta.
 * @param[in] str_jet_phi The name of the column for the jet phi.
 * @param[in] str_jet_mass The name of the column for the jet mass.
 * @param[out] str_dphi_top_b1 The name of the output column for the delta phi
 * between the reconstructed top quark and the hardest b jet.
 * @param[out] str_deta_top_sb The name of the output column for the delta eta
 * between the reconstructed top quark and the jet assigned to b quark of the
 * top quark decay.
 * @param[out] str_dphi_b1_b2 The name of the output column for the delta phi
 * between the two b jets.
 * @param[out] str_deta_lep_b1 The name of the output column for the delta eta
 * between the lepton and the hardest b jet.
 * @param[out] str_m_lep_b2 The name of the output column for the invariant mass
 * of the lepton and the second hardest b jet.
 * @param[out] str_pt_b1_b2 The name of the output column for the transverse
 * momentum of the two combined b jets.
 * @param[out] str_costhetastar The name of the output column for the
 * polarization angle.
 * @param[out] str_sumht The name of the output column for the sum of Ht.
 * @param[out] str_wolfram The name of the output column for the third
 * Fox-Wolfram moment.
 * @param[out] str_deta_topb2_b1 The name of the output column for the delta eta
 * between the top quark reconstructed from the second hardest b jet and the
 * hardest b jet.
 *
 * @return A dataframe containing the new columns.
 */
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
    const std::string &str_deta_topb2_b1) {

    auto df2 =
        df.Define(str_dphi_top_b1,
                  [](const int reco, const ROOT::Math::PtEtaPhiMVector &top,
                     const ROOT::Math::PtEtaPhiMVector &b1) {
                      if (!reco) {
                          return -10.0;
                      }
                      return ROOT::Math::VectorUtil::DeltaPhi(top, b1);
                  },
                  {str_is_reco, str_top_p4, str_bjet_p4_1});

    auto df3 =
        df2.Define(str_deta_top_sb,
                   [](const int reco, const ROOT::Math::PtEtaPhiMVector &top,
                      const ROOT::Math::PtEtaPhiMVector &sb) {
                       if (!reco) {
                           return -10.0;
                       }
                       return abs(top.Eta() - sb.Eta());
                   },
                   {str_is_reco, str_top_p4, str_sb_p4});

    auto df4 = df3.Define(
        str_dphi_b1_b2,
        [](const int reco, const ROOT::Math::PtEtaPhiMVector &b1,
           const ROOT::Math::PtEtaPhiMVector &b2,
           const ROOT::Math::PtEtaPhiMVector &nonb1) {
            if (!reco) {
                return -10.0;
            }
            if (b2.Pt() > 0)
                return ROOT::Math::VectorUtil::DeltaPhi(b1, b2);
            return ROOT::Math::VectorUtil::DeltaPhi(b1, nonb1);
        },
        {str_is_reco, str_bjet_p4_1, str_bjet_p4_2, str_nonbjet_p4_1});

    auto df5 =
        df4.Define(str_deta_lep_b1,
                   [](const int reco, const ROOT::Math::PtEtaPhiMVector &lep,
                      const ROOT::Math::PtEtaPhiMVector &b1) {
                       if (!reco) {
                           return -10.0;
                       }
                       return abs(lep.Eta() - b1.Eta());
                   },
                   {str_is_reco, str_lep_p4, str_bjet_p4_1});

    auto df6 =
        df5.Define(str_m_lep_b2,
                   [](const int reco, const ROOT::Math::PtEtaPhiMVector &lep,
                      const ROOT::Math::PtEtaPhiMVector &b2,
                      const ROOT::Math::PtEtaPhiMVector &nonb1) {
                       if (!reco) {
                           return -10.0;
                       }
                       if (b2.Pt() > 0)
                           return (lep + b2).M();
                       return (lep + nonb1).M();
                   },
                   {str_is_reco, str_lep_p4, str_bjet_p4_2, str_nonbjet_p4_1});

    auto df7 = df6.Define(
        str_pt_b1_b2,
        [](const int reco, const ROOT::Math::PtEtaPhiMVector &b1,
           const ROOT::Math::PtEtaPhiMVector &b2,
           const ROOT::Math::PtEtaPhiMVector &nonb1) {
            if (!reco) {
                return -10.0;
            }
            if (b2.Pt() > 0)
                return (b1 + b2).Pt();
            return (b1 + nonb1).Pt();
        },
        {str_is_reco, str_bjet_p4_1, str_bjet_p4_2, str_nonbjet_p4_1});

    auto df8 = df7.Define(
        str_costhetastar,
        [](const int reco, const ROOT::Math::PtEtaPhiMVector &top,
           const ROOT::Math::PtEtaPhiMVector &sb,
           const ROOT::Math::PtEtaPhiMVector &lep) {
            if (!reco) {
                return -10.0;
            }
            double costhetastar = 0;

            auto top_boost_vec = ROOT::Math::Cartesian3D(
                top.X() / top.T(), top.Y() / top.T(), top.Z() / top.T());
            ROOT::Math::Boost top_boost(top_boost_vec);
            top_boost.Invert();
            ROOT::Math::PtEtaPhiMVector lep_boosted = top_boost(lep);
            ROOT::Math::PtEtaPhiMVector sb_boosted = top_boost(sb);
            costhetastar =
                lep_boosted.Vect().Dot(sb_boosted.Vect()) /
                sqrt(lep_boosted.Vect().Mag2() * sb_boosted.Vect().Mag2());

            Logger::get("DNN_costhetastar")
                ->debug("top_boost {} {} {}", top_boost_vec.X(),
                        top_boost_vec.Y(), top_boost_vec.Z());
            Logger::get("DNN_costhetastar")
                ->debug("sb boosted {} {} {} {}", sb_boosted.Pt(),
                        sb_boosted.Eta(), sb_boosted.Phi(), sb_boosted.M());
            Logger::get("DNN_costhetastar")
                ->debug("lep boosted {} {} {} {}", lep_boosted.Pt(),
                        lep_boosted.Eta(), lep_boosted.Phi(), lep_boosted.M());
            Logger::get("DNN_costhetastar")
                ->debug("costhetastar {}", costhetastar);

            return costhetastar;
        },
        {str_is_reco, str_top_p4, str_sb_p4, str_lep_p4});

    auto df9 = df8.Define(
        str_sumht,
        [](const int reco, const ROOT::Math::PtEtaPhiMVector &b1,
           const ROOT::Math::PtEtaPhiMVector &b2,
           const ROOT::Math::PtEtaPhiMVector &lep,
           const ROOT::Math::PtEtaPhiMVector &met) {
            if (!reco) {
                return -10.0;
            }
            double Ht = lep.Pt() + met.Pt();
            if (b1.Pt() > 0)
                Ht += b1.Pt();
            if (b2.Pt() > 0)
                Ht += b2.Pt();
            return Ht;
        },
        {str_is_reco, str_bjet_p4_1, str_bjet_p4_2, str_lep_p4, str_met_p4});

    auto df10 = df9.Define(
        str_wolfram,
        [](const int reco, const ROOT::RVec<int> &jet_mask,
           const ROOT::RVec<float> &jet_pt, const ROOT::RVec<float> &jet_eta,
           const ROOT::RVec<float> &jet_phi,
           const ROOT::RVec<float> &jet_mass) {
            if (!reco) {
                return -10.0;
            }

            Logger::get("DNN_wolfram")->debug("eval sum_e...");
            Logger::get("DNN_wolfram")->debug("jet mask {}", jet_mask);

            float sum_e = 0;
            for (std::size_t index = 0; index < jet_pt.size(); index++) {
                Logger::get("DNN_wolfram")->debug("checking jet {}", index);
                if (jet_mask.at(index)) {
                    sum_e += (ROOT::Math::PtEtaPhiMVector(
                                  jet_pt.at(index, 0), jet_eta.at(index, 0),
                                  jet_phi.at(index, 0), jet_mass.at(index, 0)))
                                 .E();
                    Logger::get("DNN_wolfram")
                        ->debug("sum_e now {} after adding jet {}", sum_e,
                                index);
                }
            }

            Logger::get("DNN_wolfram")->debug("now momenta");

            double h3 = 0;

            for (std::size_t index = 0; index < jet_pt.size(); index++) {
                for (std::size_t index2 = index + 1; index2 < jet_pt.size();
                     index2++) {
                    if (jet_mask.at(index) && jet_mask.at(index2)) {

                        auto jet_a = ROOT::Math::PtEtaPhiMVector(
                            jet_pt[index], jet_eta[index], jet_phi[index],
                            jet_mass[index]);
                        auto jet_b = ROOT::Math::PtEtaPhiMVector(
                            jet_pt[index2], jet_eta[index2], jet_phi[index2],
                            jet_mass[index2]);

                        auto jet_aXYZT = ROOT::Math::XYZTVector(
                            jet_a.X(), jet_a.Y(), jet_a.Z(), jet_a.T());
                        auto jet_bXYZT = ROOT::Math::XYZTVector(
                            jet_b.X(), jet_b.Y(), jet_b.Z(), jet_b.T());

                        float costh = ROOT::Math::VectorUtil::CosTheta(
                            jet_aXYZT, jet_bXYZT);
                        float p3 = 0.5 * (5.0 * costh * costh - 3.0 * costh);
                        float pipj = jet_aXYZT.P() * jet_bXYZT.P();

                        h3 += (pipj / (sum_e * sum_e)) * p3;
                    }
                }
            }

            return h3;
        },
        {str_is_reco, str_good_jetslowpt_mask, str_jet_pt, str_jet_eta,
         str_jet_phi, str_jet_mass});

    auto df11 =
        df10.Define(str_deta_topb2_b1,
                    [](const int reco, const ROOT::Math::PtEtaPhiMVector &b1,
                       const ROOT::Math::PtEtaPhiMVector &b2,
                       const ROOT::Math::PtEtaPhiMVector &nonb1,
                       const ROOT::Math::PtEtaPhiMVector &wlep) {
                        if (!reco) {
                            return -10.0;
                        }
                        if (b2.Pt() > 0)
                            return abs((b2 + wlep).Eta() - b1.Eta());
                        return abs((nonb1 + wlep).Eta() - b1.Eta());
                    },
                    {str_is_reco, str_bjet_p4_1, str_bjet_p4_2,
                     str_nonbjet_p4_1, str_wlep_p4});

    return df11;
}

void sf_from_root_file(TH2D *h, int nbinsx, int nbinsy, float pt, float eta,
                       int var, double *sf) {
    // double sf_from_root_file(TH2D *h,float pt, float eta, int var, int
    // nbinsx, int nbinsy, int ymax) {

    Logger::get("trigger_func")->debug("complex debug point 1");

    Logger::get("trigger_func")
        ->debug("hist bins {} {}", h->GetXaxis()->GetNbins(),
                h->GetYaxis()->GetNbins());

    int xbin_index = h->GetXaxis()->FindBin(eta);
    int ybin_index = h->GetYaxis()->FindBin(pt);

    // int ybin_max_index = h->GetYaxis()->FindBin(ymax);
    Logger::get("trigger_func")->debug("complex debug point 2");

    if (xbin_index > nbinsx)
        xbin_index = nbinsx;
    if (ybin_index > nbinsy)
        ybin_index = nbinsy;
    // if(ybin_max_index>nbinsy) ybin_max_index=nbinsy;
    Logger::get("trigger_func")->debug("complex debug point 3");

    if (xbin_index == 0)
        xbin_index = 1;
    if (ybin_index == 0)
        ybin_index = 1;

    // int bin_index;
    Logger::get("trigger_func")->debug("complex debug point 4");

    const int bin_index = h->GetBin(xbin_index, ybin_index, 0);
    // const int bin_max_index = h->GetBin(xbin_index,ybin_max_index,0);
    Logger::get("trigger_func")->debug("complex debug point 5");

    *sf = 1.;

    Logger::get("trigger_func")->debug("complex debug point 6");

    if (var == 0)
        *sf = h->GetBinContent(bin_index);
    if (var == +1)
        *sf = h->GetBinContent(bin_index) + h->GetBinError(bin_index);
    if (var == -1)
        *sf = h->GetBinContent(bin_index) - h->GetBinError(bin_index);
}

/**
 * Function to add columns for various different lepton scale factors
 *
 * @param[in] df the input dataframe
 * @param[in] str_lep_pt name of the column containing the lepton pT
 * @param[in] str_lep_eta name of the column containing the lepton eta
 * @param[in] str_lep_sceta name of the column containing the lepton
 * supercluster (electrons only)
 * @param[in] str_lep_is_mu name of the column containing the flag if the event
 * contains a muon
 * @param[in] str_lep_is_el name of the column containing the flag if the event
 * contains an electron
 * @param[in] str_lep_is_iso name of the column containing the flag if the event
 * is (anti-)isolated
 * @param[out] str_lep_sf_mu_trigger_nom name of the output column for the
 * nominal muon trigger scale factor
 * @param[out] str_lep_sf_mu_trigger_up name of the output column for the
 * up-shifted muon trigger scale factor
 * @param[out] str_lep_sf_mu_trigger_down name of the output column for the
 * down-shifted muon trigger scale factor
 * @param[out] str_lep_sf_mu_iso_nom name of the output column for the nominal
 * muon isolation scale factor
 * @param[out] str_lep_sf_mu_iso_up name of the output column for the up-shifted
 * muon isolation scale factor
 * @param[out] str_lep_sf_mu_iso_down name of the output column for the
 * down-shifted muon isolation scale factor
 * @param[out] str_lep_sf_mu_id_nom name of the output column for the nominal
 * muon ID scale factor
 * @param[out] str_lep_sf_mu_id_up name of the output column for the up-shifted
 * muon ID scale factor
 * @param[out] str_lep_sf_mu_id_down name of the output column for the
 * down-shifted muon ID scale factor
 * @param[out] str_lep_sf_el_trigger_nom name of the output column for the
 * nominal electron trigger scale factor
 * @param[out] str_lep_sf_el_trigger_up name of the output column for the
 * up-shifted electron trigger scale factor
 * @param[out] str_lep_sf_el_trigger_down name of the output column for the
 * down-shifted electron trigger scale factor
 * @param[out] str_lep_sf_el_id_nom name of the output column for the nominal
 * electron ID scale factor
 * @param[out] str_lep_sf_el_id_up name of the output column for the up-shifted
 * electron ID scale factor
 * @param[out] str_lep_sf_el_id_down name of the output column for the
 * down-shifted electron ID scale factor
 * @param[out] str_lep_sf_el_reco_nom name of the output column for the nominal
 * electron reconstruction scale factor
 * @param[out] str_lep_sf_el_reco_up name of the output column for the
 * up-shifted electron reconstruction scale factor
 * @param[out] str_lep_sf_el_reco_down name of the output column for the
 * down-shifted electron reconstruction scale factor
 * @param[in] mu_sf_era name tag for the muon scale factors
 * @param[in] mu_trigger_sf_file file name for the nominal muon trigger scale
 * factor file
 * @param[in] mu_trigger_sf_file_syst file name for the shifted muon trigger
 * scale factor file
 * @param[in] mu_trigger_sf_name name of the nominal muon trigger scale factor
 * in the JSON file
 * @param[in] mu_trigger_sf_name_syst name of the shifted muon trigger scale
 * factor in the JSON file
 * @param[in] mu_iso_sf_file file name for the nominal muon isolation scale
 * factor file
 * @param[in] mu_iso_sf_file_syst file name for the shifted muon isolation scale
 * factor file
 * @param[in] mu_iso_sf_name name of the nominal muon isolation scale factor in
 * the JSON file
 * @param[in] mu_iso_sf_name_syst name of the shifted muon isolation scale
 * factor in the JSON file
 * @param[in] mu_sf_file file name for the MUO POG correction lib file
 * @param[in] mu_id_sf_name name of the muon ID scale factors
 * @param[in] el_sf_era name tag for the electron scale factors
 * @param[in] el_trigger_sf_file file name for the nominal electron trigger
 * scale factor file
 * @param[in] el_trigger_sf_file_syst file name for the shifted electron trigger
 * scale factor file
 * @param[in] el_trigger_sf_name name of the nominal electron trigger scale
 * factor in the JSON file
 * @param[in] el_trigger_sf_name_syst name of the shifted electron trigger scale
 * factor in the JSON file
 * @param[in] el_sf_file file name for the EGM POG correction lib file
 * @param[in] el_id_sf_name name of the electron ID scale factors
 *
 * @return a dataframe containing the new columns
 */
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
    const std::string &el_id_sf_name) {

    Logger::get("lepsf_muonTriggerSF")
        ->debug("Setting up functions for muon trigger sf");
    Logger::get("lepsf_muonTriggerSF")
        ->debug("Trigger - File {}", mu_trigger_sf_file);
    Logger::get("lepsf_muonTriggerSF")
        ->debug("Trigger - Name {}", mu_trigger_sf_name);
    auto evaluator_mu_trigger =
        correction::CorrectionSet::from_file(mu_trigger_sf_file)
            ->at(mu_trigger_sf_name);
    auto mu_trigger_nom =
        [evaluator_mu_trigger](const float &pt, const float &eta,
                               const int &is_mu, const int &is_el,
                               const int &is_iso) {
            double sf = 1.;
            if (is_mu == 0 || is_iso != +1)
                return sf;
            Logger::get("lepsf_muonTriggerSF")
                ->debug("Muon - pt {}, eta {}", pt, eta);
            if (pt >= 0.0)
                sf = evaluator_mu_trigger->evaluate({std::abs(eta), pt});
            Logger::get("lepsf_muonTriggerSF")->debug("Trigger - sf {}", sf);
            return sf;
        };
    Logger::get("lepsf_muonTriggerSF")
        ->debug("Trigger syst - File {}", mu_trigger_sf_file_syst);
    Logger::get("lepsf_muonTriggerSF")
        ->debug("Trigger syst - Name {}", mu_trigger_sf_name_syst);
    auto evaluator_mu_trigger_syst =
        correction::CorrectionSet::from_file(mu_trigger_sf_file_syst)
            ->at(mu_trigger_sf_name_syst);
    auto mu_trigger_up = [evaluator_mu_trigger, evaluator_mu_trigger_syst](
                             const float &pt, const float &eta,
                             const int &is_mu, const int &is_el,
                             const int &is_iso) {
        double sf = 1.;
        if (is_mu == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_mu_trigger->evaluate({std::abs(eta), pt}) +
                 evaluator_mu_trigger_syst->evaluate({std::abs(eta), pt});
        Logger::get("lepsf_muonTriggerSF")->debug("Trigger up - sf {}", sf);
        return sf;
    };
    auto mu_trigger_down = [evaluator_mu_trigger, evaluator_mu_trigger_syst](
                               const float &pt, const float &eta,
                               const int &is_mu, const int &is_el,
                               const int &is_iso) {
        double sf = 1.;
        if (is_mu == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_mu_trigger->evaluate({std::abs(eta), pt}) -
                 evaluator_mu_trigger_syst->evaluate({std::abs(eta), pt});
        Logger::get("lepsf_muonTriggerSF")->debug("Trigger down - sf {}", sf);
        return sf;
    };

    Logger::get("lepsf_muonIsoSF")
        ->debug("Setting up functions for muon iso sf");
    Logger::get("lepsf_muonIsoSF")->debug("ISO - File {}", mu_iso_sf_file);
    Logger::get("lepsf_muonIsoSF")->debug("ISO - Name {}", mu_iso_sf_name);
    auto evaluator_mu_iso = correction::CorrectionSet::from_file(mu_iso_sf_file)
                                ->at(mu_iso_sf_name);
    auto mu_iso_nom = [evaluator_mu_iso](const float &pt, const float &eta,
                                         const int &is_mu, const int &is_el,
                                         const int &is_iso) {
        double sf = 1.;
        if (is_mu == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_mu_iso->evaluate({std::abs(eta), pt});
        Logger::get("lepsf_muonIsoSF")->debug("Iso - sf {}", sf);
        return sf;
    };
    Logger::get("lepsf_muonisoSF")
        ->debug("iso syst - File {}", mu_iso_sf_file_syst);
    Logger::get("lepsf_muonisoSF")
        ->debug("iso syst - Name {}", mu_iso_sf_name_syst);
    auto evaluator_mu_iso_syst =
        correction::CorrectionSet::from_file(mu_iso_sf_file_syst)
            ->at(mu_iso_sf_name_syst);
    auto mu_iso_up = [evaluator_mu_iso, evaluator_mu_iso_syst](
                         const float &pt, const float &eta, const int &is_mu,
                         const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_mu == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_mu_iso->evaluate({std::abs(eta), pt}) +
                 evaluator_mu_iso_syst->evaluate({std::abs(eta), pt});
        Logger::get("lepsf_muonIsoSF")->debug("Iso up - sf {}", sf);
        return sf;
    };
    auto mu_iso_down = [evaluator_mu_iso, evaluator_mu_iso_syst](
                           const float &pt, const float &eta, const int &is_mu,
                           const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_mu == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_mu_iso->evaluate({std::abs(eta), pt}) -
                 evaluator_mu_iso_syst->evaluate({std::abs(eta), pt});
        Logger::get("lepsf_muonIsoSF")->debug("Iso down - sf {}", sf);
        return sf;
    };

    Logger::get("lepsf_muonIdSF")->debug("Setting up functions for muon id sf");
    Logger::get("lepsf_muonIdSF")->debug("ID - File {}", mu_sf_file);
    Logger::get("lepsf_muonIdSF")->debug("ID - Name {}", mu_id_sf_name);
    Logger::get("lepsf_muonIdSF")->debug("ID - Era {}", mu_sf_era);
    auto evaluator_mu_id =
        correction::CorrectionSet::from_file(mu_sf_file)->at(mu_id_sf_name);

    auto mu_id_nom = [evaluator_mu_id, mu_sf_era](
                         const float &pt, const float &eta, const int &is_mu,
                         const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_mu == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf =
                evaluator_mu_id->evaluate({mu_sf_era, std::abs(eta), pt, "sf"});
        Logger::get("lepsf_muonIdSF")->debug("ID - sf {}", sf);
        return sf;
    };

    auto mu_id_up = [evaluator_mu_id, mu_sf_era](
                        const float &pt, const float &eta, const int &is_mu,
                        const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_mu == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_mu_id->evaluate(
                {mu_sf_era, std::abs(eta), pt, "systup"});
        Logger::get("lepsf_muonIdSF")->debug("ID up - sf {}", sf);
        return sf;
    };

    auto mu_id_down = [evaluator_mu_id, mu_sf_era](
                          const float &pt, const float &eta, const int &is_mu,
                          const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_mu == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_mu_id->evaluate(
                {mu_sf_era, std::abs(eta), pt, "systdown"});
        Logger::get("lepsf_muonIdSF")->debug("ID down - sf {}", sf);
        return sf;
    };

    auto df1 = df.Define(str_lep_sf_mu_trigger_nom, mu_trigger_nom,
                         {str_lep_pt, str_lep_eta, str_lep_is_mu, str_lep_is_el,
                          str_lep_is_iso});
    auto df2 = df1.Define(str_lep_sf_mu_trigger_up, mu_trigger_up,
                          {str_lep_pt, str_lep_eta, str_lep_is_mu,
                           str_lep_is_el, str_lep_is_iso});
    auto df3 = df2.Define(str_lep_sf_mu_trigger_down, mu_trigger_down,
                          {str_lep_pt, str_lep_eta, str_lep_is_mu,
                           str_lep_is_el, str_lep_is_iso});
    auto df4 = df3.Define(str_lep_sf_mu_iso_nom, mu_iso_nom,
                          {str_lep_pt, str_lep_eta, str_lep_is_mu,
                           str_lep_is_el, str_lep_is_iso});
    auto df5 = df4.Define(str_lep_sf_mu_iso_up, mu_iso_up,
                          {str_lep_pt, str_lep_eta, str_lep_is_mu,
                           str_lep_is_el, str_lep_is_iso});
    auto df6 = df5.Define(str_lep_sf_mu_iso_down, mu_iso_down,
                          {str_lep_pt, str_lep_eta, str_lep_is_mu,
                           str_lep_is_el, str_lep_is_iso});
    auto df7 = df6.Define(str_lep_sf_mu_id_nom, mu_id_nom,
                          {str_lep_pt, str_lep_eta, str_lep_is_mu,
                           str_lep_is_el, str_lep_is_iso});
    auto df8 = df7.Define(str_lep_sf_mu_id_up, mu_id_up,
                          {str_lep_pt, str_lep_eta, str_lep_is_mu,
                           str_lep_is_el, str_lep_is_iso});
    auto df9 = df8.Define(str_lep_sf_mu_id_down, mu_id_down,
                          {str_lep_pt, str_lep_eta, str_lep_is_mu,
                           str_lep_is_el, str_lep_is_iso});

    Logger::get("lepsf_electronTriggerSF")
        ->debug("Setting up functions for electron trigger sf");
    Logger::get("lepsf_electronTriggerSF")
        ->debug("Trigger - File {}", el_trigger_sf_file);
    Logger::get("lepsf_electronTriggerSF")
        ->debug("Trigger - Name {}", el_trigger_sf_name);
    auto evaluator_el_trigger =
        correction::CorrectionSet::from_file(el_trigger_sf_file)
            ->at(el_trigger_sf_name);
    auto el_trigger_nom = [evaluator_el_trigger](
                              const float &pt, const float &sceta,
                              const int &is_mu, const int &is_el,
                              const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        Logger::get("lepsf_electronTriggerSF")
            ->debug("Electron - pt {}, sceta {}", pt, sceta);
        if (pt >= 0.0)
            sf = evaluator_el_trigger->evaluate({sceta, pt});
        Logger::get("lepsf_electronTriggerSF")->debug("Trigger - sf {}", sf);
        return sf;
    };

    auto evaluator_el_trigger_syst =
        correction::CorrectionSet::from_file(el_trigger_sf_file_syst)
            ->at(el_trigger_sf_name_syst);
    auto el_trigger_up = [evaluator_el_trigger, evaluator_el_trigger_syst](
                             const float &pt, const float &sceta,
                             const int &is_mu, const int &is_el,
                             const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_el_trigger->evaluate({sceta, pt}) +
                 evaluator_el_trigger_syst->evaluate({sceta, pt});
        Logger::get("lepsf_electronTriggerSF")->debug("Trigger up - sf {}", sf);
        return sf;
    };
    auto el_trigger_down = [evaluator_el_trigger, evaluator_el_trigger_syst](
                               const float &pt, const float &sceta,
                               const int &is_mu, const int &is_el,
                               const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_el_trigger->evaluate({sceta, pt}) -
                 evaluator_el_trigger_syst->evaluate({sceta, pt});
        Logger::get("lepsf_electronTriggerSF")
            ->debug("Trigger down - sf {}", sf);
        return sf;
    };

    Logger::get("lepsf_electronSF")
        ->debug("Setting up functions for electron id+reco sf");
    Logger::get("lepsf_electronSF")->debug("SF - File {}", el_sf_file);
    Logger::get("lepsf_electronSF")->debug("SF - Name {}", el_id_sf_name);
    Logger::get("lepsf_electronSF")->debug("SF - Era {}", el_sf_era);
    auto evaluator_el_sf =
        correction::CorrectionSet::from_file(el_sf_file)->at(el_id_sf_name);

    auto el_id_nom = [evaluator_el_sf, el_sf_era](
                         const float &pt, const float &sceta, const int &is_mu,
                         const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_el_sf->evaluate(
                {el_sf_era, "sf", "Tight", sceta, pt});
        Logger::get("lepsf_electronIdSF")->debug("ID - sf {}", sf);
        return sf;
    };

    auto el_id_up = [evaluator_el_sf, el_sf_era](
                        const float &pt, const float &sceta, const int &is_mu,
                        const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_el_sf->evaluate(
                {el_sf_era, "sfup", "Tight", sceta, pt});
        Logger::get("lepsf_electronIdSF")->debug("ID up - sf {}", sf);
        return sf;
    };

    auto el_id_down = [evaluator_el_sf, el_sf_era](
                          const float &pt, const float &sceta, const int &is_mu,
                          const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_el_sf->evaluate(
                {el_sf_era, "sfdown", "Tight", sceta, pt});
        Logger::get("lepsf_electronIdSF")->debug("ID down - sf {}", sf);
        return sf;
    };

    auto el_reco_nom = [evaluator_el_sf,
                        el_sf_era](const float &pt, const float &sceta,
                                   const int &is_mu, const int &is_el,
                                   const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_el_sf->evaluate(
                {el_sf_era, "sf", "RecoAbove20", sceta, pt});
        Logger::get("lepsf_electronRecoSF")->debug("RECO - sf {}", sf);
        return sf;
    };

    auto el_reco_up = [evaluator_el_sf, el_sf_era](
                          const float &pt, const float &sceta, const int &is_mu,
                          const int &is_el, const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_el_sf->evaluate(
                {el_sf_era, "sfup", "RecoAbove20", sceta, pt});
        Logger::get("lepsf_electronRecoSF")->debug("RECO up - sf {}", sf);
        return sf;
    };

    auto el_reco_down = [evaluator_el_sf,
                         el_sf_era](const float &pt, const float &sceta,
                                    const int &is_mu, const int &is_el,
                                    const int &is_iso) {
        double sf = 1.;
        if (is_el == 0 || is_iso != +1)
            return sf;
        if (pt >= 0.0)
            sf = evaluator_el_sf->evaluate(
                {el_sf_era, "sfdown", "RecoAbove20", sceta, pt});
        Logger::get("lepsf_electronRecoSF")->debug("RECO down - sf {}", sf);
        return sf;
    };

    auto df10 = df9.Define(str_lep_sf_el_trigger_nom, el_trigger_nom,
                           {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                            str_lep_is_el, str_lep_is_iso});
    auto df11 = df10.Define(str_lep_sf_el_trigger_up, el_trigger_up,
                            {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                             str_lep_is_el, str_lep_is_iso});
    auto df12 = df11.Define(str_lep_sf_el_trigger_down, el_trigger_down,
                            {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                             str_lep_is_el, str_lep_is_iso});

    auto df13 = df12.Define(str_lep_sf_el_id_nom, el_id_nom,
                            {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                             str_lep_is_el, str_lep_is_iso});
    auto df14 = df13.Define(str_lep_sf_el_id_up, el_id_up,
                            {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                             str_lep_is_el, str_lep_is_iso});
    auto df15 = df14.Define(str_lep_sf_el_id_down, el_id_down,
                            {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                             str_lep_is_el, str_lep_is_iso});

    auto df16 = df15.Define(str_lep_sf_el_reco_nom, el_reco_nom,
                            {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                             str_lep_is_el, str_lep_is_iso});
    auto df17 = df16.Define(str_lep_sf_el_reco_up, el_reco_up,
                            {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                             str_lep_is_el, str_lep_is_iso});
    auto df18 = df17.Define(str_lep_sf_el_reco_down, el_reco_down,
                            {str_lep_pt, str_lep_sceta, str_lep_is_mu,
                             str_lep_is_el, str_lep_is_iso});

    return df18;
}

/**
 * Function to add columns for fixedWP b tagging scale factors based on the
 * event category
 *
 * @param[in] df The input dataframe
 * @param[in] str_is_iso Name of the column containing the lepton isolation
 * category
 * @param[in] str_is_reco Name of the column containing the flag if the event is
 * reconstructable
 * @param[out] str_is_jj Name of the output column for the flag if the event
 * falls into the 2j0b category.
 * @param[in] str_is_jjb Name of the column containing the flag if the event
 * falls into the 2j1b category
 * @param[in] str_is_jjbb Name of the column containing the flag if the event
 * falls into the 2j2b category
 * @param[in] str_is_jjjb Name of the column containing the flag if the event
 * falls into the 3j1b category
 * @param[in] str_is_jjjbb Name of the column containing the flag if the event
 * falls into the 3j2b category
 * @param[in] str_nonbjet_pt_1 Name of the column containing the hardest
 * untagged jet pT
 * @param[in] str_nonbjet_eta_1 Name of the column containing the hardest
 * untagged jet eta
 * @param[in] str_nonbjet_btag_1 Name of the column containing the hardest
 * untagged jet b tagging value
 * @param[in] str_nonbjet_flavor_1 Name of the column containing the hardest
 * untagged jet hadron falvor
 * @param[in] str_nonbjet_pt_2 Name of the column containing the second hardest
 * untagged jet pT
 * @param[in] str_nonbjet_eta_2 Name of the column containing the second hardest
 * untagged jet eta
 * @param[in] str_nonbjet_btag_2 Name of the column containing the second
 * hardest untagged jet b tagging value
 * @param[in] str_nonbjet_flavor_2 Name of the column containing the second
 * hardest untagged jet hadron falvor
 * @param[in] str_bjet_pt_1 Name of the column containing the hardest b-tagged
 * jet pT
 * @param[in] str_bjet_eta_1 Name of the column containing the hardest b-tagged
 * jet eta
 * @param[in] str_bjet_btag_1 Name of the column containing the hardest b-tagged
 * jet b tagging value
 * @param[in] str_bjet_flavor_1 Name of the column containing the hardest
 * b-tagged jet hadron falvor
 * @param[in] str_bjet_pt_2 Name of the column containing the second hardest
 * b-tagged jet pT
 * @param[in] str_bjet_eta_2 Name of the column containing the second hardest
 * b-tagged jet eta
 * @param[in] str_bjet_btag_2 Name of the column containing the second hardest
 * b-tagged jet b tagging value
 * @param[in] str_bjet_flavor_2 Name of the column containing the second hardest
 * b-tagged jet hadron falvor
 * @param[out] str_btag_sf_vec Name of the output column for the vector
 * containing all scale factor variations
 * @param[out] str_btagw_nom Name of the output column for the nominal scale
 * factor
 * @param[out] str_btagw_HFup_corr Name of the output column for the HF
 * up-shifted scale factor, correlated fraction among the years
 * @param[out] str_btagw_HFup_uncorr Name of the output column for the HF
 * up-shifted scale factor, uncorrelated fraction among the years
 * @param[out] str_btagw_HFdown_corr Name of the output column for the HF
 * down-shifted scale factor, correlated fraction among the years
 * @param[out] str_btagw_HFdown_uncorr Name of the output column for the HF
 * down-shifted scale factor, uncorrelated fraction among the years
 * @param[out] str_btagw_LFup_corr Name of the output column for the LF
 * up-shifted scale factor, correlated fraction among the years
 * @param[out] str_btagw_LFup_uncorr Name of the output column for the LF
 * up-shifted scale factor, uncorrelated fraction among the years
 * @param[out] str_btagw_LFdown_corr Name of the output column for the LF
 * down-shifted scale factor, correlated fraction among the years
 * @param[out] str_btagw_LFdown_uncorr Name of the output column for the LF
 * down-shifted scale factor, uncorrelated fraction among the years
 * @param[in] btag_sf_file File name for the BTV POG correction lib file
 * @param[in] btag_corr_algo_HF Name of the scale factor method for HF jets
 * @param[in] btag_corr_algo_LF Name of the scale factor method for LF jets
 * @param[in] btag_eff_file File name for the b tagging efficiency maps
 * @param[in] btag_eff_type Type of process for the efficiency map
 * @param[in] btag_wp B tagging working point
 * @param[in] max_bjet_eta_sf Maximum abs eta for b jets
 *
 * @return A dataframe containing the new columns
 */
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
    const std::string &btag_wp, const float &max_bjet_eta_sf) {

    Logger::get("btagsf")->debug("Setting up functions for b tagging sf");
    Logger::get("btagsf")->debug("B tagging SF - File {}", btag_sf_file);
    Logger::get("btagsf")->debug("B tagging SF - Name HF {}",
                                 btag_corr_algo_HF);
    auto evaluator_btag_sf_HF =
        correction::CorrectionSet::from_file(btag_sf_file)
            ->at(btag_corr_algo_HF);
    Logger::get("btagsf")->debug("B tagging SF - Name LF {}",
                                 btag_corr_algo_LF);
    auto evaluator_btag_sf_LF =
        correction::CorrectionSet::from_file(btag_sf_file)
            ->at(btag_corr_algo_LF);
    Logger::get("btagsf")->debug("B tagging EFF - File {}", btag_eff_file);

    std::string btag_eff_name_b = "top_b";
    std::string btag_eff_name_c = "top_c";
    std::string btag_eff_name_udsg = "top_l";
    if (btag_eff_type == "ewk") {
        btag_eff_name_b = "ewk_b";
        btag_eff_name_c = "ewk_c";
        btag_eff_name_udsg = "ewk_l";
    }

    Logger::get("btagsf")->debug("B tagging EFF - Name {} b, {}", btag_eff_type,
                                 btag_eff_name_b);
    auto evaluator_btag_eff_b =
        correction::CorrectionSet::from_file(btag_eff_file)
            ->at(btag_eff_name_b);
    Logger::get("btagsf")->debug("B tagging EFF - Name {} c, {}", btag_eff_type,
                                 btag_eff_name_c);
    auto evaluator_btag_eff_c =
        correction::CorrectionSet::from_file(btag_eff_file)
            ->at(btag_eff_name_c);
    Logger::get("btagsf")->debug("B tagging EFF - Name {} udsg, {}",
                                 btag_eff_type, btag_eff_name_udsg);
    auto evaluator_btag_eff_udsg =
        correction::CorrectionSet::from_file(btag_eff_file)
            ->at(btag_eff_name_udsg);

    // 0 central
    // 1 HFup_correlated
    // 2 HFup_uncorrelated
    // 3 HFdown_correlated
    // 4 HFdown_uncorrelated
    // 5 LFup_correlated
    // 6 LFup_uncorrelated
    // 7 LFdown_correlated
    // 8 LFdown_uncorrelated

    const ROOT::RVec<std::string> shift_HF{
        "central",         "up_correlated",     "up_uncorrelated",
        "down_correlated", "down_uncorrelated", "central",
        "central",         "central",           "central"};

    const ROOT::RVec<std::string> shift_LF{
        "central",         "central",         "central",
        "central",         "central",         "up_correlated",
        "up_uncorrelated", "down_correlated", "down_uncorrelated"};

    auto btag_sf = [evaluator_btag_sf_HF, evaluator_btag_sf_LF,
                    btag_corr_algo_HF, btag_corr_algo_LF, evaluator_btag_eff_b,
                    evaluator_btag_eff_c, evaluator_btag_eff_udsg, btag_wp,
                    max_bjet_eta_sf, shift_HF, shift_LF](
                       const int &is_iso, const int &is_reco, const int &is_jj,
                       const int &is_jjb, const int &is_jjbb,
                       const int &is_jjjb, const int &is_jjjbb,
                       const float &nonbjet_pt_1, const float &nonbjet_eta_1,
                       const float &nonbjet_btag_1, const int &nonbjet_flavor_1,
                       const float &nonbjet_pt_2, const float &nonbjet_eta_2,
                       const float &nonbjet_btag_2, const int &nonbjet_flavor_2,
                       const float &bjet_pt_1, const float &bjet_eta_1,
                       const float &bjet_btag_1, const int &bjet_flavor_1,
                       const float &bjet_pt_2, const float &bjet_eta_2,
                       const float &bjet_btag_2, const int &bjet_flavor_2) {
        unsigned n_vars = shift_HF.size();

        ROOT::RVec<double> sf_vec(n_vars, 1.);

        if ( ( ( is_iso != +1 || is_reco == 0 ) && is_jj != 1 ) || ( is_iso != +1 && is_jj == 1 ) )
            return sf_vec;

        double P_MC = 1.;
        ROOT::RVec<double> P_data(n_vars, 1.);

        ROOT::RVec<double> sf_b1(n_vars, 1.);
        ROOT::RVec<double> sf_b2(n_vars, 1.);
        ROOT::RVec<double> sf_nonb1(n_vars, 1.);
        ROOT::RVec<double> sf_nonb2(n_vars, 1.);

        double eff_b1 = 1.;
        double eff_b2 = 1.;
        double eff_nonb1 = 1.;
        double eff_nonb2 = 1.;

        if (is_jj) {
            if (std::abs(nonbjet_eta_1) < max_bjet_eta_sf) {
                if (nonbjet_flavor_1 == 5) {
                    eff_nonb1 = evaluator_btag_eff_b->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                } else if (nonbjet_flavor_1 == 4) {
                    eff_nonb1 = evaluator_btag_eff_c->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                } else {
                    eff_nonb1 = evaluator_btag_eff_udsg->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_LF->evaluate(
                            {shift_LF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                }

                P_MC *= (1 - eff_nonb1);
                for (unsigned i = 0; i < n_vars; i++) {
                    P_data[i] *= (1 - sf_nonb1[i] * eff_nonb1);
                    Logger::get("btagsf")->debug(
                        "updating P ... eff_nonb1 {}, sf_nonb1 {}, P_MC {}, "
                        "P_data {}",
                        eff_nonb1, sf_nonb1[i], P_MC, P_data[i]);
                }
            }
            if (std::abs(nonbjet_eta_2) < max_bjet_eta_sf) {
                Logger::get("btagsf")->debug(
                    "val 1 {}, val 2 {}, < ??? {}", std::abs(nonbjet_eta_2),
                    max_bjet_eta_sf, std::abs(nonbjet_eta_2) < max_bjet_eta_sf);
                if (nonbjet_flavor_2 == 5) {
                    eff_nonb2 = evaluator_btag_eff_b->evaluate(
                        {nonbjet_eta_2, nonbjet_pt_2});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb2[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_2,
                             std::abs(nonbjet_eta_2), nonbjet_pt_2});
                } else if (nonbjet_flavor_2 == 4) {
                    eff_nonb2 = evaluator_btag_eff_c->evaluate(
                        {nonbjet_eta_2, nonbjet_pt_2});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb2[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_2,
                             std::abs(nonbjet_eta_2), nonbjet_pt_2});
                } else {
                    eff_nonb2 = evaluator_btag_eff_udsg->evaluate(
                        {nonbjet_eta_2, nonbjet_pt_2});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb2[i] = evaluator_btag_sf_LF->evaluate(
                            {shift_LF[i], btag_wp, nonbjet_flavor_2,
                             std::abs(nonbjet_eta_2), nonbjet_pt_2});
                }

                P_MC *= (1 - eff_nonb2);
                for (unsigned i = 0; i < n_vars; i++) {
                    P_data[i] *= (1 - sf_nonb2[i] * eff_nonb2);
                    Logger::get("btagsf")->debug(
                        "updating P ... eff_nonb2 {}, sf_nonb2 {}, P_MC {}, "
                        "P_data {}",
                        eff_nonb2, sf_nonb2[i], P_MC, P_data[i]);
                }
            }
        }

        if (is_jjb) {

            Logger::get("btagsf")->debug(
                "bjet 1 ---> flavor: {}, eta: {}, pt: {}", bjet_flavor_1,
                bjet_eta_1, bjet_pt_1);
            Logger::get("btagsf")->debug(
                "nonbjet 1 ---> flavor: {}, eta: {}, pt: {}", nonbjet_flavor_1,
                nonbjet_eta_1, nonbjet_pt_1);

            if (bjet_flavor_1 == 5) {
                eff_b1 =
                    evaluator_btag_eff_b->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            } else if (bjet_flavor_1 == 4) {
                eff_b1 =
                    evaluator_btag_eff_c->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            } else {
                eff_b1 =
                    evaluator_btag_eff_udsg->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_LF->evaluate(
                        {shift_LF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            }

            P_MC *= eff_b1;
            for (unsigned i = 0; i < n_vars; i++) {
                P_data[i] *= sf_b1[i] * eff_b1;
                Logger::get("btagsf")->debug(
                    "updating P ... eff_b1 {}, sf_b1 {}, P_MC {}, P_data {}",
                    eff_b1, sf_b1[i], P_MC, P_data[i]);
            }

            if (std::abs(nonbjet_eta_1) < max_bjet_eta_sf) {
                if (nonbjet_flavor_1 == 5) {
                    eff_nonb1 = evaluator_btag_eff_b->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                } else if (nonbjet_flavor_1 == 4) {
                    eff_nonb1 = evaluator_btag_eff_c->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                } else {
                    eff_nonb1 = evaluator_btag_eff_udsg->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_LF->evaluate(
                            {shift_LF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                }

                P_MC *= (1 - eff_nonb1);
                for (unsigned i = 0; i < n_vars; i++) {
                    P_data[i] *= (1 - sf_nonb1[i] * eff_nonb1);
                    Logger::get("btagsf")->debug(
                        "updating P ... eff_nonb1 {}, sf_nonb1 {}, P_MC {}, "
                        "P_data {}",
                        eff_nonb1, sf_nonb1[i], P_MC, P_data[i]);
                }
            }
        }

        if (is_jjbb) {

            Logger::get("btagsf")->debug(
                "bjet 1 ---> flavor: {}, eta: {}, pt: {}", bjet_flavor_1,
                bjet_eta_1, bjet_pt_1);
            Logger::get("btagsf")->debug(
                "bjet 2 ---> flavor: {}, eta: {}, pt: {}", bjet_flavor_2,
                bjet_eta_2, bjet_pt_2);

            if (bjet_flavor_1 == 5) {
                eff_b1 =
                    evaluator_btag_eff_b->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            } else if (bjet_flavor_1 == 4) {
                eff_b1 =
                    evaluator_btag_eff_c->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            } else {
                eff_b1 =
                    evaluator_btag_eff_udsg->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_LF->evaluate(
                        {shift_LF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            }

            P_MC *= eff_b1;
            for (unsigned i = 0; i < n_vars; i++) {
                P_data[i] *= sf_b1[i] * eff_b1;
                Logger::get("btagsf")->debug(
                    "updating P ... eff_b1 {}, sf_b1 {}, P_MC {}, P_data {}",
                    eff_b1, sf_b1[i], P_MC, P_data[i]);
            }

            if (bjet_flavor_2 == 5) {
                eff_b2 =
                    evaluator_btag_eff_b->evaluate({bjet_eta_2, bjet_pt_2});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b2[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_2,
                         std::abs(bjet_eta_2), bjet_pt_2});
            } else if (bjet_flavor_2 == 4) {
                eff_b2 =
                    evaluator_btag_eff_c->evaluate({bjet_eta_2, bjet_pt_2});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b2[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_2,
                         std::abs(bjet_eta_2), bjet_pt_2});
            } else {
                eff_b2 =
                    evaluator_btag_eff_udsg->evaluate({bjet_eta_2, bjet_pt_2});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b2[i] = evaluator_btag_sf_LF->evaluate(
                        {shift_LF[i], btag_wp, bjet_flavor_2,
                         std::abs(bjet_eta_2), bjet_pt_2});
            }

            P_MC *= eff_b2;
            for (unsigned i = 0; i < n_vars; i++) {
                P_data[i] *= sf_b2[i] * eff_b2;
                Logger::get("btagsf")->debug(
                    "updating P ... eff_b2 {}, sf_b2 {}, P_MC {}, P_data {}",
                    eff_b2, sf_b2[i], P_MC, P_data[i]);
            }
        }

        if (is_jjjb) {

            Logger::get("btagsf")->debug(
                "bjet 1 ---> flavor: {}, eta: {}, pt: {}", bjet_flavor_1,
                bjet_eta_1, bjet_pt_1);
            Logger::get("btagsf")->debug(
                "nonbjet 1 ---> flavor: {}, eta: {}, pt: {}", nonbjet_flavor_1,
                nonbjet_eta_1, nonbjet_pt_1);
            Logger::get("btagsf")->debug(
                "nonbjet 2 ---> flavor: {}, eta: {}, pt: {}", nonbjet_flavor_2,
                nonbjet_eta_2, nonbjet_pt_2);

            if (bjet_flavor_1 == 5) {
                eff_b1 =
                    evaluator_btag_eff_b->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            } else if (bjet_flavor_1 == 4) {
                eff_b1 =
                    evaluator_btag_eff_c->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            } else {
                eff_b1 =
                    evaluator_btag_eff_udsg->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_LF->evaluate(
                        {shift_LF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            }

            P_MC *= eff_b1;
            for (unsigned i = 0; i < n_vars; i++) {
                P_data[i] *= sf_b1[i] * eff_b1;
                Logger::get("btagsf")->debug(
                    "updating P ... eff_b1 {}, sf_b1 {}, P_MC {}, P_data {}",
                    eff_b1, sf_b1[i], P_MC, P_data[i]);
            }

            if (std::abs(nonbjet_eta_1) < max_bjet_eta_sf) {
                if (nonbjet_flavor_1 == 5) {
                    eff_nonb1 = evaluator_btag_eff_b->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                } else if (nonbjet_flavor_1 == 4) {
                    eff_nonb1 = evaluator_btag_eff_c->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                } else {
                    eff_nonb1 = evaluator_btag_eff_udsg->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_LF->evaluate(
                            {shift_LF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                }

                P_MC *= (1 - eff_nonb1);
                for (unsigned i = 0; i < n_vars; i++) {
                    P_data[i] *= (1 - sf_nonb1[i] * eff_nonb1);
                    Logger::get("btagsf")->debug(
                        "updating P ... eff_nonb1 {}, sf_nonb1 {}, P_MC {}, "
                        "P_data {}",
                        eff_nonb1, sf_nonb1[i], P_MC, P_data[i]);
                }
            }

            if (std::abs(nonbjet_eta_2) < max_bjet_eta_sf) {
                Logger::get("btagsf")->debug(
                    "val 1 {}, val 2 {}, < ??? {}", std::abs(nonbjet_eta_2),
                    max_bjet_eta_sf, std::abs(nonbjet_eta_2) < max_bjet_eta_sf);
                if (nonbjet_flavor_2 == 5) {
                    eff_nonb2 = evaluator_btag_eff_b->evaluate(
                        {nonbjet_eta_2, nonbjet_pt_2});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb2[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_2,
                             std::abs(nonbjet_eta_2), nonbjet_pt_2});
                } else if (nonbjet_flavor_2 == 4) {
                    eff_nonb2 = evaluator_btag_eff_c->evaluate(
                        {nonbjet_eta_2, nonbjet_pt_2});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb2[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_2,
                             std::abs(nonbjet_eta_2), nonbjet_pt_2});
                } else {
                    eff_nonb2 = evaluator_btag_eff_udsg->evaluate(
                        {nonbjet_eta_2, nonbjet_pt_2});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb2[i] = evaluator_btag_sf_LF->evaluate(
                            {shift_LF[i], btag_wp, nonbjet_flavor_2,
                             std::abs(nonbjet_eta_2), nonbjet_pt_2});
                }

                P_MC *= (1 - eff_nonb2);
                for (unsigned i = 0; i < n_vars; i++) {
                    P_data[i] *= (1 - sf_nonb2[i] * eff_nonb2);
                    Logger::get("btagsf")->debug(
                        "updating P ... eff_nonb2 {}, sf_nonb2 {}, P_MC {}, "
                        "P_data {}",
                        eff_nonb2, sf_nonb2[i], P_MC, P_data[i]);
                }
            }
        }

        if (is_jjjbb) {

            Logger::get("btagsf")->debug(
                "bjet 1 ---> flavor: {}, eta: {}, pt: {}", bjet_flavor_1,
                bjet_eta_1, bjet_pt_1);
            Logger::get("btagsf")->debug(
                "bjet 2 ---> flavor: {}, eta: {}, pt: {}", bjet_flavor_2,
                bjet_eta_2, bjet_pt_2);
            Logger::get("btagsf")->debug(
                "nonbjet 1 ---> flavor: {}, eta: {}, pt: {}", nonbjet_flavor_1,
                nonbjet_eta_1, nonbjet_pt_1);

            if (bjet_flavor_1 == 5) {
                eff_b1 =
                    evaluator_btag_eff_b->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            } else if (bjet_flavor_1 == 4) {
                eff_b1 =
                    evaluator_btag_eff_c->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            } else {
                eff_b1 =
                    evaluator_btag_eff_udsg->evaluate({bjet_eta_1, bjet_pt_1});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b1[i] = evaluator_btag_sf_LF->evaluate(
                        {shift_LF[i], btag_wp, bjet_flavor_1,
                         std::abs(bjet_eta_1), bjet_pt_1});
            }

            P_MC *= eff_b1;
            for (unsigned i = 0; i < n_vars; i++) {
                P_data[i] *= sf_b1[i] * eff_b1;
                Logger::get("btagsf")->debug(
                    "updating P ... eff_b1 {}, sf_b1 {}, P_MC {}, P_data {}",
                    eff_b1, sf_b1[i], P_MC, P_data[i]);
            }

            if (bjet_flavor_2 == 5) {
                eff_b2 =
                    evaluator_btag_eff_b->evaluate({bjet_eta_2, bjet_pt_2});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b2[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_2,
                         std::abs(bjet_eta_2), bjet_pt_2});
            } else if (bjet_flavor_2 == 4) {
                eff_b2 =
                    evaluator_btag_eff_c->evaluate({bjet_eta_2, bjet_pt_2});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b2[i] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[i], btag_wp, bjet_flavor_2,
                         std::abs(bjet_eta_2), bjet_pt_2});
            } else {
                eff_b2 =
                    evaluator_btag_eff_udsg->evaluate({bjet_eta_2, bjet_pt_2});
                for (unsigned i = 0; i < n_vars; i++)
                    sf_b2[i] = evaluator_btag_sf_LF->evaluate(
                        {shift_LF[i], btag_wp, bjet_flavor_2,
                         std::abs(bjet_eta_2), bjet_pt_2});
            }

            P_MC *= eff_b2;
            for (unsigned i = 0; i < n_vars; i++) {
                P_data[i] *= sf_b2[i] * eff_b2;
                Logger::get("btagsf")->debug(
                    "updating P ... eff_b2 {}, sf_b2 {}, P_MC {}, P_data {}",
                    eff_b2, sf_b2[i], P_MC, P_data[i]);
            }

            if (std::abs(nonbjet_eta_1) < max_bjet_eta_sf) {
                if (nonbjet_flavor_1 == 5) {
                    eff_nonb1 = evaluator_btag_eff_b->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                } else if (nonbjet_flavor_1 == 4) {
                    eff_nonb1 = evaluator_btag_eff_c->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_HF->evaluate(
                            {shift_HF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                } else {
                    eff_nonb1 = evaluator_btag_eff_udsg->evaluate(
                        {nonbjet_eta_1, nonbjet_pt_1});
                    for (unsigned i = 0; i < n_vars; i++)
                        sf_nonb1[i] = evaluator_btag_sf_LF->evaluate(
                            {shift_LF[i], btag_wp, nonbjet_flavor_1,
                             std::abs(nonbjet_eta_1), nonbjet_pt_1});
                }

                P_MC *= (1 - eff_nonb1);
                for (unsigned i = 0; i < n_vars; i++) {
                    P_data[i] *= (1 - sf_nonb1[i] * eff_nonb1);
                    Logger::get("btagsf")->debug(
                        "updating P ... eff_nonb1 {}, sf_nonb1 {}, P_MC {}, "
                        "P_data {}",
                        eff_nonb1, sf_nonb1[i], P_MC, P_data[i]);
                }
            }
        }

        for (unsigned i = 0; i < n_vars; i++) {
            sf_vec[i] = P_data[i] / P_MC;
            Logger::get("btagsf")->debug("final btag SF (var {}): {}", i,
                                         sf_vec[i]);
        }

        return sf_vec;
    };

    auto df2 = df.Define(
        str_btag_sf_vec, btag_sf,
        {str_is_iso,         str_is_reco,          str_is_jj,
         str_is_jjb,         str_is_jjbb,          str_is_jjjb,
         str_is_jjjbb,       str_nonbjet_pt_1,     str_nonbjet_eta_1,
         str_nonbjet_btag_1, str_nonbjet_flavor_1, str_nonbjet_pt_2,
         str_nonbjet_eta_2,  str_nonbjet_btag_2,   str_nonbjet_flavor_2,
         str_bjet_pt_1,      str_bjet_eta_1,       str_bjet_btag_1,
         str_bjet_flavor_1,  str_bjet_pt_2,        str_bjet_eta_2,
         str_bjet_btag_2,    str_bjet_flavor_2});

    auto df3 =
        df2.Define(str_btagw_nom,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[0]; },
                   {str_btag_sf_vec});

    auto df4 =
        df3.Define(str_btagw_HFup_corr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[1]; },
                   {str_btag_sf_vec});

    auto df5 =
        df4.Define(str_btagw_HFup_uncorr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[2]; },
                   {str_btag_sf_vec});

    auto df6 =
        df5.Define(str_btagw_HFdown_corr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[3]; },
                   {str_btag_sf_vec});

    auto df7 =
        df6.Define(str_btagw_HFdown_uncorr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[4]; },
                   {str_btag_sf_vec});

    auto df8 =
        df7.Define(str_btagw_LFup_corr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[5]; },
                   {str_btag_sf_vec});

    auto df9 =
        df8.Define(str_btagw_LFup_uncorr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[6]; },
                   {str_btag_sf_vec});

    auto df10 =
        df9.Define(str_btagw_LFdown_corr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[7]; },
                   {str_btag_sf_vec});

    auto df11 =
        df10.Define(str_btagw_LFdown_uncorr,
                    [](const ROOT::RVec<double> &sf_vec) { return sf_vec[8]; },
                    {str_btag_sf_vec});

    return df11;
}

/**
 * Function to add columns for fixedWP b tagging scale factors based on a
 * variable number of jets and b-tagged jets.
 *
 * @param[in] df The input dataframe.
 * @param[in] str_jet_pt Name of the column containing the vector of jet pT.
 * @param[in] str_jet_eta Name of the column containing the vector of jet eta.
 * @param[in] str_jet_btag Name of the column containing the vector of jet b
 * tagging value.
 * @param[in] str_jet_flavor Name of the column containing the vector of jet
 * hadron flavor.
 * @param[in] str_jet_collection Name of the column containing the jet
 * collection.
 * @param[out] str_btag_sf_vec Name of the output column for the vector
 * containing all scale factor variations.
 * @param[out] str_btagw_nom Name of the output column for the nominal scale
 * factor.
 * @param[out] str_btagw_HFup_corr Name of the output column for the HF
 * up-shifted scale factor, correlated fraction among the years.
 * @param[out] str_btagw_HFup_uncorr Name of the output column for the HF
 * up-shifted scale factor, uncorrelated fraction among the years.
 * @param[out] str_btagw_HFdown_corr Name of the output column for the HF
 * down-shifted scale factor, correlated fraction among the years.
 * @param[out] str_btagw_HFdown_uncorr Name of the output column for the HF
 * down-shifted scale factor, uncorrelated fraction among the years.
 * @param[out] str_btagw_LFup_corr Name of the output column for the LF
 * up-shifted scale factor, correlated fraction among the years.
 * @param[out] str_btagw_LFup_uncorr Name of the output column for the LF
 * up-shifted scale factor, uncorrelated fraction among the years.
 * @param[out] str_btagw_LFdown_corr Name of the output column for the LF
 * down-shifted scale factor, correlated fraction among the years.
 * @param[out] str_btagw_LFdown_uncorr Name of the output column for the LF
 * down-shifted scale factor, uncorrelated fraction among the years.
 * @param[in] btag_sf_file File name for the BTV POG correction lib file.
 * @param[in] btag_corr_algo_HF Name of the scale factor method for HF jets.
 * @param[in] btag_corr_algo_LF Name of the scale factor method for LF jets.
 * @param[in] btag_eff_file File name for the b tagging efficiency maps.
 * @param[in] btag_eff_type Type of process for the efficiency map.
 * @param[in] btag_wp B tagging working point name.
 * @param[in] btag_cut B tagging cut value.
 * @param[in] max_bjet_eta_sf Maximum abs eta for b jets.
 *
 * @return A dataframe containing the new columns.
 */
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
    const float &max_bjet_eta_sf) {

    Logger::get("btagsf")->debug("Setting up functions for b tagging sf");
    Logger::get("btagsf")->debug("B tagging SF - File {}", btag_sf_file);
    Logger::get("btagsf")->debug("B tagging SF - Name HF {}",
                                 btag_corr_algo_HF);
    auto evaluator_btag_sf_HF =
        correction::CorrectionSet::from_file(btag_sf_file)
            ->at(btag_corr_algo_HF);
    Logger::get("btagsf")->debug("B tagging SF - Name LF {}",
                                 btag_corr_algo_LF);
    auto evaluator_btag_sf_LF =
        correction::CorrectionSet::from_file(btag_sf_file)
            ->at(btag_corr_algo_LF);
    Logger::get("btagsf")->debug("B tagging EFF - File {}", btag_eff_file);

    std::string btag_eff_name_b = "top_b";
    std::string btag_eff_name_c = "top_c";
    std::string btag_eff_name_udsg = "top_l";
    if (btag_eff_type == "ewk") {
        btag_eff_name_b = "ewk_b";
        btag_eff_name_c = "ewk_c";
        btag_eff_name_udsg = "ewk_l";
    }

    Logger::get("btagsf")->debug("B tagging EFF - Name {} b, {}", btag_eff_type,
                                 btag_eff_name_b);
    auto evaluator_btag_eff_b =
        correction::CorrectionSet::from_file(btag_eff_file)
            ->at(btag_eff_name_b);
    Logger::get("btagsf")->debug("B tagging EFF - Name {} c, {}", btag_eff_type,
                                 btag_eff_name_c);
    auto evaluator_btag_eff_c =
        correction::CorrectionSet::from_file(btag_eff_file)
            ->at(btag_eff_name_c);
    Logger::get("btagsf")->debug("B tagging EFF - Name {} udsg, {}",
                                 btag_eff_type, btag_eff_name_udsg);
    auto evaluator_btag_eff_udsg =
        correction::CorrectionSet::from_file(btag_eff_file)
            ->at(btag_eff_name_udsg);

    // 0 central
    // 1 HFup_correlated
    // 2 HFup_uncorrelated
    // 3 HFdown_correlated
    // 4 HFdown_uncorrelated
    // 5 LFup_correlated
    // 6 LFup_uncorrelated
    // 7 LFdown_correlated
    // 8 LFdown_uncorrelated

    const ROOT::RVec<std::string> shift_HF{
        "central",         "up_correlated",     "up_uncorrelated",
        "down_correlated", "down_uncorrelated", "central",
        "central",         "central",           "central"};

    const ROOT::RVec<std::string> shift_LF{
        "central",         "central",         "central",
        "central",         "central",         "up_correlated",
        "up_uncorrelated", "down_correlated", "down_uncorrelated"};

    auto btag_sf = [evaluator_btag_sf_HF, evaluator_btag_sf_LF,
                    btag_corr_algo_HF, btag_corr_algo_LF, evaluator_btag_eff_b,
                    evaluator_btag_eff_c, evaluator_btag_eff_udsg, btag_wp,
                    btag_cut, max_bjet_eta_sf, shift_HF,
                    shift_LF](const ROOT::RVec<float> &jet_pt,
                              const ROOT::RVec<float> &jet_eta,
                              const ROOT::RVec<float> &jet_btag,
                              const ROOT::RVec<int> &jet_flavor,
                              const ROOT::RVec<int> &jet_collection) {
        unsigned n_vars = shift_HF.size();
        unsigned n_jets = jet_collection.size();

        ROOT::RVec<double> sf_vec(n_vars, 1.);

        double P_MC = 1.;
        ROOT::RVec<double> P_data(n_vars, 1.);

        ROOT::RVec<double> eff(n_jets, 1.);
        ROOT::RVec<ROOT::RVec<double>> sf(n_jets,
                                          ROOT::RVec<double>(n_vars, 1.));

        Logger::get("btagsf")->debug("jet collection: {}", jet_collection);
        Logger::get("btagsf")->debug("default sf vec: {}", sf);

        for (unsigned j = 0; j < n_jets; j++) {

            if (std::abs(jet_eta[jet_collection.at(j)]) >= max_bjet_eta_sf)
                continue;

            Logger::get("btagsf")->debug("looking up jet index: {}", j);
            Logger::get("btagsf")->debug("pt, eta, btag, flavor: {} {} {} {}",
                                         jet_pt[jet_collection.at(j)],
                                         jet_eta[jet_collection.at(j)],
                                         jet_btag[jet_collection.at(j)],
                                         jet_flavor[jet_collection.at(j)]);

            if (jet_flavor[jet_collection.at(j)] == 5) {
                eff[j] = evaluator_btag_eff_b->evaluate(
                    {jet_eta[jet_collection.at(j)],
                     jet_pt[jet_collection.at(j)]});
                for (unsigned v = 0; v < n_vars; v++)
                    sf[j][v] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[v], btag_wp, jet_flavor[jet_collection.at(j)],
                         std::abs(jet_eta[jet_collection.at(j)]),
                         jet_pt[jet_collection.at(j)]});
            } else if (jet_flavor[jet_collection.at(j)] == 4) {
                eff[j] = evaluator_btag_eff_c->evaluate(
                    {jet_eta[jet_collection.at(j)],
                     jet_pt[jet_collection.at(j)]});
                for (unsigned v = 0; v < n_vars; v++)
                    sf[j][v] = evaluator_btag_sf_HF->evaluate(
                        {shift_HF[v], btag_wp, jet_flavor[jet_collection.at(j)],
                         std::abs(jet_eta[jet_collection.at(j)]),
                         jet_pt[jet_collection.at(j)]});
            } else {
                eff[j] = evaluator_btag_eff_udsg->evaluate(
                    {jet_eta[jet_collection.at(j)],
                     jet_pt[jet_collection.at(j)]});
                for (unsigned v = 0; v < n_vars; v++)
                    sf[j][v] = evaluator_btag_sf_LF->evaluate(
                        {shift_LF[v], btag_wp, jet_flavor[jet_collection.at(j)],
                         std::abs(jet_eta[jet_collection.at(j)]),
                         jet_pt[jet_collection.at(j)]});
            }

            if (jet_btag[jet_collection.at(j)] >= btag_cut) {
                P_MC *= eff[j];
                for (unsigned v = 0; v < n_vars; v++)
                    P_data[v] *= sf[j][v] * eff[j];
            } else {
                P_MC *= (1 - eff[j]);
                for (unsigned v = 0; v < n_vars; v++)
                    P_data[v] *= (1 - sf[j][v] * eff[j]);
            }

        } // end jet loop

        for (unsigned v = 0; v < n_vars; v++) {
            sf_vec[v] = P_data[v] / P_MC;
            Logger::get("btagsf")->debug("final btag SF (var {}): {}", v,
                                         sf_vec[v]);
        }

        return sf_vec;
    };

    auto df2 = df.Define(str_btag_sf_vec, btag_sf,
                         {str_jet_pt, str_jet_eta, str_jet_btag, str_jet_flavor,
                          str_jet_collection});

    auto df3 =
        df2.Define(str_btagw_nom,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[0]; },
                   {str_btag_sf_vec});

    auto df4 =
        df3.Define(str_btagw_HFup_corr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[1]; },
                   {str_btag_sf_vec});

    auto df5 =
        df4.Define(str_btagw_HFup_uncorr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[2]; },
                   {str_btag_sf_vec});

    auto df6 =
        df5.Define(str_btagw_HFdown_corr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[3]; },
                   {str_btag_sf_vec});

    auto df7 =
        df6.Define(str_btagw_HFdown_uncorr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[4]; },
                   {str_btag_sf_vec});

    auto df8 =
        df7.Define(str_btagw_LFup_corr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[5]; },
                   {str_btag_sf_vec});

    auto df9 =
        df8.Define(str_btagw_LFup_uncorr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[6]; },
                   {str_btag_sf_vec});

    auto df10 =
        df9.Define(str_btagw_LFdown_corr,
                   [](const ROOT::RVec<double> &sf_vec) { return sf_vec[7]; },
                   {str_btag_sf_vec});

    auto df11 =
        df10.Define(str_btagw_LFdown_uncorr,
                    [](const ROOT::RVec<double> &sf_vec) { return sf_vec[8]; },
                    {str_btag_sf_vec});

    return df11;
}

/**
 * Function to combine lepton flavor dependent DNN output into a single column
 *
 * @param[in] df The input dataframe.
 * @param[in] str_lep_is_mu name of the column containing the flag if the event
 * contains a muon
 * @param[in] str_lep_is_el name of the column containing the flag if the event
 * contains an electron
 * @param[in] str_lep_is_iso name of the column containing the flag if the event
 * is (anti-)isolated
 * @param[in] str_is_reco The name of the column containing the flag if a top
 * quark was reconstructed.
 * @param[in] str_dnn_mu name of the column containing DNN output for muons
 * @param[in] str_dnn_el name of the column containing DNN output for electrons
 *
 * @param[out] str_dnn Name of the output column for the final DNN output column
 *
 * @return A dataframe containing the new column.
 */
ROOT::RDF::RNode
CombineDNNOutputs(ROOT::RDF::RNode df, const std::string &str_lep_is_mu,
                  const std::string &str_lep_is_el,
                  const std::string &str_lep_is_iso,
                  const std::string &str_is_reco, const std::string &str_dnn_mu,
                  const std::string &str_dnn_el, const std::string &str_dnn) {

    auto DNN_output = [](const int is_mu, const int is_el, const int is_iso,
                         const int is_reco, const double dnn_mu,
                         const double dnn_el) {
        double output = -10;

        if (!is_reco or is_iso == 0)
            return output;

        if (is_mu)
            output = dnn_mu;

        if (is_el)
            output = dnn_el;

        return output;
    };

    auto df2 = df.Define(str_dnn, DNN_output,
                         {str_lep_is_mu, str_lep_is_el, str_lep_is_iso,
                          str_is_reco, str_dnn_mu, str_dnn_el});

    return df2;
}

} // end namespace topreco

#endif /* GUARD_TOPRECO_H */
