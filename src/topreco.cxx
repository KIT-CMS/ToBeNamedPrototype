#ifndef GUARD_TOPRECO_H
#define GUARD_TOPRECO_H

#include "../include/topreco.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/utility/utility.hxx"
#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include "TMinuit.h"


ROOT::RDF::RNode LeptonSelection(ROOT::RDF::RNode df,
				 const std::string &str_n_loose_mu,
				 const std::string &str_n_loose_el,
				 const std::string &str_n_tight_mu,
				 const std::string &str_n_tight_el,
				 const std::string &str_tight_muons_mask,
				 const std::string &str_tight_electrons_mask,
				 const std::string &str_mu_pt,
				 const std::string &str_mu_eta,
				 const std::string &str_mu_phi,
				 const std::string &str_mu_mass,
				 const std::string &str_el_pt,
				 const std::string &str_el_eta,
				 const std::string &str_el_phi,
				 const std::string &str_el_mass,
				 const std::string &str_n_loose_lep,
				 const std::string &str_n_tight_lep,
				 const std::string &str_is_mu,
				 const std::string &str_is_el,
				 const std::string &str_lep_p4
				 ) {

  auto df1 = df.Define(str_n_loose_lep,
		       [](const int &n_loose_mu,
			  const int &n_loose_el) {
			 Logger::get("lepsel")->debug("size of n_loose_mu and n_loose_el: {} {}",
						      n_loose_mu, n_loose_el);
			 return n_loose_mu + n_loose_el;
		       },
		       {str_n_loose_mu, str_n_loose_el}
		       );

  auto df2 = df1.Define(str_n_tight_lep,
			[](const int &n_tight_mu,
			   const int &n_tight_el) {
			  Logger::get("lepsel")->debug("size of n_tight_mu and n_tight_el: {} {}",
						       n_tight_mu, n_tight_el);
			  return n_tight_mu + n_tight_el;
			},
			{str_n_tight_mu, str_n_tight_el}
			);


  auto df3 = df2.Filter([](const int nloose,
			   const int ntight) {
			  return (nloose == 1) && (ntight == 1);
			},
			{str_n_loose_lep, str_n_tight_lep},
			"lepton selection (exactly one tight lepton)"
			);

  auto df4 = df3.Define(str_is_mu,
			[](const int &n_tight_mu) {
			  return n_tight_mu;
			},
			{str_n_tight_mu}
			);

  auto df5 = df4.Define(str_is_el,
			[](const int &n_tight_el) {
			  return n_tight_el;
			},
			{str_n_tight_el}
			);

  auto lep_p4 = [](const int is_mu,
		   const int is_el,
		   const ROOT::RVec<int> &tight_muons_mask,
		   const ROOT::RVec<int> &tight_electrons_mask,
		   const ROOT::RVec<float> &mu_pt,
		   const ROOT::RVec<float> &mu_eta,
		   const ROOT::RVec<float> &mu_phi,
		   const ROOT::RVec<float> &mu_mass,
		   const ROOT::RVec<float> &el_pt,
		   const ROOT::RVec<float> &el_eta,
		   const ROOT::RVec<float> &el_phi,
		   const ROOT::RVec<float> &el_mass) {

    Logger::get("lep_p4")->debug("masks mu {}, el {}",
    				 tight_muons_mask, tight_electrons_mask);
    Logger::get("lep_p4")->debug("mask sizes mu {}, el {}",
    				 tight_muons_mask.size(), tight_electrons_mask.size());
    Logger::get("lep_p4")->debug("max in mu {}, el {}",
    				 ROOT::VecOps::Max(tight_muons_mask), ROOT::VecOps::Max(tight_electrons_mask));
    Logger::get("lep_p4")->debug("index of max in mu {}, el {}",
    				 ROOT::VecOps::ArgMax(tight_muons_mask), ROOT::VecOps::ArgMax(tight_electrons_mask));

    ROOT::Math::PtEtaPhiMVector lep;

    if (is_mu) {
      Logger::get("lep_p4")->debug("---> should reco mu...");
      lep = ROOT::Math::PtEtaPhiMVector(mu_pt.at(ROOT::VecOps::ArgMax(tight_muons_mask),2),
					mu_eta.at(ROOT::VecOps::ArgMax(tight_muons_mask),2),
					mu_phi.at(ROOT::VecOps::ArgMax(tight_muons_mask),2),
					mu_mass.at(ROOT::VecOps::ArgMax(tight_muons_mask),2));
    }
    else if (is_el) {
      Logger::get("lep_p4")->debug("---> should reco el...");
      lep = ROOT::Math::PtEtaPhiMVector(el_pt.at(ROOT::VecOps::ArgMax(tight_electrons_mask),3),
					el_eta.at(ROOT::VecOps::ArgMax(tight_electrons_mask),3),
					el_phi.at(ROOT::VecOps::ArgMax(tight_electrons_mask),3),
					el_mass.at(ROOT::VecOps::ArgMax(tight_electrons_mask),3));
    }
    else {
      lep = ROOT::Math::PtEtaPhiMVector(1,1,1,1);
    }

    Logger::get("final_lep")->debug("building p4 from lepton with {} {} {} {}",
				    lep.Pt(), lep.Eta(), lep.Phi(), lep.M());

    return lep;

  };

  auto df6 = df5.Define(str_lep_p4,
			lep_p4,
			{str_is_mu, str_is_el,
			    str_tight_muons_mask, str_tight_electrons_mask,
			    str_mu_pt, str_mu_eta, str_mu_phi, str_mu_mass,
			    str_el_pt, str_el_eta, str_el_phi, str_el_mass}
			);

  // ROOT::Math::PtEtaPhiMVector lep_p4 = ROOT::Math::PtEtaPhiMVector(0,0,0,0);
  // lep_p4.SetPtEtaPhiM(0,0,0,0);



  // float lep_pt = 0;
  // float lep_eta = 0;
  // float lep_phi = 0;
  // float lep_mass = 0;



  // auto df3 = df2.Define(str_lep_p4,
  // 			lep_p4,
  // 			{lep_pt, lep_eta, lep_phi, lep_mass}
  // 			);

			// [](ROOT::Math::PtEtaPhiMVector lep_p4) {
			//   return lep_p4;
			// },
			// [](const float &pt,
			//    const float &eta,
			//    const float &phi,
			//    const float &mass) {
			//   auto lep_p4 = ROOT::Math::PtEtaPhiMVector(pt,eta,phi,mass);
			//   return lep_p4;
			// },


  //


  return df6;
}



const float W_MASS = 80.377; // PDG value as of 10/22

// helper function for minimizer constraint
double rad_py(double x, double lep_px) {
  return W_MASS*W_MASS + 4*lep_px*x;
}

// the delta plus function with the py nu plus solution
double min_fplus(double *par) {
  // par[0] = x, par[1] = lep_px, par[2] = lep_py, par[3] = lep_pt, par[4] = px_miss, par[5] = py_miss
  double r = rad_py(par[0],par[1]);
  double y = 0;
  //double res = 99999;
  //if (r>=0) {
  y = (W_MASS*W_MASS * par[2] + 2 * par[1] * par[2] * par[0] + W_MASS * par[3] * sqrt(r))/(2 * par[1]*par[1]);
  double res = sqrt((par[0]-par[4])*(par[0]-par[4]) + (y-par[5])*(y-par[5]));
  // }
  // else // FIXME: proper constraint in TMinuit?
  // res = 99999;
  return res;
}

// the delta minus function with the py nu minus solution
double min_fminus(double *par) {
  // par[0] = x, par[1] = lep_px, par[2] = lep_py, par[3] = lep_pt, par[4] = px_miss, par[5] = py_miss
  double r = rad_py(par[0],par[1]);
  double y = 0;
  double res = 99999;
  if (r>=0) {
    y = (W_MASS*W_MASS * par[2] + 2 * par[1] * par[2] * par[0] - W_MASS * par[3] * sqrt(r))/(2 * par[1]*par[1]);
    res = sqrt((par[0]-par[4])*(par[0]-par[4]) + (y-par[5])*(y-par[5]));
  }
  else
    res = 99999;
  return res;
}


// TMinuit fit function for the py nu plus solution
void fcn_plus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
  f = min_fplus(par);
}

// TMinuit fit function for the py nu minus solution
void fcn_minus(Int_t &npar, Double_t *gin, Double_t &f, Double_t *par, Int_t iflag) {
f = min_fminus(par);
}

ROOT::RDF::RNode ReconstructLeptonicW(ROOT::RDF::RNode df,
				      const std::string &str_lep_p4,
				      const std::string &str_met_p4,
				      const std::string &str_wlep_p4
				      ) {

  auto leptonicW = [](const ROOT::Math::PtEtaPhiMVector lep_p4,
		      const ROOT::Math::PtEtaPhiMVector met_p4) {

    double lep_e  = lep_p4.E();
    double lep_pt = lep_p4.Pt();
    double lep_px = lep_p4.Px();
    double lep_py = lep_p4.Py();
    double lep_pz = lep_p4.Pz();

    ROOT::Math::PtEtaPhiMVector nu_p4 = met_p4;
    double nu_e = met_p4.Pt();
    double nu_px = met_p4.Px();
    double nu_py = met_p4.Py();

    ROOT::Math::PtEtaPhiMVector wlep_p4;

    bool solution_is_real;


    Logger::get("wlep")->debug("building wlep p4 from lepton with E: {} px: {} py: {} pz: {}",
			       lep_e, lep_px, lep_py, lep_pz);
    Logger::get("wlep")->debug("building wlep p4 from pTmiss with E: {} px: {} py: {} pz: ???",
			       nu_e, nu_px, nu_py);

    // definition of the constant mu in Eq. 4.5 (here called alpha to not confuse mu and nu)
    // also converting p_T and cos dphi into p_x and p_y
    double alpha = (W_MASS*W_MASS)/2 + (lep_px*nu_px) + (lep_py*nu_py);

    // for p_z,nu there is a quadratic equation with up to two solutions as shown in Eq. 4.6 and A.7
    // (NOTE: there is a 'power of two' missing in the first denominator of Eq. 4.6)
    // first, check if we have complex solution, i.e. if the radicand is negative
    double rad = ((alpha*alpha * lep_pz*lep_pz)/(lep_pt*lep_pt*lep_pt*lep_pt)) - ((lep_e*lep_e * nu_e*nu_e - alpha*alpha)/(lep_pt*lep_pt));

    if(rad < 0){
      // complex solutions, in around 30% of all cases
      //    cout << "Complex neutrino p_z" << endl;

      // assumption: p_T^miss does not directly correspond to p_T,nu
      // set m_T^W to m^W, result is a quadratic equation for p_(x,y) that depends on p_(y,x)

      // save p_x^miss and p_y^miss as we need them later to determine the better solution

      Logger::get("wlep")->debug("complex solution");

      double px_miss = nu_px;
      double py_miss = nu_py;

      Logger::get("wlep")->debug("complex debug point 1");

      // initialize TMinuit with a maximum of 6 params for py nu plus and minus solution
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

      // set strategy (0 is less accurate, 1 is default, 2 is more accurate)
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
	lower = -W_MASS*W_MASS/(4*lep_px) + 1e-5;
	start_val = lower + 1;
	//      cout << "lower: " << lower << endl;
      }
      if (lep_px < 0) {
	upper = -W_MASS*W_MASS/(4*lep_px) - 1e-5;
	start_val = upper - 1;
	//      cout << "upper: " << upper << endl;
      }
      gMinuit_plus->mnparm(0, "x", start_val, step, lower, upper, ierflg);
      gMinuit_plus->mnparm(1, "lep_px", lep_px, step, -999, 999, ierflg);
      gMinuit_plus->mnparm(2, "lep_py", lep_py, step, -999, 999, ierflg);
      gMinuit_plus->mnparm(3, "lep_pt", lep_pt, step, 0, 999, ierflg);
      gMinuit_plus->mnparm(4, "px_miss", px_miss, step, -999, 999, ierflg);
      gMinuit_plus->mnparm(5, "py_miss", py_miss, step, -999, 999, ierflg);
      gMinuit_plus->FixParameter(1);
      gMinuit_plus->FixParameter(2);
      gMinuit_plus->FixParameter(3);
      gMinuit_plus->FixParameter(4);
      gMinuit_plus->FixParameter(5);

      gMinuit_minus->mnparm(0, "x", start_val, step, lower, upper, ierflg);
      gMinuit_minus->mnparm(1, "lep_px", lep_px, step, -999, 999, ierflg);
      gMinuit_minus->mnparm(2, "lep_py", lep_py, step, -999, 999, ierflg);
      gMinuit_minus->mnparm(3, "lep_pt", lep_pt, step, 0, 999, ierflg);
      gMinuit_minus->mnparm(4, "px_miss", px_miss, step, -999, 999, ierflg);
      gMinuit_minus->mnparm(5, "py_miss", py_miss, step, -999, 999, ierflg);
      gMinuit_minus->FixParameter(1);
      gMinuit_minus->FixParameter(2);
      gMinuit_minus->FixParameter(3);
      gMinuit_minus->FixParameter(4);
      gMinuit_minus->FixParameter(5);

      Logger::get("wlep")->debug("complex debug point 6");

      // now ready for minimization step
      arglist[0] = 500000; // maximum number of iterations
      arglist[1] = 1.; // related to errors
      gMinuit_plus->mnexcm("MIGARD", arglist, 2, ierflg);
      gMinuit_minus->mnexcm("MIGRAD", arglist, 2, ierflg);

      // obtain fit results and calculate values of delta minus and delta plus functions
      // choose solution that leads to a smaller delta value
      double x_plus, x_pluserr;
      double d_plus;
      gMinuit_plus->GetParameter(0,x_plus, x_pluserr);
      double par_plus[6] = {x_plus,lep_px,lep_py,lep_pt,px_miss,py_miss};
      d_plus = min_fplus(par_plus);
      //    cout << "Fit result plus: x=" << x_plus << " " << "d(x)=" << d_plus << endl;

      double x_minus, x_minuserr;
      double d_minus;
      gMinuit_minus->GetParameter(0,x_minus, x_minuserr);
      double par_minus[6] = {x_minus,lep_px,lep_py,lep_pt,px_miss,py_miss};
      d_minus = min_fminus(par_minus);
      //    cout << "Fit result minus: x=" << x_minus << " d(x)=" << d_minus << endl;

      Logger::get("wlep")->debug("complex debug point 7");

      double nu_pxnew, nu_pynew, r_new;
      if (d_plus<d_minus){
	nu_pxnew = x_plus;
	r_new = rad_py(nu_pxnew,lep_px);
	nu_pynew = (W_MASS*W_MASS * lep_py + 2 * lep_px * lep_py * nu_pxnew + W_MASS * lep_pt * sqrt(r_new))/(2 * lep_px * lep_px);
      }
      else{
	nu_pxnew = x_minus;
	r_new = rad_py(nu_pxnew,lep_px);
	nu_pynew = (W_MASS*W_MASS * lep_py + 2 * lep_px * lep_py * nu_pxnew - W_MASS * lep_pt * sqrt(r_new))/(2 * lep_px * lep_px);
      }
      // calculate new nu pz (only one solution with fixed px and py)
      double nu_pznew = lep_pz / (lep_pt*lep_pt) * ((W_MASS*W_MASS / 2) + (lep_px * nu_pxnew) + (lep_py * nu_pynew));
      //    cout << "new nu px: " << nu_pxnew << ", new nu py: " << nu_pynew << ", new nu pz: " << nu_pznew << endl;

      // set 4 momenta of neutrino and W boson
      nu_p4.SetPxPyPzE(nu_pxnew,nu_pynew,nu_pznew,sqrt(nu_pxnew*nu_pxnew + nu_pynew*nu_pynew + nu_pznew*nu_pznew));

      Logger::get("wlep")->debug("complex debug point 8");

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
      nu_p4.SetPxPyPzE(nu_px, nu_py, nu_pz, sqrt(nu_px * nu_px + nu_py * nu_py + nu_pz * nu_pz));
      solution_is_real = true;
    }

    wlep_p4 = lep_p4 + nu_p4;

    Logger::get("wlep")->debug("final leptonic W boson pT: {} eta: {} phi: {} mass: {}",
			       wlep_p4.Pt(), wlep_p4.Eta(), wlep_p4.Phi(), wlep_p4.M());

    return wlep_p4;

  };



  return df.Define(str_wlep_p4,
		   leptonicW,
		   {str_lep_p4, str_met_p4}
		   );

}



ROOT::RDF::RNode ReconstructLeptonicW_mt(ROOT::RDF::RNode df,
					 const std::string &outputname,
					 const std::string &particle_p4) {
  auto calculate_mt = [](ROOT::Math::PtEtaPhiMVector &particle_p4) {
    return particle_p4.Mt();;
  };
  return df.Define(outputname, calculate_mt, {particle_p4});
}


#endif /* GUARD_LEPSELECTION_H */