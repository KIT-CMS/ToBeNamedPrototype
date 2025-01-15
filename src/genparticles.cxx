#ifndef GUARD_GENPARTICLES_H
#define GUARD_GENPARTICLES_H

#include "../include/utility/Logger.hxx"
#include "../include/defaults.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "bitset"
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>
#include <cmath>
typedef std::bitset<20> IntBits;

enum class GenMatchingCode : int {
    NONE = -1,
    IS_ELE_PROMPT = 1,
    IS_MUON_PROMPT = 2,
    IS_ELE_FROM_TAU = 3,
    IS_MUON_FROM_TAU = 4,
    IS_TAU_HAD_DECAY = 5,
    IS_FAKE = 6,
    IS_ELE_FROM_W = 7,
    IS_MUON_FROM_W = 8
};

enum class GenStatusFlag : int {
    NONE = -1,
    isPrompt = 0,
    isDecayedLeptonHadron = 1,
    isTauDecayProduct = 2,
    isPromptTauDecayProduct = 3,
    isDirectTauDecayProduct = 4,
    isDirectPromptTauDecayProduct = 5,
    isDirectHadronDecayProduct = 6,
    isHardProcess = 7,
    fromHardProcess = 8,
    isHardProcessTauDecayProduct = 9,
    isDirectHardProcessTauDecayProduct = 10,
    fromHardProcessBeforeFSR = 11,
    isFirstCopy = 12,
    isLastCopy = 13,
    isLastCopyBeforeFSR = 14
};

namespace genleptons {
ROOT::RDF::RNode GenLepton(ROOT::RDF::RNode df,
                           const std::string &genparticles_pt,
                           const std::string &genparticles_eta,
                           const std::string &genparticles_phi,
                           const std::string &genparticles_mass,
                           const std::string &genparticles_pdgid,
                           const std::string &genparticles_status,
                           const std::string &genparticles_statusFlag) {

    auto lambda_idx_1 = [](
        const ROOT::RVec<int> &rvec_pdgId,
        const ROOT::RVec<int> &rvec_status,
        const ROOT::RVec<UShort_t> &rvec_status_flag) {
        int idx = -99;
        for (unsigned int i = 0; i < rvec_pdgId.size(); i++) {
            int pdgid = rvec_pdgId.at(i);
            int status = rvec_status.at(i);
            int status_flag = rvec_status_flag.at(i);
            if ((pdgid == 11 || pdgid == 13) && status == 1 &&
                 IntBits(status_flag).test((int)GenStatusFlag::fromHardProcess) &&
                 IntBits(status_flag).test((int)GenStatusFlag::isPrompt) &&
                 !IntBits(status_flag).test((int)GenStatusFlag::isTauDecayProduct)) {
                if (idx > -1) {
                    Logger::get("GenLepton")->critical("ERROR: more than 1 final state lepton with pdgId {}, prev idx {}, new idx {}", pdgid, idx, i);
                    return -99;
                }
                idx = i;
            }
        }
        return idx;
    };

    auto lambda_idx_2 = [](
        const ROOT::RVec<int> &rvec_pdgId,
        const ROOT::RVec<int> &rvec_status,
        const ROOT::RVec<UShort_t> &rvec_status_flag) {
        int idx = -99;
        for (unsigned int i = 0; i < rvec_pdgId.size(); i++) {
            int pdgid = rvec_pdgId.at(i);
            int status = rvec_status.at(i);
            int status_flag = rvec_status_flag.at(i);
            if ((pdgid == -11 || pdgid == -13) && status == 1 &&
                 IntBits(status_flag).test((int)GenStatusFlag::fromHardProcess) &&
                 IntBits(status_flag).test((int)GenStatusFlag::isPrompt) &&
                 !IntBits(status_flag).test((int)GenStatusFlag::isTauDecayProduct)) {
                if (idx > -1) {
                    Logger::get("GenLepton")->critical("ERROR: more than 1 final state lepton with pdgId {}, prev idx {}, new idx {}", pdgid, idx, i);
                    return -99;
                }
                idx = i;
            }
        }
        return idx;
    };

    auto lambda_float = [](const int &idx, const ROOT::RVec<float> &col) {
        return col.at(idx, default_float);
    };
    auto lambda_int = [](const int &idx, const ROOT::RVec<int> &col) {
        return col.at(idx, default_int);
    };

    auto lambd_p4 = [](const float &pt, const float &eta, const float &phi, const float &mass) {
        ROOT::Math::PtEtaPhiMVector p4;
        if (pt > 0.) {
            p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
        } else {
            p4 = ROOT::Math::PtEtaPhiMVector(
                default_float, default_float,
                default_float, default_float
            );
        }
        return p4;
    };

    auto df1  =   df.Define("genlep_idx_1", lambda_idx_1, {genparticles_pdgid, genparticles_status, genparticles_statusFlag});
    auto df2  =  df1.Define("genlep_idx_2", lambda_idx_2, {genparticles_pdgid, genparticles_status, genparticles_statusFlag});
    auto df3  =  df2.Define("genlep_pt_1", lambda_float, {"genlep_idx_1", genparticles_pt});
    auto df4  =  df3.Define("genlep_eta_1", lambda_float, {"genlep_idx_1", genparticles_eta});
    auto df5  =  df4.Define("genlep_phi_1", lambda_float, {"genlep_idx_1", genparticles_phi});
    auto df6  =  df5.Define("genlep_mass_1", lambda_float, {"genlep_idx_1", genparticles_mass});
    auto df7  =  df6.Define("genlep_pdgId_1", lambda_int, {"genlep_idx_1", genparticles_pdgid});
    auto df8  =  df7.Define("genlep_pt_2", lambda_float, {"genlep_idx_2", genparticles_pt});
    auto df9  =  df8.Define("genlep_eta_2", lambda_float, {"genlep_idx_2", genparticles_eta});
    auto df10 =  df9.Define("genlep_phi_2", lambda_float, {"genlep_idx_2", genparticles_phi});
    auto df11 = df10.Define("genlep_mass_2", lambda_float, {"genlep_idx_2", genparticles_mass});
    auto df12 = df11.Define("genlep_pdgId_2", lambda_int, {"genlep_idx_2", genparticles_pdgid});
    auto df13 = df12.Define("genlep_p4_1", lambd_p4, {"genlep_pt_1", "genlep_eta_1", "genlep_phi_1", "genlep_mass_1"});
    auto df14 = df13.Define("genlep_p4_2", lambd_p4, {"genlep_pt_2", "genlep_eta_2", "genlep_phi_2", "genlep_mass_2"});

    return df14;
}

ROOT::RDF::RNode GenLeptonPreFSR(ROOT::RDF::RNode df,
                           const std::string &genparticles_pt,
                           const std::string &genparticles_eta,
                           const std::string &genparticles_phi,
                           const std::string &genparticles_mass,
                           const std::string &genparticles_pdgid,
                           const std::string &genparticles_status,
                           const std::string &genparticles_statusFlag) {

    auto lambda_idx_1 = [](
        const ROOT::RVec<int> &rvec_pdgId,
        const ROOT::RVec<int> &rvec_status,
        const ROOT::RVec<UShort_t> &rvec_status_flag) {
        int idx = -99;
        for (unsigned int i = 0; i < rvec_pdgId.size(); i++) {
            int pdgid = rvec_pdgId.at(i);
            int status = rvec_status.at(i);
            int status_flag = rvec_status_flag.at(i);
            if ((pdgid == 11 || pdgid == 13) &&
                 IntBits(status_flag).test((int)GenStatusFlag::isHardProcess) &&
                 IntBits(status_flag).test((int)GenStatusFlag::isPrompt) &&
                 !IntBits(status_flag).test((int)GenStatusFlag::isTauDecayProduct)) {
                if (idx > -1) {
                    Logger::get("GenLepton")->critical("ERROR: more than 1 final state lepton with pdgId {}, prev idx {}, new idx {}", pdgid, idx, i);
                    return -99;
                }
                idx = i;
            }
        }
        return idx;
    };

    auto lambda_idx_2 = [](
        const ROOT::RVec<int> &rvec_pdgId,
        const ROOT::RVec<int> &rvec_status,
        const ROOT::RVec<UShort_t> &rvec_status_flag) {
        int idx = -99;
        for (unsigned int i = 0; i < rvec_pdgId.size(); i++) {
            int pdgid = rvec_pdgId.at(i);
            int status = rvec_status.at(i);
            int status_flag = rvec_status_flag.at(i);
            if ((pdgid == -11 || pdgid == -13) &&
                 IntBits(status_flag).test((int)GenStatusFlag::isHardProcess) &&
                 IntBits(status_flag).test((int)GenStatusFlag::isPrompt) &&
                 !IntBits(status_flag).test((int)GenStatusFlag::isTauDecayProduct)) {
                if (idx > -1) {
                    Logger::get("GenLepton")->critical("ERROR: more than 1 final state lepton with pdgId {}, prev idx {}, new idx {}", pdgid, idx, i);
                    return -99;
                }
                idx = i;
            }
        }
        return idx;
    };

    auto lambda_float = [](const int &idx, const ROOT::RVec<float> &col) {
        return col.at(idx, default_float);
    };
    auto lambda_int = [](const int &idx, const ROOT::RVec<int> &col) {
        return col.at(idx, default_int);
    };

    auto lambd_p4 = [](const float &pt, const float &eta, const float &phi, const float &mass) {
        ROOT::Math::PtEtaPhiMVector p4;
        if (pt > 0.) {
            p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
        } else {
            p4 = ROOT::Math::PtEtaPhiMVector(
                default_float, default_float,
                default_float, default_float
            );
        }
        return p4;
    };

    auto df1  =   df.Define("genlepPreFSR_idx_1", lambda_idx_1, {genparticles_pdgid, genparticles_status, genparticles_statusFlag});
    auto df2  =  df1.Define("genlepPreFSR_idx_2", lambda_idx_2, {genparticles_pdgid, genparticles_status, genparticles_statusFlag});
    auto df3  =  df2.Define("genlepPreFSR_pt_1", lambda_float, {"genlepPreFSR_idx_1", genparticles_pt});
    auto df4  =  df3.Define("genlepPreFSR_eta_1", lambda_float, {"genlepPreFSR_idx_1", genparticles_eta});
    auto df5  =  df4.Define("genlepPreFSR_phi_1", lambda_float, {"genlepPreFSR_idx_1", genparticles_phi});
    auto df6  =  df5.Define("genlepPreFSR_mass_1", lambda_float, {"genlepPreFSR_idx_1", genparticles_mass});
    auto df7  =  df6.Define("genlepPreFSR_pdgId_1", lambda_int, {"genlepPreFSR_idx_1", genparticles_pdgid});
    auto df8  =  df7.Define("genlepPreFSR_pt_2", lambda_float, {"genlepPreFSR_idx_2", genparticles_pt});
    auto df9  =  df8.Define("genlepPreFSR_eta_2", lambda_float, {"genlepPreFSR_idx_2", genparticles_eta});
    auto df10 =  df9.Define("genlepPreFSR_phi_2", lambda_float, {"genlepPreFSR_idx_2", genparticles_phi});
    auto df11 = df10.Define("genlepPreFSR_mass_2", lambda_float, {"genlepPreFSR_idx_2", genparticles_mass});
    auto df12 = df11.Define("genlepPreFSR_pdgId_2", lambda_int, {"genlepPreFSR_idx_2", genparticles_pdgid});
    auto df13 = df12.Define("genlepPreFSR_p4_1", lambd_p4, {"genlepPreFSR_pt_1", "genlepPreFSR_eta_1", "genlepPreFSR_phi_1", "genlepPreFSR_mass_1"});
    auto df14 = df13.Define("genlepPreFSR_p4_2", lambd_p4, {"genlepPreFSR_pt_2", "genlepPreFSR_eta_2", "genlepPreFSR_phi_2", "genlepPreFSR_mass_2"});

    return df14;
}

ROOT::RDF::RNode GenDressedLepton(ROOT::RDF::RNode df,
                           const std::string &genparticles_pt,
                           const std::string &genparticles_eta,
                           const std::string &genparticles_phi,
                           const std::string &genparticles_mass,
                           const std::string &genparticles_pdgid,
                           const std::string &genparticles_hasTauAnc,
                           const std::string &genlep_p4_1,
                           const std::string &genlep_p4_2) {

    auto lambda_idx_1 = [](
        const ROOT::RVec<float> &rvec_pt,
        const ROOT::RVec<float> &rvec_eta,
        const ROOT::RVec<float> &rvec_phi,
        const ROOT::RVec<float> &rvec_mass,
        const ROOT::RVec<int> &rvec_pdgId,
        const ROOT::RVec<bool> &rvec_hasTauAnc,
        ROOT::Math::PtEtaPhiMVector &genlep_p4) {
        int idx = -99;
        float dR_min = 0.3;

        if (genlep_p4.pt() < 0.) {
            return idx;
        }

        for (unsigned int i = 0; i < rvec_pt.size(); i++) {
            ROOT::Math::PtEtaPhiMVector p4 = ROOT::Math::PtEtaPhiMVector(rvec_pt.at(i), rvec_eta.at(i), rvec_phi.at(i), rvec_mass.at(i));
            int pdgid = rvec_pdgId.at(i);
            bool hasTauAnc = rvec_hasTauAnc.at(i);
            float dR_tmp = (float)ROOT::Math::VectorUtil::DeltaR(p4, genlep_p4);
            if ((pdgid == 11 || pdgid == 13) && !hasTauAnc && (dR_tmp < dR_min)) {
                idx = i;
                dR_min = dR_tmp;
            }
        }
        return idx;
    };
    auto lambda_idx_2 = [](
        const ROOT::RVec<float> &rvec_pt,
        const ROOT::RVec<float> &rvec_eta,
        const ROOT::RVec<float> &rvec_phi,
        const ROOT::RVec<float> &rvec_mass,
        const ROOT::RVec<int> &rvec_pdgId,
        const ROOT::RVec<bool> &rvec_hasTauAnc,
        ROOT::Math::PtEtaPhiMVector &genlep_p4) {
        int idx = -99;
        float dR_min = 0.3;

        if (genlep_p4.pt() < 0.) {
            return idx;
        }

        for (unsigned int i = 0; i < rvec_pt.size(); i++) {
            ROOT::Math::PtEtaPhiMVector p4 = ROOT::Math::PtEtaPhiMVector(rvec_pt.at(i), rvec_eta.at(i), rvec_phi.at(i), rvec_mass.at(i));
            int pdgid = rvec_pdgId.at(i);
            bool hasTauAnc = rvec_hasTauAnc.at(i);
            float dR_tmp = (float)ROOT::Math::VectorUtil::DeltaR(p4, genlep_p4);
            if ((pdgid == -11 || pdgid == -13) && !hasTauAnc && (dR_tmp < dR_min)) {
                idx = i;
                dR_min = dR_tmp;
            }
        }
        return idx;
    };

    auto lambda_float = [](const int &idx, const ROOT::RVec<float> &col) {
        return col.at(idx, default_float);
    };
    auto lambda_int = [](const int &idx, const ROOT::RVec<int> &col) {
        return col.at(idx, default_int);
    };

    auto lambd_p4 = [](const float &pt, const float &eta, const float &phi, const float &mass) {
        ROOT::Math::PtEtaPhiMVector p4;
        if (pt > 0.) {
            p4 = ROOT::Math::PtEtaPhiMVector(pt, eta, phi, mass);
        } else {
            p4 = ROOT::Math::PtEtaPhiMVector(
                default_float, default_float,
                default_float, default_float
            );
        }
        return p4;
    };

    auto lambd_dR = [](ROOT::Math::PtEtaPhiMVector &genlep_p4, ROOT::Math::PtEtaPhiMVector &genDressed_p4) {
        if (genlep_p4.pt() < 0. || genDressed_p4.pt() < 0.) {
            return (float)99.;
        }

        return (float)ROOT::Math::VectorUtil::DeltaR(genDressed_p4, genlep_p4);
    };

    auto df1  =   df.Define("genDressed_idx_1", lambda_idx_1, {genparticles_pt, genparticles_eta, genparticles_phi, genparticles_mass, genparticles_pdgid, genparticles_hasTauAnc, genlep_p4_1});
    auto df2  =  df1.Define("genDressed_idx_2", lambda_idx_2, {genparticles_pt, genparticles_eta, genparticles_phi, genparticles_mass, genparticles_pdgid, genparticles_hasTauAnc, genlep_p4_2});
    auto df3  =  df2.Define("genDressed_pt_1", lambda_float, {"genDressed_idx_1", genparticles_pt});
    auto df4  =  df3.Define("genDressed_eta_1", lambda_float, {"genDressed_idx_1", genparticles_eta});
    auto df5  =  df4.Define("genDressed_phi_1", lambda_float, {"genDressed_idx_1", genparticles_phi});
    auto df6  =  df5.Define("genDressed_mass_1", lambda_float, {"genDressed_idx_1", genparticles_mass});
    auto df7  =  df6.Define("genDressed_pdgId_1", lambda_int, {"genDressed_idx_1", genparticles_pdgid});
    auto df8  =  df7.Define("genDressed_pt_2", lambda_float, {"genDressed_idx_2", genparticles_pt});
    auto df9  =  df8.Define("genDressed_eta_2", lambda_float, {"genDressed_idx_2", genparticles_eta});
    auto df10 =  df9.Define("genDressed_phi_2", lambda_float, {"genDressed_idx_2", genparticles_phi});
    auto df11 = df10.Define("genDressed_mass_2", lambda_float, {"genDressed_idx_2", genparticles_mass});
    auto df12 = df11.Define("genDressed_pdgId_2", lambda_int, {"genDressed_idx_2", genparticles_pdgid});
    auto df13 = df12.Define("genDressed_p4_1", lambd_p4, {"genDressed_pt_1", "genDressed_eta_1", "genDressed_phi_1", "genDressed_mass_1"});
    auto df14 = df13.Define("genDressed_p4_2", lambd_p4, {"genDressed_pt_2", "genDressed_eta_2", "genDressed_phi_2", "genDressed_mass_2"});
    auto df15 = df14.Define("genDressed_dR_1", lambd_dR, {genlep_p4_1, "genDressed_p4_1"});
    auto df16 = df15.Define("genDressed_dR_2", lambd_dR, {genlep_p4_2, "genDressed_p4_2"});

    return df16;
}

}

namespace genflag {
/**
 * @brief Function to writeout a boolean flag to select a specific DY decay mode based on gen-level PDG ID
 *
 * @param df The input dataframe
 * @param outputname The name of the output column
 * @param genparticles_pdgid The name of the column containing the pdgids of the
 genparticles
 * @param genparticles_statusFlag The name of the column containing the
 statusFlags of the genparticles
 * @param pdgId The PDG ID of decayed leptons
 * @return a dataframe with the output flag as a column named outputname
 */
ROOT::RDF::RNode DYGenFlag(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const int &pdgId) {
    auto lambda = [pdgId](const ROOT::RVec<int> &pdgids, const ROOT::RVec<UShort_t> &status_flags) {
        bool found_0 = false;
        bool found_1 = false;
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            int pdgid = pdgids.at(i);
            int status_flag = status_flags.at(i);
            if (pdgid == pdgId && IntBits(status_flag).test((int)GenStatusFlag::isHardProcess)) {
                found_0 = true;
            } else if (pdgid == -pdgId && IntBits(status_flag).test((int)GenStatusFlag::isHardProcess)) {
                found_1 = true;
            }
            if (found_0 && found_1)
                break;
        }
        return (found_0 && found_1);
    };
    auto df1 =
        df.Define(outputname, lambda,
                  {genparticles_pdgid, genparticles_statusFlag});
    return df1;
}
ROOT::RDF::RNode WGenFlag(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const int &pdgId) {
    auto lambda = [pdgId](const ROOT::RVec<int> &pdgids, const ROOT::RVec<UShort_t> &status_flags) {
        bool found_0 = false;
        bool found_1 = false;
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            int pdgid = pdgids.at(i);
            int status_flag = status_flags.at(i);
            if (pdgid == pdgId && IntBits(status_flag).test((int)GenStatusFlag::isHardProcess)) {
                found_0 = true;
            } else if (pdgid == -pdgId && IntBits(status_flag).test((int)GenStatusFlag::isHardProcess)) {
                found_1 = true;
            }
            if (found_0 || found_1)
                break;
        }
        return (found_0 || found_1);
    };
    auto df1 =
        df.Define(outputname, lambda,
                  {genparticles_pdgid, genparticles_statusFlag});
    return df1;
}
}

namespace genmatching {
namespace tau {
/**
 * @brief implementation of the genmatching used in tau analyses, based
 * on
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * implementation. The genmatches are represented as integer flags:
    Decaytype         | Value
    ------------------|-------
    NONE              | -1,
    IS_ELE_PROMPT     | 1,
    IS_MUON_PROMPT    | 2,
    IS_ELE_FROM_TAU   | 3,
    IS_MUON_FROM_TAU  | 4,
    IS_TAU_HAD_DECAY  | 5,
    IS_FAKE           | 6

    The genmatch is caluclated individually for each taupair candidate.

 *
 * @param df The dataframe to be used for the genmatching.
 * @param outputname The name of the output column.
 * @param hadronicGenTaus The name of the column containing the hadronicGenTaus
 from genmatching::tau::hadronicGenTaus
 * @param genparticles_pdgid The name of the column containing the pdgids of the
 genparticles
 * @param genparticles_statusFlag The name of the column containing the
 statusFlags of the genparticles
 * @param genparticles_pt The name of the column containing the pT of the
 genparticles
 * @param genparticles_eta The name of the column containing the eta of the
 genparticles
 * @param genparticles_phi The name of the column containing the phi of the
 genparticles
 * @param genparticles_mass The name of the column containing the mass of the
 genparticles
 * @param lepton_p4 The name of the column containing the p4 of the lepton to be
 genmatched
 * @return a dataframe with the genmatching as a column named outputname
 */
ROOT::RDF::RNode genmatching(ROOT::RDF::RNode df, const std::string &outputname,
                             const std::string &hadronicGenTaus,
                             const std::string &genparticles_pdgid,
                             const std::string &genparticles_statusFlag,
                             const std::string &genparticles_pt,
                             const std::string &genparticles_eta,
                             const std::string &genparticles_phi,
                             const std::string &genparticles_mass,
                             const std::string &lepton_p4) {
    auto match_lepton = [](const std::vector<int> &hadronicGenTaus,
                           const ROOT::RVec<int> &pdgids,
                           const ROOT::RVec<UShort_t> &status_flags,
                           const ROOT::RVec<float> &pts,
                           const ROOT::RVec<float> &etas,
                           const ROOT::RVec<float> &phis,
                           const ROOT::RVec<float> &masses,
                           const ROOT::Math::PtEtaPhiMVector &lepton_p4) {
        // find closest lepton fulfilling the requirements
        float min_delta_r = 9999;
        int closest_genparticle_index = 0;
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            int pdgid = std::abs(pdgids.at(i));
            // check
            // 1. if there is a gen electron or muon close to the lepton
            // 2. that the genparticle pt is larger than 8 GeV
            // 3. the genparticle is isPrompt (statusbit 0) or
            // isDirectPromptTauDecayProduct (statusbit 5)
            bool statusbit = (IntBits(status_flags.at(i)).test(0) ||
                              IntBits(status_flags.at(i)).test(5));
            if ((pdgid == 11 || pdgid == 13) && pts.at(i) > 8 && statusbit) {
                ROOT::Math::PtEtaPhiMVector gen_p4(pts.at(i), etas.at(i),
                                                   phis.at(i), masses.at(i));
                float delta_r =
                    ROOT::Math::VectorUtil::DeltaR(gen_p4, lepton_p4);
                if (delta_r < min_delta_r) {
                    closest_genparticle_index = i;
                    min_delta_r = delta_r;
                }
            }
        }
        Logger::get("genmatching::tau::genmatching")
            ->debug("closest genlepton {} // DeltaR {}",
                    closest_genparticle_index, min_delta_r);
        // now loop trough the gentaus and check, if they are closer to the
        // lepton than the closest lepton genparticle
        for (auto hadronicGenTau : hadronicGenTaus) {
            // check if the hadronicGenTau is closer to the lepton than the
            // closest lepton genparticle
            ROOT::Math::PtEtaPhiMVector hadronicGenTau_p4(
                pts.at(hadronicGenTau), etas.at(hadronicGenTau),
                phis.at(hadronicGenTau), masses.at(hadronicGenTau));
            float gentau_delta_r =
                ROOT::Math::VectorUtil::DeltaR(hadronicGenTau_p4, lepton_p4);
            // the decay is considered a hadronic decay (statusbit 5) if
            // 1. the hadronicGenTau pt is larger than 15 GeV
            // 2. the delta_r is smaller than 0.2
            // 3. the delta_r is smaller than the closest lepton genparticle
            // delta_r
            if (hadronicGenTau_p4.Pt() > 15 && gentau_delta_r < 0.2 &&
                gentau_delta_r < min_delta_r) {
                // statusbit 5 is hadronic tau decay
                Logger::get("genmatching::tau::genmatching")
                    ->debug(
                        "found hadronicGenTau closer than closest lepton: {}",
                        gentau_delta_r);
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_TAU_HAD_DECAY");
                return (int)GenMatchingCode::IS_TAU_HAD_DECAY;
            }
        }
        // if it is not a hadronic decay, check if the lepton is close
        // enough (deltaR < 0.2)
        int closest_pdgid = std::abs(pdgids.at(closest_genparticle_index));
        if (min_delta_r < 0.2) {
            bool prompt =
                IntBits(status_flags.at(closest_genparticle_index)).test(0);
            bool from_tau =
                IntBits(status_flags.at(closest_genparticle_index)).test(5);
            if (closest_pdgid == 11 && prompt) {
                // statusbit 1 is prompt electron
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_ELE_PROMPT");
                return (int)GenMatchingCode::IS_ELE_PROMPT;
            }
            if (closest_pdgid == 13 && prompt) {
                // statusbit 2 is prompt muon
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_MUON_PROMPT");
                return (int)GenMatchingCode::IS_MUON_PROMPT;
            }
            if (closest_pdgid == 11 && from_tau) {
                // statusbit 3 is electron from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_ELE_FROM_TAU");
                return (int)GenMatchingCode::IS_ELE_FROM_TAU;
            }
            if (closest_pdgid == 13 && from_tau) {
                // statusbit 4 is muon from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_MUON_FROM_TAU");
                return (int)GenMatchingCode::IS_MUON_FROM_TAU;
            }
        }
        // if no genlepton was found within the deltaR < 0.2, return fake
        // (statusbit 6)
        Logger::get("genmatching::tau::genmatching")->debug("IS_FAKE");
        return (int)GenMatchingCode::IS_FAKE;
    };

    auto df1 =
        df.Define(outputname, match_lepton,
                  {hadronicGenTaus, genparticles_pdgid, genparticles_statusFlag,
                   genparticles_pt, genparticles_eta, genparticles_phi,
                   genparticles_mass, lepton_p4});
    return df1;
}
/**
 * @brief implementation of the genmatching used in tau analyses, based
 * on
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * implementation. The genmatches are represented as integer flags:
    Decaytype         | Value
    ------------------|-------
    NONE              | -1,
    IS_ELE_PROMPT     | 1,
    IS_MUON_PROMPT    | 2,
    IS_ELE_FROM_TAU   | 3,
    IS_MUON_FROM_TAU  | 4,
    IS_TAU_HAD_DECAY  | 5,
    IS_FAKE           | 6
    IS_ELE_FROM_W     | 7,
    IS_MUON_FROM_W    | 8

    The genmatch is caluclated individually for each taupair candidate.

 *
 * @param df The dataframe to be used for the genmatching.
 * @param outputname The name of the output column.
 * @param hadronicGenTaus The name of the column containing the hadronicGenTaus
 from genmatching::tau::hadronicGenTaus
 * @param genparticles_pdgid The name of the column containing the pdgids of the
 genparticles
 * @param genparticles_statusFlag The name of the column containing the
 statusFlags of the genparticles
 * @param genparticles_pt The name of the column containing the pT of the
 genparticles
 * @param genparticles_eta The name of the column containing the eta of the
 genparticles
 * @param genparticles_phi The name of the column containing the phi of the
 genparticles
 * @param genparticles_mass The name of the column containing the mass of the
 genparticles
 * @param genparticle_motheridx The name of the column containing the mother
 indices of the genparticles
 * @param genparticles_status The name of the column containing the status of
the genparticles
 * @param lepton_p4 The name of the column containing the p4 of the lepton to be
 genmatched
 * @return a dataframe with the genmatching as a column named outputname
 */
ROOT::RDF::RNode genmatching_wh(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &hadronicGenTaus, const std::string &genparticles_pdgid,
    const std::string &genparticles_statusFlag,
    const std::string &genparticles_pt, const std::string &genparticles_eta,
    const std::string &genparticles_phi, const std::string &genparticles_mass,
    const std::string &genparticle_motheridx,
    const std::string &genparticles_status, const std::string &lepton_p4) {
    auto match_lepton = [](const std::vector<int> &hadronicGenTaus,
                           const ROOT::RVec<int> &pdgids,
                           const ROOT::RVec<int> &mother_idx,
                           const ROOT::RVec<UShort_t> &status_flags,
                           const ROOT::RVec<int> &status,
                           const ROOT::RVec<float> &pts,
                           const ROOT::RVec<float> &etas,
                           const ROOT::RVec<float> &phis,
                           const ROOT::RVec<float> &masses,
                           const ROOT::Math::PtEtaPhiMVector &lepton_p4) {
        // find closest lepton fulfilling the requirements
        float min_delta_r = 9999;
        int closest_genparticle_index = 0;
        int closest_genparticle_mother_pdgid = 0;
        int closest_genparticle_mother_statusFlag = 0;
        int closest_genparticle_mother_status = 0;
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            int pdgid = std::abs(pdgids.at(i));
            // check
            // 1. if there is a gen electron or muon close to the lepton
            // 2. that the genparticle pt is larger than 8 GeV
            // 3. the genparticle is isPrompt (statusbit 0) or
            // isDirectPromptTauDecayProduct (statusbit 5)
            bool statusbit = (IntBits(status_flags.at(i)).test(0) ||
                              IntBits(status_flags.at(i)).test(5));
            if ((pdgid == 11 || pdgid == 13) && pts.at(i) > 8 && statusbit) {
                ROOT::Math::PtEtaPhiMVector gen_p4(pts.at(i), etas.at(i),
                                                   phis.at(i), masses.at(i));
                float delta_r =
                    ROOT::Math::VectorUtil::DeltaR(gen_p4, lepton_p4);
                if (delta_r < min_delta_r) {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("pdgids {}, status {}, status_flags {}, "
                                "mother_idx {}",
                                pdgids, status, status_flags, mother_idx);
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("mother_idx {}, pdgids {}, status {}, "
                                "status_flags {}",
                                mother_idx.at(i), pdgids.at(i),
                                status_flags.at(i), status.at(i));
                    if (mother_idx.at(i) == -1) {
                        closest_genparticle_index = i;
                        closest_genparticle_mother_pdgid = pdgids.at(i);
                        closest_genparticle_mother_status = status.at(i);
                        closest_genparticle_mother_statusFlag =
                            status_flags.at(i);
                        min_delta_r = delta_r;
                    } else {
                        closest_genparticle_index = i;
                        closest_genparticle_mother_pdgid =
                            pdgids.at(mother_idx.at(i));
                        closest_genparticle_mother_status =
                            status.at(mother_idx.at(i));
                        closest_genparticle_mother_statusFlag =
                            status_flags.at(mother_idx.at(i));
                        min_delta_r = delta_r;
                    }
                }
            }
        }
        Logger::get("genmatching::tau::genmatching")
            ->debug("closest genlepton {} // DeltaR {}",
                    closest_genparticle_index, min_delta_r);
        // now loop through the gentaus and check, if they are closer to the
        // lepton than the closest lepton genparticle
        for (auto hadronicGenTau : hadronicGenTaus) {
            // check if the hadronicGenTau is closer to the lepton than the
            // closest lepton genparticle
            ROOT::Math::PtEtaPhiMVector hadronicGenTau_p4(
                pts.at(hadronicGenTau), etas.at(hadronicGenTau),
                phis.at(hadronicGenTau), masses.at(hadronicGenTau));
            float gentau_delta_r =
                ROOT::Math::VectorUtil::DeltaR(hadronicGenTau_p4, lepton_p4);
            // the decay is considered a hadronic decay (statusbit 5) if
            // 1. the hadronicGenTau pt is larger than 15 GeV
            // 2. the delta_r is smaller than 0.2
            // 3. the delta_r is smaller than the closest lepton genparticle
            // delta_r
            if (hadronicGenTau_p4.Pt() > 15 && gentau_delta_r < 0.2 &&
                gentau_delta_r < min_delta_r) {
                // statusbit 5 is hadronic tau decay
                Logger::get("genmatching::tau::genmatching")
                    ->debug(
                        "found hadronicGenTau closer than closest lepton: {}",
                        gentau_delta_r);
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_TAU_HAD_DECAY");
                return (int)GenMatchingCode::IS_TAU_HAD_DECAY;
            }
        }
        // if it is not a hadronic decay, check if the lepton is close
        // enough (deltaR < 0.2)
        int closest_pdgid = std::abs(pdgids.at(closest_genparticle_index));
        if (min_delta_r < 0.2) {
            bool prompt =
                IntBits(status_flags.at(closest_genparticle_index)).test(0);
            bool from_tau =
                IntBits(status_flags.at(closest_genparticle_index)).test(5);
            if (closest_pdgid == 11 && prompt) {
                // statusbit 1 is prompt electron
                if (abs(closest_genparticle_mother_pdgid) == 24) {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("IS_ELE_FROM_W");
                    return (int)GenMatchingCode::IS_ELE_FROM_W;
                } else {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("IS_ELE_PROMPT");
                    return (int)GenMatchingCode::IS_ELE_PROMPT;
                }
            }
            if (closest_pdgid == 13 && prompt) {
                // statusbit 2 is prompt muon
                if (abs(closest_genparticle_mother_pdgid) == 24) {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("IS_MUON_FROM_W");
                    return (int)GenMatchingCode::IS_MUON_FROM_W;
                } else {
                    Logger::get("genmatching::tau::genmatching")
                        ->debug("IS_MUON_PROMPT");
                    return (int)GenMatchingCode::IS_MUON_PROMPT;
                }
            }
            if (closest_pdgid == 11 && from_tau) {
                // statusbit 3 is electron from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_ELE_FROM_TAU");
                return (int)GenMatchingCode::IS_ELE_FROM_TAU;
            }
            if (closest_pdgid == 13 && from_tau) {
                // statusbit 4 is muon from tau
                Logger::get("genmatching::tau::genmatching")
                    ->debug("IS_MUON_FROM_TAU");
                return (int)GenMatchingCode::IS_MUON_FROM_TAU;
            }
        }
        // if no genlepton was found within the deltaR < 0.2, return fake
        // (statusbit 6)
        Logger::get("genmatching::tau::genmatching")->debug("IS_FAKE");
        return (int)GenMatchingCode::IS_FAKE;
    };

    auto df1 = df.Define(
        outputname, match_lepton,
        {hadronicGenTaus, genparticles_pdgid, genparticle_motheridx,
         genparticles_statusFlag, genparticles_status, genparticles_pt,
         genparticles_eta, genparticles_phi, genparticles_mass, lepton_p4});
    return df1;
}
/**
 * @brief function to find all hadronicGenTaus needed for the genmatching, based
 * on
 * https://github.com/KIT-CMS/Artus/blob/dictchanges/KappaAnalysis/src/Utility/GeneratorInfo.cc
 * implementation. Loop trough all genparticles and check, if a genparticle,
 that is prompt, without any leptonic daughters can be found. If this
 genparticle has a neutrino as daughter, the hadronicGenTau is added to the list
 of hadronicGenTaus.
 *
 * @param df The dataframe to be extended with the hadronicGenTaus
 * @param outputname The name of the column to be added to the dataframe
 * @param genparticles_pdgid The name of the column containing the pdgids of the
 genparticles
 * @param genparticles_statusFlag The name of the column containing the
 statusflags of the genparticles
 * @param genparticles_motherid The name of the column containing the motherids
 of the genparticles
 * @return a new dataframe with the hadronicGenTaus added to the dataframe

 */

ROOT::RDF::RNode hadronicGenTaus(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const std::string &genparticles_pdgid,
                                 const std::string &genparticles_statusFlag,
                                 const std::string &genparticles_motherid) {

    auto gentaus = [](const ROOT::RVec<int> &pdgids,
                      const ROOT::RVec<UShort_t> &status_flags,
                      const ROOT::RVec<int> &mother_index) {
        // set default values for the output
        std::vector<int> hadronicGenTaus;
        if (pdgids.size() == 0) {
            hadronicGenTaus.push_back(-1);
            return hadronicGenTaus;
        }
        // debug printout
        for (int i = 0; i < pdgids.size(); i++) {
            Logger::get("genmatching::tau::genpair")
                ->debug("genparticle index {}, pdgid: {}, status_flags {}, "
                        "mother_index {}",
                        i, pdgids.at(i), status_flags.at(i),
                        mother_index.at(i));
        }
        // now loop though the genparticles and find the taus
        for (unsigned int i = 0; i < pdgids.size(); i++) {
            // check if the particle is a tau
            if (std::abs(pdgids.at(i)) == 15) {
                // check if the particle is a stable one
                bool prompt = IntBits(status_flags.at(i)).test(0);
                if (prompt) {
                    Logger::get("genmatching::tau::genpair")
                        ->debug("Found prompt tau: {}", i);
                    // find all daughters of the tau by checking which particles
                    // have the tauindex i as mother
                    std::vector<int> daughters;
                    for (unsigned int j = 0; j < mother_index.size(); j++) {
                        if (mother_index.at(j) == i) {
                            daughters.push_back(j);
                            Logger::get("genmatching::tau::genpair")
                                ->debug("daughters of {} : {}", i, j);
                        }
                    }
                    // check if the tau has at least one daughter
                    if (daughters.size() > 0) {
                        // check if the tau has at least one daughter that is a
                        // tau
                        bool hasTauDaughter = false;
                        bool hasLeptonDaughter = false;
                        for (unsigned int j = 0; j < daughters.size(); j++) {
                            int daughter_pdgid =
                                std::abs(pdgids.at(daughters.at(j)));
                            Logger::get("genmatching::tau::genpair")
                                ->debug("daughter {} : pdgid {}", j,
                                        daughter_pdgid);
                            if (daughter_pdgid == 15) {
                                hasTauDaughter = true;
                            }
                            if (daughter_pdgid == 11 || daughter_pdgid == 13) {
                                hasLeptonDaughter = true;
                            }
                        }
                        // in this case, it is not the correct tau, continue
                        if (hasTauDaughter || hasLeptonDaughter) {
                            continue;
                        }
                        // now run again through the daughters and check if
                        // there is a neutrino, if this is the case, we have the
                        // correct gentau
                        for (unsigned int j = 0; j < daughters.size(); j++) {
                            int daughter_pdgid =
                                std::abs(pdgids.at(daughters.at(j)));
                            if (daughter_pdgid == 12 || daughter_pdgid == 14 ||
                                daughter_pdgid == 16) {
                                Logger::get("genmatching::tau::genpair")
                                    ->debug("gentau found: {}", i);
                                hadronicGenTaus.push_back(i);
                            }
                        }
                    }
                }
            }
        }
        Logger::get("genmatching::tau::genpair")
            ->debug("found {} hadronic hadronicGenTaus",
                    hadronicGenTaus.size());
        for (int i = 0; i < hadronicGenTaus.size(); i++) {
            Logger::get("genmatching::tau::genpair")
                ->debug("hadronicGenTaus {} : {}", i, hadronicGenTaus.at(i));
        }
        return hadronicGenTaus;
    };

    auto df1 = df.Define(
        outputname, gentaus,
        {genparticles_pdgid, genparticles_statusFlag, genparticles_motherid});
    return df1;
}

} // end namespace tau

} // end namespace genmatching

#endif /* GUARD_GENPARTICLES_H */