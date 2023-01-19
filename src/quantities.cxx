#ifndef GUARD_QUANTITIES_H
#define GUARD_QUANTITIES_H

#include "../include/basefunctions.hxx"
#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "../include/vectoroperations.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>
#include <Math/VectorUtil.h>

/// The namespace that is used to hold the functions for basic quantities that
/// are needed for every event
namespace quantities {

/// Function to calculate the pt from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode pt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &inputvector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.pt(); },
        {inputvector});
}
/// Function to calculate the eta from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the eta value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode eta(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &inputvector) {
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4) { return (float)p4.eta(); },
        {inputvector});
}
/// Function to calculate the eta from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the eta value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode phi(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &inputvector) {
    return df.Define(outputname,
                     [](const ROOT::Math::PtEtaPhiMVector &p4) {
                         if (p4.pt() <
                             0.0) // negative pt is used to mark invalid LVs
                             return default_float;
                         return (float)p4.phi();
                     },
                     {inputvector});
}
/// Function to calculate the mass from a given lorentz vector and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the mass value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode mass(ROOT::RDF::RNode df, const std::string &outputname,
                      const std::string &inputvector) {
    return df.Define(outputname,
                     [](const ROOT::Math::PtEtaPhiMVector &p4) {
                         if (p4.pt() <
                             0.0) // negative pt is used to mark invalid LVs
                             return default_float;
                         return (float)p4.mass();
                     },
                     {inputvector});
}
/// Function to writeout the dxy impact parameter from a particle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the dxy value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param dxycolumn name of the column containing the dxy values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode dxy(ROOT::RDF::RNode df, const std::string &outputname,
                     const int &position, const std::string &pairname,
                     const std::string &dxycolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<float> &dxy) {
            const int index = pair.at(position);
            return dxy.at(index, default_float);
        },
        {pairname, dxycolumn});
}
/// Function to writeout the dz impact parameter from a particle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the dz value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param dzcolumn name of the column containing the dz values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode dz(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &dzcolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<float> &dz) {
            const int index = pair.at(position);
            return dz.at(index, default_float);
        },
        {pairname, dzcolumn});
}
/// Function to writeout the charge of a particle. The particle is identified
/// via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the charge value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param chargecolumn name of the column containing the charge values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode charge(ROOT::RDF::RNode df, const std::string &outputname,
                        const int &position, const std::string &pairname,
                        const std::string &chargecolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<int> &charge) {
            const int index = pair.at(position);
            return charge.at(index, default_int);
        },
        {pairname, chargecolumn});
}
/// Function to calculate the scalar sum of pts for given lorentz vectors and add it to the
/// dataframe
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvector name of the column containing the lorentz vector
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode Lt(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &pt_1, const std::string &pt_2, const std::string &pt_3) {
    // build scalar sum of pts of 3 objects
    return df.Define(
        outputname,
        [](const float &pt_1,
           const float &pt_2, const float &pt_3) {
            if (pt_3 < 0.0 || pt_3 < 0.0 || pt_3 < 0.0)
                return default_float;
            auto const triple_lepton_pt = pt_1 + pt_2 + pt_3;
            return (float)triple_lepton_pt;
        },
        {pt_1, pt_2, pt_3});
}
/**
 * @brief function used to calculate the deltaPhi between two lorentz vectors. $\phi_1$ is from the first lorentz vector and $\phi_2$ is from the second lorentz vector.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the deltaR value
 * @param p_1_p4 first lorentz vector
 * @param p_2_p4 second lorentz vector of
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode deltaPhi(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4) {
    auto calculate_deltaPhi = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4) {
        return ROOT::Math::VectorUtil::DeltaPhi(p_1_p4, p_2_p4);
    };
    return df.Define(outputname, calculate_deltaPhi, {p_1_p4, p_2_p4});
}
/**
 * @brief function used to calculate the deltaPhi between the lepton from a W and the visible Higgs decay products. $\phi_1$ is from the first lorentz vector and $\phi_2$ is from the second lorentz vector and \phi_3$ is from the third lorentz vector.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the deltaR value
 * @param p_1_p4 first lorentz vector
 * @param p_2_p4 second lorentz vector
 * @param p_3_p4 second lorentz vector
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode deltaPhi_WH(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4, const std::string &p_3_p4) {
    auto calculate_deltaPhi = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4, ROOT::Math::PtEtaPhiMVector &p_3_p4) {
        auto const dileptonsystem = p_2_p4 + p_3_p4;
        return ROOT::Math::VectorUtil::DeltaPhi(p_1_p4, dileptonsystem);
    };
    return df.Define(outputname, calculate_deltaPhi, {p_1_p4, p_2_p4, p_3_p4});
}
/// Function to calculate the visible mass from a pair of lorentz vectors and
/// add it to the dataframe. The visible mass is calculated as the mass of the
/// lorentz vector of the dilepton system.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode m_vis(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::vector<std::string> &inputvectors) {
    // build visible mass from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
                return default_float;
            auto const dileptonsystem = p4_1 + p4_2;
            return (float)dileptonsystem.mass();
        },
        inputvectors);
}
/// Function to calculate the visible pt from a pair of lorentz vectors and
/// add it to the dataframe. The visible pt is calculated as the pt of the
/// lorentz vector of the dilepton system.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode pt_vis(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors) {
    // build visible pt from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
                return default_float;
            auto const dileptonsystem = p4_1 + p4_2;
            return (float)dileptonsystem.pt();
        },
        inputvectors);
}
/// Function to calculate the Higgs boson mass from a pair of lorentz vectors (visible Higgs system) and the neutrino 
/// four vector and
/// add it to the dataframe.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode m_tt_lt(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors) {
    // build visible pt from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2,
           const ROOT::Math::PtEtaPhiEVector &p4_3) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0 || p4_3.pt() < 0.0)
                return default_float;
            auto const higgssystem = p4_1 + p4_2 + p4_3;
            return (float)higgssystem.M();
        },
        inputvectors);
}
/// Function to calculate the pt of the composed system of all visible particles.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode pt_123(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors) {
    // build visible pt from the three particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2,
           const ROOT::Math::PtEtaPhiMVector &p4_3) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0 || p4_3.pt() < 0.0)
                return default_float;
            auto const pt_123 = p4_1 + p4_2 + p4_3;
            return (float)pt_123.Pt();
        },
        inputvectors);
}
/// Function to calculate the pt of the W from a the visible lepton fourvector, the met four vector and the neutrino four vector from the Higgs system and
/// add it to the dataframe.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode pt_W_lt(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors) {
    // build visible pt from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2,
           const ROOT::Math::PtEtaPhiEVector &p4_3) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0 || p4_3.pt() < 0.0)
                return default_float;
            auto const w_p4 = p4_1 + p4_2 - p4_3;
            return (float)w_p4.Pt();
        },
        inputvectors);
}
/// Function to calculate the visible eta from a pair of lorentz vectors and
/// add it to the dataframe. The visible eta is calculated as eta of the
/// lorentz vector of the dilepton system.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode eta_vis(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors) {
    // build visible pt from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
                return default_float;
            auto const dileptonsystem = p4_1 + p4_2;
            return (float)dileptonsystem.eta();
        },
        inputvectors);
}
/// Function to calculate the visible phi from a pair of lorentz vectors and
/// add it to the dataframe. The visible phi is calculated as phi of the
/// lorentz vector of the dilepton system.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value
/// \param inputvectors a vector of the two names of the columns containing the
/// required lorentz vectors
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode phi_vis(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors) {
    // build visible pt from the two particles
    return df.Define(
        outputname,
        [](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
                return default_float;
            auto const dileptonsystem = p4_1 + p4_2;
            return (float)dileptonsystem.phi();
        },
        inputvectors);
}
/**
 * @brief Function to calculate the quantity `pZetaMissVis` from the two leptons
 in the event + the met vector. The variable is defined as
 * \f[
     D_\zeta = p_\zeta^\text{miss} - 0.85 p_\zeta^\text{vis}
    \qquad;
    p_\zeta^\text{miss} = \vec{p}_\text{T}^\text{miss} \cdot \hat{\zeta}
    \qquad;
    p_\zeta^\text{vis} = (\vec{p}_\text{T}^{p_1} + \vec{p}_\text{T}^{p_2}) \cdot
 \hat{\zeta} \f] where \f$\vec{p}_\text{T}^{p_{1,2}}\f$ corresponds to the
 transverse momentum vector of the first (second) lepton and \f$\hat{\zeta}\f$
 to the bisectional direction between the two leptons in the transverse plane.
 For more information check

 D. Jang, “Search for MSSM Higgs decaying to
 tau pairs in pp collision at √s=1.96 TeV at CDF”. PhD thesis, Rutgers
 University, 2006. FERMILAB-THESIS-2006-11.

 * @param df the input dataframe
 * @param p_1_p4 the lorentz vector of the first particle
 * @param p_2_p4 the lorentz vector of the second particle
 * @param met the lorentz vector of the met
 * @param outputname the name of the new column containing the pZetaMissVis
 value
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode pzetamissvis(ROOT::RDF::RNode df,
                              const std::string &outputname,
                              const std::string &p_1_p4,
                              const std::string &p_2_p4,
                              const std::string &met) {
    float alpha = 0.85;
    auto calculate_pzetamissvis = [alpha](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                                          ROOT::Math::PtEtaPhiMVector &p_2_p4,
                                          ROOT::Math::PtEtaPhiMVector &met) {
        auto met_3dvec = met.Vect();
        met_3dvec.SetZ(0.0);
        // calculate zeta for the delepton system
        auto p1_norm = p_1_p4.Vect().Unit();
        auto p2_norm = p_2_p4.Vect().Unit();
        p1_norm.SetZ(0.0);
        p2_norm.SetZ(0.0);
        p1_norm = p1_norm.Unit();
        p2_norm = p2_norm.Unit();
        auto zeta = (p1_norm + p2_norm).Unit();

        auto dileptonsystem = p_1_p4.Vect() + p_2_p4.Vect();
        dileptonsystem.SetZ(0);
        auto pzetaVis = dileptonsystem.Dot(zeta);
        return met_3dvec.Dot(zeta) - (alpha * pzetaVis);
    };
    return df.Define(outputname, calculate_pzetamissvis, {p_1_p4, p_2_p4, met});
}
/**
 * @brief function used to calculate mTdileptonMET, which is the transverse mass
 * of the di-lepton system. The transverse mass is calculated using the
 * vectoroperations::calculateMT function.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the mTdileptonMET value
 * @param p_1_p4 lorentz vector of the first particle
 * @param p_2_p4 lorentz vector of the second particle
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode mTdileptonMET(ROOT::RDF::RNode df,
                               const std::string &outputname,
                               const std::string &p_1_p4,
                               const std::string &p_2_p4,
                               const std::string &met) {
    auto calculate_mTdileptonMET = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                                      ROOT::Math::PtEtaPhiMVector &p_2_p4,
                                      ROOT::Math::PtEtaPhiMVector &met) {
        ROOT::Math::PtEtaPhiMVector dilepton = p_1_p4 + p_2_p4;
        return vectoroperations::calculateMT(dilepton, met);
    };
    return df.Define(outputname, calculate_mTdileptonMET,
                     {p_1_p4, p_2_p4, met});
}

/**
 * @brief function used to calculate the deltaR between two lorentz vectors. It
 is defined as \f[ \Delta R = \sqrt{(\eta_1 - \eta_2)^2 + (\phi_1 - \phi_2)^2}
 \f$ where $\eta_1$ and $\phi_1$ are from the first lorentz vector and $\eta_2$
 and $\phi_2$ are from the second lorentz vector.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the deltaR value
 * @param p_1_p4 first lorentz vector
 * @param p_2_p4 second lorentz vector of
 * @return a new dataframe with the new column
 */
ROOT::RDF::RNode deltaR(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4) {
    auto calculate_deltaR = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4) {
        return ROOT::Math::VectorUtil::DeltaR(p_1_p4, p_2_p4);
    };
    return df.Define(outputname, calculate_deltaR, {p_1_p4, p_2_p4});
}
/**
 * @brief function used to calculate the transverse mass of a particle. The
 * transverse mass is calculated using the vectoroperations::calculateMT
 * function.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the mT value
 * @param particle_p4 lorentz vector of the particle
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode mT(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &particle_p4, const std::string &met) {
    auto calculate_mt = [](ROOT::Math::PtEtaPhiMVector &particle_p4,
                           ROOT::Math::PtEtaPhiMVector &met) {
        return vectoroperations::calculateMT(particle_p4, met);
    };
    return df.Define(outputname, calculate_mt, {particle_p4, met});
}

/**
 * @brief function used to calculate the pt of the dilepton + met system.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the pt_tt value
 * @param p_1_p4 lorentz vector of the first particle
 * @param p_2_p4 lorentz vector of the second particle
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode pt_tt(ROOT::RDF::RNode df, const std::string &outputname,
                       const std::string &p_1_p4, const std::string &p_2_p4,
                       const std::string &met) {
    auto calculate_pt_tt = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                              ROOT::Math::PtEtaPhiMVector &p_2_p4,
                              ROOT::Math::PtEtaPhiMVector &met) {
        auto dileptonmet = p_1_p4 + p_2_p4 + met;
        return (float)dileptonmet.Pt();
    };
    return df.Define(outputname, calculate_pt_tt, {p_1_p4, p_2_p4, met});
}

/**
 * @brief function used to calculate the pt of the dilepton + two leading jets +
 * met system. If the number of jets is less than 2, the quantity is set to 10
 * instead.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the pt_ttjj value
 * @param p_1_p4 lorentz vector of the first particle
 * @param p_2_p4 lorentz vector of the second particle
 * @param jet_1_p4 lorentz vector of the first jet
 * @param jet_2_p4 lorentz vector of the second jet
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode pt_ttjj(ROOT::RDF::RNode df, const std::string &outputname,
                         const std::string &p_1_p4, const std::string &p_2_p4,
                         const std::string &jet_1_p4,
                         const std::string &jet_2_p4, const std::string &met) {
    auto calculate_pt_ttjj = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                                ROOT::Math::PtEtaPhiMVector &p_2_p4,
                                ROOT::Math::PtEtaPhiMVector &jet_1_p4,
                                ROOT::Math::PtEtaPhiMVector &jet_2_p4,
                                ROOT::Math::PtEtaPhiMVector &met) {
        if (jet_1_p4.pt() < 0.0 || jet_2_p4.pt() < 0.0)
            return default_float;
        auto jetlepmet = p_1_p4 + p_2_p4 + met + jet_1_p4 + jet_2_p4;
        return (float)jetlepmet.Pt();
    };
    return df.Define(outputname, calculate_pt_ttjj,
                     {p_1_p4, p_2_p4, jet_1_p4, jet_2_p4, met});
}

/**
 * @brief function used to calculate `mt_tot`. It is defined as
 \f[
    m_{T}^{tot} = \sqrt{m_{T}^2(p_{1},E_{T}^{miss}) +
 m_{T}^2(p_{2},E_{T}^{miss}) + m_{T}^2(p_{1},p_2) } \f] where \f$ m_{T}^2 \f$ is
 the transverse mass, \f$ p_{1} \f$ and \f$ p_{2} \f$ are the lepton lorentz
 vectors and \f$ E_{T}^{miss} \f$ is the missing energy.
 *
 * @param df name of the dataframe
 * @param outputname name of the new column containing the mt_tot value
 * @param p_1_p4 lorentz vector of the first particle
 * @param p_2_p4 lorentz vector of the second particle
 * @param met lorentz vector of the met
 * @return a new dataframe with the new column
 */

ROOT::RDF::RNode mt_tot(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &p_1_p4, const std::string &p_2_p4,
                        const std::string &met) {
    auto calculate_mt_tot = [](ROOT::Math::PtEtaPhiMVector &p_1_p4,
                               ROOT::Math::PtEtaPhiMVector &p_2_p4,
                               ROOT::Math::PtEtaPhiMVector &met) {
        const float mt_1 = vectoroperations::calculateMT(p_1_p4, met);
        const float mt_2 = vectoroperations::calculateMT(p_2_p4, met);
        const float mt_mix = vectoroperations::calculateMT(p_1_p4, p_2_p4);
        return (float)sqrt(mt_1 * mt_1 + mt_2 * mt_2 + mt_mix * mt_mix);
    };
    return df.Define(outputname, calculate_mt_tot, {p_1_p4, p_2_p4, met});
}

/// Function to writeout the isolation of a particle. The particle is
/// identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the isolation value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param isolationcolumn name of the column containing the isolation
/// values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode isolation(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &isolationcolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<float> &isolation) {
                         const int index = pair.at(position);
                         return isolation.at(index, default_float);
                     },
                     {pairname, isolationcolumn});
}
/// Function to writeout the PDGID from a genparticle. The particle
/// is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the PDGID
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param pdgidcolumn name of the column containing the pdgID
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode pdgid(ROOT::RDF::RNode df, const std::string &outputname,
                       const int &position, const std::string &pairname,
                       const std::string &pdgidcolumn) {
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<int> &pdgid) {
            const int index = pair.at(position);
            return pdgid.at(index, default_pdgid);
        },
        {pairname, pdgidcolumn});
}
/// Function to determine number of good leptons
///
/// \param[in] df the input dataframe
/// \param[out] outputname the name of the produced quantity
/// \param[in] goodleptons name of the vector that contains a lepton mask of
/// good leptons, its length of non-zero values constitutes the output quantity
///
/// \return a dataframe containing the number of good leptons in an event
ROOT::RDF::RNode NumberOfGoodLeptons(ROOT::RDF::RNode df,
                                     const std::string &outputname,
                                     const std::string &goodleptons) {
    return df.Define(outputname,
                     [](const ROOT::RVec<int> &goodleptons) {
                         return (int)ROOT::VecOps::Nonzero(goodleptons).size();
                     },
                     {goodleptons});
}
/// namespace for tau specific quantities
namespace tau {
/// Function to writeout the decaymode of a tau. The particle is identified
/// via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the decaymode value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param decaymodecolumn name of the column containing the decaymode
/// values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode decaymode(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &decaymodecolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &decaymode) {
                         const int index = pair.at(position);
                         return decaymode.at(index, default_int);
                     },
                     {pairname, decaymodecolumn});
}
/// Function to writeout the genmatch of a tau. The particle is identified
/// via the index stored in the pair vector Genmatch values are defined as
/// \code
///   1 = prompt electron,
///   2 = prompt muon,
///   3 = tau->e decay,
///   4 = tau->mu decay,
///   5 = hadronic tau decay,
///   0 = unknown or unmatched
///   \endcode
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the genmatch value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param genmatchcolumn name of the column containing the genmatch values
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode genmatch(ROOT::RDF::RNode df, const std::string &outputname,
                          const int &position, const std::string &pairname,
                          const std::string &genmatchcolumn) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<UChar_t> &genmatch) {
                         const int index = pair.at(position);
                         return genmatch.at(index, default_uchar);
                     },
                     {pairname, genmatchcolumn});
}
/// Function to writeout the pt of the reco jet associated with the given
/// tau.
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the jet pt value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param taujet_index name of the column containing the association
/// between the tau and the reco jet \param jetpt_column name of the column
/// containing the recojet pt values
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode matching_jet_pt(ROOT::RDF::RNode df,
                                 const std::string &outputname,
                                 const int &position,
                                 const std::string &pairname,
                                 const std::string &taujet_index,
                                 const std::string &jetpt_column) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &taujets,
                                const ROOT::RVec<float> &jetpt) {
                         const int tauindex = pair.at(position);
                         const int jetindex = taujets.at(tauindex, -1);
                         return jetpt.at(jetindex, default_float);
                     },
                     {pairname, taujet_index, jetpt_column});
}
/// Function to writeout the pt of the gen jet associated with the reco jet,
/// which is associated with the given tau. \code
///  Tau --> recoJet --> GenJet
///   \endcode
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the jet pt value
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param taujet_index name of the column containing the association
/// between the tau and the reco jet \param genjet_index name of the column
/// containing the association between the reco jet and the gen jet \param
/// genjetpt_column name of the column containing the genJet pt values
///
/// \returns a dataframe with the new column
ROOT::RDF::RNode matching_genjet_pt(
    ROOT::RDF::RNode df, const std::string &outputname, const int &position,
    const std::string &pairname, const std::string &taujet_index,
    const std::string &genjet_index, const std::string &genjetpt_column) {
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<int> &taujets,
                                const ROOT::RVec<int> &genjets,
                                const ROOT::RVec<float> &genjetpt) {
                         const int tauindex = pair.at(position);
                         const int jetindex = taujets.at(tauindex, -1);
                         const int genjetindex = genjets.at(jetindex, -1);
                         return genjetpt.at(genjetindex, default_float);
                     },
                     {pairname, taujet_index, genjet_index, genjetpt_column});
}
/// Function to writeout a flag if a tau passes a specific tau id cut. The
/// particle is identified via the index stored in the pair vector
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the flag
/// \param position index of the position in the pair vector
/// \param pairname name of the column containing the pair vector
/// \param nameID name of the ID column in the NanoAOD
/// \param idxID bitvalue of the WP the has to be passed
///
/// \returns a dataframe with the new column

ROOT::RDF::RNode TauIDFlag(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &nameID, const int &idxID) {
    return df.Define(
        outputname,
        [position, idxID](const ROOT::RVec<int> &pair,
                          const ROOT::RVec<UChar_t> &IDs) {
            Logger::get("tauIDFlag")
                ->debug(
                    "position tau in pair {}, pair {}, id bit {}, vsjet ids {}",
                    position, pair, idxID, IDs);
            const int index = pair.at(position);
            const int ID = IDs.at(index, default_int);
            if (ID != default_int)
                return std::min(1, int(ID & 1 << (idxID - 1)));
            else
                return int(ID);
        },
        {pairname, nameID});
}
} // end namespace tau
/// namespace for muon specific quantities
namespace muon {
/**
 * @brief Function to writeout the id of a muon.
 *
 * @param df the dataframe to add the quantity to
 * @param outputname the name of the new quantity
 * @param position position of the muon in the pair vector
 * @param pairname name of the column containing the pair vector
 * @param idcolumn name of the column containing the muon id values
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &idcolumn) {
    Logger::get("muonIDflag")
                ->debug(
                    "muon ID  {}", position);
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<bool> &id) {
            const int index = pair.at(position, -1);
            return id.at(index, false);
        },
        {pairname, idcolumn});
}
/**
 * @brief Function to writeout the is global flag of a muon.
 *
 * @param df the dataframe to add the quantity to
 * @param outputname the name of the new quantity
 * @param position position of the muon in the pair vector
 * @param pairname name of the column containing the pair vector
 * @param globalflagcolumn name of the column containing the muon is global flag
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode is_global(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &globalflagcolumn) {
    Logger::get("muonIsGlobalflag")
                ->debug(
                    "is global pos {}", position);
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<bool> &globalflag) {
                         const int index = pair.at(position, -1);
                         return globalflag.at(index, false);
                     },
                     {pairname, globalflagcolumn});
}
/**
 * @brief Function to writeout the is tracker flag of a muon.
 *
 * @param df the dataframe to add the quantity to
 * @param outputname the name of the new quantity
 * @param position position of the muon in the pair vector
 * @param pairname name of the column containing the pair vector
 * @param trackerflagcolumn name of the column containing the muon is global flag
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode is_tracker(ROOT::RDF::RNode df, const std::string &outputname,
                           const int &position, const std::string &pairname,
                           const std::string &trackerflagcolumn) {
    Logger::get("muonIsTrackerflag")
                ->debug(
                    "is tracker pos {}", position);
    return df.Define(outputname,
                     [position](const ROOT::RVec<int> &pair,
                                const ROOT::RVec<bool> &trackerflag) {
                         const int index = pair.at(position, -1);
                         return trackerflag.at(index, false);
                     },
                     {pairname, trackerflagcolumn});
}
} // end namespace muon
namespace electron {
/**
 * @brief Function to writeout the id of a electron.
 *
 * @param df the dataframe to add the quantity to
 * @param outputname the name of the new quantity
 * @param position position of the muon in the pair vector
 * @param pairname name of the column containing the pair vector
 * @param idcolumn name of the column containing the muon id values
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode id(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &position, const std::string &pairname,
                    const std::string &idcolumn) {
    Logger::get("electronIDflag")
                ->debug(
                    "ele ID position {}", position);
    return df.Define(
        outputname,
        [position](const ROOT::RVec<int> &pair, const ROOT::RVec<bool> &id) {
            const int index = pair.at(position, -1);
            return id.at(index, false);
        },
        {pairname, idcolumn});
}
} // end namespace electron
} // end namespace quantities
#endif /* GUARD_QUANTITIES_H */