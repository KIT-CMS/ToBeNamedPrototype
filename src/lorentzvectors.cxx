#ifndef GUARDLVECS_H
#define GUARDLVECS_H

#include "../include/defaults.hxx"
#include "../include/utility/Logger.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include <Math/Vector4D.h>

/// Namespace used for lorentzvector operations

namespace lorentzvectors {

/// Function to build the lorentzvector from the pt, eta, phi and mass of a
/// particle. This utilizes the [PtEtaPhiMVector from
/// ROOT](https://root.cern/doc/master/namespaceROOT_1_1Math.html#a6cea5921731c7ac99dea921fb188df31)
///
/// \param df The input dataframe
/// \param quantities This vector contains the names of the columns:
///    -# pair - In this vector the particle index is stored. This index is used
///    to select the correct particle from the four particle quantity vectors.
///    -# pts - In this vector, the pts of the particle are stored.
///    -# etas - In this vector, the etas of the particle are stored.
///    -# phis - In this vector, the phis of the particle are stored.
///    -# masses - In this vector, the masses of the particle are stored.
///      This order has to be kept!
/// \param outputname The name of the output column in the new dataframe
/// \param position The position in the pair vector, which is used to store the
/// index of the particle in the particle quantity vectors.
///
/// \returns a new dataframe, which contains the new lorentz vector
ROOT::RDF::RNode buildparticle(ROOT::RDF::RNode df,
                               const std::vector<std::string> &quantities,
                               const std::string &outputname,
                               const int &position) {
    auto df1 = df.Define(
        outputname,
        [position, outputname](
            const ROOT::RVec<int> &pair, const ROOT::RVec<float> &pts,
            const ROOT::RVec<float> &etas, const ROOT::RVec<float> &phis,
            const ROOT::RVec<float> &masses) {
            // the index of the particle is stored in the pair vector
            ROOT::Math::PtEtaPhiMVector p4;
            Logger::get("lorentzvectors")
                ->debug("starting to build 4vector {}!", outputname);
            try {
                const int index = pair.at(position);
                Logger::get("lorentzvectors")->debug("pair {}", pair);
                Logger::get("lorentzvectors")->debug("pts {}", pts);
                Logger::get("lorentzvectors")->debug("etas {}", etas);
                Logger::get("lorentzvectors")->debug("phis {}", phis);
                Logger::get("lorentzvectors")->debug("masses {}", masses);
                Logger::get("lorentzvectors")->debug("Index {}", index);

                p4 = ROOT::Math::PtEtaPhiMVector(pts.at(index), etas.at(index),
                                                 phis.at(index),
                                                 masses.at(index));
            } catch (const std::out_of_range &e) {
                p4 = ROOT::Math::PtEtaPhiMVector(default_float, default_float,
                                                 default_float, default_float);
                Logger::get("lorentzvectors")
                    ->debug("Index not found, retuning dummy vector !");
            }
            Logger::get("lorentzvectors")
                ->debug("P4 - Particle {} : {}", position, p4);
            return p4;
        },
        quantities);
    return df1;
}

/**
 * @brief Function used to construct a 4-vector for a pair particle.
 *
 * @param df the input dataframe
 * @param obj_quantities a vector of strings containing the names of the
 * quantities used for the 4 vector building. As this function is a wrapper for
 * lorentzvectors::buildparticle, the order of the inputs has to match the order
 * used in that function
 * @param pairindex the index of the particle in the pair to be build
 * @param obj_p4_name name of the column containing the 4-vector
 * @return a new df, containing the thw column
 */

ROOT::RDF::RNode build(ROOT::RDF::RNode df,
                       const std::vector<std::string> &obj_quantities,
                       const int pairindex, const std::string &obj_p4_name) {
    Logger::get("lorentzvectors")->debug("Building {}", obj_p4_name);
    for (auto i : obj_quantities)
        Logger::get("lorentzvectors")->debug("Used object quantities {}", i);
    return lorentzvectors::buildparticle(df, obj_quantities, obj_p4_name,
                                         pairindex);
}

/**
 * @brief Function used to construct the missing transverse energy lorentz
 * vector.
 *
 * @param df the input dataframe
 * @param met_pt the met pt value
 * @param met_phi the met phi value
 * @param outputname name of the new column containing the missing transverse
 * energy lorentz vector.
 * @return a new df, containing the new column
 */
ROOT::RDF::RNode buildMet(ROOT::RDF::RNode df, const std::string &met_pt,
                          const std::string &met_phi,
                          const std::string &outputname) {
    auto construct_metvector = [](const float &pt, const float &phi) {
        // for Met, eta is zero
        auto met = ROOT::Math::PtEtaPhiEVector(pt, 0, phi, pt);
        // cast Met vector to a ROOT::Math::PtEtaPhiMVector to make latter
        // functions easier to use
        return (ROOT::Math::PtEtaPhiMVector)met;
    };
    return df.Define(outputname, construct_metvector, {met_pt, met_phi});
}
/**
 * @brief Function used to construct the missing transverse energy lorentz
 * vector from the dileptonsystem from the Higgs boson
 *
 * @param df the input dataframe
 * @param outputname name of the new column containing the missing transverse
 * energy lorentz vector.
 * @param inputvectors a vector of the two names of the columns containing the required lorentz vectors
 * @param p4_miss_sf the scale factor, that is applied to the di lepton system to obtain the vector of missing
 * energy
 * @return a new df, containing the new column
 */
/// Function to calculate the neutrino four vector from a pair of lorentz vectors (visible Higgs system) and
/// add it to the dataframe. see AN(HIG-19-010)
///
/// \param df the dataframe to add the quantity to
/// \param outputname name of the new column containing the pt value

///
/// \returns a dataframe with the new column
ROOT::RDF::RNode p4_miss(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvectors, const float &p4_miss_sf) {
    return df.Define(
        outputname,
        [p4_miss_sf](const ROOT::Math::PtEtaPhiMVector &p4_1,
           const ROOT::Math::PtEtaPhiMVector &p4_2) {
            if (p4_1.pt() < 0.0 || p4_2.pt() < 0.0)
                return ROOT::Math::PtEtaPhiEVector(default_float,default_float,default_float,
                default_float);;
            auto dileptonsystem = p4_1 + p4_2;
            auto dileptonsystem_scaled = p4_miss_sf*dileptonsystem;
            ROOT::Math::PtEtaPhiEVector neutrinos_p4 = ROOT::Math::PtEtaPhiEVector(dileptonsystem_scaled.Pt(),dileptonsystem_scaled.Eta(),dileptonsystem_scaled.Phi(),
                dileptonsystem_scaled.E());
            return neutrinos_p4;
        },
        inputvectors);
}
} // namespace lorentzvectors
#endif /* GUARDLVECS_H */