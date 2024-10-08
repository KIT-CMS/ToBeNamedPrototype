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
                    ->debug("Index not found, returning dummy vector !");
            }
            Logger::get("lorentzvectors")
                ->debug("P4 - Particle {} : {}", position, p4);
            return p4;
        },
        quantities);
    return df1;
}

/**
 * @brief Build a new column with a collection of Lorentz vectors per dataframe entry, which are
 * created from the pt, eta, phi and mass columns of a collection.
 *
 * For instance, this can be used to create four-vector objects from the four-vector component
 * columns of a NanoAOD collection.
 *
 * The function expects pt, eta phi and mass columns, which contain `ROOT::VecOps::RVec<float>`
 * objects. In addition, the argument `collection_index` must point to a column which contains
 * index lists of type `ROOT::VecOps::RVec<int>`. These index lists contain the indices of
 * elements, for which four-vectors are built. The output column contains collections of
 * four-vectors `ROOT::VecOps::RVec<PtEtaPhiMVector>`. The output collection only contains
 * the elements, which have been selected in the `collection_index` list.
 *
 * @param df Input dataframe
 * @param outputname name of the output column
 * @param pts column name with the pt values of the collection's objects
 * @param etas column name with the eta values of the collection's objects
 * @param phis column name with the phi values of the collection's objects
 * @param masses column name with the mass values of the collection's objects
 * @param collection_index column name for index lists of the elements for which four-vectors are built
 *
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode BuildP4Collection(ROOT::RDF::RNode df,
                                   const std::string &outputname,
                                   const std::string &pts,
                                   const std::string &etas,
                                   const std::string &phis,
                                   const std::string &masses,
                                   const std::string &collection_index) {
    auto build_collection_p4 = [](const ROOT::RVec<float> &pts,
                                  const ROOT::RVec<float> &etas,
                                  const ROOT::RVec<float> &phis,
                                  const ROOT::RVec<float> &masses,
                                  const ROOT::RVec<int> &collection_index) {
        auto pts_ = ROOT::VecOps::Take(pts, collection_index);
        auto etas_ = ROOT::VecOps::Take(etas, collection_index);
        auto phis_ = ROOT::VecOps::Take(phis, collection_index);
        auto masses_ = ROOT::VecOps::Take(masses, collection_index);
        return ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(pts_, etas_, phis_, masses_);
    };
    return df.Define(outputname, build_collection_p4, {pts, etas, phis, masses, collection_index});
}

/** 
 * @brief Build a new column with a collection of Lorentz vectors per dataframe entry, which are
 * created from the pt, eta, phi and mass columns of a collection.
 *
 * For instance, this can be used to create four-vector objects from the four-vector component
 * columns of a NanoAOD collection.
 *
 * The function expects pt, eta phi and mass columns, which contain `ROOT::VecOps::RVec<float>`
 * objects. The output column contains four-vectors `ROOT::VecOps::RVec<PtEtaPhiMVector>`.
 *
 * @param df Input dataframe
 * @param outputname name of the output column
 * @param pts column name with the pt values of the collection's objects
 * @param etas column name with the eta values of the collection's objects
 * @param phis column name with the phi values of the collection's objects
 * @param masses column name with the mass values of the collection's objects
 *
 * @return a dataframe with the new column
 */
ROOT::RDF::RNode BuildP4Collection(ROOT::RDF::RNode df,
                                          const std::string &outputname,
                                          const std::string &pts,
                                          const std::string &etas,
                                          const std::string &phis,
                                          const std::string &masses) {
    auto build_collection_p4 = [](const ROOT::RVec<float> &pts,
                                  const ROOT::RVec<float> &etas,
                                  const ROOT::RVec<float> &phis,
                                  const ROOT::RVec<float> &masses) {
        return ROOT::VecOps::Construct<ROOT::Math::PtEtaPhiMVector>(pts, etas, phis, masses);
    };
    return df.Define(outputname, build_collection_p4, {pts, etas, phis, masses});
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
 * @brief Function to scale a lorentz vector by a scale factor.
 *
 * @param df the input dataframe
 * @param outputname name of the new column containing the missing transverse
 * energy lorentz vector.
 * @param inputvector a vector of the two names of the columns containing the
 * required lorentz vectors
 * @param p4_sf the scale factor, that is applied to the input lorentz vector
 * @return a new df, containing the new column
 */
ROOT::RDF::RNode scaleP4(ROOT::RDF::RNode df, const std::string &outputname,
                         const std::vector<std::string> &inputvector,
                         const float &p4_sf) {
    return df.Define(
        outputname,
        [p4_sf](const ROOT::Math::PtEtaPhiMVector &p4) {
            if (p4.pt() < 0.0)
                return ROOT::Math::PtEtaPhiMVector(
                    default_float, default_float, default_float, default_float);
            auto p4_scaled = p4_sf * p4;
            return p4_scaled;
        },
        inputvector);
}
} // namespace lorentzvectors
#endif /* GUARDLVECS_H */