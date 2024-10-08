#ifndef GUARDLVECS_H
#define GUARDLVECS_H

#include "ROOT/RDFHelpers.hxx"
#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "utility/utility.hxx"
#include "defaults.hxx"
#include <Math/Vector4D.h>

namespace lorentzvectors {
ROOT::RDF::RNode buildparticle(ROOT::RDF::RNode df,
                               const std::vector<std::string> &quantities,
                               const std::string &outputname,
                               const int &position);
ROOT::RDF::RNode BuildP4Collection(ROOT::RDF::RNode df,
                                   const std::string &outputname,
                                   const std::string &pts,
                                   const std::string &etas,
                                   const std::string &phis,
                                   const std::string &masses,
                                   const std::string &collection_index);
ROOT::RDF::RNode BuildP4Collection(ROOT::RDF::RNode df,
                                          const std::string &outputname,
                                          const std::string &pts,
                                          const std::string &etas,
                                          const std::string &phis,
                                          const std::string &masses);
ROOT::RDF::RNode build(ROOT::RDF::RNode df,
                       const std::vector<std::string> &obj_quantities,
                       const int pairindex, const std::string &obj_p4_name);
ROOT::RDF::RNode buildMet(ROOT::RDF::RNode df, const std::string &met_pt,
                          const std::string &met_phi,
                          const std::string &outputname);

/// This function constructs a vectorial sum of n lorentz vectors. Since it uses
/// an aribtrary number of input vectors, it is implemented in the header file.
///
/// \param df The input dataframe
/// \param outputflag The name of the new column
/// \param p4s Parameter pack of column names that contain the considered
/// lorentz vectors (must be of type ROOT::Math::PtEtaPhiMVector)
///
/// \returns a dataframe containing the new column
template <class... Lorentzvectors>
ROOT::RDF::RNode CombineP4s(ROOT::RDF::RNode df, const std::string &outputflag,
                            const Lorentzvectors &...p4s) {
    std::vector<std::string> P4List;
    utility::appendParameterPackToVector(P4List, p4s...);
    const auto np4s = sizeof...(Lorentzvectors);
    using namespace ROOT::VecOps;
    return df.Define(
        outputflag,
        utility::PassAsVec<np4s, ROOT::Math::PtEtaPhiMVector>(
            [](const ROOT::RVec<ROOT::Math::PtEtaPhiMVector> &p4s) {
                for (int i = 0; i < p4s.size(); i++) {
                    auto p4 = p4s.at(i);
                    if (p4.pt() < 0.) {
                        return default_lorentzvector;
                    }
                }
                return Sum(p4s, ROOT::Math::PtEtaPhiMVector());
            }),
        P4List);
}

ROOT::RDF::RNode scaleP4(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::vector<std::string> &inputvector, const float &p4_miss_sf);   
/// namespace used for mutau lorentzvectors
namespace mutau {
ROOT::RDF::RNode build(ROOT::RDF::RNode df, const std::string &pairname,
                       const std::vector<std::string> &muon_quantities,
                       const std::vector<std::string> &tau_quantities,
                       const std::string &muon_p4_name,
                       const std::string &tau_p4_name);
} // namespace mutau

} // namespace lorentzvectors
#endif /* GUARDLVECS_H */