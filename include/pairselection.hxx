#ifndef GUARD_PAIRSELECTION_H
#define GUARD_PAIRSELECTION_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include "bitset"
#include "utility/Logger.hxx"
#include "utility/utility.hxx"

typedef std::bitset<15> StatusBits;
struct GenParticle {
    int index;
    int status;
    std::bitset<15> statusflag;
    int pdgid;
    int motherid;
};
bool check_mother(ROOT::RVec<GenParticle> genparticles, const int index,
                  const int mother_pdgid);

namespace ditau_pairselection {
ROOT::RDF::RNode buildgenpair(ROOT::RDF::RNode df, const std::string &recopair,
                              const std::string &genindex_particle1,
                              const std::string &genindex_particle2,
                              const std::string &genpair);
ROOT::RDF::RNode
buildtruegenpair(ROOT::RDF::RNode df, const std::string &statusflags,
                 const std::string &status, const std::string &pdgids,
                 const std::string &motherids, const std::string &pts,
                 const std::string &genpair, const int mother_pdgid,
                 const int daughter_1_pdgid, const int daughter_2_pdgid);
ROOT::RDF::RNode flagGoodPairs(ROOT::RDF::RNode df, const std::string &flagname,
                               const std::string &pairname);
auto compareForPairs(const ROOT::RVec<float> &lep1pt,
                     const ROOT::RVec<float> &lep1iso,
                     const ROOT::RVec<float> &lep2pt,
                     const ROOT::RVec<float> &lep2iso);
namespace semileptonic {
auto PairSelectionAlgo(const float &mindeltaR);
} // end namespace semileptonic

namespace fullhadronic {
auto PairSelectionAlgo(const float &mindeltaR);
} // end namespace fullhadronic

// namespace for full leptonic pairselection
namespace leptonic {
auto ElMuPairSelectionAlgo(const float &mindeltaR);
auto PairSelectionAlgo(const float &mindeltaR);
auto ZBosonPairSelectionAlgo(const float &mindeltaR);
} // namespace leptonic
namespace mutau {

ROOT::RDF::RNode PairSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &pairname,
                               const float &mindeltaR);
} // end namespace mutau

namespace eltau {

ROOT::RDF::RNode PairSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &pairname,
                               const float &mindeltaR);
} // end namespace eltau

namespace tautau {

ROOT::RDF::RNode PairSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &pairname,
                               const float &mindeltaR);

} // namespace tautau

namespace elmu {

ROOT::RDF::RNode PairSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &pairname,
                               const float &mindeltaR);
} // namespace elmu

namespace mumu {

ROOT::RDF::RNode PairSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &pairname,
                               const float &mindeltaR);

ROOT::RDF::RNode
ZBosonPairSelection(ROOT::RDF::RNode df,
                    const std::vector<std::string> &input_vector,
                    const std::string &pairname, const float &mindeltaR);
} // end namespace mumu
namespace elel {

ROOT::RDF::RNode PairSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &pairname,
                               const float &mindeltaR);

ROOT::RDF::RNode
ZBosonPairSelection(ROOT::RDF::RNode df,
                    const std::vector<std::string> &input_vector,
                    const std::string &pairname, const float &mindeltaR);
} // end namespace elel
} // namespace ditau_pairselection
#endif /* GUARD_PAIRSELECTION_H */