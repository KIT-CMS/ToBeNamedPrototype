#ifndef GUARD_TRIPLESELECTION_H
#define GUARD_TRIPLESELECTION_H

#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TVector2.h"
#include "bitset"
#include "utility/Logger.hxx"
#include "utility/utility.hxx"

typedef std::bitset<15> StatusBits;
bool check_mother(ROOT::RVec<GenParticle> genparticles, const int index,
                  const int mother_pdgid);


namespace whtautau_tripleselection {
auto compareForTriples(const ROOT::RVec<float> &lep1pt,
                     const ROOT::RVec<float> &lep1iso,
                     const ROOT::RVec<float> &lep2pt,
                     const ROOT::RVec<float> &lep2iso);
ROOT::RDF::RNode
buildtruegentriple(ROOT::RDF::RNode df, const std::string &statusflags,
                 const std::string &status, const std::string &pdgids,
                 const std::string &motherids, const std::string &pts,
                 const std::string &gentriple, const int mother_pdgid_1, const int mother_pdgid_23,
                 const int daughter_1_pdgid, const int daughter_2_pdgid, const int daughter_3_pdgid);
ROOT::RDF::RNode flagGoodTriples(ROOT::RDF::RNode df, const std::string &flagname,
                               const std::string &triplename);
ROOT::RDF::RNode buildgentriple(ROOT::RDF::RNode df, const std::string &recotriple,
                              const std::string &genindex_particle1,
                              const std::string &genindex_particle2,
                              const std::string &genindex_particle3,
                              const std::string &gentriple);
namespace three_flavor {
auto TripleSelectionAlgo(const float &mindeltaR_leptau, const float &mindeltaR_leplep, const std::string &triple);
} // end namespace three_flavor
namespace two_flavor {
auto TripleSelectionAlgo(const float &mindeltaR_leptau, const float &mindeltaR_leplep);
} // end namespace two_flavor
namespace lep_tautau {
auto TripleSelectionAlgo(const float &mindeltaR_leptau, const float &mindeltaR_tautau);
} //end namespace lep_tautau
namespace lep1lep1_lep2 {
auto TripleSelectionAlgo(const float &mindeltaR_lep1lep1, const float &mindeltaR_lep1lep2); //end namespace lep1lep1_lep2
}
namespace elemutau {
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_leplep);
ROOT::RDF::RNode TripleSelectionWOEle(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau);                               
} //end namespace elemutau
namespace mueletau {
    ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_leplep);
}//end namespace mueletau
namespace mumutau {
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_leplep);
}
namespace mu_tautau {
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_tautau);
}
namespace ele_tautau {
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_leptau, const float &mindeltaR_tautau);
}
namespace mumuele {
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_lep1lep1, const float &mindeltaR_lep1lep2);
}
namespace eleelemu {
ROOT::RDF::RNode TripleSelection(ROOT::RDF::RNode df,
                               const std::vector<std::string> &input_vector,
                               const std::string &triplename,
                               const float &mindeltaR_lep1lep1, const float &mindeltaR_lep1lep2);
}
}
#endif /* GUARD_TRIPLESELECTION_H */