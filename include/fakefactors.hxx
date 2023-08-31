#ifndef GUARDFAKEFACTORS_H
#define GUARDFAKEFACTORS_H

namespace fakefactors {

ROOT::RDF::RNode
raw_fakefactor_nmssm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                        const std::string &tau_pt, const std::string &njets,
                        const std::string &lep_mt, const std::string &nbtags,
                        const std::string &variation,
                        const std::string &ff_file);
ROOT::RDF::RNode
raw_fakefactor_nmssm_tt(ROOT::RDF::RNode df, const std::string &outputname,
                        const int &tau_idx, const std::string &tau_pt_1,
                        const std::string &tau_pt_2, const std::string &njets,
                        const std::string &variation,
                        const std::string &ff_file);
ROOT::RDF::RNode
fakefactor_nmssm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                    const std::string &tau_pt, const std::string &njets,
                    const std::string &lep_mt, const std::string &nbtags,
                    const std::string &lep_pt, const std::string &lep_iso,
                    const std::string &m_vis, const std::string &variation,
                    const std::string &ff_file,
                    const std::string &ff_corr_file);
ROOT::RDF::RNode
fakefactor_nmssm_tt(ROOT::RDF::RNode df, const std::string &outputname,
                    const int &tau_idx, const std::string &tau_pt_1,
                    const std::string &tau_pt_2, const std::string &njets,
                    const std::string &m_vis, const std::string &variation,
                    const std::string &ff_file,
                    const std::string &ff_corr_file);
  
ROOT::RDF::RNode
raw_fakefactor_sm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                     const std::string &qcd_outputname,
                     const std::string &wjets_outputname,
                     const std::string &ttbar_outputname,
                     const std::string &tau_pt, const std::string &njets,
                     const std::string &lep_mt, const std::string &delta_r,
                     const std::string &variation, const std::string &ff_file);
ROOT::RDF::RNode
fakefactor_sm_lt(ROOT::RDF::RNode df, const std::string &outputname,
                 const std::string &qcd_outputname,
                 const std::string &wjets_outputname,
                 const std::string &ttbar_outputname, const std::string &tau_pt,
                 const std::string &njets, const std::string &lep_mt,
                 const std::string &lep_pt, const std::string &lep_iso,
                 const std::string &m_vis, const std::string &delta_r,
                 const std::string &variation, const std::string &ff_file,
                 const std::string &ff_corr_file);
  
ROOT::RDF::RNode fakefactor_sm_lt_no_deltaR(
    ROOT::RDF::RNode df, const std::string &outputname,
    const std::string &qcd_outputname, const std::string &wjets_outputname,
    const std::string &ttbar_outputname, const std::string &tau_pt,
    const std::string &njets, const std::string &lep_mt,
    const std::string &lep_pt, const std::string &lep_iso,
    const std::string &m_vis, const std::string &variation,
    const std::string &ff_file, const std::string &ff_corr_file);
} // namespace fakefactors
#endif /* GUARDFAKEFACTORS_H */