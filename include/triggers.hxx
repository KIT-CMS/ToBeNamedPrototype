#ifndef GUARD_TRIGGERS_H
#define GUARD_TRIGGERS_H

typedef std::bitset<20> IntBits;

namespace trigger {

bool matchParticle(const ROOT::Math::PtEtaPhiMVector &particle,
                   ROOT::RVec<float> &triggerobject_pts,
                   ROOT::RVec<float> &triggerobject_etas,
                   ROOT::RVec<float> &triggerobject_phis,
                   ROOT::RVec<int> &triggerobject_bits,
                   ROOT::RVec<int> &triggerobject_ids, const float &matchDeltaR,
                   const float &pt_cut, const float &eta_cut,
                   const int &trigger_particle_id_cut,
                   const int &triggerbit_cut);

ROOT::RDF::RNode GenerateSingleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &hltpath, const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold);

ROOT::RDF::RNode GenerateDoubleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle1_p4, const std::string &particle2_p4,
    const std::string &triggerobject_bits, const std::string &triggerobject_id,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &hltpath,
    const float &p1_pt_cut, const float &p2_pt_cut, const float &p1_eta_cut,
    const float &p2_eta_cut, const int &p1_trigger_particle_id_cut,
    const int &p2_trigger_particle_id_cut, const int &p1_triggerbit_cut,
    const int &p2_triggerbit_cut, const float &DeltaR_threshold);

ROOT::RDF::RNode MatchSingleTriggerObject(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold);

ROOT::RDF::RNode MatchDoubleTriggerObject(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle1_p4, const std::string &particle2_p4,
    const std::string &triggerobject_bits, const std::string &triggerobject_id,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const float &p1_pt_cut,
    const float &p2_pt_cut, const float &p1_eta_cut, const float &p2_eta_cut,
    const int &p1_trigger_particle_id_cut,
    const int &p2_trigger_particle_id_cut, const int &p1_triggerbit_cut,
    const int &p2_triggerbit_cut, const float &DeltaR_threshold);

namespace tagandprobe {

bool matchParticle(const ROOT::Math::PtEtaPhiMVector &particle,
                   ROOT::RVec<float> &triggerobject_pts,
                   ROOT::RVec<float> &triggerobject_etas,
                   ROOT::RVec<float> &triggerobject_phis,
                   ROOT::RVec<int> &triggerobject_bits,
                   ROOT::RVec<int> &triggerobject_ids, const float &matchDeltaR,
                   const float &pt_cut, const float &eta_cut,
                   const int &trigger_particle_id_cut,
                   const int &triggerbit_cut,
                   const float &trigger_particle_pt_cut);

ROOT::RDF::RNode GenerateSingleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle_p4, const std::string &triggerobject_bits,
    const std::string &triggerobject_id, const std::string &triggerobject_pt,
    const std::string &triggerobject_eta, const std::string &triggerobject_phi,
    const std::string &hltpath, const float &pt_cut, const float &eta_cut,
    const int &trigger_particle_id_cut, const int &triggerbit_cut,
    const float &DeltaR_threshold, const float &trigger_particle_pt_cut);

ROOT::RDF::RNode GenerateDoubleTriggerFlag(
    ROOT::RDF::RNode df, const std::string &triggerflag_name,
    const std::string &particle1_p4, const std::string &particle2_p4,
    const std::string &triggerobject_bits, const std::string &triggerobject_id,
    const std::string &triggerobject_pt, const std::string &triggerobject_eta,
    const std::string &triggerobject_phi, const std::string &hltpath,
    const float &p1_pt_cut, const float &p2_pt_cut, const float &p1_eta_cut,
    const float &p2_eta_cut, const int &p1_trigger_particle_id_cut,
    const int &p2_trigger_particle_id_cut, const int &p1_triggerbit_cut,
    const int &p2_triggerbit_cut, const float &DeltaR_threshold,
    const float &p1_trigger_particle_pt_cut,
    const float &p2_trigger_particle_pt_cut);

} // end namespace tagandprobe

} // end namespace trigger
#endif /* GUARD_TRIGGERS_H */