#include "ROOT/RDataFrame.hxx"
#include "ROOT/RVec.hxx"
#include "TRandom3.h"

struct UnsignedInt {
    unsigned int mMax;
    UnsignedInt(unsigned int max) : mMax(max) {};
    unsigned int operator()() { return gRandom->Integer(mMax); }
};

struct GaussVector {
    float mMean;
    float mSigma;
    float mMin;
    GaussVector(float mean, float sigma, float min) : mMean(mean), mSigma(sigma), mMin(min) {};
    ROOT::RVec<float> operator()(unsigned int n) {
        ROOT::RVec<float> x(n);
        for (auto i=0u; i<n; i++)
            x[i] = std::max((float)gRandom->Gaus(mMean, mSigma), mMin);
        return x;
    }
};

struct BoolVector {
    float mThreshold;
    BoolVector(float threshold) : mThreshold(threshold) {};
    ROOT::RVec<bool> operator()(unsigned int n) {
        ROOT::RVec<bool> x(n);
        for (auto i=0u; i<n; i++)
            x[i] = gRandom->Rndm() > mThreshold;
        return x;
    }
};


template<class T>
struct ConstantVector {
    T mValue;
    ConstantVector(T value) : mValue(value) {};
    ROOT::RVec<T> operator()(unsigned int n) {
        ROOT::RVec<T> x(n, mValue);
        return x;
    }
};

auto constValue() { return (unsigned int) 42; }


int main() {
    gRandom->SetSeed(1234);
    const auto n_events = 100;
    ROOT::RDataFrame df(n_events);
    auto df2 = df.Define("run", constValue, {})
                 .Define("luminosityBlock", constValue, {})
                 .Define("event", constValue, {})
                 .Define("nMuon", UnsignedInt(5), {})
                 .Define("Muon_pt", GaussVector(30, 50, 0), {"nMuon"})
                 .Define("Muon_eta", GaussVector(0, 3, -999), {"nMuon"})
                 .Define("Muon_phi", GaussVector(0, 3, -999), {"nMuon"})
                 .Define("Muon_mass", GaussVector(0.1, 0.1, 0), {"nMuon"})
                 .Define("Muon_mediumId", BoolVector(0.5), {"nMuon"})
                 .Define("Muon_pfRelIso04_all", GaussVector(0.3, 0.1, 0), {"nMuon"})
                 .Define("nTau", UnsignedInt(5), {})
                 .Define("Tau_pt", GaussVector(30, 50, 0), {"nTau"})
                 .Define("Tau_eta", GaussVector(0, 3, -999), {"nTau"})
                 .Define("Tau_phi", GaussVector(0, 3, -999), {"nTau"})
                 .Define("Tau_mass", GaussVector(0.1, 0.1, 0), {"nTau"})
                 .Define("Tau_dz", GaussVector(0, 0.1, -999), {"nTau"})
                 .Define("Tau_decayMode", ConstantVector<int>(1), {"nTau"})
                 .Define("Tau_idDeepTau2017v2p1VSjet", ConstantVector<unsigned char>(1), {"nTau"})
                 .Define("Tau_idDeepTau2017v2p1VSe", ConstantVector<unsigned char>(1), {"nTau"})
                 .Define("Tau_idDeepTau2017v2p1VSmu", ConstantVector<unsigned char>(1), {"nTau"})
                 .Define("Tau_rawDeepTau2017v2p1VSjet", ConstantVector<float>(1), {"nTau"})
                 .Define("nVertices", UnsignedInt(10), {})
                 .Define("Flag_goodVertices", [](){ return true; }, {});
    const auto output_path = "sample.root";
    df2.Snapshot("Events", output_path);
}
