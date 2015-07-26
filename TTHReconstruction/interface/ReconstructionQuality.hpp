#ifndef RECONSTRUCTIONQUALITY
#define RECONSTRUCTIONQUALITY

#include "TFile.h"
#include "TLorentzVector.h"
#include "Interpretation.hpp"
#include "MECalculator.hpp"
#include "TH1F.h"
#include "TMath.h"
#include <iostream>
#include <algorithm> 

class ReconstructionQuality{
public:
  ReconstructionQuality(std::string filename="data/likelihoodhistos.root");

  float TTHChi2(float mthad, float mtlep, float mhiggs);
  float TTChi2(float mthad, float mtlep);
  float TTWHChi2(float mthad, float mtlep, float mhiggs, float mwhad);
  float TTWChi2(float mthad, float mtlep, float mwhad);
  float TTWHLikelihood(float mthad, float mtlep, float mhiggs, float mwhad);
  float TTWLikelihood(float mthad, float mtlep, float mwhad);
  float TTWHishLikelihood(float mthad, float mtlep, float mhiggs, float mwhad);
  float TTWishLikelihood(float mthad, float mtlep, float mwhad);

  float TTHChi2_tagged(float mthad, float mtlep, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2);
  float TTChi2_tagged(float mthad, float mtlep, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99.);
  float TTWHChi2_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2);
  float TTWChi2_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99.);
  float TTWHLikelihood_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2);
  float TTWLikelihood_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99.);
  float TTWHishLikelihood_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2);
  float TTWishLikelihood_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99.);

  float TTChi2_tagged_higgspt(float mthad, float mtlep, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99.,float pt =1.);
  float TTWChi2_tagged_higgspt(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99.,float pt =1.);
  float TTChi2_tagged_higgsjetpt(float mthad, float mtlep, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99., float p1=1, float pt2=1);
  float TTWChi2_tagged_higgsjetpt(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99., float p1=1, float pt2=1);


  float GetTag(std::string tag, Interpretation& i);

  float TTHChi2(Interpretation& i);
  float TTChi2(Interpretation& i);
  float TTWHChi2(Interpretation& i);
  float TTWChi2(Interpretation& i);

  float TTHChi2_tagged(Interpretation& i);
  float TTChi2_tagged(Interpretation& i, bool inclHiggsTags=true);
  float TTWHChi2_tagged(Interpretation& i);
  float TTWChi2_tagged(Interpretation& i, bool inclHiggsTags=true);

  float TTWHLikelihood(Interpretation& i);
  float TTWLikelihood(Interpretation& i);
  float TTWHLikelihood_tagged(Interpretation& i);
  float TTWLikelihood_tagged(Interpretation& i, bool inclHiggsTags=true);

  float TTWHishLikelihood(Interpretation& i);
  float TTWishLikelihood(Interpretation& i);
  float TTWHishLikelihood_tagged(Interpretation& i);
  float TTWishLikelihood_tagged(Interpretation& i, bool inclHiggsTags=true);

  float TTChi2_tagged_higgspt(Interpretation& i);
  float TTWChi2_tagged_higgspt(Interpretation& i);
  float TTChi2_tagged_higgsjetpt(Interpretation& i);
  float TTWChi2_tagged_higgsjetpt(Interpretation& i);

  float BLikelihood(float csv);
  float LLikelihood(float csv);
  float NBLikelihood(uint ntagged, uint njets, float* csvs);
  float TopHadLikelihood(float m);
  float TopHadishLikelihood(float m);
  float TopLepLikelihood(float m);
  float TopLepishLikelihood(float m);
  float WHadLikelihood(float m);
  float WHadishLikelihood(float m);
  float HiggsLikelihood(float m);
  float HiggsishLikelihood(float m);

  float TTH_ME(Interpretation& i);
  float TTBB_OFF_ME(Interpretation& i);
  float TTBB_ON_ME(Interpretation& i);
  float TTHBB_ME(Interpretation& i);
  float TTH_TTBB_ON_ME_RATIO(Interpretation& i);
  float TTHBB_TTBB_ON_ME_RATIO(Interpretation& i);
  float TTH_TTBB_OFF_ME_RATIO(Interpretation& i);
  float TTHBB_TTBB_OFF_ME_RATIO(Interpretation& i);

  float Interpolate(TH1F* histo, float value);

 
private:
  TFile* file;
  TH1F* h_CSV_b;
  TH1F* h_CSV_l_w_c;
  TH1F* h_M_Higgs_reco;
  TH1F* h_M_TopHad_reco;
  TH1F* h_M_TopLep_reco;
  TH1F* h_M_WHad_reco;
  TH1F* h_M_Higgs_best;
  TH1F* h_M_TopHad_best;
  TH1F* h_M_TopLep_best;
  TH1F* h_M_WHad_best;

  float higgs_mean;
  float tophad_mean;
  float toplep_mean;
  float whad_mean;
  float higgs_sigma;
  float tophad_sigma;
  float toplep_sigma;
  float whad_sigma;
  
  float btagcut;

  float btagbonus;
  float tiny_likelihood;

  MECalculator me;
};

#endif
