#ifndef RECONSTRUCTIONQUALITY
#define RECONSTRUCTIONQUALITY

#include "TFile.h"
#include "TLorentzVector.h"
#include "Interpretation.hpp"
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

  float TTHChi2_tagged(float mthad, float mtlep, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2);
  float TTChi2_tagged(float mthad, float mtlep, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99.);
  float TTWHChi2_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2);
  float TTWChi2_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1=-99., float tagb2=-99.);

  float GetTag(std::string tag, Interpretation& i);

  float TTHChi2(Interpretation& i);
  float TTChi2(Interpretation& i);
  float TTWHChi2(Interpretation& i);
  float TTWChi2(Interpretation& i);

  float TTHChi2_tagged(Interpretation& i);
  float TTChi2_tagged(Interpretation& i, bool inclHiggsTags=true);
  float TTWHChi2_tagged(Interpretation& i);
  float TTWChi2_tagged(Interpretation& i, bool inclHiggsTags=true);

  float BLikelihood(float csv);
  float LLikelihood(float csv);
  float NBLikelihood(uint ntagged, uint njets, float* csvs);
  float TopHadLikelihood(float m);
  float TopHadishLikelihood(float m);
  float TopLepLikelihood(float m);
  float TopLepishLikelihood(float m);
  float WHadLikelihood(float m);
  float WHadishLikelihood(float m);

 
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

};

#endif
