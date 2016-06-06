#ifndef RECONSTRUCTIONQUALITY2
#define RECONSTRUCTIONQUALITY2

#include "TFile.h"
#include "TLorentzVector.h"
#include "Interpretation.hpp"
#include "InterpretationGenerator.hpp"
#include "MECalculator.hpp"
#include "TH1F.h"
#include "TMath.h"
#include <iostream>
#include <algorithm> 

class ReconstructionQuality2{
public:
  ReconstructionQuality2(std::string filename="data/likelihoodhistos.root");

  float GetTag(std::string tag, Interpretation& i);
  // get chi2 containing 3 (4) terms: hadronic top-mass, leptonic top-mass, W-mass (, Higgs/BB-mass)
  // get likelihood for the masses of hadronic top, leptonic top, W (, Higgs/BB-mass) for a perfectly reconstructed tth-interpreatation
  float TTWHLikelihood(Interpretation& i);
  float TTWHishLikelihood(Interpretation& i);
  float TTWHLikelihoodTimesPSME(Interpretation& i);
  float TTWHishLikelihoodTimesPSME(Interpretation& i);
  float TTWHLikelihoodTimesSCME(Interpretation& i);
  float TTWHishLikelihoodTimesSCME(Interpretation& i);
  float TTWHLikelihoodTimesCSVTimesPSME(Interpretation& i);
  float TTWHLikelihoodTimesCSVTimesSCME(Interpretation& i);


  // get the b-tagger likelihoods
  float BLikelihood(float csv);
  float LLikelihood(float csv);
  float NBLikelihood(uint ntagged, uint njets, const float* csvs);

  // get likelihoods for different masses
  float TopHadLikelihood(float m);
  float TopHadishLikelihood(float m);
  float TopLepLikelihood(float m);
  float TopLepishLikelihood(float m);
  float WHadLikelihood(float m);
  float WHadishLikelihood(float m);
  float HiggsLikelihood(float m);
  float HiggsishLikelihood(float m);

  float TTX0_SC_ME(Interpretation& i);
  float TTX0_PS_ME(Interpretation& i);

  float TTX0_SC_ME(const TLorentzVector& t, const TLorentzVector& tbar, const TLorentzVector& H);
  float TTX0_PS_ME(const TLorentzVector& t, const TLorentzVector& tbar, const TLorentzVector& H);

  float TTX0_SC_ME_PDF(Interpretation& i);
  float TTX0_PS_ME_PDF(Interpretation& i);
  
  float Interpolate(TH1F* histo, float value);

 
private:

  float TTWHLikelihood(float mthad, float mtlep, float mhiggs, float mwhad);
  float TTWHishLikelihood(float mthad, float mtlep, float mhiggs, float mwhad);

  TFile* file;
  TH1F* h_CSV_b;
  TH1F* h_CSV_l_w_c;
  TH1F* h_M_Higgs_reco;
  TH1F* h_M_BB_reco;
  TH1F* h_M_TopHad_reco;
  TH1F* h_M_TopLep_reco;
  TH1F* h_M_WHad_reco;
  TH1F* h_M_Higgs_best;
  TH1F* h_M_TopHad_best;
  TH1F* h_M_TopLep_best;
  TH1F* h_M_WHad_best;

  float tiny_likelihood;

  MECalculator me;
};

#endif
