#ifndef ANALYZER_HPP
#define ANALYZER_HPP

#include "ReconstructionQuality.hpp"
#include "ReconstructionMCMatching.hpp"
#include "InterpretationGenerator.hpp"
#include "Interpretation.hpp"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include <map>

class Analyzer{
public:
  Analyzer(std::string outfilename="test");
  ~Analyzer();
  void Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
	       const TLorentzVector& lepvec, const TVector2& metvec,
	       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, const TLorentzVector& vQ2_true, 
	       const TLorentzVector& vBLep_true, const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
	       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true);

private:
  ReconstructionMCMatching mcmatcher;
  ReconstructionQuality quality;
  InterpretationGenerator generator;
  TFile* outfile;

  int nevents;
  int nselected;
  int nmatched_tth;
  int nmatched_tt;
  int nmatched_ttbb;
  int nmatched_tth2;
  int nmatched_tt2;
  int nmatched_ttbb2;

  TH1F* h_best_tth_me_likelihood;
  TH1F* h_best_ttbb_me_likelihood;
  TH1F* h_best_ratio_me_likelihood;
  TH1F* h_best_tth_likelihood;
  TH1F* h_best_ttbb_likelihood;
  TH1F* h_best_ratio_likelihood;
  TH1F* h_best_tth_me;
  TH1F* h_best_ttbb_me;
  TH1F* h_best_ratio_me;
  TH1F* h_best_m_higgs_tthreco;
  TH1F* h_best_m_higgs_ttbbreco;
  TH1F* h_best_m_higgs_ttreco;
  TH1F* h_best_m_higgs_tthreco2;
  TH1F* h_best_m_higgs_ttbbreco2;
  TH1F* h_best_m_higgs_ttreco2;
  TH1F* h_sum_tth_me_likelihood;
  TH1F* h_sum_ttbb_me_likelihood;
  TH1F* h_sum_ratio_me_likelihood;
  TH1F* h_sum_tth_likelihood;
  TH1F* h_sum_ttbb_likelihood;
  TH1F* h_sum_ratio_likelihood;
  TH1F* h_sum_tth_me;
  TH1F* h_sum_ttbb_me;
  TH1F* h_sum_ratio_me; 
};

#endif
