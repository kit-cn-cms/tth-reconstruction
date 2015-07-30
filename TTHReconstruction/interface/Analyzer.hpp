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
  Analyzer(std::string outfilename="test",float min_combo_prob_tth=0.95,float min_combo_prob_tt=0,float min_combo_prob_ttbb=0.);
  ~Analyzer();
  void Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
	       const TLorentzVector& lepvec, const TVector2& metvec,
	       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, const TLorentzVector& vQ2_true, 
	       const TLorentzVector& vBLep_true, const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
	       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true);

private:
  std::string outfilename;
  TFile* outfile;
  ReconstructionMCMatching mcmatcher;
  ReconstructionQuality quality;
  InterpretationGenerator generator;
  int nevents;
  int nselected;
  int nmatched_tth;
  int nmatched_tt;
  int nmatched_ttbb;
  float min_combo_prob_tth;
  float min_combo_prob_tt;
  float min_combo_prob_ttbb;
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
 
};

#endif
