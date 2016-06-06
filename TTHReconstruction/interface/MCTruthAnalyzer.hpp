#ifndef MCTRUTHANALYZER_HPP
#define MCTRUTHANALYZER_HPP

#include "ReconstructionQuality.hpp"
#include "ReconstructionMCMatching.hpp"
#include "InterpretationGenerator.hpp"
#include "Interpretation.hpp"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include <map>

class MCTruthAnalyzer{
public:
  MCTruthAnalyzer(std::string outfilename="test");
  ~MCTruthAnalyzer();
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

  TH1F* h_M_TopHad_reco;
  TH1F* h_M_WHad_reco;
  TH1F* h_M_Higgs_reco;
  TH1F* h_M_TopLep_reco;

  TH1F* h_M_TopHad_best;
  TH1F* h_M_WHad_best;
  TH1F* h_M_Higgs_best;
  TH1F* h_M_TopLep_best;



};

#endif
