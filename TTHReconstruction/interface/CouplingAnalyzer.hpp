#ifndef COUPLINGANALYZER_HPP
#define COUPLINGANALYZER_HPP

#include "ReconstructionQuality2.hpp"
#include "ReconstructionMCMatching.hpp"
#include "InterpretationGenerator.hpp"
#include "Interpretation.hpp"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include "TTree.h"
#include <map>

class CouplingAnalyzer{
public:
  CouplingAnalyzer(std::string outfilename="test");
  ~CouplingAnalyzer();
  void Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
		 const TLorentzVector& lepvec, const TVector2& metvec,
		 const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, const TLorentzVector& vQ2_true, 
		 const TLorentzVector& vBLep_true, const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
		 const TLorentzVector& vB1_true, const TLorentzVector& vB2_true,
		 const TLorentzVector& vBHad_boosted,const TLorentzVector& vQ1_boosted, const TLorentzVector& vQ2_boosted, 
		 const TLorentzVector& vBLep_boosted, 
		 const TLorentzVector& vH_boosted);

private:
  ReconstructionMCMatching mcmatcher;
  ReconstructionQuality2 quality;
  InterpretationGenerator generator;
  TFile* outfile;

  int nevents;
  int nselected;

  TH1F* h_ratio_me_best_likelihood;
  TH1F* h_ratio_me_best_csv_likelihood;
  TH1F* h_ratio_me_bestish_likelihood;
  TTree* tree;
  std::map<std::string,float> vars;

  void InitVar(std::string name);
  void FillVar(std::string name,float value);

};

#endif
