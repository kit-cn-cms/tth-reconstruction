#ifndef BTAGANALYZER_HPP
#define BTAGANALYZER_HPP

#include "ReconstructionQuality.hpp"
#include "ReconstructionMCMatching.hpp"
#include "InterpretationGenerator.hpp"
#include "Interpretation.hpp"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include <map>

class BTagAnalyzer{
public:
  BTagAnalyzer(std::string outfilename="test");
  ~BTagAnalyzer();
  void Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, const std::vector<float>& jetflavs);

private:
  ReconstructionMCMatching mcmatcher;
  ReconstructionQuality quality;
  InterpretationGenerator generator;
  TFile* outfile;

  int nevents;
  int nselected;

  TH1* h_r32;
  TH1* h_r42;
  TH1* h_r43;

  TH1* h_r32_4b;
  TH1* h_r42_4b;
  TH1* h_r43_4b;

  TH1* h_r32_3b;
  TH1* h_r42_3b;
  TH1* h_r43_3b;

  TH1* h_r32_2b;
  TH1* h_r42_2b;
  TH1* h_r43_2b;

  TH1* h_ntagsT;
  TH1* h_ntagsM;
  TH1* h_ntagsL;

  TH1* h_ntagsT_4b;
  TH1* h_ntagsM_4b;
  TH1* h_ntagsL_4b;

  TH1* h_ntagsT_3b;
  TH1* h_ntagsM_3b;
  TH1* h_ntagsL_3b;

  TH1* h_ntagsT_2b;
  TH1* h_ntagsM_2b;
  TH1* h_ntagsL_2b;

  TH1* h_averageTagged;
  TH1* h_averageTagged_4b;
  TH1* h_averageTagged_3b;
  TH1* h_averageTagged_2b;

};

#endif
