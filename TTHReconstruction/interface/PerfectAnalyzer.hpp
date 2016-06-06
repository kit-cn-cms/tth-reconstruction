#ifndef PERFECTANALYZER_HPP
#define PERFECTANALYZER_HPP

#include "ReconstructionQuality.hpp"
#include "ReconstructionMCMatching.hpp"
#include "InterpretationGenerator.hpp"
#include "Interpretation.hpp"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include <map>

class PerfectAnalyzer{
public:
  PerfectAnalyzer(std::string outfilename="test");
  ~PerfectAnalyzer();
  void Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
	       const TLorentzVector& lepvec, const TVector2& metvec,
	       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, const TLorentzVector& vQ2_true, 
	       const TLorentzVector& vBLep_true, const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
	       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true);

private:
  ReconstructionQuality quality;
  TFile* outfile;

  int nevents;
  int nselected;

  TH1F* h_perfect_tth_me;
  TH1F* h_perfect_tth_ps_me;
  TH1F* h_perfect_tth_sc_me;
  TH1F* h_perfect_tth_ratio_me;
  TH1F* h_perfect_tth_ratio_me_pdf;
  TH1F* h_perfect_ttbb_off_me;
  TH1F* h_perfect_tth_nd_me;
  TH1F* h_perfect_ttbb_me;
  TH1F* h_perfect_ratio_me;

  TH1F* h_perfect_h_pt;
  TH1F* h_perfect_deta_tt;
  TH1F* h_perfect_dr_th;



  TH1F* h_njets;
  TH1F* h_ntags;
  TH1F* h_HT;
  TH1F* h_mass;


};

#endif
