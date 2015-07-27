#ifndef RECONSTRUCTIONTESTER
#define RECONSTRUCTIONTESTER

#include "ReconstructionQuality.hpp"
#include "ReconstructionMCMatching.hpp"
#include "InterpretationGenerator.hpp"
#include "Interpretation.hpp"
#include "TLorentzVector.h"
#include "TH1.h"
#include "TH2.h"
#include <map>

class ReconstructionTester{
public:
  ReconstructionTester(std::vector<std::string> tags, std::string outfilename="test",bool testMEs=false);
  ~ReconstructionTester();
  void Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
	       const TLorentzVector& lepvec, const TVector2& metvec,
	       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, const TLorentzVector& vQ2_true, 
	       const TLorentzVector& vBLep_true, const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
	       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true);
  void PlotInt(std::string tag, Interpretation* i, std::string suffix);
  void InitHisto(std::string tag,std::string name,int nbins,float xmin, float xmax);
  void FillHisto(std::string tag,std::string name,float value);
  void Init2dHisto(std::string tag,std::string name,int nbinsx,float xmin, float xmax,int nbinsy,float ymin, float ymax);
  void Fill2dHisto(std::string tag,std::string name,float xvalue,float yvalue);

private:
  std::vector<std::string> tags;
  std::string outfilename;
  bool testMEs;
  TFile* outfile;
  ReconstructionMCMatching mcmatcher;
  ReconstructionQuality quality;
  InterpretationGenerator generator;
  std::map<std::string , std::map<std::string,TH1* > > allHistos;
  std::map<std::string,int> foundH;
  std::map<std::string,int> foundAll;
  std::map<std::string,int> foundWHad;

  int ntotal;
  int nmatchAll;
  int nmatchH;
  int nmatchWhad;
  int nmatchBhad;
  int nmatchBlep;
 
};

#endif
