#ifndef RECONSTRUCTIONTESTER
#define RECONSTRUCTIONTESTER

#include "ReconstructionQuality.hpp"
#include "ReconstructionMCMatching.hpp"
#include "InterpretationGenerator.hpp"
#include "Interpretation.hpp"
#include "TLorentzVector.h"
#include <map>

class ReconstructionTester{
public:
  ReconstructionTester(std::vector<std::string> tags, std::string outfilename="test");
  ~ReconstructionTester();
  void Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
	       const TLorentzVector& lepvec, const TVector2& metvec,
	       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, const TLorentzVector& vQ2_true, 
	       const TLorentzVector& vBLep_true, const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
	       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true);
  void PlotInt(std::string tag, Interpretation* i, std::string suffix);
  void InitHisto(std::string tag,std::string name,int nbins,float xmin, float xmax);
  void FillHisto(std::string tag,std::string name,float value);
private:
  std::vector<std::string> tags;
  std::string outfilename;
    
  ReconstructionMCMatching mcmatcher;
  ReconstructionQuality quality;
  InterpretationGenerator generator;
  std::map<std::string , std::map<std::string,TH1F* > > allHistos;
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
