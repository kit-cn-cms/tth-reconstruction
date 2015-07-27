#ifndef INTERPRETATIONGENERATOR
#define INTERPRETATIONGENERATOR

#include "TLorentzVector.h"
#include "TVector2.h"
#include <vector>
#include <iostream>
#include <algorithm>
#include "Interpretation.hpp"
#include "next_combination.h"
namespace IntType{enum IntType{tth,tt};};

class InterpretationGenerator{
public:
  InterpretationGenerator(IntType::IntType type=IntType::tth, int allowedMistags=1, int maxJets=10, float minMWHad=-99999,float maxMWHad=99999,float btagCut=0.89);

  std::vector<Interpretation*> GenerateTTHInterpretations(std::vector<TLorentzVector> jetvecs, std::vector<float> jetcsvs, TLorentzVector lepvec, TVector2 metvec);

  void GetNuVecs(const TLorentzVector & lepvec, const TVector2 & metvec, TLorentzVector & nu1, TLorentzVector & nu2);

private:
  bool NextSubsetPermutation(std::vector<int>& idxs, int length);

  IntType::IntType type;
  int allowedMistags;
  int maxJets;
  float minMWHad;
  float maxMWHad;
  float btagCut;
  
 
};

#endif
