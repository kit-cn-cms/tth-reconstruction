#ifndef MECALCULATOR
#define MECALCULATOR

#include "tthProcess.h"
#include "TLorentzVector.h"
#include <iostream>
#include <iomanip> 
#include <assert.h> 

class MECalculator{
public:
  MECalculator();

  float GetMEsq(const TLorentzVector & top, const TLorentzVector & tophad, const TLorentzVector & higgs);
  float test();
  

private:  
  tthProcess tthME;
 
};

#endif
