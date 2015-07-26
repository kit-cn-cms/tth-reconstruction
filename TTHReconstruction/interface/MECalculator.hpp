#ifndef MECALCULATOR
#define MECALCULATOR

#include "tthProcess.h"
#include "ttbbProcess.h"
#include "tthbbProcess.h"
#include "TLorentzVector.h"
#include <iostream>
#include <iomanip> 
#include <assert.h> 

class MECalculator{
public:
  MECalculator();

  float GetTTHMEsq(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & higgs);
  float GetTTBBMEsq_onshell(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & b, const TLorentzVector & bbar);
  float GetTTBBMEsq_offshell(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & b, const TLorentzVector & bbar);
  float GetTTHBBMEsq(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & b, const TLorentzVector & bbar);
  

private:  
  tthProcess tthME;
  ttbbProcess ttbbME;
  tthbbProcess tthbbME;
 
};

#endif
