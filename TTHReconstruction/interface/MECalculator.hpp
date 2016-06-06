#ifndef MECALCULATOR
#define MECALCULATOR

#include "tthProcess.h"
#include "LHAPDF/LHAPDF.h"
#include "ttbbProcess.h"
#include "tthbbProcess.h"
#include "ggttx0Process.h"
#include "qqttx0Process.h"
#include "TLorentzVector.h"
#include <iostream>
#include <iomanip> 
#include <assert.h> 

class MECalculator{
public:
  MECalculator();

  float GetTTHMEsq(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & higgs);
  float GetTTX0MEsq(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & higgs,bool pseudoscalar=false);
  float GetTTX0MEsqWpdf(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & higgs,bool pseudoscalar=false);

  float GetTTBBMEsq_onshell(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & b, const TLorentzVector & bbar);
  float GetTTBBMEsq_offshell(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & b, const TLorentzVector & bbar);
  float GetTTHBBMEsq(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & b, const TLorentzVector & bbar);
  

private:  
  tthProcess tthME;
  ttbbProcess ttbbME;
  tthbbProcess tthbbME;
  ggttx0Process ggttx0ME_sc;
  qqttx0Process qqttx0ME_sc;
  ggttx0Process ggttx0ME_ps;
  qqttx0Process qqttx0ME_ps;
  LHAPDF::PDF* pdf;
};

#endif
