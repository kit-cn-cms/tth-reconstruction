//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <iostream> 
#include <iomanip> 
#include "Parameters_bsm.h"

// Initialize static instance
Parameters_bsm * Parameters_bsm::instance = 0; 

// Function to get static instance - only one instance per program
/*Parameters_bsm * Parameters_bsm::getInstance()
{
  if (instance == 0)
    instance = new Parameters_bsm(); 

  return instance; 
  }*/

void Parameters_bsm::setIndependentParameters(SLHAReader& slha)
{
  // Define "zero"
  zero = 0; 
  ZERO = 0; 
  // Prepare a vector for indices
  vector<int> indices(2, 0); 
  mdl_WX2 = slha.get_block_entry("decay", 5000002, 4.070000e-03); 
  mdl_WX1 = slha.get_block_entry("decay", 5000001, 4.070000e-03); 
  mdl_WX0 = slha.get_block_entry("decay", 5000000, 4.070000e-03); 
  mdl_WW = slha.get_block_entry("decay", 24, 2.085000e+00); 
  mdl_WZ = slha.get_block_entry("decay", 23, 2.495200e+00); 
  mdl_WT = slha.get_block_entry("decay", 6, 1.508336e+00); 
  mdl_ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00); 
  mdl_ymt = slha.get_block_entry("yukawa", 6, 1.720000e+02); 
  mdl_ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00); 
  aS = slha.get_block_entry("sminputs", 3, 1.184000e-01); 
  mdl_Gf = slha.get_block_entry("sminputs", 2, 1.166370e-05); 
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.279000e+02); 
  mdl_MX2 = slha.get_block_entry("mass", 5000002, 1.250000e+02); 
  mdl_MX1 = slha.get_block_entry("mass", 5000001, 1.250000e+02); 
  mdl_MX0 = slha.get_block_entry("mass", 5000000, 1.250000e+02); 
  mdl_MZ = slha.get_block_entry("mass", 23, 9.118760e+01); 
  mdl_MTA = slha.get_block_entry("mass", 15, 1.777000e+00); 
  mdl_MT = slha.get_block_entry("mass", 6, 1.720000e+02); 
  mdl_MB = slha.get_block_entry("mass", 5, 4.700000e+00); 
  mdl_kza = slha.get_block_entry("frblock", 45, 0.000000e+00); 
  mdl_kw = slha.get_block_entry("frblock", 44, 1.000000e+00); 
  mdl_kz = slha.get_block_entry("frblock", 43, 1.000000e+00); 
  mdl_ka = slha.get_block_entry("frblock", 42, 1.000000e+00); 
  mdl_kg = slha.get_block_entry("frblock", 41, 1.000000e+00); 
  mdl_kl = slha.get_block_entry("frblock", 40, 1.000000e+00); 
  mdl_kq3 = slha.get_block_entry("frblock", 39, 1.000000e+00); 
  mdl_kq = slha.get_block_entry("frblock", 38, 1.000000e+00); 
  mdl_kz5 = slha.get_block_entry("frblock", 37, 0.000000e+00); 
  mdl_kz3 = slha.get_block_entry("frblock", 36, 1.000000e+00); 
  mdl_kz1 = slha.get_block_entry("frblock", 35, 0.000000e+00); 
  mdl_kw5 = slha.get_block_entry("frblock", 34, 0.000000e+00); 
  mdl_kw4 = slha.get_block_entry("frblock", 33, 0.000000e+00); 
  mdl_kw3 = slha.get_block_entry("frblock", 32, 0.000000e+00); 
  mdl_kw2 = slha.get_block_entry("frblock", 31, 1.000000e+00); 
  mdl_kw1 = slha.get_block_entry("frblock", 30, 1.000000e+00); 
  mdl_klb = slha.get_block_entry("frblock", 29, 1.000000e+00); 
  mdl_kla = slha.get_block_entry("frblock", 28, 1.000000e+00); 
  mdl_kqb = slha.get_block_entry("frblock", 27, 1.000000e+00); 
  mdl_kqa = slha.get_block_entry("frblock", 26, 1.000000e+00); 
  mdl_kAAgg = slha.get_block_entry("frblock", 25, 1.000000e+00); 
  mdl_kHHgg = slha.get_block_entry("frblock", 24, 1.000000e+00); 
  mdl_kHdwI = slha.get_block_entry("frblock", 23, 0.000000e+00); 
  mdl_kHdwR = slha.get_block_entry("frblock", 22, 0.000000e+00); 
  mdl_kHdz = slha.get_block_entry("frblock", 21, 0.000000e+00); 
  mdl_kHda = slha.get_block_entry("frblock", 20, 0.000000e+00); 
  mdl_kAww = slha.get_block_entry("frblock", 19, 0.000000e+00); 
  mdl_kHww = slha.get_block_entry("frblock", 18, 0.000000e+00); 
  mdl_kAzz = slha.get_block_entry("frblock", 17, 0.000000e+00); 
  mdl_kHzz = slha.get_block_entry("frblock", 16, 0.000000e+00); 
  mdl_kAgg = slha.get_block_entry("frblock", 15, 1.000000e+00); 
  mdl_kHgg = slha.get_block_entry("frblock", 14, 1.000000e+00); 
  mdl_kAza = slha.get_block_entry("frblock", 13, 1.000000e+00); 
  mdl_kHza = slha.get_block_entry("frblock", 12, 1.000000e+00); 
  mdl_kAaa = slha.get_block_entry("frblock", 11, 1.000000e+00); 
  mdl_kHaa = slha.get_block_entry("frblock", 10, 1.000000e+00); 
  mdl_kAll = slha.get_block_entry("frblock", 9, 1.000000e+00); 
  mdl_kHll = slha.get_block_entry("frblock", 8, 1.000000e+00); 
  mdl_kAbb = slha.get_block_entry("frblock", 7, 1.000000e+00); 
  mdl_kHbb = slha.get_block_entry("frblock", 6, 1.000000e+00); 
  mdl_kAtt = slha.get_block_entry("frblock", 5, 1.000000e+00); 
  mdl_kHtt = slha.get_block_entry("frblock", 4, 1.000000e+00); 
  mdl_kSM = slha.get_block_entry("frblock", 3, 1.000000e+00); 
  mdl_ca = slha.get_block_entry("frblock", 2, 1.000000e+00); 
  mdl_Lambda = slha.get_block_entry("frblock", 1, 1.000000e+03); 
  mdl_cabi = slha.get_block_entry("ckmblock", 1, 2.277360e-01); 
  mdl_cos__cabi = cos(mdl_cabi); 
  mdl_CKM1x1 = mdl_cos__cabi; 
  mdl_sin__cabi = sin(mdl_cabi); 
  mdl_CKM1x2 = mdl_sin__cabi; 
  mdl_CKM2x1 = -mdl_sin__cabi; 
  mdl_CKM2x2 = mdl_cos__cabi; 
  mdl_ca__exp__2 = pow(mdl_ca, 2.); 
  mdl_sa = sqrt(1. - mdl_ca__exp__2); 
  mdl_complexi = std::complex<double> (0., 1.); 
  mdl_kHdw = mdl_complexi * mdl_kHdwI + mdl_kHdwR; 
  mdl_MZ__exp__2 = pow(mdl_MZ, 2.); 
  mdl_MZ__exp__4 = pow(mdl_MZ, 4.); 
  mdl_sqrt__2 = sqrt(2.); 
  mdl_nb__2__exp__0_75 = pow(2., 0.75); 
  mdl_MX0__exp__2 = pow(mdl_MX0, 2.); 
  mdl_ca__exp__4 = pow(mdl_ca, 4.); 
  mdl_ca__exp__3 = pow(mdl_ca, 3.); 
  mdl_conjg__CKM1x1 = conj(mdl_CKM1x1); 
  mdl_conjg__CKM1x2 = conj(mdl_CKM1x2); 
  mdl_conjg__CKM2x1 = conj(mdl_CKM2x1); 
  mdl_conjg__CKM2x2 = conj(mdl_CKM2x2); 
  mdl_conjg__kHdw = conj(mdl_kHdw); 
  mdl_aEW = 1./aEWM1; 
  mdl_MW = sqrt(mdl_MZ__exp__2/2. + sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = sqrt(mdl_aEW); 
  mdl_ee = 2. * mdl_sqrt__aEW * sqrt(M_PI); 
  mdl_MW__exp__2 = pow(mdl_MW, 2.); 
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2; 
  mdl_cw = sqrt(1. - mdl_sw2); 
  mdl_sqrt__sw2 = sqrt(mdl_sw2); 
  mdl_sw = mdl_sqrt__sw2; 
  mdl_ad = (mdl_ee * (-0.5 + (2. * mdl_sw2)/3.))/(2. * mdl_cw * mdl_sw); 
  mdl_al = (mdl_ee * (-0.5 + 2. * mdl_sw2))/(2. * mdl_cw * mdl_sw); 
  mdl_an = mdl_ee/(4. * mdl_cw * mdl_sw); 
  mdl_au = (mdl_ee * (0.5 - (4. * mdl_sw2)/3.))/(2. * mdl_cw * mdl_sw); 
  mdl_bd = -mdl_ee/(4. * mdl_cw * mdl_sw); 
  mdl_bl = -mdl_ee/(4. * mdl_cw * mdl_sw); 
  mdl_bn = mdl_ee/(4. * mdl_cw * mdl_sw); 
  mdl_bu = mdl_ee/(4. * mdl_cw * mdl_sw); 
  mdl_gwwz = -((mdl_cw * mdl_ee)/mdl_sw); 
  mdl_g1 = mdl_ee/mdl_cw; 
  mdl_gw = mdl_ee/mdl_sw; 
  mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee; 
  mdl_ee__exp__2 = pow(mdl_ee, 2.); 
  mdl_gAaa = mdl_ee__exp__2/(3. * pow(M_PI, 2.) * mdl_vev); 
  mdl_vev__exp__2 = pow(mdl_vev, 2.); 
  mdl_cw__exp__2 = pow(mdl_cw, 2.); 
  mdl_gAza = ((-5. + 8. * mdl_cw__exp__2) * sqrt(mdl_ee__exp__2 * mdl_Gf *
      mdl_MZ__exp__2))/(6. * mdl_nb__2__exp__0_75 * pow(M_PI, 2.) * mdl_vev);
  mdl_gHaa = (47. * mdl_ee__exp__2)/(72. * pow(M_PI, 2.) * mdl_vev); 
  mdl_gHza = ((-13. + 94. * mdl_cw__exp__2) * sqrt(mdl_ee__exp__2 * mdl_Gf *
      mdl_MZ__exp__2))/(36. * mdl_nb__2__exp__0_75 * pow(M_PI, 2.) * mdl_vev);
  mdl_lam = mdl_MX0__exp__2/(2. * mdl_vev__exp__2); 
  mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_vev; 
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev; 
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev; 
  mdl_muH = sqrt(mdl_lam * mdl_vev__exp__2); 
  mdl_sw__exp__2 = pow(mdl_sw, 2.); 
}
void Parameters_bsm::setIndependentCouplings()
{
  GC_103 = -((mdl_ca * mdl_complexi * mdl_kHtt * mdl_yt)/mdl_sqrt__2); 
  GC_104 = (mdl_kAtt * mdl_sa * mdl_yt)/mdl_sqrt__2; 
}
void Parameters_bsm::setDependentParameters()
{
  mdl_sqrt__aS = sqrt(aS); 
  G = 2. * mdl_sqrt__aS * sqrt(M_PI); 
  mdl_G__exp__2 = pow(G, 2.); 
  mdl_gAAgg = mdl_G__exp__2/(8. * pow(M_PI, 2.) * mdl_vev__exp__2); 
  mdl_gAgg = mdl_G__exp__2/(8. * pow(M_PI, 2.) * mdl_vev); 
  mdl_gHgg = -mdl_G__exp__2/(12. * pow(M_PI, 2.) * mdl_vev); 
  mdl_gHHgg = mdl_G__exp__2/(12. * pow(M_PI, 2.) * mdl_vev__exp__2); 
}
void Parameters_bsm::setDependentCouplings()
{
  GC_13 = -(mdl_ca * mdl_complexi * mdl_gHHgg * mdl_kHHgg); 
  GC_6 = -G; 
  GC_11 = -(mdl_ca * G * mdl_gHgg * mdl_kHgg); 
  GC_10 = -(mdl_ca * mdl_complexi * mdl_gHgg * mdl_kHgg); 
  GC_64 = (mdl_complexi * mdl_gAgg * mdl_kAgg * mdl_sa)/8.; 
  GC_65 = (G * mdl_gAgg * mdl_kAgg * mdl_sa)/4.; 
  GC_62 = (mdl_complexi * mdl_gAAgg * mdl_kAAgg * mdl_sa)/8.; 
  GC_7 = mdl_complexi * G; 
}

// Routines for printing out parameters
void Parameters_bsm::printIndependentParameters()
{
  cout <<  "HC_UFO model parameters independent of event kinematics:" << endl; 
  cout << setw(20) <<  "mdl_WX2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WX2 << endl;
  cout << setw(20) <<  "mdl_WX1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WX1 << endl;
  cout << setw(20) <<  "mdl_WX0 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WX0 << endl;
  cout << setw(20) <<  "mdl_WW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WW << endl;
  cout << setw(20) <<  "mdl_WZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WZ << endl;
  cout << setw(20) <<  "mdl_WT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_WT << endl;
  cout << setw(20) <<  "mdl_ymtau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymtau << endl;
  cout << setw(20) <<  "mdl_ymt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymt << endl;
  cout << setw(20) <<  "mdl_ymb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ymb << endl;
  cout << setw(20) <<  "aS " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aS << endl;
  cout << setw(20) <<  "mdl_Gf " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_Gf << endl;
  cout << setw(20) <<  "aEWM1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << aEWM1 << endl;
  cout << setw(20) <<  "mdl_MX2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MX2 << endl;
  cout << setw(20) <<  "mdl_MX1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MX1 << endl;
  cout << setw(20) <<  "mdl_MX0 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MX0 << endl;
  cout << setw(20) <<  "mdl_MZ " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MZ << endl;
  cout << setw(20) <<  "mdl_MTA " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MTA << endl;
  cout << setw(20) <<  "mdl_MT " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MT << endl;
  cout << setw(20) <<  "mdl_MB " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MB << endl;
  cout << setw(20) <<  "mdl_kza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kza << endl;
  cout << setw(20) <<  "mdl_kw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kw << endl;
  cout << setw(20) <<  "mdl_kz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kz << endl;
  cout << setw(20) <<  "mdl_ka " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ka << endl;
  cout << setw(20) <<  "mdl_kg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kg << endl;
  cout << setw(20) <<  "mdl_kl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kl << endl;
  cout << setw(20) <<  "mdl_kq3 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kq3 << endl;
  cout << setw(20) <<  "mdl_kq " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kq << endl;
  cout << setw(20) <<  "mdl_kz5 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kz5 << endl;
  cout << setw(20) <<  "mdl_kz3 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kz3 << endl;
  cout << setw(20) <<  "mdl_kz1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kz1 << endl;
  cout << setw(20) <<  "mdl_kw5 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kw5 << endl;
  cout << setw(20) <<  "mdl_kw4 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kw4 << endl;
  cout << setw(20) <<  "mdl_kw3 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kw3 << endl;
  cout << setw(20) <<  "mdl_kw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kw2 << endl;
  cout << setw(20) <<  "mdl_kw1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kw1 << endl;
  cout << setw(20) <<  "mdl_klb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_klb << endl;
  cout << setw(20) <<  "mdl_kla " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kla << endl;
  cout << setw(20) <<  "mdl_kqb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kqb << endl;
  cout << setw(20) <<  "mdl_kqa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kqa << endl;
  cout << setw(20) <<  "mdl_kAAgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAAgg << endl;
  cout << setw(20) <<  "mdl_kHHgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHHgg << endl;
  cout << setw(20) <<  "mdl_kHdwI " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHdwI << endl;
  cout << setw(20) <<  "mdl_kHdwR " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHdwR << endl;
  cout << setw(20) <<  "mdl_kHdz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHdz << endl;
  cout << setw(20) <<  "mdl_kHda " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHda << endl;
  cout << setw(20) <<  "mdl_kAww " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAww << endl;
  cout << setw(20) <<  "mdl_kHww " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHww << endl;
  cout << setw(20) <<  "mdl_kAzz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAzz << endl;
  cout << setw(20) <<  "mdl_kHzz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHzz << endl;
  cout << setw(20) <<  "mdl_kAgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAgg << endl;
  cout << setw(20) <<  "mdl_kHgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHgg << endl;
  cout << setw(20) <<  "mdl_kAza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAza << endl;
  cout << setw(20) <<  "mdl_kHza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHza << endl;
  cout << setw(20) <<  "mdl_kAaa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAaa << endl;
  cout << setw(20) <<  "mdl_kHaa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHaa << endl;
  cout << setw(20) <<  "mdl_kAll " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAll << endl;
  cout << setw(20) <<  "mdl_kHll " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHll << endl;
  cout << setw(20) <<  "mdl_kAbb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAbb << endl;
  cout << setw(20) <<  "mdl_kHbb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHbb << endl;
  cout << setw(20) <<  "mdl_kAtt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kAtt << endl;
  cout << setw(20) <<  "mdl_kHtt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHtt << endl;
  cout << setw(20) <<  "mdl_kSM " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kSM << endl;
  cout << setw(20) <<  "mdl_ca " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ca << endl;
  cout << setw(20) <<  "mdl_Lambda " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_Lambda << endl;
  cout << setw(20) <<  "mdl_cabi " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cabi << endl;
  cout << setw(20) <<  "mdl_cos__cabi " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cos__cabi << endl;
  cout << setw(20) <<  "mdl_CKM1x1 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM1x1 << endl;
  cout << setw(20) <<  "mdl_sin__cabi " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sin__cabi << endl;
  cout << setw(20) <<  "mdl_CKM1x2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM1x2 << endl;
  cout << setw(20) <<  "mdl_CKM2x1 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM2x1 << endl;
  cout << setw(20) <<  "mdl_CKM2x2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_CKM2x2 << endl;
  cout << setw(20) <<  "mdl_ca__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ca__exp__2 << endl;
  cout << setw(20) <<  "mdl_sa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sa << endl;
  cout << setw(20) <<  "mdl_complexi " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_complexi << endl;
  cout << setw(20) <<  "mdl_kHdw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_kHdw << endl;
  cout << setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__2 << endl;
  cout << setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MZ__exp__4 << endl;
  cout << setw(20) <<  "mdl_sqrt__2 " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__2 << endl;
  cout << setw(20) <<  "mdl_nb__2__exp__0_75 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_nb__2__exp__0_75 << endl;
  cout << setw(20) <<  "mdl_MX0__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MX0__exp__2 << endl;
  cout << setw(20) <<  "mdl_ca__exp__4 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ca__exp__4 << endl;
  cout << setw(20) <<  "mdl_ca__exp__3 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ca__exp__3 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM1x1 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM1x1 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM1x2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM1x2 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM2x1 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM2x1 << endl;
  cout << setw(20) <<  "mdl_conjg__CKM2x2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__CKM2x2 << endl;
  cout << setw(20) <<  "mdl_conjg__kHdw " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_conjg__kHdw << endl;
  cout << setw(20) <<  "mdl_aEW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_aEW << endl;
  cout << setw(20) <<  "mdl_MW " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_MW << endl;
  cout << setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__aEW << endl;
  cout << setw(20) <<  "mdl_ee " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ee << endl;
  cout << setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_MW__exp__2 << endl;
  cout << setw(20) <<  "mdl_sw2 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw2 << endl;
  cout << setw(20) <<  "mdl_cw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_cw << endl;
  cout << setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sqrt__sw2 << endl;
  cout << setw(20) <<  "mdl_sw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_sw << endl;
  cout << setw(20) <<  "mdl_ad " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ad << endl;
  cout << setw(20) <<  "mdl_al " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_al << endl;
  cout << setw(20) <<  "mdl_an " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_an << endl;
  cout << setw(20) <<  "mdl_au " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_au << endl;
  cout << setw(20) <<  "mdl_bd " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_bd << endl;
  cout << setw(20) <<  "mdl_bl " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_bl << endl;
  cout << setw(20) <<  "mdl_bn " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_bn << endl;
  cout << setw(20) <<  "mdl_bu " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_bu << endl;
  cout << setw(20) <<  "mdl_gwwz " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gwwz << endl;
  cout << setw(20) <<  "mdl_g1 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_g1 << endl;
  cout << setw(20) <<  "mdl_gw " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gw << endl;
  cout << setw(20) <<  "mdl_vev " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_vev << endl;
  cout << setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_ee__exp__2 << endl;
  cout << setw(20) <<  "mdl_gAaa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gAaa << endl;
  cout << setw(20) <<  "mdl_vev__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_vev__exp__2 << endl;
  cout << setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_cw__exp__2 << endl;
  cout << setw(20) <<  "mdl_gAza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gAza << endl;
  cout << setw(20) <<  "mdl_gHaa " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gHaa << endl;
  cout << setw(20) <<  "mdl_gHza " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gHza << endl;
  cout << setw(20) <<  "mdl_lam " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_lam << endl;
  cout << setw(20) <<  "mdl_yb " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yb << endl;
  cout << setw(20) <<  "mdl_yt " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_yt << endl;
  cout << setw(20) <<  "mdl_ytau " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_ytau << endl;
  cout << setw(20) <<  "mdl_muH " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_muH << endl;
  cout << setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_sw__exp__2 << endl;
}
void Parameters_bsm::printIndependentCouplings()
{
  cout <<  "HC_UFO model couplings independent of event kinematics:" << endl; 
  cout << setw(20) <<  "GC_103 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_103 << endl;
  cout << setw(20) <<  "GC_104 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_104 << endl;
}
void Parameters_bsm::printDependentParameters()
{
  cout <<  "HC_UFO model parameters dependent on event kinematics:" << endl; 
  cout << setw(20) <<  "mdl_sqrt__aS " <<  "= " << setiosflags(ios::scientific)
      << setw(10) << mdl_sqrt__aS << endl;
  cout << setw(20) <<  "G " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << G << endl;
  cout << setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      setiosflags(ios::scientific) << setw(10) << mdl_G__exp__2 << endl;
  cout << setw(20) <<  "mdl_gAAgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gAAgg << endl;
  cout << setw(20) <<  "mdl_gAgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gAgg << endl;
  cout << setw(20) <<  "mdl_gHgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gHgg << endl;
  cout << setw(20) <<  "mdl_gHHgg " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << mdl_gHHgg << endl;
}
void Parameters_bsm::printDependentCouplings()
{
  cout <<  "HC_UFO model couplings dependent on event kinematics:" << endl; 
  cout << setw(20) <<  "GC_13 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_13 << endl;
  cout << setw(20) <<  "GC_6 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_6 << endl;
  cout << setw(20) <<  "GC_11 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_11 << endl;
  cout << setw(20) <<  "GC_10 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_10 << endl;
  cout << setw(20) <<  "GC_64 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_64 << endl;
  cout << setw(20) <<  "GC_65 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_65 << endl;
  cout << setw(20) <<  "GC_62 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_62 << endl;
  cout << setw(20) <<  "GC_7 " <<  "= " << setiosflags(ios::scientific) <<
      setw(10) << GC_7 << endl;
}


