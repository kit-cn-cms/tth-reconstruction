//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.4.0, 2016-05-12
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef Parameters_bsm_H
#define Parameters_bsm_H

#include <complex> 

#include "read_slha.h"
using namespace std; 

class Parameters_bsm
{
  public:

  //static Parameters_bsm * getInstance(); 

    // Define "zero"
    double zero, ZERO; 
    // Model parameters independent of aS
    double mdl_WX2, mdl_WX1, mdl_WX0, mdl_WW, mdl_WZ, mdl_WT, mdl_ymtau,
        mdl_ymt, mdl_ymb, aS, mdl_Gf, aEWM1, mdl_MX2, mdl_MX1, mdl_MX0, mdl_MZ,
        mdl_MTA, mdl_MT, mdl_MB, mdl_kza, mdl_kw, mdl_kz, mdl_ka, mdl_kg,
        mdl_kl, mdl_kq3, mdl_kq, mdl_kz5, mdl_kz3, mdl_kz1, mdl_kw5, mdl_kw4,
        mdl_kw3, mdl_kw2, mdl_kw1, mdl_klb, mdl_kla, mdl_kqb, mdl_kqa,
        mdl_kAAgg, mdl_kHHgg, mdl_kHdwI, mdl_kHdwR, mdl_kHdz, mdl_kHda,
        mdl_kAww, mdl_kHww, mdl_kAzz, mdl_kHzz, mdl_kAgg, mdl_kHgg, mdl_kAza,
        mdl_kHza, mdl_kAaa, mdl_kHaa, mdl_kAll, mdl_kHll, mdl_kAbb, mdl_kHbb,
        mdl_kAtt, mdl_kHtt, mdl_kSM, mdl_ca, mdl_Lambda, mdl_cabi,
        mdl_cos__cabi, mdl_sin__cabi, mdl_ca__exp__2, mdl_sa, mdl_MZ__exp__2,
        mdl_MZ__exp__4, mdl_sqrt__2, mdl_nb__2__exp__0_75, mdl_MX0__exp__2,
        mdl_ca__exp__4, mdl_ca__exp__3, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee,
        mdl_MW__exp__2, mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_ad, mdl_al,
        mdl_an, mdl_au, mdl_bd, mdl_bl, mdl_bn, mdl_bu, mdl_gwwz, mdl_g1,
        mdl_gw, mdl_vev, mdl_ee__exp__2, mdl_gAaa, mdl_vev__exp__2,
        mdl_cw__exp__2, mdl_gAza, mdl_gHaa, mdl_gHza, mdl_lam, mdl_yb, mdl_yt,
        mdl_ytau, mdl_muH, mdl_sw__exp__2;
    std::complex<double> mdl_CKM1x1, mdl_CKM1x2, mdl_CKM2x1, mdl_CKM2x2,
        mdl_complexi, mdl_kHdw, mdl_conjg__CKM1x1, mdl_conjg__CKM1x2,
        mdl_conjg__CKM2x1, mdl_conjg__CKM2x2, mdl_conjg__kHdw;
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G, mdl_G__exp__2, mdl_gAAgg, mdl_gAgg, mdl_gHgg,
        mdl_gHHgg;
    // Model couplings independent of aS
    std::complex<double> GC_103, GC_104; 
    // Model couplings dependent on aS
    std::complex<double> GC_13, GC_6, GC_11, GC_10, GC_64, GC_65, GC_62, GC_7; 

    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha); 
    // Set couplings that are unchanged during the run
    void setIndependentCouplings(); 
    // Set parameters that are changed event by event
    void setDependentParameters(); 
    // Set couplings that are changed event by event
    void setDependentCouplings(); 

    // Print parameters that are unchanged during the run
    void printIndependentParameters(); 
    // Print couplings that are unchanged during the run
    void printIndependentCouplings(); 
    // Print parameters that are changed event by event
    void printDependentParameters(); 
    // Print couplings that are changed event by event
    void printDependentCouplings(); 


  private:
    static Parameters_bsm * instance; 
}; 

#endif  // Parameters_bsm_H

