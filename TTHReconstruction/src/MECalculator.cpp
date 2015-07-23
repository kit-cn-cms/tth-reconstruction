#include "MECalculator.hpp"

MECalculator::MECalculator(){
  process.initProc("data/param_card.dat");
  test();
}

float MECalculator::test(){
  cout << "ME is " << GetMEsq(TLorentzVector(2.406089e+02,5.390753e+02,-2.747567e+02,6.737322e+02),TLorentzVector(-2.602682e+01,-4.937100e+02,1.510072e+02,5.451230e+02),TLorentzVector(-2.145821e+02,-4.536532e+01,1.237495e+02,2.811448e+02)) << endl;
}

float MECalculator::GetMEsq(const TLorentzVector & top, const TLorentzVector & topbar, const TLorentzVector & higgs){
  double p4_g1[4];
  double p4_g2[4];
  double p4_top[4];
  double p4_topbar[4];
  double p4_h[4];
  
  TLorentzVector tth=top+topbar+higgs;
  double cme = tth.M();
  TLorentzVector g1(0,0,cme/2,cme/2);
  TLorentzVector g2(0,0,-cme/2,cme/2);
  TVector3 boost=tth.BoostVector();
  g1.Boost(boost);
  g2.Boost(boost);
  

  p4_g1[0]=g1.E();
  p4_g1[1]=g1.Px();
  p4_g1[2]=g1.Py();
  p4_g1[3]=g1.Pz();
  p4_g2[0]=g2.E();
  p4_g2[1]=g2.Px();
  p4_g2[2]=g2.Py();
  p4_g2[3]=g2.Pz();
  p4_top[0]=top.E();
  p4_top[1]=top.Px();
  p4_top[2]=top.Py();
  p4_top[3]=top.Pz();
  p4_topbar[0]=topbar.E();
  p4_topbar[1]=topbar.Px();
  p4_topbar[2]=topbar.Py();
  p4_topbar[3]=topbar.Pz();
  p4_h[0]=higgs.E();
  p4_h[1]=higgs.Px();
  p4_h[2]=higgs.Py();
  p4_h[3]=higgs.Pz();

  vector<double*> momenta;
  momenta.push_back(p4_g1);
  momenta.push_back(p4_g2);
  momenta.push_back(p4_top);
  momenta.push_back(p4_topbar);
  momenta.push_back(p4_h);

  process.setMomenta(momenta);

  process.sigmaKin();

  const double* matrix_elements = process.getMatrixElements();

  assert(process.nprocesses==1);
  
  return matrix_elements[0];
}
