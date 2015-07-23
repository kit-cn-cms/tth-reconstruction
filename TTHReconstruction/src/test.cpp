#include "TStopwatch.h"
#include "TChain.h"
#include "TBranch.h"
#include "TFile.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include <iostream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>
#include "ReconstructionTester.hpp"

using namespace std;
typedef vector<TLorentzVector> LVs;
typedef TLorentzVector LV;

TLorentzVector getLV(float pt,float eta,float phi,float e){
  TLorentzVector v;
  v.SetPtEtaPhiE(pt,eta,phi,e);
  return v;
}

TLorentzVector getLV(float pt,float eta,float phi){
  TLorentzVector v;
  v.SetPtEtaPhiM(pt,eta,phi,0);
  return v;
}

LVs getLVs(uint n, float* pt,float* eta,float* phi,float* e){
  vector<TLorentzVector> vs;
  for(uint i=0; i<n; i++){
    TLorentzVector v= getLV(pt[i],eta[i],phi[i],e[i]);
    vs.push_back(v);
  }

  return vs;
}


void test(){
  TH1F::SetDefaultSumw2();
  
  TChain* chain = new TChain("MVATree");
  char* filenames = getenv ("FILENAMES");
  char* outfilename = getenv ("OUTFILENAME");
  int maxevents = atoi(getenv ("MAXEVENTS"));
  int skipevents = atoi(getenv ("SKIPEVENTS"));
  string buf;
  stringstream ss(filenames); 
  while (ss >> buf){
    chain->Add(buf.c_str());
  }
  chain->SetBranchStatus("*",0);

  float Weight;
  chain->SetBranchAddress("Weight",&Weight);
  int N_Jets;
  chain->SetBranchAddress("N_Jets",&N_Jets);
  int N_BTagsM;
  chain->SetBranchAddress("N_BTagsM",&N_BTagsM);
  int N_BTagsT;
  chain->SetBranchAddress("N_BTagsT",&N_BTagsT);
  int N_BTagsL;
  chain->SetBranchAddress("N_BTagsL",&N_BTagsL);
  int N_TightLeptons;
  chain->SetBranchAddress("N_TightLeptons",&N_TightLeptons);
  float* Jet_Pt = new float[20];
  chain->SetBranchAddress("Jet_Pt",Jet_Pt);
  float* Jet_Phi = new float[20];
  chain->SetBranchAddress("Jet_Phi",Jet_Phi);
  float* Jet_Eta = new float[20];
  chain->SetBranchAddress("Jet_Eta",Jet_Eta);
  float* Jet_E = new float[20];
  chain->SetBranchAddress("Jet_E",Jet_E);
  float* Jet_CSV = new float[20];
  chain->SetBranchAddress("Jet_CSV",Jet_CSV);
  float* Jet_Flav = new float[20];
  chain->SetBranchAddress("Jet_Flav",Jet_Flav);
  float Evt_Pt_PrimaryLepton;
  chain->SetBranchAddress("Evt_Pt_PrimaryLepton",&Evt_Pt_PrimaryLepton);
  float Evt_Eta_PrimaryLepton;
  chain->SetBranchAddress("Evt_Eta_PrimaryLepton",&Evt_Eta_PrimaryLepton);
  float Evt_Phi_PrimaryLepton;
  chain->SetBranchAddress("Evt_Phi_PrimaryLepton",&Evt_Phi_PrimaryLepton);
  float Evt_E_PrimaryLepton;
  chain->SetBranchAddress("Evt_E_PrimaryLepton",&Evt_E_PrimaryLepton);
  float Evt_Pt_MET;
  chain->SetBranchAddress("Evt_Pt_MET",&Evt_Pt_MET);
  float Evt_Phi_MET;
  chain->SetBranchAddress("Evt_Phi_MET",&Evt_Phi_MET);
  int GenEvt_I_TTPlusCC;
  chain->SetBranchAddress("GenEvt_I_TTPlusCC",&GenEvt_I_TTPlusCC);
  int GenEvt_I_TTPlusBB;
  chain->SetBranchAddress("GenEvt_I_TTPlusBB",&GenEvt_I_TTPlusBB);
  int N_GenTopHad;
  chain->SetBranchAddress("N_GenTopHad",&N_GenTopHad);
  float* GenTopHad_Pt = new float[20];
  chain->SetBranchAddress("GenTopHad_Pt",GenTopHad_Pt);
  float* GenTopHad_Eta = new float[20];
  chain->SetBranchAddress("GenTopHad_Eta",GenTopHad_Eta);
  float* GenTopHad_Phi = new float[20];
  chain->SetBranchAddress("GenTopHad_Phi",GenTopHad_Phi);
  float* GenTopHad_B_Pt = new float[20];
  chain->SetBranchAddress("GenTopHad_B_Pt",GenTopHad_B_Pt);
  float* GenTopHad_Q1_Pt = new float[20];
  chain->SetBranchAddress("GenTopHad_Q1_Pt",GenTopHad_Q1_Pt);
  float* GenTopHad_Q2_Pt = new float[20];
  chain->SetBranchAddress("GenTopHad_Q2_Pt",GenTopHad_Q2_Pt);
  float* GenTopHad_B_Eta = new float[20];
  chain->SetBranchAddress("GenTopHad_B_Eta",GenTopHad_B_Eta);
  float* GenTopHad_Q1_Eta = new float[20];
  chain->SetBranchAddress("GenTopHad_Q1_Eta",GenTopHad_Q1_Eta);
  float* GenTopHad_Q2_Eta = new float[20];
  chain->SetBranchAddress("GenTopHad_Q2_Eta",GenTopHad_Q2_Eta);
  float* GenTopHad_B_Phi = new float[20];
  chain->SetBranchAddress("GenTopHad_B_Phi",GenTopHad_B_Phi);
  float* GenTopHad_Q1_Phi = new float[20];
  chain->SetBranchAddress("GenTopHad_Q1_Phi",GenTopHad_Q1_Phi);
  float* GenTopHad_Q2_Phi = new float[20];
  chain->SetBranchAddress("GenTopHad_Q2_Phi",GenTopHad_Q2_Phi);
  int N_GenTopLep;
  chain->SetBranchAddress("N_GenTopLep",&N_GenTopLep);
  float* GenTopLep_Pt = new float[20];
  chain->SetBranchAddress("GenTopLep_Pt",GenTopLep_Pt);
  float* GenTopLep_Eta = new float[20];
  chain->SetBranchAddress("GenTopLep_Eta",GenTopLep_Eta);
  float* GenTopLep_Phi = new float[20];
  chain->SetBranchAddress("GenTopLep_Phi",GenTopLep_Phi);
  float* GenTopLep_B_Pt = new float[20];
  chain->SetBranchAddress("GenTopLep_B_Pt",GenTopLep_B_Pt);
  float* GenTopLep_Lep_Pt = new float[20];
  chain->SetBranchAddress("GenTopLep_Lep_Pt",GenTopLep_Lep_Pt);
  float* GenTopLep_Nu_Pt = new float[20];
  chain->SetBranchAddress("GenTopLep_Nu_Pt",GenTopLep_Nu_Pt);
  float* GenTopLep_B_Eta = new float[20];
  chain->SetBranchAddress("GenTopLep_B_Eta",GenTopLep_B_Eta);
  float* GenTopLep_Lep_Eta = new float[20];
  chain->SetBranchAddress("GenTopLep_Lep_Eta",GenTopLep_Lep_Eta);
  float* GenTopLep_Nu_Eta = new float[20];
  chain->SetBranchAddress("GenTopLep_Nu_Eta",GenTopLep_Nu_Eta);
  float* GenTopLep_B_Phi = new float[20];
  chain->SetBranchAddress("GenTopLep_B_Phi",GenTopLep_B_Phi);
  float* GenTopLep_Lep_Phi = new float[20];
  chain->SetBranchAddress("GenTopLep_Lep_Phi",GenTopLep_Lep_Phi);
  float* GenTopLep_Nu_Phi = new float[20];
  chain->SetBranchAddress("GenTopLep_Nu_Phi",GenTopLep_Nu_Phi);
  float GenHiggs_Pt;
  chain->SetBranchAddress("GenHiggs_Pt",&GenHiggs_Pt);
  float GenHiggs_Eta;
  chain->SetBranchAddress("GenHiggs_Eta",&GenHiggs_Eta);
  float GenHiggs_Phi;
  chain->SetBranchAddress("GenHiggs_Phi",&GenHiggs_Phi);
  float GenHiggs_B1_Pt;
  chain->SetBranchAddress("GenHiggs_B1_Pt",&GenHiggs_B1_Pt);
  float GenHiggs_B2_Pt;
  chain->SetBranchAddress("GenHiggs_B2_Pt",&GenHiggs_B2_Pt);
  float GenHiggs_B1_Eta;
  chain->SetBranchAddress("GenHiggs_B1_Eta",&GenHiggs_B1_Eta);
  float GenHiggs_B2_Eta;
  chain->SetBranchAddress("GenHiggs_B2_Eta",&GenHiggs_B2_Eta);
  float GenHiggs_B1_Phi;
  chain->SetBranchAddress("GenHiggs_B1_Phi",&GenHiggs_B1_Phi);
  float GenHiggs_B2_Phi;
  chain->SetBranchAddress("GenHiggs_B2_Phi",&GenHiggs_B2_Phi);
  int N_AdditionalGenBJets;
  chain->SetBranchAddress("N_AdditionalGenBJets",&N_AdditionalGenBJets);
  float* AdditionalGenBJet_Pt = new float[20];
  chain->SetBranchAddress("AdditionalGenBJet_Pt",AdditionalGenBJet_Pt);
  float* AdditionalGenBJet_Eta = new float[20];
  chain->SetBranchAddress("AdditionalGenBJet_Eta",AdditionalGenBJet_Eta);
  float* AdditionalGenBJet_Phi = new float[20];
  chain->SetBranchAddress("AdditionalGenBJet_Phi",AdditionalGenBJet_Phi);
  float* AdditionalGenBJet_E = new float[20];
  chain->SetBranchAddress("AdditionalGenBJet_E",AdditionalGenBJet_E);


  
  vector<string> tags;
  tags.push_back("TTWChi2_tagged");
  tags.push_back("TTWHChi2_tagged");
  //  tags.push_back("TTWChi2_tagged_higgspt");
  //  tags.push_back("TTWChi2_tagged_higgsjetpt");
  tags.push_back("TTWLikelihood_tagged");
  tags.push_back("TTWHLikelihood_tagged");
  tags.push_back("TTWishLikelihood_tagged");
  tags.push_back("TTWHishLikelihood_tagged");

  //  tags.push_back("TTWLikelihood");
  //  tags.push_back("TTWHLikelihood");



  ReconstructionTester tester(tags,outfilename);

  // loop
  long nentries = chain->GetEntries(); 
  cout << "total number of events: " << nentries << endl;
  TStopwatch watch;
  watch.Start();
  int nselected=0;
  for (long iEntry=0;iEntry<nentries;iEntry++) {
    if(maxevents>0&&iEntry>maxevents) break;
    if(iEntry%10000==0||iEntry<100){
      cout << "analyzing event " << iEntry << endl;
      watch.Print();
      watch.Continue();
    }

    chain->GetEntry(iEntry+skipevents); 
    // selection
    if(N_Jets<6||N_BTagsM<4||N_GenTopHad!=1||N_GenTopLep!=1) continue;
    nselected++;
    if(nselected%100==0){
      cout << "selected events " << nselected << endl;
    }
    TLorentzVector vHiggs_true=getLV(GenHiggs_Pt,GenHiggs_Eta,GenHiggs_Phi);
    TLorentzVector vTopHad_true=getLV(GenTopHad_Pt[0],GenTopHad_Eta[0],GenTopHad_Phi[0]);
    TLorentzVector vTopLep_true=getLV(GenTopLep_Pt[0],GenTopLep_Eta[0],GenTopLep_Phi[0]);
    TLorentzVector vQ1_true=getLV(GenTopHad_Q1_Pt[0],GenTopHad_Q1_Eta[0],GenTopHad_Q1_Phi[0]);
    TLorentzVector vQ2_true=getLV(GenTopHad_Q2_Pt[0],GenTopHad_Q2_Eta[0],GenTopHad_Q2_Phi[0]);
    TLorentzVector vBHad_true=getLV(GenTopHad_B_Pt[0],GenTopHad_B_Eta[0],GenTopHad_B_Phi[0]);
    TLorentzVector vNu_true=getLV(GenTopLep_Nu_Pt[0],GenTopLep_Nu_Eta[0],GenTopLep_Nu_Phi[0]);
    TLorentzVector vLep_true=getLV(GenTopLep_Lep_Pt[0],GenTopLep_Lep_Eta[0],GenTopLep_Lep_Phi[0]);
    TLorentzVector vBLep_true=getLV(GenTopLep_B_Pt[0],GenTopLep_B_Eta[0],GenTopLep_B_Phi[0]);
    TLorentzVector vB1_true;
    TLorentzVector vB2_true;
    if(GenHiggs_B1_Pt>0.1){
      vB1_true=getLV(GenHiggs_B1_Pt,GenHiggs_B1_Eta,GenHiggs_B1_Phi);
      vB2_true=getLV(GenHiggs_B2_Pt,GenHiggs_B2_Eta,GenHiggs_B2_Phi);
    }
    else if(N_AdditionalGenBJets>=2){
      vB1_true=getLV(AdditionalGenBJet_Pt[0],AdditionalGenBJet_Eta[0],AdditionalGenBJet_Phi[0],AdditionalGenBJet_E[0]);
      vB2_true=getLV(AdditionalGenBJet_Pt[1],AdditionalGenBJet_Eta[1],AdditionalGenBJet_Phi[1],AdditionalGenBJet_E[1]);
    }
    else continue;


    LVs jetvecs = getLVs(N_Jets,Jet_Pt,Jet_Eta,Jet_Phi,Jet_E);
    vector<float> jetcsvs;
    for(uint i=0; i<N_Jets; i++){
      jetcsvs.push_back(Jet_CSV[i]);
    }
    LV lepvec = getLV(Evt_Pt_PrimaryLepton,Evt_Eta_PrimaryLepton,Evt_Phi_PrimaryLepton,Evt_E_PrimaryLepton);
    TVector2 metvec;
    metvec.SetMagPhi(Evt_Pt_MET,Evt_Phi_MET);
    tester.Analyze(jetvecs,jetcsvs,lepvec,metvec,
		   vBHad_true,vQ1_true,vQ2_true,vBLep_true,
		   vLep_true,vNu_true,vB1_true,vB2_true);
  }

}

int main(){
  test();
  return 0;
}
