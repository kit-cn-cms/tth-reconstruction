#include "../interface/ReconstructionQuality2.hpp"

using namespace std;

ReconstructionQuality2::ReconstructionQuality2(string filename){
  file=TFile::Open(filename.c_str());
  h_CSV_b=(TH1F*)file->Get("CSV_b");
  h_CSV_l_w_c=(TH1F*)file->Get("CSV_l_w_c");
  h_M_Higgs_reco=(TH1F*)file->Get("M_Higgs_reco");
  h_M_TopHad_reco=(TH1F*)file->Get("M_TopHad_reco");
  h_M_TopLep_reco=(TH1F*)file->Get("M_TopLep_reco");
  h_M_WHad_reco=(TH1F*)file->Get("M_WHad_reco");
  h_M_Higgs_best=(TH1F*)file->Get("M_Higgs_best");
  h_M_TopHad_best=(TH1F*)file->Get("M_TopHad_best");
  h_M_TopLep_best=(TH1F*)file->Get("M_TopLep_best");
  h_M_WHad_best=(TH1F*)file->Get("M_WHad_best");
}

float ReconstructionQuality2::GetTag(std::string tag, Interpretation& i){
   
  if(tag=="TTWHLikelihood") return TTWHLikelihood(i);
  else if(tag=="TTWHishLikelihood") return TTWHLikelihood(i);
  else if(tag=="TTX0_SC_ME") return TTX0_SC_ME(i);
  else if(tag=="TTX0_PS_ME") return TTX0_PS_ME(i);
  else{
    cout << "tag " << tag << " unknown!" << endl;
    return -1e20;
  }
}


float ReconstructionQuality2::TTWHLikelihood(float mthad, float mtlep, float mwhad, float mhiggs){
  float llh=1;
  llh*=TopHadLikelihood(mthad);
  llh*=TopLepLikelihood(mtlep);
  llh*=WHadLikelihood(mwhad);
  llh*=HiggsLikelihood(mhiggs);
  return llh;
}
float ReconstructionQuality2::TTWHLikelihood(Interpretation& i){
  if(i.HasTag("TTWHLikelihood")) return i.GetTag("TTWHLikelihood");
  float tag=TTWHLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWHLikelihood",tag);
  return tag;
}

float ReconstructionQuality2::TTWHLikelihoodTimesSCME(Interpretation& i){
  float llh=0.;
  if(i.HasTag("TTWHLikelihoodTimesSCME")) llh= i.GetTag("TTWHLikelihoodTimesSCME");
  else llh=TTWHLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  float me=TTX0_SC_ME(i);
  i.SetTag("TTWHLikelihoodTimesSCME",llh*me);
  return llh*me;
}

float ReconstructionQuality2::TTWHLikelihoodTimesPSME(Interpretation& i){
  float llh=0.;
  if(i.HasTag("TTWHLikelihoodTimesPSME")) llh= i.GetTag("TTWHLikelihoodTimesPSME");
  else llh=TTWHLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  float me=TTX0_PS_ME(i);
  i.SetTag("TTWHLikelihoodTimesPSME",llh*me);
  return llh*me;
}

float ReconstructionQuality2::TTWHLikelihoodTimesCSVTimesPSME(Interpretation& i){
  float tag=TTWHLikelihoodTimesPSME(i);
  tag*=LLikelihood(i.Q1_CSV());
  tag*=LLikelihood(i.Q2_CSV());
  tag*=BLikelihood(i.BHad_CSV());
  tag*=BLikelihood(i.BLep_CSV());
  tag*=BLikelihood(i.B1_CSV());
  tag*=BLikelihood(i.B2_CSV());
  return tag;
  
}
float ReconstructionQuality2::TTWHLikelihoodTimesCSVTimesSCME(Interpretation& i){
  float tag=TTWHLikelihoodTimesSCME(i);
  tag*=LLikelihood(i.Q1_CSV());
  tag*=LLikelihood(i.Q2_CSV());
  tag*=BLikelihood(i.BHad_CSV());
  tag*=BLikelihood(i.BLep_CSV());
  tag*=BLikelihood(i.B1_CSV());
  tag*=BLikelihood(i.B2_CSV());
  return tag;

}


float ReconstructionQuality2::TTWHishLikelihoodTimesSCME(Interpretation& i){
  float llh=0.;
  if(i.HasTag("TTWHishLikelihoodTimesME")) llh= i.GetTag("TTWHishLikelihoodTimesME");
  else llh=TTWHishLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  float me=TTX0_SC_ME(i);
  i.SetTag("TTWHishLikelihoodTimesME",llh*me);
  return llh*me;
}

float ReconstructionQuality2::TTWHishLikelihoodTimesPSME(Interpretation& i){
  float llh=0.;
  if(i.HasTag("TTWHishLikelihoodTimesME")) llh= i.GetTag("TTWHishLikelihoodTimesME");
  else llh=TTWHishLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  float me=TTX0_PS_ME(i);
  i.SetTag("TTWHishLikelihoodTimesME",llh*me);
  return llh*me;
}


float ReconstructionQuality2::TTWHishLikelihood(float mthad, float mtlep, float mwhad, float mhiggs){
  float llh=1;
  llh*=TopHadishLikelihood(mthad);
  llh*=TopLepishLikelihood(mtlep);
  llh*=WHadishLikelihood(mwhad);
  llh*=HiggsishLikelihood(mhiggs);
  return llh;
}


float ReconstructionQuality2::TTWHishLikelihood(Interpretation& i){
  if(i.HasTag("TTWHishLikelihood")) return i.GetTag("TTWHishLikelihood");
  float tag=TTWHishLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWHishLikelihood",tag);
  return tag;
}

float ReconstructionQuality2::TTX0_PS_ME(Interpretation& i){
  if(i.HasTag("TTX0_PS_ME")) return i.GetTag("TTX0_PS_ME");
  float tag=2.5*me.GetTTX0MEsq(i.TopHad(),i.TopLep(),i.Higgs(),true);
  i.SetTag("TTX0_PS_ME",tag);
  return tag;
}

float ReconstructionQuality2::TTX0_SC_ME(Interpretation& i){
  if(i.HasTag("TTX0_SC_ME")) return i.GetTag("TTX0_SC_ME");
  float tag=me.GetTTX0MEsq(i.TopHad(),i.TopLep(),i.Higgs(),false);
  i.SetTag("TTX0_SC_ME",tag);
  return tag;
}

float ReconstructionQuality2::TTX0_SC_ME(const TLorentzVector& t, const TLorentzVector& tbar, const TLorentzVector& H){
  return me.GetTTX0MEsq(t,tbar,H,true);
}
float ReconstructionQuality2::TTX0_PS_ME(const TLorentzVector& t, const TLorentzVector& tbar, const TLorentzVector& H){
  return 2.5*me.GetTTX0MEsq(t,tbar,H,false);
}

float ReconstructionQuality2::BLikelihood(float csv){
  csv=fmax(csv,-0.099);
  float llh=Interpolate(h_CSV_b,csv);
  return llh;
}
float ReconstructionQuality2::LLikelihood(float csv){
  csv=fmax(csv,-0.099);
  float llh=Interpolate(h_CSV_l_w_c, csv);
  return llh;
}
float ReconstructionQuality2::TopHadLikelihood(float m){
  float llh=Interpolate(h_M_TopHad_reco,m);
  return llh;
}
float ReconstructionQuality2::TopHadishLikelihood(float m){
  float llh=Interpolate(h_M_TopHad_best, m);
  return llh;
}
float ReconstructionQuality2::TopLepLikelihood(float m){
  float llh=Interpolate(h_M_TopLep_reco,m);
  return llh;
}
float ReconstructionQuality2::TopLepishLikelihood(float m){
  float llh=Interpolate(h_M_TopLep_best, m);
  return llh;
}
float ReconstructionQuality2::WHadLikelihood(float m){
  float llh= Interpolate(h_M_WHad_reco,m);
  return llh;
}
float ReconstructionQuality2::WHadishLikelihood(float m){
  float llh= Interpolate(h_M_WHad_best,m);
  return llh;
}
float ReconstructionQuality2::HiggsLikelihood(float m){
  float llh= Interpolate(h_M_Higgs_reco,m);
  return llh;
}
float ReconstructionQuality2::HiggsishLikelihood(float m){
  float llh= Interpolate(h_M_Higgs_best,m);
  return llh;
}

float ReconstructionQuality2::Interpolate(TH1F* histo, float value){
  if(value>histo->GetXaxis()->GetXmax()){
    return histo->GetBinContent(histo->GetNbinsX())*tiny_likelihood;
  }
  if(value<histo->GetXaxis()->GetXmin()){
    return histo->GetBinContent(0)*tiny_likelihood;
  }
  return fmax(histo->Interpolate(value),tiny_likelihood);
}

float ReconstructionQuality2::NBLikelihood(uint ntagged, uint njets,const float* csvs){
  if(njets<ntagged) return 0;
  double llh=0.; // total likelihood
  vector<int> idxs(njets);
  for(int i=0;i<njets;i++){
    idxs[i]=i;
  }
  do{
    double l=1; //likelihood of current combination
    for(uint i=0; i<ntagged;i++){
      l*=BLikelihood(csvs[idxs[i]]);
    }
    for(uint i=ntagged; i<njets;i++){
      l*=LLikelihood(csvs[idxs[i]]);
    }
    llh+=l;     
  }while(next_combination(idxs.begin(), idxs.begin()+ntagged, idxs.end()));
  return llh;
}
