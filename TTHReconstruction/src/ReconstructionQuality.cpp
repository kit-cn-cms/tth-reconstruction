#include "../interface/ReconstructionQuality.hpp"

using namespace std;

ReconstructionQuality::ReconstructionQuality(string filename){
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

float ReconstructionQuality::GetTag(std::string tag, Interpretation& i){
   
  if(tag=="TTHChi2") return TTHChi2(i);
  else if(tag=="TTChi2") return TTChi2(i);
  else if(tag=="TTWHChi2") return TTWHChi2(i);
  else if(tag=="TTWChi2") return TTWChi2(i);
  else if(tag=="TTHChi2_tagged") return TTHChi2_tagged(i);
  else if(tag=="TTChi2_tagged") return TTChi2_tagged(i);
  else if(tag=="TTWHChi2_tagged") return TTWHChi2_tagged(i);
  else if(tag=="TTWChi2_tagged") return TTWChi2_tagged(i);
  else{
    cout << "tag " << tag << " unknown!" << endl;
    return -1e20;
  }
}

float ReconstructionQuality::TTHChi2(float mthad, float mtlep, float mhiggs){
  float chi2=0;
  chi2+=(mthad-171)*(mthad-171)/(21*21);
  chi2+=(mtlep-171)*(mtlep-171)/(27*27);
  chi2+=(mhiggs-115)*(mhiggs-115)/(18*18);
  return -chi2;
}

float ReconstructionQuality::TTChi2(float mthad, float mtlep){
  float chi2=0;
  chi2+=(mthad-171)*(mthad-171)/(21*21);
  chi2+=(mtlep-171)*(mtlep-171)/(27*27);
  return -chi2;
}

float ReconstructionQuality::TTWHChi2(float mthad, float mtlep, float mwhad, float mhiggs){
  float chi2=0;
  chi2+=(mthad-171)*(mthad-171)/(21*21);
  chi2+=(mtlep-171)*(mtlep-171)/(27*27);
  chi2+=(mhiggs-115)*(mhiggs-115)/(18*18);
  chi2+=(mwhad-83)*(mwhad-83)/(13*13);
  return -chi2;
}

float ReconstructionQuality::TTWChi2(float mthad, float mtlep, float mwhad){
  float chi2=0;
  chi2+=(mthad-171)*(mthad-171)/(21*21);
  chi2+=(mtlep-171)*(mtlep-171)/(27*27);
  chi2+=(mwhad-83)*(mwhad-83)/(13*13);
  return -chi2;
}

float ReconstructionQuality::TTHChi2_tagged(float mthad, float mtlep, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTHChi2(mthad, mtlep, mhiggs);
  int ntags=0;
  if(tagblep>0.814) ntags++;
  if(tagbhad>0.814) ntags++;
  if(tagb1>0.814) ntags++;
  if(tagb2>0.814) ntags++;
  return chi2*pow(0.001,ntags);
}

float ReconstructionQuality::TTChi2_tagged(float mthad, float mtlep, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTChi2(mthad, mtlep);
  int ntags=0;
  if(tagblep>0.814) ntags++;
  if(tagbhad>0.814) ntags++;
  if(tagb1>0.814) ntags++;
  if(tagb2>0.814) ntags++;
  return chi2*pow(0.001,ntags);
}


float ReconstructionQuality::TTWHChi2_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTWHChi2(mthad, mtlep, mwhad, mhiggs);
  int ntags=0;
  if(tagblep>0.814) ntags++;
  if(tagbhad>0.814) ntags++;
  if(tagb1>0.814) ntags++;
  if(tagb2>0.814) ntags++;
  return chi2*pow(0.001,ntags);
}

float ReconstructionQuality::TTWChi2_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTWChi2(mthad, mtlep, mwhad);
  int ntags=0;
  if(tagblep>0.814) ntags++;
  if(tagbhad>0.814) ntags++;
  if(tagb1>0.814) ntags++;
  if(tagb2>0.814) ntags++;
  return chi2*pow(0.001,ntags);
}


float ReconstructionQuality::TTHChi2_tagged(Interpretation& i){
  float tag=TTHChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTHChi2_tagged",tag);
  return tag;
}
float ReconstructionQuality::TTChi2_tagged(Interpretation& i, bool inclHiggsTags){
  float tag=-1;
  if(inclHiggsTags) tag=TTChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  else tag=TTChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.BHad_CSV(),i.BLep_CSV());
  i.SetTag("TTChi2_tagged",tag);
  return tag;

}
float ReconstructionQuality::TTWHChi2_tagged(Interpretation& i){
  float tag=TTWHChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTWHChi2_tagged",tag);
  return tag;
}
float ReconstructionQuality::TTWChi2_tagged(Interpretation& i, bool inclHiggsTags){
  float tag;
  if(inclHiggsTags) tag=TTWChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  else tag=TTWChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV());
  i.SetTag("TTWChi2_tagged",tag);
  return tag;

}

float ReconstructionQuality::TTHChi2(Interpretation& i){
  float tag=TTHChi2(i.TopHad_M(),i.TopLep_M(),i.Higgs_M());
  i.SetTag("TTHChi2",tag);
  return tag;
}
float ReconstructionQuality::TTChi2(Interpretation& i){
  float tag=TTChi2(i.TopHad_M(),i.TopLep_M());
  i.SetTag("TTChi2",tag);
  return tag;
}
float ReconstructionQuality::TTWHChi2(Interpretation& i){
  float tag=TTWHChi2(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWHChi2",tag);
  return tag;
}
float ReconstructionQuality::TTWChi2(Interpretation& i){
  float tag=TTWChi2(i.TopHad_M(),i.TopLep_M(),i.WHad_M());
  i.SetTag("TTWChi2",tag);
  return tag;
}
float ReconstructionQuality::BLikelihood(float csv){
  csv=fmax(csv,-0.099);
  float llh= h_CSV_b->Interpolate(csv);
  return llh;
}
float ReconstructionQuality::LLikelihood(float csv){
  csv=fmax(csv,-0.099);
  float llh= h_CSV_l_w_c->Interpolate(csv);
  return llh;
}
float ReconstructionQuality::TopHadLikelihood(float m){
  float llh= h_M_TopHad_reco->Interpolate(m);
  return llh;
}
float ReconstructionQuality::TopHadishLikelihood(float m){
  float llh= h_M_TopHad_best->Interpolate(m);
  return llh;
}
float ReconstructionQuality::TopLepLikelihood(float m){
  float llh= h_M_TopLep_reco->Interpolate(m);
  return llh;
}
float ReconstructionQuality::TopLepishLikelihood(float m){
  float llh= h_M_TopLep_best->Interpolate(m);
  return llh;
}
float ReconstructionQuality::WHadLikelihood(float m){
  float llh= h_M_WHad_reco->Interpolate(m);
  return llh;
}
float ReconstructionQuality::WHadishLikelihood(float m){
  float llh= h_M_WHad_best->Interpolate(m);
  return llh;
}

float ReconstructionQuality::NBLikelihood(uint ntagged, uint njets, float* csvs){
  if(njets<ntagged) return 0;
  double llh=0.; // total likelihood
  uint binary=pow(2,ntagged)-1; //e.g.2^4-1=15=00001111, only last 4 tagged
  uint lastcomb=binary*pow(2,(njets-ntagged)); //e.g. 11110000, only first 4 tagged
  while(binary<=lastcomb){
    double l=1; //likelihood of current combination
    uint b=binary; //used to copute nth digit in binary
    for(uint ijet=0; ijet<njets;ijet++){
      if(b%2==1){
	l*=BLikelihood(csvs[ijet]);
      }
      else{
	l*=LLikelihood(csvs[ijet]);
      }
      b=b/2;
    }
    llh+=l;
    // get next permutation: 00001111 -> 00010111
    uint t = (binary | (binary - 1)) + 1;  
      binary = t | ((((t & -t) / (binary & -binary)) >> 1) - 1);  
      
  }   
  return llh;
}
