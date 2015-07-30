#include "../interface/ReconstructionQuality.hpp"

using namespace std;

ReconstructionQuality::ReconstructionQuality(string filename){
  file=TFile::Open(filename.c_str());
  h_CSV_b=(TH1F*)file->Get("CSV_b");
  h_CSV_l_w_c=(TH1F*)file->Get("CSV_l_w_c");
  h_M_Higgs_reco=(TH1F*)file->Get("M_Higgs_reco");
  h_M_BB_reco=(TH1F*)file->Get("M_BB_reco");
  h_M_TopHad_reco=(TH1F*)file->Get("M_TopHad_reco");
  h_M_TopLep_reco=(TH1F*)file->Get("M_TopLep_reco");
  h_M_WHad_reco=(TH1F*)file->Get("M_WHad_reco");
  h_M_Higgs_best=(TH1F*)file->Get("M_Higgs_best");
  h_M_TopHad_best=(TH1F*)file->Get("M_TopHad_best");
  h_M_TopLep_best=(TH1F*)file->Get("M_TopLep_best");
  h_M_WHad_best=(TH1F*)file->Get("M_WHad_best");
  h_M_Higgs_all=(TH1F*)file->Get("M_Higgs_all");
  h_M_BB_all=(TH1F*)file->Get("M_BB_all");
  h_M_TopHad_all=(TH1F*)file->Get("M_TopHad_all");
  h_M_TopLep_all=(TH1F*)file->Get("M_TopLep_all");
  h_M_WHad_all=(TH1F*)file->Get("M_WHad_all");

  higgs_mean=111;
  tophad_mean=165;
  toplep_mean=168;
  whad_mean=80;
  higgs_sigma=17;
  tophad_sigma=17;
  toplep_sigma=26;
  whad_sigma=10;
  btagcut=0.89;
  btagbonus=1000.;
  tiny_likelihood=1e-4;
  bb_slope=0.0128;
}

float ReconstructionQuality::GetTag(std::string tag, Interpretation& i){
   
  if(tag=="TTHChi2") return TTHChi2(i);
  if(tag=="TTBBChi2") return TTBBChi2(i);
  else if(tag=="TTChi2") return TTChi2(i);
  else if(tag=="TTWHChi2") return TTWHChi2(i);
  else if(tag=="TTWBBChi2") return TTWBBChi2(i);
  else if(tag=="TTWChi2") return TTWChi2(i);
  else if(tag=="TTHChi2_tagged") return TTHChi2_tagged(i);
  else if(tag=="TTBBChi2_tagged") return TTBBChi2_tagged(i);
  else if(tag=="TTChi2_tagged") return TTChi2_tagged(i);
  else if(tag=="TTWHChi2_tagged") return TTWHChi2_tagged(i);
  else if(tag=="TTWBBChi2_tagged") return TTWBBChi2_tagged(i);
  else if(tag=="TTWChi2_tagged") return TTWChi2_tagged(i);
  else if(tag=="TTChi2_tagged_higgspt") return TTChi2_tagged_higgspt(i);
  else if(tag=="TTChi2_tagged_higgsjetpt") return TTChi2_tagged_higgsjetpt(i);
  else if(tag=="TTWChi2_tagged_higgspt") return TTWChi2_tagged_higgspt(i);
  else if(tag=="TTWChi2_tagged_higgsjetpt") return TTWChi2_tagged_higgsjetpt(i);

  else if(tag=="TTWHLikelihood") return TTWHLikelihood(i);
  else if(tag=="TTWBBLikelihood") return TTWBBLikelihood(i);
  else if(tag=="TTWLikelihood") return TTWLikelihood(i);

  else if(tag=="TTWHLikelihood_tagged") return TTWHLikelihood_tagged(i);
  else if(tag=="TTWBBLikelihood_tagged") return TTWBBLikelihood_tagged(i);
  else if(tag=="TTWLikelihood_tagged") return TTWLikelihood_tagged(i);

  else if(tag=="TTWHLikelihood_comb") return TTWHLikelihood_comb(i);
  else if(tag=="TTWBBLikelihood_comb") return TTWBBLikelihood_comb(i);
  else if(tag=="TTWLikelihood_comb") return TTWLikelihood_comb(i);

  else if(tag=="TTWHLikelihood_comb_ratio") return TTWHLikelihood_comb_ratio(i);
  else if(tag=="TTWBBLikelihood_comb_ratio") return TTWBBLikelihood_comb_ratio(i);
  else if(tag=="TTWLikelihood_comb_ratio") return TTWLikelihood_comb_ratio(i);

  else if(tag=="TTWHLikelihood_comb_ratio_tagged") return TTWHLikelihood_comb_ratio_tagged(i);
  else if(tag=="TTWBBLikelihood_comb_ratio_tagged") return TTWBBLikelihood_comb_ratio_tagged(i);
  else if(tag=="TTWLikelihood_comb_ratio_tagged") return TTWLikelihood_comb_ratio_tagged(i);

  else if(tag=="H_BB_Likelihoodratio") return H_BB_Likelihoodratio(i);

  else if(tag=="TTWHishLikelihood") return TTWHishLikelihood(i);
  else if(tag=="TTWishLikelihood") return TTWishLikelihood(i);
  else if(tag=="TTWHishLikelihood_tagged") return TTWHishLikelihood_tagged(i);
  else if(tag=="TTWishLikelihood_tagged") return TTWishLikelihood_tagged(i);
  else if(tag=="TTH_ME") return TTH_ME(i);
  else if(tag=="TTHBB_ME") return TTHBB_ME(i);
  else if(tag=="TTBB_ON_ME") return TTBB_ON_ME(i);
  else if(tag=="TTBB_OFF_ME") return TTBB_OFF_ME(i);
  else if(tag=="TTH_TTBB_ON_ME_RATIO") return TTH_TTBB_ON_ME_RATIO(i);
  else if(tag=="TTHBB_TTBB_ON_ME_RATIO") return TTHBB_TTBB_ON_ME_RATIO(i);
  else if(tag=="TTH_TTBB_OFF_ME_RATIO") return TTH_TTBB_OFF_ME_RATIO(i);
  else if(tag=="TTHBB_TTBB_OFF_ME_RATIO") return TTHBB_TTBB_OFF_ME_RATIO(i);
  else{
    cout << "tag " << tag << " unknown!" << endl;
    return -1e20;
  }
}

float ReconstructionQuality::TTHChi2(float mthad, float mtlep, float mhiggs){
  float chi2=0;
  chi2+=(mthad-tophad_mean)*(mthad-tophad_mean)/(tophad_sigma*tophad_sigma);
  chi2+=(mtlep-toplep_mean)*(mtlep-toplep_mean)/(toplep_sigma*toplep_sigma);
  chi2+=(mhiggs-higgs_mean)*(mhiggs-higgs_mean)/(higgs_sigma*higgs_sigma);
  return -chi2;
}

float ReconstructionQuality::TTBBChi2(float mthad, float mtlep, float mhiggs){
  float chi2=0;
  chi2+=(mthad-tophad_mean)*(mthad-tophad_mean)/(tophad_sigma*tophad_sigma);
  chi2+=(mtlep-toplep_mean)*(mtlep-toplep_mean)/(toplep_sigma*toplep_sigma);
  chi2+=2*bb_slope*mhiggs;
  return -chi2;
}

float ReconstructionQuality::TTChi2(float mthad, float mtlep){
  float chi2=0;
  chi2+=(mthad-tophad_mean)*(mthad-tophad_mean)/(tophad_sigma*tophad_sigma);
  chi2+=(mtlep-toplep_mean)*(mtlep-toplep_mean)/(toplep_sigma*toplep_sigma);
  return -chi2;
}

float ReconstructionQuality::TTWHChi2(float mthad, float mtlep, float mwhad, float mhiggs){
  float chi2=0;
  chi2+=(mthad-tophad_mean)*(mthad-tophad_mean)/(tophad_sigma*tophad_sigma);
  chi2+=(mtlep-toplep_mean)*(mtlep-toplep_mean)/(toplep_sigma*toplep_sigma);
  chi2+=(mhiggs-higgs_mean)*(mhiggs-higgs_mean)/(higgs_sigma*higgs_sigma);
  chi2+=(mwhad-whad_mean)*(mwhad-whad_mean)/(whad_sigma*whad_sigma);
  return -chi2;
}

float ReconstructionQuality::TTWBBChi2(float mthad, float mtlep, float mwhad, float mhiggs){
  float chi2=0;
  chi2+=(mthad-tophad_mean)*(mthad-tophad_mean)/(tophad_sigma*tophad_sigma);
  chi2+=(mtlep-toplep_mean)*(mtlep-toplep_mean)/(toplep_sigma*toplep_sigma);
  chi2+=2*bb_slope*mhiggs;
  chi2+=(mwhad-whad_mean)*(mwhad-whad_mean)/(whad_sigma*whad_sigma);
  return -chi2;
}
float ReconstructionQuality::TTWChi2(float mthad, float mtlep, float mwhad){
  float chi2=0;
  chi2+=(mthad-tophad_mean)*(mthad-tophad_mean)/(tophad_sigma*tophad_sigma);
  chi2+=(mtlep-toplep_mean)*(mtlep-toplep_mean)/(toplep_sigma*toplep_sigma);
  chi2+=(mwhad-whad_mean)*(mwhad-whad_mean)/(whad_sigma*whad_sigma);
  return -chi2;
}

float ReconstructionQuality::TTWHLikelihood(float mthad, float mtlep, float mwhad, float mhiggs){
  float llh=1;
  llh*=TopHadLikelihood(mthad);
  llh*=TopLepLikelihood(mtlep);
  llh*=WHadLikelihood(mwhad);
  llh*=HiggsLikelihood(mhiggs);
  return llh;
}
float ReconstructionQuality::TTWBBLikelihood(float mthad, float mtlep, float mwhad, float mhiggs){
  float llh=1;
  llh*=TopHadLikelihood(mthad);
  llh*=TopLepLikelihood(mtlep);
  llh*=WHadLikelihood(mwhad);
  llh*=BBLikelihood(mhiggs);
  return llh;
}
float ReconstructionQuality::TTWLikelihood(float mthad, float mtlep, float mwhad){
  float llh=1;
  llh*=TopHadLikelihood(mthad);
  llh*=TopLepLikelihood(mtlep);
  llh*=WHadLikelihood(mwhad);
  return llh;
}
float ReconstructionQuality::TTWHLikelihood(Interpretation& i){
  float tag=TTWHLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWHLikelihood",tag);
  return tag;
}
float ReconstructionQuality::TTWBBLikelihood(Interpretation& i){
  float tag=TTWBBLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWBBLikelihood",tag);
  return tag;
}

float ReconstructionQuality::TTWLikelihood(Interpretation& i){
  float tag=TTWLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M());
  i.SetTag("TTWLikelihood",tag);
  return tag;
}


float ReconstructionQuality::TTWHLikelihood_comb(float mthad, float mtlep, float mwhad, float mhiggs){
  float llh=1;
  llh*=TopHadLikelihood_comb(mthad);
  llh*=TopLepLikelihood_comb(mtlep);
  llh*=WHadLikelihood_comb(mwhad);
  llh*=HiggsLikelihood_comb(mhiggs);
  return llh;
}
float ReconstructionQuality::TTWBBLikelihood_comb(float mthad, float mtlep, float mwhad, float mhiggs){
  float llh=1;
  llh*=TopHadLikelihood_comb(mthad);
  llh*=TopLepLikelihood_comb(mtlep);
  llh*=WHadLikelihood_comb(mwhad);
  llh*=BBLikelihood_comb(mhiggs);
  return llh;
}
float ReconstructionQuality::TTWLikelihood_comb(float mthad, float mtlep, float mwhad){
  float llh=1;
  llh*=TopHadLikelihood_comb(mthad);
  llh*=TopLepLikelihood_comb(mtlep);
  llh*=WHadLikelihood_comb(mwhad);
  return llh;
}
float ReconstructionQuality::TTWHLikelihood_comb(Interpretation& i){
  float tag=TTWHLikelihood_comb(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWHLikelihood_comb",tag);
  return tag;
}
float ReconstructionQuality::TTWBBLikelihood_comb(Interpretation& i){
  float tag=TTWBBLikelihood_comb(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWBBLikelihood_comb",tag);
  return tag;
}

float ReconstructionQuality::TTWLikelihood_comb(Interpretation& i){
  float tag=TTWLikelihood_comb(i.TopHad_M(),i.TopLep_M(),i.WHad_M());
  i.SetTag("TTWLikelihood_comb",tag);
  return tag;
}

float ReconstructionQuality::TTWHLikelihood_comb_ratio(float mthad, float mtlep, float mwhad, float mhiggs){
  double n=1;
  double d=1;
  n*=TopHadLikelihood(mthad,false);
  d*=TopHadLikelihood_comb(mthad,false);
  n*=TopLepLikelihood(mtlep,false);
  d*=TopLepLikelihood_comb(mtlep,false);
  n*=WHadLikelihood(mwhad,false);
  d*=WHadLikelihood_comb(mwhad,false);
  n*=HiggsLikelihood(mhiggs,false);
  d*=HiggsLikelihood_comb(mhiggs,false);
  return n/(n+d);
}
float ReconstructionQuality::TTWBBLikelihood_comb_ratio(float mthad, float mtlep, float mwhad, float mhiggs){
  double n=1;
  double d=1;
  n*=TopHadLikelihood(mthad,false);
  d*=TopHadLikelihood_comb(mthad,false);
  n*=TopLepLikelihood(mtlep,false);
  d*=TopLepLikelihood_comb(mtlep,false);
  n*=WHadLikelihood(mwhad,false);
  d*=WHadLikelihood_comb(mwhad,false);
  n*=BBLikelihood(mhiggs,false);
  d*=BBLikelihood_comb(mhiggs,false);
  return n/(n+d);
}
float ReconstructionQuality::TTWLikelihood_comb_ratio(float mthad, float mtlep, float mwhad){
  double n=1;
  double d=1;
  n*=TopHadLikelihood(mthad,false);
  d*=TopHadLikelihood_comb(mthad,false);
  n*=TopLepLikelihood(mtlep,false);
  d*=TopLepLikelihood_comb(mtlep,false);
  n*=WHadLikelihood(mwhad,false);
  d*=WHadLikelihood_comb(mwhad,false);
  return n/(n+d);
}
float ReconstructionQuality::TTWHLikelihood_comb_ratio(Interpretation& i){
  float tag=TTWHLikelihood_comb_ratio(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWHLikelihood_comb_ratio",tag);
  return tag;
}
float ReconstructionQuality::TTWBBLikelihood_comb_ratio(Interpretation& i){
  float tag=TTWBBLikelihood_comb_ratio(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWBBLikelihood_comb_ratio",tag);
  return tag;
}

float ReconstructionQuality::TTWLikelihood_comb_ratio(Interpretation& i){
  float tag=TTWLikelihood_comb_ratio(i.TopHad_M(),i.TopLep_M(),i.WHad_M());
  i.SetTag("TTWLikelihood_comb_ratio",tag);
  return tag;
}


float ReconstructionQuality::TTWHLikelihood_comb_ratio_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float llh=TTWHLikelihood_comb_ratio(mthad, mtlep, mwhad, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return llh*pow(btagbonus,ntags);
}
float ReconstructionQuality::TTWBBLikelihood_comb_ratio_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float llh=TTWBBLikelihood_comb_ratio(mthad, mtlep, mwhad, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return llh*pow(btagbonus,ntags);
}
float ReconstructionQuality::TTWLikelihood_comb_ratio_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1, float tagb2){
  float llh=TTWLikelihood_comb_ratio(mthad, mtlep, mwhad);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return llh*pow(btagbonus,ntags);
}
float ReconstructionQuality::TTWHLikelihood_comb_ratio_tagged(Interpretation& i){
  float tag=TTWHLikelihood_comb_ratio_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTWHLikelihood_comb_ratio_tagged",tag);
  return tag;
}
float ReconstructionQuality::TTWBBLikelihood_comb_ratio_tagged(Interpretation& i){
  float tag=TTWBBLikelihood_comb_ratio_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTWBBLikelihood_comb_ratio_tagged",tag);
  return tag;
}
float ReconstructionQuality::TTWLikelihood_comb_ratio_tagged(Interpretation& i, bool inclHiggsTags){
  float tag;
  if(inclHiggsTags) tag=TTWLikelihood_comb_ratio_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  else tag=TTWLikelihood_comb_ratio_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV());
  i.SetTag("TTWLikelihood_comb_ratio_tagged",tag);
  return tag;
}

float ReconstructionQuality::TTWHLikelihood_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float llh=TTWHLikelihood(mthad, mtlep, mwhad, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return llh*pow(btagbonus,ntags);
}
float ReconstructionQuality::TTWBBLikelihood_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float llh=TTWBBLikelihood(mthad, mtlep, mwhad, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return llh*pow(btagbonus,ntags);
}

float ReconstructionQuality::TTWLikelihood_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1, float tagb2){
  float llh=TTWLikelihood(mthad, mtlep, mwhad);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return llh*pow(btagbonus,ntags);
}
float ReconstructionQuality::TTWHLikelihood_tagged(Interpretation& i){
  float tag=TTWHLikelihood_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTWHLikelihood_tagged",tag);
  return tag;
}
float ReconstructionQuality::TTWBBLikelihood_tagged(Interpretation& i){
  float tag=TTWBBLikelihood_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTWBBLikelihood_tagged",tag);
  return tag;
}
float ReconstructionQuality::TTWLikelihood_tagged(Interpretation& i, bool inclHiggsTags){
  float tag;
  if(inclHiggsTags) tag=TTWLikelihood_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  else tag=TTWLikelihood_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV());
  i.SetTag("TTWLikelihood_tagged",tag);
  return tag;
}

float ReconstructionQuality::H_BB_Likelihoodratio(Interpretation& i){
  float mhiggs=i.Higgs_M();
  float tag= 1/(1+BBLikelihood(mhiggs,false)/HiggsLikelihood(mhiggs,false));
  i.SetTag("H_BB_Likelihoodratio",tag);
  return tag;
}


float ReconstructionQuality::TTWHishLikelihood(float mthad, float mtlep, float mwhad, float mhiggs){
  float llh=1;
  llh*=TopHadishLikelihood(mthad);
  llh*=TopLepishLikelihood(mtlep);
  llh*=WHadishLikelihood(mwhad);
  llh*=HiggsishLikelihood(mhiggs);
  return llh;
}

float ReconstructionQuality::TTWishLikelihood(float mthad, float mtlep, float mwhad){
  float llh=1;
  llh*=TopHadishLikelihood(mthad);
  llh*=TopLepishLikelihood(mtlep);
  llh*=WHadishLikelihood(mwhad);
  return llh;
}

float ReconstructionQuality::TTWHishLikelihood_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float llh=TTWHishLikelihood(mthad, mtlep, mwhad, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return llh*pow(btagbonus,ntags);
}

float ReconstructionQuality::TTWishLikelihood_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1, float tagb2){
  float llh=TTWishLikelihood(mthad, mtlep, mwhad);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return llh*pow(btagbonus,ntags);
}


float ReconstructionQuality::TTHChi2_tagged(float mthad, float mtlep, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTHChi2(mthad, mtlep, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return chi2*pow(1./btagbonus,ntags);
}

float ReconstructionQuality::TTBBChi2_tagged(float mthad, float mtlep, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTBBChi2(mthad, mtlep, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return chi2*pow(1./btagbonus,ntags);
}

float ReconstructionQuality::TTChi2_tagged(float mthad, float mtlep, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTChi2(mthad, mtlep);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return chi2*pow(1./btagbonus,ntags);
}


float ReconstructionQuality::TTWHChi2_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTWHChi2(mthad, mtlep, mwhad, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return chi2*pow(1./btagbonus,ntags);
}
float ReconstructionQuality::TTWBBChi2_tagged(float mthad, float mtlep, float mwhad, float mhiggs, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTWBBChi2(mthad, mtlep, mwhad, mhiggs);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return chi2*pow(1./btagbonus,ntags);
}


float ReconstructionQuality::TTWChi2_tagged(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1, float tagb2){
  float chi2=TTWChi2(mthad, mtlep, mwhad);
  int ntags=0;
  if(tagblep>btagcut) ntags++;
  if(tagbhad>btagcut) ntags++;
  if(tagb1>btagcut) ntags++;
  if(tagb2>btagcut) ntags++;
  return chi2*pow(1./btagbonus,ntags);
}

float ReconstructionQuality::TTChi2_tagged_higgspt(float mthad, float mtlep, float tagbhad, float tagblep, float tagb1, float tagb2,float pt ){
  return log(1+pt/100)*TTChi2_tagged(mthad,mtlep,tagbhad,tagblep,tagb1,tagb2);
}
float ReconstructionQuality::TTWChi2_tagged_higgspt(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1, float tagb2,float pt ){
  return log(1+pt/100)*TTWChi2_tagged(mthad,mtlep,mwhad,tagbhad,tagblep,tagb1,tagb2);
}
float ReconstructionQuality::TTChi2_tagged_higgsjetpt(float mthad, float mtlep, float tagbhad, float tagblep, float tagb1, float tagb2, float pt1, float pt2){
  return log(1+pt1*pt2/10000)*TTWChi2_tagged(mthad,mtlep,tagbhad,tagblep,tagb1,tagb2);
}
float ReconstructionQuality::TTWChi2_tagged_higgsjetpt(float mthad, float mtlep, float mwhad, float tagbhad, float tagblep, float tagb1, float tagb2, float pt1, float pt2){
  return log(1+pt1*pt2/10000)*TTWChi2_tagged(mthad,mtlep,mwhad,tagbhad,tagblep,tagb1,tagb2);
}


float ReconstructionQuality::TTChi2_tagged_higgspt(Interpretation& i){
  float tag=TTChi2_tagged_higgspt(i.TopHad_M(),i.TopLep_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV(),i.Higgs().Pt());
  i.SetTag("TTChi2_tagged_higgspt",tag);
  return tag;
}
float ReconstructionQuality::TTWChi2_tagged_higgspt(Interpretation& i){
  float tag=TTWChi2_tagged_higgspt(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV(),i.Higgs().Pt());
  i.SetTag("TTWChi2_tagged_higgspt",tag);
  return tag;

}
float ReconstructionQuality::TTChi2_tagged_higgsjetpt(Interpretation& i){
  float tag=TTChi2_tagged_higgsjetpt(i.TopHad_M(),i.TopLep_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV(),i.B1().Pt(),i.B2().Pt());
  i.SetTag("TTChi2_tagged_higgsjetpt",tag);
  return tag;
}
float ReconstructionQuality::TTWChi2_tagged_higgsjetpt(Interpretation& i){
  float tag=TTWChi2_tagged_higgsjetpt(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV(),i.B1().Pt(),i.B2().Pt());
  i.SetTag("TTWChi2_tagged_higgsjetpt",tag);
  return tag;
}

float ReconstructionQuality::TTHChi2_tagged(Interpretation& i){
  float tag=TTHChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTHChi2_tagged",tag);
  return tag;
}
float ReconstructionQuality::TTBBChi2_tagged(Interpretation& i){
  float tag=TTBBChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTBBChi2_tagged",tag);
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
float ReconstructionQuality::TTWBBChi2_tagged(Interpretation& i){
  float tag=TTWBBChi2_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTWBBChi2_tagged",tag);
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
float ReconstructionQuality::TTBBChi2(Interpretation& i){
  float tag=TTBBChi2(i.TopHad_M(),i.TopLep_M(),i.Higgs_M());
  i.SetTag("TTBBChi2",tag);
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
float ReconstructionQuality::TTWBBChi2(Interpretation& i){
  float tag=TTWBBChi2(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWBBChi2",tag);
  return tag;
}

float ReconstructionQuality::TTWChi2(Interpretation& i){
  float tag=TTWChi2(i.TopHad_M(),i.TopLep_M(),i.WHad_M());
  i.SetTag("TTWChi2",tag);
  return tag;
}


float ReconstructionQuality::TTWHishLikelihood(Interpretation& i){
  float tag=TTWHishLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M());
  i.SetTag("TTWHishLikelihood",tag);
  return tag;
}
float ReconstructionQuality::TTWishLikelihood(Interpretation& i){
  float tag=TTWishLikelihood(i.TopHad_M(),i.TopLep_M(),i.WHad_M());
  i.SetTag("TTWishLikelihood",tag);
  return tag;
}
float ReconstructionQuality::TTWHishLikelihood_tagged(Interpretation& i){
  float tag=TTWHishLikelihood_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.Higgs_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  i.SetTag("TTWHishLikelihood_tagged",tag);
  return tag;
}
float ReconstructionQuality::TTWishLikelihood_tagged(Interpretation& i, bool inclHiggsTags){
  float tag;
  if(inclHiggsTags) tag=TTWishLikelihood_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV(),i.B1_CSV(),i.B2_CSV());
  else tag=TTWishLikelihood_tagged(i.TopHad_M(),i.TopLep_M(),i.WHad_M(),i.BHad_CSV(),i.BLep_CSV());
  i.SetTag("TTWishLikelihood_tagged",tag);
  return tag;
}

float ReconstructionQuality::TTH_ME(Interpretation& i){
  float tag=me.GetTTHMEsq(i.TopHad(),i.TopLep(),i.Higgs());
  i.SetTag("TTH_ME",tag);
  return tag;
}

float ReconstructionQuality::TTBB_ON_ME(Interpretation& i){
  float tag=me.GetTTBBMEsq_onshell(i.TopHad(),i.TopLep(),i.B1(),i.B2());
  i.SetTag("TTBB_ON_ME",tag);
  return tag;
}
float ReconstructionQuality::TTBB_OFF_ME(Interpretation& i){
  float tag=me.GetTTBBMEsq_offshell(i.TopHad(),i.TopLep(),i.B1(),i.B2());
  i.SetTag("TTBB_OFF_ME",tag);
  return tag;
}


float ReconstructionQuality::TTHBB_ME(Interpretation& i){
  float tag=me.GetTTHBBMEsq(i.TopHad(),i.TopLep(),i.B1(),i.B2());
  i.SetTag("TTHBB_ME",tag);
  return tag;
}

float ReconstructionQuality::TTH_TTBB_ON_ME_RATIO(Interpretation& i){
  float tag=1./(1+me.GetTTBBMEsq_onshell(i.TopHad(),i.TopLep(),i.B1(),i.B2())/me.GetTTHMEsq(i.TopHad(),i.TopLep(),i.Higgs()));
  i.SetTag("TTH_TTBB_ON_ME_RATIO",tag);
  return tag;
}
float ReconstructionQuality::TTHBB_TTBB_ON_ME_RATIO(Interpretation& i){
  float tag=1./(1+1e4*me.GetTTBBMEsq_onshell(i.TopHad(),i.TopLep(),i.B1(),i.B2())/me.GetTTHBBMEsq(i.TopHad(),i.TopLep(),i.B1(),i.B2()));
  i.SetTag("TTHBB_TTBB_ON_ME_RATIO",tag);
  return tag;
}
float ReconstructionQuality::TTH_TTBB_OFF_ME_RATIO(Interpretation& i){
  float tag=1./(1+me.GetTTBBMEsq_offshell(i.TopHad(),i.TopLep(),i.B1(),i.B2())/me.GetTTHMEsq(i.TopHad(),i.TopLep(),i.Higgs()));
  i.SetTag("TTH_TTBB_OFF_ME_RATIO",tag);
  return tag;
}
float ReconstructionQuality::TTHBB_TTBB_OFF_ME_RATIO(Interpretation& i){
  float tag=1./(1+1e4*me.GetTTBBMEsq_offshell(i.TopHad(),i.TopLep(),i.B1(),i.B2())/me.GetTTHBBMEsq(i.TopHad(),i.TopLep(),i.B1(),i.B2()));
  i.SetTag("TTHBB_TTBB_OFF_ME_RATIO",tag);
  return tag;
}


float ReconstructionQuality::BLikelihood(float csv){
  csv=fmax(csv,-0.099);
  float llh=Interpolate(h_CSV_b,csv);
  return llh;
}
float ReconstructionQuality::LLikelihood(float csv){
  csv=fmax(csv,-0.099);
  float llh=Interpolate(h_CSV_l_w_c, csv);
  return llh;
}
float ReconstructionQuality::TopHadLikelihood(float m, bool exclude_overflow){
  float llh=Interpolate(h_M_TopHad_reco,m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::TopHadLikelihood_comb(float m, bool exclude_overflow){
  float llh=Interpolate(h_M_TopHad_all,m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::TopHadishLikelihood(float m, bool exclude_overflow){
  float llh=Interpolate(h_M_TopHad_best, m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::TopLepLikelihood(float m, bool exclude_overflow){
  float llh=Interpolate(h_M_TopLep_reco,m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::TopLepLikelihood_comb(float m, bool exclude_overflow){
  float llh=Interpolate(h_M_TopLep_all,m, exclude_overflow);
  return llh;
}

float ReconstructionQuality::TopLepishLikelihood(float m, bool exclude_overflow){
  float llh=Interpolate(h_M_TopLep_best, m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::WHadLikelihood(float m, bool exclude_overflow){
  float llh= Interpolate(h_M_WHad_reco,m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::WHadLikelihood_comb(float m, bool exclude_overflow){
  float llh= Interpolate(h_M_WHad_all,m, exclude_overflow);
  return llh;
}

float ReconstructionQuality::WHadishLikelihood(float m, bool exclude_overflow){
  float llh= Interpolate(h_M_WHad_best,m, exclude_overflow);
  return llh;
}

float ReconstructionQuality::HiggsLikelihood(float m, bool exclude_overflow){
  float llh= Interpolate(h_M_Higgs_reco,m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::HiggsLikelihood_comb(float m, bool exclude_overflow){
  float llh= Interpolate(h_M_Higgs_all,m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::BBLikelihood(float m, bool exclude_overflow){
  float llh= Interpolate(h_M_BB_reco,m, exclude_overflow);
  return llh;
}
float ReconstructionQuality::BBLikelihood_comb(float m, bool exclude_overflow){
  float llh= Interpolate(h_M_BB_all,m, exclude_overflow);
  return llh;
}

float ReconstructionQuality::HiggsishLikelihood(float m, bool exclude_overflow){
  float llh= Interpolate(h_M_Higgs_best,m, exclude_overflow);
  return llh;
}

float ReconstructionQuality::Interpolate(TH1F* histo, float value, bool exclude_overflow){
  if(value>histo->GetXaxis()->GetXmax()||value<histo->GetXaxis()->GetXmin()){
    //    if(exclude_overflow)
    //      return tiny_likelihood;
    //    else
    return histo->GetBinContent(histo->FindBin(value))*tiny_likelihood;
  }
  return fmax(histo->Interpolate(value),tiny_likelihood);
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
