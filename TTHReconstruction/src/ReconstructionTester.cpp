#include "../interface/ReconstructionTester.hpp"

using namespace std;

ReconstructionTester::ReconstructionTester(std::vector<std::string> tags_, string outfilename_, bool testMEs_):
  tags(tags_),outfilename(outfilename_), testMEs(testMEs_), generator(){

  ntotal=0;

  nmatchAll=0;
  nmatchH=0;
  nmatchWhad=0;
  nmatchBhad=0;
  nmatchBlep=0;
  vector<string> tags_plus=tags;
  tags_plus.push_back("mcmatch");
  tags_plus.push_back("perfect");
  vector<string> suffixes;
  suffixes.push_back("best");
  //  suffixes.push_back("random");
  //  suffixes.push_back("second");
  suffixes.push_back("best_and_foundAll");
  suffixes.push_back("best_and_foundH");
  suffixes.push_back("best_and_foundW");
  for(uint t=0; t<tags_plus.size(); t++){
    foundH[tags_plus[t]]=0;
    foundAll[tags_plus[t]]=0;
    for(uint s=0; s<suffixes.size(); s++){
      InitHisto(tags_plus[t],"Higgs_M_"+suffixes[s],40,0,400);
      InitHisto(tags_plus[t],"TopHad_M_"+suffixes[s],60,0,600);
      InitHisto(tags_plus[t],"TopLep_M_"+suffixes[s],60,0,600);
      InitHisto(tags_plus[t],"WHad_M_"+suffixes[s],40,0,400);
      
      InitHisto(tags_plus[t],"Higgs_Eta_"+suffixes[s],60,-3,3);
      InitHisto(tags_plus[t],"TopHad_Eta_"+suffixes[s],60,-3,3);
      InitHisto(tags_plus[t],"TopLep_Eta_"+suffixes[s],60,-3,3);
      InitHisto(tags_plus[t],"WHad_Eta_"+suffixes[s],60,-3,3);
      
      InitHisto(tags_plus[t],"Higgs_Pt_"+suffixes[s],40,0,400);
      InitHisto(tags_plus[t],"TopHad_Pt_"+suffixes[s],60,0,600);
      InitHisto(tags_plus[t],"TopLep_Pt_"+suffixes[s],60,0,600);
      InitHisto(tags_plus[t],"WHad_Pt_"+suffixes[s],40,0,400);
      InitHisto(tags_plus[t],"BHad_CSV_"+suffixes[s],50,0,1);
      InitHisto(tags_plus[t],"BLep_CSV_"+suffixes[s],50,0,1);
      InitHisto(tags_plus[t],"Q1_CSV_"+suffixes[s],50,0,1);
      InitHisto(tags_plus[t],"Q2_CSV_"+suffixes[s],50,0,1);
      InitHisto(tags_plus[t],"B1_CSV_"+suffixes[s],50,0,1);
      InitHisto(tags_plus[t],"B2_CSV_"+suffixes[s],50,0,1);
      
      InitHisto(tags_plus[t],"TTWHChi2_"+suffixes[s],120,0,40);
      InitHisto(tags_plus[t],"TTWChi2_"+suffixes[s],120,0,40);
      InitHisto(tags_plus[t],"TTWBBChi2_"+suffixes[s],120,0,40);
      InitHisto(tags_plus[t],"TTWHChi2_minus_TTWBBChi2_"+suffixes[s],120,-10,30);
      InitHisto(tags_plus[t],"TTWLikelihood_"+suffixes[s],100,-30,-5);
      InitHisto(tags_plus[t],"TTWHLikelihood_"+suffixes[s],100,-30,-5);
      InitHisto(tags_plus[t],"TTWH_TTW_LikelihoodRatio_"+suffixes[s],100,0,0.5);
      
      if(testMEs){
	InitHisto(tags_plus[t],"TTH_ME_"+suffixes[s],100,-18,-8);
	InitHisto(tags_plus[t],"TTHBB_ME_"+suffixes[s],100,-15,-5);
	InitHisto(tags_plus[t],"TTBB_ON_ME_"+suffixes[s],100,-22,-10);
	InitHisto(tags_plus[t],"TTBB_OFF_ME_"+suffixes[s],100,-22,-10);
	InitHisto(tags_plus[t],"TTHBB_TTBB_ON_ME_RATIO_"+suffixes[s],120,-3,16);
	InitHisto(tags_plus[t],"TTHBB_TTBB_OFF_ME_RATIO_"+suffixes[s],120,-6,20);
	
    //    Init2dHisto(tags_plus[t],"TTH_ME_vs_CME_best",20,-20,-6,20,500,2000);
    //    Init2dHisto(tags_plus[t],"TTH_ME_lin_vs_CME_best",40,0.00,0.00004,20,500,2000);
      }
    }
  }

  outfile=new TFile((outfilename+".root").c_str(),"RECREATE");        
  
}
ReconstructionTester::~ReconstructionTester(){
  cout << "match all " << nmatchAll << " / "  << ntotal << " = " << (float)nmatchAll/ntotal << endl;
  cout << "match H " << nmatchH << " / "  << ntotal << " = " << (float)nmatchH/ntotal << endl;
  cout << "match Whad " << nmatchWhad << " / "  << ntotal << " = " << (float)nmatchWhad/ntotal << endl;
  cout << "match Bhad " << nmatchBhad << " / "  << ntotal << " = " << (float)nmatchBhad/ntotal << endl;
  cout << "match Blep " << nmatchBlep << " / "  << ntotal << " = " << (float)nmatchBlep/ntotal << endl;
  for(uint t=0; t<tags.size(); t++){
    cout << "for tag " << tags[t] << ": found H " << foundH[tags[t]] << " / "  << nmatchH << " = " << (float)foundH[tags[t]]/nmatchH << endl;
    cout << "for tag " << tags[t] << ": found All " << foundAll[tags[t]] << " / "  << nmatchAll << " = " << (float)foundAll[tags[t]]/nmatchAll << endl;
  }
  cout << "writing histos" << endl;
  TH1F* h_matches=new TH1F("matches","matches",5,0.5,5.5);
  h_matches->GetXaxis()->SetBinLabel(1,"all");
  h_matches->GetXaxis()->SetBinLabel(2,"higgs");
  h_matches->GetXaxis()->SetBinLabel(3,"whad");
  h_matches->GetXaxis()->SetBinLabel(4,"bhad");
  h_matches->GetXaxis()->SetBinLabel(5,"blep");
  h_matches->SetBinContent(1,nmatchAll);
  h_matches->SetBinContent(2,nmatchH);
  h_matches->SetBinContent(3,nmatchWhad);
  h_matches->SetBinContent(4,nmatchBhad);
  h_matches->SetBinContent(5,nmatchBlep);
  h_matches->Write();
  
  for(std::map<std::string,std::map<std::string,TH1*> >::iterator t=allHistos.begin();t!=allHistos.end(); t++){
    cout << "tag " << t->first << endl;
    for(std::map<std::string,TH1*>::iterator h=t->second.begin();h!=t->second.end(); h++){
      h->second->Write();
    }
  }
  outfile->Close();

  //  h_FoundQuarksChi2TTH->Scale(1./h_FoundQuarksChi2TTH->Integral());
  //  h_FoundQuarksChi2TTH->Write();
  //  h_FoundResonancesChi2TTH->Scale(1./h_FoundResonancesChi2TTH->Integral());
  //  h_FoundResonancesChi2TTH->Write();
}

void ReconstructionTester::Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
				   const TLorentzVector& lepvec, const TVector2& metvec,
				   const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, 
				   const TLorentzVector& vQ2_true, const TLorentzVector& vBLep_true, 
				   const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
				   const TLorentzVector& vB1_true, const TLorentzVector& vB2_true){

  mcmatcher.Setup(vBHad_true,vQ1_true,vQ2_true,vBLep_true,vLep_true,vNu_true,vB1_true,vB2_true);
  vector<Interpretation*> ints = generator.GenerateTTHInterpretations(jetvecs,jetcsvs,lepvec,metvec);
  map<string,float> best_tag;
  map<string,Interpretation*> best_int;
  bool mcmatch=false;
  bool mcmatchH=false;
  bool mcmatchWhad=false;
  bool mcmatchBhad=false;
  bool mcmatchBlep=false;
  Interpretation* best_mcmatch_int=0;
  float best_mcmatch=-1e12;
  for(uint t=0; t<tags.size(); t++){
    best_tag[tags[t]]=-1e12;
    best_int[tags[t]]=0;
  }
  ntotal++;
  for(uint i=0;i<ints.size();i++){
    if(mcmatcher.MatchTTHallQ(*(ints[i]))) mcmatch=true;
    if(mcmatcher.MatchH(*(ints[i]))) mcmatchH=true;
    if(mcmatcher.MatchWHad(*(ints[i]))) mcmatchWhad=true;
    if(mcmatcher.MatchBHad(*(ints[i]))) mcmatchBhad=true;
    if(mcmatcher.MatchBLep(*(ints[i]))) mcmatchBlep=true;
    float sumDr=-mcmatcher.SumDrTTH(*(ints[i]));
    if(sumDr>best_mcmatch){
      best_mcmatch=sumDr;
      best_mcmatch_int=ints[i];
    }
    for(uint t=0; t<tags.size(); t++){
      float tag = quality.GetTag(tags[t],*(ints[i]));
      if(tag>best_tag[tags[t]]){
	best_tag[tags[t]]=tag;
	best_int[tags[t]]=ints[i];
      }
    }
  }
  
  if(mcmatch) nmatchAll++;
  if(mcmatchH) nmatchH++;
  if(mcmatchWhad) nmatchWhad++;
  if(mcmatchBhad) nmatchBhad++;
  if(mcmatchBlep) nmatchBlep++;
  for(uint t=0; t<tags.size(); t++){
    if(best_int[tags[t]]!=0 && mcmatcher.MatchTTHallQ(*(best_int[tags[t]]))){
      foundAll[tags[t]]++;
      PlotInt(tags[t],best_int[tags[t]],"best_and_foundAll");
    }
    if(best_int[tags[t]]!=0 && mcmatcher.MatchH(*(best_int[tags[t]]))){
      foundH[tags[t]]++;
      PlotInt(tags[t],best_int[tags[t]],"best_and_foundH");
    }
    if(best_int[tags[t]]!=0 && mcmatcher.MatchWHad(*(best_int[tags[t]]))){
      PlotInt(tags[t],best_int[tags[t]],"best_and_foundW");
    }
    PlotInt(tags[t],best_int[tags[t]],"best");
  }
  if(best_mcmatch_int!=0){
    PlotInt("mcmatch",best_mcmatch_int,"best");
    if(mcmatcher.MatchWHad(*(best_mcmatch_int))){
      PlotInt("mcmatch",best_mcmatch_int,"best_and_foundW");
    }
    if(mcmatcher.MatchH(*(best_mcmatch_int))){
      PlotInt("mcmatch",best_mcmatch_int,"best_and_foundH");
    }
    if(mcmatcher.MatchTTHallQ(*(best_mcmatch_int))){
      PlotInt("mcmatch",best_mcmatch_int,"best_and_foundAll");
    }

  }

  Interpretation* perfect_int=new Interpretation(vBHad_true,1,vQ1_true,0,vQ2_true,0, vBLep_true,1, 
						 vLep_true, vNu_true, vB1_true,1, vB2_true,1);
  PlotInt("perfect",perfect_int,"best");

  for(uint i=0;i<ints.size();i++){
    delete ints[i];
  }
  ints.clear();

}

void ReconstructionTester::InitHisto(std::string tag,std::string name,int nbins,float xmin, float xmax){
  allHistos[tag][name]=new TH1F((name+"_"+tag).c_str(),(name+"_"+tag).c_str(),nbins,xmin,xmax);
}
void ReconstructionTester::FillHisto(std::string tag,std::string name,float value){
  allHistos[tag][name]->Fill(value);
}

void ReconstructionTester::Init2dHisto(std::string tag,std::string name,int nbinsx,float xmin, float xmax,int nbinsy,float ymin, float ymax){
  allHistos[tag][name]=new TH2F((name+"_"+tag).c_str(),(name+"_"+tag).c_str(),nbinsx,xmin,xmax,nbinsy,ymin,ymax);
}
void ReconstructionTester::Fill2dHisto(std::string tag,std::string name,float xvalue,float yvalue){
  ((TH2*)allHistos[tag][name])->Fill(xvalue,yvalue);
}


void ReconstructionTester::PlotInt(std::string tag, Interpretation* i, std::string suffix){
  FillHisto(tag,"Higgs_M_"+suffix,i->Higgs_M());
  FillHisto(tag,"TopHad_M_"+suffix,i->TopHad_M());
  FillHisto(tag,"TopLep_M_"+suffix,i->TopLep_M());
  FillHisto(tag,"WHad_M_"+suffix,i->WHad_M());
  
  FillHisto(tag,"Higgs_Pt_"+suffix,i->Higgs().Pt());
  FillHisto(tag,"TopHad_Pt_"+suffix,i->TopHad().Pt());
  FillHisto(tag,"TopLep_Pt_"+suffix,i->TopLep().Pt());
  FillHisto(tag,"WHad_Pt_"+suffix,i->WHad().Pt());

  FillHisto(tag,"Higgs_Eta_"+suffix,i->Higgs().Eta());
  FillHisto(tag,"TopHad_Eta_"+suffix,i->TopHad().Eta());
  FillHisto(tag,"TopLep_Eta_"+suffix,i->TopLep().Eta());
  FillHisto(tag,"WHad_Eta_"+suffix,i->WHad().Eta());

  FillHisto(tag,"BHad_CSV_"+suffix,i->BHad_CSV());
  FillHisto(tag,"BLep_CSV_"+suffix,i->BLep_CSV());
  FillHisto(tag,"Q1_CSV_"+suffix,i->Q1_CSV());
  FillHisto(tag,"Q2_CSV_"+suffix,i->Q2_CSV());
  FillHisto(tag,"B1_CSV_"+suffix,i->B1_CSV());
  FillHisto(tag,"B2_CSV_"+suffix,i->B2_CSV());

  FillHisto(tag,"TTWHChi2_"+suffix,-quality.TTWHChi2(*i));
  FillHisto(tag,"TTWHChi2_minus_TTWBBChi2_"+suffix,-quality.TTWHChi2(*i)+quality.TTWBBChi2(*i));
  FillHisto(tag,"TTWBBChi2_"+suffix,-quality.TTWHChi2(*i));
  FillHisto(tag,"TTWChi2_"+suffix,-quality.TTWChi2(*i));
  FillHisto(tag,"TTWHLikelihood_"+suffix,log(quality.TTWHLikelihood(*i)));
  FillHisto(tag,"TTWLikelihood_"+suffix,log(quality.TTWLikelihood(*i)));
  FillHisto(tag,"TTWH_TTW_LikelihoodRatio_"+suffix,quality.TTWHLikelihood(*i)/quality.TTWLikelihood(*i));


  if(testMEs){  
    FillHisto(tag,"TTH_ME_"+suffix,log(quality.TTH_ME(*i)));
    FillHisto(tag,"TTHBB_ME_"+suffix,log(quality.TTHBB_ME(*i)));
    FillHisto(tag,"TTBB_ON_ME_"+suffix,log(quality.TTBB_ON_ME(*i)));
    FillHisto(tag,"TTBB_OFF_ME_"+suffix,log(quality.TTBB_OFF_ME(*i)));
    FillHisto(tag,"TTHBB_TTBB_ON_ME_RATIO_"+suffix,log(quality.TTHBB_TTBB_ON_ME_RATIO(*i)));
    FillHisto(tag,"TTHBB_TTBB_OFF_ME_RATIO_"+suffix,log(quality.TTHBB_TTBB_OFF_ME_RATIO(*i)));
  }

  //  Fill2dHisto(tag,"TTH_ME_vs_CME_"+suffix,log(quality.TTH_ME(*i)),i->Total().M());
  //  Fill2dHisto(tag,"TTH_ME_lin_vs_CME_"+suffix,quality.TTH_ME(*i),i->Total().M());
   
}
