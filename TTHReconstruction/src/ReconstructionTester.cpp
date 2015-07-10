#include "../interface/ReconstructionTester.hpp"

using namespace std;

ReconstructionTester::ReconstructionTester(std::vector<std::string> tags_, string outfilename_):tags(tags_),outfilename(outfilename_){
  ntotal=0;

  nmatchAll=0;
  nmatchH=0;
  nmatchWhad=0;
  nmatchBhad=0;
  nmatchBlep=0;
  for(uint t=0; t<tags.size(); t++){
    foundH[tags[t]]=0;
    foundAll[tags[t]]=0;
    InitHisto(tags[t],"Higgs_M_best",40,0,400);
    InitHisto(tags[t],"TopHad_M_best",60,0,600);
    InitHisto(tags[t],"TopLep_M_best",60,0,600);
    InitHisto(tags[t],"WHad_M_best",40,0,400);
  }

  
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

  TFile* outfile=new TFile((outfilename+".root").c_str(),"RECREATE");    
  outfile->cd();
  for(std::map<std::string,std::map<std::string,TH1F*> >::iterator t=allHistos.begin();t!=allHistos.end(); t++){
    for(std::map<std::string,TH1F*>::iterator h=t->second.begin();h!=t->second.end(); h++){
      h->second->Write();
    }
  }
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
    }
    if(best_int[tags[t]]!=0 && mcmatcher.MatchH(*(best_int[tags[t]]))){
      foundH[tags[t]]++;
    }
    PlotInt(tags[t],best_int[tags[t]],"best");
  }
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


void ReconstructionTester::PlotInt(std::string tag, Interpretation* i, std::string suffix){
  FillHisto(tag,"Higgs_M_"+suffix,i->Higgs_M());
  FillHisto(tag,"TopHad_M_"+suffix,i->TopHad_M());
  FillHisto(tag,"TopLep_M_"+suffix,i->TopLep_M());
  FillHisto(tag,"WHad_M_"+suffix,i->WHad_M());

}
