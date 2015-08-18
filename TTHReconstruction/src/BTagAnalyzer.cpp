#include "../interface/BTagAnalyzer.hpp"

using namespace std;

BTagAnalyzer::BTagAnalyzer(string outfilename):generator(InterpretationGenerator(IntType::tth,0)){
  outfile=new TFile((outfilename+"_ana.root").c_str(),"RECREATE");        
  nevents=0;
  nselected=0;
  h_r32=new TH1F("r32","r32",2000,0,1);
  h_r42=new TH1F("r42","r42",2000,0,1);
  h_r43=new TH1F("r43","r43",2000,0,1);

  h_r32_4b=new TH1F("r32_4b","r32_4b",2000,0,1);
  h_r42_4b=new TH1F("r42_4b","r42_4b",2000,0,1);
  h_r43_4b=new TH1F("r43_4b","r43_4b",2000,0,1);

  h_r32_3b=new TH1F("r32_3b","r32_3b",2000,0,1);
  h_r42_3b=new TH1F("r42_3b","r42_3b",2000,0,1);
  h_r43_3b=new TH1F("r43_3b","r43_3b",2000,0,1);

  h_r32_2b=new TH1F("r32_2b","r32_2b",2000,0,1);
  h_r42_2b=new TH1F("r42_2b","r42_2b",2000,0,1);
  h_r43_2b=new TH1F("r43_2b","r43_2b",2000,0,1);

  h_ntagsT=new TH1F("ntagsT","ntagsT",6,-0.5,5.5);
  h_ntagsM=new TH1F("ntagsM","ntagsM",6,-0.5,5.5);
  h_ntagsL=new TH1F("ntagsL","ntagsL",6,-0.5,5.5);

  h_ntagsT_4b=new TH1F("ntagsT_4b","ntagsT_4b",6,-0.5,5.5);
  h_ntagsM_4b=new TH1F("ntagsM_4b","ntagsM_4b",6,-0.5,5.5);
  h_ntagsL_4b=new TH1F("ntagsL_4b","ntagsL_4b",6,-0.5,5.5);

  h_ntagsT_3b=new TH1F("ntagsT_3b","ntagsT_3b",6,-0.5,5.5);
  h_ntagsM_3b=new TH1F("ntagsM_3b","ntagsM_3b",6,-0.5,5.5);
  h_ntagsL_3b=new TH1F("ntagsL_3b","ntagsL_3b",6,-0.5,5.5);

  h_ntagsT_2b=new TH1F("ntagsT_2b","ntagsT_2b",6,-0.5,5.5);
  h_ntagsM_2b=new TH1F("ntagsM_2b","ntagsM_2b",6,-0.5,5.5);
  h_ntagsL_2b=new TH1F("ntagsL_2b","ntagsL_2b",6,-0.5,5.5);

  h_averageTagged=new TH1F("averageTagged","averageTagged",5000,0,5);
  h_averageTagged_2b=new TH1F("averageTagged_2b","averageTagged_2b",5000,0,5);
  h_averageTagged_3b=new TH1F("averageTagged_3b","averageTagged_3b",5000,0,5);
  h_averageTagged_4b=new TH1F("averageTagged_4b","averageTagged_4b",5000,0,5);

}
BTagAnalyzer::~BTagAnalyzer(){
  cout << "nevents " << nevents << endl;
  cout << "nselected " << nselected << endl;
  
  outfile->cd();
  // write histos
  h_r32->Write();
  h_r42->Write();
  h_r43->Write();

  h_r32_4b->Write();
  h_r42_4b->Write();
  h_r43_4b->Write();

  h_r32_3b->Write();
  h_r42_3b->Write();
  h_r43_3b->Write();

  h_r32_2b->Write();
  h_r42_2b->Write();
  h_r43_2b->Write();

  h_ntagsT->Write();
  h_ntagsM->Write();
  h_ntagsL->Write();

  h_ntagsT_4b->Write();
  h_ntagsM_4b->Write();
  h_ntagsL_4b->Write();

  h_ntagsT_3b->Write();
  h_ntagsM_3b->Write();
  h_ntagsL_3b->Write();

  h_ntagsT_2b->Write();
  h_ntagsM_2b->Write();
  h_ntagsL_2b->Write();

  h_averageTagged->Write();
  h_averageTagged_2b->Write();
  h_averageTagged_3b->Write();
  h_averageTagged_4b->Write();


  outfile->Close();
}

void BTagAnalyzer::Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs,const std::vector<float>& jetflavs){

  nevents++;
  nselected++;

  int nb=0;
  int ntagsT=0;
  int ntagsM=0;
  int ntagsL=0;
  float averageTagged=0;
  for(uint i=0; i< jetflavs.size(); i++){
    uint flav=abs(int(fabs(jetflavs[i])+.5));
    if(flav==5) nb++;
    if(jetcsvs[i]>.605) ntagsL++;
    if(jetcsvs[i]>.89) ntagsM++;
    if(jetcsvs[i]>.97) ntagsT++;
    if(jetcsvs[i]>.89) averageTagged+=jetcsvs[i];
  }
  averageTagged/=ntagsM;
  averageTagged+=min(4,ntagsM);

  float p_4b=quality.NBLikelihood(4,jetcsvs.size(),&(jetcsvs[0]));
  float p_3b=quality.NBLikelihood(3,jetcsvs.size(),&(jetcsvs[0]));
  float p_2b=quality.NBLikelihood(2,jetcsvs.size(),&(jetcsvs[0]));
  float r_42=p_4b/(p_4b+p_2b);
  float r_32=p_3b/(p_3b+p_2b);
  float r_43=p_4b/(p_4b+p_3b);
  h_ntagsT->Fill(ntagsT);
  h_ntagsM->Fill(ntagsM);
  h_ntagsL->Fill(ntagsL);
  h_r32->Fill(r_32);
  h_r42->Fill(r_42);
  h_r43->Fill(r_43);
  h_averageTagged->Fill(averageTagged);

  if(nb==4){
    h_ntagsT_4b->Fill(ntagsT);
    h_ntagsM_4b->Fill(ntagsM);
    h_ntagsL_4b->Fill(ntagsL);
    h_r32_4b->Fill(r_32);
    h_r42_4b->Fill(r_42);
    h_r43_4b->Fill(r_43);
    h_averageTagged_4b->Fill(averageTagged);
  }
  if(nb==3){
    h_ntagsT_3b->Fill(ntagsT);
    h_ntagsM_3b->Fill(ntagsM);
    h_ntagsL_3b->Fill(ntagsL);
    h_r32_3b->Fill(r_32);
    h_r42_3b->Fill(r_42);
    h_r43_3b->Fill(r_43);
    h_averageTagged_3b->Fill(averageTagged);
  }
  if(nb==2){
    h_ntagsT_2b->Fill(ntagsT);
    h_ntagsM_2b->Fill(ntagsM);
    h_ntagsL_2b->Fill(ntagsL);
    h_r32_2b->Fill(r_32);
    h_r42_2b->Fill(r_42);
    h_r43_2b->Fill(r_43);
    h_averageTagged_2b->Fill(averageTagged);
  }
}

