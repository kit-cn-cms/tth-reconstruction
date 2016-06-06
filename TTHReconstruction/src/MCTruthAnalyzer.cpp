#include "../interface/MCTruthAnalyzer.hpp"

using namespace std;

MCTruthAnalyzer::MCTruthAnalyzer(string outfilename):generator(InterpretationGenerator(IntType::tth,0)){
  outfile=new TFile((outfilename+"_ana.root").c_str(),"RECREATE");        
  nevents=0;
  nselected=0;

  outfile->cd();
  // create histos
  
  h_M_TopHad_reco=new TH1F("M_TopHad_reco","M_TopHad_reco",50,0,500);
  h_M_WHad_reco=new TH1F("M_WHad_reco","M_WHad_reco",30,0,300);
  h_M_Higgs_reco=new TH1F("M_Higgs_reco","M_Higgs_reco",30,0,300);
  h_M_TopLep_reco=new TH1F("M_TopLep_reco","M_TopLep_reco",50,0,500);

  h_M_TopHad_best=new TH1F("M_TopHad_best","M_TopHad_best",50,0,500);
  h_M_WHad_best=new TH1F("M_WHad_best","M_WHad_best",30,0,300);
  h_M_Higgs_best=new TH1F("M_Higgs_best","M_Higgs_best",30,0,300);
  h_M_TopLep_best=new TH1F("M_TopLep_best","M_TopLep_best",50,0,500);


}
MCTruthAnalyzer::~MCTruthAnalyzer(){
  cout << "nevents " << nevents << endl;
  cout << "nselected " << nselected << endl;

  outfile->cd();
  // write histos

  //  h_mcmatched_tth_me->Write();
  outfile->Write();
  outfile->Close();
}

void MCTruthAnalyzer::Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
		       const TLorentzVector& lepvec, const TVector2& metvec,
		       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, 
		       const TLorentzVector& vQ2_true, const TLorentzVector& vBLep_true, 
		       const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
		       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true){

  nevents++;
  nselected++;
  mcmatcher.Setup(vBHad_true,vQ1_true,vQ2_true,vBLep_true,vLep_true,vNu_true,vB1_true,vB2_true);
  vector<Interpretation*> ints = generator.GenerateTTHInterpretations(jetvecs,jetcsvs,lepvec,metvec);

  Interpretation* reco_int=0;
  Interpretation* best_int=0;

  // interpretation with highest ttH/ttbb/tt over combinatoric background likelihood
  float minDr=9999;
  for(uint i=0;i<ints.size();i++){
    float dr=mcmatcher.SumDrTTH(*(ints[i]));
    if(dr<minDr){
      minDr=dr;
      best_int=ints[i];
      if(mcmatcher.MatchNTTH(*(ints[i]))==6){
	reco_int=ints[i];
      }
    }
  }
  if(best_int!=0){
    h_M_TopHad_best->Fill(best_int->TopHad_M());
    h_M_WHad_best->Fill(best_int->WHad_M());
    h_M_Higgs_best->Fill(best_int->Higgs_M());
    h_M_TopLep_best->Fill(best_int->TopLep_M());
  }
  if(reco_int!=0){
    h_M_TopHad_reco->Fill(reco_int->TopHad_M());
    h_M_WHad_reco->Fill(reco_int->WHad_M());
    h_M_Higgs_reco->Fill(reco_int->Higgs_M());
    h_M_TopLep_reco->Fill(reco_int->TopLep_M());
  }

  for(uint i=0;i<ints.size();i++){
    delete ints[i];
  }
  ints.clear();
  
  float ht=0;
  TLorentzVector p4_total;
  for(uint i=0; i< jetvecs.size();i++){
    p4_total+=jetvecs[i];
    ht+=jetvecs[i].Pt();
  }
  int ntags=0;
  for(uint i=0; i< jetcsvs.size();i++){
    if(jetcsvs[i]>0.8) ntags++;
  }



}

