#include "../interface/Analyzer.hpp"

using namespace std;

Analyzer::Analyzer(string outfilename_,float min_combo_prob_tth_,float min_combo_prob_tt_,float min_combo_prob_ttbb_):outfilename(outfilename_),min_combo_prob_tth(min_combo_prob_tth_),min_combo_prob_tt(min_combo_prob_tt_),min_combo_prob_ttbb(min_combo_prob_ttbb_),generator(InterpretationGenerator(IntType::tth,0)){
  outfile=new TFile((outfilename+"_ana.root").c_str(),"RECREATE");        
  nevents=0;
  nselected=0;
  nmatched_tth=0;
  nmatched_ttbb=0;
  outfile->cd();
  h_best_tth_me_likelihood=new TH1F("best_tth_me_likelihood","best_tth_me_likelihood",50,-14,-9);
  h_best_ttbb_me_likelihood=new TH1F("best_ttbb_me_likelihood","best_ttbb_me_likelihood",50,-14,-9);
  h_best_ratio_me_likelihood=new TH1F("best_ratio_me_likelihood","best_ratio_me_likelihood",26,-0.1,1);
  h_best_tth_likelihood=new TH1F("best_tth_likelihood","best_tth_likelihood",40,-10,-6);
  h_best_ttbb_likelihood=new TH1F("best_ttbb_likelihood","best_ttbb_likelihood",40,-10,-6);
  h_best_ratio_likelihood=new TH1F("best_ratio_likelihood","best_ratio_likelihood",26,-0.1,1);
  h_best_tth_me=new TH1F("best_tth_me","best_tth_me",40,-6,-2);
  h_best_ttbb_me=new TH1F("best_ttbb_me","best_ttbb_me",40,-5,0);
  h_best_ratio_me=new TH1F("best_ratio_me","best_ratio_me",26,-0.1,1);
  h_best_m_higgs_tthreco=new TH1F("best_m_higgs_tthreco","best_m_higgs_tthreco",40,0,300);
  h_best_m_higgs_ttbbreco=new TH1F("best_m_higgs_ttbbreco","best_m_higgs_ttbbreco",40,0,300);
  h_best_m_higgs_ttreco=new TH1F("best_m_higgs_ttreco","best_m_higgs_ttreco",40,0,300);


}
Analyzer::~Analyzer(){
  cout << "nevents " << nevents << endl;
  cout << "nselected " << nselected << endl;
  cout << "nmatched_tth " << nmatched_tth << endl;
  cout << "nmatched_ttbb " << nmatched_ttbb << endl;
  outfile->cd();
  h_best_tth_me_likelihood->Write();
  h_best_ttbb_me_likelihood->Write();
  h_best_ratio_me_likelihood->Write();
  h_best_tth_likelihood->Write();
  h_best_ttbb_likelihood->Write();
  h_best_ratio_likelihood->Write();
  h_best_tth_me->Write();
  h_best_ttbb_me->Write();
  h_best_ratio_me->Write();
  h_best_m_higgs_tthreco->Write();
  h_best_m_higgs_ttbbreco->Write();
  h_best_m_higgs_ttreco->Write();

  outfile->Close();
}

void Analyzer::Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
				   const TLorentzVector& lepvec, const TVector2& metvec,
				   const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, 
				   const TLorentzVector& vQ2_true, const TLorentzVector& vBLep_true, 
				   const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
				   const TLorentzVector& vB1_true, const TLorentzVector& vB2_true){

  nevents++;

  mcmatcher.Setup(vBHad_true,vQ1_true,vQ2_true,vBLep_true,vLep_true,vNu_true,vB1_true,vB2_true);
  vector<Interpretation*> ints = generator.GenerateTTHInterpretations(jetvecs,jetcsvs,lepvec,metvec);
  float max_combo_prob_tth=0;
  float max_combo_prob_ttbb=0;
  float max_combo_prob_tt=0;

  Interpretation* best_combo_ttbb=0;
  Interpretation* best_combo_tth=0;
  Interpretation* best_combo_tt=0;
  vector<Interpretation*> good_tth_ints;
  vector<Interpretation*> good_ttbb_ints;
  vector<Interpretation*> good_tt_ints;
  const float min_prob=0.0;
  for(uint i=0;i<ints.size();i++){
    float combo_prob_tth=quality.TTWHLikelihood_comb_ratio(*(ints[i]));
    float combo_prob_ttbb=quality.TTWBBLikelihood_comb_ratio(*(ints[i]));
    float combo_prob_tt=quality.TTWLikelihood_comb_ratio(*(ints[i]));
    if(combo_prob_tth>max_combo_prob_tth){
      max_combo_prob_tth=combo_prob_tth;
      best_combo_tth=ints[i];
      if(combo_prob_tth>min_prob){
	good_tth_ints.push_back(ints[i]);
      }
    }
    if(combo_prob_ttbb>max_combo_prob_ttbb){
      max_combo_prob_ttbb=combo_prob_ttbb;
      best_combo_ttbb=ints[i];
      if(combo_prob_ttbb>min_prob){
	good_ttbb_ints.push_back(ints[i]);
      }
    }
    if(combo_prob_tt>max_combo_prob_tt){
      max_combo_prob_tt=combo_prob_tt;
      best_combo_tt=ints[i];
      if(combo_prob_tt>min_prob){
	good_tt_ints.push_back(ints[i]);
      }

    }
  }
  
  if(max_combo_prob_tth<min_combo_prob_tth) return;
  if(max_combo_prob_tt<min_combo_prob_tt) return;
  if(max_combo_prob_ttbb<min_combo_prob_ttbb) return;
  nselected++;

  //  cout << "intepreations " << ints.size() << endl;
  //  cout << "good tth intepreations " <<good_tth_ints.size() <<  endl;
  //  cout << "good ttbb intepreations " <<good_ttbb_ints.size() <<  endl;
  h_best_m_higgs_tthreco->Fill(best_combo_tth!=0?best_combo_tth->Higgs_M():1);
  h_best_m_higgs_ttreco->Fill(best_combo_tt!=0?best_combo_tt->Higgs_M():1);
  h_best_m_higgs_ttbbreco->Fill(best_combo_ttbb!=0?best_combo_ttbb->Higgs_M():1);

  if( best_combo_tth!=0 && mcmatcher.MatchTTHallQ(*best_combo_tth) ){
    nmatched_tth++;
  }
  if( best_combo_ttbb!=0 && mcmatcher.MatchTTHallQ(*best_combo_ttbb) ){
    nmatched_ttbb++;
  }
  if( best_combo_tt!=0 && mcmatcher.MatchTTHallQ(*best_combo_tt) ){
    nmatched_tt++;
  }

  double best_tth_me_likelihood=-1;
  double best_ttbb_me_likelihood=-1;
  double best_tth_likelihood=-1;
  double best_ttbb_likelihood=-1;
  double best_tth_me=-1;
  double best_ttbb_me=-1;

  for(uint i=0;i< good_tth_ints.size();i++){
    float tthlike= quality.TTWHLikelihood(*good_tth_ints[i]);
    float tthme = quality.TTHBB_ME(*good_tth_ints[i]);
    float tthlikeme = tthlike*tthme;
    if(best_tth_likelihood<tthlike){
      best_tth_likelihood=tthlike;
      best_tth_me=tthme;
    }
    if(best_tth_me_likelihood<tthlikeme){
      best_tth_me_likelihood=tthlikeme;
    }

  }
  for(uint i=0;i< good_ttbb_ints.size();i++){
    float ttbblike= quality.TTWBBLikelihood(*good_ttbb_ints[i]);
    float ttbbme = 1e4*quality.TTBB_ON_ME(*good_ttbb_ints[i]);
    float ttbblikeme = ttbblike*ttbbme;
    if(best_ttbb_likelihood<ttbblike){
      best_ttbb_likelihood=ttbblike;
      best_ttbb_me=ttbbme;
    }
    if(best_ttbb_me_likelihood<ttbblikeme){
      best_ttbb_me_likelihood=ttbblikeme;
    }

  }
  float best_ratio_me_likelihood=-1;
  float best_ratio_likelihood=-1;
  float best_ratio_me=-1;
  if(best_tth_me_likelihood>0&&best_ttbb_me_likelihood>0)
    best_ratio_me_likelihood = best_tth_me_likelihood/(best_ttbb_me_likelihood+best_tth_me_likelihood);
  if(best_tth_likelihood>0&&best_ttbb_likelihood>0)
    best_ratio_likelihood = best_tth_likelihood/(best_ttbb_likelihood+best_tth_likelihood);
  if(best_tth_likelihood>0&&best_ttbb_likelihood>0)
    best_ratio_me = best_tth_me/(best_ttbb_me+best_tth_me);
  h_best_tth_me_likelihood->Fill(best_tth_me_likelihood>0?log10(best_tth_me_likelihood):-1e20);
  h_best_ttbb_me_likelihood->Fill(best_ttbb_me_likelihood>0?log10(best_ttbb_me_likelihood):-1e20);
  h_best_ratio_me_likelihood->Fill(best_ratio_me_likelihood);
  h_best_tth_likelihood->Fill(best_tth_likelihood>0?log10(best_tth_likelihood):-1e20);
  h_best_ttbb_likelihood->Fill(best_ttbb_likelihood>0?log10(best_ttbb_likelihood):-1e20);
  h_best_ratio_likelihood->Fill(best_ratio_likelihood);
  h_best_tth_me->Fill(best_tth_me>0?log10(best_tth_me):-1e20);
  h_best_ttbb_me->Fill(best_ttbb_me>0?log10(best_ttbb_me):-1e20);
  h_best_ratio_me->Fill(best_ratio_me);
  
  for(uint i=0;i<ints.size();i++){
    delete ints[i];
  }
  ints.clear();

}

