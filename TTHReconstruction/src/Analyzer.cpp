#include "../interface/Analyzer.hpp"

using namespace std;

Analyzer::Analyzer(string outfilename):generator(InterpretationGenerator(IntType::tth,0)){
  outfile=new TFile((outfilename+"_ana.root").c_str(),"RECREATE");        
  nevents=0;
  nselected=0;
  nmatched_tth=0;
  nmatched_ttbb=0;
  nmatched_tth2=0;
  nmatched_ttbb2=0;

  outfile->cd();
  // create histos
  h_best_tth_me_likelihood=new TH1F("best_tth_me_likelihood","best_tth_me_likelihood",50,-20,-10);
  h_best_ttbb_me_likelihood=new TH1F("best_ttbb_me_likelihood","best_ttbb_me_likelihood",50,-20,-10);
  h_best_ratio_me_likelihood=new TH1F("best_ratio_me_likelihood","best_ratio_me_likelihood",26,-0.1,1);
  h_best_tth_likelihood=new TH1F("best_tth_likelihood","best_tth_likelihood",40,-16,-6);
  h_best_ttbb_likelihood=new TH1F("best_ttbb_likelihood","best_ttbb_likelihood",40,-16,-6);
  h_best_ratio_likelihood=new TH1F("best_ratio_likelihood","best_ratio_likelihood",26,-0.1,1);
  h_best_tth_me=new TH1F("best_tth_me","best_tth_me",40,-6,-2);
  h_best_ttbb_me=new TH1F("best_ttbb_me","best_ttbb_me",40,-8,-2);
  h_best_ratio_me=new TH1F("best_ratio_me","best_ratio_me",26,-0.1,1);

  h_best_m_higgs_tthreco=new TH1F("best_m_higgs_tthreco","best_m_higgs_tthreco",40,0,300);
  h_best_m_higgs_ttbbreco=new TH1F("best_m_higgs_ttbbreco","best_m_higgs_ttbbreco",40,0,300);
  h_best_m_higgs_ttreco=new TH1F("best_m_higgs_ttreco","best_m_higgs_ttreco",40,0,300);

  h_best_m_higgs_tthreco2=new TH1F("best_m_higgs_tthreco2","best_m_higgs_tthreco2",40,0,300);
  h_best_m_higgs_ttbbreco2=new TH1F("best_m_higgs_ttbbreco2","best_m_higgs_ttbbreco2",40,0,300);
  h_best_m_higgs_ttreco2=new TH1F("best_m_higgs_ttreco2","best_m_higgs_ttreco2",40,0,300);

  h_sum_tth_me_likelihood=new TH1F("sum_tth_me_likelihood","sum_tth_me_likelihood",50,-20,-10);
  h_sum_ttbb_me_likelihood=new TH1F("sum_ttbb_me_likelihood","sum_ttbb_me_likelihood",50,-20,-10);
  h_sum_ratio_me_likelihood=new TH1F("sum_ratio_me_likelihood","sum_ratio_me_likelihood",26,-0.1,1);
  h_sum_tth_likelihood=new TH1F("sum_tth_likelihood","sum_tth_likelihood",40,-14,-6);
  h_sum_ttbb_likelihood=new TH1F("sum_ttbb_likelihood","sum_ttbb_likelihood",40,-14,-6);
  h_sum_ratio_likelihood=new TH1F("sum_ratio_likelihood","sum_ratio_likelihood",26,-0.1,1);
  h_sum_tth_me=new TH1F("sum_tth_me","sum_tth_me",40,-6,-2);
  h_sum_ttbb_me=new TH1F("sum_ttbb_me","sum_ttbb_me",40,-8,-2);
  h_sum_ratio_me=new TH1F("sum_ratio_me","sum_ratio_me",26,-0.1,1);

}
Analyzer::~Analyzer(){
  cout << "nevents " << nevents << endl;
  cout << "nselected " << nselected << endl;
  cout << "nmatched_tth " << nmatched_tth << endl;
  cout << "nmatched_ttbb " << nmatched_ttbb << endl;
  cout << "nmatched_tth2 " << nmatched_tth2 << endl;
  cout << "nmatched_ttbb2 " << nmatched_ttbb2 << endl;

  outfile->cd();
  // write histos
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
  h_best_m_higgs_tthreco2->Write();
  h_best_m_higgs_ttbbreco2->Write();
  h_best_m_higgs_ttreco2->Write();
  h_sum_tth_me_likelihood->Write();
  h_sum_ttbb_me_likelihood->Write();
  h_sum_ratio_me_likelihood->Write();
  h_sum_tth_likelihood->Write();
  h_sum_ttbb_likelihood->Write();
  h_sum_ratio_likelihood->Write();
  h_sum_tth_me->Write();
  h_sum_ttbb_me->Write();
  h_sum_ratio_me->Write();

  outfile->Close();
}

void Analyzer::Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
		       const TLorentzVector& lepvec, const TVector2& metvec,
		       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, 
		       const TLorentzVector& vQ2_true, const TLorentzVector& vBLep_true, 
		       const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
		       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true){

  nevents++;
  nselected++;
  mcmatcher.Setup(vBHad_true,vQ1_true,vQ2_true,vBLep_true,vLep_true,vNu_true,vB1_true,vB2_true);
  vector<Interpretation*> ints = generator.GenerateTTHInterpretations(jetvecs,jetcsvs,lepvec,metvec);
  // interpretation with highest ttH/ttbb/tt likelihood
  Interpretation* best_int_ttbb_likelihood=0;
  Interpretation* best_int_tth_likelihood=0;
  Interpretation* best_int_tt_likelihood=0;
  double best_tth_likelihood=-1;
  double best_ttbb_likelihood=-1;
  double best_tt_likelihood=-1;

  // interpretation with highest ttH/ttbb/tt over combinatoric background likelihood
  Interpretation* best_int_ttbb_comb_likelihoodratio=0;
  Interpretation* best_int_tth_comb_likelihoodratio=0;
  Interpretation* best_int_tt_comb_likelihoodratio=0;
  double best_tth_comb_likelihoodratio=-1;
  double best_ttbb_comb_likelihoodratio=-1;
  double best_tt_comb_likelihoodratio=-1;
  
  // best tth*ME
  double best_tth_me_likelihood=-1;
  double best_ttbb_me_likelihood=-1;
  // ME for best likelihood
  double tth_me_best_tth_likelihood=-1;
  double ttbb_me_best_ttbb_likelihood=-1;

  // sum of likelihoods
  double sum_tth_likelihood=0;
  double sum_ttbb_likelihood=0; 
  double sum_tth_me_likelihood=0;
  double sum_ttbb_me_likelihood=0;
  double sum_tth_me=0;
  double sum_ttbb_me=0;

  // loop over all interpretations
  for(uint i=0;i<ints.size();i++){
    float tthlike= quality.TTWHLikelihood(*ints[i]);
    float ttbblike= quality.TTWBBLikelihood(*ints[i]);
    float ttlike= quality.TTWLikelihood(*ints[i]);
    float tth_comblike= quality.TTWHLikelihood_comb(*ints[i]);
    float ttbb_comblike= quality.TTWBBLikelihood_comb(*ints[i]);
    float tt_comblike= quality.TTWLikelihood_comb(*ints[i]);
    // interpretations wiht highest likelihood
    if(best_tth_likelihood<tthlike){
      best_tth_likelihood=tthlike;
      best_int_tth_likelihood=ints[i];
    }
    if(best_ttbb_likelihood<ttbblike){
      best_ttbb_likelihood=ttbblike;
      best_int_ttbb_likelihood=ints[i];
    }
    if(best_tt_likelihood<ttlike){
      best_tt_likelihood=ttlike;
      best_int_tt_likelihood=ints[i];
    }
    // interpretations wiht highest likelihood / combinatorics likelihood ratio
    float tthlikeratio=tthlike/(tthlike+tth_comblike);
    if(best_tth_comb_likelihoodratio<tthlikeratio){
      best_tth_comb_likelihoodratio=tthlike;
      best_int_tth_comb_likelihoodratio=ints[i];
    }
    float ttbblikeratio=ttbblike/(ttbblike+ttbb_comblike);
    if(best_ttbb_comb_likelihoodratio<ttbblikeratio){
      best_ttbb_comb_likelihoodratio=ttbblike;
      best_int_ttbb_comb_likelihoodratio=ints[i];
    }
    float ttlikeratio=ttlike/(ttlike+tt_comblike);
    if(best_tt_comb_likelihoodratio<ttlikeratio){
      best_tt_comb_likelihoodratio=ttlike;
      best_int_tt_comb_likelihoodratio=ints[i];
    }
    sum_tth_likelihood+=tthlike;
    sum_ttbb_likelihood+=ttbblike;
  }
  tth_me_best_tth_likelihood=quality.TTHBB_ME(*best_int_tth_likelihood);
  // scale ttbb ME up  -- at 125 GeV ttH has larger cross section
  ttbb_me_best_ttbb_likelihood=quality.TTBB_ON_ME(*best_int_ttbb_likelihood);

  // loop over selected interpretations
  for(uint i=0;i<ints.size();i++){      
    float tthlike= quality.TTWHLikelihood(*ints[i]);
    float ttbblike= quality.TTWBBLikelihood(*ints[i]);
    // skip unlikely interpretations
    //    if(tthlike<best_tth_likelihood/100) continue;
    //    if(ttbblike<best_ttbb_likelihood/100) continue;
    // calculate tth me
    float tthme = quality.TTH_ME(*ints[i]);
    float tthlikeme = tthlike*tthme;
    float ttbbme = quality.TTBB_ON_ME(*ints[i]);
    float ttbblikeme = ttbblike*ttbbme;
    if(best_tth_me_likelihood<tthlikeme){
      best_tth_me_likelihood=tthlikeme;
    }
    if(best_ttbb_me_likelihood<ttbblikeme){
      best_ttbb_me_likelihood=ttbblikeme;
    }
    sum_tth_me+=tthme;
    sum_tth_me_likelihood+=tthlikeme;
    sum_ttbb_me+=ttbbme;
    sum_ttbb_me_likelihood+=ttbblikeme;

  }  

  h_best_m_higgs_tthreco->Fill(best_int_tth_likelihood!=0?best_int_tth_likelihood->Higgs_M():1);
  h_best_m_higgs_ttreco->Fill(best_int_tt_likelihood!=0?best_int_tt_likelihood->Higgs_M():1);
  h_best_m_higgs_ttbbreco->Fill(best_int_ttbb_likelihood!=0?best_int_ttbb_likelihood->Higgs_M():1);
  h_best_m_higgs_tthreco2->Fill(best_int_tth_comb_likelihoodratio!=0?best_int_tth_comb_likelihoodratio->Higgs_M():1);
  h_best_m_higgs_ttreco2->Fill(best_int_tt_comb_likelihoodratio!=0?best_int_tt_comb_likelihoodratio->Higgs_M():1);
  h_best_m_higgs_ttbbreco2->Fill(best_int_ttbb_comb_likelihoodratio!=0?best_int_ttbb_comb_likelihoodratio->Higgs_M():1);

  if( best_int_tth_likelihood!=0 && mcmatcher.MatchTTHallQ(*best_int_tth_likelihood) ){
    nmatched_tth++;
  }
  if( best_int_ttbb_likelihood!=0 && mcmatcher.MatchTTHallQ(*best_int_ttbb_likelihood) ){
    nmatched_ttbb++;
  }
  if( best_int_tt_likelihood!=0 && mcmatcher.MatchTTHallQ(*best_int_tt_likelihood) ){
    nmatched_tt++;
  }
  if( best_int_tth_comb_likelihoodratio!=0 && mcmatcher.MatchTTHallQ(*best_int_tth_comb_likelihoodratio) ){
    nmatched_tth2++;
  }
  if( best_int_ttbb_comb_likelihoodratio!=0 && mcmatcher.MatchTTHallQ(*best_int_ttbb_comb_likelihoodratio) ){
    nmatched_ttbb2++;
  }
  if( best_int_tt_comb_likelihoodratio!=0 && mcmatcher.MatchTTHallQ(*best_int_tt_comb_likelihoodratio) ){
    nmatched_tt2++;
  }

  float best_ratio_me_likelihood=-1;
  float best_ratio_likelihood=-1;
  float best_ratio_me=-1;
  if(best_tth_me_likelihood>0&&best_ttbb_me_likelihood>0)
    best_ratio_me_likelihood = best_tth_me_likelihood/(best_ttbb_me_likelihood+best_tth_me_likelihood);
  if(best_tth_likelihood>0&&best_ttbb_likelihood>0)
    best_ratio_likelihood = best_tth_likelihood/(best_ttbb_likelihood+best_tth_likelihood);
  if(tth_me_best_tth_likelihood>0&&ttbb_me_best_ttbb_likelihood>0)
    best_ratio_me = tth_me_best_tth_likelihood/(tth_me_best_tth_likelihood+ttbb_me_best_ttbb_likelihood);


  h_best_tth_me_likelihood->Fill(best_tth_me_likelihood>0?log10(best_tth_me_likelihood):-1e20);
  h_best_ttbb_me_likelihood->Fill(best_ttbb_me_likelihood>0?log10(best_ttbb_me_likelihood):-1e20);
  h_best_ratio_me_likelihood->Fill(best_ratio_me_likelihood);
  h_best_tth_likelihood->Fill(best_tth_likelihood>0?log10(best_tth_likelihood):-1e20);
  h_best_ttbb_likelihood->Fill(best_ttbb_likelihood>0?log10(best_ttbb_likelihood):-1e20);
  h_best_ratio_likelihood->Fill(best_ratio_likelihood);
  h_best_tth_me->Fill(tth_me_best_tth_likelihood>0?log10(tth_me_best_tth_likelihood):-1e20);
  h_best_ttbb_me->Fill(ttbb_me_best_ttbb_likelihood>0?log10(ttbb_me_best_ttbb_likelihood):-1e20);
  h_best_ratio_me->Fill(best_ratio_me);

  float sum_ratio_me_likelihood=0;
  float sum_ratio_likelihood=0;
  float sum_ratio_me=0;
  if(sum_tth_me_likelihood>0&&sum_ttbb_me_likelihood>0)
    sum_ratio_me_likelihood = sum_tth_me_likelihood/(sum_ttbb_me_likelihood+sum_tth_me_likelihood);
  if(sum_tth_likelihood>0&&sum_ttbb_likelihood>0)
    sum_ratio_likelihood = sum_tth_likelihood/(sum_ttbb_likelihood+sum_tth_likelihood);
  if(sum_tth_likelihood>0&&sum_ttbb_likelihood>0)
    sum_ratio_me = sum_tth_me/(sum_ttbb_me+sum_tth_me);

  if(sum_tth_likelihood>0&&sum_ttbb_likelihood>0)

  h_sum_tth_me_likelihood->Fill(sum_tth_me_likelihood>0?log10(sum_tth_me_likelihood):-1e20);
  h_sum_ttbb_me_likelihood->Fill(sum_ttbb_me_likelihood>0?log10(sum_ttbb_me_likelihood):-1e20);
  h_sum_ratio_me_likelihood->Fill(sum_ratio_me_likelihood);
  h_sum_tth_likelihood->Fill(sum_tth_likelihood>0?log10(sum_tth_likelihood):-1e20);
  h_sum_ttbb_likelihood->Fill(sum_ttbb_likelihood>0?log10(sum_ttbb_likelihood):-1e20);
  h_sum_ratio_likelihood->Fill(sum_ratio_likelihood);
  h_sum_tth_me->Fill(sum_tth_me>0?log10(sum_tth_me):-1e20);
  h_sum_ttbb_me->Fill(sum_ttbb_me>0?log10(sum_ttbb_me):-1e20);
  h_sum_ratio_me->Fill(sum_ratio_me);
  
  for(uint i=0;i<ints.size();i++){
    delete ints[i];
  }
  ints.clear();

}

