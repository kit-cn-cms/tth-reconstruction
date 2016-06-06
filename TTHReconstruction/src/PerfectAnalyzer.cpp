#include "../interface/PerfectAnalyzer.hpp"

using namespace std;

PerfectAnalyzer::PerfectAnalyzer(string outfilename){
  outfile=new TFile((outfilename+"_perfect.root").c_str(),"RECREATE");        
  nevents=0;
  nselected=0;

  outfile->cd();
  // create histos
  h_njets = new TH1F("njets","njets",12,4.5,10.5);
  h_ntags = new TH1F("ntags","ntags",12,2.5,8.5);
  h_HT = new TH1F("HT","HT",40,0,2000);
  h_mass = new TH1F("mass","mass",60,0,3000);
  //  h_bquark_vs_tquark_pt = new TH1F("h_bquark_vs_tquark_pt","h_bquark_vs_tquark_pt",20,0,250,20,0,500);
  

  h_perfect_tth_me=new TH1F("perfect_tth_me","perfect_tth_me",40,-7,-2);
  h_perfect_tth_ps_me=new TH1F("perfect_tth_ps_me","perfect_tth_ps_me",40,-7,-2);
  h_perfect_tth_sc_me=new TH1F("perfect_tth_sc_me","perfect_tth_sc_me",40,-7,-2);
  h_perfect_tth_ratio_me=new TH1F("perfect_tth_ratio_me","perfect_tth_ratio_me",26,-0.1,1.1);
  h_perfect_tth_ratio_me_pdf=new TH1F("perfect_tth_ratio_me_pdf","perfect_tth_ratio_me_pdf",26,-0.1,1.1);
  h_perfect_ttbb_me=new TH1F("perfect_ttbb_me","perfect_ttbb_me",40,-5,-1);
  h_perfect_ratio_me=new TH1F("perfect_ratio_me","perfect_ratio_me",26,-0.1,1);
  h_perfect_tth_nd_me=new TH1F("perfect_tth_nd_me","perfect_tth_nd_me",40,-4,0);
  h_perfect_ttbb_off_me=new TH1F("perfect_ttbb_off_me","perfect_ttbb_off_me",40,-6,0);

  h_perfect_h_pt=new TH1F("perfect_h_pt","perfect_h_pt",25,0,500);
  h_perfect_deta_tt=new TH1F("perfect_deta_tt","perfect_deta_tt",25,0,5);
  h_perfect_dr_th=new TH1F("perfect_dr_th","perfect_dr_th",25,0,5);


}
PerfectAnalyzer::~PerfectAnalyzer(){
  cout << "nevents " << nevents << endl;
  cout << "nselected " << nselected << endl;

  outfile->cd();
  // write histos

  h_perfect_tth_me->Write();
  h_perfect_ttbb_me->Write();
  h_perfect_ratio_me->Write();

  h_perfect_tth_sc_me->Write();
  h_perfect_tth_ps_me->Write();
  h_perfect_tth_ratio_me->Write();
  h_perfect_tth_ratio_me_pdf->Write();

  h_perfect_tth_nd_me->Write();
  h_perfect_ttbb_off_me->Write();

  h_perfect_h_pt->Write();
  h_perfect_deta_tt->Write();
  h_perfect_dr_th->Write();


  outfile->Close();
}

void PerfectAnalyzer::Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
		       const TLorentzVector& lepvec, const TVector2& metvec,
		       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, 
		       const TLorentzVector& vQ2_true, const TLorentzVector& vBLep_true, 
		       const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
		       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true){

  nevents++;
  nselected++;

  Interpretation* perfect_int=new Interpretation(vBHad_true,1,vQ1_true,0,vQ2_true,0, vBLep_true,1, 
						 vLep_true, vNu_true, vB1_true,1, vB2_true,1);

  if(perfect_int!=0){
    h_perfect_tth_me->Fill(log10(quality.TTHBB_ME(*perfect_int)));
    h_perfect_tth_nd_me->Fill(log10(quality.TTH_ME(*perfect_int)));
    h_perfect_ttbb_me->Fill(log10(quality.TTBB_ON_ME(*perfect_int)));
    h_perfect_ttbb_off_me->Fill(log10(quality.TTBB_OFF_ME(*perfect_int)));
    h_perfect_ratio_me->Fill(quality.TTHBB_ME(*perfect_int)/(quality.TTBB_ON_ME(*perfect_int)+quality.TTHBB_ME(*perfect_int)));

    h_perfect_tth_ps_me->Fill(log10(quality.TTX0_PS_ME(*perfect_int)));
    h_perfect_tth_sc_me->Fill(log10(quality.TTX0_SC_ME(*perfect_int)));
    h_perfect_tth_ratio_me->Fill(quality.TTX0_SC_ME(*perfect_int)/(quality.TTX0_PS_ME(*perfect_int)+quality.TTX0_SC_ME(*perfect_int)));
    h_perfect_tth_ratio_me_pdf->Fill(quality.TTX0_SC_ME_PDF(*perfect_int)/(quality.TTX0_PS_ME_PDF(*perfect_int)+quality.TTX0_SC_ME_PDF(*perfect_int)));

    h_perfect_h_pt->Fill(perfect_int->Higgs().Pt());
    h_perfect_deta_tt->Fill(fabs(perfect_int->TopHad().Eta()-perfect_int->TopLep().Eta()));
    h_perfect_dr_th->Fill(fmin(perfect_int->TopHad().DeltaR(perfect_int->Higgs()),perfect_int->TopLep().DeltaR(perfect_int->Higgs())));

  }


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

  h_njets->Fill(jetvecs.size());
  h_ntags->Fill(ntags);
  h_HT->Fill(ht);
  h_mass->Fill(p4_total.M());


}

