#include "../interface/CouplingAnalyzer.hpp"

using namespace std;

CouplingAnalyzer::CouplingAnalyzer(string outfilename):generator(InterpretationGenerator(IntType::tth,0)){
  outfile=new TFile((outfilename+"_ana.root").c_str(),"RECREATE");        
  nevents=0;
  nselected=0;

  outfile->cd();
  tree = new TTree("MVATree","MVATree");
  // create histos

  h_ratio_me_best_likelihood=new TH1F("ratio_me_best_likelihood","ratio_me_best_likelihood",26,-0.1,1);
  h_ratio_me_best_csv_likelihood=new TH1F("ratio_me_best_csv_likelihood","ratio_me_best_csv_likelihood",26,-0.1,1);
  h_ratio_me_bestish_likelihood=new TH1F("ratio_me_bestish_likelihood","ratio_me_bestish_likelihood",26,-0.1,1);

  InitVar("Evt_MERatio_Best");
  InitVar("Evt_MERatioCSV_Best");
  InitVar("Evt_MaxPt_TaggedJetPair");
  InitVar("Evt_MaxDeta_TaggedJets");
  InitVar("Evt_Dr_Lepton_Jet2");
  InitVar("Evt_BoostedMERatio");

}
CouplingAnalyzer::~CouplingAnalyzer(){
  cout << "nevents " << nevents << endl;
  //  tree->Write();
  outfile->Write();
  
  outfile->Close();
}

void CouplingAnalyzer::Analyze(const std::vector<TLorentzVector>& jetvecs, const std::vector<float>& jetcsvs, 
		       const TLorentzVector& lepvec, const TVector2& metvec,
		       const TLorentzVector& vBHad_true,const TLorentzVector& vQ1_true, 
		       const TLorentzVector& vQ2_true, const TLorentzVector& vBLep_true, 
		       const TLorentzVector& vLep_true, const TLorentzVector& vNu_true, 
		       const TLorentzVector& vB1_true, const TLorentzVector& vB2_true,
		       const TLorentzVector& vBHad_boosted,const TLorentzVector& vQ1_boosted, 
		       const TLorentzVector& vQ2_boosted, const TLorentzVector& vBLep_boosted, 
		       const TLorentzVector& vH_boosted){

  nevents++;
  nselected++;
  mcmatcher.Setup(vBHad_true,vQ1_true,vQ2_true,vBLep_true,vLep_true,vNu_true,vB1_true,vB2_true);
  vector<Interpretation*> ints = generator.GenerateTTHInterpretations(jetvecs,jetcsvs,lepvec,metvec);
  // interpretation with highest ttH/ttbb/tt likelihood
  Interpretation* best_int_likelihood_psme=0;
  Interpretation* best_int_likelihood_scme=0;
  Interpretation* bestish_int_likelihood_psme=0;
  Interpretation* bestish_int_likelihood_scme=0;
  Interpretation* best_int_likelihood_csv_psme=0;
  Interpretation* best_int_likelihood_csv_scme=0;
  float pslike_best=-999.;
  float sclike_best=-999.;
  float pslikish_best=-999.;
  float sclikish_best=-999.;
  float pscsvlike_best=-999.;
  float sccsvlike_best=-999.;
  
  for(uint i=0;i<ints.size();i++){
    float pslike= quality.TTWHLikelihoodTimesPSME(*ints[i]);
    float sclike= quality.TTWHLikelihoodTimesSCME(*ints[i]);
    float pslikish= quality.TTWHishLikelihoodTimesPSME(*ints[i]);
    float sclikish= quality.TTWHishLikelihoodTimesSCME(*ints[i]);
    float pscsvlike= quality.TTWHLikelihoodTimesCSVTimesPSME(*ints[i]);
    float sccsvlike= quality.TTWHLikelihoodTimesCSVTimesSCME(*ints[i]);

    if(pslike_best<pslike){
      pslike_best=pslike;
      best_int_likelihood_psme=ints[i];
    }
    if(sclike_best<sclike){
      sclike_best=sclike;
      best_int_likelihood_scme=ints[i];
    }
    if(pslikish_best<pslikish){
      pslikish_best=pslikish;
      bestish_int_likelihood_psme=ints[i];
    }
    if(sclikish_best<sclikish){
      sclikish_best=sclikish;
      bestish_int_likelihood_scme=ints[i];
    }
    if(pscsvlike_best<pscsvlike){
      pscsvlike_best=pscsvlike;
      best_int_likelihood_csv_psme=ints[i];
    }
    if(sccsvlike_best<sccsvlike){
      sccsvlike_best=sccsvlike;
      best_int_likelihood_csv_scme=ints[i];
    }

  }
  float ratio1=-0.01;
  float ratio2=-0.01;
  float ratio3=-0.01;
  if(best_int_likelihood_scme!=0  && best_int_likelihood_psme!=0){
    float sc=quality.TTWHLikelihoodTimesSCME(*best_int_likelihood_scme);
    float ps=quality.TTWHLikelihoodTimesPSME(*best_int_likelihood_psme);
    if(ps>0 && sc>0){
      ratio1=sc/(ps+sc);
    }
    else if(sc>0){
      ratio1=1;
    }
    else if(ps>0){
      ratio1=0;
    }

  }
  if(bestish_int_likelihood_scme!=0  && bestish_int_likelihood_psme!=0){
    float sc=quality.TTWHishLikelihoodTimesSCME(*best_int_likelihood_scme);
    float ps=quality.TTWHishLikelihoodTimesPSME(*best_int_likelihood_psme);
    if(ps>0 && sc>0){
      ratio2=sc/(ps+sc);
    }
    else if(sc>0){
      ratio1=1;
    }
    else if(ps>0){
      ratio1=0;
    }

  }
  if(best_int_likelihood_csv_scme!=0  && best_int_likelihood_csv_psme!=0){
    float sc=quality.TTWHLikelihoodTimesCSVTimesSCME(*best_int_likelihood_csv_scme);
    float ps=quality.TTWHLikelihoodTimesCSVTimesPSME(*best_int_likelihood_csv_psme);
    if(ps>0 && sc>0){
      ratio3=sc/(ps+sc);
    }
    else if(sc>0){
      ratio3=1;
    }
    else if(ps>0){
      ratio3=0;
    }

  }

  h_ratio_me_best_likelihood->Fill(ratio1);
  h_ratio_me_bestish_likelihood->Fill(ratio2);
  h_ratio_me_best_csv_likelihood->Fill(ratio3);
  cout << ratio1 << endl;
  cout << ratio2 << endl;
  cout << "<<<<" << endl;
  for(uint i=0;i<ints.size();i++){
    delete ints[i];
  }
  ints.clear();
  FillVar("Evt_MERatio_Best",ratio1);
  FillVar("Evt_MERatioCSV_Best",ratio3);

  vector<TLorentzVector> tagged_jets;
  for(uint i=0; i< jetvecs.size(); i++){
    if(jetcsvs[i]<0.8) tagged_jets.push_back(jetvecs[i]);
  }
  float maxpt=-1;
  float maxdeta=-1;
  float mindr2=9;
  for(uint i=0; i< tagged_jets.size(); i++){
    for(uint j=i+1; j< tagged_jets.size(); j++){
      float pt=(tagged_jets[i]+tagged_jets[j]).Pt();
      if(pt>maxpt){
	maxpt=pt;
      }

      float deta=fabs(tagged_jets[i].Eta()-tagged_jets[j].Eta());
      if(deta>maxdeta){
	maxdeta=deta;
      }

      float dr1=lepvec.DeltaR(tagged_jets[i]);
      float dr2=lepvec.DeltaR(tagged_jets[j]);
      if(dr1>dr2 && dr1 < mindr2){
	mindr2=dr1;
      }
      else if(dr1<dr2 && dr2 < mindr2){
	mindr2=dr2;
      }
	
    }
  }
  
  FillVar("Evt_MaxPt_TaggedJetPair",maxpt);
  FillVar("Evt_MaxDeta_TaggedJets",maxdeta);
  FillVar("Evt_Dr_Lepton_Jet2",mindr2);
  
  TLorentzVector thad=vBHad_boosted+vQ2_boosted+vQ1_boosted;
  TLorentzVector h=vH_boosted;

  TLorentzVector nuvec1;
  TLorentzVector nuvec2; 
  InterpretationGenerator::GetNuVecs(lepvec,metvec,nuvec1,nuvec2);
  TLorentzVector tlep1=nuvec1+lepvec+vBLep_boosted;
  TLorentzVector tlep2=nuvec2+lepvec+vBLep_boosted;
  TLorentzVector tlep;
  if(fabs(tlep1.M()-173.)<fabs(tlep2.M()-173.)){
    tlep=tlep1;
  }
  else {
    tlep=tlep2;
  }
  cout << thad.Pt() << " : " << tlep.Pt() << " : " << h.Pt() << endl;
  float scme=quality.TTX0_SC_ME(thad,tlep,h);
  float psme=quality.TTX0_PS_ME(thad,tlep,h);
  FillVar("Evt_BoostedMERatio",scme/(scme+psme));
  tree->Fill();

}

void CouplingAnalyzer::InitVar(string name){
  vars[name]=-99.;
  tree->Branch(name.c_str(), &(vars[name]), (name+"/F").c_str() );
}
void CouplingAnalyzer::FillVar(string name,float value){
  vars[name]=value;
}
