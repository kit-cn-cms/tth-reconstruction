#include "../interface/Interpretation.hpp"
#include <iostream>
using namespace std;

Interpretation::Interpretation(TLorentzVector b_had_, float csv_b_had_, TLorentzVector q1_, float csv_q1_, TLorentzVector q2_, float csv_q2_, TLorentzVector b_lep_, float csv_b_lep_, TLorentzVector lep_, TLorentzVector nu_, TLorentzVector b1_, float csv_b1_, TLorentzVector b2_, float csv_b2_){
  b1=b1_;
  b2=b2_;
  b_had=b_had_;
  q1=q1_;
  q2=q2_;
  b_lep=b_lep_;
  lep=lep_;
  nu=nu_;
  csv_b1=csv_b1_;
  csv_b2=csv_b2_;
  csv_b_had=csv_b_had_;
  csv_b_lep=csv_b_lep_;
  csv_q1=csv_q1_;
  csv_q2=csv_q2_;
  //dervived
  w_had=q1+q2;
  t_had=b_had+w_had;
  w_lep=lep+nu;
  t_lep=b_lep+w_lep;
  higgs=b1+b2;

}

Interpretation::Interpretation(TLorentzVector b_had_, float csv_b_had_, TLorentzVector q1_, float csv_q1_, TLorentzVector q2_, float csv_q2_, TLorentzVector b_lep_, float csv_b_lep_, TLorentzVector lep_, TLorentzVector nu_){
  b_had=b_had_;
  q1=q1_;
  q2=q2_;
  b_lep=b_lep_;
  lep=lep_;
  nu=nu_;
  csv_b1=-2;
  csv_b2=-2;
  csv_b_had=csv_b_had_;
  csv_b_lep=csv_b_lep_;
  csv_q1=csv_q1_;
  csv_q2=csv_q2_;
  //dervived
  w_had=q1+q2;
  t_had=b_had+w_had;
  w_lep=lep+nu;
  t_lep=b_lep+w_lep;
}


TLorentzVector Interpretation::BHad(){
  return b_had;
}
TLorentzVector Interpretation::Q1(){
  return q1;
}
TLorentzVector Interpretation::Q2(){
  return q2;
}
TLorentzVector Interpretation::BLep(){
  return b_lep;
}
TLorentzVector Interpretation::Lep(){
  return lep;
}
TLorentzVector Interpretation::Nu(){
  return nu;
}
TLorentzVector Interpretation::TopHad(){
  return t_had;
}
TLorentzVector Interpretation::TopLep(){
  return t_lep;
}
TLorentzVector Interpretation::B1(){
  return b1;
}
TLorentzVector Interpretation::B2(){
  return b2;
}

TLorentzVector Interpretation::Higgs(){
  return higgs;
}
float Interpretation::TopHad_M(){
  return t_had.M();
}
float Interpretation::TopLep_M(){
  return t_lep.M();
}
float Interpretation::WHad_M(){
  return w_had.M();
}
float Interpretation::WLep_M(){
  return w_lep.M();
}
float Interpretation::Higgs_M(){
  return higgs.M();
}
float Interpretation::BHad_CSV(){
  return csv_b_had;
}
float Interpretation::BLep_CSV(){
  return csv_b_lep;
}
float Interpretation::Q1_CSV(){
  return csv_q1;
}
float Interpretation::Q2_CSV(){
  return csv_q2;
}
float Interpretation::B1_CSV(){
  return csv_b1;
}
float Interpretation::B2_CSV(){
  return csv_b2;
}
void Interpretation::SetTag(std::string tag, float value){
  tags[tag]=value;
}
float Interpretation::GetValue(std::string tag){
  return tags[tag];
}

void Interpretation::GetNuVecs(const TLorentzVector & lepvec, const TVector2 & metvec, TLorentzVector & nu1, TLorentzVector & nu2){
  double metvec2 = metvec.Px()*metvec.Px() + metvec.Py()*metvec.Py();
  double mu = (80.4*80.4)/2 + metvec.Px()*lepvec.Px() + metvec.Py()*lepvec.Py();
  double a = (mu*lepvec.Pz())/(lepvec.E()*lepvec.E() - lepvec.Pz()*lepvec.Pz());
  double a2 = TMath::Power(a, 2);
  double b = (TMath::Power(lepvec.E(), 2.)*metvec2 - TMath::Power(mu, 2.)) / (TMath::Power(lepvec.E(), 2)- TMath::Power(lepvec.Pz(), 2));
  float pz1,pz2;
  if (a2-b < 0) { 
    pz1 = a;
    pz2 = a;
  } else {
    double root = sqrt(a2-b);
    pz1 = a + root;
    pz2 = a - root;
  }
  nu1.SetPxPyPzE(metvec.Px(),metvec.Py(),pz1,sqrt(metvec.Mod2()+pz1*pz1));
  nu2.SetPxPyPzE(metvec.Px(),metvec.Py(),pz2,sqrt(metvec.Mod2()+pz2*pz2));
}


std::vector<Interpretation*> Interpretation::GenerateTTHInterpretations(std::vector<TLorentzVector> jetvecs_in, std::vector<float> jetcsvs_in, TLorentzVector lepvec, TVector2 metvec){

  std::vector<Interpretation*> interpretations;

  // Neutrino reconstruction
  TLorentzVector nuvec1;
  TLorentzVector nuvec2;
  Interpretation::GetNuVecs(lepvec,metvec,nuvec1,nuvec2);
  uint nNu=2;
  if(fabs(nuvec1.Pz()-nuvec2.Pz()<1)) nNu=1;

  // skim jets
  std::vector<TLorentzVector> jetvecs;
  std::vector<float> jetcsvs;
  for(uint i=0; i<jetvecs_in.size(); i++){
    jetvecs.push_back(jetvecs_in[i]);
    jetcsvs.push_back(jetcsvs_in[i]);
  }
  int ntags=0;
  for(uint i=0; i<jetcsvs.size(); i++){
    if(jetcsvs[i]>0.814) ntags++;
  }
  for(uint iQ1=0; iQ1<jetvecs.size();iQ1++){
    for(uint iQ2=0; iQ2<jetvecs.size();iQ2++){
      if(iQ2<iQ1){
	for(uint iBHad=0; iBHad<jetvecs.size();iBHad++){    
	  if(iBHad!=iQ1 && iBHad!=iQ2){
	    for(uint iBLep=0; iBLep<jetvecs.size();iBLep++){
	      if(iBLep!=iQ1 && iBLep!=iQ2 && iBLep!=iBHad) {
		for(uint iB1=0; iB1<jetvecs.size();iB1++){
		  if(iB1!=iQ1 && iB1!=iQ2 && iB1!=iBHad && iB1!=iBLep) {
		    for(uint iB2=0; iB2<jetvecs.size();iB2++){
		      if(iB2!=iQ1 && iB2!=iQ2 && iB2!=iBHad && iB2!=iBLep && iB2<iB1){
			for(uint iNu=1; iNu<=nNu;iNu++){
			  
			  if(int(jetcsvs[iBHad]>0.814)+int(jetcsvs[iBLep]>0.814)+int(jetcsvs[iB1]>0.814)+int(jetcsvs[iB2]>0.814)<ntags-1) continue;

			  interpretations.push_back(new Interpretation(jetvecs[iBHad],jetcsvs[iBHad],
								       jetvecs[iQ1],jetcsvs[iQ1],
								       jetvecs[iQ2],jetcsvs[iQ2],
								       jetvecs[iBLep],jetcsvs[iBLep],
								       lepvec,iNu==1?nuvec1:nuvec2,
								       jetvecs[iB1],jetcsvs[iB1],
								       jetvecs[iB2],jetcsvs[iB2]));
						    
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }

  return interpretations;
}

std::vector<Interpretation*> Interpretation::GenerateTTInterpretations(std::vector<TLorentzVector> jetvecs_in, std::vector<float> jetcsvs_in, TLorentzVector lepvec, TVector2 metvec){

  std::vector<Interpretation*> interpretations;

  // Neutrino reconstruction
  TLorentzVector nuvec1;
  TLorentzVector nuvec2;
  Interpretation::GetNuVecs(lepvec,metvec,nuvec1,nuvec2);
  uint nNu=2;
  if(fabs(nuvec1.Pz()-nuvec2.Pz()<1)) nNu=1;

  // skim jets
  std::vector<TLorentzVector> jetvecs;
  std::vector<float> jetcsvs;
  for(uint i=0; i<jetvecs_in.size(); i++){
    jetvecs.push_back(jetvecs_in[i]);
    jetcsvs.push_back(jetcsvs_in[i]);
  }
  int ntags=0;
  for(uint i=0; i<jetcsvs.size(); i++){
    if(jetcsvs[i]>0.814) ntags++;
  }
  int nints=0;
  for(uint iQ1=0; iQ1<jetvecs.size();iQ1++){
    for(uint iQ2=0; iQ2<jetvecs.size();iQ2++){
      if(iQ2!=iQ1){
	for(uint iBHad=0; iBHad<jetvecs.size();iBHad++){    
	  if(iBHad!=iQ1 && iBHad!=iQ2){
	    for(uint iBLep=0; iBLep<jetvecs.size();iBLep++){
	      if(iBLep!=iQ1 && iBLep!=iQ2 && iBLep!=iBHad) {
		for(uint iNu=1; iNu<=nNu;iNu++){
		  if(int(jetcsvs[iBHad]>0.814)+int(jetcsvs[iBLep]>0.814) < 1 ) continue;
		  interpretations.push_back(new Interpretation(jetvecs[iBHad],jetcsvs[iBHad],
							       jetvecs[iQ1],jetcsvs[iQ1],
							       jetvecs[iQ2],jetcsvs[iQ2],
							       jetvecs[iBLep],jetcsvs[iBLep],
							       lepvec,iNu==1?nuvec1:nuvec2));
		  nints++;

		  //		  cout << iBHad << iQ1 << iQ2 << iBLep << iNu << endl;
		  //		  cout << "q1 pt "<< jetvecs[iQ1].Pt() << " q2 pt " << jetvecs[iQ2].Pt() << " blep pt " << jetvecs[iBLep].Pt()  << " bhad pt " << jetvecs[iBHad].Pt() << " nu pz " << (iNu==1?nuvec1.Pz():nuvec2.Pz())  << endl;
		}
	      }
	    }
	  }
	}
      }
    }
  }
  //  std::cout << "njets "  << jetvecs.size() << " ints: "<< nints  <<std::endl;
  return interpretations;
}
