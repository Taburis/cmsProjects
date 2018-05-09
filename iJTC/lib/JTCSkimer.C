
#include "JTCSkimer.h"
#include "../test/corrTable/bJTCTrkCorr_pp.h"
#include "TF1.h"

float trketamaxcut = 2.;
float trkPtcut = 1.0;
float jetetacut = 1.6;
float jetptcut = 120;
float vzcut = 15;

double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
double ppPthatEntries[9] = {0,0,272902,377559,467823,447683,259111,234347,50942};

auto ppTrkCorrection = new bJTCTrkCorr_pp("../test/corrTable/PYTHIA_noDCA_mulWeighted_1stMay.root");

void JTCSkimer::eventWeight(){
		weights=1;
		if(!isHI){
				int ibin=0;
				while(pthat>pthatbins[ibin+1]) ibin++;
				weights=(xsecs[ibin]-xsecs[ibin+1])/ppPthatEntries[ibin];
				weights=weights*(1.0/ppVz->Eval(vz));
		}
}

bool JTCSkimer::trkCuts(int i){
		// return 1 for pass the cuts, otherwise return 0;
		if(trkPt->at(i)<=trkPtcut || trkPt->at(i) > 400) return false;
		if(abs(trkEta->at(i))>=trketamaxcut) return false;
		if(!highPurity->at(i)) return false;
		if(fabs(trkPtError->at(i)/trkPt->at(i))>0.3) return false;
		if(TMath::Abs(trkDz->at(i)/trkDzError->at(i))>=3.0 ||
						TMath::Abs(trkDxy->at(i)/trkDxyError->at(i))>=3.0) return false;
		//if(!ispp && (float)trkChi2/(float)trkNdof/(float)trkNlayer > 0.15) return false;
		//        //if(!ispp && trkNHit<11 && trkPt > 0.7) return false;
		if((float)trkChi2->at(i)/(float)trkNdof->at(i)/(float)trkNlayer->at(i) > 0.15) return false;
		if(trkNHit->at(i)<11 && trkPt->at(i) > 0.7) return false;

		float Et = (pfHcal->at(i)+pfEcal->at(i))/TMath::CosH(trkEta->at(i));
		if(!(trkPt->at(i)<20 || Et > 0.5*trkPt->at(i))) return false;
		return 1;
}

bool JTCSkimer::eventCuts(){
		if(TMath::Abs(vz)>vzcut) return 0;
		if(isMC && pthat<80) return 0;
		if(!isHI){
				if(!collisionEventSelection || !HBHEFilter) return 0;
		}
		return 1;
}

double JTCSkimer::findDr(double eta1, double phi1, double eta2, double phi2){
		double dphi = phi1 - phi2;
		while (dphi > M_PI) dphi-= 2*M_PI;
		while (dphi <= -M_PI) dphi += 2*M_PI;
		return sqrt(pow(dphi,2)+pow(eta1-eta2,2));
}

bool JTCSkimer::recoJetCut(int j){
		if(jteta->at(j)>jetetacut) return false;
		if(jtpt->at(j)<jetptcut) return false;
		return true;
}

bool JTCSkimer::btagger(int j){
		if(discr_csv->at(j) > 0.9) return true;
		else return false;
}

void JTCSkimer::trackLoop(){
		if(doGenTrack){
				for(int jtk=0; jtk< pt->size(); ++jtk){
						if(chg->at(jtk) == 0) continue;
						if(TMath::Abs(eta->at(jtk)) > trketamaxcut) continue;
						if(pt->at(jtk) < trkPtcut) continue;
						selectedTrack_pt =pt ->at(jtk);
						selectedTrack_eta=eta->at(jtk);
						selectedTrack_phi=phi->at(jtk);
						fillHist();
				}
		}
		else{
				for(int jtk=0; jtk< trkPt->size(); ++jtk){
						if(!trkCuts(jtk)) continue;
						selectedTrack_pt =trkPt ->at(jtk);
						selectedTrack_eta=trkEta->at(jtk);
						selectedTrack_phi=trkPhi->at(jtk);
						trkCorr = ppTrkCorrection->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, 1);
						fillHist();
				}
		}
}

void JTCSkimer::loop(){
		configHist();
		Long64_t nentries = t->GetEntriesFast();
		if(doGenJet )	std::cout<<"staring Gen-jet loop..."<<std::endl;
		else std::cout<<"staring Reco-jet loop..."<<std::endl;
		for(Long64_t jentry=0; jentry<nentries; ++jentry){
				if(jentry % 10000 ==0 ) std::cout<<"processing "<<jentry<<"th event..."<<std::endl;
				t->GetEntry(jentry);
				if(!eventCuts()) continue;
				eventWeight();
				fillEventInfo();
				if(doGenJet){
						for(int jjet = 0; jjet<gen_jtpt->size(); ++jjet){
								foundTrueBJet = 0, foundTaggedBJet=0, foundInclJet=0;
								if(gen_jtpt->at(jjet) < 120) continue;
								if(TMath::Abs(gen_jteta->at(jjet)) > jetetacut) continue;
								foundInclJet = 1; 
								double dr = 0.4;
								int recoIndex = -1;
								for(int recoj = 0; recoj<jtpt->size(); ++recoj){
										if(findDr(jteta->at(recoj), jtphi->at(recoj), gen_jteta->at(jjet), gen_jtphi->at(jjet))>dr)continue;
										dr =findDr(jteta->at(recoj), jtphi->at(recoj), gen_jteta->at(jjet), gen_jtphi->at(jjet)); 
										recoIndex = recoj;
								}
								if( recoIndex < 0) continue;
								if(btagger(recoIndex)) foundTaggedBJet = 1;
								if(fabs(flavorForB->at(recoIndex))== 5) foundTrueBJet =1;
								selectedJt_pt = gen_jtpt->at(jjet);			
								selectedJt_eta = gen_jteta->at(jjet);			
								selectedJt_phi = gen_jtphi->at(jjet);	
								fillAllJetInfo();
								trackLoop();
						}
				}
				else{
						for(int jjet = 0; jjet<jtpt->size(); ++jjet){
								foundTrueBJet = 0, foundTaggedBJet=0, foundInclJet=0;
								if(!recoJetCut(jjet)) continue;
								foundInclJet = 1; 
								if(btagger(jjet)) foundTaggedBJet = 1;
								if(fabs(flavorForB->at(jjet))== 5) foundTrueBJet =1;
								selectedJt_pt = jtpt->at(jjet);			
								selectedJt_eta = jteta->at(jjet);			
								selectedJt_phi = jtphi->at(jjet);	
								fillAllJetInfo();
								trackLoop();
						}
				}
				//end of reco jet loop
		}
}
