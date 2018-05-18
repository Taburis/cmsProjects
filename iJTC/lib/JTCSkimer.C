
#include "JTCSkimer.h"
//#include "../corrTable/bJTCTrkCorr_pp.h"
#include "TF1.h"

float trketamaxcut = 2.;
float trkPtcut = 1.0;
float jetetacut = 1.6;
float jetptcut = 120;
float vzcut = 15;

double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
double ppPthatEntries[9] = {0,0,272902,377559,467823,447683,259111,234347,50942};

void JTCSkimer::configMixing(){
		std::cout<<"configuring the trees for mixing..."<<std::endl;
		ncent_mix = isHI ? 40 : 1;
		nvz_mix = 30;
		//saving the mixing event index for different hibin and vz bin
		// here are 30 vz bins with equal bin width 1 and the index labeled based on: floor(vz+15)
		// 40 hiBin bins with equal bin width 5 and the index labeled based on: floor(float(hiBin)/5)
		// (the vz = 15 and the hiBin = 200 has to be dropped as it would lead to crash sinec index for hiBin=200 is 40)
		mixTable = new std::vector<unsigned int>*[nvz_mix*ncent_mix]; 
		for(int i=0; i<nvz_mix; ++i){
				for(int j=0; j<ncent_mix; ++j){
						mixTable[i+nvz_mix*j] = new std::vector<unsigned int >();
				}
		}

		mt->init(isMC, isHI); //build the mixing tree after added the mixing files
		Long64_t nentries = mt->t->GetEntries();
		std::cout<<"total events in mixing = "<<nentries<<std::endl;
		for(Long64_t jevt = 0; jevt<nentries; ++jevt){
				mt->t->GetEntry(jevt);
				if(jevt%10000 == 0 ) std::cout<<"scan "<<jevt<<" events in mixing tree..."<<std::endl;
				//				cout<<mt->vz<<endl;
				//need to implement the PbPb case event cut
				if(!(mt->eventCuts())) continue; 
				//if(mt->vz>=15 || mt->vz<=-15) continue;
				int ivz = floor(mt->vz+15);
				int ihibin = isHI ? floor(float(mt->hiBin)/5) : 0;
				//				cout<<"vz = "<<mt->vz<<", ivz = "<<ivz<<endl;
				if(mixTable[ivz+nvz_mix*ihibin]->size()< 100)mixTable[ivz+nvz_mix*ihibin]->push_back(jevt);
		}
		std::cout<<"vz statitcs: "<<endl;
		for(int i= 0; i<nvz_mix; ++i){
				for(int j= 0; j<ncent_mix; ++j){
						std::cout<<-15+i<<" < vz < "<<-14+i<<"; "
								<<j*5<<" < cent < "<<j*5+5<<": "<<mixTable[i+nvz_mix*j]->size()<<std::endl;
				}
		}
}

void JTCSkimer::eventWeight(float pthat, float vz, float hibin){
		//need the PbPb weight 
		weights=1;
		if(!isHI){
				int ibin=0;
				while(pthat>pthatbins[ibin+1]) ibin++;
				weights=(xsecs[ibin]-xsecs[ibin+1])/ppPthatEntries[ibin];
				weights=weights*(1.0/ppVz->Eval(vz));
		}
}

int JTCSkimer::trkCuts(float trkpt, float trkpterror, float trketa, 
				float trkchi2, int highpurity, int trknhit, float trkndof, float trknlayer,
				float pfhcal, float pfecal){
		// return 1 for pass the cuts, otherwise return 0;

		if(trkpt<=trkPtcut || trkpt > 400) return 2;
		if(TMath::Abs(trketa) >=trketamaxcut) return 3;
		if(!highpurity) return 4;
		if(fabs(trkpterror/trkpt)>0.3) return 5;
		//		if(TMath::Abs(trkDz->at(i)/trkDzError->at(i))>=3.0 ||
		//						TMath::Abs(trkDxy->at(i)/trkDxyError->at(i))>=3.0) return false;
		if(trknhit<11 ) return 6;
		if(trkchi2/trkndof/trknlayer > 0.15) return 7;

		float Et = (pfhcal+pfecal)/TMath::CosH(trketa);
		if(!(trkpt<20 || Et > 0.5*trkpt)) return 8;
		return 0;
}

bool JTCSkimer::eventCuts(){
		//need to implement the PbPb case event cut
		if(TMath::Abs(vz)>=vzcut) return 0;
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

float JTCSkimer::findDrmin(float eta, float phi, vector<float> *jetpt, vector<float> *jeteta, vector<float> *jetphi){
		float dr = 999;
		for(int i=0 ;i< jetpt->size(); ++i){
				if(jetpt->at(i) < 80) continue;
				float dr0 = findDr(eta, phi, jeteta->at(i), jetphi->at(i));
				if(dr> dr0) dr=dr0;
		}
		return dr;
}

bool JTCSkimer::btagger(int j){
		if(discr_csv->at(j) > 0.9) return true;
		else return false;
}

void JTCSkimer::trackLoop(){
		if(doGenTrack){
				for(size_t jtk=0; jtk< pt->size(); ++jtk){
						if(chg->at(jtk) == 0) continue;
						if(TMath::Abs(eta->at(jtk)) > trketamaxcut) continue;
						if(pt->at(jtk) < trkPtcut) continue;
						selectedTrack_pt =pt ->at(jtk);
						selectedTrack_eta=eta->at(jtk);
						selectedTrack_phi=phi->at(jtk);
						if(selectedTrack_pt>399) selectedTrack_pt = 399; // for setting the maximum bin content
						fillHist();
				}
		}
		else{
				for(size_t jtk=0; jtk< trkPt->size(); ++jtk){
					  int code = trkCuts(trkPt->at(jtk), trkPtError->at(jtk), trkEta->at(jtk),
												(float)trkChi2->at(jtk), highPurity->at(jtk), trkNHit->at(jtk),
												(float)trkNdof->at(jtk), (float)trkNlayer->at(jtk),
												pfHcal->at(jtk), pfEcal->at(jtk));
					//	cout<<"exit code = "<<code<<endl;
						if(code ) continue;
						float dcaCut = 0;
						if(TMath::Abs(trkDz->at(jtk)/trkDzError->at(jtk))>=3.0 ||
										TMath::Abs(trkDxy->at(jtk)/trkDxyError->at(jtk))>=3.0) dcaCut = 1;
						selectedTrack_pt =trkPt ->at(jtk);
						selectedTrack_eta=trkEta->at(jtk);
						selectedTrack_phi=trkPhi->at(jtk);
						float drmin = findDrmin(selectedTrack_eta, selectedTrack_phi, jtpt, jteta, jtphi);
						trkCorr = ppTrkCorrection->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, drmin);
						//trkCorr = ppTrkCorrection->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, 1);
						fillHist();
				}
		}
}

void JTCSkimer::loop(){
		configHist();
		Long64_t nentries = t->GetEntriesFast();
		if(doGenJet )	std::cout<<"staring Gen-jet ";
		else std::cout<<"staring Reco-jet ";
		if(doGenTrack )	std::cout<<"and Gen-track correlation... "<<std::endl;
		else std::cout<<"and Reco-track correlation... "<<std::endl;
		for(Long64_t jentry=0; jentry<nentries; ++jentry){
				if(jentry % 1000 ==0 ) std::cout<<"processing "<<jentry<<"th event..."<<std::endl;
				t->GetEntry(jentry);
				if(!eventCuts()) continue;
				eventWeight(pthat, vz, hiBin);
				fillEventInfo();
				if(doGenJet){
						for(size_t jjet = 0; jjet<gen_jtpt->size(); ++jjet){
								foundTrueBJet = 0, foundTaggedBJet=0, foundInclJet=0;
								if(gen_jtpt->at(jjet) < 120) continue;
								if(TMath::Abs(gen_jteta->at(jjet)) > jetetacut) continue;
								foundInclJet = 1; 
								double dr = 0.4;
								int recoIndex = -1;
								for(size_t recoj = 0; recoj<jtpt->size(); ++recoj){
										float dr0 =findDr(jteta->at(recoj), jtphi->at(recoj), gen_jteta->at(jjet), gen_jtphi->at(jjet)); 
										if(dr0>dr)continue;
										dr = dr0;
										recoIndex = recoj;
								}
								if( recoIndex >= 0){
										if(btagger(recoIndex)) foundTaggedBJet = 1;
										if(fabs(flavorForB->at(recoIndex))== 5) foundTrueBJet =1;
								}
								selectedJt_pt = gen_jtpt->at(jjet);			
								selectedJt_eta = gen_jteta->at(jjet);			
								selectedJt_phi = gen_jtphi->at(jjet);	
								fillAllJetInfo();
								trackLoop();
								if(!do_mixing) continue;
								mixingSection(jentry);
						}
				}
				else{
						for(size_t jjet = 0; jjet<jtpt->size(); ++jjet){
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

								if(!do_mixing) continue;
								mixingSection(jentry);
						}
				}
				//end of reco jet loop
		}
}

void JTCSkimer::mixingSection(Long64_t voidIndex){
		int kevt = 0; 
		int vzIndex = floor(vz+15);
		int centIndex = isHI ? floor(float(hiBin)/5) : 0;
		for(int kmix = 0; kmix<nPerTrigger; ++kmix){
				//				cout<<kmix<<endl;
				if(kevt == mixTable[vzIndex+centIndex*ncent_mix]->size()) kevt = 0;
//				cout<<"vz = "<<vz<<", index = "<<vzIndex+centIndex*ncent_mix<<endl;
//				cout<<"index = "<<mixTable[vzIndex+centIndex*ncent_mix]->at(kevt)<<endl;
				Long64_t index = mixTable[vzIndex+centIndex*ncent_mix]->at(kevt);
				if(index == voidIndex) continue; // to avoid the auto correlation in any case
				(mt->t)->GetEntry(index);
				trackLoopForMixing();
				kevt++;
		}
}

void JTCSkimer::trackLoopForMixing(){
		if(doGenTrack){
				for(size_t jtk=0; jtk< (mt->pt)->size(); ++jtk){
						if((mt->chg)->at(jtk) == 0) continue;
						if(TMath::Abs((mt->eta)->at(jtk)) > trketamaxcut) continue;
						if((mt->pt)->at(jtk) < trkPtcut) continue;
						selectedTrack_pt =(mt->pt )->at(jtk);
						selectedTrack_eta=(mt->eta)->at(jtk);
						selectedTrack_phi=(mt->phi)->at(jtk);
						if(selectedTrack_pt>399) selectedTrack_pt = 399; // for setting the maximum bin content
						fillHist(1);
				}
		}	else{
				for(size_t jtk=0; jtk< (mt->trkPt)->size(); ++jtk){
						//						cout<<jtk<<endl;
//						cout<<"=========================================="<<endl;
//						cout<<"ndof  ="<<(mt->trkNdof)->at(jtk)<<endl;
//						cout<<"nlayer="<<(mt->trkNlayer)->at(jtk)<<endl;
//						cout<<"chi2  ="<<(mt->trkChi2)->at(jtk)<<endl;
//						cout<<"here:"<<(mt->trkChi2)->at(jtk)/(float)(mt->trkNdof)->at(jtk)/(float)(mt->trkNlayer)->at(jtk)<<endl;
						int code = trkCuts((mt->trkPt)->at(jtk), (mt->trkPtError)->at(jtk), (mt->trkEta)->at(jtk),
												(float)(mt->trkChi2)->at(jtk), (mt->highPurity)->at(jtk), (mt->trkNHit)->at(jtk),
												(float)(mt->trkNdof)->at(jtk), (float)(mt->trkNlayer)->at(jtk),
												(mt->pfHcal)->at(jtk), (mt->pfEcal)->at(jtk));
//						cout<<"code = "<<code<<endl;
						if(code) continue;
//						float dcaCut = 0;
//						if(TMath::Abs((mt->trkDz)->at(jtk)/(mt->trkDzError)->at(jtk))>=3.0 ||
//										TMath::Abs((mt->trkDxy)->at(jtk)/(mt->trkDxyError)->at(jtk))>=3.0) dcaCut = 1;
						selectedTrack_pt =(mt->trkPt )->at(jtk);
						selectedTrack_eta=(mt->trkEta)->at(jtk);
						selectedTrack_phi=(mt->trkPhi)->at(jtk);
						float drmin = findDrmin(selectedTrack_eta, selectedTrack_phi, mt->jtpt, mt->jteta, mt->jtphi);
						trkCorr = ppTrkCorrection->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, drmin);
						fillHist(1);
				}
		}

}
