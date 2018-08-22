
#include "JTCSkimer.h"
#include "TF1.h"
//#include "../corrTable/bJTCTrkCorr_pp.h"

#include "TrkCorr_July22_Iterative_pp_eta2p4/getTrkCorr.h"
TrkCorr* trkCorr_OBJ = new TrkCorr("TrkCorr_July22_Iterative_pp_eta2p4/"); //must have the '/' after the subdirectory name

double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
double ppPthatEntries[9] = {0,0,272902,377559,467823,447683,259111,234347,50942};
double pthatEntries[9] = {0, 0, 0, 2.56444e+06, 2.84656e+06, 2.67143e+06, 2.88286e+06, 779403, 168046}; //cymbal tune

double pp_bMC_pthat_entries[4] = {352315., 399215., 645016., 595293.};
double xsecs_bMC[5] = {4.412E-04, 1.511E-04, 6.147E-05, 1.018E-05, 0};
double pthatbins_bMC[5] = {80, 100, 120, 170, 9999};


float JTCSkimer::eff_weight(float pt){
	int n = efficiency_hist->FindBin(pt);
	return 1.0/efficiency_hist->GetBinContent(n);
}

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

	mt->init(isMC, isHI, doKurtSkim); //build the mixing tree after added the mixing files
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
			if(mixTable[i+nvz_mix*j]->size()<2) std::cout<<-15+i<<" < vz < "<<-14+i<<"; "
				<<j*5<<" < cent < "<<j*5+5<<": "<<mixTable[i+nvz_mix*j]->size()<<std::endl;
		}
	}
}

void JTCSkimer::eventWeight(float pthat, float vz, float hibin){
	//need the PbPb weight 
	weights=1;
	if(!isHI){
		int ibin=0;
		if(isbMC){
			while(pthat>pthatbins_bMC[ibin+1]) ibin++;
			weights=(xsecs_bMC[ibin]-xsecs_bMC[ibin+1])/pp_bMC_pthat_entries[ibin];
			weights=weights*(pp_vz_bMC->Eval(vz));
		}else {
			while(pthat>pthatbins[ibin+1]) ibin++;
			weights=(xsecs[ibin]-xsecs[ibin+1])/ppPthatEntries[ibin];
			weights=weights*(1.0/ppVz->Eval(vz));
		}
	} else {
		int ibin=0;
		while(pthat>pthatbins[ibin+1]) ibin++;
		weights=(xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];
		float wcen = (hiBin<194) ? fit_cen->Eval(1.*hiBin) : 1.;
		weights=weights*fit_vz_cymbal->Eval(vz)*wcen;
	}
// this section for the pre-approval homework
}

int JTCSkimer::trkCuts(float trkpt, float trkpterror, float trketa, 
		float trkchi2, int highpurity, int trknhit, float trkndof, float trknlayer,
		float pfhcal, float pfecal, float trkdz, float trkdxy, float trkdzerror, float trkdxyerror){
	// return 0 for pass the cuts;
	/* standard cuts
	*/
	if(trkpt<=trkPtcut || trkpt > 400) return 2;
	if(TMath::Abs(trketa) >=trketamaxcut) return 3;
	if(!highpurity) return 4;
	float Et = (pfhcal+pfecal)/TMath::CosH(trketa);
	if(!(trkpt<20 || Et > 0.5*trkpt)) return 5;

	if(doDCATkCut) {
		if(TMath::Abs(trkdz/trkdzerror)>=3.0 ||
				TMath::Abs(trkdxy/trkdxyerror)>=3.0) return 1;
	} else {
		if(trknhit<11 ) return 6;
		if(trkchi2/trkndof/trknlayer > 0.15) return 7;
	}

	/* varing cuts
	//	if(trkpt<=trkPtcut || trkpt > 400) return 2;
	//	if(TMath::Abs(trketa) >=trketamaxcut) return 3;
	//	if(!highpurity) return 4;
	//	if(fabs(trkpterror/trkpt)>0.5) return 5;
	//	if(TMath::Abs(trkdz/trkdzerror)>=3.0 ||
	//			TMath::Abs(trkdxy/trkdxyerror)>=3.0) return false;
	//	if(trknhit<11 ) return 6;
	//	if(trkchi2/trkndof/trknlayer > 0.15) return 7;

	//float Et = (pfhcal+pfecal)/TMath::CosH(trketa);
	//if(!(trkpt<20 || Et > 0.5*trkpt)) return 8;
	*/
	return 0;
}

bool JTCSkimer::eventCuts(){
	//need to implement the PbPb case event cut
	if(TMath::Abs(vz)>=vzcut) return 0;
	if(isMC ) if( pthat<80) return 0;
	if(!isHI){
		if(!pvFilter || !HBHEFilter) return 0;
	} else {
		if(!phfCoincFilter || !collisionEventSelection || !pvFilter || !HBHEFilter) return 0;
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
	int hibin;
	//ncs corr
	if(TMath::Abs(jteta->at(j))>jetetacut) return false;
	if(!isHI) hibin = 0;
	//	if((ncscorr->getCorrection(!isHI, ncs->at(j), hibin, jtpt->at(j), jteta->at(j)))<jetptcut) return false;
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
	if(discr_csv->at(j) > csvCut) return true;
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
			if(isHI) gen_sube = sube->at(jtk);
			fillHist();
		}
	}
	else{
		for(size_t jtk=0; jtk< trkPt->size(); ++jtk){
			int code = trkCuts(trkPt->at(jtk), trkPtError->at(jtk), trkEta->at(jtk),
					(float)trkChi2->at(jtk), highPurity->at(jtk), trkNHit->at(jtk),
					(float)trkNdof->at(jtk), (float)trkNlayer->at(jtk),
					pfHcal->at(jtk), pfEcal->at(jtk), 
					trkDz->at(jtk),trkDxy->at(jtk),trkDzError->at(jtk),trkDxyError->at(jtk));
			//	cout<<"exit code = "<<code<<endl;
			if(code) continue;
			dcaCut = 0;
			if(TMath::Abs(trkDz->at(jtk)/trkDzError->at(jtk))>=3.0 ||
					TMath::Abs(trkDxy->at(jtk)/trkDxyError->at(jtk))>=3.0) dcaCut = 1;
			selectedTrack_pt =trkPt ->at(jtk);
			selectedTrack_eta=trkEta->at(jtk);
			selectedTrack_phi=trkPhi->at(jtk);
			//			float rmin = findDrmin(selectedTrack_eta, selectedTrack_phi, jtpt, jteta, jtphi);
			float rmin = 999;
			//	if(!isHI) trkCorr = ppTrkCorrection->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, rmin);
			trkCorr = trkCorr_OBJ->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, 1, rmin);
			//trkCorr = ppNominalTrkCorrection->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, 1, rmin);
			//trkCorr = ppTrkCorrection->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, 1);
			fillHist();
		}
	}
}


void JTCSkimer::loop(){
	configHist();
	std::cout<<"working config---------------------------------------"<<endl;
	if(doCaloJet) std::cout<<"	CaloJet"<<endl;
	else std::cout<<"	PFJet"<<endl;
	std::cout<<"	CSV cut at "<<csvCut<<endl;
	if(doDCATkCut )	std::cout<<"	with DCA tracking cut "<<endl;
	else std::cout<<"	no DCA tracking cut ";
	if(doWTAaxis )	std::cout<<"	using WTA axis "<<endl;
	else std::cout<<"	using anti-kT jet axis ";
	if(doGenJet )	std::cout<<"	Gen-jet ";
	else std::cout<<"	Reco-jet ";
	if(doGenTrack )	std::cout<<"and Gen-track correlation... "<<std::endl;
	else std::cout<<"and Reco-track correlation... "<<std::endl;
	std::cout<<"----------------end of config------------------------"<<endl;

	if(do_mixing) configMixing();
	Long64_t nentries = t->GetEntriesFast();
	for(Long64_t jentry=0; jentry<nentries; ++jentry){
		if(jentry % 1000 ==0 ) std::cout<<"processing "<<jentry<<"th event..."<<std::endl;
		t->GetEntry(jentry);
		if(!eventCuts()) continue;
		if(isMC) eventWeight(pthat, vz, hiBin);
		else weights = 1;
		fillEventInfo();
		if(doGenJet){
			for(size_t jjet = 0; jjet<gen_jtpt->size(); ++jjet){
				foundTrueBJet = 0, foundTaggedBJet=0, foundInclJet=0;
				if(TMath::Abs(gen_jteta->at(jjet)) > jetetacut) continue;
				float smear_pt=0;
//smearing ==========================
				if(doSmear ){
					smear_pt = rand3.Gaus(0, jer_sigma);
				}
				if(gen_jtpt->at(jjet)+smear_pt < jetptcut) continue;
//-----------------------------------
				foundInclJet = 1; 
				double dr = 0.4;
				int recoIndex = -1;
				for(size_t recoj = 0; recoj<jtpt->size(); ++recoj){
					if(jtpt->at(recoj) < 80) continue;
					float dr0 =findDr(jteta->at(recoj), jtphi->at(recoj), gen_jteta->at(jjet), gen_jtphi->at(jjet)); 
					if(dr0>dr)continue;
					dr = dr0;
					recoIndex = recoj;
				}
				if( recoIndex >= 0){
					if(btagger(recoIndex)) foundTaggedBJet = 1;
					if(fabs(flavorForB->at(recoIndex))== 5) foundTrueBJet =1;
				}
				selectedJt_pt = gen_jtpt->at(jjet) + smear_pt;
				selectedJt_eta= gen_jteta->at(jjet);			
				selectedJt_phi= gen_jtphi->at(jjet);	
				if(doWTAaxis) findWTAaxis();
				fillAllJetInfo();
				trackLoop();
				if(!do_mixing) continue;
				mixingSection(jentry);
			}
		}
		else{
			//cout<<"size "<<endl;
			//cout<<jtpt->size()<<endl;
			//cout<<ncs->size()<<endl;
			for(size_t jjet = 0; jjet<jtpt->size(); ++jjet){
				foundTrueBJet = 0, foundTaggedBJet=0, foundInclJet=0;
				if(!recoJetCut(jjet)) continue;
				//matching to gen jet
				if(doJetAxisResolution){
					int genIndex = -1;
					double dr = 0.4;
					for(size_t genj = 0; genj<gen_jtpt->size(); ++genj){
						if(gen_jtpt->at(genj) < 80) continue;
						float dr0 =findDr(jteta->at(jjet), jtphi->at(jjet), gen_jteta->at(genj), gen_jtphi->at(genj)); 
						if(dr0>dr)continue;
						dr = dr0;
						genIndex = genj;
					}
					if(genIndex !=-1){
						genJetAxis_pt  = gen_jtpt ->at(genIndex);			
						genJetAxis_eta = gen_jteta->at(genIndex);			
						genJetAxis_phi = gen_jtphi->at(genIndex);	
					} else{
						continue;
					}
				}
				//matching end
				foundInclJet = 1; 
				if(btagger(jjet)) foundTaggedBJet = 1;
				if(isMC) if(fabs(flavorForB->at(jjet))== 5) foundTrueBJet =1;
				if(doJetAxisResolution){
					selectedJt_pt  = genJetAxis_pt;     		
					selectedJt_eta = genJetAxis_eta;		
					selectedJt_phi = genJetAxis_phi;
				}else {
					selectedJt_pt = jtpt->at(jjet);			
					selectedJt_eta = jteta->at(jjet);			
					selectedJt_phi = jtphi->at(jjet);	
				}
				if(doWTAaxis) findWTAaxis();
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
	if(mixTable[vzIndex+centIndex*nvz_mix]->size()<2){
		if(vzIndex == 29) vzIndex = 28;// shift mixing vz from 14, 15 to 13, 14; 
		else if(vzIndex == 0) vzIndex = 1;// shift mixing vz from -14,-15 to -13, -14; 
		else if(centIndex>0) centIndex=centIndex-1;  // shift the most prepheral to next prepheral
	}
	if(mixTable[vzIndex+centIndex*nvz_mix]->size()==0) return; 
	for(int kmix = 0; kmix<nPerTrigger; ++kmix){
		//				cout<<kmix<<endl;
		if(kevt == mixTable[vzIndex+centIndex*nvz_mix]->size()) kevt = 0;
		//				cout<<"vz = "<<vz<<", index = "<<vzIndex+centIndex*ncent_mix<<endl;
		//				cout<<"index = "<<mixTable[vzIndex+centIndex*ncent_mix]->at(kevt)<<endl;
		Long64_t index = mixTable[vzIndex+centIndex*nvz_mix]->at(kevt);
		if(index == voidIndex) continue; // to avoid the auto correlation in any case
		(mt->t)->GetEntry(index);
		//cout<<"current vz: "<<vz<<", mix vz: "<<mt->vz<<endl;
		//if(isHI)cout<<"current hi: "<<hiBin<<", mix hi: "<<mt->hiBin<<endl;
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
			if(selectedTrack_pt>399) selectedTrack_pt = 399; //setting the maximum bin content, only for gen particle
			fillHist(1);
		}
	}	else{
		for(size_t jtk=0; jtk< (mt->trkPt)->size(); ++jtk){
			//cout<<jtk<<endl;
			//cout<<"=========================================="<<endl;
			//cout<<"ndof  ="<<(mt->trkNdof)->at(jtk)<<endl;
			//cout<<"nlayer="<<(mt->trkNlayer)->at(jtk)<<endl;
			//cout<<"chi2  ="<<(mt->trkChi2)->at(jtk)<<endl;
			//cout<<"here:"<<(mt->trkChi2)->at(jtk)/(float)(mt->trkNdof)->at(jtk)/(float)(mt->trkNlayer)->at(jtk)<<endl;
			int code = trkCuts((mt->trkPt)->at(jtk), (mt->trkPtError)->at(jtk), (mt->trkEta)->at(jtk),
					(float)(mt->trkChi2)->at(jtk), (mt->highPurity)->at(jtk), (mt->trkNHit)->at(jtk),
					(float)(mt->trkNdof)->at(jtk), (float)(mt->trkNlayer)->at(jtk),
					(mt->pfHcal)->at(jtk), (mt->pfEcal)->at(jtk),
					(mt->trkDz)->at(jtk), (mt->trkDxy)->at(jtk), (mt->trkDzError)->at(jtk), (mt->trkDxyError)->at(jtk));
			//cout<<"code = "<<code<<endl;
			if(code) continue;
			dcaCut = 0;
			if(TMath::Abs((mt->trkDz)->at(jtk)/(mt->trkDzError)->at(jtk))>=3.0 ||
					TMath::Abs((mt->trkDxy)->at(jtk)/(mt->trkDxyError)->at(jtk))>=3.0) dcaCut = 1;
			selectedTrack_pt =(mt->trkPt )->at(jtk);
			selectedTrack_eta=(mt->trkEta)->at(jtk);
			selectedTrack_phi=(mt->trkPhi)->at(jtk);
			float drmin = findDrmin(selectedTrack_eta, selectedTrack_phi, mt->jtpt, mt->jteta, mt->jtphi);
			//if(!isHI) trkCorr = ppTrkCorrection->getTrkCorr(selectedTrack_pt, selectedTrack_eta, selectedTrack_phi, drmin);
			trkCorr=1;
			fillHist(1);
		}
	}

}

void JTCSkimer::findWTAaxis(){
//finding the winner-take-all jet axis by using pick up the axis of the tracks with largest pt passing the track cuts inside the jet-cone.
//notice that here we need to take the gen particle as doGenJet 
	if(doGenJet){
		float initialPt = 0;
		for(size_t jtk=0; jtk< pt->size(); ++jtk){
			if(chg->at(jtk) == 0) continue;
			if(TMath::Abs(eta->at(jtk)) > trketamaxcut) continue;
			if(pt->at(jtk) < trkPtcut) continue;
			if(isHI) gen_sube = sube->at(jtk);
			if(pt->at(jtk)<initialPt) continue;
			if(findDr(selectedJt_eta, selectedJt_phi, eta->at(jtk), phi->at(jtk))<0.4){
				wta_eta = eta->at(jtk);
				wta_phi = phi->at(jtk);
				wta_pt  = pt->at(jtk);
				initialPt = pt->at(jtk);
			}
		}
	}
	else{
		float initialPt = 0;
		for(size_t jtk=0; jtk< trkPt->size(); ++jtk){
			int code = trkCuts(trkPt->at(jtk), trkPtError->at(jtk), trkEta->at(jtk),
					(float)trkChi2->at(jtk), highPurity->at(jtk), trkNHit->at(jtk),
					(float)trkNdof->at(jtk), (float)trkNlayer->at(jtk),
					pfHcal->at(jtk), pfEcal->at(jtk), 
					trkDz->at(jtk),trkDxy->at(jtk),trkDzError->at(jtk),trkDxyError->at(jtk));
			if(code) continue;
			if(trkPt->at(jtk)<initialPt) continue;
			if(findDr(selectedJt_eta, selectedJt_phi, trkEta->at(jtk), trkPhi->at(jtk))<0.4){
				wta_pt  = trkPt->at(jtk);
				wta_eta = trkEta->at(jtk);
				wta_phi = trkPhi->at(jtk);
				initialPt = trkPt->at(jtk);
			}
		}
	}
}
