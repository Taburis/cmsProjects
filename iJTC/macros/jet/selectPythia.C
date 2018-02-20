
#include "pythia_jet_tree.h"
double findDR(double eta1, double phi1, double eta2, double phi2){
		phi1 = TMath::ACos(TMath::Cos(phi1));
		// phi2 = TMath::ACos(TMath::Cos(phi2));
		// return sqrt(pow(phi1-phi2,2)+pow(eta1-eta2,2));
		double dphi = fabs(phi1-phi2);
		if( dphi> TMath::Pi()) dphi = 2*TMath::Pi()-dphi;
		return sqrt(pow(dphi,2)+pow(eta1-eta2,2));
}

void selectPythia(){

		auto fit_vz = new TF1("newppVz","gaus",-15,15);
		fit_vz->SetParameter(0,1.10477);
		fit_vz->SetParameter(1,2.52738);
		fit_vz->SetParameter(2,1.30296e1);
//		double ppPthatEntries[9] = {0,0,272902,377559,467823,447683,259111,234347,50942};
		double ppPthatEntries[9] = {0,0,271007,374878,464554,444518,257311,232698,50612};

		auto *t = new pythia_jet_tree();
		t->fChain->SetBranchStatus("*",0);
		t->fChain->SetBranchStatus("pthat",1);
		t->fChain->SetBranchStatus("vz",1);
		t->fChain->SetBranchStatus("genpt",1);
		t->fChain->SetBranchStatus("geneta",1);
		t->fChain->SetBranchStatus("genphi",1);
		t->fChain->SetBranchStatus("calo_refpt",1);
		t->fChain->SetBranchStatus("calo_corrpt",1);
		t->fChain->SetBranchStatus("calo_jteta",1);
		t->fChain->SetBranchStatus("calo_jtphi",1);
		t->fChain->SetBranchStatus("calo_jtpt",1);
		t->fChain->SetBranchStatus("calo_trackMax",1);
		t->fChain->SetBranchStatus("calo_discr_csvV1",1);
		t->fChain->SetBranchStatus("calo_refparton_flavorForB",1);

		TH1D* hvz = new TH1D("hvz", "",100, -15, 15);

		TH1D* hvz_weighted = new TH1D("hvz_weighted", "",100, -15, 15);
		TH1D* hpthat_weighted = new TH1D("hpthat_weighted", "",200, 80, 380 );

		const int ncsv=12;
		TH1D* htagged_b[ncsv][2], *htrue_b[2], *htagged_true_b[ncsv][2];
		TH1D* hincl[2];
		TH1D* heta[2], *hphi[2];
		TH1I * hcounter[ncsv][2];

		float csvBin[12] = {0.1, 0.15, 0.3, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99};
		const Double_t newbin [21] = {100, 120, 136, 152, 168, 184, 200, 216, 232, 248, 264, 280, 296, 312, 328, 344, 360,        
				                380, 400, 432, 500};

		for(int i=0; i<2; ++i){
				htrue_b        [i] = new TH1D(Form("trueB_%d",i) , "", 100, 100, 500);
				hincl        [i] = new TH1D(Form("incl_%d",i) , "", 20, newbin);
				for(int j=0; j<ncsv; ++j){
						hcounter       [j][i] = new TH1I(Form("counter%d_%d",j,i), "", 6, 0, 6);
						htagged_b      [j][i] = new TH1D(Form("taggedB_%d_%d",j,i), "", 100, 100, 500);
						htagged_true_b [j][i] = new TH1D(Form("tagged_trueB_%d_%d",j,i) , "", 100, 100, 500);
						//heta   [i] = new TH1D(Form("heta_%d",i), "", 500, -2.5, 2.5);
						//hphi   [i] = new TH1D(Form("hphi_%d",i), "", 720, -TMath::Pi(), TMath::Pi());
				}
		}


		TFile *wf = new TFile("pythiaSpectra2.root" , "recreate");
		Long64_t ne = t->fChain->GetEntriesFast();
		float wvz=1, wpthat=1, w;
		Double_t pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
		TH1D* hpthat = new TH1D("hpthat", "",9, pthatbins );
		double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};

		//		ne=10000;
		TF1 *jes_cor_f = new TF1("jes_cor_f", "[0]*exp(-pow(log([2]*x-[2]*[1]),2)/pow([3],2)/2)/[2]/(x-[1])/[3]+[4]+[5]*(x-100)",100, 600);
		jes_cor_f->SetParameters(0.1, 100, 0.005, 1.12, 0.9, 0.000023);
		Long64_t nmulti=0, nunmatch=0, jetcounter=0;
		for(Long64_t jentry = 0; jentry< ne; ++jentry){
				t->GetEntry(jentry);
				//if(t->pthat<80) continue;
				wvz=1.0/fit_vz->Eval(t->vz);
				int ibin =0;
				while(t->pthat>pthatbins[ibin+1]) ibin++;
				wpthat= (xsecs[ibin]-xsecs[ibin+1])/ppPthatEntries[ibin];
				w= wvz*wpthat;
				if(jentry %100000 == 0) cout<<"processed "<<jentry<<" events..."<<endl;
				hvz->Fill(t->vz); hpthat->Fill(t->pthat);
				hvz_weighted->Fill(t->vz,w); hpthat_weighted->Fill(t->pthat,w);
				int i=0;
				for(int j=0;j < t->genpt->size(); ++j){
						if(t->genpt->at(j)<120) continue;
						if(fabs(t->geneta->at(j))>1.6) continue;
						hincl[i]->Fill(t->genpt->at(j),w);
				}
				for(int j=0;j < t->calo_corrpt->size(); ++j){
						if((t->calo_corrpt->at(j))<60) continue;
						//						if((t->calo_trackMax->at(j))/t->calo_corrpt->at(j)<=0.01) continue;

						int imatch = -1;
						double dr = 4;
						int nmatch = 0;
						int foundGen = 0;
						for(unsigned int igen=0; igen<t->genpt->size(); igen++){
							//	if( fabs(t->genpt->at(igen)-t->calo_refpt->at(j))< 1) {
										foundGen = 1;
										double dr1 = findDR(t->geneta->at(igen), t->genphi->at(igen), t->calo_jteta->at(j),t->calo_jtphi->at(j));
										if(dr> dr1){
												dr = dr1;
												imatch = igen;
										}
							//	}
								//if( fabs(t->genpt->at(igen)-t->calo_refpt->at(j))< 1) {}
						}
						//jetcounter =jetcounter+foundGen;
						if(nmatch>1) {
								nmulti = nmulti+nmatch;
								//	cout<<"final dr= "<<dr<<endl;
						}
						if(imatch<0 ) {
//								nunmatch++;
								continue;
						}
						jetcounter++;
						float matchedPt = t->genpt ->at(imatch);
						float matchedEta= t->geneta->at(imatch);
						float matchedPhi= t->genphi->at(imatch);

						if(fabs(matchedPt-t->calo_refpt->at(j))>.1){
								nunmatch++;
						}

						/*
						   float matchedPt = t->calo_jtpt ->at(j);
						   float matchedEta= t->calo_jteta->at(j);
						   float matchedPhi= t->calo_jtphi->at(j);
						   */
						//float matchedPt = t->calo_corrpt->at(j);
						if(matchedPt<120) continue;
						if(fabs(matchedEta)>1.6) continue;
						//if(matchedPt<80) continue;
						//if(fabs(matchedEta)>2) continue;

						for(int k=0; k<ncsv; k++){
								if(t->calo_discr_csvV1->at(j)>=csvBin[k]) {
										htagged_b[k][i]->Fill(matchedPt , w);
										hcounter [k][i]->Fill(1);
										if(fabs(t->calo_refparton_flavorForB->at(j)) ==5){
												hcounter [k][i]->Fill(3);
												htagged_true_b [k][i]->Fill(matchedPt , w);
										}
								}
						}
						/*
						*/
						if(fabs(t->calo_refparton_flavorForB->at(j)) !=5) continue;
						hcounter[0][i]->Fill(5);
						//cout<<"filling"<<endl;
						htrue_b [i]->Fill(matchedPt , w);
						//heta   [i]->Fill(matchedEta, w);
						//hphi   [i]->Fill(matchedPhi, w);
				}

		}
		cout<<"Total jets: "<<jetcounter<<endl;
		cout<<"Unmatched jets: "<<nunmatch<<" ("<<float(100*nunmatch)/float(jetcounter)<<")%"<<endl;
		cout<<"Multi jets: "<<nmulti<<" ("<<float(100*nmulti)/float(jetcounter)<<")%"<<endl;
		for(int i=0; i<2; ++i){
				htrue_b[i]->Write();
				hincl  [i]->Write();
				for(int j=0; j<ncsv; j++){
						htagged_b     [j][i]->Write();
						htagged_true_b[j][i]->Write();
						hcounter      [j][i]->Write();
						//heta   [i]->Write();
						//hphi   [i]->Write();
				}
		}
		hpthat->Write();
		hvz_weighted->Write();
		hpthat_weighted->Write();
}
