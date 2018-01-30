
#include "jetTree.h"
double findDR(double eta1, double phi1, double eta2, double phi2){
		phi1 = TMath::ACos(TMath::Cos(phi1));
		// phi2 = TMath::ACos(TMath::Cos(phi2));
		// return sqrt(pow(phi1-phi2,2)+pow(eta1-eta2,2));
		double dphi = fabs(phi1-phi2);
		if( dphi> TMath::Pi()) dphi = 2*TMath::Pi()-dphi;
		return sqrt(pow(dphi,2)+pow(eta1-eta2,2));
}
void mcSelector(){


		jetTree *t = new jetTree();
		t->fChain->SetBranchStatus("*",0);
		t->fChain->SetBranchStatus("pthat",1);
		t->fChain->SetBranchStatus("hiBin",1);
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

		TH1D* hvz    = new TH1D("hvz", "",100, -15, 15);
		TH1D* hcen   = new TH1D("hcen", "",100, 0, 200);
		TH1D* hpthat = new TH1D("hpthat", "",200, 80, 380 );

		TH1D* hvz_weighted    = new TH1D("hvz_weighted", "",100, -15, 15);
		TH1D* hcen_weighted   = new TH1D("hcen_weighted", "",100, 0, 200);
		TH1D* hpthat_weighted = new TH1D("hpthat_weighted", "",200, 80, 380 );

		TFile *fcorr= TFile::Open("test.root");
		TH1D* corr =(TH1D*) fcorr->Get("corr");
		const int ncsv=12;
		float csvBin[12] = {0.1, 0.15, 0.3, 0.5, 0.6, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95, 0.99};
		TH1D* htagged_b[ncsv][2], *htrue_b[2], *htagged_true_b[ncsv][2];
		TH1D* heta[2], *hphi[2];
		TH1I * hcounter[ncsv][2];
		TH1D* hincl[2];

		for(int i=0; i<2; ++i){
				htrue_b        [i] = new TH1D(Form("trueB_%d",i) , "", 100, 100, 500);
				hincl        [i] = new TH1D(Form("incl_%d",i) , "", 100, 100, 500);
				for(int j=0; j<ncsv; ++j){
						hcounter       [j][i] = new TH1I(Form("counter%d_%d",j,i), "", 6, 0, 6);
						htagged_b      [j][i] = new TH1D(Form("taggedB_%d_%d",j,i), "", 100, 100, 500);
						htagged_true_b [j][i] = new TH1D(Form("tagged_trueB_%d_%d",j,i) , "", 100, 100, 500);
				}
		}

		TFile *wf = new TFile("phSpectra.root" , "recreate");
		Long64_t nev = t->fChain->GetEntriesFast();
		float wvz=1, whi=1, wpthat=1, w;
		double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
		double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
		double pthatEntries[9] = {0, 0, 0, 2.56444e+06, 2.84656e+06, 2.67143e+06, 2.88286e+06, 779403, 168046}; //cymbal tune

		auto fit_vz = new TF1("fit_vz","gaus(0)/gaus(3)",-30,30);
		fit_vz->SetParameters(0.08,0.44,5.12,0.08,3.25,5.23);
		auto fit_cen = new TF1("fcent1","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[7]*exp([5]+[6]*x)",0,180);
		fit_cen->SetParameters(4.40810, -7.75301e-02, 4.91953e-04, -1.34961e-06, 1.44407e-09, -160, 1, 3.68078e-15);

		TF1 *jes_cor_f = new TF1("jes_cor_f", "[0]*exp(-pow(log([2]*x-[2]*[1]),2)/pow([3],2)/2)/[2]/(x-[1])/[3]+[4]+[5]*(x-100)",100, 600);
		jes_cor_f->SetParameters(0.1, 100, 0.005, 1.12, 0.9, 0.000023);
		Long64_t jetcounter= 0; 
		//nev=300000;
		Long64_t nmulti=0, nunmatch=0;
		for(Long64_t jentry = 0; jentry< nev; ++jentry){
				t->GetEntry(jentry);
				wvz=fit_vz->Eval(t->vz);
				whi = (t->hiBin<194) ? fit_cen->Eval(1.*t->hiBin) : 1.;
				int ibin =0;
				while(t->pthat>pthatbins[ibin+1]) ibin++;
				wpthat= (xsecs[ibin]-xsecs[ibin+1])/pthatEntries[ibin];
				w= wvz*whi*wpthat;
				if(jentry %100000 == 0) cout<<"processed "<<jentry<<" events..."<<endl;
				hvz->Fill(t->vz); hcen->Fill(t->hiBin); hpthat->Fill(t->pthat);
				hvz_weighted->Fill(t->vz,w); hcen_weighted->Fill(t->hiBin,w); hpthat_weighted->Fill(t->pthat,w);
				int i=0;
				if(t->hiBin > 60) i=1;
				for(int j=0;j < t->genpt->size(); ++j){
						if(t->genpt->at(j)<120) continue;
						if(fabs(t->geneta->at(j))>1.6) continue;
						hincl[i]->Fill(t->genpt->at(j),w);
				}
				for(int j=0;j < t->calo_refpt->size(); ++j){
						//						if(fabs(t->calo_jteta->at(j))>1.6) continue;
						if((t->calo_refpt->at(j))<80) continue;

						//						if((t->calo_trackMax->at(j))/t->calo_corrpt->at(j)<=0.01) continue;

						//if(t->calo_discr_csvV1->at(j)<0.9) continue;
						int imatch = -1;
						double dr = 5;
						int nmatch = 0;
						for(unsigned int igen=0; igen<t->genpt->size(); igen++){
								if( fabs(t->genpt->at(igen)-t->calo_refpt->at(j))< .1) {
										double dr1 = findDR(t->geneta->at(igen), t->genphi->at(igen), t->calo_jteta->at(j),t->calo_jtphi->at(j));
										if(dr> dr1){
												dr = dr1;
												imatch = igen;
										}
								}
						}
						if(nmatch>1) {
								nmulti = nmulti+nmatch;
								//	cout<<"final dr= "<<dr<<endl;
						}
						if(imatch<0 ) {
								nunmatch++;
								continue;
						}

						float matchedPt = t->genpt ->at(imatch);
						float matchedEta= t->geneta->at(imatch);
						float matchedPhi= t->genphi->at(imatch);
						/*
						*/
						//float matchedPt = t->calo_corrpt->at(j);
						if(matchedPt<120) continue;
						if(fabs(matchedEta)>1.6) continue;

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
						jetcounter ++;
						hcounter[0][i]->Fill(5);
						//cout<<"filling"<<endl;
						htrue_b [i]->Fill(matchedPt , w);
						//heta   [i]->Fill(matchedEta, w);
						//hphi   [i]->Fill(matchedPhi, w);
				}
		}
		cout<<"Total jets: "<<jetcounter<<endl;
		cout<<"Unmatched jets: "<<nunmatch<<" ("<<float(nunmatch)/float(jetcounter)<<")%"<<endl;
		cout<<"Multi jets: "<<nmulti<<" ("<<float(nmulti)/float(jetcounter)<<")%"<<endl;
		for(int i=0; i<2; ++i){
				htrue_b       [i]->Write();
				hincl       [i]->Write();
				for(int j=0; j<ncsv; j++){
						htagged_b     [j][i]->Write();
						htagged_true_b[j][i]->Write();
						hcounter      [j][i]->Write();
				}
		}
		hvz->Write();
		hcen->Write();
		hpthat->Write();
		hvz_weighted->Write();
		hcen_weighted->Write();
		hpthat_weighted->Write();
}
