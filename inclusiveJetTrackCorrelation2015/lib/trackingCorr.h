
#ifndef config_H
#include "config_v1.h"
#endif

#ifndef inputTree_H
#include "inputTree.h"
#endif

#ifndef xthf4_H
#include "xthf4.h"
#endif
#ifndef trackingCorr_H
#define trackingCorr_H
using namespace jetTrack;

class trackingCorr  {
	public: 
		trackingCorr(inputTree* tt):
			t(tt),
			ntrkpt(trackingCorrConfig::ntrkpt_in),
			ncent (trackingCorrConfig::ncent_in),
			trkpt (trackingCorrConfig::trkpt_in),
			cent  (trackingCorrConfig::cent_in),
			nx    (trackingCorrConfig::neta),
			ny    (trackingCorrConfig::nphi),
			xmin  (trackingCorrConfig::etamin),
			xmax  (trackingCorrConfig::etamax),
			ymin  (trackingCorrConfig::phimin),
			ymax  (trackingCorrConfig::phimax)
	{
		readonly =0;
		corr=NULL;
		gen = new xthf4("gen","gen", nx,xmin, xmax , ny, ymin, ymax,ntrkpt, trkpt,ncent, cent);
		rec = new xthf4("rec","rec", nx,xmin, xmax , ny, ymin, ymax,ntrkpt, trkpt,ncent, cent);
	};
		trackingCorr(TFile *ff ):
			ntrkpt(trackingCorrConfig::ntrkpt_in),
			ncent (trackingCorrConfig::ncent_in),
			trkpt (trackingCorrConfig::trkpt_in),
			cent  (trackingCorrConfig::cent_in),
			nx    (trackingCorrConfig::neta),
			ny    (trackingCorrConfig::nphi),
			xmin  (trackingCorrConfig::etamin),
			xmax  (trackingCorrConfig::etamax),
			ymin  (trackingCorrConfig::phimin),
			ymax  (trackingCorrConfig::phimax),
			f(ff)
	{
		readonly = 1;
		corr=NULL;
		gen = new xthf4();
		rec = new xthf4();
	}

		void runScan();
		void Read();
		void getCorr(TString file);
		void showCorr(int i=-1, int j=-1);

	public: 
		inputTree * t;
		int ntrkpt;
		int ncent;
		int nx, ny;
		float *trkpt, *cent;
		float xmin, xmax, ymin, ymax;
		xthf4 * gen;
		xthf4 * rec;
		bool readonly;
		TFile* f;
		TH2F** corr;
		int ntrkpt_out  ;
		int ncent_out   ;
		float *trkpt_out;
		float *cent_out ;
};

void trackingCorr::runScan(){
	Long64_t nentries = t->fChain->GetEntriesFast();
	for(Long64_t jentry = 0; jentry <nentries ; ++jentry){
		t->GetEntry(jentry);
		float w_vz=1;
		if(eventCuts(t)) continue;
		w_vz = weight::vz_weight(t->vz);
		for(int j=0;j<t->trkPt->size(); ++j){
			if(trackQualityCuts(t,j))continue;
			rec->Fill(t->trkEta->at(j), t->trkPhi->at(j), t->trkPt->at(j), t->hiBin, w_vz);
		}
		for(int j=0;j<t->pt->size(); ++j){
			if(genParticleCuts(t,j))continue;
			gen->Fill(t->eta->at(j), t->phi->at(j), t->pt->at(j), t->hiBin, w_vz);
		}
	}
	gen->Write();
	rec->Write();
}

void trackingCorr::Read(){
	if(!readonly) {
		std::cout<<"already filled by obj, can't load anymore!"<<std::endl;
		return ;
	}
	gen->Read("gen", "gen", f, ntrkpt, trkpt, ncent, cent);
	rec->Read("rec", "rec", f, ntrkpt, trkpt, ncent, cent);
	return ;
}

void trackingCorr::getCorr(TString file){
	TFile* wf = TFile::Open(file, "recreate");
	ntrkpt_out = trackingCorrConfig::ntrkpt_out;
	ncent_out = trackingCorrConfig::ncent_out;
	trkpt_out = trackingCorrConfig::trkpt_out;
	cent_out = trackingCorrConfig::cent_out;
	gen->RebinZ(ntrkpt_out, trkpt_out);
	rec->RebinZ(ntrkpt_out, trkpt_out);
	gen->RebinW(ncent_out, cent_out);
	rec->RebinW(ncent_out, cent_out);

	int nhist = ntrkpt_out*ncent_out;
	corr= new TH2F*[nhist];
	for(int i=0; i<ntrkpt_out; ++i){
		for(int j=0; j<ncent_out; ++j){
			corr[i+j*ntrkpt_out]= (TH2F*) gen->hist(i,j)->Clone(Form("corr_%d_%d",i,j));
			corr[i+j*ntrkpt_out]->Divide( (TH2F*) rec->hist(i,j));
			corr[i+j*ntrkpt_out]->Write();
		}
	}
//	wf->Close();
}

void trackingCorr::showCorr(int ih, int jh){
	int nc= ceil(float(ncent_out)/10);
	if( corr==NULL){
		std::cout<<"no correction loaded!"<<std::endl;
		return;
	}
	TCanvas **cc;
	//show in 4 by 5 subpads
	if(jh==-1){
		if(ih==-1){
			cc= new TCanvas*[int(nc*ntrkpt_out)];
			for(int j=0; j<ntrkpt_out; ++j){
				for(int ic=0; ic<nc; ++ic){
					cc[ic+j*nc]=new TCanvas(Form("ccorr_%d",ic+j*nc),"", 2000,1600);
					cc[ic+j*nc]->Divide(5,4,0,0);
					for(int k=0;k<ncent_out && k<10; ++k){
						cc[ic+j*nc]->cd(k+1);
						corr[j+k*ntrkpt_out]->Draw("colz");
					}
				}
			}
		}
		else if(ih< ntrkpt_out) {
			cc= new TCanvas*[nc];
			for(int ic=0; ic<nc; ++ic){
				cc[ic]=new TCanvas(Form("ccorr_%d",ic),"", 2000,1600);
				cc[ic]->Divide(5,4,0,0);
				for(int k=0;k<ncent_out && k<10; ++k){
					cc[ic]->cd(k+1);
					cout<<k<<endl;
					corr[ih+k*ntrkpt_out]->Draw("colz");
				}
			}
		}
	}
	else if(ih< ntrkpt_out && jh<ncent_out){
		auto ccc= new TCanvas("cccc","", 600,500);
		corr[ih+jh*ntrkpt_out]->Draw("colz");
	}
}
#endif
