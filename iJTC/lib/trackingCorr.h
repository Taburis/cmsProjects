
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

#include "../Austin_TrkCorr/getTrkCorr.h"
using namespace jetTrack;

//TrkCorr * Austin_corr = new TrkCorr("../Austin_TrkCorr/");

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
		void fixHole(TH2F* h);
		void loadDataTracks(TString datafile);
		void symmetrization(TH2F* h, TH2F* dh);
		void symmetrization2(TH2F* h); // copy the right to the left
		void patchCorr(TH2F* h, float pt, float cent);

	public: 
		inputTree * t;
		int ntrkpt;
		int ncent;
		int nx, ny;
		float *trkpt, *cent;
		float xmin, xmax, ymin, ymax;
		xthf4 * gen;
		xthf4 * rec;
		xthf4 * data_cre;
		xthf4 * data_rec;
		bool readonly;
		TFile* f;
		TFile* data_f;
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

void trackingCorr::loadDataTracks(TString file){
	doSymmetrization = true;
	data_f = TFile::Open(file);
	data_cre= new xthf4(); 
	data_rec= new xthf4(); 
	data_cre->Read("cre", "cre", data_f, ntrkpt, trkpt, ncent, cent);
	data_rec->Read("rec", "rec", data_f, ntrkpt, trkpt, ncent, cent);
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
	trkpt_out = trackingCorrConfig::trkpt_out;
	ncent_out = trackingCorrConfig::ncent_out;
	cent_out = trackingCorrConfig::cent_out;

	/*
	//for gettign aux file to fix the hole problem
	ncent_out = 4;
	float fcent[5] = {0,20, 60, 100, 200};
	cent_out = fcent;
	*/

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
			if(i<2) fixHole(corr[i+j*ntrkpt_out]);
			corr[i+j*ntrkpt_out]->Write();
		}
	}
	//here didn't delete the histogram so the potential memory leak is expected!! needs to fix
	int nptSymm = 2;
	
/*
	if(!doSymmetrization) return;
	TFile* wf2 = TFile::Open("symmetrized_"+file, "recreate");
	// symmetrization only apply to the 3 lowest pt bins
	for(int i=0; i<ntrkpt_out; ++i){
		for(int j=0; j<ncent_out; ++j){
			corr[i+j*ntrkpt_out]= (TH2F*) gen->hist(i,j)->Clone(Form("corr_%d_%d",i,j));
			corr[i+j*ntrkpt_out]->Divide( (TH2F*) rec->hist(i,j));
			if(i<2) fixHole(corr[i+j*ntrkpt_out]);
			if(i<nptSymm) symmetrization((TH2F*)data_cre->hist(i,j), corr[i+j*ntrkpt_out]);
			corr[i+j*ntrkpt_out]->Write();
		}
	}
*/
	TFile* wf2 = TFile::Open("eta_symmetry_"+file, "recreate");
	for(int i=0; i<ntrkpt_out; ++i){
		for(int j=0; j<ncent_out; ++j){
			corr[i+j*ntrkpt_out]= (TH2F*) gen->hist(i,j)->Clone(Form("corr_%d_%d",i,j));
			corr[i+j*ntrkpt_out]->Divide( (TH2F*) rec->hist(i,j));
			if(i<2) fixHole(corr[i+j*ntrkpt_out]);
			if(i<nptSymm) symmetrization2(corr[i+j*ntrkpt_out]);
			corr[i+j*ntrkpt_out]->Write();
		}
	}
	/*
	TFile* wf3 = TFile::Open("patched_"+file, "recreate");

	for(int i=0; i<ntrkpt_out; ++i){
		for(int j=0; j<ncent_out; ++j){
			corr[i+j*ntrkpt_out]= (TH2F*) gen->hist(i,j)->Clone(Form("corr_%d_%d",i,j));
			corr[i+j*ntrkpt_out]->Divide( (TH2F*) rec->hist(i,j));
			if(i<2) fixHole(corr[i+j*ntrkpt_out]);
			if(i<nptSymm) patchCorr(corr[i+j*ntrkpt_out],(trkpt_out[i]+trkpt_out[i+1])/2, (cent_out[j]+cent_out[j+1])/2);
			corr[i+j*ntrkpt_out]->Write();
		}
	}
	*/

}

void trackingCorr::showCorr(int ih, int jh){
	int nc= ceil(float(ncent_out)/20);
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
					for(int k=0;k<ncent_out && k<20; ++k){
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
				for(int k=0;k<ncent_out && k<20; ++k){
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

//for fixing the hole in the tracker in 2015 PbPb data
void trackingCorr::fixHole(TH2F* h){
	int xh=36, yh=65;
	int radius=3;
	float peak=0;
	int ip, jp;
	for(int i=xh-radius; i<xh+radius+1;i++){
		for(int j=yh-radius; j<yh+radius+1;j++){
			//	if(h->GetBinContent(i,j)!=0) continue;
			if(peak< h->GetBinContent(i,j)){
				peak = h->GetBinContent(i,j);
				ip = i;
				jp = j;
			}
		}
	}
	int rr=2;
	float boundary=float(h->GetBinContent(ip,jp+radius)+h->GetBinContent(ip,jp-radius)
			+h->GetBinContent(ip+radius,jp)+h->GetBinContent(ip-radius,jp))/4;
	for(int i=ip-radius; i<ip+radius+1;i++){
		for(int j=jp-radius; j<jp+radius+1;j++){
			//cout<<"!"<<i<<", "<<j<<endl;
			if(h->GetBinContent(i,j)==0){
				float bin = peak- (peak-boundary)/5*pow(pow((j-jp),2)+pow((i-ip),2),0.5);
				if(bin<0){ 
					cout<<"!"<<endl;
					bin = boundary;
				}
				h->SetBinContent(i, j, bin);
			}
		}
	}
}

void trackingCorr::symmetrization2(TH2F* h){
	int eta_low=34, eta_up=41, phi_low=63, phi_up=67;
	for(int i=eta_low; i<eta_up+1; ++i){
		for(int j=phi_low; j<phi_up+1; ++j){
			float content = h->GetBinContent(101-i, j);
			cout<<" eta="<<h->GetXaxis()->GetBinCenter(i)<<" reflected eta="<<h->GetXaxis()->GetBinCenter(101-i)<<endl;
//			cout<<content<<endl;
			h->SetBinContent(i,j,content );
		}
	}
}
void trackingCorr::symmetrization(TH2F* h, TH2F* corr){
	float cont ;
	int nx = h->GetNbinsX();
	int ny = h->GetNbinsY();
	TH2F* temp = (TH2F*)h->Clone();
	for(int j=0; j<ny+1;j++){
		for(int i=1;i<nx+1;i++){
//			if( j ==72)cout<<j<<endl;
			if( float(ny)/2<j && h->GetBinContent(i,j) !=0){
				cont = h->GetBinContent(i,j);
				cont = h->GetBinContent(i,ny-j+1)/cont;
				//                      if(cont == 0){
				//                              cout<<"den: "<<h->GetBinContent(i,j)<<endl;
				//                              cout<<"num: "<<h->GetBinContent(i,ny-j+1)<<" at x= "<<\
				//                                      h->GetXaxis()->GetBinCenter(i)<<"; y= "<<\
				//                                      h->GetYaxis()->GetBinCenter(ny-j+1)<<endl;
				//                              cout<<"with i= "<<i<<"; "<<"j= "<<ny-j+1<<endl;
				//                      }
				temp->SetBinContent(i,j, cont);
			}
			//      else if( h->GetBinContent(i,j) !=0) temp->SetBinContent(i,j, 1);
			else temp->SetBinContent(i,j, 1);
		}
	}
	corr->Multiply(temp);
	delete temp;
}

void trackingCorr::patchCorr(TH2F* h, float pt, float cent){
	int eta_low=32, eta_up=43, phi_low=64, phi_up=67;
	for(int i=eta_low; i<eta_up+1; ++i){
		for(int j=phi_low; j<phi_up+1; ++j){
			float eta=h->GetXaxis()->GetBinCenter(i);
			float phi=h->GetYaxis()->GetBinCenter(j);
			float content = Austin_corr->getTrkCorr(pt, eta, phi, cent, 999);
			h->SetBinContent(i,j,content );
		}
	}
}

#endif
