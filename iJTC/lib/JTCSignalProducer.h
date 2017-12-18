
#ifndef JTCSignalProducer_H
#define JTCSignalProducer_H

#ifndef signalFactoryBase_H
#include "signalFactoryBase.h"
#endif

const Double_t etabin[22] ={-3, -2.5,-2.,-1.5, -1., -0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.5,2.,2.5, 3};
const float drbin [15] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};
class JTCSignalProducer :public signalFactoryBase {
		public :
				JTCSignalProducer() : signalFactoryBase() {};
				JTCSignalProducer(TH2D* h, TH2D* hmix = NULL) : signalFactoryBase() {
						raw_sig = h; mix = hmix;
				};
				TH2D* getSignal(TString name);
				void read(TFile *f, TString name);
				TH1* projectionX(TH2D* h, float x, float y);
				TH1* projX(bool doRebin, TH2D* h, float x, float y); // for general projection 
				TH1* signal_X(bool doRebin=1){ return projX( doRebin, sig, -1, 1);}
				TH1* raw_X(bool doRebin=1){ return projX( doRebin, raw_sig, -1, 1);}
				TH1D* doDrIntegral(TString name);
				void WriteTH2();
		public :
				TH2D* raw_sig =NULL;
				TH2D* sig =NULL;
				TH2D* mix =NULL;
				TH2D* mix_normalized =NULL;
				TH2D* bkg =NULL;
				TH1D* dr_integral = NULL;
				bool doSmoothME = true;
				float sideMin = 1.5, sideMax = 2.5;
};

TH2D* JTCSignalProducer::getSignal(TString name){
	   	raw_sig->SetName("raw_"+name);
		mix->SetName("mixing_"+name);
		sig = (TH2D*) raw_sig->Clone("signal_"+name);
		sig->GetXaxis()->SetTitle("d#eta");
		sig->GetYaxis()->SetTitle("d#phi");
		sig->GetXaxis()->CenterTitle();
		sig->GetYaxis()->CenterTitle();
		sig->Scale(1.0/sig->GetXaxis()->GetBinWidth(1)/sig->GetYaxis()->GetBinWidth(1));
		mix_normalized= mixingTableMaker(mix, doSmoothME);
		mix_normalized->SetName("smoothed_mixing_"+name);
		sig->Divide(mix_normalized);
		bkg = (TH2D*) getV2Bkg(sig,sideMin , sideMax );
		bkg->SetName("bkg_"+name);
		sig->Add(sig, bkg, 1, -1);
		return sig;
}

TH1D* JTCSignalProducer::doDrIntegral(TString name ){
		TString title = sig->GetTitle();
		TString hname = "dr_dist_"+name;
		dr_integral = drDistMaker(sig, hname, title, 14, drbin );
		dr_integral->GetXaxis()->SetTitle("#Deltar");
		return dr_integral;
}

void JTCSignalProducer::WriteTH2(){
		raw_sig->Write();
		sig->Write();
		mix_normalized->Write();
		bkg->Write();
}

void JTCSignalProducer::read(TFile *f , TString name){
		TString tmp;
		tmp = "signal_"+name;
//		cout<<tmp<<endl;
		sig=(TH2D*) f->Get(tmp);
		tmp = "smoothed_mixing_"+name;
		mix_normalized=(TH2D*) f->Get(tmp);
		tmp = "raw_"+name;
		raw_sig=(TH2D*) f->Get(tmp);
		tmp = "bkg_"+name;
		bkg=(TH2D*) f->Get(tmp);
}

TH1* JTCSignalProducer::projectionX(TH2D* h, float x, float y){
		int xbin = h->GetYaxis()->FindBin(x);
		int ybin = h->GetYaxis()->FindBin(y)-1;
		//cout<<xbin<<endl;
		//cout<<ybin<<endl;
		TString name = h->GetName(); name = "projX_"+name;
		return h->ProjectionX(name, xbin, ybin );
}

TH1* JTCSignalProducer::projX(bool doRebin, TH2D*h2, float x, float y){
		TH1* h=projectionX(h2, x, y);
		h->Scale(sig->GetYaxis()->GetBinWidth(1));  // as signal is invariant in default
		if(doRebin){
				TString name = h->GetName();
				name = "rebined_"+name;
			   	TH1* hh=h;
				h=invariantRebin(h,name, 21, etabin);
				delete hh;
		}
		return h;
}

#endif 
