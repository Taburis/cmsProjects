
#ifndef JTCSignalProducer_H
#define JTCSignalProducer_H

#ifndef signalFactoryBase_H
#include "signalFactoryBase.h"
#endif

const Double_t etabin[22] ={-3, -2.5,-2.,-1.5, -1., -0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.5,2.,2.5, 3};
const Double_t phibin[18] ={-1.50796, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531,1.50796};
const float drbin [15] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};

class JTCSignalProducer :public signalFactoryBase {
		public :
				JTCSignalProducer() : signalFactoryBase() {};
				JTCSignalProducer(TH2D* h, TH2D* hmix = NULL) : signalFactoryBase() {
						raw_sig = h; mix = hmix;
				};
				TH2D* getSignal(TString name);
				void read(TFile *f, TString name);
				void read1D(TFile *f, TString name);
				TH1* projectionX(TH2D* h, float x, float y);
				TH1* projectionY(TH2D* h, float x, float y);
				TH1* projX(bool doRebin, TH2D* h, float x, float y); // for general projection 
				TH1* projY(bool doRebin, TH2D* h, float x, float y); // for general projection 
				TH1* signal_X(bool doRebin=1){ return projX( doRebin, sig, -1, 1);}
				TH1* raw_X(bool doRebin=1){ return projX( doRebin, raw_sig, -1, 1);}
				TH1D* doDrIntegral(TString name);
				void WriteTH2();
				void WriteTH1();
				void getAllProj(TString name);
				void drawBkgCheck(bool doX = 1);
				void drawSideBandCheck();
				void histStyle(TH1* h);
		public :
				// 3 steps to get the final signal:
				// 1: get the raw_sig by skiming
				// 2: correct the raw_sig by the mixing table
				// 3: subtract the bkg to get the sig
				TH2D* raw_sig =NULL;
				TH2D* sig =NULL;
				TH2D* sig_step2 =NULL;  // right after the mixing correction
				TH2D* mix =NULL;
				TH2D* mix_normalized =NULL;
				TH2D* bkg =NULL;
				TH1D* dr_integral = NULL;
				bool doSmoothME = true;
				float sideMin = 1.5, sideMax = 2.5;
				TH1D *sig_deta=0, *sig_dphi=0, *sig_step2_deta=0, *sig_step2_dphi=0,*bkg_deta=0, *bkg_dphi=0,*side_deta=0, *side_deta_mix=0;
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
		sig_step2 = (TH2D*) sig->Clone("sig_mix_corrected_"+name);
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
		sig_step2->Write();
}

void JTCSignalProducer::read(TFile *f , TString name){
		TString tmp;
		tmp = "signal_"+name;
		//cout<<tmp<<endl;
		sig=(TH2D*) f->Get(tmp);
		tmp = "smoothed_mixing_"+name;
		mix_normalized=(TH2D*) f->Get(tmp);
		tmp = "raw_"+name;
		raw_sig=(TH2D*) f->Get(tmp);
		tmp = "bkg_"+name;
		bkg=(TH2D*) f->Get(tmp);
		tmp = "sig_mix_corrected_"+name;
		sig_step2=(TH2D*) f->Get(tmp);
}

TH1* JTCSignalProducer::projectionX(TH2D* h, float x, float y){
		int xbin = h->GetYaxis()->FindBin(x);
		int ybin = h->GetYaxis()->FindBin(y)-1;
		TString name = h->GetName(); name = "projX_"+name;
		return h->ProjectionX(name, xbin, ybin );
}

TH1* JTCSignalProducer::projectionY(TH2D* h, float x, float y){
		int xbin = h->GetXaxis()->FindBin(x);
		int ybin = h->GetXaxis()->FindBin(y)-1;
		TString name = h->GetName(); name = "projY_"+name;
		return h->ProjectionY(name, xbin, ybin );
}

TH1* JTCSignalProducer::projX(bool doRebin, TH2D*h2, float x, float y){
		// here h2 needs to be invariant 
		TH1* h=projectionX(h2, x, y);
		h->Scale(h2->GetYaxis()->GetBinWidth(1));  
		if(doRebin){
				TString name = h->GetName();
				name = "rebined_"+name;
			   	TH1* hh=h;
				h=invariantRebin(h,name, 21, etabin);
				delete hh;
		}
		return h;
}

TH1* JTCSignalProducer::projY(bool doRebin, TH2D*h2, float x, float y){
		// here h2 needs to be invariant 
		TH1* h=projectionY(h2, x, y);
		h->Scale(h2->GetXaxis()->GetBinWidth(1));  
		if(doRebin){
				TString name = h->GetName();
				name = "rebined_"+name;
			   	TH1* hh=h;
				h=invariantRebin(h,name, 17, phibin);
				delete hh;
		}
		return h;
}

void JTCSignalProducer::getAllProj(TString name){
		TString tmp;
		tmp = "signal_deta_"+name;
		sig_deta = (TH1D*) projX(1, sig, -1, 1); sig_deta->SetName(tmp);
		tmp = "signal_dphi_"+name;
		sig_dphi = (TH1D*) projY(1, sig, -1, 1); sig_dphi->SetName(tmp);

		tmp = "sig_step2_deta_"+name;
		sig_step2_deta = (TH1D*) projX(1, sig_step2, -1, 1); sig_step2_deta->SetName(tmp);
		tmp = "sig_step2_dphi_"+name;
		sig_step2_dphi = (TH1D*) projY(1, sig_step2, -1, 1); sig_step2_dphi->SetName(tmp);

		tmp = "bkg_deta_"+name;
		bkg_deta = (TH1D*) projX(1, bkg, -1, 1); bkg_deta->SetName(tmp);
		tmp = "bkg_dphi_"+name;
		bkg_dphi = (TH1D*) projY(1, bkg, -1, 1); bkg_dphi->SetName(tmp);

		tmp = "side_deta_"+name;
		side_deta = (TH1D*) projX(1, sig, 1.2, 2.2); side_deta->SetName(tmp);

		tmp = "side_deta_mix_"+name;
		side_deta_mix = (TH1D*) projX(1, mix_normalized, 1.2, 2.2); side_deta_mix->SetName(tmp);
}


void JTCSignalProducer::WriteTH1(){
		sig_deta->Write(); sig_step2_deta->Write(); bkg_deta->Write();
		sig_dphi->Write(); sig_step2_dphi->Write(); bkg_dphi->Write();
		side_deta->Write();
		side_deta_mix->Write();
}

void JTCSignalProducer::read1D(TFile *f , TString name){
		TString tmp;
		tmp = "signal_deta_"+name; sig_deta=(TH1D*) f->Get(tmp);
		tmp = "signal_dphi_"+name; sig_dphi=(TH1D*) f->Get(tmp);
		//cout<<tmp<<endl;

		tmp = "sig_step2_deta_"+name; sig_step2_deta=(TH1D*) f->Get(tmp);
		tmp = "sig_step2_dphi_"+name; sig_step2_dphi=(TH1D*) f->Get(tmp);

		tmp = "bkg_deta_"+name; bkg_deta=(TH1D*) f->Get(tmp);
		tmp = "bkg_dphi_"+name; bkg_dphi=(TH1D*) f->Get(tmp);

		tmp = "side_deta_"+name; side_deta=(TH1D*) f->Get(tmp);
		tmp = "side_deta_mix_"+name; side_deta_mix=(TH1D*) f->Get(tmp);
}

void JTCSignalProducer::drawBkgCheck(bool doX){
		TH1 *h1, *h2;
		if(doX) {
				h1 = sig_step2_deta; 
				h2 = bkg_deta; 
		}
		else {
				h1 =sig_step2_dphi; 
				h2 =bkg_dphi; 
		}
		histStyle(h1);
		float mean = h2->GetBinContent(h1->FindBin(0.5));
		float dvt =  h2->GetBinError(h1->FindBin(0.5));
		h1->SetAxisRange(mean-4*dvt, mean+6*dvt, "Y");
		h1->SetLineWidth(2);
		h2->SetLineWidth(2);
		h1->SetLineColor(kBlack);
		h2->SetLineColor(kRed);
		h1->Draw("same");
		h2->Draw("same");
}


void JTCSignalProducer::drawSideBandCheck(){
		TH1 *h1, *h2;
		h1 = sig_deta; 
		h2 = side_deta; 
		cout<<h2->GetName()<<endl;
		histStyle(h1);
		float mean = h2->GetBinContent(h1->FindBin(0));
		float dvt =  h2->GetBinError(h1->FindBin(0));
		h1->SetAxisRange(mean-6*dvt, mean+10*dvt, "Y");
		h1->SetLineWidth(2);
		h2->SetLineWidth(2);
		h1->SetLineColor(kBlack);
		h2->SetLineColor(kOrange+7);
		h1->Draw("same");
		h2->Draw("same");
}

void JTCSignalProducer::histStyle(TH1* h){
		h->SetTitle("");
		h->GetXaxis()->SetLabelSize(0.08);
		h->GetYaxis()->SetLabelSize(0.08);
		h->GetXaxis()->CenterTitle();
		h->GetXaxis()->SetTitleOffset(1);
		h->GetXaxis()->SetTitleSize(0.08);
		h->GetYaxis()->SetTitleSize(0.08);
}

#endif 
