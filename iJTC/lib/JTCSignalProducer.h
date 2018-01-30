
#ifndef JTCSignalProducer_H
#define JTCSignalProducer_H

#ifndef signalFactoryBase_H
#include "signalFactoryBase.h"
#endif

const Double_t etabin[22] ={-3, -2.5,-2.,-1.5, -1., -0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.5,2.,2.5, 3};
const Double_t phibin[18] ={-1.50796, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531,1.50796};
const float drbin [16] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1., 1.2};


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
				void getAllProj(TString name, bool rebin = 1);
				void drawBkgCheck(bool doX = 1);
				void drawSideBandCheck();
				void histStyle(TH1* h);
				void gausFit();
				float getBkgError();
				TH1D* getDrBkgErr();
				
		public :
				// 3 steps to get the final signal:
				// 1: get the raw_sig by skiming
				// 2: correct the raw_sig by the mixing table
				// 3: subtract the bkg to get the sig
				float sidebandmin=1.2 , sidebandmax=2.2;
				//float sidebandmin=-TMath::Pi()/2 , sidebandmax=-1.2;
				bool doSideBandMixing = 0;
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
				TF1* fdeta=0, *fdphi=0;
				TH1D *bkg_est=0;
				float bkgErr, mixErr;
};

TH2D* JTCSignalProducer::getSignal(TString name){
	   	raw_sig->SetName("raw_"+name);
		mix->SetName("mixing_"+name);
		sig = (TH2D*) raw_sig->Clone("signal_"+name);
		sig->GetXaxis()->SetTitle("d#eta");
		sig->GetYaxis()->SetTitle("d#phi");
		sig->GetXaxis()->CenterTitle();
		sig->GetYaxis()->CenterTitle();
		sig->Scale(1.0/sig->GetXaxis()->GetBinWidth(1)/sig->GetYaxis()->GetBinWidth(1)); //make the h2 invariant
		if(doSideBandMixing) {
				mix_normalized= sideBandMixingTableMaker(sig, sidebandmin, sidebandmax);
				//mix_normalized->SetName("sideBand_mixing_"+name);
		}
		else {
				mix_normalized= mixingTableMaker(mix, doSmoothME);
		}
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
		int ybin = h->GetYaxis()->FindBin(y);
		ybin = h->GetYaxis()->FindBin(y-h->GetYaxis()->GetBinWidth(ybin));
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

void JTCSignalProducer::getAllProj(TString name, bool rebin){
		TString tmp;
		tmp = "signal_deta_"+name;
		sig_deta = (TH1D*) projX(rebin, sig, -1, 1); sig_deta->SetName(tmp);
		tmp = "signal_dphi_"+name;
		sig_dphi = (TH1D*) projY(rebin, sig, -1, 1); sig_dphi->SetName(tmp);

		tmp = "sig_step2_deta_"+name;
		sig_step2_deta = (TH1D*) projX(rebin, sig_step2, -1, 1); sig_step2_deta->SetName(tmp);
		tmp = "sig_step2_dphi_"+name;
		sig_step2_dphi = (TH1D*) projY(rebin, sig_step2, -1, 1); sig_step2_dphi->SetName(tmp);

		tmp = "bkg_deta_"+name;
		bkg_deta = (TH1D*) projX(rebin, bkg, -1, 1); bkg_deta->SetName(tmp);
		tmp = "bkg_dphi_"+name;
		bkg_dphi = (TH1D*) projY(rebin, bkg, -1, 1); bkg_dphi->SetName(tmp);

		tmp = "side_deta_"+name;
		side_deta = (TH1D*) projX(rebin, sig, sidebandmin, sidebandmax); side_deta->SetName(tmp);

		tmp = "side_deta_mix_"+name;
		side_deta_mix = (TH1D*) projX(rebin, mix_normalized, 1.2, 2.2); side_deta_mix->SetName(tmp);
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
//		float dvt =  fabs(0.05*h2->GetMaximum());
		h2->SetAxisRange(h2->GetMinimum()-5*dvt, h2->GetMaximum()+6*dvt, "Y");
		h1->SetLineWidth(2);
		h2->SetLineWidth(2);
		h1->SetLineColor(kBlack);
		h2->SetLineColor(kRed);
		h2->Draw("same");
		h1->Draw("same");
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

void JTCSignalProducer::gausFit(){
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

float JTCSignalProducer::getBkgError(){
		if(bkg_est==0){
			   	bkg_est = (TH1D*) projX(1, sig, -TMath::Pi()/2, TMath::Pi()/2);
		}
		int in_left = bkg_est->GetXaxis()->FindBin(-1.99);
		int in_right= bkg_est->GetXaxis()->FindBin(1.99);
		int out_left = bkg_est->GetXaxis()->FindBin(-2.49);
		int out_right= bkg_est->GetXaxis()->FindBin(2.49);
		float in_sum = bkg_est->GetBinContent(in_left)+bkg_est->GetBinContent(in_right);
//		cout<<"in sum: "<<in_sum<<endl;
		float out_sum = bkg_est->GetBinContent(out_left)+bkg_est->GetBinContent(out_right);
//		cout<<"out sum: "<<out_sum<<endl;
		float err = fabs(in_sum) >fabs(out_sum) ? in_sum : out_sum;
		float left_sum = bkg_est->GetBinContent(in_left)+bkg_est->GetBinContent(out_left);
		float right_sum = bkg_est->GetBinContent(in_right)+bkg_est->GetBinContent(out_right);
		bkgErr = err/2;
		err = fabs(left_sum) >fabs(right_sum) ? left_sum : right_sum;
		mixErr = err/2;
		float ave = (in_sum+out_sum)/4;
		return ave;
//		return bkgErr;
}

TH1D* JTCSignalProducer::getDrBkgErr(){
		TString name = dr_integral->GetName();
		name = name + "_bkgSystErr";
		TH1D* bkg = (TH1D*) dr_integral ->Clone(name);
		for(int k=1; k<15; k++){
		float ring =TMath::Pi()*(pow(bkg->GetBinLowEdge(k)+bkg->GetBinWidth(k),2)-
						pow(bkg->GetBinLowEdge(k),2))/2/bkg->GetBinWidth(k);		
		float err = pow(pow(bkgErr,2) + pow(mixErr, 2 ), 0.5);
						bkg->SetBinError (k, err*ring);
						bkg->SetBinContent(k, 1);
		}
		return bkg;
}

#endif 
