#ifndef JTCSignalProducer_H
#define JTCSignalProducer_H

#ifndef signalFactoryBase_H
#include "signalFactoryBase.h"
#endif

const Double_t etabin[24] ={-3.5, -3, -2.5,-2.,-1.5, -1., -0.8, -0.6, -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6, 0.8, 1., 1.5,2.,2.5, 3, 3.5};
const Double_t phibin[18] ={-1.50796, -1.00531,-0.879646, -.75398, -0.628319,-0.502655, -0.376991, -0.251327, -0.125664, 0.125664, 0.251327, 0.376991, 0.502655, 0.628319,.75398, 0.879646, 1.00531,1.50796};
const float drbin [16] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1., 1.2};

//Double_t mix_pol2(Double_t *x, Double_t *par){
//		if(x[0] >-.3 && x[0] < .299) {
//				TF1::RejectPoint();
//				return 0;
//		}
//		return par[0]+x[0]*par[1]+pow(x[0],2)*par[2];
////		return par[0]+pow(x[0],2)*par[1];
//}

Double_t mix_pol2(Double_t *x, Double_t *par){
		//fitting function is: a0+a1*|x-x0|+a2*x^2+a3*x;
		if(x[0] < 0.5 && x[0] > -0.5) return par[0];
//		else if( x[0]<-0.3 ) return -par[1]*(x[0]+0.5)+par[0]+par[2]*x[0]*x[0]+par[3]*x[0];
//		else return par[1]*(x[0]-0.5)+par[0]+par[2]*x[0]*x[0]+par[3]*x[0];
		else if( x[0]<-0.3 ) return par[0]+par[1]*x[0]*x[0]+par[2]*x[0];
		else return par[0]+par[1]*x[0]*x[0]+par[2]*x[0];
//		else if( x[0]<-0.5 ) return -par[1]*(x[0]+0.5)+par[0]+par[2]*x[0]+par[3]*x[0]*x[0];
//		else return par[1]*(x[0]-0.5)+par[0]+par[2]*x[0]+par[3]*x[0]*x[0];
	//	else if( x[0]<-0.5 ) return par[0]+par[1]*pow(x[0]+0.5,2);
	//	else return par[1]*(x[0]-0.5)+par[0]+par[1]*pow(x[0]-0.5,2);
}

class JTCSignalProducer :public signalFactoryBase {
		public :
				JTCSignalProducer() : signalFactoryBase() {};
				JTCSignalProducer(TH2D* h, TH2D* hmix = NULL) : signalFactoryBase() {
						raw_sig = h; mix = hmix;
				};
				~JTCSignalProducer();
				TH2D* getSignal(TString name, bool doSeagullCorr = 0);
				void pullSeagull(TH2D* hsig);
				TH2D* mixingCorr(TString name, TH2D* , TH2D*, bool );
				TH2D*  doBkgSubtraction(TString name, TH2D* h, float sideMin = 1.5, float sideMax = 2.5);
				void read(TFile *f, TString name);
				void read1D(TFile *f, TString name);
				TH1* projX(bool doRebin, TH2D* h, float x, float y, TString opt=""); // for general projection 
				TH1* projY(bool doRebin, TH2D* h, float x, float y, TString opt=""); // for general projection 
				TH1* signal_X(bool doRebin=1){ return projX( doRebin, sig, -1, 1);}
				TH1* raw_X(bool doRebin=1){ return projX( doRebin, raw_sig, -1, 1);}
				TH1D* doDrIntegral(TString name);
				float getMean(TH1* h, float x1, float x2);
				void WriteTH2();
				void WriteTH1();
				void getAllProj(TString name, bool rebin = 1);
				void drawBkgCheck(bool doX = 1);
				void drawSliceSideBand(TString name);
				void drawSideBandCheck();
				TH1D* getSignal_phiSideBand(TString name);
				TH1D* getSignal_dEta(TString name);
				void histStyle(TH1* h);
				void gausFit();
				float getBkgError();
				TH1D* getDrBkgErr();

		public :
				// 3 steps to get the final signal:
				// 1: get the raw_sig by skiming
				// 2: correct the raw_sig by the mixing table
				// 3: subtract the bkg to get the sig
				float sidebandmin=1.4 , sidebandmax=1.8;
				//float sidebandmin=-TMath::Pi()/2 , sidebandmax=-1.2;
				bool doSideBandMixing = 0, shiftSignal=0;
				TH2D* raw_sig =0;
				TH2D* sig =0;
				TH2D* sig_step2 =0;  // right after the mixing correction
				TH2D* mix =0;
				TH2D* mix_normalized =0;
				TH2D* bkg =0;
				TH1D* dr_integral = 0;
				bool doSmoothME = true;
				float sideMin = 1.5, sideMax = 2.5;
				TH1D *sig_deta=0, *sig_dphi=0, *sig_step2_deta=0, *sig_step2_dphi=0,*bkg_deta=0, *bkg_dphi=0,*side_deta=0, *side_deta_mix=0;
				TF1* fdeta=0, *fdphi=0;
				TH1D *bkg_est=0;
				float bkgErr, mixErr;
};

void JTCSignalProducer::pullSeagull(TH2D* hsig){
		TF1 *func = new TF1("func", mix_pol2, -3., 3., 3);
	//	TF1 *func = new TF1("func", "pol2", -2.5, 2.5);
		int n1 = hsig->GetYaxis()->FindBin(1.4);
		int n2 = hsig->GetYaxis()->FindBin(1.799);
		TH1D *htm = (TH1D*) hsig->ProjectionX("side_for_fitting", n1, n2); 
		htm->Scale(1.0/(n2-n1));
		htm->Rebin(4); htm->Scale(0.25);
//		TH1D *htm1=(TProfile*)invariantRebin((TH1*)htm, "tet_fitted", 21, etabin);
		float mean = getMean(htm, -0.5, 0.5);
		func->SetParameters(mean, 0, 0);
		htm->Fit(func, "", "", -3., 2.99);
		mean = func->GetParameter(0);
		//htm1->Fit(func);
		auto line = new TLine(); line->DrawLine(-3, mean, 3, mean);
		for(int i=1; i<hsig->GetNbinsX()+1; ++i){
				float x = hsig->GetXaxis()->GetBinCenter(i);
//				if(fabs(x)<0.3) continue;
				float corr = mean/func->Eval(x);
				for(int j=1; j<hsig->GetNbinsY()+1; ++j){
						if(hsig->GetBinContent(i,j) == 0) continue;
						float cont= hsig->GetBinContent(i,j);
						hsig->SetBinContent(i,j, cont*corr); 
			}
		}

//		delete htm;
//		return 0;
}

float JTCSignalProducer::getMean(TH1* h, float x1, float x2){
		int n1 = h->FindBin(x1);
		int n2 = h->FindBin(x2);
		return h->Integral(n1, n2)/(n2-n1+1);
}

TH2D* JTCSignalProducer::doBkgSubtraction(TString name, TH2D* h, float sideMin, float sideMax){
		TH2D* hh =(TH2D*) h->Clone(name);
		TH2D* htmp = (TH2D*) getV2Bkg(hh,sideMin , sideMax );
//		htmp->SetName("bkg_"+name);
		hh->Add(hh, htmp, 1, -1);
		delete htmp;
		return hh;
}

TH2D* JTCSignalProducer::mixingCorr(TString name, TH2D* h, TH2D* hmix, bool doseagull){
		TH2D* hsig =(TH2D*) h->Clone(name);
		hsig->Scale(1.0/h->GetXaxis()->GetBinWidth(1)/h->GetYaxis()->GetBinWidth(1)); //make the h2 invariant
		hsig->GetXaxis()->SetTitle("d#eta");
		hsig->GetYaxis()->SetTitle("d#phi");
		hsig->GetXaxis()->CenterTitle();
		hsig->GetYaxis()->CenterTitle();
		TH2D* mix_tmp = mixingTableMaker(hmix, doSmoothME);
		hsig->Divide(mix_tmp);
		if(doseagull) pullSeagull(hsig);
		return hsig;
}

TH2D* JTCSignalProducer::getSignal(TString name, bool doSeagullCorr){
	//	cout<<mix->GetName()<<endl;
	//	cout<<raw_sig->GetName()<<endl;
		raw_sig=(TH2D*) raw_sig->Clone("raw_"+name);
		raw_sig->Scale(1.0/raw_sig->GetXaxis()->GetBinWidth(1)/raw_sig->GetYaxis()->GetBinWidth(1)); //make the h2 invariant
		//	mix->SetName("mixing_"+name);
		sig = (TH2D*) raw_sig->Clone("signal_"+name);
		sig->GetXaxis()->SetTitle("d#eta");
		sig->GetYaxis()->SetTitle("d#phi");
		sig->GetXaxis()->CenterTitle();
		sig->GetYaxis()->CenterTitle();
		if(mix == 0) {
				cout<<"no mixing table has been assigned!"<<endl;
				return 0;
		}
		if(doSideBandMixing) {
				mix_normalized= sideBandMixingTableMaker(sig, sidebandmin, sidebandmax);
				//mix_normalized->SetName("sideBand_mixing_"+name);
		}
		else {
				mix_normalized= mixingTableMaker(mix, doSmoothME);
		}
		mix_normalized->SetName("smoothed_mixing_"+name);
		sig->Divide(mix_normalized);
		if(doSeagullCorr)	pullSeagull(sig);
		sig_step2 = (TH2D*) sig->Clone("sig_mix_corrected_"+name);
		bkg = (TH2D*) getV2Bkg(sig,sideMin , sideMax );
		bkg->SetName("bkg_"+name);
		sig->Add(sig, bkg, 1, -1);
		if(shiftSignal){
			   	getSignal_phiSideBand("sideBand");
				float mean = getMean(side_deta, -0.5, 0.5);
				cout<<"mean = "<<mean<<endl;
				cout<<"name : "<<name<<endl;
				shiftTH2D(sig, -mean);
		}
		return sig;
}

TH1D* JTCSignalProducer::doDrIntegral(TString name ){
		TString title = sig->GetTitle();
		TString hname = "dr_"+name;
		dr_integral = drDistMaker(sig, hname, title, 14, drbin );
		dr_integral->GetXaxis()->SetTitle("#Deltar");
		return dr_integral;
}

void JTCSignalProducer::WriteTH2(){
		raw_sig->Write();
		sig->Write();
		mix_normalized->Write();
		//mix->Write();
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

TH1* JTCSignalProducer::projX(bool doRebin, TH2D*h2, float x, float y, TString opt){
		// here h2 needs to be invariant 
		TH1* h=projectionX(h2, x, y, opt);
		h->Scale(h2->GetYaxis()->GetBinWidth(1));  
		if(doRebin){
				TString name = h->GetName();
				name = "rebined_"+name;
				TH1* hh=h;
				h=invariantRebin(h,name, 23, etabin);
				delete hh;
		}
		return h;
}

TH1* JTCSignalProducer::projY(bool doRebin, TH2D*h2, float x, float y, TString opt){
		// here h2 needs to be invariant 
		TH1* h=projectionY(h2, x, y, opt);
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
		sig_deta = (TH1D*) projX(rebin, sig, -1, 1, "e"); sig_deta->SetName(tmp);
		tmp = "signal_dphi_"+name;
		sig_dphi = (TH1D*) projY(rebin, sig, -1, 1, "e"); sig_dphi->SetName(tmp);

		tmp = "sig_step2_deta_"+name;
		sig_step2_deta = (TH1D*) projX(rebin, sig_step2, -1, 1, "e"); sig_step2_deta->SetName(tmp);
		tmp = "sig_step2_dphi_"+name;
		sig_step2_dphi = (TH1D*) projY(rebin, sig_step2, -1, 1, "e"); sig_step2_dphi->SetName(tmp);

		tmp = "bkg_deta_"+name;
		bkg_deta = (TH1D*) projX(rebin, bkg, -1, 1, "e"); bkg_deta->SetName(tmp);
		tmp = "bkg_dphi_"+name;
		bkg_dphi = (TH1D*) projY(rebin, bkg, -1, 1, "e"); bkg_dphi->SetName(tmp);

		tmp = "side_deta_"+name;
		side_deta = (TH1D*) projX(rebin, sig, sidebandmin, sidebandmax, "e"); side_deta->SetName(tmp);

		tmp = "side_deta_mix_"+name;
		side_deta_mix = (TH1D*) projX(rebin, mix_normalized, 1.2, 2.2, "e"); side_deta_mix->SetName(tmp);

		tmp = name;
		dr_integral = doDrIntegral(tmp); 
}


void JTCSignalProducer::WriteTH1(){
		sig_deta->Write(); sig_dphi->Write(); 
		sig_step2_deta->Write(); sig_step2_dphi->Write(); 
	   	bkg_deta->Write(); bkg_dphi->Write();
		side_deta->Write();	
		side_deta_mix->Write();
		dr_integral->Write();
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
		histStyle(h2);
		float mean = h2->GetBinContent(h1->FindBin(0.5));
		float dvt =  h2->GetBinError(h1->FindBin(0.5));
		//		float dvt =  fabs(0.05*h2->GetMaximum());
		h2->SetAxisRange(h2->GetMinimum()-5*dvt, h2->GetMaximum()+6*dvt, "Y");
		h1->SetLineWidth(1);
		h2->SetLineWidth(1);
		h1->SetLineColor(kBlack);
		h2->SetLineColor(kRed);
		h2->SetAxisRange(-2.5, 2.499, "X");
		h2->Draw("same");
		h1->Draw("same");
}


void JTCSignalProducer::drawSideBandCheck(){
		cout<<"drawing.."<<endl;
		TH1 *h1, *h2;
		h1 = sig_deta; 
		h2 = side_deta; 
		cout<<h2->GetName()<<endl;
		histStyle(h2);
		float mean = h2->GetBinContent(h1->FindBin(0));
		float dvt =  h2->GetBinError(h1->FindBin(0));
		h2->SetAxisRange(mean-6*dvt, mean+10*dvt, "Y");
		h1->SetLineWidth(1);
		h2->SetLineWidth(1);
		h1->SetLineColor(kBlack);
		h2->SetLineColor(kOrange+7);
		h2->SetAxisRange(-3., 2.99, "X");
		h2->Draw("same");
		h1->Draw("same");
		TLine* l = new TLine(); l->SetLineStyle(2);
		l->DrawLine(-2.5, 0, 2.5, 0);
}

void JTCSignalProducer::drawSliceSideBand(TString name){
		cout<<"drawing slice sideband ... "<<endl;
	TH1D* h[3];
	TString tmp = "side_deta1_"+name;
	h[0] = (TH1D*) projX(1, sig, 1.2, 1.4, "e"); h[0]->SetName(tmp);
	cout<<sig->GetName()<<endl;
	cout<<h[0]->GetName()<<endl;
	tmp = "side_deta2_"+name;
	h[1] = (TH1D*) projX(1, sig, 1.4, 1.6, "e"); h[1]->SetName(tmp);
	tmp = "side_deta3_"+name;
	h[2] = (TH1D*) projX(1, sig, 1.6, 1.8, "e"); h[2]->SetName(tmp);
	h[0]->SetLineColor(kAzure+7);
	h[1]->SetLineColor(kGreen+3);
	h[2]->SetLineColor(kOrange+7);
	h[0]->SetMarkerColor(kAzure+7);
	h[1]->SetMarkerColor(kGreen+3);
	h[2]->SetMarkerColor(kOrange+7);
	float mean = h[0]->GetBinContent(h[0]->FindBin(0));
	float dvt =  h[0]->GetBinError(h[0]->FindBin(0));
	h[0]->SetAxisRange(mean-10*dvt, mean+10*dvt, "Y");
	h[0]->SetAxisRange(-2.5, 2.499, "X");
	for(int i=0; i<3; ++i){
		histStyle(h[i]);
		h[i]->SetLineWidth(1);
		h[i]->Draw("same");
	}
		cout<<"finish drawing "<<endl;
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

TH1D* JTCSignalProducer::getSignal_phiSideBand(TString name){
		TString tmp = "side_deta_"+name;
		side_deta = (TH1D*) projX(1, sig, sidebandmin, sidebandmax, "e"); side_deta->SetName(tmp);
		return side_deta;
}

TH1D* JTCSignalProducer::getSignal_dEta(TString name){
		TString tmp = "signal_deta_"+name;
		sig_deta = (TH1D*) projX(1, sig, -1, 1, "e"); sig_deta->SetName(tmp);
		return sig_deta;
}

JTCSignalProducer::~JTCSignalProducer(){
		if(sig!=0){
				delete sig;
				delete sig_step2;  // right after the mixing correction
				delete mix_normalized;
				delete bkg;
		}
		if( dr_integral!= 0)delete dr_integral;
		if(sig_deta != 0 ){
				delete sig_deta;
				delete sig_dphi;
				delete sig_step2_deta;
				delete sig_step2_dphi;
				delete bkg_deta;
				delete bkg_dphi;
				delete side_deta;
				delete side_deta_mix;
		}
};
#endif 
