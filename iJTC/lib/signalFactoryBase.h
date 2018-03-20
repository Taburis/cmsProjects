


#ifndef signalFactoryBase_H
#define signalFactoryBase_H

#ifndef ROOT_TH1D
#include "TH1D.h"
#endif

#ifndef ROOT_TH2D
#include "TH2D.h"
#endif

#ifndef doublePad_H
#include "doublePad.h"
#endif


class signalFactoryBase {
		public : 
				signalFactoryBase(){}
				TH2D* mixingTableMaker(TH2D* mix, bool doSmooth = true);
				TH2D* getV2Bkg(TH2D* signal, float sideMin=1.5 ,float sideMax = 2.5);
				TH2D* signalMaker(TH2D* signalH2, TH2D* meH2, float sideMin, float sideMax, bool doSmoothME = true);
				void drIntegral(TH2D* signal, TH1D* drDist, bool isStatErr = true);
				TH1D* drDistMaker(TH2D* signal, TString name, TString title, int nbin, const float * bins);
				TH1D* drGeoTest(TH2D* signal, TH1D* drDist);
				void  doDrPhaseCorrection(TH2D* signal, TH1D* drDist);
				TH1* invariantRebin(TH1* h, TString name , int n, const Double_t * bins);
				TH2D* sideBandMixingTableMaker(TH2D* h2, float sidemin, float sidemax);
				TH1* projectionX(TH2D* h, float x, float y, TString opt="");
				TH1* projectionY(TH2D* h, float x, float y, TString opt="");
				void averageTH1(TH1* h, int n);
				TH1* meanProjectionX(TH2D* h, float x, float y);
				TH1* meanProjectionY(TH2D* h, float x, float y);
				void makeInvariant(TH2D*);

		public : 
				float midLeft = -0.15;
				float midRight=  0.15;
				TString dump_path = "";
				TH1D* hgeo = NULL;
};

TH1* signalFactoryBase::projectionX(TH2D* h, float x, float y, TString opt){
		int xbin = h->GetYaxis()->FindBin(x);
		int ybin = h->GetYaxis()->FindBin(y);
		ybin = h->GetYaxis()->FindBin(y-h->GetYaxis()->GetBinWidth(ybin));
		TString name = h->GetName(); name = "projX_"+name;
		return h->ProjectionX(name, xbin, ybin , opt);
}

TH1* signalFactoryBase::projectionY(TH2D* h, float x, float y, TString opt){
		int xbin = h->GetXaxis()->FindBin(x);
		int ybin = h->GetXaxis()->FindBin(y)-1;
		ybin = h->GetXaxis()->FindBin(y-h->GetXaxis()->GetBinWidth(ybin));
		TString name = h->GetName(); name = "projY_"+name;
		return h->ProjectionY(name, xbin, ybin ,opt);
}

TH2D* signalFactoryBase::mixingTableMaker(TH2D* mix, bool doSmooth){
		//make the mixing table invariant
		mix->Scale(1.0/mix->Integral()/mix->GetXaxis()->GetBinWidth(1)/mix->GetYaxis()->GetBinWidth(1));
		TH1D* temp = (TH1D*)mix->ProjectionX("_eta");
		TString name = mix->GetName();
		if( doSmooth) name = name+"_smoothed";
		else name = name+"_noEmptyBin";
		TH2D* ME = (TH2D*) mix->Clone(name);
		float mean=0;
		int binLeft = temp->FindBin(midLeft);
		int binRight= temp->FindBin(midRight)+1;
		for(int i=binLeft;i<binRight; i++){
				mean += temp->GetBinContent(i);
		}
		mean = mean /(temp->FindBin(midRight)-temp->FindBin(midLeft)+1)/mix->GetNbinsY();
		temp->Scale(1.0/mean);
		if(doSmooth){
				for(int ix=1; ix<mix->GetNbinsX()+1; ix++){
						for(int iy=1; iy<mix->GetNbinsY()+1; iy++){
								if( ix< binRight && ix>= binLeft){
										ME->SetBinContent(ix, iy, 1);
										ME->SetBinError(ix, iy, 0);
								}
								else{
										ME->SetBinContent(ix, iy, temp->GetBinContent(ix)/mix->GetNbinsY());
										ME->SetBinError(ix, iy, temp->GetBinError(ix)/sqrt(mix->GetNbinsY()));
								}
						}
				}
		}
		else { //only fill the empty bin of the mixing by the content we have
				ME->Scale(1.0/mean);
				for(int ix=1; ix<mix->GetNbinsX()+1; ix++){
						for(int iy=1; iy<mix->GetNbinsY()+1; iy++){
								if( ix< binRight && ix>= binLeft){
										ME->SetBinContent(ix, iy, 1);
										ME->SetBinError(ix, iy, 0);
								}
								else if(ME->GetBinContent(ix,iy)==0){
										ME->SetBinContent(ix, iy, temp->GetBinContent(ix)/mix->GetNbinsY());
										ME->SetBinError(ix, iy, temp->GetBinError(ix)/sqrt(mix->GetNbinsY()));
								}
						}
				}
		}
		return ME;
}
TH2D* signalFactoryBase::getV2Bkg(TH2D* signal, float sideMin,float sideMax){
		TString stemp = signal->GetName();
		TH2D* bkg = (TH2D*)signal->Clone(stemp+"_bkg");
		int outterRight= signal->GetXaxis()->FindBin(sideMax); 	
		int innerRight = signal->GetXaxis()->FindBin(sideMin); 	
		int outterLeft = signal->GetXaxis()->FindBin(-signal->GetXaxis()->GetBinCenter(outterRight)); 	
		int innerLeft  = signal->GetXaxis()->FindBin(-signal->GetXaxis()->GetBinCenter(innerRight)); 
		TH1D* temp = (TH1D*)bkg->ProjectionY("bkg_sideBandProj", innerRight, outterRight, "e");
		temp->Add(bkg->ProjectionY("", outterLeft, innerLeft, "e"));

		TH2D* aux_bkg = (TH2D*)signal->Clone("aux_bkg");
		for(int iy=1; iy<bkg->GetNbinsY()+1; iy++){
				for(int ix=1; ix<bkg->GetNbinsX()+1; ix++){
						aux_bkg->SetBinContent(ix, iy, 1);
				}
		}
		TH1D* aux_proj = (TH1D*)aux_bkg->ProjectionY("aux_proj", innerRight, outterRight);
		aux_proj->Add(aux_bkg->ProjectionY("", outterLeft, innerLeft));
		for(int iy=1; iy<bkg->GetNbinsY()+1; iy++){
				int nbin = int(aux_proj->GetBinContent(iy));
				for(int ix=1; ix<bkg->GetNbinsX()+1; ix++){
						bkg->SetBinContent(ix, iy, temp->GetBinContent(iy)/nbin);
						bkg->SetBinError(ix, iy, temp->GetBinError(iy)/sqrt(nbin));
				}
		}
		temp->Delete();
		aux_bkg->Delete();
		aux_proj->Delete();
		return bkg;
}

TH2D* signalFactoryBase::signalMaker(TH2D* signalH2, TH2D* meH2, float sideMin, float sideMax, bool doSmoothME){
		// the signal input has to be per-jet normalized before put it into this function
		// it should return the invariant signal after ME correction and bkg subtraction
		//
		TString hname = signalH2->GetName();
		TH2D* signal = (TH2D*)signalH2->Clone(hname+"_signal");
		signal->GetXaxis()->SetTitle("d#eta");
		signal->GetYaxis()->SetTitle("d#phi");
		//make the input 2D histogram invariant
		signal->Scale(1.0/signal->GetXaxis()->GetBinWidth(1)/signal->GetYaxis()->GetBinWidth(1));
		TH2D* mix= mixingTableMaker(meH2, doSmoothME);
		signal->Divide(mix);
		delete mix;
		mix = NULL;
		TH2D* bkg = (TH2D*) getV2Bkg(signal,sideMin , sideMax );
		signal->Add(signal, bkg, 1, -1);
		return signal;
}

void signalFactoryBase::drIntegral(TH2D* signal, TH1D* drDist, bool isStatError){
		// this should return the invariant dr distribution (dr bin width will be divided out ) 
		float content;
		float error;
		float width;
		float dr;
		float xwidth, ywidth;
		drDist->Sumw2();
		for(int i=1; i<drDist->GetNbinsX()+1;i++){
				drDist->SetBinContent(i, 0);
				drDist->SetBinError(i, 0);
		}
		for(int jx=0; jx<signal->GetNbinsX(); jx++){
				for(int jy=0; jy<signal->GetNbinsY(); jy++){
						dr = sqrt( pow(signal->GetXaxis()->GetBinCenter(jx),2) +\
										pow(signal->GetYaxis()->GetBinCenter(jy),2));
						xwidth = signal->GetXaxis()->GetBinWidth(jx);
						ywidth = signal->GetYaxis()->GetBinWidth(jy);
						// integrand f(x,y)dxdy
						content = signal->GetBinContent(jx,jy)*xwidth*ywidth;
						if( content ) {
								if(isStatError){
										error = sqrt(pow(drDist->GetBinError(drDist->FindBin(dr)),2)+\
														pow(signal->GetBinError(jx,jy)*xwidth*ywidth,2));
								}
								else {
										float err = signal->GetBinError(jx,jy)*xwidth*ywidth/content;
										//	if(err > 1) {
										//		error =drDist->GetBinError(drDist->FindBin(dr));
										//	}
										//	else {	
										error = sqrt(pow(drDist->GetBinError(drDist->FindBin(dr)),2)+\
														pow(err,2));
										//	}
								}
								drDist->Fill(dr, content);
								drDist->SetBinError(drDist->FindBin(dr), error);
						}
				}
		}
		// make the histogram invariant 
		for(int i=1; i<drDist->GetNbinsX()+1;i++){
				content = drDist->GetBinContent(i);
				error= drDist->GetBinError(i);
				width= drDist->GetBinWidth(i);
				drDist->SetBinContent(i, content/width);
				if( isStatError) drDist->SetBinError(i, error/width);
				else drDist->SetBinError(i, content*error/width);
		}
		return;
}

TH1D* signalFactoryBase::drDistMaker(TH2D* signal, TString name, TString title, int nbin, const float * bins){
		TH1D* drDist = new TH1D(name, title, nbin, bins);
		drIntegral(signal, drDist);
		return drDist;
}

TH1D* signalFactoryBase::drGeoTest(TH2D* signal, TH1D* drDist){
		TString temp = drDist->GetName();
		TH1D* drCounts = (TH1D*) drDist->Clone(temp+"_counts");
		TH1D* area_ideal = (TH1D*) drDist->Clone("ideal_phase");
		float dr;
		float da;
		for(int jx=1; jx<signal->GetNbinsX()+1; jx++){
				for(int jy=1; jy<signal->GetNbinsY()+1; jy++){
						dr = sqrt( pow(signal->GetXaxis()->GetBinCenter(jx),2) +\
										pow(signal->GetYaxis()->GetBinCenter(jy),2));
						da = (signal->GetXaxis()->GetBinWidth(jx))*(signal->GetYaxis()->GetBinWidth(jy));
						drCounts->Fill(dr,da); // get aera in the dr region
				}
		}
		for(int i=1; i<drCounts->GetNbinsX()+1; ++i){
				float rlow = area_ideal->GetBinLowEdge(i);
				float width = area_ideal->GetBinWidth(i);
				float rup = rlow + width;
				cout<<rup<<endl;
				//		drCounts->SetBinContent(i, drCounts->GetBinContent(i));
				area_ideal->SetBinContent(i, TMath::Pi()*(pow(rlow+width,2)-pow(rlow,2))); 
				area_ideal->SetBinError(i, 0); 
		}
		TH1D* geoCorr = (TH1D*) area_ideal ->Clone("geoCorr");
		geoCorr->Divide(drCounts);

		//plot the report map;
		TCanvas* c = new TCanvas("c_geoCorrection", "", 600,800);
		auto dp = new doublePad((TPad*)gPad);
		dp->addRatioPair(area_ideal, drCounts, kBlue, "2#pir", "grids/#Deltar");
		drCounts->SetMarkerStyle(20);
		drCounts->SetMarkerSize(0.8);
		area_ideal->SetLineColor(kRed);

		auto tl = new TLine();
		tl->SetLineStyle(2);
		dp->y2min=0.71; dp->y2max=1.29;
		//	dp->xmin=0; dp->xmax=0.99;
		dp->setXtitle("#Deltar");
		dp->setYtitle("#DeltaS/#Deltar");
		dp->Draw();
		dp->drawLegend();
		dp->putText("phase-space adjustment",0.4, 0.9, 2);
		dp->hpad->cd(2);
		tl->DrawLine(0, 1, 1, 1);
		dp->h1frame->GetYaxis()->SetLabelSize(0.06);
		dp->h1frame->GetYaxis()->CenterTitle();
		dp->h1frame->GetYaxis()->SetTitleSize(0.08);
		dp->h1frame->GetYaxis()->SetTitleOffset(0.5);
		dp->h2frame->GetYaxis()->SetLabelSize(0.05);
		dp->h2frame->GetXaxis()->SetLabelSize(0.05);
		dp->h2frame->GetXaxis()->CenterTitle();
		dp->h2frame->GetXaxis()->SetTitleSize(0.05);
		c->SaveAs(dump_path+"phaseSpaceCorrection.eps");
		return geoCorr;
}	

void signalFactoryBase::doDrPhaseCorrection(TH2D* signal, TH1D* h1){
		//const TArrayD* a = h1->GetXaxis()->GetXbins();
		//int nbin = a->GetSize();
		//Double_t* rbin = new Double_t(nbin);
		if( hgeo ==NULL ) {
				TH1D* htmp = (TH1D*) h1->Clone("_count");
				for(int k=1; k< htmp->GetNbinsX()+1; ++k){
						htmp->SetBinContent(k, 0);
						htmp->SetBinError(k, 0);
				}
				hgeo = (TH1D*) htmp->Clone("_phase_correction");
				float dr, da;
				for(int jx=1; jx<signal->GetNbinsX()+1; jx++){
						for(int jy=1; jy<signal->GetNbinsY()+1; jy++){
								dr = sqrt( pow(signal->GetXaxis()->GetBinCenter(jx),2) +\
												pow(signal->GetYaxis()->GetBinCenter(jy),2));
								da = (signal->GetXaxis()->GetBinWidth(jx))*(signal->GetYaxis()->GetBinWidth(jy));
								htmp->Fill(dr,da); // get aera in the dr region
						}
				}
				for(int i=1; i<hgeo->GetNbinsX()+1; ++i){
						float rlow = hgeo->GetBinLowEdge(i);
						float width = hgeo->GetBinWidth(i);
						float rup = rlow + width;
						//cout<<rup<<endl;
						//		hgeo->SetBinContent(i, hgeo->GetBinContent(i));
						hgeo->SetBinContent(i, TMath::Pi()*(pow(rlow+width,2)-pow(rlow,2))); 
				}
				hgeo->Divide(htmp);
				for(int i=1; i<hgeo->GetNbinsX()+1; ++i){
						hgeo->SetBinError(i, 0); 
				}
				delete htmp;
		}
		/*
		*/
		h1->Multiply(hgeo);
		return ;
}

TH1* signalFactoryBase::invariantRebin(TH1* h1, TString name , int n, const Double_t * bins){
		TH1* h=(TH1*) h1->Clone("tmp");
		//input h needs to be invariant 
		for(int i=1; i<h->GetNbinsX()+1; ++i){
				Double_t wd = h->GetBinWidth(i);
				h->SetBinContent(i, h->GetBinContent(i)*wd);
				h->SetBinError  (i, h->GetBinError(i)*wd);
		}
		TH1* hh = h->Rebin(n, name, bins);
		for(int i=1; i<hh->GetNbinsX()+1; ++i){
				Double_t wd = hh->GetBinWidth(i);
				hh->SetBinContent(i, hh->GetBinContent(i)/wd);
				hh->SetBinError  (i, hh->GetBinError(i)/wd);
		}
		delete h;
		return hh;
}

TH2D* signalFactoryBase::sideBandMixingTableMaker(TH2D* h2, float sidemin, float sidemax){
		int s1 = h2->GetYaxis()->FindBin(sidemin);
		int s2 = h2->GetYaxis()->FindBin(sidemax-0.001);
		int dbin = s2-s1;
		TH1D* temp = (TH1D*) h2->ProjectionX("_sideMix_deta", s1, s2);
		TString name = h2->GetName();
		TH2D* ME = (TH2D*) h2->Clone(name);
		float mean=0;
		int binLeft = temp->FindBin(midLeft);
		int binRight= temp->FindBin(midRight)+1;
		for(int i=binLeft;i<binRight; i++){
				mean += temp->GetBinContent(i);
		}
		mean = mean /(temp->FindBin(midRight)-temp->FindBin(midLeft)+1)/h2->GetNbinsY();
		temp->Scale(1.0/mean);
		for(int ix=1; ix<h2->GetNbinsX()+1; ix++){
				for(int iy=1; iy<h2->GetNbinsY()+1; iy++){
						if( ix< binRight && ix>= binLeft){
								ME->SetBinContent(ix, iy, 1);
								ME->SetBinError(ix, iy, 0);
						}
						else{
								ME->SetBinContent(ix, iy, temp->GetBinContent(ix)/dbin);
								ME->SetBinError(ix, iy, temp->GetBinError(ix)/sqrt(dbin));
						}
				}
		}
		return ME;
}

void signalFactoryBase::averageTH1(TH1* h, int nbin){
		h->Scale(1.0/nbin); float nn = pow(nbin,0.5);
		for(int i=1; i<h->GetNbinsX()+1; ++i){
				h->SetBinError(i, h->GetBinError(i)*nn);
		}
}

TH1* signalFactoryBase::meanProjectionX(TH2D* h, float x, float y){
		int xbin = h->GetYaxis()->FindBin(x);
		int ybin = h->GetYaxis()->FindBin(y);
		ybin = h->GetYaxis()->FindBin(y-h->GetYaxis()->GetBinWidth(ybin));
		TString name = h->GetName(); name = "meanProjX_"+name;
		TH1D* hh = h->ProjectionX(name, xbin, ybin , "e");
		averageTH1(hh, ybin-xbin+1);
		return hh;
}

TH1* signalFactoryBase::meanProjectionY(TH2D* h, float x, float y){
		int xbin = h->GetXaxis()->FindBin(x);
		int ybin = h->GetXaxis()->FindBin(y)-1;
		ybin = h->GetXaxis()->FindBin(y-h->GetXaxis()->GetBinWidth(ybin));
		TString name = h->GetName(); name = "meanProjY_"+name;
		TH1D* hh = h->ProjectionY(name, xbin, ybin , "e");
		averageTH1(hh, ybin-xbin+1);
		return hh;
}

void signalFactoryBase::makeInvariant(TH2D* h){
		h->Scale(1.0/h->GetXaxis()->GetBinWidth(1)/h->GetYaxis()->GetBinWidth(1)); //make the h2 invariant
}
#endif
