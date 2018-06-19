
#ifndef standardSequence_H
#define standardSequence_H
#include "config.h"

void pullSignal(TString fname, TH2D** raw_sig, TH2D** mixing, float sidemin=1.5, float sidemax=2.5, bool doSeagull=0){
		cout<<"pulling singal for "<<fname;
		JTCSignalProducer *sp1[nPt*nCent];
		TString trkbin [] = {"1", "2", "3", "4", "8", "12", "16", "400"};
		TString centbin [] = {"0", "30", "100"};
		cout<<".";
		TString tmp;
		for(int i=0; i<nPt ; ++i){
				cout<<".";
				for(int j=0; j<nCent ; ++j){
						tmp = "track pt in ["+trkbin[i]+", "+trkbin[i+1]\
							   +"), cent in ["+centbin[j]+", "+centbin[j+1]+"), jet pt in [120, 1000]";
						raw_sig[i+nPt*j]->SetTitle("");
						//raw_sig[i+nPt*j]->SetTitle(tmp);
						//cout<<raw_sig[i+nPt*j]->GetName()<<endl;
						mixing [i+nPt*j]->SetTitle(tmp);
						sp1[i+nPt*j] = new JTCSignalProducer(raw_sig[i+nPt*j], mixing[i+nPt*j]);
						sp1[i+nPt*j]->sideMin=sidemin; sp1[i+nPt*j]->sideMax=sidemax;
						sp1[i+nPt*j]->getSignal(fname+Form("_%d_%d", i, j), doSeagull);
						sp1[i+nPt*j]->WriteTH2();
						delete sp1[i+nPt*j];
						//cout<<sp2[i+nPt*j]->sig->GetName()<<endl;
				}
		}
		cout<<". Done"<<endl;;
}

void pullSignal(TString fname, TString name, histCase &h, float sidemin=1.5, float sidemax=2.5, bool doSeagull=0){
		TFile * wf = new TFile(dataDumpPath+fname,"recreate");
		pullSignal(name, h.sig, h.mixing, sidemin, sidemax, doSeagull);
		pullSignal(name+"_pTweighted", h.sig_pTweighted, h.mixing, sidemin, sidemax, doSeagull);
		pullSignal(name+"_noCorr", h.sig_raw, h.mixing_raw, sidemin, sidemax, doSeagull);
		pullSignal(name+"_pTweighted_noCorr", h.sig_pTweighted_raw, h.mixing_raw, sidemin, sidemax, doSeagull);
		wf->Close();
		cout<<"signal dumped to "<<dataDumpPath+fname<<endl;
}


TF1* fit0(TString name, TH1D* h){
		TF1 *f = new TF1(name, "x*[1]+[0]", 0, 1);
		h->Fit(f, "", "", 0, 1);
		return f;
}


TF1* fit1(TString name, TH1D* h){
		// a+bx  0.3  c+dx+ex^2
		TF1 *f = new TF1(name, "(x<0.3)*([0]+x*[1])+ (x>0.3)*([2]+[3]*x+[4]*x*x)", 0, 1);
		h->Fit(f, "", "", 0, 1);
		return f;
}

TF1* fit2(TString name, TH1D* h){
		// a+bx  0.3  c+dx
		TF1 *f = new TF1(name, "(x<0.3)*([0]+x*[1])+ (x>0.3)*([2]+[3]*x)", 0, 1);
		h->Fit(f, "", "", 0, 1);
		return f;
}

TF1* fit3(TString name, TH1D* h){
		// a+bx+cx^2  0.4  d
		TF1 *f = new TF1(name, "(x<0.4)*([0]+x*[1]+x*x*[2])+ (x>0.4)*([3])", 0, 1);
		h->Fit(f, "", "", 0, 1);
		return f;
}
TF1* fit4(TString name, TH1D* h){
		// a+bx+cx^2  0.3  d
		TF1 *f = new TF1(name, "(x<0.3)*([0]+x*[1]+x*x*[2])+ (x>0.3)*([3])", 0, 1);
		h->Fit(f, "", "", 0, 1);
		return f;
}
TF1* fit5(TString name, TH1D* h){
		// a+bx  0.3  c
		TF1 *f = new TF1(name, "(x<0.3)*([0]+x*[1])+ (x>0.3)*([3])", 0, 1);
		h->Fit(f, "", "", 0, 1);
		return f;
}

TF1** doFitting(TString name, TH1D** h){
		TF1** tf = new TF1*[nPt*nCent];
		auto tx = new TLatex();  tx->SetTextSize(.08);
		auto mc = new mCanvasLoose("c_"+name, "", 2, 3);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						mc->drawHist(h[i+nPt*j], i+1);
						mc->cd(i+1);
						TString tmp = track_label[i];
						tx->DrawLatexNDC(0.02,0.93, tmp); 
						if(i==2)
								tf[i+nPt*j] = fit4(name+Form("_%d_%d",i,j), h[i+nPt*j]);
						else if(i==3) 
								tf[i+nPt*j] = fit5(name+Form("_%d_%d",i,j), h[i+nPt*j]);
						else if(i==4) 
								tf[i+nPt*j] = fit4(name+Form("_%d_%d",i,j), h[i+nPt*j]);
						else if(i<4)
								tf[i+nPt*j] = fit4(name+Form("_%d_%d",i,j), h[i+nPt*j]);
						else if(i<5)
								tf[i+nPt*j] = fit3(name+Form("_%d_%d",i,j), h[i+nPt*j]);
						else 
								tf[i+nPt*j] = fit3(name+Form("_%d_%d",i,j), h[i+nPt*j]);
				}
		}
		mc->SaveAs(figDumpPath+name+"plot_"+name+".gif");
		return tf;
}

void divideTF1(TH1D** h, TF1** f){
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						//						cout<<"i = "<<i<<", j = "<<j<<endl;
						h[i+nPt*j]->Divide(f[i+nPt*j]);
				}
		}
}

void applyBiasCorrection(TH1D** h, bool isNumber =1, bool hist=0){
		cout<<"adding tagger bias correction..."<<endl;
		auto tagf = TFile::Open(dataDumpPath+"ppCSVv1taggerBias2.root");
		TF1 **tagCorr;
		TH1D **tagCorr_hist;
		if(hist){
				if(isNumber)
						tagCorr	= read_flatten<TF1>(tagf, "taggedBias");
				else 
						tagCorr	= read_flatten<TF1>(tagf, "taggedBias_pTweighted");
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								cout<<i<<", "<<j<<endl;
								h[i+nPt*j]->Divide(tagCorr[i+nPt*j]);
						}
				}
		} else {
				if(isNumber)
						tagCorr_hist = read_flatten<TH1D>(tagf, "taggerBias_hist");
				else 
						tagCorr_hist = read_flatten<TH1D>(tagf, "taggerBias_hist_pTweighted");
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								cout<<i<<", "<<j<<endl;
								h[i+nPt*j]->Divide(tagCorr_hist[i+nPt*j]);
						}
				}
		}
//cout<<tagCorr[0]->GetName()<<endl;
//divideTF1(h, tagCorr);
//	tagf->Close();
return;
}

TH1D** applyTkCorrection(TString name, TH1D** h, bool isNumber =1){
		auto tkf = TFile::Open(dataDumpPath+"ppCSVv1TkCorr.root");
		TH1D **tkCorr;
		if(isNumber) tkCorr= read_flatten<TH1D>(tkf, "tkCorr");
		else tkCorr= read_flatten<TH1D>(tkf, "tkCorr_pTweighted");
		cout<<"adding tk correction..."<<endl;
		auto trk = binary_operation<TH1D>(name, h, tkCorr, "ratio");
		//	tkf->Close();
		return trk;
}

TH1D** applyResidual(TString name, TH1D** h, bool isNumber =1){
		cout<<"adding residual correction..."<<endl;
		auto resf = TFile::Open(dataDumpPath+"ppCSVv1Residual.root");
		if(isNumber) resCorr= read_flatten<TH1D>(resf, "residual");
		else resCorr= read_flatten<TH1D>(resf, "residual_pTweighted");
		auto dr = binary_operation<TH1D>(name, h, resCorr, "ratio");
		//	resf->Close();
		return dr;
}

TH2D** decontamination(inputSet &iset, float purity, bool isNumber = 1){
		cout<<"decontamination processing.."<<endl;
		readInput(iset);
		TString endname = tpname(iset);
		TH2D** rs_tg, **rs_in;
		if(isNumber){
				rs_tg = doMixingCorr("step2_sig_taggedB", iset.tg.sig_raw, iset.tg.mixing_raw, 0);
				rs_in = doMixingCorr("step2_sig_incl", iset.in.sig_raw, iset.in.mixing_raw, 1);
		}
		else {
				rs_tg = doMixingCorr("step2_sig_taggedB_pTweighted", iset.tg.sig_pTweighted_raw, iset.tg.mixing_raw, 0);
				rs_in = doMixingCorr("step2_sig_incl_pTweighted", iset.in.sig_pTweighted_raw, iset.in.mixing_raw, 1);
		}

		scale<TH2D>(rs_in, 1-purity);
		//auto rs_pur = binary_operation<TH2D>("raw_sig_pur", rs_tg, rs_in, "diff");
		TH2D** rs_pur = new TH2D*[nPt*nCent]; 
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<i<<", "<<j<<endl;
						rs_pur[i+nPt*j]=(TH2D*) rs_tg[i+nPt*j]->Clone(Form("raw_sig_pur_%d_%d",i,j));
						if(i<nPt-1){ rs_pur[i+nPt*j]->Add(rs_in[i+nPt*j], -1);
								rs_pur[i+nPt*j]->Scale(1.0/purity);
						}
				}
		}

		//		scale(rs_pur, 1.0/purity);

		if(isNumber) endname+="_pTweighted";
		auto sig_pur = doBkgSub("signal_"+endname, rs_pur);

		cout<<"decontamination done"<<endl;
		return sig_pur;
}
#endif

