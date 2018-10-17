
#ifndef standardSequence_H
#define standardSequence_H
#include "config.h"

void addDr(TString name, TString prefix, TString path){
		auto f = TFile::Open(path, "update");
		auto sig = read_flatten<TH2D>(f, "signal_"+prefix+"_noCorr");
		auto sig_pt = read_flatten<TH2D>(f, "signal_"+prefix+"_pTweighted_noCorr");
		auto rs  = read_flatten<TH2D>(f, "sig_mix_corrected_"+prefix+"_noCorr");
		auto rs_pt = read_flatten<TH2D>(f, "sig_mix_corrected_"+prefix+"_pTweighted_noCorr");
		auto dr = getDr(name+"_sig_nc", sig);
		auto dr_pt = getDr(name+"_sig_pTweighted_nc", sig_pt);
		auto dr_rs_pt = getDr(name+"_rs_pTweighted_nc", rs_pt);
		auto dr_rs  = getDr(name+"_rs_nc", rs);
		dumpHists<TH1D>(dr);
		dumpHists<TH1D>(dr_pt);
		dumpHists<TH1D>(dr_rs_pt);
		dumpHists<TH1D>(dr_rs);
		f->Close();
}

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


TF1* fit0(TString name, TH1D* h, float xmin, float xmax){
		TF1 *f = new TF1(name, "[0]*x+[1]", 0, 1);
		h->Fit(f, "", "", xmin, xmax);
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
						tf[i+nPt*j] = fit0(name+Form("_%d_%d",i,j), h[i+nPt*j], 0, 0.2);
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
		if(!hist){
				if(isNumber)
						tagCorr	= read_flatten<TF1>(tagf, "taggedBias");
				else 
						tagCorr	= read_flatten<TF1>(tagf, "taggedBias_pTweighted");
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								cout<<i<<", "<<j<<": "<<h[i+nPt*j]->GetName()<<endl;
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
								//cout<<i<<", "<<j<<endl;
								cout<<i<<", "<<j<<": "<<h[i+nPt*j]->GetName()<<endl;
								h[i+nPt*j]->Divide(tagCorr_hist[i+nPt*j]);
						}
				}
		}
		//cout<<tagCorr[0]->GetName()<<endl;
		//divideTF1(h, tagCorr);
		//	tagf->Close();
		return;
}

TH1D** decontamination(TString name, TH1D** sig, TH1D** in, float pur){
		auto dr_co= copy<TH1D>("dr_co", in);
		scale(dr_co, 1-pur);
		TH1D** dr = new TH1D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						dr[i+nPt*j]=(TH1D*)sig[i+nPt*j]->Clone(name+Form("_%d_%d", i,j));
						//if(i==nPt-1) continue;
						dr[i+nPt*j]->Add(dr_co[i+nPt*j], -1);
						dr[i+nPt*j]->Scale(1.0/pur);
				}
		}
		free<TH1D>(dr_co);
		return dr;
}
TH2D** decontamination(TString name, TH2D** _sig, TH2D** in, float pur){
		auto sig_co= copy<TH2D>("sig_co", in);
		scale(sig_co, 1-pur);
		TH2D** sig = new TH2D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						sig[i+nPt*j]=(TH2D*)_sig[i+nPt*j]->Clone(name+Form("_%d_%d", i,j));
//						if(i==nPt-1) continue;
						sig[i+nPt*j]->Add(sig_co[i+nPt*j], -1);
						sig[i+nPt*j]->Scale(1.0/pur);
				}
		}
		free<TH2D>(sig_co);
		return sig;
}
void applyTkCorrection(TString name, TString tkEffFile){
		auto tkf = TFile::Open(dataDumpPath+tkEffFile);
		auto file = TFile::Open(dataDumpPath+name, "update");
		auto tkCorr= read_flatten<TH1D>(tkf, "trackEff_smoothed");
		auto h = read_flatten<TH1D>(file, "dr_residualCorrected");
		auto hpt = read_flatten<TH1D>(file, "dr_pTweighted_residualCorrected");
		cout<<"applying tk correction..."<<endl;
		auto h_corr= binary_operation<TH1D>("dr_final", h, tkCorr, "ratio");
		auto h_pt_corr= binary_operation<TH1D>("dr_final_pTweighted", hpt, tkCorr, "ratio");
		file->cd();
		dumpHists<TH1D>(h_corr);
		dumpHists<TH1D>(h_pt_corr);
		file->Close();
}

void applyResidual(TString name, TString resfname ){
		auto resf = TFile::Open(dataDumpPath+resfname);
		auto file = TFile::Open(dataDumpPath+name, "update");
		cout<<"adding residual correction..."<<endl;
		auto resCorr= read_flatten<TH1D>(resf, "bias");
		auto resCorr_pt= read_flatten<TH1D>(resf, "bias_pt");
		auto h = read_flatten<TH1D>(file, "dr_pur");
		auto hpt = read_flatten<TH1D>(file, "dr_pur_pTweighted");
		auto dr = binary_operation<TH1D>("dr_residualCorrected", h, resCorr, "ratio");
		auto dr_pt = binary_operation<TH1D>("dr_pTweighted_residualCorrected", h, resCorr_pt, "ratio");
		file->cd();
		dumpHists<TH1D>(dr);
		dumpHists<TH1D>(dr_pt);
		file->Close(); 
}

TH1D* get_eff(TString name, TFile *f ){
		TH1D* eff = (TH1D*) f->Get("taggedTrueBJet_corrpt_0")->Clone(name);
		TH1D* tmp = (TH1D*) f->Get("trueBJet_corrpt_0");
		eff->Divide(eff, tmp, 1, 1, "B");
		return eff;
}

TH1D* get_pur(TString name, TFile *f ){
		TH1D* eff = (TH1D*) f->Get("taggedTrueBJet_corrpt_0")->Clone(name);
		TH1D* tmp = (TH1D*) f->Get("taggedBJet_corrpt_0");
		eff->Divide(eff, tmp, 1, 1, "B");
		return eff;
}

TH1D* get_cont(TString name, TFile *f ){
		TH1D* eff = (TH1D*) f->Get("contJet_corrpt_0")->Clone(name);
		TH1D* tmp = (TH1D*) f->Get("inclJet_corrpt_0");
		eff->Divide(tmp);
		return eff;
}
#endif

