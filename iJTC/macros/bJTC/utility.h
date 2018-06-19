
//#include "../../lib/import_pf_config.h"
#include "../../lib/JTCSignalProducer.h"
#include "../../lib/stackHist.h"
//#include "../../lib/JTCSkimer.h"
#ifndef xPlotStyle_H
#include "../../lib/xPlotStyle.h"
#endif
#ifndef utility_H
#define utility_H
#endif

TString dataDumpPath = "/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/myCorrelation/trunk/";
TString figDumpPath  = "/Users/tabris/cmsProjects/iJTC/macros/bJTC/fig_newWay/";
TString track_label[] = {"1 < p_{T}^{track} < 2 GeV",
		"2 < p_{T}^{track} < 3 GeV", "3 < p_{T}^{track} < 4 GeV",
		"4 < p_{T}^{track} < 8 GeV", "8 < p_{T}^{track} < 12 GeV", 
		"p_{T}^{track} > 12 GeV"};
TString cent_label[]={"Cent. 0-30%", "Cent. 30-100%"};
Color_t color_vec[6] = {kBlue+1, kRed+1, kGreen+2, kAzure+7, kMagenta+2, kBlack};
const int nPt = 6;
const int nCent=2;


//TLatex* tex = new TLatex(); 
//tex->SetTextSize(.08);


struct histCase{
		//if want to add any hist, need to add it in the quickRegHist as well, and add the filling in the fillCase
		TH2D** sig;
		TH2D** sig_pTweighted;
		TH2D** mixing;
		TH2D** sig_raw;
		TH2D** sig_pTweighted_raw;
		TH2D** mixing_raw;
		TH1D** jet_corrpt;
		TH1D** jet_eta;
		TH1D** jet_phi;
};

void initHistCase(histCase &hc){
		hc.sig=new TH2D*[nPt*nCent];
		hc.sig_pTweighted=new TH2D*[nPt*nCent];
		hc.mixing=new TH2D*[nPt*nCent];
		hc.sig_raw=new TH2D*[nPt*nCent];
		hc.sig_pTweighted_raw=new TH2D*[nPt*nCent];
		hc.mixing_raw=new TH2D*[nPt*nCent];
		hc.jet_corrpt=new TH1D*[nCent];
		hc.jet_eta=new TH1D*[nCent];
		hc.jet_phi=new TH1D*[nCent];
}

histCase readHistCase(TString name, TString cap, TFile* f){
		histCase hc;
		initHistCase(hc);
		for(int j=0; j<nCent; ++j){
				hc.jet_corrpt[j] =(TH1D*)f->Get(cap+Form("_corrpt_%d", j));
				hc.jet_eta[j] =(TH1D*)f->Get(cap+Form("_eta_%d", j));
				hc.jet_phi[j] =(TH1D*)f->Get(cap+Form("_phi_%d", j));
		}
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<name+Form("_%d_%d", i, j)<<endl;
						hc.sig[i+j*nPt] =(TH2D*)f->Get(name+Form("_%d_%d", i, j));
						hc.sig_pTweighted[i+j*nPt] =(TH2D*)f->Get(name+Form("_pTweighted_%d_%d", i, j));
						hc.mixing[i+j*nPt] =(TH2D*)f->Get(name+Form("_mixing_%d_%d", i, j));

						hc.sig_raw[i+j*nPt] =(TH2D*)f->Get(name+Form("_noCorr_%d_%d", i, j));
						hc.sig_pTweighted_raw[i+j*nPt] =(TH2D*)f->Get(name+Form("_pTweighted_noCorr_%d_%d", i, j));
						hc.mixing_raw[i+j*nPt] =(TH2D*)f->Get(name+Form("_mixing_noCorr_%d_%d", i, j));

						hc.sig[i+j*nPt]->Scale(1.0/hc.jet_corrpt[j]->Integral());
						hc.sig_pTweighted[i+j*nPt]->Scale(1.0/hc.jet_corrpt[j]->Integral());
						hc.sig_raw[i+j*nPt]->Scale(1.0/hc.jet_corrpt[j]->Integral());
						hc.sig_pTweighted_raw[i+j*nPt]->Scale(1.0/hc.jet_corrpt[j]->Integral());
				}
		}
		return hc;
}

template<typename T>
T** flatten(T* h[nPt][nCent]){
		T** hreturn = new T*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						hreturn[i+nPt*j] = h[i][j];
				}
		}
		return hreturn;
}

template<typename T>
T* index(T** h, int ptb, int centb){
		return h[ptb+nPt*centb];
}

template<typename T>
T** read_flatten(TFile* f, TString hcap){
		T** h = new T*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h[i+nPt*j] = (T*) f->Get(hcap+Form("_%d_%d", i, j));
						//cout<<hcap+Form("_%d_%d", i, j)<<endl;;
						//cout<<h[i+nPt*j] ->GetName()<<endl;;
				}
		}
		return h;
}

template<typename T>
T** read_flatten(TString name, TString hcap){
		TFile* f = new TFile(name);
		return read_flatten<T>(f, hcap);
}

template<typename T>
void scale(T** h, float c){
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h[i+nPt*j]->Scale(c);
				}
		}
}

template<typename T>
T** binary_operation(TString name, T** h1, T** h2, TString opt, bool ispp=0){
		int ncent = ispp ? 1:nCent;
		T** h = new T*[nPt*nCent]; 
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
						cout<<i<<", "<<j<<endl;
						h[i+nPt*j]=(T*) h1[i+nPt*j]->Clone(name + Form("_%d_%d",i,j));
						if(opt == "ratio") 
								h[i+nPt*j]->Divide(h2[i+nPt*j]);
						else if( opt== "add") 
								h[i+nPt*j]->Add(h2[i+nPt*j]);
						else if( opt== "multiply") 
								h[i+nPt*j]->Multiply(h2[i+nPt*j]);
						else if( opt== "diff") {
								h[i+nPt*j]->Add(h2[i+nPt*j], -1);
						}
						else if( opt== "binomialRatio") {
								h[i+nPt*j]->Divide(h[i+nPt*j], h2[i+nPt*j], 1, 1, "B");
						}
						else cout<<"no defined operation: "<<opt<<endl;
				}
		}
		return h;
}

TH1D** getDr(TString name, TH2D** sig){
		JTCSignalProducer sp;
		TH1D** dr = new TH1D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<sig[i+nPt*j]->GetName()<<endl;
						//									cout<<i+j*nPt<<endl;
						sp.sig=sig[i+nPt*j];
						dr[i+j*nPt]= sp.doDrIntegral(name+Form("_%d_%d",i, j));
						sp.doDrPhaseCorrection(sig[i+nPt*j], dr[i+j*nPt]);
						sp.dr_integral=0;
						sp.sig=0;
						//cout<<dr[i+j*nPt]->GetName()<<endl;
						//cout<<&dr[i+j*nPt]<<endl;;
				}
		}
		return dr;
}

void setXtitle(TH1** h, TString title){
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h[i+nPt*j]->GetXaxis()->SetTitle(title);
				}
		}
}


TH1D** projX(TString name, TH2D** sig){
		JTCSignalProducer sp;
		TH1D** projx = new TH1D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<i<<", "<<j<<endl;
						projx[i+nPt*j] =(TH1D*) sp.projX(1, sig[i+nPt*j], -1, 1);
						//					projx[i+nPt*j]->Rebin(4);
						//					projx[i+nPt*j]->Scale(0.25);
						//cout<<sig[i+nPt*j]->GetName()<<endl;
						//cout<<projx[i+nPt*j]->GetName()<<endl;
				}
		}
		return projx;
}

int index(int i, int j){
		return i+j*nPt;
}

void showColz(TString name, TH2D** sig ){
		auto cm = new mCanvasLoose("c_"+name, "", 2, 3, 1000, 1000);
		int j=0; 
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<nCent ; ++j){
						sig[i+nPt*j]->SetAxisRange(-1, 0.99, "X");
						sig[i+nPt*j]->SetAxisRange(-1, 0.99, "Y");
						cm->drawHist(sig[i+nPt*j], i+1, "colz");
				}
		}
		cm->SaveAs(figDumpPath+name+".gif");
}

TCanvas* showPlot(TString name, bool isHI , float line, float x1, float x2, float y1, float y2, int n_args, ...){
		cout<<"start drawing...."<<endl;
		va_list ap;
		va_start(ap, n_args);
		int ncol = 3, nrow = 2;
		if( isHI ) {ncol= 6; nrow = 2;}
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol);
		auto tx = new TLatex();  tx->SetTextSize(.06);
		auto tl = new TLine(); tl->SetLineStyle(2);
		int ncent = isHI ? 2: 1;
		TH1D** h;
		TString tmp;
		float min[nPt*nCent];
		float max[nPt*nCent];
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<ncent; ++j){
						min[i+nPt*j] = 0;
						max[i+nPt*j] = 0;
				}
		}
		for(int k=0; k<n_args; ++k){
				h = va_arg(ap, TH1D**);
				for(int i=0; i<nPt ; ++i){
						for(int j=0; j<ncent; ++j){
								float holder = h[index(i,j)]->GetMaximum();
								if(max[i+j*nPt]< holder)  max[i+j*nPt] = holder ;
								holder = h[index(i,j)]->GetMinimum();
								if(min[i+j*nPt]> holder)  min[i+j*nPt] = holder ;
						}
				}
		}

		va_start(ap, n_args);
		for(int k=0; k<n_args; ++k){
				h = va_arg(ap, TH1D**);
				for(int i=0; i<nPt ; ++i){
						for(int j=0; j<ncent; ++j){
								//cout<<i<<", "<<j<<endl;
								float grid = (max[i+j*nPt]-min[i+j*nPt])/20;

								cout<<index(i,j)<<endl; h[index(i,j)]->SetTitle("");
								h[index(i,j)]->SetLineColor(color_vec[k]);
								h[index(i,j)]->GetXaxis()->SetNdivisions(505);
								h[index(i,j)]->SetAxisRange(x1, x2,"X");
								if(y1<y2) h[index(i,j)]->SetAxisRange(y1, y2,"Y");
								else h[index(i,j)]->SetAxisRange(min[i+j*nPt]-grid, max[i+j*nPt]+grid,"Y");
								if(!isHI ) { 
										cm->drawHist(h[index(i,j)], i+1);
										cm->cd(i+1);
										//cout<<i+1<<endl;
										tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
										tl->DrawLine(x1, line, x2, line);
								}
								else { 
										cm->drawHist(h[index(i,j)], j+1, i+1);
										cout<<"i = "<<i<<", j = "<<j<<": "<<h[index(i,j)]->GetName()<<endl;
										cm->CD(j+1, i+1);
										//gPad->SetLogy();
										tmp = track_label[i]+"; "+cent_label[j];
										tx->DrawLatexNDC(0.1, 0.93, tmp);
								}
						}
				}
		}
		cm->SaveAs(figDumpPath+name);
		return (TCanvas*) cm;
}


void checkBkg(TString name, int ispb = 0){
		ispb++;
		TString tmp =dataDumpPath+name+"_JTCSignal.root";
		TFile *f = TFile::Open(tmp);
		JTCSignalProducer *sp1[8][2];
		JTCSignalProducer *sp2[8][2];
		auto tx = new TLatex();  tx->SetTextSize(.08);
		auto cp1 = new mCanvasLoose("bkgCheck_"+name, "", nPt, ispb);
		auto cp2 = new mCanvasLoose("bkgCheck_JS_"+name, "", nPt, ispb);
		auto cp3 = new mCanvasLoose("bkgCheck_sideBand"+name, "", nPt, 2*ispb);
		auto cp4 = new mCanvasLoose("bkgCheck_sideBand_JS"+name, "", nPt, 2*ispb);
		int ncol = 2*ispb;
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j< ispb; ++j){
						if(ispb-1)  tmp = track_label[i]+", "+cent_label[j];
						else	tmp = track_label[i];
						sp1[i][j] = new JTCSignalProducer();
						sp2[i][j] = new JTCSignalProducer();
						sp1[i][j]->read(f, name+Form("_%d_%d", i, j));
						sp2[i][j]->read(f, name+Form("_pTweighted_%d_%d", i, j));
						sp1[i][j]->getAllProj(name+Form("_%d_%d", i, j), 1);
						sp2[i][j]->getAllProj(name+Form("_pTweighted_%d_%d", i, j), 1);
						cp1->CD(i+1, ispb-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp1[i][j]->drawBkgCheck();
						tx->DrawLatexNDC(0.02,0.93, tmp); 
						/*
						   cp1->CD(i+1, ncol-j);
						   gPad->SetTopMargin(0.1);
						   gPad->SetBottomMargin(0.18);
						   sp1[i][j]->drawBkgCheck(0);
						   tx->DrawLatexNDC(0.02,0.93, tmp); 
						   */

						cp2->CD(i+1, ispb-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp2[i][j]->drawBkgCheck();
						tx->DrawLatexNDC(0.02,0.93, tmp); 
						/*
						   cp2->CD(i+1, ncol-j);
						   gPad->SetTopMargin(0.1);
						   gPad->SetBottomMargin(0.18);
						   sp2[i][j]->drawBkgCheck(0);
						   tx->DrawLatexNDC(0.02,0.93, tmp); 
						   */
						cp3->CD(i+1, ispb-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp1[i][j]->drawSideBandCheck();
						tx->DrawLatexNDC(0.02,0.93, tmp); 

						cp3->CD(i+1, ncol-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp1[i][j]->histStyle(sp1[i][j]->sig_deta);
						sp1[i][j]->sig_deta->SetAxisRange(-2.5, 2.499, "X");
						sp1[i][j]->sig_deta->Draw();
						sp1[i][j]->side_deta->Draw("same");
						tx->DrawLatexNDC(0.02,0.93, tmp); 

						cp4->CD(i+1, ispb-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp2[i][j]->drawSideBandCheck();
						tx->DrawLatexNDC(0.02,0.95, tmp); 

						cp4->CD(i+1, ncol-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp2[i][j]->histStyle(sp2[i][j]->sig_deta);
						sp2[i][j]->sig_deta->SetAxisRange(-2.5, 2.499, "X");
						sp2[i][j]->sig_deta->Draw();
						sp2[i][j]->side_deta->Draw("same");
						tx->DrawLatexNDC(0.02,0.95, tmp); 
				}
		}
		cp1->SaveAs(figDumpPath+"bkgCheck_"+name+".pdf");
		cp2->SaveAs(figDumpPath+"bkgCheck_JS_"+name+".pdf");
		cp3->SaveAs(figDumpPath+"bkgCheck_sideBand_"+name+".pdf");
		cp4->SaveAs(figDumpPath+"bkgCheck_sideBand_JS_"+name+".pdf");
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<ispb ; ++j){
						delete sp1[i][j];
						delete sp2[i][j];
				}
		}
}

template <typename T>
void dumpHists(T** h, bool ispp = 0){
		int ncent = ispp ? 1 : nCent;
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<ncent ; ++j){
						h[i+nPt*j]->Write();
				}
		}
}

template <typename T>
void dumpHists(TString path, T** h){
		auto wf = new TFile(path, "recreate");
		wf->cd();	
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<nCent ; ++j){
						h[i+nPt*j]->Write();
				}
		}
		wf->Close();
}

template <typename T>
void free(T** h){
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<nCent ; ++j){
						delete h[i+nPt*j];
				}
		}
		h=0;
}

template <typename T>
T** copy(TString name, T** h){
		T** hh = new T*[nPt*nCent];
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<nCent ; ++j){
						hh[i+nPt*j]=(T*)h[i+nPt*j]->Clone(name+Form("_%d_%d",i,j));
				}
		}
		return hh;
}

TH2D** doMixingCorr(TString name, TH2D** h, TH2D** hmix, bool doseagull){
		TH2D** sig_corr = new TH2D*[nPt*nCent];
		JTCSignalProducer sp;
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<h[i+j*nPt]->GetName()<<endl;
						sig_corr[i+j*nPt]=sp.mixingCorr(name+Form("_%d_%d",i,j), h[i+j*nPt], hmix[i+j*nPt], doseagull);
				}
		}
		return sig_corr;
}

TH2D** doBkgSub(TString name, TH2D** h){
		JTCSignalProducer sp;
		TH2D** sig_corr = new TH2D*[nPt*nCent];
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<h[i+j*nPt]->GetName()<<endl;
						sig_corr[i+j*nPt]=sp.doBkgSubtraction(name+Form("_%d_%d",i,j), h[i+j*nPt]);
						cout<<sig_corr[i+j*nPt]->GetName()<<endl;
				}
		}
		return sig_corr;
}

TH2D** readSig(TFile* f, TString cap, TString cap0){
		auto hh = read_flatten<TH2D>(f, cap);
		for(int j=0; j<nCent; ++j){
				TH1D* hjet = (TH1D*) f->Get(Form(cap0+"_corrpt_%d",j));
				for(int i=0; i<nPt; ++i){
						hh[i+j*nPt]->Scale(1.0/hjet->Integral());
				}
		}
		return hh;
}

void checkSideBand(TString name, TH2D** hsig, int ispb = 0){
		JTCSignalProducer* sp[nPt][nCent];
		auto tx = new TLatex();  tx->SetTextSize(.08);
		int ncol = 3, nrow = 2;
		auto cp1 = new mCanvasLoose("SideBandCheck_"+name, "", nrow, ncol);
		auto cp2 = new mCanvasLoose("SigdEta_"+name, "", nrow, ncol);
		TString tmp;
		auto tl = new TLine(); tl->SetLineStyle(2);
		int ncent = ispb ? 4 : 1;
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j< ncent; ++j){
						if(ispb-1)  tmp = track_label[i]+", "+cent_label[j];
						else	tmp = track_label[i];
						//cout<<tmp<<endl;
						sp[i][j] = new JTCSignalProducer();
						sp[i][j]->sig=hsig[i+nPt*j];
						sp[i][j]->getSignal_phiSideBand(name+Form("_phiSide_%d_%d",i,j));
						sp[i][j]->getSignal_dEta(name+Form("_sigDeta_%d_%d",i,j));
						cp1->cd(i+1);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp[i][j]->drawSideBandCheck();
						tx->DrawLatexNDC(0.02,0.93, tmp); 

						cp2->drawHist(sp[i][j]->sig_deta, i+1);
						cp2->cd(i+1);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						tl->DrawLine(-2.5, 0, 2.499, 0);
						tx->DrawLatexNDC(0.02,0.93, tmp); 

				}
		}
		cp1->SaveAs(figDumpPath+"bkgCheck_"+name+".pdf");
		cp2->SaveAs(figDumpPath+"sig_deta_"+name+".pdf");
}

TH2D** fracCorr(TString name, TH2D** ha, TH2D** hb){
		// return 2-ha/hb; 
		TH2D** ratio = binary_operation<TH2D>("ratio_temp", ha, hb, "ratio");
		TH2D** h = new TH2D*[nPt*nCent];
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<nCent ; ++j){
						h[i+nPt*j]=(TH2D*)hb[i+j*nPt]->Clone(name+Form("_%d_%d",i,j));
						for(int k=1; k<h[i+nPt*j]->GetNbinsX()+1; ++k){
								for(int l=1; l<h[i+nPt*j]->GetNbinsY()+1; ++l){
										float content =  ratio[i+nPt*j]->GetBinContent(k,l);
										if(content == 0) 
												h[i+nPt*j]->SetBinContent(k,l,1);
										else {
												h[i+nPt*j]->SetBinContent(k,l,2-content);
												h[i+nPt*j]->SetBinError  (k,l,ratio[i+nPt*j]->GetBinError(k,l));
										}
								}
						}
				}
		}
		free<TH2D>(ratio);
		return h;
}

TCanvas* showStack(TString name, bool isHI , float line, float x1, float x2, float y1, float y2, int n_args, ...){
		cout<<"start drawing...."<<endl;
		va_list ap;
		va_start(ap, n_args);
		int ncol = 3, nrow = 2;
		if( isHI ) {ncol= 4; nrow = 9;}
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		auto tl = new TLine(); tl->SetLineStyle(2);
		int ncent = isHI ? 2: 1;
		TH1D** h;
		TString tmp;
		THStack *hs[nPt*nCent]; 
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<ncent; ++j){
						hs[i+nPt*j]=new THStack(Form("hs%d_%d",i, j), "");
				}
		}

		va_start(ap, n_args);
		for(int k=0; k<n_args; ++k){
				h = va_arg(ap, TH1D**);
				for(int i=0; i<nPt ; ++i){
						for(int j=0; j<ncent; ++j){
								cout<<index(i,j)<<endl; h[index(i,j)]->SetTitle("");
								if(y1<y2) h[index(i,j)]->SetAxisRange(y1, y2,"Y");
								h[index(i,j)]->SetAxisRange(x1, x2,"X");
								h[index(i,j)]->SetLineColor(color_vec[k]);
								h[index(i,j)]->SetFillColor(color_vec[k]);
								hs[i+nPt*j]->Add(h[index(i,j)]); 
						}
				}
		}
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<ncent; ++j){
						if(!isHI ) { 
								cm->cd(i+1);
								hs[i+nPt*j]->Draw(); 
								//cout<<i+1<<endl;
								tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
								//tl->DrawLine(x1, line, x2, line);
						}
						else { 
								cout<<"i = "<<i<<", j = "<<j<<": "<<h[index(i,j)]->GetName()<<endl;
								cm->CD(i+1, 4-j);
								gPad->SetLogy();
								tmp = track_label[i]+" "+cent_label[j];
								//tx->DrawLatexNDC(0.1, 0.93, tmp);
						}
						hs[index(i,j)]->GetXaxis()->SetNdivisions(505);
						//						hs[i+nPt*j]->Draw(); 
				}
		}
		cm->SaveAs(figDumpPath+name+".gif");
		return (TCanvas*) cm;
}

void showClosure(TString name, float x1, float x2, float y1, float y2, TString l1, TString l2, TH1D** h1, TH1D** h2, TString opt="binomialRatio", bool isHI = 0){
		float cmax;
		float cmin;
		int ncol = 3, nrow = 2;
		if( isHI ) {ncol= 6; nrow = 2;}
		auto tl = new TLine(); tl->SetLineStyle(2);
		auto leg = new TLegend(.55, .57, .95, .77); leg->SetLineColor(0);
		auto ratio = binary_operation<TH1D>(name+"_ratio", h1, h2, opt);
		auto tx = new TLatex();  tx->SetTextSize(.07);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cmax=h1[index(i,j)]->GetMaximum();
						float holder = h2[index(i,j)]->GetMaximum();
						if(cmax< holder)  cmax = holder ;
						cmin=h1[index(i,j)]->GetMinimum();
						holder = h2[index(i,j)]->GetMinimum();
						if(cmin > holder)  cmin = holder ;
						float grid = (cmax-cmin)/20;
						h1[index(i,j)]->SetAxisRange(cmin-2.*grid, cmax+4*grid, "Y");
						h1[index(i,j)]->SetAxisRange(x1, x2, "X");
				}
		}
		auto df = new doublePanelFig(name+"_df", "", nrow, ncol, 0.4);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h1[i+nPt*j]->SetLineColor(kAzure+2);
						h1[i+nPt*j]->SetMarkerColor(kAzure+2);
						h2[i+nPt*j]->SetLineColor(kRed+1);
						h2[i+nPt*j]->SetMarkerColor(kRed+1);
						h1[i+nPt*j]->SetMarkerStyle(20);
						h2[i+nPt*j]->SetMarkerStyle(20);
						h1[i+nPt*j]->SetMarkerSize(0.5);
						h2[i+nPt*j]->SetMarkerSize(0.5);
						df->addHist(h1[i+nPt*j], i/3+1,i%3+1 );
						df->addHist(h2[i+nPt*j], i/3+1,i%3+1 );

						ratio[i+nPt*j]->GetXaxis()->SetNdivisions(505);
						ratio[i+nPt*j]->GetYaxis()->SetNdivisions(505);
						ratio[i+nPt*j]->SetMarkerStyle(20);
						ratio[i+nPt*j]->SetMarkerSize(0.5);
						ratio[i+nPt*j]->SetLineColor(kBlue+3);
						ratio[i+nPt*j]->SetMarkerColor(kBlue+3);
						ratio[i+nPt*j]->SetAxisRange(y1, y2, "Y");
						ratio[i+nPt*j]->SetAxisRange(x1, x2, "X");
						df->addHist(ratio[i+nPt*j], i/3+1,i%3+1 , 1);
						df->CD(i/3+1,i%3+1, 1);
						tl->DrawLine(x1, 1, x2, 1);
						df->CD(i/3+1,i%3+1, 0);
						tx->DrawLatexNDC(0.2, 0.87, track_label[i]);
				}
		}
		leg->AddEntry(h1[0], l1 );
		leg->AddEntry(h2[0], l2 );
		df->CD(1,1 , 0);
		leg->Draw();
		df->SaveAs(figDumpPath+name);
}

void showClosure_syst(TString name, bool isHI , float x1, float x2, float y1, float y2, TString l1, TString l2, TH1D** h1, TH1D** h1err, TH1D** h2, TH1D** h2err, TString opt="binomialRatio"){
		float cmax;
		float cmin;
		int ncol = 3, nrow = 2;
		if( isHI ) {ncol= 6; nrow = 2;}
		auto tl = new TLine(); tl->SetLineStyle(2);
		auto leg = new TLegend(.55, .57, .95, .77); leg->SetLineColor(0);
		auto ratio = binary_operation<TH1D>(name+"_ratio", h1, h2, opt);
		auto ratio_err = binary_operation<TH1D>(name+"_ratio_err", h1err, h2err, "ratio");
		auto tx = new TLatex();  tx->SetTextSize(.07);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cmax=h1[index(i,j)]->GetMaximum();
						float holder = h2[index(i,j)]->GetMaximum();
						if(cmax< holder)  cmax = holder ;
						cmin=h1[index(i,j)]->GetMinimum();
						holder = h2[index(i,j)]->GetMinimum();
						if(cmin > holder)  cmin = holder ;
						float grid = (cmax-cmin)/20;
						h1err[index(i,j)]->SetAxisRange(cmin-2.*grid, cmax+4*grid, "Y");
						h1err[index(i,j)]->SetAxisRange(x1, x2, "X");
						h1err[index(i,j)]->SetFillStyle(3345);
						h2err[index(i,j)]->SetFillStyle(3345);
						h1err[index(i,j)]->SetFillColor(kAzure+6);
						h2err[index(i,j)]->SetFillColor(kRed-7);
						ratio_err[index(i,j)]->SetFillStyle(3345);
						ratio_err[index(i,j)]->SetFillColor(kGray+1);
				}
		}
		auto df = new doublePanelFig(name+"_df", "", nrow, ncol, 0.4);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h1[i+nPt*j]->SetLineColor(kAzure+2);
						h1[i+nPt*j]->SetMarkerColor(kAzure+2);
						h2[i+nPt*j]->SetLineColor(kRed+1);
						h2[i+nPt*j]->SetMarkerColor(kRed+1);
						h1[i+nPt*j]->SetMarkerStyle(20);
						h2[i+nPt*j]->SetMarkerStyle(20);
						h1[i+nPt*j]->SetMarkerSize(0.5);
						h2[i+nPt*j]->SetMarkerSize(0.5);

						ratio_err[i+nPt*j]->GetXaxis()->SetNdivisions(505);
						ratio_err[i+nPt*j]->GetYaxis()->SetNdivisions(505);
						ratio[i+nPt*j]->SetMarkerStyle(20);
						ratio[i+nPt*j]->SetMarkerSize(0.5);
						ratio[i+nPt*j]->SetLineColor(kBlue+3);
						ratio[i+nPt*j]->SetMarkerColor(kBlue+3);
						ratio_err[i+nPt*j]->SetAxisRange(y1, y2, "Y");
						ratio_err[i+nPt*j]->SetAxisRange(x1, x2, "X");
						if(!isHI) {
								df->addHist(h1err[i+nPt*j], i/3+1,i%3+1, 0, "e2");
								df->addHist(h2err[i+nPt*j], i/3+1,i%3+1, 0, "e2");
								df->addHist(h1[i+nPt*j], i/3+1,i%3+1 );
								df->addHist(h2[i+nPt*j], i/3+1,i%3+1 );
								df->addHist(ratio_err[i+nPt*j], i/3+1,i%3+1 , 1, "e2");
								df->addHist(ratio[i+nPt*j], i/3+1,i%3+1 , 1);
								df->CD(i/3+1,i%3+1, 1);
								tl->DrawLine(x1, 1, x2, 1);
								df->CD(i/3+1,i%3+1, 0);
								tx->DrawLatexNDC(0.2, 0.87, track_label[i]);
						} else {
								df->addHist(h1err[i+nPt*j], j+1,i+1, 0, "e2");
								df->addHist(h2err[i+nPt*j], j+1,i+1, 0, "e2");
								df->addHist(h1[i+nPt*j], j+1,i+1 );
								df->addHist(h2[i+nPt*j], j+1,i+1 );
								df->addHist(ratio_err[i+nPt*j], j+1,i+1 , 1, "e2");
								df->addHist(ratio[i+nPt*j], j+1,i+1 , 1);
								df->CD(j+1,i+1, 1);
								tl->DrawLine(x1, 1, x2, 1);
								df->CD(j+1,i+1, 0);
								tx->DrawLatexNDC(0.2, 0.87, track_label[i]+"; "+cent_label[j]);
						}
				}
		}
		leg->AddEntry(h1[0], l1 );
		leg->AddEntry(h2[0], l2 );
		df->CD(1,1 , 0);
		leg->Draw();
		df->SaveAs(figDumpPath+name);
}

