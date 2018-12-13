
//#include "../../lib/import_pf_config.h"
#include "JTCSignalProducer.h"
//#include "../../lib/JTCSkimer.h"
#ifndef xPlotStyle_H
#include "xPlotStyle.h"
#endif
#ifndef jtc_macro_utility_H
#define jtc_macro_utility_H
#include "THStack.h"

int nPt = 6;
int nCent=2;
int nrow1= 2, ncol1= 3; // for pp
int nrow2= 4, ncol2 = 6; // for pb

TString dataDumpPath = "/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/myCorrelation/trunk/";
TString figDumpPath  = "/Users/tabris/cmsProjects/iJTC/macros/bJTC/fig_newWay/";
//TString *track_label;
TString *cent_label;
TString track_label[] = {"1 < p_{T}^{track} < 2 GeV",
		"2 < p_{T}^{track} < 3 GeV", "3 < p_{T}^{track} < 4 GeV",
		"4 < p_{T}^{track} < 8 GeV", "8 < p_{T}^{track} < 12 GeV", 
		"p_{T}^{track} > 12 GeV"};
//TString cent_label[]={"Cent. 0-30%", "Cent. 30-100%"};
Color_t color_vec[6] = {kBlue+1, kRed+1, kGreen+2, kAzure+7, kMagenta+2, kBlack};


//TLatex* tex = new TLatex(); 
//tex->SetTextSize(.08);
using namespace std;

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

histCase readHistCase(TString name, TString cap, TFile* f, bool sigOnly){
		histCase hc;
		initHistCase(hc);
		for(int j=0; j<nCent; ++j){
				hc.jet_corrpt[j] =(TH1D*)f->Get(cap+Form("_pt_%d", j));
				hc.jet_eta[j] =(TH1D*)f->Get(cap+Form("_eta_%d", j));
				hc.jet_phi[j] =(TH1D*)f->Get(cap+Form("_phi_%d", j));
		}
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<name+Form("_%d_%d", i, j)<<endl;
						hc.sig[i+j*nPt] =(TH2D*)f->Get(name+Form("_%d_%d", i, j));
						hc.sig_pTweighted[i+j*nPt] =(TH2D*)f->Get(name+Form("_pTweighted_%d_%d", i, j));
						hc.sig[i+j*nPt]->Scale(1.0/hc.jet_corrpt[j]->Integral());
						hc.sig_pTweighted[i+j*nPt]->Scale(1.0/hc.jet_corrpt[j]->Integral());
						if(sigOnly) continue;
						hc.sig_raw[i+j*nPt] =(TH2D*)f->Get(name+Form("_noCorr_%d_%d", i, j));
						hc.sig_pTweighted_raw[i+j*nPt] =(TH2D*)f->Get(name+Form("_pTweighted_noCorr_%d_%d", i, j));
						hc.mixing_raw[i+j*nPt] =(TH2D*)f->Get(name+Form("_mixing_noCorr_%d_%d", i, j));
						hc.mixing[i+j*nPt] =(TH2D*)f->Get(name+Form("_mixing_%d_%d", i, j));

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
T** read_flatten_duplicate(TFile* f, TString hcap, TString endcap = ""){
		T** h = new T*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h[i+nPt*j] = (T*) f->Get(hcap+Form("_%d", i)+endcap)->Clone(hcap+Form("_dup_%d_%d", i, j));
				}
		}
		return h;
}

template<typename T>
T** read_flatten_nocent(TFile* f, TString hcap){
		T** h = new T*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h[i+nPt*j] = (T*) f->Get(hcap+Form("_%d", i))->Clone(hcap+Form("_dup_%d_%d", i, j));
				}
		}
		return h;
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
T** binary_operation(TString name, T** h1, T** h2, TString opt, int nn = nPt*nCent){
		T** h = new T*[nn]; 
		for(int i=0; i<nn; ++i){
				h[i]=(T*) h1[i]->Clone(name + Form("_%d",i));
				if(opt == "ratio") 
						h[i]->Divide(h2[i]);
				else if( opt== "add") 
						h[i]->Add(h2[i]);
				else if( opt== "multiply") 
						h[i]->Multiply(h2[i]);
				else if( opt== "diff") {
						h[i]->Add(h2[i], -1);
				}
				else if( opt== "binomialRatio") {
						h[i]->Divide(h[i], h2[i], 1, 1, "B");
				}
				else cout<<"no defined operation: "<<opt<<endl;
		}
		if(nn == nPt*nCent){
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								h[i+j*nPt]->SetName(name + Form("_%d_%d",i,j));
						}
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
						sp.sig=(TH2D*)sig[i+nPt*j];
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


TH2D** getBkg(TString name, TH2D** h2, float min=1.5, float max = 2.5){
		JTCSignalProducer sp;
		TH2D** sig = new TH2D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						sig[i+nPt*j] = sp.getV2Bkg(h2[i+nPt*j], min, max);
						sp.sig=0;
				}
		}
		return sig;
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
TH1D** projY(TString name, TH2D** sig){
		JTCSignalProducer sp;
		TH1D** projy = new TH1D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<i<<", "<<j<<endl;
						projy[i+nPt*j] =(TH1D*) sp.projY(1, sig[i+nPt*j], -1, 1);
				}
		}
		return projy;
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
		int ncol = ncol1, nrow = nrow1;
		if( isHI ) {ncol= ncol2; nrow = nrow2;}
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol);
		auto tx = new TLatex();  tx->SetTextSize(.06);
		auto tl = new TLine(); tl->SetLineStyle(2);
		int ncent = isHI ? ncol2: 1;
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

		h = va_arg(ap, TH1D**);
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<nCent; ++j){
						max[i+j*nPt] = h[index(i,j)]->GetMaximum();
						min[i+j*nPt] = h[index(i,j)]->GetMinimum();
				}
		}
		for(int k=1; k<n_args; ++k){
				h = va_arg(ap, TH1D**);
				for(int i=0; i<nPt ; ++i){
						for(int j=0; j<nCent; ++j){
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
//						cout<<i<<", "<<j<<endl;
								float grid = (max[i+j*nPt]-min[i+j*nPt])/16;

								//cout<<index(i,j)<<endl; h[index(i,j)]->SetTitle("");
								h[index(i,j)]->SetLineColor(color_vec[k]);
								h[index(i,j)]->SetMarkerStyle(20);
								h[index(i,j)]->SetMarkerSize(0.8);
								if(k==3) h[index(i,j)]->SetMarkerSize(0.55);
								h[index(i,j)]->SetMarkerColor(color_vec[k]);
								h[index(i,j)]->GetXaxis()->SetNdivisions(505);
								h[index(i,j)]->SetAxisRange(x1, x2,"X");
								if(y1<y2) h[index(i,j)]->SetAxisRange(y1, y2,"Y");
								else h[index(i,j)]->SetAxisRange(min[i+j*nPt]-1.5*grid, max[i+j*nPt]+1.5*grid,"Y");

//for logy------------
//								if( k==0) h[index(i,j)]->SetMarkerSize(1.2);
//								if( index(i,j) < 4) {
//									h[index(i,j)]->SetAxisRange(.0001, .008,"Y");
//								}else {
//									h[index(i,j)]->SetAxisRange(0.00001, 0.03,"Y");
//								}
//--------------------
								if(!isHI) { 
										cm->drawHist(h[index(i,j)], i+1);
//										gPad->SetLogy();
										//cout<<i+1<<endl;
										tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
										cm->cd(i+1);
										tl->DrawLine(x1, line, x2, line);
								} else { 
										cm->drawHist(h[index(i,j)], i+1, 4-j);
		//								cout<<"i = "<<i<<", j = "<<j<<": "<<h[index(i,j)]->GetName()<<endl;
										cm->CD(i+1, 4-j); //gPad->SetLogy();
										//gPad->SetLogy();
										tmp = track_label[i]+"; "+cent_label[j];
										tx->DrawLatexNDC(0.1, 0.93, tmp);
										if(min[i+j*nPt]-grid< line && max[i+j*nPt]+grid> line) tl->DrawLine(x1, line, x2, line);
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

template <typename T>
T** copy_duplicate(TString name, T** h){
        T** hh = new T*[nPt*nCent];
        for(int i=0; i<nPt ; ++i){
                for(int j=0; j<nCent ; ++j){
                        hh[i+nPt*j]=(T*)h[i]->Clone(name+Form("_%d_%d",i,j));
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
								//cout<<index(i,j)<<endl; h[index(i,j)]->SetTitle("");
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
								//cout<<"i = "<<i<<", j = "<<j<<": "<<h[index(i,j)]->GetName()<<endl;
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

void showClosure(TString name, bool isHI , float x1, float x2, float y1, float y2, TString l1, TString l2, TH1D** h1, TH1D** h2, TString opt="binomialRatio"){
		float cmax;
		float cmin;
		int ncol = 3, nrow = 2;
		if( isHI ) {ncol= 6; nrow = 2;}
		int ncent = isHI ? nCent: 1;
		auto tl = new TLine(); tl->SetLineStyle(2);
		auto leg = new TLegend(.55, .57, .95, .77); leg->SetLineColor(0);
		auto ratio = binary_operation<TH1D>(name+"_ratio", h1, h2, opt);
		auto tx = new TLatex();  tx->SetTextSize(.07);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
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
				for(int j=0; j<ncent; ++j){
						h1[i+nPt*j]->SetLineColor  (kGreen);
						h1[i+nPt*j]->SetMarkerColor(kGreen);
						//h1[i+nPt*j]->SetLineColor  (kAzure+2);
						//h1[i+nPt*j]->SetMarkerColor(kAzure+2);
						h2[i+nPt*j]->SetLineColor  (kRed+1);
						h2[i+nPt*j]->SetMarkerColor(kRed+1);
						//h2[i+nPt*j]->SetLineColor  (kAzure+7);
						//h2[i+nPt*j]->SetMarkerColor(kAzure+7);
						//h1[i+nPt*j]->SetLineColor  (kMagenta+2);
						//h1[i+nPt*j]->SetMarkerColor(kMagenta+2);

						h1[i+nPt*j]->SetMarkerStyle(20);
						h2[i+nPt*j]->SetMarkerStyle(20);
						h1[i+nPt*j]->SetMarkerSize(0.8);
						h2[i+nPt*j]->SetMarkerSize(0.8);
						df->addHist(h1[i+nPt*j], i/3+1,i%3+1 );
						df->addHist(h2[i+nPt*j], i/3+1,i%3+1 );

						ratio[i+nPt*j]->GetXaxis()->SetNdivisions(505);
						ratio[i+nPt*j]->GetYaxis()->SetNdivisions(505);
						ratio[i+nPt*j]->SetMarkerStyle(20);
						ratio[i+nPt*j]->SetMarkerSize(0.8);
						ratio[i+nPt*j]->SetLineColor(kBlue+3);
						ratio[i+nPt*j]->SetMarkerColor(kBlue+3);
						ratio[i+nPt*j]->SetAxisRange(y1, y2, "Y");
						ratio[i+nPt*j]->SetAxisRange(x1, x2, "X");
						if(isHI){
								df->addHist(h1[i+nPt*j], j+1,i+1 );
								df->addHist(h2[i+nPt*j], j+1,i+1 );
								df->addHist(ratio[i+nPt*j], j+1,i+1 , 1);
								df->CD(j+1,i+1, 1);
								tl->DrawLine(x1, 1, x2, 1);
								df->CD(j+1,i+1, 0);
								tx->DrawLatexNDC(0.2, 0.87, track_label[i]+"; "+cent_label[j]);
						} else {
								df->addHist(h1[i+nPt*j], i/3+1,i%3+1 );
								df->addHist(h2[i+nPt*j], i/3+1,i%3+1 );
								df->addHist(ratio[i+nPt*j], i/3+1,i%3+1 , 1);
								df->CD(i/3+1,i%3+1, 1);
								tl->DrawLine(x1, 1, x2, 1);
								df->CD(i/3+1,i%3+1, 0);
								tx->DrawLatexNDC(0.2, 0.87, track_label[i]);
						}
				}
		}
		leg->AddEntry(h1[0], l1 );
		leg->AddEntry(h2[0], l2 );
		df->CD(1,1 , 0);
		leg->Draw();
		df->SaveAs(figDumpPath+name);
}


void show_plot_syst(TString name, bool isHI ,float l,  float x1, float x2, float y1, float y2, TH1D** ratio, TH1D** ratio_err){
		float cmax;
		float cmin;
		int ncent = isHI ? 2 : 1;
		int ncol = 3, nrow = 2;
		if( isHI ) {ncol= 6; nrow = 2;}
		auto tl = new TLine(); tl->SetLineStyle(2);
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol);
		auto tx = new TLatex();  tx->SetTextSize(.07);

		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
						ratio_err[index(i,j)]->SetFillStyle(1001);
						ratio_err[index(i,j)]->SetFillColorAlpha(kGray+2, 0.6);
						ratio[i+nPt*j]->SetLineColor(kAzure+2);
						ratio[i+nPt*j]->SetMarkerColor(kAzure+2);
						ratio[i+nPt*j]->SetLineColor(kRed+1);
						ratio[i+nPt*j]->SetMarkerColor(kRed+1);
						ratio[i+nPt*j]->SetMarkerStyle(20);
						ratio[i+nPt*j]->SetMarkerStyle(20);
						ratio[i+nPt*j]->SetMarkerSize(0.5);
						ratio[i+nPt*j]->SetMarkerSize(0.5);

						ratio_err[i+nPt*j]->GetXaxis()->SetNdivisions(505);
						ratio_err[i+nPt*j]->GetYaxis()->SetNdivisions(505);
						ratio[i+nPt*j]->SetMarkerStyle(20);
						ratio[i+nPt*j]->SetMarkerSize(0.5);
						ratio[i+nPt*j]->SetLineColor(kBlue+3);
						ratio[i+nPt*j]->SetMarkerColor(kBlue+3);
						ratio_err[i+nPt*j]->SetAxisRange(y1, y2, "Y");
						ratio_err[i+nPt*j]->SetAxisRange(x1, x2, "X");
						ratio_err[i+nPt*j]->SetTitle("");
						if(!isHI) {
								cm->drawHist(ratio_err[i+nPt*j], i/3+1,i%3+1 , "e2");
								cm->drawHist(ratio[i+nPt*j], i/3+1,i%3+1 );
								cm->CD(i/3+1,i%3+1);
								tl->DrawLine(x1, l, x2, l);
								tx->DrawLatexNDC(0.2, 0.92, track_label[i]);
						} else {
								cm->drawHist(ratio_err[i+nPt*j], j+1,i+1 ,  "e2");
								cm->drawHist(ratio[i+nPt*j], j+1,i+1 );
								cm->CD(j+1,i+1);
								tl->DrawLine(x1, l, x2, l);
								tx->DrawLatexNDC(0.2, 0.92, track_label[i]+"; "+cent_label[j]);
						}
				}
		}
		cm->SaveAs(figDumpPath+name);
}

TCanvas* showClosure_syst(TString name, bool isHI , float x1, float x2, float y1, float y2, TString l1, TString l2, TH1D** h1, TH1D** h1err, TH1D** h2, TH1D** h2err, TString opt="binomialRatio", TString opt2= "ratio"){
		float cmax;
		float cmin;
		int ncol = 3, nrow = 2;
		if( isHI ) {ncol= 6; nrow = 2;}
		auto tl = new TLine(); tl->SetLineStyle(2);
		auto leg = new TLegend(.45, .5, .95, .77); leg->SetLineColor(0);
		auto ratio = binary_operation<TH1D>(name+"_ratio", h1, h2, opt);
		auto ratio_err = binary_operation<TH1D>(name+"_ratio_err", h1err, h2err, opt2);
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
		return (TCanvas*) df;
}

TH1D* invariant_rebin(TH1D* h, const int nbin, Double_t *newbin){
		//	Double_t newdrbin[n_new_drbin+1] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.45,0.6,0.8,1. };
		for(int i=1; i<h->GetNbinsX()+1; ++i){
				h->SetBinContent(i, h->GetBinContent(i)*h->GetBinWidth(i));
				h->SetBinError(i, h->GetBinError(i)*h->GetBinWidth(i));
		}
		TString name = h->GetName();
		name+="_rebin";
		TH1D* htm = (TH1D*)h->Rebin(nbin, "", newbin);
		//		delete h;
		for(int i=1; i<htm->GetNbinsX()+1; ++i){
				//				cout<<h->GetBinWidth(i)<<endl;
				htm->SetBinContent(i, htm->GetBinContent(i)/htm->GetBinWidth(i));
				htm->SetBinError(i, htm->GetBinError(i)/htm->GetBinWidth(i));
		}
		return htm;
}

TH1D** rebin_all(TString name, TH1D** h, const int nbin, Double_t *newbin){
		TH1D** dr = new TH1D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						dr[i+nPt*j]=invariant_rebin( h[i+nPt*j], nbin, newbin);
				}
		}
		return dr;
}

void normalization(TH2D** h){
		float xw = h[0]->GetXaxis()->GetBinWidth(1);
		float yw = h[0]->GetYaxis()->GetBinWidth(1);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h[i+nPt*j]->Scale(1.0/xw/yw);
				}
		}
}

void normal2unit(TH1D** h){
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						float cc =  h[i+nPt*j]->Integral();
						h[i+nPt*j]->Scale(1.0/cc);
				}
		}
}

void smooth(TH1D** h)
{
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						h[i+nPt*j]->Smooth();
				}
		}
}

TH1D** projectionX(TString name, TH2D** h){
		TH1D** projx = new TH1D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						projx[i+nPt*j] =(TH1D*) h[i+nPt*j]->ProjectionX();
				}
		}
		return projx;

}

void set_bin0(TH1D** dr, TH1D** syst, float cc = 0){
		// if the bin content is less than the error, set the bin as 0 with error 0;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){	
						for(int k=1; k<dr[i+j*nPt]->GetNbinsX()+1; k++){
								//								cout<<fabs(dr[i+j*nPt]->GetBinContent(k))<<" err = "<<syst[i+j*nPt]->GetBinError(k)<<endl;
								if(fabs(dr[i+j*nPt]->GetBinContent(k)-cc)>pow(pow(syst[i+j*nPt]->GetBinError(k),2)+pow(dr[i+j*nPt]->GetBinError(k),2), 0.5))continue;
								//								if(dr[i+j*nPt]->GetBinContent(k)>0.06)continue;
								dr[i+j*nPt]->SetBinContent(k, 0);
								dr[i+j*nPt]->SetBinError  (k, 0);
						}
				}
		}
}

void set_bin_with_err(TH1D** dr, TH1D** syst, float x0, float cc = 0, float err = 0.02){
		// if the bin content larger than x0, is less than the error, set the bin as 0 with error=err, do both for dr and syst;
		// to make sure the bin is visible, set the error not equal 0;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){	
						for(int k=1; k<dr[i+j*nPt]->GetNbinsX()+1; k++){
								if(dr[i+j*nPt]->GetBinCenter(k)<x0 || fabs(dr[i+j*nPt]->GetBinContent(k)-cc)>pow(pow(syst[i+j*nPt]->GetBinError(k),2)+pow(dr[i+j*nPt]->GetBinError(k),2), 0.5))continue;
								dr[i+j*nPt]->SetBinContent(k, cc);
								dr[i+j*nPt]->SetBinError  (k, .3*err);
								syst[i+j*nPt]->SetBinContent(k, cc);
								syst[i+j*nPt]->SetBinError  (k, err);
						}
				}
		}
}

void swip_hist(TH1D** dr, float xx, float cc, float err = 0){
		// set the hist bin content at x< xx to be the cc with error 0;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){	
						for(int k=1; k<dr[i+j*nPt]->GetNbinsX()+1; k++){
								if(dr[i+j*nPt]->GetBinCenter(k)<= xx)continue;
								dr[i+j*nPt]->SetBinContent(k, cc);
								dr[i+j*nPt]->SetBinError  (k, err);
						}
				}
		}
}

TH1D**  sum_pt(TString name, TH1D** h){
		TH1D** hh = new TH1D*[nCent];
		for(int j=0; j<nCent; ++j){	
				hh[j]= (TH1D*)h[j*nPt]->Clone(name+Form("_sum_%d",j));
				for(int i=1; i<nPt; ++i){
						hh[j]->Add(h[i+j*nPt]);
				}
		}
		return hh;
}

TCanvas* showPanel_syst(TString name, bool isHI , float x1, float x2, float y1, float y2, TString l1, TString l2, TH1D** h1, TH1D** h1err, TH1D** h2, TH1D** h2err, TString opt="binomialRatio", TString opt2= "ratio"){
		float cmax;
		float cmin;
		int ncol = 1, nrow = 1;
		if( isHI ) {ncol= 2; nrow = 1;}
		auto tl = new TLine(); tl->SetLineStyle(2);
		auto leg = new TLegend(.45, .5, .95, .77); leg->SetLineColor(0);
		auto ratio = binary_operation<TH1D>(name+"_ratio", h1, h2, opt, 2);
		auto ratio_err = binary_operation<TH1D>(name+"_ratio_err", h1err, h2err, opt2, 2);
		auto tx = new TLatex();  tx->SetTextSize(.07);
		int j=0;
		cout<<"here"<<endl;
		for(int i=0; i<ncol; ++i){
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
		auto df = new doublePanelFig(name+"_df", "", nrow, ncol, 0.4);
		for(int i=0; i<ncol; ++i){
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
						tx->DrawLatexNDC(0.2, 0.87, "p_{T}^{track} > 1GeV");
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
		leg->AddEntry(h1[0], l1 );
		leg->AddEntry(h2[0], l2 );
		df->CD(1,1 , 0);
		leg->Draw();
		df->SaveAs(figDumpPath+name);
		return (TCanvas*) df;

}

TH1D* integral_yield(TString name, TH1D** hh ){
		float ptbin[7] = {1, 2, 3, 4, 8, 12, 20};
		TString labels[7] = {"1", "2", "3", "4", "8", "12", "18"};
		TH1D* h = new TH1D(name, "", 6, ptbin);
		h->Sumw2();
		for(int i=0; i<6; i++){
				float yield = hh[i]->Integral();
				h->SetBinContent(i+1, yield/(ptbin[i+1]-ptbin[i]));
				h->SetBinError(i+1, 1/yield/(ptbin[i+1]-ptbin[i]));
				h->GetXaxis()->SetBinLabel(i+1, labels[i]);
		}
		h->GetXaxis()->SetTitle("p_{T}^{trk}");
		return h;
}

#endif
