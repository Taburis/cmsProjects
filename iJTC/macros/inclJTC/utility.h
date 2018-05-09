
#include "../../lib/JTCSignalProducer.h" 
#include "../../lib/stackHist.h" 
#ifndef xPlotStyle_H
#include "../../lib/xPlotStyle.h" 
#endif  
#include "read.h" 
#include "stdarg.h"

#define nPt 9

TString tmp;

TString track_label[] = {"0.7 < p_{T}^{track} < 1 GeV","1 < p_{T}^{track} < 2 GeV",
		"2 < p_{T}^{track} < 3 GeV", "3 < p_{T}^{track} < 4 GeV",
		"4 < p_{T}^{track} < 8 GeV", "8 < p_{T}^{track} < 12 GeV", 
		"12 < p_{T}^{track} < 16 GeV", "16 < p_{T}^{track} < 20 GeV", "p_{T}^{track} > 20 GeV"};
TString cent_label[] = {"Cent. 0-10%", "Cent. 10-30%", "Cent. 30-50%","Cent. 50-100%"};

Color_t color_vec[6] = {kBlue, kRed+1, kGreen+2, kAzure+7, kMagenta+2, kBlack};
//		"12 < p_{T}^{track} < 16 GeV", "p_{T}^{track} > 16 Gev"};
void clean_input(){
		for(int i=0; i<9; ++i){
				for(int j=0; j<5; ++j){
						delete hraw_sig[i][j];
						delete hraw_sig_pTweighted[i][j];
						delete hmixing[i][j];
				}
		}
		for(int i=0; i<8; ++i){
				for(int j=0; j<1; ++j){
						delete raw_sig[i][j];
						delete raw_sig_pTweighted[i][j];
						delete mixing[i][j];
				}
		}
}

void batch_operation(TString name,TH1D* h[9][4], TH1D* h1[9][4], TH1D* h2[9][4],  TString opt, bool ispp){
		cout<<"oprating on "<<name<<"..."<<endl;
		int ncent = ispp ? 1 : 4;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
						h[i][j]=(TH1D*) h1[i][j]->Clone(name + Form("_%d_%d",i,j));
						if(opt == "ratio") 
								h[i][j]->Divide(h2[i][j]);
						else if( opt== "add") 
								h[i][j]->Add(h2[i][j]);
						else if( opt== "multiply") 
								h[i][j]->Multiply(h2[i][j]);
						else if( opt== "diff") {
								h[i][j]->Add(h2[i][j], -1);
						}
						else if( opt== "binomialRatio") {
								h[i][j]->Divide(h[i][j], h2[i][j], 1, 1, "B");
						}
						else cout<<"no defined operation: "<<opt<<endl;
				}
		}
}

void batch_operation2(TString name,TH1D* h[9][4], TH1D* h1[9][4], TH1D* h2[9][4],  TString opt){
		cout<<"oprating on "<<name<<"..."<<endl;
		int ncent = 4;//ispp ? 1 : 4;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
						h[i][j]=(TH1D*) h1[i][j]->Clone(name + Form("_%d_%d",i,j));
						if(opt == "ratio") 
								h[i][j]->Divide(h2[i][0]);
						else if( opt== "add") 
								h[i][j]->Add(h2[i][0]);
						else if( opt== "multiply") 
								h[i][j]->Multiply(h2[i][0]);
						else if( opt== "diff") {
								h[i][j]->Add(h2[i][0], -1);
						}
						else if( opt== "binomialRatio") {
								h[i][j]->Divide(h[i][j], h2[i][0], 1, 1, "B");
						}
						else cout<<"no defined operation: "<<opt<<endl;
				}
		}
}

void batch_operation2D(TString name,TH2D* h[9][4], TH2D* h1[9][4], TH2D* h2[9][4],  TString opt, bool ispp){
		cout<<"oprating on "<<name<<"..."<<endl;
		int ncent = ispp ? 1 : 4;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
						h[i][j]=(TH2D*) h1[i][j]->Clone(name + Form("_%d_%d",i,j));
						if(opt == "ratio") 
								h[i][j]->Divide(h2[i][j]);
						else if( opt== "add") 
								h[i][j]->Add(h2[i][j]);
						else if( opt== "multiply") 
								h[i][j]->Multiply(h2[i][j]);
						else if( opt== "diff") {
								h[i][j]->Add(h2[i][j], -1);
						}
						else if( opt== "binomialRatio") {
								h[i][j]->Divide(h[i][j], h2[i][j], 1, 1, "B");
						}
						else cout<<"no defined operation: "<<opt<<endl;
				}
		}
}

void batch_dump(TString name, TString fname, bool ispp, TH1D* h1[9][4]){
		int ncent = ispp ? 1 : 4;
		auto wf = new TFile(dataDumpPath+fname+".root", "recreate");
		wf->cd();
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
						h1[i][j]->SetName(name+Form("_%d_%d", i, j));
						h1[i][j]->Write();
				}
		}
		wf->Close();
}

void batch_read1D(TString name, TString fname,  bool ispp, TH1D* h[9][4]){
	int ncent = ispp ? 1 : 4;
		TFile *f = new TFile(dataDumpPath +fname+".root");
		for(int i=0; i<9; ++i){
				for(int j=0; j<ncent; ++j){
						h[i][j]=(TH1D*) f->Get(name+Form("_%d_%d",i,j));
						//				cout<<h[i][j]->GetName()<<endl;
				}
		}
}

void batch_ratio(TString name, TH1D* h1[9][4], TH1D* h2[9][4], TH1D* h[9][4], bool ispp){
		int ncent = ispp ? 1 : 4;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
						h[i][j]=(TH1D*) h1[i][j]->Clone(name + Form("_%d_%d",i,j));
						h[i][j]->Divide(h2[i][j]);
						h[i][j]->SetAxisRange(0,2, "Y");
				}
		}
}

void batchGet2D(TString name,TString name2, TH2D* h[9][4], bool isNumber=0, bool ispp = 1){
		int ncent = ispp ? 1 : 4;
		TFile *f = new TFile(dataDumpPath +name2+"_JTCSignal.root");
		for(int i=0; i<9; ++i){
				for(int j=0; j<ncent; ++j){
	//					cout<<name+Form("_%d_%d",i,j)<<endl;
						if(isNumber) h[i][j]=(TH2D*) f->Get(name+Form("_%d_%d",i,j));
						else h[i][j]=(TH2D*) f->Get(name+Form("_pTweighted_%d_%d",i,j));
	//									cout<<h[i][j]->GetName()<<endl;
				}
		}
}

TH1D** flattenTH1D_94(TH1D* h[9][4] ){
		TH1D** hh;
		hh=new TH1D*[36];
		for(int i=0; i<4; ++i){
				for(int j=0; j<9; ++j){
						hh[j+i*9] = h[j][i];
				}
		}
		return hh;
}

int index_94(int i, int j){
		return i+j*9;
}

void showPlot(TString name, bool ispp , float x1, float x2, int n_args, ...){
		va_list ap;
		va_start(ap, n_args);
		int ncol = 3, nrow = 3;
		if( !ispp ) {ncol= 4; nrow = 9;}
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		int ncent = ispp ? 1: 4;
		TH1D** h;
		for(int k=0; k<n_args; ++k){
				h = va_arg(ap, TH1D**);
				for(int i=0; i<9 ; ++i){
						for(int j=0; j<ncent; ++j){
						//		cout<<index_94(i,j)<<endl;
								h[index_94(i,j)]->SetTitle("");
								h[index_94(i,j)]->SetLineColor(color_vec[k]);
								h[index_94(i,j)]->GetXaxis()->SetNdivisions(505);
								h[index_94(i,j)]->SetAxisRange(x1, x2,"X");
								if(ispp ) { cm->drawHist(h[index_94(i,j)], i+1);
										cm->cd(i+1);
										tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
								}
								else { 
										cm->drawHist(h[index_94(i,j)], i+1, 4-j);
										cout<<"i = "<<i<<", j = "<<j<<": "<<h[index_94(i,j)]->GetName()<<endl;
										cm->CD(i+1, 4-j);
										gPad->SetLogy();
										tmp = track_label[i]+" "+cent_label[j];
										tx->DrawLatexNDC(0.1, 0.93, tmp);
								}
						}
				}
		}
		cm->SaveAs(FigDumpPath+name+".gif");
}

void getSig(TString name, TH2D* h[9][4], bool isNumber=0, bool ispp = 1){
		batchGet2D("signal_"+name,name, h, isNumber, ispp);
}

void getRawSig(TString name, TH2D* h[9][4], bool isNumber=0, bool ispp= 1){
		batchGet2D("sig_mix_corrected_"+name,name, h, isNumber, ispp);
}

void pullSig(TString fname, float sidemin=1.5, float sidemax=2.5, bool doSeagull=0, bool ispp =1){
		cout<<"pulling singal for "<<fname;
		TString trkbin [] = {"0.7", "1", "2", "3", "4", "8", "12", "16", "20", "999"};
		TString centbin [] = {"0", "30", "100"};
		JTCSignalProducer *sp1[9][4];
		JTCSignalProducer *sp2[9][4];
		TFile * wf = new TFile(dataDumpPath+fname+"_JTCSignal.root","recreate");
		wf->cd();
		cout<<".";
		int nCent= ispp ? 1 : 4;
		for(int i=0; i<9 ; ++i){
				cout<<".";
				for(int j=0; j<nCent ; ++j){
						tmp = "track pt in ["+trkbin[i]+", "+trkbin[i+1]\
							   +"), cent in ["+centbin[j]+", "+centbin[j+1]+"), jet pt in [120, 1000]";
						raw_sig[i][j]->SetTitle(tmp);
						raw_sig_pTweighted[i][j]->SetTitle(tmp);
						//cout<<raw_sig[i][j]->GetName()<<endl;
						mixing [i][j]->SetTitle(tmp);
						sp1[i][j] = new JTCSignalProducer(raw_sig[i][j], mixing[i][j]);
						sp1[i][j]->sideMin=sidemin; sp1[i][j]->sideMax=sidemax;
						sp1[i][j]->getSignal(fname+Form("_%d_%d", i, j), doSeagull);
						sp1[i][j]->WriteTH2();
						sp2[i][j] = new JTCSignalProducer(raw_sig_pTweighted[i][j], mixing[i][j]);
						sp2[i][j]->sideMin=sidemin; sp2[i][j]->sideMax=sidemax;
						sp2[i][j]->getSignal(fname+Form("_pTweighted_%d_%d", i, j), doSeagull);
						sp2[i][j]->WriteTH2();
						delete sp1[i][j];
						delete sp2[i][j];
						//								cout<<sp2[i][j]->sig->GetName()<<endl;
				}
		}
		cout<<". Done"<<endl;;
		wf->Close();
		//clearInput();
		cout<<"signal dumped to "<<dataDumpPath+fname+"_JTCSignal.root"<<endl;

		clean_input();
}

void checkBkg(TString name, int ispp = 1){
		tmp =dataDumpPath+name+"_JTCSignal.root";
		TFile *f = TFile::Open(tmp);
		JTCSignalProducer *sp1[9][4];
		JTCSignalProducer *sp2[9][4];
		int nrow = 9, ncol=4;
		int ncent=4;
		if(ispp) { nrow = 3; ncol = 3; ncent =1;
		}
		cout<<"cent "<<ncent<<endl;
		auto tx = new TLatex();  tx->SetTextSize(.08);
		auto cp1 = new mCanvasLoose("sig_deta_"+name,            "", nrow, ncol);
		auto cp2 = new mCanvasLoose("sig_deta_JS_"+name,         "", nrow, ncol);
		auto cp3 = new mCanvasLoose("bkgCheck_sideBand"+name,    "", nrow, ncol);
		auto cp4 = new mCanvasLoose("bkgCheck_sideBand_JS"+name, "", nrow, ncol);
		auto cp5 = new mCanvasLoose("bkgCheck_sliceSB_JS"+name, "", nrow, ncol);
		for(int i=0; i<9 ; ++i){
				for(int j=0; j<ncent; ++j){
						cout<<i<<", "<<j<<endl;
						tmp = track_label[i];
						sp1[i][j] = new JTCSignalProducer();
						sp2[i][j] = new JTCSignalProducer();
						sp1[i][j]->sig=(TH2D*) f->Get("signal_"+name+Form("_%d_%d", i, j));
						sp2[i][j]->sig=(TH2D*) f->Get("signal_"+name+Form("_pTweighted_%d_%d", i, j));
						//cout<<sp1[i][j]->sig->GetName()<<endl;
						sp1[i][j]->getSignal_phiSideBand(Form("sideBand_%d_%d",i,j));
						sp2[i][j]->getSignal_phiSideBand(Form("sideBand_pTweighted_%d_%d",i,j));
						sp1[i][j]->getSignal_dEta(Form("sig_deta_%d_%d",i,j));
						sp2[i][j]->getSignal_dEta(Form("sig_deta_pTweighted_%d_%d",i,j));

						if( ispp) cp5->cd(i+1);
						else cp5->CD(i+1, 4-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp2[i][j]->drawSliceSideBand(name+Form("_slice_%d_%d",i,j));
						tx->DrawLatexNDC(0.02,0.95, tmp); 

						if( ispp) cp3->cd(i+1);
						else cp3->CD(i+1, 4-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp1[i][j]->drawSideBandCheck();
						tx->DrawLatexNDC(0.02,0.93, tmp); 

						if( ispp) cp1->cd(i+1);
						else cp1->CD(i+1, 4-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp1[i][j]->histStyle(sp1[i][j]->sig_deta);
						sp1[i][j]->sig_deta->SetAxisRange(-2.5, 2.499, "X");
						sp1[i][j]->sig_deta->Draw();
						sp1[i][j]->side_deta->Draw("same");
						tx->DrawLatexNDC(0.02,0.93, tmp); 

						if( ispp) cp4->cd(i+1);
						else cp4->CD(i+1, 4-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp2[i][j]->drawSideBandCheck();
						tx->DrawLatexNDC(0.02,0.95, tmp); 

						if( ispp) cp2->cd(i+1);
						else cp2->CD(i+1, 4-j);
						gPad->SetTopMargin(0.1);
						gPad->SetBottomMargin(0.18);
						sp2[i][j]->histStyle(sp2[i][j]->sig_deta);
						sp2[i][j]->sig_deta->SetAxisRange(-2.5, 2.499, "X");
						sp2[i][j]->sig_deta->Draw();
						sp2[i][j]->side_deta->Draw("same");
						tx->DrawLatexNDC(0.02,0.95, tmp); 
				}
		}
		cp1->SaveAs(FigDumpPath+"bkgCheck_"+name+".pdf");
		cp2->SaveAs(FigDumpPath+"bkgCheck_JS_"+name+".pdf");
		cp3->SaveAs(FigDumpPath+"bkgCheck_sideBand_"+name+".pdf");
		cp4->SaveAs(FigDumpPath+"bkgCheck_sideBand_JS_"+name+".pdf");
		cp5->SaveAs(FigDumpPath+"bkgCheck_sliceSB_JS_"+name+".pdf");
		for(int i=0; i<8 ; ++i){
				for(int j=0; j<1 ; ++j){
						delete sp1[i][j];
						delete sp2[i][j];
				}
		}
}


void batch_sig_subtraction(TString name, TString name1, TString name2, bool isNumber = 0, bool ispp = 1){
		int ncent = ispp ? 1: 4;
		TH2D* sig[9][4], *sub[9][4];
		float trkbin [] = {0.3, 1, 1, 1, 4, 4, 4, 4, 1};
		getSig(name1, sig, isNumber, ispp);
		getSig(name2, sub, isNumber, ispp);
		TFile *wf = new TFile(dataDumpPath+name+Form("_JTCSignal.root"), "recreate");
		wf->cd();
		TH2D *corr[9][4];
		for(int i=0; i<9; ++i){
				for(int j=0; j<ncent; ++j){
						cout<<sig[i][j]->GetName()<<endl;
						if(isNumber) corr[i][j]=(TH2D*)sig[i][j]->Clone("signal_"+name+Form("_%d_%d",i,j));
						else corr[i][j]=(TH2D*)sig[i][j]->Clone("signal_"+name+Form("_pTweighted_%d_%d",i,j));
						corr[i][j]->Add(sub[i][j],-1);
						if(isNumber)corr[i][j]->Scale(1.0/trkbin[i]);
						corr[i][j]->Write();
				}
		}
		wf->Close();
}

void batch_subtraction(TString name, TString name1, bool ispp1, TString name2, bool ispp2){
		TH2D* sig[9][4], *sub[9][4];
		bool isNumber =0;
		float trkbin [] = {0.3, 1, 1, 1, 4, 4, 4, 4, 1};
		getSig(name1, sig, isNumber, ispp1);
		getSig(name2, sub, isNumber, ispp2);
		TFile *wf = new TFile(dataDumpPath+name+Form("_JTCSignal.root"), "recreate");
		wf->cd();
		TH2D *corr[9][4];
		for(int i=0; i<9; ++i){
				for(int j=0; j<4; ++j){
						int k1= ispp1 ? 0 : j;
						int k2= ispp2 ? 0 : j;
						//						cout<<sig[i][k1]->GetName()<<endl;
						if(isNumber) corr[i][j]=(TH2D*)sig[i][k1]->Clone("signal_"+name+Form("_%d_%d",i,j));
						else corr[i][j]=(TH2D*)sig[i][k1]->Clone("signal_"+name+Form("_pTweighted_%d_%d",i,j));
						corr[i][j]->Add(sub[i][k2],-1);
						if(isNumber)corr[i][j]->Scale(1.0/trkbin[i]);
						corr[i][j]->Write();
				}
		}
		wf->Close();

}

void batch_sig_pp_sub_pb(TString name, TString pp, TString pb, bool isNumber){
		TH2D* sig[9][4], *sub[9][4];
		float trkbin [] = {0.3, 1, 1, 1, 4, 4, 4, 4, 1};
		getSig(pp, sig, isNumber, 1);
		getSig(pb, sub, isNumber, 0);
		TH2D *corr[9][4];
		TFile *wf = new TFile(dataDumpPath+name+Form("_JTCSignal.root"), "recreate");
		wf->cd();
		for(int i=0; i<9; ++i){
				for(int j=0; j<4; ++j){
						if(isNumber) corr[i][j]=(TH2D*)sig[i][0]->Clone("signal_"+name+Form("_%d_%d",i,j));
						else corr[i][j]=(TH2D*)sig[i][0]->Clone("signal_"+name+Form("_pTweighted_%d_%d",i,j));
						corr[i][j]->Add(sub[i][j],-1);
						if(isNumber)corr[i][j]->Scale(1.0/trkbin[i]);
						corr[i][j]->Write();
				}
		}
		wf->Close();
}	

mCanvasLoose * show_1D_pp(TString name, TH1D* h[9][4], float x1=0, float x2=-1, bool ispp = 1, float y1=0, float y2 = -1){
		int ncol = 3, nrow = 3;
		if( !ispp ) {ncol= 4; nrow = 9;}
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		int ncent = ispp ? 1: 4;
		for(int i=0; i<9 ; ++i){
				for(int j=0; j<ncent; ++j){
						h[i][0]->SetTitle("");
						if(ispp ) { cm->drawHist(h[i][0], i+1);
								h[i][0]->GetXaxis()->SetNdivisions(505);
								if( x1< x2) h[i][0]->SetAxisRange(x1, x2,"X");
								if( y1< y2) h[i][0]->SetAxisRange(y1, y2,"Y");
								cm->cd(i+1);
								tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
						}
						else { cm->drawHist(h[i][j], i+1, 4-j);
								h[i][j]->GetXaxis()->SetNdivisions(505);
								if( x1< x2) h[i][j]->SetAxisRange(x1, x2,"X");
								if( y1< y2) h[i][j]->SetAxisRange(y1, y2,"Y");
								cm->CD(i+1, 4-j);
								tmp = track_label[i]+" "+cent_label[j];
								tx->DrawLatexNDC(0.1, 0.93, tmp);
						}
				}
		}
		cm->SaveAs(FigDumpPath+name+".gif");
		return cm;
}


mCanvasLoose * show_1D_overlay(TString name, TH1D* h[9][4], TH1D* h2[9][4], TString lab1, TString lab2, float x1=0, float x2=-1, bool ispp = 1){
		int ncol = 3, nrow = 3, ncent = 1;
		if( !ispp ) {ncol= 4; nrow = 9; ncent = 4;}
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		auto line = new TLine();  line ->SetLineStyle(2);
		for(int i=0; i<9 ; ++i){
				for(int j=0; j<ncent ; ++j){
						float ymin = min(h2[i][j]->GetMinimum(), h[i][j]->GetMinimum());
						float ymax = max(h2[i][j]->GetMaximum(), h[i][j]->GetMaximum());
						float range = ymax-ymin;
						h[i][j]->SetAxisRange(ymin-0.05*range, ymax+0.05*range, "Y");
						h[i][j]->SetTitle("");
						h[i][j]->GetXaxis()->SetNdivisions(505);
						h2[i][j]->SetLineColor(kRed);
						if( x1< x2) h[i][j]->SetAxisRange(x1, x2,"X");
						if(ispp) {
								cm->drawHist(h[i][j], i+1);
								cm->drawHist(h2[i][j], i+1);
								cm->cd(i+1);
								line->DrawLine(x1, 0, x2, 0);
								tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
						}else {
								cm->drawHist(h[i][j], i+1, 4-j);
								cm->drawHist(h2[i][j], i+1, 4-j);
								cm->CD(i+1, 4-j);
								line->DrawLine(x1, 0, x2, 0);
								tmp = track_label[i]+" "+cent_label[j];
								tx->DrawLatexNDC(0.1, 0.93, tmp);
						}
				}
		}
		auto tl = new TLegend(0.65, 0.6, 0.95, 0.87); tl->SetLineColor(0);
		//auto tl = new TLegend(0.6, 0.2, 0.95, 0.5); tl->SetLineColor(0);
		tl->AddEntry(h[0][0], lab1);
		tl->AddEntry(h2[0][0], lab2);
		cm->cd(1); tl->Draw();
		cm->SaveAs(FigDumpPath+name+".gif");
		return cm;
}

mCanvasLoose * show_1D_sig_deta(TString name, TString name1, bool isNumber = 0){
		TH2D* sig[9][4];
		TH1D* h[9][4];
		auto sp = new JTCSignalProducer();
		getSig(name1, sig, isNumber);
		for(int i=0; i<9 ; ++i){
				sp->sig=sig[i][0];
				h[i][0]= sp->getSignal_dEta(name1+Form("_%d_%d",i,0));
				h[i][0]->SetAxisRange(-1., 0.99,"X");
		}
		return show_1D_pp(name, h);
}

void show_2D_colz(TString name, TH2D* h[9][4], bool doLogz= 0, bool ispp = 1){
		int ncent = ispp ? 1: 4;
		int nrow = 9, ncol = 4;
		if(ispp ) { nrow =3; ncol =3;}
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol, 1000, 1000);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		for(int i=0; i<9 ; ++i){
				for(int j=0; j<ncent; ++j){
						h[i][j]->SetAxisRange(-1, 0.99, "X");
						h[i][j]->SetAxisRange(-1, 0.99, "Y");
						h[i][j]->GetXaxis()->SetNdivisions(505);
						h[i][j]->SetTitle("");
						if(ispp) cm->drawHist(h[i][j], i+1, "colz");
						else cm->drawHist(h[i][j], i+1, 4-j, "colz");
						if(ispp) cm->cd(i+1);
						else cm->CD(i+1, 4-j);
						if(doLogz)	gPad->SetLogz();
						tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
				}
		}
		cm->SaveAs(FigDumpPath+name+".gif");
}

void show_2D_sig(TString name, TString name1,bool doLogz=0,  bool isNumber = 0, bool ispp = 1){
		TH2D* sig[9][4];
		getSig(name1, sig, isNumber, ispp);
		show_2D_colz(name, sig, doLogz, ispp );
}

void show_2D_raw_sig(TString name, TString name1,bool doLogz=0, bool isNumber = 0){
		TH2D* sig[9][4];
		getRawSig(name1, sig, isNumber);
		show_2D_colz(name, sig, doLogz);
}

void getdEta(TString name, TH1D* h[9][4], bool isNumber = 0, bool ispp = 1){
		TH2D* sig[9][4];
		auto sp = new JTCSignalProducer();
		getSig(name, sig, isNumber, ispp);
		int ncent = ispp ? 1 : 4;
		for(int i=0; i<9 ; ++i){
				for(int j=0;j<ncent ; ++j){
						sp->sig=sig[i][j];
						h[i][j]= sp->getSignal_dEta(name+Form("_%d_%d",i, j));
				}
		}
}

void getRawDr(TString name, TH1D* h[9][4], bool isNumber = 0, bool ispp = 1){
		TH2D* sig[9][4];
		auto sp = new JTCSignalProducer();
		getRawSig(name, sig, isNumber, ispp);
		int ncent = ispp ? 1 : 4;
		for(int i=0; i<9 ; ++i){
				for(int j=0;j<ncent ; ++j){
						sp->sig=sig[i][j];
						h[i][j]= sp->doDrIntegral(name+Form("_raw_%d_%d",i, j));
						sp->doDrPhaseCorrection(sig[i][j], h[i][j]);
						//cout<<h[i][j]->GetName()<<endl;
				}
		}
}

void getDr(TString name, TH1D* h[9][4], bool isNumber = 0, bool ispp = 1){
		TH2D* sig[9][4];
		auto sp = new JTCSignalProducer();
		getSig(name, sig, isNumber, ispp);
		int ncent = ispp ? 1 : 4; for(int i=0; i<9 ; ++i){
				for(int j=0;j<ncent ; ++j){
						sp->sig=sig[i][j];
						h[i][j]= sp->doDrIntegral(name+Form("_%d_%d",i, j));
						sp->doDrPhaseCorrection(sig[i][j], h[i][j]);
						//cout<<h[i][j]->GetName()<<endl;
				}
		}
		//		delete sp;
}

mCanvasLoose *drawDrDist(TString name, TString name1, bool isNumber = 0, float x1=0, float x2 = -1){
		TH1D* h[9][4];
		getDr(name1, h, isNumber);
		return show_1D_pp(name, h, x1, x2);
}

void drawDrOverlay(TString name, TString name1, TString name2, TString lab1, TString lab2, float x1 =0, float x2=-1, bool isNumber = 0, bool ispp =1){
		TH1D* h[9][4], *h2[9][4];
		getDr(name1, h, isNumber, ispp);
		getDr(name2, h2, isNumber, ispp);
		show_1D_overlay(name, h, h2, lab1, lab2, x1, x2, ispp);
}


void drawdEtaOverlay(TString name, TString name1, TString name2, TString lab1, TString lab2, float x1 =-1.5, float x2=1.499, bool isNumber = 0, bool ispp = 1){
		TH1D* h[9][4], *h2[9][4];
		getdEta(name1, h, isNumber , ispp);
		getdEta(name2, h2, isNumber, ispp);
		show_1D_overlay(name, h, h2, lab1, lab2, x1, x2, ispp);
		TH1D* ratio[9][4]; 
		batch_ratio(name+"raito", h, h2, ratio, ispp);
		show_1D_pp("detaRatio_"+name, ratio);
}

void drRatio_pbOverpp(TString name, TString pbf, TString ppf, TH1D* ratio[9][4]){
		TH1D* ppdr[9][4], *pbdr[9][4];
		getDr(pbf, pbdr, 0, 0);
		getDr(ppf, ppdr, 0, 1);
		for(int i=0; i<nPt; ++i){
				for(int j=0; i<4; ++j){
						ratio[i][j]=(TH1D*)pbdr[i][j]->Clone(name+Form("_%d_%d", i,j ));
						ratio[i][j]->Divide(ppdr[i][0]);
				}
		}
}

mCanvasLoose *drRatio(TString name, TString name1, TString name2, bool isNumber = 0){
		TH1D* dr1[9][4], *dr2[9][4], *ratio[9][4];
		getDr(name1, dr1, isNumber);
		getDr(name2, dr2, isNumber);
		for(int i=0; i<9 ; ++i){
				ratio[i][0]=(TH1D*) dr1[i][0]->Clone(Form("ratio_%d",i));
				ratio[i][0]->Divide(ratio[i][0], dr2[i][0], 1, 1 , "B");
				//				ratio[i][0]->SetAxisRange(0.5, 1.5, "Y");
				ratio[i][0]->SetAxisRange(0, 0.99,"X");
		}
		return show_1D_pp(name, ratio);
}

float calDr(int i, int j, TH2D* h){
		float deta = h->GetXaxis()->GetBinCenter(i);
		float dphi = h->GetYaxis()->GetBinCenter(j);
		return sqrt(deta*deta+dphi*dphi);
}

void applyCorr(TString name, TString dataname, TString corrname, bool isNumber =0, bool ispp = 1, float drCut = 0.5){
		cout<<"applying correction.. "<<endl;
		TH2D* data[9][4], *corr[9][4], *h[9][4];
		getSig(dataname, data, isNumber, ispp);
		getSig(corrname, corr, isNumber, ispp);
		int ncent = ispp ? 1: 4;
		auto wf = new TFile(dataDumpPath+name+"_JTCSignal.root","recreate"); wf->cd();
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<ncent ; ++j){
						if( isNumber ) 	h[i][j]=(TH2D*) data[i][j]->Clone("signal_"+name+Form("_%d_%d",i,j));
						else h[i][j]=(TH2D*) data[i][j]->Clone("signal_"+name+Form("_pTweighted_%d_%d",i,j));
						for(int k=1; k<corr[i][j]->GetNbinsX()+1; ++k){
								for(int l=1; l<corr[i][j]->GetNbinsY()+1; ++l){
										if(calDr(k,l, corr[i][j] ) >drCut ) {
												corr[i][j]->SetBinContent(k,l,0);	
												corr[i][j]->SetBinError(k,l,0);	
										}
								}
						}
						h[i][j]->Add(corr[i][j],-1);
						h[i][j]->Write();
				}
		}
		wf->Close();
}

void applySpillOver(TString name, TString dataname, TString corrname, bool isNumber =0, bool ispp = 1){
		TH2D* data[9][4], *corr[9][4], *h[9][4];
		getSig(dataname, data, isNumber, ispp);
		getSig(corrname, corr, isNumber, ispp);
		int ncent = ispp ? 1: 4;
		auto wf = new TFile(dataDumpPath+name+"_JTCSignal.root","recreate"); wf->cd();
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<ncent ; ++j){
						if( isNumber ) 	h[i][j]=(TH2D*) data[i][j]->Clone("signal_"+name+Form("_%d_%d",i,j));
						else h[i][j]=(TH2D*) data[i][j]->Clone("signal_"+name+Form("_pTweighted_%d_%d",i,j));
						if(i<4) h[i][j]->Add(corr[i][j],-1);
						h[i][j]->Write();
				}
		}
		wf->Close();
}

void dumpDr_pp(TString name, TString dataname, bool isNumber =0 ){
		TH1D* dr[9][4];
		TH1D* dr_total[4];
		getDr(dataname, dr, isNumber);
		auto wf = new TFile(dataDumpPath+name+"_JTCProj.root", "recreate");
		wf->cd();
		int ncent=1;
		for(int j=0;j<ncent; ++j){
				dr_total[j]=(TH1D*) dr[0][j]->Clone(Form("dr_total_%d",j));
				for(int i=1; i<9; ++i){
						dr_total[j]->Add(dr[i][j]);
				}
				dr_total[j]->Write();
		}
		for(int i=0; i<9; ++i){
				dr[i][0]->Write();
		}
}

void JSRatio(TString name, TString ppname, TString pbname){
		TH1D* ppdr[9][4], *pbdr[9][4], *ratio[9][4];
		getDr(ppname, ppdr, 0);
		getDr(pbname, pbdr, 0, 0);

		for(int i=0; i<9; ++i){
				for(int j=0;j<4; ++j){
						cout<<i<<", "<<j << endl;
						ratio[i][j]=(TH1D*) pbdr[i][j]->Clone(Form("ratio_%d_%d", i,j));
						ratio[i][j]->Divide(ppdr[i][0]);
				}
		}
		show_1D_pp(name, ratio, 0, 0.99, 0);
}

TH1D* sumPt(TString name ,TH1D* h[9][4], TH1D* sum, int j=0){
		sum=(TH1D*) h[0][j]->Clone(name);
		for(int i=1; i<9; ++i){
				sum->Add(h[i][j]);
				cout<<"summing "<<h[i][j]->GetName()<<endl;
		}
		return sum;
}


void getJSClosure(TString name, TString ppname, TString pbname, TH1D* ratio[4]){
		TH1D* ppdr[9][4], *pbdr[9][4], *pp;
		getDr(ppname, ppdr, 0);
		getDr(pbname, pbdr, 0, 0);
		pp=sumPt(name+"_pp", ppdr, pp, 0);
		for(int i=0; i<4; ++i){
				ratio[i]=sumPt(name+Form("_pb_%d",i), pbdr, ratio[i], i);
				ratio[i]->Divide(pp);
				ratio[i]->SetAxisRange(0.,3.4, "Y");
				ratio[i]->SetNdivisions(505);
		}
}

void JSRatioTotal(TString name, TString ppname, TString pbname){
		TH1D* ratio[4];
		getJSClosure("ratio", ppname, pbname, ratio);
		auto cm  = new mCanvasLoose("cm_"+name, "", 1, 4);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		auto tl = new TLine(); tl->SetLineStyle(2);
		for(int i=0; i<4; ++i){
				cm->drawHist(ratio[i], 1, 4-i);
				cm->CD(1, 4-i);
				tl->DrawLine(0, 1, 1, 1);
				tx->DrawLatexNDC(0.1, 0.93, cent_label[i]);
		}
		cm->SaveAs(FigDumpPath+name+".gif");
}

void subComp_pp_pb1D(TString name, TString ppst, TString pbst){
		bool isNumber =0;
		TH1D* ppdr[9][4], *pbdr[9][4], *ppdeta[9][4], *pbdeta[9][4];
		getDr(ppst, ppdr, isNumber, 1);
		getDr(pbst, pbdr, isNumber, 0);
		getdEta(ppst, ppdeta, isNumber, 1);
		getdEta(pbst, pbdeta, isNumber, 0);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<4; ++j){
						pbdr[i][j]->Add(ppdr[i][0],-1);
						pbdr[i][j]->SetAxisRange(0., 0.99, "X");
						pbdeta[i][j]->Add(ppdeta[i][0],-1);
						pbdeta[i][j]->SetAxisRange(-1.5,1.499, "X");
				}
		}
		show_1D_pp("dr_"+name, pbdr, 0, 0.99, 0);
		show_1D_pp("deta_"+name, pbdeta, -1.5, 1.499, 0);
}

void subComp_pp_pb2D(TString name, TString ppst, TString pbst){
		bool isNumber =0;
		TH2D *pp[9][4], *pb[9][4];
		getSig(ppst, pp, isNumber, 1);
		getSig(pbst, pb, isNumber, 0);
		for(int i=0 ;i<nPt; ++i){
				for(int j=0; j<4; ++j){
						pb[i][j]->Add(pp[i][0],-1);
				}
		}
		show_2D_colz("logz_"+name, pb, 1, 0);
		show_2D_colz(name, pb, 0, 0);
}


mCanvasLoose * show_1D_pp_pb_overlay(TString name, TH1D* h[9][4], TH1D* h2[9][4], TString lab1, TString lab2, float x1=0, float x2=-1){
		auto cm = new mCanvasLoose("c_"+name, "", 9, 4);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		for(int i=0; i<9 ; ++i){
				for(int j=0; j<4 ; ++j){
						//						cout<<i<<", "<<j<<endl;
						h[i][0]->SetTitle("");
						h[i][0]->GetXaxis()->SetNdivisions(505);
						float ymin = min(h2[i][j]->GetMinimum(), h[i][0]->GetMinimum());
						float ymax = max(h2[i][j]->GetMaximum(), h[i][0]->GetMaximum());
						float range = ymax-ymin;
						h[i][0]->SetAxisRange(ymin-0.05*range, ymax+0.05*range, "Y");
						h2[i][j]->SetLineColor(kRed);
						if( x1< x2) h[i][0]->SetAxisRange(x1, x2,"X");
						cm->drawHist(h[i][0], i+1, 4-j);
						cm->drawHist(h2[i][j], i+1, 4-j);
						cm->CD(i+1, 4-j);
						//						gPad->SetLogy();
						tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
				}
		}
		//		auto tl = new TLegend(0.6, 0.6, 0.95, 0.87); tl->SetLineColor(0);
		auto tl = new TLegend(0.6, 0.2, 0.95, 0.5); tl->SetLineColor(0);
		//auto tl = new TLegend(0.6, 0.6, 0.95, 0.87); tl->SetLineColor(0);
		tl->AddEntry(h[0][0], lab1);
		tl->AddEntry(h2[0][0], lab2);
		cm->cd(1); tl->Draw();
		cm->SaveAs(FigDumpPath+name+".gif");
		return cm;
}

void show_1D_pp_pb_overlayRatio(TString name, TH1D* h[9][4], TH1D* h2[9][4], TString lab1, TString lab2, float x1=0, float x2=-1){
		auto cm = new doublePanelFig("c_"+name, "", 9, 4);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		TH1D* ratio[9][4];
		for(int i=0; i<9 ; ++i){
				for(int j=0; j<4 ; ++j){
						float ymin = min(h2[i][j]->GetMinimum(), h[i][0]->GetMinimum());
						float ymax = max(h2[i][j]->GetMaximum(), h[i][0]->GetMaximum());
						float range = ymax-ymin;
						//						cout<<i<<", "<<j<<endl;
						ratio[i][j]=(TH1D*) h2[i][j]->Clone(name+Form("ratio_%d_%d",i,j));
						ratio[i][j]->Divide(h[i][0]);
						ratio[i][j]->SetAxisRange(x1, x2, "X");
						//						ratio[i][j]->SetAxisRange(0, 5, "Y");
						h[i][0]->SetTitle("");
						h[i][0]->GetXaxis()->SetNdivisions(505);
						h[i][0]->SetAxisRange(ymin-0.05*range, ymax+0.05*range, "Y");
						h2[i][j]->SetLineColor(kRed);
						if( x1< x2) h[i][0]->SetAxisRange(x1, x2,"X");
						cm->addHist(h[i][0], i+1, 4-j);
						cm->addHist(h2[i][j], i+1, 4-j);
						cm->addHist(ratio[i][j], i+1, 4-j, 1);
						cm->CD(i+1, 4-j, 0);
						//						gPad->SetLogy();
						tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
						cm->CD(i+1, 4-j, 1);
						tx->DrawLatexNDC(0.2, 0.9, lab2+" / "+lab1);
				}
		}
		//		auto tl = new TLegend(0.6, 0.6, 0.95, 0.87); tl->SetLineColor(0);
		auto tl = new TLegend(0.6, 0.2, 0.95, 0.5); tl->SetLineColor(0);
		//auto tl = new TLegend(0.6, 0.6, 0.95, 0.87); tl->SetLineColor(0);
		tl->AddEntry(h[0][0], lab1);
		tl->AddEntry(h2[0][0], lab2);
		cm->CD(1, 1, 0); tl->Draw();
		cm->SaveAs(FigDumpPath+name+".gif");
}

void overlay_pb_pp(TString name, TString name1, TString name2, TString lab1, TString lab2, float x1=-1.5, float x2=1.499){
		TH1D* h[9][4], *h2[9][4];
		getdEta(name1, h, 0, 0);
		getdEta(name2, h2, 0, 1);
		show_1D_pp_pb_overlayRatio("dEta_"+name, h2, h, lab1, lab2, x1, x2);
		TH1D* dr1[9][4], *dr2[9][4];
		getDr(name1, dr1, 0, 0);
		getDr(name2, dr2, 0, 1);
		show_1D_pp_pb_overlayRatio("dr_"+name, dr2, dr1, lab2, lab1, 0, 0.99);
}

void mixTables(TString name, TString name1, TString name2, bool ispp=1){
		int isNumber = 0;
		TH2D* h1[9][4], *h2[9][4];
		getSig(name1, h1, isNumber, ispp);
		getSig(name2, h2, isNumber, ispp);
		int ncent = ispp ? 1 : 4;
		auto wf = new TFile(dataDumpPath+name+"_JTCSignal.root", "recreate"); wf->cd();
		TH2D* h;
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<ncent; ++j){
						if(i>6){ 
								h = h2[i][j];
						}
						else {
								h=h1[i][j];
						}
						if(!isNumber) tmp = "signal_"+name+Form("_pTweighted_%d_%d",i,j);
						else "signal_"+name+Form("_%d_%d",i,j);
						h->SetName(tmp);
						h->Write();
				}
		}
		wf->Close();
}

void fig5Stack(TString name, TH1D* h[3][4], bool isNumber =0 ){
		TH1D* dr[9][4];
		getDr(name, dr, isNumber , 0);		
		int bin[3] = {0, 2, 5};
		for(int i=0; i<3; ++i){
				for(int k=0; k<4; ++k){
						h[i][k]=(TH1D*) dr[bin[i]][k]->Clone(name+Form("_stacked_%d_%d",i,k));
						for(int j=bin[i]+1; j<9; ++j){
								h[i][k]->Add(dr[j][k]);
						}
				}
		}
}

void addStackLabel(mCanvasLoose* cp){
		TString label[]={"p_{T}^{track}>0.7", "p_{T}^{track}>2", "p_{T}^{track}>4"};
		TString cent[]={"Cent. 0-10%", "Cent. 10-30%", "Cent. 30-50%", "Cent. 50-100%"};
		auto tx = new TLatex();  tx->SetTextSize(.08);
		for(int j=0; j<4; ++j){
				cp->CD(1, 4-j);
				tx->DrawLatexNDC(0.6, 0.8, cent[j]);
		}
		for(int j=0; j<3; ++j){
				cp->CD(j+1, 1);
				tx->DrawLatexNDC(0.2, 0.8, label[j]);
		}
}

void ratioStack(TString name, TString name1, TString name2, bool isNumber = 0 ){
		TH1D* h1[3][4], *h2[3][4], *ratio[3][4];
		fig5Stack(name1, h1, isNumber);
		fig5Stack(name2, h2, isNumber);
		auto cp = new mCanvasLoose(name+"_st" , "", 3, 4);
		for(int j=0; j<3; ++j){
				for(int i=0; i<4; ++i){
						ratio[j][i]=(TH1D*) h1[j][i]->Clone(Form(name+"_ratio_%d_%d",i , j));
						ratio[j][i]->Divide(h2[j][i]);
						cp->drawHist(ratio[j][i], j+1, 4-i);
						cp->CD(j+1, 4-i);
						gPad->SetLogy();
				}
		}	   
		cp->SaveAs(FigDumpPath+name+".gif");
}

void readJSRes(TH1D* pp[9][4], bool ispp = 1){
		auto reff= new TFile("/Users/tabris/cmsProjects/iJTC/dataSet/inclJTC/nominalJSref.root");
		for(int i=0; i<9; ++i){
				pp[i][0]=(TH1D*)reff->Get(Form("JS_pp_%d",i));
		}
}

void showPlot2(TString name, bool ispp , float x1, float x2, int n_args, ...){
		va_list ap;
		va_start(ap, n_args);
		int ncol = 3, nrow = 3;
		if( !ispp ) {ncol= 4; nrow = 4;}
		auto cm = new mCanvasLoose("c_"+name, "", nrow, ncol);
		auto tx = new TLatex();  tx->SetTextSize(.08);
		int ncent = ispp ? 1: 4;
		TH1D** h;
		for(int k=0; k<n_args; ++k){
				h = va_arg(ap, TH1D**);
				for(int i=5; i<9 ; ++i){
int ii = i-5;
						for(int j=0; j<ncent; ++j){
						//		cout<<index_94(i,j)<<endl;
								h[index_94(i,j)]->SetTitle("");
								h[index_94(i,j)]->SetLineColor(color_vec[k]);
								h[index_94(i,j)]->GetXaxis()->SetNdivisions(505);
								h[index_94(i,j)]->SetAxisRange(x1, x2,"X");
								if(ispp ) { cm->drawHist(h[index_94(i,j)], ii+1);
										cm->cd(i+1);
										tx->DrawLatexNDC(0.2, 0.93, track_label[i]);
								}
								else { cm->drawHist(h[index_94(i,j)], ii+1, 4-j);
										cm->CD(i+1, 4-j);
										tmp = track_label[i]+" "+cent_label[j];
										tx->DrawLatexNDC(0.1, 0.93, tmp);
								}
						}
				}
		}
		cm->SaveAs(FigDumpPath+name+".gif");
}
