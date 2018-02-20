
#include "../../lib/import_pf_config.h"
#include "../../lib/JTCSignalProducer.h"
#ifndef xPlotStyle_H
#include "../../lib/xPlotStyle.h"
#endif

TString dataDumpPath = "/Users/tabris/cmsProjects/iJTC/dataSet/correlation/pfJetCorrelation/";
TString FigDumpPath  = "/Users/tabris/cmsProjects/iJTC/macros/fig_pfJets/";

TString track_label[] = {"0.7 < p_{T}^{track} < 1 GeV","1 < p_{T}^{track} < 2 GeV",
		"2 < p_{T}^{track} < 3 GeV", "3 < p_{T}^{track} < 4 GeV",
		"4 < p_{T}^{track} < 8 GeV", "8 < p_{T}^{track} < 12 GeV", 
		"p_{T}^{track} > 12 GeV"};
//		"12 < p_{T}^{track} < 16 GeV", "p_{T}^{track} > 16 Gev"};
TString cent_label[]={"Cent. 0-30%", "Cent. 30-100%"};
TString tmp;

const int nPt = 7;
const int nCent=2;
#include "anaFig_recoCheck.h"

//JTCSignalProducer*** batchRead1D(TString name, JTCSignalProducer* sp[8][2], bool isNumber = 1){
JTCSignalProducer*** batchRead1D(TString name,bool isNumber = 0){
		TFile *f = TFile::Open(dataDumpPath+name+"_JTCProj.root");
		JTCSignalProducer*** sp = new JTCSignalProducer**[8];
		for(int i=0; i<8; ++i){
				sp[i]=new JTCSignalProducer*[2];
				for(int j=0; j<2; ++j){
						sp[i][j]=new JTCSignalProducer();
						if(isNumber) sp[i][j]->read1D(f, name+Form("_%d_%d", i, j));
						else sp[i][j]->read1D(f, name+Form("_pTweighted_%d_%d", i, j));
				}
		}
		return sp;
}

void getSig(TString name, TH2D* h[8][2], bool isNumber=1){
		TFile *f = TFile::Open(dataDumpPath +name+"_JTCSignal.root");
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						h[i][j]=(TH2D*) f->Get("signal_"+name+Form("_%d_%d",i,j));
				}
		}
}

void getSig(TString name, TH2D* h[8][2], TH2D* h2[8][2]){
		TFile *f = TFile::Open(dataDumpPath +name+"_JTCSignal.root");
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						h[i][j]=(TH2D*) f->Get("signal_"+name+Form("_%d_%d",i,j));
						//h2[i][j]=(TH2D*) f->Get("signal_"+name+Form("_pTweighted_%d_%d",i,j));
	//					cout<<h[i][j]->GetName()<<endl;
				}
		}
}

void pullSig(TString fname, float sidemin=1.5, float sidemax=2.5, bool doSeagull=0){
		cout<<"pulling singal for "<<fname;
		TString trkbin [] = {"0.7", "1", "2", "3", "4", "8", "12", "16", "999"};
		TString centbin [] = {"0", "30", "100"};
		JTCSignalProducer *sp1[8][2];
		JTCSignalProducer *sp2[8][2];
		TFile * wf = new TFile(dataDumpPath+fname+"_JTCSignal.root","recreate");
		tmp;
		cout<<".";
		for(int i=0; i<nPt ; ++i){
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

		clean_pf_input();
}

void pull1D(TString fname, TFile *f, bool rebin=1){
		cout<<"pulling 1D histograms for "<<fname;
		JTCSignalProducer *sp1[8][2];
		JTCSignalProducer *sp2[8][2];
		TFile * wf = new TFile(dataDumpPath+fname+"_JTCProj.root","recreate");
		cout<<".";
		for(int i=0; i<nPt ; ++i){
				cout<<".";
				for(int j=0; j<nCent ; ++j){
						sp1[i][j] = new JTCSignalProducer();
						sp2[i][j] = new JTCSignalProducer();
						sp1[i][j]->read(f, fname+Form("_%d_%d", i, j));
						sp1[i][j]->getAllProj(fname+Form("_%d_%d", i, j), rebin);
						sp1[i][j]->WriteTH1();
						sp2[i][j]->read(f, fname+Form("_pTweighted_%d_%d", i, j));
						sp2[i][j]->getAllProj(fname+Form("_pTweighted_%d_%d", i, j), rebin);
						sp2[i][j]->WriteTH1();
						delete sp1[i][j];
						delete sp2[i][j];
				}
		}
		cout<<". Done"<<endl;;
		wf->Close();
		//clearInput();
		cout<<"signal dumped to "<<dataDumpPath+fname+"_JTCProj.root"<<endl;
}


void pull1D(TString fname, bool rebin=1){
		TFile *f = TFile::Open(dataDumpPath+fname+"_JTCSignal.root");
		cout<<f->GetName()<<endl;
		pull1D(fname, f, rebin);
}

void checkBkg(TString name, int ispb = 0){
		ispb++;
		tmp =dataDumpPath+name+"_JTCProj.root";
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
						sp1[i][j]->read1D(f, name+Form("_%d_%d", i, j));
						sp2[i][j]->read1D(f, name+Form("_pTweighted_%d_%d", i, j));
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
		cp1->SaveAs(FigDumpPath+"bkgCheck_"+name+".pdf");
		cp2->SaveAs(FigDumpPath+"bkgCheck_JS_"+name+".pdf");
		cp3->SaveAs(FigDumpPath+"bkgCheck_sideBand_"+name+".pdf");
		cp4->SaveAs(FigDumpPath+"bkgCheck_sideBand_JS_"+name+".pdf");
		for(int i=0; i<nPt ; ++i){
				for(int j=0; j<ispb ; ++j){
						delete sp1[i][j];
						delete sp2[i][j];
				}
		}
}

void spectraRatio(TString name, TString name1, TString name2, TFile *f1, TFile *f2, bool ispp = 0){
		auto *cp = new doublePanelFig("c_"+name1+"_"+name2, "", 1, nCent );
		auto tx = new TLatex();  tx->SetTextSize(.08);
		auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
		TString tmp;
		TH1D* hr[nCent];
		TH1D* h[2][nCent];
		auto sp = new JTCSignalProducer();
		int nbin = 21;
		const Double_t newbin [22] = {100, 120, 136, 152, 168, 184, 200, 216, 232, 248, 264, 280, 296, 312, 328, 344, 360,
				380, 400, 432, 500, 600};
		getSpectra(f1, name1, h[0], "sp1");
		getSpectra(f2, name2, h[1], "sp2");
		for(int j=0; j<nCent; ++j){
				h[0][j]->GetXaxis()->SetTitle("p_{T}^{jet}"); h[j][0]->SetTitle("");
				h[0][j]->SetMarkerColor(kBlue+2);
				h[0][j]->SetLineWidth(1);
				hr[j]=(TH1D*)h[0][j]->Clone(Form("hr_%d",j));
				h[1][j]->SetAxisRange(100, 500 ,"X");
				cp->CD(1, 2-j, 0); gPad->SetLogy();
				h[1][j]->SetLineColor(kRed);  h[1][j]->SetMarkerColor(kRed);
				h[1][j]->SetAxisRange(1e-10, 1e-5 ,"Y");
				cp->addHist(h[1][j], 1, 2-j); 
				cp->addHist(h[0][j], 1, 2-j);
				hr[j]->Divide(h[1][j]); hr[j]->GetYaxis()->SetNdivisions(505);
				hr[j]->SetAxisRange(100, 500, "X");
				hr[j]->GetXaxis()->SetNdivisions(505);	
				hr[j]->SetAxisRange(0., 1, "Y");	
				cp->addHist(hr[j], 1, 2-j, 1);
				cp->CD(1, 2-j, 0); tl->DrawLine(-2.5, 0, 2.5, 0);
				if( ispp ) 
						tmp = "";
				else
						tmp = "Jet spectrum: "+cent_label[j];
				tx->DrawLatexNDC(0.15,0.87, tmp); 
				cp->CD(1, 2-j, 1); tl->DrawLine(-2.5, 1, 2.5, 1);
				float frac = h[0][j]->Integral()/h[1][j]->Integral();
				cout<<"fraction: "<<frac*100<<"%"<<endl;
		}
		cp->SaveAs(FigDumpPath+name+"_jetSpectra.gif");
}

void drawRecoCheck(TString name1, TString name2, bool ver=1){
		recCheck_pb("recoCheck_"+name1+"_over_"+name2, name1, name2, ver);
}
/*
void getHist(TString name, TH1* h[8][2]){
		TString tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/"+name+"_drIntegral.root";
		TFile*f = TFile::Open(tmp);
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						h[i][j]=(TH1*) f->Get(name+Form("_%d_%d",i,j));
				}
		}
}
*/

void getDr(TString name, TH1D* h[8][2], bool isNumber = 1){
		TString tmp = dataDumpPath+name+"_JTCProj.root";
		TFile *f = TFile::Open(tmp);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<2; ++j){
						if(isNumber) tmp = "dr_"+name+Form("_%d_%d", i,j);
						else tmp = "dr_"+name+Form("_pTweighted_%d_%d", i,j);
						h[i][j]=(TH1D*)f->Get(tmp)->Clone(tmp);
				}
		}
//		f->Close();
}

void JSSubtraction(TString name1, TString name2, TString wname,  float frac){
		TString tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/pfJetCorrelation/"+name1+"_JTCSignal.root";
		TFile *f1 = TFile::Open(tmp);
		tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/"+name2+"_JTCSignal.root";
		TFile *f2 = TFile::Open(tmp);
		JTCSignalProducer* sp1[8][2];
		JTCSignalProducer* sp2[8][2];
		JTCSignalProducer* sp3[8][2];
		JTCSignalProducer* sp4[8][2];
		TH1* hjs[8][2];
		TH1* hpy[8][2];
		TH1* hpy_ref[8][2];
		TH1* hjs_ref[8][2];
		float comp = 1-frac;
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						sp1[i][j] = new JTCSignalProducer();
						sp2[i][j] = new JTCSignalProducer();
						sp1[i][j]->read(f1, name1+Form("_%d_%d", i,j));
						sp2[i][j]->read(f1, name1+Form("_pTweighted_%d_%d", i,j));
						sp1[i][j]->doDrIntegral(name1+Form("_%d_%d", i,j));
						sp2[i][j]->doDrIntegral(name1+Form("_pTweighted_%d_%d",i,j));

						sp3[i][j] = new JTCSignalProducer();
						sp4[i][j] = new JTCSignalProducer();
						sp3[i][j]->read(f2, name2+Form("_%d_%d", i,j));
						sp4[i][j]->read(f2, name2+Form("_pTweighted_%d_%d", i,j));
						sp3[i][j]->doDrIntegral(name2+Form("_%d_%d", i,j));
						sp4[i][j]->doDrIntegral(name2+Form("_pTweighted_%d_%d",i,j));

						hjs_ref[i][j]=(TH1D*)sp2[i][j]->dr_integral->Clone(Form("JS_ref_%d_%d",i,j));
						hpy_ref[i][j]=(TH1D*)sp1[i][j]->dr_integral->Clone(Form("PY_ref_%d_%d",i,j));
						sp3[i][j]->dr_integral->Scale(comp);
						sp4[i][j]->dr_integral->Scale(comp);
						hpy[i][j]=(TH1D*) sp1[i][j]->dr_integral->Clone(wname+Form("_%d_%d",i,j));
						hjs[i][j]=(TH1D*) sp2[i][j]->dr_integral->Clone(wname+Form("_pTweighted_%d_%d",i,j));
						hpy[i][j]->Add(sp3[i][j]->dr_integral,-1);
						hpy[i][j]->Scale(1.0/frac);
						hjs[i][j]->Add(sp4[i][j]->dr_integral,-1);
						hjs[i][j]->Scale(1.0/frac);
				}
		}
		TFile* wf = new TFile(dataDumpPath+wname+"_drIntegral.root", "recreate");

		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						hpy[i][j]->Write();
						hjs[i][j]->Write();
				}
		}
		wf->Close();

		//quickRatio(wname, hjs, hjs_ref, 0.0, .99);
		/*
		   cpdr = new mCanvasLoose("cpdr",  "", 7, 2);
		   for(int i=1; i<8; ++i){
		   for(int j=0; j<2; ++j){
		   cpdr->CD(i, 2-j); 
		   }
		   }
		   */
}


void adjustYrange(TH1* h ){
		float range = h->GetMaximum()-h->GetMinimum();
		h->SetAxisRange(h->GetMinimum()-0.1*range, h->GetMaximum()+0.12*range, "Y");
}


