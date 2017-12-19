
#include "../lib/import_config.h"
#ifndef JTCSignalProducer_H
#include "../lib/JTCSignalProducer.h"
#endif
#ifndef multiPad_H
#include "../lib/multiPad.h"
#endif
#ifndef xPlotStyle_H
#include "../lib/xPlotStyle.h"
#endif

//const int nPt = 8; const int nCent =2 ;
TString dataDumpPath = "/Users/tabris/cmsProjects/iJTC/dataSet/correlation/";

TString trk_tag[] = {"0.7 < p_{T}^{track} < 1 GeV","1 < p_{T}^{track} < 2 GeV",
		"2 < p_{T}^{track} < 3 GeV", "3 < p_{T}^{track} < 4 GeV",
		"4 < p_{T}^{track} < 8 GeV", "8 < p_{T}^{track} < 12 GeV", 
		"12 < p_{T}^{track} < 16 GeV", "p_{T}^{track} > 16 Gev"};
TString cent_label[]={"Cent. 0-30%", "Cent. 30-100%"};
namespace input_raw2D{
		void pullSig(TString fname){
				cout<<"pulling singal for "<<fname;
				TString trkbin [] = {"0.7", "1", "2", "3", "4", "8", "12", "16", "999"};
				TString centbin [] = {"0", "30", "100"};
				JTCSignalProducer *sp1[8][2];
				JTCSignalProducer *sp2[8][2];
				TFile * wf = new TFile(dataDumpPath+fname+"_JTCSignal.root","recreate");
				TString title;
				cout<<".";
				for(int i=0; i<nPt ; ++i){
						cout<<".";
						for(int j=0; j<nCent ; ++j){
								title = "track pt in ["+trkbin[i]+", "+trkbin[i+1]\
										 +"), cent in ["+centbin[j]+", "+centbin[j+1]+"), jet pt in [120, 1000]";
								raw_sig[i][j]->SetTitle(title);
								raw_sig_pTweighted[i][j]->SetTitle(title);
								mixing [i][j]->SetTitle(title);
								sp1[i][j] = new JTCSignalProducer(raw_sig[i][j], mixing[i][j]);
								sp1[i][j]->getSignal(fname+Form("_%d_%d", i, j));
								sp1[i][j]->WriteTH2();
								sp2[i][j] = new JTCSignalProducer(raw_sig_pTweighted[i][j], mixing[i][j]);
								sp2[i][j]->getSignal(fname+Form("_pTweighted_%d_%d", i, j));
								sp2[i][j]->WriteTH2();
//								cout<<sp2[i][j]->sig->GetName()<<endl;
						}
				}
				cout<<". Done"<<endl;;
				wf->Close();
				//clearInput();
				cout<<"signal dumped to "<<dataDumpPath+fname+"_JTCSignal.root"<<endl;
		}
		void showSpectra(TString name1, TString name2, TFile *f1, TFile *f2){
				auto *cp = new doublePanelFig("c_"+name1+"_"+name2, "", 1, nCent );
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				TString tmp;
				TH1D* hr[nCent];
				TH1D* h[nCent][2];
				for(int j=0; j<nCent; ++j){
						tmp = name1+"_all_bjets_corrpT"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt120_Pt1000";
						cout<<tmp<<endl;
						h[j][0] =(TH1D*) f1->Get(tmp); h[j][0]->GetXaxis()->SetTitle("#Delta#eta"); h[j][0]->SetTitle("");
						h[j][0]->SetMarkerColor(kBlue+2);
						h[j][0]->SetLineWidth(1);
						hr[j]=(TH1D*)h[j][0]->Clone(Form("hr_%d",j));
						h[j][0]->SetAxisRange(h[j][0]->GetMinimum()-h[j][0]->GetMaximum()*0.1, h[j][0]->GetMaximum()*1.5, "Y");
						h[j][0]->SetAxisRange(100, 500 ,"X");
						h[j][0]->SetAxisRange(1e-13, 1e-4,"Y");
						cp->addHist(h[j][0], 1, 2-j);
						cp->CD(1, 2-j, 0); gPad->SetLogy();
						tmp = name2+"_all_bjets_corrpT"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt120_Pt1000";
						h[j][1] =(TH1D*) f2->Get(tmp); h[j][1]->SetLineColor(kRed);  h[j][1]->SetMarkerColor(kRed);
						cp->addHist(h[j][1], 1, 2-j); 
						hr[j]->Divide(h[j][1]); hr[j]->GetYaxis()->SetNdivisions(505);
						hr[j]->SetAxisRange(100, 500, "X");  
						cp->addHist(hr[j], 1, 2-j, 1);
						cp->CD(1, 2-j, 0); tl->DrawLine(-2.5, 0, 2.5, 0);
						tmp = "Jet spectra: "+cent_label[j];
						tx->DrawLatexNDC(0.15,0.87, tmp); 
						cp->CD(1, 2-j, 1); tl->DrawLine(-2.5, 1, 2.5, 1);
						cp->draw95Area(1,2-j, 100, 500);
				}
				cp->SaveAs(name1+"_"+name2+"_jetSpectra.gif");
		}

}

namespace signal2D{
		void drawTableWithRatio(TString name1, TString name2, TFile *f1, TFile *f2, bool isNumber = 1){
				cout<<"drawing the ratio of "<<name1<<" over "<<name2<<endl;	
				JTCSignalProducer *sp1 [nPt][nCent];
				JTCSignalProducer *sp2 [nPt][nCent];
				auto *cp = new doublePanelFig("cp_"+name1+"_"+name2, "", nPt,nCent );
				auto *cp2 = new doublePanelFig("cp2_"+name1+"_"+name2, "", nPt,nCent );
				TH1D* hr[nPt][nCent];
				TH1D* hr2[nPt][nCent];
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				TString tmp;
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								cout<<i<<", "<<j<<endl;
								sp1[i][j]= new JTCSignalProducer();
								sp2[i][j]= new JTCSignalProducer();
								if(!isNumber) {sp1[i][j]->read(f1, name1+Form("_pTweighted_%d_%d", i,j));
										sp2[i][j]->read(f2, name2+Form("_pTweighted_%d_%d", i,j));}
								else {  sp1[i][j]->read(f1, name1+Form("_%d_%d", i,j));
										sp2[i][j]->read(f2, name2+Form("_%d_%d", i,j));}
								TH1* h=sp1[i][j]->projX(1,sp1[i][j]->raw_sig,-1,1); h->GetXaxis()->SetTitle("#Delta#eta"); h->SetTitle("");
								TH1* h2 = sp1[i][j]->signal_X(); h2->GetXaxis()->SetTitle("#Delta#eta"); h2->SetTitle("");
								if(i==0 ){
										h->SetAxisRange(h->GetMinimum()-h->GetMaximum()*0.1, h->GetMaximum()*1.5, "Y");
										h2->SetAxisRange(h2->GetMinimum()-h2->GetMaximum()*0.1, h2->GetMaximum()*1.5, "Y");
								}
								else {
										h->SetAxisRange(h->GetMinimum()-h->GetMaximum()*0.1, h->GetMaximum()*1.4, "Y");
										h2->SetAxisRange(h2->GetMinimum()-h2->GetMaximum()*0.1, h2->GetMaximum()*1.4, "Y");
								}
								h->SetMarkerColor(kBlue+2);
								h->SetLineWidth(1);
								h->SetAxisRange(-2, 1.99, "X");
								h2->SetAxisRange(-2, 1.99, "X");
								hr[i][j]=(TH1D*)h->Clone(Form("hr_%d_%d",i,j));
								hr2[i][j]=(TH1D*)h2->Clone(Form("hr2_%d_%d",i,j));
								cp->addHist(h, i+1, 2-j);
								cp2->addHist(h2, i+1, 2-j);
								h=sp2[i][j]->projX(1,sp2[i][j]->raw_sig,-1,1); h->SetLineColor(kRed);  h->SetMarkerColor(kRed);
								h2 = sp2[i][j]->signal_X(); h2->SetLineColor(kRed);  h2->SetMarkerColor(kRed);
								cp->addHist(h, i+1, 2-j); 
								cp2->addHist(h2, i+1, 2-j); 
								hr[i][j]->Divide(h); hr[i][j]->GetYaxis()->SetNdivisions(505);
								hr[i][j]->SetAxisRange(0, 2, "Y");  cp->addHist(hr[i][j], i+1, 2-j, 1);
								hr2[i][j]->Divide(h2); hr2[i][j]->GetYaxis()->SetNdivisions(505);
								hr2[i][j]->SetAxisRange(0, 2, "Y");  cp2->addHist(hr2[i][j], i+1, 2-j, 1);
								cp->CD(i+1, 2-j, 0); tl->DrawLine(-2., 0, 2., 0);
								cp2->CD(i+1, 2-j, 0); tl->DrawLine(-2., 0, 2., 0);
								tmp = trk_tag[i]+", "+cent_label[j];
								tx->DrawLatexNDC(0.15,0.87, tmp); 
								cp->CD(i+1, 2-j, 1); tl->DrawLine(-2., 1, 2., 1);
								cp->draw95Area(i+1,2-j, -2, 2);
								cp2->CD(i+1, 2-j, 1); tl->DrawLine(-2., 1, 2., 1);
								cp2->draw95Area(i+1,2-j, -2, 2);
						}
				}
				if(isNumber) {
						cp2->SaveAs(name1+"_"+name2+"_overlay.gif");
						cp->SaveAs(name1+"_"+name2+"_raw_overlay.gif");
				}
				else {
						cp2->SaveAs(name1+"_"+name2+"_pTweighted_overlay.gif");
						cp->SaveAs(name1+"_"+name2+"_pTweighted_raw_overlay.gif");
				}
		}

		void drawTable(TString name1, TString name2, TFile *f1, TFile *f2){
				JTCSignalProducer *sp1 [nPt][nCent];
				JTCSignalProducer *sp2 [nPt][nCent];
				multiPad *cp = new multiPad("cp_"+name1+"_"+name2, "", nPt,nCent );
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								sp1[i][j]= new JTCSignalProducer();
								sp1[i][j]->read(f1, name1+Form("_%d_%d", i,j));
								sp2[i][j]= new JTCSignalProducer();
								sp2[i][j]->read(f2, name2+Form("_%d_%d", i,j));
								TH1* h = sp1[i][j]->signal_X(); h->GetXaxis()->SetTitle("#Delta#eta"); h->SetTitle("");
								if(i==0 )h->SetAxisRange(h->GetMinimum()-h->GetMaximum()*0.1, h->GetMaximum()*1.5, "Y");
								else h->SetAxisRange(h->GetMinimum()-h->GetMaximum()*0.1, h->GetMaximum()*1.4, "Y");
								cp->addHist(h, i+1, 2-j);
								h = sp2[i][j]->signal_X(); h->SetLineColor(kRed);
								cp->addHist(h, i+1, 2-j);
						}
				}
				cp->addLine(-2.5, 0, 2.5, 0);
				cp->tl->SetLineStyle(3);
				cp->draw();
				auto tx = new TLatex();  tx->SetTextSize(.08);
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								cp->cd(i*2+2-j);
								if( j==1) tx->DrawLatexNDC(0.2,0.87, trk_tag[i]);
								if( i==0) tx->DrawLatexNDC(0.6,0.87, cent_label[j]);
						}
				}
				cp->SaveAs(name1+"_"+name2+"_overlay.gif");
		}

		void drawJetShapeRatio(TString name1, TString name2, TFile *f1, TFile *f2){
				float binwidth [] = {0.3, 1, 1,1, 4, 4, 4, 1};
				auto *cp = new doublePanelFig("cj_"+name1+"_"+name2, "",nCent, nPt);
				JTCSignalProducer *sp1 [nPt][nCent];
				JTCSignalProducer *sp2 [nPt][nCent];
				TH1D* hr[nPt][nCent];
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				TString tmp;
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								sp1[i][j]= new JTCSignalProducer();
								sp2[i][j]= new JTCSignalProducer();
								sp1[i][j]->read(f1, name1+Form("_pTweighted_%d_%d", i,j));
								sp2[i][j]->read(f2, name2+Form("_pTweighted_%d_%d", i,j));
								sp1[i][j]->sig->Scale(1.0/binwidth[i]);
								sp2[i][j]->sig->Scale(1.0/binwidth[i]);
								sp1[i][j]->doDrIntegral(name1+Form("_%d_%d", i, j));
								sp2[i][j]->doDrIntegral(name2+Form("_%d_%d", i, j));
								TH1* h=sp1[i][j]->dr_integral; // h->SetTitle("");
								h->SetAxisRange(h->GetMinimum()-h->GetMaximum()*0.1, h->GetMaximum()*1.5, "Y");
								h->SetMarkerColor(kBlue+2);
								h->SetLineWidth(1);
								hr[i][j]=(TH1D*)h->Clone(Form("hr_%d_%d",i,j)); 
								cp->addHist(h, j+1, i+1);
								h=sp2[i][j]->dr_integral; h->SetLineColor(kRed);  h->SetMarkerColor(kRed);
								cp->addHist(h, j+1, i+1); 
								hr[i][j]->Divide(h); hr[i][j]->GetYaxis()->SetNdivisions(505);
								hr[i][j]->GetXaxis()->SetTitle("#Deltar");
								hr[i][j]->GetXaxis()->SetTitleColor(kBlack);
								hr[i][j]->SetAxisRange(0, 2, "Y");  cp->addHist(hr[i][j], j+1,i+1, 1);
								cp->CD(j+1, i+1, 0); tl->DrawLine(0, 0, 1, 0);
								tmp = trk_tag[i]+", "+cent_label[j];
								tx->DrawLatexNDC(0.15,0.87, tmp); 
								cp->CD(j+1, i+1, 1); tl->DrawLine(0, 1, 1, 1);
								cp->draw95Area(j+1, i+1, 0, 1);
						}
				}
				cp->SaveAs("jet_shape_"+name1+"_"+name2+"_overlay.gif");
		}
}
