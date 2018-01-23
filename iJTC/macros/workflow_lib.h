
#include "../lib/import_config.h"
#include "../lib/import_inclusive.h"
#ifndef JTCSignalProducer_H
#include "../lib/JTCSignalProducer.h"
#endif
#ifndef multiPad_H
#include "../lib/multiPad.h"
#endif
#ifndef xPlotStyle_H
#include "../lib/xPlotStyle.h"
#endif

#include "anaFig_stack.h"
#include "utility.h"
// namespace :
//  utility
//  input_raw2D
//  signal2D
//  signal1D
//  inclusive_input
//  ana_fig
//const int nPt = 8; const int nCent =2 ;



namespace ana_fig{
		void closure(TString name, TString name1, TString name2, TFile *f1, TFile *f2, bool isNumber = 1){
				// need to add the systematic uncertainty
				auto cp1 = new mCanvasLoose("cp1", "", 7, 2, 300, 125);
				TH1D * hr [8][2];
				TH1D * hr2[8][2];
				TString tmp;
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				JTCSignalProducer *sp1[8][2];
				JTCSignalProducer *sp2[8][2];
				for(int i=1; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								sp1[i][j] = new JTCSignalProducer();
								sp2[i][j] = new JTCSignalProducer();
								if(!isNumber) {sp1[i][j]->read(f1, name1+Form("_pTweighted_%d_%d", i,j));
										sp2[i][j]->read(f2, name2+Form("_pTweighted_%d_%d", i,j));}
								else {  sp1[i][j]->read(f1, name1+Form("_%d_%d", i,j));
										sp2[i][j]->read(f2, name2+Form("_%d_%d", i,j));}
								sp1[i][j]->doDrIntegral(name1+Form("_%d_%d", i, j));
								sp2[i][j]->doDrIntegral(name2+Form("_%d_%d", i, j));
								hr[i][j]=(TH1D*)sp1[i][j]->dr_integral->Clone(Form("ratio_%d_%d",i,j));
								hr[i][j]->Divide(sp2[i][j]->dr_integral);
								hr2[i][j]=(TH1D*)hr[i][j]->Clone(Form("ratio2_%d_%d",i,j));

								tmp = trk_tag[i]+", "+cent_label[j];
								cp1->CD(i, 2-j);
								hr[i][j]->SetTitle("");
								hr[i][j]->Draw();
								tx->DrawLatexNDC(0.15,0.87, tmp); 
						}
				}
				cp1->SaveAs(name+"_closure.gif");
		}
}

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
		void showSpectra(TString name, TString name1, TString name2, TFile *f1, TFile *f2){
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
			//			h[j][0]->Scale(1.0/h[j][0]->Integral("width"));
						h[j][0]->SetMarkerColor(kBlue+2);
						h[j][0]->SetLineWidth(1);
						hr[j]=(TH1D*)h[j][0]->Clone(Form("hr_%d",j));
//						h[j][0]->SetAxisRange(1e-6, h[j][0]->GetMaximum()*10, "Y");
						h[j][0]->SetAxisRange(100, 500 ,"X");
						//						h[j][0]->SetAxisRange(1e-13, 1e-4,"Y");
						cp->addHist(h[j][0], 1, 2-j);
						cp->CD(1, 2-j, 0); gPad->SetLogy();
						tmp = name2+"_all_bjets_corrpT"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt120_Pt1000";
						h[j][1] =(TH1D*) f2->Get(tmp); h[j][1]->SetLineColor(kRed);  h[j][1]->SetMarkerColor(kRed);
			//			h[j][1]->Scale(1.0/h[j][1]->Integral("width"));
						cp->addHist(h[j][1], 1, 2-j); 
						hr[j]->Divide(h[j][1]); hr[j]->GetYaxis()->SetNdivisions(505);
						hr[j]->SetAxisRange(100, 500, "X");  
						cp->addHist(hr[j], 1, 2-j, 1);
						cp->CD(1, 2-j, 0); tl->DrawLine(-2.5, 0, 2.5, 0);
						tmp = "Jet spectra: "+cent_label[j];
						tx->DrawLatexNDC(0.15,0.87, tmp); 
						cp->CD(1, 2-j, 1); tl->DrawLine(-2.5, 1, 2.5, 1);
			//			cp->draw95Area(1,2-j, 100, 500);
				}
				cp->SaveAs(FigDumpPath+name+"_jetSpectra.gif");
		}

}

namespace signal2D{
		void drawJSratio(TString name1, TString name2, TString name, bool ispp=0){
				TString tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/"+name1+"_JTCSignal.root";
				TFile *f1 = TFile::Open(tmp);
				tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/"+name2+"_JTCSignal.root";
				TFile *f2 = TFile::Open(tmp);
				JTCSignalProducer *sp1[8][2];
				JTCSignalProducer *sp2[8][2];
				TFile* wf = TFile::Open("debug.root","recreate");
				for(int i=0; i<8; ++i){
						for(int j=0; j<2; ++j){
								sp1[i][j]=new JTCSignalProducer();
								sp2[i][j]=new JTCSignalProducer();
								sp1[i][j]->read(f1, name1+Form("_pTweighted_%d_%d", i,j));
								sp2[i][j]->read(f2, name2+Form("_pTweighted_%d_%d", i,j));
								sp1[i][j]->doDrIntegral(Form("_1_%d_%d", i, j));
								sp2[i][j]->doDrIntegral(Form("_2_%d_%d", i, j));
						}
				}
				utility::quickJSratio(name, sp1, sp2, 0, ispp);
		}

		void drawStackJSDiff(TString name1, TString name2, TString name, bool isNumber = 1){
				TFile *f1 = TFile::Open(dataDumpPath+name1+"_JTCSignal.root");
				TFile *f2 = TFile::Open(dataDumpPath+name2+"_JTCSignal.root");
				JTCSignalProducer *sp1[8][2];
				JTCSignalProducer *sp2[8][2];
				TString tmp;
				TH1* hpb[8][2];
				TH1* hpp[8];
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								sp1[i][j] = new JTCSignalProducer();
								sp2[i][j] = new JTCSignalProducer();
								if(!isNumber) {sp1[i][j]->read(f1, name1+Form("_pTweighted_%d_%d", i,j));
										sp2[i][j]->read(f2, name2+Form("_pTweighted_%d_%d", i,j));}
								else {  sp1[i][j]->read(f1, name1+Form("_%d_%d", i,j));
										sp2[i][j]->read(f2, name2+Form("_%d_%d", i,j));}
								sp1[i][j]->doDrIntegral(name1+Form("_%d_%d", i, j));
								sp2[i][j]->doDrIntegral(name2+Form("_%d_%d", i, j));
								hpb[i][j]=sp1[i][j]->dr_integral;
								if(j==0)	hpp[i]=sp2[i][j]->dr_integral;
						}
				}
				tmp = FigDumpPath+"stack_plot_"+name+".pdf";
				stackPlot_diff(hpb, hpp, tmp);
		}

		void pull1D(TString fname, TFile *f){
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
								sp1[i][j]->getAllProj(fname+Form("_%d_%d", i, j));
								sp1[i][j]->WriteTH1();
								sp2[i][j]->read(f, fname+Form("_pTweighted_%d_%d", i, j));
								sp2[i][j]->getAllProj(fname+Form("_pTweighted_%d_%d", i, j));
								sp2[i][j]->WriteTH1();
						}
				}
				cout<<". Done"<<endl;;
				wf->Close();
				//clearInput();
				cout<<"signal dumped to "<<dataDumpPath+fname+"_JTCProj.root"<<endl;
		}


		void pull1D(TString fname){
				TFile *f = TFile::Open(dataDumpPath+fname+"_JTCSignal.root");
				cout<<f->GetName()<<endl;
				pull1D(fname, f);
		}

}

namespace signal1D {
		void checkSide(TString name, int ncent = 2){
				TString tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/"+name+"_JTCProj.root";
				TFile *f = TFile::Open(tmp);
				JTCSignalProducer *sp1[8][2];
				TH1* h1[8][2], *h2[8][2];
				for(int i=0; i<nPt ; ++i){
						for(int j=0; j<ncent ; ++j){
								sp1[i][j] = new JTCSignalProducer();
								sp1[i][j]->read1D(f, name+Form("_%d_%d", i, j));
								h1[i][j]= (TH1*)sp1[i][j]->side_deta;
								h2[i][j]= (TH1*)sp1[i][j]->side_deta_mix;
								h1[i][j]->Scale(1.0/h1[i][j]->Integral("width"));
								h2[i][j]->Scale(1.0/h1[i][j]->Integral("width"));
								h2[i][j]->SetLineColor(kRed);
						}
				}
				utility::quickShow(name+"_sideCheck", h1, h2, ncent);
		}
		void checkBkg(TString name, int ncent=2){
				TString tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/"+name+"_JTCProj.root";
				TFile *f = TFile::Open(tmp);
				JTCSignalProducer *sp1[8][2];
				JTCSignalProducer *sp2[8][2];
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto cp1 = new mCanvasLoose("bkgCheck_"+name, "", nPt, 2*ncent);
				auto cp2 = new mCanvasLoose("bkgCheck_JS_"+name, "", nPt, 2*ncent);
				auto cp3 = new mCanvasLoose("bkgCheck_sideBand"+name, "", nPt, 2*ncent);
				for(int i=0; i<nPt ; ++i){
						for(int j=0; j<ncent ; ++j){
								if(ncent == 2 )  tmp = trk_tag[i]+", "+cent_label[j];
								else	tmp = trk_tag[i];
								sp1[i][j] = new JTCSignalProducer();
								sp2[i][j] = new JTCSignalProducer();
								sp1[i][j]->read1D(f, name+Form("_%d_%d", i, j));
								sp2[i][j]->read1D(f, name+Form("_pTweighted_%d_%d", i, j));
								cp1->CD(i+1, ncent-j);
								sp1[i][j]->drawBkgCheck();
								tx->DrawLatexNDC(0.02,0.98, tmp); 
								cp1->CD(i+1, 2*ncent-j);
								sp1[i][j]->drawBkgCheck(0);
								tx->DrawLatexNDC(0.02,0.98, tmp); 
								cp2->CD(i+1, ncent-j);
								sp2[i][j]->drawBkgCheck();
								tx->DrawLatexNDC(0.02,0.98, tmp); 
								cp2->CD(i+1, 2*ncent-j);
								sp2[i][j]->drawBkgCheck(0);
								tx->DrawLatexNDC(0.02,0.98, tmp); 

								cp3->CD(i+1, ncent-j);
								sp1[i][j]->drawSideBandCheck();
								tx->DrawLatexNDC(0.02,0.98, tmp); 
								cp3->CD(i+1, 2*ncent-j);
								sp2[i][j]->drawSideBandCheck();
								tx->DrawLatexNDC(0.02,0.98, tmp); 
						}
				}
				cp1->SaveAs(FigDumpPath+"bkgCheck_"+name+".gif");
				cp2->SaveAs(FigDumpPath+"bkgCheck_JS_"+name+".gif");
				cp3->SaveAs(FigDumpPath+"bkgCheck_sideBand_"+name+".gif");
		}
}


namespace inclusive_input{
		TFile* gengen_pb_sub0_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/incl_gen_gen_sub0_JTCSignal.root");
		void drawRatio_sub0(TString name1,  TString name2,TString name, bool doBinomial=0){
				TString tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/"+name1+"_JTCSignal.root";
				TFile *f1 = TFile::Open(tmp);
				tmp ="/Users/tabris/cmsProjects/iJTC/dataSet/correlation/"+name2+"_JTCSignal.root";
				TFile *f2 = TFile::Open(tmp);
				JTCSignalProducer *sp1[8][2];
				JTCSignalProducer *sp2[8][2];
				TFile* wf = TFile::Open("debug.root","recreate");
				for(int i=0; i<8; ++i){
						for(int j=0; j<2; ++j){
								sp1[i][j]=new JTCSignalProducer();
								sp2[i][j]=new JTCSignalProducer();
								sp1[i][j]->read(f1, name1+Form("_pTweighted_%d_%d", i,j));
								sp2[i][j]->read(f2, name2+Form("_pTweighted_%d_%d", i,j));
								sp1[i][j]->doDrIntegral(Form("sub0_1_%d_%d", i, j));
								sp2[i][j]->doDrIntegral(Form("sub0_2_%d_%d", i, j));
						}
				}
				utility::quickJSratio(name, sp1, sp2, doBinomial);
		}
		void sub0Ratio(TString fname){
				cout<<"pulling singal for "<<fname;
				TString trkbin [] = {"0.7", "1", "2", "3", "4", "8", "12", "16", "20", "999"};
				TString centbin [] = {"0", "30", "100"};
				JTCSignalProducer *sp1[8][2];
				JTCSignalProducer *sp2[8][2];
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
								sp2[i][j] = new JTCSignalProducer(raw_sig_pTweighted[i][j], mixing[i][j]);
								sp1[i][j]->sig=raw_sig[i][j];
								sp2[i][j]->sig=raw_sig[i][j];
						}
				}
				cout<<". Done"<<endl;;
				//clearInput();
		}
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

		void showSpectra(TString name, TString name1, TString name2, TFile *f1, TFile *f2){
				int ncent =2;
				if( !f1->IsOpen()) cout<<f1->GetName()<<" haven't been open"<<endl;
				getH2(name2, f2);
				//comparing the b-tagged jets with the inclusive jets
				// the second input has to come from the inclusive jets
				auto *cp = new doublePanelFig("c_"+name, "", 1, nCent );
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				TString tmp;
				TH1D* hr[nCent];
				TH1D* h[nCent][2];
				for(int j=0; j<ncent; ++j){
						tmp = name1+"_all_bjets_corrpT"+input_raw2D::cent_tag[j]+"_"+input_raw2D::cent_tag[j+1]+"_Pt120_Pt1000";
						cout<<tmp<<endl;
						h[j][0] =(TH1D*) f1->Get(tmp); h[j][0]->GetXaxis()->SetTitle("#Delta#eta"); h[j][0]->SetTitle("");
		//				h[j][0]->Scale(1.0/h[j][0]->Integral("width"));
						h[j][0]->SetMarkerColor(kBlue+2);
						h[j][0]->SetLineWidth(1);
						hr[j]=(TH1D*)h[j][0]->Clone(Form("hr_%d",j));
						h[j][0]->SetAxisRange(1e-6, h[j][0]->GetMaximum()*10, "Y");
						h[j][0]->SetAxisRange(100, 500 ,"X");
						//						h[j][0]->SetAxisRange(1e-13, 1e-4,"Y");
						cp->addHist(h[j][0], 1, ncent-j);
						cp->CD(1, 2-j, 0); gPad->SetLogy();
						h[j][1] = hjet[j]; h[j][1]->SetLineColor(kRed);  h[j][1]->SetMarkerColor(kRed);
		//				h[j][1]->Scale(1.0/h[j][1]->Integral("width"));
						cp->addHist(h[j][1], 1, ncent-j); 
						hr[j]->Divide(h[j][1]); hr[j]->GetYaxis()->SetNdivisions(505);
						hr[j]->SetAxisRange(100, 500, "X");  
						cp->addHist(hr[j], 1, ncent-j, 1);
						cp->CD(1, 2-j, 0); tl->DrawLine(-2.5, 0, 2.5, 0);
						tmp = "Jet spectra: "+cent_label[j];
						tx->DrawLatexNDC(0.15,0.87, tmp); 
						cp->CD(1, 2-j, 1); tl->DrawLine(-2.5, 1, 2.5, 1);
				}
				cp->SaveAs(FigDumpPath+name+"_jetSpectra.gif");
		}
}

