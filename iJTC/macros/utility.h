

TString dataDumpPath = "/Users/tabris/cmsProjects/iJTC/dataSet/correlation/";
TString FigDumpPath = "/Users/tabris/cmsProjects/iJTC/macros/figures/";

TString trk_tag[] = {"0.7 < p_{T}^{track} < 1 GeV","1 < p_{T}^{track} < 2 GeV",
		"2 < p_{T}^{track} < 3 GeV", "3 < p_{T}^{track} < 4 GeV",
		"4 < p_{T}^{track} < 8 GeV", "8 < p_{T}^{track} < 12 GeV", 
		"12 < p_{T}^{track} < 16 GeV", "p_{T}^{track} > 16 Gev"};
TString cent_label[]={"Cent. 0-30%", "Cent. 30-100%"};

namespace utility{
		void quickLabels_loosCanvase(mCanvasLoose *c){
				auto tx = new TLatex();  
				tx->SetTextSize(.08);
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								TString tmp = trk_tag[i]+", "+cent_label[j];
								c->CD(i+1, 2-j);
								tx->DrawLatexNDC(0.02,0.98, tmp); 
						}
				}
		}

		void quickRatio(TString name, TH1* h1[8][2], TH1* h2[8][2], float x1, float x2, bool ispp=1){
				int ncent = ispp ? 1: 2; int drop =1;
				auto *cp = new doublePanelFig("c_"+name, "", nPt-drop,ncent );
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				TH1D* hr[nPt][ncent];
				TString tmp;
				for(int i=drop; i<nPt; ++i){
						int ic=0;
						for(int j=0; j<ncent; ++j){
								h1[i][j]->SetMarkerColor(kBlue+2);
								h1[i][j]->SetLineColor(kBlue+2);
								h1[i][j]->SetLineWidth(2);
								h1[i][j]->SetAxisRange(x1, x2, "X");
								hr[i][j]=(TH1D*)h1[i][j]->Clone(Form("hr_%d_%d",i,j));
								cp->addHist(h1[i][j], ic+1, ncent-j);
								h2[i][j]->SetLineColor(kRed);  
								h2[i][j]->SetMarkerColor(kRed);
								h2[i][j]->SetLineWidth(2);
								cp->addHist(h2[i][j], ic+1, ncent-j); 
								hr[i][j]->Divide(h2[i][j]); hr[i][j]->GetYaxis()->SetNdivisions(505);
								hr[i][j]->SetAxisRange(0, 2, "Y");  cp->addHist(hr[i][j], ic+1, ncent-j, 1);
								cp->CD(ic+1, ncent-j, 0); tl->DrawLine(x1, 0, x2, 0);
								if(ispp) tmp = trk_tag[i];
								else	tmp = trk_tag[i]+", "+cent_label[j];
								tx->DrawLatexNDC(0.15,0.87, tmp); 
								cp->CD(ic+1, ncent-j, 1); tl->DrawLine(x1, 1,x2, 1);
								ic++;
						}
				}
				cp->SaveAs(FigDumpPath+"ratio_"+name+".gif");
		}

		void quickJSOverlay(TString name1, TString name2, TFile *f1, TFile *f2, float x1, float x2){
				cout<<"drawing the overlay for "<<name1<<" and "<<name2<<endl;	
				JTCSignalProducer *sp1 [nPt][nCent];
				JTCSignalProducer *sp2 [nPt][nCent];
				auto *cp = new doublePanelFig("c_"+name1+"_"+name2, "", nPt,nCent );
				TH1D* hr[nPt][nCent];
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				TString tmp;
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								cout<<i<<", "<<j<<endl;
								sp1[i][j]= new JTCSignalProducer();
								sp2[i][j]= new JTCSignalProducer();
								sp1[i][j]->read(f1, name1+Form("_%d_%d", i,j));
								sp2[i][j]->read(f2, name2+Form("_%d_%d", i,j));
								sp1[i][j]->doDrIntegral(name1+Form("_%d_%d", i, j));
								sp2[i][j]->doDrIntegral(name2+Form("_%d_%d", i, j));
								TH1* h = sp1[i][j]->dr_integral; h->SetTitle("");
								h->SetMarkerColor(kBlue+2);
								h->SetLineColor(kBlue+2);
								h->SetLineWidth(2);
								h->SetAxisRange(x1, x2, "X");
								hr[i][j]=(TH1D*)h->Clone(Form("hr_%d_%d",i,j));
								cp->addHist(h, i+1, 2-j);
								h=sp2[i][j]->dr_integral; h->SetLineColor(kRed);  h->SetMarkerColor(kRed);
								h->SetLineWidth(2);
								cp->addHist(h, i+1, 2-j); 
								hr[i][j]->Divide(h); hr[i][j]->GetYaxis()->SetNdivisions(505);
								hr[i][j]->SetAxisRange(0, 2, "Y");  cp->addHist(hr[i][j], i+1, 2-j, 1);
								cp->CD(i+1, 2-j, 0); tl->DrawLine(x1, 0, x2, 0);
								tmp = trk_tag[i]+", "+cent_label[j];
								tx->DrawLatexNDC(0.15,0.87, tmp); 
								cp->CD(i+1, 2-j, 1); tl->DrawLine(x1, 1,x2, 1);
								cp->draw95Area(i+1,2-j, x1, x2);
						}
				}
				cp->SaveAs(FigDumpPath+"quickLook_"+name1+"_"+name2+"_overlay.gif");
		}
		void quickJSratio(TString name, JTCSignalProducer* sp1[8][2], JTCSignalProducer* sp2[8][2], bool doBinomial = 0, bool ispp=0){
				int ncent = ispp ? 1:2;
				cout<<ncent<<endl;
				float x1=0, x2=0.99;
				auto *cp = new doublePanelFig("c_"+name, "", nPt,ncent );
				TH1D* hr[nPt][nCent];
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				TString tmp;
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<ncent; ++j){
								TH1* h = sp1[i][j]->dr_integral; h->SetTitle("");
								h->SetMarkerColor(kBlue+2);
								h->SetLineColor(kBlue+2);
								h->SetLineWidth(2);
								h->SetAxisRange(x1, x2, "X");
								h->SetAxisRange(0, 1.4*h->GetMaximum(), "Y");
								hr[i][j]=(TH1D*)h->Clone(Form("hr_%d_%d",i,j));
								cp->addHist(h, i+1, ncent-j);
								h=sp2[i][j]->dr_integral; h->SetLineColor(kRed);  h->SetMarkerColor(kRed);
								h->SetLineWidth(2);
								cp->addHist(h, i+1, ncent-j); 
								if( doBinomial) hr[i][j]->Divide(hr[i][j], h, 1,1, "B");
								else hr[i][j]->Divide(h); 
								hr[i][j]->GetYaxis()->SetNdivisions(505);
								hr[i][j]->SetAxisRange(0.5, 1.5, "Y");  cp->addHist(hr[i][j], i+1, ncent-j, 1);
								cp->CD(i+1, ncent-j, 0); tl->DrawLine(x1, 0, x2, 0);
								if(ispp) tmp = trk_tag[i];
								else tmp = trk_tag[i]+", "+cent_label[j];
								tx->DrawLatexNDC(0.15,0.87, tmp); 
								cp->CD(i+1, ncent-j, 1); tl->DrawLine(x1, 1,x2, 1);
								//								cp->draw95Area(i+1,2-j, x1, x2);
						}
				}
				cp->SaveAs(FigDumpPath+"quickLook_"+name+"_ratio.gif");
		}

		void quickShow(TString name, TH1* h1[8][2], TH1* h2[8][2], int ncent=2){
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				auto c = new mCanvasLoose("quick_"+name, "", nPt, ncent);
				TString tmp;
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<ncent; ++j){
								if( ncent ==2 )	tmp = trk_tag[i]+", "+cent_label[j];
								else tmp = trk_tag[i];
								c->CD(i+1, 2-j);
								h1[i][j]->Draw();
								h2[i][j]->Draw("same");
								tx->DrawLatexNDC(0.02,0.98, tmp); 
						}
				}
				c->SaveAs(FigDumpPath+"quickLook_"+name+".gif");
		}
/*
		void quickTwoPanelFig(TString name, TH1* h1[8][2], TH1* h2[8][2], TH1* hh[8][2], float y1, float y2){
				float x1=0, x2=0.99;
				auto *cp = new doublePanelFig("c_"+name, "", nPt,nCent );
				auto tx = new TLatex();  tx->SetTextSize(.08);
				auto tl = new TLine();   tl->SetLineStyle(2); tl->SetLineWidth(2);
				TString tmp;
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								TH1* h = sp1[i][j]->dr_integral; h->SetTitle("");
								h->SetMarkerColor(kBlue+2);
								h->SetLineColor(kBlue+2);
								h->SetLineWidth(2);
								h->SetAxisRange(x1, x2, "X");
								h->SetAxisRange(0, 1.4*h->GetMaximum(), "Y");
								hr[i][j]=(TH1D*)h->Clone(Form("hr_%d_%d",i,j));
								cp->addHist(h, i+1, 2-j);
								h=sp2[i][j]->dr_integral; h->SetLineColor(kRed);  h->SetMarkerColor(kRed);
								h->SetLineWidth(2);
								cp->addHist(h, i+1, 2-j); 
								if( doBinomial) hr[i][j]->Divide(hr[i][j], h, 1,1, "B");
								else hr[i][j]->Divide(h); 
								hr[i][j]->GetYaxis()->SetNdivisions(505);
								for(int i=0; i<nPt; ++i){
										for(int j=0; j<nCent; ++j){

										}
								}
						}
				}
		}
		*/
}
