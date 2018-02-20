

void pull_one_from_another(TString name, TString name1, TString name2, float frac){
		TH2D* sig[8][2];
		TH2D* sub[8][2];

		getSig(name1, sig);
		getSig(name2, sub);
		float comp = 1-frac;
		TFile *wf = new TFile(dataDumpPath+name+"_JTCSignal.root","recreate");
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						sub[i][j]->Scale(comp);
						sig[i][j]->Add(sub[i][j],-1);
						sig[i][j]->Scale(1.0/frac);
						sig[i][j]->SetName("signal_"+name+Form("_%d_%d", i,j));
						sig[i][j]->Write();
				}
		}
		wf->Close();
}

void plot_closure_and_overlay(TString name, TString label1, TString label2, TH1D* heta1[8][2], TH1D* heta2[8][2],float x1, float x2, bool ispp, bool doDiff = 0){
		auto tx = new TLatex(); 
		if(ispp){
				auto tl = new TLegend(0.2,0.65,0.6, 0.83 ); tl->SetLineColor(0);
				tx->SetTextSize(0.058);
				auto cm = new doublePanelFig("c_"+name, "", 2, 3);
				auto line = new TLine(); line ->SetLineStyle(2);
				TH1D* hratio[8][2];
				int n=1;
				int j=0;
				for(int i=1; i<7; ++i){
						//								cm->CD(i%3, ceil(float(i)/3)+1);
						heta1[i][j]->SetAxisRange(x1, x2);
						hratio[i][j]=(TH1D*) heta1[i][j]->Clone(Form("ratio_%d_%d", i, j));
						if(!doDiff) hratio[i][j]->Divide(heta2[i][j]);
						else hratio[i][j]->Add(heta2[i][j],-1);
						hratio[i][j]->SetMarkerStyle(20);
						heta1[i][j]->SetMarkerStyle(20);
						heta2[i][j]->SetMarkerStyle(20);
						hratio[i][j]->SetMarkerSize(.5);
						heta1[i][j] ->SetMarkerSize(.5);
						heta2[i][j] ->SetMarkerSize(.5);
						adjustYrange(heta1[i][j]);
						hratio[i][j]->SetMarkerColor(kBlack);
						heta1[i][j] ->SetMarkerColor(kAzure+2);
						heta2[i][j] ->SetMarkerColor(kRed-4);
						heta1[i][j]->SetLineColor(kAzure+2);
						heta2[i][j]->SetLineColor(kRed-4);
						hratio[i][j]->SetLineColor(kBlack);
						hratio[i][j]->GetXaxis()->SetNdivisions(505);
						hratio[i][j]->GetYaxis()->SetNdivisions(505);
						if(!doDiff) hratio[i][j]->SetAxisRange(.2,1.8, "Y");
						cm->addHist(heta1[i][0], n/4+1, n%4);
						cm->addHist(heta2[i][0], n/4+1, n%4);
						cm->addHist(hratio[i][0], n/4+1, n%4, 1);
						//						cout<<n%4<<", "<<ceil(float(n)/3)<<endl;
						cm->CD(n/4+1, n%4,0 );
						line->DrawLine(-1, 0, 1,0 );
						tx->DrawLatexNDC(0.62, 0.85, track_label[i]);
						cm->CD(n/4+1, n%4,1 );
						if( doDiff) line->DrawLine(-1, 0, 1,0 );
						else line->DrawLine(-1, 1, 1,1 );
						n++;
						if(n==4 ) n++;
				}
				tl->AddEntry(heta1[1][0], label1);
				tl->AddEntry(heta2[1][0], label2);
				cm->CD(1,1, 0); tl->Draw();
				cm->CD(1,1, 1); 
				tx->SetTextSize(0.08);
				tx->DrawLatexNDC(0.2, 0.85, label1+" / "+label2);

				cm->SaveAs(name+".pdf");
		}
}

void dEtaClosure(TString name1, TString name2, TString lg1, TString lg2, bool ispp=1, bool doDiff=0){
		float trkbin[] = {1, 1, 1, 4, 4, 1}; 
		auto sp1 = batchRead1D(name1,0);
		auto sp2 = batchRead1D(name2,0);
		TH1D* heta1[8][2];
		TH1D* heta2[8][2];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<2; ++j){
						heta1[i][j]=sp1[i][j]->sig_deta;
						heta2[i][j]=sp2[i][j]->sig_deta;
						heta1[i][j]->Scale(trkbin[i]);
						heta2[i][j]->Scale(trkbin[i]);
				}
		}
		TString name = "closure_"+name1+"_"+name2;
		if(doDiff) name = "Diff_"+name1+"_"+name2;
		plot_closure_and_overlay(FigDumpPath+name, lg1, lg2, heta1, heta2,-1, 0.99, ispp, doDiff);
}

void overlay_dr_pp(TString name, TH1D* h1[8][2], TH1D* h2[8][2]){
		auto cpdr = new mCanvasLoose("overlay_"+name,  "", 2, 4); 
		auto tx = new TLatex();  tx->SetTextSize(.06); 
		for(int i=1; i<nPt; ++i){
				cout<<h2[i][0]->GetName()<<endl;
				h2[i][0]->SetLineColor(kGreen);
				h1[i][0]->SetTitle("");
				cpdr->drawHist(h1[i][0],ceil(float(i)/8), i);
				cpdr->drawHist(h2[i][0],ceil(float(i)/8), i);
				tx->DrawLatexNDC(0.45, 0.85, track_label[i]);
		}
		cpdr->SaveAs(FigDumpPath+"overlay_"+name+".gif");
}

void drRatioPlot(TString cname, TH1D* h1[8][2], TH1D* h2[8][2]){
		TH1* hratio[8][2];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<2; ++j){
						hratio[i][j]=(TH1*) h1[i][j]->Clone(Form("ratio_%d_%d",i,j));
						hratio[i][j]->Divide(h2[i][j]);
				}
		}
		int ncol=1;	
		auto cpdr = new mCanvasLoose("cpdr_"+cname,  "", 2, 3);
		TString tmp;
		auto line = new TLine(); line->SetLineStyle(2);
		auto tx = new TLatex();  tx->SetTextSize(.06);
		for(int i=1; i<nPt; ++i){
				for(int j=0; j<ncol; ++j){
						///cpdr->CD(i,ncol-j);
						cpdr->cd(i);
						if(ncol == 2 )  tmp = track_label[i]+", "+cent_label[j];
						else    tmp = track_label[i];
						hratio[i][j]->SetTitle("");
						hratio[i][j]->SetAxisRange(0.5, 1.5, "Y");
						TH1* h = hratio[i][j];
						h->SetMarkerStyle(20);
						h->SetMarkerSize(0.8);
						h->SetMarkerColor(kBlack);
						h->SetLineWidth(2);
						cpdr->drawHist(hratio[i][j], i/4+1, i%3+int((i%3)==0)*3);
						line->DrawLine(0,1, 1, 1);
						tx->DrawLatexNDC(0.2,0.88, tmp);
						//						tmp = "tagged b(remove incl)/tagged & true b";
						//				tx->DrawLatexNDC(0.1,0.98, tmp);
				}
		}
		cpdr->SaveAs(FigDumpPath + cname+".gif");
}

void purify_dr(TString name, TString name1, TString name2, TString name3, float frac, bool isNumber =1 ){
		TH1D* dr1[8][2];
		TH1D* dr2[8][2];
		TH1D* dr3[8][2];
		getDr(name1, dr1, isNumber);
		getDr(name2, dr2, isNumber);
		getDr(name3, dr3, isNumber);
		float comp = 1-frac;
		auto cc = new mCanvasLoose("cc", "", 7, 1);
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<2; ++j){
						//				cout<<dr1[i][j]->GetName()<<endl;
						dr2[i][j]->Scale(comp);
						dr1[i][j]->Add(dr2[i][j],-1);
						dr1[i][j]->Scale(1.0/frac);
						//						for(int n= 1; n<dr2[i][j]->GetNbinsX()+1; ++n){
						//								dr2[i][j]->SetBinContent(n , 1);
						//						}
						//						dr1[i][j]->Write();
				}
		}
		TFile *wf = new TFile(dataDumpPath+name+"_JTCProj.root", "recreate");
		TH1D* dr4[8][2];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<2; ++j){
						dr4[i][j]=(TH1D*) dr1[i][j]->Clone("dr_"+name+Form("_%d_%d", i, j));
						dr4[i][j]->Write();
				}
		}
		overlay_dr_pp("check", dr1, dr3);
		drRatioPlot(name, dr1, dr3);
}

void drRatio(TString cname, TString name1, TString name2){
		TH1D* h1[8][2];
		getDr(name1, h1);
		TH1D* h2[8][2];
		getDr(name2, h2);
		drRatioPlot(cname, h1, h2);
		overlay_dr_pp(cname, h1, h2);
}



