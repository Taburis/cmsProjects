
void pull_one_from_another(TString name, TString name1, TString name2, float frac){
		TH2D* sig[8][2];
		TH2D* sub[8][2];

		getSig(name1, sig);
		getSig(name2, sub);
		float comp = 1-frac;
		TFile *wf = new TFile(dataDumpPath+name+"_JTCSignal.root","recreate");
		for(int i=0; i<nPt; ++i){
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
				auto tl = new TLegend(0.15,0.6,0.51, 0.8 ); tl->SetLineColor(0);
				tx->SetTextSize(0.058);
				auto cm = new doublePanelFig("c_"+name, "", 2, 3);
				auto line = new TLine(); line ->SetLineStyle(2);
				TH1D* hratio[8][2];
				int n=1;
				int j=0;
				for(int i=1; i<7; ++i){
						//								cm->CD(i%3, ceil(float(i)/3)+1);
						heta1[i][j]->SetAxisRange(x1, x2);
//						symmetrize(heta1[i][0], -1, 1);
//						symmetrize(heta2[i][0], -1, 1);
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
						if(!doDiff) hratio[i][j]->SetAxisRange(.6,1.4, "Y");
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
				if( doDiff)
						tx->DrawLatexNDC(0.2, 0.85, label1+" - "+label2);
				else
						tx->DrawLatexNDC(0.2, 0.85, label1+" / "+label2);

				cm->SaveAs(name+".pdf");
		}
}

void dEtaRawClosure(TString name1, TString name2, TString lg1, TString lg2, bool ispp=1, bool doDiff=0, bool isNumber =0 ){
		TH2D* h1[8][2];
		TH2D* h2[8][2];
		getRawSig(name1, h1, isNumber);
		getRawSig(name2, h2, isNumber);

		TH1D* heta1[8][2];
		TH1D* heta2[8][2];
		//cout<<h1[1][0]->GetName()<<endl;
		float x1= h1[1][0]->GetYaxis()->FindBin(-1);
		float x2= h1[1][0]->GetYaxis()->FindBin(0.99);
		float width= h1[1][0]->GetYaxis()->GetBinWidth(1);
		for(int i=1; i<nPt; ++i){
				heta1[i][0]=(TH1D*) h1[i][0]->ProjectionX(Form("deta1_%d",i), x1, x2, "E");
				heta2[i][0]=(TH1D*) h2[i][0]->ProjectionX(Form("deta2_%d",i), x1, x2, "E");
				heta1[i][0]->Rebin(1);
				heta2[i][0]->Rebin(1);
				heta1[i][0]->Scale(.5*width);
				heta2[i][0]->Scale(.5*width);
				heta1[i][0]->GetXaxis()->SetTitle("d#eta");
		}
		TString name = "raw_dEta_closure_"+name1+name2;
		plot_closure_and_overlay(FigDumpPath+name, lg1, lg2, heta1, heta2,-1, 0.99, ispp, doDiff);
}

void dEtaClosure(TString name1, TString name2, TString lg1, TString lg2, bool ispp=1, bool doDiff=0, bool isNumber =0 ){
		TH2D* h1[8][2];
		TH2D* h2[8][2];
		getSig(name1, h1, isNumber);
		getSig(name2, h2, isNumber);

		float trkbin[] = {1, 1, 1, 4, 4, 1}; 
		TH1D* heta1[8][2];
		TH1D* heta2[8][2];
		//cout<<h1[1][0]->GetName()<<endl;
		float x1= h1[1][0]->GetYaxis()->FindBin(-1);
		float x2= h1[1][0]->GetYaxis()->FindBin(0.99);
		float width= h1[1][0]->GetYaxis()->GetBinWidth(1);
		for(int i=1; i<nPt; ++i){
				heta1[i][0]=(TH1D*) h1[i][0]->ProjectionX(Form("deta1_%d",i), x1, x2, "B");
				heta2[i][0]=(TH1D*) h2[i][0]->ProjectionX(Form("deta2_%d",i), x1, x2, "B");
				heta1[i][0]->Rebin(4);
				heta2[i][0]->Rebin(4);
				heta1[i][0]->Scale(.25*width);
				heta2[i][0]->Scale(.25*width);
				if(isNumber) heta1[i][0]->Scale(1.0/trkbin[i]);
				if(isNumber) heta2[i][0]->Scale(1.0/trkbin[i]);
		}
		TString name = "dEta_closure_"+name1+"_"+name2;
		if(doDiff) name = "dEta_diff_"+name1+"_"+name2;
		plot_closure_and_overlay(FigDumpPath+name, lg1, lg2, heta1, heta2,-1, 0.99, ispp, doDiff);
}

void dPhiClosure(TString name1, TString name2, TString lg1, TString lg2, bool ispp=1, bool doDiff=0, bool isNumber =0, float eta1=-1, float eta2= 0.99){
		TH2D* h1[8][2];
		TH2D* h2[8][2];
		getSig(name1, h1, isNumber);
		getSig(name2, h2, isNumber);

		float trkbin[] = {1, 1, 1, 4, 4, 1}; 
		TH1D* hphi1[8][2];
		TH1D* hphi2[8][2];
		//cout<<h1[1][0]->GetName()<<endl;
		float x1= h1[1][0]->GetXaxis()->FindBin(eta1);
		float x2= h1[1][0]->GetXaxis()->FindBin(eta2);
		float width= h1[1][0]->GetXaxis()->GetBinWidth(1);
		for(int i=1; i<nPt; ++i){
				hphi1[i][0]=(TH1D*) h1[i][0]->ProjectionY(Form("dphi1_%d",i), x1, x2, "B");
				hphi2[i][0]=(TH1D*) h2[i][0]->ProjectionY(Form("dphi2_%d",i), x1, x2, "B");
				hphi1[i][0]->Rebin(2);
				hphi2[i][0]->Rebin(2);
				hphi1[i][0]->Scale(.5*width);
				hphi2[i][0]->Scale(.5*width);
				if(isNumber) hphi1[i][0]->Scale(1.0/trkbin[i]);
				if(isNumber) hphi2[i][0]->Scale(1.0/trkbin[i]);
		}
		TString name = "dPhi_closure_"+name1+"_"+name2;
		if(doDiff) name = "dPhi_diff_"+name1+"_"+name2;
		plot_closure_and_overlay(FigDumpPath+name, lg1, lg2, hphi1, hphi2,-1, 0.99, ispp, doDiff);
}

void overlay_dr_pp(TString name, TH1D* h1[8][2], TH1D* h2[8][2]){
		auto cpdr = new mCanvasLoose("overlay_"+name,  "", 2, 3); 
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

void drRatioPlot(TString cname, TH1D* h1[8][2], TH1D* h2[8][2], TString opt = ""){
		TH1* hratio[8][2];
		for(int i=1; i<nPt; ++i){
				for(int j=0; j<1; ++j){
						hratio[i][j]=(TH1*) h1[i][j]->Clone(Form("ratio_%d_%d",i,j));
						hratio[i][j]->Divide(hratio[i][j], h2[i][j], 1, 1, opt);
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
						h->GetXaxis()->SetNdivisions(505);
						h->SetMarkerColor(kBlack);
						//h->SetLineWidth(2);
						cpdr->drawHist(hratio[i][j], i/4+1, i%3+int((i%3)==0)*3);
						line->DrawLine(0,1, 1, 1);
						tx->DrawLatexNDC(0.2,0.93, tmp);
						//						tmp = "tagged b(remove incl)/tagged & true b";
						//				tx->DrawLatexNDC(0.1,0.98, tmp);
				}
		}
		cpdr->SaveAs(FigDumpPath + cname+".pdf");
}

void purify_dr(TString name, TString name1, TString name2, TString name3, float frac, bool isNumber =0 ){
		TH1D* dr1[8][2];
		TH1D* dr2[8][2];
		TH1D* dr3[8][2];
		getDr(name3, dr3, isNumber);

		TH2D* h1[8][2];
		TH2D* h2[8][2];
		getSig(name1, h1, isNumber);
		getSig(name2, h2, isNumber);
		float comp = 1-frac;
		auto cc = new mCanvasLoose("cc", "", 7, 1);
		for(int i=1; i<nPt; ++i){
				for(int j=0; j<2; ++j){
						//				cout<<dr1[i][j]->GetName()<<endl;
						h2[i][j]->Scale(comp);
						h1[i][j]->Add(h2[i][j],-1);
						h1[i][j]->Scale(1.0/frac);
						//						for(int n= 1; n<dr2[i][j]->GetNbinsX()+1; ++n){
						//								dr2[i][j]->SetBinContent(n , 1);
						//						}
						//						dr1[i][j]->Write();
				}
		}
		TFile *wf = new TFile(dataDumpPath+name+"_JTCSignal.root", "recreate");
		TH1D* dr4[8][2];
		JTCSignalProducer* sp[8][2];
		for(int i=1; i<nPt; ++i){
				for(int j=0; j<2; ++j){
				//		cout<<"here i= "<<i<<", j= "<<j<<endl;
						sp[i][j]=new JTCSignalProducer();
						if(!isNumber ) {
								sp[i][j]->sig=(TH2D*) h1[i][j]->Clone("signal_"+name+Form("_pTweighted_%d_%d", i, j));
								sp[i][j]->doDrIntegral(name+Form("_pTweighted_%d_%d", i, j));
				  		}
						else {
								sp[i][j]->sig=(TH2D*) h1[i][j]->Clone("signal_"+name+Form("_%d_%d", i, j));
								sp[i][j]->doDrIntegral(name+Form("_%d_%d", i, j));
						}
						dr4[i][j]=(TH1D*) sp[i][j]->dr_integral;
						sp[i][j]->sig->Write();
				}
		}
		overlay_dr_pp("check", dr4, dr3);
		drRatioPlot(name, dr4, dr3, "B");
		wf->Close();
}

void drRatio(TString cname, TString name1, TString name2, TString opt= ""){
		TH1D* h1[8][2];
		getDr(name1, h1, 0);
		TH1D* h2[8][2];
		getDr(name2, h2, 0);
		drRatioPlot(cname, h1, h2, opt);
		overlay_dr_pp(cname, h1, h2 );
}

void drRawRatio(TString cname, TString name1, TString name2){
		int isNumber = 0;
		TH2D* h1[8][2];
		TH2D* h2[8][2];
		getRawSig(name1, h1, isNumber);
		getRawSig(name2, h2, isNumber);
		JTCSignalProducer* sp[8][2];
		TH1D* dr1[8][2], *dr2[8][2];
		for(int i=1; i<nPt; ++i){
				sp[i][0]=new JTCSignalProducer();
				sp[i][1]=new JTCSignalProducer();
				sp[i][0]->sig=h1[i][0];
				sp[i][1]->sig=h2[i][0];
				sp[i][0]->doDrIntegral(Form("js1_%d",i));
				sp[i][1]->doDrIntegral(Form("js2_%d",i));
				dr1[i][0]=sp[i][0]->dr_integral;
				dr2[i][0]=sp[i][1]->dr_integral;
				cout<<dr2[i][0]->GetName()<<endl;
		}
		cout<<"plotting"<<endl;
		drRatioPlot(cname, dr1, dr2);
}

void draw_pp_SideCheck(TString name, JTCSignalProducer *sp[8][2], TH1D* ref[8]){
		/*
		   TFile *f3 = TFile::Open(dataDumpPath+"gen_gen_PYTHIA_tagged_trueB_JTCSignal.root");
		   TFile *f4 = TFile::Open(dataDumpPath+"gen_gen_PYTHIA_trueB_JTCSignal.root");
		   JTCSignalProducer *sp3[8][2];
		   JTCSignalProducer *sp4[8][2];
		   */
		auto cm = new mCanvasLoose("pp_side_check_"+name, "", 2,3);
		auto tx = new TLatex; tx->SetTextSize(0.06);
		auto line = new TLine(); line->SetLineStyle(2);
		for(int i=1; i<nPt ; ++i){
				tmp = track_label[i];
				sp[i][0]->getAllProj(name+Form("_%d_%d",i,0), 1);
				sp[i][0]->sig_deta->SetTitle("");
				sp[i][0]->sig_deta->SetAxisRange(-2.5, 2.49, "X");
				cm->drawHist(ref[i], i/4+1, i%3+int((i%3)==0)*3 );
				cm->drawHist(sp[i][0]->sig_deta, i/4+1, i%3+int((i%3)==0)*3 );
				sp[i][0]->side_deta->SetLineColor(kOrange+7);
				//cm->drawHist(sp[i][0]->side_deta, i/4+1, i%3+int((i%3)==0)*3 );

				cm->CD(i/4+1, i%3+int((i%3)==0)*3 );
				tx->DrawLatexNDC(0.2, 0.93, tmp);		
				line->DrawLine(-2.5,0, 2.5, 0);
		}
		/*
		   TH1D* diff[8][2];
		   for(int i=1; i<nPt ; ++i){
		   sp3[i][0] = new JTCSignalProducer();
		   sp4[i][0] = new JTCSignalProducer();
		   sp3[i][0]->read(f3, Form("gen_gen_PYTHIA_tagged_trueB_pTweighted_%d_%d", i, 0));
		   sp4[i][0]->read(f4, Form("gen_gen_PYTHIA_trueB_pTweighted_%d_%d", i, 0));
		   sp3[i][0]->getAllProj("gen_gen_PYTHIA_tagged_trueB");
		   sp4[i][0]->getAllProj("gen_gen_PYTHIA_trueB");
		   diff[i][0]=(TH1D*) sp3[i][0]->sig_deta->Clone(Form("diff_%d",i));
		   diff[i][0]->Add(sp4[i][0]->sig_deta, -1);
		   diff[i][0]->SetLineColor(kRed);
		   cm->drawHist(diff[i][0], i/4+1, i%3+int((i%3)==0)*3 );
		   }
		   */
		cout<<"drawing"<<endl;
		cm->SaveAs(FigDumpPath+name+".pdf");
}

void pullCorrection(TString name,TString fname1, TString fname2, float sidemin=1.5, float sidemax=2.5, bool doSeagull=0){

		TFile *f1 = TFile::Open(dataDumpPath+fname1+"_JTCSignal.root");
		TFile *f2 = TFile::Open(dataDumpPath+fname2+"_JTCSignal.root");
		JTCSignalProducer *sp1[8][2];
		JTCSignalProducer *sp2[8][2];
		JTCSignalProducer *sp3[8][2];
		TH1D* ref[8];
		TFile * wf = new TFile(dataDumpPath+name+"_JTCSignal.root","recreate");
		for(int i=1; i<nPt ; ++i){
				//				cout<<".";
				for(int j=0; j<1 ; ++j){
						sp1[i][j] = new JTCSignalProducer();
						sp2[i][j] = new JTCSignalProducer();
						sp3[i][j] = new JTCSignalProducer();
						sp1[i][j]->read(f1, fname1+Form("_pTweighted_%d_%d", i, j));
						sp2[i][j]->read(f2, fname2+Form("_pTweighted_%d_%d", i, j));
						//						sp1[i][j]->getAllProj(fname1+Form("_pTweighted_%d_%d", i, j));
						//						sp2[i][j]->getAllProj(fname2+Form("_pTweighted_%d_%d", i, j));
						ref[i]=(TH1D*) sp1[i][0]->sig_deta->Clone(Form("ref_%d",i));
						ref[i]->Add(sp2[i][0]->sig_deta, -1);
						ref[i]->SetLineColor(kRed);
						cout<<ref[i]->GetName()<<endl;
						sp3[i][j]->sig_step2 = (TH2D*)sp1[i][j]->sig_step2->Clone("sig_mix_corrected_"+name+Form("_%d_%d",i,j));
						sp3[i][j]->sig_step2->Add(sp2[i][j]->sig_step2, -1);
						sp3[i][j]->sig=(TH2D*) sp3[i][j]->sig_step2->Clone("signal_"+name+Form("_%d_%d",i,j));
						//						cout<<sp3[i][j]->sig->GetName()<<endl;
						sp3[i][j]->doBkgSubtraction(sidemin, sidemax);
						sp3[i][j]->sig_step2->Write();
						sp3[i][j]->sig->Write();

						delete sp1[i][j];
						delete sp2[i][j];
				}
		}
		draw_pp_SideCheck("fig_"+name, sp3, ref);
}

void addSignal(TString name, TString name1, TString name2){
		TH2D *hsig[8][2], *hsub[8][2];
		//		getSig(name1, , 0);
}

void drawSideCheck_pp(TString name, TH1D* sig_deta[8][2], TH1D* side_deta[8][2]){
		auto cm = new mCanvasLoose("pp_side_check_"+name, "", 2,3);
		auto tx = new TLatex; tx->SetTextSize(0.06);
		auto line = new TLine(); line->SetLineStyle(2);
		for(int i=1; i<nPt ; ++i){
				tmp = track_label[i];
				sig_deta[i][0]->SetTitle("");
				sig_deta[i][0]->SetAxisRange(-2.5, 2.49, "X");
				sig_deta[i][0]->GetXaxis()->SetTitle("#Delta#eta");
				//	cm->drawHist(ref[i], i/4+1, i%3+int((i%3)==0)*3 );
				cm->drawHist(sig_deta[i][0], i/4+1, i%3+int((i%3)==0)*3 );
				side_deta[i][0]->SetLineColor(kOrange+7);
				//			cm->drawHist(side_deta[i][0], i/4+1, i%3+int((i%3)==0)*3 );

				cm->CD(i/4+1, i%3+int((i%3)==0)*3 );
				tx->DrawLatexNDC(0.2, 0.93, tmp);		
				line->DrawLine(-2.5,0, 2.5, 0);
		}
		cm->SaveAs(FigDumpPath+name+".pdf");
}

void taggingBiasCorrection(){
		TH2D *hrawsig[8][2];
		TH2D *hrawsig_pTweighted[8][2];
		TH2D *hmix[8][2];
		TH2D *hmix2[8][2];
		get2DInput_GenGen_PYTHIA_taggedAndTrueB();	
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						hrawsig[i][j]=(TH2D*)raw_sig[i][j]->Clone(Form("taggingBias_%d_%d",i,j));
						hrawsig_pTweighted[i][j]=(TH2D*)raw_sig_pTweighted[i][j]->Clone(Form("taggingBias_pTweighted_%d_%d",i,j));
				}
		}
		get2DInput_GenGen_PYTHIA_trueB();	
		JTCSignalProducer * sp1[8][2];
		JTCSignalProducer * sp2[8][2];
		TFile* wf = new TFile(dataDumpPath + "taggingBias_JTCSignal.root","recreate");
		TH1D* side_deta[8][2];
		TH1D* sig_deta[8][2];
		TH1D* side_deta2[8][2];
		TH1D* sig_deta2[8][2];
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						hrawsig[i][j]->Add(raw_sig[i][j],-1);
						hrawsig_pTweighted[i][j]->Add(raw_sig_pTweighted[i][j],-1);
						//					hmix2[i][j]=mixing[i][j];
						sp1[i][j]=new JTCSignalProducer();
						//					sp1[i][j]=sig_step2=(TH2D*) hrawsig[i][j]->Clone(Form("sig_mix_corrected_taggingBias_%d_%d",i,j));
						sp1[i][j]->sig=hrawsig[i][j];
						sp1[i][j]->makeInvariant(sp1[i][j]->sig);
						sp1[i][j]->doBkgSubtraction(Form("taggingBias_%d_%d",i,j));
						sp1[i][j]->sig->Write();
						side_deta[i][j]=sp1[i][j]->getSignal_phiSideBand(Form("taggingBias_%d_%d",i,j));
						sig_deta[i][j]=sp1[i][j]->getSignal_dEta(Form("taggingBias_%d_%d",i,j));
						//hrawsig[i][j]->Write();
						sp2[i][j]=new JTCSignalProducer();
						sp2[i][j]->sig=hrawsig_pTweighted[i][j];
						sp2[i][j]->makeInvariant(sp2[i][j]->sig);
						sp2[i][j]->doBkgSubtraction(Form("taggingBias_%d_%d",i,j));
						sp2[i][j]->sig->Write();
						side_deta2[i][j]=sp2[i][j]->getSignal_phiSideBand(Form("taggingBias_%d_%d",i,j));
						sig_deta2[i][j]=sp2[i][j]->getSignal_dEta(Form("taggingBias_%d_%d",i,j));
				}
		}
		drawSideCheck_pp("tagging_bias_correction", sig_deta, side_deta);
		drawSideCheck_pp("tagging_bias_JS_correction", sig_deta2, side_deta2);
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						delete sp1[i][j];
						delete sp2[i][j];
				}
		}
}

void addBiasCorrection(){
		bool isNumber =0;
		TH1D *dr1[8][2];
		TH1D *dr_ref[8][2];
		getDr("purified_gen_gen_PYTHIA_taggedB", dr1, isNumber);
		getDr("gen_gen_PYTHIA_trueB", dr_ref, isNumber);
		JTCSignalProducer* sp[8];
		TFile *f = TFile::Open(dataDumpPath+"taggingBias_JTCSignal.root");
		for(int i=0; i<nPt; ++i){
				sp[i]= new JTCSignalProducer();
				sp[i]->sig = (TH2D*) f->Get(Form("taggingBias_pTweighted_%d_0",i));
				sp[i]->doDrIntegral(Form("taggingBias_pTweighted_%d_0",i));
				dr1[i][0]->Add(sp[i]->dr_integral,-1);
		}
		drRatioPlot("closure_taggingBiasCorrected", dr1, dr_ref, "B");
}

void slicedEtaRawClosure(TString name1, TString name2, TString lg1, TString lg2, bool doDiff=0, bool isNumber =0 ){
		TH2D* h1[8][2];
		TH2D* h2[8][2];
		getRawSig(name1, h1, isNumber);
		getRawSig(name2, h2, isNumber);

		TH1D* heta1[8][2];
		TH1D* heta2[8][2];
		int nn = 20;
		float v1 = -1, deta= 2.0/nn;
		//cout<<h1[1][0]->GetName()<<endl;
		for(int j=0; j<nn; j++){
				float x1= h1[1][0]->GetYaxis()->FindBin(v1);
				float x2= h1[1][0]->GetYaxis()->FindBin(v1+deta);
				v1+=deta;
				float width= h1[1][0]->GetYaxis()->GetBinWidth(1);
				for(int i=1; i<nPt; ++i){
						heta1[i][0]=(TH1D*) h1[i][0]->ProjectionX(Form("deta1_%d",i), x1, x2, "B");
						heta2[i][0]=(TH1D*) h2[i][0]->ProjectionX(Form("deta2_%d",i), x1, x2, "B");
						heta1[i][0]->Rebin(4);
						heta2[i][0]->Rebin(4);
						heta1[i][0]->Scale(.25*width);
						heta2[i][0]->Scale(.25*width);
				}
				TString name = Form("raw_dEtaSlice_%d",j)+name1+name2;
				plot_closure_and_overlay(FigDumpPath+name, lg1, lg2, heta1, heta2,-1, 0.99, 1, doDiff);
		}
}


void dPhiRawClosure(TString name1, TString name2, TString lg1, TString lg2, bool ispp=1, bool doDiff=0, bool isNumber =0, float eta1=-1, float eta2= 0.99){
		TH2D* h1[8][2];
		TH2D* h2[8][2];
		getRawSig(name1, h1, isNumber);
		getRawSig(name2, h2, isNumber);

/* =======================quick test =====================*/ 
/*=============================================*/

		TH1D* hphi1[8][2];
		TH1D* hphi2[8][2];
		//cout<<h1[1][0]->GetName()<<endl;
		float x1= h1[1][0]->GetXaxis()->FindBin(eta1);
		float x2= h1[1][0]->GetXaxis()->FindBin(eta2);
		float width= h1[1][0]->GetYaxis()->GetBinWidth(1);
		for(int i=1; i<nPt; ++i){
				//hphi1[i][0]=(TH1D*) h1[i][0]->ProjectionY(Form("dphi1_%d",i), 1, -1, "B");
				//hphi2[i][0]=(TH1D*) h2[i][0]->ProjectionY(Form("dphi2_%d",i), 1, -1, "B");
				hphi1[i][0]=(TH1D*) h1[i][0]->ProjectionY(Form("dphi1_%d",i), x1, x2, "B");
				hphi2[i][0]=(TH1D*) h2[i][0]->ProjectionY(Form("dphi2_%d",i), x1, x2, "B");
				cout<<h1[i][0]->GetName()<<endl;
				cout<<h2[i][0]->GetName()<<endl;
				hphi1[i][0]->Rebin(2);
				hphi2[i][0]->Rebin(2);
				hphi1[i][0]->Scale(.5*width);
				hphi2[i][0]->Scale(.5*width);
				hphi1[i][0]->GetXaxis()->SetTitle("#Delta#phi");
		}
		TString name = "raw_signal_closure_dphi_"+name1+name2;
		plot_closure_and_overlay(FigDumpPath+name, lg1, lg2, hphi1, hphi2,-1, 0.99, ispp, doDiff);
}

Double_t func(Double_t x,Double_t y,Double_t *par)
{
		Double_t value=( par[0]*TMath::Exp(-pow((x-par[1])/par[2], 2)-pow((y-par[3])/par[4],2)));
		//		cout<<" value = "<<value<<endl;
		//		cout<<" x = "<<x<<endl;
		//		cout<<" y = "<<y<<endl;
		return value;
}

#include "fitter.h"
void fit2D(bool isNumber =0){
		TFile *f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/pfJetCorrelation/taggingBias_JTCSignal.root");
		TH2D* h2[8][2];
		TH2D* hfit[8][2];
		TF2* f2= new TF2("f2","[0]*TMath::Exp(-pow((x-[1])/[2],2)-pow((y-[3])/[4],2))", -1, 1, -1., 1);
		for(int i=0; i<8; ++i){
				for(int j=0; j<2; ++j){
						h2[i][j]= (TH2D*)f->Get(Form("taggingBias_%d_0", i));
						fitter::fitter(h2[i][j], f2);
						hfit[i][j]=fitter::evalToHist(Form("fitHist_%d",i), f2, h2[i][j]);
				}
		}
		TFile *wf = new TFile(dataDumpPath+"tagBias2DFit_JTCSignal.root", "recreate");
		wf->cd();
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<2; ++j){
						if(isNumber) hfit[i][j]->SetName(Form("signal_tagBias2DFit_%d_%d",i,j));
						else hfit[i][j]->SetName(Form("signal_tagBias2DFit_pTweighted_%d_%d",i,j));
						hfit[i][j]->Write();
				}
		}
		wf->Close();
		TH1D* heta1[8][2];
		TH1D* heta2[8][2];
		float x1= h2[1][0]->GetYaxis()->FindBin(-1);
		float x2= h2[1][0]->GetYaxis()->FindBin(0.99);
		float width= h2[1][0]->GetYaxis()->GetBinWidth(1);
		for(int i=1; i<nPt; ++i){
				heta1[i][0]=(TH1D*) h2[i][0]->ProjectionX(Form("deta1_%d",i), x1, x2, "B");
				heta2[i][0]=(TH1D*) hfit[i][0]->ProjectionX(Form("deta2_%d",i), x1, x2, "B");
				heta1[i][0]->Rebin(4);
				heta2[i][0]->Rebin(1);
				heta1[i][0]->Scale(0.25*width);
				heta2[i][0]->Scale(1*width);
				heta1[i][0]->GetXaxis()->SetTitle("d#eta");
		}
		TString lg1 = "bias histo";
		TString lg2 = "fit func";
		TString name = "fittingCheck";
		plot_closure_and_overlay(FigDumpPath+name, lg1, lg2, heta1, heta2,-1, 0.99, 1, 1);
//		hfit->Draw("surf2");

}

void applyCorrection2D(TString name, TString name1,TString corr_name, bool isNumber = 0){
		TH2D* h0[8][2];
		TH2D* h2[8][2];
		getSig(name1, h0, isNumber);
		TH2D* corr[8][2];
		getSig(corr_name, corr, isNumber);
		auto wf = new TFile(dataDumpPath+name+"_JTCSignal.root", "recreate");
		wf->cd();
		for(int i=1; i<nPt; ++i){
				for(int j=0; j<2; ++j){
						if(isNumber) h2[i][j]=(TH2D*) h0[i][j]->Clone("signal_"+name+Form("_%d_%d",i,j));
						else h2[i][j]=(TH2D*) h0[i][j]->Clone("signal_"+name+Form("_pTweighted_%d_%d",i,j));
//						cout<<corr[i][j]->GetName()<<endl;
						h2[i][j]->Add(corr[i][j],-1);
						h2[i][j]->Write();
				}
		}
}

void show_2D_colz(TString name, TH2D* h[8][2]){
		auto cm = new mCanvasLoose("c_"+name, "", 2, 3, 1000, 1000);
		for(int i=1; i<nPt ; ++i){
				h[i][0]->SetAxisRange(-1, 0.99, "X");
				h[i][0]->SetAxisRange(-1, 0.99, "Y");
				cm->drawHist(h[i][0], i/4+1, i%3+int((i%3)==0)*3, "colz");
				cm->CD(i/4+1, i%3+int((i%3)==0)*3);
				gPad->SetLogz();
		}
		cm->SaveAs(FigDumpPath+name+".gif");
}

void rawRatio2D(TString name, TString name1, TString name2, bool isNumber = 0){
		TH2D* h1[8][2];
		TH2D* h2[8][2];
		getRawSig(name1, h1, isNumber);
		getRawSig(name2, h2, isNumber);
		TH2D* ratio[8][2];
		for(int i=1; i<nPt ; ++i){
				ratio[i][0] = (TH2D*) h1[i][0]->Clone(Form("ratio_%d",i));
				//ratio[i][0]->Divide(h2[i][0]);
				ratio[i][0]->Add(h2[i][0], -1);
				ratio[i][0]->GetXaxis()->SetTitle("d#eta");
				ratio[i][0]->GetYaxis()->SetTitle("d#phi");
				ratio[i][0]->SetTitle(track_label[i]);
		}
		show_2D_colz(name, ratio);
}

void getTaggingBias(){
				/*
		TH2D *rawsig[8][2];
		TH2D *rawsig_pTweighted[8][2];
		TH2D *rawsig2[8][2];
		TH2D *rawsig2_pTweighted[8][2];
		getRawSig("gen_gen_PYTHIA_tagged_trueB", rawsig_pTweighted, 0);
		getRawSig("gen_gen_PYTHIA_trueB", rawsig2_pTweighted, 0);
		TH1D* deta[8][2];
		TH1D* dphi[8][2];
		int x1= rawsig_pTweighted[1][0]->GetYaxis()->FindBin(-1);
		int x2= rawsig_pTweighted[1][0]->GetYaxis()->FindBin(0.99);
		TH1D* h;
		JTCUtility *su[8];
		for(int i=1; i<nPt; ++i){
				su[i]=new JTCUtility();
				
				deta[i][0]=(TH1D*) rawsig_pTweighted[i][0]->ProjectionX(Form("deta_%d",i),x1, x2,"E");
				h=(TH1D*) rawsig2_pTweighted[i][0]->ProjectionX(Form("h_%d",i),x1, x2, "E");
				deta[i][0]->Divide(h);
				delete h;
		}	
				*/
}
