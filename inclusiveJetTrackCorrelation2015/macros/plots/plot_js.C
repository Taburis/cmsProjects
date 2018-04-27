
#include "../../lib/stackHist.h"
#include "../../lib/plottingInput.h"
#include "../../lib/labels_PAS.h"
#include "../../lib/xCanvas.h"
using namespace plottingInput;

TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-100%"};
float rbin[16] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,.985, .99 };
const float drbin [16] = {0.,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1., 1.2};
//only for check
TH1D* changeBin(TString name, TH1D* href){
		TH1D* h=new TH1D(name, "", 14, drbin);
		for(int i =1; i<15; ++i){
				h->SetBinContent(i, href->GetBinContent(i));
				h->SetBinError(i, href->GetBinError(i));
		}
		return h;
}

TH1D* auxi_hist(TH1D* hh){
		TString name = hh->GetName();
		name = name+"_auxi_band";
		TH1D* h= new TH1D(name, "", 15, rbin);
		for(int i =1; i<15; ++i){
				h->SetBinContent(i, hh->GetBinContent(i));
				h->SetBinError(i, hh->GetBinError(i));
		}
		h->SetBinError(15, hh->GetBinError(14));
		h->SetBinContent(15, hh->GetBinContent(14));
		//	h->SetBinError(16, hh->GetBinError(14));
		return h;
}
void plot_js(){

		prePlot_js_dr();

		stackHist *st[5];
		auto tl = new TLatex();

		TH1D* ratio[4];
		TLine *line = new TLine();

		TBox * box = new TBox();
		for(int i=0; i<5; ++i){
				js_dr_err_all[i]->SetFillStyle(1001);
				js_dr_err_all[i]->SetFillColorAlpha(kGray+3, .4);
				js_dr_err_all[i]->SetMarkerStyle(20);
				js_dr_err_all[i]->SetMarkerSize(1.6);
				js_dr_err_all[i]->SetMarkerColor(kBlack);

				st[i]= new stackHist(Form("st_%d",i));
				st[i]->setRange(0., 0.99, "x");
				st[i]->setRange(0.5, 3000, "y");
				for(int j=0; j<9; ++j){
						st[i]->addHist((TH1*) js_dr[j][i]);
				}
		}

		for(int i=0 ;i<4; ++i){
				ratio[i]=(TH1D*)js_dr_err_all[i]->Clone(Form("js_ratio_%d",i));
				ratio[i]->Divide(js_dr_err_all[4]);
				for(int j=0; j<9; ++j){
						//js_dr_err[j][i]->Divide(js_dr_err[j][4]);
				}
		}


		//Color_t co[4] = {kBlue-4, kOrange-3, kGreen+1, kRed+1};
		Color_t co[4] = {kBlue-9, kOrange+1, kGreen+3, kRed+2};
		TString leg[4] = {"p_{T}^{trk} > 0.7 GeV", "p_{T}^{trk} > 2 GeV", "p_{T}^{trk} > 4 GeV", " 12 < p_{T}^{trk}< 16"};
		TLegend *tll = new TLegend(0.25, 0.6, 0.7, 0.98);
		tll->SetTextSize(0.07);
		tll->SetLineColor(kWhite);
		tll->SetFillColor(kWhite);


		TH1D* sub_ratio [4][5];
		TH1D* sub_ratio2[4][5];
		int bin[3] = {0, 2, 5};
		for(int j=0; j<5; ++j){
				for(int i=0; i<3; ++i){
						sub_ratio [i][j] = (TH1D*) js_dr_err[bin[i]][j]->Clone(Form("sub_ratio_%d_%d",i,j));
						sub_ratio2[i][j] = (TH1D*) js_dr    [bin[i]][j]->Clone(Form("sub_ratio2_%d_%d",i,j));
						cout<<i<<", "<<j<<endl;
						cout<<js_dr    [bin[0]][j]->GetName()<<endl;
						for(int k=bin[i]+1; k<9; k++){
								//								cout<<js_dr_err[k][j]->GetName()<<endl;
								sub_ratio [i][j] ->Add(js_dr_err[k][j]);
								sub_ratio2[i][j] ->Add(js_dr    [k][j]);
						}
				}
		}



		auto wf = new TFile("~/cmsProjects/iJTC/dataSet/inclJTC/JS5TeV_HIN_16_020.root","recreate");
		TH1D* js_sub[3][5], *js_pp[9], *jspb[9][4], *jserrpp[9], *jserrpb[9][4];
		wf->cd();
		//	for(int i=0; i<3; ++i){
		//			js_sub[i][4] = changeBin(Form("Rho_pp_%d",i), sub_ratio2[i][4]);
		//			js_sub[i][4]->Write();
		//			for(int j=0; j<4; ++j){
		//					js_sub[i][j] = changeBin(Form("Rho_pb_%d_%d",i,j), sub_ratio2[i][j]);
		//					js_sub[i][j]->Write();
		//			}
		//	}
		TString pts[] = {"0.7", "1", "2", "3", "4", "8", "12", "16", "20", "999"};
		TString cents[] = {"0", "10", "30", "50", "100"};
		for(int i=0; i<9; ++i){
				js_pp[i] = changeBin(Form("JS_pp_%d",i), js_dr[i][4]);
				js_pp[i]->SetTitle(pts[i]+" < pT < "+pts[i+1]);
				js_pp[i]->Write();
				jserrpp[i]=changeBin(Form("SystError_pp_%d", i), js_dr_err[i][4]);
				jserrpp[i]->SetTitle(pts[i]+" < pT < "+pts[i+1]);
				jserrpp[i]->Write();
		}
		for(int i=0; i<9; ++i){
				for(int j=0; j<4; ++j){
						jspb[i][j] = changeBin(Form("JS_pb_%d_%d",i, j), js_dr[i][j]);
						jspb[i][j]->SetTitle(pts[i]+" < pT <"+pts[i+1]+", Cent: "+cents[j]+"-"+cents[j+1]);
						jspb[i][j]->Write();
						jserrpb[i][j]=changeBin(Form("SystError_pb_%d_%d", i, j), js_dr_err[i][j]);
						jserrpb[i][j]->SetTitle(pts[i]+" < pT < "+pts[i+1]+", Cent: "+cents[j]+"-"+cents[j+1]);
						jserrpb[i][j]->Write();
				}
		}
		wf->Close();

		for(int i=0; i<3; ++i){
				for(int j=0; j<4; ++j){
						//cout<<i<<", "<<j<<endl;
						sub_ratio [i][j]->Divide(sub_ratio [i][4]);
						sub_ratio2[i][j]->Divide(sub_ratio2[i][4]);
						sub_ratio2[i][j]->SetAxisRange(0, .99, "X");
						//combind the syst and stat err togther
						/*
						   for(int k=1 ;k<sub_ratio[i][j]->GetNbinsX()+1; ++k){
						   float err1 = sub_ratio[i][j]->GetBinError(k);
						   float err2 = sub_ratio2[i][j]->GetBinError(k);
						   float err = pow(pow(err1, 2)+pow(err2,2), .5);
						//	sub_ratio[i][j]->SetBinError(k, err);
						}
						*/
						sub_ratio [i][j]= auxi_hist(sub_ratio[i][j]);
						sub_ratio [i][j]->SetLineColor(co[i]);
						sub_ratio2[i][j]->SetLineColor(co[i]);
						sub_ratio2[i][j]->SetMarkerColor(co[i]);
						sub_ratio2[i][j]->SetMarkerStyle(20);
						sub_ratio2[i][j]->SetMarkerSize(1);
						sub_ratio[i][j]->SetFillStyle(1001);
						sub_ratio[i][j]->SetFillColorAlpha(co[i],0.5);
						sub_ratio[i][j]->SetMarkerColorAlpha(co[i],0.3);
				}
				tll->AddEntry(sub_ratio[i][0], leg[i]);
		}

		cout<<sub_ratio[0][0]->GetName()<<endl;
		cout<<"err = "<<sub_ratio[0][0]->GetBinError(sub_ratio[0][0]->FindBin(0.95))<<endl;;
		cout<<js_dr_err_all[4]->Integral("width")/js_dr_err_all[0]->Integral("width")<<endl;

		auto c = new auxi_canvas("c", "", 2500, 2000);
		c->SetMargin(0.06, 0.01, 0.08, 0.02);
		c->divide(3,4);

		for(int i=0; i<4; ++i){
				c->CD(8-i);
				gPad->SetLogy();
				st[i]->drawStack();
				st[i]->hst->GetYaxis()->SetNdivisions(505);
				st[i]->hst->GetYaxis()->SetLabelSize(0.08);
				st[i]->hst->GetYaxis()->SetTitleOffset(.9);
				st[i]->hst->GetYaxis()->SetTitleSize(0.1);
				st[i]->hst->GetYaxis()->SetTitle("#Rho(#Deltar)");
				st[i]->hst->Draw();
				js_dr_err_all[i]->Draw("same e2");
				if(i==3 ){
						tl->SetTextFont(22);
						tl->SetTextSize(.085);
						tl->DrawLatexNDC(.35, 0.9, "PbPb");
						tl->DrawLatexNDC(.35, .82, cent_lab[i]);
				}
				else {
						tl->SetTextFont(22);
						tl->SetTextSize(.09);
						tl->DrawLatexNDC(.19, 0.9, "PbPb");
						tl->DrawLatexNDC(.19, .82, cent_lab[i]);
				}
				c->CD(12-i);
				//ratio[i]->GetYaxis()->SetNdivisions(505);
				sub_ratio[0][i]->GetXaxis()->SetTitleOffset(1.1);
				sub_ratio[0][i]->GetXaxis()->SetTitle("#Deltar");
				sub_ratio[0][i]->SetAxisRange(0., 3.2, "Y");
				sub_ratio[0][i]->SetAxisRange(0., .99, "X");
				if( i<3 )  {
						sub_ratio[0][i]->GetXaxis()->SetTitleOffset(0.7);
						sub_ratio[0][i]->GetXaxis()->SetTitleSize(0.11);
						sub_ratio[0][i]->GetXaxis()->SetNdivisions(505);
						sub_ratio[0][i]->GetXaxis()->SetLabelSize(0.08);
						sub_ratio[0][i]->GetXaxis()->SetLabelOffset(0.004);
						sub_ratio[0][i]->GetXaxis()->SetNdivisions(505);
				}
				if(i==3 ){
						sub_ratio[0][i]->GetXaxis()->SetTitleOffset(0.92);
						sub_ratio[0][i]->GetXaxis()->SetTitleSize(0.086);
						sub_ratio[0][i]->GetXaxis()->SetNdivisions(505);
						sub_ratio[0][i]->GetXaxis()->SetLabelOffset(0.02);
						sub_ratio[0][i]->GetXaxis()->SetLabelSize(0.07);
						sub_ratio[0][i]->GetYaxis()->SetNdivisions(505);
						sub_ratio[0][i]->GetYaxis()->SetLabelSize(0.07);
						sub_ratio[0][i]->GetYaxis()->SetLabelOffset(0.02);
						sub_ratio[0][i]->GetYaxis()->SetTitleOffset(0.9);
						sub_ratio[0][i]->GetYaxis()->SetTitleSize(0.08);
						sub_ratio[0][i]->GetYaxis()->SetTitle("#Rho(#Deltar)_{PbPb}/#Rho(#Deltar)_{pp}");
						tl->SetTextSize(0.073);
						//tl->DrawLatexNDC(.25, .92, "PbPb - pp");
						//tl->DrawLatexNDC(.25, .84, cent_lab[i]);
				}
				else{
						//tl->SetTextSize(0.09);
						//tl->DrawLatexNDC(.05, .92, "PbPb - pp");
						//tl->DrawLatexNDC(.05, .84, cent_lab[i]);
				}
				sub_ratio[0][i]->GetYaxis()->CenterTitle();
				sub_ratio[0][i]->GetXaxis()->CenterTitle();
				sub_ratio[0][i]->SetStats(0);
				sub_ratio[0][i]->Draw("same ce5");
				sub_ratio2[0][i]->Draw("same");
				line->SetLineStyle(2);
				line->DrawLine(0, 1, 1, 1);
				for(int j=1; j<3; ++j){
						sub_ratio[j][i]->Draw("same ce5");
						sub_ratio2[j][i]->Draw("same");
				}
		}
		c->CD(9);
		tll->Draw();
		c->CD(1);
		gPad->SetLogy();
		st[4]->drawStack();
		st[4]->hst->GetYaxis()->SetNdivisions(505);
		st[4]->hst->GetYaxis()->SetLabelSize(0.08);
		st[4]->hst->GetYaxis()->SetTitleOffset(.9);
		st[4]->hst->GetYaxis()->SetTitleSize(0.1);
		st[4]->hst->GetYaxis()->SetTitle("#Rho(#Deltar)");
		tl->SetTextSize(.085);
		tl->DrawLatexNDC(.35, .88, "pp reference");
		js_dr_err_all[4]->Draw("same e2");

		TLegend* lt1 = new TLegend(0.01,0.1,1.,0.5);
		TLegend* lt2 = new TLegend(0.0,0.1,1.,0.5);
		TLegend* lt3 = new TLegend(0.0,0.1,1 ,0.5);
		TLegend* lt5 = new TLegend(0.01,0.5,1 ,.68);
		TLegend* lt4 = new TLegend(0.25,0.75,.9,0.95);
		lt1->SetTextSize(0.07);
		lt1->SetLineColor(kWhite);
		lt1->SetFillColor(kWhite);
		lt2->SetTextSize(0.07);
		lt2->SetLineColor(kWhite);
		lt2->SetFillColor(kWhite);
		lt3->SetTextSize(0.07);
		lt3->SetLineColor(kWhite);
		lt3->SetFillColor(kWhite);
		lt4->SetTextSize(0.07);
		lt4->SetLineColor(kWhite);
		lt4->SetFillColor(kWhite);
		lt5->SetTextSize(0.07);
		lt5->SetLineColor(kWhite);
		lt5->SetFillColor(kWhite);

		lt1->AddEntry(st[4]->hist_trunk.at(0), "0.7 < p_{T}^{trk}< 1 GeV","f");
		lt1->AddEntry(st[4]->hist_trunk.at(1), "1 < p_{T}^{trk}< 2 GeV","f");
		lt1->AddEntry(st[4]->hist_trunk.at(2), "2 < p_{T}^{trk}< 3 GeV","f");

		lt2->AddEntry(st[4]->hist_trunk.at(3), "3 < p_{T}^{trk}< 4 GeV","f");
		lt2->AddEntry(st[4]->hist_trunk.at(4), "4 < p_{T}^{trk}< 8 GeV","f");
		lt2->AddEntry(st[4]->hist_trunk.at(5), "8 < p_{T}^{trk}< 12 GeV","f");

		lt3->AddEntry(st[4]->hist_trunk.at(6), "12 < p_{T}^{trk}< 16 GeV","f");
		lt3->AddEntry(st[4]->hist_trunk.at(7), "16 < p_{T}^{trk}< 20 GeV","f");
		lt3->AddEntry(st[4]->hist_trunk.at(8), "20 < p_{T}^{trk}< 300 GeV","f");

		c->CD(1);
		TLegend* lt6 = new TLegend(0.3,0.7,0.98,0.85);
		lt6->SetTextSize(0.07);
		lt6->SetLineColor(kWhite);
		lt6->SetFillColor(kWhite);
		lt6->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
		lt6->Draw();

		c->CD(2);
		//	lt5->AddEntry(js_dr_err_all[0], "0.7 < p_{T}^{track}< 20 GeV","lpfe");
		//	lt5->Draw();
		line->SetLineStyle(1);
		line->DrawLineNDC(0, 0, 0, 1);
		lt1->Draw();
		c->CD(3);
		lt2->Draw();
		c->CD(4);
		lt3->Draw();

		c->cd(0);
		tl->SetTextFont(62);
		tl->SetTextSize(0.035);
		tl->DrawLatexNDC(0.45, 0.94, "CMS");

		tl->SetTextFont(42);
		tl->SetTextSize(0.035);
		tl->DrawLatexNDC(0.55, 0.94, "Radial momentum distribution");

		tl->SetTextSize(0.03);
		tl->DrawLatexNDC(.4, .9, "pp 27.4 pb^{-1} (5.02 TeV)  PbPb 404 #mub^{-1} (5.02 TeV)");
		tl->SetTextSize(0.025);
		tl->DrawLatexNDC(.45, .86, "anti-k_{T} R=0.4 jets, p_{T}> 120 GeV, |#eta_{jet}|<1.6");
		box->SetFillColor(kWhite);
		c->cd(0);
		box->DrawBox(0.285,.047, 0.3, 0.072);
		tl->SetTextSize(.0268);
		box->DrawBox(0.518,.047, 0.533, 0.072);
		box->DrawBox(0.75,.047, 0.765, 0.072);

		tl->DrawLatex(0.29,  0.0533, "0");
		tl->DrawLatex(0.523, 0.0533, "0");
		tl->DrawLatex(0.755, 0.0533, "0");
		tl->DrawLatex(0.985, 0.0533, "1");

		c->SaveAs("js_dr_new.eps");
		c->SaveAs("js_dr_new.pdf");
		//c->SaveAs("js_dr_v3.eps");
		//c->SaveAs("js_dr_v3.pdf");

}


