
#include "../../lib/stackHist.h"
#include "../../lib/plottingInput.h"
#include "../../lib/labels_PAS.h"
#include "../../lib/xCanvas.h"
using namespace plottingInput;

TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-100%"};
void plot_js_normal(){
	prePlot_js_dr();
	stackHist *st[5];
	labels_PAS* lb = new labels_PAS();
	TBox * box = new TBox();
	auto tl = new TLatex();

	TH1D* ratio[4];
	TH1D* ratio_auxi[4];
	TLine *line = new TLine();

	for(int i=0; i<5; ++i){
		js_dr_err_all[i]->SetFillStyle(1001);
		js_dr_err_all[i]->SetFillColorAlpha(kGray+3,.4);
		js_dr_err_all[i]->SetMarkerStyle(20);
		js_dr_err_all[i]->SetMarkerSize(1.6);
		js_dr_err_all[i]->SetMarkerColor(kBlack);
		js_dr_err_all[i]->SetLineColor(kBlack);

		st[i]= new stackHist(Form("st_%d",i));
		st[i]->setRange(0., 0.99, "x");
		st[i]->setRange(0.005, 30, "y");
		float fac = js_dr_all[i]->Integral("width");
		cout<<fac<<endl;
		for(int j=0; j<9; ++j){
			js_dr[j][i]->Scale(1.0/fac);
			st[i]->addHist((TH1*) js_dr[j][i]);
		}
		js_dr_err_all[i]->Scale(1.0/fac);
		js_dr_all[i]->Scale(1.0/fac);
		js_dr_all[i]->SetMarkerColor(0);
	}

	for(int i=0 ;i<4; ++i){
		ratio[i]=(TH1D*)js_dr_all[i]->Clone(Form("js_ratio_%d",i));
		ratio[i]->Divide(js_dr_all[4]);
		ratio_auxi[i]=(TH1D*) js_dr_err_all[i]->Clone(Form("js_syst_ratio_%d",i));
		ratio_auxi[i]->Divide(js_dr_err_all[4]);
		ratio_auxi[i]->SetFillStyle(1001);
		ratio_auxi[i]->SetFillColorAlpha(kGray+3, .4);
		ratio_auxi[i]->SetMarkerColor(0);
//		ratio_auxi[i]->SetLineColor(kBlack);
	}
	auto c = new auxi_canvas("c", "", 2500, 2000);
	c->SetMargin(0.06, 0.01, 0.08, 0.02);
	c->divide(3,4);

	cout<<ratio_auxi[0]->GetName()<<endl;
	cout<<"err = "<<ratio_auxi[0]->GetBinError(ratio_auxi[0]->FindBin(0.95))<<endl;;
	
	for(int i=0; i<4; ++i){
		c->CD(8-i);
		gPad->SetLogy();
		st[i]->drawStack();
		st[i]->hst->GetYaxis()->SetNdivisions(505);
		st[i]->hst->GetYaxis()->SetLabelSize(0.08);
		st[i]->hst->GetYaxis()->SetTitleOffset(.9);
		st[i]->hst->GetYaxis()->SetTitleSize(0.1);
		st[i]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
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
		ratio[i]->GetXaxis()->SetTitleOffset(1.1);
		ratio[i]->GetXaxis()->SetTitle("#Deltar");
		ratio[i]->SetAxisRange(0., 3.2, "Y");
		ratio[i]->SetAxisRange(0, .99, "X");
		if( i<3 )  {
			ratio[i]->GetXaxis()->SetTitleOffset(0.7);
			ratio[i]->GetXaxis()->SetTitleSize(0.11);
			ratio[i]->GetXaxis()->SetNdivisions(505);
			ratio[i]->GetXaxis()->SetLabelSize(0.08);
			ratio[i]->GetYaxis()->SetNdivisions(505);
		}
		if(i==3 ){
			ratio[i]->GetXaxis()->SetTitleOffset(0.92);
			ratio[i]->GetXaxis()->SetTitleSize(0.086);
			ratio[i]->GetXaxis()->SetNdivisions(505);
			ratio[i]->GetXaxis()->SetLabelOffset(0.016);
			ratio[i]->GetXaxis()->SetLabelSize(0.064);
			ratio[i]->GetYaxis()->SetNdivisions(505);
			ratio[i]->GetYaxis()->SetLabelOffset(0.02);
			ratio[i]->GetYaxis()->SetLabelSize(0.07);
			ratio[i]->GetYaxis()->SetTitleOffset(0.9);
			ratio[i]->GetYaxis()->SetTitleSize(0.08);
			ratio[i]->GetYaxis()->SetTitle("#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
			tl->SetTextSize(0.073);
			//tl->DrawLatexNDC(.25, .92, "PbPb - pp");
			//tl->DrawLatexNDC(.25, .84, cent_lab[i]);
		}
		else{
			//tl->SetTextSize(0.09);
			//tl->DrawLatexNDC(.05, .92, "PbPb - pp");
			//tl->DrawLatexNDC(.05, .84, cent_lab[i]);
		}
		ratio[i]->GetYaxis()->CenterTitle();
		ratio[i]->GetXaxis()->CenterTitle();
		ratio[i]->SetStats(0);
		ratio[i]->Draw("same");
		ratio_auxi[i]->Draw("same e2");
		line->SetLineStyle(2);
		line->DrawLine(0, 1, 1, 1);
	}
	c->CD(1);
	gPad->SetLogy();
	st[4]->drawStack();
	st[4]->hst->GetYaxis()->SetNdivisions(505);
	st[4]->hst->GetYaxis()->SetLabelSize(0.08);
	st[4]->hst->GetYaxis()->SetTitleOffset(.9);
	st[4]->hst->GetYaxis()->SetTitleSize(0.1);
	st[4]->hst->GetYaxis()->SetTitle("#rho(#Deltar)");
	tl->SetTextSize(.085);
	tl->DrawLatexNDC(.35, .88, "pp reference");
	js_dr_err_all[4]->Draw("same e2");

	TLegend* lt1 = new TLegend(0.01,0.1,1.,0.5);
	TLegend* lt2 = new TLegend(0.0,0.1,1.,0.5);
	TLegend* lt3 = new TLegend(0.0,0.1,1 ,0.5);
	TLegend* lt5 = new TLegend(0.01,0.5,1 ,.68);
	TLegend* lt4 = new TLegend(0.25,0.85,.85,0.95);
	lt1->SetTextSize(0.07);
	lt1->SetLineColor(kWhite);
	lt1->SetFillColor(kWhite);
	lt2->SetTextSize(0.07);
	lt2->SetLineColor(kWhite);
	lt2->SetFillColor(kWhite);
	lt3->SetTextSize(0.07);
	lt3->SetLineColor(kWhite);
	lt3->SetFillColor(kWhite);
	lt4->SetTextSize(0.06);
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

	//lt4->AddEntry(ratio_auxi[0], "#rho(#Deltar)_{PbPb}/#rho(#Deltar)_{pp}");
	lt4->AddEntry(ratio_auxi[0], "0.7 < p_{T}^{trk}< 300 GeV","lpfe");
	c->CD(9);
	lt4->Draw();

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
	tl->DrawLatexNDC(0.5, 0.94, "CMS");
	
	tl->SetTextFont(42);
	tl->SetTextSize(0.035);
	tl->DrawLatexNDC(0.6, 0.94, "Jet shapes");

	tl->SetTextSize(0.03);
	tl->DrawLatexNDC(.4, .9, "pp 27.4 pb^{-1} (5.02 TeV)  PbPb 404 #mub^{-1} (5.02 TeV)");
	tl->SetTextSize(0.025);
	tl->DrawLatexNDC(.45, .86, "anti-k_{T} R=0.4 jets, p_{T}> 120 GeV, |#eta_{jet}|<1.6");
	//	lb->drawText("(p_{T}> 120 GeV, |#eta_{jet}|<1.6)", 0.2, 0.25, 4);
	box->SetFillColor(kWhite);
	c->cd(0);
	box->DrawBox(0.285,.047, 0.3, 0.072);
	tl->SetTextSize(.025);
	tl->DrawLatex(0.29, 0.055, "0");
	box->DrawBox(0.518,.047, 0.533, 0.072);
	box->DrawBox(0.75,.047, 0.765, 0.072);

	tl->DrawLatex(0.523, 0.055, "0");
	tl->DrawLatex(0.755, 0.055, "0");
	tl->DrawLatex(0.985, 0.055, "1");

	c->SaveAs("js_dr_normal_new.eps");
	c->SaveAs("js_dr_normal_new.pdf");
	//c->SaveAs("js_dr_normal_v3.eps");
	//c->SaveAs("js_dr_normal_v3.pdf");
}
