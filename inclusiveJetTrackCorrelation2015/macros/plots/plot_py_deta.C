
#include "../../lib/stackHist.h"
#include "../../lib/plottingInput.h"
#include "../../lib/labels_PAS.h"
#include "../../lib/xCanvas.h"
using namespace plottingInput;

TString cent_lab[4] = {"0-10%", "10-30%", "30-50%", "50-100%"};
void plot_py_deta(){
		prePlot_py_deta();
	stackHist *st[5];
	stackHist *diff[5];
	labels_PAS* lb = new labels_PAS();
	auto tl = new TLatex();



	for(int i=0; i<5; ++i){

		py_deta_err_all[i]->SetFillStyle(1001);
		py_deta_err_all[i]->SetFillColorAlpha(kGray+3,.4);
		py_deta_err_all[i]->SetMarkerStyle(20);
		py_deta_err_all[i]->SetMarkerSize(1.1);
		py_deta_err_all[i]->SetMarkerColor(kWhite);

		st[i]= new stackHist(Form("st_%d",i));
		st[i]->setRange(-1.5, 1.5, "x");
		st[i]->setRange(-2, 32, "y");
		for(int j=0; j<8; ++j){
			st[i]->addHist((TH1*) py_deta[j][i]);
		}
	}
	py_deta_all[4]->SetLineColor(kBlack);
	py_deta_all[4]->SetLineWidth(1);
	py_deta_all[4]->SetMarkerStyle(20);
	py_deta_all[4]->SetMarkerColor(kWhite);
	for(int i=0 ;i<4; ++i){
		py_deta_all[i]->SetLineColor(kBlack);
		py_deta_all[i]->SetLineWidth(1);
		py_deta_all[i]->SetMarkerStyle(20);
		py_deta_all[i]->SetMarkerColor(kWhite);
		py_deta_all[i]->SetMarkerSize(1);

		py_deta_diff_all[i]->SetLineColor(kBlack);
		py_deta_diff_all[i]->SetLineWidth(1);
		py_deta_diff_all[i]->SetMarkerStyle(20);
		py_deta_diff_all[i]->SetMarkerColor(kWhite);
		py_deta_diff_all[i]->SetMarkerSize(1);

		py_deta_diff_err[i]->SetFillStyle(1001);
		py_deta_diff_err[i]->SetFillColorAlpha(kGray+3, 0.4);
		py_deta_diff_err[i]->SetMarkerStyle(20);
		py_deta_diff_err[i]->SetMarkerSize(1.1);
		py_deta_diff_err[i]->SetMarkerColor(kWhite);
	}
	//add diff
	TH1D** hdif;
	for(int i=0; i<4; ++i){
		diff[i]= new stackHist(Form("dif_%d",i));
		diff[i]->setRange(-1.5, 1.5, "xd");
		diff[i]->setRange(-0.9, 11, "yd");
		hdif = new TH1D*[8];
		for(int j=0; j<8; ++j){
			hdif[j] = (TH1D*) py_deta[j][i]->Clone(Form("hdif_%d_%d",i, j));
			hdif[j]->Add(py_deta[j][4], -1);
		}
		diff[i]->addDiff((TH1**)hdif, 8);
	}
	auto c = new auxi_canvas("c", "", 2500, 2000);
	c->SetMargin(0.06, 0.01, 0.08, 0.02);
	c->divide(3,4);

	for(int i=0; i<4; ++i){
		c->CD(8-i);
		st[i]->drawStack("r");
		st[i]->hst->GetYaxis()->SetNdivisions(505);
		st[i]->hst->GetYaxis()->SetLabelSize(0.08);
		st[i]->hst->GetYaxis()->SetTitleOffset(0.8);
		st[i]->hst->GetYaxis()->SetTitleSize(0.1);
		st[i]->hst->GetYaxis()->SetTitle("Y=#frac{1}{N_{jet}}#frac{dN}{d#Delta#eta}");
		st[i]->hst->Draw();
	//	py_deta_all[i]->Draw("same");
		py_deta_err_all[i]->Draw("same e2");
		if(i==3 ){
			tl->SetTextFont(22);
			tl->SetTextSize(.085);
			tl->DrawLatexNDC(.25, 0.85, "PbPb");
			tl->DrawLatexNDC(.25, .72, cent_lab[i]);
		}
		else {
			tl->SetTextFont(22);
			tl->SetTextSize(.09);
			tl->DrawLatexNDC(.07, 0.85, "PbPb");
			tl->DrawLatexNDC(.07, .72, cent_lab[i]);
		}
		c->CD(12-i);
		diff[i]->drawDiff("r");
		diff[i]->hst_up->GetXaxis()->SetTitleOffset(1);
		diff[i]->hst_up->GetXaxis()->SetTitle("#Delta#eta");
		if( i<3 )  {
			diff[i]->hst_up->GetXaxis()->SetTitleOffset(0.8);
			diff[i]->hst_up->GetXaxis()->SetTitleSize(0.09);
			diff[i]->hst_up->GetXaxis()->SetLabelSize(0.08);
			diff[i]->hst_up->GetXaxis()->SetLabelOffset(0.001);
		}
		diff[i]->hst_up->Draw();
		diff[i]->hst_down->Draw("same");
		if(i==3 ){
			diff[i]->hst_up->GetXaxis()->SetTitleSize(.075);
			diff[i]->hst_up->GetXaxis()->SetTitleOffset(.95);
			diff[i]->hst_up->GetYaxis()->SetTitleOffset(0.8);
			diff[i]->hst_up->GetYaxis()->SetTitleSize(0.08);
			diff[i]->hst_up->GetXaxis()->SetLabelSize(0.063);
			diff[i]->hst_up->GetXaxis()->SetLabelOffset(0.013);
			diff[i]->hst_up->GetYaxis()->SetTitle("Y_{PbPb} - Y_{pp}");
			tl->SetTextSize(0.073);
			tl->DrawLatexNDC(.25, .85, "PbPb - pp");
			tl->DrawLatexNDC(.25, .72, cent_lab[i]);
		}
		else{
			tl->SetTextSize(0.09);
			tl->DrawLatexNDC(.05, .85, "PbPb - pp");
			tl->DrawLatexNDC(.05, .72, cent_lab[i]);
		}
		py_deta_diff_all[i]->Draw("same");
		py_deta_diff_err[i]->Draw("same e2");
	}
	c->CD(1);
	st[4]->drawStack("r");
	st[4]->hst->GetYaxis()->SetNdivisions(505);
	st[4]->hst->GetYaxis()->SetLabelSize(0.08);
	st[4]->hst->GetYaxis()->SetTitleOffset(0.8);
	st[4]->hst->GetYaxis()->SetTitleSize(0.1);
	st[4]->hst->GetYaxis()->SetTitle("Y=#frac{1}{N_{jet}}#frac{dN}{d#Delta#eta}");
	tl->SetTextSize(.085);
	tl->DrawLatexNDC(.25, .85, "pp reference");
	py_deta_err_all[4]->Draw("same e2");

	TLegend* lt1 = new TLegend(0.1,0.1,1.,0.5);
	TLegend* lt2 = new TLegend(0.0,0.1,1.,0.5);
	TLegend* lt3 = new TLegend(0.0,0.1,1 ,0.5);
	lt1->SetTextSize(0.07);
	lt1->SetLineColor(kWhite);
	lt1->SetFillColor(kWhite);
	lt2->SetTextSize(0.07);
	lt2->SetLineColor(kWhite);
	lt2->SetFillColor(kWhite);
	lt3->SetTextSize(0.07);
	lt3->SetLineColor(kWhite);
	lt3->SetFillColor(kWhite);

	lt1->AddEntry(st[4]->hist_trunk.at(0), "0.7 < p_{T}^{trk}< 1 GeV","f");
	lt1->AddEntry(st[4]->hist_trunk.at(1), "1 < p_{T}^{trk}< 2 GeV","f");
	lt1->AddEntry(st[4]->hist_trunk.at(2), "2 < p_{T}^{trk}< 3 GeV","f");

	lt2->AddEntry(st[4]->hist_trunk.at(3), "3 < p_{T}^{trk}< 4 GeV","f");
	lt2->AddEntry(st[4]->hist_trunk.at(4), "4 < p_{T}^{trk}< 8 GeV","f");
	lt2->AddEntry(st[4]->hist_trunk.at(5), "8 < p_{T}^{trk}< 12 GeV","f");

	lt3->AddEntry(st[4]->hist_trunk.at(6), "12 < p_{T}^{trk}< 16 GeV","f");
	lt3->AddEntry(st[4]->hist_trunk.at(7), "16 < p_{T}^{trk}< 20 GeV","f");
	lt3->AddEntry(py_deta_err_all[0], "0.7 < p_{T}^{trk}< 20 GeV","lpfe");
	c->CD(2);

	TLine *line = new TLine();
	line->DrawLineNDC(0, 0, 0, 1);
	lt1->Draw();
	c->CD(3);
	lt2->Draw();
	c->CD(4);
	lt3->Draw();

	c->cd(0);
	lb->drawText("CMS", 0.45, 0.94, 0);
	lb->drawText("Particle yield vs. #Delta#eta", 0.55, 0.94, 3);

	lb->drawText("pp 27.4 pb^{-1} (5.02 TeV)  PbPb 404 #mub^{-1} (5.02 TeV)", 0.38, 0.9, 4);
	lb->drawText("anti-k_{T} R=0.4 jets, p_{T}> 120 GeV, |#eta_{jet}|<1.6", 0.45, 0.86, 5);
	//	lb->drawText("(p_{T}> 120 GeV, |#eta_{jet}|<1.6)", 0.2, 0.25, 4);
	c->SaveAs("py_deta_new.eps");
	c->SaveAs("py_deta_new.pdf");
	//auto cc = new TCanvas("cc", "", 600, 500);
	//cc->cd();
	//py_deta_all[1]->Draw();
}
