
#include "../lib/stackHist.h"
void addLabels(TCanvas* c){
		TString cent_label[]={"Cent. 0-30%", "Cent. 30-100%"};
		c->cd(3);
		auto tx = new TLatex();  tx->SetTextSize(.08);
}

void stackPlot_diff(TH1* pb[8][2], TH1* pp[8], TString name){
		float ydmin = -9.9, ydmax= 17;
		const int npt =7; int drop = 1;
		TH1* h[8][3];
		TH1* hdiff[8][2];
		for(int j=drop; j<npt; ++j){
				h[j][2] = pp[j];
				for(int i=0; i<2; ++i){
						h[j][i]= pb[j][i];
						hdiff[j][i] = (TH1*)pb[j][i]->Clone(Form("hdiff_%d_%d",j,i));
						hdiff[j][i]->Add(pp[j],-1);
				}
		}

		stackHist *st[3];
		stackHist *diff[2];
		for(int i=0; i<3; ++i){
				st[i]= new stackHist(Form("st_%d",i));
				st[i]->setRange(0, 0.99, "x");
				st[i]->setRange(-5, 40, "y");
				for(int j=drop; j<npt; ++j){
						st[i]->addHist((TH1*) h[j][i]);
				}
		}
		for(int i=0; i<2; ++i){
				TH1** htm = new TH1*[8];
				diff[i]= new stackHist(Form("diff_%d", i));
				diff[i]->setRange(0, 0.99, "xd");
				diff[i]->setRange(ydmin, ydmax, "yd");
				for(int j=0; j<npt-drop; ++j){
					   	htm[j] = hdiff[j+drop][i];}
				diff[i]->addDiff(htm, npt-drop);
		}

		auto c = new TCanvas("c", "", 1400, 1000);
		auto tx = new TLatex(); 
		c->SetMargin(0.17, 0.05, 0.15, 0.05);
		c->Divide(3,2, 0, 0);
		c->cd(1);
		gPad->SetLeftMargin(0.2);
		for(int i=0; i<2; ++i){
				c->cd(3-i);
				st[i]->drawStack("r");

				c->cd(6-i);
				gPad->SetBottomMargin(0.18);
				diff[i]->drawDiff("r");
				diff[i]->hst_up->GetXaxis()->SetTitle("#Deltar");
				diff[i]->hst_up->GetXaxis()->CenterTitle();
				diff[i]->hst_up->GetXaxis()->SetTitleOffset(0.98);
				diff[i]->hst_up->GetXaxis()->SetTitleSize(0.07);
				diff[i]->hst_up->GetXaxis()->SetLabelOffset(0.01);
		}
		c->cd(1);
		st[2]->drawStack("r");
		c->cd(4);
		TLegend* lt = new TLegend(0.1,0.18,1.,0.75);
		lt->SetTextSize(0.07);
		lt->SetLineColor(kWhite);
		lt->SetFillColor(kWhite);
		lt->AddEntry(st[2]->hist_trunk.at(0), "1 < p_{T}^{trk}< 2 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(1), "2 < p_{T}^{trk}< 3 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(2), "3 < p_{T}^{trk}< 4 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(3), "4 < p_{T}^{trk}< 8 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(4), "8 < p_{T}^{trk}< 12 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(5), "12 < p_{T}^{trk}< 16 GeV","f");
		lt->Draw();
		/*axis setting*/
		c->cd(1);
		st[2]->hst->GetYaxis()->SetTitle("Y=#frac{1}{N_{jet}}#frac{dN}{d#Deltar}");
		st[2]->hst->GetYaxis()->SetNdivisions(505);
		st[2]->hst->GetYaxis()->SetTitleOffset(1.05);
		st[2]->hst->Draw();
		auto axis = new TGaxis();
		c->cd(4);
		gPad->SetBottomMargin(0.18);
		axis->SetLabelFont(42);
		axis->SetTitleFont(42);
		axis->SetTitle("Y_{PbPb} - Y_{pp}");
		axis->SetTitleSize(0.07);
		axis->SetTitleOffset(0.98);
		axis->CenterTitle();
		axis->SetLabelSize(0.07); axis->SetLabelOffset(0.01);
		axis->DrawAxis(1, 0.18, 1, 1, ydmin, ydmax, 510);

		axis->SetTitle("#Deltar");
		axis->DrawAxis(0.2, 1, 1, 1, 0, 1, 505);


		/*draw caption*/
		c->cd(1);
		tx->SetTextSize(0.07);
		tx->SetTextFont(22);
		tx->DrawLatexNDC(0.26, 0.86, "pp reference");
		c->cd(2);
		tx->DrawLatexNDC(0.1, 0.85, "PbPb Cent. 30-100%");
		c->cd(3);
		tx->DrawLatexNDC(0.1, 0.85, "PbPb Cent. 0-30%");
		c->cd(5);
		tx->DrawLatexNDC(0.1, 0.85, "PbPb - pp");


		/*
		lt->AddEntry(st[2]->hist_trunk.at(0), "0.7 < p_{T}^{trk}< 1 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(1), "1 < p_{T}^{trk}< 2 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(2), "2 < p_{T}^{trk}< 3 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(3), "3 < p_{T}^{trk}< 4 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(4), "4 < p_{T}^{trk}< 8 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(5), "8 < p_{T}^{trk}< 12 GeV","f");
		lt->AddEntry(st[2]->hist_trunk.at(6), "12 < p_{T}^{trk}< 16 GeV","f");
		*/
		c->SaveAs(name);
}
