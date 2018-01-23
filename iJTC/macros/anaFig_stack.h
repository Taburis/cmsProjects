
#include "../lib/stackHist.h"
void addLabels(TCanvas* c){
		TString cent_label[]={"Cent. 0-30%", "Cent. 30-100%"};
		c->cd(3);
		auto tx = new TLatex();  tx->SetTextSize(.08);
}

void stackPlot_diff(TH1* pb[8][2], TH1* pp[8], TString name){
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
				diff[i]->setRange(-10, 10, "yd");
				for(int j=0; j<npt-drop; ++j){
					   	htm[j] = hdiff[j+drop][i];}
				diff[i]->addDiff(htm, npt-drop);
		}

		auto c = new TCanvas("c", "", 1400, 1000);
		auto tx = new TLatex(); 
		c->SetMargin(0.17, 0.05, 0.15, 0.05);
		c->Divide(3,2, 0, 0);
		for(int i=0; i<2; ++i){
				c->cd(3-i);
				st[i]->drawStack("r");

				c->cd(6-i);
				diff[i]->drawDiff("r");
		}
		c->cd(1);
		st[2]->drawStack("r");
		c->cd(4);
		TLegend* lt = new TLegend(0.1,0.1,1.,0.85);
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

		/*draw caption*/
		c->cd(1);
		tx->SetTextSize(0.07);
		tx->SetTextFont(22);
		tx->DrawLatexNDC(0.24, 0.85, "pp reference");
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
