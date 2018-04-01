
#include "../lib/import_config.h"
#include "../lib/JTCSystError.h"
#include "../lib/systErrorTable.h"
#include "../lib/read_residual_corr.h"
using namespace hallie_reg;

void errSummary_js(){
		TString trk_tag[] = {"0.7 < p_{T}^{track} GeV < 1 ", "1 < p_{T}^{track} GeVTrkPt < 2", 
				"2 < p_{T}^{track} GeVTrkPt < 3", "3 < p_{T}^{track} GeV < 4", 
				"4 < p_{T}^{track} GeV < 8", 
				"8 < p_{T}^{track} GeV < 12", "12 < p_{T}^{track} GeV < 16","16 < p_{T}^{track} GeV < 20", 
				"p_{T}^{track} GeV > 20"};

		py_f = new TFile("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/Jet_Shapes_pTweighted_kurt_pp.root");	
		getPYErr_dr();
		TLatex tx;
		tx.SetTextSize(.1);
		readResidualCorr();
		JTCSystError* pb_dr_err[9][4];
		JTCSystError* pp_dr_err[9];
		auto cpb_dr = new JTCSystCanvas("cpb_dr", "", 9, 5);
		for(int i=0; i<9; ++i){
				//dr 
						cout<<endl;
						cout<<py_dr_syst_err[4][i][4]->GetName()<<endl;
						pp_dr_err[i]=new JTCSystError();
						pp_dr_err[i]->bkg_err=pp_pTweighted_err[i];		
						pp_dr_err[i]->jff_err= jff_err;		
						pp_dr_err[i]->spill_err= 1;	
						pp_dr_err[i]->rela_err= rela_err_pp;	
						pp_dr_err[i]->getAllError(py_dr_syst_err[0][i][4], py_dr_syst_err[4][i][4], NULL, 1);
						cpb_dr->addError(pp_dr_err[i],i+1, 1);

				for(int j=0; j<4; ++j){
						pb_dr_err[i][j]=new JTCSystError();
						pb_dr_err[i][j]->bkg_err=pb_pTweighted_err[i][j];		
						pb_dr_err[i][j]->jff_err= jff_err;		
						pb_dr_err[i][j]->spill_err= spill_err[i];	
						pb_dr_err[i][j]->rela_err= rela_err_pb;	
						pb_dr_err[i][j]->getAllError(py_dr_syst_err[0][i][j],py_dr_syst_err[4][i][j], py_dr_syst_err[3][i][j], 1);
						cpb_dr->addError(pb_dr_err[i][j],i+1, 5-j);
				}
		}
		cpb_dr->drawSummary(0, 0.99);
		for(int i=0; i<9; ++i){
				cpb_dr->cd(i*5+1);
				tx.DrawLatexNDC(0.2, 0.8, trk_tag[i]);
		}

		TFile *wf = TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/js_syst_err.root","recreate");
		for(int i=0; i<9; ++i){
				for(int j=0; j<4; ++j){
						pb_dr_err[i][j]->error->SetName(Form("js_dr_Pb_Syst_Error_%d_%d",i,j));
						pb_dr_err[i][j]->error->Write();
				}
				pp_dr_err[i]->error->SetName(Form("js_dr_pp_Syst_Error_%d",i));
				pp_dr_err[i]->error->Write();
		}
		/*
		*/
		js_f  = TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/Jet_Shapes_pTweighted_v3.root");
		getJS();
		for(int i=0; i<9; ++i){
				for(int j=0; j<5; ++j){
						for(int k=1; k<js_dr_err[i][j]->GetNbinsX()+1; ++k){
								js_dr_err[i][j]->SetBinContent(k, js_dr_err[i][j]->GetBinError(k));
								js_dr_err[i][j]->SetBinError(k, 0);
						}
						cpb_dr->cd(i*5+5-j);
						js_dr_err[i][j]->SetLineColor(kRed-2);
						js_dr_err[i][j]->SetLineStyle(2);
						js_dr_err[i][j]->SetLineWidth(2);
						js_dr_err[i][j]->Draw("same");
				}
		}
		TLegend * l1  = new TLegend(0.4, 0.2, 0.9, 0.7);
		l1->SetFillColor(0);
		l1->SetLineColor(0);
		TLegend * l2  = new TLegend(0.4, 0.3, 0.9, 0.8);
		l2->SetFillColor(0);
		l2->SetLineColor(0);
		l1->AddEntry(cpb_dr->getErr(1,1)->error, "New Syst. Error");
		l1->AddEntry(js_dr_err[1][1], "Old Syst. Error");
		l2->AddEntry(cpb_dr->getErr(1,1)->sig, "Rela. Error");
		l2->AddEntry(cpb_dr->getErr(1,1)->bkg, "ME+BG");
		l2->AddEntry(cpb_dr->getErr(1,1)->jff, "Jff residual");
		l2->AddEntry(cpb_dr->getErr(1,1)->spill, "Spill-over");
		cpb_dr->cd(1);
		l1->Draw();
		cpb_dr->cd(2);
		l2->Draw();
		cpb_dr->SaveAs("js_dr_err_summary_overlay.pdf");

}
