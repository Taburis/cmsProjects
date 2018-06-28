
#include "standardSequence.C"

void tkCorrection(bool isHI = 0){
		TString prefix = isHI ? "PH_": ""; 
		TString cap1 = "GenGen", cap = "GenJet_GenTrack";
		auto sig_tb_gg = read_flatten<TH2D>(dataDumpPath+prefix+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_noCorr");
		auto sig_tb_gg_pt = read_flatten<TH2D>(dataDumpPath+prefix+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_pTweighted_noCorr");
		cap1 = "GenRec", cap = "GenJet_RecoTrack";
		auto sig_tb_gr = read_flatten<TH2D>(dataDumpPath+prefix+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_noCorr");
		auto sig_tb_gr_pt = read_flatten<TH2D>(dataDumpPath+prefix+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_pTweighted_noCorr");
		auto dr_gg = getDr("dr_gg", sig_tb_gg);
		auto dr_gr = getDr("dr_gr", sig_tb_gr);
		auto dr_gg_pt = getDr("dr_gg_pt", sig_tb_gg_pt);
		auto dr_gr_pt = getDr("dr_gr_pt", sig_tb_gr_pt);

		cap1 = "RecGen", cap = "RecoJet_GenTrack";
		auto sig_tb_rg = read_flatten<TH2D>(dataDumpPath+prefix+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_noCorr");
		auto sig_tb_rg_pt = read_flatten<TH2D>(dataDumpPath+prefix+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_pTweighted_noCorr");
		cap1 = "RecRec", cap = "RecoJet_RecoTrack";
		auto sig_tb_rr = read_flatten<TH2D>(dataDumpPath+prefix+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_noCorr");
		auto sig_tb_rr_pt = read_flatten<TH2D>(dataDumpPath+prefix+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_pTweighted_noCorr");
		auto dr_rg = getDr("dr_rg", sig_tb_rg);
		auto dr_rr = getDr("dr_rr", sig_tb_rr);
		auto dr_rg_pt = getDr("dr_gg_pt", sig_tb_rg_pt);
		auto dr_rr_pt = getDr("dr_gr_pt", sig_tb_rr_pt);

		/* testing if tk correction could be pull from signal before bkg subtraction
		cap1 = "GenGen", cap = "GenJet_GenTrack";
		auto fgg = TFile::Open(str_gg_v1);
		auto rs_tb_gg = readSig(fgg, "trueBJet_"+cap+"_noCorr", "trueBJet");
		cap1 = "GenRec", cap = "GenJet_RecoTrack";
		auto fgr = TFile::Open(str_gr_v1);
		auto rs_tb_gr = readSig(fgr, "trueBJet_"+cap+"_noCorr", "trueBJet");
		auto dr_rs_gg = getDr("dr_gg", rs_tb_gg);
		auto dr_rs_gr = getDr("dr_gr", rs_tb_gr);
		auto ratio3 = binary_operation<TH1D>("ratio3", dr_rs_gr, dr_rs_gg, "binomialRatio");
		*/

		auto ratio1 = binary_operation<TH1D>("tkCorr", dr_gr, dr_gg, "binomialRatio");
		auto ratio2 = binary_operation<TH1D>("ratio2", dr_rr, dr_rg, "binomialRatio");
		auto ratio3 = binary_operation<TH1D>("tkCorr_pTweighted", dr_gr_pt, dr_gg_pt, "binomialRatio");
		auto ratio4 = binary_operation<TH1D>("ratio4", dr_rr_pt, dr_rg_pt, "binomialRatio");
		auto dr_cor= binary_operation<TH1D>("dr_corr", dr_rr, ratio1, "ratio");
		auto dr_cor_pt= binary_operation<TH1D>("dr_corr_pt", dr_rr_pt, ratio3, "ratio");
		auto ratio22 = binary_operation<TH1D>("ratio22", dr_cor, dr_rg, "binomialRatio");
		auto ratio44 = binary_operation<TH1D>("ratio44", dr_cor_pt, dr_rg_pt, "binomialRatio");

		auto ratio_js_err = binary_operation<TH1D>("tkError_pTweighted", dr_cor, dr_rg, "binomialRatio");
		auto ratio_py_err = binary_operation<TH1D>("tkError", dr_cor_pt, dr_rg_pt, "binomialRatio");

		showPlot("tkRecoEff.pdf", 0, 1, 0, 0.99, 0.4, 1.2, 2, ratio1, ratio2);
		showPlot("tkRecoEff_pTweighted.pdf", 0, 1, 0, 0.99, 0.4, 1.2, 2, ratio3, ratio4);
		showClosure("tkPeformance.pdf", 0, .99, 0.5, 1.5, "Reco-Corr. ", "Reco-Gen", dr_cor, dr_rg);
		showClosure("tkPeformance_pTweighted.pdf", 0, .99, 0.5, 1.5, "Reco-Corr. ", "Reco-Gen", dr_cor_pt, dr_rg_pt);
		auto wf = new TFile(dataDumpPath+prefix+"ppCSVv1TkCorr.root", "recreate");
		wf->cd();
		dumpHists<TH1D>(ratio1);
		dumpHists<TH1D>(ratio3);
		dumpHists<TH1D>(ratio_js_err);
		dumpHists<TH1D>(ratio_py_err);
		wf->Close();
}
