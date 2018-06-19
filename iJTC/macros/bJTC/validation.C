
#include "standardSequence.C"
#include "syst_utility.h"

void validation(){
		config();
		bool isnumber = 1;
		TString type = isnumber ? "": "_pTweighted";

// decontamination
//		auto sig_tg_rr = read_flatten<TH2D>(dataDumpPath+"taggedBJet_RecRec_JTCSignal.root", "signal_taggedBJet_RecoJet_RecoTrack"+type+"_noCorr");
		auto rs = decontamination(p6rr, 0.7, isnumber);
		auto rs_deta = projX("rs_deta", rs);
		auto rs_dr = getDr("rs_dr", rs);
		auto rs_dr_err = getDrErr("rs_dr_err", rs_dr, rs_deta, 1);
		//ref for decontamination
		auto sig_tt_rr = read_flatten<TH2D>(dataDumpPath+"taggedTrueBJet_RecRec_JTCSignal.root", "signal_taggedTrueBJet_RecoJet_RecoTrack"+type+"_noCorr");
		auto js_tt_rr = getDr("js_tt_rr", sig_tt_rr);
		auto sig_tt_rr_deta = projX("sig_tt_rr_deta", sig_tt_rr);
		auto js_tt_rr_err = getDrErr("js_tt_rr_err", js_tt_rr, sig_tt_rr_deta, 1);

		add_frac_error(rs_dr_err, 0.1);// decontamination non-closure
		showClosure_syst("validation/dr_afterDecont_RecRec.pdf", 0, .99, 0.5, 1.5, "decont.ite.(RecRec)", "t&t (RecRec)", rs_dr, rs_dr_err, js_tt_rr, js_tt_rr_err);

		auto sig_dr = copy<TH1D>("sig_dr", rs_dr);
		applyBiasCorrection(sig_dr, isnumber);
		auto sig_dr_err = copy<TH1D>("sig_dr_err", sig_dr);
		sign_err(sig_dr_err, rs_dr_err);

		//ref for tagger bias correciton
		auto sig_tb_rr = read_flatten<TH2D>(dataDumpPath+"trueBJet_RecRec_JTCSignal.root", "signal_trueBJet_RecoJet_RecoTrack"+type+"_noCorr");
		auto js_tb_rr = getDr("js_tt_rr", sig_tb_rr);
		auto sig_tb_rr_deta = projX("sig_tb_rr_deta", sig_tb_rr);
		auto js_tb_rr_err = getDrErr("js_tb_rr_err", js_tb_rr, sig_tb_rr_deta, 1);
		showClosure_syst("validation/dr_afterBias_RecRec.pdf", 0, .99, 0.5, 1.5, "bias.ite.(RecRec)", "true b (RecRec)", sig_dr, sig_dr_err, js_tb_rr, js_tb_rr_err);

		auto drCorr = applyTkCorrection("drtkCorr", sig_dr, isnumber);
		auto drCorrErr = applyTkCorrection("drtkCorr", sig_dr_err, isnumber);
		auto tkErr = read_flatten<TH1D>(dataDumpPath+"ppCSVv1TkCorr.root", "tkError"+type);
		add_frac_error(drCorrErr , tkErr);

		//ref for tracking correciton
		auto sig_tb_rg = read_flatten<TH2D>(dataDumpPath+"trueBJet_RecGen_JTCSignal.root", "signal_trueBJet_RecoJet_GenTrack"+type+"_noCorr");
		auto js_tb_rg = getDr("js_tb_rg", sig_tb_rg);
		auto sig_tb_rg_deta = projX("sig_tb_rg_deta", sig_tb_rg);
		auto js_tb_rg_err = getDrErr("js_tb_rg_err", js_tb_rg, sig_tb_rg_deta, 1);
		showClosure_syst("validation/dr_afterTk_RecRec.pdf", 0, .99, 0.5, 1.5, "tk.ite.(RecRec)", "true b (RecGen)", drCorr, drCorrErr, js_tb_rg, js_tb_rg_err);
		
		auto drfinal = applyResidual("dr_final", drCorr, isnumber);
		auto drfinal_err = copy<TH1D>("dr_final_err", drfinal);
		sign_err(drfinal_err, drCorrErr);

		//ref for residual correciton
		auto sig_tb_gg = read_flatten<TH2D>(dataDumpPath+"trueBJet_GenGen_JTCSignal.root", "signal_trueBJet_GenJet_GenTrack"+type+"_noCorr");
		auto js_tb_gg = getDr("js_tb_rg", sig_tb_gg);
		auto sig_tb_gg_deta = projX("sig_tb_gg_deta", sig_tb_gg);
		auto js_tb_gg_err = getDrErr("js_tb_gg_err", js_tb_gg, sig_tb_gg_deta, 1);
		showClosure_syst("validation/dr_full_RecRec.pdf", 0, .99, 0.5, 1.5, "full ite.(RecRec)", "true b (GenGen)", drfinal, drfinal_err, js_tb_gg, js_tb_gg_err);
//		showClosure("validation/dr_afterDecont_RecRec.pdf", 0, .99, 0.5, 1.5, "decont.(RecRec)", "t&t (RecRec)", rs_dr, js_tt_rr);

/*
*/
}
