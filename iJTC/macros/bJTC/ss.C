
#ifndef standardSequence_H
#include "standardSequence.C"
#endif
//#include "decontamination.C"
#include "taggerBias.C"
#include "get_jet_spectra.C"
#include "pullSig.C"
//#include "SSLib.h"
void ss_v0(TString name, inputSet &iset, TString capname, bool dohist = 0){
		float purity = 0.48;
		auto sig = decontamination(iset, purity);
		auto dr = getDr("dr_tmp", sig);
		applyBiasCorrection(dr, 0, dohist);
		auto drCorr = applyTkCorrection(name+"_tkCorr", dr);
		auto drfinal = applyResidual(name, drCorr);

		auto sig2 = decontamination(iset, purity, 0);
		auto dr2 = getDr("dr_tmp2", sig2);
		applyBiasCorrection(dr2, 0, dohist);
		auto drCorr2 = applyTkCorrection(name+"_pTweighted_tkCorr", dr2, 0);
		auto drfinal2 = applyResidual(name+"_pTweighted", drCorr2, 0);
		TString label = tpname(iset);
		auto wf = TFile::Open(dataDumpPath+capname+"_"+label+".root", "recreate");
		wf->cd();
		dumpHists<TH1D>(drfinal);
		dumpHists<TH1D>(drfinal2);
		wf->Close();
}

void ss(){
		config();
		pullSig(p6rrCSV85, "_CSV85");
		//ss_v0("JS_final_Data_CSV85", p6rrCSV85, "bJTC_fullCorrected_85");
//		ss_v0("JS_final_Data", ppData, "bJTC_fullCorrected");
//		ss_v0("JS_final_Data", ppData, "bJTC_fullCorrected");
/*
		pullSig(p6gg);
		pullSig(p6gr);
		pullSig(p6rg);
		pullSig(ppData);
		pullSig(phgg);
		pullSig(p6rr);

//		ss_v0("JS_final_RecRec", p6rr,1);

		readInput(ppData);
		readInput(p6gg);
		auto sig_py = doMixingCorr("mix_corrected_py", ppData.tg.sig_raw , ppData.tg.mixing_raw, 0);
		auto sig_js = doMixingCorr("mix_corrected_js", ppData.tg.sig_pTweighted_raw , ppData.tg.mixing_raw, 0);
		auto sig_gg_py = doMixingCorr("mix_corrected_py", p6gg.tg.sig_raw ,            p6gg.tg.mixing_raw, 0);
		auto sig_gg_js = doMixingCorr("mix_corrected_js", p6gg.tg.sig_pTweighted_raw , p6gg.tg.mixing_raw, 0);
		auto dr_py_gg = getDr("py_gg", p6gg.tb.sig_raw);
		auto dr_js_gg = getDr("js_gg", p6gg.tb.sig_pTweighted_raw);
		auto dr_py = getDr("py", ppData.tg.sig_raw);
		auto dr_js = getDr("js", ppData.tg.sig_pTweighted_raw);
		auto ratio1 = binary_operation<TH1D>("ratio1", dr_py, dr_js, "ratio");
		auto ratio2 = binary_operation<TH1D>("ratio2", dr_py_gg, dr_js_gg, "ratio");
		showPlot("pt_distribution_check.gif", 0, 1, 0, 0.99, 0.0, -1., 2, ratio1, ratio2);
		*/
		/*
		auto sig_b_gg = read_flatten<TH2D>(dataDumpPath+"trueBJet_GenRec_JTCSignal.root", "signal_trueBJet_GenJet_RecoTrack_pTweighted_noCorr");
		auto sig_in_gg = read_flatten<TH2D>(dataDumpPath+"inclJet_GenRec_JTCSignal.root", "signal_inclJet_GenJet_RecoTrack_pTweighted_noCorr");
		//auto sig_b_gg = read_flatten<TH2D>(dataDumpPath+"trueBJet_GenGen_JTCSignal.root", "signal_trueBJet_GenJet_GenTrack_pTweighted_noCorr");
		//auto sig_in_gg = read_flatten<TH2D>(dataDumpPath+"inclJet_GenGen_JTCSignal.root", "signal_inclJet_GenJet_GenTrack_pTweighted_noCorr");
		auto js_in_gg = getDr("js_in_gg", sig_in_gg);
		auto js_b_gg  = getDr("js_b_gg", sig_b_gg);
//		auto ratio1 = binary_operation<TH1D>("ratio1", js_b_gg, js_in_gg, "ratio");
//		auto ratio2 = binary_operation<TH1D>("ratio2", js_b_data, js_incl_data, "ratio");
//		auto ratio2 = binary_operation<TH1D>("ratio2", js_rr, js_incl, "ratio");
		showClosure("fullClosure_validation.pdf", 0, .99, 0.5, 1.5, "Corr. tag.", "b-jet GenGen", js_b_gg, js_in_gg);
		*/
//		showClosure("drRatio_validation.pdf", 0, .99, 0.5, 1.5, "tagged GenGen", "incl GenGen", js_gg, js_in_gg);
		//showClosure("fullClosure_validation.pdf", 0, .99, 0.5, 1.5, "b-jet", "incl jet", js_data, js_incl);
//		cout<<js_gr[0]->GetName()<<endl;
//		showPlot("dr_final_data.pdf", 0, 0, 0, 0.99, 0.5, -1.5, 1, js_data);

}
