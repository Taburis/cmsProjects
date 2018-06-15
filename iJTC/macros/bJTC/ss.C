
#ifndef standardSequence_H
#include "standardSequence.C"
#endif
//#include "decontamination.C"
#include "taggerBias.C"
#include "get_jet_spectra.C"
#include "pullSig.C"
//#include "SSLib.h"
TH1D** ss_v0(TString name, inputSet &iset){
		auto dr = decontamination(iset, 0.7);
		applyBiasCorrection(dr);
		auto drCorr = applyTkCorrection(name+"_tkCorr", dr);
		auto drfinal = applyResidual(name, dr);
		if(iset.isRecoTk) return drCorr;
		if(iset.isRecoJet) return drfinal;
		else return dr;
}
void ss(){
		config();

//		pullSig(p6rg);
//		pullSig(p6rr);
//		pullSig(p6rg);
//		auto js_gr = ss_v0("JS_final_RecRec", p6rr);
TString name = "test";
		auto dr = decontamination(p6rr, 0.7);
		applyBiasCorrection(dr);
		auto drCorr = applyTkCorrection(name+"_tkCorr", dr);
		auto drfinal = applyResidual(name, drCorr);
		auto js_gr= drfinal;

		TString tyn = tpname(p6rr);
		TString cap1 = "GenGen", cap = "GenJet_GenTrack";
		auto sig_tb = read_flatten<TH2D>(dataDumpPath+"trueBJet_"+cap1+"_JTCSignal.root", "signal_trueBJet_"+cap+"_noCorr");
		auto js_tb = getDr("js_gg_tb", sig_tb);
//		
		showClosure("fullClosure_"+tyn+"_over_"+cap1+".pdf", 0, .99, 0.5, 1.5, "corr.", "true b", js_gr, js_tb);
		showClosure("test.pdf", 0, .99, 0.5, 1.5, "corr.", "true b", drCorr, drfinal);
//decontamination2(p6rr);
/*
		TString cap1 = "RecRec", cap = "RecoJet_RecoTrack";
		auto sig_in = read_flatten<TH2D>(dataDumpPath+"inclJet_"+cap1+"_JTCSignal.root", "signal_inclJet_"+cap+"_noCorr");
		auto sig_co = read_flatten<TH2D>(dataDumpPath+"contJet_"+cap1+"_JTCSignal.root", "signal_contJet_"+cap+"_noCorr");
		auto dr_in = getDr("dr_incl", sig_in);
		auto dr_co = getDr("dr_cont", sig_co);
		showClosure("contamination_"+cap1+".pdf", 0, .99, 0.5, 1.5, "corr.", "true b", dr_in, dr_co);
		*/
}
