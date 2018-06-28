
#ifndef standardSequence_H
#include "standardSequence.C"
#endif
//#include "decontamination.C"
#include "taggerBias.C"
#include "get_jet_spectra.C"
#include "pullSig.C"
//#include "SSLib.h"
void ss_v0(TString name, inputSet &iset, bool dohist = 0){
		auto sig = decontamination(iset, 0.7);
		auto dr = getDr("dr_tmp", sig);
		applyBiasCorrection(dr, 0, dohist);
		auto drCorr = applyTkCorrection(name+"_tkCorr", dr);
		auto drfinal = applyResidual(name, drCorr);

		auto sig2 = decontamination(iset, 0.7, 0);
		auto dr2 = getDr("dr_tmp2", sig2);
		applyBiasCorrection(dr2, 0, dohist);
		auto drCorr2 = applyTkCorrection(name+"_pTweighted_tkCorr", dr2, 0);
		auto drfinal2 = applyResidual(name+"_pTweighted", drCorr2, 0);
		TString label = tpname(iset);
		auto wf = TFile::Open(dataDumpPath+"bJTC_fullCorrected_"+label+".root", "recreate");
		wf->cd();
		dumpHists<TH1D>(drfinal);
		dumpHists<TH1D>(drfinal2);
		wf->Close();
}

void ss(){
		config();

		pullSig(p6gg);
		pullSig(p6gr);
		pullSig(p6rr);
		pullSig(p6rg);
		pullSig(ppData);
//		pullSig(phgg);

//		ss_v0("JS_final_RecRec", p6rr,1);
//		ss_v0("JS_final_Data", ppData);

		//auto js_data = read_flatten<TH1D>(dataDumpPath+"bJTC_fullCorrected_Data.root", "JS_final_Data_pTweighted");
//		auto js_b_data = read_flatten<TH1D>(dataDumpPath+"bJTC_fullCorrected_Data.root", "JS_final_Data_pTweighted");
		//auto js_incl_data = read_flatten<TH1D>(dataDumpPath+"incl_pp_referenceForbJTC.root", "inclJTC_pp_Data_pTweighted");
//		auto js_rr = read_flatten<TH1D>(dataDumpPath+"bJTC_fullCorrected_RecRec.root", "JS_final_RecRec_pTweighted");
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
//		showPlot("drOverlay_data_vs_GenGen.pdf", 0, 1, 0, 0.99, 0.2, 1.8, 2, ratio1, ratio2);
		showClosure("fullClosure_validation.pdf", 0, .99, 0.5, 1.5, "Corr. tag.", "b-jet GenGen", js_b_gg, js_in_gg);
		*/
//		showClosure("drRatio_validation.pdf", 0, .99, 0.5, 1.5, "tagged GenGen", "incl GenGen", js_gg, js_in_gg);
		//showClosure("fullClosure_validation.pdf", 0, .99, 0.5, 1.5, "b-jet", "incl jet", js_data, js_incl);
//		cout<<js_gr[0]->GetName()<<endl;
//		showPlot("dr_final_data.pdf", 0, 0, 0, 0.99, 0.5, -1.5, 1, js_data);

}
