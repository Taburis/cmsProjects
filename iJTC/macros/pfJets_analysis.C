
#include "./analyzer/utility_pf.h"
#include "../lib/import_inclusive.h"
#include "./analyzer/incl_utility.h"
#include "./analyzer/sub_macros.h"

void pfJets_analysis(){

/* ================================ reconstruction validation =================================*/

//		getTaggingBias();
//fit2D();		
//		applyCorrection2D("tagCorrectedB", "purified_gen_gen_PYTHIA_taggedB", "tagBias2DFit");
//		dEtaClosure("tagCorrectedB","gen_gen_PYTHIA_trueB", "tag b", "true b");
//		dPhiClosure("tagCorrectedB","gen_gen_PYTHIA_trueB", "tag b", "true b");
		drRatio("dr_closure_tagCorrB_over_trueB", "tagCorrectedB","gen_gen_PYTHIA_trueB");
		//get2DInput_GenGen_PYTHIA_trueB();
		//pullSig("Kurt");
		//rawRatio2D("Colz2D_Diff_officialCorr_gengen_over_genrec","genrec_incl_off","gengen_incl_off");
		//rawRatio2D("Colz2D_Diff_gengen_over_genrec","genrec_incl","gengen_incl");
//		rawRatio2D("Colz2D_Diff_gengen_taggedB_over_trueB","gen_gen_PYTHIA_taggedB","gen_gen_PYTHIA_trueB");
//		rawRatio2D("Colz2D_Diff_gengen_tagged_trueB_over_trueB","gen_gen_PYTHIA_tagged_trueB","gen_gen_PYTHIA_trueB");
//		get2DInput_GenGen_PYTHIA_trueB_noMixing();
//		pullSig("gengen_trueB_nomixing");
//		get2DInput_GenGen_PYTHIA_trueB_noMixing(0);
//		pullSig("gengen_incl_off");
//		get2DInput_GenRec_PYTHIA_trueB_noMixing();
//		pullSig("genrec_trueB_nomixing");
//		get2DInput_GenRec_PYTHIA_trueB_noMixing(0);
//		pullSig("genrec_incl_off");
//		//dPhiRawClosure("Mine","Kurt", "Mine", "Kurt's");
//		dPhiRawClosure("genrec_trueB_nomixing","gengen_trueB_nomixing", "gen-rec", "gen-gen", 1, 0, 0, 0.1, 0.3);
//		dPhiRawClosure("genrec_incl_off","gengen_incl_off", "gen-rec", "gen-gen");
//		dEtaRawClosure("genrec_incl_off","gengen_incl_off", "gen-rec", "gen-gen");
//		dEtaRawClosure("genrec_trueB_nomixing","gengen_trueB_nomixing", "gen-rec", "gen-gen");
//		drRawRatio("raw_dr_closure_genrec_true_over_gen_gen_true", "genrec_trueB_nomixing","gengen_trueB_nomixing");
//		dPhiRawClosure("genrec_incl","gengen_incl", "gen-rec", "gen-gen");
//		dEtaRawClosure("genrec_incl","gengen_incl", "gen-rec", "gen-gen");
//		drRawRatio("raw_dr_closure_genrec_incl_gen_gen_incl", "genrec_incl","gengen_incl");
//		slicedEtaRawClosure("rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB", "reco-reco", "reco-gen");
//		dPhiRawClosure("rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB", "reco-reco", "reco-gen");
//		dEtaRawClosure("rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB", "reco-reco", "reco-gen");
//		dEtaClosure("rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB", "reco-reco", "reco-gen");
//		drRatio("dr_taggedBclosure_recrec_over_recgen", "rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB");
//		drRawRatio("raw_dr_closure_recrec_tagged_over_rec_gen_tagged", "rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB");
//		dPhiRawClosure("gen_gen_PYTHIA_taggedB","gen_gen_PYTHIA_trueB", "tagged b", "true b");
//		drRawRatio("raw_dr_closure_gengen_taggedB_over_gengen_trueB", "gen_gen_PYTHIA_taggedB","gen_gen_PYTHIA_trueB");
//		drRatio("dr_closure_gengen_taggedB_over_gengen_trueB", "gen_gen_PYTHIA_taggedB","gen_gen_PYTHIA_trueB");
//		drRatio("dr_closure_gengen_tAtB_over_gengen_trueB", "gen_gen_PYTHIA_tagged_trueB","gen_gen_PYTHIA_trueB");
//		drRatio("dr_closure_gengen_purified_over_gengen_trueB", "purified_gen_gen_PYTHIA_taggedB","gen_gen_PYTHIA_trueB");
//		drRatio("dr_closure_gengen_purified_over_gengen_tAtB", "purified_gen_gen_PYTHIA_taggedB","gen_gen_PYTHIA_tagged_trueB");
//taggingBiasCorrection();
//addBiasCorrection();
//		pullCorrection("taggingBias", "gen_gen_PYTHIA_tagged_trueB", "gen_gen_PYTHIA_trueB");
//		batch_sig_subtraction("tagging_bias","gen_gen_PYTHIA_tagged_trueB", "gen_gen_PYTHIA_trueB");
//		dEtaClosure("gen_gen_PYTHIA_tagged_trueB","gen_gen_PYTHIA_trueB", "tagged & true b", "true b", 1, 1);
//		dEtaClosure("rec_gen_PYTHIA_taggedB","gen_gen_PYTHIA_taggedB", "reco jet", "gen jet");
//		dEtaClosure("rec_gen_PYTHIA_taggedB","gen_gen_PYTHIA_taggedB", "reco jet", "gen jet", 1, 1);
//		purify_dr("purified_gen_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_taggedB", "pf_incl_gen_gen_PYTHIA", "gen_gen_PYTHIA_tagged_trueB", 0.699);
//		drRatio("recgen_taggedB_over_gengen_taggedB", "rec_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_taggedB");
//		drRatio("recgen_taggedB_over_gengen_taggedB", "rec_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_taggedB");
		//drRatio("gengen_trueB_over_incl_gengen", "gen_gen_PYTHIA_trueB", "pf_incl_gen_gen_PYTHIA");
//		drRatio("validation_recrec_over_recgen", "rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB");
		/*
		drRatio("validation_recrec_over_recgen", "rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB");
		dEtaClosure("rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB", "reco tracks", "gen tracks");
		dEtaClosure("rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB", "reco tracks", "gen tracks");
		drRatio("validation_gengen_over_recgen", "gen_gen_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB");
		dEtaClosure("rec_gen_PYTHIA_taggedB","gen_gen_PYTHIA_taggedB", "reco jet", "gen jet");
		*/
/* contanmination  
		get2DInput_GenGen_PYTHIA_contan();
		pullSig("gen_gen_PYTHIA_contan");
		pull1D("gen_gen_PYTHIA_contan");
		checkBkg("gen_gen_PYTHIA_contan");
		*/
//		spectraRatio("pf_tt_spectra", "GenJet_GenTrack", "GenJet_GenTrack", gengen_pythia_tagged_trueB_f, gengen_pythia_taggedB_f, 1);
//		spectraRatio("pf_eff_spectra", "GenJet_GenTrack", "GenJet_GenTrack", gengen_pythia_tagged_trueB_f, gengen_pythia_trueB_f, 1);
	//	drawRecoCheck("gen_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_tagged_trueB");
	//	drawRecoCheck("gen_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_trueB");
	//	drawRecoCheck("gen_gen_PYTHIA_tagged_trueB", "gen_gen_PYTHIA_trueB");
	//	drawRecoCheck("gen_gen_PYTHIA_taggedB", "incl_gen_gen_PYTHIA");
//		purify_dr("purified_gen_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_taggedB", "pf_incl_gen_gen_PYTHIA", "gen_gen_PYTHIA_tagged_trueB", 0.699);
//		purify_dr("purified_gen_gen_PYTHIA_taggedB2", "gen_gen_PYTHIA_taggedB", "calo_incl_gen_gen_PYTHIA", "gen_gen_PYTHIA_tagged_trueB", 0.7);
//		drRatio("remove_incl", "purified_gen_gen_PYTHIA_taggedB2", "gen_gen_PYTHIA_tagged_trueB");
		//drRatio("inclRatio_pf_over_calo", "pf_incl_gen_gen_PYTHIA_test", "incl_gen_gen_PYTHIA");
		//drRatio("inclRatio_pf_over_calo", "pf_incl_gen_gen_PYTHIA_test", "calo_incl_gen_gen_PYTHIA");
//		drRatio("inclRatio_pf_test_over_pf", "pf_incl_gen_gen_PYTHIA_test", "pf_incl_gen_gen_PYTHIA");
//		drRatio("inclRatio_pf_over_calo", "pf_incl_gen_gen_PYTHIA", "incl_gen_gen_PYTHIA");
//		drRatio("inclRatio_mine_over_offical", "incl_gen_gen_PYTHIA", "calo_incl_gen_gen_PYTHIA");
//		drRatio("contan_over_incl", "gen_gen_PYTHIA_contan", "incl_gen_gen_PYTHIA");
//		drRatio("dr_ratio_gen_gen_taggedB_over_tagged_trueB", "gen_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_tagged_trueB");
//		drRatio("dr_ratio_gen_gen_taggedB_over_tagged_trueB", "gen_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_tagged_trueB");
//		dEtaClosure("gen_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_tagged_trueB", "tagged b", "tagged & true b");
//		dEtaClosure("gen_gen_PYTHIA_tagged_trueB", "gen_gen_PYTHIA_trueB", "tagged & true b", "true b");
//		dEtaClosure("purified_gen_gen_PYTHIA_tagged_trueB", "gen_gen_PYTHIA_trueB", "tagged & true b", "true b");
//		dEtaClosure("pf_incl_gen_gen_PYTHIA", "incl_gen_gen_PYTHIA", "pf_incl", "calo_incl");
}
