
#include "./analyzer/utility_pf.h"
#include "../lib/import_inclusive.h"
#include "./analyzer/incl_utility.h"
#include "./analyzer/sub_macros.h"

void pfJets_analysis(){

/* ================================ reconstruction validation =================================*/
		dEtaClosure("rec_gen_PYTHIA_taggedB","gen_gen_PYTHIA_taggedB", "reco jet", "gen jet");
		dEtaClosure("rec_gen_PYTHIA_taggedB","gen_gen_PYTHIA_taggedB", "reco jet", "gen jet", 1, 1);
//		drRatio("recgen_taggedB_over_gengen_taggedB", "rec_gen_PYTHIA_taggedB", "gen_gen_PYTHIA_taggedB");
		/*
		drRatio("validation_recrec_over_recgen", "rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB");
		dEtaClosure("rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB", "reco tracks", "gen tracks");
		drRatio("validation_recrec_over_recgen", "rec_rec_PYTHIA_taggedB","rec_gen_PYTHIA_taggedB");
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
