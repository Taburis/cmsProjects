
#include "./analyzer/utility_pf.h"
#include "../lib/import_inclusive.h"
#include "./analyzer/incl_utility.h"
#include "./analyzer/sub_macros.h"

void pullSignal_pfJets(){
/* ================================ Gen Gen ===============================*/
/* inclusive 
		inclusive_input::getPYTHIA_input("GenJet_GenTrack", inclusive_input::GenGen_pythia_f);
		inclusive_input::getMix("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_f);
		inclusive_input::pullSig("incl_gen_gen_PYTHIA", 1);
		pull1D  ("incl_gen_gen_PYTHIA");
		checkBkg("incl_gen_gen_PYTHIA");
		*/
		/*
		get2DInput_GenGen_PYTHIA_incl();
		//getMix("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_f);
		pullSig( "pf_incl_gen_gen_PYTHIA", 1.5, 2.5, 1);
		pull1D  ("pf_incl_gen_gen_PYTHIA");
		checkBkg("pf_incl_gen_gen_PYTHIA");
//		purify_dr("purified_gen_gen_PYTHIA_test", "gen_gen_PYTHIA_taggedB", "pf_incl_gen_gen_PYTHIA_test", "gen_gen_PYTHIA_tagged_trueB", 0.7);
		dEtaClosure("gen_gen_PYTHIA_trueB", "pf_incl_gen_gen_PYTHIA_test", "true b", "incl");
		drRatio("inclRatio_pf_seagullCheck", "pf_incl_gen_gen_PYTHIA_test", "pf_incl_gen_gen_PYTHIA");
		drRatio("ratio_trueB_over_incl", "gen_gen_PYTHIA_trueB", "pf_incl_gen_gen_PYTHIA_test");
		*/

/* tagged & true B 
		get2DInput_GenGen_PYTHIA_taggedAndTrueB();
		pullSig("gen_gen_PYTHIA_tagged_trueB");
		pull1D("gen_gen_PYTHIA_tagged_trueB");
		checkBkg("gen_gen_PYTHIA_tagged_trueB");
		*/
/* tagged B 
		get2DInput_GenGen_PYTHIA_taggedB();
		pullSig("gen_gen_PYTHIA_taggedB", 1.5, 2.5, 0);
		pull1D("gen_gen_PYTHIA_taggedB");
		checkBkg("gen_gen_PYTHIA_taggedB");

		get2DInput_RecRec_PYTHIA_taggedB();
		pullSig("rec_rec_PYTHIA_taggedB", 1.5, 2.5, 0);
		pull1D("rec_rec_PYTHIA_taggedB");
		checkBkg("rec_rec_PYTHIA_taggedB");

*/
/* true B 
		get2DInput_GenGen_PYTHIA_trueB();
		pullSig("gen_gen_PYTHIA_trueB");
		pull1D("gen_gen_PYTHIA_trueB");
		checkBkg("gen_gen_PYTHIA_trueB");
 */
/* ======================================= Gen - Reco =================================*/
		get2DInput_RecGen_PYTHIA_taggedB();
		pullSig("rec_gen_PYTHIA_taggedB", 1.5, 2.5, 0);
		pull1D("rec_gen_PYTHIA_taggedB");
		checkBkg("rec_gen_PYTHIA_taggedB");
}
