
#include "workflow_lib.h"

void runAnalysis(){
		signal2D::loadFile();
		// pull the signal for MC
		//input_raw2D::get2DInput_GenGen(); input_raw2D::pullSig("gen_gen");
		//input_raw2D::get2DInput_RecGen(); input_raw2D::pullSig("rec_gen");
		//input_raw2D::get2DInput_RecRec(); input_raw2D::pullSig("rec_rec");
		//input_raw2D::get2DInput_GenRec(); input_raw2D::pullSig("gen_rec");
		//input_raw2D::get2DInput_GenRec(); input_raw2D::pullSig("gen_rec");
/* P+H gen-gen sube0 trueB*/
		//input_raw2D::get2DInput_GenGen_sub0_trueB(); input_raw2D::getMix_GenGen_nsube0(); input_raw2D::pullSig("gen_gen_sub0_trueB");
		//input_raw2D::get2DInput_GenGen_sub0_trueB(); input_raw2D::getMix_GenGen_nsube0(); input_raw2D::pullSig("gen_gen_sub0_trueB", 1.7, 2.7);
		//input_raw2D::get2DInput_GenGen_sub0_trueB(); input_raw2D::pullSig("gen_gen_sub0_trueB");
		//signal2D::pull1D("gen_gen_sub0_trueB", 1);
		//signal1D::checkBkg("gen_gen_sub0_trueB");
		
		//input_raw2D::get2DInput_GenGen_sub0_trueB(); input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0_trueB_fineBin");
		//signal2D::pull1D("gen_gen_sub0_trueB_fineBin", 0);
/* P+H gen-gen sube0 tagged trueB*/
		//input_raw2D::get2DInput_GenGen_sub0_tagged_trueB(); input_raw2D::getMix_GenGen_nsube0(); input_raw2D::pullSig("gen_gen_sub0_tagged_trueB");
		//input_raw2D::get2DInput_GenGen_sub0_tagged_trueB(); input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0_tagged_trueB");
		//input_raw2D::get2DInput_GenGen_sub0_tagged_trueB(); input_raw2D::pullSig("gen_gen_sub0_tagged_trueB");
		//signal2D::pull1D("gen_gen_sub0_tagged_trueB",1);
		//signal1D::checkBkg("gen_gen_sub0_tagged_trueB");
		
/* P+H gen-gen sube0 */
		//input_raw2D::get2DInput_GenGen_sub0(); input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0");
		//input_raw2D::get2DInput_GenGen_sub0(); input_raw2D::getMix_GenGen_nsube0(); input_raw2D::pullSig("gen_gen_sub0");
		//input_raw2D::get2DInput_GenGen_sub0(); input_raw2D::pullSig("gen_gen_sub0", 1.7, 2.7);
		//signal2D::pull1D("gen_gen_sub0",1);
		//signal1D::checkBkg("gen_gen_sub0");

		//input_raw2D::get2DInput_GenGen_sub0(); input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0_fineBin");
		//signal2D::pull1D("gen_gen_sub0_fineBin",0);
		//signal1D::checkBkg("gen_gen_sub0_fineBin");

/* P+H gen-gen nsube0 */
		input_raw2D::get2DInput_GenGen_nsube0(); input_raw2D::pullSig("gen_gen_nsub0", 1.3, 2.);
		//signal2D::pull1D("gen_gen_nsub0");

/* P+H rec-gen sube0 */
		//input_raw2D::get2DInput_RecGen_sub0(); input_raw2D::pullSig("rec_gen_sub0");
		//signal2D::pull1D("rec_gen_sub0");
		//signal1D::checkBkg("rec_gen_sub0");

/* P+H rec-gen nsube0 */
		//input_raw2D::get2DInput_RecGen_nsube0(); input_raw2D::pullSig("rec_gen_nsub0");
		//signal2D::pull1D("rec_gen_nsub0");
		//signal1D::checkBkg("rec_gen_nsub0");
/* data */
		//input_raw2D::get2DInput_Data(); input_raw2D::pullSig("data_pb");
		//signal2D::pull1D("data_pb");
		//signal1D::checkBkg("data_pb");
		//signal2D::loadFile();
		//

//pull signal from inclusive jet
		//inclusive_input::getMix("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_f);
		//inclusive_input::getMix("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_nsube0_f);
		//inclusive_input::getH2("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_sube0_f);
		//inclusive_input::pullSig("incl_gen_gen_sub0", 1.5, 2.5, 1);
		//signal2D::pull1D("incl_gen_gen_sub0", 1);
		//signal1D::checkBkg("incl_gen_gen_sub0");
	

/* 4 pumutation stacks in particle yield */
	//tracking 
		//ana_fig::drawRecoCheck("gen_gen_sub0", "gen_rec");
	//tagging 
		ana_fig::drawRecoCheck("gen_gen_sub0", "incl_gen_gen_sub0");
		//ana_fig::drawRecoCheck("gen_gen_sub0_trueB", "incl_gen_gen_sub0");
		//ana_fig::drawRecoCheck("rec_gen_sub0", "gen_gen_sub0_trueB");
		//ana_fig::drawRecoCheck("rec_gen_sub0", "gen_gen_sub0");
		//ana_fig::drawRecoCheck("rec_gen_sub0", "gen_gen_sub0_trueB");
		//ana_fig::drawRecoCheck("gen_gen_sub0", "gen_gen_sub0_tagged_trueB");
		//ana_fig::drawRecoCheck("gen_gen_sub0_tagged_trueB", "gen_gen_sub0_trueB");
		ana_fig::drawRecoCheck("gen_gen_sub0", "gen_gen_sub0_trueB");
		//ana_fig::drawRecoCheck("gen_rec", "gen_gen");
		//ana_fig::drawTrackRecoCheck("rec_jet_track_test","rec_gen", "rec_rec");
		//signal2D::drawStackJSDiff("rec_gen","rec_gen_pythia", "rec_gen_stackClosure", 1);
		//signal2D::drawStackJSDiff("gen_gen","gen_gen_pythia", "gen_gen_stackClosure", 1);
		//signal2D::drawStackJSDiff("data_pb","gen_gen_pythia", "pre", 1);
		//signal2D::drawStackJSDiff("data_pb","rec_gen_pythia", "pre2", 1);

	// pythia:	
		//input_raw2D::get2DPythiaInput_GenGen(); input_raw2D::pullSig("gen_gen_pythia");
		//input_raw2D::get2DPythiaInput_RecGen(); input_raw2D::pullSig("rec_gen_pythia");
		//input_raw2D::get2DPythiaInput_RecRec(); input_raw2D::pullSig("rec_rec_pythia");
		//input_raw2D::get2DPythiaInput_GenRec(); input_raw2D::pullSig("gen_rec_pythia");
		//
		//input_raw2D::get2DInput_GenGen_sub0_trueB();  input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0_trueB");
		//input_raw2D::get2DInput_GenGen_sub0_tagged_trueB();  input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0_tagged_trueB");
		//input_raw2D::get2DInput_GenGen_sub0();  input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0");
		//input_raw2D::get2DInput_RecGen_nsub0(); input_raw2D::pullSig("rec_gen_nsub0");
		//input_raw2D::get2DInput_RecGen_sub0();  input_raw2D::pullSig("rec_gen_sub0");


// get 1D histogram for checking 
		//TFile *f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/gen_gen_sub0_mix_JTCSignal.root");
		//TFile *f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/gen_gen_pythia_JTCSignal.root");
		//signal2D::pull1D("gen_gen_pythia", f);
		//signal2D::pull1D("rec_rec_pythia");
		//signal2D::pull1D("rec_gen_pythia");
		//signal2D::pull1D("gen_gen_sub0_mix", f);
		//signal2D::pull1D("gen_gen_sub0_mix", f);
		//signal2D::pull1D("gen_rec", signal2D::genrec_pb_f);
		//signal2D::pull1D("gen_gen_sub0", signal2D::gengen_pb_sub0_f);
		//signal2D::pull1D("incl_gen_gen_sub0", inclusive_input::gengen_pb_sub0_f);
		//signal2D::drawStackJSDiff("rec_rec","gen_gen_pythia", "recrec_test", 1);
		//signal2D::drawStackJSDiff("gen_gen","gen_gen_pythia", "gengen_test", 1);
		//signal2D::drawStackJSDiff("gen_gen_sub0","gen_gen_pythia", "gengen_sub0", 1);
		//signal2D::drawStackJSDiff("gen_gen","gen_gen_pythia", "gengen", 1);
//		signal2D::drawJSratio("rec_gen","gen_gen", "recgen_over_gengen", 0);
		//signal2D::drawJSratio("gen_gen_pythia","rec_gen_pythia", "gengen_over_recgen", 1);


// bkg checking 
		//signal1D::checkBkg("rec_gen_pythia", 1);
		//signal1D::checkBkg("rec_rec_pythia", 1);
		//signal1D::checkBkg("gen_rec_pythia", 1);
		//signal1D::checkBkg("gen_gen_pythia", 1);
		//signal1D::checkBkg("gen_gen_sub0_mix");
		//signal1D::checkSide("gen_gen_sub0");
		//signal1D::checkBkg("incl_gen_gen_sub0");
		//signal1D::checkSide("incl_gen_gen_sub0");
		//input_raw2D::showSpectra("sube0_tagged_trueB_over_trueB","GenJet_GenTrack", "GenJet_GenTrack", input_raw2D::gengen_pb_sube0_tagged_trueB_f,input_raw2D::gengen_pb_sube0_trueB_f);
		//input_raw2D::showSpectra("sube0_taggedB_over_trueB","GenJet_GenTrack", "GenJet_GenTrack", input_raw2D::gengen_pb_sube0_f,input_raw2D::gengen_pb_sube0_trueB_f);
		//input_raw2D::showSpectra("sube0_tagged_trueB_over_tagged_B","GenJet_GenTrack", "GenJet_GenTrack", input_raw2D::gengen_pb_sube0_tagged_trueB_f,input_raw2D::gengen_pb_sube0_f);
		//input_raw2D::showSpectra("sube0_trueB_over_tagged_trueB","GenJet_GenTrack", "GenJet_GenTrack",input_raw2D::gengen_pb_sube0_trueB_f, input_raw2D::gengen_pb_sub0_tagged_trueB_f);
		//input_raw2D::showSpectra("sube0_taggedB_over_tagged_trueB","GenJet_GenTrack", "GenJet_GenTrack",input_raw2D::gengen_pb_sube0_f, input_raw2D::gengen_pb_sub0_tagged_trueB_f);
		//inclusive_input::showSpectra("sube0_taggedB_over_incl","GenJet_GenTrack", "GenJet_GenTrack",input_raw2D::gengen_pb_sube0_f, inclusive_input::GenGen_MC_pb_sube0_f);
		//inclusive_input::drawRatio_sub0("gen_gen_sub0", "incl_gen_gen_sub0", "gengen_sub0");
		//inclusive_input::drawRatio_sub0("gen_gen_sub0_trueB", "incl_gen_gen_sub0", "gengen_sub0_trueB_over_incl");
		//inclusive_input::drawRatio_sub0("gen_gen_sub0", "incl_gen_gen_sub0", "gengen_sub0_tagged_over_incl");
		//inclusive_input::drawRatio_sub0("gen_gen_sub0_trueB", "gen_gen_sub0_tagged_trueB", "gengen_sub0_trueB_over_tagged_trueB", 1);
		//inclusive_input::drawRatio_sub0("gen_gen_sub0", "gen_gen_sub0_tagged_trueB", "gengen_sub0_tagged_trueBratio");
	/*
		input_raw2D::get2DInput_GenGen_sub0();
		inclusive_input::getH2("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_sube0_f);
						*/



		//inclusive_input::getH2("RecoJet_RecoTrack", inclusive_input::RecRec_MC_pb_f);  inclusive_input::pullSig("incl_rec_rec");
		//input_raw2D::showSpectra("tagged_trueB_over_trueB","GenJet_GenTrack", "GenJet_GenTrack", input_raw2D::gengen_pb_sube0_tagged_trueB_f,input_raw2D::gengen_pb_sube0_trueB_f);
		//input_raw2D::showSpectra("trueB_over_tagged_trueB","GenJet_GenTrack", "GenJet_GenTrack",input_raw2D::gengen_pb_sube0_trueB_f, input_raw2D::gengen_pb_sube0_tagged_trueB_f);
		//input_raw2D::showSpectra("RecoJet_GenTrack", "RecoJet_RecoTrack",input_raw2D::recgen_pb_f, input_raw2D::recrec_pb_f);
// TH2D sector: deal within all the th2 histograms pulling from raw inputs
		//input_raw2D::get2DInput_GenGen();
		//inclusive_input::test(input_raw2D::raw_sig);
		/*
		input_raw2D::get2DInput_GenGen_sub0();
		inclusive_input::getH2("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_sube0_f);
		inclusive_input::drawRatio_sub0(input_raw2D::raw_sig_pTweighted,input_raw2D::mixing, inclusive_input::raw_sig_pTweighted,
						inclusive_input::mixing);
						*/
//		ana_fig::closure("track","gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::getDr("gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::pull1D("gen_gen", signal2D::gengen_pb_f);
		//signal1D::checkBkg("gen_gen");
		//getting the 4-commuted tabels for jet and track reco validation check 
}
