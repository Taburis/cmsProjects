
#include "workflow_lib.h"

void runAnalysis(){
		// pull the signal for MC
		//input_raw2D::get2DInput_GenGen(); input_raw2D::pullSig("gen_gen");
		//input_raw2D::get2DInput_RecGen(); input_raw2D::pullSig("rec_gen");
		//input_raw2D::get2DInput_RecRec(); input_raw2D::pullSig("rec_rec");
		//input_raw2D::get2DInput_GenRec(); input_raw2D::pullSig("gen_rec");
		//input_raw2D::get2DInput_GenGen_sub0(); input_raw2D::pullSig("gen_gen_sub0_mix");
		
		//input_raw2D::get2DPythiaInput_GenGen(); input_raw2D::pullSig("gen_gen_pythia");
		//input_raw2D::get2DInput_GenGen_sub0_trueB();  input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0_trueB");
		//input_raw2D::get2DInput_GenGen_sub0_tagged_trueB();  input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0_tagged_trueB");
		//input_raw2D::get2DInput_GenGen_sub0();  input_raw2D::getMix_GenGen(); input_raw2D::pullSig("gen_gen_sub0");
		//input_raw2D::get2DInput_RecGen_nsub0(); input_raw2D::pullSig("rec_gen_nsub0");
		//input_raw2D::get2DInput_RecGen_sub0();  input_raw2D::pullSig("rec_gen_sub0");

//pull signal from inclusive jet
		//inclusive_input::getH2("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_sub0_f);
		//inclusive_input::getMix("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_f);
		//inclusive_input::pullSig("incl_gen_gen_sub0");
	

// get 1D histogram for checking 
		signal2D::loadFile();
		//TFile *f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/gen_gen_sub0_mix_JTCSignal.root");
		//TFile *f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/gen_gen_pythia_JTCSignal.root");
		//signal2D::pull1D("gen_gen_pythia", f);
		//signal2D::pull1D("gen_gen_sub0_mix", f);
		//signal2D::pull1D("gen_gen_sub0_mix", f);
		//signal2D::pull1D("gen_rec", signal2D::genrec_pb_f);
		//signal2D::pull1D("gen_gen_sub0", signal2D::gengen_pb_sub0_f);
		//signal2D::pull1D("incl_gen_gen_sub0", inclusive_input::gengen_pb_sub0_f);
		signal2D::drawStackJSDiff("gen_gen_sub0","gen_gen_pythia", "gengen_sub0", 1);
		signal2D::drawStackJSDiff("gen_gen","gen_gen_pythia", "gengen", 1);


// bkg checking 
		//signal1D::checkBkg("gen_gen_pythia");
		//signal1D::checkBkg("gen_gen_sub0_mix");
		//signal1D::checkSide("gen_gen_sub0");
		//signal1D::checkBkg("incl_gen_gen_sub0");
		//signal1D::checkSide("incl_gen_gen_sub0");
		//inclusive_input::drawRatio_sub0("gen_gen_sub0_trueB", "incl_gen_gen_sub0", "gengen_sub0_trueB");
		//inclusive_input::drawRatio_sub0("gen_gen_sub0", "incl_gen_gen_sub0", "gengen_sub0");
		//inclusive_input::drawRatio_sub0("gen_gen_sub0_trueB", "incl_gen_gen_sub0", "gengen_sub0_trueB_over_incl");
		//inclusive_input::drawRatio_sub0("gen_gen_sub0", "incl_gen_gen_sub0", "gengen_sub0_tagged_over_incl");
		//inclusive_input::drawRatio_sub0("gen_gen_sub0_trueB", "gen_gen_sub0_tagged_trueB", "gengen_sub0_trueB_over_tagged_trueB", 1);
		//inclusive_input::drawRatio_sub0("gen_gen_sub0", "gen_gen_sub0_tagged_trueB", "gengen_sub0_tagged_trueBratio");
	/*
		input_raw2D::get2DInput_GenGen_sub0();
		inclusive_input::getH2("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_sub0_f);
						*/



		//inclusive_input::getH2("RecoJet_RecoTrack", inclusive_input::RecRec_MC_pb_f);  inclusive_input::pullSig("incl_rec_rec");
		//input_raw2D::showSpectra("GenJet_RecoTrack", "RecoJet_RecoTrack",input_raw2D::genrec_pb_f, input_raw2D::recrec_pb_f);
		//input_raw2D::showSpectra("GenJet_GenTrack", "RecoJet_GenTrack",input_raw2D::gengen_pb_f, input_raw2D::recgen_pb_f);
		//input_raw2D::showSpectra("RecoJet_GenTrack", "RecoJet_RecoTrack",input_raw2D::recgen_pb_f, input_raw2D::recrec_pb_f);
// TH2D sector: deal within all the th2 histograms pulling from raw inputs
		//input_raw2D::get2DInput_GenGen();
		//inclusive_input::test(input_raw2D::raw_sig);
		/*
		input_raw2D::get2DInput_GenGen_sub0();
		inclusive_input::getH2("GenJet_GenTrack", inclusive_input::GenGen_MC_pb_sub0_f);
		inclusive_input::drawRatio_sub0(input_raw2D::raw_sig_pTweighted,input_raw2D::mixing, inclusive_input::raw_sig_pTweighted,
						inclusive_input::mixing);
						*/
//		ana_fig::closure("track","gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::getDr("gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::pull1D("gen_gen", signal2D::gengen_pb_f);
		//signal1D::checkBkg("gen_gen");
		//getting the 4-commuted tabels for jet and track reco validation check 
}
