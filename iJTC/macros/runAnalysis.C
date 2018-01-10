
#include "workflow_lib.h"

void runAnalysis(){
		// pull the signal for MC
		//input_raw2D::get2DInput_GenGen(); input_raw2D::pullSig("gen_gen");
		//input_raw2D::get2DInput_RecGen(); input_raw2D::pullSig("rec_gen");
		//input_raw2D::get2DInput_RecRec(); input_raw2D::pullSig("rec_rec");
		//input_raw2D::get2DInput_GenRec(); input_raw2D::pullSig("gen_rec");
		input_raw2D::get2DInput_GenGen_sub0(); input_raw2D::pullSig("gen_gen_sub0");
		//input_raw2D::get2DInput_RecGen_nsub0(); input_raw2D::pullSig("rec_gen_nsub0");
		//input_raw2D::get2DInput_RecGen_sub0();  input_raw2D::pullSig("rec_gen_sub0");

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
//		signal2D::loadFile();
//		ana_fig::closure("track","gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::getDr("gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::pull1D("gen_gen", signal2D::gengen_pb_f);
		//signal1D::checkBkg("gen_gen");
		//signal2D::pull1D("gen_gen", signal2D::gengen_pb_f);
		//signal1D::checkBkg("gen_gen");
		//getting the 4-commuted tabels for jet and track reco validation check 
		//signal2D::drawTableWithRatio("gen_gen", "rec_rec", signal2D::gengen_pb_f, signal2D::recrec_pb_f);
		//signal2D::drawJetShapeRatio("gen_gen", "rec_rec", signal2D::gengen_pb_f, signal2D::recrec_pb_f);
		//signal2D::drawTableWithRatio("rec_gen", "rec_rec", signal2D::recgen_pb_f, signal2D::recrec_pb_f);
		//signal2D::drawJetShapeRatio("rec_gen", "rec_rec", signal2D::recgen_pb_f, signal2D::recrec_pb_f);
		//signal2D::drawJetShapeRatio("gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::drawTableWithRatio("gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::drawJetShapeRatio("gen_gen", "rec_gen", signal2D::gengen_pb_f, signal2D::recgen_pb_f);
		//signal2D::drawTableWithRatio("gen_gen", "rec_gen", signal2D::gengen_pb_f, signal2D::recgen_pb_f);
		//signal2D::drawJetShapeRatio("gen_rec", "rec_rec", signal2D::genrec_pb_f, signal2D::recrec_pb_f);
		//signal2D::drawTableWithRatio("gen_rec", "rec_rec", signal2D::genrec_pb_f, signal2D::recrec_pb_f);
		/*
		signal2D::drawTable("gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		signal2D::drawTable("rec_gen", "rec_rec", signal2D::recgen_pb_f, signal2D::recrec_pb_f);
		signal2D::drawTable("gen_gen", "rec_gen", signal2D::gengen_pb_f, signal2D::recgen_pb_f);
		signal2D::drawTable("gen_rec", "rec_rec", signal2D::genrec_pb_f, signal2D::recrec_pb_f);
		*/
		//signal2D::drawTableWithRatio("gen_gen", "gen_rec", signal2D::gengen_pb_f, signal2D::genrec_pb_f);
		//signal2D::drawTableWithRatio("rec_gen", "rec_rec", signal2D::recgen_pb_f, signal2D::recrec_pb_f);
		//signal2D::drawTableWithRatio("gen_rec", "rec_rec", signal2D::genrec_pb_f, signal2D::recrec_pb_f);
		//signal2D::drawTableWithRatio("rec_gen_nsub0", "rec_rec", signal2D::recgen_pb_nsub0_f, signal2D::recrec_pb_f);
		//signal2D::drawTableWithRatio("rec_gen_sub0", "rec_gen_nsub0", signal2D::recgen_pb_sub0_f, signal2D::recgen_pb_nsub0_f);
}
