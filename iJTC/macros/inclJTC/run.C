
#include "utility.h" 
#include "submacros.h" 
void run(){
//		get2DInput_RecGen_PYTHIA();
//		pullSig("rec_gen_PYTHIA", 1.5, 2.5, 1);
//		checkBkg("rec_gen_PYTHIA", 0);
//		get2DInput_GenGen_PYTHIA();
//		pullSig("gen_gen_PYTHIA", 1.5, 2.5, 1);
//		checkBkg("gen_gen_PYTHIA", 0);
		/*

		batch_sig_subtraction("JFF_PYTHIA", "rec_gen_PYTHIA", "gen_gen_PYTHIA");
		show_2D_sig("check2D_JFFCorr", "JFF_PYTHIA");
		show_2D_sig("check2D_rec_gen_sig", "rec_gen_PYTHIA", 1);
		show_2D_sig("check2D_gen_gen_sig", "gen_gen_PYTHIA", 1);
		show_2D_raw_sig("check2D_rec_gen_rawSig", "rec_gen_PYTHIA", 1);
		show_2D_raw_sig("check2D_gen_gen_rawSig", "gen_gen_PYTHIA", 1);
		drawDrDist("jff_drProj", "JFF_PYTHIA");
		get2DInput_pp_data();
		pullSig("Data_pp", 1.5, 2.5, 0);
		checkBkg("Data_pp");
		*/
//		dump_nominal_res("data_pb_nominal", data_pb_f, 0);
//		dump_nominal_res("data_pp_nominal", data_pp_f, 1);
//		checkBkg("data_pb_nominal", 1);
	//	JSRatio("JS_ratio_nominal", "data_pp_nominal", "data_pb_nominal");
//		applyJFF("data_pp_Res_nominal","data_pp_nominal", "JFF_PYTHIA_nominal");
//		applyJFF("data_pp_Res_new","Data_pp", "JFF_PYTHIA");
//		drRatio("dr_JS_newOverNominal","data_pp_Res_new","data_pp_Res_nominal");
//		drawDrOverlay("overlay_data_pp_newAndNominal","data_pp_Res_new","data_pp_Res_nominal",  "new", "nominal",0., 0.99);
		
//		drRatio("dr_signal_newOverNominal","Data_pp","data_pp_nominal");
//		drawDrOverlay("overlay_data_pp_newAndNominal","Data_pp","data_pp_nominal",  "new", "nominal",0., 0.99);
//		drawDrDist("JS_dr_Data_pp", "Data_pp");
	//	show_1D_sig_deta("jff_detaProj", "JFF_PYTHIA");
	//	drRatio("dr_jff_over_gengen","JFF_PYTHIA", "gen_gen_PYTHIA");
		//dumpDr_pp("data_pp_Res_new", "data_pp_Res_new");

//		drawDrDist("jff_nominal_drProj", "JFF_PYTHIA_nominal");
//		drRatio("dr_jffRatio","JFF_PYTHIA", "JFF_PYTHIA_nominal");
		//applyJFF("data_new_jffCorr","Data_pp", "JFF_PYTHIA");
		//applyJFF("data_norm_jffCorr","Data_pp", "JFF_PYTHIA_nominal");
//		drRatio("dr_newOverNorm","data_new_jffCorr", "data_norm_jffCorr");
		//drawDrOverlay("overlay_data_pp", "data_new_jffCorr", "data_norm_jffCorr", "new", "nominal");
//		drawDrOverlay("overlay_data_pp", "data_new_jffCorr", "data_norm_jffCorr", "new", "nominal",0.5, 0.99);
//		drawDrOverlay("jff_overlay_PYTHIA", "JFF_PYTHIA", "JFF_PYTHIA_nominal", "jff_new", "jff_nominal", 0.5, 0.99);
//		checkBkg("rec_gen_sube0_nominal", 1);
//		dump_nominal_res("rec_gen_PYTHIA_nominal", recgen_pythia_norm_f);
//
/* reproduce the nominal resutl */
		/*
		dump_nominal_res("rec_gen_PYTHIA_nominal", recgen_pythia_norm_f);
		dump_nominal_res("gen_gen_PYTHIA_nominal", gengen_pythia_norm_f);
		batch_sig_subtraction("JFF_PYTHIA_nominal", "rec_gen_PYTHIA_nominal", "gen_gen_PYTHIA_nominal");
		dump_nominal_res("data_pp_nominal", data_pp_f);
		dump_nominal_res("data_pb_nominal", data_pb_f, 0);
		//no correction result
		// getting the jff from hydjet
		dump_nominal_res("rec_gen_sube0_nominal", recgen_hydjet_sube0_nom_f, 0);
		dump_nominal_res("gen_gen_sube0_nominal", gengen_hydjet_sube0_nom_f, 0);
		batch_sig_subtraction("jff_hydjet", "rec_gen_sube0_nominal", "gen_gen_sube0_nominal", 0, 0);
		JSRatioTotal("JS_ratio_noCorr", "data_pp_nominal", "data_pb_nominal");
		//apply jff
		dump_nominal_res("data_pb_nominal", data_pb_f, 0);
		cout<<"adding jff for pp.."<<endl;
		applyCorr("data_pp_Res_nominal","data_pp_nominal", "JFF_PYTHIA_nominal");
		applyCorr("data_pb_jffOn","data_pb_nominal", "jff_hydjet", 0, 0);

		//apply spill+jff
		applySpillOver("data_pb_fullCorr","data_pb_jffOn", "spillOver", 0, 0);
		JSRatioTotal("JS_ratio_jffOn", "data_pp_Res_nominal", "data_pb_jffOn");
		JSRatioTotal("JS_ratio_fullCorr", "data_pp_Res_nominal", "data_pb_fullCorr");
		tripleTogether("nominalResultCutAt0p5", "data_pp_nominal", "data_pb_nominal", "data_pp_Res_nominal", "data_pb_jffOn", "data_pp_Res_nominal", "data_pb_fullCorr");
//		batch_sig_subtraction("comp_rec_gen", "rec_gen_sube0_nominal", "rec_gen_PYTHIA", 0, 1);
//		show_2D_sig("colz_recgen_hydjet", "comp_rec_gen",0, 0,  1);
//		batch_sig_subtraction("comp_gen_gen", "gen_gen_sube0_nominal", "gen_gen_PYTHIA", 0, 1);
//		show_2D_sig("colz_gengen_hydjet", "comp_gen_gen",0, 0,  1);

		dump_nominal_res("gen_gen_hydjet_nominal", gengen_hydjet_norm_f, 0);
		dump_nominal_res("rec_gen_hydjet_nominal", recgen_hydjet_norm_f, 0);
		batch_sig_subtraction("diff_hydjet", "rec_gen_hydjet_nominal", "gen_gen_hydjet_nominal", 0, 0);
		drawDrOverlay("overlay_diff","jff_hydjet","diff_hydjet",  "jff", "recgen-gengen",0., 0.99, 0, 0);
		drawdEtaOverlay("overlay_diff","jff_hydjet","diff_hydjet",  "jff", "recgen-gengen", -1.5, 1.499, 0, 0);
		*/
}
