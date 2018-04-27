
#include "utility.h" 
#include "submacros.h" 
void run(){

		TH1D* gengenpbsube0[9][4], *recgenpb[9][4], *jff[9][4];
		getDr("gen_gen_sube0_nominal", gengenpbsube0, 0, 0);
		getDr("rec_gen_hydjet_nominal", recgenpb, 0, 0);
		batch_operation("jff",jff, recgenpb, gengenpbsube0, "diff", 0);

		TH1D* pb[9][4], *pp[9][4], *ppRes[9][4], *pbRes[9][4];
		getDr("data_pp_nominal", pp, 0);
		getDr("data_pb_ref", pb, 0, 0);
		getDr("data_pp_Res_nominal", ppRes, 0);
		getDr("data_pb_final_nominal", pbRes, 0, 0);

		TH1D* ratio1[9][4], *ratio2[9][4], *ratio3[9][4], *ratio4[9][4];

		batch_operation("newPb", pbRes, pb, jff, "diff", 0);
		batch_operation2("newPbJff",ratio1, pbRes, ppRes, "ratio");
		show_1D_pp("ratio_newPbJff", ratio1, 0, 0.49, 0);
		show_1D_pp("newPbjff", jff, 0, 0.49, 0);
		/*
		batch_operation2("noJFF",ratio1, pb, pp, "ratio");
		batch_operation2("ppJFF",ratio2, pb, ppRes, "ratio");
		batch_operation2("pbJFF",ratio3, pbRes, pp, "ratio");
		batch_operation2("allJff",ratio4, pbRes, ppRes, "ratio");

		show_1D_pp("ratio_noJff", ratio1, 0, 0.49, 0);
		show_1D_pp("ratio_ppJff", ratio2, 0, 0.49, 0);
		show_1D_pp("ratio_pbJff", ratio3, 0, 0.49, 0);
		show_1D_pp("ratio_allJff", ratio4, 0, 0.49, 0);

		dump_nominal_res("gen_gen_hydjet_raw_nominal", gengen_hydjet_norm_f,0);
		TH1D* gengenPb[9][4], *gengenPP[9][4];
		getDr("gen_gen_raw_PYTHIA_nominal", gengenPP, 0);
		//getDr("gen_gen_hydjet_raw_nominal", gengenPb, 0, 0);
		getDr("gen_gen_sube0_raw_nominal", gengenPb, 0, 0);

		TH1D* recrecPb[9][4], *recrecPP[9][4] , *recgenPP[9][4];
		getDr("rec_rec_PYTHIA_nominal", recrecPP, 0);

		TH1D* jffpp[9][4], *ratioPP[9][4], *jff[9][4];
		getDr("JFF_PYTHIA_nominal", jffpp, 0);
		TH1D* recrecPP_corrected[9][4];
		batch_operation("recPP",recrecPP_corrected, recrecPP, jffpp, "diff", 1);
		batch_operation("closurePP", ratioPP, recrecPP_corrected, gengenPP, "binomialRatio", 1);
		show_1D_pp("PYTHIA_Closure", ratioPP, 0, 0.49, 1);

		TH1D* ratio[9][4];
		batch_operation2("ratio",ratio, gengenPb, gengenPP, "ratio");
		TH1D** pp1 = flattenTH1D_94(recrecPP);
		TH1D** pp2 = flattenTH1D_94(recrecPP_corrected);
		TH1D** pp3 = flattenTH1D_94(gengenPP);
		showPlot("comparison", 1, 0 ,0.49,3, pp3, pp1, pp2);

		//show_1D_pp("PYTHIA_Closure", ratioPP, 0, 0.49, 1);
		show_1D_pp("ratio_raw_sube0_PbOverPP", ratio, 0, 0.49, 0);
		*/
}
