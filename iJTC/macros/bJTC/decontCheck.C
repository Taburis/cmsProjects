
#ifndef standardSequence_H
#include "standardSequence.C"
#include "syst_utility.h"
#endif

void doCheck(inputSet &iset, float purity, bool isNumber = 1){

		TString prefix = iset.isHi ? "PH_": ""; 
		readInput(iset);
		TH2D** rs_tg, **rs_tt, **rs_in;
		if(isNumber){
				rs_tg = doMixingCorr("step2_sig_taggedB", iset.tg.sig_raw, iset.tg.mixing_raw, 0);
				rs_in = doMixingCorr("step2_sig_incl", iset.in.sig_raw, iset.in.mixing_raw, 1);
				rs_tt = doMixingCorr("step2_sig_ttB", iset.tt.sig_raw, iset.tt.mixing_raw, 0);
		} else {
				rs_tg = doMixingCorr("step2_sig_taggedB", iset.tg.sig_pTweighted_raw, iset.tg.mixing_raw, 0);
				rs_in = doMixingCorr("step2_sig_incl", iset.in.sig_pTweighted_raw, iset.in.mixing_raw, 1);
				rs_tt = doMixingCorr("step2_sig_ttB", iset.tt.sig_pTweighted_raw, iset.tt.mixing_raw, 0);
		}
		TString label = tpname(iset), endname = tpname(iset);
		if(!isNumber) label+="_pTweighted";
		auto deta_rs_tt = projX("deta_rs_tt", rs_tt);
		auto dr_tt = getDr("dr_tt", rs_tt);
		auto dr_tt_err = getDrErr("dr_tt_err", dr_tt, deta_rs_tt, 0);
		auto dr_before= getDr("dr_before",rs_tg);
		showClosure("decontCheck/"+prefix+"dr_beforeDecont_"+label+".pdf", iset.isHi, 0, .99, 0.5, 1.5, "tag.("+endname+")", "t&t ("+endname+")", dr_before, dr_tt);

		scale<TH2D>(rs_in, 1-purity);
		//auto rs_pur = binary_operation<TH2D>("raw_sig_pur", rs_tg, rs_in, "diff");
		TH2D** rs_pur = new TH2D*[nPt*nCent]; 
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						cout<<i<<", "<<j<<endl;
						rs_pur[i+nPt*j]=(TH2D*) rs_tg[i+nPt*j]->Clone(Form("raw_sig_pur_%d_%d",i,j));
						if(i<nPt-2){ rs_pur[i+nPt*j]->Add(rs_in[i+nPt*j], -1);
								rs_pur[i+nPt*j]->Scale(1.0/purity);
						}
				}
		}
		auto deta_rs_pur = projX("deta_rs_pur", rs_pur);

		auto dr_pur= getDr("dr_pur",rs_pur);
		auto dr_pur_err= getDrErr("dr_pur_err",dr_pur, deta_rs_pur , 0);

		//showClosure("decontCheck/"+prefix+"dr_afterDecont_"+label+".pdf", iset.isHi, 0, .99, 0.5, 1.5, "decont.("+endname+")", "t&t ("+endname+")", dr_pur, dr_tt);
		showClosure_syst("decontCheck/"+prefix+"dr_afterDecont_"+label+".pdf", iset.isHi, 0, .99, 0.5, 1.5, "decont.("+endname+")", "t&t ("+endname+")", dr_pur, dr_pur_err,dr_tt, dr_tt_err, "binomialRatio", "binomialRatio");
		free<TH1D>(dr_tt);
		free<TH1D>(deta_rs_tt);
		free<TH1D>(deta_rs_pur);
		free<TH1D>(dr_tt_err);
		free<TH1D>(dr_pur);
		free<TH1D>(dr_pur_err);
		free<TH1D>(dr_before);
		free<TH2D>(rs_tg);
		free<TH2D>(rs_in);
		free<TH2D>(rs_tt);
		free<TH2D>(rs_pur);
}

void decontCheck(){
		config();
	//	readInput(p6gg);
	//	readInput(p6gr);
	//	readInput(p6rr);
	//	readInput(p6rg);
	/*
		doCheck(p6gg, 0.7, 0);
		doCheck(p6gg, 0.7, 1);
		doCheck(phgg, 0.7, 0);
		doCheck(phgg, 0.7, 1);
		*/
		doCheck(p6gr, 0.7, 0);
		doCheck(p6gr, 0.7, 1);
		doCheck(p6rg, 0.7, 0);
		doCheck(p6rg, 0.7, 1);
		doCheck(p6rr, 0.7, 0);
		doCheck(p6rr, 0.7, 1);
}

TH1D** decontamination2(inputSet &iset){

		readInput(iset);
		TString endname = tpname(iset);

		//		auto hc_tb = readHistCase("trueBJet_"+cap0, "trueBJet",f);

		auto rs_tg = doMixingCorr("step2_sig_taggedB", iset.tg.sig_raw, iset.tg.mixing_raw, 0);
		auto rs_in = doMixingCorr("step2_sig_incl", iset.in.sig_raw, iset.in.mixing_raw, 1);
		auto rs_tt = doMixingCorr("step2_sig_ttB", iset.tt.sig_raw, iset.tt.mixing_raw, 0);
		//		auto rs_tb = doMixingCorr("step2_sig_tB", hc_tb.sig_raw, hc_tb.mixing_raw, 0);
		//checkSideBand("taggedBgengen", sig_tg_gg);
		//checkSideBand("taggedTrueBgengen", sig_tt_gg);

		auto dr_tt = getDr("dr_tt", rs_tt);
		auto dr_pur= getDr("dr_pur",rs_tg);
		//		auto dr_tb= getDr("dr_tb",rs_tb);
		auto ratio1 = binary_operation<TH1D>("ratio1", dr_pur, dr_tt, "binomialRatio");
		showClosure("dr_tag_tt_"+endname+".pdf", iset.isHi, 0, .99, 0.5, 1.5, "tag.("+endname+")", "t&t ("+endname+")", dr_pur, dr_tt);
		//		showPlot("drRatio_before_decont", 0, 1, 0, 0.99, 0.5, 1.5, 1, ratio);

		float purity = 0.7;
		scale<TH2D>(rs_in, 1-purity);
		auto rs_pur = binary_operation<TH2D>("raw_sig_pur", rs_tg, rs_in, "diff");
		scale(rs_pur, 1.0/purity);

		free<TH1D>(dr_tt);
		free<TH1D>(dr_pur);
		dr_tt = getDr("dr_tt", rs_tt);
		dr_pur= getDr("dr_pur",rs_pur);
		auto ratio2 = binary_operation<TH1D>("ratio2", dr_pur, dr_tt, "binomialRatio");
		showPlot("drRatio_before_after_decont_"+endname+".gif", 0, 1, 0, 0.99, 0.5, 1.5, 2, ratio1, ratio2);

		showClosure("dr_afterDecont_"+endname+".pdf", iset.isHi, 0, .99, 0.5, 1.5, "decont.("+endname+")", "t&t ("+endname+")", dr_pur, dr_tt);
		//showPlot("drOverlay_decont_ttb", 0, 0, 0, 0.99, 0.5, -1.5, 2, dr_tt, dr_pur);
		free<TH1D>(dr_tt);
		free<TH1D>(dr_pur);
		auto sig_tt = doBkgSub("sig_ttB", rs_tt);
		auto sig_pur= doBkgSub("sig_pur", rs_pur);
		dr_tt = getDr("dr_tt2", sig_tt);
		dr_pur= getDr("dr_pur2",sig_pur);
		showPlot("drOverlay_bkgSubtracted_"+endname+".gif", 0, 0, 0, 0.99, 0.5, -1.5, 2, dr_tt, dr_pur);

		return dr_pur;
}

void checkContaminationBias(inputSet &iset){
		auto f = TFile::Open(iset.path);
		TString cap0; 
		if(iset.isRecoJet) cap0 = "RecoJet_";
		else cap0 = "GenJet_";
		if(iset.isRecoTk) cap0+= "RecoTrack";
		else cap0+= "GenTrack";
		//		cap0+="_noCorr";
		auto rs_tg = read_flatten<TH2D>(iset.path, "taggedBJet_"+cap0+"_noCorr"); 
		auto rs_tt = read_flatten<TH2D>(iset.path, "taggedTrueBJet_"+cap0+"_noCorr"); 
		auto rs_in = read_flatten<TH2D>(iset.path, "inclJet_"+cap0+"_noCorr");	
		auto mix_tg = read_flatten<TH2D>(iset.path, "taggedBJet_"+cap0+"_mixing_noCorr"); 
		auto mix_tt = read_flatten<TH2D>(iset.path, "taggedTrueBJet_"+cap0+"_mixing_noCorr"); 
		auto mix_in = read_flatten<TH2D>(iset.path, "inclJet_"+cap0+"_mixing_noCorr");	

		/*
		   auto rs_tg = read_flatten<TH2D>(iset.path, "taggedBJet_"+cap0); 
		   auto rs_tt = read_flatten<TH2D>(iset.path, "taggedTrueBJet_"+cap0); 
		   auto rs_in = read_flatten<TH2D>(iset.path, "inclJet_"+cap0);	
		   auto mix_tg = read_flatten<TH2D>(iset.path, "taggedBJet_"+cap0); 
		   auto mix_tt = read_flatten<TH2D>(iset.path, "taggedTrueBJet_"+cap0); 
		   auto mix_in = read_flatten<TH2D>(iset.path, "inclJet_"+cap0);	
		   */

		auto sig_tg = doMixingCorr("step2_sig_cont", rs_tg, mix_tg, 0);
		auto sig_tt = doMixingCorr("step2_sig_cont", rs_tt, mix_tt, 0);
		auto sig_cont = binary_operation<TH2D>("rs_cont", sig_tg, sig_tt, "diff");
		auto sig_in = doMixingCorr("step2_sig_incl", rs_in, mix_in, 1);

		auto dr_raw0 = getDr("dr_raw0", sig_cont);
		auto dr_raw1 = getDr("dr_raw1", sig_in);
		auto dr_raw2 = getDr("dr_raw2", sig_tt);

		auto jt_tt = (TH1D*) f->Get("taggedTrueBJet_corrpt_0");
		auto jt_tg = (TH1D*) f->Get("taggedBJet_corrpt_0");
		auto jt_in = (TH1D*) f->Get("inclJet_corrpt_0");
		float num_cont = jt_tg->Integral()-jt_tt->Integral();
		float num_incl = jt_in->Integral();
		scale(sig_cont, 1.0/num_cont);
		scale(sig_tt, 1.0/jt_tt->Integral());
		scale(sig_in, 1.0/num_incl);
		//auto signal_cont = doBkgSub("signal_cont", sig_cont);
		//auto signal_in = doBkgSub("signal_in", sig_in);
		auto dr0 = getDr("dr0", sig_cont);
		auto dr1 = getDr("dr1", sig_in);
		auto dr2 = getDr("dr2", sig_tt);
		showPlot("drOverlay_cont_incl.gif", 0, 0, 0, 0.99, 0.5, -1.5, 2, dr0, dr1);
		auto ratio = binary_operation<TH1D>("contaminationBias", dr0, dr1, "ratio");
		showPlot("drRatio_contOverIncl.gif", 0, 1, 0, 0.99, 0.5, 1.5, 1, ratio);
		showStack("drStack_contOverIncl", 0, 1, 0, 0.99, 0.5, 1.5, 2, dr_raw0, dr_raw2);
}
