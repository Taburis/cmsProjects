
#include "utility.h"

struct inputSet {
		bool isRecoJet;
		bool isRecoTk;
		bool isHi = 0;
		bool ismc = 1;
		bool open = 0;
		bool hasCont = 1;
		TString path;
		TH2D** sig_tb, **sig_tt, **sig_in, **sig_tg;
		histCase tb, tt, tg, in, co;
};


TString dataPath = "/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/myCorrelation/";

TString str_pp_data_v1 = dataPath + "bJTC_Data_5TeV_25mix_csvV1_drCorr_15June18.root";
TString str_pp_data_csvV1_85 = dataPath + "bJTC_Data_5TeV_25mix_csvV1_noTkCorr_CSV85.root";
TString str_pp_rr_csvV1_85 = dataPath + "bJTC_RecRec_5TeV_25mix_csvV1_drCorr_CSV85.root";

TString str_gg_v1 = dataPath + "bJTC_GenGen_5TeV_25mix_csvV1_drCorr_CSV90.root";
TString str_gr_v1 = dataPath + "bJTC_GenRec_5TeV_25mix_2p4eta_drBinCorr_25May18.root";
TString str_rg_v1 = dataPath + "bJTC_RecGen_5TeV_25mix_csvV1_drCorr_12June18.root";
TString str_rr_v1 = dataPath + "bJTC_RecRec_5TeV_25mix_csvV1_drCorr_12June18.root";
TString str_gg_v2 = dataPath + "CSVv2/bJTC_GenGen_5TeV_25mix_csvV2_drCorr_12June18.root";
TString str_rg_v2 = dataPath + "CSVv2/bJTC_RecGen_PFJet_noDCA_offTkCorr_csvV2p9_4Aug18.root";
TString str_gr_v2 = dataPath + "CSVv2/bJTC_GenRec_PFJet_noDCA_offTkCorr_csvV2p9_6Aug18.root";
TString str_rr_v2 = dataPath + "CSVv2/bJTC_RecRec_PFJet_noDCA_offTkCorr_csvV2p9_6Aug18.root";
TString str_gg_HF_v2 = dataPath + "HFMC/bJTC_GenGen_HFMC_PFJet_noDCA_offTkCorr_csvV2p9_31July18.root";
TString str_gr_HF_v2 = dataPath + "HFMC/bJTC_GenRec_HFMC_PFJet_noDCA_offTkCorr_csvV2p9_4Aug18.root";
TString str_rg_HF_v2 = dataPath + "HFMC/bJTC_RecGen_HFMC_PFJet_noDCA_offTkCorr_csvV2p9_1Aug18.root";
TString str_rr_HF_v2 = dataPath + "HFMC/bJTC_RecRec_HFMC_PFJet_noDCA_offTkCorr_csvV2p9_1Aug18.root";
TString str_data_v2 = dataPath + "CSVv2/bJTC_Data_PFJet_5TeV_withDCAcut_nominalTkCorr_csvV2p9_10Aug18.root";

TString str_rg_genjt_axis_v1 = dataPath + "bJTC_RecGen_5TeV_25mix_csvV1_noTkCorr_genJetAxis.root";

inputSet p6gg, p6gr, p6rg, p6rr, p6gg2, ppData, p6v2rr, p6v2rg, p6v2gr, p6HFv2gg, p6HFv2rg, p6HFv2rr, p6HFv2gr, p6v2gg;

inputSet ppDatav2;

inputSet p6rg_genjt_axis, ppDataCSV85, p6rrCSV85;


inputSet ppcheck;
inputSet p6rr_calojet;
 
//TF1 **tagCorr;
//TH1D **tkCorr, **resCorr;

void config(){
		ppDatav2.path= str_data_v2;
		ppDatav2.ismc = 0;
		ppDatav2.hasCont = 0;

		p6v2gr.path= str_gr_v2;
		p6v2gr.isRecoJet = 0;
		p6v2gr.isRecoTk = 1;

		p6v2gg.path= str_gg_v2;
		p6v2gg.isRecoJet = 0;
		p6v2gg.isRecoTk = 0;

		p6HFv2gg.path= str_gg_HF_v2;
		p6HFv2gg.isRecoJet = 0;
		p6HFv2gg.isRecoTk = 0;

		p6HFv2rg.path= str_rg_HF_v2;
		p6HFv2rg.isRecoJet = 1;
		p6HFv2rg.isRecoTk = 0;

		p6HFv2gr.path= str_gr_HF_v2;
		p6HFv2gr.isRecoJet = 0;
		p6HFv2gr.isRecoTk = 1;

		p6HFv2rr.path= str_rr_HF_v2;
		p6HFv2rr.isRecoJet = 1;
		p6HFv2rr.isRecoTk = 1;
		
		p6v2rr.path = str_rr_v2;
		p6v2rr.isRecoJet = 1;
		p6v2rr.isRecoTk = 1;

		p6v2rg.path = str_rg_v2;
		p6v2rg.isRecoJet = 1;
		p6v2rg.isRecoTk = 0;

		p6rr_calojet.path = dataPath+"bJTC_RecRec_CaloJet_5TeV_tkDrCorr_csvV1p9_16July18.root";
		p6rr_calojet.isRecoJet = 1;
		p6rr_calojet.isRecoTk = 1;

		ppcheck.path = dataPath + "bJTC_GenGen_5TeV_25mix_csvV1_drCorr_CSV90.root";
		ppcheck.isRecoJet = 0;
		ppcheck.isRecoTk = 0;

		ppData.path = str_pp_data_v1;
		ppData.ismc = 0;

		p6rrCSV85.path = str_pp_rr_csvV1_85;
		p6rrCSV85.isRecoJet = 1;
		p6rrCSV85.isRecoTk = 1;

		ppDataCSV85.path = str_pp_data_csvV1_85;
		ppDataCSV85.ismc = 0;

		p6gg.path = str_gg_v1;
		p6gg.isRecoJet = 0;
		p6gg.isRecoTk = 0;
		p6gg.hasCont = 0;

		p6gr.path = str_gr_v1;
		p6gr.isRecoJet = 0;
		p6gr.isRecoTk = 1;
		p6rg.hasCont = 0;

		p6rg.path = str_rg_v1;
		p6rg.isRecoJet = 1;
		p6rg.isRecoTk = 0;

		p6rg_genjt_axis.path = str_rg_genjt_axis_v1;
		p6rg_genjt_axis.isRecoJet = 1;
		p6rg_genjt_axis.isRecoTk = 0;

		p6rr.path = str_rr_v1;
		p6rr.isRecoJet = 1;
		p6rr.isRecoTk = 1;

}

TString tpname(inputSet &iset){
		TString cap0; 
		if(iset.ismc){
			if(iset.isRecoJet) cap0 = "Rec";
			else cap0 = "Gen";
			if(iset.isRecoTk) cap0+= "Rec";
			else cap0+= "Gen";
		} else cap0 = "Data";
		return cap0;
}

TString histname(inputSet &iset){
		TString cap0; 
		if(iset.ismc){ 
				if(iset.isRecoJet) cap0 = "RecoJet_";
				else cap0 = "GenJet_";
				if(iset.isRecoTk) cap0+= "RecoTrack";
				else cap0+= "GenTrack";
		}else  cap0 = "Data";
		return cap0;
}
void readInput(inputSet &iset){
		if(iset.open) return;
		iset.open = 1;

		TFile *f = TFile::Open(iset.path);
		TString cap0; 
		if(iset.ismc){
		if(iset.isRecoJet) cap0 = "RecoJet_";
		else cap0 = "GenJet_";
		if(iset.isRecoTk) cap0+= "RecoTrack";
		else cap0+= "GenTrack";
		} else cap0 = "Data";

		iset.in = readHistCase("inclJet_"+cap0, "inclJet", f); 
		iset.tg = readHistCase("taggedBJet_"+cap0, "taggedBJet", f );
		if(iset.ismc){	
				iset.tt = readHistCase("taggedTrueBJet_"+cap0, "taggedTrueBJet",f);
				iset.tb = readHistCase("trueBJet_"+cap0, "trueBJet",f);
		}
		if(iset.hasCont){
				iset.co = readHistCase("contJet_"+cap0, "contJet",f);
		}
}



