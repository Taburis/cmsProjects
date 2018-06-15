
#include "utility.h"

struct inputSet {
		bool isRecoJet;
		bool isRecoTk;
		bool isHi = 0;
		bool ismc = 1;
		bool open = 0;
		TString path;
		histCase tb, tt, tg, in;
};

TString dataPath = "/Users/tabris/cmsProjects/iJTC/dataSet/bJTC/myCorrelation/";

TString str_gg_v1 = dataPath + "bJTC_GenGen_5TeV_25mix_2p4eta_drBinCorr_25May18.root";
TString str_gr_v1 = dataPath + "bJTC_GenRec_5TeV_25mix_2p4eta_drBinCorr_25May18.root";
TString str_rg_v1 = dataPath + "bJTC_RecGen_5TeV_25mix_csvV1_drCorr_12June18.root";
TString str_rr_v1 = dataPath + "bJTC_RecRec_5TeV_25mix_csvV1_drCorr_12June18.root";
TString str_gg_v2 = dataPath + "bJTC_GenGen_5TeV_25mix_csvV2_drCorr_12June18.root";
TString str_gr_v2 = dataPath + "bJTC_GenRec_5TeV_25mix_csvV2_drCorr_8June18.root";

inputSet p6gg, p6gr, p6rg, p6rr, p6gg2, p6gr2;

void config(){
		p6gg.path = str_gg_v1;
		p6gg.isRecoJet = 0;
		p6gg.isRecoTk = 0;
		p6gr.path = str_gr_v1;
		p6gr.isRecoJet = 0;
		p6gr.isRecoTk = 1;
		
		p6rg.path = str_rg_v1;
		p6rg.isRecoJet = 1;
		p6rg.isRecoTk = 0;

		p6rr.path = str_rr_v1;
		p6rr.isRecoJet = 1;
		p6rr.isRecoTk = 1;

		p6gg2.path = str_gg_v2;
		p6gg2.isRecoJet = 0;
		p6gg2.isRecoTk = 0;
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
		}else {
				cap0 = "Data";
		}
		return cap0;
}
void readInput(inputSet &iset){
		if(iset.open) return;
		iset.open = 1;
		
		TFile *f = TFile::Open(iset.path);
		TString cap0; 
		if(iset.isRecoJet) cap0 = "RecoJet_";
		else cap0 = "GenJet_";
		if(iset.isRecoTk) cap0+= "RecoTrack";
		else cap0+= "GenTrack";

		iset.in = readHistCase("inclJet_"+cap0, "inclJet", f); 
		iset.tg = readHistCase("taggedBJet_"+cap0, "taggedBJet", f );
		if(iset.ismc){	iset.tt = readHistCase("taggedTrueBJet_"+cap0, "taggedTrueBJet",f);
				iset.tb = readHistCase("trueBJet_"+cap0, "trueBJet",f);
		}
}

//auto gengen_inclf = TFile::Open(dataPath + "PYTHIA_GenGen_inclusiveJTC.root");
//auto genrecf_dcaCut = TFile::Open(dataPath + "bJTC_GenRec_5TeV_25mix_dcaCutIncl_drBinCorr_31May18.root");

TString tkCorrFile = dataDumpPath+"bJet_tkCorr.root";
// config for the input 

bool isData      = 0;

bool doDecontamination = 1;
bool doBiasCorrection  = 1;
bool doTkCorr          = 1;
bool isHi              = 0;

bool isValidation  = 1;

// setting up the variables will be used according to the config

