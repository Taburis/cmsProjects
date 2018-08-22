
#include <iostream>
//#include "corrTable/bJTCTrkCorr_pp_drBin.h"
//#include "jffcorr/nCScorr.h"
//#include "corrTable/bJTCTrkCorr_pp_type2_mulBin.h"
//#include "corrTable/bJTCTrkCorr_pp_mulBin.h"
//auto ppTrkCorrection = new bJTCTrkCorr_pp_drBin("corrTable/PYTHIA_tkCorr_drBin_15May18.root");
//nCScorr *ncscorr;
//auto ppNominalTrkCorrection = new TrkCorr("TrkCorr_Jun7_Iterative_PbPb_etaLT2p4/");
#include "TH1.h"
#include "lib/JTCSkimer.C"
void ReadFileList(std::vector<std::string> &my_file_names, TString file_of_names, bool debug);

void skim(bool doCrab = 0, int jobID=0, int globalCode=0, int startfile=0, int nFiles = 1){
// config here 
	bool doKurtSkim = 0;
	bool usebMC = 1; // HF dedicated sample will change the weight function
	bool doGenJets = 1, doGenTracks = 1;
	float csvCut = 0.9;
	bool isHI =0 , isMC = 1;
	bool do_mixing = 0;
	//-------------------------------------
	bool doCaloJet = 0;
	bool doDCATkCut= 0;
	bool doWTAaxis = 0;
	bool doGenJetAxis = 0;
	float jetetacut = 1.6;
	float jetptcut = 115;

	bool doSmear = 0;
	float smear_mean = 125, smear_sigma = 12;
//-------------------------------------------------------

	TString redirector_xrd = "root://cms-xrd-global.cern.ch//";
	TString filename;

	if(doKurtSkim) usebMC = 0; 
	std::vector<std::string> file_names; 
	TString mixingFileName;
	if(doCrab){
		ReadFileList( file_names, Form("job_input_file_list_%d.txt",jobID), true);
		filename  = file_names.at(0);
		if(doKurtSkim && !isHI && isMC) mixingFileName = "root://cmseos.fnal.gov//eos/uscms/store/user/xiaowang/PYTHIA6_5TeV_forestSkim_finalState_12June18/PYTHIA_5TeV_hiforestSkim_looseMerge/crab_PYTHIA6_5TeV_forestSkim_finalState_12June18/180612_162512/0000/PYTHIA6_1.root";
		else if(doKurtSkim && !isMC && !isHI) mixingFileName = "/mnt/hadoop/store/user/kjung/ppData_5TeV_inclJet_reskimnCS/crab_pp_5TeV_Data_skimJetTrack_reskimnCS_Aug2017/170807_200435/0000/Data_pp_68.root";
		else if(!doKurtSkim && !isMC && !isHI) mixingFileName = "/mnt/hadoop/store/user/wangx/bJTC_5TeV_pp/ppData_5TeV_withCSVv2Skim/CRAB_UserFiles/crab_ppData_5TeV_Skim_withCSVv2_6Aug18/180806_194117/0000/pp_5TeV_Data_100.root";
		else if(!doKurtSkim && !isHI && isMC && !usebMC) mixingFileName = "/mnt/hadoop/store/user/wangx/bJTC_5TeV_pp/PYTHIA6_5TeV_withCSVv2Skim/CRAB_UserFiles/crab_PYTHIA6_5TeV_forestSkim_finalState_1Aug18/180802_033700/0000/PYTHIA6_12.root";
		else if(doKurtSkim && isMC && !isHI)  mixingFileName = "/mnt/hadoop/store/user/kjung/pp_Pythia6MC_5TeV_Skim_Aug2017/pp_Pythia6_5TeV_ppReco_skim/crab_pp_Pythia6MC_5TeV_Skim_Aug2017/170807_150502/0000/Pythia15_1.root";
		else if(doKurtSkim && !isMC && !isHI) mixingFileName = "/mnt/hadoop/store/user/kjung/ppData_5TeV_inclJet_reskimnCS/crab_pp_5TeV_Data_skimJetTrack_reskimnCS_Aug2017/170807_200435/0000/Data_pp_68.root";
		else if(doKurtSkim && isHI && !isMC) mixingFileName = "/mnt/hadoop/store/user/kjung/PbPb_5TeV_MinBiasSkim/Data2015_finalTrkCut_1Mevts.root";
		else if(doKurtSkim && isHI && isMC)mixingFileName = "/mnt/hadoop/store/user/kjung/PbPbMC_Py6H_skim_looseTrkCuts_finalJFFs_lowPtGenPartCut_CymbalTune/crab_PbPb_Pythia6Hydjet_MC_JetTrackSkim_finalizedJFFs_CymbalTune/mergedMixFile/Pythia6Hydjet_PbPbMC_cymbalTune_mixMerged.root";
		if(usebMC){
//			mixingFileName = redirector_xrd+"/store/user/wangx/Pythia6_ppMCskim_bJetforBJTC/list_PYTHIA6_bjet80/CRAB_UserFiles/crab_list_PYTHIA6_bjet80_addGenParticles/180713_213600/0000/PYTHIA6_10.root";
			mixingFileName = "/mnt/hadoop/store/user/wangx/Pythia6_ppMCskim_bJetforBJTC/list_PYTHIA6_bjet120/CRAB_UserFiles/crab_list_PYTHIA6_bjet120_addGenParticles/180814_051856/0000/PYTHIA6_19.root";
			filename  = file_names.at(0);
		}
	}
	else {
		//if(doKurtSkim && !isHI) filename ="/uscms_data/d3/xiaowang/sampleSet/Data_pp_68.root" ;
		if(doKurtSkim&& !isHI) filename ="/uscms_data/d3/xiaowang/sampleSet/Pythia15_22.root" ;
//		else if(!isHI) filename ="root://131.225.204.161:1094//store/user/xiaowang/PYTHIA6_5TeV_forestSkim/PYTHIA_5TeV_hiforestSkim_looseMerge/crab_PYTHIA6_5TeV_forestSkim/180606_044203/0000/PYTHIA6_1.root";
		else if(doKurtSkim && isHI && !isMC) filename = "/uscms_data/d3/xiaowang/sampleSet/Data2015_12.root";
		else if(!doKurtSkim && !isHI && !isMC) filename = "/uscms_data/d3/xiaowang/sampleSet/pp_5TeV_Data_96.root";
		else if(!doKurtSkim && !isHI && isMC) filename = "/uscms_data/d3/xiaowang/sampleSet/PYTHIA6_5TeV_bJetSampe_pthat50.root";
//		filename = "/uscms_data/d3/xiaowang/sampleSet/pp_5TeV_Data_96.root";
		mixingFileName = filename;
	}
cout<<filename<<endl;

	config();
	auto sk = new JTCSkimer();
	TFile *f = TFile::Open(filename);
	TTree* t = (TTree*) f->Get("mixing_tree");

//config setup here---------------------
	sk->doGenJet=doGenJets; sk->doGenTrack=doGenTracks; sk->doCaloJet = doCaloJet;
	sk->csvCut = csvCut;
	sk->connectTree(t, isMC, isHI, doKurtSkim);
	sk->do_mixing = do_mixing;
	sk->doDCATkCut = doDCATkCut;
	sk->trketamaxcut = 2.4;
	sk->doSmear = doSmear;
	sk->jer_mean = smear_mean;
	sk->jer_sigma = smear_sigma;
	sk->jetetacut = jetetacut;
	sk->jetptcut = jetptcut;
	sk->doJetAxisResolution = doGenJetAxis;
	sk->doWTAaxis = doWTAaxis;
	sk->isbMC = usebMC; //using the HF dedicated sample, in this case, the weight function will changed.
//-------------------------------------

//	ncscorr = new nCScorr(!isHI, 1); // implemented when skimming from hiForest
	sk->addMixingTree(mixingFileName);
	//sk->configMixing();
	sk->loop();
	sk->write("correlation.root");
}

void ReadFileList(std::vector<std::string> &my_file_names, TString file_of_names, bool debug)
{
	ifstream file_stream(file_of_names);
	std::string line;
	my_file_names.clear();
	if( debug ) std::cout << "Open file " << file_of_names << " to extract files to run over" << std::endl;
	if( file_stream.is_open() ) {
		if( debug ) std::cout << "Opened " << file_of_names << " for reading" << std::endl;
		int line_num = 0;
		while( file_stream >> line ) {
			if( debug ) std::cout << line_num << ": " << line << std::endl;
			std::string tstring_line(line);
			tstring_line.erase(std::remove(tstring_line.begin(), tstring_line.end(), '"'), tstring_line.end());
			tstring_line.erase(std::remove(tstring_line.begin(), tstring_line.end(), ','), tstring_line.end());
			tstring_line.erase(std::remove(tstring_line.begin(), tstring_line.end(), '['), tstring_line.end());
			tstring_line.erase(std::remove(tstring_line.begin(), tstring_line.end(), ']'), tstring_line.end());
			if( tstring_line != "" ) my_file_names.push_back(tstring_line);
			line_num++;
		}
	} else {
		std::cout << "Error, could not open " << file_of_names << " for reading" << std::endl;
		assert(0);
	}
}
