
#ifndef ROOT_TH2D
#include "TH2D.h"
#endif

TH2D* pb_sig[9][5];
TH2D* pb_me [9][5];
TH1D* jet_pt[5];
namespace kurt_reg{
		TFile *GenGen_MC_pp_f ;
		TFile *RecRec_MC_pp_f ;
		TFile *GenGen_MC_pp_f_calo5jet ;
		TFile *GenGen_MC_pp_f_calo3jet ;
		TFile *GenGen_MC_pb_f ;
		TFile *RecRec_MC_pb_f ;
		TFile *GenGen_MC_pb_f_sub0 ;
		TFile *RecGen_MC_pb_f_sub0 ;
		TFile *data_pb_f ;
		void loadConfig(){
				std::cout<<"loading kurt's configuration.."<<std::endl;
				//pp input 
				data_pb_f =TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/PbPb_5TeVData_fixMix_fineBinTrkCorrs_withTrkCorrEtaSymmV2_finalJFFs.root");	

				GenGen_MC_pb_f_sub0 =TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/PbPb_5TeVMC_withMix_GenGen_CymbalTune_fullCymbalCorrs_inclJetBinning.root");	

				GenGen_MC_pp_f =TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/pp_Pythia6MC_GenGen_withMix_inclJetBinning_finalJFFs.root");	
				RecRec_MC_pp_f =TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/pp_Pythia6MC_RecoReco_withMix_inclJetBinning_finalJFFs.root");	
				// pb mc input 
				GenGen_MC_pb_f =TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/PbPb_5TeVMC_withMix_GenGen_CymbalTune_fullCymbalCorrs_inclJetBinning.root");	
				RecRec_MC_pb_f =TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/PbPb_5TeVMC_withMix_RecoReco_CymbalTune_fullCymbalCorrs_inclJetBinning.root");	
		}
		TString trk_tag[] = {"TrkPt0p7", "TrkPt1", "TrkPt2","TrkPt3","TrkPt4","TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt999"};
		TString cent_tag[]= {"Cent0","Cent10","Cent30","Cent50","Cent70","Cent100"};
		TH1D* GenJet_pp;	
		TH2D* GenGen_MC_pp_sig[9];	
		TH2D* GenGen_MC_pp_mix[9];	
		TH1D* GenJet_pb[5];	
		TH2D* GenGen_MC_pb_sig[9][5];	
		TH2D* GenGen_MC_pb_mix[9][5];	
		TH2D* RecRec_MC_pb_sig[9][5];
		TH2D* RecRec_MC_pb_mix[9][5];

		TH1D* hjet[5];
		TH2D* mix[9][5];
		TH2D* sig[9][5];

		void getRaw2D(TFile *f, TString cap){
				float trkBinWidth[] = {0.3, 1, 1,1, 4, 4, 4, 4, 50}; 
				TString stmp;
				for(int j=0; j<5; ++j){
						stmp = cap+"_all_jets_corrpT"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000";
						hjet[j]= (TH1D*)f->Get(stmp);
						for(int i=0;i<9; ++i){
								stmp = cap+"_hJetTrackME"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000_"
										+trk_tag[i]+"_"+trk_tag[i+1];
								mix[i][j]=(TH2D*)f->Get(stmp);
								stmp = cap+"_hJetTrackSignalBackground"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000_"
										+trk_tag[i]+"_"+trk_tag[i+1];
								sig[i][j]=(TH2D*)f->Get(stmp);
								sig[i][j]->Scale(1./trkBinWidth[i]);
						}
				}
				hjet[3]->Add(GenJet_pb[4]);
				for(int i=0;i<9; ++i){
						//combined the last two bins together
						sig[i][3]->Add(sig[i][4]);
						for(int j=0; j<4; ++j){
								//adding the cent50-70 and cent70-100 together:
								sig[i][j]->Scale(1.0/hjet[j]->Integral());
						}
				}
		}

		void getGenGen_sube0();

		void getGenGen_pb(){
				float trkBinWidth[] = {0.3, 1, 1,1, 4, 4, 4, 4, 50}; 
				TString stmp;
				for(int j=0; j<5; ++j){
						stmp = "GenJet_GenTrack_all_jets_corrpT"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000";
						GenJet_pb[j]= (TH1D*)GenGen_MC_pb_f->Get(stmp);
						for(int i=0;i<9; ++i){
								stmp = "GenJet_GenTrack_hJetTrackME"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000_"
										+trk_tag[i]+"_"+trk_tag[i+1];
								GenGen_MC_pb_mix[i][j]=(TH2D*)GenGen_MC_pb_f->Get(stmp);
								stmp = "GenJet_GenTrack_hJetTrackSignalBackground"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000_"
										+trk_tag[i]+"_"+trk_tag[i+1];
								GenGen_MC_pb_sig[i][j]=(TH2D*)GenGen_MC_pb_f->Get(stmp);
								GenGen_MC_pb_sig[i][j]->Scale(1./trkBinWidth[i]);
						}
				}
				//combined the last two bins together
				GenJet_pb[3]->Add(GenJet_pb[4]);
				for(int i=0;i<9; ++i){
						//combined the last two bins together
						GenGen_MC_pb_sig[i][3]->Add(GenGen_MC_pb_sig[i][4]);
						for(int j=0; j<4; ++j){
								//adding the cent50-70 and cent70-100 together:
								GenGen_MC_pb_sig[i][j]->Scale(1.0/GenJet_pb[j]->Integral());
						}
				}
		}

		void getGenGen_pb_sub0(){
				float trkBinWidth[] = {0.3, 1, 1,1, 4, 4, 4, 4, 50}; 
				TString stmp;
				for(int j=0; j<5; ++j){
						stmp = "GenJet_GenTrack_all_jets_corrpT"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000";
						GenJet_pb[j]= (TH1D*)GenGen_MC_pb_f_sub0->Get(stmp);
						for(int i=0;i<9; ++i){
								stmp = "GenJet_GenTrack_hJetTrackME"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000_"
										+trk_tag[i]+"_"+trk_tag[i+1];
								GenGen_MC_pb_mix[i][j]=(TH2D*)GenGen_MC_pp_f->Get(stmp);
								stmp = "GenJet_GenTrack_hJetTrackSignalBackground"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000_"
										+trk_tag[i]+"_"+trk_tag[i+1];
								GenGen_MC_pb_sig[i][j]=(TH2D*)GenGen_MC_pb_f_sub0->Get(stmp);
								GenGen_MC_pb_sig[i][j]->Scale(1./trkBinWidth[i]);
						}
				}
				GenJet_pb[3]->Add(GenJet_pb[4]);
				for(int i=0;i<9; ++i){
						//combined the last two bins together
						GenGen_MC_pb_sig[i][3]->Add(GenGen_MC_pb_sig[i][4]);
						for(int j=0; j<4; ++j){
								//adding the cent50-70 and cent70-100 together:
								GenGen_MC_pb_sig[i][j]->Scale(1.0/GenJet_pb[j]->Integral());
						}
				}
		}

		void getGenGen_pp(){
				TString stmp;
				float trkBinWidth[] = {0.3, 1, 1,1, 4, 4, 4, 4, 50}; 
				stmp = "GenJet_GenTrack_all_jets_corrpT"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000";
				GenJet_pp= (TH1D*)GenGen_MC_pp_f->Get(stmp);
				for(int i=0;i<9; ++i){
						stmp = "GenJet_GenTrack_hJetTrackME"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000_"
								+trk_tag[i]+"_"+trk_tag[i+1];
						GenGen_MC_pp_mix[i]=(TH2D*)GenGen_MC_pp_f->Get(stmp);
						stmp = "GenJet_GenTrack_hJetTrackSignalBackground"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000_"
								+trk_tag[i]+"_"+trk_tag[i+1];
						GenGen_MC_pp_sig[i]=(TH2D*)GenGen_MC_pp_f->Get(stmp);
						GenGen_MC_pp_sig[i]->Scale(1.0/GenJet_pp->Integral()/trkBinWidth[i]);
				}
		}

		void getGenGen_pp_calo5jet(){
				TString stmp;
				float trkBinWidth[] = {0.3, 1, 1,1, 4, 4, 4, 4, 50}; 
				stmp = "GenJet_GenTrack_all_jets_corrpT"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000";
				GenJet_pp= (TH1D*)GenGen_MC_pp_f_calo5jet->Get(stmp);
				for(int i=0;i<9; ++i){
						stmp = "GenJet_GenTrack_hJetTrackME"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000_"
								+trk_tag[i]+"_"+trk_tag[i+1];
						GenGen_MC_pp_mix[i]=(TH2D*)GenGen_MC_pp_f->Get(stmp);
						stmp = "GenJet_GenTrack_hJetTrackSignalBackground"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000_"
								+trk_tag[i]+"_"+trk_tag[i+1];
						//cout<<stmp<<endl;
						GenGen_MC_pp_sig[i]=(TH2D*)GenGen_MC_pp_f_calo5jet->Get(stmp);
						GenGen_MC_pp_sig[i]->Scale(1.0/GenJet_pp->Integral()/trkBinWidth[i]);
				}
		}

		void getGenGen_pp_calo3jet(){
				TString stmp;
				float trkBinWidth[] = {0.3, 1, 1,1, 4, 4, 4, 4, 50}; 
				stmp = "GenJet_GenTrack_all_jets_corrpT"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000";
				GenJet_pp= (TH1D*)GenGen_MC_pp_f_calo3jet->Get(stmp);
				for(int i=0;i<9; ++i){
						stmp = "GenJet_GenTrack_hJetTrackME"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000_"
								+trk_tag[i]+"_"+trk_tag[i+1];
						GenGen_MC_pp_mix[i]=(TH2D*)GenGen_MC_pp_f->Get(stmp);
						stmp = "GenJet_GenTrack_hJetTrackSignalBackground"+cent_tag[0]+"_"+cent_tag[1]+"_Pt100_Pt1000_"
								+trk_tag[i]+"_"+trk_tag[i+1];
						//			std::cout<<stmp<<std::endl;
						GenGen_MC_pp_sig[i]=(TH2D*)GenGen_MC_pp_f_calo3jet->Get(stmp);
						GenGen_MC_pp_sig[i]->Scale(1.0/GenJet_pp->Integral()/trkBinWidth[i]);
				}
		}

		void getPbData(){
				TString stmp;
				float trkBinWidth[] = {0.3, 1, 1,1, 4, 4, 4, 4, 50}; 
				for(int j=0; j<5; ++j){
						stmp = "Data_all_jets_corrpT"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000";
						jet_pt[j]= (TH1D*)data_pb_f->Get(stmp);
						for(int i=0;i<9; ++i){
								stmp = "Data_hJetTrackME"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000_"
										+trk_tag[i]+"_"+trk_tag[i+1];
								pb_me[i][j]=(TH2D*)data_pb_f->Get(stmp);
								stmp = "Data_hJetTrackSignalBackground"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt1000_"
										+trk_tag[i]+"_"+trk_tag[i+1];
								pb_sig[i][j]=(TH2D*)data_pb_f->Get(stmp);
								pb_sig[i][j]->Scale(1./trkBinWidth[i]);
						}
				}

				jet_pt[3]->Add(jet_pt[4]);
				for(int i=0;i<9; ++i){
						//combined the last two bins together
						pb_sig[i][3]->Add(pb_sig[i][4]);
						for(int j=0; j<4; ++j){
								//adding the cent50-70 and cent70-100 together:
								pb_sig[i][j]->Scale(1.0/jet_pt[j]->Integral());
						}
				}
		}

}

namespace hallie_reg{
		//TFile * py_f = new TFile("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/Jet_Shapes_hallie.root");
		TFile * py_f = new TFile("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/Jet_Shapes_kurt_pp.root");
		TFile * py_proj_f = new TFile("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/Particle_Yields_combined.root");
		//TFile * py_proj_f = new TFile("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/Particle_Yields_fixErr.root");
		TFile * js_f = new TFile("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/Jet_Shapes_pTweighted.root");
		TString trk_tag[] = {"TrkPt07", "TrkPt1", "TrkPt2","TrkPt3","TrkPt4","TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt300"};
		TString cent_tag[]= {"Cent0","Cent10","Cent30","Cent50","Cent100"};
		TString errType[] = {"_rel", "_bkg","_me", "_spill", "_jff"};

		TH1D* spill_deta[9][4];
		TH1D* jff_deta[9][4];
		TH1D* spill_dphi[9][4];
		TH1D* jff_dphi[9][4];

		TFile * spill_py_f= TFile::Open("/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/spill_over/Inclusive_Hydjet_SpillOvers.root");
		TFile * spill_js_f= TFile::Open("/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/spill_over/Inclusive_Hydjet_SpillOvers_pTweighted.root");
		TFile * jff_pb_py_f= TFile::Open("/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/jff_residual/Inclusive_Hydjet_JFFResiduals.root");
		TFile * jff_pb_js_f= TFile::Open("/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/jff_residual/Inclusive_Hydjet_JFFResiduals_pTweighted.root");
		void getCorr(bool isjs = 0){
				TFile * sf, *jf;
				int npt;
				if(isjs) {
						sf = spill_js_f; jf = jff_pb_js_f; npt =9;
				}
				else {
						sf = spill_py_f; jf = jff_pb_py_f; npt =8;
				}
				for(int i=0; i<npt; ++i){
						for(int j=0; j<4; ++j){
								TString tmp ="Eta_SpillOver_Points_"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt300_" + trk_tag[i] + "_" + trk_tag[i+1];
								spill_deta[i][j]=(TH1D*) sf->Get(tmp);
								tmp ="Phi_SpillOver_Points_"+cent_tag[j]+"_"+cent_tag[j+1]+"_Pt100_Pt300_" + trk_tag[i] + "_" + trk_tag[i+1];
								spill_dphi[i][j]=(TH1D*) sf->Get(tmp);
								if( isjs) tmp ="Raw_JFF_Residual_Eta_pTweighted"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i] + "_" + trk_tag[i+1];
								else tmp ="Raw_JFF_Residual_Eta_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i] + "_" + trk_tag[i+1];
								jff_deta[i][j]=(TH1D*) jf->Get(tmp);
								if( isjs) tmp ="Raw_JFF_Residual_Phi_pTweighted"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i] + "_" + trk_tag[i+1];
								else tmp ="Raw_JFF_Residual_Phi_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i] + "_" + trk_tag[i+1];
								jff_dphi[i][j]=(TH1D*) jf->Get(tmp);
						}
				}
		}
		TH1D* dphi_syst_err[5][9][4];
		TH1D* deta_syst_err[5][9][4];
		void getAllErr_proj(int npt){
				TString tmp;
				for(int i=0; i<npt; ++i){
						for(int j=0; j<4; ++j){
								for(int k=0; k<5; ++k){
										tmp = "dPhi_Syst_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1]+
												errType[k];
										dphi_syst_err[k][i][j]=(TH1D*) py_proj_f->Get(tmp);
			//							cout<<tmp<<endl;
										tmp = "dEta_Syst_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1]+
												errType[k];
										deta_syst_err[k][i][j]=(TH1D*) py_proj_f->Get(tmp);
								}
						}
				}
		}

		float trkPtBin[] = {0.7, 1, 2, 3, 4, 8, 12, 16,20,300};
		float centBin [] = {100, 50, 30, 10, 0};
		float errBin  []  = {0, 1, 2, 3, 4, 5};
//		auto py_err_table = new xTHD13("err","", 9, trkPtBin, 4, centBin, 5, errBin);
		TH1D* py_dr_syst_err[5][9][5];
		void getPYErr_dr(){
				TString tmp;
				for(int i=0; i<9; ++i){
						for(int j=0; j<4; ++j){
								for(int k=0; k<5; ++k){
										tmp = "dR_Syst_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1]+errType[k];
										py_dr_syst_err[k][i][j]=(TH1D*)py_f->Get(tmp);
										if( j==0){
												tmp = "dR_Syst_pp_"+trk_tag[i]+"_"+trk_tag[i+1]+errType[k];
												py_dr_syst_err[k][i][4]=(TH1D*)py_f->Get(tmp);
												cout<<i<<", "<<j<<endl;
												cout<<py_dr_syst_err[k][i][4]->GetName()<<endl;
										}
								}
						}
				}
		}
		TH1D* integral[5];
		TH1D* integral_syst[5];
		TH1D* diff[4];
		TH1D* diff_syst[4];
		TString tmp;

		TH1D* js_dr[9][5];
		TH1D* js_dr_err[9][5];
		TH1D* js_dr_all[5];
		TH1D* js_dr_err_all[5];

		TH1D* py_deta    [8][5];
		TH1D* py_dphi    [8][5];
		TH1D* py_dr      [8][5];
		TH1D* py_dr_err  [8][5];
		TH1D* py_deta_err[8][5];
		TH1D* py_dphi_err[8][5];

		TH1D* py_deta_all[5];
		TH1D* py_dphi_all[5];
		TH1D* py_dr_all   [5];
		TH1D* py_deta_err_all[5];
		TH1D* py_dphi_err_all[5];
		TH1D* py_dr_err_all[5];
		void getPY_dr(){
				for(int i=0;i<8; ++i){
						tmp = "JetShape2_Yield_BkgSub_Inclusive_pp_"+trk_tag[i]+"_"+trk_tag[i+1];
						py_dr[i][4]=(TH1D*)py_f->Get(tmp);
						tmp = "Jet_Shape_SystErr_pp_"+trk_tag[i]+"_"+trk_tag[i+1];
						py_dr_err[i][4]=(TH1D*)py_f->Get(tmp);
						for(int j=0; j<4; ++j){
								tmp = "JetShape2_Yield_BkgSub_Inclusive_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								cout<<tmp<<endl;
								py_dr[i][j]=(TH1D*)py_f->Get(tmp);
								tmp = "Jet_Shape_SystErr_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								py_dr_err[i][j]=(TH1D*)py_f->Get(tmp);
						}
				}
				py_dr_err_all[4] = (TH1D*) py_dr_err[0][4]->Clone("PY_pp_syst_err_all");
				py_dr_all[4] = (TH1D*) py_dr[0][4]->Clone("PY_pp_all");
				for(int i=1; i<8; ++i){
						py_dr_err_all[4]->Add(py_dr_err[i][4]);
						py_dr_all[4]->Add(py_dr[i][4]);
				}

/*
				*/
				for(int j=0; j<4; ++j){
						py_dr_err_all[j]=(TH1D*) py_dr_err[0][j]->Clone(Form("PY_err_all_%d",j));
						py_dr_all[j] = (TH1D*) py_dr[0][j]->Clone(Form("PY_all_%d",j));
						for(int i=1;i<8; ++i){
								py_dr_all[j]->Add(py_dr[i][j]);
								py_dr_err_all[j]->Add(py_dr_err[i][j]);
						}
				}

				//get hallie's total
				tmp = "Jet_Shape_SystErr_pp_TrkPt300_";
				py_dr_err_all[4] = (TH1D*) py_f->Get(tmp);
				for(int i=0; i<4; ++i){
						tmp = "Jet_Shape_SystErr_"+cent_tag[i]+"_"+cent_tag[i+1]+"_TrkPt300_";
						py_dr_err_all[i] = (TH1D*) py_f->Get(tmp);
				}
		}

		void getPY_proj(){
				for(int i=0;i<8; ++i){
						tmp = "Proj_dEta_pp_"+trk_tag[i]+"_"+trk_tag[i+1]+"_Rebin";
						py_deta[i][4]=(TH1D*)py_proj_f->Get(tmp);
						tmp = "Proj_dPhi_pp_"+trk_tag[i]+"_"+trk_tag[i+1]+"_Rebin";
						py_dphi[i][4]=(TH1D*)py_proj_f->Get(tmp);
						tmp = "dEta_Syst_pp_Cent0_Cent10_"+trk_tag[i]+"_"+trk_tag[i+1];
						py_deta_err[i][4]=(TH1D*)py_proj_f->Get(tmp);
						tmp = "dPhi_Syst_pp_Cent0_Cent10_"+trk_tag[i]+"_"+trk_tag[i+1];
						//cout<<tmp<<endl;
						py_dphi_err[i][4]=(TH1D*)py_proj_f->Get(tmp);
						for(int j=0; j<4; ++j){
								tmp = "Proj_dEta_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1]+"_Rebin";
								py_deta[i][j]=(TH1D*)py_proj_f->Get(tmp);
								cout<<py_deta[i][j]->GetName()<<endl;
								tmp = "Proj_dPhi_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1]+"_Rebin";
								py_dphi[i][j]=(TH1D*)py_proj_f->Get(tmp);
								tmp = "dEta_Syst_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								py_deta_err[i][j]=(TH1D*)py_proj_f->Get(tmp);
								tmp = "dPhi_Syst_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								py_dphi_err[i][j]=(TH1D*)py_proj_f->Get(tmp);
						}
				}
				for(int j=0; j<5; ++j){
						py_deta_all[j]=(TH1D*) py_deta[0][j]->Clone(Form("py_deta_all_%d",j));
						py_dphi_all[j]=(TH1D*) py_dphi[0][j]->Clone(Form("py_dphi_all_%d",j));
						py_deta_err_all[j]=(TH1D*) py_deta_err[0][j]->Clone(Form("py_deta_err_all_%d",j));
						py_dphi_err_all[j]=(TH1D*) py_dphi_err[0][j]->Clone(Form("py_dphi_err_all_%d",j));
						for(int i=1;i<8; ++i){
								py_deta_all[j]->Add(py_deta[i][j]);
								py_dphi_all[j]->Add(py_dphi[i][j]);
								py_deta_err_all[j]->Add(py_deta_err[i][j]);
								py_dphi_err_all[j]->Add(py_dphi_err[i][j]);
						}
				}
		}

		void getPY_proj2(){
				for(int i=0;i<8; ++i){
						tmp = "Proj_dEta_pp_"+trk_tag[i]+"_"+trk_tag[i+1];
						py_deta[i][4]=(TH1D*)py_proj_f->Get(tmp);
						tmp = "Proj_dPhi_pp_"+trk_tag[i]+"_"+trk_tag[i+1];
						py_dphi[i][4]=(TH1D*)py_proj_f->Get(tmp);
						tmp = "dEta_Syst_pp_Cent0_Cent10_"+trk_tag[i]+"_"+trk_tag[i+1];
						py_deta_err[i][4]=(TH1D*)py_proj_f->Get(tmp);
						tmp = "dPhi_Syst_pp_Cent0_Cent10_"+trk_tag[i]+"_"+trk_tag[i+1];
						//cout<<tmp<<endl;
						py_dphi_err[i][4]=(TH1D*)py_proj_f->Get(tmp);
						for(int j=0; j<4; ++j){
								tmp = "Proj_dEta_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								py_deta[i][j]=(TH1D*)py_proj_f->Get(tmp);
								cout<<py_deta[i][j]->GetName()<<endl;
								tmp = "Proj_dPhi_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								py_dphi[i][j]=(TH1D*)py_proj_f->Get(tmp);
								tmp = "dEta_Syst_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								py_deta_err[i][j]=(TH1D*)py_proj_f->Get(tmp);
								tmp = "dPhi_Syst_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								py_dphi_err[i][j]=(TH1D*)py_proj_f->Get(tmp);
						}
				}
				for(int j=0; j<5; ++j){
						py_deta_all[j]=(TH1D*) py_deta[0][j]->Clone(Form("py_deta_all_%d",j));
						py_dphi_all[j]=(TH1D*) py_dphi[0][j]->Clone(Form("py_dphi_all_%d",j));
						py_deta_err_all[j]=(TH1D*) py_deta_err[0][j]->Clone(Form("py_deta_err_all_%d",j));
						py_dphi_err_all[j]=(TH1D*) py_dphi_err[0][j]->Clone(Form("py_dphi_err_all_%d",j));
						for(int i=1;i<8; ++i){
								py_deta_all[j]->Add(py_deta[i][j]);
								py_dphi_all[j]->Add(py_dphi[i][j]);
								py_deta_err_all[j]->Add(py_deta_err[i][j]);
								py_dphi_err_all[j]->Add(py_dphi_err[i][j]);
						}
				}
		}
		void getJS(){
				for(int i=0;i<9; ++i){
						tmp = "JetShape2_Yield_BkgSub_pTweightedInclusive_pp_"+
								trk_tag[i]+"_"+trk_tag[i+1];
						js_dr[i][4]=(TH1D*)js_f->Get(tmp);
						tmp = "Jet_Shape_SystErr_pp_"+
								trk_tag[i]+"_"+trk_tag[i+1];
						js_dr_err[i][4]=(TH1D*)js_f->Get(tmp);
						for(int j=0;j<4; ++j){
								tmp = "JetShape2_Yield_BkgSub_pTweightedInclusive_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+
										trk_tag[i]+"_"+trk_tag[i+1];
								js_dr[i][j]=(TH1D*)js_f->Get(tmp);
								tmp = "Jet_Shape_SystErr_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+
										trk_tag[i]+"_"+trk_tag[i+1];
								js_dr_err[i][j]=(TH1D*)js_f->Get(tmp);
						}
				}
						tmp = "Jet_Shape_SystErr_pp_TrkPt300_";
						js_dr_err_all[4]=(TH1D*)js_f->Get(tmp);
				for(int i=0;i<4; ++i){
						tmp = "JetShape2_Yield_BkgSub_pTweightedInclusive_"+cent_tag[i]+"_"+cent_tag[i+1]+"_"+
								trk_tag[9]+"_";
						js_dr_all[i]=(TH1D*)js_f->Get(tmp);

						tmp = "Jet_Shape_SystErr_"+cent_tag[i]+"_"+cent_tag[i+1]+"_"+
								trk_tag[9]+"_";
						js_dr_err_all[i]=(TH1D*)js_f->Get(tmp);
						//js_dr_err_all[i]=(TH1D*)js_dr_err[0][i]->Clone(Form("JetShape_PbPb_Total_Syst_err_%d",i));
				}
				/*
				tmp = "JetShape2_Yield_BkgSub_pTweightedInclusive_pp_"+
						trk_tag[9]+"_";
				js_dr_all[4]=(TH1D*)js_f->Get(tmp);

				js_dr_err_all[4]=(TH1D*)js_dr_err[0][4]->Clone("JetShape_pp_Total_Syst_err");
				for(int j=0; j<5; ++j){
						for(int i=1; i<9; ++i){
								js_dr_err_all[j]->Add(js_dr_err[i][j]);
						}
				}
				*/
		}


		TFile * bkgSub_f= TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/nominal/PbPb_Inclusive_Correlations.root");
		TH2D* py[9][4];
		TH2D* js[9][4];
		TH2D* raw_js[9][4];
		TH2D* raw_py[9][4];
		TH2D* mix_py[9][4];
		TH2D* mix_js[9][4];
		TH2D* bkg_js[9][4];
		TH2D* bkg_py[9][4];
		TH2D* py_nobkgSub[9][4];
		TH2D* js_nobkgSub[9][4];
		void get2DHist(bool ispp = false){
				float ncent = 4;
				if(ispp) ncent =1;
				for(int i=0;i<ncent; ++i){
						for(int j=0;j<8; ++j){
								tmp = "Yield_BkgSub_"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								py[j][i]=(TH2D*)bkgSub_f->Get(tmp);
								tmp = "Yield_BkgSub_pTweighted"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								js[j][i]=(TH2D*)bkgSub_f->Get(tmp);
								tmp = "Raw_Yield_"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								raw_py[j][i]=(TH2D*)bkgSub_f->Get(tmp);
								tmp = "Raw_Yield_pTweighted"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								raw_js[j][i]=(TH2D*)bkgSub_f->Get(tmp);
								tmp = "Mixed_Event_"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								mix_py[j][i]=(TH2D*)bkgSub_f->Get(tmp);
								tmp = "SummedBkg_"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								bkg_py[j][i]=(TH2D*)bkgSub_f->Get(tmp);
								tmp = "SummedBkg_pTweighted"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								bkg_js[j][i]=(TH2D*)bkgSub_f->Get(tmp);
								tmp = "Yield_PbPb_"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								if( ispp){ tmp = "Yield_pp_"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+\
										trk_tag[j]+"_"+trk_tag[j+1];
								}
								py_nobkgSub[j][i]=(TH2D*)bkgSub_f->Get(tmp);
								tmp = "Yield_pTweighted_PbPb_"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+
										trk_tag[j]+"_"+trk_tag[j+1];
								if( ispp ){ tmp = "Yield_pTweighted_pp_"+cent_tag[i]+"_"+cent_tag[i+1]+"_Pt100_Pt1000_"+\
										trk_tag[j]+"_"+trk_tag[j+1];
								}
								js_nobkgSub[j][i]=(TH2D*)bkgSub_f->Get(tmp);
						}
				}
		}
		TH1D* ppyield;
		TH1D* pbyield[4];
		TH1D* yield_diff[4];
		TH1D* ppyield_err;
		TH1D* pbyield_err[4];
		TH1D* yield_diff_err[4];

		Double_t bc [2][9][5];
		Double_t err[2][9][5];
		void prePlot_yield_integral(){
				getPY_proj();
				int x1 = py_deta[0][0]->GetXaxis()->FindBin(-1);
				int x2 = py_deta[0][0]->GetXaxis()->FindBin(1)-1;
				Double_t bins[] = {0., 0.7, 1, 2, 3, 4, 8, 12, 16, 20};
				ppyield = new TH1D("pp_yield", "",  9, bins);
				ppyield_err = new TH1D("pp_yield_err", "", 9, bins);
				for(int k=0; k<8; ++k){
						bc[0][k][4] = py_deta[k][4]->IntegralAndError(x1, x2, err[0][k][4], "width");
						ppyield->SetBinContent(k+2, bc [0][k][4]);
						ppyield->SetBinError  (k+2, err[0][k][4]);
						bc[1][k][4] = py_deta_err[k][4]->IntegralAndError(x1, x2, err[1][k][4], "width");
						ppyield_err->SetBinContent(k+2, bc [1][k][4]);
						ppyield_err->SetBinError  (k+2, err[1][k][4]);
				}
				for(int i=0; i<4; i++){
						pbyield[i] = (TH1D*) ppyield->Clone(Form("pb_yield_%d",i));
						pbyield_err[i] = (TH1D*) ppyield->Clone(Form("pb_yield_err_%d",i));
						for(int k=0; k<8; ++k){
								bc[0][k][i] = py_deta[k][i]->IntegralAndError(x1, x2, err[0][k][i], "width");
								pbyield[i]->SetBinContent(k+2, bc[0][k][i]);
								pbyield[i]->SetBinError  (k+2, err[0][k][i]);
								bc[1][k][i] = py_deta[k][i]->IntegralAndError(x1, x2, err[1][k][i], "width");
								pbyield_err[i]->SetBinContent(k+2, bc [1][k][i]);
								pbyield_err[i]->SetBinError  (k+2, err[1][k][i]);
						}
						yield_diff[i]=(TH1D*) pbyield[i]->Clone(Form("yield_diff_%d",i));
						yield_diff[i]->Add(ppyield,-1);
						yield_diff_err[i]=(TH1D*) pbyield_err[i]->Clone(Form("yield_diff_err_%d",i));
						yield_diff_err[i]->Add(ppyield_err,-1);
				}
		}
}

