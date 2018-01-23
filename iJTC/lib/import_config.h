
const int nPt = 8;
const int nCent = 2;
namespace input_raw2D{
		TString trk_tag[] = {"TrkPt07", "TrkPt1", "TrkPt2","TrkPt3","TrkPt4","TrkPt8","TrkPt12","TrkPt16","TrkPt999"};
		TString cent_tag[]= {"Cent0","Cent30","Cent100"};

		TFile *gengen_pb_sube0_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/GenGen_PbPb_5TeV_bJTC_csvCut_sube0.root");
		TFile *gengen_pb_sube0_trueB_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/GenGen_PbPb_5TeV_bJTC_sub0_trueB.root");
		TFile *recgen_pb_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/RecGen_PbPb_5TeV_bJTC_HPOn.root");
		TFile *recrec_pb_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/RecRec_PbPb_5TeV_bJTC_HPOn.root");
		TFile *gengen_pb_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/GenGen_PbPb_5TeV_bJTC_HPOn_CSVcutOnly.root");
		TFile *genrec_pb_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/GenRec_PbPb_5TeV_bJTC_HPOn_CSVcutOnly.root");
		TFile *gengen_pb_f_sub0_tagged_trueB_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/GenGen_PbPb_5TeV_bJTC_bTagged_trueB_sub0.root");

		TFile *gengen_pythia_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/GenGen_5TeV_bJTC_pythia.root");
		TH2D* raw_sig[8][2];
		TH2D* raw_sig_pTweighted[8][2];
		TH2D* mixing [8][2];

		void get2DInput(TFile *f, TString cap){
				float trkbinWidth[] = {0.3, 1, 1,1 , 4, 4, 4, 1};
				TH1D* h1; TString tmp; double njet;
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								tmp = cap+"_all_bjets_corrpT"+cent_tag[j]+"_"+cent_tag[j+1]\
									   +"_Pt120_Pt1000";
								h1 = (TH1D*) f->Get(tmp);
								njet = h1->Integral();
//								cout<<tmp<<endl;
								tmp = cap+"_hbJetTrackSignalBackground_"+cent_tag[j]+"_"+cent_tag[j+1]\
									   +"_Pt120_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
								raw_sig[i][j] = (TH2D*) f->Get(tmp);
								raw_sig[i][j]->Scale(1.0/njet/trkbinWidth[i]);
								tmp = cap+"_hbJetTrackSignalBackground_pTweighted"+cent_tag[j]+"_"+cent_tag[j+1]\
									   +"_Pt120_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
								raw_sig_pTweighted[i][j] = (TH2D*) f->Get(tmp);
								raw_sig_pTweighted[i][j]->Scale(1.0/njet/trkbinWidth[i]);
								tmp = cap+"_hbJetTrackME"+cent_tag[j]+"_"+cent_tag[j+1]\
									   +"_Pt120_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
						//		cout<<mixing[i][j]->GetName()<<endl;
								mixing[i][j] = (TH2D*) f->Get(tmp);
						}
				}
		}

		void clearInput(){
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								delete raw_sig[i][j];
								delete mixing[i][j];
						}
				}
		}

		void getMix(TFile *f, TString cap){
				TH1D* h1; TString tmp; 
				for(int i=0; i<nPt; ++i){
						for(int j=0; j<nCent; ++j){
								tmp = cap+"_hbJetTrackME"+cent_tag[j]+"_"+cent_tag[j+1]\
									   +"_Pt120_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
								mixing[i][j] = (TH2D*) f->Get(tmp);
						}
				}
		}

		void get2DInput_RecGen(){ get2DInput(recgen_pb_f, "RecoJet_GenTrack"); }
		void get2DInput_RecRec(){ get2DInput(recrec_pb_f, "RecoJet_RecoTrack"); }
		void get2DInput_GenGen(){ get2DInput(gengen_pb_f, "GenJet_GenTrack"); }
		void get2DInput_GenRec(){ get2DInput(genrec_pb_f, "GenJet_RecoTrack"); }
		void get2DInput_GenGen_sub0(){ get2DInput(gengen_pb_sube0_f, "GenJet_GenTrack"); }
		void get2DInput_GenGen_sub0_tagged_trueB(){ get2DInput(gengen_pb_f_sub0_tagged_trueB_f, "GenJet_GenTrack"); }
		void get2DInput_GenGen_sub0_trueB(){ get2DInput(gengen_pb_sube0_trueB_f, "GenJet_GenTrack"); }
		void getMix_GenGen(){ getMix(gengen_pb_f, "GenJet_GenTrack"); }
		void get2DPythiaInput_GenGen(){ get2DInput(gengen_pythia_f, "GenJet_GenTrack"); }
//		void get2DInput_RecGen_sub0(){ get2DInput(recgen_pb_sube0_f, "RecoJet_GenTrack"); }
}

namespace signal2D {
		TH2D* sig[nPt][nCent];
		TH2D* raw[nPt][nCent];
		TH2D* mixing[nPt][nCent];
		TH2D* bkg[nPt][nCent];

		TFile *recgen_pb_f;
		TFile *recrec_pb_f;
		TFile *gengen_pb_f;
		TFile *genrec_pb_f;
		TFile *recgen_pb_sub0_f;
		TFile *gengen_pb_sub0_f;
		TFile *recgen_pb_nsub0_f;

		void loadFile(){
				cout<<"reading the file.. "<<endl;
				recgen_pb_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/rec_gen_JTCSignal.root");
				recrec_pb_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/rec_rec_JTCSignal.root");
				gengen_pb_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/gen_gen_JTCSignal.root");
				genrec_pb_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/gen_rec_JTCSignal.root");
				recgen_pb_sub0_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/rec_gen_sub0_JTCSignal.root");
				gengen_pb_sub0_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/gen_gen_sub0_JTCSignal.root");
				recgen_pb_nsub0_f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/correlation/rec_gen_nsub0_JTCSignal.root");
		}
		void getHist(TFile *f, TString name){
				TString tmp;
				for(int i=0; i<7; ++i){
						for(int j=0; j<2; ++j){
								tmp = "signal_"+name+Form("_%d_%d",i,j);
								//cout<<tmp<<endl;
								sig[i][j]=(TH2D*) f->Get(tmp);
								tmp = "smoothed_mixing_"+name+Form("_%d_%d",i,j);
								mixing[i][j]=(TH2D*) f->Get(tmp);
								tmp = "raw_"+name+Form("_%d_%d",i,j);
								raw[i][j]=(TH2D*) f->Get(tmp);
								tmp = "bkg_"+name+Form("_%d_%d",i,j);
								bkg[i][j]=(TH2D*) f->Get(tmp);
						}
				}
		}
		void getRecRec(){ getHist(recrec_pb_f, "rec_rec");}
		void getRecGen(){ getHist(recgen_pb_f, "rec_gen");}
}

