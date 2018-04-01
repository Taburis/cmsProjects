
TString dataDumpPath = "/Users/tabris/cmsProjects/iJTC/dataSet/inclJTC/";
TString FigDumpPath  = "figs/";

TString trk_tag[] = {"TrkPt0p7", "TrkPt1", "TrkPt2","TrkPt3","TrkPt4","TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt999"};
TString trk2_tag[] = {"TrkPt07", "TrkPt1", "TrkPt2","TrkPt3","TrkPt4","TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt300"};
TString cent_tag[]= {"Cent0","Cent10", "Cent30", "Cent50", "Cent70","Cent100"};
TString cent2_tag[]= {"Cent0","Cent10", "Cent30", "Cent50", "Cent100"};

TString gengen_pythia_f = "/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/pp_Pythia6MC_GenGen_withMix_inclJetBinning_finalJFFs.root";
TString recgen_pythia_f = "/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/pp_Pythia6MC_RecoGen_withMix_inclJetBinning_finalJFFs.root";
TString recrec_pythia_f = "/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/pp_Pythia6MC_RecoReco_withMix_inclJetBinning_finalJFFs.root";
TString pp_f = "/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/pp_Data_withMix_inclJetBinning_finalJFFs.root";

TString gengen_pythia_norm_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/Pythia_GenJet_GenTrack_Inclusive_Correlations.root";
TString recgen_pythia_norm_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/Pythia_RecoJet_GenTrack_Inclusive_Correlations.root";
TString recrec_pythia_norm_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/Pythia_RecoJet_RecoTrack_Inclusive_Correlations.root";
//TString data_pp_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/pp_Inclusive_Correlations.root";
TString data_pp_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/pp_Inclusive_Correlations_kurtCorr.root";
TString data_pb_f2 = "/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/nominal/PbPb_Inclusive_Correlations.root";
TString data_pb_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/PbPb_Inclusive_Correlations.root";

TString gengen_hydjet_sube0_f = "/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/PbPb_5TeVMC_withMix_GenGen_sube0_CymbalTune_fullCymbalCorrs_inclJetBinning.root";
TString recgen_hydjet_sube0_f = "/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/PbPb_5TeVMC_withMix_RecoGen_sube0_CymbalTune_fullCymbalCorrs_inclJetBinning.root";

TString gengen_hydjet_sube0_nom_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/Hydjet_GenJet_GenTrack_Sube0_Inclusive_Correlations.root";
TString recgen_hydjet_sube0_nom_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/Hydjet_RecoJet_GenTrack_Sube0_Inclusive_Correlations.root";
TString recrec_hydjet_nom_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/Hydjet_RecoJet_RecoTrack_Inclusive_Correlations.root";

TString gengen_hydjet_norm_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/Hydjet_GenJet_GenTrack_Inclusive_Correlations.root";
TString recgen_hydjet_norm_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/Hydjet_RecoJet_GenTrack_Inclusive_Correlations.root";
TString pb_norm_f = "/Users/tabris/Research/HIN_5.02TeV/JetTrack2016/me_correct/PbPb_Inclusive_Correlations.root";

TH2D* hraw_sig[9][5];
TH2D* hraw_sig_pTweighted[9][5];
TH2D* hmixing[9][5];

TH1D* hjet_tmp[5];

TH1D* hjet[4];
TH2D* raw_sig[9][4];
TH2D* raw_sig_pTweighted[9][4];
TH2D* mixing [9][4];

void getSig_hallie(TFile *f, TH2D* h[9][4], TString cap, int ncent = 4 ) {
		TString tmp;
		for(int i=0; i<9; ++i){
				for(int j=0; j<ncent; ++j){
						tmp = cap+cent2_tag[j]+"_"+cent2_tag[j+1]+"_Pt100_Pt1000_"+trk2_tag[i]+"_"+trk2_tag[i+1];

						h[i][j] = (TH2D*) f->Get(tmp);
						h[i][j]->Scale(1.0/h[i][j]->GetXaxis()->GetBinWidth(1)/h[i][j]->GetYaxis()->GetBinWidth(1));
						cout<<h[i][j]->GetName()<<endl;
				}
		}
}

void dump_nominal_res(TString name, TString fname, bool ispp = 1, bool israw = 1){
		int ncent = ispp ? 1: 4;
		TH2D* h1[9][4];
		TH2D* h2[9][4];
		TH2D* h3[9][4];
		TFile *f = TFile::Open(fname);
		if(israw){
				getSig_hallie(f, h1, "Raw_Yield_", ncent);
				getSig_hallie(f, h2, "Raw_Yield_pTweighted", ncent);
		}else {
				getSig_hallie(f, h1, "Yield_BkgSub_", ncent);
				getSig_hallie(f, h2, "Yield_BkgSub_pTweighted", ncent);
		}
		getSig_hallie(f, h3, "SummedBkg_pTweighted_", ncent);
		auto wf = new TFile(dataDumpPath+name+"_JTCSignal.root", "recreate");
		wf->cd();
		for(int i=0; i<9; ++i){
				for(int j=0; j<ncent; ++j){
						h1[i][j]->SetName("signal_"+name+Form("_%d_%d",i,j));
						h1[i][j]->Write();
						h2[i][j]->SetName("signal_"+name+Form("_pTweighted_%d_%d",i,j));
						h2[i][j]->Write();
						h3[i][j]->SetName("bkg_"+name+Form("_pTweighted_%d_%d",i,j));
						h3[i][j]->Write();
				}
		}
		wf->Close();
}

void get2DInput_OneCent(TFile *f, TString cap, int j){
		TString tmp;
		TString jetname, jetname2;
		jetname = "_all_jets_corrpT";
		jetname2= "_hJet";
		tmp = cap+jetname+cent_tag[j]+"_"+cent_tag[j+1]\
			  +"_Pt100_Pt1000";
		hjet_tmp[j] = (TH1D*) f->Get(tmp);
		for(int i=0; i<9; ++i){
				tmp = cap+jetname2+"TrackSignalBackground"+cent_tag[j]+"_"+cent_tag[j+1]\
					  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
				//  cout<<tmp<<endl;
				hraw_sig[i][j] = (TH2D*) f->Get(tmp);
				tmp = cap+jetname2+"TrackSignalBackground_pTweighted"+cent_tag[j]+"_"+cent_tag[j+1]\
					  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
				hraw_sig_pTweighted[i][j] = (TH2D*) f->Get(tmp);
				tmp = cap+jetname2+"TrackME"+cent_tag[j]+"_"+cent_tag[j+1]\
					  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
				//      cout<<mixing[i][j]->GetName()<<endl;
				hmixing[i][j] = (TH2D*) f->Get(tmp);
		}
}


void transfer_tmp(){
		hjet[0] = (TH1D*) hjet_tmp[0]->Clone("jet_spectra_0");
		hjet[1] = (TH1D*) hjet_tmp[2]->Clone("jet_spectra_1");

		for(int i=0; i<9; ++i){
				raw_sig[i][0]=(TH2D*) hraw_sig[i][0]->Clone(Form("raw_sig_input_%d_0", i));
				raw_sig[i][1]=(TH2D*) hraw_sig[i][2]->Clone(Form("raw_sig_input_%d_1", i));
				raw_sig_pTweighted[i][0]=(TH2D*) hraw_sig_pTweighted[i][0]->Clone(Form("raw_sig_input_pTweighted_%d_0", i));
				raw_sig_pTweighted[i][1]=(TH2D*) hraw_sig_pTweighted[i][2]->Clone(Form("raw_sig_input_pTweighted_%d_1", i));
				double njet1 = hjet[0]->Integral();
				double njet2 = hjet[1]->Integral();
				raw_sig[i][0]->Scale(1.0/njet1);
				raw_sig[i][1]->Scale(1.0/njet2);
				raw_sig_pTweighted[i][0]->Scale(1.0/njet1);
				raw_sig_pTweighted[i][1]->Scale(1.0/njet2);
				mixing[i][0]=(TH2D*) hmixing[i][0]->Clone(Form("mixing_input_%d_0", i));
				mixing[i][1]=(TH2D*) hmixing[i][2]->Clone(Form("mixing_input_%d_1", i));
		}
}

void get2DInput_pp(TFile *f , TString cap ){
		for(int j=0; j<5; ++j){
				//cout<<hmixing[i][j]->GetName()<<endl;
				get2DInput_OneCent(f, cap, j);
		}
		cout<<"transfering"<<endl;
		transfer_tmp();
}

void get2DInput_Pb(TFile *f, TString cap){
		TString tmp;
		TString jetname, jetname2;
		jetname = "_all_jets_corrpT";
		jetname2= "_hJet";
		for(int i=0; i<9; ++i){
				for(int j=0; j<4; ++j){
						tmp = cap+jetname2+"TrackSignalBackground"+cent_tag[j]+"_"+cent_tag[j+1]\
							  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
						//  cout<<tmp<<endl;
						raw_sig[i][j] = (TH2D*) f->Get(tmp);
						tmp = cap+jetname2+"TrackSignalBackground_pTweighted"+cent_tag[j]+"_"+cent_tag[j+1]\
							  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
						raw_sig_pTweighted[i][j] = (TH2D*) f->Get(tmp);
						tmp = cap+jetname2+"TrackME"+cent_tag[j]+"_"+cent_tag[j+1]\
							  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
						//cout<<tmp<<endl;
						mixing[i][j] = (TH2D*) f->Get(tmp);
						//cout<<mixing[i][j]->GetName()<<endl;
				}
				int j=4;
				tmp = cap+jetname2+"TrackSignalBackground"+cent_tag[j]+"_"+cent_tag[j+1]\
					  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
				TH2D* htmp = (TH2D*) f->Get(tmp);
				cout<<htmp->GetName()<<endl;
				raw_sig[i][3]->Add(htmp);
				tmp = cap+jetname2+"TrackSignalBackground_pTweighted"+cent_tag[j]+"_"+cent_tag[j+1]\
					  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
				htmp = (TH2D*) f->Get(tmp);
				raw_sig_pTweighted[i][3]->Add(htmp);
				tmp = cap+jetname2+"TrackME"+cent_tag[j]+"_"+cent_tag[j+1]\
					  +"_Pt100_Pt1000_"+trk_tag[i]+"_"+trk_tag[i+1];
				h = (TH2D*) f->Get(tmp);
				mixing[i][3]->Add(htmp);
		}
		for(int j=0; j<5; ++j){
				tmp = cap+jetname+cent_tag[j]+"_"+cent_tag[j+1]\
					  +"_Pt100_Pt1000";
				cout<<tmp<<endl;
				hjet_tmp[j] = (TH1D*) f->Get(tmp);
				//cout<<hjet_tmp[j]->Integral()<<endl;
		}
		hjet_tmp[3]->Add(hjet_tmp[4]);
		for(int j=0; j<4; ++j){
				float njet = hjet_tmp[j]->Integral();
				for(int i=0; i<9; ++i){
//						cout<<raw_sig[i][j]->GetName()<<endl;
						raw_sig[i][j]->Scale(1.0/njet);
						raw_sig_pTweighted[i][j]->Scale(1.0/njet);
				}
		}
}

void get2DInput_GenGen_PYTHIA(){
		TFile *f = TFile::Open(gengen_pythia_f);
		get2DInput_pp(f, "GenJet_GenTrack");
}

void get2DInput_RecGen_PYTHIA(){
		TFile *f = TFile::Open(recgen_pythia_f);
		get2DInput_pp(f, "RecoJet_GenTrack");
}
void get2DInput_RecRec_PYTHIA(){
		TFile *f = TFile::Open(recrec_pythia_f);
		get2DInput_pp(f, "RecoJet_RecoTrack");
}
void get2DInput_pp_data(){
		TFile *f = TFile::Open(pp_f);
		get2DInput_pp(f, "Data");
}
void get2DInput_RecGen_sube0_HYDJET(){
		TFile *f = TFile::Open(recgen_hydjet_sube0_f);
		get2DInput_Pb(f, "RecoJet_GenTrack");
}
void get2DInput_GenGen_sube0_HYDJET(){
		TFile *f = TFile::Open(gengen_hydjet_sube0_f);
		get2DInput_Pb(f, "GenJet_GenTrack");
}


