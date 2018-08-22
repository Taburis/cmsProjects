
#ifndef mixingTree_H
#define mixingTree_H
#include "TChain.h"
#include "TMath.h"
class mixingTree{
	public : mixingTree(TString name){t=new TChain(name);};
		 void Loop(TString);
		 void make_jet_miniTree(TString, TString);
		 void init(bool ismc, bool ishi, bool isKurtSkim);
		 bool eventCuts();

		 TChain *t;
		 bool isMC, isHI;

		 float pthat;
		 float vz ;
		 int hiBin;
		 int b_mode;

		 vector<float> *refpt=0, *jtpt=0, *jteta=0, *jtphi=0, *discr_csv = 0, *trackMax = 0;// need jet info. for tk correction
		 vector<int> *pdg=0, *status=0;
		 vector<float>  *gen_jtpt=0, *gen_jteta=0, *gen_jtphi=0;
		 vector<float> *pt=0, *eta=0, *phi=0, *chg=0;
		 vector<float> *trkPt=0, *trkEta=0, *trkPhi=0, *trkChi2=0,  *trkPtError=0;
		 vector<float> *trkDz=0, *trkDzError=0, *trkDxy=0, *trkDxyError=0, *pfEcal=0, *pfHcal=0;
		 vector<int> *trkNlayer=0, *trkNdof=0, *trkNHit=0, *flavorForB=0;
		 vector<bool> *highPurity=0; 
		 int  HBHEFilter=0, collisionEventSelection = 0, pvFilter=0, phfCoincFilter=0;

		 // the parameters used below are the configuration for the loop and scan, not for the correlation
		 bool doCaloJet = 0;
		 bool doJFF = 0;
		 bool doJetScan = 0;
		 bool doTkScan = 0;
		 int npthat = 5;
		 float *bin_pthat; //pthat bin need to be specify by hands
		 int nvz = 200;
		 float vz_min = -20, vz_max = 20;
		 int njtpt = 200; // has to be fine to get the JES resolution
		 float jtpt_min = 50, jtpt_max = 450;//50 for the jet without matching
		 int njtphi = 36, njteta=  300;
		 float jteta_min = -3, jteta_max = 3, jtphi_min = -TMath::Pi(), jtphi_max = TMath::Pi();
		 int ntkpt = 30, ntketa = 300, ntkphi = 36;
		 float tkpt_min = 0, tkpt_max = 30, tketa_min = 3, tketa_max = 3, tkphi_min = -TMath::Pi(), tkphi_max = TMath::Pi();

};

void mixingTree::init(bool ismc, bool ishi, bool isKurtSkim){
	isMC = ismc; isHI = ishi;

	t->SetBranchStatus("*", 0);
	if(ishi){
		t->SetBranchStatus("hiBin",1);
		t->SetBranchAddress("hiBin",&hiBin);
		t->SetBranchStatus("pprimaryVertexFilter",1);
		t->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
		t->SetBranchStatus("pcollisionEventSelection",1);
		t->SetBranchAddress("pcollisionEventSelection",&collisionEventSelection);
		t->SetBranchStatus("phfCoincFilter3",1);
		t->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
	}
	else {
		if(isKurtSkim){
			t->SetBranchStatus("pprimaryVertexFilter",1);
			t->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
		}
		else {
			t->SetBranchStatus("pPAprimaryVertexFilter",1);
			t->SetBranchAddress("pPAprimaryVertexFilter",&pvFilter);
		}
	}
	t->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
	t->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHEFilter);
	t->SetBranchStatus("vz",1);
	t->SetBranchAddress("vz",&vz);

	if(ismc){
		t->SetBranchStatus("flavor_b_mode",1);
		t->SetBranchAddress("flavor_b_mode",&b_mode);
		t->SetBranchStatus("pthat",1);
		t->SetBranchAddress("pthat",&pthat);
		t->SetBranchStatus("pt",1);
		t->SetBranchAddress("pt",&pt);
		t->SetBranchStatus("eta",1);
		t->SetBranchAddress("eta",&eta);
		t->SetBranchStatus("phi",1);
		t->SetBranchAddress("phi",&phi);
		t->SetBranchStatus("chg",1);
		t->SetBranchAddress("chg",&chg);
		t->SetBranchStatus("genpt",1);
		t->SetBranchAddress("genpt",&gen_jtpt);
		t->SetBranchStatus("geneta",1);
		t->SetBranchAddress("geneta",&gen_jteta);
		t->SetBranchStatus("genphi",1);
		t->SetBranchAddress("genphi",&gen_jtphi);
		t->SetBranchStatus("pdg",1);
		t->SetBranchAddress("pdg",&pdg);
		t->SetBranchStatus("status",1);
		t->SetBranchAddress("status",&status);
	}
	if(doCaloJet){
		if(doJFF){	
			t->SetBranchStatus ("calo_corrpt",1);
			t->SetBranchAddress("calo_corrpt",&jtpt);
		}else {
			t->SetBranchStatus ("calo_jtpt",1);
			t->SetBranchAddress("calo_jtpt",&jtpt);
		}
		t->SetBranchStatus ("calo_refpt",1);
		t->SetBranchAddress("calo_refpt",&refpt);
		t->SetBranchStatus ("calo_jteta",1);
		t->SetBranchAddress("calo_jteta" ,&jteta);
		t->SetBranchStatus ("calo_jtphi",1);
		t->SetBranchAddress("calo_jtphi" ,&jtphi);
		t->SetBranchStatus ("calo_discr_csvV2" ,1);
		t->SetBranchAddress("calo_discr_csvV2" ,&discr_csv);
		t->SetBranchStatus ("calo_trackMax" ,1);
		t->SetBranchAddress("calo_trackMax" ,&trackMax);
		if(ismc){
			t->SetBranchStatus ("calo_refparton_flavorForB" ,1);
			t->SetBranchAddress("calo_refparton_flavorForB" ,&flavorForB);
		}
	}else {
		if(doJFF){	
			t->SetBranchStatus ("pf_corrpt",1);
			t->SetBranchAddress("pf_corrpt",&jtpt);
		} else {
			t->SetBranchStatus ("pf_jtpt",1);
			t->SetBranchAddress("pf_jtpt",&jtpt);
		}
		t->SetBranchStatus ("pf_refpt",1);
		t->SetBranchAddress("pf_refpt",&refpt);
		t->SetBranchStatus ("pf_jteta",1);
		t->SetBranchAddress("pf_jteta" ,&jteta);
		t->SetBranchStatus ("pf_jtphi",1);
		t->SetBranchAddress("pf_jtphi" ,&jtphi);
		t->SetBranchStatus ("pf_discr_csvV2" ,1);
		t->SetBranchAddress("pf_discr_csvV2" ,&discr_csv);
		t->SetBranchStatus ("pf_trackMax" ,1);
		t->SetBranchAddress("pf_trackMax" ,&trackMax);
		if(ismc){
			t->SetBranchStatus ("pf_refparton_flavorForB" ,1);
			t->SetBranchAddress("pf_refparton_flavorForB" ,&flavorForB);
		}
	}
	t->SetBranchStatus("highPurity",1);
	t->SetBranchAddress("highPurity",&highPurity);
	t->SetBranchStatus("trkPt",1);
	t->SetBranchAddress("trkPt",&trkPt);
	t->SetBranchStatus("trkEta",1);
	t->SetBranchAddress("trkEta",&trkEta);
	t->SetBranchStatus("trkPhi",1);
	t->SetBranchAddress("trkPhi",&trkPhi);
	t->SetBranchStatus("trkDz",1);
	t->SetBranchAddress("trkDz",&trkDz);
	t->SetBranchStatus("trkDxy",1);
	t->SetBranchAddress("trkDxy",&trkDxy);
	t->SetBranchStatus("trkDxyError",1);
	t->SetBranchAddress("trkDxyError",&trkDxyError);
	t->SetBranchStatus("trkDzError",1);
	t->SetBranchAddress("trkDzError",&trkDzError);
	t->SetBranchStatus("trkPtError",1);
	t->SetBranchAddress("trkPtError",&trkPtError);
	t->SetBranchStatus("trkNdof",1);
	t->SetBranchAddress("trkNdof",&trkNdof);
	t->SetBranchStatus("trkNHit",1);
	t->SetBranchAddress("trkNHit",&trkNHit);
	t->SetBranchStatus("trkNlayer",1);
	t->SetBranchAddress("trkNlayer",&trkNlayer);
	t->SetBranchStatus("trkChi2",1);
	t->SetBranchAddress("trkChi2",&trkChi2);
	t->SetBranchStatus("pfEcal",1);
	t->SetBranchAddress("pfEcal",&pfEcal);
	t->SetBranchStatus("pfHcal",1);
	t->SetBranchAddress("pfHcal",&pfHcal);
}

bool mixingTree::eventCuts(){
	//	cout<<"vz = "<<vz<<", pthat = "<<pthat<<", HBHEFilter="<<HBHEFilter<<", eventSelection= "<<collisionEventSelection<<endl;
	if(TMath::Abs(vz)>=15) return 0;
	if(isMC && pthat<80) return 0;
	if(!isHI){
		if(!pvFilter || !HBHEFilter){
			//cout<<"!"<<endl;
			return 0;
		}
	}else {
		if(!phfCoincFilter || !collisionEventSelection || !pvFilter || !HBHEFilter) return 0;
	}
	return 1;
}

#endif
