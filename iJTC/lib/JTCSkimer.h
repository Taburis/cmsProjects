#include <stdlib.h>
#include <fstream>
#include <vector>
#include "TTree.h"
#include "histManager.h"
#include "TFile.h"
#include "TMath.h"
#include "TH1.h"
#include "TH2.h"
#include "TF1.h"
#include "TString.h"
#include "xAxis.h"
#include "mixingTree.h"
#include "TRandom3.h"

struct histCase{
	//if want to add any hist, need to add it in the quickHistReg as well, and add the filling in the fillCase
	TH2D** sig;
	TH2D** sig_raw;
	TH2D** sig_pTweighted;
	TH2D** sig_pTweighted_raw;
	TH2D** mixing;
	TH2D** mixing_raw;
	TH1D** jet_corrpt;
	TH1D** jet_eta;
	TH1D** jet_phi;
	TH1D** jt_wta_pt;
	TH1D** jt_wta_deta;
	TH1D** jt_wta_dphi;
};

const int nPt = 6;
const int nCent = 2;
TF1 *ppVz = new TF1("newppVz","gaus",-15,15);
TF1 *pp_vz_bMC = new TF1("pp_vz_bMC","[0]*exp(([1]*x*x+[2]*x+[3])*0.5)",-15,15);
TF1* fit_cen = new TF1("fcent1","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[7]*exp([5]+[6]*x)",0,180);
TF1* fit_vz_cymbal = new TF1("fit_vz_cymbal","gaus(0)/gaus(3)",-30,30);
void config(){
	pp_vz_bMC->SetParameters(0.922211, 0.00525714, -0.0287517, 0.012718);
	ppVz->SetParameter(0,1.10477);
	ppVz->SetParameter(1,2.52738);
	ppVz->SetParameter(2,1.30296e1);
	fit_cen->SetParameters(4.40810, -7.75301e-02, 4.91953e-04, -1.34961e-06, 1.44407e-09, -160, 1, 3.68078e-15);
	fit_vz_cymbal->SetParameters(0.08,0.44,5.12,0.08,3.25,5.23);
}

class JTCSkimer{
	public : 
		JTCSkimer(){ 
		};
		bool doGenJet=1, doGenTrack=1;
		void connectTree(TTree* ct, bool isMC, bool isHI, bool doKurtSkim);
		void initTree(TTree* ct, bool isMC, bool isHI);
		void addMixingTree(TString fname);
		void configMixing();
		bool isHI, isMC, doKurtSkim, domixing = 0;
		void eventWeight(float pthat, float vz, float hibin);
		void fillEventInfo();
		void findWTAaxis();
		void addGenJet   (TTree* t);
		void addRecoJet  (TTree* t);
		void addGenTrack (TTree* t);
		void addRecoTrack(TTree* t);
		int trkCuts(float trkPt, float trkPtError, float trkEta, float trkChi2, int highPurity, int trkNHit, float trkNdof, float trkNlayer, float pfHcal, float pfEcal, float trkdz,float trkdxy, float trkdzerror, float trkdxyerror);
		bool eventCuts();
		bool recoJetCut(int j);
		bool btagger(int j);
		void quickHistReg(histManager *h, TString cap, histCase & hc);
		void configHist();
		void loop();
		void write(TString name);
		double findDr(double eta1, double phi1, double eta2, double phi2);
		int jetMatch(int j);
		void trackLoop();
		void trackLoopForMixing();
		void fillCase(histCase &hc, bool ismixing = 0, float weight = 1);
		void fillCaseMix(histCase &hc);
		void fillHist(bool isMixing = 0);
		void fillAllJetInfo();
		void fillJetInfo(histCase &hc, float weight);
		float findDrmin(float eta, float phi, vector<float> *jetpt, vector<float> *jeteta, vector<float> *jetphi);
		void mixingSection(Long64_t voidIndex);
	public: 
		TRandom3 rand3;
		float jer_mean=0, jer_sigma=1;
		histManager *hm;
		//control config ----------------------------------
		bool doSmear = 0;
		bool isbMC = 0;
		bool doCaloJet = 0;
		bool doWTAaxis = 0;
		bool doDCATkCut= 0;
		float csvCut = 0.9;
		//mixing config------------------------------
		float trketamaxcut = 2.;
		float trkPtcut = 1.0;
		float jetetacut= 1.6;
		float jetptcut = 120;
		float vzcut = 15;
		bool do_mixing=0;
		int nPerTrigger=25;
		int nvz_mix, ncent_mix;
		//-----------------------------------------
		float ptbins[nPt+1] = {1, 2, 3, 4, 8, 12, 400};
		TString ptLabel[nPt+1] = {"Trk1", "Trk2", "Trk3", "Trk4", "Trk8", "trk12", "Trk400"};
		float centbins[nCent+1] = {0, 60, 200};
		TString centLabel[nCent+1] = {"Cent0", "Cent30", "Cent100"};
		TTree* t;
		// event branch
		float pthat;
		float vz ;
		int hiBin=0, b_mode = 0;
		mixingTree *mt = 0;
		// jet branch
		vector<float> *gen_jtpt=0, *gen_jteta=0, *gen_jtphi=0, *jtpt=0, *jteta=0, *jtphi=0;
		vector<float> *discr_csv=0;
		vector<int>	*flavorForB=0;
		// gen jet branch

		//gen particle branch
		vector<float> *pt=0, *eta=0, *phi=0, *chg=0;
		vector<int>	*sube=0;

		//track branch
		vector<float> *trkPt=0, *trkEta=0, *trkPhi=0, *trkChi2=0, *trkPtError=0;
		vector<float> *trkDz=0, *trkDzError=0, *trkDxy=0, *trkDxyError=0, *pfEcal=0, *pfHcal=0;
		vector<int> *trkNlayer=0, *trkNdof=0, *trkNHit=0, *ncs=0;
		vector<bool> *highPurity=0; 

		int  HBHEFilter=0, collisionEventSelection = 0, phfCoincFilter=0, pvFilter=0;

		//variables using in the correlations
		bool dcaCut;
		float selectedJt_pt;
		float selectedJt_eta;
		float selectedJt_phi;
		float selectedTrack_pt;
		float selectedTrack_eta;
		float selectedTrack_phi;
		float genJetAxis_pt ;			
		float genJetAxis_eta;			
		float genJetAxis_phi;	
		float weights;
		bool foundTrueBJet, foundTaggedBJet, foundInclJet, doJetAxisResolution = 0;
		float trkCorr=1, JESCorr=1;
		int gen_sube;
		xAxis *ptax, *centax;
		float wta_pt, wta_eta, wta_phi;
		// end of the variables declaration	
		std::vector<unsigned int >** mixTable;

		TH1D* hvz, *hpthat, *hcent;
		// if add a new histCase here, need to add it in the fillJetInfo function, register in the configHist
		// and add it into the fillHist
		histCase inclCase;
		histCase contCase;
		histCase taggedBCase;
		histCase trueBCase;
		histCase taggedTrueBCase;
};

void JTCSkimer::addGenJet(TTree *t){
	t->SetBranchStatus("genpt",1);
	t->SetBranchStatus("geneta",1);
	t->SetBranchStatus("genphi",1);
	t->SetBranchAddress("genpt",&gen_jtpt);
	t->SetBranchAddress("geneta",&gen_jteta);
	t->SetBranchAddress("genphi",&gen_jtphi);
}

void JTCSkimer::addRecoJet(TTree *t){
	if( doCaloJet){
		t->SetBranchStatus ("calo_corrpt",1);
		t->SetBranchAddress("calo_corrpt",&jtpt);
		t->SetBranchStatus ("calo_jteta",1);
		t->SetBranchAddress("calo_jteta" ,&jteta);
		t->SetBranchStatus ("calo_jtphi",1);
		t->SetBranchAddress("calo_jtphi" ,&jtphi);
		t->SetBranchStatus ("calo_discr_csvV1" ,1); t->SetBranchAddress("calo_discr_csvV1" ,&discr_csv);
		//		t->SetBranchStatus ("calo_discr_csvV2" ,1);
		//		t->SetBranchAddress("calo_discr_csvV2" ,&discr_csv);
	}else {
		t->SetBranchStatus ("pf_jtpt",1);
		t->SetBranchAddress("pf_jtpt",&jtpt);
		//t->SetBranchStatus ("pf_corrpt",1);
		//t->SetBranchAddress("pf_corrpt",&jtpt);
		t->SetBranchStatus("pf_jteta",1);
		t->SetBranchAddress("pf_jteta" ,&jteta);
		t->SetBranchStatus("pf_jtphi",1);
		t->SetBranchAddress("pf_jtphi" ,&jtphi);
		if(doKurtSkim){
			t->SetBranchStatus("pf_discr_csvV1" ,1);
			t->SetBranchAddress("pf_discr_csvV1" ,&discr_csv);
		} else {
			t->SetBranchStatus("pf_discr_csvV2" ,1);
			t->SetBranchAddress("pf_discr_csvV2" ,&discr_csv);
		}
	}
}

void JTCSkimer::addGenTrack(TTree *t){
	t->SetBranchStatus("pt",1);
	t->SetBranchAddress("pt",&pt);
	t->SetBranchStatus("eta",1);
	t->SetBranchAddress("eta",&eta);
	t->SetBranchStatus("phi",1);
	t->SetBranchAddress("phi",&phi);
	t->SetBranchStatus("chg",1);
	t->SetBranchAddress("chg",&chg);
}
void JTCSkimer::addRecoTrack(TTree *t){
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


void JTCSkimer::initTree(TTree *ct, bool ismc, bool ishi){
	ct->SetBranchStatus("*", 0);
	if(ismc){
		ct->SetBranchStatus("pthat",1);
		ct->SetBranchAddress("pthat",&pthat);
		ct->SetBranchStatus("flavor_b_mode",1);
		ct->SetBranchAddress("flavor_b_mode",&b_mode);
		addGenJet(ct);
		addGenTrack(ct);
		if(ishi) {
			ct->SetBranchStatus("sube",1);
			ct->SetBranchAddress("sube",&sube);
		}
		if(doCaloJet){
			t->SetBranchStatus ("calo_refparton_flavorForB" ,1);
			t->SetBranchAddress("calo_refparton_flavorForB" ,&flavorForB);
		} else { 
			t->SetBranchStatus("pf_refparton_flavorForB" ,1);
			t->SetBranchAddress("pf_refparton_flavorForB" ,&flavorForB);
		}
	}
	if(ishi){
		ct->SetBranchStatus("hiBin",1);
		ct->SetBranchAddress("hiBin",&hiBin);
		ct->SetBranchStatus("phfCoincFilter3",1);
		ct->SetBranchAddress("phfCoincFilter3",&phfCoincFilter);
		ct->SetBranchStatus("pcollisionEventSelection",1);
		ct->SetBranchAddress("pcollisionEventSelection",&collisionEventSelection);
		ct->SetBranchStatus("pprimaryVertexFilter",1);
		ct->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
		ct->SetBranchStatus("pprimaryVertexFilter",1);
		ct->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
	}
	else {
		if(doKurtSkim){
			ct->SetBranchStatus("pprimaryVertexFilter",1);
			ct->SetBranchAddress("pprimaryVertexFilter",&pvFilter);
		} else {
			ct->SetBranchStatus("pPAprimaryVertexFilter",1);
			ct->SetBranchAddress("pPAprimaryVertexFilter",&pvFilter);
		}
	}
	//	ct->SetBranchStatus("pf_nPFpartGT2_id1",1);
	//	ct->SetBranchAddress("pf_nPFpartGT2_id1",&ncs);
	ct->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
	ct->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHEFilter);
	ct->SetBranchStatus("vz",1);
	ct->SetBranchAddress("vz",&vz);
	addRecoJet  (ct);
	addRecoTrack(ct);
}

void JTCSkimer::connectTree(TTree *ct, bool ismc, bool ishi, bool dokurtSkim){
	isMC = ismc; isHI = ishi, doKurtSkim = dokurtSkim;
	t=ct;
	if( !isMC){ doGenJet = 0; doGenTrack=0;}
	initTree(t, isMC, isHI);
}

void JTCSkimer::quickHistReg(histManager *h, TString cap, histCase &hc){
	int nHistoBinsX = 500;
	int nHistoBinsY = 200;
	hc.sig= new TH2D*[nPt*nCent];
	hc.sig_pTweighted= new TH2D*[nPt*nCent];
	hc.mixing= new TH2D*[nPt*nCent];
	hc.sig_raw= new TH2D*[nPt*nCent];
	hc.sig_pTweighted_raw= new TH2D*[nPt*nCent];
	hc.mixing_raw= new TH2D*[nPt*nCent];
	TString tmp, name;
	if(isMC){
		if(doGenJet ) name= cap+"_GenJet_";
		else name=cap+"_RecoJet_";
		if(doGenTrack) name= name+"GenTrack";
		else name = name +"RecoTrack";
	} else name = cap+"_Data";
	for(int i=0; i<nPt; ++i){
		for(int j=0; j<nCent; ++j){
			tmp = centLabel[j]+"_"+centLabel[j+1]+"_"+ptLabel[i]+"_"+ptLabel[i+1];
			hc.sig[i+j*nPt] = hm->regHist<TH2D>(name+Form("_%d_%d",i, j), tmp,
					nHistoBinsX,-5,5,nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
			hc.sig_pTweighted[i+j*nPt] = hm->regHist<TH2D>(name+Form("_pTweighted_%d_%d",i, j), tmp,
					nHistoBinsX,-5,5, nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
			hc.mixing[i+j*nPt] = hm->regHist<TH2D>(name+Form("_mixing_%d_%d",i, j), tmp,
					nHistoBinsX,-5,5, nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
			hc.sig_raw[i+j*nPt] = hm->regHist<TH2D>(name+Form("_noCorr_%d_%d",i, j), tmp,
					nHistoBinsX,-5,5,nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
			hc.sig_pTweighted_raw[i+j*nPt] = hm->regHist<TH2D>(name+Form("_pTweighted_noCorr_%d_%d",i, j), tmp,
					nHistoBinsX,-5,5, nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
			hc.mixing_raw[i+j*nPt] = hm->regHist<TH2D>(name+Form("_mixing_noCorr_%d_%d",i, j), tmp,
					nHistoBinsX,-5,5, nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
		}
	}

	const float newbin [21] = {110, 120, 136, 152, 168, 184, 200, 216, 232, 248, 264, 280, 296, 312, 328, 344, 360,
		380, 400, 432, 500};
	int nbin = 20;
	name=cap;
	hc.jt_wta_pt = new TH1D*[nCent];
	hc.jt_wta_deta = new TH1D*[nCent];
	hc.jt_wta_dphi = new TH1D*[nCent];
	hc.jet_corrpt = new TH1D*[nCent];
	hc.jet_eta = new TH1D*[nCent];
	hc.jet_phi = new TH1D*[nCent];
	const float wta_ptbin [17] = {1, 2, 3, 4,5 ,6,7 ,8, 10, 12, 14, 16, 18, 20, 50, 100, 500};
	for(int j=0; j<nCent; ++j){
		tmp = centLabel[j]+" to "+centLabel[j+1];
		//				hc.jet_corrpt[j] = hm->regHist<TH1D>(name+Form("_corrpt_%d",j), tmp, 50, 0, 500.0);
		hc.jet_corrpt[j] = hm->regHist<TH1D>(name+Form("_corrpt_%d",j), tmp, nbin, newbin);
		hc.jet_eta[j] = hm->regHist<TH1D>(name+Form("_eta_%d",j), tmp, 100, -2.0, 2.0);
		hc.jet_phi[j] = hm->regHist<TH1D>(name+Form("_phi_%d",j), tmp, 72, -TMath::Pi(), TMath::Pi());
		if(doWTAaxis){
			hc.jt_wta_pt[j] = hm->regHist<TH1D>(name+Form("_wta_pt_%d",j), tmp, 16, wta_ptbin);
			hc.jt_wta_deta[j] = hm->regHist<TH1D>(name+Form("_wta_deta_%d",j), tmp, 20, -0.5, .5);
			hc.jt_wta_dphi[j] = hm->regHist<TH1D>(name+Form("_wta_dphi_%d",j), tmp, 20, -TMath::Pi()/4, TMath::Pi()/4);
		}
	}
}

void JTCSkimer::configHist(){
	hm = new histManager();

	ptax = new xAxis(nPt, ptbins);
	centax= new xAxis(nCent, centbins);
	hvz = hm->regHist<TH1D>("vzInfo", "", 200, -20, 20);
	hcent = hm->regHist<TH1D>("centInfo","",  50, 0, 200);
	if(isMC) hpthat = hm->regHist<TH1D>("pthatInfo", "", 100, 0, 400);
	quickHistReg(hm, "inclJet", inclCase);
	quickHistReg(hm, "taggedBJet", taggedBCase);
	if(isMC) {
		quickHistReg(hm, "contJet", contCase);
		quickHistReg(hm, "trueBJet", trueBCase);
		quickHistReg(hm, "taggedTrueBJet", taggedTrueBCase);
	}
	hm->sumw2();
}

void JTCSkimer::write(TString name){
	TFile *wf = new TFile(name, "recreate");
	wf->cd();
	hm->write();
	wf->Close();
}

void JTCSkimer::fillAllJetInfo(){
	float extra = 1;
	fillJetInfo(inclCase, weights*extra);
	if(foundTaggedBJet)fillJetInfo(taggedBCase, weights*extra);
	if(isMC) {
		if(foundTaggedBJet && !foundTrueBJet) fillJetInfo(contCase, weights*extra);
		if(foundTrueBJet)fillJetInfo(trueBCase, weights);
		if(foundTaggedBJet&&foundTrueBJet)fillJetInfo(taggedTrueBCase, weights*extra);
	}
}

void JTCSkimer::fillJetInfo(histCase &hc, float weight){
	if(isHI){
		int j = centax->findBin(hiBin);
		hc.jet_corrpt[j]->Fill(selectedJt_pt, weight);
		hc.jet_eta[j]->Fill(selectedJt_eta,   weight);
		hc.jet_phi[j]->Fill(selectedJt_phi,   weight);
	}else {
		for(int j=0; j<nCent; ++j){
			hc.jet_corrpt[j]->Fill(selectedJt_pt, weight);
			hc.jet_eta[j]->Fill(selectedJt_eta,   weight);
			hc.jet_phi[j]->Fill(selectedJt_phi,   weight);
			if(doWTAaxis){
				hc.jt_wta_pt  [j]->Fill(wta_pt, weight);
				hc.jt_wta_deta[j]->Fill(wta_eta-selectedJt_eta, weight);
				hc.jt_wta_dphi[j]->Fill(wta_phi-selectedJt_phi, weight);
			}
		}
	}
}

void JTCSkimer::fillCase(histCase &hc, bool ismixing, float weight){
	float deta;
	float dphi;
	if(doWTAaxis){
		deta = wta_eta - selectedTrack_eta;
		dphi = wta_phi - selectedTrack_phi;
	}else {
		deta = selectedJt_eta - selectedTrack_eta;
		dphi = selectedJt_phi - selectedTrack_phi;
	}
	while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
	while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}

	int i = ptax->findBin(selectedTrack_pt);
	if(!isHI){ 
		for(int j=0; j<nCent; ++j){
			if(ismixing){
				//cout<<"filling mixing"<<endl;
				hc.mixing[i+j*nPt]->Fill(deta, dphi, weight*trkCorr);
				hc.mixing_raw[i+j*nPt]->Fill(deta, dphi, weight);
			} else {
				hc.sig[i+j*nPt]->Fill(deta, dphi, weight*trkCorr);
				hc.sig_pTweighted[i+j*nPt]->Fill(deta, dphi, weight*selectedTrack_pt*trkCorr);
				hc.sig_raw[i+j*nPt]->Fill(deta, dphi, weight);
				hc.sig_pTweighted_raw[i+j*nPt]->Fill(deta, dphi, weight*selectedTrack_pt);
			}
		}
	} else {
		//special case, noCorr will be fill sube0 and ordinary will be filled by nsube0 for gen
		int j = centax->findBin(hiBin);
		//cout<<j<<", hiBin: "<<hiBin<<endl;
		if(ismixing){
			//cout<<"filling mixing"<<endl;
			if(doGenTrack && gen_sube != 0 )hc.mixing[i+j*nPt]->Fill(deta, dphi, weight);
			else if(doGenTrack) hc.mixing_raw[i+j*nPt]->Fill(deta, dphi, weight);
			else hc.mixing_raw[i+j*nPt]->Fill(deta, dphi, weight);
		} else {
			if(doGenTrack){
				if(gen_sube != 0){
					hc.sig[i+j*nPt]->Fill(deta, dphi, weight);
					hc.sig_pTweighted[i+j*nPt]->Fill(deta, dphi, weight*selectedTrack_pt);
				} else {
					hc.sig_raw[i+j*nPt]->Fill(deta, dphi, weight);
					hc.sig_pTweighted_raw[i+j*nPt]->Fill(deta, dphi, weight*selectedTrack_pt);
				}
			} else {
				hc.sig_raw[i+j*nPt]->Fill(deta, dphi, weight);
				hc.sig_pTweighted_raw[i+j*nPt]->Fill(deta, dphi, weight*selectedTrack_pt);
			}
		}
	}
}

void JTCSkimer::fillEventInfo(){
	hvz->Fill(vz, weights);
	if(isMC) hpthat->Fill(pthat, weights);
}

void JTCSkimer::addMixingTree(TString name){
	if(mt==0) mt = new mixingTree(t->GetName());
	mt->doCaloJet = doCaloJet;
	cout<<t->GetName()<<" has been created"<<endl;
	mt->t->Add(name);
}

void JTCSkimer::fillHist(bool isMixing){
	float extra= 1;
	fillCase(inclCase, isMixing, weights*extra);
	//if(!dcaCut) fillCase(inclCase, isMixing);
	if(foundTaggedBJet)fillCase(taggedBCase, isMixing, weights*extra);
	if(isMC) {
		if(foundTrueBJet)fillCase(trueBCase, isMixing, weights);
		if(foundTaggedBJet && foundTrueBJet)fillCase(taggedTrueBCase, isMixing, weights*extra);
		if(foundTaggedBJet && !foundTrueBJet)fillCase(contCase, isMixing, weights*extra);
	}
}

