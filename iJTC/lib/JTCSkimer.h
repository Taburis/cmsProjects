
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

struct histCase{
		//if want to add any hist, need to add it in the quickRegHist as well, and add the filling in the fillCase
		TH2D** sig;
		TH2D** sig_pTweighted;
		TH2D** mixing;
		TH1D** jet_corrpt;
		TH1D** jet_eta;
		TH1D** jet_phi;
};

const int nPt = 6;
const int nCent = 2;
TF1 *ppVz = new TF1("newppVz","gaus",-15,15);
void config(){
		ppVz->SetParameter(0,1.10477);
		ppVz->SetParameter(1,2.52738);
		ppVz->SetParameter(2,1.30296e1);
}

class JTCSkimer{
		public : 
				JTCSkimer(){ 
				};
				bool doGenJet=1, doGenTrack=1;
				void connectTree(TTree* ct, bool isMC, bool isHI);
				void initTree(TTree* ct, bool isMC, bool isHI);
				void addMixingTree(TString fname);
				void configMixing();
				bool isHI, isMC, domixing = 0;
				void eventWeight(float pthat, float vz, float hibin);
				void fillEventInfo();
				void addGenJet   (TTree* t);
				void addRecoJet  (TTree* t);
				void addGenTrack (TTree* t);
				void addRecoTrack(TTree* t);
				int trkCuts(float trkPt, float trkPtError, float trkEta, float trkChi2, int highPurity, int trkNHit, float trkNdof, float trkNlayer, float pfHcal, float pfEcal);
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
				void fillCase(histCase &hc, bool ismixing = 0);
				void fillCaseMix(histCase &hc);
				void fillHist(bool isMixing = 0);
				void fillAllJetInfo();
				void fillJetInfo(histCase &hc);
				float findDrmin(float eta, float phi, vector<float> *jetpt, vector<float> *jeteta, vector<float> *jetphi);
				void mixingSection(Long64_t voidIndex);
		public: 
				histManager *hm;
				//mixing config------------------------------
				bool do_mixing=0;
				int nPerTrigger=25;
				int nvz_mix, ncent_mix;
				//-----------------------------------------
				float ptbins[nPt+1] = {1, 2, 3, 4, 8, 12, 400};
				TString ptLabel[nPt+1] = {"Trk1", "Trk2", "Trk3", "Trk4", "Trk8", "trk12", "Trk400"};
				float centbins[nCent+1] = {0, 30, 100};
				TString centLabel[nCent+1] = {"Cent0", "Cent30", "Cent100"};
				TTree* t;
				// event branch
				float pthat;
				float vz ;
				int hiBin=0;
				mixingTree *mt = 0;
				// jet branch
				vector<float> *gen_jtpt=0, *gen_jteta=0, *gen_jtphi=0, *jtpt=0, *jteta=0, *jtphi=0;
				vector<float> *discr_csv=0;
				vector<int>	*flavorForB=0;
				// gen jet branch

				//gen particle branch
				vector<float> *pt=0, *eta=0, *phi=0, *chg=0;

				//track branch
				vector<float> *trkPt=0, *trkEta=0, *trkPhi=0, *trkChi2=0, *trkPtError=0;
				vector<float> *trkDz=0, *trkDzError=0, *trkDxy=0, *trkDxyError=0, *pfEcal=0, *pfHcal=0;
				vector<int> *trkNlayer=0, *trkNdof=0, *trkNHit=0;
				vector<bool> *highPurity=0; 

				int  HBHEFilter=0, collisionEventSelection = 0;

				//variables using in the correlations
				float selectedJt_pt;
				float selectedJt_eta;
				float selectedJt_phi;
				float selectedTrack_pt;
				float selectedTrack_eta;
				float selectedTrack_phi;
				float weights;
				bool foundTrueBJet, foundTaggedBJet, foundInclJet;
				float trkCorr=1, JESCorr=1;
				xAxis *ptax, *centax;
				// end of the variables declaration	
				std::vector<unsigned int >** mixTable;

				TH1D* hvz, *hpthat, *hcent;
				// if add a new histCase here, need to add it in the quickHistReg function as well
				histCase inclCase;
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
		t->SetBranchStatus("pf_corrpt",1);
		t->SetBranchAddress("pf_corrpt",&jtpt);
		t->SetBranchStatus("pf_jteta",1);
		t->SetBranchAddress("pf_jteta" ,&jteta);
		t->SetBranchStatus("pf_jtphi",1);
		t->SetBranchAddress("pf_jtphi" ,&jtphi);
		t->SetBranchStatus("pf_refparton_flavorForB" ,1);
		t->SetBranchAddress("pf_refparton_flavorForB" ,&flavorForB);
		t->SetBranchStatus("pf_discr_csvV1" ,1);
		t->SetBranchAddress("pf_discr_csvV1" ,&discr_csv);
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
				addGenJet(ct);
				addGenTrack(ct);
		}
		if(ishi){
				ct->SetBranchStatus("hiBin",1);
				ct->SetBranchAddress("hiBin",&hiBin);
		}
		else {
				ct->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
				ct->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHEFilter);
		}
		ct->SetBranchStatus("pprimaryVertexFilter",1);
		ct->SetBranchAddress("pprimaryVertexFilter",&collisionEventSelection);
		ct->SetBranchStatus("vz",1);
		ct->SetBranchAddress("vz",&vz);
		addRecoJet  (ct);
		addRecoTrack(ct);
}

void JTCSkimer::connectTree(TTree *ct, bool ismc, bool ishi){
		isMC = ismc; isHI = ishi;
		t=ct;
		initTree(t, isMC, isHI);
}

void JTCSkimer::quickHistReg(histManager *h, TString cap, histCase &hc){
		int nHistoBinsX = 500;
		int nHistoBinsY = 200;
		hc.sig= new TH2D*[nPt*nCent];
		hc.sig_pTweighted= new TH2D*[nPt*nCent];
		hc.mixing= new TH2D*[nPt*nCent];
		TString tmp, name;
		if(doGenJet ) name= cap+"_GenJet_";
		else name=cap+"RecoJet_";
		if(doGenTrack) name= name+"GenTrack";
		else name = name +"RecoTrack";
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						tmp = centLabel[j]+"_"+centLabel[j+1]+"_"+ptLabel[i]+"_"+ptLabel[i+1];
						hc.sig[i+j*nPt] = hm->regHist<TH2D>(name+Form("_%d_%d",i, j), tmp,
										nHistoBinsX,-5,5,nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
						hc.sig_pTweighted[i+j*nPt] = hm->regHist<TH2D>(name+Form("_pTweighted_%d_%d",i, j), tmp,
										nHistoBinsX,-5,5, nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
						hc.mixing[i+j*nPt] = hm->regHist<TH2D>(name+Form("_mixing_%d_%d",i, j), tmp,
										nHistoBinsX,-5,5, nHistoBinsY,-TMath::Pi()/2,3*TMath::Pi()/2);
				}
		}

		if(doGenJet ) name=cap+ "_GenJet";
		else name=cap+ "_RecoJet";
		hc.jet_corrpt = new TH1D*[nCent];
		hc.jet_eta = new TH1D*[nCent];
		hc.jet_phi = new TH1D*[nCent];
		for(int j=0; j<nCent; ++j){
				tmp = centLabel[j]+" to "+centLabel[j+1];
				hc.jet_corrpt[j] = hm->regHist<TH1D>(name+Form("_corrpt_%d",j), tmp, 50, 0, 500.0);
				hc.jet_eta[j] = hm->regHist<TH1D>(name+Form("_eta_%d",j), tmp, 100, -2.0, 2.0);
				hc.jet_phi[j] = hm->regHist<TH1D>(name+Form("_phi_%d",j), tmp, 72, -TMath::Pi(), TMath::Pi());
		}
}

void JTCSkimer::configHist(){
		hm = new histManager();

		ptax = new xAxis(nPt, ptbins);
		centax= new xAxis(nCent, centbins);
		hvz = hm->regHist<TH1D>("vzInfo", "", 40, -20, 20);
		hcent = hm->regHist<TH1D>("centInfo","",  50, 0, 200);
		hpthat = hm->regHist<TH1D>("pthatInfo", "", 100, 0, 400);
		quickHistReg(hm, "inclJet", inclCase);
		quickHistReg(hm, "taggedBJet", taggedBCase);
		quickHistReg(hm, "trueBJet", trueBCase);
		quickHistReg(hm, "taggedTrueBJet", taggedTrueBCase);
		hm->sumw2();
}

void JTCSkimer::write(TString name){
		TFile *wf = new TFile(name, "recreate");
		wf->cd();
		hm->write();
		wf->Close();
}

void JTCSkimer::fillAllJetInfo(){
		fillJetInfo(inclCase);
		if(foundTrueBJet)fillJetInfo(trueBCase);
		if(foundTaggedBJet)fillJetInfo(taggedBCase);
}

void JTCSkimer::fillJetInfo(histCase &hc){
		for(int j=0; j<nCent; ++j){
				hc.jet_corrpt[j]->Fill(selectedJt_pt, weights);
				hc.jet_eta[j]->Fill(selectedJt_eta, weights);
				hc.jet_phi[j]->Fill(selectedJt_phi, weights);
		}
}

void JTCSkimer::fillCase(histCase &hc, bool ismixing){
		float deta = selectedJt_eta - selectedTrack_eta;
		float dphi = selectedJt_phi - selectedTrack_phi;
		while(dphi>(1.5*TMath::Pi())){dphi+= -2*TMath::Pi();}
		while(dphi<(-0.5*TMath::Pi())){dphi+= 2*TMath::Pi();}

		int i = ptax->findBin(selectedTrack_pt);
		if(!isHI){ 
				for(int j=0; j<nCent; ++j){
						if(ismixing){
								//								cout<<"filling mixing"<<endl;
								hc.mixing[i+j*nPt]->Fill(deta, dphi, weights*trkCorr);
						} else {
								hc.sig[i+j*nPt]->Fill(deta, dphi, weights*trkCorr);
								hc.sig_pTweighted[i+j*nPt]->Fill(deta, dphi, weights*selectedTrack_pt*trkCorr);
						}
				}
		}
}

void JTCSkimer::fillEventInfo(){
		hvz->Fill(vz, weights);
		hpthat->Fill(pthat, weights);
}

void JTCSkimer::addMixingTree(TString name){
		if(mt==0) mt = new mixingTree(t->GetName());
		cout<<t->GetName()<<" has been created"<<endl;
		mt->t->Add(name);
}

void JTCSkimer::fillHist(bool isMixing){
		fillCase(inclCase, isMixing);
		if(foundTrueBJet)fillCase(trueBCase, isMixing);
		if(foundTaggedBJet)fillCase(taggedBCase, isMixing);
		if(foundTaggedBJet&&foundTrueBJet)fillCase(taggedTrueBCase, isMixing);
}

