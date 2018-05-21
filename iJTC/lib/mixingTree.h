
#ifndef mixingTree_H
#define mixingTree_H
#include "TChain.h"

class mixingTree{
		public : mixingTree(TString name){t=new TChain(name);};
				 void init(bool ismc, bool ishi);
				 bool eventCuts();

				 TChain *t;
				 bool isMC, isHI;

				 float pthat;
				 float vz ;
				 int hiBin;

				 vector<float>  *jtpt=0, *jteta=0, *jtphi=0;// need jet info. for tk correction
				 vector<float> *pt=0, *eta=0, *phi=0, *chg=0;
				 vector<float> *trkPt=0, *trkEta=0, *trkPhi=0, *trkChi2=0,  *trkPtError=0;
				 vector<float> *trkDz=0, *trkDzError=0, *trkDxy=0, *trkDxyError=0, *pfEcal=0, *pfHcal=0;
				 vector<int> *trkNlayer=0, *trkNdof=0, *trkNHit=0;
				 vector<bool> *highPurity=0; 
				 int  HBHEFilter=0, collisionEventSelection = 0;
};

void mixingTree::init(bool ismc, bool ishi){
		isMC = ismc; isHI = ishi;

		t->SetBranchStatus("*", 0);
		if(ishi){
				t->SetBranchStatus("hiBin",1);
				t->SetBranchAddress("hiBin",&hiBin);
		}
		else {
				t->SetBranchStatus("HBHENoiseFilterResultRun2Loose",1);
				t->SetBranchAddress("HBHENoiseFilterResultRun2Loose",&HBHEFilter);
		}
		t->SetBranchStatus("pprimaryVertexFilter",1);
		t->SetBranchAddress("pprimaryVertexFilter",&collisionEventSelection);
		t->SetBranchStatus("vz",1);
		t->SetBranchAddress("vz",&vz);

		if(ismc){
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
		}
		t->SetBranchStatus("pf_corrpt",1);
		t->SetBranchAddress("pf_corrpt",&jtpt);
		t->SetBranchStatus("pf_jteta",1);
		t->SetBranchAddress("pf_jteta" ,&jteta);
		t->SetBranchStatus("pf_jtphi",1);
		t->SetBranchAddress("pf_jtphi" ,&jtphi);
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
		if(!collisionEventSelection || !HBHEFilter){
				//cout<<"!"<<endl;
			   	return 0;
		}
	}
	return 1;
}

#endif
