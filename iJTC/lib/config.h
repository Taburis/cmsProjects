
#ifndef config_H
#define config_H
#include "TF1.h"

float trkPtcut = 1;
float trketamaxcut = 2.4;

TF1 *ppVz = new TF1("newppVz","gaus",-15,15);
double xsecs[10] = {5.335E-01, 3.378E-02, 3.778E-03, 4.412E-04, 6.147E-05, 1.018E-05, 2.477E-06, 6.160E-07, 1.088E-07, 0};
double pthatbins[10] = {15,30,50,80,120,170,220,280,370,9999};
double ppPthatEntries[9] = {0,0,272902,377559,467823,447683,259111,234347,50942};

void config(){
	ppVz->SetParameter(0,1.10477);
	ppVz->SetParameter(1,2.52738);
	ppVz->SetParameter(2,1.30296e1);
}

bool eventCut(float vz,float pthat, int HBHEFilter, int collisionEvtSel,  bool isHI, bool isMC){
		// for pp the collisionEvtSel should be filled by pprimaryVertexFilter
		if(TMath::Abs(vz)>15) return 1;
		if(isMC && pthat<80) return 1;
		if(!isHI){
				if(!collisionEvtSel||!HBHEFilter) return 1;
		}
		return 0;
}

Double_t wpthat(float pthat, bool isHI){
	Double_t weights=1;
	if(!isHI){
		int ibin=0;
		while(pthat>pthatbins[ibin+1]) ibin++;
		weights=(xsecs[ibin]-xsecs[ibin+1])/ppPthatEntries[ibin];
	}
	return weights;
}

float wvz(float vz, bool isHI){
		if(!isHI)
				return 1.0/ppVz->Eval(vz);
		else return 1;
}

int trkCuts(bool doDCA, float trkpt, float trkpterror, float trketa, float trkchi2, int highpurity, int trknhit, float trkndof, float trknlayer, float pfhcal, float pfecal, float trkdz,float trkdxy, float trkdzerror, float trkdxyerror){
		//return 0 if it passed cuts;
	if(trkpt<=trkPtcut || trkpt > 400) return 2;
	if(TMath::Abs(trketa) >=trketamaxcut) return 3;
	if(!highpurity) return 4;
	if(fabs(trkpterror/trkpt)>0.3) return 5;
	if(trknhit<11 ) return 6;
	if(trkchi2/trkndof/trknlayer > 0.15) return 7;

	float Et = (pfhcal+pfecal)/TMath::CosH(trketa);
	if(!(trkpt<20 || Et > 0.5*trkpt)) return 8;
	if(doDCA) {
			if(TMath::Abs(trkdz/trkdzerror)>=3.0 ||
					TMath::Abs(trkdxy/trkdxyerror)>=3.0) return 1;
	}
	return 0;
}

#endif

