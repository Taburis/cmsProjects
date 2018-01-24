
#define config_H
#include "TROOT.h"
#include "string"
#include <iostream>
#include "TMath.h"
#include "inputTree.h"
#include "dataTree.h"
#include "xthf4.h"
#include "xTree.h"
#include "TString.h"
// old cymbal correction
//#include "../dataSet/corrTableCymbal/xiaoTrkCorr.h"
// fine bin cymbal correction
#include "xiaoTrkCorr.h" //new correction
#include "xiaoTrkCorr_pp.h" //new correction
#include "TF1.h"

namespace jetTrack{
		namespace ioconfig{
				TString trackingOutput = "../hist/";
		}
		namespace weight {
				TF1* fvz   = new TF1("fvz","pol6",-15,15);
				TF1* fppvz   = new TF1("fppvz","gaus",-15,15);
				TF1* fcent1= new TF1("fcent1","[0]+[1]*x+[2]*x^2+[3]*x^3+[4]*x^4+[7]*exp([5]+[6]*x)",0,180);
				float cent_weight(int hiBin ){
						return (hiBin < 194) ? fcent1->Eval(hiBin) : 1;
				}
				float vz_weight(float vz){
						return fvz->Eval(vz);
				}
				float pp_vz_weight(float vz){
						return 1.0/fppvz->Eval(vz);
				}
		}
		namespace correction{
				// no dca cut correction
				xiaoTrkCorr* trkc= new xiaoTrkCorr("../../corrTable/patched_corrTable_cymbal_noDCAcuts.root");
				xiaoTrkCorr_pp* trkc_pp= new xiaoTrkCorr_pp("../../corrTable/corrTable_pp_noDCA.root");
				// broder bin correction
				//xiaoTrkCorr* trkc= new xiaoTrkCorr("../dataSet/corrTableCymbal/inputCorr_cymbalTune.root");
				float trk_corr(inputTree *t, int j ){
						return trkc->getTrkCorr(t->trkPt->at(j), t->trkEta->at(j), t->trkPhi->at(j),t->hiBin);	
				}
				float trk_corr_pp(inputTree *t, int j ){
						return trkc_pp->getTrkCorr(t->trkPt->at(j), t->trkEta->at(j), t->trkPhi->at(j),0,0);	
				}
				float trk_corr(dataTree *t, int j ){
						return trkc->getTrkCorr(t->trkPt->at(j), t->trkEta->at(j), t->trkPhi->at(j),t->hiBin);	
				}
				float trk_corr(xTree *t, int j ){
						return trkc->getTrkCorr(t->trkPt->at(j), t->trkEta->at(j), t->trkPhi->at(j),t->hiBin);	
				}
		}
		/*
		*/

		namespace trackingCorrConfig{
				//the bin config for initial scan 
				int ntrkpt_in = 22;
				float trkpt_in_c[23]= {0.7,0.8, 0.9, 1, 1.2, 1.4, 1.6, 1.8, 2,
						2.5, 3, 3.5, 4,5, 6, 7, 8, 10, 12, 16, 20, 50, 999};
				int ncent_in = 20; 
				float cent_in_c[21]; 
				int neta= 100;
				int nphi= 72;
				float etamin = -5, etamax =5;
				float phimin = -TMath::Pi(), phimax = TMath::Pi() ;
				//the bin config fo the correction table scan 
				int ncent_out = 17; 
				float cent_out_c [18]; 
				int ntrkpt_out = ntrkpt_in;
				float* trkpt_out_c= trkpt_in_c;
				//
				//----------------------------------------------------------------
				//
				float * cent_in = cent_in_c;
				float * cent_out = cent_out_c;
				float* trkpt_in= trkpt_in_c;
				float *trkpt_out = trkpt_out_c;
		}

		namespace trackingClosureConfig{
				//bin scheme of closure table we want to build, usually the same as corr output;
				int ntrkpt = trackingCorrConfig::ntrkpt_out; 
				float* trkpt_c = trackingCorrConfig::trkpt_in_c;
				int ncent = trackingCorrConfig::ncent_out; 
				float* cent= trackingCorrConfig::cent_out;
				int neta= 100;
				int nphi= 72;
				float* trkpt = trkpt_c;
				//bin scheme for the pt distribution
				int npt = 3993;
				float ptmin = 0.7;
				float ptmax = 400;
				float etamin = -5, etamax =5;
				float phimin = -TMath::Pi(), phimax = TMath::Pi() ;
				int ncent_pt = 4;
				Double_t cent_pt[5] = {0, 20, 60, 100, 200};
		}
		namespace correlationConfig{
				int ntrkpt = 9;
				float trkpt[10]= {0.7, 1, 2, 3, 4, 8,12, 16,20,999};
				int ncent = 4; 
				float cent[5] = {0, 20, 60, 100, 200}; 
				int njetpt=3;
				float jetpt[4]={120,150,200,999};
				int ndeta= 500;
				int ndphi= 200;
				Double_t detalow = -5, detaup =5;
				//Double_t dphilow = -TMath::Pi()/2, dphiup = 3*TMath::Pi()/2 ;
				Double_t dphilow = -TMath::Pi()/2, dphiup = 4.7 ;
		}

		namespace anaConfig{
				int ntrkpt = 9;
				float trkpt_c[10]= {0.7,1, 2, 3, 4, 8,12, 16, 20,999};
				int ncent = 4; 
				float cent[5] = {0, 20, 60, 100, 200}; 
				float *trkpt= trkpt_c;
		}

		void loadConfig(){
				weight::fvz->SetParameters(1.18472,-0.132675,0.00857998,-0.000326085,-1.48786e-06,4.68665e-07,-7.32942e-09 );	
				weight::fcent1->SetParameters(4.40810, -7.75301e-02, 4.91953e-04, -1.34961e-06, 1.44407e-09, -160, 1, 3.68078e-15);
				weight::fppvz->SetParameter(0,1.10477);
				weight::fppvz->SetParameter(1,2.52738);
				weight::fppvz->SetParameter(2,1.30296e1);
				/* initiate the cent bin for tracking */
				for(int i=0;i<trackingCorrConfig::ncent_in+1;++i) trackingCorrConfig::cent_in_c[i] = i*10;	
				for(int i=0;i<trackingCorrConfig::ncent_out;++i) trackingCorrConfig::cent_out_c[i] = i*10;
				trackingCorrConfig::cent_out_c[trackingCorrConfig::ncent_out] = 200;
		}


		bool trackQualityCutsImp (bool highPurity, float pfHcal, float pfEcal, float trkPt, float trkPtError, float trkEta, float trkDz, float trkDzError, float trkDxy, float trkDxyError, int trkNHit, float trkChi2, int trkNdof, int trkNlayer){
				if (!highPurity) return 1;
				if ( TMath::Abs(trkEta)>=2.4) return 1;
				if ( trkNHit< 11 ) return 1;
				if ( trkPt<= 0.7 || trkPt > 400) return 1;
				if ( trkPtError/trkPt>=0.1) return 1;
				//if ( trkPtError/trkPt>=0.3) return 1;
				//		if (TMath::Abs(trkEta)>1.6) return 1; // check mid-rapidity
				//if ( TMath::Abs(trkDz/trkDzError) >=3) return 1;
				//if ( TMath::Abs(trkDxy/trkDxyError) >=3) return 1;
				if (float(trkChi2)/float(trkNdof)/float(trkNlayer)>0.15) return 1;
				float Et = (pfHcal+pfEcal)/TMath::CosH(trkEta);
				if (!(trkPt<20 || (Et>0.5*trkPt))) return 1;
				/*
				*/
				return 0;
		}

		bool trackQualityCuts(inputTree * t , int j){
				return trackQualityCutsImp(t->highPurity->at(j), t->pfHcal->at(j), t->pfEcal->at(j),
								t->trkPt->at(j), t->trkPtError->at(j), t->trkEta->at(j), t->trkDz->at(j), t->trkDzError->at(j),
								t->trkDxy->at(j), t->trkDxyError->at(j), t->trkNHit->at(j), t->trkChi2->at(j),
								t->trkNdof->at(j), t->trkNlayer->at(j));
		}
		bool trackQualityCuts(dataTree * t , int j){
				return trackQualityCutsImp(t->highPurity->at(j), t->pfHcal->at(j), t->pfEcal->at(j),
								t->trkPt->at(j), t->trkPtError->at(j), t->trkEta->at(j), t->trkDz->at(j), t->trkDzError->at(j),
								t->trkDxy->at(j), t->trkDxyError->at(j), t->trkNHit->at(j), t->trkChi2->at(j),
								t->trkNdof->at(j), t->trkNlayer->at(j));
		}
		bool trackQualityCuts(xTree * t , int j){
				return trackQualityCutsImp(t->highPurity->at(j), t->pfHcal->at(j), t->pfEcal->at(j),
								t->trkPt->at(j), t->trkPtError->at(j), t->trkEta->at(j), t->trkDz->at(j), t->trkDzError->at(j),
								t->trkDxy->at(j), t->trkDxyError->at(j), t->trkNHit->at(j), t->trkChi2->at(j),
								t->trkNdof->at(j), t->trkNlayer->at(j));
		}
		bool eventCutsImp(int HBHENoiseFilterResultRun2Loose, int pcollisionEventSelection, int pprimaryVertexFilter, int phfCoincFilter3, float vz){
				if ( HBHENoiseFilterResultRun2Loose ==0) return 1;
				if ( pcollisionEventSelection ==0) return 1;
				if ( pprimaryVertexFilter ==0) return 1;
				if ( phfCoincFilter3 ==0) return 1;
				if ( vz>15 || vz<-15) return 1;
				return 0;
		}
		bool eventCutsImp_pp(int HBHENoiseFilterResultRun2Loose, int pprimaryVertexFilter, float vz){
				if ( HBHENoiseFilterResultRun2Loose ==0) return 1;
				if ( pprimaryVertexFilter ==0) return 1;
				if ( vz>15 || vz<-15) return 1;
				return 0;
		}
		bool eventCuts(inputTree *t){
if( t->ispp)
				return eventCutsImp_pp(t->HBHENoiseFilterResultRun2Loose, t->pprimaryVertexFilter, t->vz);
else
				return eventCutsImp(t->HBHENoiseFilterResultRun2Loose, t->pcollisionEventSelection, t->pprimaryVertexFilter, t->phfCoincFilter3, t->vz);
		}
		bool eventCuts(dataTree *t){
				return eventCutsImp(t->HBHENoiseFilterResultRun2Loose, t->pcollisionEventSelection, t->pprimaryVertexFilter, t->phfCoincFilter3, t->vz);
		}
		bool eventCuts(xTree *t){
				if(t->isMC && !(t->pthat<80)) return 1;
				if( fabs(t->vz)>15.) return 1;
				if(t->ispp){
						if(!(t->collisionEventSelection)) return 1;
						if(!(t->HBHEFilter)) return 1;
				}
				else {
						if(!(t->HBHEFilter)) return 1;
						if(!(t->collisionEventSelection)) return 1;
						if(!(t->phfCoincFilter)) return 1;
						if(!(t->pprimaryVertexFilter)) return 1;
				}
				return 0;
		}

		bool genParticleCuts(inputTree *t,int j){
				if (t->chg->at(j)==0) return 1;
				//		if (TMath::Abs(t->eta->at(j))>1.6) return 1; // check mid-rapidity
				return 0;
		}
		bool jetCutImp(float pt,float trackMax){
				if( pt< 120) return 1; 
				if( trackMax/pt<=0.01) return 1;
				return 0;
		}
		bool jetCut(xTree *t,int j){
				return jetCutImp(t->corrpt->at(j), t->trackMax->at(j));
		}
}

