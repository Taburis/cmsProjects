
#ifndef config_H
#include "config.h"
#endif

#ifndef inputTree_H
#include "inputTree.h"
#endif
#ifndef dataTree_H
#include "dataTree.h"
#endif

#ifndef xthf4_H
#include "xthf4.h"
#endif

#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif

#ifndef ROOT_TLegend
#include "TLegend.h"
#endif

#ifndef ROOT_TLine
#include "TLine.h"
#endif

#ifndef xUtility
#include "utility.h"
#endif

#ifndef trackingClosure_H
#define trackingClosure_H
using namespace jetTrack;

class trackingClosure  {
		public: 
				trackingClosure(TString file):
						ntrkpt(trackingClosureConfig::ntrkpt),
						ncent(trackingClosureConfig::ncent),
						nx(trackingClosureConfig::neta),
						ny(trackingClosureConfig::nphi),
						etamin(trackingClosureConfig::etamin),
						etamax(trackingClosureConfig::etamax),
						phimin(trackingClosureConfig::phimin),
						phimax(trackingClosureConfig::phimax),
						trkpt(trackingClosureConfig::trkpt),
						cent(trackingClosureConfig::cent)
		{
				i_t=NULL;
				d_t=NULL;
				isread = 1;
				gen = new xthf4();
				rec = new xthf4();
				cre = new xthf4();
				f= TFile::Open(file);
		};
				trackingClosure(inputTree* tt):
						i_t(tt),
						ntrkpt(trackingClosureConfig::ntrkpt),
						ncent(trackingClosureConfig::ncent),
						nx(trackingClosureConfig::neta),
						ny(trackingClosureConfig::nphi),
						etamin(trackingClosureConfig::etamin),
						etamax(trackingClosureConfig::etamax),
						phimin(trackingClosureConfig::phimin),
						phimax(trackingClosureConfig::phimax),
						trkpt(trackingClosureConfig::trkpt),
						cent(trackingClosureConfig::cent)
		{
				int   npt  = trackingClosureConfig::npt;
				float ptmin= trackingClosureConfig::ptmin;
				float ptmax= trackingClosureConfig::ptmax;
				int   ncent_pt   =trackingClosureConfig::ncent_pt;
				Double_t* cent_pt=trackingClosureConfig::cent_pt;
				d_t=NULL;
				f=NULL;
				isread = 0;
				if(i_t->ispp){
						fvz_weight = &(weight::pp_vz_weight);
						ncent = 1;
						float ppcent[] = {0,1};
						Double_t ppcent_pt[] = {0,1};
						gen = new xthf4("gen","gen", nx,etamin, etamax , ny, phimin, phimax,ntrkpt, trkpt,ncent, ppcent);
						rec = new xthf4("rec","rec", nx,etamin, etamax , ny, phimin, phimax,ntrkpt, trkpt,ncent, ppcent);
						cre = new xthf4("cre","cre", nx,etamin, etamax , ny, phimin, phimax,ntrkpt, trkpt,ncent, ppcent);
						genpt= new TH2F("genpt", "gen pt distriubtion"        ,npt, ptmin, ptmax, ncent, ppcent_pt);
						recpt= new TH2F("recpt", "track pt distriubtion"      ,npt, ptmin, ptmax, ncent, ppcent_pt);
						crept= new TH2F("crept", "corr. track pt distriubtion",npt, ptmin, ptmax, ncent, ppcent_pt);
				}
				else{
						fvz_weight = &(weight::vz_weight);
						gen = new xthf4("gen","gen", nx,etamin, etamax , ny, phimin, phimax,ntrkpt, trkpt,ncent, cent);
						rec = new xthf4("rec","rec", nx,etamin, etamax , ny, phimin, phimax,ntrkpt, trkpt,ncent, cent);
						cre = new xthf4("cre","cre", nx,etamin, etamax , ny, phimin, phimax,ntrkpt, trkpt,ncent, cent);
						genpt= new TH2F("genpt", "gen pt distriubtion"        ,npt, ptmin, ptmax, ncent_pt, cent_pt);
						recpt= new TH2F("recpt", "track pt distriubtion"      ,npt, ptmin, ptmax, ncent_pt, cent_pt);
						crept= new TH2F("crept", "corr. track pt distriubtion",npt, ptmin, ptmax, ncent_pt, cent_pt);
				}
				genpt->Sumw2();
				recpt->Sumw2();
				crept->Sumw2();
		};

				trackingClosure(dataTree* tt):
						d_t(tt),
						ntrkpt(trackingClosureConfig::ntrkpt),
						ncent(trackingClosureConfig::ncent),
						nx(trackingClosureConfig::neta),
						ny(trackingClosureConfig::nphi),
						etamin(trackingClosureConfig::etamin),
						etamax(trackingClosureConfig::etamax),
						phimin(trackingClosureConfig::phimin),
						phimax(trackingClosureConfig::phimax),
						trkpt(trackingClosureConfig::trkpt),
						cent(trackingClosureConfig::cent)
		{
				int   npt  = trackingClosureConfig::npt;
				float ptmin= trackingClosureConfig::ptmin;
				float ptmax= trackingClosureConfig::ptmax;
				int   ncent_pt   =trackingClosureConfig::ncent_pt;
				Double_t* cent_pt=trackingClosureConfig::cent_pt;
				f=NULL;
				i_t=NULL;
				isread = 0;
				rec = new xthf4("rec","rec", nx,etamin, etamax , ny, phimin, phimax,ntrkpt, trkpt,ncent, cent);
				cre = new xthf4("cre","cre", nx,etamin, etamax , ny, phimin, phimax,ntrkpt, trkpt,ncent, cent);
				recpt= new TH2F("recpt", "track pt distriubtion"      ,npt, ptmin, ptmax, ncent_pt, cent_pt);
				crept= new TH2F("crept", "corr. track pt distriubtion",npt, ptmin, ptmax, ncent_pt, cent_pt);
				recpt->Sumw2();
				crept->Sumw2();
		};

				void runScan();
				void runMCScan(inputTree* t);
				void runDataScan(dataTree *t);
				void Read();
				void DrawClosure(TString name="", TString type= "");
				TPad* ratioPlot(TH1F* hden, TH1F* hnum1, TH1F* hnum2, 
								float xmin, float xmax,float ymin,float ymax, bool doLegend, bool dofit=0);
				~trackingClosure();

		public: 
				float (*fvz_weight)(float );
				TFile *f;
				inputTree * i_t;
				dataTree * d_t;
				int ntrkpt;
				int ncent;
				int nx, ny;
				float *trkpt, *cent;
				float etamin, etamax, phimin, phimax;
				xthf4 * gen;
				xthf4 * rec;
				xthf4 * cre;
				TH2F * genpt;
				TH2F * recpt;
				TH2F * crept;
				bool isread;
				std::vector<TH1F*> ratioVec;
				bool ispp;
};

void trackingClosure::runScan(){
		if(d_t ==NULL) runMCScan(i_t);
		else runDataScan(d_t);
}

void trackingClosure::runMCScan(inputTree *t){
		std::cout<<"prcessing the MC tree.."<<std::endl;
		loadConfig();
		Long64_t nentries = t->fChain->GetEntriesFast();
		cout<<"here"<<endl;
		for(Long64_t jentry = 0; jentry <nentries ; ++jentry){
				if(jentry %1000 == 0) std::cout<<"processing "<<jentry<<"th event.."<<std::endl;
				t->GetEntry(jentry);
				float w_vz=1;
				float w_cent=1;
				if(eventCuts(t)) continue;
				w_vz = fvz_weight(t->vz);
				if(i_t->ispp) w_cent=1;
				else	w_cent = weight::cent_weight(t->hiBin);
				for(int j=0;j<t->trkPt->size(); ++j){
						if(trackQualityCuts(t,j))continue;
						float trkcorr = 1;
						if( i_t->ispp){
								trkcorr= correction::trk_corr_pp(t, j);
						}
						else trkcorr= correction::trk_corr(t, j);
						if(i_t->ispp) {
								rec->Fill(t->trkEta->at(j), t->trkPhi->at(j), t->trkPt->at(j), 0, w_vz);
								cre->Fill(t->trkEta->at(j), t->trkPhi->at(j), t->trkPt->at(j), 0, w_vz*trkcorr);
								recpt->Fill(t->trkPt->at(j), 0.0, w_cent*w_vz);
								crept->Fill(t->trkPt->at(j), 0.0, w_cent*w_vz*trkcorr);
						}
						else {
								rec->Fill(t->trkEta->at(j), t->trkPhi->at(j), t->trkPt->at(j), t->hiBin, w_vz);
								cre->Fill(t->trkEta->at(j), t->trkPhi->at(j), t->trkPt->at(j), t->hiBin, w_vz*trkcorr);
								recpt->Fill(t->trkPt->at(j), t->hiBin, w_cent*w_vz);
								crept->Fill(t->trkPt->at(j), t->hiBin, w_cent*w_vz*trkcorr);
						}
				}
				for(int j=0;j<t->pt->size(); ++j){
						if(genParticleCuts(t,j))continue;
						if(i_t->ispp) {
								gen->Fill(t->eta->at(j), t->phi->at(j), t->pt->at(j), 0.0, w_vz);
								genpt->Fill(t->pt->at(j), 0.0, w_cent*w_vz);
						}
						else{
								gen->Fill(t->eta->at(j), t->phi->at(j), t->pt->at(j), t->hiBin, w_vz);
								genpt->Fill(t->pt->at(j), t->hiBin, w_cent*w_vz);
						}
				}
		}
		genpt->Scale(1.0/nentries);
		recpt->Scale(1.0/nentries);
		crept->Scale(1.0/nentries);
		gen->Write();
		rec->Write();
		cre->Write();
		genpt->Write();
		recpt->Write();
		crept->Write();
}

void trackingClosure::runDataScan(dataTree* t){
		std::cout<<"prcessing the data tree.."<<std::endl;
		loadConfig();
		Long64_t nentries = t->fChain->GetEntriesFast();
		for(Long64_t jentry = 0; jentry <nentries ; ++jentry){
				if(jentry %1000 == 0) std::cout<<"processing "<<jentry<<"th event.."<<std::endl;
				t->GetEntry(jentry);
				float w_vz=1;
				float w_cent=1;
				if(eventCuts(t)) continue;
				for(int j=0;j<t->trkPt->size(); ++j){
						if(trackQualityCuts(t,j))continue;
						float trkcorr = correction::trk_corr(t, j);
						rec->Fill(t->trkEta->at(j), t->trkPhi->at(j), t->trkPt->at(j), t->hiBin);
						cre->Fill(t->trkEta->at(j), t->trkPhi->at(j), t->trkPt->at(j), t->hiBin, trkcorr);
						recpt->Fill(t->trkPt->at(j), t->hiBin);
						crept->Fill(t->trkPt->at(j), t->hiBin, trkcorr);
				}
		}
		recpt->Scale(1.0/nentries);
		crept->Scale(1.0/nentries);
		rec->Write();
		cre->Write();
		recpt->Write();
		crept->Write();
}

void trackingClosure::Read(){
		if(f==NULL) {
				std::cout<<"trackingClosure contained histograms, can't read anymore"<<std::endl;
				return ;
		}
		int ana_ntrkpt=anaConfig::ntrkpt;
		int ana_ncent=anaConfig::ncent;
		float* ana_trkpt = anaConfig::trkpt;
		float* ana_cent = anaConfig::cent;
		if(ispp){
				float ppcent[2] = {0,1 };
				ncent =1;
				cent = ppcent;
		}
		TFile* wf = TFile::Open("closure_Output.root", "recreate");
		gen->Read("gen", "gen", f, ntrkpt, trkpt, ncent, cent);
		rec->Read("rec", "rec", f, ntrkpt, trkpt, ncent, cent);
		cre->Read("cre", "cre", f, ntrkpt, trkpt, ncent, cent);
		gen->RebinZ(ana_ntrkpt, ana_trkpt);
		rec->RebinZ(ana_ntrkpt, ana_trkpt);
		cre->RebinZ(ana_ntrkpt, ana_trkpt);
		if( !ispp){
				gen->RebinW(ana_ncent, ana_cent);
				rec->RebinW(ana_ncent, ana_cent);
				cre->RebinW(ana_ncent, ana_cent);
		}
		genpt = (TH2F*) f->Get("genpt");
		recpt = (TH2F*) f->Get("recpt");
		crept = (TH2F*) f->Get("crept");
		gen->Write();
		rec->Write();
		cre->Write();
}

trackingClosure::~trackingClosure(){
		delete gen;
		delete rec; 
		delete cre; 
		for(auto &it : ratioVec) delete it;
		ratioVec.clear();
}

void trackingClosure::DrawClosure(TString name, TString type){
		float xmin =-2.7 , xmax=2.7 , ymin = 0.0, ymax= 2;
		int ana_ntrkpt=anaConfig::ntrkpt;
		int ana_ncent=anaConfig::ncent;
		float* ana_trkpt = anaConfig::trkpt;
		float* ana_cent = anaConfig::cent;
		if(ispp) {
				ana_ncent=1;
		}
		ratioVec.reserve(ana_ntrkpt*ana_ncent);
		TCanvas * c[ana_ntrkpt][2];
		TH1F *hgen, *hrec, *hcre;
		TLine *ll = new TLine();
		ll->SetLineStyle(8);


		auto ut= new padPolisher();
		ut->tx->SetTextSize(0.07);
		TPad* pd;
		float fitRange0 = xmin, fitRange1 = xmax;
		TF1* p0 = new TF1("p0", "pol0", fitRange0, fitRange1);
		for(int i=0;i<ana_ntrkpt ; ++i){
				c[i][0]= new TCanvas(Form("c_%d_0",i),"", 1200,600);
				c[i][0]->Divide(4,1,0,0);
				for(int j=0;j<ana_ncent; ++j){
						c[i][0]->cd(ana_ncent-j);
						hgen= (TH1F*)gen->hist(i,j)->ProjectionX();	
						hrec= (TH1F*)rec->hist(i,j)->ProjectionX();	
						hcre= (TH1F*)cre->hist(i,j)->ProjectionX();	
						hgen->Rebin();
						hrec->Rebin(); hcre->Rebin();
						hcre->GetXaxis()->SetTitle("#eta_{track}");
						pd=ratioPlot(hgen, hcre, hrec,xmin, xmax, ymin, ymax, !(ana_ncent-1-j), 1);
						ll->DrawLine(xmin, 1, xmax,1);
						ll->DrawLine(xmin, 1.05, xmax,1.05);
						ll->DrawLine(xmin, .95, xmax,.95);
						if(j==2){
								pd->cd(2);
								ut->ptLabel(i,ana_trkpt);
						}
						pd->cd(1);
						ut->centLabel(j);
				}		
				if(name !="") c[i][0]->SaveAs("TrackingEfficiency_Eta_"+name+Form("_%d",i)+"."+type);
		}
		xmin =-3.2 , xmax=3.2 , ymin = 0.0, ymax= 2;
		for(int i=0;i<ana_ntrkpt ; ++i){
				c[i][1]= new TCanvas(Form("c_%d_1",i),"", 1200,600);
				c[i][1]->Divide(4,1,0,0);
				for(int j=0;j<ana_ncent; ++j){
						c[i][1]->cd(ana_ncent-j);
						hgen= (TH1F*)gen->hist(i,j)->ProjectionY();	
						hrec= (TH1F*)rec->hist(i,j)->ProjectionY();	
						hcre= (TH1F*)cre->hist(i,j)->ProjectionY();	
						hgen->Rebin();
						hrec->Rebin();
						hcre->Rebin();
						hcre->GetXaxis()->SetTitle("#phi_{track}");
						pd=ratioPlot(hgen, hcre, hrec,xmin, xmax, ymin, ymax, !(ana_ncent-1-j));
						ll->DrawLine(xmin, 1, xmax,1);
						ll->DrawLine(xmin, 1.05, xmax,1.05);
						ll->DrawLine(xmin, .95, xmax,.95);
						if(j==2){
								pd->cd(2);
								ut->ptLabel(i,ana_trkpt);
						}
						pd->cd(1);
						ut->centLabel(j);
				}		
				//ut->cent4Labels(c[i][1]);
				if(name !="") c[i][1]->SaveAs("TrackingEfficiency_Phi_"+name+Form("_%d",i)+"."+type);
		}

		TH1F* closure[4];
		if( name =="" || ana_ncent>4) return;
		TCanvas *cpt = new TCanvas("cpt", "", 1200,340);
		cpt->Divide(ana_ncent,1,0,0);
		for(int ic=0; ic<ana_ncent; ++ic){
				cpt->cd(ana_ncent-ic);
				std::stringstream stream;
				stream<<"cent : "<<fixed<<setprecision(0)<<ana_cent[ic]/2<<"-"<<fixed<<setprecision(0)<<ana_cent[ic+1]/2<<"%";
				TString title = stream.str();
				closure[ic]= new TH1F(Form("closure_pt_%d",ic),title , ana_ntrkpt, ana_trkpt);
				for(int i=0;i<ana_ntrkpt;++i){
						cout<<i+ic*ana_ncent<<endl;
						ratioVec.at(i+ic*ana_ncent)->Fit(p0, "", "", fitRange0, fitRange1);
						closure[ic]->SetBinContent(i+1, p0->GetParameter(0));
						closure[ic]->SetBinError(i+1, p0->GetParError(0));
				}
				closure[ic]->Draw();
				gPad->SetLogx();
		}
		cpt->SaveAs(name+"_closure_pt"+"."+type);
		return;
}

TPad* trackingClosure::ratioPlot(TH1F* hden, TH1F* hnum1, TH1F* hnum2, 
				float xmin, float xmax,float ymin,float ymax, bool doLegend, bool dofit)
{
		TH1F* hratio1, *hratio2;
		hden->SetAxisRange(0, 1.8*hden->GetMaximum(), "y");
		hden->SetMarkerSize(1.4);
		hden->SetMarkerStyle(21);
		hden->SetMarkerColor(kBlue-5);
		hden->SetLabelSize(0.05, "Y");
		hden->SetLabelSize(0.05, "X");
		hnum1->SetMarkerSize(0.9);
		hnum1->SetMarkerStyle(20);
		hnum1->SetMarkerColor(kRed+1);
		hnum1->SetTitle("");
		hnum1->SetLabelSize(0.05, "Y");
		hnum1->SetLabelSize(0.05, "X");
		hnum1->GetXaxis()->SetTitleOffset(0.8);
		hnum1->GetXaxis()->CenterTitle();
		hnum1->GetXaxis()->SetTitleSize(0.06);
		hratio1=(TH1F*)hnum1->Clone();
		hden->SetAxisRange(xmin, xmax, "X");
		hratio1->SetAxisRange(xmin, xmax, "X");
		hratio1->Divide(hnum1, hden, 1, 1, "B");
		hnum2->SetMarkerSize(0.9);
		hnum2->SetMarkerStyle(20);
		hnum2->SetTitle("");
		hnum2->SetMarkerColor(kBlue-3);
		hratio2=(TH1F*)hnum2->Clone();
		hratio2->Divide(hnum2, hden, 1 ,1, "B");
		TPad* pad = (TPad*) gPad;
		pad->Divide(1, 2, 0,0);
		pad->cd(1);
		hden->SetStats(0);
		hden->SetTitle("");
		hden->Draw();
		gPad->SetTicks(2,2);
		hnum1->Draw("same");
		hnum2->Draw("same");
		TLegend *ltemp;
		if(doLegend){
				ltemp= new TLegend(0.15,0.67,0.6,0.97);
				ltemp->SetLineColor(kWhite);
				ltemp->SetFillColorAlpha(kWhite, 0);
				ltemp->AddEntry(hden, "gen. particles");
				ltemp->AddEntry(hnum1, "corr. tracks");
				ltemp->AddEntry(hnum2, "reco. tracks");
				ltemp->Draw();
		}
		pad->cd(2);
		hratio1->SetStats(0);
		hratio1->SetAxisRange(ymin, ymax , "Y");
		hratio1->Draw();
		hratio2->Draw("same");
		gPad->SetTicks(2,2);
		if(doLegend){
				TLegend *ltemp2 = new TLegend(0.15,0.72,0.8,0.97);
				ltemp2->SetLineColor(kWhite);
				ltemp2->SetFillColorAlpha(kWhite, 0);
				ltemp2->AddEntry(hratio1, "reco.Tracks/gen.Particles");
				ltemp2->AddEntry(hratio2, "corr.Tracks/gen.Particles");
				ltemp2->Draw();
		}
		ratioVec.push_back(hratio1);
		return (TPad*)pad;
}

#endif
