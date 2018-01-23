
#ifndef auxi_canvas_H
#define auxi_canvas_H
class auxi_canvas: public TCanvas {
		public :
				auxi_canvas(TString cname, TString title, float w, float h): 
						name(cname), cw(w), ch(h),TCanvas(cname, title, w, h){
						};
				void divide(int, int);
				void CD(int);
		public :
				TString name;
				TPad **pad;
				float cw, ch;
};

void auxi_canvas::divide(int nx, int ny){
		float ml = this->GetLeftMargin();
		float mr = this->GetRightMargin();
		float mt = this->GetTopMargin();
		float mb = this->GetBottomMargin();
		this->SetMargin(0,0,0,0);
		float h = (1-mt-mb)/nx;
		float w = (1-ml-mr)/ny;
		float absl = cw*ml;
		float absb = ch*mb;
		pad = new TPad*[nx*ny];
		for(int i=0; i<nx; ++i){
				for(int j=0; j<ny; ++j){
						float sl=0, sb =0;
						if(j==0) sl = ml;
						if(i==nx-1) sb = mb;
						pad[i*ny+j] = new TPad(name, "", w*j+ml-sl, h*(nx-i-1)+mb-sb, ml+w*(j+1), h*(nx-i)+mb);
						if(j==0) sl = absl/(absl+w*cw);
						if(i==nx-1) sb = absb/(absb+h*ch);
						pad[i*ny+j]->SetMargin(sl,0,sb,0);
						pad[i*ny+j]->Draw();
				}
		}
}

void auxi_canvas::CD(int i){
		pad[i-1]->cd();
}
#endif

#include "signalFactoryBase.h"
namespace plottingInput{

		TString trk_tag[] = {"TrkPt07", "TrkPt1", "TrkPt2","TrkPt3","TrkPt4","TrkPt8","TrkPt12","TrkPt16","TrkPt20","TrkPt300"};
		TString cent_tag[]= {"Cent0","Cent10","Cent30","Cent50","Cent100"};
		TString path = "/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/";
		// nominal plotting input files 
		//loading syst error
		TFile *js_syst_f = TFile::Open(path+"js_syst_err.root");
		TFile *yield_proj_syst_f = TFile::Open(path+"yield_proj_syst_err.root");
		TFile *run1_ref_f= TFile::Open(path + "Inclusive_Data_AllPlots.root");
		//TFile *js_dr_pb_f = TFile::Open(path+"Jet_Shapes_pTweighted_fixErr.root");
		TFile *js_dr_pp_f = TFile::Open(path+"Jet_Shapes_pTweighted_final.root");
		TFile *js_dr_pb_f = TFile::Open(path+"Jet_Shapes_pTweighted_final.root");
//		TFile *js_dr_pp_f = TFile::Open(path+"Jet_Shapes_pTweighted_kurt_pp.root");
//		TFile *js_dr_pp_f = TFile::Open(path+"Jet_Shapes_pTweighted_test.root");
		TFile *yield_dr_pb_f = TFile::Open(path+"Jet_Shapes_fixErr.root");
		TFile *yield_dr_pp_f = TFile::Open(path+"Jet_Shapes_kurt_pp.root");
//		TFile *yield_proj_pb_f = TFile::Open(path+"Particle_Yields_combined.root");
//		TFile *yield_proj_pp_f = TFile::Open(path+"Particle_Yields_for_pp_and_stat_err.root");
		TFile *yield_proj_pb_f = TFile::Open(path+"Particle_Yields_final.root");
		TFile *yield_proj_pp_f = TFile::Open(path+"Particle_Yields_final.root");

		void changeErr(TH1D* h , TH1D* err){
				for(int k=1; k<h->GetNbinsX()+1;++k){
						//						cout<<err->GetBinContent(k)<<endl;
						h->SetBinError(k, err->GetBinContent(k));
				}
		}

		TH1D* etaBin (TH1D* h ){
				Double_t etabin[] = {-2.5,-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.,2.5};
				//Double_t etabin[] = {-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.};
				//Double_t etabin[] = {-2.5,-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.,2.5};
				//Double_t etabin[] = {-3, -2.5,-2.,-1.5,  -1.,  -0.8,  -0.6,  -0.4, -0.3, -0.2, -0.1, 0.1, 0.2, 0.3, 0.4, 0.6,  0.8, 1., 1.5,2.,2.5, 3};
				TString tmp = h->GetName();
				tmp = tmp+"_rebin";
				TH1D* hn = new TH1D(tmp, "", sizeof(etabin) / sizeof(Double_t)-1, etabin);
				for(int k=1; k<hn->GetNbinsX()+1; ++k){
						int bin = h->FindBin(hn->GetBinCenter(k));
						hn->SetBinContent(k, h->GetBinContent(bin));
						hn->SetBinError(k, h->GetBinError(bin));
				}
				return (TH1D*) hn;
		}

		TH1D* drBin (TH1D* h ){
				Double_t drbin[] = {0, 0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.45,0.5,0.6,0.7,0.8,1.};
				TString tmp = h->GetName();
				tmp = tmp+"_rebin";
				TH1D* hn = new TH1D(tmp, "", 14, drbin);
				for(int k=1; k<hn->GetNbinsX()+1; ++k){
						hn->SetBinContent(k, h->GetBinContent(k+1));
						hn->SetBinError(k, h->GetBinError(k+1));
				}
				return (TH1D*) hn;
		}

		void symmetrization(TH1D* h){
				int stop; 
				if( (h->FindBin(0)-1)%2== 0) stop = h->FindBin(0);
				else stop = h->FindBin(0)+1;
				//cout<<stop<<endl;
				//cout<<h->GetNbinsX()<<endl;
				for(int k=1; k<stop+1 ; ++k){
						//		cout<<h->GetBinCenter(k)<<endl;
						float c = h->GetBinContent(k);
						float err = pow(h->GetBinError(k),2);
						float x = h->GetBinCenter(k);
						c=c+ h->GetBinContent(h->FindBin(-x));
						err = pow((err+ pow(h->GetBinError(h->FindBin(-x)),2))/2,0.5);
						h->SetBinContent(k, c/2);
						h->SetBinError(k, err);
						h->SetBinContent(h->FindBin(-x), c/2);
						h->SetBinError(h->FindBin(-x), err);
				}
				return ;
		}

		TH1D* py_dphi_diff_err[5];
		TH1D* py_dphi_diff[5];
		TH1D* py_dphi_err_all[8];
		TH1D* py_dphi_all[5]; 
		TH1D* py_dphi    [8][5];
		TH1D* py_dphi_err[8][5];

		void prePlot_py_dphi(){
				TString tmp;
				TH1D* h;
				for(int i=0;i<8; ++i){
						tmp = "Proj_dPhi_pp_"+trk_tag[i]+"_"+trk_tag[i+1]+"_Rebin";
						py_dphi[i][4]=(TH1D*)yield_proj_pp_f->Get(tmp);
						symmetrization(py_dphi[i][4]);

						py_dphi_err[i][4]=(TH1D*)py_dphi[i][4]->Clone(Form("py_pp_err_dphi_%d",i));
						h=(TH1D*)yield_proj_syst_f->Get(Form("Yield_dPhi_pp_Syst_Error_%d",i));//pp error
						changeErr(py_dphi_err[i][4], h);

						for(int j=0; j<4; ++j){
								tmp = "Proj_dPhi_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1]+"_Rebin";
								py_dphi[i][j]=(TH1D*)yield_proj_pb_f->Get(tmp);
								py_dphi_err[i][j]=(TH1D*)yield_proj_pb_f->Get(tmp)->Clone(Form("Yield_dPhi_pb_syst_Error_%d_%d",i,j));
								h=(TH1D*)yield_proj_syst_f->Get(Form("Yield_dPhi_Pb_Syst_Error_%d_%d",i,j));
						//		h->Print("ALL");
								changeErr(py_dphi_err[i][j], h);
						}
				}
				//sum them to get total 
				for(int j=0; j<5; ++j){
						//cout<<py_dphi[0][j]->GetName()<<endl;
						py_dphi_all[j]=(TH1D*) py_dphi[0][j]->Clone(Form("py_dphi_all_%d",j));
						py_dphi_err_all[j] = (TH1D*) py_dphi_err[0][j]->Clone(Form("py_dphi_err_all_%d",j)); 
						for(int i=1;i<8; ++i){
								py_dphi_all[j]->Add(py_dphi[i][j]);
								py_dphi_err_all[j]->Add(py_dphi_err[i][j]);
						}
				}

				for(int i=0 ;i<4; ++i){
						py_dphi_diff[i]=(TH1D*)py_dphi_all[i]->Clone(Form("py_dphi_diff_%d",i));
						py_dphi_diff[i]->Add(py_dphi_all[4],-1);
						py_dphi_diff_err[i] =(TH1D*)py_dphi_err_all[i]->Clone(Form("py_dphi_diff_syst_%d",i));
						py_dphi_diff_err[i]->Add(   py_dphi_err_all[4], -1);
				}
		}

		TH1D* py_deta_diff_err[5];
		TH1D* py_deta_diff_all[5];
		TH1D* py_deta_err_all[5];
		TH1D* py_deta_all[5]; 
		TH1D* py_deta    [8][5];
		TH1D* py_deta_err[8][5];
		void prePlot_py_deta(){
				TString tmp;
				TH1D* h;
				for(int i=0;i<8; ++i){
						//cout<<i<<endl;
						tmp = "Proj_dEta_pp_"+trk_tag[i]+"_"+trk_tag[i+1]+"_Rebin";
						py_deta[i][4]=(TH1D*)yield_proj_pp_f->Get(tmp);
						py_deta[i][4]=(TH1D*) etaBin(py_deta[i][4]);
						//cout<<py_deta[i][4]->GetNbinsX()<<endl;
						//						symmetrization(py_deta[i][4]);

						py_deta_err[i][4]=(TH1D*)py_deta[i][4]->Clone(Form("py_pp_err_deta_%d",i));
						h=(TH1D*)yield_proj_syst_f->Get(Form("Yield_dEta_pp_Syst_Error_%d",i));//pp error
						h=(TH1D*) etaBin(h);
						//cout<<h->GetNbinsX()<<endl;
						changeErr(py_deta_err[i][4], h);

						for(int j=0; j<4; ++j){
								tmp = "Proj_dEta_PbPb_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1]+"_Rebin";
								py_deta[i][j]=(TH1D*)yield_proj_pb_f->Get(tmp);
								py_deta[i][j]=(TH1D*) etaBin(py_deta[i][j]);
								py_deta_err[i][j]=(TH1D*)yield_proj_pb_f->Get(tmp)->Clone(Form("Yield_dEta_pb_syst_Error_%d_%d",i,j));
								h=(TH1D*)yield_proj_syst_f->Get(Form("Yield_dEta_Pb_Syst_Error_%d_%d",i,j));
								h=(TH1D*) etaBin(h);
								changeErr(py_deta_err[i][j], h);
						}
				}
				for(int j=0; j<5; ++j){
						py_deta_all[j]=(TH1D*) py_deta[0][j]->Clone(Form("py_deta_all_%d",j));
						py_deta_err_all[j] = (TH1D*) py_deta_err[0][j]->Clone(Form("py_deta_err_all_%d",j)); 
						for(int i=1;i<8; ++i){
								py_deta_all[j]->Add(py_deta[i][j]);
								py_deta_err_all[j]->Add(py_deta_err[i][j]);
						}
				}
				for(int i=0 ;i<4; ++i){
						py_deta_diff_all[i]=(TH1D*)py_deta_all[i]->Clone(Form("py_deta_diff_%d",i));
						py_deta_diff_all[i]->Add(py_deta_all[4],-1);
						py_deta_diff_err[i] =(TH1D*)py_deta_err_all[i]->Clone(Form("py_deta_diff_syst_%d",i));
						py_deta_diff_err[i]->Add(   py_deta_err_all[4], -1);
				}
		}

		TH1D* py_dr_diff_err[5];
		TH1D* py_dr_diff[5];
		TH1D* py_dr_err_all[5];
		TH1D* py_dr_all[5]; 
		TH1D* py_dr    [8][5];
		TH1D* py_dr_err[8][5];
		void prePlot_py_dr(){
				auto s = new signalFactoryBase();
				TString tmp;
				TH1D* h;
				TFile * f= TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/nominal/PbPb_Inclusive_Correlations.root");
				TH2D* h2 = (TH2D*) f->Get("Yield_BkgSub_Cent0_Cent10_Pt100_Pt1000_TrkPt07_TrkPt1");
				for(int i=0;i<8; ++i){
						//cout<<i<<endl;
						tmp = "JetShape2_Yield_BkgSub_Inclusive_pp_"+trk_tag[i]+"_"+trk_tag[i+1];
						py_dr[i][4]=(TH1D*)yield_dr_pp_f->Get(tmp);
						s->doDrPhaseCorrection(h2, py_dr[i][4]);
						py_dr_err[i][4]=(TH1D*)py_dr[i][4]->Clone(Form("py_pp_err_dr_%d",i));
						h=(TH1D*)yield_proj_syst_f->Get(Form("Yield_dr_pp_Syst_Error_%d",i));//pp error
						changeErr(py_dr_err[i][4], h);
						//						py_dr[i][4] = drBin(py_dr[i][4]);
						//						py_dr_err[i][4] = drBin(py_dr_err[i][4]);
						//		cout<<py_deta[i][4]->GetNbinsX()<<endl;

						for(int j=0; j<4; ++j){
								tmp = "JetShape2_Yield_BkgSub_Inclusive_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								py_dr[i][j]=(TH1D*)yield_dr_pb_f->Get(tmp);
								s->doDrPhaseCorrection(h2, py_dr[i][j]);
								py_dr_err[i][j]=(TH1D*)py_dr[i][j]->Clone(Form("Yield_dR_pb_syst_Error_%d_%d",i,j));
								h=(TH1D*)yield_proj_syst_f->Get(Form("Yield_dr_Pb_Syst_Error_%d_%d",i,j));
								changeErr(py_dr_err[i][j], h);
						}
				}
				for(int j=0; j<5; ++j){
						py_dr_all[j]=(TH1D*) py_dr[0][j]->Clone(Form("py_dr_all_%d",j));
						py_dr_err_all[j] = (TH1D*) py_dr_err[0][j]->Clone(Form("py_dr_err_all_%d",j)); 
						for(int i=1;i<8; ++i){
								py_dr_all[j]->Add(py_dr[i][j]);
								py_dr_err_all[j]->Add(py_dr_err[i][j]);
						}
				}
				for(int i=0 ;i<4; ++i){
						py_dr_diff[i]=(TH1D*)py_dr_all[i]->Clone(Form("py_dr_diff_%d",i));
						py_dr_diff[i]->Add(py_dr_all[4],-1);
						py_dr_diff_err[i]= (TH1D*)py_dr_err_all[i]->Clone(Form("py_dr_diff_syst_%d",i));
						py_dr_diff_err[i]->Add(py_dr_err_all[4], -1);
				}
		}

		TH1D* js_dr_diff_err[5];
		TH1D* js_dr_diff[5];
		TH1D* js_dr_err_all[5];
		TH1D* js_dr_all[5]; 
		TH1D* js_dr    [9][5];
		TH1D* js_dr_err[9][5];

		void prePlot_js_dr(){
				auto s = new signalFactoryBase();
				TString tmp;
				TH1D* h;
				TFile * f= TFile::Open("/Users/tabris/cmsProjects/inclusiveJetTrackCorrelation2015/dataSet/correlation/nominal/PbPb_Inclusive_Correlations.root");
				TH2D* h2 = (TH2D*) f->Get("Yield_BkgSub_Cent0_Cent10_Pt100_Pt1000_TrkPt07_TrkPt1");
				for(int i=0;i<9; ++i){
						tmp = "JetShape2_Yield_BkgSub_pTweightedInclusive_pp_"+
								trk_tag[i]+"_"+trk_tag[i+1];
						js_dr[i][4]=(TH1D*)js_dr_pp_f->Get(tmp);
//						s->doDrPhaseCorrection(h2, js_dr[i][4]);
						js_dr_err[i][4]=(TH1D*)js_dr[i][4]->Clone(Form("JS_pp_err_dr_%d",i));
						h=(TH1D*)js_syst_f->Get(Form("js_dr_pp_Syst_Error_%d",i));//pp error
						changeErr(js_dr_err[i][4], h);
						for(int j=0; j<4; ++j){
								tmp = "JetShape2_Yield_BkgSub_pTweightedInclusive_"+cent_tag[j]+"_"+cent_tag[j+1]+"_"+trk_tag[i]+"_"+trk_tag[i+1];
								js_dr[i][j]=(TH1D*)js_dr_pb_f->Get(tmp);
//								s->doDrPhaseCorrection(h2, js_dr[i][j]);
								js_dr_err[i][j]=(TH1D*)js_dr[i][j]->Clone(Form("js_dR_pb_Error_%d_%d",i,j));
								h=(TH1D*)js_syst_f->Get(Form("js_dr_Pb_Syst_Error_%d_%d",i,j));
								changeErr(js_dr_err[i][j], h);
						}
				}
				for(int j=0; j<5; ++j){
						js_dr_all[j]=(TH1D*) js_dr[0][j]->Clone(Form("JS_dr_all_%d",j));
						js_dr_err_all[j] = (TH1D*) js_dr_err[0][j]->Clone(Form("JS_dr_err_all_%d",j)); 
						for(int i=1;i<9; ++i){
								js_dr_all[j]->Add(js_dr[i][j]);
								js_dr_err_all[j]->Add(js_dr_err[i][j]);
						}
						//cout<<js_dr_all[j]->GetName()<<endl;
				}
		}


		TH1D* ppyield;
		TH1D* pbyield[4];
		TH1D* yield_diff[4];
		TH1D* ppyield_err;
		TH1D* pbyield_err[4];
		TH1D* yield_diff_err[4];

		void prePlot_yield_integral(){
				prePlot_py_dr();
				int x1 = (py_dr[0][0]->GetXaxis())->FindBin(0.0);
				int x2 = (py_dr[0][0]->GetXaxis())->FindBin(1);
				Double_t bins[] = {0., 0.7, 1, 2, 3, 4, 8, 12, 16, 20 , 20.3};
				ppyield = new TH1D("pp_yield", "",  10, bins);
				ppyield_err = new TH1D("pp_yield_err", "", 10, bins);
				for(int k=0; k<8; ++k){
						Double_t err;
						Double_t bc = py_dr[k][4]->IntegralAndError(x1, x2, err, "width");
						ppyield->SetBinContent(k+2, bc);
						ppyield->SetBinError  (k+2, err);
						bc = py_dr_err[k][4]->IntegralAndError(x1, x2, err, "width");
						err=0; 
						for(int l=x1; l<x2+1; ++l){
								err =err+ py_dr_err[k][4]->GetBinError(l)*py_dr_err[k][4]->GetBinWidth(l);
						}
						ppyield_err->SetBinContent(k+2, bc);
						ppyield_err->SetBinError  (k+2, err);
				}
				for(int i=0; i<4; i++){
						pbyield[i] = (TH1D*) ppyield->Clone(Form("pb_yield_%d",i));
						pbyield_err[i] = (TH1D*) ppyield->Clone(Form("pb_yield_err_%d",i));
						for(int k=0; k<8; ++k){
								Double_t err;
								Double_t bc = py_dr[k][i]->IntegralAndError(x1, x2, err, "width");
								pbyield[i]->SetBinContent(k+2, bc);
								pbyield[i]->SetBinError  (k+2, err);
								bc = py_dr[k][i]->IntegralAndError(x1, x2, err, "width");
								err=0; 
								for(int l=x1; l<x2+1; ++l){
										err= err+py_dr_err[k][i]->GetBinError(l)*py_dr_err[k][i]->GetBinWidth(l);
								}
								pbyield_err[i]->SetBinContent(k+2, bc);
								pbyield_err[i]->SetBinError  (k+2, err);
						}
						yield_diff[i]=(TH1D*) pbyield[i]->Clone(Form("yield_diff_%d",i));
						yield_diff[i]->Add(ppyield,-1);
						yield_diff_err[i]=(TH1D*) pbyield_err[i]->Clone(Form("yield_diff_err_%d",i));
						yield_diff_err[i]->Add(ppyield_err,-1);
				}
		}

		TGraphAsymmErrors* py_diff_err_ref[4];
		TH1D* py_diff_err[4];
		void getReference(){
				TString trk2_tag[] = {"TrkPt05", "TrkPt1", "TrkPt2","TrkPt3","TrkPt4","TrkPt8","TrkPt300"};
				TH1D* h ;
				TString tmp;
				Double_t bins[] = {1, 2, 3, 4, 8};
				for(int j=0; j<4; j++){
						tmp = "Integrated_Yield_Phi_Syst_"+cent_tag[j]+"_"+cent_tag[j+1];
						py_diff_err_ref[j] = (TGraphAsymmErrors*)run1_ref_f->Get(tmp);
				}

		}

}
