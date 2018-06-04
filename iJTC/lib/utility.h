
#ifndef xUtility
#define xUtility

#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif

#ifndef ROOT_TString
#include "TString.h"
#endif

#ifndef ROOT_TLegend
#include "TLegend.h"
#endif

#ifndef ROOT_TLine
#include "TLine.h"
#endif

#ifndef ROOT_TLatex
#include "TLatex.h"
#endif

#ifndef ROOT_TString
#include "TString.h"
#endif
#include <iomanip> // setprecision
#include <sstream>
#include <stdio.h>

class padPolisher{
	public : 
		padPolisher(){
			tl=new TLegend();
			tx=new TLatex();
			line=new TLine();
		};
		void cent4Labels(TCanvas *, float x=0.65, float y=0.8);
		void centLabel(int k, TPad *pad=NULL, float x=0.65, float y=0.8);
		void ptLabel( int k, float *ptbins, TString name="p_{T}", TPad *pad = NULL, float x=0.6, float y=0.8);
		~padPolisher(){
			delete tl;
			delete tx;
			delete line;
		}

	public :
		TLegend* tl;
		TLatex* tx;
		TLine * line;
		std::stringstream stream;
		TString str;
};

void padPolisher::cent4Labels(TCanvas *c, float x, float y){
	TString l[4]= {"Cent. 0-10%", "Cent. 10-30%","Cent. 30-50%","Cent. 50-100%"};
	for(int i=0;i<4;++i){
		c->cd(4-i);
		tx->DrawLatexNDC(x,y,l[i]);
	}
}

void padPolisher::ptLabel(int k, float *ptbins, TString name, TPad *pad, float x, float y){
	if(pad!=NULL) pad->cd();
	int np1=0, np2=0;
	if(ceil(ptbins[k])!=floor(ptbins[k])) np1 = 1;
	if(ceil(10*ptbins[k])!=floor(10*ptbins[k])) np1 = 2;
	if(ceil(ptbins[k+1])!=floor(ptbins[k+1])) np2 = 1;
	if(ceil(10*ptbins[k+1])!=floor(10*ptbins[k+1])) np2 = 2;
	stream.str("");
	stream<<fixed<<setprecision(np1)<<ptbins[k]<<"< "<<name<<"< "<<fixed<<setprecision(np2)<<ptbins[k+1];
	str=stream.str();
	tx->DrawLatexNDC(x,y,str);	
}

void padPolisher::centLabel(int k, TPad *pad, float x, float y){
	TString l[4]= {"0-10%", "10-30%","30-50%","50-100%"};
	if(pad!=NULL) pad->cd();
	str="cent:"+l[k];
	tx->DrawLatexNDC(x,y,str);
}

TH1D* xRebin(TH1D* h, int nbin, const float * bin){
	TString hname = h->GetName();
	TString htitle = h->GetTitle();
	TH1D* hnew = new TH1D(hname+"_temp",htitle, nbin, bin);
	std::vector<int > node; // the position of the new bin node in the old bin scheme
	int Nbin = h->GetNbinsX();
	//test if new bin edge is contained in the old binning scheme
	if(bin[0] != h->GetBinLowEdge(1) || bin[nbin] != h->GetBinLowEdge(Nbin+1)){
		std::cout<<"Error in Rebin: rebin range changed!"<<std::endl;
		return NULL;
	}
	node.push_back(1);
	for(int j=1; j<nbin ; ++j){
		for(int k=node.back()+1 ; k<Nbin; ++k){
			if( h->GetBinLowEdge(k) == bin[j]){
				node.push_back(k);
				break;
			}
		}
	}
	node.push_back(Nbin+1);
	for(int j=0; j<node.size()-1; ++j){
//		cout<<node.at(j)<<endl;
		float var = 0;
		float cont =0;
		float width=0;
		for(int k=node.at(j); k<node.at(j+1); ++k){
			cont += h->GetBinContent(k)*h->GetBinWidth(k);
			var  += pow(h->GetBinError(k)*h->GetBinWidth(k), 2);
		}
		width = hnew->GetBinLowEdge(j+2)-hnew->GetBinLowEdge(j+1);
		cont=cont/width;
		var=var/pow(width,2);
		hnew->SetBinContent(j+1, cont);
		hnew->SetBinError(j+1, pow(var, 0.5));
	} 
	hnew ->SetName(hname);
	*h=*hnew;
	return hnew;
}

void h1Symm(TH1D &h, float x , float y){
	int nl = h.FindBin(x);
	int nr = h.FindBin(y);
	int midl = floor(float(nr+nl)/2); 
	int midr =  ceil(float(nr+nl)/2); 
	float err, cont;
	for(int i=0; i<midl; ++i){
		cont = (h.GetBinContent(midl-i)+h.GetBinContent(i+midr))/2;
		err  = (h.GetBinError(midl-i)  +h.GetBinError(midr+i))/2;
		h.SetBinContent(midl-i, cont);
		h.SetBinContent(midr+i, cont);
		h.SetBinError(midl-i, err);
		h.SetBinError(midr+i, err);
	}
	return;
}

double findDr(double eta1, double phi1, double eta2, double phi2){
	double dphi = phi1 - phi2;
	while (dphi > M_PI) dphi-= 2*M_PI;
	while (dphi <= -M_PI) dphi += 2*M_PI;
	return sqrt(pow(dphi,2)+pow(eta1-eta2,2));
}

#endif
