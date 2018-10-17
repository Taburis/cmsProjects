
#ifndef doublePad_H
#define doublePad_H
#include "TPad.h"
#include "TLegend.h"
#include "TLatex.h"

class doublePad{
	public:
		doublePad(TPad* pad){
			hpad=pad;
			hpad->Divide(1,2, 0, 0);
			hpad->SetMargin(0.125, 0.05, 0.1, 0.075);
			xmin=1; xmax=0;
			y1min=1; y1max=0;y2min =1; y2max=0;
			l1=new TLegend(0.38, 0.7, 0.65,0.95);
			l2=new TLegend(0.38, 0.8, 0.9,0.95);
			tx=new TLatex();
		};
		void addRatioPair(TH1* hnum, TH1* hden, Color_t color, TString leg1, TString leg2, TString legr="");
		void Draw();
		~doublePad();
		void setRange(float, float , float y1=1, float y2=0, float y3=1, float y4=0);
		void drawLegend();
		void setXtitle(TString title){ xtitle = title;};
		void setYtitle(TString title1, TString title2=""){ y1title = title1; y2title = title2;};
		void putText(TString, float x=0.6, float=0.7, int n=1);
	public:
		TPad* hpad;
		std::vector<TH1*> hratio;
		std::vector<TString> hr_leg;
		std::vector<TString> hleg;
		std::vector<TH1*> hist;
		TLegend* l1;
		TLegend* l2;
		float xmin,xmax;
		float y1min, y1max, y2min, y2max;
		TString xtitle, y1title, y2title;
		TH1 *h1frame, * h2frame;	
		TLatex* tx;
};

void doublePad::addRatioPair(TH1* hnum, TH1* hden, Color_t color, TString leg1, TString leg2, TString legr){
	int cache1=0, cache2=0;
	for(auto &it : hist){
		if(it == hnum) cache1++; //test if cached
		if(it == hden) cache2++; //test if cached
	}
	if( !cache2 ) {
		hist.push_back(hden);
		hleg.push_back(leg2);
	}
	if( !cache1 ) {
		hist.push_back(hnum);
		hleg.push_back(leg1);
	}
	if(legr !="") hr_leg.push_back(legr);
	else hr_leg.push_back(leg1+"/"+leg2);
	hratio.push_back((TH1*)hnum->Clone(Form("dc_ratio_%d",int(hratio.size()))));
	hratio.back()->Divide(hden);
	hratio.back()->SetLineColor(color);
	hratio.back()->SetMarkerStyle(20);
	hratio.back()->SetMarkerColor(color);
	hratio.back()->SetTitle("");
}

doublePad::~doublePad(){
	for(auto &it: hist) delete it;
	for(auto &it: hratio) delete it;
	hist.clear();
	hratio.clear();
	delete l1;
	delete l2;
}

void doublePad::Draw(){
	if( hist.size()==0) return;
	hpad->cd(1);
	if(xmin<xmax){
		hist[0]->SetAxisRange(xmin, xmax, "X");
		hratio[0]->SetAxisRange(xmin, xmax, "X");
	}
	if(y1min<y1max) hist[0]->SetAxisRange(y1min, y1max, "Y");
	if(y2min<y2max)  hratio[0]->SetAxisRange(y2min, y2max, "Y");
	h1frame = hist[0];
	h2frame = hratio[0];
	hist[0]->GetYaxis()->SetTitle(y1title);
	hist[0]->SetStats(0);
	for(auto &it: hist) it->Draw("same");
	gPad->SetTicks(1,1);
	hpad->cd(2);
	hratio[0]->SetStats(0);
	hratio[0]->GetYaxis()->SetTitle(y2title);
	hratio[0]->GetXaxis()->SetTitle(xtitle);
	for(auto &it: hratio) it->Draw("same");
	gPad->SetTicks(1,1);
}

void doublePad::setRange(float x1, float x2, float y3, float y4, float y1, float y2 ){
	xmin= x1; xmax=x2; y1min =y1;y1max=y2; y2min =y3; y2max=y4;
}
void doublePad::drawLegend(){
	l1->SetFillColorAlpha(kWhite,0);
	l1->SetLineColorAlpha(kWhite,0);
	l2->SetFillColorAlpha(kWhite,0);
	l2->SetLineColorAlpha(kWhite,0);
	for(int j=0 ; j<hist.size(); ++j){
		l1->AddEntry(hist[j], hleg[j]);
	}
	for(int j=0 ; j<hratio.size(); ++j){
		l2->AddEntry(hratio[j], hr_leg[j]);
	}
	hpad->cd(1);
	l1->Draw();
	hpad->cd(2);
//	l2->Draw();
}
void doublePad::putText(TString txt, float x, float y, int n){
	hpad->cd(n);
	tx->DrawLatexNDC(x, y, txt);
}
#endif
