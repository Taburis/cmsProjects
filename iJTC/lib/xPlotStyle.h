

#ifndef xPlotStyle
#define xPlotStyle

#ifndef ROOT_TCanvas
#include "TCanvas.h"
#endif

class mCanvasBase : public TCanvas{
		public : 
				mCanvasBase(const char * name, const char *title, int x, int y, int w, int h, int wm = 100, int hm=80);
		public :
				int ncol, nrow;
};

mCanvasBase::mCanvasBase(const char * name, const char *title, int x, int y, int w, int h, int wm, int hm) : 
		nrow(x), ncol(y), TCanvas(name, title , y*w+wm, x*h+hm)
{}

class doublePanelFig : public mCanvasBase {
		public :
				doublePanelFig(const char *name, const char *title, int n, int m, float r=0.37);
				int at(int i, int j, int n=0){ return 2*(i-1)*ncol+2*j+n-2;}
				void addHist(TH1* h, int i, int j, int n=0);
				TStyle * style();
				void histStyle(TH1*, int n);
				void CD(int i, int j, int n) { pad[at(i,j,n)]->cd();}
				void drawShadowArea(int i, int j, float x1, float x2, float y1, float y2);
				void draw95Area(int i, int j, float x1, float x2){return drawShadowArea(i,j, x1, x2,0.5, 1.5);}
		public : 
				TLine *tl; TBox *box;
				TPad** pad;
				TH1* hr;
};

doublePanelFig::doublePanelFig(const char * name, const char *title, int x, int y, float r):
		mCanvasBase(name, title, x, y, 300, 400)
{
		style()->cd();
		this->SetMargin(0.5, 0.5, 0.2, 0.10);
		this->Divide(ncol,nrow);
		tl = new TLine();
		pad= new TPad*[2*ncol*nrow];
		TPad* tp; TString sname;
		box = new TBox(); box->SetFillColorAlpha(kGray+1, 0.5);
		for(int i=1; i<nrow+1; ++i){
				for(int j=1; j<ncol+1; ++j){
						//cout<<i<<", "<<j<<endl;
						sname = (this->cd(ncol*(i-1)+j))->GetName();
						pad[at(i,j, 0)] = new  TPad(sname+"_0", "", 0.0, r, 1, 1);
						pad[at(i,j, 1)] = new  TPad(sname+"_1", "", 0.0, 0.0, 1, r);
						pad[at(i,j, 0)]->SetBottomMargin(0);
						pad[at(i,j, 1)]->SetTopMargin(0);
						pad[at(i,j, 0)]->SetRightMargin(0.05); pad[at(i,j, 0)]->SetTopMargin(0.03);
						pad[at(i,j, 1)]->SetRightMargin(0.05); 
						pad[at(i,j,0)]->Draw(); pad[at(i,j,1)]->Draw();
						//cout<<pad[at(i,j,1)]->GetName()<<endl;
				}
		}
}

void doublePanelFig::addHist(TH1* h, int i, int j, int n){
		pad[at(i,j, n)]->cd();
		histStyle(h, n);
		h->Draw("same");
}
void doublePanelFig::histStyle(TH1* h, int n){
		if(n==1){
				h->GetXaxis()->SetLabelSize(0.12);
				h->GetYaxis()->SetLabelSize(0.12);
				h->GetXaxis()->CenterTitle();
				h->GetXaxis()->SetTitleSize(0.15);
		}
		else{
		h->GetYaxis()->SetLabelSize(0.08);
		}
}
void doublePanelFig::drawShadowArea(int i, int j, float x1, float x2, float y1, float y2){
		CD(i,j,1);
		box->DrawBox(x1, y1, x2, y2);
}

TStyle * doublePanelFig::style(){
		auto st = new TStyle(*gStyle);
		st->SetOptStat(0);
		st->SetStatStyle(0);
		st->SetFrameFillColor(0);
//		st->SetFrameFillStyle(1001);
		st->SetFrameLineWidth(0);
		st->SetCanvasBorderMode(0);
		st->SetPadBorderMode(0); st->SetPadColor(0); st->SetPadBorderSize(0);
		st->SetOptTitle(0);
		st->SetPadTickX(1); st->SetPadTickY(1);
		st->SetTitleX(2);
		return st;
//		st->SetStats(0);
}
#endif
