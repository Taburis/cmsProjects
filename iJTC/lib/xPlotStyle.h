

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

class mCanvasLoose : public mCanvasBase {
		public: 
				mCanvasLoose(const char *name, const char *title , int n, int m , float width = 300, float height=275);
				void CD(int i, int j){this->cd((i-1)*ncol+j);}
				void drawHist(TH1* h, int i, int j, TString opt ="");
				void drawHist(TH1* h, int i, TString opt ="");
				void histStyle(TH1* h);
		public:
				TLine *tl; TBox *box;
};

mCanvasLoose::mCanvasLoose(const char *name, const char *title , int x, int y , float width, float height):
		mCanvasBase(name, title, x, y, width, height, 0, 0)
{
		//		style()->cd();
		this->SetMargin(0.5, 0.5, 0.4, 0.40);
		gStyle->SetOptStat(0);
		this->Divide(ncol,nrow);
		tl = new TLine();
		box = new TBox(); box->SetFillColorAlpha(kGray+1, 0.5);
}

void mCanvasLoose::drawHist(TH1* h, int i, int j, TString opt){
		CD(i,j);
		gPad->SetMargin(0.15,0.04, 0.15, 0.08);
		histStyle(h);
		h->Draw(opt+"same");
}

void mCanvasLoose::drawHist(TH1* h, int i, TString opt){
		cd(i);
		gPad->SetMargin(0.15,0.04, 0.15, 0.08);
		histStyle(h);
		h->Draw(opt+"same");
}

void mCanvasLoose::histStyle(TH1 *h){
		h->GetXaxis()->SetLabelSize(0.07);
		h->GetYaxis()->SetLabelSize(0.07);
		h->GetXaxis()->CenterTitle();
		h->GetXaxis()->SetTitleOffset(0.8);
		h->GetXaxis()->SetTitleSize(0.085);
		h->GetYaxis()->SetTitleSize(0.09);
}

class doublePanelFig : public mCanvasBase {
		public :
				doublePanelFig(const char *name, const char *title, int n, int m, float r=0.37);
				int at(int i, int j, int n=0){ return 2*(i-1)*ncol+2*j+n-2;}
				void addHist(TH1* h, int i, int j, int n=0, TString opt = "");
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
		mCanvasBase(name, title, x, y, 350, 450, 0, 0)
{
		style()->cd();
		//		this->SetMargin(0.5, 0.5, 0.4, 0.40);
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
						pad[at(i,j, 0)]->SetTopMargin(0.05);
						pad[at(i,j, 0)]->SetLeftMargin(0.14);
						pad[at(i,j, 1)]->SetLeftMargin(0.14);
						pad[at(i,j, 0)]->SetBottomMargin(0);
						pad[at(i,j, 1)]->SetTopMargin(0);
						pad[at(i,j, 0)]->SetRightMargin(0.02); pad[at(i,j, 0)]->SetTopMargin(0.03);
						pad[at(i,j, 1)]->SetRightMargin(0.02); 
						pad[at(i,j, 1)]->SetBottomMargin(0.25); 
						pad[at(i,j,0)]->Draw(); pad[at(i,j,1)]->Draw();
						//cout<<pad[at(i,j,1)]->GetName()<<endl;
				}
		}
}

void doublePanelFig::addHist(TH1* h, int i, int j, int n, TString opt){
		pad[at(i,j, n)]->cd();
		histStyle(h, n);
		h->Draw("same"+opt);
}


void doublePanelFig::histStyle(TH1* h, int n){
		if(n==1){
				h->GetXaxis()->SetLabelSize(0.11);
				h->GetYaxis()->SetLabelSize(0.12);
				h->GetXaxis()->CenterTitle();
				h->GetXaxis()->SetTitleSize(0.15);
				h->GetXaxis()->SetTitleOffset(0.7);
				h->GetYaxis()->SetTitleOffset(0.5);
				h->GetYaxis()->SetTitleSize(0.12);
				h->GetYaxis()->CenterTitle();
		}
		else{
				h->GetYaxis()->SetLabelSize(0.08);
				h->GetXaxis()->SetLabelSize(0.08);
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
