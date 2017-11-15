

#ifndef xCanvas_H
#define xCanvas_H

class xCanvas : public TCanvas{
	public :
		xCanvas(TString cname, TString title, float w, float h): 
			name(cname), cw(w), ch(h),TCanvas(cname, title, w, h){
			}
		void divide(int, int);
		void CD(int);
	public :
		TString name;
		TPad **pad;
		float cw, ch;
};

void xCanvas::divide(int nx, int ny){
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

void xCanvas::CD(int i){
	pad[i-1]->cd();
}

#endif
