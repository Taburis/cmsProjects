
#ifndef TH2DBlock_H
#define TH2DBlock_H
#include "histPlayer.h"

class TH2DBlock : public histPlayer2D<TH2D>{
		public: TH2DBlock(TString name, int n1, int n2, int nx, double x1, double x2, int ny, double y1, double y2);
				TH2DBlock();
				TH2DBlock* binaryOperation(TString name, TH2DBlock* th2, TString opt);
				~TH2DBlock(){}
};

TH2DBlock::TH2DBlock(TString hname, int n1, int n2, int nx, double x1, double x2, int ny, double y1, double y2):
		histPlayer2D<TH2D>()
{
		name = hname;
		hist = new TH2D*[n1*n2];
		for(int i=0; i<n1; ++i){
				for(int j=0;j<n2;++j){
					   	hist[flatten(i, j)] = new TH2D(name+Form("_%d_%d",i,j), "", nx, x1, x2, ny, y1, y2);
						hist[flatten(i, j)]->Sumw2();
				}
		}

}

TH2DBlock::TH2DBlock(): histPlayer2D<TH2D>(){
}

TH2DBlock* TH2DBlock::binaryOperation(TString name, TH2DBlock* th2, TString opt){
		auto hr = new TH2DBlock();
		hr->hist = binary_operation(name, hist, th2->hist, opt);
		hr->name=name;
		hr->n1 = n1, hr->n2=n2;
		return hr;
}

#endif
