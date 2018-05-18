
#ifndef histPlayer_H
#define histPlayer_H
#include "xAxis.h"

template <typename T>
class histPlayer2D {
		public : histPlayer2D(){};
				 histPlayer2D(T** h, int x, int y): n1(x), n2(y){ hist = h;};
				 ~histPlayer2D();
				 void read(TFile* f, TString hname, int x, int y);
				 void saveHist();
				 int* rebinIndex(int nin, float *inbin, int nout, float *outbin);
				 T** rebin2D(TString name, int nx0, float *x0bin, int ny0, float *y0bin, int nx, float *xbin, int ny, float *ybin);
				 void selfRebin2D(TString name, int nx0, float *x0bin, int ny0, float *y0bin, int nx, float *xbin, int ny, float *ybin);
				 T** binary_operation(TString name, T** h1, T** h2, TString opt);
				 unsigned int flatten(int x, int y){ return x+y*n1;};
				 T* at(int x, int y){ return hist[x+y*n1]; };
		public :
				 TString name; 
				 int n1, n2;
				 T** hist;
				 TFile* hf;
};

template <typename T>
void histPlayer2D<T>::read(TFile *f, TString hname, int x, int y){
		name = hname;
		n1 = x; n2 = y;
		hist = new T*[n1*n2];
		hf=f;
		for(int i=0; i<n1; ++i){
				for(int j=0;j<n2;++j) hist[flatten(i, j)] = (T*)f->Get(name+Form("_%d_%d",i,j));
		}
}

template<typename T>
T** histPlayer2D<T>::binary_operation(TString name, T** h1, T** h2, TString opt){
		T** h = new T*[n1*n2]; 
		for(int i=0; i<n1; ++i){
				for(int j=0; j<n2; ++j){
						//cout<<i<<", "<<j<<endl;
						h[i+n1*j]=(T*) h1[i+n1*j]->Clone(name + Form("_%d_%d",i,j));
						if(opt == "ratio") 
								h[i+n1*j]->Divide(h2[i+n1*j]);
						else if( opt== "add") 
								h[i+n1*j]->Add(h2[i+n1*j]);
						else if( opt== "multiply") 
								h[i+n1*j]->Multiply(h2[i+n1*j]);
						else if( opt== "diff") {
								h[i+n1*j]->Add(h2[i+n1*j], -1);
						}
						else if( opt== "binomialRatio") {
								h[i+n1*j]->Divide(h[i+n1*j], h2[i+n1*j], 1, 1, "B");
						}
						else cout<<"no defined operation: "<<opt<<endl;
				}
		}
		return h;
}

template <typename T>
void histPlayer2D<T>::saveHist(){
		for(int i=0; i<n1; ++i){
				for(int j=0;j<n2;++j) hist[flatten(i, j)]->Write();
		}
}

template <typename T>
int* histPlayer2D<T>::rebinIndex(const int nin, float *inbin, const int nout, float *outbin){
		auto ina = xAxis(nin, inbin);
		auto oua = xAxis(nout, outbin);
		int* rebinIndx= new int[nout+1];
		for(int i=0; i<nout+1; ++i){
				if(inbin[ina.findBin(outbin[i])] == outbin[i]){ 
						rebinIndx[i]=ina.findBin(outbin[i]);
						continue;
				}
				std::cout<<"WARNING: the bin scheme isn't match!"<<std::endl;
				return 0;
		}
		return rebinIndx;
}

template <typename T>
T** histPlayer2D<T>::rebin2D(TString name, int nx0, float *x0bin, int ny0, float *y0bin, int nx, float *xbin, int ny, float *ybin){
		T** newh = new T*[nx*ny];
		int* xindex = rebinIndex(nx0, x0bin, nx, xbin);
		int* yindex = rebinIndex(ny0, y0bin, ny, ybin);
		if(xindex == 0 || yindex ==0 ){
				std::cout<<"WARNING: Rein abort!"<<std::endl;
				return 0;
		}

		for(int j=0; j<ny; ++j){
				for(int i=0; i<nx; ++i){
						for(int k=xindex[i]; k<xindex[i+1]; ++k){
								for(int l=yindex[j]; l<yindex[j+1]; ++l){
				//						cout<<"i= "<<i<<", j= "<<j<<", k="<<k<<", l="<<l<<endl;
										if(k==xindex[i] && l==yindex[j])
											   	newh[i+j*nx]=(T*)hist[xindex[i]+nx0*yindex[j]]->Clone(name+Form("_%d_%d", i, j));
										else 
												newh[i+j*nx]->Add(hist[k+l*nx0]);
								}
						}
				}
		}
		return newh;
}

template <typename T>
void histPlayer2D<T>::selfRebin2D(TString hname, int nx0, float *x0bin, int ny0, float *y0bin, int nx, float *xbin, int ny, float *ybin){
		T** htmp = rebin2D(hname, nx0, x0bin, ny0, y0bin,  nx, xbin, ny, ybin);
		n1 = nx; n2 = ny;
		name = hname;
		for(int i=0; i<n1*n2; ++i){
				delete (T*)hist[i];
		}
		hist = htmp;
}

template <typename T>
histPlayer2D<T>::~histPlayer2D(){
		for(int i=0; i<n1*n2; ++i){
				delete (T*)hist[i];
		}
		delete[] hist;
		hist = 0;
}

#endif
