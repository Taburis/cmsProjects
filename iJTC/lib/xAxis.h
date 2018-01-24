

#ifndef xAxis_H
#define xAxis_H

#ifndef xAlgo_H
#include "xAlgo.h"
#endif

class xAxis{
	public : xAxis (){
		 };
	 	 xAxis (int nbin, float *bin);
//	 	 xAxis (int nbin, const float *bin);
	 	 xAxis (int nbin, float xmin, float xmax);
		 int findBin(float x);
	
	public : 
		 float * bins;
		 int     nbin;
		 float   xmax;
		 float   xmin;
		 char*   name;
};

xAxis::xAxis(int n, float *bin){
	bins = new float[n+1];
	nbin = n;
	for(int i=0; i<n+1; ++i) bins[i]=bin[i];
}
//xAxis::xAxis(int n, const float *bin){
//	bins = new float[n+1];
//	nbin = n;
//	for(int i=0; i<n+1; ++i) bins[i]=bin[i];
//}

xAxis::xAxis(int n, float xmin, float xmax){
	float width =(xmax-xmin)/n;
	bins = new float[n+1];
	nbin = n;
	for(int i=0; i<n+1; ++i) bins[i]=xmin+i*width;
}

int xAxis::findBin(float x){
	return xAlgo::BinarySearch(nbin+1, bins, float(x));
}

#endif

