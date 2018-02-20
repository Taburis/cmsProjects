
#ifndef JTCUtility_H
#define JTCUtility_H
#ifndef signalFactoryBase_H
#include "signalFactoryBase.h"
#endif
#ifndef JTCSignalProducer_H
#include "JTCSignalProducer.h"
#endif

class JTCUtility :public signalFactoryBase{
		public : 
				JTCUtility() : signalFactoryBase(){};
				void getJFF(TH2D* hrec, TH2D* hgen);
		public : 
				TF1 *fjff_deta=0, *fjff_dphi=0;
				TH1D* hjff_deta=0, *hjff_dphi=0;
				JTCSignalProducer *sp=0;
};

void JTCUtility::getJFF(TH2D* hrec, TH2D* hgen){
	int n1 = hrec->GetYaxis()->FindBin(-1);	
	int n2 = hrec->GetYaxis()->FindBin(1);	
	hjff_deta = hrec->ProjectionX("jff_deta", n1, n2);	
	TH1D* gen = hgen->ProjectionX("gen_deta", n1, n2);
	hjff_deta->Add(gen,-1);
	hjff_deta->Rebin(4); hjff_deta->Scale(0.25);
	n1 = hjff_deta->FindBin(-2.5);
	n2 = hjff_deta->FindBin(-1.49);
	int nsum = n2-n1;
	float p0 = hjff_deta->Integral(n1, n2);
	n1 = hjff_deta->FindBin(1.5);
	n2 = hjff_deta->FindBin(2.49);
	nsum = n2-n1+nsum;
	p0 = (p0+hjff_deta->Integral(n1,n2))/nsum;
	n1 = hjff_deta->FindBin(-.15);
	n2 = hjff_deta->FindBin(.149);
	float p1 = hjff_deta->Integral(n1,n2)/(n2-n1);
	fjff_deta = new TF1("jff_deta", "[0]+[1]*TMath::Exp([2]*x*x)", -2.5, 2.5);
	fjff_deta->SetParameter(0, p0);
	fjff_deta->SetParameter(1, p1);
	hjff_deta->Fit(fjff_deta,"","", -2.5, 2.49);

	n1 = hrec->GetXaxis()->FindBin(-1);	
	n2 = hrec->GetXaxis()->FindBin(1);	
	hjff_dphi = hrec->ProjectionY("jff_dphi", n1, n2);	
	gen = hgen->ProjectionY("gen_dphi", n1, n2);
	hjff_dphi->Rebin(4); hjff_dphi->Scale(0.25);
	hjff_dphi->Add(gen,-1);
	n1 = hjff_dphi->FindBin(-1.3);
	n2 = hjff_dphi->FindBin(-1.);
	nsum = n2-n1;
	p0 = hjff_dphi->Integral(n1, n2);
	n1 = hjff_dphi->FindBin(1.);
	n2 = hjff_dphi->FindBin(1.29);
	nsum = n2-n1+nsum;
	p0 = (p0+hjff_dphi->Integral(n1,n2))/nsum;
	n1 = hjff_dphi->FindBin(-.15);
	n2 = hjff_dphi->FindBin(.149);
	p1 = hjff_dphi->Integral(n1,n2)/(n2-n1);
	fjff_dphi = new TF1("jff_dphi", "[0]+[1]*TMath::Exp([2]*x*x)", -1.5, 1.5);
	fjff_dphi->SetParameter(0, p0);
	fjff_dphi->SetParameter(1, p1);
	hjff_dphi->Fit(fjff_dphi,"0","", -1.5, 1.49);
}

#endif
