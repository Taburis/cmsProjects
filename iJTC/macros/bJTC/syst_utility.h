void add_frac_error(TH1D* h , TH1D* herr){
		for(int i=1; i<h->GetNbinsX()+1; ++i){
				float err0 = h->GetBinError(i);
				float err1 = h->GetBinContent(i)*TMath::Abs(1-herr->GetBinContent(i));
				//float err1 = h->GetBinContent(i)*TMath::Abs(1-herr->GetBinContent(i))/herr->GetBinContent(i);
				float err = pow(err0*err0+err1*err1, .5);
				h->SetBinError(i, err);
		}
}
void scale_error(TH1D* h , float c){
		for(int i=1; i<h->GetNbinsX()+1; ++i){
				float err0 = h->GetBinError(i);
				//float err1 = h->GetBinContent(i)*TMath::Abs(1-herr->GetBinContent(i))/herr->GetBinContent(i);
				h->SetBinError(i, err0*c);
		}
}
void add_frac_error(TH1D* h , float c){
		for(int i=1; i<h->GetNbinsX()+1; ++i){
				float err0 = h->GetBinError(i);
				float err1 = h->GetBinContent(i)*c;
				//float err1 = h->GetBinContent(i)*TMath::Abs(1-herr->GetBinContent(i))/herr->GetBinContent(i);
				float err = pow(err0*err0+err1*err1, .5);
				h->SetBinError(i, err);
		}
}

void sign_err(TH1D** h, TH1D** herr){
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						for(int k=1; k<h[i+j*nPt]->GetNbinsX()+1; ++k){
								h[i+j*nPt]->SetBinError(k, herr[i+j*nPt]->GetBinError(k));
						}
				}
		}
}


void add_frac_error(TH1D** h , float c){
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						add_frac_error(h[i+nPt*j], c);
				}
		}
}

void add_frac_error(TH1D** h , TH1D** herr){
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						add_frac_error(h[i+nPt*j], herr[i+nPt*j]);
				}
		}
}

TH1D* get_bkg_dr_err(TString name, TH1D* h1,TH1D* dr, bool withbg = 1){
		JTCSignalProducer sp;
		auto herr = new TH1D(name, "", ndr, drbin); herr->Sumw2();
		int l1p5 = h1->FindBin(-1.4999)-1;
		int l2   = h1->FindBin(-2.0001);
		int r2   = h1->FindBin( 2.0);
		int r1p5 = h1->FindBin( 1.5);
		float mean = 0;
		float left_ave  = (h1->GetBinContent(l1p5)+h1->GetBinContent(l2))/2;
		float right_ave = (h1->GetBinContent(r1p5)+h1->GetBinContent(r2))/2;
		float in_ave = (h1->GetBinContent(l1p5)+h1->GetBinContent(r1p5))/2;
		float out_ave= (h1->GetBinContent(l2)+h1->GetBinContent(r2))/2;
		float me_err = max(fabs(left_ave-mean), fabs(right_ave-mean));
		float bg_err = max(fabs(in_ave-mean), fabs(out_ave-mean));
		for(int k=1; k< ndr+1; ++k){
				float ring =TMath::Pi()*(pow(herr->GetBinLowEdge(k)+herr->GetBinWidth(k),2)-
								pow(herr->GetBinLowEdge(k),2))/2/herr->GetBinWidth(k);
				float err;
				if(withbg) err = ring*pow(bg_err*bg_err+me_err*me_err, 0.5);
				else err	= ring*me_err;
				herr->SetBinError(k, err);
				herr->SetBinContent(k, dr->GetBinContent(k));
		}
		herr->GetXaxis()->SetTitle(dr->GetXaxis()->GetTitle());
		return herr;
}

TH1D** getDrErr(TString name,TH1D** dr, TH1D** h1, bool withbg = 1){
		TH1D** drerr = new TH1D*[nPt*nCent];
		for(int i=0; i<nPt; ++i){
				for(int j=0; j<nCent; ++j){
						drerr[i+nPt*j]=get_bkg_dr_err(name+Form("_%d_%d",i, j), h1[i+nPt*j], dr[i+nPt*j], withbg);
				}
		}
		return drerr;
}
