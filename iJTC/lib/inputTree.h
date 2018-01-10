//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Fri Apr 21 18:19:32 2017 by ROOT version 6.08/05
// from TTree mixing_tree/
// found on file: root://cmsxrootd.fnal.gov//store/user/kjung/PbPb_Pythia6Hydjet_MC_JetTrackSkim_looseTrkCutsForB/PbPb_Pythia6Hydjet_5TeVSkim_Gendxyz_March16/crab_PbPb_Pythia6Hydjet_MC_JetTrackSkim_looseTrkCutsForB/170316_164800/0000/HydJet15_9.root
//////////////////////////////////////////////////////////

#ifndef inputTree_h
#define inputTree_h

//#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include "vector"
using namespace std;
class inputTree {
		public :
				TTree          *fChain;   //!pointer to the analyzed TTree or TChain
				Int_t           fCurrent; //!current Tree number in a TChain

				// Fixed size dimensions of array or collections stored in the TTree if any.

				bool ispp;

				Int_t           HLT_ak4CaloJet80;
				Int_t           HLT_ak4CaloJet100;
				Int_t           HLT_ak4CaloJet80_Prescale;
				Int_t           HLT_ak4CaloJet100_Prescale;
				Int_t           HBHENoiseFilterResult;
				Int_t           HBHENoiseFilterResultRun1;
				Int_t           HBHENoiseFilterResultRun2Loose;
				Int_t           HBHENoiseFilterResultRun2Tight;
				Int_t           HBHEIsoNoiseFilterResult;

				//Pb 
				Int_t           pcollisionEventSelection;
				Int_t           pclusterCompatibilityFilter;
				Int_t           pprimaryVertexFilter;
				Int_t           phfCoincFilter1;
				Int_t           phfCoincFilter2;
				Int_t           phfCoincFilter3;
				Int_t           phfCoincFilter4;
				Int_t           phfCoincFilter5;

				Int_t           pPAprimaryVertexFilter;// pp
				// Declaration of leaf types

				vector<float>   *trkPt;
				vector<float>   *trkEta;
				vector<float>   *trkPhi;
				vector<int>     *trkAlgo;
				vector<bool>    *highPurity;
				Float_t         vz;
				Int_t           hiBin;
				//   Float_t         evtPlane_HF2;
				Float_t         pthat;
				vector<float>   *trkDxy;
				vector<float>   *trkDxyError;
				vector<float>   *trkDz;
				vector<float>   *trkDzError;
				vector<float>   *trkPtError;
				vector<float>   *trkChi2;
				vector<int>     *trkNdof;
				vector<int>     *trkNHit;
				vector<int>     *trkNlayer;
				vector<float>   *pfEcal;
				vector<float>   *pfHcal;
				vector<float>   *trkMVALoose;
				vector<float>   *trkMVATight;
				vector<float>   *pt;
				vector<float>   *phi;
				vector<float>   *eta;
				vector<int>     *chg;
				vector<int>     *sube;
				vector<int>     *pdg;
//				vector<float>   *dxy;
//				vector<float>   *dz;

				// List of branches
				TBranch        *b_HLT_ak4CaloJet80;   //!
				TBranch        *b_HLT_ak4CaloJet100;   //!
				TBranch        *b_HLT_ak4CaloJet80_Prescale;   //!
				TBranch        *b_HLT_ak4CaloJet100_Prescale;   //!
				TBranch        *b_HBHENoiseFilterResult;   //!
				TBranch        *b_HBHENoiseFilterResultRun1;   //!
				TBranch        *b_HBHENoiseFilterResultRun2Loose;   //!
				TBranch        *b_HBHENoiseFilterResultRun2Tight;   //!
				TBranch        *b_HBHEIsoNoiseFilterResult;   //!
				TBranch        *b_pcollisionEventSelection;   //!
				TBranch        *b_pprimaryVertexFilter;   //!
				TBranch        *b_pPAprimaryVertexFilter;   //!
				TBranch        *b_pclusterCompatibilityFilter;   //!
				TBranch        *b_phfCoincFilter1;   //!
				TBranch        *b_phfCoincFilter2;   //!
				TBranch        *b_phfCoincFilter3;   //!
				TBranch        *b_phfCoincFilter4;   //!
				TBranch        *b_phfCoincFilter5;   //!
				TBranch        *b_trkPt;   //!
				TBranch        *b_trkEta;   //!
				TBranch        *b_trkPhi;   //!
				TBranch        *b_trkAlgo;   //!
				TBranch        *b_highPurity;   //!
				TBranch        *b_vz;   //!
				TBranch        *b_hiBin;   //!
				//   TBranch        *b_evtPlane_HF2;   //!
				TBranch        *b_pthat;   //!
				TBranch        *b_trkDxy;   //!
				TBranch        *b_trkDxyError;   //!
				TBranch        *b_trkDz;   //!
				TBranch        *b_trkDzError;   //!
				TBranch        *b_trkPtError;   //!
				TBranch        *b_trkChi2;   //!
				TBranch        *b_trkNdof;   //!
				TBranch        *b_trkNHit;   //!
				TBranch        *b_trkNlayer;   //!
				TBranch        *b_pfEcal;   //!
				TBranch        *b_pfHcal;   //!
				TBranch        *b_trkMVALoose;   //!
				TBranch        *b_trkMVATight;   //!
				TBranch        *b_pt;   //!
				TBranch        *b_phi;   //!
				TBranch        *b_eta;   //!
				TBranch        *b_chg;   //!
				TBranch        *b_sube;   //!
				TBranch        *b_pdg;   //!
//				TBranch        *b_dxy;   //!
//				TBranch        *b_dz;   //!

				inputTree(TTree *tree=0);
				virtual ~inputTree();
				virtual Int_t    Cut(Long64_t entry);
				virtual Int_t    GetEntry(Long64_t entry);
				virtual Long64_t LoadTree(Long64_t entry);
				virtual void     Init(TTree *tree);
				virtual Bool_t   Notify();
				virtual void     Show(Long64_t entry = -1);
};


inputTree::inputTree(TTree *tree) : fChain(0) 
{
		// if parameter tree is not specified (or zero), connect the file
		// used to generate this class and read the Tree.
	//	Init(tree);
}

inputTree::~inputTree()
{
		if (!fChain) return;
		delete fChain->GetCurrentFile();
}

Int_t inputTree::GetEntry(Long64_t entry)
{
		// Read contents of entry.
		if (!fChain) return 0;
		return fChain->GetEntry(entry);
}
Long64_t inputTree::LoadTree(Long64_t entry)
{
		// Set the environment to read one entry
		if (!fChain) return -5;
		Long64_t centry = fChain->LoadTree(entry);
		if (centry < 0) return centry;
		if (fChain->GetTreeNumber() != fCurrent) {
				fCurrent = fChain->GetTreeNumber();
				Notify();
		}
		return centry;
}

void inputTree::Init(TTree *tree)
{
		// The Init() function is called when the selector needs to initialize
		// a new tree or chain. Typically here the branch addresses and branch
		// pointers of the tree will be set.
		// It is normally not necessary to make changes to the generated
		// code, but the routine can be extended by the user if needed.
		// Init() will be called many times when running on PROOF
		// (once per file to be processed).

		// Set object pointer
		trkPt = 0;
		trkEta = 0;
		trkPhi = 0;
		trkAlgo = 0;
		highPurity = 0;
		trkDxy = 0;
		trkDxyError = 0;
		trkDz = 0;
		trkDzError = 0;
		trkPtError = 0;
		trkChi2 = 0;
		trkNdof = 0;
		trkNHit = 0;
		trkNlayer = 0;
		pfEcal = 0;
		pfHcal = 0;
		trkMVALoose = 0;
		trkMVATight = 0;
		pt = 0;
		phi = 0;
		eta = 0;
		chg = 0;
		sube = 0;
		pdg = 0;
//		dxy = 0;
//		dz = 0;
		// Set branch addresses and branch pointers
		if (!tree) return;
		fChain = tree;
		fCurrent = -1;
		fChain->SetMakeClass(1);

		fChain->SetBranchAddress("HLT_ak4CaloJet80", &HLT_ak4CaloJet80, &b_HLT_ak4CaloJet80);
		fChain->SetBranchAddress("HLT_ak4CaloJet100", &HLT_ak4CaloJet100, &b_HLT_ak4CaloJet100);
		fChain->SetBranchAddress("HLT_ak4CaloJet80_Prescale", &HLT_ak4CaloJet80_Prescale, &b_HLT_ak4CaloJet80_Prescale);
		fChain->SetBranchAddress("HLT_ak4CaloJet100_Prescale", &HLT_ak4CaloJet100_Prescale, &b_HLT_ak4CaloJet100_Prescale);
		fChain->SetBranchAddress("HBHENoiseFilterResult", &HBHENoiseFilterResult, &b_HBHENoiseFilterResult);
		fChain->SetBranchAddress("HBHENoiseFilterResultRun1", &HBHENoiseFilterResultRun1, &b_HBHENoiseFilterResultRun1);
		fChain->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose, &b_HBHENoiseFilterResultRun2Loose);
		fChain->SetBranchAddress("HBHENoiseFilterResultRun2Tight", &HBHENoiseFilterResultRun2Tight, &b_HBHENoiseFilterResultRun2Tight);
		fChain->SetBranchAddress("HBHEIsoNoiseFilterResult", &HBHEIsoNoiseFilterResult, &b_HBHEIsoNoiseFilterResult);

		if( ispp) fChain->SetBranchAddress("pPAprimaryVertexFilter", &pPAprimaryVertexFilter, &b_pPAprimaryVertexFilter);
		else {
				fChain->SetBranchAddress("hiBin", &hiBin, &b_hiBin);
				fChain->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection, &b_pcollisionEventSelection);
				fChain->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter, &b_pprimaryVertexFilter);
				fChain->SetBranchAddress("pclusterCompatibilityFilter", &pclusterCompatibilityFilter, &b_pclusterCompatibilityFilter);
				fChain->SetBranchAddress("phfCoincFilter1", &phfCoincFilter1, &b_phfCoincFilter1);
				fChain->SetBranchAddress("phfCoincFilter2", &phfCoincFilter2, &b_phfCoincFilter2);
				fChain->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3, &b_phfCoincFilter3);
				fChain->SetBranchAddress("phfCoincFilter4", &phfCoincFilter4, &b_phfCoincFilter4);
				fChain->SetBranchAddress("phfCoincFilter5", &phfCoincFilter5, &b_phfCoincFilter5);
		}

		fChain->SetBranchAddress("trkPt", &trkPt, &b_trkPt);
		fChain->SetBranchAddress("trkEta", &trkEta, &b_trkEta);
		fChain->SetBranchAddress("trkPhi", &trkPhi, &b_trkPhi);
		fChain->SetBranchAddress("trkAlgo", &trkAlgo, &b_trkAlgo);
		fChain->SetBranchAddress("highPurity", &highPurity, &b_highPurity);
		fChain->SetBranchAddress("vz", &vz, &b_vz);
		//   fChain->SetBranchAddress("evtPlane_HF2", &evtPlane_HF2, &b_evtPlane_HF2);
		fChain->SetBranchAddress("pthat", &pthat, &b_pthat);
		fChain->SetBranchAddress("trkDxy", &trkDxy, &b_trkDxy);
		fChain->SetBranchAddress("trkDxyError", &trkDxyError, &b_trkDxyError);
		fChain->SetBranchAddress("trkDz", &trkDz, &b_trkDz);
		fChain->SetBranchAddress("trkDzError", &trkDzError, &b_trkDzError);
		fChain->SetBranchAddress("trkPtError", &trkPtError, &b_trkPtError);
		fChain->SetBranchAddress("trkChi2", &trkChi2, &b_trkChi2);
		fChain->SetBranchAddress("trkNdof", &trkNdof, &b_trkNdof);
		fChain->SetBranchAddress("trkNHit", &trkNHit, &b_trkNHit);
		fChain->SetBranchAddress("trkNlayer", &trkNlayer, &b_trkNlayer);
		fChain->SetBranchAddress("pfEcal", &pfEcal, &b_pfEcal);
		fChain->SetBranchAddress("pfHcal", &pfHcal, &b_pfHcal);
		fChain->SetBranchAddress("trkMVALoose", &trkMVALoose, &b_trkMVALoose);
		fChain->SetBranchAddress("trkMVATight", &trkMVATight, &b_trkMVATight);
		fChain->SetBranchAddress("pt", &pt, &b_pt);
		fChain->SetBranchAddress("phi", &phi, &b_phi);
		fChain->SetBranchAddress("eta", &eta, &b_eta);
		fChain->SetBranchAddress("chg", &chg, &b_chg);
		fChain->SetBranchAddress("sube", &sube, &b_sube);
		fChain->SetBranchAddress("pdg", &pdg, &b_pdg);
//		fChain->SetBranchAddress("dxy", &dxy, &b_dxy);
//		fChain->SetBranchAddress("dz", &dz, &b_dz);
		Notify();
}

Bool_t inputTree::Notify()
{
		// The Notify() function is called when a new file is opened. This
		// can be either for a new TTree in a TChain or when when a new TTree
		// is started when using PROOF. It is normally not necessary to make changes
		// to the generated code, but the routine can be extended by the
		// user if needed. The return value is currently not used.

		return kTRUE;
}

void inputTree::Show(Long64_t entry)
{
		// Print contents of entry.
		// If entry is not specified, print current entry
		if (!fChain) return;
		fChain->Show(entry);
}
Int_t inputTree::Cut(Long64_t entry)
{
		// This function may be called from Loop.
		// returns  1 if entry is accepted.
		// returns -1 otherwise.
		return 1;
}
#endif // #ifdef inputTree_cxx
