

void minJetTree(){
		TFile * f = TFile::Open("/Users/tabris/cmsProjects/iJTC/dataSet/Hydjet_cymbalTune_5TeV_1.root");
		TTree *t;
		t = (TTree*) f->Get("mixing_tree");
		t->SetBranchStatus("*",0);

		Float_t vz, pthat;
		Int_t hiBin, HBHENoiseFilterResultRun2Loose, pcollisionEventSelection, pprimaryVertexFilter, phfCoincFilter3;
		vector<float> *jteta=0, *jtphi=0, *jtpt=0, *corrpt=0, *rawpt=0;
		vector<float> *geneta=0, *genphi=0, *genpt=0;
		vector<float> *trackMax=0;
		vector<int > *flavorForB=0;
		vector<float> *discr_csvV1=0;

		Float_t ovz, opthat;
		Int_t ohiBin;
		vector<float> *ojteta=0, *ojtphi=0, *ojtpt=0, *ocorrpt=0, *orawpt=0;
		vector<float> *ogeneta=0, *ogenphi=0, *ogenpt=0;
		vector<float> *otrackMax=0;
		vector<int > *oflavorForB=0;
		vector<float> *odiscr_csvV1=0;

		t->SetBranchStatus("HBHENoiseFilterResultRun2Loose", 1);
		t->SetBranchStatus("pprimaryVertexFilter", 1);
		t->SetBranchStatus("pcollisionEventSelection", 1);
		t->SetBranchStatus("phfCoincFilter3", 1);
		t->SetBranchAddress("HBHENoiseFilterResultRun2Loose", &HBHENoiseFilterResultRun2Loose);
		t->SetBranchAddress("pprimaryVertexFilter", &pprimaryVertexFilter);
		t->SetBranchAddress("pcollisionEventSelection", &pcollisionEventSelection);
		t->SetBranchAddress("phfCoincFilter3", &phfCoincFilter3);

		t->SetBranchStatus("vz", 1);
		t->SetBranchStatus("hiBin",1);
		t->SetBranchAddress("hiBin", &hiBin);
		t->SetBranchAddress("vz", &vz);

		t->SetBranchStatus("calo_jteta" , 1);
		t->SetBranchStatus("calo_jtphi" , 1);
		t->SetBranchStatus("calo_jtpt"  , 1);
		t->SetBranchStatus("calo_rawpt" , 1);
		t->SetBranchStatus("calo_corrpt", 1);
		t->SetBranchAddress("calo_jteta", &jteta);
		t->SetBranchAddress("calo_jtphi", &jtphi);
		t->SetBranchAddress("calo_jtpt", &jtpt);
		t->SetBranchAddress("calo_rawpt", &rawpt);
		t->SetBranchAddress("calo_corrpt", &corrpt);
		t->SetBranchStatus("calo_trackMax", 1);
		t->SetBranchAddress("calo_trackMax", &trackMax);

		t->SetBranchStatus("genpt" , 1);
		t->SetBranchStatus("geneta" , 1);
		t->SetBranchStatus("genphi"  , 1);
		t->SetBranchAddress("genpt" , &genpt);
		t->SetBranchAddress("geneta", &geneta);
		t->SetBranchAddress("genphi", &genphi);

		t->SetBranchStatus("calo_refparton_flavorForB", 1);
		t->SetBranchAddress("calo_refparton_flavorForB", &flavorForB);

		t->SetBranchStatus("calo_discr_csvV1", 1);
		t->SetBranchAddress("calo_discr_csvV1",&discr_csvV1);

		TFile *wf = TFile::Open("output.root", "recreate");
		TTree * outTree = new TTree("miniTree", "");
		outTree->Branch("vz", &ovz);
		outTree->Branch("hiBin", &ohiBin);
		outTree->Branch("calo_jteta", &ojteta);
		outTree->Branch("calo_jtphi", &ojtphi);
		outTree->Branch("calo_jtpt", &ojtpt);
		outTree->Branch("calo_rawpt", &orawpt);
		outTree->Branch("calo_corrpt", &ocorrpt);

		outTree->Branch("genpt", &ogenpt);
		outTree->Branch("geneta", &ogeneta);
		outTree->Branch("genphi", &ogenphi);
		outTree->Branch("calo_trackMax", &otrackMax);
		outTree->Branch("calo_refparton_flavorForB", &oflavorForB);
		outTree->Branch("calo_discr_csvV1",&odiscr_csvV1);

		Long64_t nentries = t->GetEntriesFast();
		for(Long64_t jentry = 0; jentry<nentries; ++jentry){
				if(jentry %1000==0) std::cout<<"processing event "<<jentry<<std::endl;
				t->GetEntry(jentry);
				if ( HBHENoiseFilterResultRun2Loose ==0) continue;
				if ( pcollisionEventSelection ==0) continue;
				if ( pprimaryVertexFilter ==0) continue;
				if ( phfCoincFilter3 ==0) continue;
				if ( vz>15 || vz<-15) continue;
				ohiBin = hiBin;
				ovz = vz;
				opthat = pthat;
				for( int i=0; i<int(jtpt->size()); ++i){
						ojtpt->push_back    (jtpt->at(i));	
						ojteta->push_back   (jteta->at(i));	
						ojtphi->push_back   (jtphi->at(i));	
						orawpt->push_back   (rawpt->at(i));	
						otrackMax->push_back(trackMax->at(i));
						odiscr_csvV1->push_back(discr_csvV1->at(i));
				}
				for( int i=0; i<int(genpt->size()); ++i){
						ogenpt->push_back (genpt->at(i));
						ogeneta->push_back(geneta->at(i));
						ogenphi->push_back(genphi->at(i));
				}
				outTree->Fill();
				ojtpt    ->clear();
				ojteta   ->clear();
				ojtphi   ->clear();
				orawpt   ->clear();
				otrackMax->clear();
				ogenpt ->clear();
				ogeneta->clear();
				ogenphi->clear();
				odiscr_csvV1->clear();
		}
		outTree->Write();
}
