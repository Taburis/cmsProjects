
#include "../../lib/config.h"
#include "../../lib/trackingCorr2.h"
using namespace jetTrack;
void corrMaker(){
	loadConfig();
	TFile *f = TFile::Open("../../dataSet/corrScan_pp_5TeV.root");
	auto a= new trackingCorr2(f);
	a->is_pp=1;
	a->Read();
	//do symmetrization;
//	a->loadDataTracks("../../dataSet/tracking/data_clsoure.root");
	a->getCorr("corrTable_pp_noDCA.root");
	//a->getCorr("cymbalCorr_aux.root");
	a->showCorr(0,1);
}
