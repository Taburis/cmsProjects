
#include <iostream>
#include "edmJtcProducer.h"
//#include "edmMCBase.h"
#include "ppjtc_cfg.h"
#include "TTree.h"

class mcCompare : public edmJtcProducer{
		public : mcCompare(ParaSet &s): edmJtcProducer(s){}
				 virtual void handleGenParticle(const char * name){
						 gen_particle_t = handle(name);
						 gen_particle_t->SetBranchAddress("pt" , &pt);
						 gen_particle_t->SetBranchAddress("eta", &eta);
						 gen_particle_t->SetBranchAddress("phi", &phi);
						 gen_particle_t->SetBranchAddress("chg", &chg);
						 gen_particle_t->SetBranchAddress("pdg", &pdg);
				 }
				 virtual void beginJob(){
						 std::cout<<"starting job"<<std::endl;
						 handleGenParticle("genParticleTree");
						 auto jes4 = handleJetSet("ak4GenJetTree_ESchemeWTA", 1);
						 regHist("ak4Escheme", hc, jes4);
						 auto t = jes4->jtTree;
						 t->SetBranchAddress("pthat", &pthat);
						 t->SetBranchAddress("weight", &weight);
						 t->SetBranchAddress("genJtPt",  &jes4->jt_pt);
						 t->SetBranchAddress("genJtPhi", &jes4->jt_phi);
						 t->SetBranchAddress("genJtEta", &jes4->jt_eta);
				 }
};

class myMC : public edmJtcProducer{
		public : myMC(ParaSet &s) : edmJtcProducer(s){}
				 virtual void beginJob(){
						 handleGenParticle("particle_tree");
						 auto jw3 = handleJetSet("ak3WTA");
						 regHist("ak3WTA", hc3, jw3);
						 auto jw4 = handleJetSet("ak4WTA");
						 regHist("ak4WTA", hc4, jw4);
						 auto jw5 = handleJetSet("ak5WTA");
						 regHist("ak5WTA", hc5, jw5);
				 }
				 virtual void analyzer(){
						 jtc_wf001(hc3);
						 jtc_wf001(hc4);
						 jtc_wf001(hc5);
				 }
		public : 
				 histCase hc3, hc4, hc5;
};
/*
class wtaAnalyzer : public edmJtcProducer{
		public : wtaAnalyzer(ParaSet &s): edmJtcProducer(s) { 
						 hm = new histManager();
				 };
				 virtual void beginJob(){
						 handleGenParticle("particle_tree");
						 hcwta.set = handleJetSet("ak4WTA", 1);
						 hcesm.set = handleJetSet("ak4ES", 1);
						 hcmix.set = handleJetSet("ak4ES_WTA");
						 regHist("ak4WTA", hcwta);
						 regHist("ak4ESM", hcesm);
						 regHist("ak4MIX", hcmix);
						 hm->sumw2();
				 }

				 bool eventCut(){
						 // return 1 indicate to pass this event 
						 if(pthat < 80) return 1;
						 bool hasESMJet = 0 ;
						 bool hasWTAJet = 0 ;
						 for(int i=0 ; i<hcwta.set->jt_pt->size(); ++i){
								 if(hcwta.set->jt_pt->at(i) > 80) {
										 if(fabs(hcwta.set->jt_eta->at(i))<2){
												 hasWTAJet = 1;
												 break;
										 }
								 }
						 }
						 for(int i=0 ; i<hcesm.set->jt_pt->size(); ++i){
								 if(hcesm.set->jt_pt->at(i) > 80) {
										 if(fabs(hcesm.set->jt_eta->at(i))<2){
												 hasESMJet = 1;
												 break;
										 }
								 }
						 }
						 bool res = (!hasESMJet) || (!hasWTAJet);
						 return res;
				 }
				 unsigned int loadLeading(jetSet *js, int j){
						 float pt = 0; 
						 unsigned int index = 0;
						 unsigned int upper = j < js->jt_pt->size()-1 ? js->pfcand_indx->at(j+1): js->pfcand_pt->size(); 
						 if( upper > js->pfcand_pt->size()){
								 std::cout<<"index is wrong: "<<js->pfcand_pt->size()<<" : "<<upper<<std::endl;
						 }
						 for(int k = js->pfcand_indx->at(j); k<upper; ++k){
								 if( pt > js->pfcand_pt->at(k)) continue;
								 pt = js->pfcand_pt->at(k);
								 index = k;
						 }

						 return index;
				 }
				 bool jetAcceptanceCut(histCase &hc){
						 if(hc.jtpt<80) return 1;
						 if(fabs(hc.jteta)>2) return 1;
						 return 0;
				 }

				 void analyzer(){
						 if (eventCut()) return ;

						 for(unsigned int j=0; j<hcesm.set->jt_pt->size(); ++j){
								 loadJet(hcesm, j);
							 if(jetAcceptanceCut(hcesm)) continue;
						 }
				 };
				 virtual void endJob(){
						 auto f = TFile::Open("test.root","recreate");
						 f->cd();
						 hm->write();
						 f->Close();
				 }

		public : histManager* hm;
				 TH2D* h2D;
				 histCase hcwta, hcesm, hcmix;
				 histCase hcwta1, hcesm1; //for leading pfcandidate
};
*/

int main(){
		auto ps = makePSet_edmJtcDefault();
		rootEDMFrame * ff = new rootEDMFrame();

		//ps.setPara<TString>("output_name", "output_my.root");
		ps.setPara<TString>("output_name", "output_wta.root");
		ff->open("/Users/tabris/cmsProjects/myGenerator/test_CUETP8M1_NEvt100000.root");
//		auto a = new edmJtcProducer(ps);
		auto a = new myMC(ps);

		//ps.setPara<TString>("output_name", "output_chris.root");
		//ff->open("/Applications/physPlayGround/Generators/PYTHIA8/bin/output/20181211/forJetAxis_CUETP8M1_NEvt100000_DoFlatPthatTrue_FlatPower4p50_PthatMin100p0_PthatMax999999p0_20181211.root");
		//auto a = new mcCompare(ps);

		ff->EventMax=-1;
		ff->addAnalyzer((rootEDMAnalyzer*) a);
		ff->eventLoop();
		delete ff;
		return 0;
}
