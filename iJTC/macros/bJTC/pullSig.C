

void pullSig(inputSet &iset){
		readInput(iset);
		TString prefix = tpname(iset);
		TString cap = histname(iset);
		TString header = iset.isHi ? "HI_" : "";
		bool isppincl = !(iset.isHi);
		pullSignal(header+"taggedBJet_"+prefix+"_JTCSignal.root", "taggedBJet_"+cap, iset.tg, 1.5, 2.5, isppincl);
		pullSignal(header+"inclJet_"+prefix+"_JTCSignal.root", "inclJet_"+cap, iset.in, 1.5, 2.5, isppincl);
		auto sig_tg = read_flatten<TH2D>(dataDumpPath+header+"taggedBJet_"+prefix+"_JTCSignal.root", "signal_taggedBJet_"+cap+"_noCorr");
		auto sig_in = read_flatten<TH2D>(dataDumpPath+header+"inclJet_"+prefix+"_JTCSignal.root", "signal_inclJet_"+cap+"_noCorr");
		auto sig_co = read_flatten<TH2D>(dataDumpPath+header+"contJet_"+prefix+"_JTCSignal.root", "signal_contJet_"+cap+"_noCorr");
		checkSideBand(header+"taggedB_"+prefix, sig_tg);
		checkSideBand(header+"incl_"+prefix, sig_in);
		if(iset.ismc){
				pullSignal(header+"taggedTrueBJet_"+prefix+"_JTCSignal.root", "taggedTrueBJet_"+cap, iset.tt, 1.5, 2.5, 0);
				pullSignal(header+"trueBJet_"+prefix+"_JTCSignal.root", "trueBJet_"+cap, iset.tb, 1.5, 2.5, 0);
				pullSignal(header+"contJet_"+prefix+"_JTCSignal.root", "contJet_"+cap, iset.tb, 1.5, 2.5, 0);
				auto sig_co = read_flatten<TH2D>(dataDumpPath+header+"contJet_"+prefix+"_JTCSignal.root", "signal_contJet_"+cap+"_noCorr");
				auto sig_tt = read_flatten<TH2D>(dataDumpPath+header+"taggedTrueBJet_"+prefix+"_JTCSignal.root", "signal_taggedTrueBJet_"+cap+"_noCorr");
				auto sig_tb = read_flatten<TH2D>(dataDumpPath+header+"trueBJet_"+prefix+"_JTCSignal.root", "signal_trueBJet_"+cap+"_noCorr");
				checkSideBand(header+"trueB_"+prefix, sig_tb);
				checkSideBand(header+"taggedTrueB_"+prefix, sig_tt);
		}
}
