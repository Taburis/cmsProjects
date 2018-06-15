

void pullSig(inputSet &iset){
		readInput(iset);
		TString prefix = tpname(iset);
		TString cap = histname(iset);
		bool isppincl = !(iset.isHi);
		pullSignal("taggedBJet_"+prefix+"_JTCSignal.root", "taggedBJet_"+cap, iset.tg, 1.5, 2.5, 0);
		pullSignal("inclJet_"+prefix+"_JTCSignal.root", "inclJet_"+cap, iset.tg, 1.5, 2.5, isppincl);
		auto sig_tg = read_flatten<TH2D>(dataDumpPath+"taggedBJet_"+prefix+"_JTCSignal.root", "signal_taggedBJet_"+cap+"_noCorr");
		auto sig_in = read_flatten<TH2D>(dataDumpPath+"inclJet_"+prefix+"_JTCSignal.root", "signal_inclJet_"+cap+"_noCorr");
		auto sig_co = read_flatten<TH2D>(dataDumpPath+"contJet_"+prefix+"_JTCSignal.root", "signal_contJet_"+cap+"_noCorr");
		checkSideBand("taggedB_"+prefix, sig_tg);
		checkSideBand("incl_"+prefix, sig_in);
		if(iset.ismc){
				pullSignal("taggedTrueBJet_"+prefix+"_JTCSignal.root", "taggedTrueBJet_"+cap, iset.tt, 1.5, 2.5, 0);
				pullSignal("trueBJet_"+prefix+"_JTCSignal.root", "trueBJet_"+cap, iset.tb, 1.5, 2.5, 0);
				pullSignal("contJet_"+prefix+"_JTCSignal.root", "contJet_"+cap, iset.tb, 1.5, 2.5, 0);
				auto sig_co = read_flatten<TH2D>(dataDumpPath+"contJet_"+prefix+"_JTCSignal.root", "signal_contJet_"+cap+"_noCorr");
				auto sig_tt = read_flatten<TH2D>(dataDumpPath+"taggedTrueBJet_"+prefix+"_JTCSignal.root", "signal_taggedTrueBJet_"+cap+"_noCorr");
				auto sig_tb = read_flatten<TH2D>(dataDumpPath+"trueBJet_"+prefix+"_JTCSignal.root", "signal_trueBJet_"+cap+"_noCorr");
				checkSideBand("trueB_"+prefix, sig_tb);
				checkSideBand("taggedTrueB_"+prefix, sig_tt);
		}
}
