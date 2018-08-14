

void pullSig(inputSet &iset, TString label = "", bool doCont = 0){
		readInput(iset);
		TString prefix = tpname(iset);
		prefix=prefix+label;
		TString cap = histname(iset);
		TString header = iset.isHi ? "HI_" : "";
		bool isppincl = !(iset.isHi);
		//pullSignal(header+"taggedBJet_"+prefix+"_JTCSignal.root", "taggedBJet_"+cap, iset.tg, 1.5, 2.5);
		//pullSignal(header+"inclJet_"+prefix+"_JTCSignal.root", "inclJet_"+cap, iset.in, 1.5, 2.5);
		pullSignal(header+"taggedBJet_"+prefix+"_JTCSignal.root", "taggedBJet_"+cap, iset.tg, 1.5, 2.5, isppincl);
		pullSignal(header+"inclJet_"+prefix+"_JTCSignal.root", "inclJet_"+cap, iset.in, 1.5, 2.5, isppincl);
		auto sig_tg = read_flatten<TH2D>(dataDumpPath+header+"taggedBJet_"+prefix+"_JTCSignal.root", "signal_taggedBJet_"+cap+"_noCorr");
		auto sig_in = read_flatten<TH2D>(dataDumpPath+header+"inclJet_"+prefix+"_JTCSignal.root", "signal_inclJet_"+cap+"_noCorr");
	//	auto sig_co = read_flatten<TH2D>(dataDumpPath+header+"contJet_"+prefix+"_JTCSignal.root", "signal_contJet_"+cap+"_noCorr");
		checkSideBand(header+"taggedB_"+prefix, sig_tg);
		checkSideBand(header+"incl_"+prefix, sig_in);
		if(iset.ismc){
				pullSignal(header+"taggedTrueBJet_"+prefix+"_JTCSignal.root", "taggedTrueBJet_"+cap, iset.tt, 1.5, 2.5, 0);
				pullSignal(header+"trueBJet_"+prefix+"_JTCSignal.root", "trueBJet_"+cap, iset.tb, 1.5, 2.5, 0);
				if(doCont&&iset.hasCont){
					   	pullSignal(header+"contJet_"+prefix+"_JTCSignal.root", "contJet_"+cap, iset.co, 1.5, 2.5, 0);
						addDr("contJet_"+tpname(iset), "contJet_"+cap, dataDumpPath+header+"contJet_"+prefix+"_JTCSignal.root");
				}
	//			auto sig_co = read_flatten<TH2D>(dataDumpPath+header+"contJet_"+prefix+"_JTCSignal.root", "signal_contJet_"+cap+"_noCorr");
				auto sig_tt = read_flatten<TH2D>(dataDumpPath+header+"taggedTrueBJet_"+prefix+"_JTCSignal.root", "signal_taggedTrueBJet_"+cap+"_noCorr");
				auto sig_tb = read_flatten<TH2D>(dataDumpPath+header+"trueBJet_"+prefix+"_JTCSignal.root", "signal_trueBJet_"+cap+"_noCorr");
				checkSideBand(header+"trueB_"+prefix, sig_tb);
				checkSideBand(header+"taggedTrueB_"+prefix, sig_tt);
				addDr("trueBJet_"+tpname(iset), "trueBJet_"+cap, dataDumpPath+header+"trueBJet_"+prefix+"_JTCSignal.root");
				addDr("taggedTrueBJet_"+tpname(iset), "taggedTrueBJet_"+cap, dataDumpPath+header+"taggedTrueBJet_"+prefix+"_JTCSignal.root");
		}
		addDr("taggedBJet_"+tpname(iset), "taggedBJet_"+cap, dataDumpPath+header+"taggedBJet_"+prefix+"_JTCSignal.root");
		addDr("inclJet_"+tpname(iset), "inclJet_"+cap, dataDumpPath+header+"inclJet_"+prefix+"_JTCSignal.root");
}
