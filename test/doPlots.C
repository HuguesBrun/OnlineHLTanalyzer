TFile *theOutFile = new TFile("theHistosFile_MC.root","RECREATE");
TChain *chain = new TChain("eventsTree");


void preparePlots(TString name, TString variable, TString cuts, int nBins, float lowBin, float highBin){
    cout << "creating " << name << endl;
    TH1F *histo = new TH1F("histo","",nBins,lowBin,highBin);
    chain->Draw(variable+">>histo",cuts);
    
    theOutFile->cd();
    histo->Write(name);
    
    delete histo;
}




void doPlots(){
    chain->Add("/tmp/hbrun/MCminBias/muonL2tree_*.root");
    
    
    preparePlots("ptL2","T_L2muon_Pt","1",100,0,200);
    preparePlots("etaL2","T_L2muon_Eta","1",20,-2.5,2.5);
    preparePlots("etaL2highBins","T_L2muon_Eta","1",100,-2.5,2.5);
    preparePlots("phiL2","T_L2muon_Phi","1",10,-3.1416,3.1416);
    preparePlots("phiL2highBin","T_L2muon_Phi","1",100,-3.1416,3.1416);
    preparePlots("dxyL2","T_L2muon_dxy","1",10,-200,200);
    preparePlots("dxyL2highBin","T_L2muon_dxy","1",100,-200,200);
    preparePlots("dzL2","T_L2muon_dz","1",10,-400,400);
    preparePlots("dzL2highBin","T_L2muon_dz","1",100,-400,400);
    preparePlots("nbStationsWithHits","T_L2muon_nbStationWithHits","1",10,0,10);
    preparePlots("nbStationsWithHitsCSC","T_L2muon_nbStationWithHitsCSC","1",10,0,10);
    preparePlots("nbStationsWithHitsDT","T_L2muon_nbStationWithHitsDT","1",10,0,10);
    preparePlots("nbValidHits","T_L2muon_nbValidHits","1",12,0,24);
    
    preparePlots("ptL3","T_L3muon_Pt","1",100,0,200);
    preparePlots("etaL3","T_L3muon_Eta","1",20,-2.5,2.5);
    preparePlots("etaL3highBins","T_L3muon_Eta","1",100,-2.5,2.5);
    preparePlots("phiL3","T_L3muon_Phi","1",10,-3.1416,3.1416);
    preparePlots("phiL3highBins","T_L3muon_Phi","1",100,-3.1416,3.1416);
    preparePlots("drL3","T_L3muon_dr","1",10,0,0.5);
    preparePlots("drL3highBins","T_L3muon_dr","1",100,0,0.5);
    preparePlots("dzL3","T_L3muon_dz","1",20,0,50);
    preparePlots("dzL3highBins","T_L3muon_dz","1",100,0,50);
    preparePlots("dxyL3","T_L3muon_dxyBS","1",20,0,0.5);
    preparePlots("dxyL3highBins","T_L3muon_dxyBS","1",100,0,0.5);
    preparePlots("Chi2L3","T_L3muon_Chi2","1",20,0,100);
    preparePlots("Chi2L3highBins","T_L3muon_Chi2","1",100,0,100);

    TH2F *L2etaPhi = new TH2F("L2etaPhi","",100,-3.1416,3.1416,100,-2.5,2.5);
    chain->Draw("T_L2muon_Eta:T_L2muon_Phi>>L2etaPhi");
    theOutFile->cd();
    L2etaPhi->Write();
    
    TH2F *L3etaPhi = new TH2F("L3etaPhi","",100,-3.1416,3.1416,100,-2.5,2.5);
    chain->Draw("T_L3muon_Eta:T_L3muon_Phi>>L3etaPhi");
    theOutFile->cd();
    L3etaPhi->Write();

    theOutFile->Close();
}