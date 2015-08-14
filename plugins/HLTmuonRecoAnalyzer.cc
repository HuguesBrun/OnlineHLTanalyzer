#include "HLTmuonRecoAnalyzer.h"

//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
HLTmuonRecoAnalyzer::HLTmuonRecoAnalyzer(const edm::ParameterSet& iConfig)

{
   //now do what ever initialization is needed
    isMC_       = iConfig.getParameter<bool>("isMC");
    pathsToSave_ = iConfig.getParameter<std::vector<std::string> >("pathsToSave");
    L2muonsLabel_=  iConfig.getParameter<edm::InputTag>("L2muonsCollection");
    L3muonsLabel_=  iConfig.getParameter<edm::InputTag>("L3muonsCollection");
    TKmuonsLabel_=  iConfig.getParameter<edm::InputTag>("TKmuonsCollection");
    beamSpotLabel_=  iConfig.getParameter<edm::InputTag>("BeamSpot");
    theMuonRecHitBuilderName_ = iConfig.getParameter<std::string>("MuonRecHitBuilder");
    outputFile_     = iConfig.getParameter<std::string>("outputFile");
    rootFile_       = TFile::Open(outputFile_.c_str(),"RECREATE");
}


HLTmuonRecoAnalyzer::~HLTmuonRecoAnalyzer()
{
 
   // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
HLTmuonRecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    using namespace edm;
    using namespace std;
    beginEvent();
    
   // edm::ESHandle<TransientTrackingRecHitBuilder> theMuonRecHitBuilder;
    //iSetup.get<TransientRecHitRecord>().get(theMuonRecHitBuilderName_,theMuonRecHitBuilder);
   
    edm::ESHandle<GlobalTrackingGeometry> trackingGeometry; 
    iSetup.get<GlobalTrackingGeometryRecord>().get(trackingGeometry); 

	edm::ESHandle<MagneticField> magField;
    iSetup.get<IdealMagneticFieldRecord>().get(magField);
 
   /* edm::Handle<std::vector<reco::Track> > L2Tracks;
    iEvent.getByLabel(L2muonsLabel_, L2Tracks);*/
    
    edm::InputTag trackingParticlesTag = edm::InputTag("mix","MergedTrackTruth");
    edm::Handle<TrackingParticleCollection>  TPCollectionH ;
    TrackingParticleCollection tPC;
    
    if (isMC_){
        iEvent.getByLabel(trackingParticlesTag,TPCollectionH);
        if (TPCollectionH.isValid()) tPC   = *(TPCollectionH.product());
        else cout << "not found tracking particles collection" << endl;
    }
    
    Handle<reco::RecoChargedCandidateCollection> allMuons;
    iEvent.getByLabel(L2muonsLabel_, allMuons);
    
    Handle<reco::RecoChargedCandidateCollection> allL3Muons;
    iEvent.getByLabel(L3muonsLabel_, allL3Muons);
    
    
    edm::Handle<reco::MuonCollection> allTKMuons;
    iEvent.getByLabel(TKmuonsLabel_, allTKMuons);
    
    Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel(beamSpotLabel_, beamSpotHandle);
    const reco::BeamSpot& theBeamSpot = *beamSpotHandle;
    reco::BeamSpot::Point beamSpot = beamSpotHandle->position();
    
    edm::Handle<trigger::TriggerFilterObjectWithRefs> L3muonFilter;
    iEvent.getByLabel(edm::InputTag("hltL3fL1sMu16L1f0L2f10QL3Filtered20Q","","HLTfancy"), L3muonFilter);
    
    edm::Handle<reco::RecoChargedCandidateIsolationMap> IsoCaloMapECAL;
    iEvent.getByLabel(edm::InputTag("hltMuonEcalPFClusterIsoForMuons","","HLTfancy"),IsoCaloMapECAL);
    
    edm::Handle<reco::RecoChargedCandidateIsolationMap> IsoCaloMapHCAL;
    iEvent.getByLabel(edm::InputTag("hltMuonHcalPFClusterIsoForMuons","","HLTfancy"),IsoCaloMapHCAL);
    
    edm::Handle<reco::IsoDepositMap> IsoTkMap;
    iEvent.getByLabel(edm::InputTag("hltMuonTkRelIsolationCut0p09Map","trkIsoDeposits","HLTfancy"),IsoTkMap);
    
    // get rho
    edm::Handle<double> rhoHandle;
    iEvent.getByLabel(edm::InputTag("hltFixedGridRhoFastjetAllCaloForMuons","","HLTfancy"), rhoHandle);
    
    
    
     edm::Handle<LumiScalersCollection> lumiScalers;
    if (!isMC_){
        iEvent.getByLabel( InputTag("scalersRawToDigi", "", "HLTfancy"), lumiScalers );
        if (lumiScalers.isValid()) {
            LumiScalersCollection::const_iterator it3 = lumiScalers->begin();
           // T_Event_LuminosityBlock = it3->sectionNumber();
           // std::cout << "section number" << T_Event_LuminosityBlock << std::endl;
            T_Event_InstLumi = it3->instantLumi();
            cout << "hello the inst. lumi is " << T_Event_InstLumi << endl;
        }
    }
    

    
    edm::ESHandle<TrackerGeometry> geom;
    iSetup.get<TrackerDigiGeometryRecord>().get( geom );
    const TrackerGeometry& theTracker( *geom );
    
    bool changedConfig = false;
    if (!hltConfig.init(iEvent.getRun(), iSetup, "HLTfancy", changedConfig)) {
        cout << "Initialization of HLTConfigProvider failed!!" << endl;
        return;
    }
    if (changedConfig){
        unsigned int nbPaths = pathsToSave_.size();
        std::cout << "the curent menu is " << hltConfig.tableName() << std::endl;
        for (unsigned int itePath=0 ; itePath<nbPaths ; itePath++){
            for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
                if (TString(hltConfig.triggerNames()[j]).Contains(pathsToSave_.at(itePath))){
                    triggerBits_.push_back(j);
                    cout << "found the path " << pathsToSave_.at(itePath) << endl;
                }
                
            }
        }
        if (triggerBits_.size() <nbPaths) cout << "an HLT paths is not found ! ! " << endl;
        
    }
    
    bool changedConfigData = false;
    if (!hltConfigData.init(iEvent.getRun(), iSetup, "HLT", changedConfigData)) {
        cout << "Initialization of HLTConfigProvider failed!!" << endl;
        return;
    }
    if (changedConfigData){
        unsigned int nbPaths = pathsToSave_.size();
        std::cout << "the curent menu is " << hltConfigData.tableName() << std::endl;
        for (unsigned int itePath=0 ; itePath<nbPaths ; itePath++){
            for (size_t j = 0; j < hltConfigData.triggerNames().size(); j++) {
                if (TString(hltConfigData.triggerNames()[j]).Contains(pathsToSave_.at(itePath))){
                    triggerBitsData_.push_back(j);
                    cout << "found the path " << pathsToSave_.at(itePath) << endl;
                }
                
            }
        }
        if (triggerBits_.size() <nbPaths) cout << "an HLT paths is not found ! ! " << endl;
        
    }
    
    
    edm::InputTag triggerResultsLabel = edm::InputTag("TriggerResults", "", "HLTfancy");
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResultsLabel, triggerResults);
    
    edm::InputTag triggerResultsDataLabel = edm::InputTag("TriggerResults", "", "HLT");
    edm::Handle<edm::TriggerResults> triggerResultsData;
    iEvent.getByLabel(triggerResultsDataLabel, triggerResultsData);
    
    
    if (isMC_){
        edm::Handle<GenEventInfoProduct> genEvent;
        iEvent.getByLabel("generator", genEvent);
        
        
        if (genEvent->binningValues().size()>0)  T_Event_PUptHat2 = genEvent->binningValues()[0]; else T_Event_PUptHat2=-1;
        T_Event_PUptHat=genEvent->qScale();
        
        float truePu=0.;
        Handle<std::vector< PileupSummaryInfo > > puInfo;
        try {
            iEvent.getByLabel("addPileupInfo",puInfo);
            std::vector<PileupSummaryInfo>::const_iterator PVI;
            //The in-time crossing is getBunchCrossing = 0; negative ones are early, positive ones are late.
            for(PVI = puInfo->begin(); PVI != puInfo->end(); ++PVI) {
                
                //std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
                if(PVI->getBunchCrossing()==0){
                    T_Event_nPU =PVI->getPU_NumInteractions();
                    T_Event_nTruePU=PVI->getTrueNumInteractions();
                    float pu_pT_hat_max = -1;
                    for(const auto& pu_pT_hat :  PVI->getPU_pT_hats()){
                        if (pu_pT_hat>pu_pT_hat_max) pu_pT_hat_max = pu_pT_hat;
                    }
                    T_Event_PUeventPtHat = pu_pT_hat_max;
                    if (T_Event_PUeventPtHat>genEvent->qScale()) T_Event_PUdominated = 1; else T_Event_PUdominated=0;
                    
                }
                
                else if(PVI->getBunchCrossing()==-1){
                    T_Event_nPUm=PVI->getPU_NumInteractions();
                }
                else if(PVI->getBunchCrossing()==1){
                    T_Event_nPUp=PVI->getPU_NumInteractions();
                }
                truePu += PVI->getTrueNumInteractions();
            }
        } catch (...) {}
        T_Event_AveNTruePU=truePu;

        
    }
    
    
  //  std::cout << "passing the trigger=" << triggerResults->accept(triggerBit) << std::endl;
   // T_Event_PassL3muon = triggerResults->accept(triggerBit);
    T_Event_RunNumber =  iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    T_Event_LuminosityBlock = iEvent.id().luminosityBlock();
    
    for (unsigned int itePath = 0 ; itePath < triggerBits_.size() ; itePath++){
        if (triggerResults->accept(triggerBits_.at(itePath))) {
            edm::LogVerbatim("HLTmuonRecoAnalyzer") << "passing path " << itePath << endl;
            T_Event_pathsFired->push_back(1);
        }
        else T_Event_pathsFired->push_back(0);
    }
    
    for (unsigned int itePath = 0 ; itePath < triggerBitsData_.size() ; itePath++){
        if (triggerResultsData->accept(triggerBitsData_.at(itePath))) {
            edm::LogVerbatim("HLTmuonRecoAnalyzer") << "passing path " << itePath << endl;
            T_Event_pathsFiredData->push_back(1);
        }
        else T_Event_pathsFiredData->push_back(0);
    }
    
    BSx = theBeamSpot.x0();
    BSy = theBeamSpot.y0();
    BSz = theBeamSpot.z0();
    
    hltRho = -1;
    if (rhoHandle.isValid()) hltRho = *(rhoHandle.product());

    
    if (allMuons.isValid()){
    for(reco::RecoChargedCandidateCollection::const_iterator cand=allMuons->begin(); cand!=allMuons->end(); cand++){
        reco::TrackRef mu = cand->get<reco::TrackRef>();
        //const reco::Track theTrack = (*L2Tracks)[i];
        T_L2muon_Pt->push_back(mu->pt());
        T_L2muon_Eta->push_back(mu->eta());
        T_L2muon_Phi->push_back(mu->phi());
        T_L2muon_dz->push_back(mu->dz());
        T_L2muon_dxy->push_back(mu->dxy());
        T_L2muon_dzBS->push_back(mu->dz(beamSpot));
        T_L2muon_dxyBS->push_back(mu->dxy(beamSpot));
        T_L2muon_nbStationWithHits->push_back(mu->hitPattern().muonStationsWithAnyHits());
        T_L2muon_nbStationWithHitsDT->push_back(mu->hitPattern().dtStationsWithAnyHits());
        T_L2muon_nbStationWithHitsCSC->push_back(mu->hitPattern().cscStationsWithAnyHits());
        T_L2muon_nbValidHits->push_back(mu->numberOfValidHits());
        double abspar0 = std::abs(mu->parameter(0));
        double ptLx = mu->pt();
      // convert 50% efficiency threshold to 90% efficiency threshold
        if(abspar0 > 0) ptLx += 1.*mu->error(0)/abspar0*mu->pt();
        T_L2muon_ptLx->push_back(ptLx);
    }
    }
    
    std::vector<reco::RecoChargedCandidateRef> passL3filterMuon;
    
    if(L3muonFilter.isValid()){
        L3muonFilter->getObjects(trigger::TriggerMuon, passL3filterMuon);
    }
    
    if (allL3Muons.isValid()){
        int m = 0;
        for(reco::RecoChargedCandidateCollection::const_iterator cand=allL3Muons->begin(); cand!=allL3Muons->end(); cand++){
            reco::TrackRef l3track = cand->track();
            T_L3muon_Pt->push_back(cand->pt());
            T_L3muon_Px->push_back(cand->px());
            T_L3muon_Py->push_back(cand->py());
            T_L3muon_Pz->push_back(cand->pz());
            T_L3muon_Energy->push_back(cand->energy());
            T_L3muon_Eta->push_back(cand->eta());
            T_L3muon_Phi->push_back(cand->phi());
            T_L3muon_dr->push_back(std::abs( (- (cand->vx()-theBeamSpot.x0()) * cand->py() + (cand->vy()-theBeamSpot.y0()) * cand->px() ) / cand->pt() ));
            T_L3muon_dz->push_back(std::abs((cand->vz()-theBeamSpot.z0()) - ((cand->vx()-theBeamSpot.x0())*cand->px()+(cand->vy()-theBeamSpot.y0())*cand->py())/cand->pt() * cand->pz()/cand->pt()));
            T_L3muon_dxyBS->push_back(std::abs(l3track->dxy(beamSpot)));
            T_L3muon_Chi2->push_back(l3track->normalizedChi2());
            const reco::RecoChargedCandidateRef refCand(allL3Muons,m);
            bool isPassingTheFilter=false;
            reco::RecoChargedCandidateRef ref;
            for (size_t j = 0 ; j < passL3filterMuon.size() ; j++){
                ref = passL3filterMuon[j];
                if (ref==refCand) isPassingTheFilter = true;
            }
            T_L3muon_PassingL3filter->push_back(isPassingTheFilter);
            if (isPassingTheFilter&&IsoCaloMapECAL.isValid()) T_L3muon_pfEcal->push_back((*IsoCaloMapECAL)[refCand]);
            else T_L3muon_pfEcal->push_back(-1);
            if (isPassingTheFilter&&IsoCaloMapHCAL.isValid()) T_L3muon_pfHcal->push_back((*IsoCaloMapHCAL)[refCand]);
            else T_L3muon_pfHcal->push_back(-1);
            if (isPassingTheFilter&&IsoTkMap.isValid()){
                reco::IsoDeposit theTkIsolation = (*IsoTkMap)[ref];
                T_L3muon_trkIso->push_back(theTkIsolation.depositWithin(0.3));
            }
            else T_L3muon_trkIso->push_back(-1);
            m++;
        }
    }
    
    muon::SelectionType m_trkMuonId = muon::SelectionType(0);
    
    if (allTKMuons.isValid()){
        //int m = 0;
        for(reco::MuonCollection::const_iterator muons=allTKMuons->begin(); muons!=allTKMuons->end(); muons++){
            //reco::TrackRef tktrack = cand->track();
            
            T_Tkmuon_Pt->push_back(muons->pt());
            T_Tkmuon_Px->push_back(muons->px());
            T_Tkmuon_Py->push_back(muons->py());
            T_Tkmuon_Pz->push_back(muons->pz());
            T_Tkmuon_Energy->push_back(muons->energy());
            T_Tkmuon_Eta->push_back(muons->eta());
            T_Tkmuon_Phi->push_back(muons->phi());
            T_Tkmuon_Phi->push_back(muons->phi());
            T_Tkmuon_type->push_back( muons->type());
            T_Tkmuon_matchedStation->push_back(muons->numberOfMatchedStations());
            if ( !muons->innerTrack().isNull() ) T_Tkmuon_validHit->push_back( muons->innerTrack()->numberOfValidHits());
            else T_Tkmuon_validHit->push_back(-1);
            
            if ( !muons->globalTrack().isNull() ) {
                T_Tkmuon_normChi2->push_back(muons->globalTrack()->normalizedChi2());
                T_Tkmuon_validMuonHit->push_back(muons->globalTrack()->hitPattern().numberOfValidMuonHits());
            }
            else {
                T_Tkmuon_normChi2->push_back(-1);
                T_Tkmuon_validMuonHit->push_back(-1);
            }
            T_Tkmuon_goodTrackerMuon->push_back(( muons->isTrackerMuon() && !muon::isGoodMuon(*muons,m_trkMuonId) ));
  
        }
        
    }
    
    bool needToVeto = false;
    if (isMC_){
        for (TrackingParticleCollection::size_type i=0; i<tPC.size(); i++) {
            TrackingParticleRef trpart(TPCollectionH, i);
            
            /*cout << "tp eta=" << trpart->eta() << endl;
             cout << "tp phi=" << trpart->phi() << endl;
             
             cout << "tp pT=" << trpart->pt() << "pdgID=" << trpart->pdgId() << " vx=" << trpart->vx() << " vy=" << trpart->vy() << endl;*/
            if (fabs(trpart->pdgId())!=13) continue;
            //  cout << "pdg" << endl;
            float rhoVtx = sqrt(pow(trpart->vx(),2)+pow(trpart->vy(),2));
            float zVtx = fabs(trpart->vz());
            //  cout << "rhoVtx=" << rhoVtx << endl;
            //cout << "zVtx=" << zVtx << endl;
            if ((trpart->pt()>5)&&(fabs(trpart->eta())<2.5)&&(rhoVtx<2000)&&(zVtx<4000)) needToVeto = true;
            
        }
    }
    T_Event_antiMuEnrichementVeto = needToVeto;
    
    /*edm::Handle<L3MuonTrajectorySeedCollection> L3seedCollection;
    edm::InputTag L3seedCollectionTag = edm::InputTag("hltL3TrajSeedOIStateNoVtx");
    iEvent.getByLabel(L3seedCollectionTag, L3seedCollection);
    cout << L3seedCollection.isValid() << endl;
    
    vector<L3MuonTrajectorySeed> theTrajSeeds = *(L3seedCollection);
    int sizeTrajSeeds = theTrajSeeds.size();
    cout << "sizeTrajSeeds" << sizeTrajSeeds << endl;
    
    for(L3MuonTrajectorySeedCollection::const_iterator seed = L3seedCollection->begin(); seed!=L3seedCollection->end(); ++seed){
			PTrajectoryStateOnDet state = seed->startingState();
			DetId seedDetId( state.detId() );
			const GeomDet* gdet = trackingGeometry->idToDet( seedDetId );
			TrajectoryStateOnSurface tsos = trajectoryStateTransform::transientState(state, &(gdet->surface()), magField.product());
					float pt = tsos.globalMomentum().perp();
					float eta = tsos.globalMomentum().eta();
                    float phi = tsos.globalMomentum().phi();
                    float xSeed = tsos.globalPosition().x();
                    float ySeed = tsos.globalPosition().y();
                    float zSeed = tsos.globalPosition().z();
        cout << "pt="<< pt  << " eta=" << eta << " phi=" << phi << endl;
        cout << "xSeed="<< xSeed  << " ySeed=" << ySeed << " zSeed=" << zSeed << endl;
        T_L3seed_Pt->push_back(pt);
        T_L3seed_Eta->push_back(eta);
        T_L3seed_Phi->push_back(phi);
        T_L3seed_X->push_back(xSeed);
        T_L3seed_Y->push_back(ySeed);
        T_L3seed_Z->push_back(zSeed);
    }
 
   // if (T_Event_PassL3muon) mytree_->Fill();*/
    mytree_->Fill();
    endEvent();

}


// ------------ method called once each job just before starting event loop  ------------
void 
HLTmuonRecoAnalyzer::beginJob()
{
    mytree_ = new TTree("eventsTree","");
    
    mytree_->Branch("T_Event_RunNumber", &T_Event_RunNumber, "T_Event_RunNumber/I");
    mytree_->Branch("T_Event_EventNumber", &T_Event_EventNumber, "T_Event_EventNumber/I");
    mytree_->Branch("T_Event_LuminosityBlock", &T_Event_LuminosityBlock, "T_Event_LuminosityBlock/I");
    mytree_->Branch("T_Event_InstLumi", &T_Event_InstLumi, "T_Event_InstLumi/F");
    mytree_->Branch("T_nbStripClusters", &T_nbStripClusters, "T_nbStripClusters/I");
    mytree_->Branch("T_Event_PassL3muon", &T_Event_PassL3muon, "T_Event_PassL3muon/I");
    
    mytree_->Branch("T_Event_nPU", &T_Event_nPU, "T_Event_nPU/I");
    mytree_->Branch("T_Event_nTruePU", &T_Event_nTruePU, "T_Event_nTruePU/I");
    mytree_->Branch("T_Event_PUeventPtHat", &T_Event_PUeventPtHat, "T_Event_PUeventPtHat/F");
    mytree_->Branch("T_Event_PUdominated", &T_Event_PUdominated, "T_Event_PUdominated/I");
    mytree_->Branch("T_Event_AveNTruePU", &T_Event_AveNTruePU, "T_Event_AveNTruePU/I");
    mytree_->Branch("T_Event_PUptHat2", &T_Event_PUptHat2, "T_Event_PUptHat2/F");
    mytree_->Branch("T_Event_PUptHat", &T_Event_PUptHat, "T_Event_PUptHat/F");
    mytree_->Branch("T_Event_nPUm", &T_Event_nPUm, "T_Event_nPUm/I");
    mytree_->Branch("T_Event_nPUp", &T_Event_nPUp, "T_Event_nPUp/I");

    
    
    mytree_->Branch("T_Event_antiMuEnrichementVeto", &T_Event_antiMuEnrichementVeto, "T_Event_antiMuEnrichementVeto/I");
    
    mytree_->Branch("T_Event_pathsFired", "std::vector<int>", &T_Event_pathsFired);
    mytree_->Branch("T_Event_pathsFiredData", "std::vector<int>", &T_Event_pathsFiredData);
    
    mytree_->Branch("BSx", &BSx, "BSx/F");
    mytree_->Branch("BSy", &BSy, "BSy/F");
    mytree_->Branch("BSz", &BSz, "BSz/F");
    
    mytree_->Branch("hltRho", &hltRho, "hltRho/F");

    
    mytree_->Branch("T_L2muon_Pt", "std::vector<float>", &T_L2muon_Pt);
    mytree_->Branch("T_L2muon_Eta", "std::vector<float>", &T_L2muon_Eta);
    mytree_->Branch("T_L2muon_Phi", "std::vector<float>", &T_L2muon_Phi);
    mytree_->Branch("T_L2muon_dxy", "std::vector<float>", &T_L2muon_dxy);
    mytree_->Branch("T_L2muon_dz", "std::vector<float>", &T_L2muon_dz);
    mytree_->Branch("T_L2muon_dxyBS", "std::vector<float>", &T_L2muon_dxyBS);
    mytree_->Branch("T_L2muon_dzBS", "std::vector<float>", &T_L2muon_dzBS);
    mytree_->Branch("T_L2muon_ptLx", "std::vector<float>", &T_L2muon_ptLx);
    mytree_->Branch("T_L2muon_nbStationWithHits", "std::vector<int>", &T_L2muon_nbStationWithHits);
    mytree_->Branch("T_L2muon_nbStationWithHitsDT", "std::vector<int>", &T_L2muon_nbStationWithHitsDT);
    mytree_->Branch("T_L2muon_nbStationWithHitsCSC", "std::vector<int>", &T_L2muon_nbStationWithHitsCSC);
    mytree_->Branch("T_L2muon_nbValidHits", "std::vector<int>", &T_L2muon_nbValidHits);
    
    mytree_->Branch("T_L3muon_PassingL3filter", "std::vector<int>", &T_L3muon_PassingL3filter);
    mytree_->Branch("T_L3muon_Pt", "std::vector<float>", &T_L3muon_Pt);
    mytree_->Branch("T_L3muon_Px", "std::vector<float>", &T_L3muon_Px);
    mytree_->Branch("T_L3muon_Py", "std::vector<float>", &T_L3muon_Py);
    mytree_->Branch("T_L3muon_Pz", "std::vector<float>", &T_L3muon_Pz);
    mytree_->Branch("T_L3muon_Energy", "std::vector<float>", &T_L3muon_Energy);
    mytree_->Branch("T_L3muon_Eta", "std::vector<float>", &T_L3muon_Eta);
    mytree_->Branch("T_L3muon_Phi", "std::vector<float>", &T_L3muon_Phi);
    mytree_->Branch("T_L3muon_dr", "std::vector<float>", &T_L3muon_dr);
    mytree_->Branch("T_L3muon_dz", "std::vector<float>", &T_L3muon_dz);
    mytree_->Branch("T_L3muon_dxyBS", "std::vector<float>", &T_L3muon_dxyBS);
    mytree_->Branch("T_L3muon_Chi2", "std::vector<float>", &T_L3muon_Chi2);
    mytree_->Branch("T_L3muon_pfEcal", "std::vector<float>", &T_L3muon_pfEcal);
    mytree_->Branch("T_L3muon_pfHcal", "std::vector<float>", &T_L3muon_pfHcal);
    mytree_->Branch("T_L3muon_trkIso", "std::vector<float>", &T_L3muon_trkIso);

    mytree_->Branch("T_Tkmuon_Pt", "std::vector<float>", &T_Tkmuon_Pt);
    mytree_->Branch("T_Tkmuon_Px", "std::vector<float>", &T_Tkmuon_Px);
    mytree_->Branch("T_Tkmuon_Py", "std::vector<float>", &T_Tkmuon_Py);
    mytree_->Branch("T_Tkmuon_Pz", "std::vector<float>", &T_Tkmuon_Pz);
    mytree_->Branch("T_Tkmuon_Energy", "std::vector<float>", &T_Tkmuon_Energy);
    mytree_->Branch("T_Tkmuon_Eta", "std::vector<float>", &T_Tkmuon_Eta);
    mytree_->Branch("T_Tkmuon_Phi", "std::vector<float>", &T_Tkmuon_Phi);
    mytree_->Branch("T_Tkmuon_matchedStation", "std::vector<int>", &T_Tkmuon_matchedStation);
    mytree_->Branch("T_Tkmuon_normChi2", "std::vector<float>", &T_Tkmuon_normChi2);
    mytree_->Branch("T_Tkmuon_validMuonHit", "std::vector<int>", &T_Tkmuon_validMuonHit);
    mytree_->Branch("T_Tkmuon_validHit", "std::vector<int>", &T_Tkmuon_validHit);
    mytree_->Branch("T_Tkmuon_type", "std::vector<int>", &T_Tkmuon_type);
    mytree_->Branch("T_Tkmuon_goodTrackerMuon", "std::vector<int>", &T_Tkmuon_goodTrackerMuon);


}

// ------------ method called once each job just after ending the event loop  ------------
void
HLTmuonRecoAnalyzer::endJob() 
{
    rootFile_->Write();
    rootFile_->Close();
}

void
HLTmuonRecoAnalyzer::beginEvent()
{
    T_Event_pathsFired = new std::vector<int>;
    T_Event_pathsFiredData = new std::vector<int>;
    
    T_L2muon_Pt = new std::vector<float>;
    T_L2muon_Eta = new std::vector<float>;
    T_L2muon_Phi = new std::vector<float>;
    T_L2muon_dxy = new std::vector<float>;
    T_L2muon_dz = new std::vector<float>;
    T_L2muon_dxyBS = new std::vector<float>;
    T_L2muon_dzBS = new std::vector<float>;
    T_L2muon_ptLx = new std::vector<float>;
    T_L2muon_nbStationWithHits = new std::vector<int>;
    T_L2muon_nbStationWithHitsDT = new std::vector<int>;
    T_L2muon_nbStationWithHitsCSC = new std::vector<int>;
    T_L2muon_nbValidHits = new std::vector<int>;
    
    T_L3muon_PassingL3filter = new std::vector<int>;
    T_L3muon_Pt = new std::vector<float>;
    T_L3muon_Px = new std::vector<float>;
    T_L3muon_Py = new std::vector<float>;
    T_L3muon_Pz = new std::vector<float>;
    T_L3muon_Energy = new std::vector<float>;
    T_L3muon_Eta = new std::vector<float>;
    T_L3muon_Phi = new std::vector<float>;
    T_L3muon_dr = new std::vector<float>;
    T_L3muon_dz = new std::vector<float>;
    T_L3muon_dxyBS = new std::vector<float>;
    T_L3muon_Chi2 = new std::vector<float>;
    T_L3muon_pfEcal = new std::vector<float>;
    T_L3muon_pfHcal = new std::vector<float>;
    T_L3muon_trkIso = new std::vector<float>;

    T_Tkmuon_Pt = new std::vector<float>;
    T_Tkmuon_Px = new std::vector<float>;
    T_Tkmuon_Py = new std::vector<float>;
    T_Tkmuon_Pz = new std::vector<float>;
    T_Tkmuon_Energy = new std::vector<float>;
    T_Tkmuon_Eta = new std::vector<float>;
    T_Tkmuon_Phi = new std::vector<float>;
    T_Tkmuon_matchedStation = new std::vector<int>;
    T_Tkmuon_normChi2 = new std::vector<float>;
    T_Tkmuon_validMuonHit = new std::vector<int>;
    T_Tkmuon_validHit = new std::vector<int>;
    T_Tkmuon_type = new std::vector<int>;
    T_Tkmuon_goodTrackerMuon = new std::vector<int>;

}

void
HLTmuonRecoAnalyzer::endEvent()
{
    delete T_Event_pathsFired;
    delete T_Event_pathsFiredData;
    
    delete T_L2muon_Pt;
    delete T_L2muon_Eta;
    delete T_L2muon_Phi;
    delete T_L2muon_dxy;
    delete T_L2muon_dz;
    delete T_L2muon_dxyBS;
    delete T_L2muon_dzBS;
    delete T_L2muon_ptLx;
    delete T_L2muon_nbStationWithHits;
    delete T_L2muon_nbStationWithHitsDT;
    delete T_L2muon_nbStationWithHitsCSC;
    delete T_L2muon_nbValidHits;

    
    delete T_L3muon_PassingL3filter;
    delete T_L3muon_Pt;
    delete T_L3muon_Px;
    delete T_L3muon_Py;
    delete T_L3muon_Pz;
    delete T_L3muon_Energy;
    delete T_L3muon_Eta;
    delete T_L3muon_Phi;
    delete T_L3muon_dr;
    delete T_L3muon_dz;
    delete T_L3muon_dxyBS;
    delete T_L3muon_Chi2;
    delete T_L3muon_pfEcal;
    delete T_L3muon_pfHcal;
    delete T_L3muon_trkIso;
    
    delete T_Tkmuon_Pt;
    delete T_Tkmuon_Px;
    delete T_Tkmuon_Py;
    delete T_Tkmuon_Pz;
    delete T_Tkmuon_Energy;
    delete T_Tkmuon_Eta;
    delete T_Tkmuon_Phi;
    delete T_Tkmuon_matchedStation;
    delete T_Tkmuon_normChi2;
    delete T_Tkmuon_validMuonHit;
    delete T_Tkmuon_validHit;
    delete T_Tkmuon_type;
    delete T_Tkmuon_goodTrackerMuon;


}

// ------------ method called when starting to processes a run  ------------
/*
void 
HLTmuonRecoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a run  ------------
/*
void 
HLTmuonRecoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when starting to processes a luminosity block  ------------
/*
void 
HLTmuonRecoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method called when ending the processing of a luminosity block  ------------
/*
void 
HLTmuonRecoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}
*/

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
HLTmuonRecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(HLTmuonRecoAnalyzer);
