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
    L2muonsLabel_=  iConfig.getParameter<edm::InputTag>("L2muonsCollection");
    L3muonsLabel_=  iConfig.getParameter<edm::InputTag>("L3muonsCollection");
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
    
    Handle<reco::RecoChargedCandidateCollection> allMuons;
    iEvent.getByLabel(L2muonsLabel_, allMuons);
    
    Handle<reco::RecoChargedCandidateCollection> allL3Muons;
    iEvent.getByLabel(L3muonsLabel_, allL3Muons);
    
    Handle<reco::BeamSpot> beamSpotHandle;
    iEvent.getByLabel(beamSpotLabel_, beamSpotHandle);
    const reco::BeamSpot& theBeamSpot = *beamSpotHandle;
    reco::BeamSpot::Point beamSpot = beamSpotHandle->position();
    
    
    edm::ESHandle<TrackerGeometry> geom;
    iSetup.get<TrackerDigiGeometryRecord>().get( geom );
    const TrackerGeometry& theTracker( *geom );
    
    bool changedConfig = false;
    if (!hltConfig.init(iEvent.getRun(), iSetup, "HLTfancy", changedConfig)) {
        cout << "Initialization of HLTConfigProvider failed!!" << endl;
        return;
    }
    if (changedConfig){
        std::cout << "the curent menu is " << hltConfig.tableName() << std::endl;
        triggerBit = -1;
        for (size_t j = 0; j < hltConfig.triggerNames().size(); j++) {
            if (TString(hltConfig.triggerNames()[j]).Contains("HLT_L3Mu0_NoFilterNoVtx_v1")) triggerBit = j;
        }
        if (triggerBit == -1) cout << "HLT path not found" << endl;
        
    }
    
    
    
    edm::InputTag triggerResultsLabel = edm::InputTag("TriggerResults", "", "HLTfancy");
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResultsLabel, triggerResults);
    
  //  std::cout << "passing the trigger=" << triggerResults->accept(triggerBit) << std::endl;
   // T_Event_PassL3muon = triggerResults->accept(triggerBit);
    T_Event_RunNumber =  iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    
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
    if (allL3Muons.isValid()){
    for(reco::RecoChargedCandidateCollection::const_iterator cand=allL3Muons->begin(); cand!=allL3Muons->end(); cand++){
        reco::TrackRef l3track = cand->track();
        T_L3muon_Pt->push_back(cand->pt());
        T_L3muon_Eta->push_back(cand->eta());
        T_L3muon_Phi->push_back(cand->phi());
        T_L3muon_dr->push_back(std::abs( (- (cand->vx()-theBeamSpot.x0()) * cand->py() + (cand->vy()-theBeamSpot.y0()) * cand->px() ) / cand->pt() ));
        T_L3muon_dz->push_back(std::abs((cand->vz()-theBeamSpot.z0()) - ((cand->vx()-theBeamSpot.x0())*cand->px()+(cand->vy()-theBeamSpot.y0())*cand->py())/cand->pt() * cand->pz()/cand->pt()));
        T_L3muon_dxyBS->push_back(std::abs(l3track->dxy(beamSpot)));
        T_L3muon_Chi2->push_back(l3track->normalizedChi2());
        
        
    }
    }
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
    mytree_->Branch("T_nbStripClusters", &T_nbStripClusters, "T_nbStripClusters/I");
    mytree_->Branch("T_Event_PassL3muon", &T_Event_PassL3muon, "T_Event_PassL3muon/I");
    
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
    
    mytree_->Branch("T_L3muon_Pt", "std::vector<float>", &T_L3muon_Pt);
    mytree_->Branch("T_L3muon_Eta", "std::vector<float>", &T_L3muon_Eta);
    mytree_->Branch("T_L3muon_Phi", "std::vector<float>", &T_L3muon_Phi);
    mytree_->Branch("T_L3muon_dr", "std::vector<float>", &T_L3muon_dr);
    mytree_->Branch("T_L3muon_dz", "std::vector<float>", &T_L3muon_dz);
    mytree_->Branch("T_L3muon_dxyBS", "std::vector<float>", &T_L3muon_dxyBS);
    mytree_->Branch("T_L3muon_Chi2", "std::vector<float>", &T_L3muon_Chi2);



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
    
    T_L3muon_Pt = new std::vector<float>;
    T_L3muon_Eta = new std::vector<float>;
    T_L3muon_Phi = new std::vector<float>;
    T_L3muon_dr = new std::vector<float>;
    T_L3muon_dz = new std::vector<float>;
    T_L3muon_dxyBS = new std::vector<float>;
    T_L3muon_Chi2 = new std::vector<float>;


}

void
HLTmuonRecoAnalyzer::endEvent()
{
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

    
    delete T_L3muon_Pt;
    delete T_L3muon_Eta;
    delete T_L3muon_Phi;
    delete T_L3muon_dr;
    delete T_L3muon_dz;
    delete T_L3muon_dxyBS;
    delete T_L3muon_Chi2;
    


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
