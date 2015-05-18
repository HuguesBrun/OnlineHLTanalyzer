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

	edm::ESHandle<MagneticField> magField;				iSetup.get<IdealMagneticFieldRecord>().get(magField);
 
    edm::Handle<std::vector<reco::Track> > L2Tracks;
    iEvent.getByLabel(L2muonsLabel_, L2Tracks);
    
    edm::ESHandle<TrackerGeometry> geom;
    iSetup.get<TrackerDigiGeometryRecord>().get( geom );
    const TrackerGeometry& theTracker( *geom );
    
    bool changedConfig = false;
    if (!hltConfig.init(iEvent.getRun(), iSetup, "HLTX", changedConfig)) {
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
    
    
    
    edm::InputTag triggerResultsLabel = edm::InputTag("TriggerResults", "", "HLTX");
    edm::Handle<edm::TriggerResults> triggerResults;
    iEvent.getByLabel(triggerResultsLabel, triggerResults);
    
    std::cout << "passing the trigger=" << triggerResults->accept(triggerBit) << std::endl;
    T_Event_PassL3muon = triggerResults->accept(triggerBit);
    T_Event_RunNumber =  iEvent.id().run();
    T_Event_EventNumber = iEvent.id().event();
    
    int nbL2muons = (*L2Tracks).size();
    if (nbL2muons==0) return;
    
    for (int i=0 ; i < nbL2muons ; i++){
        const reco::Track theTrack = (*L2Tracks)[i];
        T_L2muon_Pt->push_back(theTrack.pt());
        T_L2muon_Eta->push_back(theTrack.innerMomentum().eta());
        T_L2muon_Phi->push_back(theTrack.innerMomentum().phi());
        T_L2muon_dz->push_back(theTrack.dz());
        T_L2muon_dxy->push_back(theTrack.dxy());
        int count = 0;
       //for (trackingRecHit_iterator hit = theTrack.recHitsBegin(); hit != theTrack.recHitsEnd(); ++hit) {
         //   if((*hit)->isValid()) {
           //     count++;
                // TransientTrackingRecHit::ConstRecHitPointer ttrh(theMuonRecHitBuilder->build(*hit));
           // }
       // }
        T_L2muon_Hits->push_back(count);
    }
 
    edm::InputTag clusterProducerStripLabel_ = edm::InputTag("hltSiStripRawToClustersFacility");
    //edm::InputTag clusterProducerStripLabel_ = edm::InputTag("siStripClusters");
     edm::Handle< edmNew::DetSetVector<SiStripCluster> > cluster_detsetvektor;
     iEvent.getByLabel(clusterProducerStripLabel_, cluster_detsetvektor);
    T_nbStripClusters = -1;
    if (cluster_detsetvektor.isValid()){
	const edmNew::DetSetVector<SiStripCluster> * StrC= cluster_detsetvektor.product();
        int NStripClusters= StrC->data().size();
        T_nbStripClusters = NStripClusters;
            }
    
    
    for (edmNew::DetSetVector<SiStripCluster>::const_iterator clustSet = cluster_detsetvektor->begin(); clustSet!=cluster_detsetvektor->end(); ++clustSet) {
        edmNew::DetSet<SiStripCluster>::const_iterator clustIt = clustSet->begin();
        edmNew::DetSet<SiStripCluster>::const_iterator end     = clustSet->end();
        
        DetId detIdObject( clustSet->detId() );
    //    edmNew::DetSetVector<SiStripCluster>::FastFiller spc(*output, detIdObject.rawId());
        const StripGeomDetUnit* theGeomDet = dynamic_cast<const StripGeomDetUnit*> (theTracker.idToDet(detIdObject) );
        const StripTopology * topol = dynamic_cast<const StripTopology*>(&(theGeomDet->specificTopology()));
        
        for(; clustIt!=end;++clustIt) {
            LocalPoint lpclust = topol->localPosition(clustIt->barycenter());
            GlobalPoint GPclust = theGeomDet->surface().toGlobal(Local3DPoint(lpclust.x(),lpclust.y(),lpclust.z()));
            double clustX = GPclust.x();
            double clustY = GPclust.y();
            double clustZ = GPclust.z();
            
            T_StripCluster_x->push_back(clustX);
            T_StripCluster_y->push_back(clustY);
            T_StripCluster_z->push_back(clustZ);
        }
    }
 
    edm::Handle<L3MuonTrajectorySeedCollection> L3seedCollection;
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
 
   // if (T_Event_PassL3muon) mytree_->Fill();
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
    mytree_->Branch("T_L2muon_Hits", "std::vector<int>", &T_L2muon_Hits);
    
    mytree_->Branch("T_L3seed_Pt", "std::vector<float>", &T_L3seed_Pt);
    mytree_->Branch("T_L3seed_Eta", "std::vector<float>", &T_L3seed_Eta);
    mytree_->Branch("T_L3seed_Phi", "std::vector<float>", &T_L3seed_Phi);
    mytree_->Branch("T_L3seed_X", "std::vector<float>", &T_L3seed_X);
    mytree_->Branch("T_L3seed_Y", "std::vector<float>", &T_L3seed_Y);
    mytree_->Branch("T_L3seed_Z", "std::vector<float>", &T_L3seed_Z);
    
    
    mytree_->Branch("T_StripCluster_x", "std::vector<float>", &T_StripCluster_x);
    mytree_->Branch("T_StripCluster_y", "std::vector<float>", &T_StripCluster_y);
    mytree_->Branch("T_StripCluster_z", "std::vector<float>", &T_StripCluster_z);


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
    T_L2muon_Hits = new std::vector<int>;
    
    T_L3seed_Pt = new std::vector<float>;
    T_L3seed_Eta = new std::vector<float>;
    T_L3seed_Phi = new std::vector<float>;
    T_L3seed_X = new std::vector<float>;
    T_L3seed_Y = new std::vector<float>;
    T_L3seed_Z = new std::vector<float>;

    T_StripCluster_x = new std::vector<float>;
    T_StripCluster_y = new std::vector<float>;
    T_StripCluster_z = new std::vector<float>;

}

void
HLTmuonRecoAnalyzer::endEvent()
{
    delete T_L2muon_Pt;
    delete T_L2muon_Eta;
    delete T_L2muon_Phi;
    delete T_L2muon_dxy;
    delete T_L2muon_dz;
    delete T_L2muon_Hits;
    
    delete T_L3seed_Pt;
    delete T_L3seed_Eta;
    delete T_L3seed_Phi;
    delete T_L3seed_X;
    delete T_L3seed_Y;
    delete T_L3seed_Z;
    
    delete T_StripCluster_x;
    delete T_StripCluster_y;
    delete T_StripCluster_z;


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
