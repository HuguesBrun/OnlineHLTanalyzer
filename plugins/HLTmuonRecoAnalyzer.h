// -*- C++ -*-
//
// Package:    hugues/HLTmuonRecoAnalyzer
// Class:      HLTmuonRecoAnalyzer
// 
/**\class HLTmuonRecoAnalyzer HLTmuonRecoAnalyzer.cc hugues/HLTmuonRecoAnalyzer/plugins/HLTmuonRecoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Hugues Louis Brun
//         Created:  Sat, 14 Feb 2015 10:36:11 GMT
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"

#include "DataFormats/TrackReco/interface/TrackFwd.h"
#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/Common/interface/RefToBase.h"

#include "DataFormats/RecoCandidate/interface/RecoChargedCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoChargedCandidateFwd.h"

#include "DataFormats/BeamSpot/interface/BeamSpot.h"

#include "TrackingTools/TransientTrack/interface/TransientTrack.h"
#include "TrackingTools/Records/interface/TransientRecHitRecord.h"
#include "TrackingTools/TransientTrackingRecHit/interface/TransientTrackingRecHitBuilder.h"
 #include "CalibTracker/Records/interface/SiStripDetCablingRcd.h"
 #include "CalibFormats/SiStripObjects/interface/SiStripDetCabling.h"
 #include "CondFormats/DataRecord/interface/SiStripCondDataRecords.h"
 #include "CondFormats/SiStripObjects/interface/SiStripNoises.h"
 #include "CalibFormats/SiStripObjects/interface/SiStripGain.h"
 #include "CalibFormats/SiStripObjects/interface/SiStripQuality.h"
 #include "DataFormats/SiStripCluster/interface/SiStripClusterCollection.h"
 #include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
// #include "DQM/SiStripCommon/interface/SiStripFolderOrganizer.h"
// #include "DQM/SiStripCommon/interface/SiStripHistoId.h"
 //#include "DQM/SiStripMonitorCluster/interface/SiStripMonitorCluster.h"
 //#include "DQMServices/Core/interface/DQMStore.h"
 //#include "DQMServices/Core/interface/MonitorElement.h"
 #include "DataFormats/SiStripDetId/interface/StripSubdetector.h"
 #include "DataFormats/SiStripDetId/interface/SiStripSubStructure.h"
#include "CalibTracker/SiStripCommon/interface/SiStripDCSStatus.h"
 #include "DataFormats/SiPixelCluster/interface/SiPixelCluster.h" 
 #include "DPGAnalysis/SiStripTools/interface/APVCyclePhaseCollection.h"
 #include "DPGAnalysis/SiStripTools/interface/EventWithHistory.h"
 
#include "CommonTools/TriggerUtils/interface/GenericTriggerEventFlag.h"

#include "TrackingTools/TrajectoryState/interface/TrajectoryStateTransform.h"
#include "Geometry/Records/interface/GlobalTrackingGeometryRecord.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeed.h"
#include "DataFormats/MuonSeed/interface/L3MuonTrajectorySeedCollection.h"

#include "MagneticField/Records/interface/IdealMagneticFieldRecord.h"
#include "MagneticField/Engine/interface/MagneticField.h"

#include "DataFormats/SiStripCluster/interface/SiStripCluster.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/StripGeomDetUnit.h"
#include "Geometry/CommonTopologies/interface/StripTopology.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Common/interface/DetSetVectorNew.h"

#include "HLTrigger/HLTcore/interface/HLTConfigProvider.h"
#include "DataFormats/Common/interface/TriggerResults.h"

#include "DataFormats/HLTReco/interface/TriggerEvent.h"
#include "DataFormats/HLTReco/interface/TriggerObject.h"

// root stuff 
#include "TH1D.h"
#include <map>
#include "TFile.h"
#include <math.h>
#include <sstream>
#include <string>
#include <stdlib.h>
#include <string.h>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TString.h"
#include "TTree.h"

//
// class declaration
//

class HLTmuonRecoAnalyzer : public edm::EDAnalyzer {
   public:
      explicit HLTmuonRecoAnalyzer(const edm::ParameterSet&);
      ~HLTmuonRecoAnalyzer();

      static void fillDescriptions(edm::ConfigurationDescriptions& descriptions);


   private:
      virtual void beginJob() override;
      virtual void analyze(const edm::Event&, const edm::EventSetup&) override;
      virtual void endJob() override;
    
    virtual void beginEvent();
    virtual void endEvent();
    

      //virtual void beginRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void endRun(edm::Run const&, edm::EventSetup const&) override;
      //virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;
      //virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&) override;

      // ----------member data ---------------------------
    edm::InputTag L2muonsLabel_;
    edm::InputTag L3muonsLabel_;
    edm::InputTag beamSpotLabel_;
    std::string theMuonRecHitBuilderName_;

    
    // Tree and outfile
    // root file to store histograms
    TFile*  rootFile_;
    std::string outputFile_; // output file name
    
    //Tree
    TTree* mytree_;
    
    int T_Event_RunNumber;
    int T_Event_EventNumber;
    int T_Event_LuminosityBlock;
    int T_Event_PassL3muon;
    
    std::vector<float>* T_L2muon_Pt;
    std::vector<float>* T_L2muon_Eta;
    std::vector<float>* T_L2muon_Phi;
    std::vector<float>* T_L2muon_dxy;
    std::vector<float>* T_L2muon_dz;
    std::vector<float>* T_L2muon_dxyBS;
    std::vector<float>* T_L2muon_dzBS;
    std::vector<float>* T_L2muon_ptLx;
    std::vector<int>* T_L2muon_nbStationWithHits;
    std::vector<int>* T_L2muon_nbStationWithHitsDT;
    std::vector<int>* T_L2muon_nbStationWithHitsCSC;
    std::vector<int>* T_L2muon_nbValidHits;
    
    
    
    std::vector<float>* T_L3muon_Pt;
    std::vector<float>* T_L3muon_Eta;
    std::vector<float>* T_L3muon_Phi;
    std::vector<float>* T_L3muon_dr;
    std::vector<float>* T_L3muon_dz;
    std::vector<float>* T_L3muon_dxyBS;
    std::vector<float>* T_L3muon_Chi2;


    
    
    int T_nbStripClusters;

    
    HLTConfigProvider hltConfig;
    int triggerBit;
    
};

//
// constants, enums and typedefs
//

//
// static data member definitions
//

