import FWCore.ParameterSet.Config as cms

process = cms.Process("ALZ")

process.load("FWCore.MessageService.MessageLogger_cfi")
#process.load("Configuration.StandardSequences.GeometryRecoDB_cff")
process.load("Configuration.StandardSequences.MagneticField_0T_cff")
process.load("Configuration.StandardSequences.Services_cff")
##iprocess.load("Configuration.StandardSequences.Reconstruction_cff")
process.load("Configuration.StandardSequences.FrontierConditions_GlobalTag_cff")

process.GlobalTag.globaltag = 'GR_H_V43A'
if 'GlobalTag' in process.__dict__:
    from Configuration.AlCa.GlobalTag_condDBv2 import GlobalTag as customiseGlobalTag
    process.GlobalTag = customiseGlobalTag(process.GlobalTag, globaltag = 'GR_P_V49')
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_CONDITIONS'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    for pset in process.GlobalTag.toGet.value():
        pset.connect = pset.connect.value().replace('frontier://FrontierProd/', 'frontier://FrontierProd/')
    # fix for multi-run processing
    process.GlobalTag.RefreshEachRun = cms.untracked.bool( False )
    process.GlobalTag.ReconnectEachRun = cms.untracked.bool( False )
#process.load("SimMuon.MCTruth.MuonAssociatorByHits_cfi")
#process.load("SimTracker.TrackAssociation.TrackAssociatorByPosition_cff")
#process.load("SimMuon.MCTruth.MuonAssociatorByHitsESProducer_NoSimHits_cfi")

process.hltESPGlobalTrackingGeometryESProducer = cms.ESProducer( "GlobalTrackingGeometryESProducer" )
#process.hltESPGlobalDetLayerGeometry = cms.ESProducer( "GlobalDetLayerGeometryESProducer",
#  ComponentName = cms.string( "hltESPGlobalDetLayerGeometry" )
#)

process.MessageLogger.cerr.FwkReport.reportEvery = 1
process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )

process.source = cms.Source("PoolSource",
                            # replace 'myfile.root' with the source file you want to use
                            fileNames = cms.untracked.vstring(
                                                              #'/store/data/Commissioning2015/Cosmics/RECO/PromptReco-v1/000/234/390/00000/78284AF3-D9B6-E411-BA96-02163E012120.root'
                                                              #'/store/data/Commissioning2015/Cosmics/RECO/PromptReco-v1/000/233/238/00000/A8E21980-8CA9-E411-B51A-02163E0125DE.root'
                                                              #'file:/afs/cern.ch/work/h/hbrun/CMSSW_7_3_1_patch2_MWGR/src/testHLT/output/outputA_new.root'
                                                              'file:/afs/cern.ch/work/h/hbrun/CMSSW_7_3_1_patch2_MWGR/src/testHLT/testAgain/outputA.root'
                                                              #'file:/afs/cern.ch/work/h/hbrun/CMSSW_7_3_1_patch2_MWGR/src/testHLT/outputA.root'
                                                              ),
                            
#                            secondaryFileNames = cms.untracked.vstring('file:/afs/cern.ch/work/h/hbrun/CMSSW_7_3_1_patch2_MWGR/src/testHLT/output/outputA_new.root')
)

process.load("FWCore.MessageService.MessageLogger_cfi")



process.runL2seed = cms.EDAnalyzer('HLTmuonRecoAnalyzer',
                                L2muonsCollection = cms.InputTag("hltL2Muons",""),
                                MuonRecHitBuilder = cms.string("MuonRecHitBuilder"),
                                outputFile = cms.string("muonL2tree.root")
                                   )


process.p = cms.Path(process.runL2seed)




