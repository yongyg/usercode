# Auto generated configuration file
# using: 
# Revision: 1.149 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: promptCollisionReco -s RAW2DIGI,L1Reco,RECO,DQM,ALCA:SiStripCalZeroBias --datatier RECO --eventcontent RECO --conditions CRAFT09_R_V4::All --scenario pp --no_exec --data --magField AutoFromDBCurrent -n 100
import FWCore.ParameterSet.Config as cms

process = cms.Process('RERECO')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 10000 # print one event every

process.options = cms.untracked.PSet(
  wantSummary = cms.untracked.bool(True)
)
process.load('Configuration/StandardSequences/GeometryExtended_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
process.load('Configuration/StandardSequences/L1Reco_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('DQMOffline/Configuration/DQMOffline_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/data/Run2010B/AlCaP0/ALCARECO/preproduction_3110_ALCARECOEcalCalEtaCalib-v1/0062/8C9AA317-CC33-E011-994B-002481E0DC7C.root'
    )
)


#process.GlobalTag.globaltag = 'GR10_E_V3::All'
process.GlobalTag.globaltag = 'GR_R_311_V1A::All'


process.recoAnalyzer = cms.EDAnalyzer("RecoAnalyzer",
                                      outputFile = cms.string('recoAnalyzer.clusters.AlCaP0Run2010B-preproduction_3110_ALCARECOEcalCalEtaCalib-v1ALCARECO.root'),
                                      
                                      dataformat = cms.int32(10), ### 10 data   ### 0 MC
                                      L1GtTmLInputTag = cms.InputTag("l1GtTriggerMenuLite"),
                                      
                                      alcaL1trigNames = cms.vstring("L1_SingleIsoEG5","L1_SingleIsoEG8","L1_SingleIsoEG10","L1_SingleIsoEG12","L1_SingleIsoEG15","L1_SingleEG2","L1_SingleEG5","L1_SingleEG8","L1_SingleEG10","L1_SingleEG12","L1_SingleEG15","L1_SingleEG20","L1_SingleJet6U","L1_SingleJet10U","L1_SingleJet20U","L1_SingleJet30U","L1_SingleJet40U","L1_SingleJet50U","L1_DoubleJet30U","L1_DoubleEG5","L1_DoubleEG2"),

                                      barrelHits = cms.InputTag('ecalPi0Corrected','pi0EcalRecHitsEB'),
                                      endcapHits = cms.InputTag('ecalPi0Corrected','pi0EcalRecHitsEE'),
                                      preshRecHitProducer = cms.InputTag( 'hltAlCaPi0RecHitsFilter','pi0EcalRecHitsES' ),
                                      
                                      barrelHitsEta = cms.InputTag( 'ecalEtaCorrected','etaEcalRecHitsEB' ),
                                      endcapHitsEta = cms.InputTag( 'ecalEtaCorrected','etaEcalRecHitsEE' ),
                                      preshRecHitProducerEta = cms.InputTag( 'hltAlCaEtaRecHitsFilter','etaEcalRecHitsES' ),
 
                                      clusSeedThr = cms.double( 0.5 ),
                                      clusSeedThrEndCap = cms.double( 1.0),
                                      clusEtaSize = cms.int32( 3 ),
                                      clusPhiSize = cms.int32( 3 ),
                                      
                                      L1GtRecordInputTag = cms.InputTag('hltGtDigis'),
                                      TriggerResultsTag = cms.InputTag('TriggerResults','','HLT'),
                                      
                                      doSelForPi0Barrel = cms.bool( False ),
                                      selePtGamma = cms.double( 1.3 ),
                                      selePtPi0 = cms.double( 2.6 ),
                                      seleMinvMaxPi0 = cms.double( 0.23 ),
                                      seleMinvMinPi0 = cms.double( 0.04 ),
                                      seleS4S9Gamma = cms.double( 0.83 ),
                                      selePi0BeltDR = cms.double( 0.2 ),
                                      selePi0BeltDeta = cms.double( 0.05 ),
                                      ptMinForIsolation = cms.double(0.5),
                                      selePi0Iso = cms.double( 0.5 ),
                                      
                                      doSelForPi0Endcap = cms.bool( False ),
                                      region1_Pi0EndCap = cms.double( 2 ),
                                      selePtGammaPi0EndCap_region1 = cms.double( 0.6 ),
                                      selePtPi0EndCap_region1 = cms.double( 2.5 ),
                                      region2_Pi0EndCap = cms.double( 2.5 ),
                                      selePtGammaPi0EndCap_region2 = cms.double( 0.6 ),
                                      selePtPi0EndCap_region2 = cms.double( 2.5 ),
                                      selePtGammaPi0EndCap_region3 = cms.double( 0.5 ),
                                      selePtPi0EndCap_region3 = cms.double( 1.0 ),
                                      seleS4S9GammaEndCap = cms.double( 0.90 ),
                                      seleMinvMaxPi0EndCap = cms.double( 0.3 ),
                                      seleMinvMinPi0EndCap = cms.double( 0.05 ),
                                      selePi0BeltDREndCap = cms.double( 0.2 ),
                                      selePi0BeltDetaEndCap = cms.double( 0.05 ),
                                      ptMinForIsolationEndCap = cms.double(0.5),
                                      selePi0IsoEndCap = cms.double( 0.5 ),
                                      
                                      removeSpike = cms.untracked.int32(1),
                                      doSelForEtaBarrel = cms.bool( True ),
                                      removePi0CandidatesForEta = cms.bool( True ),
                                      massLowPi0Cand = cms.double( 0.084 ),
                                      massHighPi0Cand = cms.double( 0.156 ),
                                      selePtGammaEta = cms.double( 1.2 ),
                                      selePtEta = cms.double( 4.0 ),
                                      seleMinvMaxEta = cms.double( 0.8 ),
                                      seleMinvMinEta = cms.double( 0.3 ),
                                      seleS4S9GammaEta = cms.double( 0.87 ),
                                      seleS9S25GammaEta = cms.double( 0.8 ),
                                      seleEtaBeltDR = cms.double( 0.3 ),
                                      seleEtaBeltDeta = cms.double( 0.1 ),
                                      ptMinForIsolationEta = cms.double(0.5),
                                      seleEtaIso = cms.double( 0.5 ),
                                      
                                      
                                      doSelForEtaEndcap = cms.bool( True ),
                                      region1_EtaEndCap = cms.double( 2 ),
                                      selePtGammaEtaEndCap_region1 = cms.double( 1.0),
                                      selePtEtaEndCap_region1 = cms.double( 3.0 ),
                                      region2_EtaEndCap = cms.double( 2.5 ),
                                      selePtGammaEtaEndCap_region2 = cms.double( 1.0 ),
                                      selePtEtaEndCap_region2 = cms.double( 3.0 ),
                                      selePtGammaEtaEndCap_region3 = cms.double( 0.7 ),
                                      selePtEtaEndCap_region3 = cms.double( 3.0 ),
                                      seleS4S9GammaEtaEndCap = cms.double( 0.90 ),
                                      seleS9S25GammaEtaEndCap = cms.double( 0.85 ),
                                      seleMinvMaxEtaEndCap = cms.double( 0.9 ),
                                      seleMinvMinEtaEndCap = cms.double( 0.2 ),
                                      seleEtaBeltDREndCap = cms.double( 0.3 ),
                                      seleEtaBeltDetaEndCap = cms.double( 0.1 ),
                                      ptMinForIsolationEtaEndCap = cms.double(0.5),
                                      seleEtaIsoEndCap = cms.double( 0.5 ),
                                      
                                      
                                      storeRecHitES = cms.bool( True ),
                                      preshNclust = cms.int32( 4 ),
                                      preshClusterEnergyCut = cms.double( 0.0 ),
                                      preshStripEnergyCut = cms.double( 0.0 ),
                                      preshSeededNstrip = cms.int32( 15 ),
                                      preshCalibPlaneX = cms.double( 1.0 ),
                                      preshCalibPlaneY = cms.double( 0.7 ),
                                      preshCalibGamma = cms.double( 0.024 ),
                                      preshCalibMIP = cms.double( 9.0E-5 ),
                                      debugLevelES = cms.string( " " ),
                                      posCalcParameters = cms.PSet( T0_barl      = cms.double(7.4),
                                                                    T0_endc      = cms.double(3.1),
                                                                    T0_endcPresh = cms.double(1.2),
                                                                    LogWeighted  = cms.bool(True),
                                                                    W0           = cms.double(4.2),
                                                                    X0           = cms.double(0.89)
                                                                    )          ,
                                      debugLevel = cms.int32(0)
                                      
                                      )

process.out_step = cms.EndPath(process.recoAnalyzer )

process.schedule = cms.Schedule(process.out_step)
