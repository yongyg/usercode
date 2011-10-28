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
process.load('Configuration/StandardSequences/GeometryDB_cff')
process.load('Configuration/StandardSequences/MagneticField_AutoFromDBCurrent_cff')
process.load('Configuration/StandardSequences/RawToDigi_Data_cff')
process.load('Configuration/StandardSequences/L1Reco_cff')
process.load('Configuration/StandardSequences/Reconstruction_cff')
process.load('DQMOffline/Configuration/DQMOffline_cff')
process.load('Configuration/StandardSequences/EndOfProcess_cff')
process.load('Configuration/StandardSequences/FrontierConditions_GlobalTag_cff')
process.load('Configuration/EventContent/EventContent_cff')
process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(5000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/data/Run2011A/AlCaP0/RAW/v1/000/163/869/EEBEE372-8B75-E011-9780-001D09F24D4E.root'
    '/store/mc/Summer11/ZJetToMuMu_Pt-15to3000_TuneD6T_Flat_7TeV_pythia6/GEN-SIM-RECO/PU_S3_START42_V11-v2/0000/82F5509B-4281-E011-9D46-00215E21DA2C.root'
    )
)

#process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT")

# Other statements

#process.GlobalTag.globaltag = 'GR10_E_V3::All'
#process.GlobalTag.globaltag = 'GR_R_42_V14A::All'
process.GlobalTag.globaltag = 'START42_V12::All'

#process.GlobalTag.toGet = cms.VPSet(
#    cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
             #tag = cms.string("EcalLaserAPDPNRatios_2011fit_noVPT_nolim_online"),
#             tag = cms.string("EcalLaserAPDPNRatios_2011V3_online"),
#             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
##             )
#    )


##


##from IOMC.RandomEngine.IOMC_cff import *
from SimCalorimetry.EcalSimProducers.ecalSimParameterMap_cff import *
from SimCalorimetry.EcalSimProducers.ecalElectronicsSim_cff import *
from SimCalorimetry.EcalSimProducers.ecalNotContainmentSim_cff import *


process.recoAnalyzer = cms.EDAnalyzer("RecoAnalyzer",
                                          ecal_electronics_sim,
                                          ecal_sim_parameter_map,
                                          ecal_notCont_sim,
                                          hitsProducer = cms.string('g4SimHits'),
                                      
                                      #outputFile = cms.string('recoAnalyzer.clusters.AlCaP0Run2011A-v1RAW2011V3_online.GR_R_42_V14A.laser2011fit_noVPT_nolim_online_run163476.root'),
                                      outputFile = cms.string('pi0Analyzer.QCD_Pt-15to30_TuneZ2_7TeV_pythia6Summer11-PU_S3_START42_V11-v2GEN-SIM-RECO.root'),

                                      dataformat = cms.int32(0), ### 10 data   ### 0 MC
                                      L1GtTmLInputTag = cms.InputTag("l1GtTriggerMenuLite"),
                                      
                                      alcaL1trigNames = cms.vstring("L1_DoubleEG10","L1_DoubleEG2_FwdVeto","L1_DoubleEG3","L1_DoubleEG5","L1_DoubleEG5_HTT50","L1_DoubleEG5_HTT75","L1_DoubleEG8","L1_DoubleEG_12_5","L1_DoubleForJet20_EtaOpp","L1_DoubleForJet36_EtaOpp","L1_DoubleIsoEG10","L1_DoubleIsoEG5","L1_DoubleIsoEG8","L1_DoubleJet36_Central","L1_DoubleJet52","L1_EG5_HTT100","L1_EG5_HTT125","L1_EG5_HTT75","L1_SingleEG12","L1_SingleEG12_Eta2p17","L1_SingleEG15","L1_SingleEG20","L1_SingleEG30","L1_SingleEG5","L1_SingleIsoEG10","L1_SingleIsoEG12","L1_SingleIsoEG12_Eta2p17","L1_SingleIsoEG15","L1_SingleJet128","L1_SingleJet16","L1_SingleJet36","L1_SingleJet36_FwdVeto","L1_SingleJet52","L1_SingleJet68","L1_SingleJet80_Central","L1_SingleJet92","L1_TripleEG5","L1_TripleEG7","L1_TripleJet28"),
                                          fullRECO = cms.int32(0),

                                      #barrelHits = cms.InputTag('ecalPi0Corrected','pi0EcalRecHitsEB'),
                                      barrelHits = cms.InputTag('ecalRecHit','EcalRecHitsEB'),

                                      #endcapHits = cms.InputTag('ecalPi0Corrected','pi0EcalRecHitsEE'),
                                      
                                      endcapHits = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
                                      #preshRecHitProducer = cms.InputTag( 'hltAlCaPi0RecHitsFilter','pi0EcalRecHitsES' ),
                                      
                                      preshRecHitProducer = cms.InputTag('ecalPreshowerRecHit','EcalRecHitsES'),


                                      barrelHitsEta = cms.InputTag('ecalRecHit','EcalRecHitsEB'),
                                      endcapHitsEta = cms.InputTag('ecalRecHit','EcalRecHitsEE'),
                                      preshRecHitProducerEta = cms.InputTag('ecalPreshowerRecHit','EcalRecHitsES'),
                                                                            

                                      #barrelHitsEta = cms.InputTag( 'ecalEtaCorrected','etaEcalRecHitsEB' ),
                                      #endcapHitsEta = cms.InputTag( 'ecalEtaCorrected','etaEcalRecHitsEE' ),
                                      #preshRecHitProducerEta = cms.InputTag( 'hltAlCaEtaRecHitsFilter','etaEcalRecHitsES' ),
                                      
                                      saveEGObj = cms.bool (False),
                                      saveAllRecHitEB = cms.bool(False),
                                      saveAllRecHitEE = cms.bool(False),
                                      
                                      clusSeedThr = cms.double( 0.5 ),
                                      clusSeedThrEndCap = cms.double( 1.0),
                                      clusEtaSize = cms.int32( 3 ),
                                      clusPhiSize = cms.int32( 3 ),
                                      
                                      #L1GtRecordInputTag = cms.InputTag('hltGtDigis'),
                                      L1GtRecordInputTag = cms.InputTag('gtDigis'),
                                      
                                      ecalDigiReFit = cms.bool(False),
                                      EBdigiCollection = cms.InputTag("ecalDigis","ebDigis"),
                                      EEdigiCollection = cms.InputTag("ecalDigis","eeDigis"),

                                      TriggerResultsTag = cms.InputTag('TriggerResults','','HLT'),
                                      HLTPaths = cms.vstring('HLT_PhysicsDeclared'),               # provide the list of HLT paths (or patterns) you want
                                      HLTPathsPrescales = cms.vuint32(1),      # provide the list of prescales correseponding to the paths and patterns
                                      HLTOverallPrescale = cms.uint32(1),     # privide the overall prescale used on top of the final result
                                      eventSetupPathsKey = cms.string(''),    # not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
                                      andOr = cms.bool(True),                 # how to deal with multiple triggers: True (OR) accept if ANY is true, False (AND) accept if ALL are true
                                      throw = cms.bool(True),                  # throw exception on unknown path names
                                      
                                      beamSpotInputTag = cms.InputTag('offlineBeamSpot'),
                                      
                                      ebDigiCollectionv1 = cms.InputTag("simEcalUnsuppressedDigis","myebdigiv1"),
                                      reRunPixelRecHits = cms.untracked.bool(True),
                                      
                                      seleXtalMinEnergy  = cms.double(-10),
                                      seleXtalMinEnergyEndCap  = cms.double(-10),
                                      
                                      maxNumberofSeeds = cms.int32(70000),
                                      maxNumberofClusters = cms.int32(2000),
                                      useRecoFlag = cms.bool( False ),
                                      flagLevelRecHitsToUse = cms.int32( 1 ),
                                      useDBStatus = cms.bool( True ),
                                      statusLevelRecHitsToUse = cms.int32( 1 ),
                                      
                                      RegionalMatch = cms.bool(False),
                                      ptMinEMObj = cms.double(2.0),
                                      EMregionEtaMargin = cms.double(0.25),
                                      EMregionPhiMargin = cms.double(0.4),
                                      l1IsolatedTag = cms.InputTag("l1extraPi0s","Isolated"),
                                      l1NonIsolatedTag = cms.InputTag("l1extraPi0s","NonIsolated"),
                                      
                                      nMinRecHitsSel1stCluster = cms.int32( 1 ),
                                      nMinRecHitsSel2ndCluster = cms.int32( 1 ),
                                      
                                      saveAllPhotonBarrel = cms.bool(False),
                                      saveAllPhotonEndcap = cms.bool(False),
                                      
                                      doSelForPi0Barrel = cms.bool( True ),
                                      selePtGamma = cms.double( 1 ),
                                      selePtPi0 = cms.double( 2. ),
                                      seleMinvMaxPi0 = cms.double( 0.3 ),
                                      seleMinvMinPi0 = cms.double( 0.0 ),
                                      seleS4S9Gamma = cms.double( 0.8 ),
                                      selePi0BeltDR = cms.double( 0.2 ),
                                      selePi0BeltDeta = cms.double( 0.05 ),
                                      ptMinForIsolation = cms.double(0.5),
                                      selePi0Iso = cms.double( 999 ),
                                      
                                      doSelForPi0Endcap = cms.bool( True ),
                                      region1_Pi0EndCap = cms.double( 2 ),
                                      selePtGammaPi0EndCap_region1 = cms.double( 0.5 ),
                                      selePtPi0EndCap_region1 = cms.double( 1.0 ),
                                      region2_Pi0EndCap = cms.double( 2.5 ),
                                      selePtGammaPi0EndCap_region2 = cms.double( 0.5 ),
                                      selePtPi0EndCap_region2 = cms.double( 1.0 ),
                                      selePtGammaPi0EndCap_region3 = cms.double( 0.5 ),
                                      selePtPi0EndCap_region3 = cms.double( 1.0 ),
                                      selePtPi0MaxEndCap_region3 = cms.double(9999 ),
                                      seleS4S9GammaEndCap = cms.double( 0.80 ),
                                      seleMinvMaxPi0EndCap = cms.double( 0.3 ),
                                      seleMinvMinPi0EndCap = cms.double( 0.0 ),
                                      selePi0BeltDREndCap = cms.double( 0.2 ),
                                      selePi0BeltDetaEndCap = cms.double( 0.05 ),
                                      ptMinForIsolationEndCap = cms.double(0.5),
                                      selePi0IsoEndCap = cms.double( 999 ),
                                      
                                      
                                      doSelForEtaBarrel = cms.bool( True ),
                                      removePi0CandidatesForEta = cms.bool( True ),
                                      massLowPi0Cand = cms.double( 0.084 ),
                                      massHighPi0Cand = cms.double( 0.156 ),
                                      selePtGammaEta = cms.double( 1 ),
                                      selePtEta = cms.double( 3.0 ),
                                      seleMinvMaxEta = cms.double( 0.8 ),
                                      seleMinvMinEta = cms.double( 0.3 ),
                                      seleS4S9GammaEta = cms.double( 0.8 ),
                                      seleS9S25GammaEta = cms.double( 0.7 ),
                                      seleEtaBeltDR = cms.double( 0.3 ),
                                      seleEtaBeltDeta = cms.double( 0.1 ),
                                      ptMinForIsolationEta = cms.double(0.5),
                                      seleEtaIso = cms.double( 999 ),
                                                                            
                                      doSelForEtaEndcap = cms.bool( True ),
                                      region1_EtaEndCap = cms.double( 2 ),
                                      selePtGammaEtaEndCap_region1 = cms.double( 0.7),
                                      selePtEtaEndCap_region1 = cms.double( 2.0 ),
                                      region2_EtaEndCap = cms.double( 2.5 ),
                                      selePtGammaEtaEndCap_region2 = cms.double( 0.7 ),
                                      selePtEtaEndCap_region2 = cms.double( 2.0 ),
                                      selePtGammaEtaEndCap_region3 = cms.double( 0.7 ),
                                      selePtEtaEndCap_region3 = cms.double( 2.0 ),
                                      seleS4S9GammaEtaEndCap = cms.double( 0.50 ),
                                      seleS9S25GammaEtaEndCap = cms.double( 0.50 ),
                                      seleMinvMaxEtaEndCap = cms.double( 0.9 ),
                                      seleMinvMinEtaEndCap = cms.double( 0.2 ),
                                      seleEtaBeltDREndCap = cms.double( 0.3 ),
                                      seleEtaBeltDetaEndCap = cms.double( 0.1 ),
                                      ptMinForIsolationEtaEndCap = cms.double(0.5),
                                      seleEtaIsoEndCap = cms.double(999 ),
                                      
                                      
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
#process.schedule = cms.Schedule(process.pathALCARECOEcalCalPi0Calib,process.pathALCARECOEcalCalEtaCalib,process.out_step)
process.schedule = cms.Schedule(process.out_step)

#process.schedule = cms.Schedule(process.raw2digi_step,process.L1GtObjectMap_step,process.L1Reco_step,process.recoEcal,process.out_step)
