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
#process.load('Configuration/StandardSequences/GeometryExtended_cff')
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
    input = cms.untracked.int32(1000)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
#    sys.argv[2]
    #'/store/data/BeamCommissioning09/MinimumBias/RECO/Jan29ReReco-v2/0016/204C63BA-280D-DF11-84E2-001A92811746.root'
    #'/store/data/Commissioning10/MinimumBias/RECO/v7/000/132/442/66C60D63-273C-DF11-B790-0030487CD718.root'
    #'/store/mc/Spring10/MinBias/GEN-SIM-RECO/START3X_V25B_356ReReco-v1/0004/FEE68846-F73B-DF11-98C8-002618943959.root'
    #'/store/data/Run2010A/AlCaP0/ALCARECO/394preproduction_ALCARECOEcalCalPi0Calib-v1/0020/7EF3940B-09F9-DF11-9424-0030487D8633.root'
    #'/store/data/Run2011A/MinimumBias/RECO/May10ReReco-v2/10003/FE8CCB04-4389-E011-B563-003048678FF8.root'
    ##'/store/data/Run2011A/MinimumBias/RECO/May10ReReco-v2/10003//FEB522FB-D98A-E011-BD6D-001A92971B0C.root'
    #'/store/data/Run2011A/AlCaP0/ALCARECO/ALCARECOEcalCalPi0Calib-May10ReReco-v3/0002/6AEE8D69-6B80-E011-833B-00261894396D.root'
#    '/store/data/Run2011A/AlCaP0/RAW/v1/000/163/869/EEBEE372-8B75-E011-9780-001D09F24D4E.root'
    '/store/data/Run2011B/AlCaP0/RAW/v1/000/176/868/4E8EEAE6-8FE4-E011-AD19-BCAEC5364CFB.root'

    #'/store/data/Commissioning10/AlCaP0/RAW/v4/000/135/575/42FE739C-CB61-DF11-BC7D-000423D99F1E.root'
#    '/store/express/Commissioning10/ExpressPhysics/FEVT/v4/000/130/910/FE71291D-662F-DF11-B7A8-001D09F24493.root'
    ###'/store/data/BeamCommissioning09/MinimumBias/RECO/Feb9ReReco_v2/0026/8CF1AC61-0916-DF11-8760-003048678B94.root'
#'file:/data/withvertex732.root'
#'/store/express/BeamCommissioning09/ExpressPhysics/FEVT/v2/000/123/615/38379AF1-B4E2-DE11-BB10-001617C3B706.root'
#'rfio:/castor.cern.ch/cms/store/data/BeamCommissioning09/castor/MinimumBias/RAW/v1/000/122/314/CC89C4BC-DE11-B365-0030487D0D3A.root'
    )
)

#process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT")

# Other statements

#process.GlobalTag.globaltag = 'GR10_E_V3::All'
#process.GlobalTag.globaltag = 'GR_R_42_V21::All'
process.GlobalTag.globaltag = 'GR_P_V22::All'

process.GlobalTag.toGet = cms.VPSet(
    cms.PSet(record = cms.string("EcalLaserAPDPNRatiosRcd"),
             #tag = cms.string("EcalLaserAPDPNRatios_2011fit_noVPT_nolim_online"),
             tag = cms.string("EcalLaserAPDPNRatios_data_20111010_158851_177954"),
             connect = cms.untracked.string("frontier://FrontierPrep/CMS_COND_ECAL")
             )
    )


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
                                      
                                      #outputFile = cms.string('recoAnalyzer.clusters.AlCaP0Run2011B-v1RAWdata_20111010_158851_177954.GR_R_42_V14A.laser2011fit_noVPT_nolim_online_run169985.root'),
                                      outputFile = cms.string('recoAnalyzer.clusters.AlCaP0Run2011B-v1RAWdata_20111010_158851_177954.GR_R_42_V14A_run169985.root'),

                                      dataformat = cms.int32(10), ### 10 data   ### 0 MC
                                      L1GtTmLInputTag = cms.InputTag("l1GtTriggerMenuLite"),
                                      
                                      alcaL1trigNames = cms.vstring("L1_DoubleEG10","L1_DoubleEG2_FwdVeto","L1_DoubleEG3","L1_DoubleEG5","L1_DoubleEG5_HTT50","L1_DoubleEG5_HTT75","L1_DoubleEG8","L1_DoubleEG_12_5","L1_DoubleForJet20_EtaOpp","L1_DoubleForJet36_EtaOpp","L1_DoubleIsoEG10","L1_DoubleIsoEG5","L1_DoubleIsoEG8","L1_DoubleJet36_Central","L1_DoubleJet52","L1_EG5_HTT100","L1_EG5_HTT125","L1_EG5_HTT75","L1_SingleEG12","L1_SingleEG12_Eta2p17","L1_SingleEG15","L1_SingleEG20","L1_SingleEG30","L1_SingleEG5","L1_SingleIsoEG10","L1_SingleIsoEG12","L1_SingleIsoEG12_Eta2p17","L1_SingleIsoEG15","L1_SingleJet128","L1_SingleJet16","L1_SingleJet36","L1_SingleJet36_FwdVeto","L1_SingleJet52","L1_SingleJet68","L1_SingleJet80_Central","L1_SingleJet92","L1_TripleEG5","L1_TripleEG7","L1_TripleJet28"),
                                          fullRECO = cms.int32(0),

                                      barrelHits = cms.InputTag('ecalPi0Corrected','pi0EcalRecHitsEB'),
                                      endcapHits = cms.InputTag('ecalPi0Corrected','pi0EcalRecHitsEE'),
                                      preshRecHitProducer = cms.InputTag( 'hltAlCaPi0RecHitsFilter','pi0EcalRecHitsES' ),
                                      
                                      barrelHitsEta = cms.InputTag( 'ecalEtaCorrected','etaEcalRecHitsEB' ),
                                      endcapHitsEta = cms.InputTag( 'ecalEtaCorrected','etaEcalRecHitsEE' ),
                                      preshRecHitProducerEta = cms.InputTag( 'hltAlCaEtaRecHitsFilter','etaEcalRecHitsES' ),
                                      
                                      saveEGObj = cms.bool (False),
                                      saveAllRecHitEB = cms.bool(False),
                                      saveAllRecHitEE = cms.bool(False),
                                      
                                      clusSeedThr = cms.double( 0.5 ),
                                      clusSeedThrEndCap = cms.double( 1.0),
                                      clusEtaSize = cms.int32( 3 ),
                                      clusPhiSize = cms.int32( 3 ),
                                      
                                      L1GtRecordInputTag = cms.InputTag('hltGtDigis'),
                                      
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
                                          selePtGamma = cms.double( 1.3 ),
                                          selePtPi0 = cms.double( 2.6 ),
                                          seleMinvMaxPi0 = cms.double( 0.23 ),
                                          seleMinvMinPi0 = cms.double( 0.04 ),
                                          seleS4S9Gamma = cms.double( 0.83 ),
                                          selePi0BeltDR = cms.double( 0.2 ),
                                          selePi0BeltDeta = cms.double( 0.05 ),
                                          ptMinForIsolation = cms.double(0.5),
                                          selePi0Iso = cms.double( 0.5 ),
                                      
                                          doSelForPi0Endcap = cms.bool( True ),
                                          region1_Pi0EndCap = cms.double( 2 ),
                                          selePtGammaPi0EndCap_region1 = cms.double( 0.6 ),
                                          selePtPi0EndCap_region1 = cms.double( 2.5 ),
                                          region2_Pi0EndCap = cms.double( 2.5 ),
                                          selePtGammaPi0EndCap_region2 = cms.double( 0.6 ),
                                          selePtPi0EndCap_region2 = cms.double( 2.5 ),
                                          selePtGammaPi0EndCap_region3 = cms.double( 0.5 ),
                                          selePtPi0EndCap_region3 = cms.double( 1.0 ),
                                          selePtPi0MaxEndCap_region3 = cms.double(2.5 ),
                                          seleS4S9GammaEndCap = cms.double( 0.90 ),
                                          seleMinvMaxPi0EndCap = cms.double( 0.3 ),
                                          seleMinvMinPi0EndCap = cms.double( 0.05 ),
                                          selePi0BeltDREndCap = cms.double( 0.2 ),
                                          selePi0BeltDetaEndCap = cms.double( 0.05 ),
                                          ptMinForIsolationEndCap = cms.double(0.5),
                                          selePi0IsoEndCap = cms.double( 0.5 ),


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

process.beamspot = cms.Path(process.offlineBeamSpot)
process.out_step = cms.EndPath(process.recoAnalyzer )
process.schedule = cms.Schedule(process.beamspot,process.pathALCARECOEcalCalPi0Calib,process.pathALCARECOEcalCalEtaCalib,process.out_step)

#process.schedule = cms.Schedule(process.raw2digi_step,process.L1GtObjectMap_step,process.L1Reco_step,process.recoEcal,process.out_step)
