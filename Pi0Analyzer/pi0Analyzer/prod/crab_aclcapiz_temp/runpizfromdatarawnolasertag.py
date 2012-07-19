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

#process.load('Configuration.StandardSequences.AlCaRecoStreams_cff')

import HLTrigger.HLTfilters.hltHighLevel_cfi
import RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi


process.ecalpi0ebCalibHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    HLTPaths = ['AlCa_EcalPi0EBonly_v*'],
    #eventSetupPathsKey='EcalCalPi0CalibEB',
    throw = False
    )

process.ecalPi0CorrectedEB =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
                doEnergyScale = cms.bool(True),
                doIntercalib = cms.bool(True),
                EERecHitCollection = cms.InputTag("",""),
                EBRecHitCollection = cms.InputTag("hltAlCaPi0EBUncalibrator","pi0EcalRecHitsEB"),
                doLaserCorrections = cms.bool(True),
                EBRecalibRecHitCollection = cms.string('pi0EcalRecHitsEB'),
                EERecalibRecHitCollection = cms.string('pi0EcalRecHitsEE')
                )

process.ecalpi0eeCalibHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
        # HLTPaths = ['AlCa_EcalPi0','AlCa_EcalEta'],
    HLTPaths = ['AlCa_EcalPi0EEonly_v*'],
    #eventSetupPathsKey='EcalCalPi0CalibEE',
    throw = False
            )

process.ecalPi0CorrectedEE =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
                doEnergyScale = cms.bool(True),
                doIntercalib = cms.bool(True),
                EBRecHitCollection = cms.InputTag("",""),
                EERecHitCollection = cms.InputTag("hltAlCaPi0EEUncalibrator","pi0EcalRecHitsEE"),
                doLaserCorrections = cms.bool(True),
                EBRecalibRecHitCollection = cms.string('pi0EcalRecHitsEB'),
                EERecalibRecHitCollection = cms.string('pi0EcalRecHitsEE')
                )

process.ecalEtaCorrectedEB =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
                doEnergyScale = cms.bool(True),
                doIntercalib = cms.bool(True),
                EERecHitCollection = cms.InputTag("",""),
                EBRecHitCollection = cms.InputTag("hltAlCaEtaEBUncalibrator","etaEcalRecHitsEB"),
                doLaserCorrections = cms.bool(True),
                EBRecalibRecHitCollection = cms.string('etaEcalRecHitsEB'),
                EERecalibRecHitCollection = cms.string('etaEcalRecHitsEE')
                )


process.ecalEtaCorrectedEE =  RecoLocalCalo.EcalRecProducers.ecalRecalibRecHit_cfi.ecalRecHit.clone(
                doEnergyScale = cms.bool(True),
                doIntercalib = cms.bool(True),
                EBRecHitCollection = cms.InputTag("",""),
                EERecHitCollection = cms.InputTag("hltAlCaEtaEEUncalibrator","etaEcalRecHitsEE"),
                doLaserCorrections = cms.bool(True),
                EBRecalibRecHitCollection = cms.string('etaEcalRecHitsEB'),
                EERecalibRecHitCollection = cms.string('etaEcalRecHitsEE')
                )


process.ecaletaebCalibHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    # HLTPaths = ['AlCa_EcalPi0','AlCa_EcalEta'],
    HLTPaths = ['AlCa_EcalEtaEBonly_v*'],
    #eventSetupPathsKey='EcalCalEtaCalibEB',
    throw = False
    )

process.ecaletaeeCalibHLT = HLTrigger.HLTfilters.hltHighLevel_cfi.hltHighLevel.clone(
    # HLTPaths = ['AlCa_EcalPi0','AlCa_EcalEta'],
    #HLTPaths = ['AlCa_EcalEtaEEonly_v3'],
    HLTPaths = ['AlCa_EcalEtaEEonly_v*'],
    #TriggerResultsTag = cms.InputTag("TriggerResults","","TEST"),
    #eventSetupPathsKey='EcalCalEtaCalibEE',
    throw = False
    )

process.maxEvents = cms.untracked.PSet(
        input = cms.untracked.int32(200000)
        )

# Input source
process.source = cms.Source("PoolSource",
                                fileNames = cms.untracked.vstring(
        '/store/data/Run2012A/AlCaP0/RAW/v1/000/190/949/005F9643-2D84-E111-A7BE-001D09F2910A.root',
            )
                            )
#process.GlobalTag.globaltag = 'GR_P_V32::All'
process.GlobalTag.globaltag = 'GRTagName::All'


process.pi0Analyzer = cms.EDAnalyzer("Pi0Analyzer",
                                     outputFile = cms.string('pi0Analyzer.datasetname.GRTagName.runxxxxxx.root'),
                                     usegtDigis = cms.untracked.bool (True),
                                     l1GtRecordInputTag = cms.untracked.InputTag("hltGtDigis"),
                                     beamSpotInputTag = cms.untracked.InputTag("offlineBeamSpot"),
                                     doSelForPi0Barrel = cms.bool( True ),
                                     barrelHits = cms.untracked.InputTag( 'ecalPi0CorrectedEB','pi0EcalRecHitsEB' ),
                                     doSelForPi0Endcap = cms.bool( True ),
                                     endcapHits = cms.untracked.InputTag( 'ecalPi0CorrectedEE','pi0EcalRecHitsEE' ),
                                     preshRecHitProducer = cms.untracked.InputTag( 'hltAlCaPi0RecHitsFilterEEonly','pi0EcalRecHitsES' ),
                                     doSelForEtaBarrel = cms.bool( True ),
                                     barrelHitsEta = cms.untracked.InputTag( 'ecalEtaCorrectedEB','etaEcalRecHitsEB' ),
                                     doSelForEtaEndcap = cms.bool( True ),
                                     endcapHitsEta = cms.untracked.InputTag( 'ecalEtaCorrectedEE','etaEcalRecHitsEE' ),
                                     preshRecHitProducerEta = cms.untracked.InputTag( 'hltAlCaEtaRecHitsFilterEEonly','etaEcalRecHitsES' ),
                                     
                                     posCalcParameters = cms.untracked.PSet( T0_barl      = cms.double(7.4),
                                                                             T0_endc      = cms.double(3.1),
                                                                             T0_endcPresh = cms.double(1.2),
                                                                             LogWeighted  = cms.bool(True),
                                                                             W0           = cms.double(4.2),
                                                                             X0           = cms.double(0.89)
                                                                             )          ,
                                     debugLevel = cms.int32(0)
                                     )


process.beamspot = cms.Path(process.offlineBeamSpot)
process.out_step = cms.EndPath(process.pi0Analyzer )

process.pathALCARECOEcalCalPi0CalibEB = cms.Path(process.ecalpi0ebCalibHLT*process.ecalPi0CorrectedEB)
process.pathALCARECOEcalCalPi0CalibEE = cms.Path(process.ecalpi0eeCalibHLT*process.ecalPi0CorrectedEE)
process.pathALCARECOEcalCalEtaCalibEB = cms.Path(process.ecaletaebCalibHLT*process.ecalEtaCorrectedEB)
process.pathALCARECOEcalCalEtaCalibEE = cms.Path(process.ecaletaeeCalibHLT*process.ecalEtaCorrectedEE)

process.schedule = cms.Schedule(process.beamspot,process.pathALCARECOEcalCalPi0CalibEB,process.pathALCARECOEcalCalPi0CalibEE,process.pathALCARECOEcalCalEtaCalibEB,process.pathALCARECOEcalCalEtaCalibEE,process.out_step)


