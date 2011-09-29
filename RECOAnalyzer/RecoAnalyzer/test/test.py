# Auto generated configuration file
# using: 
# Revision: 1.149 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: promptCollisionReco -s RAW2DIGI,L1Reco,RECO,DQM,ALCA:SiStripCalZeroBias --datatier RECO --eventcontent RECO --conditions CRAFT09_R_V4::All --scenario pp --no_exec --data --magField AutoFromDBCurrent -n 100
import FWCore.ParameterSet.Config as cms

process = cms.Process('RecoAnalysis')

# import of standard configurations
process.load('Configuration/StandardSequences/Services_cff')
process.load('FWCore/MessageService/MessageLogger_cfi')

process.MessageLogger.cerr.FwkReport.reportEvery = 1000 # print one event every

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


process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)


# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    #'/store/mc/Fall10/QCD_Pt_1000to1400_TuneZ2_7TeV_pythia6/GEN-SIM-RECO/START38_V12-v1/0005/D85FC8E6-8DCC-DF11-9479-0024E8769B2C.root'
    #'/store/mc/Fall10/QCD_Pt_15to3000_TuneZ2_Flat_7TeV_pythia6/GEN-SIM-RECO/START38_V12-v1/0005/F24E17D7-F3CD-DF11-89A4-00215E221938.root'
    #'/store/mc/Spring10/ZprimeSSMToMuMu_M-1000_7TeV-pythia6/GEN-SIM-RECO/START3X_V26-v1/0002/90B1AD26-EF55-DF11-ACCA-001F29C9C492.root'
    #'/store/data/Run2010A/Mu/RECO/Nov4ReReco_v1/0011/82BC70AF-50EC-DF11-AB25-00237DA1FD7C.root'
    #'/store/relval/CMSSW_4_2_0/Electron/RECO/GR_R_42_V8_RelVal_wzEG2010B-v1/0055/FE29BC99-F05E-E011-8DED-001A92971B90.root'
    #'file:step2_RAW2DIGI_L1Reco_RECO_DQM_98_1_G9J.root'
    #'file:../../CMSSW_4_2_3/src/AEE6B4AB-F87F-E011-BC0F-485B39800BB1.GluGluToHToGG_M-115_7TeV-powheg-pythia6AODSIMPU_S3_START42_V11-v2.root'

    "/store/data/Run2011A/Photon/AOD/05Jul2011ReReco-ECAL-v1/0001/526B9884-27A8-E011-BD16-003048673E86.root"
    #"/store/data/Run2011A/Photon/AOD/05Jul2011ReReco-ECAL-v1/0001/A068ED1E-24A8-E011-BF4A-00E08178C0E7.root",
    #"/store/data/Run2011A/Photon/AOD/05Jul2011ReReco-ECAL-v1/0001/06E7A6B6-25A8-E011-A5C2-0025B3E05DFE.root",
    #"/store/data/Run2011A/Photon/AOD/05Jul2011ReReco-ECAL-v1/0001/3A1B910D-4EA9-E011-AA51-001A64789470.root"

    #'/store/data/Run2010B/Photon/RECO/Dec22ReReco_v1/0024/FA52A8DD-E00F-E011-AA1B-003048678BB8.root'
    #'/store/data/Run2011A/Photon/RECO/May10ReReco-v1/0005/029FBE15-1D7C-E011-BC71-001A64789E44.root'

    #'/store/relval/CMSSW_4_2_3/RelValTTbar/GEN-SIM-RECO/START42_V12-v2/0062/728877FF-717B-E011-9989-00261894395B.root'
    ) 
#  skipEvents=cms.untracked.uint32(1)

)

#process.source.inputCommands = cms.untracked.vstring("keep *", "drop *_MEtoEDMConverter_*_*", "drop L1GlobalTriggerObjectMapRecord_hltL1GtObjectMap__HLT")

# Other statements


#process.GlobalTag.globaltag = 'GR10_P_V6::All'
#process.GlobalTag.globaltag = 'GR_R_38X_V9::All'
##GR_R_38X_V9
process.GlobalTag.globaltag = 'FT_R_42_V17A::All'
#process.GlobalTag.globaltag = 'START38_V12::All'
#process.GlobalTag.globaltag = 'START3X_V26::All'
#process.GlobalTag.globaltag = 'FT_R_38X_V14A::All'
#process.GlobalTag.globaltag = 'GR_R_42_V12::All'
#process.GlobalTag.globaltag = 'FT_R_42_V13A::All'

# KH, adding lazy-download

###this is not Good all, dead runn
#process.AdaptorConfig = cms.Service("AdaptorConfig",
#                                                                        stats = cms.untracked.bool(True),
#                                                                        enable = cms.untracked.bool(True),
#                                                                        #                                    tempDir = cms.untracked.string(""),
#                                                                        cacheHint = cms.untracked.string("lazy-download"),
#                                                                        readHint = cms.untracked.string("auto-detect"))


#####################################################################################################
####
####  Top level replaces for handling strange scenarios of early collisions
####

## TRACKING:
## Skip events with HV off


###
###  end of top level replacements
###
###############################################################################################


#process.out_step = cms.EndPath(process.FEVT)

process.hltTrigReport = cms.EDAnalyzer( 'HLTrigReport',
     HLTriggerResults = cms.InputTag( 'TriggerResults','','HLT' )
)

## Geometry and Detector Conditions (needed for a few patTuple production steps)
#process.GlobalTag.globaltag = cms.string( 'GR_R_42_V12::All' )
##-------------------- Import the JEC services -----------------------

#process.load('JetMETCorrections/Configuration/DefaultJEC_cff')
##-------------------- Import the Jet RECO modules -----------------------

#process.load('RecoJets/Configuration/RecoPFJets_cff')
#process.kt6PFJets.doRhoFastjet = True
#process.kt6PFJets.Rho_EtaMax = cms.double(2.5)


##-------------------- Turn-on the FastJet jet area calculation for your favorite algorithm -----------------------
#process.ak5PFJets.doAreaFastjet = True


from RecoJets.Configuration.RecoPFJets_cff import *
process.kt6PFJetsForRhoComputationEtaMaxDefault = kt6PFJets.clone(doRhoFastjet = True)

process.kt6PFJetsForRhoComputationEtaMax25 = kt6PFJets.clone(Rho_EtaMax = 2.5,doRhoFastjet = True)

##.............................................................................................
##-------------------- User analyzer  --------------------------------
#process.MyAnalyzer  = cms.EDAnalyzer('MyAnalyzer',
##...............................................................................................
#JetCorrectionService = cms.string('ak5PFL1FastL2L3')
##...............................................................................................
##-------------------- Include in the path the jet reconstruction  --------------------------------
#process.p = cms.Path(process.kt6PFJets * process.ak5PFJets * process.MyAnalyzer)
#)

#process.load('RecoJets.JetProducers.kt4PFJets_cfi')
#process.kt6PFJets = process.kt4PFJets.clone( rParam = 0.6, doRhoFastjet = True )
#process.kt6PFJets.Rho_EtaMax = cms.double(2.5)

process.recoAnalyzer = cms.EDAnalyzer("RecoAnalyzer",
                                      outputFile = cms.string('analysis.PhotonRun2011A-05Jul2011ReReco-ECAL-v1AOD.root'),
                                      ##L1GtRecordInputTag = cms.InputTag('gtDigis'),

                                      trigResultTag = cms.untracked.InputTag('TriggerResults','','HLT'),
                                      trigSummaryTag =  cms.untracked.InputTag('hltTriggerSummaryAOD','','HLT'),
                                      
                                      ##                                           dataformat = cms.untracked.int32(0), ### >= 10 data   ### 0 MC
                                      ##                                        rhoCorrection = cms.untracked.InputTag()
                                      ###beamSpotInputTag = cms.InputTag('offlineBeamSpot'),
                                      
                                      rhoCorrection    = cms.untracked.InputTag('kt6PFJetsForRhoComputationEtaMax25','rho'),
                                      rhoCorrection2    = cms.untracked.InputTag('kt6PFJetsForRhoComputationEtaMaxDefault','rho'),
                                      
                                      debugLevel = cms.int32(0)
                                      )

process.out_step = cms.EndPath(process.recoAnalyzer )
                                 
#process.rhopath1 = cms.Path(process.kt6PFJets)
process.rhocomputeDefaultEtaMax = cms.Sequence(process.kt6PFJetsForRhoComputationEtaMax25*process.kt6PFJetsForRhoComputationEtaMaxDefault)
process.rhopath2 = cms.Path(process.rhocomputeDefaultEtaMax)
#process.schedule = cms.Schedule(process.rhopath1,process.rhopath2,process.out_step)

process.schedule = cms.Schedule(process.rhopath2,process.out_step)
#process.schedule = cms.Schedule(process.triggerReport)

#process.schedule = cms.Schedule(process.raw2digi_step,process.L1GtObjectMap_step,process.L1Reco_step,process.recoEcal,process.out_step)
