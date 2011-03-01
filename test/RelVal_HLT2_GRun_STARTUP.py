# Auto generated configuration file
# using: 
# Revision: 1.284 
# Source: /cvs_server/repositories/CMSSW/CMSSW/Configuration/PyReleaseValidation/python/ConfigBuilder.py,v 
# with command line options: RelVal --step=HLT:GRun --conditions=auto:startup --filein=file:RelVal_HLT_GRun_STARTUP.root --custom_conditions= --fileout=RelVal_HLT2_GRun_STARTUP.root --number=100 --mc --no_exec --datatier RAW-HLT --eventcontent=FEVTDEBUGHLT --customise=HLTrigger/Configuration/CustomConfigs.L1THLT --python_filename=RelVal_HLT2_GRun_STARTUP.py --processName=HLT2
import FWCore.ParameterSet.Config as cms

process = cms.Process('HLT2')

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('SimGeneral.HepPDTESSource.pythiapdt_cfi')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('SimGeneral.MixingModule.mixNoPU_cfi')
process.load('Configuration.StandardSequences.GeometryDB_cff')
process.load('Configuration.StandardSequences.MagneticField_38T_cff')
process.load('HLTrigger.Configuration.HLT_GRun_cff')
process.load('Configuration.StandardSequences.EndOfProcess_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')
process.load('Configuration.EventContent.EventContent_cff')

process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32(-1)
)

# Input source
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    '/store/relval/CMSSW_3_11_0/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/START311_V1_64bit-v2/0011/2CCA681C-042B-E011-8A05-001D09F23F2A.root',
    '/store/relval/CMSSW_3_11_0/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/START311_V1_64bit-v2/0011/C2AFB799-022B-E011-BA3C-0030487C8CB6.root',
    '/store/relval/CMSSW_3_11_0/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/START311_V1_64bit-v2/0011/DACE2A5D-032B-E011-A90D-003048F118DE.root',
    '/store/relval/CMSSW_3_11_0/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/START311_V1_64bit-v2/0011/FABFCD5D-032B-E011-A9C1-003048F1C420.root',
    '/store/relval/CMSSW_3_11_0/RelValMinBias/GEN-SIM-DIGI-RAW-HLTDEBUG/START311_V1_64bit-v2/0012/0C1E411B-B72B-E011-BB5E-003048F11C5C.root'
                                      )
)

process.options = cms.untracked.PSet(

)

# Production Info
process.configurationMetadata = cms.untracked.PSet(
    version = cms.untracked.string('$Revision: 1.284 $'),
    annotation = cms.untracked.string('RelVal nevts:100'),
    name = cms.untracked.string('PyReleaseValidation')
)

# Output definition


# Additional output definition

# Other statements
process.GlobalTag.globaltag = 'START311_V1::All'

# Path and EndPath definitions
#process.endjob_step = cms.EndPath(process.endOfProcess)
#process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)

process.FEVTDEBUGHLToutput = cms.OutputModule("PoolOutputModule",
    splitLevel = cms.untracked.int32(0),
    #outputCommands = process.FEVTDEBUGHLTEventContent.outputCommands,
    fileName = cms.untracked.string('RelVal_HLT2_MB_GRun_STARTUP.root'),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'AlCa_EcalPi0','AlCa_EcalEta' ) ),
    outputCommands = cms.untracked.vstring( 'drop *',
    'keep L1GlobalTriggerReadoutRecord_hltGtDigis_*_HLT2',
    'keep *_hltAlCaPi0RecHitsFilter_*_HLT2' ,
    'keep *_hltAlCaEtaRecHitsFilter_*_HLT2',
    'keep *_TriggerResults_*_HLT2',
    'keep *_hltTriggerSummaryAOD_*_HLT2',                                        
    'keep *_hltOfflineBeamSpot_*_HLT2'                                     
    ),                                             
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string(''),
        dataTier = cms.untracked.string('RAW-HLT')
    )
)


process.EcalPi0Output = cms.OutputModule("EventStreamFileWriter",
  fileName = cms.untracked.string( "AlCaRawEcalPi0Stream_MB_noBeamSpot.dat" ),
  compression_level = cms.untracked.int32( 1 ),
  use_compression = cms.untracked.bool(True),
  SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'AlCa_EcalPi0','AlCa_EcalEta' ) ),
  outputCommands = cms.untracked.vstring( 'drop *',
    'keep L1GlobalTriggerReadoutRecord_hltGtDigis_*_HLT2',
    'keep *_hltAlCaPi0RecHitsFilter_*_HLT2' ,
    'keep *_hltAlCaEtaRecHitsFilter_*_HLT2',
    'keep *_TriggerResults_*_HLT2',
    'keep *_hltTriggerSummaryAOD_*_HLT2'                                        
    ##'keep *_hltOfflineBeamSpot_*_HLT2'        
  )
)

process.endjob_step = cms.EndPath(process.endOfProcess)
process.AlCaOutput = cms.EndPath(process.EcalPi0Output)
#process.FEVTDEBUGHLToutput_step = cms.EndPath(process.FEVTDEBUGHLToutput)


# Schedule definition
process.schedule = cms.Schedule()
process.schedule.extend(process.HLTSchedule)
#process.schedule.extend([process.endjob_step,process.FEVTDEBUGHLToutput_step])
process.schedule.extend([process.endjob_step,process.AlCaOutput])


# customisation of the process

# Automatic addition of the customisation function from HLTrigger.Configuration.CustomConfigs

from L1Trigger.Configuration import patchToRerunL1Emulator


def Base(process):
#   default modifications

    process.options.wantSummary = cms.untracked.bool(True)

    process.MessageLogger.categories.append('TriggerSummaryProducerAOD')
    process.MessageLogger.categories.append('L1GtTrigReport')
    process.MessageLogger.categories.append('HLTrigReport')

    if 'hltTrigReport' in process.__dict__:
        process.hltTrigReport.HLTriggerResults = cms.InputTag( 'TriggerResults','',process.name_() )

    if 'hltDQMHLTScalers' in process.__dict__:
        process.hltDQMHLTScalers.triggerResults = cms.InputTag( 'TriggerResults','',process.name_() )

    if 'hltDQML1SeedLogicScalers' in process.__dict__:
        process.hltDQML1SeedLogicScalers.processname = process.name_()

    return(process)


def L1T(process):
#   modifications when running L1T only

    process.load('L1Trigger.GlobalTriggerAnalyzer.l1GtTrigReport_cfi')
    process.l1GtTrigReport.L1GtRecordInputTag = cms.InputTag( "simGtDigis" )

    process.L1AnalyzerEndpath = cms.EndPath( process.l1GtTrigReport )
    process.schedule.append(process.L1AnalyzerEndpath)

    process=Base(process)

    return(process)


def L1THLT(process):
#   modifications when running L1T+HLT

    process=Base(process)

    return(process)


def L1THLT2(process):
#   modifications when re-running L1T+HLT    

#   run trigger primitive generation on unpacked digis, then central L1

    process.load("L1Trigger.Configuration.CaloTriggerPrimitives_cff")
    process.simEcalTriggerPrimitiveDigis.Label = 'ecalDigis'
    process.simHcalTriggerPrimitiveDigis.inputLabel = ('hcalDigis', 'hcalDigis')

#   patch the process to use 'sim*Digis' from the L1 emulator
#   instead of 'hlt*Digis' from the RAW data

    patchToRerunL1Emulator.switchToSimGtDigis( process )

    process=Base(process)

    return(process)


def HLTDropPrevious(process):
#   drop on input the previous HLT results
    process.source.inputCommands = cms.untracked.vstring (
        'keep *',
        'drop *_hltL1GtObjectMap_*_*',
        'drop *_TriggerResults_*_*',
        'drop *_hltTriggerSummaryAOD_*_*',
    )

    process=Base(process)
    
    return(process)


process = L1THLT(process)


# End of customisation functions
