# /users/mgataull/hltpi0_2012/V27 (CMSSW_5_2_0_HLT1)

import FWCore.ParameterSet.Config as cms

process = cms.Process( "HLT" )
process.load("setup_cff")

process.HLTConfigVersion = cms.PSet(
  tableName = cms.string('/users/mgataull/hltpi0_2012/V27')
)

process.streams = cms.PSet(  ALCAP0 = cms.vstring( 'AlCaP0' ) )
process.datasets = cms.PSet(  AlCaP0 = cms.vstring( 'AlCa_EcalEtaEBonly_v1',
  'AlCa_EcalEtaEEonly_v1',
  'AlCa_EcalPi0EBonly_v1',
  'AlCa_EcalPi0EEonly_v1' ) )

process.hltTriggerType = cms.EDFilter( "HLTTriggerTypeFilter",
    SelectedTriggerType = cms.int32( 1 )
)
process.hltGtDigis = cms.EDProducer( "L1GlobalTriggerRawToDigi",
    DaqGtFedId = cms.untracked.int32( 813 ),
    DaqGtInputTag = cms.InputTag( "rawDataCollector" ),
    UnpackBxInEvent = cms.int32( 5 ),
    ActiveBoardsMask = cms.uint32( 0xffff )
)
process.hltGctDigis = cms.EDProducer( "GctRawToDigi",
    unpackSharedRegions = cms.bool( False ),
    numberOfGctSamplesToUnpack = cms.uint32( 1 ),
    verbose = cms.untracked.bool( False ),
    numberOfRctSamplesToUnpack = cms.uint32( 1 ),
    inputLabel = cms.InputTag( "rawDataCollector" ),
    unpackerVersion = cms.uint32( 0 ),
    gctFedId = cms.untracked.int32( 745 ),
    hltMode = cms.bool( True )
)
process.hltL1GtObjectMap = cms.EDProducer( "L1GlobalTrigger",
    TechnicalTriggersUnprescaled = cms.bool( True ),
    ProduceL1GtObjectMapRecord = cms.bool( True ),
    AlgorithmTriggersUnmasked = cms.bool( False ),
    EmulateBxInEvent = cms.int32( 1 ),
    AlgorithmTriggersUnprescaled = cms.bool( True ),
    ProduceL1GtDaqRecord = cms.bool( False ),
    ReadTechnicalTriggerRecords = cms.bool( True ),
    RecordLength = cms.vint32( 3, 0 ),
    TechnicalTriggersUnmasked = cms.bool( False ),
    ProduceL1GtEvmRecord = cms.bool( False ),
    GmtInputTag = cms.InputTag( "hltGtDigis" ),
    TechnicalTriggersVetoUnmasked = cms.bool( True ),
    AlternativeNrBxBoardEvm = cms.uint32( 0 ),
    TechnicalTriggersInputTags = cms.VInputTag( 'simBscDigis' ),
    CastorInputTag = cms.InputTag( "castorL1Digis" ),
    GctInputTag = cms.InputTag( "hltGctDigis" ),
    AlternativeNrBxBoardDaq = cms.uint32( 0 ),
    WritePsbL1GtDaqRecord = cms.bool( False ),
    BstLengthBytes = cms.int32( -1 )
)
process.hltL1extraParticles = cms.EDProducer( "L1ExtraParticlesProd",
    tauJetSource = cms.InputTag( 'hltGctDigis','tauJets' ),
    etHadSource = cms.InputTag( "hltGctDigis" ),
    etTotalSource = cms.InputTag( "hltGctDigis" ),
    centralBxOnly = cms.bool( True ),
    centralJetSource = cms.InputTag( 'hltGctDigis','cenJets' ),
    etMissSource = cms.InputTag( "hltGctDigis" ),
    hfRingEtSumsSource = cms.InputTag( "hltGctDigis" ),
    produceMuonParticles = cms.bool( True ),
    forwardJetSource = cms.InputTag( 'hltGctDigis','forJets' ),
    ignoreHtMiss = cms.bool( False ),
    htMissSource = cms.InputTag( "hltGctDigis" ),
    produceCaloParticles = cms.bool( True ),
    muonSource = cms.InputTag( "hltGtDigis" ),
    isolatedEmSource = cms.InputTag( 'hltGctDigis','isoEm' ),
    nonIsolatedEmSource = cms.InputTag( 'hltGctDigis','nonIsoEm' ),
    hfRingBitCountsSource = cms.InputTag( "hltGctDigis" )
)
process.hltScalersRawToDigi = cms.EDProducer( "ScalersRawToDigi",
    scalersInputTag = cms.InputTag( "rawDataCollector" )
)
process.hltOnlineBeamSpot = cms.EDProducer( "BeamSpotOnlineProducer",
    maxZ = cms.double( 40.0 ),
    src = cms.InputTag( "hltScalersRawToDigi" ),
    gtEvmLabel = cms.InputTag( "" ),
    changeToCMSCoordinates = cms.bool( False ),
    setSigmaZ = cms.double( 0.0 ),
    maxRadius = cms.double( 2.0 )
)
process.hltOfflineBeamSpot = cms.EDProducer( "BeamSpotProducer" )
process.hltL1sAlCaEcalPi0Eta = cms.EDFilter( "HLTLevel1GTSeed",
    saveTags = cms.bool( True ),
    L1SeedsLogicalExpression = cms.string( "L1_SingleEG5 OR L1_SingleEG7 OR L1_SingleEG12 OR L1_SingleEG20 OR L1_SingleEG22 OR L1_SingleEG24 OR L1_SingleEG30 OR L1_DoubleEG_13_7 OR L1_TripleEG7 OR L1_TripleEG_12_7_5 OR L1_DoubleEG5 OR L1_TripleJet_64_44_24_VBF OR L1_TripleJet_64_48_28_VBF OR  L1_TripleJetC_52_28_28 OR L1_QuadJetC32 OR L1_QuadJetC36 OR L1_QuadJetC40  OR L1_DoubleEG6_HTT100 OR  L1_DoubleEG6_HTT125 OR L1_EG8_DoubleJetC20 OR L1_Mu12_EG7 OR L1_MuOpen_EG12 OR L1_DoubleMu3p5_EG5 OR L1_DoubleMu5_EG5 OR L1_Mu12_EG7 OR L1_Mu5_DoubleEG5 OR L1_Mu5_DoubleEG6 OR L1_MuOpen_EG5" ),
    L1MuonCollectionTag = cms.InputTag( "hltL1extraParticles" ),
    L1UseL1TriggerObjectMaps = cms.bool( True ),
    L1UseAliasesForSeeding = cms.bool( True ),
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    L1CollectionsTag = cms.InputTag( "hltL1extraParticles" ),
    L1NrBxInEvent = cms.int32( 3 ),
    L1GtObjectMapTag = cms.InputTag( "hltL1GtObjectMap" ),
    L1TechTriggerSeeding = cms.bool( False )
)
process.hltPreAlCaEcalEtaEEonly = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltESRawToRecHitFacility = cms.EDProducer( "EcalRawToRecHitFacility",
    sourceTag = cms.InputTag( "rawDataCollector" ),
    workerName = cms.string( "hltESPESUnpackerWorker" )
)
process.hltEcalRawToRecHitFacility = cms.EDProducer( "EcalRawToRecHitFacility",
    sourceTag = cms.InputTag( "rawDataCollector" ),
    workerName = cms.string( "" )
)
process.hltEcalRegionalPi0EtaFEDs = cms.EDProducer( "EcalRawToRecHitRoI",
    JetJobPSet = cms.VPSet( 
    ),
    sourceTag_es = cms.InputTag( "hltESRawToRecHitFacility" ),
    doES = cms.bool( True ),
    type = cms.string( "egamma" ),
    sourceTag = cms.InputTag( "hltEcalRawToRecHitFacility" ),
    EmJobPSet = cms.VPSet( 
      cms.PSet(  regionEtaMargin = cms.double( 0.25 ),
        regionPhiMargin = cms.double( 0.4 ),
        Ptmin = cms.double( 2.0 ),
        Source = cms.InputTag( 'hltL1extraParticles','Isolated' )
      ),
      cms.PSet(  regionEtaMargin = cms.double( 0.25 ),
        regionPhiMargin = cms.double( 0.4 ),
        Ptmin = cms.double( 2.0 ),
        Source = cms.InputTag( 'hltL1extraParticles','NonIsolated' )
      )
    ),
    CandJobPSet = cms.VPSet( 
    ),
    MuonJobPSet = cms.PSet(  ),
    esInstance = cms.untracked.string( "es" ),
    MuJobPSet = cms.PSet(  )
)
process.hltESRegionalPi0EtaRecHit = cms.EDProducer( "EcalRawToRecHitProducer",
    splitOutput = cms.bool( False ),
    rechitCollection = cms.string( "EcalRecHitsES" ),
    EErechitCollection = cms.string( "" ),
    EBrechitCollection = cms.string( "" ),
    sourceTag = cms.InputTag( 'hltEcalRegionalPi0EtaFEDs','es' ),
    lazyGetterTag = cms.InputTag( "hltESRawToRecHitFacility" )
)
process.hltEcalRegionalPi0EtaRecHit = cms.EDProducer( "EcalRawToRecHitProducer",
    splitOutput = cms.bool( True ),
    rechitCollection = cms.string( "NotNeededsplitOutputTrue" ),
    EErechitCollection = cms.string( "EcalRecHitsEE" ),
    EBrechitCollection = cms.string( "EcalRecHitsEB" ),
    sourceTag = cms.InputTag( "hltEcalRegionalPi0EtaFEDs" ),
    lazyGetterTag = cms.InputTag( "hltEcalRawToRecHitFacility" )
)
process.hltSimple3x3Clusters = cms.EDProducer( "EgammaHLTNxNClusterProducer",
    statusLevelRecHitsToUse = cms.int32( 1 ),
    barrelClusterCollection = cms.string( "Simple3x3ClustersBarrel" ),
    flagLevelRecHitsToUse = cms.int32( 1 ),
    maxNumberofClusters = cms.int32( 38 ),
    clusPhiSize = cms.int32( 3 ),
    posCalcParameters = cms.PSet( 
      T0_barl = cms.double( 7.4 ),
      LogWeighted = cms.bool( True ),
      T0_endc = cms.double( 3.1 ),
      T0_endcPresh = cms.double( 1.2 ),
      W0 = cms.double( 4.2 ),
      X0 = cms.double( 0.89 )
    ),
    clusEtaSize = cms.int32( 3 ),
    useRecoFlag = cms.bool( False ),
    endcapHitProducer = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEE' ),
    maxNumberofSeeds = cms.int32( 250 ),
    endcapClusterCollection = cms.string( "Simple3x3ClustersEndcap" ),
    debugLevel = cms.int32( 0 ),
    barrelHitProducer = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEB' ),
    clusSeedThrEndCap = cms.double( 1.0 ),
    doBarrel = cms.bool( True ),
    doEndcaps = cms.bool( True ),
    useDBStatus = cms.bool( True ),
    clusSeedThr = cms.double( 0.5 )
)
process.hltAlCaEtaRecHitsFilterEEonly = cms.EDFilter( "HLTEcalResonanceFilter",
    barrelSelection = cms.PSet( 
      massLowPi0Cand = cms.double( 0.084 ),
      selePtGamma = cms.double( 1.2 ),
      seleMinvMaxBarrel = cms.double( 0.8 ),
      selePtPair = cms.double( 4.0 ),
      seleMinvMinBarrel = cms.double( 0.3 ),
      seleS4S9Gamma = cms.double( 0.87 ),
      seleS9S25Gamma = cms.double( 0.8 ),
      seleIso = cms.double( 0.5 ),
      seleBeltDR = cms.double( 0.3 ),
      ptMinForIsolation = cms.double( 0.5 ),
      store5x5RecHitEB = cms.bool( True ),
      seleBeltDeta = cms.double( 0.1 ),
      removePi0CandidatesForEta = cms.bool( True ),
      barrelHitCollection = cms.string( "etaEcalRecHitsEB" ),
      massHighPi0Cand = cms.double( 0.156 )
    ),
    statusLevelRecHitsToUse = cms.int32( 1 ),
    endcapHits = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEE' ),
    doSelBarrel = cms.bool( False ),
    flagLevelRecHitsToUse = cms.int32( 1 ),
    preshRecHitProducer = cms.InputTag( 'hltESRegionalPi0EtaRecHit','EcalRecHitsES' ),
    doSelEndcap = cms.bool( True ),
    storeRecHitES = cms.bool( True ),
    endcapClusters = cms.InputTag( 'hltSimple3x3Clusters','Simple3x3ClustersEndcap' ),
    barrelHits = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEB' ),
    useRecoFlag = cms.bool( False ),
    barrelClusters = cms.InputTag( 'hltSimple3x3Clusters','Simple3x3ClustersBarrel' ),
    debugLevel = cms.int32( 0 ),
    endcapSelection = cms.PSet( 
      selePtGammaEndCap_region1 = cms.double( 1.0 ),
      region2_EndCap = cms.double( 2.5 ),
      selePtGammaEndCap_region2 = cms.double( 1.0 ),
      ptMinForIsolationEndCap = cms.double( 0.5 ),
      region1_EndCap = cms.double( 2.0 ),
      selePtGammaEndCap_region3 = cms.double( 0.7 ),
      selePtPairMaxEndCap_region3 = cms.double( 9999.0 ),
      seleMinvMinEndCap = cms.double( 0.2 ),
      seleS4S9GammaEndCap = cms.double( 0.9 ),
      seleS9S25GammaEndCap = cms.double( 0.85 ),
      selePtPairEndCap_region1 = cms.double( 3.0 ),
      seleBeltDREndCap = cms.double( 0.3 ),
      selePtPairEndCap_region3 = cms.double( 3.0 ),
      selePtPairEndCap_region2 = cms.double( 3.0 ),
      seleIsoEndCap = cms.double( 0.5 ),
      seleMinvMaxEndCap = cms.double( 0.9 ),
      endcapHitCollection = cms.string( "etaEcalRecHitsEE" ),
      seleBeltDetaEndCap = cms.double( 0.1 ),
      store5x5RecHitEE = cms.bool( True )
    ),
    preshowerSelection = cms.PSet( 
      preshCalibGamma = cms.double( 0.024 ),
      preshStripEnergyCut = cms.double( 0.0 ),
      debugLevelES = cms.string( "" ),
      preshCalibPlaneY = cms.double( 0.7 ),
      preshCalibPlaneX = cms.double( 1.0 ),
      preshCalibMIP = cms.double( 9.0E-5 ),
      ESCollection = cms.string( "etaEcalRecHitsES" ),
      preshNclust = cms.int32( 4 ),
      preshClusterEnergyCut = cms.double( 0.0 ),
      preshSeededNstrip = cms.int32( 15 )
    ),
    useDBStatus = cms.bool( True )
)





process.hltBoolEnd = cms.EDFilter( "HLTBool",
    result = cms.bool( True )
)
process.hltPreAlCaEcalEtaEBonly = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltAlCaEtaRecHitsFilterEBonly = cms.EDFilter( "HLTEcalResonanceFilter",
    barrelSelection = cms.PSet( 
      massLowPi0Cand = cms.double( 0.084 ),
      selePtGamma = cms.double( 1.2 ),
      seleMinvMaxBarrel = cms.double( 0.8 ),
      selePtPair = cms.double( 4.0 ),
      seleMinvMinBarrel = cms.double( 0.3 ),
      seleS4S9Gamma = cms.double( 0.87 ),
      seleS9S25Gamma = cms.double( 0.8 ),
      seleIso = cms.double( 0.5 ),
      seleBeltDR = cms.double( 0.3 ),
      ptMinForIsolation = cms.double( 0.5 ),
      store5x5RecHitEB = cms.bool( True ),
      seleBeltDeta = cms.double( 0.1 ),
      removePi0CandidatesForEta = cms.bool( True ),
      barrelHitCollection = cms.string( "etaEcalRecHitsEB" ),
      massHighPi0Cand = cms.double( 0.156 )
    ),
    statusLevelRecHitsToUse = cms.int32( 1 ),
    endcapHits = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEE' ),
    doSelBarrel = cms.bool( True ),
    flagLevelRecHitsToUse = cms.int32( 1 ),
    preshRecHitProducer = cms.InputTag( 'hltESRegionalPi0EtaRecHit','EcalRecHitsES' ),
    doSelEndcap = cms.bool( False ),
    storeRecHitES = cms.bool( True ),
    endcapClusters = cms.InputTag( 'hltSimple3x3Clusters','Simple3x3ClustersEndcap' ),
    barrelHits = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEB' ),
    useRecoFlag = cms.bool( False ),
    barrelClusters = cms.InputTag( 'hltSimple3x3Clusters','Simple3x3ClustersBarrel' ),
    debugLevel = cms.int32( 0 ),
    endcapSelection = cms.PSet( 
      selePtGammaEndCap_region1 = cms.double( 1.0 ),
      region2_EndCap = cms.double( 2.5 ),
      selePtGammaEndCap_region2 = cms.double( 1.0 ),
      ptMinForIsolationEndCap = cms.double( 0.5 ),
      region1_EndCap = cms.double( 2.0 ),
      selePtGammaEndCap_region3 = cms.double( 0.7 ),
      selePtPairMaxEndCap_region3 = cms.double( 9999.0 ),
      seleMinvMinEndCap = cms.double( 0.2 ),
      seleS4S9GammaEndCap = cms.double( 0.9 ),
      seleS9S25GammaEndCap = cms.double( 0.85 ),
      selePtPairEndCap_region1 = cms.double( 3.0 ),
      seleBeltDREndCap = cms.double( 0.3 ),
      selePtPairEndCap_region3 = cms.double( 3.0 ),
      selePtPairEndCap_region2 = cms.double( 3.0 ),
      seleIsoEndCap = cms.double( 0.5 ),
      seleMinvMaxEndCap = cms.double( 0.9 ),
      endcapHitCollection = cms.string( "etaEcalRecHitsEE" ),
      seleBeltDetaEndCap = cms.double( 0.1 ),
      store5x5RecHitEE = cms.bool( True )
    ),
    preshowerSelection = cms.PSet( 
      preshCalibGamma = cms.double( 0.024 ),
      preshStripEnergyCut = cms.double( 0.0 ),
      debugLevelES = cms.string( "" ),
      preshCalibPlaneY = cms.double( 0.7 ),
      preshCalibPlaneX = cms.double( 1.0 ),
      preshCalibMIP = cms.double( 9.0E-5 ),
      ESCollection = cms.string( "etaEcalRecHitsES" ),
      preshNclust = cms.int32( 4 ),
      preshClusterEnergyCut = cms.double( 0.0 ),
      preshSeededNstrip = cms.int32( 15 )
    ),
    useDBStatus = cms.bool( True )
)
process.hltPreAlCaEcalPi0EBonly = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltAlCaPi0RecHitsFilterEBonly = cms.EDFilter( "HLTEcalResonanceFilter",
    barrelSelection = cms.PSet( 
      massLowPi0Cand = cms.double( 0.084 ),
      selePtGamma = cms.double( 1.3 ),
      seleMinvMaxBarrel = cms.double( 0.23 ),
      selePtPair = cms.double( 2.6 ),
      seleMinvMinBarrel = cms.double( 0.04 ),
      seleS4S9Gamma = cms.double( 0.83 ),
      seleS9S25Gamma = cms.double( 0.0 ),
      seleIso = cms.double( 0.5 ),
      seleBeltDR = cms.double( 0.2 ),
      ptMinForIsolation = cms.double( 0.5 ),
      store5x5RecHitEB = cms.bool( False ),
      seleBeltDeta = cms.double( 0.05 ),
      removePi0CandidatesForEta = cms.bool( False ),
      barrelHitCollection = cms.string( "pi0EcalRecHitsEB" ),
      massHighPi0Cand = cms.double( 0.156 )
    ),
    statusLevelRecHitsToUse = cms.int32( 1 ),
    endcapHits = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEE' ),
    doSelBarrel = cms.bool( True ),
    flagLevelRecHitsToUse = cms.int32( 1 ),
    preshRecHitProducer = cms.InputTag( 'hltESRegionalPi0EtaRecHit','EcalRecHitsES' ),
    doSelEndcap = cms.bool( False ),
    storeRecHitES = cms.bool( True ),
    endcapClusters = cms.InputTag( 'hltSimple3x3Clusters','Simple3x3ClustersEndcap' ),
    barrelHits = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEB' ),
    useRecoFlag = cms.bool( False ),
    barrelClusters = cms.InputTag( 'hltSimple3x3Clusters','Simple3x3ClustersBarrel' ),
    debugLevel = cms.int32( 0 ),
    endcapSelection = cms.PSet( 
      selePtGammaEndCap_region1 = cms.double( 0.6 ),
      region2_EndCap = cms.double( 2.5 ),
      selePtGammaEndCap_region2 = cms.double( 0.6 ),
      ptMinForIsolationEndCap = cms.double( 0.5 ),
      region1_EndCap = cms.double( 2.0 ),
      selePtGammaEndCap_region3 = cms.double( 0.5 ),
      selePtPairMaxEndCap_region3 = cms.double( 2.5 ),
      seleMinvMinEndCap = cms.double( 0.05 ),
      seleS4S9GammaEndCap = cms.double( 0.9 ),
      seleS9S25GammaEndCap = cms.double( 0.0 ),
      selePtPairEndCap_region1 = cms.double( 2.5 ),
      seleBeltDREndCap = cms.double( 0.2 ),
      selePtPairEndCap_region3 = cms.double( 99.0 ),
      selePtPairEndCap_region2 = cms.double( 1.5 ),
      seleIsoEndCap = cms.double( 0.5 ),
      seleMinvMaxEndCap = cms.double( 0.3 ),
      endcapHitCollection = cms.string( "pi0EcalRecHitsEE" ),
      seleBeltDetaEndCap = cms.double( 0.05 ),
      store5x5RecHitEE = cms.bool( False )
    ),
    preshowerSelection = cms.PSet( 
      preshCalibGamma = cms.double( 0.024 ),
      preshStripEnergyCut = cms.double( 0.0 ),
      debugLevelES = cms.string( "" ),
      preshCalibPlaneY = cms.double( 0.7 ),
      preshCalibPlaneX = cms.double( 1.0 ),
      preshCalibMIP = cms.double( 9.0E-5 ),
      ESCollection = cms.string( "pi0EcalRecHitsES" ),
      preshNclust = cms.int32( 4 ),
      preshClusterEnergyCut = cms.double( 0.0 ),
      preshSeededNstrip = cms.int32( 15 )
    ),
    useDBStatus = cms.bool( True )
)
process.hltPreAlCaEcalPi0EEonly = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)
process.hltAlCaPi0RecHitsFilterEEonly = cms.EDFilter( "HLTEcalResonanceFilter",
    barrelSelection = cms.PSet( 
      massLowPi0Cand = cms.double( 0.084 ),
      selePtGamma = cms.double( 1.3 ),
      seleMinvMaxBarrel = cms.double( 0.23 ),
      selePtPair = cms.double( 2.6 ),
      seleMinvMinBarrel = cms.double( 0.04 ),
      seleS4S9Gamma = cms.double( 0.83 ),
      seleS9S25Gamma = cms.double( 0.0 ),
      seleIso = cms.double( 0.5 ),
      seleBeltDR = cms.double( 0.2 ),
      ptMinForIsolation = cms.double( 0.5 ),
      store5x5RecHitEB = cms.bool( False ),
      seleBeltDeta = cms.double( 0.05 ),
      removePi0CandidatesForEta = cms.bool( False ),
      barrelHitCollection = cms.string( "pi0EcalRecHitsEB" ),
      massHighPi0Cand = cms.double( 0.156 )
    ),
    statusLevelRecHitsToUse = cms.int32( 1 ),
    endcapHits = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEE' ),
    doSelBarrel = cms.bool( False ),
    flagLevelRecHitsToUse = cms.int32( 1 ),
    preshRecHitProducer = cms.InputTag( 'hltESRegionalPi0EtaRecHit','EcalRecHitsES' ),
    doSelEndcap = cms.bool( True ),
    storeRecHitES = cms.bool( True ),
    endcapClusters = cms.InputTag( 'hltSimple3x3Clusters','Simple3x3ClustersEndcap' ),
    barrelHits = cms.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEB' ),
    useRecoFlag = cms.bool( False ),
    barrelClusters = cms.InputTag( 'hltSimple3x3Clusters','Simple3x3ClustersBarrel' ),
    debugLevel = cms.int32( 0 ),
    endcapSelection = cms.PSet( 
      selePtGammaEndCap_region1 = cms.double( 0.6 ),
      region2_EndCap = cms.double( 2.5 ),
      selePtGammaEndCap_region2 = cms.double( 0.6 ),
      ptMinForIsolationEndCap = cms.double( 0.5 ),
      region1_EndCap = cms.double( 2.0 ),
      selePtGammaEndCap_region3 = cms.double( 0.5 ),
      selePtPairMaxEndCap_region3 = cms.double( 2.5 ),
      seleMinvMinEndCap = cms.double( 0.05 ),
      seleS4S9GammaEndCap = cms.double( 0.9 ),
      seleS9S25GammaEndCap = cms.double( 0.0 ),
      selePtPairEndCap_region1 = cms.double( 2.5 ),
      seleBeltDREndCap = cms.double( 0.2 ),
      selePtPairEndCap_region3 = cms.double( 99.0 ),
      selePtPairEndCap_region2 = cms.double( 1.5 ),
      seleIsoEndCap = cms.double( 0.5 ),
      seleMinvMaxEndCap = cms.double( 0.3 ),
      endcapHitCollection = cms.string( "pi0EcalRecHitsEE" ),
      seleBeltDetaEndCap = cms.double( 0.05 ),
      store5x5RecHitEE = cms.bool( False )
    ),
    preshowerSelection = cms.PSet( 
      preshCalibGamma = cms.double( 0.024 ),
      preshStripEnergyCut = cms.double( 0.0 ),
      debugLevelES = cms.string( "" ),
      preshCalibPlaneY = cms.double( 0.7 ),
      preshCalibPlaneX = cms.double( 1.0 ),
      preshCalibMIP = cms.double( 9.0E-5 ),
      ESCollection = cms.string( "pi0EcalRecHitsES" ),
      preshNclust = cms.int32( 4 ),
      preshClusterEnergyCut = cms.double( 0.0 ),
      preshSeededNstrip = cms.int32( 15 )
    ),
    useDBStatus = cms.bool( True )
)
process.hltPreALCAP0Output = cms.EDFilter( "HLTPrescaler",
    L1GtReadoutRecordTag = cms.InputTag( "hltGtDigis" ),
    offset = cms.uint32( 0 )
)


#process.hltOutputALCAP0 = cms.OutputModule( "PoolOutputModule",
#    fileName = cms.untracked.string( "outputALCAP0.root" ),
#    fastCloning = cms.untracked.bool( False ),
#    dataset = cms.untracked.PSet(
#        filterName = cms.untracked.string( "" ),
#        dataTier = cms.untracked.string( "RAW" )
#    ),
#    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'AlCa_EcalEtaEBonly_v1',
#  'AlCa_EcalEtaEEonly_v1',
#  'AlCa_EcalPi0EBonly_v1',
###  'AlCa_EcalPi0EEonly_v1' ) ),
#    outputCommands = cms.untracked.vstring( 'drop *',
#      'keep *_hltAlCaEtaRecHitsFilterEBonly_*_*',
#     'keep *_hltAlCaEtaRecHitsFilterEEonly_*_*',
#      'keep *_hltAlCaPi0RecHitsFilterEBonly_*_*',
#      'keep *_hltAlCaPi0RecHitsFilterEEonly_*_*',
#      'keep L1GlobalTriggerReadoutRecord_hltGtDigis_*_*',
#      'keep edmTriggerResults_*_*_*' )
#)


process.hltOutputALCAP0 = cms.OutputModule( "PoolOutputModule",
    fileName = cms.untracked.string( "outputALCAP0.root" ),
    fastCloning = cms.untracked.bool( False ),
    dataset = cms.untracked.PSet(
        filterName = cms.untracked.string( "" ),
        dataTier = cms.untracked.string( "RAW" )
    ),
    SelectEvents = cms.untracked.PSet(  SelectEvents = cms.vstring( 'AlCa_EcalEtaEBonly_v1',
  'AlCa_EcalEtaEEonly_v1',
  'AlCa_EcalPi0EBonly_v1',
  'AlCa_EcalPi0EEonly_v1' ) ),
    outputCommands = cms.untracked.vstring( 'drop *',
      'keep *_hltAlCaEtaEEUncalibrator_*_*',
      'keep *_hltAlCaEtaEBUncalibrator_*_*',
      'keep *_hltAlCaPi0EEUncalibrator_*_*',
      'keep *_hltAlCaPi0EBUncalibrator_*_*',
      'keep L1GlobalTriggerReadoutRecord_hltGtDigis_*_*',
      'keep edmTriggerResults_*_*_*' )
)

process.pi0Analyzer = cms.EDAnalyzer("Pi0Analyzer",
                                     outputFile = cms.string('clusters_MinimumBiasRun2011B-v1RAW_180072.root'),
                                     getHLT = cms.untracked.bool (True),
                                     usegtDigis = cms.untracked.bool (True),
                                     l1GtRecordInputTag = cms.untracked.InputTag("simGtDigis"),
                                     beamSpotInputTag = cms.untracked.InputTag("hltOfflineBeamSpot"),
                                     saveAllPhotonBarrel = cms.untracked.bool (True),
                                     saveAllPhotonEndcap = cms.untracked.bool (True),
                                     doSelForPi0Barrel = cms.bool( False ),
                                     barrelHits = cms.untracked.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEB' ),
                                     doSelForPi0Endcap = cms.bool( False ),
                                     endcapHits = cms.untracked.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEE' ),
                                     
                                     doSelForEtaBarrel = cms.bool( False ),
                                     barrelHitsEta = cms.untracked.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEB' ),
                                     doSelForEtaEndcap = cms.bool( False ),
                                     barrelHitsEtaEndcap = cms.untracked.InputTag( 'hltEcalRegionalPi0EtaRecHit','EcalRecHitsEE' ),
                                     
                                     posCalcParameters = cms.untracked.PSet( T0_barl      = cms.double(7.4),
                                                                             T0_endc      = cms.double(3.1),
                                                                             T0_endcPresh = cms.double(1.2),
                                                                             LogWeighted  = cms.bool(True),
                                                                             W0           = cms.double(4.2),
                                                                             X0           = cms.double(0.89)
                                                                             )          ,
                                     debugLevel = cms.int32(0)
                                     )

process.HLTL1UnpackerSequence = cms.Sequence( process.hltGtDigis + process.hltGctDigis + process.hltL1GtObjectMap + process.hltL1extraParticles )
process.HLTBeamSpot = cms.Sequence( process.hltScalersRawToDigi + process.hltOnlineBeamSpot + process.hltOfflineBeamSpot )
process.HLTBeginSequence = cms.Sequence( process.hltTriggerType + process.HLTL1UnpackerSequence + process.HLTBeamSpot )
process.HLTDoRegionalPi0EtaSequence = cms.Sequence( process.hltESRawToRecHitFacility + process.hltEcalRawToRecHitFacility + process.hltEcalRegionalPi0EtaFEDs + process.hltESRegionalPi0EtaRecHit + process.hltEcalRegionalPi0EtaRecHit )
process.HLTEndSequence = cms.Sequence( process.hltBoolEnd )

#process.AlCa_EcalEtaEEonly_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sAlCaEcalPi0Eta + process.hltPreAlCaEcalEtaEEonly + process.HLTDoRegionalPi0EtaSequence + process.hltSimple3x3Clusters + process.hltAlCaEtaRecHitsFilterEEonly + process.hltAlCaEtaEEUncalibrator + process.HLTEndSequence )

#process.AlCa_EcalEtaEBonly_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sAlCaEcalPi0Eta + process.hltPreAlCaEcalEtaEBonly + process.HLTDoRegionalPi0EtaSequence + process.hltSimple3x3Clusters + process.hltAlCaEtaRecHitsFilterEBonly + process.hltAlCaEtaEBUncalibrator + process.HLTEndSequence )

#process.AlCa_EcalPi0EBonly_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sAlCaEcalPi0Eta + process.hltPreAlCaEcalPi0EBonly + process.HLTDoRegionalPi0EtaSequence + process.hltSimple3x3Clusters + process.hltAlCaPi0RecHitsFilterEBonly + process.hltAlCaPi0EBUncalibrator + process.HLTEndSequence )

#process.AlCa_EcalPi0EEonly_v1 = cms.Path( process.HLTBeginSequence + process.hltL1sAlCaEcalPi0Eta + process.hltPreAlCaEcalPi0EEonly + process.HLTDoRegionalPi0EtaSequence + process.hltSimple3x3Clusters + process.hltAlCaPi0RecHitsFilterEEonly + process.hltAlCaPi0EEUncalibrator + process.HLTEndSequence )

#process.ALCAP0Output = cms.EndPath( process.hltPreALCAP0Output + process.hltOutputALCAP0 )

process.regionalpath = cms.Path(process.HLTBeginSequence  + process.HLTDoRegionalPi0EtaSequence + process.pi0Analyzer )
##process.out_step = cms.EndPath(process.pi0Analyzer )


process.source = cms.Source( "PoolSource",
    fileNames = cms.untracked.vstring(
        '/store/data/Run2011B/MinimumBias/RAW/v1/000/178/160/1A6B0678-0DF3-E011-9197-002481E0D646.root'
    ),
    secondaryFileNames = cms.untracked.vstring(
    ),
    inputCommands = cms.untracked.vstring(
        'keep *'
    )
)

# En-able HF Noise filters in GRun menu
if 'hltHfreco' in process.__dict__:
    process.hltHfreco.setNoiseFlags = cms.bool( True )

# remove the HLT prescales
if 'PrescaleService' in process.__dict__:
    process.PrescaleService.lvl1DefaultLabel = cms.string( '0' )
    process.PrescaleService.lvl1Labels       = cms.vstring( '0', '1', '2', '3', '4', '5', '6', '7', '8', '9' )
    process.PrescaleService.prescaleTable    = cms.VPSet( )

# CMSSW version specific customizations
import os
cmsswVersion = os.environ['CMSSW_VERSION']

# override the process name
process.setName_('TEST')

# adapt HLT modules to the correct process name
if 'hltTrigReport' in process.__dict__:
    process.hltTrigReport.HLTriggerResults                    = cms.InputTag( 'TriggerResults', '', 'TEST' )

if 'hltPreExpressCosmicsOutputSmart' in process.__dict__:
    process.hltPreExpressCosmicsOutputSmart.TriggerResultsTag = cms.InputTag( 'TriggerResults', '', 'TEST' )

if 'hltPreExpressOutputSmart' in process.__dict__:
    process.hltPreExpressOutputSmart.TriggerResultsTag        = cms.InputTag( 'TriggerResults', '', 'TEST' )

if 'hltPreDQMForHIOutputSmart' in process.__dict__:
    process.hltPreDQMForHIOutputSmart.TriggerResultsTag       = cms.InputTag( 'TriggerResults', '', 'TEST' )

if 'hltPreDQMForPPOutputSmart' in process.__dict__:
    process.hltPreDQMForPPOutputSmart.TriggerResultsTag       = cms.InputTag( 'TriggerResults', '', 'TEST' )

if 'hltPreHLTDQMResultsOutputSmart' in process.__dict__:
    process.hltPreHLTDQMResultsOutputSmart.TriggerResultsTag  = cms.InputTag( 'TriggerResults', '', 'TEST' )

if 'hltPreHLTDQMOutputSmart' in process.__dict__:
    process.hltPreHLTDQMOutputSmart.TriggerResultsTag         = cms.InputTag( 'TriggerResults', '', 'TEST' )

if 'hltPreHLTMONOutputSmart' in process.__dict__:
    process.hltPreHLTMONOutputSmart.TriggerResultsTag         = cms.InputTag( 'TriggerResults', '', 'TEST' )

if 'hltDQMHLTScalers' in process.__dict__:
    process.hltDQMHLTScalers.triggerResults                   = cms.InputTag( 'TriggerResults', '', 'TEST' )
    process.hltDQMHLTScalers.processname                      = 'TEST'

if 'hltDQML1SeedLogicScalers' in process.__dict__:
    process.hltDQML1SeedLogicScalers.processname              = 'TEST'

# limit the number of events to be processed
process.maxEvents = cms.untracked.PSet(
    input = cms.untracked.int32( 100 )
)

# enable the TrigReport and TimeReport
process.options = cms.untracked.PSet(
    wantSummary = cms.untracked.bool( True )
)

# override the GlobalTag, connection string and pfnPrefix
if 'GlobalTag' in process.__dict__:
    process.GlobalTag.connect   = 'frontier://FrontierProd/CMS_COND_31X_GLOBALTAG'
    process.GlobalTag.pfnPrefix = cms.untracked.string('frontier://FrontierProd/')
    from Configuration.AlCa.autoCond import autoCond
    process.GlobalTag.globaltag = autoCond['hltonline']

# override the L1 menu
if 'GlobalTag' in process.__dict__:
    process.GlobalTag.toGet.append(
        cms.PSet(
            record  = cms.string( 'L1GtTriggerMenuRcd' ),
            tag     = cms.string( 'L1GtTriggerMenu_L1Menu_Collisions2012_v0_mc' ),
            label   = cms.untracked.string( '' ),
            connect = cms.untracked.string( 'frontier://FrontierProd/CMS_COND_31X_L1T' )
        )
    )

# customize the L1 emulator to run customiseL1GtEmulatorFromRaw with HLT to switchToSimGtDigis
process.load( 'Configuration.StandardSequences.RawToDigi_Data_cff' )
process.load( 'Configuration.StandardSequences.SimL1Emulator_cff' )
import L1Trigger.Configuration.L1Trigger_custom
process = L1Trigger.Configuration.L1Trigger_custom.customiseL1GtEmulatorFromRaw( process )
process = L1Trigger.Configuration.L1Trigger_custom.customiseResetPrescalesAndMasks( process )

# customize the HLT to use the emulated results
import HLTrigger.Configuration.customizeHLTforL1Emulator
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToL1Emulator( process )
process = HLTrigger.Configuration.customizeHLTforL1Emulator.switchToSimGtDigis( process )

if 'MessageLogger' in process.__dict__:
    process.MessageLogger.categories.append('TriggerSummaryProducerAOD')
    process.MessageLogger.categories.append('L1GtTrigReport')
    process.MessageLogger.categories.append('HLTrigReport')

