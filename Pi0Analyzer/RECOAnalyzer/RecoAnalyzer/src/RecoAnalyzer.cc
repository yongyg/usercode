#include "RECOAnalyzer/RecoAnalyzer/interface/RecoAnalyzer.h"
#include "DataFormats/EgammaReco/interface/BasicCluster.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "DataFormats/EcalDetId/interface/ESDetId.h"
#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/Math/interface/Point3D.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/Common/interface/Handle.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"

#include "DataFormats/RecoCandidate/interface/RecoEcalCandidate.h"
#include "DataFormats/RecoCandidate/interface/RecoEcalCandidateFwd.h"

#include "DataFormats/GeometryVector/interface/LocalPoint.h"

//Level 1 Trigger
#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "Geometry/Records/interface/CaloGeometryRecord.h"

/// EgammaCoreTools
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalEtaPhiRegion.h"


// Ecal Mapping 
#include "DataFormats/EcalRawData/interface/EcalListOfFEDS.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"
#include <DataFormats/FEDRawData/interface/FEDNumbering.h>

// Jets stuff
#include "DataFormats/L1Trigger/interface/L1JetParticle.h"
#include "DataFormats/L1Trigger/interface/L1JetParticleFwd.h"

//// Ecal Electrons Id
#include "DataFormats/EcalDetId/interface/EcalElectronicsId.h"
#include "DataFormats/EcalDetId/interface/EcalTriggerElectronicsId.h"


// ES stuff
#include "Geometry/EcalAlgo/interface/EcalPreshowerGeometry.h"
#include "Geometry/CaloTopology/interface/EcalPreshowerTopology.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"
#include "DataFormats/EgammaReco/interface/PreshowerClusterFwd.h"

//Ecal status
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"

#include "CondFormats/L1TObjects/interface/L1GtTriggerMenuFwd.h"
#include "CondFormats/L1TObjects/interface/L1GtTriggerMenu.h"
#include "CondFormats/DataRecord/interface/L1GtTriggerMenuRcd.h"


#include "CondFormats/EcalObjects/interface/EcalPedestals.h"
#include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h"

#include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h"
#include "CondFormats/EcalObjects/interface/EcalGainRatios.h"
#include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h"

#include "CondFormats/EcalObjects/interface/EcalIntercalibConstantsMC.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsMCRcd.h"

#include "TVector3.h"

#include "DataFormats/TrackerRecHit2D/interface/SiPixelRecHitCollection.h"

#include "CLHEP/Matrix/Vector.h"
#include "CLHEP/Matrix/SymMatrix.h"

#include "FWCore/ServiceRegistry/interface/Service.h"
#include "FWCore/Utilities/interface/RandomNumberGenerator.h"
#include "CLHEP/Random/RandGaussQ.h"

#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"
#include "Geometry/Records/interface/TrackerDigiGeometryRecord.h"
#include "Geometry/TrackerGeometryBuilder/interface/PixelGeomDetUnit.h"
#include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h"
///#include "Geometry/TrackerTopology/interface/RectangularPixelTopology.h"
#include "Geometry/TrackerGeometryBuilder/interface/RectangularPixelTopology.h"

// needed for trigger bits from EventSetup as in ALCARECO paths
#include "CondFormats/HLTObjects/interface/AlCaRecoTriggerBits.h"
#include "CondFormats/DataRecord/interface/AlCaRecoTriggerBitsRcd.h"


//#include "HLTrigger/HLTfilters/interface/HLTHighLevelDev.h"
#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"

#include "DataFormats/HcalDetId/interface/HcalDetId.h"
#include "DataFormats/HcalRecHit/interface/HcalRecHitCollections.h"
#include "DataFormats/CaloTowers/interface/CaloTower.h"
#include "DataFormats/CaloTowers/interface/CaloTowerFwd.h"

#include "DataFormats/TrackReco/interface/Track.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h" 
#include "DataFormats/TrackCandidate/interface/TrackCandidate.h" 
#include "DataFormats/TrackReco/interface/Track.h" 
#include "DataFormats/TrackReco/interface/TrackFwd.h" 

#include "PhysicsTools/FWLite/interface/CommandLineParser.h" 
#include "DataFormats/EgammaCandidates/interface/Conversion.h"


#define TWOPI 6.283185308

using namespace l1extra;
using namespace edm;
using namespace std;
using namespace reco;
//using namespace trigger;




class PtSorter {
public:
  template <class T> bool operator() ( const T& a, const T& b ) {
    return ( a.pt() > b.pt() );
  }
};



bool sort_pred(const std::pair<double,int> & left, const std::pair<double,int>& right)
{
  return left.first < right.first;
}


RecoAnalyzer::RecoAnalyzer(const edm::ParameterSet& iConfig):
  
  //  triggerNamesID_(),
  triggerNames_ (),
  watchAlCaRecoTriggerBitsRcd_(0),
  HLTPrescalesExpanded_(), 
  HLTPathsByName_(),
  HLTPathsByIndex_()
  
{

  cout<<"gettting parameters.." <<endl; 
  clusSeedThr_ = iConfig.getParameter<double> ("clusSeedThr");
  clusSeedThrEndCap_ = iConfig.getParameter<double> ("clusSeedThrEndCap");

  clusEtaSize_ = iConfig.getParameter<int> ("clusEtaSize");
  clusPhiSize_ = iConfig.getParameter<int> ("clusPhiSize");
  if ( clusPhiSize_ % 2 == 0 ||  clusEtaSize_ % 2 == 0) {
    edm::LogError("AlCaPi0RecHitsProducerError") << "Size of eta/phi for simple clustering should be odd numbers, reset to be 3 ";
    clusPhiSize_ = 3; 
    clusEtaSize_ = 3; 
  }
  
  fullRECO_  = iConfig.getParameter<int> ("fullRECO");
  
  
  RegionalMatch_ = iConfig.getParameter<bool>("RegionalMatch");
  ptMinEMObj_ = iConfig.getParameter<double>("ptMinEMObj");
  EMregionEtaMargin_ = iConfig.getParameter<double>("EMregionEtaMargin"); ///0.25, //
  EMregionPhiMargin_ = iConfig.getParameter<double>("EMregionPhiMargin"); ///0.4 
  l1IsolatedTag_ = iConfig.getParameter< edm::InputTag > ("l1IsolatedTag");
  l1NonIsolatedTag_ = iConfig.getParameter< edm::InputTag > ("l1NonIsolatedTag");
  

  //  seleNRHMax_ = iConfig.getParameter<int> ("seleNRHMax");
  seleXtalMinEnergy_ = iConfig.getParameter<double>("seleXtalMinEnergy");
  seleXtalMinEnergyEndCap_ = iConfig.getParameter<double>("seleXtalMinEnergyEndCap");

  saveAllPhotonBarrel_ = iConfig.getParameter<bool>("saveAllPhotonBarrel");
  saveAllPhotonEndcap_ = iConfig.getParameter<bool>("saveAllPhotonEndcap");
  
  srcVertex_= iConfig.getUntrackedParameter<string>("srcVertex",""); 
  
  getBeamSpotOnly_ = iConfig.getUntrackedParameter<bool>("getBeamSpotOnly",false);
  
  m_vertexSrc =  iConfig.getUntrackedParameter<edm::InputTag>("vertex",edm::InputTag("offlinePrimaryVerticesWithBS"));
  
  m_vertexSrc2 =  iConfig.getUntrackedParameter<edm::InputTag>("vertex2",edm::InputTag("offlinePrimaryVertices"));
  
  reRunPixelRecHits_ =  iConfig.getUntrackedParameter<bool>("reRunPixelRecHits",false);
  
  srcPixels_ = iConfig.getUntrackedParameter<string>("pixelRecHits","siPixelRecHits"); 
  m_tracksSrc  = iConfig.getUntrackedParameter<edm::InputTag>("tracks",edm::InputTag("generalTracks"));
  
  
  useRecoFlag_ = iConfig.getParameter<bool>("useRecoFlag");
  flagLevelRecHitsToUse_ = iConfig.getParameter<int>("flagLevelRecHitsToUse"); 
  
  useDBStatus_ = iConfig.getParameter<bool>("useDBStatus");
  statusLevelRecHitsToUse_ = iConfig.getParameter<int>("statusLevelRecHitsToUse"); 
  
  nMinRecHitsSel1stCluster_ = iConfig.getParameter<int>("nMinRecHitsSel1stCluster");
  nMinRecHitsSel2ndCluster_ = iConfig.getParameter<int>("nMinRecHitsSel2ndCluster");
  
  maxNumberofSeeds_    = iConfig.getParameter<int> ("maxNumberofSeeds");
  maxNumberofClusters_ = iConfig.getParameter<int> ("maxNumberofClusters");

  removeSpike_ = iConfig.getUntrackedParameter<int>("removeSpike",0);
  

  cout <<"maxcluster: "<< MAX3x3ClusEB <<" "<<  maxNumberofClusters_ <<endl; 

  srcTowers_ = iConfig.getUntrackedParameter<string>("caloTowers","towerMaker"); 

  
  doSelForPi0Barrel_ = iConfig.getParameter<bool> ("doSelForPi0Barrel");  
  barrelHits_ = iConfig.getParameter< edm::InputTag > ("barrelHits");
  if(doSelForPi0Barrel_){

    ///for Pi0 barrel selection
    selePtGamma_ = iConfig.getParameter<double> ("selePtGamma");  
    selePtPi0_ = iConfig.getParameter<double> ("selePtPi0");  
    seleMinvMaxPi0_ = iConfig.getParameter<double> ("seleMinvMaxPi0");  
    seleMinvMinPi0_ = iConfig.getParameter<double> ("seleMinvMinPi0");  
    seleS4S9Gamma_ = iConfig.getParameter<double> ("seleS4S9Gamma");  
    selePi0BeltDR_ = iConfig.getParameter<double> ("selePi0BeltDR");  
    selePi0BeltDeta_ = iConfig.getParameter<double> ("selePi0BeltDeta");  
    selePi0Iso_ = iConfig.getParameter<double> ("selePi0Iso");  
    ptMinForIsolation_  = iConfig.getParameter<double>("ptMinForIsolation");

  }
  
  doSelForPi0Endcap_ = iConfig.getParameter<bool>("doSelForPi0Endcap");  
  endcapHits_ = iConfig.getParameter< edm::InputTag > ("endcapHits");
  if(doSelForPi0Endcap_){
    

    ///for Pi0 endcap selection
    ///    selePtGammaEndCap_ = iConfig.getParameter<double> ("selePtGammaEndCap");  
    /// selePtPi0EndCap_ = iConfig.getParameter<double> ("selePtPi0EndCap");   


    ///try to divide endcap region into 3 parts
    /// eta< 2 ; eta>2 && eta<2.5 ; eta>2.5; 
    region1_Pi0EndCap_ = iConfig.getParameter<double> ("region1_Pi0EndCap");
    selePtGammaPi0EndCap_region1_ = iConfig.getParameter<double> ("selePtGammaPi0EndCap_region1");  
    selePtPi0EndCap_region1_ = iConfig.getParameter<double> ("selePtPi0EndCap_region1");   
    
    region2_Pi0EndCap_ = iConfig.getParameter<double> ("region2_Pi0EndCap");
    selePtGammaPi0EndCap_region2_ = iConfig.getParameter<double> ("selePtGammaPi0EndCap_region2");  
    selePtPi0EndCap_region2_ = iConfig.getParameter<double> ("selePtPi0EndCap_region2");   
    
    selePtGammaPi0EndCap_region3_ = iConfig.getParameter<double> ("selePtGammaPi0EndCap_region3");  
    selePtPi0EndCap_region3_ = iConfig.getParameter<double> ("selePtPi0EndCap_region3"); 
    selePtPi0MaxEndCap_region3_ = iConfig.getParameter<double> ("selePtPi0MaxEndCap_region3"); 
    
    
    seleS4S9GammaEndCap_ = iConfig.getParameter<double> ("seleS4S9GammaEndCap");  
    seleMinvMaxPi0EndCap_ = iConfig.getParameter<double> ("seleMinvMaxPi0EndCap");  
    seleMinvMinPi0EndCap_ = iConfig.getParameter<double> ("seleMinvMinPi0EndCap");  
    selePi0BeltDREndCap_ = iConfig.getParameter<double> ("selePi0BeltDREndCap");  
    selePi0BeltDetaEndCap_ = iConfig.getParameter<double> ("selePi0BeltDetaEndCap");  

    selePi0IsoEndCap_ = iConfig.getParameter<double> ("selePi0IsoEndCap");  
  
    preshHitProducer_   = iConfig.getParameter<edm::InputTag>("preshRecHitProducer");
    ptMinForIsolationEndCap_  = iConfig.getParameter<double>("ptMinForIsolationEndCap");
    
  }
  
  doSelForEtaBarrel_ = iConfig.getParameter<bool> ("doSelForEtaBarrel");  
  barrelHitsEta_ = iConfig.getParameter< edm::InputTag > ("barrelHitsEta");
  if(doSelForEtaBarrel_){

    selePtGammaEta_ = iConfig.getParameter<double> ("selePtGammaEta");  
    selePtEta_ = iConfig.getParameter<double> ("selePtEta");  
    seleMinvMaxEta_ = iConfig.getParameter<double> ("seleMinvMaxEta");  
    seleMinvMinEta_ = iConfig.getParameter<double> ("seleMinvMinEta");  
    seleS4S9GammaEta_ = iConfig.getParameter<double> ("seleS4S9GammaEta");  
    seleS9S25GammaEta_ = iConfig.getParameter<double> ("seleS9S25GammaEta");  
    
    seleEtaBeltDR_ = iConfig.getParameter<double> ("seleEtaBeltDR");  
    seleEtaBeltDeta_ = iConfig.getParameter<double> ("seleEtaBeltDeta");  
    
    seleEtaIso_ = iConfig.getParameter<double> ("seleEtaIso");  
    ptMinForIsolationEta_  = iConfig.getParameter<double>("ptMinForIsolationEta");
    
    removePi0CandidatesForEta_ = iConfig.getParameter<bool>("removePi0CandidatesForEta");
    massLowPi0Cand_ = iConfig.getParameter<double> ("massLowPi0Cand");
    massHighPi0Cand_ = iConfig.getParameter<double> ("massHighPi0Cand");
    
  }
  
  
  doSelForEtaEndcap_ = iConfig.getParameter<bool>("doSelForEtaEndcap");  
  endcapHitsEta_ = iConfig.getParameter< edm::InputTag > ("endcapHitsEta");
  if(doSelForEtaEndcap_){


    ///for Eta endcap selection
    ///    selePtGammaEndCap_ = iConfig.getParameter<double> ("selePtGammaEndCap");  
    /// selePtEtaEndCap_ = iConfig.getParameter<double> ("selePtEtaEndCap");   


    ///try to divide endcap region into 3 parts
    /// eta< 2 ; eta>2 && eta<2.5 ; eta>2.5; 
    region1_EtaEndCap_ = iConfig.getParameter<double> ("region1_EtaEndCap");
    selePtGammaEtaEndCap_region1_ = iConfig.getParameter<double> ("selePtGammaEtaEndCap_region1");  
    selePtEtaEndCap_region1_ = iConfig.getParameter<double> ("selePtEtaEndCap_region1");   
    
    region2_EtaEndCap_ = iConfig.getParameter<double> ("region2_EtaEndCap");
    selePtGammaEtaEndCap_region2_ = iConfig.getParameter<double> ("selePtGammaEtaEndCap_region2");  
    selePtEtaEndCap_region2_ = iConfig.getParameter<double> ("selePtEtaEndCap_region2");   
    
    selePtGammaEtaEndCap_region3_ = iConfig.getParameter<double> ("selePtGammaEtaEndCap_region3");  
    selePtEtaEndCap_region3_ = iConfig.getParameter<double> ("selePtEtaEndCap_region3"); 
        
    
    seleS4S9GammaEtaEndCap_ = iConfig.getParameter<double> ("seleS4S9GammaEtaEndCap");  
    seleS9S25GammaEtaEndCap_ = iConfig.getParameter<double> ("seleS9S25GammaEtaEndCap");  
    
    seleMinvMaxEtaEndCap_ = iConfig.getParameter<double> ("seleMinvMaxEtaEndCap");  
    seleMinvMinEtaEndCap_ = iConfig.getParameter<double> ("seleMinvMinEtaEndCap");  
    seleEtaBeltDREndCap_ = iConfig.getParameter<double> ("seleEtaBeltDREndCap");  
    seleEtaBeltDetaEndCap_ = iConfig.getParameter<double> ("seleEtaBeltDetaEndCap");  
    ptMinForIsolationEtaEndCap_  = iConfig.getParameter<double>("ptMinForIsolationEtaEndCap");
    
    seleEtaIsoEndCap_ = iConfig.getParameter<double> ("seleEtaIsoEndCap");  
    preshHitProducerEta_   = iConfig.getParameter<edm::InputTag>("preshRecHitProducerEta");
  }
  

  // input tag for L1GtTriggerMenuLite
  m_l1GtTmLInputTag = iConfig.getParameter<edm::InputTag> ("L1GtTmLInputTag"); 


  
  
  saveAllRecHitEB_ = iConfig.getParameter<bool> ("saveAllRecHitEB");
  saveAllRecHitEE_ = iConfig.getParameter<bool> ("saveAllRecHitEE");
  saveEGObj_= iConfig.getParameter<bool>("saveEGObj");
    
  
  
  
  ///for storing rechits ES for each selected EE clusters.
  storeRecHitES_ = iConfig.getParameter<bool>("storeRecHitES");  
  if(storeRecHitES_){
    // maximum number of matched ES clusters (in each ES layer) to each BC
    preshNclust_             = iConfig.getParameter<int>("preshNclust");
    // min energy of ES clusters
    preshClustECut = iConfig.getParameter<double>("preshClusterEnergyCut");
    // algo params
    double preshStripECut = iConfig.getParameter<double>("preshStripEnergyCut");
    int preshSeededNst = iConfig.getParameter<int>("preshSeededNstrip");
    // calibration parameters:
    calib_planeX_ = iConfig.getParameter<double>("preshCalibPlaneX");
    calib_planeY_ = iConfig.getParameter<double>("preshCalibPlaneY");
    gamma_        = iConfig.getParameter<double>("preshCalibGamma");
    mip_          = iConfig.getParameter<double>("preshCalibMIP");

    // The debug level
    std::string debugString = iConfig.getParameter<std::string>("debugLevelES");
    if      (debugString == "DEBUG")   debugL = PreshowerClusterAlgo::pDEBUG;
    else if (debugString == "INFO")    debugL = PreshowerClusterAlgo::pINFO;
    else                               debugL = PreshowerClusterAlgo::pERROR;
    // ES algo constructor:
    presh_algo = new PreshowerClusterAlgo(preshStripECut,preshClustECut,preshSeededNst,debugL);

    
  }
  
  InputDataFormat_ = iConfig.getParameter<int>("dataformat");  ///defautl input data
  
  
  //ParameterLogWeighted_ = iConfig.getParameter<bool> ("ParameterLogWeighted");
  //ParameterX0_ = iConfig.getParameter<double> ("ParameterX0");
  //  ParameterT0_barl_ = iConfig.getParameter<double> ("ParameterT0_barl");
  //ParameterT0_endc_ = iConfig.getParameter<double> ("ParameterT0_endc");
  //ParameterT0_endcPresh_ = iConfig.getParameter<double> ("ParameterT0_endcPresh");
  //ParameterW0_ = iConfig.getParameter<double> ("ParameterW0");
  
  
  theMinBunch = -10; 
  theMaxBunch = 10; 
  
  
  convLabel = iConfig.getUntrackedParameter<edm::InputTag>("convTrack",edm::InputTag("trackerOnlyConversions"));
  
  

  //SimHits
  hitsProducer_ = iConfig.getParameter<std::string>("hitsProducer");
  // initialize the default valuer for hardcoded parameters and the EB/EE shape

  bool addNoise = iConfig.getParameter<bool>("doNoise"); 
  bool applyConstantTerm = iConfig.getParameter<bool>("applyConstantTerm");
  double rmsConstantTerm = iConfig.getParameter<double> ("ConstantTerm");
  
  cout<<" applyNosie_ ConstantTerm "<<addNoise<<" "<< applyConstantTerm <<" "<< rmsConstantTerm <<endl;
  
  double simHitToPhotoelectronsBarrel = iConfig.getParameter<double>("simHitToPhotoelectronsBarrel");
  double simHitToPhotoelectronsEndcap = iConfig.getParameter<double>("simHitToPhotoelectronsEndcap");
  double photoelectronsToAnalogBarrel = iConfig.getParameter<double>("photoelectronsToAnalogBarrel");
  double photoelectronsToAnalogEndcap = iConfig.getParameter<double>("photoelectronsToAnalogEndcap");
  double samplingFactor = iConfig.getParameter<double>("samplingFactor");
  double timePhase = iConfig.getParameter<double>("timePhase");
  int readoutFrameSize = iConfig.getParameter<int>("readoutFrameSize");
  int binOfMaximum = iConfig.getParameter<int>("binOfMaximum");
  bool doPhotostatistics = iConfig.getParameter<bool>("doPhotostatistics");
  bool syncPhase = iConfig.getParameter<bool>("syncPhase");

  
  theParameterMap = new EcalSimParameterMap(simHitToPhotoelectronsBarrel, simHitToPhotoelectronsEndcap, 
					    photoelectronsToAnalogBarrel, photoelectronsToAnalogEndcap, 
					    samplingFactor, timePhase, readoutFrameSize, binOfMaximum,
                                            doPhotostatistics, syncPhase);
  
  //  theEcalShape = new EcalShape(timePhase);

  //theEcalResponse = new CaloHitResponse(theParameterMap, theEcalShape);
  
  
  double EBs25notCont = iConfig.getParameter<double>("EBs25notContainment");
  
  //  double pedrms12 = iConfig.getParameter<double>("pedWidth_gain12");
  
  cout<<"EBs25notCont: "<< EBs25notCont <<" "<< endl; 
  
  m_l1GtRecordInputTag = iConfig.getParameter<edm::InputTag>("L1GtRecordInputTag");
  
  ///from Raw only
  
  ecalDigiReFit_ = iConfig.getParameter<bool>("ecalDigiReFit");
  if( ecalDigiReFit_){
    ebDigiCollection_ = iConfig.getParameter<edm::InputTag>("EBdigiCollection");
    eeDigiCollection_ = iConfig.getParameter<edm::InputTag>("EEdigiCollection");
  }
  
  ebDigiCollectionv1_ = iConfig.getParameter<edm::InputTag>("ebDigiCollectionv1");
  
  
  debug_ = iConfig.getParameter<int> ("debugLevel");
  
  TheMapping = new EcalElectronicsMapping();
  first_ = true;
  

//   providedParameters.insert(std::make_pair("LogWeighted",ParameterLogWeighted_));
//   providedParameters.insert(std::make_pair("X0",ParameterX0_));
//   providedParameters.insert(std::make_pair("T0_barl",ParameterT0_barl_));
//   providedParameters.insert(std::make_pair("T0_endc",ParameterT0_endc_));
//   providedParameters.insert(std::make_pair("T0_endcPresh",ParameterT0_endcPresh_));
//   providedParameters.insert(std::make_pair("W0",ParameterW0_));
//   posCalculator_ = PositionCalc(providedParameters);
    
  posCalculator_ = PositionCalc( iConfig.getParameter<edm::ParameterSet>("posCalcParameters") );
  
  
  
  hfETowerh_ = iConfig.getUntrackedParameter<double>("hfETowerThreshold",3.); 
  
  
  ///HLT all bits, and physics decarled
  
  inputTag_    = iConfig.getParameter<edm::InputTag> ("TriggerResultsTag"); 
  andOr_        = iConfig.getParameter<bool> ("andOr"); 
  throw_        = iConfig.getParameter<bool> ("throw"); 
  eventSetupPathsKey_ = iConfig.getParameter<std::string>("eventSetupPathsKey"); 
  HLTPatterns_  = iConfig.getParameter<std::vector<std::string> >("HLTPaths"); 
  HLTPrescales_ = iConfig.getParameter<std::vector<uint32_t> >("HLTPathsPrescales"); 
  HLTOverallPrescale_ = iConfig.getParameter<uint32_t>("HLTOverallPrescale"); 
  
  if (eventSetupPathsKey_.size()) {
    // If paths come from eventsetup, we must watch for IOV changes.
    if (!HLTPatterns_.empty()) {
      // We do not want double trigger path setting, so throw!
      throw cms::Exception("Configuration")
        << " HLTHighLevelDev instance: "<< iConfig.getParameter<std::string>("@module_label")
        << "\n configured with " << HLTPatterns_.size() << " HLTPaths and\n"
        << " eventSetupPathsKey " << eventSetupPathsKey_ << ", choose either of them.";
    }
    watchAlCaRecoTriggerBitsRcd_ = new edm::ESWatcher<AlCaRecoTriggerBitsRcd>;
  }
  
  
  
  beamSpotInputTag_ = iConfig.getParameter<edm::InputTag> ("beamSpotInputTag"); 
  
  

  //register your products
  // doBarrel = true; 
 //  if(doSelForPi0Barrel_ && doSelForEtaBarrel_){
//     BarrelHits_ = pi0BarrelHits_;
//   }
//   else if(!doSelForPi0Barrel_ && doSelForEtaBarrel_){
//     BarrelHits_ = etaBarrelHits_;
//   }
//   else if(doSelForPi0Barrel_ && !doSelForEtaBarrel_){
//     BarrelHits_ = pi0BarrelHits_;
//   }else{
//     doBarrel = false; 
//   }
  
  
  
  // doEndcap = true; 

//   if(doSelForPi0Endcap_ && doSelForEtaEndcap_){
//     EndcapHits_ = pi0EndcapHits_;
//     ESHits_ = pi0ESHits_;
//   }else if(!doSelForPi0Endcap_ && doSelForEtaEndcap_){
//     EndcapHits_ = etaEndcapHits_;
//     ESHits_ = etaESHits_;
//   }else if(doSelForPi0Endcap_ && !doSelForEtaEndcap_){
//     EndcapHits_ = pi0EndcapHits_;
//     ESHits_ = pi0ESHits_;
//   }else{
//     doEndcap = false; 
//   }
  
  outputFile_   = iConfig.getParameter<std::string>("outputFile");
  rootFile_ = new TFile(outputFile_.c_str(),"RECREATE"); // open output file to store histograms
  
  std::cout<<"RecoAnalyzer initilized .." <<" "<<outputFile_.c_str()<<" opened."<<endl;
  
  std::cout<<"dataformat: "<< InputDataFormat_ <<" regional "<<RegionalMatch_ <<endl;
  alcaL1trigNames_ = iConfig.getParameter<vstring>("alcaL1trigNames");

  
  
  nEventsProcessed = 0; 
  
  nErrorPrinted = 0; 

  nPassedBSC = 0; 
  nPassedBSC_noBeamHalo = 0; 
  
  for(int j=0; j< 200; j++){
    nL1bits_fired[j] = 0; 
    nL1bitsTech_fired[j] = 0; 
  }

  curLumiBlock = -1; 
  
}


RecoAnalyzer::~RecoAnalyzer()
{
 
  //delete TheMapping;
  delete watchAlCaRecoTriggerBitsRcd_; // safe on null pointer...
  
  if(storeRecHitES_){
    delete presh_algo;
  }
  
  cout<<"totalProcessed: "<< nEventsProcessed<< " BSC4041 "<< nPassedBSC  <<" noBeamHalo "<< nPassedBSC_noBeamHalo <<endl; 
  //  TimingReport::current()->dump(std::cout);
  
}


// ------------ method called to produce the data  ------------
void
RecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
//RecoAnalyzer::analyze(edm::Event& iEvent, const edm::EventSetup& iSetup)
{

  
  if( nEventsProcessed ==0){
    l1algoName->clear();
    edm::ESHandle<L1GtTriggerMenu> menuRcd;
    iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
    const L1GtTriggerMenu* menu = menuRcd.product();
    cout <<" L1 Menu Name: "<< menuRcd->gtTriggerMenuName()<<endl;
    for (CItAlgo algo = menu->gtAlgorithmMap().begin(); algo!=menu->gtAlgorithmMap().end(); ++algo) {
      string tmp = (algo->second).algoName();
      l1algoName->push_back(tmp);
      
      cout<< "L1MenuPhy" << "Name: " << (algo->second).algoName() << " Alias: " << (algo->second).algoAlias() << std::endl;
    }
    for (CItAlgo techTrig = menu->gtTechnicalTriggerMap().begin(); techTrig != menu->gtTechnicalTriggerMap().end(); ++techTrig) {
      cout << "L1MenuTech" << (techTrig->second).algoName() << std::endl;
    }
    
    
    
    cout<<"getting ecalGainRatiosRcd.."<<endl;
    // Gain Ratios
    LogDebug("EcalUncalibRecHitDebug") << "fetching gainRatios....";
    iSetup.get<EcalGainRatiosRcd>().get(pRatio);
    cout<<"getting ecalGainRatiosRcd done"<<endl;
    
    cout<<"getting ecalPedestalRcd.."<<endl;
    // fetch the pedestals from the cond DB via EventSetup
    iSetup.get<EcalPedestalsRcd>().get( pedHandle );
    cout<<"getting ecalPedestalRcd done"<<endl;
    
    
    ///printout the icalconst
    edm::ESHandle<EcalIntercalibConstants> pIcal;
    const EcalIntercalibConstants *iical = 0;
    //print startup31x-v2 EB smeared constants
    iSetup.get<EcalIntercalibConstantsRcd>().get(pIcal);
    iical = pIcal.product();
    icalp = pIcal.product();
    
    const EcalIntercalibConstantMap &icalMap = iical->getMap();
    
    //   // Ecal Intercalibration Constants MC 
    //     edm::ESHandle<EcalIntercalibConstantsMC> pIcalmc;
    //     iSetup.get<EcalIntercalibConstantsMCRcd>().get(pIcalmc);
    //     const EcalIntercalibConstantsMC *icalmc = pIcalmc.product();
    //     const EcalIntercalibConstantMCMap &icalMapmc = icalmc->getMap();
    
    
    cout.setf(ios::fixed,ios::floatfield);
    cout.precision(8);
    for( int x = -85; x<= -85 ; x++){
      if( x ==0) continue; 
      for( int y = 1; y<= 3; y++){
	try{
	  EBDetId det = EBDetId(x,y,EBDetId::ETAPHIMODE);
	  EcalIntercalibConstant icalconst = 1.;
	  
	  EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(det);
	  if( icalit!=icalMap.end() ){
	    icalconst = (*icalit);
	    //  cout<<" ieb "<<x<<" "<<y<<" "<<icalconst<<endl; 
	  } else {
	    cout << "No intercalib const found for xtal " << det << "! something wrong with EcalIntercalibConstants in your DB? ";
	  }
	  EcalIntercalibConstantMC thisconstmc = 1.;
	  //   EcalIntercalibConstantMCMap::const_iterator icalitmc = icalMapmc.find(det);
	  // 	  if( icalitmc!=icalMapmc.end() ){
	  // 	    thisconstmc = (*icalitmc);
	  // 	  }
	  
	  ////cout<<"ebietaiphi "<< x <<" "<<y <<" "<<det.tower_ieta()<<" "<<det.tower_iphi()<<" "<<det.ism()<<endl; 
	  ////for idealtag, icalconst = icalconstmc =icalcosnt_preCalib
	  cout<<"ieb "<<x <<" "<< y <<" "<< icalconst <<" "<< thisconstmc <<endl;
	  
	}catch(...){}
      }
    }
    
    
    for( int j=0; j<0; j++){
      int iz = -1; 
      if( j==1) iz = 1; 
      
      for( int ix = 1; ix <=100; ix++){
	for( int iy =1; iy <=100; iy++){
	  try{
	    EEDetId det = EEDetId(ix,iy,iz,EEDetId::XYMODE);  
	    EcalIntercalibConstant icalconst = 1.;
	    EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(det);
	    if( icalit!=icalMap.end() ){
	      icalconst = (*icalit);
	      //  cout<<" ieb "<<x<<" "<<y<<" "<<icalconst<<endl; 
	    } else {
	      cout << "No intercalib const found for xtal " << det << "! something wrong with EcalIntercalibConstants in your DB? ";
	    }

	    EcalIntercalibConstantMC thisconstmc = 1.;
	    // EcalIntercalibConstantMCMap::const_iterator icalitmc = icalMapmc.find(det);
	    //if( icalitmc !=icalMapmc.end() ){
	    //  thisconstmc = (*icalitmc);
	    // }
	    
	    cout<<"iee: "<< iz <<" "<< ix <<" "<< iy <<" "<< icalconst <<" "<< endl; 
	    
	  }catch(...){}
	  
	}
      }
    }
    

    
    
    iSetup.get<EcalTimeCalibConstantsRcd>().get(itime);
    iSetup.get<EcalADCToGeVConstantRcd>().get(agc);
    agcp = agc.product();
    cout<<"agc: "<< agc->getEBValue() << " "<< agc->getEEValue() <<endl;
    
    
    
    
  }  
  
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  ////needed for each evetn 
  ///edm::ESHandle<EcalLaserDbService> laser;
  iSetup.get<EcalLaserDbRecord>().get(laser);
  
  
  edm::ESHandle<CaloGeometry> geoHandle;
  iSetup.get<CaloGeometryRecord>().get(geoHandle);
  

  if (first_) {
    edm::ESHandle< EcalElectronicsMapping > ecalmapping;
    iSetup.get< EcalMappingRcd >().get(ecalmapping);
    const EcalElectronicsMapping* TheMapping_ = ecalmapping.product();
    *TheMapping = *TheMapping_;
    first_ = false;

    
    //    const CaloGeometry * pGeometry = &*geoHandle;
    // theEcalResponse->setGeometry(pGeometry);
    
  }
  

//   ///for each detid, save the energy and time of simHits of barrel
//   std::map<DetId, vector<float> > ebId_simEn_map;
//   std::map<DetId, vector<float> > ebId_simTime_map;
  
//   //int nhit_sim[180][360];
//   // float en_sim[180][360][20]; 
//   // float time_sim[180][360][20]; 
//   for(int j=0; j< 180; j++){
//     for(int k=0; k<360; k++){
//       nhit_simeb[j][k] = 0; 
//     }
//   }
  
//   if(debug_ >=2) cout<<"getting CrossingFrame<PCaloHit> " <<endl;
//   // Get input
//   edm::Handle<CrossingFrame<PCaloHit> > crossingFrame;
//   if(debug_ >=2) cout<<"getting SimHitsEB " <<endl;
//   // test access to SimHits
//   const std::string barrelHitsName(hitsProducer_+"EcalHitsEB");
//   const std::string endcapHitsName(hitsProducer_+"EcalHitsEE");
  ///bool isEBHitValid = true;

  //iEvent.getByLabel("mix",barrelHitsName,crossingFrame);
  // MixCollection<PCaloHit> * EBHits = 0 ;
  //if (crossingFrame.isValid()) { 
    // EBHits = new MixCollection<PCaloHit>(crossingFrame.product());
    //std::auto_ptr<MixCollection<PCaloHit> >  barrelHits( EBHits );
    
    // if(debug_ >=2) cout<<"run ecalResponse barre.."<<endl;
    //theEcalResponse->run(*barrelHits);
    // if(debug_ >=2) cout<<"done ecalResponse barre.."<<endl;
    
  //    if(debug_>=1) cout<<"barreHits1: "<< barrelHits->size() <<endl;
  //    for(MixCollection<PCaloHit>::MixItr hitItr = barrelHits->begin();
  //	hitItr != barrelHits->end(); ++hitItr)
  //      {
	// check the bunch crossing range
  //	if ( hitItr.bunch() < theMinBunch || hitItr.bunch() > theMaxBunch )
  //	  { continue; }
  //	EBDetId ebid((*hitItr).id());
  //	int ieta = ebid.ieta(); 
  //	int iphi = ebid.iphi();
  //	convxtalid(iphi,ieta);
  //	if( nhit_simeb[ieta+85][iphi] < 20){
  //	  en_simeb[ieta+85][iphi][nhit_simeb[ieta+85][iphi]] = (*hitItr).energy(); 
  //	  time_simeb[ieta+85][iphi][nhit_simeb[ieta+85][iphi]] = (*hitItr).time(); 
  //	  nhit_simeb[ieta+85][iphi] ++; 
  //	}else{
  //	  cout<<"warning_nhit_simeb_bigger_than"<< nhit_simeb[ieta+85][iphi]<<" " << ieta <<" "<< iphi<<" "<< (*hitItr).energy()<<" "<< (*hitItr).time() <<endl;
  //	}
	
  //	if( (ebid.ieta()== -3 && ebid.iphi() == 70)   ||(ebid.ieta()== -4 && ebid.iphi() == 72) || (ebid.ieta()== -8 && ebid.iphi() == 65)
  //            || (ebid.ieta()== -11 && ebid.iphi() == 68)
  //            || (ebid.ieta()== -19 && ebid.iphi() == 115) ){
  //	  if(debug_>=2) cout<< "CaloSamples_ecaleb_checkHitEnergy: "<<ebid.ieta()<<" "<<ebid.iphi()<<" "<< (*hitItr).energy() <<endl;
  //	}
	
  //      }
  //  }
  //  else { 
  //  cout<< "Error! can't get the product " << barrelHitsName.c_str() <<endl;
  //  isEBHitValid = false;
  // }

  //  bool isEEHitValid = true;
  // iEvent.getByLabel("mix",endcapHitsName,crossingFrame);
  //MixCollection<PCaloHit> * EEHits = 0 ;
  //if (crossingFrame.isValid()) {
  // EEHits = new MixCollection<PCaloHit>(crossingFrame.product());
  //}
  //else {
  // cout<< "Error! can't get the product " << endcapHitsName.c_str() <<endl;
  // isEEHitValid = false;
  //}
  
  
//   Handle<EBDigiv1Collection> pEBDigisv1;
//   alldetIdEBDigiv1.clear();
  
//   if(ecalDigiReFit_){
//     try{
//       iEvent.getByLabel(ebDigiCollectionv1_, pEBDigisv1);
//       if ( pEBDigisv1.isValid() ) {
// 	ebDigisv1 = pEBDigisv1.product(); // get a ptr to the produc
// 	if(debug_ >=2) cout<< "EcalUncalibRecHitInfov1 " << "total # ebDigisv1: " << ebDigisv1->size() <<endl;
// 	for(EBDigiv1Collection::const_iterator itdg = ebDigisv1->begin(); itdg != ebDigisv1->end(); ++itdg) {
// 	  alldetIdEBDigiv1.push_back(itdg->id());
// 	}
	
//       }else{
// 	nErrorPrinted++;
// 	if(nErrorPrinted< maxErrorToPrint) cout<<"ebDigiCollectionv1_notValid "<<endl;
//       }
      
//     }catch(std::exception& ex ){
//       nErrorPrinted++;
//       if(nErrorPrinted< maxErrorToPrint) cout<<"ebDigiCollectionv1_NA.."<<endl;
//     }
    
//   }
  

//   Handle< EBDigiCollection > pEBDigis;
//   Handle< EEDigiCollection > pEEDigis;
    
//   alldetIdEBDigi.clear();
//   alldetIdEEDigi.clear();


  ///double adc_samp[20];
  //// double adc_sampv1[20];
  
 //  nEB = 0;
//   nEE = 0; 
  
  
//   int findEbRecHit = 0; 
//   const EcalRecHitCollection *hitCollection_p = 0 ; 
//   Handle<EBRecHitCollection> barrelRecHitsHandle;
//   try{
//     iEvent.getByLabel(barrelHits_,barrelRecHitsHandle);
//     hitCollection_p = barrelRecHitsHandle.product();
//     findEbRecHit  = 1; 
//   }catch(std::exception& ex ){
//     nErrorPrinted++;
//     if(nErrorPrinted< maxErrorToPrint) cout<<"barrelRecHits NA.."<<endl;
//   }
  
  
//   int findEbRecHitEta = 0; 
//   const EcalRecHitCollection *hitCollection_pEta = 0 ; 
//   Handle<EBRecHitCollection> barrelRecHitsHandleEta;
//   try{
//     iEvent.getByLabel(barrelHitsEta_,barrelRecHitsHandleEta);
//     hitCollection_pEta = barrelRecHitsHandleEta.product();
//     findEbRecHitEta  = 1; 
//   }catch(std::exception& ex ){
//     nErrorPrinted++;
//     if(nErrorPrinted< maxErrorToPrint) cout<<"barrelRecHitsEta NA.."<<endl;
//   }
    
  
//   int findEeRecHit = 0; 
//   Handle<EERecHitCollection> endcapRecHitsHandle;
//   const EcalRecHitCollection *hitCollection_e = 0 ; 
//   try{
//     iEvent.getByLabel(endcapHits_,endcapRecHitsHandle);
//     hitCollection_e = endcapRecHitsHandle.product();
//     findEeRecHit = 1; 
//   }catch(std::exception& ex ){
//     nErrorPrinted++;
//     if(nErrorPrinted< maxErrorToPrint) cout<<"endcapRecHits NA.."<<endl;
//   }
  
  
  

  
   nMCpart =0;

   
   if( InputDataFormat_ < 10){
     
     for(int i=0; i<MAXMC; i++){ ///MC
     pxMCpart[i] = 0;
     pyMCpart[i] = 0;
     pzMCpart[i] = 0;
     ptotMCpart[i] = 0;
     ptMCpart[i]=0; 
     etMCpart[i] = 0;
     eMCpart[i] = 0;
     etaMCpart[i] = -10;
     phiMCpart[i] = -10;
     pidMCpart[i] = 0;
     pidmomMCpart[i] = 0;
     vtxXMCpart[i] =0; 
     vtxYMCpart[i] =0; 
     vtxZMCpart[i] =0;
     barcodemomMCpart[i]=-1; 
     statusMCpart[i] =-1; 
     
     nDauMCpart[i] =0; 
     for( int j=0; j<MAXDAU; j++){
       barcodeDauMCpart[i][j] = -1; 
     }
    
     convPhtMCpart[i][0] = 0; 
     convPhtMCpart[i][1] = -1; 
     convPhtMCpart[i][2] = -1; 
     
     }
     
   }

   nGenpi0 = 0; 
   nGeneta = 0; 
   
   procID = -1; 
   ptHAT = -1; 
   
   npizallgen = 0; 
   netaallgen = 0; 
   
   nChaGen = 0; 
   
   
   

  runNumber      = iEvent.id().run();
  evtNumber    = iEvent.id().event();
  lumiBlock = iEvent.luminosityBlock();
  bunchX          = iEvent.bunchCrossing();
  orbitNumber       = iEvent.orbitNumber();
  const edm::Timestamp jtime = iEvent.time();
  evtTime = jtime.value() >> 32;
  
  isNewLumiBlock = false; 
  
  if( curLumiBlock != lumiBlock){
    curLumiBlock = lumiBlock; 
    isNewLumiBlock = true; 
  }

  if( getBeamSpotOnly_) { ///new lumiBlock
    if( !isNewLumiBlock ) return; 
  }
  
  
  ///let's get the gen pi0/eta
  if( InputDataFormat_ < 10){
    
    Handle<GenEventInfoProduct> geninfos; 
    try{
      iEvent.getByLabel("generator",geninfos);
      
      procID = int( geninfos->signalProcessID());
      ///this is q, not pthat.
      ///      ptHAT = geninfos->qScale();
      ptHAT = geninfos->binningValues()[0]; //this is ptHAT
      
      
    }catch(std::exception& ex ){
      nErrorPrinted++;
      if(nErrorPrinted< maxErrorToPrint) cout<<"generator not working.."<<endl;
    }

    if(debug_>=2) cout<<"fill MCpart: "<<endl; 
    

    Handle<GenParticleCollection> genPart;
    try {
      //iEvent.getByLabel( "genParticleCandidates", genParticles );
      iEvent.getByLabel( "genParticles", genPart);
      ///const HepMC::GenEvent* ievt = EvtHandle->GetEvent() ;
      nMCpart = 0; 
      for(int i=0; i< int(genPart->size())&& nMCpart<MAXMC; i++){
	const GenParticle & p = (*genPart)[i];
	eMCpart[nMCpart] = p.energy();
	ptotMCpart[nMCpart] = p.p();
	pxMCpart[nMCpart] =p.px();
	pyMCpart[nMCpart] =p.py();
	pzMCpart[nMCpart] =p.pz();
	ptMCpart[nMCpart] =p.pt();
	etMCpart[nMCpart] = p.et();
	etaMCpart[nMCpart] = p.eta();
	phiMCpart[nMCpart] = p.phi();
	mMCpart[nMCpart] = p.mass();
	pidMCpart[nMCpart] = p.pdgId();
	chargeMCpart[nMCpart] = p.charge();
	statusMCpart[nMCpart] = p.status();
	nDauMCpart[nMCpart] = p.numberOfDaughters();
	vtxXMCpart[nMCpart] = p.vx();
	vtxYMCpart[nMCpart] = p.vy();
	vtxZMCpart[nMCpart] = p.vz();

	

	if(p.charge() !=0 && p.status() ==1 ){
	  //	  cout<<"charge: "<< nMCpart <<" "<< p.pt()<<" "<< p.eta()<<" "<<p.phi()<<" "<<p.pdgId()<<" "<< p.charge()<<" "<< p.status()<<endl; 
	  ptChaGen[nChaGen] = p.pt();
	  etaChaGen[nChaGen] = p.eta();
	  pidChaGen[nChaGen] = p.pdgId();
	  chaChaGen[nChaGen] = p.charge();
	  staChaGen[nChaGen] = p.status();
	  nChaGen ++; 
	} 
	

	nMCpart++;
      }
      if(debug_>=1) cout<<"n total MCpart: "<<nMCpart<<endl;
      
        ///let's try to find the index of mother and daughters.
      for(int i=0; i< int(genPart->size())&& nMCpart<MAXMC; i++){
	const Candidate & p = (*genPart)[i];
	///don't look if it's not photon
	//	if( p.pdgId() != 22) continue; 
	int pid = p.pdgId();
        if( pid != 22 && pid != 111 && pid !=  221 ) continue;
	const Candidate * mom = p.mother();
	if( mom !=NULL){
	  pidmomMCpart[i] = mom->pdgId();
	  barcodemomMCpart[i] = indexofParticle(mom->px(),mom->pz(),mom->status());
	}
      }
      
      
    
    } catch( std::exception& ex){
      nErrorPrinted++;
      if(nErrorPrinted< maxErrorToPrint) cout<<"genparticle  not working.."<<endl;
      
    }
      // Simulated tracks (i.e. GEANT particles).
    Handle<edm::SimTrackContainer> psimtracks; 
    // Get the associated vertices
    Handle<SimVertexContainer> simvertices;
    try{
      iEvent.getByLabel("g4SimHits", simvertices);
      
    }catch( std::exception& ex){
      nErrorPrinted++;
      if(nErrorPrinted< maxErrorToPrint) cout<<" SimVertexContainer not working.."<<endl;
    
    }
  
  
    try{
      iEvent.getByLabel("g4SimHits", psimtracks);
    
      
      if(debug_ >=1) cout<<"get g4SimHits..."<<endl;
      
      //const edm::SimTrackContainer simTracks = *(psimtracks.product());
    
      // cout<<"g4simHitSize: "<<psimtracks->size()<<endl;
    
      // Need to check that SimTrackContainer is sorted; otherwise, copy and sort :-(
      std::auto_ptr<SimTrackContainer> simtracksTmp;

      if(debug_ >=1) cout<<"get g4SimHits..get1."<<endl;
      const SimTrackContainer * simtracksSorted = &* psimtracks;

      if(debug_ >=1) cout<<"get g4SimHits..get2."<<endl;
 

      if (!__gnu_cxx::is_sorted(psimtracks->begin(), psimtracks->end(), LessById())) {
	simtracksTmp.reset(new SimTrackContainer(*psimtracks));
	std::sort(simtracksTmp->begin(), simtracksTmp->end(), LessById());
	if(debug_ >=1) cout<<"get g4SimHits..sorted."<<endl;
	
	simtracksSorted = &* simtracksTmp;
      }
            
      
      int nMCpartGen = nMCpart; 
      
      
      map<unsigned int,int> map_simtrackid; 
    
      // loop through all simParticles
      for (SimTrackContainer::const_iterator iM = psimtracks->begin(); 
	   iM != psimtracks->end(); ++iM) {
      
	// Skip PYTHIA tracks.
	if (iM->genpartIndex() != -1) continue; 
      
	/// skip if it is not electron
	if( iM->type() != 11 &&  iM->type() != -11) continue; 
      
	
	///int indp = iM->genpartIndex()-1; 
      
	
	if( nMCpart >= MAXMC) break; 
	
	map_simtrackid.insert(make_pair(iM->trackId(),nMCpart));
      
	pxMCpart[nMCpart] = iM->momentum().px();
	pyMCpart[nMCpart] = iM->momentum().py();
	pzMCpart[nMCpart] = iM->momentum().pz();
	ptMCpart[nMCpart] = iM->momentum().pt();
	etaMCpart[nMCpart] = iM->momentum().eta();
	phiMCpart[nMCpart] = iM->momentum().phi();
	pidMCpart[nMCpart] = iM->type();
      
	chargeMCpart[nMCpart] = int(iM->charge());
	//      pidFDauMCpart[nMCpart] = 0; 
      
	eMCpart[nMCpart] = sqrt(0.000511* 0.000511 + iM->momentum().px()*iM->momentum().px() + iM->momentum().py()*iM->momentum().py() + iM->momentum().pz()*iM->momentum().pz()); 
	      
	statusMCpart[nMCpart] = -99; 
      
	vtxXMCpart[nMCpart] = 0;
	vtxYMCpart[nMCpart] = 0;
	vtxZMCpart[nMCpart] = 0;
      
	if (! iM->noVertex()) {
	  const SimVertex &vtx = (*simvertices)[iM->vertIndex()];
	  vtxXMCpart[nMCpart] = vtx.position().x(); 
	  vtxYMCpart[nMCpart] = vtx.position().y(); 
	  vtxZMCpart[nMCpart] = vtx.position().z(); 
	}
	
	barcodemomMCpart[nMCpart] = -2; 
	nMCpart ++; 
      
      }
      
      
      vector<int> inde1;  ///+ /
      vector<int> inde2; ///- 
      vector<int> indpht; 
      
      
      
      int nSimTrack = 0; 
    
      ///loop over again
      // loop through all simParticles
      for (SimTrackContainer::const_iterator iM = psimtracks->begin(); 
	   iM != psimtracks->end(); ++iM) {
      
	// Skip PYTHIA tracks.
	if (iM->genpartIndex() != -1) continue; 

	/// skip if it is not electron
	if( iM->type() != 11 &&  iM->type() != -11) continue; 
		
      
	if( nMCpartGen + nSimTrack >= MAXMC) break; 
	
	
	barcodemomMCpart[nMCpartGen+nSimTrack] = -2; 
      
	
	// find simtrack that has a genParticle match to its parent
	// Look at the production vertex. If there is no vertex, I can do nothing...
	if (! iM->noVertex()) {
	  // Pick the vertex (isimtrk.vertIndex() is really an index)
	  const SimVertex &vtx = (*simvertices)[iM->vertIndex()];
	
	  // Now note that vtx.parentIndex() is NOT an index, it's a track id, so I have to search for it 
	  unsigned int idx = vtx.parentIndex();
	  SimTrackContainer::const_iterator it = std::lower_bound(simtracksSorted->begin(), simtracksSorted->end(), idx, LessById());
	
	  if ((it != simtracksSorted->end()) && (it->trackId() == idx)) { //it is the parent sim track
	    if (it->genpartIndex() != -1) {
	    
	      barcodemomMCpart[nMCpartGen+nSimTrack] = it->genpartIndex() -1; 
	    
	      convPhtMCpart[it->genpartIndex() -1][0] = 1; ///this gen photon converts
	    

	      vector<int>::iterator iit = find(indpht.begin(),indpht.end(),it->genpartIndex() -1);
	      if( iit == indpht.end()){
		indpht.push_back(it->genpartIndex() -1); 
	      }
	      if( iM->type() == 11) inde1.push_back(nMCpartGen+nSimTrack); 
	      else inde2.push_back(nMCpartGen+nSimTrack); 
	    

	      ////  cout<<"found efromg: "<<" "<<evtNumber<<" "<<iM->type()<<" e/eta/phi " <<eMCpart[nMCpartGen+nSimTrack]<<" "<<ptMCpart[nMCpartGen+nSimTrack]<<" "<<iM->momentum().eta()<<" "<<iM->momentum().phi()<<" mom_ind: "<<eMCpart[it->genpartIndex()-1]<<" "<<etaMCpart[it->genpartIndex()-1]<<" "<<phiMCpart[it->genpartIndex()-1]<<endl;
	    
	    }else{
	      ///	    cout<<"no genPartIndexOfParent_of_nogenpartIndex.._"<<endl;
	      map<unsigned int,int>::iterator intt = map_simtrackid.find(idx);
	      if( intt != map_simtrackid.end()){
	      
		barcodemomMCpart[nMCpartGen+nSimTrack]  = intt->second;
	      }else{
		barcodemomMCpart[nMCpartGen+nSimTrack] = -3; ////check this if any
		//	      cout<<"warninig.. noSimIndexOfParent_of_nogenpartIndex.."<<endl;
	      }
	    
	    }
	  } ///found parent in simtrack
	
	}
	
	nSimTrack ++; 
	
      
      } ///end of 2nd loop 
      

      ///now let's try to find for each indpht, it's correponding two electrons index
      vector<int> pht;
      vector<int> e1;
      vector<int> e2;
      

      if(debug_ >= 2) cout<<" indpht_conv: "<< int( indpht.size()) <<endl; 

      for(int j =0; j< int(indpht.size());j++){
	int k = indpht[j];
      
	int found1 = -1;
	int found2 = -1;
	vector<int>::iterator it1= inde1.begin();
	for(int j1 =0; j1< int(inde1.size()); j1++){
	  if(k ==   barcodemomMCpart[inde1[j1]]){
	    found1 = inde1[j1];
	    inde1.erase(it1+j1);
	    break;
	  }
	}
	vector<int>::iterator it2= inde2.begin();
	for(int j2 =0; j2< int(inde2.size()); j2++){
	  if(k ==   barcodemomMCpart[inde2[j2]]){
	    found2 = inde2[j2];
	    inde2.erase(it2+j2);
	    break;
	  }
	}
      
	convPhtMCpart[k][1] = found1; 
	convPhtMCpart[k][2] = found2; 
	
	if( debug_ >=2) cout<<"convPhtMCpart: "<< k<<" " << found1 <<" "<< found2<<" "<< eMCpart[k] << endl; 
	
	if( found1 > 0 && found2 >0) convPhtMCpart[k][0] = 2; /// both electron found in sim
	else if( found1 >0 && found2 < 0) convPhtMCpart[k][0] = 3; /// only one found in sim 
	else if( found1 <0 && found2 > 0) convPhtMCpart[k][0] = 3; /// only one found in sim 
	else if( found1 <0 && found2 < 0) convPhtMCpart[k][0] = -10; /// neither found in sim 
            
      }
      
    
    }catch( std::exception& ex){
      nErrorPrinted++;
      if(nErrorPrinted< maxErrorToPrint) cout<<"SimTrackContainer  not working.."<<endl;
      
    }
  
    findgenpi0eta();
    
    
  } ///If for MC samples
  
  
  
  // get vertex collection offlinePrimaryVerticesWithBS
  
  if(fullRECO_==1){
    if(debug_ >=1) {
      cout<<" check vertex.."<<endl;
    }
    
    
    Handle<VertexCollection> vertices;
    try{
      iEvent.getByLabel(m_vertexSrc, vertices);
      nVertex = 0; 
      for (int j = 0; j < int(vertices->size()) && j < MAXVX; j++){
	vertexx[nVertex] = (*vertices)[j].x();
	vertexy[nVertex] = (*vertices)[j].y();
	vertexz[nVertex] = (*vertices)[j].z();
	vertexchi2[nVertex] = (*vertices)[j].chi2();
	vertexndof[nVertex] = (*vertices)[j].ndof(); 
	vertexnormalizedChi2[nVertex] = (*vertices)[j].normalizedChi2(); 
	vertextrackSize[nVertex] = (*vertices)[j].tracksSize(); 
	vertexisFake[nVertex] = (*vertices)[j].isFake();
	vertexisValid[nVertex] = (*vertices)[j].isValid();
	nVertex ++;
      }
    }catch( std::exception& ex){
      nErrorPrinted++;
      if(nErrorPrinted< maxErrorToPrint) cout<<" vertex not working.."<<endl;
    }
    
  }
  
  
  ///now save beamSpot for each pair as well
  ///but the beamSpot is the same within each lumiBlock of each run
  /// so to save space, we should save only the beamspot information when a new lumiBlock or a new run come
  
  //  if( getBeamSpotOnly_){

  if( isNewLumiBlock){
    reco::TrackBase::Point beamPoint(0,0, 0);
    vBeamSpot[0] = vBeamSpot[1] = vBeamSpot[2] = 0; 
    reco::BeamSpot beamSpot;
    edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
    try{
      iEvent.getByLabel(beamSpotInputTag_,recoBeamSpotHandle);
      beamSpot = *recoBeamSpotHandle;
      vBeamSpot[0] = beamSpot.position().X();
      vBeamSpot[1] = beamSpot.position().Y();
      vBeamSpot[2] = beamSpot.position().Z();
      if(debug_ >=2) cout<<"beamPoint: "<< beamSpot.position().X() <<" "<< beamSpot.position().Y() <<" "<< beamSpot.position().Z() <<endl; 
    }catch(std::exception& ex ){
      nErrorPrinted++;
      if(nErrorPrinted< maxErrorToPrint) cout<<"beamSpot not working.."<<endl;
    }
    
    mytree_bs->Fill(); 
    
    
    if( getBeamSpotOnly_){
      mytree_evtInfo->Fill();
      nEventsProcessed ++; 
      return; 
    }
    
  }
  
  
  
  int    ePxHits      = 0;
  double eClusVtxQual = 0;
  double eClusVtxDiff = 0;
  int    ePxbHits     = 0;
  
  nHfTowersP      = 0; ///save into root tree
  nHfTowersN      = 0; ///save into root tree
  
  double sumEsubEpPos    = 0;
  double sumEaddEpPos    = 0;
  double sumHfEsubEpPos  = 0;
  double sumHfEaddEpPos  = 0;
  double sumEsubEpNeg    = 0;
  double sumEaddEpNeg    = 0;
  double sumHfEsubEpNeg  = 0;
  double sumHfEaddEpNeg  = 0;
  //double eHPTrkFrac      = 0;
  
  
  
  

  
  
  // get L1GlobalTriggerReadoutRecord
  edm::Handle<L1GlobalTriggerReadoutRecord> gtReadoutRecord;

  m_l1GtUtils.retrieveL1EventSetup(iSetup);
  m_l1GtUtils.retrieveL1GtTriggerMenuLite(iEvent, m_l1GtTmLInputTag);

  if(debug_ >=1) {
    cout<<" check L1.."<<endl;
  }
  
  ///    std::string m_nameAlgTechTrig = "L1_SingleEG2";
  ///vector<string> listOfAlCaL1Trigs; 
  int iErrorCode = -1;
  nL1Alca = int(alcaL1trigNames_.size());
  for(int j=0; j< nL1Alca; j++){
    string m_nameAlgTechTrig = alcaL1trigNames_[j];
    bool decisionBeforeMaskAlgTechTrig = m_l1GtUtils.decisionBeforeMask(iEvent, m_nameAlgTechTrig, iErrorCode);
    if (iErrorCode == 0) {
      // do something
      L1Alca[j] = int(decisionBeforeMaskAlgTechTrig); 
    } else if (iErrorCode == 1) {
      L1Alca[j] = -1; 
      // algorithm / technical trigger  does not exist in the L1 menu
      // do something
    } else {
      L1Alca[j] = -2; 
      // error - see error code
      // do whatever needed
    }
    
  }
  
  if(debug_ >=1) {
    cout<<" check HLT.."<<endl;
  }
  

  if (fullRECO_){
    bool changed = false;
    hltConfig_.init(iEvent.getRun(),iSetup,inputTag_.process(),changed);
    edm::Handle<edm::TriggerResults> h_triggerResults_HLT1;
    iEvent.getByLabel(inputTag_, h_triggerResults_HLT1);
    hlt_pathName->clear();
    hlt_bitFired->clear();
    if (h_triggerResults_HLT1.isValid()) {
      for (size_t i = 0; i < hltConfig_.size(); ++i){
	hlt_pathName->push_back(hltConfig_.triggerName(i));
	if(h_triggerResults_HLT1->accept(i)){
	  hlt_bitFired->push_back((unsigned short)(i));
	}
      }
    }
  }    
  
  edm::ESHandle<L1GtTriggerMenu> menuRcd;
  iSetup.get<L1GtTriggerMenuRcd>().get(menuRcd) ;
  const L1GtTriggerMenu* menu = menuRcd.product();
  //l1bitFired = new std::vector<unsigned short>; l1bitFired->clear();
  //l1algoName = new std::vector<std::string>; l1algoName->clear();
  l1bitFired->clear();
  
  for(unsigned int j=0; j< l1algoName->size(); j++){
    string m_nameAlgTechTrig = l1algoName->at(j);
    bool decisionBeforeMaskAlgTechTrig = m_l1GtUtils.decisionBeforeMask(iEvent, m_nameAlgTechTrig, iErrorCode);
    if (iErrorCode == 0) {
      if(decisionBeforeMaskAlgTechTrig){
	l1bitFired->push_back(j);
      }
    }else if (iErrorCode == 1) {
      //L1bits[j] = -1;
    }else {
      ///L1bits[j] = -2;
    }
  }
  if(debug_>=2) cout<<"l1algo " << l1bitFired->size()<<"/"<<l1algoName->size()<<endl;
  
  if( nEventsProcessed % 10000 ==0) {
    cout<<"nEventProcessed: "<< nEventsProcessed<< endl; 
  }
  
  nEventsProcessed ++; 
  
  //phyDeclared = filter(iEvent, iSetup);  //save this into root tree
  
  phyDeclared = 1; 
  
  
  
  ///

  if( doSelForPi0Barrel_ || saveAllPhotonBarrel_ ){
    
  try{
    Handle<EBRecHitCollection> barrelRecHitsHandle;
    iEvent.getByLabel(barrelHits_,barrelRecHitsHandle);
    const EcalRecHitCollection *hitCollection_p = barrelRecHitsHandle.product();
        
    if(debug_>=2) cout<<"making cluster pi0 barrel.." <<" "<< hitCollection_p->size()<<endl; 

    makeNxNClusters(iEvent,iSetup,hitCollection_p, reco::CaloID::DET_ECAL_BARREL);

    if( doSelForPi0Barrel_){
      if(debug_>=2) cout<<"do selection pi0 barrel.." <<endl; 
      
      doSelectionAndFillTree(iEvent,iSetup,reco::CaloID::DET_ECAL_BARREL,false,PIZ,seleS4S9Gamma_,-999,seleMinvMinPi0_,seleMinvMaxPi0_,selePi0BeltDeta_,selePi0BeltDR_,ptMinForIsolation_,selePi0Iso_);
    }
    
  }catch(std::exception& ex ){
    nErrorPrinted++;
    if(nErrorPrinted< maxErrorToPrint) cout<<"barrelRecHits NA.."<<endl;
  }
  
  }
  
  
  if(doSelForEtaBarrel_){
    try{
    Handle<EBRecHitCollection> barrelRecHitsHandle;
    iEvent.getByLabel(barrelHitsEta_,barrelRecHitsHandle);
    const EcalRecHitCollection *hitCollection_p = barrelRecHitsHandle.product();
    
    if(debug_>=2) cout<<"making cluster eta barrel.."<<" "<< hitCollection_p->size()<<endl; 

    makeNxNClusters(iEvent,iSetup,hitCollection_p, reco::CaloID::DET_ECAL_BARREL);
    
    if(debug_>=2) cout<<"do selection eta barrel.." <<endl; 

    doSelectionAndFillTree(iEvent,iSetup,reco::CaloID::DET_ECAL_BARREL,removePi0CandidatesForEta_,ETA,seleS4S9GammaEta_,seleS9S25GammaEta_,seleMinvMinEta_,seleMinvMaxEta_,seleEtaBeltDeta_,seleEtaBeltDR_,ptMinForIsolationEta_,seleEtaIso_);
    
    
  }catch(std::exception& ex ){
    nErrorPrinted++;
    if(nErrorPrinted< maxErrorToPrint) cout<<"barrelRecHitsEta NA.."<<endl;
  }
  }
  
  
  if(doSelForPi0Endcap_ || saveAllPhotonBarrel_ ){
  try{
    Handle<EERecHitCollection> endcapRecHitsHandle;
    iEvent.getByLabel(endcapHits_,endcapRecHitsHandle);
    const EcalRecHitCollection *hitCollection_p = endcapRecHitsHandle.product();
    
    makeNxNClusters(iEvent,iSetup,hitCollection_p, reco::CaloID::DET_ECAL_ENDCAP);
    
    if(debug_>=2) cout<<"do selection pi0  endcap.." <<endl; 
    if(doSelForPi0Endcap_){
      doSelectionAndFillTree(iEvent,iSetup,reco::CaloID::DET_ECAL_ENDCAP,false,PIZ,seleS4S9GammaEndCap_,-999,seleMinvMinPi0EndCap_,seleMinvMaxPi0EndCap_,selePi0BeltDetaEndCap_,selePi0BeltDREndCap_,ptMinForIsolationEndCap_,selePi0IsoEndCap_);
    }
    
  }catch(std::exception& ex ){
    nErrorPrinted++;
    if(nErrorPrinted< maxErrorToPrint) cout<<"encapRecHits NA.."<<endl;
  }
  }
  
  if(doSelForEtaEndcap_){
  try{
    Handle<EERecHitCollection> endcapRecHitsHandle;
    iEvent.getByLabel(endcapHitsEta_,endcapRecHitsHandle);
    const EcalRecHitCollection *hitCollection_p = endcapRecHitsHandle.product();
    
    makeNxNClusters(iEvent,iSetup,hitCollection_p, reco::CaloID::DET_ECAL_ENDCAP);
    
    if(debug_>=2) cout<<"do selection eta  endcap.." <<endl; 
    
    
    doSelectionAndFillTree(iEvent,iSetup,reco::CaloID::DET_ECAL_ENDCAP,false,ETA,seleS4S9GammaEtaEndCap_,seleS9S25GammaEtaEndCap_,seleMinvMinEtaEndCap_,seleMinvMaxEtaEndCap_,seleEtaBeltDetaEndCap_,seleEtaBeltDREndCap_,ptMinForIsolationEtaEndCap_,seleEtaIsoEndCap_);
    
    
  }catch(std::exception& ex ){
    nErrorPrinted++;
    if(nErrorPrinted< maxErrorToPrint) cout<<"endcapRecHitsEta NA.."<<endl;
  }
  }
  
  
  if( fullRECO_   ){
    mytree_evtInfo->Fill();
  }
  
}



///void 
//RecoAnalyzer::beginJob(const edm::EventSetup& iSetup)
void RecoAnalyzer::beginJob()
///void RecoAnalyzer::beginRun(edm::Run const & iRun, edm::EventSetup const& iSetup)
{  
  rootFile_->cd();

  
  
  l1bitFired = new std::vector<unsigned short>; l1bitFired->clear();                         
  l1algoName = new std::vector<std::string>; l1algoName->clear();   
  hlt_bitFired = new std::vector<unsigned short>; hlt_bitFired->clear();
  hlt_pathName = new std::vector<std::string>; hlt_pathName->clear();
  
  
  

  mytree_bs = new TTree("evtInfobs","some info from reco");
  mytree_bs->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
  mytree_bs->Branch("runNumber",&runNumber,"runNumber/I");
  mytree_bs->Branch("evtNumber",&evtNumber,"evtNumber/I");
  mytree_bs->Branch("evtTime",&evtTime,"evtTime/I");
  mytree_bs->Branch("bunchX",&bunchX,"bunchX/I");                                            
  mytree_bs->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
  mytree_bs->Branch("vBeamSpot",vBeamSpot,"vBeamSpot[3]/F");
  mytree_bs->Branch("l1algoName","std::vector<std::string>",&l1algoName);
    
  if( fullRECO_    || getBeamSpotOnly_ ){
    mytree_evtInfo = new TTree("evtInfo","some info from reco");
    
    mytree_evtInfo->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
    mytree_evtInfo->Branch("runNumber",&runNumber,"runNumber/I");
    mytree_evtInfo->Branch("evtNumber",&evtNumber,"evtNumber/I");
    mytree_evtInfo->Branch("evtTime",&evtTime,"evtTime/I");
    mytree_evtInfo->Branch("bunchX",&bunchX,"bunchX/I");                                            
    mytree_evtInfo->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
    mytree_evtInfo->Branch("vBeamSpot",vBeamSpot,"vBeamSpot[3]/F");
    
    if( getBeamSpotOnly_) return; 
    
    
    mytree_evtInfo->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
    mytree_evtInfo->Branch("l1algoName","std::vector<std::string>",&l1algoName);
    
    mytree_evtInfo->Branch("hlt_bitFired", "std::vector<unsigned short>", &hlt_bitFired); 
    ///vector  of those bit which fired the trigger                                             
    mytree_evtInfo->Branch("hlt_pathName", "std::vector<std::string>", &hlt_pathName);/// vector of       string for all trigger path       
    
    
    //mytree_evtInfo->Branch("mpairv1",&mpairv1,"mpairv1/F"); ///using PV 
    ///vertex all svaed
    mytree_evtInfo->Branch("nVertex",&nVertex,"nVertex/I");
    mytree_evtInfo->Branch("vertexx",vertexx,"vertexx[nVertex]/F");
    mytree_evtInfo->Branch("vertexy",vertexy,"vertexy[nVertex]/F");
    mytree_evtInfo->Branch("vertexz",vertexz,"vertexz[nVertex]/F");
    mytree_evtInfo->Branch("vertexndof",vertexndof,"vertexndof[nVertex]/F");
    mytree_evtInfo->Branch("vertexisFake",vertexisFake,"vertexisFake[nVertex]/I");
    mytree_evtInfo->Branch("vertexisValid",vertexisValid,"vertexisValid[nVertex]/I");
    
    
    ///clusters
    if( saveAllPhotonBarrel_){
      
      mytree_evtInfo->Branch("n3x3ClusEB",&n3x3ClusEB,"n3x3ClusEB/I");
      mytree_evtInfo->Branch("e3x3ClusEB",e3x3ClusEB,"e3x3ClusEB[n3x3ClusEB]/F");
      //mytree_evtInfo->Branch("eta3x3ClusEB",eta3x3ClusEB,"eta3x3ClusEB[n3x3ClusEB]/F");
      //mytree_evtInfo->Branch("phi3x3ClusEB",phi3x3ClusEB,"phi3x3ClusEB[n3x3ClusEB]/F");
      mytree_evtInfo->Branch("x3x3ClusEB",x3x3ClusEB,"x3x3ClusEB[n3x3ClusEB]/F");
      mytree_evtInfo->Branch("y3x3ClusEB",y3x3ClusEB,"y3x3ClusEB[n3x3ClusEB]/F");
      mytree_evtInfo->Branch("z3x3ClusEB",z3x3ClusEB,"z3x3ClusEB[n3x3ClusEB]/F");
      mytree_evtInfo->Branch("nXt3x3ClusEB",nXt3x3ClusEB,"nXt3x3ClusEB[n3x3ClusEB]/I");
      mytree_evtInfo->Branch("eXt3x3ClusEB",eXt3x3ClusEB,"eXt3x3ClusEB[n3x3ClusEB][9]/F");
      mytree_evtInfo->Branch("tXt3x3ClusEB",tXt3x3ClusEB,"tXt3x3ClusEB[n3x3ClusEB][9]/F");
      mytree_evtInfo->Branch("ietaXt3x3ClusEB",ietaXt3x3ClusEB,"ietaXt3x3ClusEB[n3x3ClusEB][9]/I");
      mytree_evtInfo->Branch("iphiXt3x3ClusEB",iphiXt3x3ClusEB,"iphiXt3x3ClusEB[n3x3ClusEB][9]/I");
      mytree_evtInfo->Branch("laserCorr3x3ClusEB",laserCorr3x3ClusEB,"laserCorr3x3ClusEB[n3x3ClusEB][9]/F");
      ///mytree_evtInfo->Branch("nXt5x5ClusEB",nXt5x5ClusEB,"nXt5x5ClusEB[n3x3ClusEB]/I"); ///additional 5x5 
      
      mytree_evtInfo->Branch("s4s93x3ClusEB",s4s93x3ClusEB,"s4s93x3ClusEB[n3x3ClusEB]/F");
      mytree_evtInfo->Branch("s9s253x3ClusEB",s9s253x3ClusEB,"s9s253x3ClusEB[n3x3ClusEB]/F");
      
      
    }
    
    if( saveAllPhotonEndcap_){
      mytree_evtInfo->Branch("n3x3ClusEE",&n3x3ClusEE,"n3x3ClusEE/I");
      mytree_evtInfo->Branch("e3x3ClusEE",e3x3ClusEE,"e3x3ClusEE[n3x3ClusEE]/F");
      //   mytree_evtInfo->Branch("eta3x3ClusEE",eta3x3ClusEE,"eta3x3ClusEE[n3x3ClusEE]/F");
      //       mytree_evtInfo->Branch("phi3x3ClusEE",phi3x3ClusEE,"phi3x3ClusEE[n3x3ClusEE]/F");
      mytree_evtInfo->Branch("x3x3ClusEE",x3x3ClusEE,"x3x3ClusEE[n3x3ClusEE]/F");
      mytree_evtInfo->Branch("y3x3ClusEE",y3x3ClusEE,"y3x3ClusEE[n3x3ClusEE]/F");
      mytree_evtInfo->Branch("z3x3ClusEE",z3x3ClusEE,"z3x3ClusEE[n3x3ClusEE]/F");
      mytree_evtInfo->Branch("nXt3x3ClusEE",nXt3x3ClusEE,"nXt3x3ClusEE[n3x3ClusEE]/I");
      mytree_evtInfo->Branch("eXt3x3ClusEE",eXt3x3ClusEE,"eXt3x3ClusEE[n3x3ClusEE][9]/F");
      mytree_evtInfo->Branch("tXt3x3ClusEE",tXt3x3ClusEE,"tXt3x3ClusEE[n3x3ClusEE][9]/F");
      mytree_evtInfo->Branch("ixXt3x3ClusEE",ixXt3x3ClusEE,"ixXt3x3ClusEE[n3x3ClusEE][9]/I");
      mytree_evtInfo->Branch("iyXt3x3ClusEE",iyXt3x3ClusEE,"iyXt3x3ClusEE[n3x3ClusEE][9]/I"); 
      //mytree_evtInfo->Branch("izXt3x3ClusEE",izXt3x3ClusEE,"izXt3x3ClusEE[n3x3ClusEE][9]/I"); 
      mytree_evtInfo->Branch("laserCorr3x3ClusEE",laserCorr3x3ClusEE,"laserCorr3x3ClusEE[n3x3ClusEE][9]/F");      
      mytree_evtInfo->Branch("s4s93x3ClusEE",s4s93x3ClusEE,"s4s93x3ClusEE[n3x3ClusEE]/F");
      mytree_evtInfo->Branch("s9s253x3ClusEE",s9s253x3ClusEE,"s9s253x3ClusEE[n3x3ClusEE]/F");
      
      
    }
    
    
  }
  
  


  
    if( doSelForPi0Barrel_ ){

    mytree_pizeb = new TTree("pizSelb","Reco Simple Analysis");

    

    mytree_pizeb->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
    mytree_pizeb->Branch("runNumber",&runNumber,"runNumber/I");
    mytree_pizeb->Branch("evtNumber",&evtNumber,"evtNumber/I");

    //   bunchX ,, those are never used. ??
    //mytree_pizeb->Branch("bunchX",&bunchX,"bunchX/I");
    //mytree_pizeb->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
    mytree_pizeb->Branch("evtTime",&evtTime,"evtTime/I");
    
    mytree_pizeb->Branch("nL1Alca",&nL1Alca,"nL1Alca/I");
    mytree_pizeb->Branch("L1Alca",L1Alca,"L1Alca[nL1Alca]/I");
    
    mytree_pizeb->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
    //mytree_pizeb->Branch("l1algoName","std::vector<std::string>",&l1algoName);
    //algoName should be same for each run at least 
    
    mytree_pizeb->Branch("mpair",&mpair,"mpair/F");
    mytree_pizeb->Branch("ptpair",&ptpair,"ptpair/F");
    mytree_pizeb->Branch("etapair",&etapair,"etapair/F");
    ///mytree_pizeb->Branch("phipair",&phipair,"phipair/F");
    mytree_pizeb->Branch("ptmin",&ptmin,"ptmin/F");
    mytree_pizeb->Branch("isolation",&isolation,"isolation/F");
    
    mytree_pizeb->Branch("s4s9min",&s4s9min,"s4s9min/F");
    //mytree_pizeb->Branch("s9s25min",&s9s25min,"s9s25min/F");
    
  
  //   if( fullRECO_ ){
      
//       //mytree_pizeb->Branch("mpairv1",&mpairv1,"mpairv1/F"); ///using PV 
//       mytree_pizeb->Branch("vBeamSpot",vBeamSpot,"vBeamSpot[3]/F");
//       ///vertex all svaed
//       mytree_pizeb->Branch("nVertex",&nVertex,"nVertex/I");
//       mytree_pizeb->Branch("vertexx",vertexx,"vertexx[nVertex]/F");
//       mytree_pizeb->Branch("vertexy",vertexy,"vertexy[nVertex]/F");
//       mytree_pizeb->Branch("vertexz",vertexz,"vertexz[nVertex]/F");
//       mytree_pizeb->Branch("vertexndof",vertexndof,"vertexndof[nVertex]/F");
//       mytree_pizeb->Branch("vertexisFake",vertexisFake,"vertexisFake[nVertex]/I");
//       mytree_pizeb->Branch("vertexisValid",vertexisValid,"vertexisValid[nVertex]/I");
      

//       mytree_pizeb->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
//       mytree_pizeb->Branch("l1algoName","std::vector<std::string>",&l1algoName);
      
//       mytree_pizeb->Branch("hlt_bitFired", "std::vector<unsigned short>", &hlt_bitFired);  ///vector of those bit which fired the trigger 
//       mytree_pizeb->Branch("hlt_pathName", "std::vector<std::string>", &hlt_pathName);/// vector of string for all trigger path
      
//     }
    
    ///now add x,y,z of each Clus, mainly for history plot with vertex                                                                                                                 
    mytree_pizeb->Branch("xClus1",&xClus1,"xClus1/F");
    mytree_pizeb->Branch("yClus1",&yClus1,"yClus1/F");
    mytree_pizeb->Branch("zClus1",&zClus1,"zClus1/F");
    mytree_pizeb->Branch("xClus2",&xClus2,"xClus2/F");
    mytree_pizeb->Branch("yClus2",&yClus2,"yClus2/F");
    mytree_pizeb->Branch("zClus2",&zClus2,"zClus2/F");
    
    
    
    mytree_pizeb->Branch("nxtClus1",&nxtClus1,"nxtClus1/I");
    mytree_pizeb->Branch("eXtalClus1",eXtalClus1,"eXtalClus1[nxtClus1]/F");
    mytree_pizeb->Branch("laserCorrXtalClus1",laserCorrXtalClus1,"laserCorrXtalClus1[nxtClus1]/F");
    mytree_pizeb->Branch("ietaXtalClus1",ietaXtalClus1,"ietaXtalClus1[nxtClus1]/I");
    mytree_pizeb->Branch("iphiXtalClus1",iphiXtalClus1,"iphiXtalClus1[nxtClus1]/I");
    mytree_pizeb->Branch("tXtalClus1",tXtalClus1,"tXtalClus1[nxtClus1]/F"); //timing
  
    mytree_pizeb->Branch("nxtClus2",&nxtClus2,"nxtClus2/I");
    mytree_pizeb->Branch("laserCorrXtalClus2",laserCorrXtalClus2,"laserCorrXtalClus2[nxtClus2]/F");
    mytree_pizeb->Branch("eXtalClus2",eXtalClus2,"eXtalClus2[nxtClus2]/F");
    mytree_pizeb->Branch("ietaXtalClus2",ietaXtalClus2,"ietaXtalClus2[nxtClus2]/I");
    mytree_pizeb->Branch("iphiXtalClus2",iphiXtalClus2,"iphiXtalClus2[nxtClus2]/I");
    mytree_pizeb->Branch("tXtalClus2",tXtalClus2,"tXtalClus2[nxtClus2]/F"); //timing          
   
    ///gen info
    if( InputDataFormat_ <10){
      mytree_pizeb->Branch("genMatched",&genMatched,"genMatched/I");
      mytree_pizeb->Branch("geninfo",geninfo,"geninfo[9]/F");
      mytree_pizeb->Branch("convinfo",convinfo,"convinfo[2]/I");
      //      mytree_pizeb->Branch("ptHAT",&ptHAT,"ptHAT/F");
      //mytree_pizeb->Branch("procID",&procID,"procID/I");
    }
    
  }
  
  
  
  if( doSelForPi0Endcap_ ){
    
    mytree_pizee = new TTree("pizSele","Reco Simple Analysis");
    mytree_pizee->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
    mytree_pizee->Branch("runNumber",&runNumber,"runNumber/I");
    mytree_pizee->Branch("evtNumber",&evtNumber,"evtNumber/I");
    //   bunchX
    //mytree_pizee->Branch("bunchX",&bunchX,"bunchX/I");
    // mytree_pizee->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
    mytree_pizee->Branch("evtTime",&evtTime,"evtTime/I");
    

    
//     mytree_pizee->Branch("nL1Alca",&nL1Alca,"nL1Alca/I");
//     mytree_pizee->Branch("L1Alca",L1Alca,"L1Alca[nL1Alca]/I");
    
    
    mytree_pizee->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
    //mytree_pizee->Branch("l1algoName","std::vector<std::string>",&l1algoName);
    

    mytree_pizee->Branch("mpair",&mpair,"mpair/F");
    mytree_pizee->Branch("ptpair",&ptpair,"ptpair/F");
    mytree_pizee->Branch("etapair",&etapair,"etapair/F");
    ////mytree_pizee->Branch("phipair",&phipair,"phipair/F");
    mytree_pizee->Branch("ptmin",&ptmin,"ptmin/F");
    mytree_pizee->Branch("isolation",&isolation,"isolation/F");
    
    mytree_pizee->Branch("s4s9min",&s4s9min,"s4s9min/F");
    //mytree_pizee->Branch("s9s25min",&s9s25min,"s9s25min/F");
    
 //    if( fullRECO_ ){
      
//       mytree_pizee->Branch("vBeamSpot",vBeamSpot,"vBeamSpot[3]/F");
//       mytree_pizee->Branch("nVertex",&nVertex,"nVertex/I");
//       mytree_pizee->Branch("vertexx",vertexx,"vertexx[nVertex]/F");
//       mytree_pizee->Branch("vertexy",vertexy,"vertexy[nVertex]/F");
//       mytree_pizee->Branch("vertexz",vertexz,"vertexz[nVertex]/F");
//       mytree_pizee->Branch("vertexndof",vertexndof,"vertexndof[nVertex]/F");
//       mytree_pizee->Branch("vertexisFake",vertexisFake,"vertexisFake[nVertex]/I");
//       mytree_pizee->Branch("vertexisValid",vertexisValid,"vertexisValid[nVertex]/I");
      
      
//       mytree_pizee->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
//       mytree_pizee->Branch("l1algoName","std::vector<std::string>",&l1algoName);
//       mytree_pizee->Branch("hlt_bitFired", "std::vector<unsigned short>", &hlt_bitFired);  ///vector of those bit which fired the trigger 
//       mytree_pizee->Branch("hlt_pathName", "std::vector<std::string>", &hlt_pathName);/// vector of string for all trigger path
      
//     }
    
    mytree_pizee->Branch("nxtClus1",&nxtClus1,"nxtClus1/I");

    ///now add x,y,z of each Clus, mainly for history plot with vertex 
    mytree_pizee->Branch("xClus1",&xClus1,"xClus1/F");
    mytree_pizee->Branch("yClus1",&yClus1,"yClus1/F");
    mytree_pizee->Branch("zClus1",&zClus1,"zClus1/F");
    mytree_pizee->Branch("xClus2",&xClus2,"xClus2/F");
    mytree_pizee->Branch("yClus2",&yClus2,"yClus2/F");
    mytree_pizee->Branch("zClus2",&zClus2,"zClus2/F");
    


    mytree_pizee->Branch("eXtalClus1",eXtalClus1,"eXtalClus1[nxtClus1]/F");
    mytree_pizee->Branch("laserCorrXtalClus1",laserCorrXtalClus1,"laserCorrXtalClus1[nxtClus1]/F");
    mytree_pizee->Branch("ietaXtalClus1",ietaXtalClus1,"ietaXtalClus1[nxtClus1]/I");
    mytree_pizee->Branch("iphiXtalClus1",iphiXtalClus1,"iphiXtalClus1[nxtClus1]/I");
    mytree_pizee->Branch("izXtalClus1",&izXtalClus1,"izXtalClus1/I");
    mytree_pizee->Branch("tXtalClus1",tXtalClus1,"tXtalClus1[nxtClus1]/F"); //timing
    
    mytree_pizee->Branch("nxtClus2",&nxtClus2,"nxtClus2/I");
    mytree_pizee->Branch("laserCorrXtalClus2",laserCorrXtalClus2,"laserCorrXtalClus2[nxtClus2]/F");
    mytree_pizee->Branch("eXtalClus2",eXtalClus2,"eXtalClus2[nxtClus2]/F");
    mytree_pizee->Branch("izXtalClus2",&izXtalClus2,"izXtalClus2/I");
    mytree_pizee->Branch("ietaXtalClus2",ietaXtalClus2,"ietaXtalClus2[nxtClus2]/I");
    mytree_pizee->Branch("iphiXtalClus2",iphiXtalClus2,"iphiXtalClus2[nxtClus2]/I");
    mytree_pizee->Branch("tXtalClus2",tXtalClus2,"tXtalClus2[nxtClus2]/F"); //timing          
    

    //PreshowerS infoESX[2][8],  infoESY[2][8], 
    mytree_pizee->Branch("infoESX",infoESX,"infoESX[2][8]/F"); //timing   
    mytree_pizee->Branch("infoESY",infoESY,"infoESY[2][8]/F"); //timing   
    
    
    ///gen info
    if( InputDataFormat_ <10){
      mytree_pizee->Branch("genMatched",&genMatched,"genMatched/I");
      mytree_pizee->Branch("geninfo",geninfo,"geninfo[9]/F");
      mytree_pizee->Branch("convinfo",convinfo,"convinfo[2]/I");
      //mytree_pizee->Branch("ptHAT",&ptHAT,"ptHAT/F");
      // mytree_pizee->Branch("procID",&procID,"procID/I");
    }
    
  }
  
  
  
  
  if( doSelForEtaBarrel_ ){
    
    mytree_etaeb = new TTree("etaSelb","Reco Simple Analysis");
    mytree_etaeb->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
    mytree_etaeb->Branch("runNumber",&runNumber,"runNumber/I");
    mytree_etaeb->Branch("evtNumber",&evtNumber,"evtNumber/I");
    //   bunchX

    //mytree_etaeb->Branch("bunchX",&bunchX,"bunchX/I");
    // mytree_etaeb->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
    mytree_etaeb->Branch("evtTime",&evtTime,"evtTime/I");
    
    
//     mytree_etaeb->Branch("nL1Alca",&nL1Alca,"nL1Alca/I");
//     mytree_etaeb->Branch("L1Alca",L1Alca,"L1Alca[nL1Alca]/I");
    
    
    mytree_etaeb->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
    ///mytree_etaeb->Branch("l1algoName","std::vector<std::string>",&l1algoName);
    
    
    mytree_etaeb->Branch("mpair",&mpair,"mpair/F");
    mytree_etaeb->Branch("ptpair",&ptpair,"ptpair/F");
    mytree_etaeb->Branch("etapair",&etapair,"etapair/F");
    ///mytree_etaeb->Branch("phipair",&phipair,"phipair/F");
    mytree_etaeb->Branch("ptmin",&ptmin,"ptmin/F");
    mytree_etaeb->Branch("isolation",&isolation,"isolation/F");
    
    mytree_etaeb->Branch("s4s9min",&s4s9min,"s4s9min/F");
    mytree_etaeb->Branch("s9s25min",&s9s25min,"s9s25min/F");
    
    
//     if( fullRECO_ ){
      
//       mytree_etaeb->Branch("vBeamSpot",vBeamSpot,"vBeamSpot[3]/F");
//       mytree_etaeb->Branch("nVertex",&nVertex,"nVertex/I");
//       mytree_etaeb->Branch("vertexx",vertexx,"vertexx[nVertex]/F");
//       mytree_etaeb->Branch("vertexy",vertexy,"vertexy[nVertex]/F");
//       mytree_etaeb->Branch("vertexz",vertexz,"vertexz[nVertex]/F");
//       mytree_etaeb->Branch("vertexndof",vertexndof,"vertexndof[nVertex]/F");
//       mytree_etaeb->Branch("vertexisFake",vertexisFake,"vertexisFake[nVertex]/I");
//       mytree_etaeb->Branch("vertexisValid",vertexisValid,"vertexisValid[nVertex]/I");
      

//       mytree_etaeb->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
//       mytree_etaeb->Branch("l1algoName","std::vector<std::string>",&l1algoName);
//       mytree_etaeb->Branch("hlt_bitFired", "std::vector<unsigned short>", &hlt_bitFired);  ///vector of those bit which fired the trigger 
//       mytree_etaeb->Branch("hlt_pathName", "std::vector<std::string>", &hlt_pathName);/// vector of string for all trigger path
      
      
//     }

    
    mytree_etaeb->Branch("xClus1",&xClus1,"xClus1/F");
    mytree_etaeb->Branch("yClus1",&yClus1,"yClus1/F");
    mytree_etaeb->Branch("zClus1",&zClus1,"zClus1/F");
    mytree_etaeb->Branch("xClus2",&xClus2,"xClus2/F");
    mytree_etaeb->Branch("yClus2",&yClus2,"yClus2/F");
    mytree_etaeb->Branch("zClus2",&zClus2,"zClus2/F");


    
    mytree_etaeb->Branch("nxtClus1",&nxtClus1,"nxtClus1/I");
    mytree_etaeb->Branch("laserCorrXtalClus1",laserCorrXtalClus1,"laserCorrXtalClus1[nxtClus1]/F");
    mytree_etaeb->Branch("eXtalClus1",eXtalClus1,"eXtalClus1[nxtClus1]/F");
    mytree_etaeb->Branch("ietaXtalClus1",ietaXtalClus1,"ietaXtalClus1[nxtClus1]/I");
    mytree_etaeb->Branch("iphiXtalClus1",iphiXtalClus1,"iphiXtalClus1[nxtClus1]/I");
    mytree_etaeb->Branch("tXtalClus1",tXtalClus1,"tXtalClus1[nxtClus1]/F"); //timing
    
    mytree_etaeb->Branch("nxtClus2",&nxtClus2,"nxtClus2/I");
    mytree_etaeb->Branch("laserCorrXtalClus2",laserCorrXtalClus2,"laserCorrXtalClus2[nxtClus2]/F");
    mytree_etaeb->Branch("eXtalClus2",eXtalClus2,"eXtalClus2[nxtClus2]/F");
    mytree_etaeb->Branch("ietaXtalClus2",ietaXtalClus2,"ietaXtalClus2[nxtClus2]/I");
    mytree_etaeb->Branch("iphiXtalClus2",iphiXtalClus2,"iphiXtalClus2[nxtClus2]/I");
    mytree_etaeb->Branch("tXtalClus2",tXtalClus2,"tXtalClus2[nxtClus2]/F"); //timing      
    
    ///5x5
    mytree_etaeb->Branch("nxt5x5Clus1",&nxt5x5Clus1,"nxt5x5Clus1/I");
    mytree_etaeb->Branch("eXtal5x5Clus1",eXtal5x5Clus1,"eXtal5x5Clus1[nxt5x5Clus1]/F");
    mytree_etaeb->Branch("ietaXtal5x5Clus1",ietaXtal5x5Clus1,"ietaXtal5x5Clus1[nxt5x5Clus1]/I");
    mytree_etaeb->Branch("iphiXtal5x5Clus1",iphiXtal5x5Clus1,"iphiXtal5x5Clus1[nxt5x5Clus1]/I");
    mytree_etaeb->Branch("tXtal5x5Clus1",tXtal5x5Clus1,"tXtal5x5Clus1[nxt5x5Clus1]/F"); //timing
  
    mytree_etaeb->Branch("nxt5x5Clus2",&nxt5x5Clus2,"nxt5x5Clus2/I");
    mytree_etaeb->Branch("eXtal5x5Clus2",eXtal5x5Clus2,"eXtal5x5Clus2[nxt5x5Clus2]/F");
    mytree_etaeb->Branch("ietaXtal5x5Clus2",ietaXtal5x5Clus2,"ietaXtal5x5Clus2[nxt5x5Clus2]/I");
    mytree_etaeb->Branch("iphiXtal5x5Clus2",iphiXtal5x5Clus2,"iphiXtal5x5Clus2[nxt5x5Clus2]/I");
    mytree_etaeb->Branch("tXtal5x5Clus2",tXtal5x5Clus2,"tXtal5x5Clus2[nxt5x5Clus2]/F"); //timing   
        
    
    ///gen info
    if( InputDataFormat_ <10){
      mytree_etaeb->Branch("genMatched",&genMatched,"genMatched/I");
      mytree_etaeb->Branch("geninfo",geninfo,"geninfo[9]/F");
      mytree_etaeb->Branch("convinfo",convinfo,"convinfo[2]/I");
      
      //mytree_etaeb->Branch("ptHAT",&ptHAT,"ptHAT/F");
      //mytree_etaeb->Branch("procID",&procID,"procID/I");
    }
    
  }
  
  
  if( doSelForEtaEndcap_ ){
    
    mytree_etaee = new TTree("etaSele","Reco Simple Analysis");
    mytree_etaee->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
    mytree_etaee->Branch("runNumber",&runNumber,"runNumber/I");
    mytree_etaee->Branch("evtNumber",&evtNumber,"evtNumber/I");
    //   bunchX
    //mytree_etaee->Branch("bunchX",&bunchX,"bunchX/I");
    // mytree_etaee->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
    mytree_etaee->Branch("evtTime",&evtTime,"evtTime/I");
    
    
//     mytree_etaee->Branch("nL1Alca",&nL1Alca,"nL1Alca/I");
//     mytree_etaee->Branch("L1Alca",L1Alca,"L1Alca[nL1Alca]/I");
    
    
    mytree_etaee->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
    //mytree_etaee->Branch("l1algoName","std::vector<std::string>",&l1algoName);
    
    
    mytree_etaee->Branch("mpair",&mpair,"mpair/F");
    mytree_etaee->Branch("ptpair",&ptpair,"ptpair/F");
    mytree_etaee->Branch("etapair",&etapair,"etapair/F");
    ///mytree_etaee->Branch("phipair",&phipair,"phipair/F");
    mytree_etaee->Branch("ptmin",&ptmin,"ptmin/F");
    mytree_etaee->Branch("isolation",&isolation,"isolation/F");
    
    mytree_etaee->Branch("s4s9min",&s4s9min,"s4s9min/F");
    mytree_etaee->Branch("s9s25min",&s9s25min,"s9s25min/F");
    ///mytree_etaee->Branch("s9s25min",&s9s25min,"s9s25min/F");
    
 //    if( fullRECO_ ){
      
//       mytree_etaee->Branch("vBeamSpot",vBeamSpot,"vBeamSpot[3]/F");
//       mytree_etaee->Branch("nVertex",&nVertex,"nVertex/I");
//       mytree_etaee->Branch("vertexx",vertexx,"vertexx[nVertex]/F");
//       mytree_etaee->Branch("vertexy",vertexy,"vertexy[nVertex]/F");
//       mytree_etaee->Branch("vertexz",vertexz,"vertexz[nVertex]/F");
//       mytree_etaee->Branch("vertexndof",vertexndof,"vertexndof[nVertex]/F");
//       mytree_etaee->Branch("vertexisFake",vertexisFake,"vertexisFake[nVertex]/I");
//       mytree_etaee->Branch("vertexisValid",vertexisValid,"vertexisValid[nVertex]/I");
      

//       mytree_etaee->Branch("l1bitFired","std::vector<unsigned short>",&l1bitFired);
//       mytree_etaee->Branch("l1algoName","std::vector<std::string>",&l1algoName);

//       mytree_etaee->Branch("hlt_bitFired", "std::vector<unsigned short>", &hlt_bitFired);  ///vector of those bit which fired the trigger 
//       mytree_etaee->Branch("hlt_pathName", "std::vector<std::string>", &hlt_pathName);/// vector of string for all trigger path
      
//     }
    
    
    
    mytree_etaee->Branch("xClus1",&xClus1,"xClus1/F");
    mytree_etaee->Branch("yClus1",&yClus1,"yClus1/F");
    mytree_etaee->Branch("zClus1",&zClus1,"zClus1/F");
    mytree_etaee->Branch("xClus2",&xClus2,"xClus2/F");
    mytree_etaee->Branch("yClus2",&yClus2,"yClus2/F");
    mytree_etaee->Branch("zClus2",&zClus2,"zClus2/F");
    
    
    mytree_etaee->Branch("nxtClus1",&nxtClus1,"nxtClus1/I");
    mytree_etaee->Branch("laserCorrXtalClus1",laserCorrXtalClus1,"laserCorrXtalClus1[nxtClus1]/F");
    mytree_etaee->Branch("eXtalClus1",eXtalClus1,"eXtalClus1[nxtClus1]/F");
    mytree_etaee->Branch("ietaXtalClus1",ietaXtalClus1,"ietaXtalClus1[nxtClus1]/I");
    mytree_etaee->Branch("iphiXtalClus1",iphiXtalClus1,"iphiXtalClus1[nxtClus1]/I");
    mytree_etaee->Branch("tXtalClus1",tXtalClus1,"tXtalClus1[nxtClus1]/F"); //timing
    mytree_etaee->Branch("izXtalClus1",&izXtalClus1,"izXtalClus1/I");
    
    mytree_etaee->Branch("nxtClus2",&nxtClus2,"nxtClus2/I");
    mytree_etaee->Branch("laserCorrXtalClus2",laserCorrXtalClus2,"laserCorrXtalClus2[nxtClus2]/F");
    mytree_etaee->Branch("eXtalClus2",eXtalClus2,"eXtalClus2[nxtClus2]/F");
    mytree_etaee->Branch("ietaXtalClus2",ietaXtalClus2,"ietaXtalClus2[nxtClus2]/I");
    mytree_etaee->Branch("iphiXtalClus2",iphiXtalClus2,"iphiXtalClus2[nxtClus2]/I");
    mytree_etaee->Branch("tXtalClus2",tXtalClus2,"tXtalClus2[nxtClus2]/F"); //timing      
    mytree_etaee->Branch("izXtalClus2",&izXtalClus2,"izXtalClus2/I");
    

    ///5x5
    mytree_etaee->Branch("nxt5x5Clus1",&nxt5x5Clus1,"nxt5x5Clus1/I");
    mytree_etaee->Branch("eXtal5x5Clus1",eXtal5x5Clus1,"eXtal5x5Clus1[nxt5x5Clus1]/F");
    mytree_etaee->Branch("ietaXtal5x5Clus1",ietaXtal5x5Clus1,"ietaXtal5x5Clus1[nxt5x5Clus1]/I");
    mytree_etaee->Branch("iphiXtal5x5Clus1",iphiXtal5x5Clus1,"iphiXtal5x5Clus1[nxt5x5Clus1]/I");
    ////mytree_etaee->Branch("iphiXtal5x5Clus1",iphiXtal5x5Clus1,"izXtal5x5Clus1[nxt5x5Clus1]/I"); for 5x5 don't save iz, the same as 
    mytree_etaee->Branch("tXtal5x5Clus1",tXtal5x5Clus1,"tXtal5x5Clus1[nxt5x5Clus1]/F"); //timing
  
    mytree_etaee->Branch("nxt5x5Clus2",&nxt5x5Clus2,"nxt5x5Clus2/I");
    mytree_etaee->Branch("eXtal5x5Clus2",eXtal5x5Clus2,"eXtal5x5Clus2[nxt5x5Clus2]/F");
    mytree_etaee->Branch("ietaXtal5x5Clus2",ietaXtal5x5Clus2,"ietaXtal5x5Clus2[nxt5x5Clus2]/I");
    mytree_etaee->Branch("iphiXtal5x5Clus2",iphiXtal5x5Clus2,"iphiXtal5x5Clus2[nxt5x5Clus2]/I");
    mytree_etaee->Branch("tXtal5x5Clus2",tXtal5x5Clus2,"tXtal5x5Clus2[nxt5x5Clus2]/F"); //timing   
    
    //PreshowerS infoESX[2][8],  infoESY[2][8], 
    mytree_etaee->Branch("infoESX",infoESX,"infoESX[2][8]/F"); //timing   
    mytree_etaee->Branch("infoESY",infoESY,"infoESY[2][8]/F"); //timing   
    
    ///gen info
    if( InputDataFormat_ <10){
      mytree_etaee->Branch("genMatched",&genMatched,"genMatched/I");
      mytree_etaee->Branch("geninfo",geninfo,"geninfo[9]/F");
      mytree_etaee->Branch("convinfo",convinfo,"convinfo[2]/I");
      //mytree_etaee->Branch("ptHAT",&ptHAT,"ptHAT/F");
      //mytree_etaee->Branch("procID",&procID,"procID/I");
    }
    
  }
  
  
  
  
}



// ------------ method called once each job just after ending the event loop  ------------
void 
RecoAnalyzer::endJob() {
  

  rootFile_->cd();

  
  if( fullRECO_   ||  getBeamSpotOnly_){
    mytree_evtInfo->Write();
  }
  

  // hh_mul_ieta->Write();
  // hh_mul_ietaSeed->Write();
  // hh_mul_ietaSeedClus->Write();
  // hh_mul_ietaSeedSel->Write();
  // hh_mul_iphi->Write();
  // hh_mul_iphiSeed->Write();
  // hh_mul_iphiSeedClus->Write();
  // hh_mul_iphiSeedSel->Write();


  if(! getBeamSpotOnly_) {

    mytree_bs->Write(); ///now save also beamSpot info per lumiBlock into a new tree
    
  
  if(doSelForPi0Barrel_){
    mytree_pizeb->Write(); 
  }
  if(doSelForPi0Endcap_){
    mytree_pizee->Write();
  } 
  
  if(doSelForEtaBarrel_){
    mytree_etaeb->Write(); 
  }
  if(doSelForEtaEndcap_){
    mytree_etaee->Write();
  } 
  }
  
  
  rootFile_->Write() ;
  rootFile_->Close() ;
  
  
 
  
  
}

bool RecoAnalyzer::checkStatusOfEcalRecHit(const EcalChannelStatus &channelStatus,const EcalRecHit &rh){
  
  if(useRecoFlag_ ){ ///from recoFlag()
    int flag = rh.recoFlag();
    if( flagLevelRecHitsToUse_ ==0){ ///good 
      if( flag != 0) return false; 
    }
    else if( flagLevelRecHitsToUse_ ==1){ ///good || PoorCalib 
      if( flag !=0 && flag != 4 ) return false; 
    }
    else if( flagLevelRecHitsToUse_ ==2){ ///good || PoorCalib || LeadingEdgeRecovered || kNeighboursRecovered,
      if( flag !=0 && flag != 4 && flag != 6 && flag != 7) return false; 
    }
  }
  if ( useDBStatus_){ //// from DB
    int status =  int(channelStatus[rh.id().rawId()].getStatusCode()); 
    if ( status > statusLevelRecHitsToUse_ ) return false; 
  }
  
  return true; 
}



// void RecoAnalyzer::makeClusterES(double x, double y, double z,const CaloSubdetectorGeometry*& geometry_es,
// 					CaloSubdetectorTopology*& topology_es
// 					){
  
  
//   ///get assosicated ES clusters of this endcap cluster
//   const GlobalPoint point(x,y,z);
//   DetId tmp1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_es))->getClosestCellInPlane(point, 1);
//   DetId tmp2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_es))->getClosestCellInPlane(point, 2);
//   ESDetId strip1 = (tmp1 == DetId(0)) ? ESDetId(0) : ESDetId(tmp1);
//   ESDetId strip2 = (tmp2 == DetId(0)) ? ESDetId(0) : ESDetId(tmp2);     
  
//   // Get ES clusters (found by the PreshSeeded algorithm) associated with a given EE cluster.           
//   for (int i2=0; i2<preshNclust_; i2++) {
//     reco::PreshowerCluster cl1 = presh_algo->makeOneCluster(strip1,&used_strips,&esrechits_map,geometry_es,topology_es);   
//     reco::PreshowerCluster cl2 = presh_algo->makeOneCluster(strip2,&used_strips,&esrechits_map,geometry_es,topology_es); 
//   } // end of cycle over ES clusters
    
  
// }


// //////FED list this is obsolete 
// std::vector<int> RecoAnalyzer::ListOfFEDS(double etaLow, double etaHigh, double phiLow, 
// 					 double phiHigh, double etamargin, double phimargin)
// {
  
// 	std::vector<int> FEDs;

// 	if (phimargin > Geom::pi()) phimargin =  Geom::pi() ;

	
// 	if (debug_>=2) std::cout << " etaLow etaHigh phiLow phiHigh " << etaLow << " " << 
// 			etaHigh << " " << phiLow << " " << phiHigh << std::endl;

//         etaLow -= etamargin;
//         etaHigh += etamargin;
//         double phiMinus = phiLow - phimargin;
//         double phiPlus = phiHigh + phimargin;

//         bool all = false;
//         double dd = fabs(phiPlus-phiMinus);
// 	if (debug_>=2) std::cout << " dd = " << dd << std::endl;
//         if (dd > 2.*Geom::pi() ) all = true;

//         while (phiPlus > Geom::pi()) { phiPlus -= 2.*Geom::pi() ; }
//         while (phiMinus < 0) { phiMinus += 2.*Geom::pi() ; }
//         if ( phiMinus > Geom::pi()) phiMinus -= 2.*Geom::pi() ;

//         double dphi = phiPlus - phiMinus;
//         if (dphi < 0) dphi += 2.*Geom::pi() ;
// 	if (debug_>=2) std::cout << "dphi = " << dphi << std::endl;
//         if (dphi > Geom::pi()) {
//                 int fed_low1 = TheMapping -> GetFED(etaLow,phiMinus*180./Geom::pi());
//                 int fed_low2 = TheMapping -> GetFED(etaLow,phiPlus*180./Geom::pi());
// 		if (debug_>=2) std::cout << "fed_low1 fed_low2 " << fed_low1 << " " << fed_low2 << std::endl;
//                 if (fed_low1 == fed_low2) all = true;
//                 int fed_hi1 = TheMapping -> GetFED(etaHigh,phiMinus*180./Geom::pi());
//                 int fed_hi2 = TheMapping -> GetFED(etaHigh,phiPlus*180./Geom::pi());
// 		if (debug_>=2) std::cout << "fed_hi1 fed_hi2 " << fed_hi1 << " " << fed_hi2 << std::endl;
//                 if (fed_hi1 == fed_hi2) all = true;
//         }

// 	if (all) {
// 		if (debug_>=2) std::cout << " unpack everything in phi ! " << std::endl;
// 		phiMinus = -20 * Geom::pi() / 180.;  // -20 deg
// 		phiPlus = -40 * Geom::pi() / 180.;  // -20 deg
// 	}

//         if (debug_>=2) std::cout << " with margins : " << etaLow << " " << etaHigh << " " << 
// 			phiMinus << " " << phiPlus << std::endl;


//         const EcalEtaPhiRegion ecalregion(etaLow,etaHigh,phiMinus,phiPlus);

//         FEDs = TheMapping -> GetListofFEDs(ecalregion);
	
// /*
// 	if (debug_) {
//            int nn = (int)FEDs.size();
//            for (int ii=0; ii < nn; ii++) {
// 	   std::cout << "unpack fed " << FEDs[ii] << std::endl;
//            }
//    	   }
// */
	
//         return FEDs;
	
// }


////already existing , int EcalElectronicsMapping::DCCid(const EBDetId& id)
///obsolete
// int RecoAnalyzer::convertSmToFedNumbBarrel(int ieta, int smId){
    
//   if( ieta<=-1) return smId - 9; 
//   else return smId + 27; 
  
  
// }


void RecoAnalyzer::convxtalid(int &nphi,int &neta)
{
  // Barrel only
  // Output nphi 0...359; neta 0...84; nside=+1 (for eta>0), or 0 (for eta<0).
  // neta will be [-85,-1] , or [0,84], the minus sign indicates the z<0 side.
  
  if(neta > 0) neta -= 1;
  if(nphi > 359) nphi=nphi-360;
  
  // final check
  if(nphi >359 || nphi <0 || neta< -85 || neta > 84)
    {
      std::cout <<" unexpected fatal error in RecoAnalyzer::convxtalid "<<  nphi <<  " " << neta <<  " " <<std::endl;
      //exit(1);
    }
} //end of convxtalid




int RecoAnalyzer::diff_neta_s(int neta1, int neta2){
  int mdiff;
  mdiff=(neta1-neta2);
  return mdiff;
}

// Calculate the distance in xtals taking into account the periodicity of the Barrel
int RecoAnalyzer::diff_nphi_s(int nphi1,int nphi2) {
   int mdiff;
   if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) {
     mdiff=nphi1-nphi2;
   }
   else {
   mdiff=360-abs(nphi1-nphi2);
   if(nphi1>nphi2) mdiff=-mdiff;
   }
   return mdiff;
}




double RecoAnalyzer::DeltaPhi(double v1, double v2)
{ // Computes the correctly normalized phi difference
  // v1, v2 = phi of object 1 and 2
  double diff = fabs(v2 - v1);
  double corr = 2*acos(-1.) - diff; //// 2*pi - diff;
  if (diff < acos(-1.)){ return diff;} else { return corr;}
  
}


double RecoAnalyzer::GetDeltaR(double eta1, double eta2, double phi1, double phi2){ 
  // Computes the DeltaR of two objects from their eta and phi values
  
  return sqrt( (eta1-eta2)*(eta1-eta2) 
	       + DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
  
}

double RecoAnalyzer::getcosd(double eta1, double phi1, double eta2, double phi2) {
  double theta1 = 2*atan(exp(-eta1));
  double theta2 = 2*atan(exp(-eta2));
  double cosd;
  double dphi = DeltaPhi(phi1,phi2);
  cosd = cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(dphi);  //opening angle 
  return cosd;
}

/// from EcalElectronicsSim:
/// input is 
void RecoAnalyzer::amplify(CaloSamples & clf) const
{
  clf *= theParameterMap->simParameters(clf.id()).photoelectronsToAnalog();
  if (applyConstantTerm_) {
    clf *= (1.+constantTerm());
  }
}

double RecoAnalyzer::constantTerm() const
{
  edm::Service<edm::RandomNumberGenerator> rng;
  if ( ! rng.isAvailable()) {
    throw cms::Exception("Configuration")
      << "EcalElectroncSim requires the RandomNumberGeneratorService\n"
      "which is not present in the configuration file.  You must add the service\n"
      "in the configuration file or remove the modules that require it.";
  }
  
  double thisCT = rmsConstantTerm_;
  CLHEP::RandGaussQ gaussQDistribution(rng->getEngine(), 0.0, thisCT);
  return gaussQDistribution.fire();
}









int RecoAnalyzer::indexofParticle(float px, float pz, int status){
  
  float err = 0.00001; 
  
  // int ind = -10; 
  
  vector<int> nn; 

  for( int j=0; j<nMCpart; j++){
    if(fabs(px-pxMCpart[j])<err && fabs(pz-pzMCpart[j])<err
       && status == statusMCpart[j]){
      nn.push_back(j);
    }
  }
  //   if(int(nn.size()) !=1){ ///diquarks or gluon has two copy sometimesss!!! 
  //     cout<<"wrong! indexofParticle: n: "<<int(nn.size())<<endl;
  //     for(int j=0; j<int(nn.size()); j++){
  //       cout<<"  dup j: "<<j<<" ind "<<nn[j]<<endl;
  //     }
  //   }
  return nn[0]; 
  
}


int RecoAnalyzer::getMotherIndex(int j){
  if(j<0) {
    cout<<"no mother. input -1!!!!!!!!"<<endl;
    //  cout<<"entry: "<<entry<<" input j: "<<j<<endl;
    // exit(1);
    return -1; 
  }
  
  int indmom = barcodemomMCpart[j];
  if(indmom>=0){
    int st = statusMCpart[indmom]; 
    int indf = indmom;
    if(st==3 && pidMCpart[indmom]==pidMCpart[j]){ //if ==3  && it's itself.  trace back more.
      indf = barcodemomMCpart[indmom];
      if(indf>=0){
	st = statusMCpart[indf];
      }
      else cout<<"warning! ptcl has no mother!"<<endl;
    }
    return indf; 
    
  }
  else{
    return -1; 
  }
  
}




void RecoAnalyzer::findgenpi0eta(){
  
  indpi0Gen.clear();
  indetaGen.clear();
  indpht1pi0Gen.clear();
  indpht2pi0Gen.clear();
  indpht1etaGen.clear();
  indpht2etaGen.clear();
  

  int pidPi0 = 111; 
  int pidPht = 22; 
  int pidEta = 221; 
  
  ////float eta_max = 1000000; 
  //  float et_min = 2.0; 
  
  ////float epht_min = 0.1; 
  
  ///epht_min = 0.; 
    

  ////txtout<<"entry: "<<entry<<endl;
  
  vector<int> indphtpi0Gen; 
  vector<int> indphtetaGen; 
  
  int npi0all = 0; 
  int netaall = 0; 
  
  
  for( int j=0; j<nMCpart; j++){
    

    //if( pidMCpart[j] == pidPi0 && eMCpart[j]>2*epht_min && fabs(etaMCpart[j])<eta_max ) npi0all++; 
    if( pidMCpart[j] == pidPi0 ) npi0all++; 
    if( pidMCpart[j] == pidEta ) netaall++; 
    
    ///  if( pidMCpart[j] == pidPi0 && eMCpart[j]>2*epht_min ) npi0all++; 
    
    //tesme
    ///float ecalEta(float EtaParticle ,float Zvertex, float RhoVertex){
    ////photon from pi0 has status =1
     //    if( statusMCpart[j] ==1 && pidMCpart[j] == pidPht && eMCpart[j]>epht_min){
    if( statusMCpart[j] ==1 && pidMCpart[j] == pidPht  ){
      
      //float eta = ecalEta(etaMCpart[j],vtzMCpart[j]/10.0,sqrt(vtxMCpart[j]*vtxMCpart[j]+vtyMCpart[j]*vtyMCpart[j])/10.0);
      //float eta = etaMCpart[j];
      int indmom = getMotherIndex(j);
      ////vtx of photon is not exactly the same as vtx of it's mother pi0 and two photons not exactly the same either...
      ///	cout<<"vtz/R: "<<vtzMCpart[j]<<"/"<<sqrt(vtxMCpart[j]*vtxMCpart[j]+vtyMCpart[j]*vtyMCpart[j])<<" "<<"eta/etacorr: "<<etaMCpart[j]<<"/"<<eta<<" vtxmomz: "<<vtzMCpart[indmom]<<endl;
      
      if(pidMCpart[indmom] == pidPi0){
	indphtpi0Gen.push_back(j);
	////txtout<<"pi0->gamma: "<<indmom<<" "<<j<<endl;
      }
      if(pidMCpart[indmom] == pidEta){
	indphtetaGen.push_back(j);
      }
    } //end of found one photon
    
     //     if(pidMCpart[j] == pidPi0 && eMCpart[j]>2*epht_min ){
//       int nda = nDauMCpart[j]; 
//       txtout<<" pi0: "<<nda<<" e: "<<eMCpart[j]<<" eta: "<<etaMCpart[j]<<endl;
//       //   for( int n =0; n<nda; n++){ ///16X doesn't work
//       // 	int indd = barcodeDauMCpart[j][n];
//       // 	txtout<<" n : "<<n<<" "<<pidMCpart[indd]<<" e: "<<eMCpart[indd]<<endl; 
//       //       }
//     }
    
    
    /////131 barcodeDauMCpart[j][n] for pi0 doesn't work
     
  }
  ///txtout<<"nph: "<<int(indphtpi0Gen.size())<<endl;
  
  vector< std::pair<double,int> > indpiz_eta; 
  vector< std::pair<double,int> > indeta_eta; 
  
  int npi0_sel = 0; 
  int neta_sel = 0; 
  
  
  for( int j=0; j<int(indphtpi0Gen.size());j++){
    int ind = getMotherIndex(indphtpi0Gen[j]);
    int nsame = 0; 
    vector<int> isame; 
    for(int k=0; k<int(indphtpi0Gen.size());k++){
      if( k!=j && getMotherIndex(indphtpi0Gen[k]) == ind) {
	isame.push_back(k);
	nsame++; 
      }
    }
    if(nsame==1){
      // txtout<<"j isame0: "<<j<<" "<<isame[0]<<endl;
      
      int ind1 =indphtpi0Gen[j];
      int ind2 = indphtpi0Gen[isame[0]];

      indpht1pi0Gen.push_back(ind1);
      indpht2pi0Gen.push_back(ind2);
     
      indpi0Gen.push_back(ind);
      indphtpi0Gen.erase(indphtpi0Gen.begin()+j);
      indphtpi0Gen.erase(indphtpi0Gen.begin()+isame[0]-1);

      indpiz_eta.push_back( make_pair(fabs(etaMCpart[ind]), npi0_sel));
      npi0_sel ++; 
      
      
      j = j-1; ////start from previous j again becasue size(indphtpi0Gen) decreased by 2. 
    }
    
  }
  /////eta
  for( int j=0; j<int(indphtetaGen.size());j++){
    int ind = getMotherIndex(indphtetaGen[j]);
    int nsame = 0; 
    vector<int> isame; 
    for(int k=0; k<int(indphtetaGen.size());k++){
      if( k!=j && getMotherIndex(indphtetaGen[k]) == ind) {
	isame.push_back(k);
	nsame++; 
      }
    }
    if(nsame==1){
      // txtout<<"j isame0: "<<j<<" "<<isame[0]<<endl;
      
      int ind1 =indphtetaGen[j];
      int ind2 = indphtetaGen[isame[0]];
      
      indpht1etaGen.push_back(ind1);
      indpht2etaGen.push_back(ind2);
     
      indetaGen.push_back(ind);
      indphtetaGen.erase(indphtetaGen.begin()+j);
      indphtetaGen.erase(indphtetaGen.begin()+isame[0]-1);
      
      indeta_eta.push_back( make_pair(fabs(etaMCpart[ind]), neta_sel) );
      neta_sel ++; 
      
      j = j-1; ////start from previous j again becasue size(indphtpi0Gen) decreased by 2. 
    }
    
  }
  
  
  //fillthe tree
  
  
  //  vector<int> indpi0GenEtaSorted; 
  // vector<int> indetaGenEtaSorted; 
  
  
  sort(indpiz_eta.begin(),indpiz_eta.end(),sort_pred);
  sort(indeta_eta.begin(),indeta_eta.end(),sort_pred);
  
  //// small to large
  // sortGenEnergy(indpi0Gen,indpi0GenESorted);
  //sortGenEnergy(indetaGen,indetaGenESorted);
  
  ///eta smaller to large
  ///  sortGenEta(indpi0Gen,indpi0GenEtaSorted);
  /// sortGenEta(indetaGen,indetaGenEtaSorted);
  
  
  nGenpi0 = 0; 
  
  
  if(debug_>=1) cout<<"ngenpi0: "<< npi0all <<" "<< int(indpi0Gen.size())<<endl;
  

  //  for( int j=0; j<int(indpi0GenEtaSorted.size()) && nGenpi0<MAXGenPIZ; j++){
  for( int j=0; j<int(indpiz_eta.size()) && nGenpi0<MAXGenPIZ; j++){
    int k = indpiz_eta[j].second; 
    
    // int k = indpi0GenEtaSorted[j];
    
    
    int indmc = indpi0Gen[k];
    int indph1 = indpht1pi0Gen[k];
    int indph2 = indpht2pi0Gen[k];


    // if(debug_ >=2) cout<<"  indpiz_eta: "<< indmc<<" "<< etaMCpart[indmc]<<" "<< indph1<<" "<<indph2<<endl;
    
    
      ///mGenpi0[nGenpi0] = getMasspairMCpart(indph1,indph2);
      
      eGenpi0[nGenpi0] = eMCpart[indmc];
      //etaGenpi0[nGenpi0] = ecalEta(etaMCpart[indmc],vtzMCpart[indmc]/10.0,sqrt(vtxMCpart[indmc]*vtxMCpart[indmc]+vtyMCpart[indmc]*vtyMCpart[indmc])/10.0);
      phiGenpi0[nGenpi0] = phiMCpart[indmc];
      etaGenpi0[nGenpi0] = etaMCpart[indmc];
      //etGenpi0[nGenpi0] = eGenpi0[nGenpi0]*sin(2*atan(exp(-etaGenpi0[nGenpi0])));
    
      ////txtout<<"entry: "<<entry<<" "<<j<<" pi0 et/eta/phi: "<<etGenpi0[nGenpi0]<<" "<< etaGenpi0[nGenpi0]<<" " <<phiGenpi0[nGenpi0]<<endl;
      
      
      //let's do those later.keep the orignal and vtx info
       
      //  if(debug_ >=2) cout<<"  indpiz_eta1: "<< indmc<<" "<< etaMCpart[indmc]<<" "<< indph1<<" "<<indph2<<endl;


      vtxGenpi0[nGenpi0][0] = vtxXMCpart[indmc];
      vtxGenpi0[nGenpi0][1] = vtxYMCpart[indmc];
      vtxGenpi0[nGenpi0][2] = vtxZMCpart[indmc];
      
      barcodeMomGenpi0[nGenpi0]= getMotherIndex(indmc);
      pidMomGenpi0[nGenpi0]= pidMCpart[barcodeMomGenpi0[nGenpi0]];
      

      ///a littel difference
      vtxPhtGenpi0[nGenpi0][0] = vtxXMCpart[indph1];
      vtxPhtGenpi0[nGenpi0][1] = vtxYMCpart[indph1];
      vtxPhtGenpi0[nGenpi0][2] = vtxZMCpart[indph1];

      //  if(debug_ >=2) cout<<"  indpiz_eta2: "<< indmc<<" "<< etaMCpart[indmc]<<" "<< indph1<<" "<<indph2<<endl;

      ePhtGenpi0[nGenpi0][0] = eMCpart[indph1];
      ePhtGenpi0[nGenpi0][1] = eMCpart[indph2];
      etaPhtGenpi0[nGenpi0][0] = etaMCpart[indph1];
      etaPhtGenpi0[nGenpi0][1] = etaMCpart[indph2]; 
      phiPhtGenpi0[nGenpi0][0] = phiMCpart[indph1];
      phiPhtGenpi0[nGenpi0][1] = phiMCpart[indph2];
      
      //      if(debug_ >=2) cout<<"  indpiz_eta3: "<< indmc<<" "<< etaMCpart[indmc]<<" "<< indph1<<" "<<indph2<<endl;
      
      // if(debug_ >=2) cout<<"  indpiz_eta4: "<< indmc<<" "<< etaMCpart[indmc]<<endl;
      
      isConvPhtGenpi0[nGenpi0][0] = convPhtMCpart[indph1][0];
      isConvPhtGenpi0[nGenpi0][1] = convPhtMCpart[indph2][0];

      int e1 = convPhtMCpart[indph1][1]; 
      int e2 = convPhtMCpart[indph1][2]; 
      
      if( e1 >0){
	convPht1Genpi0[nGenpi0][0] = eMCpart[e1];
	convPht1Genpi0[nGenpi0][1] = etaMCpart[e1];
	convPht1Genpi0[nGenpi0][2] = phiMCpart[e1];
	
	
      }
      if( e2 >0){
	convPht1Genpi0[nGenpi0][3] = eMCpart[e2];
	convPht1Genpi0[nGenpi0][4] = etaMCpart[e2];
	convPht1Genpi0[nGenpi0][5] = phiMCpart[e2];
      }

      for(int jj=0; jj<6; jj++){
	convVtxPhtGenpi0[nGenpi0][jj] = 0; 
      }
      if( e1 >0){
	convVtxPhtGenpi0[nGenpi0][0] = vtxXMCpart[e1];
	convVtxPhtGenpi0[nGenpi0][1] = vtxYMCpart[e1];
	convVtxPhtGenpi0[nGenpi0][2] = vtxZMCpart[e1];

	if(debug_ >=2) cout<<"convVtx pht1 : "<<convVtxPhtGenpi0[nGenpi0][0] <<" "<< convVtxPhtGenpi0[nGenpi0][1] <<" "<< convVtxPhtGenpi0[nGenpi0][2] <<endl; 
	
      }else if( e2 >0){
	convVtxPhtGenpi0[nGenpi0][0] = vtxXMCpart[e2];
	convVtxPhtGenpi0[nGenpi0][1] = vtxYMCpart[e2];
	convVtxPhtGenpi0[nGenpi0][2] = vtxZMCpart[e2];
      }else{
	convVtxPhtGenpi0[nGenpi0][0] = 0; 
	convVtxPhtGenpi0[nGenpi0][1] = 0; 
	convVtxPhtGenpi0[nGenpi0][2] = 0; 
      }
      


      e1 = convPhtMCpart[indph2][1]; 
      e2 = convPhtMCpart[indph2][2]; 
      if( e1 >0){
	convPht2Genpi0[nGenpi0][0] = eMCpart[e1];
	convPht2Genpi0[nGenpi0][1] = etaMCpart[e1];
	convPht2Genpi0[nGenpi0][2] = phiMCpart[e1];
      }
      if( e2 >0){
	convPht2Genpi0[nGenpi0][3] = eMCpart[e2];
	convPht2Genpi0[nGenpi0][4] = etaMCpart[e2];
	convPht2Genpi0[nGenpi0][5] = phiMCpart[e2];
      }
      if( e1 >0){
	convVtxPhtGenpi0[nGenpi0][3] = vtxXMCpart[e1];
	convVtxPhtGenpi0[nGenpi0][4] = vtxYMCpart[e1];
	convVtxPhtGenpi0[nGenpi0][5] = vtxZMCpart[e1];
      }else if( e2 >0){
	convVtxPhtGenpi0[nGenpi0][3] = vtxXMCpart[e2];
	convVtxPhtGenpi0[nGenpi0][4] = vtxYMCpart[e2];
	convVtxPhtGenpi0[nGenpi0][5] = vtxZMCpart[e2];
      }else {
	convVtxPhtGenpi0[nGenpi0][3] = 0; 
	convVtxPhtGenpi0[nGenpi0][4] = 0; 
	convVtxPhtGenpi0[nGenpi0][5] = 0;
      }
      
      
      /// dr2phGenpi0[nGenpi0] = dr; 
      
      if(debug_>=2){
	TLorentzVector v[2]; 
	int indp[2]={indph1,indph2};
	for( int n =0; n<2; n++){
	  float px = pxMCpart[indp[n]];
	  float py = pyMCpart[indp[n]];
	  float pz = pzMCpart[indp[n]];
	  float e = eMCpart[indp[n]];
	  v[n].SetXYZT(px,py,pz,e);
	}
	////
	TLorentzVector vv = v[0]+v[1];
	cout<<" pi0: "<<pidMCpart[indmc]<<" "<<etaMCpart[indmc]<<" "<<phiMCpart[indmc]<<" "<<" "<<vv.M()<<" "<<vv.Eta()<<" "<<vv.Phi()<<endl;
      }
      

      if(debug_ >=2) cout<<"checkpi0conv: "<< nGenpi0 <<" "<<indph1<<" "<<indph2<<" "<< eGenpi0[nGenpi0] <<" "<< ePhtGenpi0[nGenpi0][0] <<" "<< isConvPhtGenpi0[nGenpi0][0] <<" "<<ePhtGenpi0[nGenpi0][1] <<" "<< isConvPhtGenpi0[nGenpi0][1] <<endl;
      
      
      nGenpi0++; 
      
      
  }
  
  
  

  if(debug_>=1) cout<<"ngenpi0gg: "<< npi0all <<" "<< int(indpi0Gen.size())<<" "<< nGenpi0 <<endl;  
  
  

  nGeneta = 0;
  
  
    
  /// for( int j=0; j<int(indetaGenEtaSorted.size()) && nGeneta<MAXGenPIZ; j++){
  for( int j=0; j<int(indeta_eta.size()) && nGeneta<MAXGenPIZ; j++){
    int k = indeta_eta[j].second; 
    
    //      int k = indetaGenEtaSorted[j];
    
      int indmc = indetaGen[k];
      int indph1 = indpht1etaGen[k];
      int indph2 = indpht2etaGen[k];
    
      eGeneta[nGeneta] = eMCpart[indmc];
      //etaGeneta[nGeneta] = ecalEta(etaMCpart[indmc],vtzMCpart[indmc]/10.0,sqrt(vtxMCpart[indmc]*vtxMCpart[indmc]+vtyMCpart[indmc]*vtyMCpart[indmc])/10.0);
      phiGeneta[nGeneta] = phiMCpart[indmc];
      etaGeneta[nGeneta] = etaMCpart[indmc];
      //etGeneta[nGeneta] = eGeneta[nGeneta]*sin(2*atan(exp(-etaGeneta[nGeneta])));
    
      ////txtout<<"entry: "<<entry<<" "<<j<<" eta et/eta/phi: "<<etGeneta[nGeneta]<<" "<< etaGeneta[nGeneta]<<" " <<phiGeneta[nGeneta]<<endl;
    
    
      //let's do those later.keep the orignal and vtx info
      
    
      vtxGeneta[nGeneta][0] = vtxXMCpart[indmc];
      vtxGeneta[nGeneta][1] = vtxYMCpart[indmc];
      vtxGeneta[nGeneta][2] = vtxZMCpart[indmc];
      vtxPhtGeneta[nGeneta][0] = vtxXMCpart[indph1];
      vtxPhtGeneta[nGeneta][1] = vtxYMCpart[indph1];
      vtxPhtGeneta[nGeneta][2] = vtxZMCpart[indph1];

      
      barcodeMomGeneta[nGeneta] = getMotherIndex(indmc);
      pidMomGeneta[nGeneta]= pidMCpart[barcodeMomGeneta[nGeneta]];
      
      

      ePhtGeneta[nGeneta][0] = eMCpart[indph1];
      ePhtGeneta[nGeneta][1] = eMCpart[indph2];
      etaPhtGeneta[nGeneta][0] = etaMCpart[indph1];
      etaPhtGeneta[nGeneta][1] = etaMCpart[indph2]; 
      phiPhtGeneta[nGeneta][0] = phiMCpart[indph1];
      phiPhtGeneta[nGeneta][1] = phiMCpart[indph2];
    
      
      
      isConvPhtGeneta[nGeneta][0] = convPhtMCpart[indph1][0];
      isConvPhtGeneta[nGeneta][1] = convPhtMCpart[indph2][0];

      int e1 = convPhtMCpart[indph1][1]; 
      int e2 = convPhtMCpart[indph1][2]; 
      
      if( e1 >0){
	convPht1Geneta[nGeneta][0] = eMCpart[e1];
	convPht1Geneta[nGeneta][1] = etaMCpart[e1];
	convPht1Geneta[nGeneta][2] = phiMCpart[e1];
      }else{
	convPht1Geneta[nGeneta][0] = 0; 
	convPht1Geneta[nGeneta][1] = 0; 
	convPht1Geneta[nGeneta][2] = 0; 
      }
      if( e2 >0){
	convPht1Geneta[nGeneta][3] = eMCpart[e2];
	convPht1Geneta[nGeneta][4] = etaMCpart[e2];
	convPht1Geneta[nGeneta][5] = phiMCpart[e2];
      }else{
	convPht1Geneta[nGeneta][3] = 0; 
	convPht1Geneta[nGeneta][4] = 0; 
	convPht1Geneta[nGeneta][5] = 0; 
      }

      
      if( e1 >0){
	convVtxPhtGeneta[nGeneta][0] = vtxXMCpart[e1];
	convVtxPhtGeneta[nGeneta][1] = vtxYMCpart[e1];
	convVtxPhtGeneta[nGeneta][2] = vtxZMCpart[e1];
      }else if( e2 >0){
	convVtxPhtGeneta[nGeneta][0] = vtxXMCpart[e2];
	convVtxPhtGeneta[nGeneta][1] = vtxYMCpart[e2];
	convVtxPhtGeneta[nGeneta][2] = vtxZMCpart[e2];
      }else{
	convVtxPhtGeneta[nGeneta][0] = 0; 
	convVtxPhtGeneta[nGeneta][1] = 0; 
	convVtxPhtGeneta[nGeneta][2] = 0; 
      }


      e1 = convPhtMCpart[indph2][1]; 
      e2 = convPhtMCpart[indph2][2]; 
      if( e1 >0){
	convPht2Geneta[nGeneta][0] = eMCpart[e1];
	convPht2Geneta[nGeneta][1] = etaMCpart[e1];
	convPht2Geneta[nGeneta][2] = phiMCpart[e1];
      }else{
	convPht2Geneta[nGeneta][0] = 0; 
	convPht2Geneta[nGeneta][1] = 0; 
	convPht2Geneta[nGeneta][2] = 0; 
      }
      if( e2 >0){
	convPht2Geneta[nGeneta][3] = eMCpart[e2];
	convPht2Geneta[nGeneta][4] = etaMCpart[e2];
	convPht2Geneta[nGeneta][5] = phiMCpart[e2];
      }else{
	convPht2Geneta[nGeneta][3] = 0; 
	convPht2Geneta[nGeneta][4] = 0; 
	convPht2Geneta[nGeneta][5] = 0; 
      }
      
      
      if( e1 >0){
	convVtxPhtGeneta[nGeneta][3] = vtxXMCpart[e1];
	convVtxPhtGeneta[nGeneta][4] = vtxYMCpart[e1];
	convVtxPhtGeneta[nGeneta][5] = vtxZMCpart[e1];
      }else if( e2 >0){
	convVtxPhtGeneta[nGeneta][3] = vtxXMCpart[e2];
	convVtxPhtGeneta[nGeneta][4] = vtxYMCpart[e2];
	convVtxPhtGeneta[nGeneta][5] = vtxZMCpart[e2];
      }else{
	convVtxPhtGeneta[nGeneta][3] = 0; 
	convVtxPhtGeneta[nGeneta][4] = 0; 
	convVtxPhtGeneta[nGeneta][5] = 0; 
      }



      ////dr2phGeneta[nGeneta] = dr; 
      
      
      nGeneta++; 
      
  }

  
  if(debug_>=1) cout<<"ngenetagg: "<< netaall <<" "<< int(indetaGen.size())<<" "<< nGeneta <<endl; 
  
  
  npizallgen = npi0all; 
  netaallgen = netaall; 
  
  
}


////transform eta ( z, pho), to eta at ecal ( w.r.t 0,0,0,)
double RecoAnalyzer::ecalEta(double EtaParticle ,double Zvertex, double RhoVertex){
  
  
  //  const Double_t PI    = 3.1415927;
  double PI    = acos(-1);
  
  //---Definitions for ECAL
  double R_ECAL           = 136.5;
  double Z_Endcap         = 328.0;
  double etaBarrelEndcap  = 1.479; 

  if (EtaParticle!= 0.)
    {
      double Theta = 0.0  ;
      double ZEcal = (R_ECAL-RhoVertex)*sinh(EtaParticle)+Zvertex;
      
      if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
      if(Theta<0.0) Theta = Theta+PI;

      double ETA = - log(tan(0.5*Theta));
      
      if( fabs(ETA) > etaBarrelEndcap )
	{
	  double Zend = Z_Endcap ;
	  if(EtaParticle<0.0 )  Zend = -Zend ;
	  double Zlen = Zend - Zvertex ;
	  double RR = Zlen/sinh(EtaParticle);
	  Theta = atan((RR+RhoVertex)/Zend);
	  if(Theta<0.0) Theta = Theta+PI;
	  ETA = - log(tan(0.5*Theta));
	}
      return ETA;
    }
  else
    {
      return EtaParticle;
    }
}



double  RecoAnalyzer::ecalPhi(double phi,double x0,double y0){
  
  //double R_ECAL = 136.5; ///cm 
  double r = 136.5; 

  double r0 = sqrt(x0*x0 + y0*y0);
  
  if(r0<1E-5) return phi; 
  
  if( r0 >= r){
    cout<<"warning. ecalPhi vtx outside ecal return input" << r0 <<" "<< r <<endl;
    return phi; 
  }
  
  double theta0 ;
  if(fabs(y0)>0) theta0= y0/fabs(y0) * acos(x0/r0);
  else theta0 = acos(x0/r0);
  
  ///  cout<<theta0<<" "<<phi<<endl;
  
  double theta = phi + asin( r0/r *sin(theta0-phi));

  //phinorm2(theta);
  double PI    = acos(-1);
  while ( theta < -PI) theta += PI; 
  while ( theta > PI) theta -= PI; 
  
  return theta; 
  
  
}



void RecoAnalyzer::calcPairPhoton(float en[],float eta[],float phi[],float res[]){
  
  TLorentzVector vpht[2];
  
  TLorentzVector vpair; 
  
  
  for( int j= 0; j<2; j++){
    
    float e = en[j];
    float px = e * sin ( 2*atan(exp(-eta[j]))) * cos(phi[j]);
    float py = e * sin ( 2*atan(exp(-eta[j]))) * sin(phi[j]);
    float pz = e * cos ( 2*atan(exp(-eta[j]))) ;
    
    vpht[j].SetXYZT(px,py,pz,e);
    
  }
  
  vpair = vpht[0] + vpht[1];

  res[0] = vpair.M();
  res[1] = vpair.Eta();
  res[2] = vpair.Phi();
  res[3] = vpair.Pt();
  
  
  
  
}






void RecoAnalyzer::matchGenpi0Pht(float etapair, float phipair,float eta1,float phi1,float eta2,float phi2,float res[]){
  
  float eta_2pht[2] ={eta1,eta2};
  float phi_2pht[2] ={phi1,phi2};


  float zpht;
  float rpht; 
  float xpht; 
  float ypht; 
  
  for( int j=0; j<30; j++){
    res[j] = 10;
  }
  
  float drsum_twopht_min = 1; 
  int indmin = -1; 
  int indga1 = 0; ////the index of gen-photon for eta1,phi1
  float dr_twopht[2]={1,1};
  
  for( int j=0; j<nGenpi0; j++){
    
    float drsum = 0; 
    ///first of all see which photon is from which gen phton
    ///2-2 match,deltaR 
    float allDR[4]; 
    int kk = 0; 
    
    
    for( int k1 =0; k1<2; k1++){
      
      float eta = eta_2pht[k1];
      float phi = phi_2pht[k1];
      
      for( int k2 = 0; k2<2; k2++){

	zpht = vtxPhtGenpi0[j][2]; 
	rpht = sqrt( vtxPhtGenpi0[j][0] * vtxPhtGenpi0[j][0] + vtxPhtGenpi0[j][1] * vtxPhtGenpi0[j][1]); 
	xpht = vtxPhtGenpi0[j][0]; 
	ypht = vtxPhtGenpi0[j][1]; 
	float etag = ecalEta(etaPhtGenpi0[j][k2],zpht,rpht);
	float phig = ecalPhi(phiPhtGenpi0[j][k2],xpht,ypht);
	allDR[kk] = GetDeltaR(eta,etag,phi,phig);
	kk++;
      }
    }
    float drsum1 = allDR[0] + allDR[3];
    float drsum2 = allDR[1] + allDR[2];
    drsum = drsum1 < drsum2? drsum1: drsum2; 
    if( drsum < drsum_twopht_min ){
      drsum_twopht_min = drsum; 
      indmin = j; 
      if(drsum1<drsum2){
	dr_twopht[0] = allDR[0];
	dr_twopht[1] = allDR[3];
	indga1 = 0;
      }else{
	dr_twopht[0] = allDR[1];
	dr_twopht[1] = allDR[2];
	indga1 = 1; 
      }
    }
  }
  
  ///dr two photon
  res[0] = dr_twopht[0];
  res[1] = dr_twopht[1];
  //dr piz
  res[2] = 1; 
  
  if( indmin>=0){
    
    ///gen z, r 
    zpht = vtxPhtGenpi0[indmin][2];; 
    rpht = sqrt( vtxPhtGenpi0[indmin][0] * vtxPhtGenpi0[indmin][0] + vtxPhtGenpi0[indmin][1] * vtxPhtGenpi0[indmin][1]); 
    float men[2] ={ePhtGenpi0[indmin][indga1],ePhtGenpi0[indmin][1-indga1] };
    float meta[2];
    float mphi[2];
    meta[0] = ecalEta(etaPhtGenpi0[indmin][indga1],zpht,rpht);
    meta[1] = ecalEta(etaPhtGenpi0[indmin][1-indga1],zpht,rpht);
    mphi[0] = ecalPhi(phiPhtGenpi0[indmin][indga1],vtxPhtGenpi0[indmin][0],vtxPhtGenpi0[indmin][1]);
    mphi[1] = ecalPhi(phiPhtGenpi0[indmin][1-indga1],vtxPhtGenpi0[indmin][0],vtxPhtGenpi0[indmin][1]);
    float mres[10];
    calcPairPhoton(men,meta,mphi,mres);
    res[2] = GetDeltaR(mres[1],etapair,mres[2],phipair);
    res[3] = etaPhtGenpi0[indmin][indga1];
    res[4] = phiPhtGenpi0[indmin][indga1];
    res[5] = ePhtGenpi0[indmin][indga1];
    res[6] = etaPhtGenpi0[indmin][1-indga1];
    res[7] = phiPhtGenpi0[indmin][1-indga1];
    res[8] = ePhtGenpi0[indmin][1-indga1];
    res[9] = vtxPhtGenpi0[indmin][0]; ///x
    res[10] = vtxPhtGenpi0[indmin][1]; ///y
    res[11] = vtxPhtGenpi0[indmin][2]; ///z
    res[12] = pidMomGenpi0[indmin];
    res[13] = barcodeMomGenpi0[indmin];


    res[17] = indga1; //// the first cluster's gen photon order 
    res[18] = indmin;  ///matched to gen pi0
  }else{
    res[18] = -1; 
  }
  

}



void RecoAnalyzer::matchGenetaPht(float etapair, float phipair,float eta1,float phi1,float eta2,float phi2,float res[]){
  
  float eta_2pht[2] ={eta1,eta2};
  float phi_2pht[2] ={phi1,phi2};


  float zpht;
  float rpht; 
  float xpht; 
  float ypht; 
  
  for( int j=0; j<30; j++){
    res[j] = 10;
  }
  
  float drsum_twopht_min = 1; 
  int indmin = -1; 
  int indga1 = 0; ////the index of gen-photon for eta1,phi1
  float dr_twopht[2]={1,1};
  
  for( int j=0; j<nGeneta; j++){
    
    float drsum = 0; 
    ///first of all see which photon is from which gen phton
    ///2-2 match,deltaR 
    float allDR[4]; 
    int kk = 0; 
    
    
    for( int k1 =0; k1<2; k1++){
      
      float eta = eta_2pht[k1];
      float phi = phi_2pht[k1];
      
      for( int k2 = 0; k2<2; k2++){

	zpht = vtxPhtGeneta[j][2]; 
	rpht = sqrt( vtxPhtGeneta[j][0] * vtxPhtGeneta[j][0] + vtxPhtGeneta[j][1] * vtxPhtGeneta[j][1]); 
	xpht = vtxPhtGeneta[j][0]; 
	ypht = vtxPhtGeneta[j][1]; 
	float etag = ecalEta(etaPhtGeneta[j][k2],zpht,rpht);
	float phig = ecalPhi(phiPhtGeneta[j][k2],xpht,ypht);
	allDR[kk] = GetDeltaR(eta,etag,phi,phig);
	kk++;
      }
    }
    float drsum1 = allDR[0] + allDR[3];
    float drsum2 = allDR[1] + allDR[2];
    drsum = drsum1 < drsum2? drsum1: drsum2; 
    if( drsum < drsum_twopht_min ){
      drsum_twopht_min = drsum; 
      indmin = j; 
      if(drsum1<drsum2){
	dr_twopht[0] = allDR[0];
	dr_twopht[1] = allDR[3];
	indga1 = 0;
      }else{
	dr_twopht[0] = allDR[1];
	dr_twopht[1] = allDR[2];
	indga1 = 1; 
      }
    }
  }
  
  ///dr two photon
  res[0] = dr_twopht[0];
  res[1] = dr_twopht[1];
  //dr piz
  res[2] = 1; 
  
  if( indmin>=0){
    
    ///gen z, r 
    zpht = vtxPhtGeneta[indmin][2];; 
    rpht = sqrt( vtxPhtGeneta[indmin][0] * vtxPhtGeneta[indmin][0] + vtxPhtGeneta[indmin][1] * vtxPhtGeneta[indmin][1]); 
    float men[2] ={ePhtGeneta[indmin][indga1],ePhtGeneta[indmin][1-indga1] };
    float meta[2];
    float mphi[2];
    meta[0] = ecalEta(etaPhtGeneta[indmin][indga1],zpht,rpht);
    meta[1] = ecalEta(etaPhtGeneta[indmin][1-indga1],zpht,rpht);
    mphi[0] = ecalPhi(phiPhtGeneta[indmin][indga1],vtxPhtGeneta[indmin][0],vtxPhtGeneta[indmin][1]);
    mphi[1] = ecalPhi(phiPhtGeneta[indmin][1-indga1],vtxPhtGeneta[indmin][0],vtxPhtGeneta[indmin][1]);
    float mres[10];
    calcPairPhoton(men,meta,mphi,mres);
    res[2] = GetDeltaR(mres[1],etapair,mres[2],phipair);
    res[3] = etaPhtGeneta[indmin][indga1];
    res[4] = phiPhtGeneta[indmin][indga1];
    res[5] = ePhtGeneta[indmin][indga1];
    res[6] = etaPhtGeneta[indmin][1-indga1];
    res[7] = phiPhtGeneta[indmin][1-indga1];
    res[8] = ePhtGeneta[indmin][1-indga1];
    res[9] = vtxPhtGeneta[indmin][0]; ///x
    res[10] = vtxPhtGeneta[indmin][1]; ///y
    res[11] = vtxPhtGeneta[indmin][2]; ///z
    res[12] = pidMomGeneta[indmin];
    res[13] = barcodeMomGeneta[indmin];
    
    res[17] = indga1; ///the first phton's gen order
    res[18] = indmin;  ///matched to gen pi0
  }else{
    res[18] = -1; 
  }
  
  

}


void RecoAnalyzer::makeEcalUncalibRecHitEB(const EBDataFrame &dataFrame, const double* pedestals,
					 const double* gainRatios,
					 double res[]
					 ){
  
  
    
  
  double chi2_(-1.);
  
  //  double Gain12Equivalent[4]={0,1,2,12};
  double frame[EcalDataFrame::MAXSAMPLES];// will contain the ADC values
  double pedestal =0;     // carries pedestal for gain12 i.e. gainId==1
  
  int gainId0 = 1;        // expected gainId at the beginning of dataFrame 
  int iGainSwitch = 0;    // flags whether there's any gainId other than gainId0
  int GainId= 0;          // stores gainId at every sample
  double maxsample(-1);   // ADC value of maximal ped-subtracted sample
  int imax(-1);           // sample number of maximal ped-subtracted sample
  bool external_pede = false;
  bool isSaturated = 0;   // flag reporting whether gain0 has been found


  if(debug_ >=2) cout<<" EcalDataFrame::MAXSAMPLES "<< EcalDataFrame::MAXSAMPLES <<endl;
  
  // Get time samples checking for Gain Switch and pedestals
  if(pedestals){
    external_pede = true;
    if(dyn_pedestal) { pedestal = (double(dataFrame.sample(0).adc()) + double(dataFrame.sample(1).adc()))/2.;}
    else{ pedestal  = pedestals[0];}
    for(int iSample = 0; iSample < EcalDataFrame::MAXSAMPLES; iSample++) {
	//create frame in adc gain 12 equivalent
	GainId = dataFrame.sample(iSample).gainId();

	// FIX-ME: warning: the vector pedestal is supposed to have in the order G12, G6 and G1
        // if GainId is zero treat it as 3 temporarily to protect against undefined
	// frame will be set to ~max of gain1
	if ( GainId == 0 )
	  { 
	    GainId = 3;
	    isSaturated = 1;
	  }

	if (GainId != gainId0) iGainSwitch = 1;

	if(GainId==gainId0){frame[iSample] = double(dataFrame.sample(iSample).adc())-pedestal ;}
	else {frame[iSample] = (double(dataFrame.sample(iSample).adc())-pedestals[GainId-1])*gainRatios[GainId-1];}
	if( frame[iSample]>maxsample ) {
          maxsample = frame[iSample];
          imax = iSample;
	}
	if( debug_ >=2){
	  cout<<" makeEcalUncalibRecHitEB df: "<< iSample<<" "<< dataFrame.sample(iSample).adc() <<" "<< frame[iSample] <<" "<< GainId <<" "<< imax<<endl;
	}
    }
  }
  else {// pedestal from pre-sample
    external_pede = false;
    pedestal = (double(dataFrame.sample(0).adc()) + double(dataFrame.sample(1).adc()))/2.;

    for(int iSample = 0; iSample < EcalDataFrame::MAXSAMPLES; iSample++) {
      //create frame in adc gain 12 equivalent
      GainId = dataFrame.sample(iSample).gainId();
      //no gain switch forseen if there is no external pedestal
      if ( GainId == 0 ) 
	{
	  GainId = 3;
	  isSaturated = 1;
	}

      frame[iSample] = double(dataFrame.sample(iSample).adc())-pedestal ;
      // if gain has switched but no pedestals are available, no much good you can do...
      if (GainId > gainId0) iGainSwitch = 1;
      if( frame[iSample]>maxsample ) {
	maxsample = frame[iSample];
	imax = iSample;
      }
    } 
  }

  if( (iGainSwitch==1 && external_pede==false) ||  // ... thus you return dummy rechit
      imax ==-1 ){                                 // protect against all frames being <-1
    // return EcalUncalibratedRecHit( dataFrame.id(), -1., -100., -1. , -1.);
    res[0] = -1; 
    res[1] = -100; 
    res[2] = -1; 
    res[3] = -1; 
        
    
    return; 
    
  }
  
  InitFitParameters(frame, imax);
  

  if(debug_ >=2) cout<<"PerformAnalyticFit now.."<<endl;

  chi2_ = PerformAnalyticFit(frame,imax);
  uint32_t flags = 0;
  if (isSaturated) flags = EcalUncalibratedRecHit::kSaturated;


  res[0] = fAmp_max_; 
  res[1] = pedestal+fPed_max_; 
  res[2] = fTim_max_ - 5; 
  res[3] = chi2_; 
  res[4] = flags; 
  
  

  /*    std::cout << "separate fits\nA: " << fAmp_max_  << ", ResidualPed: " <<  fPed_max_
              <<", pedestal: "<<pedestal << ", tPeak " << fTim_max_ << std::endl;
  */
  //return EcalUncalibratedRecHit( dataFrame.id(),fAmp_max_, pedestal+fPed_max_, fTim_max_ - 5 , chi2_, flags );
}



// void RecoAnalyzer::getDigiSample_v1(DetId detid, double adc[]){
//   EBDetId det = (EBDetId)detid;

//   /////if(debug_ >=2) cout<<"getDigiSample_v1; find.. "<< det.ieta()<<"_" << det.iphi()<<endl;
  
//   for(int j=0; j< 15; j++){
//     adc[j] = 0; 
//   }
  
//   std::vector<EBDetId>::const_iterator idd = find(alldetIdEBDigiv1.begin(),alldetIdEBDigiv1.end(),det);
//   if( idd != alldetIdEBDigiv1.end()){
//     int nn = int(idd - alldetIdEBDigiv1.begin());
//     EBDigiv1Collection::const_iterator itdg = ebDigisv1->begin() + nn;
//     const EcalDataFramev1 &dataFrame = *itdg;
//     ///    if( debug_ >=2) cout<<"getDigiSample_v1 "<< det.ieta()<<"_" << det.iphi()<<" sample_v1: "; 
//     for(int iSample = 0; iSample < 10; iSample++) {
//       //create frame in adc gain 12 equivalent
//       //GainId = dataFrame.sample(iSample).gainId();
//       adc[iSample] = dataFrame.sample(iSample).adc(); 
//       ///if(debug_>=2) cout<< adc[iSample]<<" "; 
//     }
//     //// if( debug_ >=2) cout<<endl;
    
//   }else{
//     cout<<"warning_getDigiSample_v1_findnothing.. "<< det.ieta()<<" "<< det.iphi()<<" "<< evtNumber<<" "<< runNumber<<" "<< lumiBlock<<endl;
//   }
  
// }


///return 10 samples of each digi

void RecoAnalyzer::getDigiSample(DetId detid, double adc[]){
  
  EBDetId det = (EBDetId)detid; 
  std::vector<EBDetId>::const_iterator idd = find(alldetIdEBDigi.begin(),alldetIdEBDigi.end(),det);
  
  double gainRatios[3];
  double pedVec[3];
  double pedestal =0;     // carries pedestal for gain12 i.e. gainId==1
  
  int gainId0 = 1;        // expected gainId at the beginning of dataFrame 
  int iGainSwitch = 0;    // flags whether there's any gainId other than gainId0
  int GainId= 0;          // stores gainId at every sample
  double frame[EcalDataFrame::MAXSAMPLES];// will contain the ADC values
  bool isSaturated = 0;

  const EcalGainRatioMap& gainMap = pRatio.product()->getMap(); // map of gain ratios
  EcalGainRatioMap::const_iterator gainIter; // gain iterator
  EcalMGPAGainRatio aGain; // gain object for a single xtal

  const EcalPedestalsMap & pedMap = pedHandle.product()->getMap(); // map of pedestals
  EcalPedestalsMapIterator pedIter; // pedestal iterator
  EcalPedestals::Item aped; // pedestal object for a single xtal
  

  if( idd != alldetIdEBDigi.end()){
    int nn = int(idd - alldetIdEBDigi.begin());
    EBDigiCollection::const_iterator itdg = ebDigis->begin() + nn; 
      gainIter = gainMap.find(itdg->id());
      if( gainIter != gainMap.end() ) {
	aGain = (*gainIter);
      } else {
	cout<< "EcalUncalibRecHitWorkerFixedAlphaBetaFit" << "error!! could not find gain ratios for channel: ";
	return ; 
      }
      gainRatios[0] = 1.;
      gainRatios[1] = aGain.gain12Over6();
      gainRatios[2] = aGain.gain6Over1()*aGain.gain12Over6();
      
      ///ped from pedMap
      pedIter = pedMap.find(itdg->id());
      if( pedIter != pedMap.end() ) {
	aped = (*pedIter); ///will crash here
      }else {
	cout<< "EcalUncalibRecHitWorkerFixedAlphaBetaFit" << "error!! could not find ped map for channel: ";
	return ; 
      }
      pedVec[0] = aped.mean_x12;
      pedVec[1] = aped.mean_x6;
      pedVec[2] = aped.mean_x1;
      
      const EBDataFrame &dataFrame = *itdg; 
      ///dataframe
      if(dyn_pedestal) { pedestal = (double(dataFrame.sample(0).adc()) + double(dataFrame.sample(1).adc()))/2.;}
      //else{ pedestal  = pedestals[0];}
      else{ pedestal  = pedVec[0];}
      
      for(int iSample = 0; iSample < EcalDataFrame::MAXSAMPLES; iSample++) {
	//create frame in adc gain 12 equivalent
	GainId = dataFrame.sample(iSample).gainId();
	if ( GainId == 0 ){ 
	  GainId = 3;
	  isSaturated = 1;
	}
	if (GainId != gainId0) iGainSwitch = 1;

	//if(GainId==gainId0){frame[iSample] = double(dataFrame.sample(iSample).adc())-pedestal ;}
	if(GainId==gainId0){frame[iSample] = double(dataFrame.sample(iSample).adc());}
	else {frame[iSample] = (double(dataFrame.sample(iSample).adc())-pedVec[GainId-1])*gainRatios[GainId-1];}
	
	adc[iSample] = frame[iSample];
	
      }
  } 
  adc[10] = int(isSaturated); //check if saturated
  adc[11] = iGainSwitch;
  
}


///return 10 samples of each digi

void RecoAnalyzer::getDigiSample_EndCap(DetId detid, double adc[]){
  
  EEDetId det = (EEDetId)detid; 
  std::vector<EEDetId>::const_iterator idd = find(alldetIdEEDigi.begin(),alldetIdEEDigi.end(),det);
  
  double gainRatios[3];
  double pedVec[3];
  double pedestal =0;     // carries pedestal for gain12 i.e. gainId==1
  
  int gainId0 = 1;        // expected gainId at the beginning of dataFrame 
  int iGainSwitch = 0;    // flags whether there's any gainId other than gainId0
  int GainId= 0;          // stores gainId at every sample
  double frame[EcalDataFrame::MAXSAMPLES];// will contain the ADC values
  bool isSaturated = 0;

  const EcalGainRatioMap& gainMap = pRatio.product()->getMap(); // map of gain ratios
  EcalGainRatioMap::const_iterator gainIter; // gain iterator
  EcalMGPAGainRatio aGain; // gain object for a single xtal

  const EcalPedestalsMap & pedMap = pedHandle.product()->getMap(); // map of pedestals
  EcalPedestalsMapIterator pedIter; // pedestal iterator
  EcalPedestals::Item aped; // pedestal object for a single xtal
  
  

  if( idd != alldetIdEEDigi.end()){
    int nn = int(idd - alldetIdEEDigi.begin());
    EEDigiCollection::const_iterator itdg = eeDigis->begin() + nn; 
    gainIter = gainMap.find(itdg->id());
    if( gainIter != gainMap.end() ) {
      aGain = (*gainIter);
    } else {
      cout<< "EcalUncalibRecHitWorkerFixedAlphaBetaFit" << "error!! could not find gain ratios for channel: ";
      return ; 
    }
    gainRatios[0] = 1.;
    gainRatios[1] = aGain.gain12Over6();
    gainRatios[2] = aGain.gain6Over1()*aGain.gain12Over6();
    
    ///ped from pedMap
    pedIter = pedMap.find(itdg->id());
    if( pedIter != pedMap.end() ) {
      aped = (*pedIter); ///will crash here
    }else {
      cout<< "EcalUncalibRecHitWorkerFixedAlphaBetaFit" << "error!! could not find ped map for channel: ";
      return ; 
    }
    pedVec[0] = aped.mean_x12;
    pedVec[1] = aped.mean_x6;
    pedVec[2] = aped.mean_x1;
    
    const EEDataFrame &dataFrame = *itdg; 
    ///dataframe
    if(dyn_pedestal) { pedestal = (double(dataFrame.sample(0).adc()) + double(dataFrame.sample(1).adc()))/2.;}
    //else{ pedestal  = pedestals[0];}
    else{ pedestal  = pedVec[0];}
    
    for(int iSample = 0; iSample < EcalDataFrame::MAXSAMPLES; iSample++) {
      //create frame in adc gain 12 equivalent
      GainId = dataFrame.sample(iSample).gainId();
      if ( GainId == 0 ){ 
	GainId = 3;
	isSaturated = 1;
      }
      if (GainId != gainId0) iGainSwitch = 1;
      //if(GainId==gainId0){frame[iSample] = double(dataFrame.sample(iSample).adc())-pedestal ;}
      if(GainId==gainId0){frame[iSample] = double(dataFrame.sample(iSample).adc()) ;}
      else {frame[iSample] = (double(dataFrame.sample(iSample).adc())-pedVec[GainId-1])*gainRatios[GainId-1];}
      
      adc[iSample] = frame[iSample];
      
    }
  }
  
  adc[10] = int(isSaturated); //check if saturated
  adc[11] = iGainSwitch; 
  
  
}





void RecoAnalyzer::runEcalUncalibRecHitFixedAlphaBetaFitCluster(vector<DetId> id_clus, double res[]){


  double frame[EcalDataFrame::MAXSAMPLES];// will contain the ADC values
  double pedestal =0;     // carries pedestal for gain12 i.e. gainId==1
  
  int gainId0 = 1;        // expected gainId at the beginning of dataFrame 
  int iGainSwitch = 0;    // flags whether there's any gainId other than gainId0
  int GainId= 0;          // stores gainId at every sample

  const EcalGainRatioMap& gainMap = pRatio.product()->getMap(); // map of gain ratios
  EcalGainRatioMap::const_iterator gainIter; // gain iterator
  EcalMGPAGainRatio aGain; // gain object for a single xtal

  const EcalPedestalsMap & pedMap = pedHandle.product()->getMap(); // map of pedestals
  EcalPedestalsMapIterator pedIter; // pedestal iterator
  EcalPedestals::Item aped; // pedestal object for a single xtal
  

  double gainRatios[3];
  
  double frame_tot[EcalDataFrame::MAXSAMPLES]; 
  for(int j=0; j< EcalDataFrame::MAXSAMPLES;j++){
    frame_tot[j] = 0; 
  }
  
  double pedVec[3];
  const EcalIntercalibConstantMap& icalMap = icalp->getMap(); 
  bool isSaturated = 0;
  
  for(int j=0; j< int( id_clus.size());j++){
    EBDetId det = (EBDetId)id_clus[j];
    std::vector<EBDetId>::const_iterator idd = find(alldetIdEBDigi.begin(),alldetIdEBDigi.end(),det);
    if( idd != alldetIdEBDigi.end()){
      int nn = int(idd - alldetIdEBDigi.begin());
      EBDigiCollection::const_iterator itdg = ebDigis->begin() + nn; 
      gainIter = gainMap.find(itdg->id());
      if( gainIter != gainMap.end() ) {
	aGain = (*gainIter);
      } else {
	cout<< "EcalUncalibRecHitWorkerFixedAlphaBetaFit" << "error!! could not find gain ratios for channel: ";
	return ; 
      }
      gainRatios[0] = 1.;
      gainRatios[1] = aGain.gain12Over6();
      gainRatios[2] = aGain.gain6Over1()*aGain.gain12Over6();


      if(debug_>=2) cout<<" runEcalUncalibRecHitFixedAlphaBetaFitClusterGainRatio : " << j <<" "<< gainRatios[1] <<" "<< gainRatios[2] <<endl;

      ///ped from pedMap
      pedIter = pedMap.find(itdg->id());
      if( pedIter != pedMap.end() ) {
	aped = (*pedIter); ///will crash here
      }else {
	cout<< "EcalUncalibRecHitWorkerFixedAlphaBetaFit" << "error!! could not find ped map for channel: ";
	return ; 
      }
      pedVec[0] = aped.mean_x12; ///200
      pedVec[1] = aped.mean_x6; ///200 
      pedVec[2] = aped.mean_x1; ///200
      
      if(debug_>=2) cout<<" runEcalUncalibRecHitFixedAlphaBetaFitClusterPED : " << j <<" "<< pedVec[0]<<" "<< pedVec[1] <<" "<<pedVec[2] <<endl;
      
      
      
      const EBDataFrame &dataFrame = *itdg; 
      ///dataframe
      if(dyn_pedestal) { pedestal = (double(dataFrame.sample(0).adc()) + double(dataFrame.sample(1).adc()))/2.;}
      //else{ pedestal  = pedestals[0];}
      else{ pedestal  = pedVec[0];}

      DetId detid(itdg->id());
      // first intercalibration constants
      EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detid);
      EcalIntercalibConstant icalconst = 1;
      if( icalit!=icalMap.end() ) {
	icalconst = (*icalit);
	if(debug_ >=2) cout<<"icalconst: "<< icalconst <<endl;

	if( icalconst <0 || icalconst >5){
	  cout<<"warning... icalconst..." << icalconst <<endl;
	}
      } else {
	cout<<"warning... no inter-calib found.."<<endl; 
      }
      
      if(debug_ >=2) cout<<"frame_cluster: "; 
      
      for(int iSample = 0; iSample < EcalDataFrame::MAXSAMPLES; iSample++) {
	//create frame in adc gain 12 equivalent
	GainId = dataFrame.sample(iSample).gainId();
	if ( GainId == 0 ){ 
	  GainId = 3;
	  isSaturated = 1;
	}
	if (GainId != gainId0) iGainSwitch = 1;
	if(GainId==gainId0){frame[iSample] = double(dataFrame.sample(iSample).adc())-pedestal ;}
	//	else {frame[iSample] = (double(dataFrame.sample(iSample).adc())-pedestals[GainId-1])*gainRatios[GainId-1];}
	else {frame[iSample] = (double(dataFrame.sample(iSample).adc())-pedVec[GainId-1])*gainRatios[GainId-1];}
	
	frame_tot[iSample] += frame[iSample] * icalconst; 
	
	if(debug_ >=2) cout<< frame[iSample]<<" "; 
	
      }///end of loop over sample
      if(debug_ >=2) cout<<endl; 
      
    }
    
  }///loop over clustered digi
    
  SetAlphaBeta(alphaEB_, betaEB_);
  SetMinAmpl(8);
  int imax = -1; 
  double maxsample(-1); 

  if(dyn_pedestal) { pedestal = (double(frame_tot[0]) + double(frame_tot[1]))/2.;}
  else{ pedestal  = pedVec[0];}
  if(debug_ >=2) cout<<"frame_cluster_tot: ";
  for(int iSample = 0; iSample < EcalDataFrame::MAXSAMPLES; iSample++) {
    if(debug_ >=2) cout<< frame_tot[iSample]<<" "; 
    if( frame_tot[iSample]>maxsample ) {
      maxsample = frame_tot[iSample];
      imax = iSample;
    }
  }
  if(debug_ >=2) cout<<endl; 
  
  
  InitFitParameters(frame_tot, imax);
  if(debug_ >=2) cout<<"PerformAnalyticFit of cluster now.."<<endl;
  double chi2_ = PerformAnalyticFit(frame_tot,imax);
  uint32_t flags = 0;
  if (isSaturated) flags = EcalUncalibratedRecHit::kSaturated;
  

  
  float adctoGev = 0.03; 
  if ( id_clus[0].subdetId() == EcalEndcap ) {
    //    rechitMaker_->setADCToGeVConstant( float(agc->getEEValue()) );
    adctoGev = float (agc->getEEValue()); 
  } else {
    //    rechitMaker_->setADCToGeVConstant( float(agc->getEBValue()) );
    adctoGev = float (agc->getEBValue()); 
    if(debug_>=2) cout<<"agc_eb: "<< agc->getEBValue() <<endl;
  }
  
  
  
  res[0] = fAmp_max_ * adctoGev; 
  res[1] = (pedestal+fPed_max_) * adctoGev; 
  res[2] = (fTim_max_ - 5) * 25; 
  res[3] = chi2_; 
  res[4] = flags; 
  
  

  
}

///return amplitude of ADC count
void RecoAnalyzer::runEcalUncalibRecHitFixedAlphaBetaFit(const EcalDigiCollection::const_iterator & itdg,double res[]){
  
  const EcalGainRatioMap& gainMap = pRatio.product()->getMap(); // map of gain ratios
  EcalGainRatioMap::const_iterator gainIter; // gain iterator
  EcalMGPAGainRatio aGain; // gain object for a single xtal
  
  const EcalPedestalsMap & pedMap = pedHandle.product()->getMap(); // map of pedestals
  EcalPedestalsMapIterator pedIter; // pedestal iterator
  
  EcalPedestals::Item aped; // pedestal object for a single xtal
  
  DetId detid( itdg->id() );

  if( debug_ >=2) cout<<"now  in runEcalUncalibRecHitFixedAlphaBetaFit " <<endl;
  
  res[0] = -1; 
  res[1] = -100; 
  res[2] = -1; 
  res[3] = -1; 


  pedIter = pedMap.find(itdg->id());

  if( debug_ >=2) cout<<"now  in runEcalUncalibRecHitFixedAlphaBetaFit 1" <<endl;

  if( pedIter != pedMap.end() ) {
    
    if( debug_ >=2) cout<<"now  in runEcalUncalibRecHitFixedAlphaBetaFit 11" <<endl;

    aped = (*pedIter); ///will crash here

    if( debug_ >=2) cout<<"now  in runEcalUncalibRecHitFixedAlphaBetaFit 2" <<endl;

  } else {
    cout <<"EcalUncalibRecHitWorkerFixedAlphaBetaFit" << "error!! could not find pedestals for channel ! " <<endl;
    if ( detid.subdetId() == EcalBarrel ) {
      cout << EBDetId( detid )<<endl;
    }else{
      cout << EEDetId( detid ) <<endl;
    }
    return ; 
  }

  if( debug_ >=2) cout<<"now  in runEcalUncalibRecHitFixedAlphaBetaFit 3" <<endl;

  double pedVec[3];
  pedVec[0] = aped.mean_x12;
  pedVec[1] = aped.mean_x6;
  pedVec[2] = aped.mean_x1;
  
  //  pedVec[0] = 200;
  //  pedVec[1] = 200;
  // pedVec[2] = 200;
  
  
  if( debug_ >=2) cout<<"pedVec: "<< pedVec[0]<<" "<< pedVec[1]<<" "<< pedVec[2]<<endl;
  

  gainIter = gainMap.find(itdg->id());
  if( gainIter != gainMap.end() ) {
    aGain = (*gainIter);
  } else {
    cout<< "EcalUncalibRecHitWorkerFixedAlphaBetaFit" << "error!! could not find gain ratios for channel: ";
    if ( detid.subdetId() == EcalBarrel ) {
      cout<<"EcalUncalibRecHitWorkerFixedAlphaBetaFitbarrel" << EBDetId( detid ) <<endl;
    } else {
      cout<<"EcalUncalibRecHitWorkerFixedAlphaBetaFitendcap" << EEDetId( detid ) <<endl;
    }
    return ; 
  }
  double gainRatios[3];
  gainRatios[0] = 1.;
  gainRatios[1] = aGain.gain12Over6();
  gainRatios[2] = aGain.gain6Over1()*aGain.gain12Over6();

  if( debug_ >=2) cout<<"gainVec: "<< gainRatios[0] <<" "<<gainRatios[1] <<" "<<gainRatios[2]<<endl;
  
  
  if ( detid.subdetId() == EcalBarrel ) {
    SetAlphaBeta(alphaEB_, betaEB_);
    SetMinAmpl(8);

    
    if( debug_ >=2) cout<<"now makeEcalUncalibRecHitEB " <<endl;
    
    makeEcalUncalibRecHitEB(*itdg, pedVec, gainRatios,res);
  }else{
    
    SetAlphaBeta(alphaEE_,betaEE_);
    SetMinAmpl(14);
  }
    

  



}











void RecoAnalyzer::SetAlphaBeta( double alpha, double beta){
  fAlpha_ = alpha;
  fBeta_=  beta;
  alfabeta_ = fAlpha_*fBeta_;
}

void RecoAnalyzer::SetMinAmpl( double ampl){
  MinAmpl_ = ampl;
}

double RecoAnalyzer::pulseShapeFunction(double t){
  if( alfabeta_ <= 0 ) return((double)0.);
  double dtsbeta,variable,puiss;
  double dt = t-fTim_max_ ;
  if(dt > -alfabeta_)  {
    dtsbeta=dt/fBeta_ ;
    variable=1.+dt/alfabeta_ ;
    puiss=pow(variable,fAlpha_);
    return fAmp_max_*puiss*exp(-dtsbeta)+fPed_max_ ;
  }
  return  fPed_max_ ;
}


void RecoAnalyzer::InitFitParameters(double* samples, int max_sample){

  // in a first attempt just use the value of the maximum sample 
  fAmp_max_ = samples[max_sample];
  fTim_max_ = max_sample;
  fPed_max_ = 0;

  // amplitude too low for fit to converge
  // timing set correctly is assumed
  if(fAmp_max_ <  MinAmpl_){
    fAmp_max_      = samples[5];
    double sumA    = samples[5]+samples[4]+samples[6];
    if(sumA != 0) { fTim_max_ = 5+(samples[6]-samples[4])/sumA; }
    else{ fTim_max_ = -993; }//-999+6
    doFit_ = false;
  }
  // if timing very badly off, that just use max sample
  else if(max_sample <1 || max_sample > 7)
    {    doFit_ = false;}
  else{
    //y=a*(x-xM)^2+b*(x-xM)+c
    doFit_ = true;
    float a = float(samples[max_sample-1]+samples[max_sample+1]-2*samples[max_sample])/2.;
    if(a==0){doFit_ =false; return;}
    float b = float(samples[max_sample+1]-samples[max_sample-1])/2.;
  
    fTim_max_ = max_sample - b/(2*a);
    fAmp_max_ =  samples[max_sample] - b*b/(4*a); 
  }  
  
} 


float RecoAnalyzer::PerformAnalyticFit(double* samples, int max_sample){
  
  //int fValue_tim_max = max_sample;  
  //! fit electronic function from simulation
  //! parameters fAlpha_ and fBeta_ are fixed and fit is providing the 3 following parameters
  //! the maximum amplitude ( fAmp_max_ ) 
  //! the time of the maximum  ( fTim_max_)
  //| the pedestal (fPed_max_) 	
  
  double chi2=-1 , db[3] ;
  
  if(debug_ >=2) cout<<"before_fit.."<< fAmp_max_ <<" "<<fTim_max_ <<" "<<fPed_max_ <<endl;
  

  
  CLHEP::HepSymMatrix DDD_(3) ;
  CLHEP::HepVector ddd_(3) ;
  
  
  int num_fit_min =(int)(max_sample - fNum_samp_bef_max_ ) ;
  int num_fit_max =(int)(max_sample + fNum_samp_after_max_) ;

  if (num_fit_min<0) num_fit_min = 0 ; 
  //if (num_fit_max>=fNsamples-1) num_fit_max = fNsamples-2 ;
  if (num_fit_max>= EcalDataFrame::MAXSAMPLES) {num_fit_max = EcalDataFrame::MAXSAMPLES - 1 ;}
  
  if(! doFit_ ) {
    LogDebug("EcalUncalibRecHitFixedAlphaBetaAlgo")<<"No fit performed. The amplitude of sample 5 will be used";
    return -1;
  }

  double func,delta ;
  double variation_func_max = 0. ;double variation_tim_max = 0. ; double variation_ped_max = 0. ;


  
  
  //!          Loop on iterations
  for (int iter=0 ; iter < fNb_iter_ ; iter ++) {
    //!          initialization inside iteration loop !
    chi2 = 0. ; //PROD.Zero() ;  DM1.Zero() ;

     for(int i1=0 ; i1<3 ; i1++) {
       ddd_[i1]=0;
	for(int j1=i1 ; j1<3 ; j1++) { 
	  DDD_.fast(j1+1,i1+1) = 0; }
      }

    fAmp_max_ += variation_func_max ;
    fTim_max_ += variation_tim_max ;
    fPed_max_ += variation_ped_max ;
    
    //! Then we loop on samples to be fitted
    for( int i = num_fit_min ; i <= num_fit_max ; i++) {
      //if(i>fsamp_edge_fit && i<num_fit_min) continue ; // remove front edge samples
      //! calculate function to be fitted
      func = pulseShapeFunction( (double)i  ) ;
      //! then calculate derivatives of function to be fitted
      double dt =(double)i - fTim_max_ ;
      if(dt > -alfabeta_)  {      
	double dt_sur_beta = dt/fBeta_ ;
	double variable = (double)1. + dt/alfabeta_ ;
	double expo = exp(-dt_sur_beta) ;	 
	double puissance = pow(variable,fAlpha_) ;
	
	db[0]=un_sur_sigma*puissance*expo ;
	db[1]=fAmp_max_*db[0]*dt_sur_beta/(alfabeta_*variable) ;
      }
      else {
	db[0]=0. ; db[1]=0. ; 
      }
      db[2]=un_sur_sigma ;
      //! compute matrix elements DM1
      for(int i1=0 ; i1<3 ; i1++) {
	for(int j1=i1 ; j1<3 ; j1++) { 
	  //double & fast(int row, int col);
	  DDD_.fast(j1+1,i1+1) += db[i1]*db[j1]; }
      }
      //! compute delta
      delta = (samples[i]-func)*un_sur_sigma ;
      //! compute vector elements PROD
      for(int ii=0 ; ii<3 ;ii++) {ddd_[ii] += delta*db[ii] ;}
      chi2 += delta * delta ;
    }//! end of loop on samples 
    
    int fail=0 ;
    DDD_.invert(fail) ;

    
    if(fail != 0.) {
      //just a guess from the value of the parameters in the previous interaction;
      //printf("wH4PulseFitWithFunction =====> determinant error --> No Fit Provided !\n") ;
      InitFitParameters(samples,max_sample);
      return -101 ;
    }
/*     for(int i1=0 ; i1<3 ; i1++) { */
/*       for(int j1=0 ; j1<3 ; j1++) {  */
/* 	//double & fast(int row, int col); */
/* 	std::cout<<"inverted: "<<DM1[j1][i1]<<std::endl;;} */
/*     } */
/*     std::cout<<"vector temp: "<< temp[0]<<" "<<temp[1]<<" "<<temp[2]<<std::endl; */
    //! compute variations of parameters fAmp_max and fTim_max 
    CLHEP::HepVector PROD = DDD_*ddd_ ;
    //    std::cout<<"vector PROD: "<< PROD[0]<<" "<<PROD[1]<<" "<<PROD[2]<<std::endl;

    // Probably the fastest way to protect against
    // +-inf value in the matrix DDD_ after inversion
    // (which is nevertheless flagged as successfull...)
    if ( isnan( PROD[0] ) ) {
            InitFitParameters(samples,max_sample);
            return -103 ;
    }

    variation_func_max = PROD[0] ;
    variation_tim_max = PROD[1] ;
    variation_ped_max = PROD[2] ;


    if(debug_ >=2) cout<<"fit iteration "<< iter <<" "<< db[0] <<" "<< db[1]<<" "<< db[2] <<" "<< PROD[0] <<" "<< PROD[1] <<" "<< PROD[2] <<" "<< chi2 <<endl;
    
    //chi2 = chi2/((double)nsamp_used - 3.) ;
  }//!end of loop on iterations       
  

  //!   protection again diverging/unstable fit 
  if( variation_func_max > 2000. || variation_func_max < -1000. ) {
    InitFitParameters(samples,max_sample);
    return -102;
  }
  

  //!      results of the fit are calculated
  fAmp_max_ += variation_func_max ;
  fTim_max_ += variation_tim_max ;
  fPed_max_ += variation_ped_max ;


  // protection against unphysical results:
  // ampli mismatched to MaxSample, ampli largely negative, time off preselected range
  if( fAmp_max_>2*samples[max_sample]  ||  fAmp_max_<-100 ||  fTim_max_<0  ||  9<fTim_max_ ) {
    InitFitParameters(samples,max_sample);
    return -104;
  }
  

  //std::cout <<"chi2: "<<chi2<<" ampl: "<<fAmp_max_<<" time: "<<fTim_max_<<" pede: "<<fPed_max_<<std::endl;
  return chi2;
}


// void RecoAnalyzer::makeEcalRecHit(const edm::Event & evt, DetId detid, double res_uncalib[],double res_calib[]){
  
  
//   float adctoGev = 0.03; 
  
//   const EcalIntercalibConstantMap& icalMap = icalp->getMap(); 
  
 
//   if ( detid.subdetId() == EcalEndcap ) {
//     //    rechitMaker_->setADCToGeVConstant( float(agc->getEEValue()) );
//     adctoGev = float (agc->getEEValue()); 
//   } else {
//     //    rechitMaker_->setADCToGeVConstant( float(agc->getEBValue()) );
//     adctoGev = float (agc->getEBValue()); 

//     if(debug_>=2) cout<<"agc_eb: "<< agc->getEBValue() <<endl;
    
//   }
//   // first intercalibration constants
//   EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detid);
//   EcalIntercalibConstant icalconst = 1;
//   if( icalit!=icalMap.end() ) {
//     icalconst = (*icalit);

//     if(debug_ >=2) cout<<"icalconst: "<< icalconst <<endl;
    
//   } else {
//     cout<<"warning... no inter-calib found.."<<endl; 
//   }
//   // get laser coefficient
//   float lasercalib = laser->getLaserCorrection( detid, evt.time());
  
//   // get time calibration coefficient
//   const EcalTimeCalibConstantMap & itimeMap = itime->getMap();  
//   EcalTimeCalibConstantMap::const_iterator itime = itimeMap.find(detid);
//   EcalTimeCalibConstant itimeconst = 0;
//   if( icalit!=icalMap.end() ) {
//     itimeconst = (*itime);
//   } else {
//     cout<<"warning... no time-calib found.."<<endl; 
//   }
//   // make the rechit and put in the output collection
  

//   if(debug_ >=2) cout<<"makeEcalRecHit " << res_uncalib[0] <<" "<< adctoGev <<" "<< icalconst << " "<< lasercalib  <<" "<< itimeconst <<endl; 
  
//   ///energy 
//   res_calib[0] = res_uncalib[0] * icalconst * lasercalib * adctoGev; 
//   //time 
//   res_calib[1] = res_uncalib[2]  * 25  + itimeconst; 
  
  
  
// }


//--------------------------------------------------------------------------------------------------
int RecoAnalyzer::getContainedHits(const std::vector<VertexHit> &hits, double z0, double &chi) 
{
  // Calculate number of hits contained in v-shaped window in cluster y-width vs. z-position.

  int n = 0;
  chi   = 0.;

  for(vector<VertexHit>::const_iterator hit = hits.begin(); hit!= hits.end(); hit++) {
    double p = 2 * fabs(hit->z - z0)/hit->r + 0.5; // FIXME
    if(TMath::Abs(p - hit->w) <= 1.) { 
      chi += fabs(p - hit->w);
      n++;
    }
  }
  return n;
}


//
// member functions
//

// Initialize the internal trigger path representation (names and indices) from the 
// patterns specified in the configuration.
// This needs to be called once at startup, whenever the trigger table has changed
// or in case of paths from eventsetup and IOV changed
//void RecoAnalyzer::init(const edm::TriggerResults & result, const edm::EventSetup& iSetup)




void RecoAnalyzer::makeNxNClusters(const edm::Event &evt, const edm::EventSetup &es,
				   const EcalRecHitCollection *hits, const reco::CaloID::Detectors detector){
  

  ///get status from DB
  edm::ESHandle<EcalChannelStatus> csHandle;
  if ( useDBStatus_ ) es.get<EcalChannelStatusRcd>().get(csHandle);
  const EcalChannelStatus &channelStatus = *csHandle; 
    
  nClus = 0; 
  
  eClus.clear();
  etaClus.clear();
  phiClus.clear();
  thetaClus.clear();
  etClus.clear();
  s4s9Clus.clear();
  s9s25Clus.clear();
  xClus.clear();
  yClus.clear();
  zClus.clear();
  

  RecHitsCluster.clear();
  RecHitsCluster5x5.clear();
  
  
  std::vector<EcalRecHit> seeds;
  
  double clusterSeedThreshold ; 
  if (detector == reco::CaloID::DET_ECAL_BARREL){
    clusterSeedThreshold = clusSeedThr_;
  }else{
    clusterSeedThreshold = clusSeedThrEndCap_; 
  }
  

  for(EcalRecHitCollection::const_iterator itt = hits->begin(); itt != hits->end(); itt++){
    double energy = itt->energy();
    
    if( detector == reco::CaloID::DET_ECAL_BARREL) {
      EBDetId ebd = EBDetId(itt->id());
      // hh_mul_ieta ->Fill(ebd.ieta());
      // hh_mul_iphi ->Fill(ebd.iphi());
      if (energy > clusterSeedThreshold ){
	// hh_mul_ietaSeed ->Fill(ebd.ieta());
	// hh_mul_iphiSeed ->Fill(ebd.iphi());
      }
    }
    
    if( ! checkStatusOfEcalRecHit(channelStatus, *itt) ) continue; 
    if (energy > clusterSeedThreshold ) seeds.push_back(*itt);
  }
  
  
  // get the geometry and topology from the event setup:
  edm::ESHandle<CaloGeometry> geoHandle;
  es.get<CaloGeometryRecord>().get(geoHandle);
  
  const CaloSubdetectorGeometry *geometry_p;
  CaloSubdetectorTopology *topology_p;
  if (detector == reco::CaloID::DET_ECAL_BARREL) {
    geometry_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalBarrel);
    topology_p = new EcalBarrelTopology(geoHandle);
  }else {
    geometry_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalEndcap);
    topology_p = new EcalEndcapTopology(geoHandle); 
  }
  
  const CaloSubdetectorGeometry *geometryES_p;
  geometryES_p = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
    
  
  std::vector<reco::BasicCluster> clusters;
  std::vector<DetId> usedXtals;
  // sort seed according to Energy
  sort(seeds.begin(), seeds.end(), ecalRecHitLess());
  
  
  for (std::vector<EcalRecHit>::iterator itseed=seeds.begin(); itseed!=seeds.end(); itseed++) {
    DetId seed_id = itseed->id();
    std::vector<DetId>::const_iterator usedIds;
    
    ///sed not yet used
    std::vector<DetId>::iterator  itdet = find(usedXtals.begin(),usedXtals.end(),seed_id);
    if(itdet != usedXtals.end()) continue; 
    
    std::vector<DetId> clus_v = topology_p->getWindow(seed_id,clusEtaSize_,clusPhiSize_);	
    float clus_energy = 0; 
    vector<EcalRecHit> RecHitsInWindow;
    vector<EcalRecHit> RecHitsInWindow5x5;
    vector<DetId> DetIdInWindow; 
    std::vector<std::pair<DetId, float> > clus_used;
    
    for (std::vector<DetId>::iterator det=clus_v.begin(); det!=clus_v.end(); det++) {
	DetId detid = *det;
	
	//not yet used 
	std::vector<DetId>::iterator itdet = find(usedXtals.begin(),usedXtals.end(),detid);
	if(itdet != usedXtals.end()) continue; 
	//inside the collection
	EcalRecHitCollection::const_iterator hit  = hits->find(detid); 
	if( hit == hits->end()) continue; 
	
	if( ! checkStatusOfEcalRecHit(channelStatus, *hit) ) continue; 
	
	usedXtals.push_back(detid);
	clus_used.push_back(std::pair<DetId, float>(detid, 1.) );
	clus_energy += hit->energy();
	
	
	RecHitsInWindow.push_back(*hit);
	DetIdInWindow.push_back(detid);
	
	
    } //// end of making one nxn simple cluster
    
    if( clus_energy <= 0 ) continue; 
        
      
    math::XYZPoint clus_pos = posCalculator_.Calculate_Location(clus_used,hits,geometry_p,geometryES_p);
    
    sort(RecHitsInWindow.begin(), RecHitsInWindow.end(), ecalRecHitLess());
    
    if (debug_>=2 ) LogDebug("")<<"nxn_cluster in run "<< evt.id().run()<<" event "<<evt.id().event()<<" energy: "<<clus_energy <<" eta: "<< clus_pos.Eta()<<" phi: "<< clus_pos.Phi()<<" nRecHits: "<< RecHitsInWindow.size()<<endl; 
    
    
    int seedx; /// = seed_id.ieta();
    int seedy; /// = seed_id.iphi();
    
    if (detector == reco::CaloID::DET_ECAL_BARREL) {
      EBDetId ebd = EBDetId(seed_id);
      seedx = ebd.ieta();
      seedy = ebd.iphi();
      convxtalid(seedy,seedx);
      
      // hh_mul_ietaSeedClus->Fill(ebd.ieta());
      // hh_mul_iphiSeedClus->Fill(ebd.iphi());
      
    }else{
      EEDetId eed = EEDetId(seed_id);
      seedx = eed.ix();
      seedy = eed.iy();
    }
    
    
    ///    convxtalid( seed_iphi,seed_ieta);
    float e2x2[4]={0};
    float e3x3 = 0; 
    float e5x5 = 0; 
        
    std::vector<DetId> clus_v5x5; 
    clus_v5x5 = topology_p->getWindow(seed_id,5,5);
    
    for( std::vector<DetId>::const_iterator idItr = clus_v5x5.begin(); idItr != clus_v5x5.end(); idItr++){
      DetId ed = *idItr; 
      EcalRecHitCollection::const_iterator rit =  hits ->find( ed );
      if ( rit == hits ->end() ) continue; 
      if( ! checkStatusOfEcalRecHit(channelStatus, *rit) ) continue; 
      float energy = (*rit).energy();
      e5x5 += energy; 
      
      //check if already in 3x3 , so no need to push_back here
      std::vector<DetId>::const_iterator already3x3 =  find(DetIdInWindow.begin(),DetIdInWindow.end(),ed);
      if( already3x3 != DetIdInWindow.end() ) continue; 
      RecHitsInWindow5x5.push_back( *rit);
    }
    
    
    for(unsigned int j=0; j<RecHitsInWindow.size();j++){
      DetId det = RecHitsInWindow[j].id(); 
      int dx,dy; 
      if (detector == reco::CaloID::DET_ECAL_BARREL) {
	EBDetId ebd = EBDetId(det); 
	int x = ebd.ieta();
	int y = ebd.iphi();
	convxtalid(y,x);
	dx = diff_neta_s(x,seedx);
	dy = diff_nphi_s(y,seedy);
      }else{
	EEDetId eed = EEDetId(det); 
	int x = eed.ix();
	int y = eed.iy();
	dx = x-seedx;
	dy = y-seedy; 
      }
      float energy = RecHitsInWindow[j].energy();
      if(abs(dx)<=1 && abs(dy)<=1){
	if(dx <= 0 && dy <=0) e2x2[0] += energy; 
	if(dx >= 0 && dy <=0) e2x2[1] += energy; 
	if(dx <= 0 && dy >=0) e2x2[2] += energy; 
	if(dx >= 0 && dy >=0) e2x2[3] += energy; 
        e3x3 += energy; 
      }
    }
        
    if (detector == reco::CaloID::DET_ECAL_BARREL ) {
      if( removeSpike_==1){ ///check E9 >2GeV && Emax/E9 > 0.95
	if( clus_energy >2 && e3x3>0 && RecHitsInWindow[0].energy() / e3x3 > 0.95) continue; 
      }
    }
    
    
    if( clusEtaSize_ ==3 && fabs( e3x3 - clus_energy)>0.001){
      cout<<"warning... check e3x2: "<< e3x3 <<" "<< clus_energy<<endl; 
    }
    
    ///Now fill the vector of everything
    eClus.push_back( clus_energy);
    etaClus.push_back( clus_pos.Eta());
    phiClus.push_back( clus_pos.Phi());
    thetaClus.push_back(2. * atan(exp(-clus_pos.Eta())));
    etClus.push_back( clus_energy *  sin(2. * atan(exp(-clus_pos.Eta()))));
    
    float s4max = *max_element( e2x2,e2x2+4); 
    s4s9Clus.push_back( s4max/ e3x3);
    float s9s25 = e5x5 !=0 ? e3x3/e5x5:-1; 
    s9s25Clus.push_back(s9s25);
    xClus.push_back( clus_pos.X());
    yClus.push_back( clus_pos.Y());
    zClus.push_back( clus_pos.Z());
    
    
    RecHitsCluster.push_back(RecHitsInWindow);
    RecHitsCluster5x5.push_back(RecHitsInWindow5x5);
    
    nClus ++; 
    
    
  }
  
  
  if( saveAllPhotonBarrel_ && detector == reco::CaloID::DET_ECAL_BARREL){
    ///save alll 3x3 clusters into the events, and all gen pi0/etas
    n3x3ClusEB = nClus < MAX3x3ClusEB ?  nClus : MAX3x3ClusEB ; 
    for(int j=0; j< n3x3ClusEB; j++){
      
      e3x3ClusEB[j] = eClus[j]; 
      eta3x3ClusEB[j] = etaClus[j];
      phi3x3ClusEB[j] = phiClus[j];
      x3x3ClusEB[j] = xClus[j];
      y3x3ClusEB[j] = yClus[j];
      z3x3ClusEB[j] = zClus[j]; 
      nXt3x3ClusEB[j] = int(RecHitsCluster[j].size()); 
      s4s93x3ClusEB[j] = s4s9Clus[j];
      s9s253x3ClusEB[j] = s9s25Clus[j];

      for( int k =0 ; k< 9; k++){
	if( k<nXt3x3ClusEB[j]){
	  eXt3x3ClusEB[j][k] = RecHitsCluster[j][k].energy();
	  tXt3x3ClusEB[j][k] = RecHitsCluster[j][k].time();
	  EBDetId det = (EBDetId)RecHitsCluster[j][k].id(); 
	  ietaXt3x3ClusEB[j][k] = det.ieta();
	  iphiXt3x3ClusEB[j][k] = det.iphi();
	  laserCorr3x3ClusEB[j][k] = laser->getLaserCorrection( RecHitsCluster[j][k].id(), evt.time());
	}else{
	  eXt3x3ClusEB[j][k]  = 0; 
	  tXt3x3ClusEB[j][k]  = -99; 
	  ietaXt3x3ClusEB[j][k] = -99;
	  iphiXt3x3ClusEB[j][k] = -99; 
	  laserCorr3x3ClusEB[j][k] = -99;
	}
      }
    }
  }
  
  if( saveAllPhotonEndcap_ && detector == reco::CaloID::DET_ECAL_ENDCAP){
    n3x3ClusEE = nClus < MAX3x3ClusEE ? nClus : MAX3x3ClusEE ; 
    for(int j=0; j< n3x3ClusEE; j++){
      e3x3ClusEE[j] = eClus[j]; 
      eta3x3ClusEE[j] = etaClus[j];
      phi3x3ClusEE[j] = phiClus[j];
      x3x3ClusEE[j] = xClus[j];
      y3x3ClusEE[j] = yClus[j];
      z3x3ClusEE[j] = zClus[j]; 
      nXt3x3ClusEE[j] = int(RecHitsCluster[j].size()); 
      s4s93x3ClusEE[j] = s4s9Clus[j];
      //s6s93x3ClusEE[j] = s6s9Clus[j];
      s9s253x3ClusEE[j] = s9s25Clus[j];
      
      for( int k =0 ; k< 9; k++){
	if( k<nXt3x3ClusEE[j]){
	  eXt3x3ClusEE[j][k] = RecHitsCluster[j][k].energy();
	  tXt3x3ClusEE[j][k] = RecHitsCluster[j][k].time();
	  EEDetId det = (EEDetId)RecHitsCluster[j][k].id(); 
	  iyXt3x3ClusEE[j][k] = det.ix();
	  ixXt3x3ClusEE[j][k] = det.iy();
	  izXt3x3ClusEE[j][k] = det.zside();
	  laserCorr3x3ClusEE[j][k] = laser->getLaserCorrection( RecHitsCluster[j][k].id(), evt.time());
	}else{
	  eXt3x3ClusEE[j][k]  = 0; 
	  tXt3x3ClusEE[j][k]  = -99; 
	  ixXt3x3ClusEE[j][k] = -99; 
	  iyXt3x3ClusEE[j][k] = -99; 
	  izXt3x3ClusEE[j][k] = -99; 
	  laserCorr3x3ClusEE[j][k] = -99;
	}
      }
    }
    
  }
  
  
}



void RecoAnalyzer::doSelectionAndFillTree(const edm::Event &evt, const edm::EventSetup &es, const reco::CaloID::Detectors detector, bool rmPiz, int pizEta, double s4s9Cut, double s9s25Cut, double minvMinCut, double minvMaxCut, double isodetaCut, double isodrCut, double isoptCut, double isoCut){
  
  
  //

  double ptminCut = 0 ; 
  double ptpairCut = 0 ;///  double s4s9Cut, double s9s25Cut, double minvMinCut, double minvMaxCut, double isodetaCut, double isodrCut, double isoptCut, double isoCut; 
  double ptminCut_endcap_region1 = 0;
  double ptminCut_endcap_region2 = 0 ;
  double ptminCut_endcap_region3 = 0 ;
  double ptpairCut_endcap_region1 = 0 ;
  double ptpairCut_endcap_region2 = 0 ;
  double ptpairCut_endcap_region3 = 0 ;
  double ptpairCutMax_endcap_region3 = 99999 ;
  
  double eta_endcap_region1 = 0 ;
  double eta_endcap_region2 = 0 ;
  
  
  if( detector == reco::CaloID::DET_ECAL_ENDCAP){
    
    if(pizEta == PIZ){
      
      ptminCut_endcap_region1 = selePtGammaPi0EndCap_region1_; 
      ptminCut_endcap_region2 = selePtGammaPi0EndCap_region2_; 
      ptminCut_endcap_region3 = selePtGammaPi0EndCap_region3_; 

      ptpairCut_endcap_region1 = selePtPi0EndCap_region1_;
      ptpairCut_endcap_region2 = selePtPi0EndCap_region2_;
      ptpairCut_endcap_region3 = selePtPi0EndCap_region3_;
      ptpairCutMax_endcap_region3 = selePtPi0MaxEndCap_region3_;
      
      eta_endcap_region1 = region1_Pi0EndCap_; 
      eta_endcap_region2 = region2_Pi0EndCap_; 
      
    }else{
      
      ptminCut_endcap_region1 = selePtGammaEtaEndCap_region1_; 
      ptminCut_endcap_region2 = selePtGammaEtaEndCap_region2_; 
      ptminCut_endcap_region3 = selePtGammaEtaEndCap_region3_; 

      ptpairCut_endcap_region1 = selePtEtaEndCap_region1_;
      ptpairCut_endcap_region2 = selePtEtaEndCap_region2_;
      ptpairCut_endcap_region3 = selePtEtaEndCap_region3_;
      
      eta_endcap_region1 = region1_EtaEndCap_; 
      eta_endcap_region2 = region2_EtaEndCap_; 
      
    }
    
  }else{
    if(pizEta == PIZ){
      ptminCut = selePtGamma_; 
      ptpairCut = selePtPi0_; 
    }else{
      ptminCut = selePtGammaEta_; 
      ptpairCut = selePtEta_; 
    }
  }
  
  //cout<<"checkme: "<< detector <<" "<< pizEta <<" "<< ptminCut <<" "<< ptpairCut <<" "<< s4s9Cut<<" "<< s9s25Cut <<" "<< isoCut <<" "<< minvMinCut <<" "<< minvMaxCut <<" "<< isodetaCut <<" "<< isodrCut<<" "<< isoptCut<<endl; 
    
  const CaloSubdetectorGeometry *geometry_es; 
  CaloSubdetectorTopology *topology_es=0;
  if( detector == reco::CaloID::DET_ECAL_ENDCAP){
    edm::ESHandle<CaloGeometry> geoHandle;
    es.get<CaloGeometryRecord>().get(geoHandle);
    geometry_es = geoHandle->getSubdetectorGeometry(DetId::Ecal, EcalPreshower);
    topology_es  = new EcalPreshowerTopology(geoHandle);
  }
  
  esrechits_map.clear();
  
  
  if( detector == reco::CaloID::DET_ECAL_ENDCAP && storeRecHitES_){
    Handle<ESRecHitCollection> esRecHitsHandle;
    try{
      if(pizEta==1){
	evt.getByLabel(preshHitProducer_,esRecHitsHandle);
      }else{
	evt.getByLabel(preshHitProducerEta_,esRecHitsHandle);
      }
      
      ///const EcalRecHitCollection* hitCollection_es = esRecHitsHandle.product();
      EcalRecHitCollection::const_iterator iter;
      for (iter = esRecHitsHandle->begin(); iter != esRecHitsHandle->end(); iter++) {
	//Make the map of DetID, EcalRecHit pairs
	esrechits_map.insert(std::make_pair(iter->id(), *iter));   
      }
      
    }catch(std::exception& ex ){
      storeRecHitES_ = false; 
      cout<<"preshowerRecHit not working.."<<endl;
    }
    
  }
  

  vector<int> indClusPi0Candidates;  ///those clusters identified as pi0s
  if( detector == reco::CaloID::DET_ECAL_BARREL && rmPiz){
    for(int i=0 ; i<nClus ; i++){
      for(int j=i+1 ; j<nClus ; j++){
	
	double p0x = etClus[i] * cos(phiClus[i]);
	double p1x = etClus[j] * cos(phiClus[j]);
	double p0y = etClus[i] * sin(phiClus[i]);
	double p1y = etClus[j] * sin(phiClus[j]);
	double p0z = eClus[i] * cos(thetaClus[i]);
	double p1z = eClus[j] * cos(thetaClus[j]);
      
	double m_inv = sqrt ( (eClus[i] + eClus[j])*(eClus[i] + eClus[j]) - (p0x+p1x)*(p0x+p1x) - (p0y+p1y)*(p0y+p1y) - (p0z+p1z)*(p0z+p1z) );  
      
	if( m_inv > massLowPi0Cand_  &&   m_inv <  massHighPi0Cand_){

	  int indtmp[2] = {i,j};
	  for(int k=0; k<2; k++){
	    std::vector<int>::iterator it = find(indClusPi0Candidates.begin(),indClusPi0Candidates.end(),indtmp[k]);
	    if( it == indClusPi0Candidates.end()) indClusPi0Candidates.push_back(indtmp[k]);
	  }
	  
	}
	
      }
      
    }
  }
  

  used_strips.clear();
  ///we need to make sure for every endcap cluster we have the same preshower cluster
  std::map<int, vector<float > > infoES_saved; 
  vector<int> indClusSelected;
  

  ///after selection do NxN matching with MC
  vector<int> indPairCand;  ///pair of cluster passed the cuts    
  vector<float> mPairCand; /// mpair of of cluster passed the cuts 
  vector<float> mPairv1Cand; /// mpair of of cluster passed the cuts   
    
  vector<float> etaPairCand; 
  vector<float> phiPairCand; 
  vector<float> ptPairCand; 
  vector<float> isoPairCand; 


  

  
  
  for(int i=0 ; i<nClus ; i++){
    
    if( detector == reco::CaloID::DET_ECAL_BARREL && rmPiz){
      std::vector<int>::iterator it = find(indClusPi0Candidates.begin(),indClusPi0Candidates.end(),i);
      if( it != indClusPi0Candidates.end()) continue; 
    }
    
    if ( s4s9Clus[i] < s4s9Cut || s9s25Clus[i] < s9s25Cut ) continue; 
    
    for(int j=i+1 ; j<nClus ; j++){
      
      if( detector == reco::CaloID::DET_ECAL_BARREL && rmPiz){
	std::vector<int>::iterator it = find(indClusPi0Candidates.begin(),indClusPi0Candidates.end(),j);
	if( it != indClusPi0Candidates.end()) continue; 
      }
      
      if ( s4s9Clus[j] < s4s9Cut || s9s25Clus[j] < s9s25Cut ) continue; 
      
      double p0x = etClus[i] * cos(phiClus[i]);
      double p1x = etClus[j] * cos(phiClus[j]);
      double p0y = etClus[i] * sin(phiClus[i]);
      double p1y = etClus[j] * sin(phiClus[j]);
      double p0z = eClus[i] * cos(thetaClus[i]);
      double p1z = eClus[j] * cos(thetaClus[j]);
      ptpair = sqrt( (p0x+p1x)*(p0x+p1x) + (p0y+p1y)*(p0y+p1y));
      if (ptpair < ptpairCut )continue;
      TVector3 pairVect = TVector3((p0x+p1x), (p0y+p1y), (p0z+p1z));
      etapair = pairVect.Eta();
      phipair = pairVect.Phi();
      
      ptmin = etClus[i]< etClus[j] ?  etClus[i] : etClus[j] ; 
      if( detector == reco::CaloID::DET_ECAL_ENDCAP){
	if( fabs(etapair) < eta_endcap_region1){
	  if(ptmin < ptminCut_endcap_region1 || ptpair < ptpairCut_endcap_region1) continue; 
	}else if (fabs(etapair)< eta_endcap_region2){
	  if(ptmin < ptminCut_endcap_region2 || ptpair < ptpairCut_endcap_region2) continue; 
	}else{
	  if(ptmin < ptminCut_endcap_region3 || ptpair < ptpairCut_endcap_region3) continue; 
	  if(ptpair > ptpairCutMax_endcap_region3 ) continue; 
	}
      }else{
	
	if(ptmin < ptminCut || ptpair < ptpairCut) continue; 
	
      }
      

      mpair = sqrt ( (eClus[i] + eClus[j])*(eClus[i] + eClus[j]) - (p0x+p1x)*(p0x+p1x) - (p0y+p1y)*(p0y+p1y) - (p0z+p1z)*(p0z+p1z) );  

      
      if( mpair < minvMinCut || mpair > minvMaxCut ) continue; 
            
      if(debug_>=2) cout<<"selpi0pair1: "<< ptpair <<" "<< etClus[i]<<" "<< etClus[j]<<" "<< s4s9Clus[i]<<" "<< s4s9Clus[j]<<" "<<mpair<<endl;
      
      
      ///Isolation
      isolation = 0; 
      for(int k=0 ; k<nClus ; k++){
	if(k==i || k==j)continue;
	if( etClus[k] < isoptCut ) continue; 
	
	TVector3 clusVect = TVector3(etClus[k] *cos(phiClus[k]), etClus[k] * sin(phiClus[k]) , eClus[k] * cos(thetaClus[k]));
	double dretacl = fabs(etaClus[k] - pairVect.Eta());
	double drcl = clusVect.DeltaR(pairVect);
	
	if( drcl < isodrCut && dretacl < isodetaCut ){
	  isolation += etClus[k];
	}
      }
      isolation /= ptpair; 
      
      
      if( isolation > isoCut ) continue; 
      

      TVector3 vtmp0( xClus[i] - xOffPVwithBS, yClus[i] - yOffPVwithBS,  zClus[i] - zOffPVwithBS); 
      TVector3 vtmp1( xClus[j] - xOffPVwithBS, yClus[j] - yOffPVwithBS,  zClus[j] - zOffPVwithBS); 
      
      float cosd = getcosd( vtmp0.Eta(),vtmp0.Phi(),vtmp1.Eta(),vtmp1.Phi());
      mpairv1 = sqrt( 2 * eClus[i] * eClus[j] * (1- cosd));
      
      indPairCand.push_back(i * nClus + j);
      mPairCand.push_back(mpair);
      mPairv1Cand.push_back(mpairv1);
      ptPairCand.push_back(ptpair);
      etaPairCand.push_back(etapair);
      phiPairCand.push_back(phipair);
      isoPairCand.push_back(isolation);
                  
    }
  }
    
  
  vector<int> pairRecoMatchedtoGenPiz; 
  vector<int> pairRecoMatchedtoGen1stPhtInd; //the 1st photon int the pair matched to 0 /1 genpht
  vector<int> genpi0MatchedtoReco;  //for pi0/eta both
  
  float drMatch = 0.05; 
  if( detector == reco::CaloID::DET_ECAL_BARREL){
    drMatch  = 0.03; 
  }
  
  if( InputDataFormat_ <10){
    if(pizEta == PIZ){
      NtoNmatchingPi0(drMatch, indPairCand,pairRecoMatchedtoGenPiz,pairRecoMatchedtoGen1stPhtInd,genpi0MatchedtoReco);
    }else{
      NtoNmatchingEta(drMatch, indPairCand,pairRecoMatchedtoGenPiz,pairRecoMatchedtoGen1stPhtInd,genpi0MatchedtoReco);
    }
    
  }
  
  
  ///Now Fill what selected
  for(int k=0;  k< int(indPairCand.size()); k++){
    
    mpair = mPairCand[k];
    mpairv1 = mPairv1Cand[k];
    
    int i = indPairCand[k] / nClus; 
    int j = indPairCand[k] % nClus; 
    
    ptpair = ptPairCand[k];
    etapair = etaPairCand[k];
    phipair = phiPairCand[k];
    isolation = isoPairCand[k];
    ptmin = etClus[i] < etClus[j]? etClus[i] : etClus[j]; 

    s4s9min = s4s9Clus[i] < s4s9Clus[j] ? s4s9Clus[i] : s4s9Clus[j]; 
    s9s25min = s9s25Clus[i] < s9s25Clus[j] ? s9s25Clus[i] : s9s25Clus[j]; 
        

    if( InputDataFormat_ <10){
      for(int kk =0; kk<9; kk++) geninfo[kk]= -99;
      genMatched = 0; 
      int indgen = pairRecoMatchedtoGenPiz[k]; 
      if(indgen>=0){
	genMatched = 1; 
	int indga1 = pairRecoMatchedtoGen1stPhtInd[k]; 
	if(pizEta == PIZ){
	  geninfo[0] = ePhtGenpi0[indgen][indga1];
	  geninfo[1] = etaPhtGenpi0[indgen][indga1];
	  geninfo[2] = phiPhtGenpi0[indgen][indga1];
	  geninfo[3] = ePhtGenpi0[indgen][1-indga1];
	  geninfo[4] = etaPhtGenpi0[indgen][1-indga1];
	  geninfo[5] = phiPhtGenpi0[indgen][1-indga1];
	  geninfo[6] = vtxPhtGenpi0[indgen][0];
	  geninfo[7] = vtxPhtGenpi0[indgen][1];
	  geninfo[8] = vtxPhtGenpi0[indgen][2];
	  
	  convinfo[0] = isConvPhtGenpi0[indgen][indga1];
	  convinfo[1] = isConvPhtGenpi0[indgen][1-indga1];
	  
	}else{
	  geninfo[0] = ePhtGeneta[indgen][indga1];
	  geninfo[1] = etaPhtGeneta[indgen][indga1];
	  geninfo[2] = phiPhtGeneta[indgen][indga1];
	  geninfo[3] = ePhtGeneta[indgen][1-indga1];
	  geninfo[4] = etaPhtGeneta[indgen][1-indga1];
	  geninfo[5] = phiPhtGeneta[indgen][1-indga1];
	  geninfo[6] = vtxPhtGeneta[indgen][0];
	  geninfo[7] = vtxPhtGeneta[indgen][1];
	  geninfo[8] = vtxPhtGeneta[indgen][2];
	  
	  convinfo[0] = isConvPhtGeneta[indgen][indga1];
	  convinfo[1] = isConvPhtGeneta[indgen][1-indga1];
	  
	}
	
      }
    }
    

    int indtmp[2]={i,j};
        
    
    izXtalClus1 = 0; 
    izXtalClus2 = 0; 
    
    xClus1 = xClus[i];
    yClus1 = yClus[i];
    zClus1 = zClus[i];
    
    
    xClus2 = xClus[j];
    yClus2 = yClus[j];
    zClus2 = zClus[j];
    
    

    
    for(int jj =0; jj<2; jj++){
      int ind = indtmp[jj];
      
      if( jj==0) nxtClus1 = int(RecHitsCluster[ind].size()); 
      else nxtClus2 = int(RecHitsCluster[ind].size()); 
            
      
      for(unsigned int Rec=0;Rec<RecHitsCluster[ind].size();Rec++) {
	int ietathis;
	int iphithis; 
	
	if( detector == reco::CaloID::DET_ECAL_BARREL){
	  EBDetId det = (EBDetId)RecHitsCluster[ind][Rec].id(); 
	  ietathis = det.ieta();
	  iphithis = det.iphi();
	  
	  if(Rec==0){
	    // hh_mul_ietaSeedSel->Fill(ietathis);
	    // hh_mul_iphiSeedSel->Fill(iphithis);
	  }
	  
	}else{
	  EEDetId det = (EEDetId)RecHitsCluster[ind][Rec].id(); 
	  ietathis = det.ix();
	  iphithis = det.iy();
	  if(jj==0) 	    izXtalClus1 = det.zside();
	  else   izXtalClus2 = det.zside();
	}
	float lasercalib = laser->getLaserCorrection( RecHitsCluster[ind][Rec].id(), evt.time());
	if( jj==0){
	  eXtalClus1[Rec] = RecHitsCluster[ind][Rec].energy();
	  ietaXtalClus1[Rec] = ietathis; 
	  iphiXtalClus1[Rec] = iphithis; 
	  
	  laserCorrXtalClus1[Rec] = lasercalib; 
	  
	  fXtalClus1[Rec] = RecHitsCluster[ind][Rec].recoFlag(); 
	  tXtalClus1[Rec] = RecHitsCluster[ind][Rec].time();
	}else{
	  ietaXtalClus2[Rec] = ietathis; 
	  iphiXtalClus2[Rec] = iphithis; 
	  laserCorrXtalClus2[Rec] = lasercalib; 
	  fXtalClus2[Rec] = RecHitsCluster[ind][Rec].recoFlag(); 
	  eXtalClus2[Rec] = RecHitsCluster[ind][Rec].energy();
	  tXtalClus2[Rec] = RecHitsCluster[ind][Rec].time();
	}
      }
    }

    if(pizEta==2){ //save also additional 5x5 around the 3x3 alredy saved
	
      for(int jj =0; jj<2; jj++){
	int ind = indtmp[jj];
	  
	if( jj==0) nxt5x5Clus1 = int(RecHitsCluster5x5[ind].size()); 
	else nxt5x5Clus2 = int(RecHitsCluster5x5[ind].size()); 
	
	for(unsigned int Rec=0;Rec<RecHitsCluster5x5[ind].size();Rec++) {
	  int ietathis;
	  int iphithis; 
	  if( detector == reco::CaloID::DET_ECAL_BARREL){
	    EBDetId det = (EBDetId)RecHitsCluster5x5[ind][Rec].id(); 
	    ietathis = det.ieta();
	    iphithis = det.iphi();
	  }else{
	    EEDetId det = (EEDetId)RecHitsCluster5x5[ind][Rec].id(); 
	    ietathis = det.ix();
	    iphithis = det.iy();
	  }
	  	  
	  if( jj==0){
	    eXtal5x5Clus1[Rec] = RecHitsCluster5x5[ind][Rec].energy();
	    ietaXtal5x5Clus1[Rec] = ietathis; 
	    iphiXtal5x5Clus1[Rec] = iphithis; 
	    fXtal5x5Clus1[Rec] = RecHitsCluster5x5[ind][Rec].recoFlag(); 
	    tXtal5x5Clus1[Rec] = RecHitsCluster5x5[ind][Rec].time();
	  }else{
	    ietaXtal5x5Clus2[Rec] = ietathis; 
	    iphiXtal5x5Clus2[Rec] = iphithis; 
	    fXtal5x5Clus2[Rec] = RecHitsCluster5x5[ind][Rec].recoFlag(); 
	    eXtal5x5Clus2[Rec] = RecHitsCluster5x5[ind][Rec].energy();
	    tXtal5x5Clus2[Rec] = RecHitsCluster5x5[ind][Rec].time();
	  }
	}
      }
    }
      
    ///Now finally fill the trees
      
    if( detector == reco::CaloID::DET_ECAL_BARREL){
      if(pizEta ==1) mytree_pizeb->Fill();
      else  mytree_etaeb->Fill();
    }else{
      
      //adding preshower information as well
      if( storeRecHitES_){
	

	for(int jj =0; jj<2; jj++){
	  for(int n=0; n<8; n++) {
	    infoESX[jj][n] = 0; 
	    infoESY[jj][n] = 0; 
	  }
	  int ind = indtmp[jj];
	  std::vector<int>::iterator it = find(indClusSelected.begin(),indClusSelected.end(),ind);
	  if( it == indClusSelected.end()){
	    indClusSelected.push_back(ind);
	    
	    double e1=0;
	    double e2=0;
	    double deltaE=0;
	    double deltaE1=0;
	    double deltaE2=0;
	    
	    const GlobalPoint point(xClus[ind],yClus[ind],zClus[ind]);
	    
	    
	    DetId tmp1 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_es))->getClosestCellInPlane(point, 1);
	    DetId tmp2 = (dynamic_cast<const EcalPreshowerGeometry*>(geometry_es))->getClosestCellInPlane(point, 2);
	    ESDetId strip1 = (tmp1 == DetId(0)) ? ESDetId(0) : ESDetId(tmp1);
	    ESDetId strip2 = (tmp2 == DetId(0)) ? ESDetId(0) : ESDetId(tmp2);     
	      
	    std::set<DetId>  used_strips_before = used_strips; 
	    // Get ES clusters (found by the PreshSeeded algorithm) associated with a given EE cluster.     
	    /// Energy weight cluster 1 
	    double wtposEs[2][3]={{0}};
	    double mposEs[2][3]={{0}}; //energy highest
	    double max_es[2]={0};
	      
	    for (int i2=0; i2<preshNclust_; i2++) {
	      reco::PreshowerCluster cl1 = presh_algo->makeOneCluster(strip1,&used_strips,&esrechits_map,geometry_es,topology_es);  
	      if ( cl1.energy() > preshClustECut) {
		e1 += cl1.energy();   
		
		wtposEs[0][0] += cl1.energy() * cl1.position().X();
		wtposEs[0][1] += cl1.energy() * cl1.position().Y();
		wtposEs[0][2] += cl1.energy() * cl1.position().Z();
		  
		if( cl1.energy() > max_es[0]){
		  max_es[0] = cl1.energy() ; 
		  mposEs[0][0] = cl1.position().X();
		  mposEs[0][1] = cl1.position().Y();
		  mposEs[0][2] = cl1.position().Z();
		}
		  
	      }
	      reco::PreshowerCluster cl2 = presh_algo->makeOneCluster(strip2,&used_strips,&esrechits_map,geometry_es,topology_es); 
	      if ( cl2.energy() > preshClustECut) {
		e2 += cl2.energy();       
		wtposEs[1][0] +=  cl2.energy() * cl2.position().X();
		wtposEs[1][1] +=  cl2.energy() * cl2.position().Y();
		wtposEs[1][2] +=  cl2.energy() * cl2.position().Z();
		if( cl2.energy() > max_es[1]){
		  max_es[1] = cl2.energy(); 
		  mposEs[1][0] = cl2.position().X();
		  mposEs[1][1] = cl2.position().Y();
		  mposEs[1][2] = cl2.position().Z();
		}
	      }
	    } // end of cycle over ES clusters
	      
	    if(e1+e2 > 1.0E-10) {
		
	      if( e1 > 0){
		for(int n=0; n<3; n++){
		  wtposEs[0][n] /= e1; 
		}
	      }
	      if( e2 > 0){
		for(int n=0; n<3; n++){
		  wtposEs[1][n] /= e2; 
		}
	      }
	      // GeV to #MIPs
	      e1 = e1 / mip_;
	      e2 = e2 / mip_;
	      deltaE = gamma_*(calib_planeX_*e1+calib_planeY_*e2);
	      deltaE1 = gamma_*(calib_planeX_*e1);
	      deltaE2 = gamma_*(calib_planeY_*e2);
	      
	      max_es[0] *= gamma_*(calib_planeX_)/ mip_; 
	      max_es[1] *= gamma_*(calib_planeY_)/ mip_; 
	      
	      if(debug_>=2){
		cout<<"e1e1max: "<< deltaE1 <<" "<< max_es[0]<<" "<< deltaE2 <<" "<<max_es[1]<<endl;
	      }
	      
	      ///x plane
	      infoESX[jj][0] = deltaE1;
	      for(int n=0; n<3; n++){
		infoESX[jj][n+1] = wtposEs[0][n]; 
	      }
	      infoESX[jj][4] = max_es[0];
	      for(int n=0; n<3; n++){
		infoESX[jj][n+5] = mposEs[0][n]; 
	      }
	      //y plane
	      infoESY[jj][0] = deltaE2;
	      for(int n=0; n<3; n++){
		infoESY[jj][n+1] = wtposEs[1][n]; 
	      }
	      infoESY[jj][4] = max_es[1];
	      for(int n=0; n<3; n++){
		infoESY[jj][n+5] = mposEs[1][n]; 
	      }
		
	    } 
	    //now push_back to the map of ES information
	    vector<float> infoES_v; 
	    for(int n=0; n<8; n++) {
	      infoES_v.push_back(infoESX[jj][n]) ;
	    }
	    for(int n=0; n<8; n++) {
	      infoES_v.push_back(infoESY[jj][n]) ;
	    }
	    infoES_saved.insert( make_pair(ind, infoES_v) ); 
	    	  
	  
	  }else{
	    std::map<int, vector<float> >::iterator itmap = infoES_saved.find(ind);
	    ///double es = (itmap->second).energy();
	    if( itmap != infoES_saved.end()){
	      for(int n=0;n<8; n++) infoESX[jj][n] = (itmap->second)[n];
	      for(int n=0;n<8; n++) infoESY[jj][n] = (itmap->second)[n+8];
	    }
	  }
	  
	}
	
      }
      
      if(pizEta ==1) mytree_pizee->Fill();
      else  mytree_etaee->Fill();
    } /// endcap 
    
    
    
  } ///end of filling tree
  
  
  
  
  
  
  delete topology_es;
  
  
  
}

void RecoAnalyzer::NtoNmatchingPi0(float drMatch, vector<int> pairReco ,vector<int> &pairRecoMatchedtoGen, vector<int> &pairRecoMatchedtoGen1stPhtInd,vector<int> &genpi0MatchedtoReco){
  
  float zpht;
  float rpht; 
  float xpht; 
  float ypht; 
  
  //a list of pair reco ( direction w.r.t (0,0,0))
  
  ///matched to a list of gen 
  
  
  ////based on the smallest dr of sum of dr_photon1+ dr_photon2 
  vector< std::pair<double,int> > dr_jk; 
  int n =0; 
  
  ////try to use x,y,Z of cluster w.r.t to gen-Level pi0/eta
  ///bool useGenVtxforRecoGenMatch = true; 
  int indga1 = 0;
  
  vector< std::pair<int,int> > indga1_jk; 
  

  for(int j=0; j< int(pairReco.size()); j++){
    
    int j1 = pairReco[j]/nClus; 
    int j2 = pairReco[j]%nClus; 

    int jtmp[2] = {j1,j2};
    
    ////eta ,phi from (0,0,0);
    float eta_2pht[2] ={etaClus[j1],etaClus[j2]};
    float phi_2pht[2] ={phiClus[j1],phiClus[j2]};
    for(int tt=0;tt<2; tt++){
      TVector3 vv(xClus[jtmp[tt]],yClus[jtmp[tt]],zClus[jtmp[tt]]);
      eta_2pht[tt] = vv.Eta();
      phi_2pht[tt] = vv.Phi();
    }
    
    
    for(int k=0; k< nGenpi0; k++){
      
      float dr_twopht[2]={1,1};
      float drsum = 0; 
      float allDR[4]; 
    
      int kk = 0; 
      for( int k1 =0; k1<2; k1++){
	float eta = eta_2pht[k1];
	float phi = phi_2pht[k1];

	for( int k2 = 0; k2<2; k2++){

	  zpht = vtxPhtGenpi0[k][2]; 
	  rpht = sqrt( vtxPhtGenpi0[k][0] * vtxPhtGenpi0[k][0] + vtxPhtGenpi0[k][1] * vtxPhtGenpi0[k][1]); 
	  xpht = vtxPhtGenpi0[k][0]; 
	  ypht = vtxPhtGenpi0[k][1]; 

	  float etag = etaPhtGenpi0[k][k2]; 
	  float phig = phiPhtGenpi0[k][k2]; 

	  etag = ecalEta(etaPhtGenpi0[k][k2],zpht,rpht);
	  phig = ecalPhi(phiPhtGenpi0[k][k2],xpht,ypht);
	  
	  
	  
	  allDR[kk] = GetDeltaR(eta,etag,phi,phig);
	  kk++;
	  
	}
      }
      
      float drsum1 = allDR[0] + allDR[3];
      float drsum2 = allDR[1] + allDR[2];
      drsum = drsum1 < drsum2? drsum1: drsum2; 
      
      if(drsum1<drsum2){
	dr_twopht[0] = allDR[0];
	dr_twopht[1] = allDR[3];
	indga1 = 0;
      }else{
	dr_twopht[0] = allDR[1];
	dr_twopht[1] = allDR[2];
	indga1 = 1;
      }
      
      if( dr_twopht[0]< drMatch && dr_twopht[1] < drMatch){
	dr_jk.push_back( make_pair ( drsum, n));
	
	indga1_jk.push_back( make_pair ( indga1, n)); 
      }
      
      n++; 
      
    }
  }
  
  
  sort(dr_jk.begin(),dr_jk.end(),sort_pred);
  
  vector<int> nn_matched;
  vector<int> j_used; 
  vector<int> k_used; 
  vector<int>::iterator it;
  vector<int>::iterator it2;
  
  int nn = int( dr_jk.size());

  for(int n=0; n< nn; n++){

    int nn0 = dr_jk[n].second;
    int j = nn0 /  nGenpi0; 
    int k = nn0 %  nGenpi0; 
    
    it = find(j_used.begin(),j_used.end(),j);
    it2 = find(k_used.begin(),k_used.end(),k);
    
    if( it == j_used.end() && it2 == k_used.end()){
      j_used.push_back(j);
      k_used.push_back(k);
      nn_matched.push_back(nn0); /// this matters a lot if to use more than once when matching...
      
      ////cout<<"nn_matched: "<<entry<<" "<< nn0 <<" "<<j<<" "<<k <<endl; 
      
      // }else{
      //cout<<"entry: "<<j<<" "<<k<<" "<< "used already.." <<entry<<endl; 
      // exit(1);
      ////entry_drjk: 2364 0.00712235 0.0422643 0 1 1 eta/phi: -0.135937 -0.233336 -0.983332 -1.16227 gen: -0.217764 -0.146435 -1.13747 -0.990584
      /////entry_drjk: 2364 0.00712235 0.00420043 1 1 11 eta/phi: -0.135937 -0.205423 -0.983332 -1.1358 gen: -0.217764 -0.146435 -1.13747 -0.990584
      
    }
    
  }
  
  
  ///nn_matched are the final matched pair reco <==> gen pi0 
  for(int j=0; j< int(pairReco.size()); j++){
    int matched = -1;  //if matched, be the index in the array of genpi0
    int nn0 = -1; 
    for(int n=0; n< int( nn_matched.size()); n++){
      int j1 = nn_matched[n]/ nGenpi0; 
      int k = nn_matched[n] % nGenpi0; 
      if( j1 == j){
	matched = k; 
	nn0 = nn_matched[n];
	break; 
      }
    }
    
    ////cout<<"pairREco Mathced: "<< j<<" "<< matched <<endl; 

    pairRecoMatchedtoGen.push_back(matched); 
    indga1 = 0; 
    for(int n =0; n< int( indga1_jk.size()); n++){
      if( indga1_jk[n].second == nn0){
	indga1 = indga1_jk[n].first; 
	n = int( indga1_jk.size()+1); 
      }
    }
    pairRecoMatchedtoGen1stPhtInd.push_back(indga1); 
    
  }
  
  
//   for(int k=0; k< nGenpi0; k++){
//     int matched = -1;  //if matched, be the index in the array of pairReco
    
//     for(int n=0; n< int( nn_matched.size()); n++){
//       int j = nn_matched[n]/ nGenpi0; 
//       int k1 = nn_matched[n] % nGenpi0; 
//       if( k1 == k){
// 	matched = j; 
// 	break; 
//       }
//     }
//     genpi0MatchedtoReco.push_back(matched);
//   }
  
  
}

void RecoAnalyzer::NtoNmatchingEta(float drMatch, vector<int> pairReco ,vector<int> &pairRecoMatchedtoGen, vector<int> &pairRecoMatchedtoGen1stPhtInd, vector<int> &genetaMatchedtoReco){
  
  
  float zpht;
  float rpht; 
  float xpht; 
  float ypht; 
  
  //a list of pair reco ( direction w.r.t (0,0,0))
  
  ///matched to a list of gen 
  int indga1 = 0;
  vector< std::pair<int,int> > indga1_jk; 
  
  
  
  ////based on the smallest dr of sum of dr_photon1+ dr_photon2 
  vector< std::pair<double,int> > dr_jk; 
  int n =0; 
  
  ////try to use x,y,Z of cluster w.r.t to gen-Level eta/eta
  ///bool useGenVtxforRecoGenMatch = true; 
  

  for(int j=0; j< int(pairReco.size()); j++){
    
    int j1 = pairReco[j]/nClus ; 
    int j2 = pairReco[j]%nClus; 

    int jtmp[2] = {j1,j2};
    
    ////eta ,phi from (0,0,0);
    float eta_2pht[2] ={etaClus[j1],etaClus[j2]};
    float phi_2pht[2] ={phiClus[j1],phiClus[j2]};
    for(int tt=0;tt<2; tt++){
      TVector3 vv(xClus[jtmp[tt]],yClus[jtmp[tt]],zClus[jtmp[tt]]);
      eta_2pht[tt] = vv.Eta();
      phi_2pht[tt] = vv.Phi();
    }
    

    for(int k=0; k< nGeneta; k++){
      
      float dr_twopht[2]={1,1};
      float drsum = 0; 
      float allDR[4]; 
    
      int kk = 0; 
      for( int k1 =0; k1<2; k1++){
	float eta = eta_2pht[k1];
	float phi = phi_2pht[k1];

	for( int k2 = 0; k2<2; k2++){
	  
	  zpht = vtxPhtGeneta[k][2]; 
	  rpht = sqrt( vtxPhtGeneta[k][0] * vtxPhtGeneta[k][0] + vtxPhtGeneta[k][1] * vtxPhtGeneta[k][1]); 
	  xpht = vtxPhtGeneta[k][0]; 
	  ypht = vtxPhtGeneta[k][1]; 
	  
	  float etag = etaPhtGeneta[k][k2]; 
	  float phig = phiPhtGeneta[k][k2]; 
	  etag = ecalEta(etaPhtGeneta[k][k2],zpht,rpht);
	  phig = ecalPhi(phiPhtGeneta[k][k2],xpht,ypht);
	  
	  allDR[kk] = GetDeltaR(eta,etag,phi,phig);
	  kk++;
	  
	}
      }
      
      float drsum1 = allDR[0] + allDR[3];
      float drsum2 = allDR[1] + allDR[2];
      drsum = drsum1 < drsum2? drsum1: drsum2; 
      
      if(drsum1<drsum2){
	dr_twopht[0] = allDR[0];
	dr_twopht[1] = allDR[3];
      }else{
	dr_twopht[0] = allDR[1];
	dr_twopht[1] = allDR[2];
      }
      
      if( dr_twopht[0]< drMatch && dr_twopht[1] <drMatch ){
	dr_jk.push_back( make_pair ( drsum, n));
	indga1_jk.push_back( make_pair ( indga1, n)); 
	
	
      }
      
      n++; 
      
    }
  }
  
  
  sort(dr_jk.begin(),dr_jk.end(),sort_pred);
  
  vector<int> nn_matched;
  vector<int> j_used; 
  vector<int> k_used; 
  vector<int>::iterator it;
  vector<int>::iterator it2;
  
  int nn = int( dr_jk.size());

  for(int n=0; n< nn; n++){

    int nn0 = dr_jk[n].second;
    int j = nn0 /  nGeneta; 
    int k = nn0 %  nGeneta; 
    
    it = find(j_used.begin(),j_used.end(),j);
    it2 = find(k_used.begin(),k_used.end(),k);
    
    if( it == j_used.end() && it2 == k_used.end()){
      j_used.push_back(j);
      k_used.push_back(k);
      nn_matched.push_back(nn0); /// this matters a lot if to use more than once when matching...
      
      ////cout<<"nn_matched: "<<entry<<" "<< nn0 <<" "<<j<<" "<<k <<endl; 

      // }else{
      //cout<<"entry: "<<j<<" "<<k<<" "<< "used already.." <<entry<<endl; 
      // exit(1);
      ////entry_drjk: 2364 0.00712235 0.0422643 0 1 1 eta/phi: -0.135937 -0.233336 -0.983332 -1.16227 gen: -0.217764 -0.146435 -1.13747 -0.990584
      /////entry_drjk: 2364 0.00712235 0.00420043 1 1 11 eta/phi: -0.135937 -0.205423 -0.983332 -1.1358 gen: -0.217764 -0.146435 -1.13747 -0.990584

    }
    
  }
  
  
  ///nn_matched are the final matched pair reco <==> gen eta 
  for(int j=0; j< int(pairReco.size()); j++){
    int matched = -1;  //if matched, be the index in the array of geneta

    int nn0 = -1; 
    
    for(int n=0; n< int( nn_matched.size()); n++){
      int j1 = nn_matched[n]/ nGeneta; 
      int k = nn_matched[n] % nGeneta; 
      if( j1 == j){
	matched = k; 
	
	nn0 = nn_matched[n]; 
	
	break; 
      }
    }

    ////cout<<"pairREco Mathced: "<< j<<" "<< matched <<endl; 
    
    pairRecoMatchedtoGen.push_back(matched); 
    indga1 = 0; 
    for(int n =0; n< int( indga1_jk.size()); n++){
      if( indga1_jk[n].second == nn0){
	indga1 = indga1_jk[n].first; 
	n = int( indga1_jk.size()+1);  ///break; 
      }
    }
    pairRecoMatchedtoGen1stPhtInd.push_back(indga1); 
    
    
  }
  
  
//   for(int k=0; k< nGeneta; k++){
//     int matched = -1;  //if matched, be the index in the array of pairReco
    
//     for(int n=0; n< int( nn_matched.size()); n++){
//       int j = nn_matched[n]/ nGeneta; 
//       int k1 = nn_matched[n] % nGeneta; 
//       if( k1 == k){
// 	matched = j; 
// 	break; 
//       }
//     }
//     genetaMatchedtoReco.push_back(matched);
//   }
  
  
}



//define this as a plug-in
DEFINE_FWK_MODULE(RecoAnalyzer);
