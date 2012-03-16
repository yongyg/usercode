/* // system include files */
/* #include <memory> */

/* // user include files */
/* #include "FWCore/Framework/interface/Frameworkfwd.h" */
/* //#include "FWCore/Framework/interface/EDProducer.h" */
/* #include "FWCore/Framework/interface/EDAnalyzer.h" */

/* #include "FWCore/Framework/interface/Event.h" */
/* #include "FWCore/Framework/interface/MakerMacros.h" */

/* #include "FWCore/ParameterSet/interface/ParameterSet.h" */
/* //#include "FWCore/ParameterSet/interface/InputTag.h" */
/* #include "FWCore/Utilities/interface/InputTag.h" */
/* #include "FWCore/Framework/interface/EventSetup.h" */
/* #include "FWCore/Utilities/interface/Exception.h" */
/* #include "FWCore/Framework/interface/ESHandle.h" */


/* #include "DataFormats/VertexReco/interface/Vertex.h" */
/* #include "DataFormats/VertexReco/interface/VertexFwd.h" */

/* #include "DataFormats/EcalRecHit/interface/EcalRecHit.h" */
/* #include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h" */
/* #include "DataFormats/EcalRecHit/interface/EcalUncalibratedRecHit.h" */

/* //#include "TrackingTools/TrackAssociator/interface/TimerStack.h" */
/* #include "Utilities/Timing/interface/TimerStack.h" */

/* #include "HLTrigger/HLTcore/interface/HLTFilter.h" */
/* // Geometry */
/* #include "Geometry/CaloEventSetup/interface/CaloTopologyRecord.h" */
/* #include "Geometry/CaloTopology/interface/CaloTopology.h" */
/* #include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h" */
/* #include "Geometry/Records/interface/IdealGeometryRecord.h" */
/* #include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h" */
/* #include "Geometry/CaloGeometry/interface/CaloCellGeometry.h" */
/* #include "Geometry/CaloGeometry/interface/CaloGeometry.h" */
/* #include "Geometry/CaloTopology/interface/EcalEndcapTopology.h" */
/* #include "Geometry/CaloTopology/interface/EcalBarrelTopology.h" */


/* #include <vector> */
/* #include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h" */


/* #include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h" */
/* #include "DataFormats/Math/interface/Point3D.h" */

/* // ES stuff */
/* #include "DataFormats/EgammaReco/interface/PreshowerCluster.h" */
/* #include "RecoEcal/EgammaClusterAlgos/interface/PreshowerClusterAlgo.h" */



/* #include "DataFormats/L1Trigger/interface/L1EmParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h" */

/* //L1 Trigger */
/* #include "DataFormats/L1Trigger/interface/L1EmParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1JetParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1MuonParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1EtMissParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h" */
/* #include "DataFormats/L1Trigger/interface/L1ParticleMap.h" */
/* #include "L1Trigger/L1ExtraFromDigis/interface/L1ExtraParticleMapProd.h" */
/* #include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h" */

/* #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" */
/* #include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" */
/* #include "DataFormats/HepMCCandidate/interface/GenParticle.h" */


/* #include "DataFormats/BeamSpot/interface/BeamSpot.h" */

/* #include "DataFormats/Common/interface/Ref.h" */
/* #include "DataFormats/Common/interface/RefVector.h" */

/* #include "Geometry/TrackerGeometryBuilder/interface/TrackerGeometry.h" */
/* #include "Geometry/CommonDetUnit/interface/GeomDet.h" */
/* #include "Geometry/Records/interface/TrackerDigiGeometryRecord.h" */


/* #include "SimDataFormats/Track/interface/SimTrack.h" */
/* #include "SimDataFormats/Track/interface/SimTrackContainer.h" */
/* #include "SimDataFormats/Vertex/interface/SimVertexContainer.h" */


/* #include "DataFormats/EcalDigi/interface/EcalDigiCollections.h" */

/* //#include "DataFormats/EcalDigi/interface/EcalDigiv1Collections.h" */

/* #include <set> */

/* #include "CondFormats/EcalObjects/interface/EcalGainRatios.h" */
/* #include "CondFormats/EcalObjects/interface/EcalPedestals.h" */
/* #include "CondFormats/DataRecord/interface/EcalPedestalsRcd.h" */
/* #include "CondFormats/DataRecord/interface/EcalGainRatiosRcd.h" */
/* #include "CondFormats/EcalObjects/interface/EcalMGPAGainRatio.h" */


/* #include "DataFormats/EcalDigi/interface/EEDataFrame.h" */
/* #include "DataFormats/EcalDigi/interface/EBDataFrame.h" */
/* #include "DataFormats/EcalDigi/interface/ESDataFrame.h" */
/* #include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h" */
/* #include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h" */
/* #include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h" */
/* #include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h" */
/* #include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h" */
/* #include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h" */

/* #include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h" */
/* #include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h" */
/* #include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h" */
/* #include "CondFormats/EcalObjects/interface/EcalChannelStatus.h" */

/* #include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h" */


/* // class declaration */
/* // */
/* #include "DataFormats/Common/interface/Handle.h" */
/* #include "FWCore/Utilities/interface/RegexMatch.h" */

/* //ROOT  */
/* #include "TROOT.h" */
/* #include "TFile.h" */
/* #include "TLorentzVector.h" */
/* #include "TTree.h" */
/* #include "TBranch.h" */
/* #include "TH1.h" */
/* #include "TH2.h" */
/* #include "TVector3.h" */

/* #include <vector> */
/* #include <boost/foreach.hpp> */

/* #include <ext/algorithm> */
/* #include "CLHEP/Matrix/Vector.h" */
/* #include "CLHEP/Matrix/SymMatrix.h" */

/* #include "SimCalorimetry/EcalSimAlgos/interface/EcalSimParameterMap.h" */
/* #include "SimCalorimetry/CaloSimAlgos/interface/CaloHitResponse.h" */
/* //#include "SimCalorimetry/EcalSimAlgos/interface/EcalShape.h" */
/* #include "SimDataFormats/CrossingFrame/interface/MixCollection.h" */
/* #include "SimDataFormats/CaloHit/interface/PCaloHit.h" */

/* #include "CalibFormats/CaloObjects/interface/CaloSamples.h" */
/* #include "DataFormats/Provenance/interface/ParameterSetID.h" */
/* #include "FWCore/Framework/interface/ESWatcher.h" */
/* #include "HLTrigger/HLTcore/interface/HLTFilter.h" */

/* #include "DataFormats/Common/interface/TriggerResults.h" */

/* ///#include "FWCore/Framework/interface/TriggerNames.h" */
/* #include "FWCore/Common/interface/TriggerNames.h" */

/* #include "L1Trigger/GlobalTriggerAnalyzer/interface/L1GtUtils.h" */


/* #include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" */



/* // forward declarations */
/* namespace edm { */
/*   class TriggerResults; */
/* } */

/* class AlCaRecoTriggerBitsRcd; */


/* typedef std::map<DetId, EcalRecHit> RecHitsMap; */

/* // Less than operator for sorting EcalRecHits according to energy. */
/* class ecalRecHitLess : public std::binary_function<EcalRecHit, EcalRecHit, bool>  */
/* { */
/*  public: */
/*   bool operator()(EcalRecHit x, EcalRecHit y)  */
/*     {  */
/*       return (x.energy() > y.energy());  */
/*     } */
/* }; */

/* struct VertexHit { */
/*   float z; */
/*   float r; */
/*   float w; */
/* }; */




/* using namespace std; */


/* class Pi0Analyzer : public edm::EDAnalyzer { */
/*    public: */
/*       explicit Pi0Analyzer(const edm::ParameterSet&); */
/*       ~Pi0Analyzer(); */
      
      
/*       int filter(const edm::Event&, const edm::EventSetup&); */
      
/*       int convertSmToFedNumbBarrel(int, int);  */
/*       void convxtalid(int & , int &); */
/*       int diff_neta_s(int,int); */
/*       int diff_nphi_s(int,int); */
      
/*       void makeClusterES(double x, double y, double z,const CaloSubdetectorGeometry*&iSubGeom,CaloSubdetectorTopology*& topology_p); */
      
/*       std::vector<int> ListOfFEDS(double etaLow, double etaHigh, double phiLow, */
/*                                     double phiHigh, double etamargin, double phimargin); */
      
      
      
/*       double DeltaPhi(double v1, double v2);  */
/*       double GetDeltaR(double eta1, double eta2, double phi1, double phi2);  */
/*       double getcosd(double eta1, double phi1, double eta2, double phi2);  */
/*       int getMotherIndex(int ); */
/*       int indexofParticle(float px, float pz, int status); */
/*       void findgenpi0eta(); */
/*       //      virtual bool filter(edm::Event &, const edm::EventSetup&); */
/*       void matchGenpi0Pht(float etapair, float phipair,float eta1,float phi1,float eta2,float phi2,float res[]);  */
/*       void matchGenetaPht(float etapair, float phipair,float eta1,float phi1,float eta2,float phi2,float res[]);  */
      
/*       void calcPairPhoton(float en[],float eta[],float phi[],float res[]);  */
      
/*       double ecalEta(double EtaParticle ,double Zvertex, double RhoVertex);  */
/*       double ecalPhi(double phi,double x0,double y0);  */
      
/*       void NtoNmatchingEta(float drMatch, vector<int> pairReco ,vector<int> &pairRecoMatchedtoGen, vector<int> &pairRecoMatchedtoGen1stPhtInd, vector<int> &genetaMatchedtoReco); */
/*       void NtoNmatchingPi0(float drMatch, vector<int> pairReco ,vector<int> &pairRecoMatchedtoGen, vector<int> &pairRecoMatchedtoGen1stPhtInd, vector<int> &genetaMatchedtoReco); */
      
/*       void doSelectionAndFillTree(const edm::Event &evt, const edm::EventSetup &es, const reco::CaloID::Detectors detector, bool rmPiz, int pizEta,double s4s9Cut, double s9s25Cut, double minvMinCut, double minvMaxCut, */
/* 				  double isodetaCut, double isodrCut, double isoptCut, double isoCut );  */
      
      
/*  private: */
/*       virtual void beginJob() ; */
/*       virtual void analyze(const edm::Event&, const edm::EventSetup&); */
      
/*       virtual void endJob() ; */

//      static const int PIZ = 1; 
//      static const int ETA = 2; 
      
/*       L1GtUtils m_l1GtUtils; */
      
/*       typedef std::vector<std::string> vstring; */
/*       bool checkStatusOfEcalRecHit(const EcalChannelStatus &channelStatus, const EcalRecHit &rh); */
      

      
      std::vector<std::string> alcaL1trigNames_;
      
	
      vector<string> listOfAlCaL1Trigs;
      
/*       void makeNxNClusters(const edm::Event &evt, const edm::EventSetup &es, */
/* 			   const EcalRecHitCollection *hits, const reco::CaloID::Detectors detector); */


/*       /// initialize the trigger conditions (call this if the trigger paths have changed) */
/*       //void init(const edm::TriggerResults & results, const edm::EventSetup &iSetup); */
/*       int init(const edm::TriggerResults & results, const edm::EventSetup &iSetup); */
      
/*       edm::InputTag  beamSpotInputTag_;  */
      

/*       edm::InputTag m_l1GtTmLInputTag;  */
      
/*       edm::InputTag inputTag_;  */
/*       edm::TriggerNames triggerNames_; */
      /// false = and-mode (all requested triggers), true = or-mode (at least one)
      bool andOr_;
      /// throw on any requested trigger being unknown
      bool throw_;
      /// not empty => use read paths from AlCaRecoTriggerBitsRcd via this key
      std::string eventSetupPathsKey_;
/*       /// Watcher to be created and used if 'eventSetupPathsKey_' non empty: */
/*       edm::ESWatcher<AlCaRecoTriggerBitsRcd> *watchAlCaRecoTriggerBitsRcd_; */
      /// input patterns that will be expanded into trigger names
      std::vector<std::string>  HLTPatterns_;
      /// for each trigger/pattern define a prescale.
/*       /// Only a subset of the events per trigger will be selected */
/*       std::vector<uint32_t> HLTPrescales_; */
/*       /// This variable stores the prescales after pattern matching and expansion */
/*       std::vector<uint32_t> HLTPrescalesExpanded_;  */
/*       /// Scalers for the prescaling, L1 style. */
/*       std::vector<uint32_t> HLTPrescalesScalers;  */
/*       /// You can define a prescale after the filter logic has been applied */
      /// this is different from the single prescales. It can be applied on
      /// top of the single trigger prescales
/*       uint32_t HLTOverallPrescale_ ; */
/*       uint32_t HLTOverallPrescalesScaler_; */
      /// list of required HLT triggers by HLT name
      std::vector<std::string> HLTPathsByName_;
      /// list of required HLT triggers by HLT index
      std::vector<unsigned int> HLTPathsByIndex_;
      
      

      

      double MinAmpl_;
      bool dyn_pedestal;
      double fAlpha_;//parameter of the shape
      double fBeta_;//parameter of the shape
      double fAmp_max_;// peak amplitude 
      double fTim_max_ ;// time of the peak (in 25ns units)
      double fPed_max_;// pedestal value
      double alfabeta_; 
      
  
      int fNb_iter_;
      int  fNum_samp_bef_max_ ;
      int fNum_samp_after_max_;
  
      float fSigma_ped; 
      double un_sur_sigma;
      //temporary solution for different alpha and beta
      //  float alpha_table_[36][1701];
      // float beta_table_[36][1701];
      
      bool doFit_;
      
      double alphaEB_; 
      double betaEB_; 
      double alphaEE_; 
      double betaEE_; 

      
/*       const EcalSimParameterMap * theParameterMap; */
/*       const CaloVShape * theEcalShape; */

/*       CaloHitResponse * theEcalResponse; */

      bool applyConstantTerm_;
      double rmsConstantTerm_;
      
      std::string srcVertex_; 
      
      std::string hitsProducer_;
      
/*       edm::ESHandle<EcalIntercalibConstants> ical; */

/*       const EcalIntercalibConstants *icalp;  */
/*       //const EcalIntercalibConstantMap &iicalMap;  */
/*       const EcalADCToGeVConstant *agcp;  */
      
/*       const EBDigiCollection* ebDigis ; */
/*       const EEDigiCollection* eeDigis ; */
      
/*       //      const EBDigiv1Collection* ebDigisv1 ; */
/*       // const EEDigiv1Collection* eeDigisv1 ; */
      

      
/*       edm::ESHandle<EcalTimeCalibConstants> itime; */
/*       edm::ESHandle<EcalADCToGeVConstant> agc; */
/*       edm::ESHandle<EcalChannelStatus> chStatus; */
/*       std::vector<int> v_chstatus_; */
/*       edm::ESHandle<EcalLaserDbService> laser; */
      

/*       EcalElectronicsMapping* TheMapping; */
/*       edm::InputTag srcPixels_;  */
      
/*       edm::InputTag  srcTowers_;  */
      

/*       edm::InputTag barrelHits_; */
/*       edm::InputTag endcapHits_; */
      
      
/*       edm::InputTag barrelHitsEta_; */
/*       edm::InputTag endcapHitsEta_; */
      
      

      int nhit_simeb[180][360];
      float en_simeb[180][360][20];
      float time_simeb[180][360][20];
      


      std::string pi0BarrelHits_;
      std::string pi0EndcapHits_;
      std::string etaBarrelHits_;
      std::string etaEndcapHits_;
      
      std::string pi0ESHits_;
      std::string etaESHits_;

      ///interal use
      std::string BarrelHits_;
      std::string EndcapHits_;
      std::string ESHits_;
      
/*       struct LessById { */
/* 	bool operator()(const SimTrack &tk1, const SimTrack &tk2) const { return tk1.trackId() < tk2.trackId(); } */
/* 	bool operator()(const SimTrack &tk1, unsigned int    id ) const { return tk1.trackId() < id;            } */
/* 	bool operator()(unsigned int     id, const SimTrack &tk2) const { return id            < tk2.trackId(); } */
/*       }; */
      
      
      

      ////gen-level pi0/eta
      vector<int> indpi0Gen;
      vector<int> indetaGen;
      vector<int> indpht1pi0Gen;
      vector<int> indpht2pi0Gen;
      vector<int> indpht1etaGen;
      vector<int> indpht2etaGen;
      
      
      int bunchX; 
      int orbitNumber; 
      int evtTime; 
      
      

      int InputDataFormat_; 
      int procID; 
      float ptHAT; 
      

      int maxNumberofSeeds_; 
      int maxNumberofClusters_; 
      
      int gammaCandEtaSize_;
      int gammaCandPhiSize_;
      
      double clusSeedThr_;
      double clusSeedThrEndCap_;
      
      int clusEtaSize_;
      int clusPhiSize_;
      double seleXtalMinEnergy_;
      double seleXtalMinEnergyEndCap_;
      int seleNRHMax_;

      //// for pi0->gg barrel 
      bool doSelForPi0Barrel_; 
      double selePtGamma_;
      double selePtPi0_;
      double seleMinvMaxPi0_;
      double seleMinvMinPi0_;
      double seleS4S9Gamma_;
      double selePi0BeltDR_;
      double selePi0BeltDeta_;
      double selePi0Iso_;
      double ptMinForIsolation_; 
      bool storeIsoClusRecHitPi0EB_; 
      
      
      ///for pi0->gg endcap
      bool doSelForPi0Endcap_; 
      double selePtGammaEndCap_;
      double selePtPi0EndCap_;
      double region1_Pi0EndCap_;
      double selePtGammaPi0EndCap_region1_; 
      double selePtPi0EndCap_region1_;
      double region2_Pi0EndCap_;
      double selePtGammaPi0EndCap_region2_; 
      double selePtPi0EndCap_region2_;
      double selePtGammaPi0EndCap_region3_; 
      double selePtPi0EndCap_region3_;
      double seleMinvMaxPi0EndCap_;
      double seleMinvMinPi0EndCap_;
      double seleS4S9GammaEndCap_;
      double selePi0IsoEndCap_;
      double selePi0BeltDREndCap_;
      double selePi0BeltDetaEndCap_;
      double ptMinForIsolationEndCap_; 
      bool storeIsoClusRecHitPi0EE_;    
      double selePtPi0MaxEndCap_region3_;
      
      
      ///for eta->gg barrel
      bool doSelForEtaBarrel_; 
      double selePtGammaEta_;
      double selePtEta_;
      double seleS4S9GammaEta_; 
      double seleS9S25GammaEta_; 
      double seleMinvMaxEta_; 
      double seleMinvMinEta_; 
      double ptMinForIsolationEta_; 
      double seleEtaIso_; 
      double seleEtaBeltDR_; 
      double seleEtaBeltDeta_; 
      bool removePi0CandidatesForEta_; 
      double massLowPi0Cand_; 
      double massHighPi0Cand_; 
      bool store5x5RecHitEtaEB_;
      bool store5x5IsoClusRecHitEtaEB_;
      bool storeIsoClusRecHitEtaEB_;


      ///for eta->gg endcap
      bool doSelForEtaEndcap_; 
      double selePtGammaEtaEndCap_;
      double seleS4S9GammaEtaEndCap_;
      double seleS9S25GammaEtaEndCap_;
      double selePtEtaEndCap_;
      double region1_EtaEndCap_;
      double selePtGammaEtaEndCap_region1_; 
      double selePtEtaEndCap_region1_;
      double region2_EtaEndCap_;
      double selePtGammaEtaEndCap_region2_; 
      double selePtEtaEndCap_region2_;
      double selePtGammaEtaEndCap_region3_; 
      double selePtEtaEndCap_region3_;
      double seleMinvMaxEtaEndCap_;
      double seleMinvMinEtaEndCap_;
      double ptMinForIsolationEtaEndCap_;
      double seleEtaIsoEndCap_;
      double seleEtaBeltDREndCap_;
      double seleEtaBeltDetaEndCap_;
      bool storeIsoClusRecHitEtaEE_;
      bool store5x5RecHitEtaEE_;
      bool store5x5IsoClusRecHitEtaEE_;
      
      
      bool doBarrel; 
      bool doEndcap; 
    
      
      
      //pre-scale factors for different eta regions in Endcap, to suppress rate from high eta regions.
      ///region1: |eta|<2; region2: 2<|eta|<2.5; region3: |eta|>2.5;
      int preScale_endcapPi0_region1_; 
      int preScale_endcapPi0_region2_; 
      int preScale_endcapPi0_region3_; 
      
      int preScale_endcapEta_region1_; 
      int preScale_endcapEta_region2_; 
      int preScale_endcapEta_region3_; 
      
      ///global counts for selected 
      long int selected_endcapPi0_region1; 
      long int selected_endcapPi0_region2; 
      long int selected_endcapPi0_region3; 
      
      long int selected_endcapEta_region1; 
      long int selected_endcapEta_region2; 
      long int selected_endcapEta_region3; 
            

      bool storeRecHitES_; 

/*       edm::InputTag preshHitProducer_;         // name of module/plugin/producer producing hits */
/*       edm::InputTag preshHitProducerEta_;         // name of module/plugin/producer producing hits */
      
      ///whether or not to add ES energy when do selections
      bool addESEnergyToEECluster_; 
      bool saveAllPhotonBarrel_;
      bool saveAllPhotonEndcap_;
       
      bool doSelPiz_; 
      bool doSelEta_; 
      
      
      int preshNclust_;
      double preshClustECut;
      double etThresh_;
      double calib_planeX_;
      double calib_planeY_;
      double mip_;
      double gamma_;

/*       PreshowerClusterAlgo * presh_algo; // algorithm doing the real work */

      ///PreshowerClusterAlgo::DebugLevel debugL;  
      // name out output ES cluster collections
      std::string preshClusterCollectionX_;  
      std::string preshClusterCollectionY_;  
     

      bool ParameterLogWeighted_;
      double ParameterX0_;
      double ParameterT0_barl_;
      double ParameterT0_endc_;
      double ParameterT0_endcPresh_;
      double ParameterW0_;
      

/*       edm::InputTag l1IsolatedTag_; */
/*       edm::InputTag l1NonIsolatedTag_; */
/*       edm::InputTag l1SeedFilterTag_; */


      /// std::map<DetId, EcalRecHit> *recHitsEB_map;
      ///replace by two vectors. 

/*       std::vector<EBDetId> detIdEBRecHits;  */
/*       std::vector<EcalRecHit> EBRecHits;  */
 
  
/*       std::vector<EEDetId> detIdEERecHits;  */
/*       std::vector<EcalRecHit> EERecHits;  */

      
/*       std::vector<EBDetId> detIdEBRecHitsEta;  */
/*       std::vector<EcalRecHit> EBRecHitsEta;  */
 
  
/*       std::vector<EEDetId> detIdEERecHitsEta;  */
/*       std::vector<EcalRecHit> EERecHitsEta;  */
      

      
/*       std::vector<EBDetId> alldetIdEBDigi;  */
/*       std::vector<EEDetId> alldetIdEEDigi;  */
      
/*       //      std::vector<EBDetId> alldetIdEBDigiv1; */
/*       //      std::vector<EEDetId> alldetIdEEDigiv1; */


/*       //std::vector<EcalDigi> allEBDigi;  */
/*       // std::vector<dataFrame> allEBDigi;  */
      
/*       edm::ESHandle<EcalGainRatios> pRatio; */
/*       edm::ESHandle<EcalPedestals> pedHandle; */
      
 
      bool Jets_; 

 
/*       edm::InputTag CentralSource_; */
/*       edm::InputTag ForwardSource_; */
/*       edm::InputTag TauSource_; */

      bool JETSdoCentral_ ;
      bool JETSdoForward_ ;
      bool JETSdoTau_ ;
      double Ptmin_jets_; 
      double Ptmin_taujets_; 
      double JETSregionEtaMargin_;
      double JETSregionPhiMargin_;
 

      int theMinBunch; 
      int theMaxBunch; 

      int debug_; 
      bool first_; 
      double EMregionEtaMargin_;
      double EMregionPhiMargin_;
 
      //parameter which decide which level of flags of EcalRecHit used for clustering.


      bool useRecoFlag_; 
      bool useDBStatus_; 
      int flagLevelRecHitsToUse_; 
      int statusLevelRecHitsToUse_;
      //number of minimal rechits for selected clusters.
      int nMinRecHitsSel1stCluster_; 
      int nMinRecHitsSel2ndCluster_; 
      
/*       edm::InputTag m_l1GtRecordInputTag; */
      
/*       edm::InputTag m_vertexSrc;  //offlinePViwthBS */
/*       edm::InputTag m_vertexSrc2; //offlinwPV */
      
      
      
      int isFakeOffPVwithBS; 
      int isValidOffPVwithBS; 
      float xOffPVwithBS;
      float yOffPVwithBS;
      float zOffPVwithBS;
      
      float xErrOffPVwithBS;
      float yErrOffPVwithBS;
      float zErrOffPVwithBS;
      
      int NtrkOffPVwithBS;
      float chi2OffPVwithBS;
      int ndofOffPVwithBS;
      float nchi2OffPVwithBS;
      
      //without BS
      
      int isFakeOffPV; 
      int isValidOffPV; 
      float xOffPV;
      float yOffPV;
      float zOffPV;
      
      float xErrOffPV;
      float yErrOffPV;
      float zErrOffPV;
      
      int NtrkOffPV;
      float chi2OffPV;
      int ndofOffPV;
      float nchi2OffPV;
      
      
      
      
      std::map<std::string,double> providedParameters;  
      
      
      
      std::vector<int> FEDListUsed; ///by regional objects.  ( em, jet, etc)

      std::vector<int> FEDListUsedBarrel; 
      std::vector<int> FEDListUsedEndcap; 

      bool RegionalMatch_;
 
      
      double ptMinEMObj_ ; 
 
/*       std::map<DetId, EcalRecHit> esrechits_map; */
/*       std::set<DetId> used_strips; */
      
/*       edm::InputTag ebDigiCollection_; // collection of EB digis */
/*       edm::InputTag eeDigiCollection_; // collection of EE digis */
      
/*       edm::InputTag ebDigiCollectionv1_; // collection of EB digis */
/*       edm::InputTag eeDigiCollectionv1_; // collection of EE digis */

      
      int nErrorPrinted; 
      
      static const int maxErrorToPrint = 100 ; 
      
      
      bool ecalDigiReFit_;
      

/*       const CaloGeometry * theGeometry; */
      
      //  const CaloSubdetectorGeometry *geometry_eb;
      // const CaloSubdetectorGeometry *geometry_ee;
      // const CaloSubdetectorGeometry *geometry_es;
      //const CaloSubdetectorTopology *topology_eb;
      // const CaloSubdetectorTopology *topology_ee;
      //CaloSubdetectorTopology *topology_es;
      
 
/*       PositionCalc posCalculator_; */
 
      static const int MAXCLUS = 2000;
      static const int MAXPI0S = 200;

      int nEventsProcessed ; 

      int nPassedBSC; 
      int nPassedBSC_noBeamHalo; 
      
      
      std::string outputFile_; 
      

      
      TFile* rootFile_; 

      
      ///mostly vertex information if run on RECO 
      TTree* mytree_evtInfo; 
      TTree* mytree_clusters; 

      TTree* mytree_bs; 
      
      
      TTree* mytree_pizeb;
      TTree* mytree_pizee;

      
      TTree* mytree_etaeb;
      TTree* mytree_etaee;
      
      
      
      TTree* mytree_pizbe;
      

      TTree* mytree_hiteb; 
      
      TTree* mytree_clus; 
      
      
      ///variable
      float mpair_3x3; 
      float mpair_3x3v1; 
      float mpair_3x3v2; 
      float mpair_3x3v3; 
      float mpair_3x3test;  //Z - 60 cm
            
      
      float etapair_3x3; 
      float phipair_3x3; 
      float ptpair_3x3; 
      
      float isolationPiz[5];/// 0, 0.3, 0.5, 0.7,1; 
      float isolationEta[5];/// 0, 0.3, 0.5, 0.7,1; 
      
      ///no belt cone dr <0.3, 0.4, 0.5; 
      float isolationv1[3];
      
      float drpair_3x3; 
      
      float gs4s9min;
      float gs4s6min;
      float gs6s9min;
      float gs9s25min;
      float openangle_3x3; 
      
      int runNumber; 
      int evtNumber; 
      int lumiBlock; 
      
      static const int MAXL1bits = 500;
      static const int MAXHLTbits = 500;
      int nL1bits;
      int L1bits[MAXL1bits];
      
      int nL1bitsTech;
      int L1bitsTech[MAXL1bits];

      int nL1Alca; 
      int L1Alca[MAXL1bits];
      

      int nHLTbits;
      int HLTbits[MAXHLTbits];
      
      ///for each of rechits, now add the sim energy and time
      int nhitSimClus1[25];
      float ehitSimClus1[25][20];
      float thitSimClus1[25][20];
      
      int nhitSimClus2[25];
      float ehitSimClus2[25][20];
      float thitSimClus2[25][20];
      
      ///for each of rechit, add the added sim energy of each time , ( after pulse shape, )
      /// in 
      float sXtalClus1[25][10];
      float sXtalClus2[25][10];
      
      
      ///this is what clustered
      int nxtClus1; 
      float eXtalClus1[25];
      int ietaXtalClus1[25];
      int iphiXtalClus1[25];
      float tXtalClus1[25];
      int fXtalClus1[25];
      int nxtClus2; 
      float eXtalClus2[25];
      int ietaXtalClus2[25];
      int iphiXtalClus2[25];
      float tXtalClus2[25];
      int fXtalClus2[25];

      float laserCorrXtalClus1[25];
      float laserCorrXtalClus2[25];
      

      float xClus1; 
      float yClus1; 
      float zClus1; 
      
      float xClus2;
      float yClus2;
      float zClus2;
      
      
      

      ///this is the 25- 
      
      int nxt5x5Clus1; 
      float eXtal5x5Clus1[25];
      int ietaXtal5x5Clus1[25];
      int iphiXtal5x5Clus1[25];
      float tXtal5x5Clus1[25];
      int fXtal5x5Clus1[25];
      int nxt5x5Clus2; 
      float eXtal5x5Clus2[25];
      int ietaXtal5x5Clus2[25];
      int iphiXtal5x5Clus2[25];
      float tXtal5x5Clus2[25];
      int fXtal5x5Clus2[25];
      
      //the same for all rechits in the cluster. 
      int izXtalClus1; 
      int izXtalClus2; 
            
      


      //digi
      float dXtalClus1[25][11];
      float dXtalClus2[25][11];
      


      float gptmin; 
      
      
      
      int  nxt9v1[2]; 
      
      float E25_3x3[2];
      float E9_3x3[2];
      float E6_3x3[2];
      float E4_3x3[2];
      float Et_3x3[2];
      float pos_3x3[2][3];
      
      
      ///ES cluster X/Y
      float infoESX[2][8]; //[2] for the first / second cluster, [5]; es cluster x,y,z, phi, highest energy strip,eta,phi
      float infoESY[2][8]; //[2] for the first / second cluster, [5]; es cluster e,eta, phi, highest energy strip,eta,phi
      
      
      static const int MAXMC = 100000;
      int nMCpart;
      float pxMCpart[MAXMC];
      float pyMCpart[MAXMC];
      float pzMCpart[MAXMC];
      float ptotMCpart[MAXMC];
      float ptMCpart[MAXMC];
      float eMCpart[MAXMC];
      float etMCpart[MAXMC];
      float etaMCpart[MAXMC];
      float phiMCpart[MAXMC];
      float mMCpart[MAXMC];
      int chargeMCpart[MAXMC];
      int pidMCpart[MAXMC];
      int pidmomMCpart[MAXMC];
      int barcodemomMCpart[MAXMC];
      int nDauMCpart[MAXMC];
      static const int MAXDAU = 100;
      int barcodeDauMCpart[MAXMC][MAXDAU];
      int statusMCpart[MAXMC];
      float vtxXMCpart[MAXMC];
      float vtxYMCpart[MAXMC];
      float vtxZMCpart[MAXMC];
      
      int convPhtMCpart[MAXMC][3]; //only for status==1 Gen photon, and the two index to the two electrons 

      
      int npizallgen; 
      int netaallgen; 
      
      
      static const int MAXGenPIZ = 1000; 

      int nGenpi0; 
      int nGeneta; 
      float etaGenpi0[MAXGenPIZ];
      float phiGenpi0[MAXGenPIZ];
      float etGenpi0[MAXGenPIZ];
      float eGenpi0[MAXGenPIZ];
      float mGenpi0[MAXGenPIZ];
      float dr2phGenpi0[MAXGenPIZ];
      float dr2phGeneta[MAXGenPIZ];
      float etaGeneta[MAXGenPIZ];
      float phiGeneta[MAXGenPIZ];
      float etGeneta[MAXGenPIZ];
      float eGeneta[MAXGenPIZ];
      float ePhtGenpi0[MAXGenPIZ][2];
      float etaPhtGenpi0[MAXGenPIZ][2];
      float phiPhtGenpi0[MAXGenPIZ][2];
      float ePhtGeneta[MAXGenPIZ][2];
      float etaPhtGeneta[MAXGenPIZ][2];
      float phiPhtGeneta[MAXGenPIZ][2];
      float vtxGenpi0[MAXGenPIZ][3];
      float vtxPhtGenpi0[MAXGenPIZ][3];
      float vtxGeneta[MAXGenPIZ][3];
      float vtxPhtGeneta[MAXGenPIZ][3];


      int pidMomGenpi0[MAXGenPIZ];
      int barcodeMomGenpi0[MAXGenPIZ];
      
      int pidMomGeneta[MAXGenPIZ];
      int barcodeMomGeneta[MAXGenPIZ];

      int isConvPhtGeneta[MAXGenPIZ][2];
      float convPht1Geneta[MAXGenPIZ][6]; ///e, eta, phi, e, eta, phi
      float convPht2Geneta[MAXGenPIZ][6]; ///e, eta, phi, e, eta, phi

      int isConvPhtGenpi0[MAXGenPIZ][2];
      float convPht1Genpi0[MAXGenPIZ][6]; ///e, eta, phi, e, eta, phi
      float convPht2Genpi0[MAXGenPIZ][6]; ///e, eta, phi, e, eta, phi
      float convVtxPhtGenpi0[MAXGenPIZ][6];
      float convVtxPhtGeneta[MAXGenPIZ][6];
      
      

/*       float geninfo[14]; */
/*       float geninfoeta[14]; */
      
      
/*       float geninfo[14]; */
/*       float geninfoeta[14]; */
      
      float geninfo[9];
      int convinfo[2];
      

      
      int conv[2];
      
      // L1 objects
      static const int MAXL1OBJ = 20;
      //     static const int MAXL1BITS = 128;
      // L1 EM isolated objects
      int nL1EMIso;
      float L1EMIso_e[MAXL1OBJ];
      float L1EMIso_et[MAXL1OBJ];
      float L1EMIso_eta[MAXL1OBJ];
      float L1EMIso_phi[MAXL1OBJ];
      // L1 EM non-isolated objects
      int nL1EMnonIso;
      float L1EMnonIso_e[MAXL1OBJ];
      float L1EMnonIso_et[MAXL1OBJ];
      float L1EMnonIso_eta[MAXL1OBJ];
      float L1EMnonIso_phi[MAXL1OBJ];
      
      
      bool saveEGObj_;
      bool saveAllRecHitEB_; 
      bool saveAllRecHitEE_; 
     
     
      static const int MAXEB = 61200; 
      int nEB; 
      int ietaEB[MAXEB];
      int iphiEB[MAXEB];
      int rEB[MAXEB]; 
      int digiEB[MAXEB][12];
      float simEB[MAXEB][10];
      float eEB[MAXEB];
     
      static const int MAXEE = 61200; 
      int nEE; 
      int ietaEE[MAXEE];
      int iphiEE[MAXEE];
      int izEE[MAXEE];
      int rEE[MAXEE];
      int digiEE[MAXEE][12];
      float eEE[MAXEE];
      
      int nSeeds; 
      int nClusters; 
            
      
      static const int MAX3x3ClusEB = 2000;
      int n3x3ClusEB;
      float e3x3ClusEB[MAX3x3ClusEB];
      float eta3x3ClusEB[MAX3x3ClusEB];
      float phi3x3ClusEB[MAX3x3ClusEB];
      float laserCorr3x3ClusEB[MAX3x3ClusEB][9];
      ///tmp by simplelog
      float leta3x3ClusEB[MAX3x3ClusEB];
      float lphi3x3ClusEB[MAX3x3ClusEB];
      
      float x3x3ClusEB[MAX3x3ClusEB];
      float y3x3ClusEB[MAX3x3ClusEB];
      float z3x3ClusEB[MAX3x3ClusEB];
      
      int nXt3x3ClusEB[MAX3x3ClusEB];
      float eXt3x3ClusEB[MAX3x3ClusEB][9];
      float tXt3x3ClusEB[MAX3x3ClusEB][9];
      int ietaXt3x3ClusEB[MAX3x3ClusEB][9]; ///index from -85,0,84; ,
      int iphiXt3x3ClusEB[MAX3x3ClusEB][9]; ///phi, 0,359
      float etaXt3x3ClusEB[MAX3x3ClusEB][9];
      float phiXt3x3ClusEB[MAX3x3ClusEB][9];
      
      
      float xXt3x3ClusEB[MAX3x3ClusEB][9];
      float yXt3x3ClusEB[MAX3x3ClusEB][9];
      float zXt3x3ClusEB[MAX3x3ClusEB][9];
      float S43x3ClusEB[MAX3x3ClusEB];
      float S63x3ClusEB[MAX3x3ClusEB];
      float S253x3ClusEB[MAX3x3ClusEB];
      
      float s4s93x3ClusEB[MAX3x3ClusEB];
      float s6s93x3ClusEB[MAX3x3ClusEB];
      float s9s253x3ClusEB[MAX3x3ClusEB];
      
      
    
      static const int MAX3x3ClusEE = 2000;
      int n3x3ClusEE;
      float e3x3ClusEE[MAX3x3ClusEE];
      float eta3x3ClusEE[MAX3x3ClusEE];
      float phi3x3ClusEE[MAX3x3ClusEE];

      ///tmp by simplelog
      float leta3x3ClusEE[MAX3x3ClusEE];
      float lphi3x3ClusEE[MAX3x3ClusEE];
      

      float x3x3ClusEE[MAX3x3ClusEE];
      float y3x3ClusEE[MAX3x3ClusEE];
      float z3x3ClusEE[MAX3x3ClusEE];

      int nXt3x3ClusEE[MAX3x3ClusEE];
      float eXt3x3ClusEE[MAX3x3ClusEE][9];
      float tXt3x3ClusEE[MAX3x3ClusEE][9];
      int ixXt3x3ClusEE[MAX3x3ClusEE][9]; ///index from -85,0,84; ,
      int iyXt3x3ClusEE[MAX3x3ClusEE][9]; ///phi, 0,359
      int izXt3x3ClusEE[MAX3x3ClusEE][9]; ///phi, 0,359
      float laserCorr3x3ClusEE[MAX3x3ClusEE][9];
      
      float etaXt3x3ClusEE[MAX3x3ClusEE][9];
      float phiXt3x3ClusEE[MAX3x3ClusEE][9];
      
      
      float xXt3x3ClusEE[MAX3x3ClusEE][9];
      float yXt3x3ClusEE[MAX3x3ClusEE][9];
      float zXt3x3ClusEE[MAX3x3ClusEE][9];
      float S43x3ClusEE[MAX3x3ClusEE];
      float S63x3ClusEE[MAX3x3ClusEE];
      float S253x3ClusEE[MAX3x3ClusEE];
      
      float s4s93x3ClusEE[MAX3x3ClusEE];
      float s6s93x3ClusEE[MAX3x3ClusEE];
      float s9s253x3ClusEE[MAX3x3ClusEE];
      

      double hfETowerh_; 
      
      int nPxlHits ; 
      float clusVtxQual; 
      
      int    nHfTowersP     ;
      int    nHfTowersN     ; 

      float sumHfEsubEpPlus; 
      float sumHfEsubEpMinus; 
      
      float sumHfEaddEpPlus; 
      float sumHfEaddEpMinus; 
      
      
      int phyDeclared; 
      
      float highPurityTrackFrac; 
      
//      edm::InputTag m_tracksSrc;
      static const int MAXTRK = 5000;
      int nTracks;
      float pxTracks[MAXTRK];
      float pyTracks[MAXTRK];
      float pzTracks[MAXTRK];
      float pTracks[MAXTRK];
      float ptTracks[MAXTRK];
      float etaTracks[MAXTRK];
      float phiTracks[MAXTRK];
      int pidTracks[MAXTRK];
      int chargeTracks[MAXTRK]; 
      float vxTracks[MAXTRK];
      float vyTracks[MAXTRK];
      float vzTracks[MAXTRK];
      int nhitsTracks[MAXTRK];
      float nChi2Tracks[MAXTRK];
      int algoTracks[MAXTRK];
      int ndofTracks[MAXTRK];
      int nValidhitsTracks[MAXTRK];
      int nValidpixelhitsTracks[MAXTRK];
      int nValidstriphitsTracks[MAXTRK];
      int qualityFlagTracks[MAXTRK];
      
      float vBeamSpot[3];
      std::vector<float> *beamspotv;
      bool isNewLumiBlock; 
            
      ///generated ch
      int nChaGen; 
      float ptChaGen[MAXMC];
      float etaChaGen[MAXMC];
      int pidChaGen[MAXMC];
      int chaChaGen[MAXMC];
      int staChaGen[MAXMC];
      
      
      bool reRunPixelRecHits_; 
      

      ////peak after some selection. 
      TH1F *hh_mpair[10]; 
      
      
      TH1F *hh_vtx_mc[10];
      
      TH1F *hh_L1bitFired;
      TH1F *hh_L1bitTechFired;
      
      int nL1bits_fired[200];
      int nL1bitsTech_fired[200];
      

      ///trackerOnlyConversion
      
      static const int nMaxTrkOnlyConv = 1000; 
      int nTrkOnlyConv; 
      float vtxLTrkOnlyConv[nMaxTrkOnlyConv][3];
      float vtxRTrkOnlyConv[nMaxTrkOnlyConv][3];
      float nChi2LTrkOnlyConv[nMaxTrkOnlyConv];
      float nChi2RTrkOnlyConv[nMaxTrkOnlyConv];
      int rechitSizeLTrkOnlyConv[nMaxTrkOnlyConv];
      int rechitSizeRTrkOnlyConv[nMaxTrkOnlyConv];
      float zPVTrkOnlyConv[nMaxTrkOnlyConv];
      
      float nChi2vtxTrkOnlyConv[nMaxTrkOnlyConv];
      int ndofvtxTrkOnlyConv[nMaxTrkOnlyConv];
      

      float ptLTrkOnlyConv[nMaxTrkOnlyConv];
      float etaLTrkOnlyConv[nMaxTrkOnlyConv];
      float phiLTrkOnlyConv[nMaxTrkOnlyConv];
            
      float ptRTrkOnlyConv[nMaxTrkOnlyConv];
      float etaRTrkOnlyConv[nMaxTrkOnlyConv];
      float phiRTrkOnlyConv[nMaxTrkOnlyConv];
      

      
      float vtxTrkOnlyConv[nMaxTrkOnlyConv][3];
      float dmApproachTrkOnlyConv[nMaxTrkOnlyConv];
      float dphiTrackatVtxTrkOnlyConv[nMaxTrkOnlyConv];
      float pairPtTrkOnlyConv[nMaxTrkOnlyConv];
      float pairEtaTrkOnlyConv[nMaxTrkOnlyConv];
      float pairPhiTrkOnlyConv[nMaxTrkOnlyConv];
      float pairMTrkOnlyConv[nMaxTrkOnlyConv];
      

      float pairCotThetaTrkOnlyConv[nMaxTrkOnlyConv];
      float dPhiTracksAtEcalTrkOnlyConv[nMaxTrkOnlyConv];
      float dEtaTracksAtEcalTrkOnlyConv[nMaxTrkOnlyConv];

///      edm::InputTag convLabel; 
      

      int nClus; 
      
      ///this is for pi0->gg barrel 
      vector<float> eClus;
      vector<float> eClusv1;
      vector<float> etClus;
      vector<float> etaClus;
      vector<float> thetaClus;
      vector<float> phiClus;
//      vector< vector<EcalRecHit> > RecHitsCluster;
//     vector< vector<EcalRecHit> > RecHitsCluster5x5;
      vector<float> s4s9Clus;
      vector<float> s4s6Clus;
      vector<float> s6s9Clus;
      vector<float> s9s25Clus;
      vector<float> e4Clus; 
      vector<float> e6Clus; 
      vector<float> e9Clus; 
      vector<float> e25Clus; 
      vector<float> xClus;
      vector<float> yClus;
      vector<float> zClus;
      vector<int> nXtv1Clus; ///number of Xtal , no matter if it is clustered or not
  
      
      float mpair; 
      float mpairv1; 
      float ptpair; 
      float etapair; 
      float phipair; 
      float ptmin; 
      float isolation;
      
      int genMatched; 
      
      float s4s9min; 
      float s9s25min;
      
      ///adding L1 trigger info
      int L1_SingleIsoEG5;
      int L1_SingleIsoEG8;
      int L1_SingleIsoEG10;
      int L1_SingleIsoEG12;
      int L1_SingleIsoEG15;
      int L1_SingleEG2;
      int L1_SingleEG5;
      int L1_SingleEG8;
      int L1_SingleEG10;
      int L1_SingleEG12;
      int L1_SingleEG15;
      int L1_SingleEG20;
      int L1_SingleJet6U;
      int L1_SingleJet10U;
      int L1_SingleJet20U;
      int L1_SingleJet30U;
      int L1_SingleJet40U;
      int L1_SingleJet50U;
      int L1_DoubleJet30U;
      int L1_DoubleEG5;
      

      int fullRECO_; 
      
      bool removeSpike_;
      
      
      
      TH1I *hh_mul_ieta; 
      TH1I *hh_mul_iphi; 
      
      
      TH1I *hh_mul_ietaSeed; 
      TH1I *hh_mul_iphiSeed; 
      
      
      TH1I *hh_mul_ietaSeedClus; 
      TH1I *hh_mul_iphiSeedClus; 
      

      TH1I *hh_mul_ietaSeedSel; 
      TH1I *hh_mul_iphiSeedSel; 

/*       HLTConfigProvider hltConfig_;        // to get configuration for L1s/Pre */

      std::vector<std::string>  hlNames_; 

      bool getBeamSpotOnly_;
      int curLumiBlock;

      //vertex 
      int nVertex; 
      static const int MAXVX = 100;
      float vertexx[MAXVX];
      float vertexy[MAXVX];
      float vertexz[MAXVX];
      float vertexchi2[MAXVX];
      float vertexndof[MAXVX];
      float vertexnormalizedChi2[MAXVX];
      int vertextrackSize[MAXVX];
      int vertexisFake[MAXVX];
      int vertexisValid[MAXVX];
      float vertexsumtrackspt[MAXVX];
      float vertexsumtracksptSquare[MAXVX];
      float vertexsumRefittedtrackspt[MAXVX];
      float vertexsumRefittedtracksptSquare[MAXVX];

      std::vector<unsigned short>* hlt_bitFired;
      std::vector<std::string> *hlt_pathName;
      std::vector<unsigned short>* l1bitFired;
      std::vector<std::string> *l1algoName;
      std::vector<std::string> *l1algoNameRun;

      int isRealData;

int nSeedsEB;
int nSeedsEE;

      
/* }; */
