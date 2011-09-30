
/* // system include files */
/* #include <memory> */

/* // user include files */
/* #include "FWCore/Framework/interface/Frameworkfwd.h" */
/* #include "FWCore/Framework/interface/EDAnalyzer.h" */

/* #include "FWCore/Framework/interface/Event.h" */
/* #include "FWCore/Framework/interface/MakerMacros.h" */

/* #include "FWCore/ParameterSet/interface/ParameterSet.h" */
/* // */
/* // class declaration */
/* // */


/* //GEN */
/* #include "SimDataFormats/GeneratorProducts/interface/GenEventInfoProduct.h" */
/* #include "SimDataFormats/GeneratorProducts/interface/HepMCProduct.h" */
/* #include "DataFormats/HepMCCandidate/interface/GenParticle.h" */
/* #include "SimDataFormats/GeneratorProducts/interface/PdfInfo.h" */


/* ////SimTrack */
/* #include "SimDataFormats/Track/interface/SimTrack.h" */
/* #include "SimDataFormats/Track/interface/SimTrackContainer.h" */
/* #include "SimDataFormats/Vertex/interface/SimVertexContainer.h" */


/* //BEamSpot */
/* #include "DataFormats/BeamSpot/interface/BeamSpot.h" */


/* //Track */
/* #include "DataFormats/TrackReco/interface/Track.h" */
/* #include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"  */
/* #include "DataFormats/TrackCandidate/interface/TrackCandidate.h"  */
/* #include "DataFormats/TrackReco/interface/Track.h"  */
/* #include "DataFormats/TrackReco/interface/TrackFwd.h"  */

/* //HLT */
/* #include "DataFormats/HLTReco/interface/TriggerEvent.h" */
/* #include "DataFormats/HLTReco/interface/TriggerObject.h" */
/* #include "DataFormats/Common/interface/TriggerResults.h" */
/* #include "DataFormats/HLTReco/interface/TriggerTypeDefs.h" */
/* #include "DataFormats/HLTReco/interface/TriggerEventWithRefs.h" */
/* #include "HLTrigger/HLTcore/interface/HLTConfigProvider.h" */


/* //L1 Trigger */
/* #include "DataFormats/L1Trigger/interface/L1EmParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1JetParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1MuonParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1EtMissParticle.h" */
/* #include "DataFormats/L1Trigger/interface/L1ParticleMapFwd.h" */
/* #include "DataFormats/L1Trigger/interface/L1ParticleMap.h" */
/* #include "L1Trigger/L1ExtraFromDigis/interface/L1ExtraParticleMapProd.h" */
/* #include "DataFormats/L1GlobalTrigger/interface/L1GlobalTriggerReadoutRecord.h" */
/* #include "DataFormats/L1GlobalTrigger/interface/L1GtFdlWord.h" */



/* ///MET */
/* #include "DataFormats/METReco/interface/CaloMET.h" */
/* #include "DataFormats/METReco/interface/CaloMETFwd.h" */
/* #include "DataFormats/METReco/interface/MET.h" */
/* #include "DataFormats/METReco/interface/METFwd.h" */
/* #include "DataFormats/METReco/interface/PFMET.h" */
/* #include "DataFormats/METReco/interface/PFMETFwd.h" */

/* ///muon */
/* #include "DataFormats/MuonReco/interface/MuonCocktails.h" */
/* #include "RecoMuon/TrackingTools/interface/MuonServiceProxy.h" */
/* #include "DataFormats/MuonReco/interface/Muon.h" */
/* #include "DataFormats/MuonReco/interface/MuonFwd.h" */
/* #include "DataFormats/MuonReco/interface/MuonSelectors.h" */
/* #include "DataFormats/MuonReco/interface/MuonTime.h" */
/* #include "DataFormats/MuonReco/interface/MuonTimeExtra.h" */
/* #include "DataFormats/MuonReco/interface/MuonTimeExtraMap.h" */
/* #include "DataFormats/MuonReco/interface/MuonCosmicCompatibility.h" */

/* ///ECAL */
/* #include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h" */
/* #include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h" */
/* #include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h" */
/* #include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h" */
/* #include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h" */
/* #include "CondFormats/EcalObjects/interface/EcalIntercalibConstants.h" */
/* #include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h" */
/* #include "CondFormats/EcalObjects/interface/EcalADCToGeVConstant.h" */
/* #include "CondFormats/EcalObjects/interface/EcalChannelStatus.h" */
/* #include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbService.h" */
/* #include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h" */
/* #include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgoRcd.h" */
/* #include "RecoEcal/EgammaCoreTools/interface/EcalTools.h" */
/* #include "DataFormats/CaloRecHit/interface/CaloCluster.h" */


/* //IsoDep */
/* #include "DataFormats/RecoCandidate/interface/IsoDeposit.h" */
/* #include "DataFormats/RecoCandidate/interface/IsoDepositFwd.h" */
/* #include "DataFormats/RecoCandidate/interface/IsoDepositDirection.h" */


/* ///vertex */
/* #include "DataFormats/VertexReco/interface/Vertex.h" */
/* #include "DataFormats/VertexReco/interface/VertexFwd.h" */


/* ///caloCluster */
/* #include "DataFormats/CaloRecHit/interface/CaloCluster.h" */
/* #include "DataFormats/CaloRecHit/interface/CaloClusterFwd.h" */


/* //covnersion */
/* #include "DataFormats/EgammaCandidates/interface/Conversion.h" */
/* #include "DataFormats/EgammaCandidates/interface/ConversionFwd.h" */
/* #include "RecoEgamma/EgammaTools/interface/ConversionTools.h" */


/* ///Egamma */
/* #include "DataFormats/EgammaCandidates/interface/PhotonFwd.h" */
/* #include "DataFormats/EgammaCandidates/interface/Photon.h" */
/* #include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h" */
/* #include "DataFormats/EgammaCandidates/interface/GsfElectron.h" */
/* #include "DataFormats/GsfTrackReco/interface/GsfTrackFwd.h" */
/* #include "DataFormats/GsfTrackReco/interface/GsfTrack.h" */


/* //PU */
/* #include "SimDataFormats/PileupSummaryInfo/interface/PileupSummaryInfo.h"  */

/* //Jet */
/* #include "DataFormats/JetReco/interface/Jet.h" */
/* #include "DataFormats/JetReco/interface/CaloJetCollection.h" */
/* #include "DataFormats/JetReco/interface/GenJet.h" */
/* #include "RecoJets/JetProducers/interface/JetIDHelper.h" */
/* #include "JetMETCorrections/Objects/interface/JetCorrector.h" */
/* #include "DataFormats/JetReco/interface/PFJetCollection.h" */


/* ///CommonTools */
/* #include "CommonTools/Statistics/interface/ChiSquaredProbability.h" */

/* ///ROOT  */
/* #include "TROOT.h" */
/* #include "TFile.h" */
/* #include "TTree.h" */
/* #include "TBranch.h" */
/* #include "TLorentzVector.h" */
/* #include "TClonesArray.h" */

/* #include "TMatrixDSym.h" */
/* #include "TMatrixDSymEigen.h" */
/* #include "TVectorD.h" */

/* ///c++ */
/* #include <boost/foreach.hpp> */
/* #include <vector> */
/* #include <ext/algorithm> */
/* #include <numeric> */
/* #include <iterator> */
/* #include <algorithm> */


/* using namespace std; */
/* using namespace edm; */


/* class RecoAnalyzer : public edm::EDAnalyzer { */
/*    public: */
/*       explicit RecoAnalyzer(const edm::ParameterSet&); */
/*       ~RecoAnalyzer(); */

/*       static void fillDescriptions(edm::ConfigurationDescriptions& descriptions); */
     
      
      
  
/*    private: */
/*       virtual void beginJob() ; */
/*       virtual void analyze(const edm::Event&, const edm::EventSetup&); */
/*       virtual void endJob() ; */

/*       virtual void beginRun(edm::Run const&, edm::EventSetup const&); */
/*       virtual void endRun(edm::Run const&, edm::EventSetup const&); */
/*       virtual void beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&); */
/*       virtual void endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&); */

/*       // ----------member data --------------------------- */
/*       double DeltaPhi(double v1, double v2);  */
/*       double GetDeltaR(double eta1, double eta2, double phi1, double phi2);  */
/*       int indexofParticle(float px, float pz, int status); */
/*       int getMotherIndex(int ); */
/*       void MatchToGenMuon(float eta, float phi, float res[]); */
/*       void MatchToGenElectron(float eta, float phi, float res[]); */
/*       void MatchToGenPhoton(float eta, float phi, float res[]); */
/*       int partonMatchingAlgo(float,float,float); */
/*       int partonMatchingPhys(float,float); */
/*       int findIndexofSC(float en,float eta,float phi); */
/*       int findIndexofTRK(float, float, float);  */
/*       double mye2overe9( const DetId id, const EcalRecHitCollection & recHits);  */
/*       float recHitE( const  DetId id,  const EcalRecHitCollection &recHits ); */
/*       float recHitE( const DetId id, const EcalRecHitCollection & recHits,int di, int dj );  */
/*       float recHitApproxEt(  const DetId id,  const EcalRecHitCollection &recHits ); */
/*       TLorentzVector photonp4wrtvertex(int, int); */
/*       double etaTransformation(  float EtaParticle , float Zvertex);  */
/*       int matchPhotonToConversion( int lpho);  */
/*       void higgsVertexAnalysis(int, int); */
/*       double HggVertexFromConversionsvtxZ(int j); */
/*       double HggVertexFromConversionsvtxdZ(int j); */
/*       bool filter(const edm::Event& iEvent, const edm::EventSetup& iSetup); */
      
/*       struct LessById { */
/* 	bool operator()(const SimTrack &tk1, const SimTrack &tk2) const { return tk1.trackId() < tk2.trackId(); } */
/* 	bool operator()(const SimTrack &tk1, unsigned int    id ) const { return tk1.trackId() < id;            } */
/* 	bool operator()(unsigned int     id, const SimTrack &tk2) const { return id            < tk2.trackId(); } */
/*       }; */
      
      
/*       edm::InputTag trigResultTag_;  */
/*       edm::InputTag trigSummaryTag_; */
/*       int debug_;  */
/*       HLTConfigProvider configProvider; */
      
  
  
/*   ///all inuptTag start with "m_xxx" */
/*   edm::InputTag m_vertexSrc;  //offlinePV */
/*   edm::InputTag m_vertexSrc2; //offlinwPVwithBS */
/*   edm::InputTag m_tracksSrc; ///general track */
/*   edm::InputTag m_beamSpot;  ///beamSpot */
/*   edm::InputTag m_muonSrc;  ///muon */
/*   edm::InputTag m_ecalRecHitBarrel;  //ecalrechit barrel */
/*   edm::InputTag m_ecalRecHitEndcap;  //endcap */
  
/*   edm::InputTag m_barrelSC;  */
/*   edm::InputTag m_endcapSC;  */
/*   edm::InputTag m_photonSrc; */
/*   edm::InputTag m_electronSrc; */
/*   edm::InputTag m_PileupSrc;  */
/*   edm::InputTag m_rhoCollection;  */
/*   edm::InputTag m_rhoCollection2;  */
/*   edm::InputTag m_allConversionsColl;  */
/*   edm::InputTag m_l1GtRecordInputTag;  */
/*   edm::InputTag m_l1IsolatedTag;  */
/*   edm::InputTag m_l1NonIsolatedTag;  */
/*   edm::InputTag m_l1MuTag;  */
  
  
/*   InputTag muIsoTrkMap_; */
/*   InputTag muIsoEcalMap_; */
/*   InputTag muIsoHcalMap_;  */
/*   InputTag muIsoHoMap_; */
  
/*   bool m_dataFormat;  */
  
  
/*   //OUTPUT */
/*   TFile* rootFile_;  */
/*   TTree *Analysis;  */
/*   std::string outputFile_;  */
  
  
  
  
  //Event info
  int lumiBlock; 
  int runNumber; 
  int evtNumber; 
  int bunchX; 
  int orbitNumber; 
  int evtTime; 
  int isRealData;
    
  //HLT
  static const int MAXHLTbits = 500;
  
  
  static const int MAX_HLT = 5000;
  
  std::vector<unsigned short>* hlt_bitFired;
  
  // MAX_HLT hlt objects pt/eta/phi/energy

  //Int_t hlt_n;
  //float hlt_pt[MAX_HLT];
  //float hlt_eta[MAX_HLT];
  //float hlt_phi[MAX_HLT];
  //float hlt_en[MAX_HLT];
  
  vector<float> *hlt_pt; 
  vector<float> *hlt_eta; 
  vector<float> *hlt_phi; 
  vector<float> *hlt_en; 
  
  
  
  std::vector<std::vector<unsigned short> >* hlt_candpath;
  std::vector<std::vector<unsigned short> >* hltfilter_candpath;
  std::vector<std::string> *hlt_pathName;
  std::vector<std::string> *hlt_path_names_HLT2;
///std::vector<edm::InputTag> theHLTLabels;
  TClonesArray* hlt_p4;
  
  


  std::vector<std::string> *hlt_label_names;
  
  static const int MAXSC= 300; 
  int nSC; 
  float ptSC[MAXSC];
  float eSC[MAXSC];
  float eRawSC[MAXSC];
  float etaSC[MAXSC];
  float phiSC[MAXSC];
  int flagSC[MAXSC];
    
  ///anlaysis objects
  
  static const int nPhotonMAX = 100;
  int nPhoton;
  float photonsigmaIetaIeta[nPhotonMAX];
  int photonhasPixelSeed[nPhotonMAX];
  float photonenergy[nPhotonMAX];
  float photonpt[nPhotonMAX];
  float photoneta[nPhotonMAX];
  float photonphi[nPhotonMAX];
  float photonvertexx[nPhotonMAX];
  float photonvertexy[nPhotonMAX];
  float photonvertexz[nPhotonMAX];
  int photonhasConversionTracks[nPhotonMAX];
  float photonscrawEnergy[nPhotonMAX];
  float photonsceta[nPhotonMAX];
  float photonscphi[nPhotonMAX];
  float photoncaloPositionx[nPhotonMAX];
  float photoncaloPositiony[nPhotonMAX];
  float photoncaloPositionz[nPhotonMAX];
  float photonscenergy[nPhotonMAX];
  float photone3x3[nPhotonMAX]; 
  float photone1x5[nPhotonMAX]; 
  float photone2x5[nPhotonMAX]; 
  float    photone5x5[nPhotonMAX]; 
  float   photonmaxEnergyXtal[nPhotonMAX]; 
  float   photonr9[nPhotonMAX]; 

  

  //more 
  float photonscpreshowerEnergy[nPhotonMAX]; 
  float photonscphiWidth[nPhotonMAX]; 
  float photonscetaWidth[nPhotonMAX]; 
  int photonscclusterSize[nPhotonMAX]; 
  std::vector<std::vector<float> >* photonscbclusterenergy;
  
  int photonscnhits[nPhotonMAX]; 
  float photoneLeft[nPhotonMAX];
  float photoneRight[nPhotonMAX];
  float photoneBottom[nPhotonMAX];
  float photoneTop[nPhotonMAX];

float photone1x3[nPhotonMAX];
float photone3x1[nPhotonMAX];
float photone2x2[nPhotonMAX];
float photone3x2[nPhotonMAX];
float photone4x4[nPhotonMAX];
float photone2x5Right[nPhotonMAX];
float photone2x5Left[nPhotonMAX];
float photone2x5Top[nPhotonMAX];
float photone2x5Bottom[nPhotonMAX];
float photone2x5Max[nPhotonMAX];
///float photonenergyBasketFractionEta[nPhotonMAX];
//float photonenergyBasketFractionPhi[nPhotonMAX];

float photonlat[nPhotonMAX][3];
float photonCovEtaEta[nPhotonMAX];
float photonCovEtaPhi[nPhotonMAX];
float photonCovPhiPhi[nPhotonMAX];

float photonCovIEtaIEta[nPhotonMAX];
float photonCovIEtaIPhi[nPhotonMAX];
float photonCovIPhiIPhi[nPhotonMAX];
float photonscCovIEtaIEta[nPhotonMAX];
float photonscCovIEtaIPhi[nPhotonMAX];
float photonscCovIPhiIPhi[nPhotonMAX];
float photonzernike20[nPhotonMAX];
float photonzernike42[nPhotonMAX];





  //isolation
  float    photonhadronicOverEm[nPhotonMAX]; 
  float    photonecalRecHitSumEtConeDR03[nPhotonMAX]; 
  float    photonhcalDepth1TowerSumEtConeDR03[nPhotonMAX]; 
  float    photonhcalDepth2TowerSumEtConeDR03[nPhotonMAX]; 
  float    photonhcalTowerSumEtConeDR03[nPhotonMAX]; 
  float   photontrkSumPtHollowConeDR03[nPhotonMAX]; 
  float   photontrkSumPtSolidConeDR03[nPhotonMAX]; 
  int    photonnTrkHollowConeDR03[nPhotonMAX]; 
  int   photonnTrkSolidConeDR03[nPhotonMAX]; 
  float   photonecalRecHitSumEtConeDR04[nPhotonMAX]; 
  float   photonhcalDepth1TowerSumEtConeDR04[nPhotonMAX]; 
  float   photonhcalDepth2TowerSumEtConeDR04[nPhotonMAX]; 
  float   photonhcalTowerSumEtConeDR04[nPhotonMAX]; 
  float photontrkSumPtHollowConeDR04[nPhotonMAX]; 
  float   photontrkSumPtSolidConeDR04[nPhotonMAX]; 
  int   photonnTrkHollowConeDR04[nPhotonMAX]; 
  int   photonnTrkSolidConeDR04[nPhotonMAX]; 
  
  int photonconversionsize[nPhotonMAX]; 
  float photonconversionVertexx[nPhotonMAX]; 
  float photonconversionVertexy[nPhotonMAX]; 
  float photonconversionVertexz[nPhotonMAX]; 
  float photonconversionrefittedPairMomentumx[nPhotonMAX];
  float photonconversionrefittedPairMomentumy[nPhotonMAX];
  float photonconversionrefittedPairMomentumz[nPhotonMAX];

  float photonconversionpairInvariantMass[nPhotonMAX];
  float photonconversionpairCotThetaSeparation[nPhotonMAX];
  float photonconversionEoverPrefittedTracks[nPhotonMAX];
  float photonconversionzOfPrimaryVertexFromTracks[nPhotonMAX];
  float photonconversiondistOfMinimumApproach[nPhotonMAX];
  float photonconversiondPhiTracksAtVtx[nPhotonMAX];
  float photonconversiondPhiTracksAtEcal[nPhotonMAX];
  float photonconversiondEtaTracksAtEcal[nPhotonMAX];
  int photonconversionnTracks[nPhotonMAX];
  float photonconversionMVAout[nPhotonMAX];
  int photonconversionVertexisValid[nPhotonMAX];
  float photonconversionVertexchi2[nPhotonMAX];
  float photonconversionChiSquaredProbability[nPhotonMAX];
  float photonconversion_track1_dz[nPhotonMAX];
float photonconversion_track1_dzError[nPhotonMAX];
int photonconversion_track1_charge[nPhotonMAX];
float photonconversion_track1_d0[nPhotonMAX];
float photonconversion_track1_tracksPout[nPhotonMAX];
float photonconversion_track1_tracksPin[nPhotonMAX];
int photonconversion_track1_algo[nPhotonMAX];
float photonconversion_track2_dz[nPhotonMAX];
float photonconversion_track2_dzError[nPhotonMAX];
int photonconversion_track2_charge[nPhotonMAX];
int photonconversion_track2_algo[nPhotonMAX];
float photonconversion_track2_d0[nPhotonMAX];
float photonconversion_track2_tracksPout[nPhotonMAX];
float photonconversion_track2_tracksPin[nPhotonMAX];



  float photonseedtime[nPhotonMAX];
  float photonseedoutOfTimeChi2[nPhotonMAX]; 
  float photonseedchi2[nPhotonMAX]; 
  int photonseedrecoFlag[nPhotonMAX]; 
  int photonseedseverityLevel[nPhotonMAX]; 
  int photonfiducialFlag[nPhotonMAX]; 
  int photonscindex[nPhotonMAX]; 
  float photonswissCross[nPhotonMAX];
  int photonieta[nPhotonMAX];
  int photoniphi[nPhotonMAX];
  
  float photonE2overE9[nPhotonMAX];
  int photonmatchToallConv[nPhotonMAX];
  
  float photongenphtmatch[nPhotonMAX][3];
  float photongenelematch[nPhotonMAX][3];
  float photongenphtconv[nPhotonMAX][4];
  float photonPartonmatch[nPhotonMAX][6];
  float photonPartonmatchp[nPhotonMAX][3];
  
  int photonhasMatchedPromptElectron[nPhotonMAX];


  TClonesArray *photonp4;  
  
  int nEventsProcessed; 
  
  float genhiggsm;
  float genhiggspt;
  float genhiggseta;
  float genhiggsphi;
  int genhiggsstatus;
  float genhiggsvx;
  float genhiggsvy;
  float genhiggsvz;
  
  
  static const int nConvMAX = 1000;
  int nConv; 
  float convrefittedPair4Momentumeta[nConvMAX]; 
  float convrefittedPair4Momentumphi[nConvMAX]; 
  float convrefittedPair4Momentumpt[nConvMAX]; 
  float convrefittedPair4Momentumenergy[nConvMAX]; 
  int convnTracks[nConvMAX];
  int convcaloClustersize[nConvMAX]; 
  float convcaloCluster0eta[nConvMAX]; 
  float convcaloCluster0x[nConvMAX]; 
  float convcaloCluster0y[nConvMAX]; 
  float convcaloCluster0z[nConvMAX]; 
  float convcaloCluster0phi[nConvMAX]; 
  int convconversionVertexisValid[nConvMAX];
  float convpairMomentumx[nConvMAX]; 
  float convpairMomentumy[nConvMAX]; 
  float convpairMomentumz[nConvMAX]; 
  float convrefittedPairMomentumx[nConvMAX]; 
  float convrefittedPairMomentumy[nConvMAX]; 
  float convrefittedPairMomentumz[nConvMAX]; 
  float convconversionVertexx[nConvMAX]; 
  float convconversionVertexy[nConvMAX]; 
  float convconversionVertexz[nConvMAX]; 
  float convconversionVertexchi2[nConvMAX]; 
  float convconversionVertexChiSquaredProbability[nConvMAX]; 
  float    convconversionVertexxError[nConvMAX]; 
  float    convconversionVertexyError[nConvMAX]; 
  float    convconversionVertexzError[nConvMAX]; 
  float    convconversionVertexnTracks[nConvMAX]; 
  float   convconversionVertexMVAout[nConvMAX]; 
  float conv_track1_dz[nConvMAX]; 
  float     conv_track1_dzError[nConvMAX]; 
  int     conv_track1_charge[nConvMAX]; 
  int     conv_track1_algo[nConvMAX]; 

  float     conv_track2_dz[nConvMAX]; 
  float     conv_track2_dzError[nConvMAX]; 
  int     conv_track2_charge[nConvMAX]; 
  int     conv_track2_algo[nConvMAX]; 
  
  float  convpairInvariantMass[nConvMAX]; 
  float   convpairCotThetaSeparation[nConvMAX]; 
    // will work in 420 conv_eoverp[nConvMAX]=localConv.EoverPrefittedTracks();
  float   convEoverPrefittedTracks[nConvMAX]; 
  float  convzOfPrimaryVertexFromTracks[nConvMAX]; 
  float   convdistOfMinimumApproach[nConvMAX]; 
  float   convdPhiTracksAtVtx[nConvMAX]; 
  float   convdPhiTracksAtEcal[nConvMAX]; 
  float   convdEtaTracksAtEcal[nConvMAX]; 
  
  int convnSharedHits[nConvMAX]; 
  float  conv_track1_d0[nConvMAX]; 
 float    conv_track1_pout[nConvMAX]; 
 float     conv_track1_pin[nConvMAX]; 
 float conv_track2_d0[nConvMAX]; 
 float conv_track2_pout[nConvMAX]; 
 float conv_track2_pin[nConvMAX]; 

 float conv_track1_pt[nConvMAX]; 
 float conv_track1_eta[nConvMAX]; 
 float conv_track1_phi[nConvMAX]; 
 float conv_track2_pt[nConvMAX]; 
 float conv_track2_eta[nConvMAX]; 
 float conv_track2_phi[nConvMAX]; 
 


 std::vector<std::vector<unsigned short> >* convnHitsBeforeVtx;
 
  
  int nElectron;
  static const int nElectronMAX = 100;
 
float electronpt[nElectronMAX]; 
float electroneta[nElectronMAX];
float electronphi[nElectronMAX];
int electroncharge[nElectronMAX];
float electronvertexx[nElectronMAX];
float electronvertexy[nElectronMAX];
float electronvertexz[nElectronMAX];
float electronscrawEnergy[nElectronMAX];
float electronsceta[nElectronMAX];
 float electronscphi[nElectronMAX];
 float electroncaloPositionx[nElectronMAX];
 float electroncaloPositiony[nElectronMAX];
 float electroncaloPositionz[nElectronMAX];
 
 float electronscpreshowerEnergy[nElectronMAX]; 
 float electronscphiWidth[nElectronMAX]; 
 float electronscetaWidth[nElectronMAX]; 
 int electronscclusterSize[nElectronMAX]; 
 std::vector<std::vector<float> >* electronscbclusterenergy;

std::vector<std::vector<float> >* electronscrechitenergy;
std::vector<std::vector<float> >* electronscrechitfraction;
 std::vector<std::vector<float> >* electronscrechitlaserCorr;
 std::vector<std::vector<short> >* electronscrechitieta;
 std::vector<std::vector<short> >* electronscrechitiphi;

 
 int electronscnhits[nElectronMAX]; 
 float electroneLeft[nElectronMAX];
 float electroneRight[nElectronMAX];
 float electroneBottom[nElectronMAX];
 float electroneTop[nElectronMAX];

float electrone1x3[nElectronMAX];
float electrone3x1[nElectronMAX];
float electrone2x2[nElectronMAX];
float electrone3x2[nElectronMAX];
float electrone4x4[nElectronMAX];
float electrone2x5Right[nElectronMAX];
float electrone2x5Left[nElectronMAX];
float electrone2x5Top[nElectronMAX];
float electrone2x5Bottom[nElectronMAX];
///float electrone2x5Max[nElectronMAX];
///float electronenergyBasketFractionEta[nElectronMAX];
//float electronenergyBasketFractionPhi[nElectronMAX];
float electronlat[nElectronMAX][3];
float electronCovEtaEta[nElectronMAX];
float electronCovEtaPhi[nElectronMAX];
float electronCovPhiPhi[nElectronMAX];

float electronCovIEtaIEta[nElectronMAX];
float electronCovIEtaIPhi[nElectronMAX];
float electronCovIPhiIPhi[nElectronMAX];
float electronscCovIEtaIEta[nElectronMAX];
float electronscCovIEtaIPhi[nElectronMAX];
float electronscCovIPhiIPhi[nElectronMAX];
float electronzernike20[nElectronMAX];
float electronzernike42[nElectronMAX];



float electronscenergy[nElectronMAX];
int electronfiduficalFlag[nElectronMAX];
int electrontrackerDrivenSeed[nElectronMAX];
int electronecalDrivenSeed[nElectronMAX];
float electronfbrem[nElectronMAX];
int electronnumberOfBrems[nElectronMAX];
float electrondr03TkSumPt[nElectronMAX];
float electrondr03EcalRecHitSumEt[nElectronMAX];
float electrondr03HcalDepth1TowerSumEt[nElectronMAX];
float electrondr03HcalDepth2TowerSumEt[nElectronMAX];
float electrondr03HcalTowerSumEt[nElectronMAX];
float electrondr04EcalRecHitSumEt[nElectronMAX];
float electrondr04HcalDepth1TowerSumEt[nElectronMAX];
float electrondr04HcalDepth2TowerSumEt[nElectronMAX];
float electrondr04HcalTowerSumEt[nElectronMAX];
float electrondr04TkSumPt[nElectronMAX];
float electronhcalDepth1OverEcal[nElectronMAX];
float electronhcalDepth2OverEcal[nElectronMAX];
float electronhcalOverEcal[nElectronMAX];
float electroneSuperClusterOverP[nElectronMAX];
float electroneSeedClusterOverP[nElectronMAX];
float electroneSeedClusterOverPout[nElectronMAX];
float electroneEleClusterOverPout[nElectronMAX];
float electrondeltaEtaSuperClusterTrackAtVtx[nElectronMAX];
float electrondeltaEtaSeedClusterTrackAtCalo[nElectronMAX];
float electrondeltaEtaEleClusterTrackAtCalo[nElectronMAX];
float electrondeltaPhiSuperClusterTrackAtVtx[nElectronMAX];
float electrondeltaPhiSeedClusterTrackAtCalo[nElectronMAX];
float electrondeltaPhiEleClusterTrackAtCalo[nElectronMAX];
int electronclassification[nElectronMAX];
float electronmva[nElectronMAX];
int electronnumberOfTracks[nElectronMAX];
float electronconvDist[nElectronMAX];
float electronconvDcot[nElectronMAX];
float electronconvRadius[nElectronMAX];
int electronconvFlags[nElectronMAX];

// int electronExpectednumberOfHits[nElectronMAX];
 int electronExpectedHitsInnernumberOfHits[nElectronMAX];
 int electronExpectedHitsOuternumberOfHits[nElectronMAX];
 
float electrongsfTrackvx[nElectronMAX];
float electrongsfTrackvy[nElectronMAX];
float electrongsfTrackvz[nElectronMAX];
float electrongsfTracknormalizedChi2[nElectronMAX];
float electrongsfTrackdxybeamSpot[nElectronMAX];
float electrongsfTrackdzbeamSpot[nElectronMAX];
int electrongsfTracknumberOfValidHits[nElectronMAX];
int electrongsfTracknumberOfLostHits[nElectronMAX];
int electrongsfTracknumberOfValidPixelHits[nElectronMAX];
int electrongsfTracknumberOfValidTrackerHits[nElectronMAX];
float electrongsfTrackpt[nElectronMAX];
float electrongsfTracketa[nElectronMAX];
float electrongsfTrackphi[nElectronMAX];

float electronseedtime[nElectronMAX];
float electronseedoutOfTimeChi2[nElectronMAX];
float electronseedchi2[nElectronMAX];
int electronseedrecoFlag[nElectronMAX];
int electronseedseverityLevel[nElectronMAX];
float electronseedlaserCorr[nElectronMAX];
  

float electrone1x5[nElectronMAX];
 float electrone2x5Max[nElectronMAX];
 float electrone5x5[nElectronMAX];
 float electronsigmaIetaIeta[nElectronMAX];
 float electrone3x3[nElectronMAX];
float electroneMax[nElectronMAX];
 int electronscindex[nElectronMAX];
  float electronswissCross[nElectronMAX];
  int electronieta[nElectronMAX];
  int electroniphi[nElectronMAX];
  float electronE2overE9[nElectronMAX];
  float electrongenelematch[nElectronMAX][4];
  float electrongentrkmatch[nElectronMAX][4];
  float electronsimelematch[nElectronMAX][7];
  float electronPartonmatch[nElectronMAX][6];


  static const int nMuonMAX = 50;
  ///muon 
  int nMuon; 
  TClonesArray *muonp4;  
  int  muonRecoAlgo[nMuonMAX];
  int muonisGoodMuon[nMuonMAX];
  
  float muonpt[nMuonMAX];
  int muoncharge[nMuonMAX];
  float muoneta[nMuonMAX];
  float muonphi[nMuonMAX];
  
  float muonvx[nMuonMAX];
  float muonvy[nMuonMAX];
  float muonvz[nMuonMAX];
  int muonnumberOfMatches[nMuonMAX];
  int muonnumberOfChambers[nMuonMAX];
  float muonglobalTrackvx[nMuonMAX];
  float muonglobalTrackvy[nMuonMAX];
  float muonglobalTrackvz[nMuonMAX];
  float muonglobalTrackdzbeamSpot[nMuonMAX];
  float muonglobalTrackdxybeamSpot[nMuonMAX];
  int muonglobalTracknumberOfValidPixelHits[nMuonMAX];
  int muonglobalTracknumberOfValidTrackerHits[nMuonMAX];
  float muonglobalTracknormalizedChi2[nMuonMAX];
  int muonglobalTracknumberOfValidMuonHits[nMuonMAX];
  float muonisolationR03sumPt[nMuonMAX];
  float muonisolationR03emEt[nMuonMAX];
  float muonisolationR03hadEt[nMuonMAX];
  float muonisolationR03hoEt[nMuonMAX];
  float muonisolationR04sumPt[nMuonMAX];
  float muonisolationR04emEt[nMuonMAX];
  float muonisolationR04hadEt[nMuonMAX];
  float muonisolationR04hoEt[nMuonMAX];
  float muonisolationR05sumPt[nMuonMAX];
  float muonisolationR05emEt[nMuonMAX];
  float muonisolationR05hadEt[nMuonMAX];
  float muonisolationR05hoEt[nMuonMAX];
  float muonGenmumatch[nMuonMAX][4];
  float muonGentrkmatch[nMuonMAX][4];
  float muonSimmumatch[nMuonMAX][7];
  float muonPartonmatch[nMuonMAX][6];
  
  int muoninnerTracknumberOfValidPixelHits[nMuonMAX];
  int muoninnerTracknumberOfValidTrackerHits[nMuonMAX];
  float muoninnerTrackpt[nMuonMAX];
  float muoninnerTracketa[nMuonMAX];
  float muoninnerTrackphi[nMuonMAX];
  
  int muonouterTracknumberOfValidMuonHits[nMuonMAX];
  float muonouterTrackpt[nMuonMAX];
  float muonouterTracketa[nMuonMAX];
  float muonouterTrackphi[nMuonMAX];
  
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
  


  
  
  int nVertexNoBS; 
  ///  static const int MAXVX = 100;
  float vertexNoBSx[MAXVX];
  float vertexNoBSy[MAXVX];
  float vertexNoBSz[MAXVX];
  float vertexNoBSchi2[MAXVX];
  float vertexNoBSndof[MAXVX];
  float vertexNoBSnormalizedChi2[MAXVX];
  int vertexNoBStrackSize[MAXVX];
  int vertexNoBSisFake[MAXVX];
  int vertexNoBSisValid[MAXVX];


  int phyDeclared; 
  


  float rho; 
  float rhoEtaMax44; 
  
  
  static const int nTrackMAX = 50000;
  //track 
  int nTrack; 
  float    trackhighPurityFraction; 
  float    trackd0[nTrackMAX]; 
  float    trackd0Error[nTrackMAX]; 
  float    trackdz[nTrackMAX]; 
  float    trackdzError[nTrackMAX]; 
  int    trackcharge[nTrackMAX]; 
  float    trackpt[nTrackMAX]; 
  float    trackpx[nTrackMAX]; 
  float    trackpy[nTrackMAX]; 
  float    trackpz[nTrackMAX]; 
  float    trackptError[nTrackMAX]; 
  float    tracketa[nTrackMAX]; 
  float trackphi[nTrackMAX]; 
  float    trackvx[nTrackMAX]; 
  float    trackvy[nTrackMAX]; 
  float    trackvz[nTrackMAX]; 
  float    tracknormalizedChi2[nTrackMAX]; 
  int    tracknumberOfValidHits[nTrackMAX]; 
  int trackndof[nTrackMAX]; 
  int    tracknumberOfValidPixelHits[nTrackMAX]; 
  int    tracknumberOfValidStripHits[nTrackMAX]; 
  int    tracknumberOfValidTrackerHits[nTrackMAX]; 
  int    trackalgo[nTrackMAX]; 
  int    trackqualityFlagTracks[nTrackMAX]; 
    

  ///beamspont
  float beamSpotX; 
  float beamSpotY; 
  float beamSpotZ; 
  
  
  ////MET
  float caloMETet; 
  float caloMETsumEt;
  float caloMETphi; 
  float caloMETeta; 
  float caloMETsig; 
            
      
  float muCorrMETsumEt;
  float muCorrMETet; 
  float muCorrMETphi; 
  float muCorrMETeta; 
  float muCorrMETsig; 
      
      
  float tcMETsumEt;
  float tcMETet; 
  float tcMETphi; 
  float tcMETeta; 
  float tcMETsig; 
      
      
  float pfMETsumEt;
  float pfMETet; 
  float pfMETphi; 
  float pfMETeta; 
  float pfMETsig; 

  

  std::vector<std::vector<short> >* vertex_trkind;
  std::vector<std::vector<float> >* vertex_trkWeight; 
  
  



  static const int nCJetMAX = 50;
  int nak5CaloJet;
  
  float ak5CaloJetet[nCJetMAX];
  float ak5CaloJeteta[nCJetMAX];
  float ak5CaloJetphi[nCJetMAX];
  float ak5CaloJetenergyFractionHadronic[nCJetMAX];
  float ak5CaloJetn90[nCJetMAX];
  float ak5CaloJetn60[nCJetMAX];
  float ak5CaloJettowersArea[nCJetMAX];
  int ak5CaloJetCalosize[nCJetMAX];
  float ak5CaloJetfHPD[nCJetMAX];
  float ak5CaloJetfRBX[nCJetMAX];
  float ak5CaloJetn90Hits[nCJetMAX];
  float ak5CaloJetrestrictedEMF[nCJetMAX];
  float ak5CaloJetcorrection[nCJetMAX];
  float ak5CaloJetPartonmatch[nCJetMAX][6];
  float ak5CaloJetmaxEInEmTowers[nCJetMAX];
  float ak5CaloJetmaxEInHadTowers[nCJetMAX];

  
  int nak7CaloJet;
  
  float ak7CaloJetet[nCJetMAX];
  float ak7CaloJeteta[nCJetMAX];
  float ak7CaloJetphi[nCJetMAX];
  float ak7CaloJetenergyFractionHadronic[nCJetMAX];
  float ak7CaloJetn90[nCJetMAX];
  float ak7CaloJetn60[nCJetMAX];
  float ak7CaloJettowersArea[nCJetMAX];
  int ak7CaloJetCalosize[nCJetMAX];
  float ak7CaloJetfHPD[nCJetMAX];
  float ak7CaloJetfRBX[nCJetMAX];
  float ak7CaloJetn90Hits[nCJetMAX];
  float ak7CaloJetrestrictedEMF[nCJetMAX];
  float ak7CaloJetcorrection[nCJetMAX];
  float ak7CaloJetPartonmatch[nCJetMAX][6];
  float ak7CaloJetmaxEInEmTowers[nCJetMAX];
  float ak7CaloJetmaxEInHadTowers[nCJetMAX];

///  reco::helper::JetIDHelper jetIDHelper;

  int nak5PFJet;
  
  float ak5PFJetet[nCJetMAX];
  float ak5PFJeteta[nCJetMAX];
  float ak5PFJetphi[nCJetMAX];
  
  float ak5PFJetchargedHadronEnergyFraction[nCJetMAX];
  float ak5PFJetchargedMuEnergyFraction[nCJetMAX];
  float ak5PFJetneutralEmEnergyFraction[nCJetMAX];
  int ak5PFJetchargedMultiplicity[nCJetMAX];
  int ak5PFJetneutralMultiplicity[nCJetMAX];
  int ak5PFJetPFsize[nCJetMAX];
  float ak5PFJetfHPD[nCJetMAX];
  float ak5PFJetfRBX[nCJetMAX];
  float ak5PFJetn90Hits[nCJetMAX];
  float ak5PFJetrestrictedEMF[nCJetMAX];
  float ak5PFJetcorrection[nCJetMAX];
  float ak5PFJetPartonmatch[nCJetMAX][6];
  

  int nak7PFJet;
  
  float ak7PFJetet[nCJetMAX];
  float ak7PFJeteta[nCJetMAX];
  float ak7PFJetphi[nCJetMAX];
  
  float ak7PFJetchargedHadronEnergyFraction[nCJetMAX];
  float ak7PFJetchargedMuEnergyFraction[nCJetMAX];
  float ak7PFJetneutralEmEnergyFraction[nCJetMAX];
  int ak7PFJetchargedMultiplicity[nCJetMAX];
  int ak7PFJetneutralMultiplicity[nCJetMAX];
  int ak7PFJetPFsize[nCJetMAX];
  float ak7PFJetfHPD[nCJetMAX];
  float ak7PFJetfRBX[nCJetMAX];
  float ak7PFJetn90Hits[nCJetMAX];
  float ak7PFJetrestrictedEMF[nCJetMAX];
  float ak7PFJetcorrection[nCJetMAX];
  float ak7PFJetPartonmatch[nCJetMAX][6];
  ///below only for MC
  
  vector<short> *pileupBunchX; 
  vector<short> *pileupNInteraction; 
  
  
  //PDF
  int pdfidfirst; 
  int pdfidsecond; 
  float pdfxfirst; 
  float pdfxsecond; 
  float xPDFfirst; 
  float xPDFsecond; 
  float pdfscalePDF; 

  ///signalID qscale
  int signalProcessID; 
  float qScale; 
  
  
  
  static const int MAXGenSaved = 1000;
  //gen-leve phton
  int nGenPht;
  float etaGenPht[MAXGenSaved];
  float phiGenPht[MAXGenSaved];
  float ptGenPht[MAXGenSaved];
  int pidmomGenPht[MAXGenSaved];
  int pidmom2GenPht[MAXGenSaved];
  int pidmom3GenPht[MAXGenSaved];
  int statusGenPht[MAXGenSaved];
  float vxGenPht[MAXGenSaved];
  float vyGenPht[MAXGenSaved];
  float vzGenPht[MAXGenSaved];
  
  //gen Muon 
  int nGenMu; 
  float etaGenMu[MAXGenSaved];
  float phiGenMu[MAXGenSaved];
  float ptGenMu[MAXGenSaved];
  int pidmomGenMu[MAXGenSaved];
  int pidmom2GenMu[MAXGenSaved];
  int pidmom3GenMu[MAXGenSaved];
  int statusGenMu[MAXGenSaved];
  int pidGenMu[MAXGenSaved];
  float vxGenMu[MAXGenSaved];
  float vyGenMu[MAXGenSaved];
  float vzGenMu[MAXGenSaved];
  
  //gen electron
  int nGenEle; 
  float etaGenEle[MAXGenSaved];
  float phiGenEle[MAXGenSaved];
  float ptGenEle[MAXGenSaved];
  int pidmomGenEle[MAXGenSaved];
  int pidmom2GenEle[MAXGenSaved];
  int pidmom3GenEle[MAXGenSaved];
  int statusGenEle[MAXGenSaved];
  int pidGenEle[MAXGenSaved];
  float vxGenEle[MAXGenSaved];
  float vyGenEle[MAXGenSaved];
  float vzGenEle[MAXGenSaved];
  
  
  
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
  int pidFDauMCpart[MAXMC]; ///pid of first daughter, mainly for checking if is 91 or 92. 
  
  
  
  static const int MAXDAU = 100;
  int barcodeDauMCpart[MAXMC][MAXDAU];
  int statusMCpart[MAXMC];
  float vtxXMCpart[MAXMC];
  float vtxYMCpart[MAXMC];
  float vtxZMCpart[MAXMC];
  int convPhtMCpart[MAXMC][3]; //only for status==1 Gen photon, and the two index to the two electrons 
  
  vector<int> partonList; 
  
  
   // L1 objects
      static const int MAXL1OBJ = 100;
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
      
      int nL1Muon;
      float L1Muon_e[MAXL1OBJ];
      float L1Muon_et[MAXL1OBJ];
      float L1Muon_pt[MAXL1OBJ];
      float L1Muon_eta[MAXL1OBJ];
      float L1Muon_phi[MAXL1OBJ];

      bool saveTrk_;
      bool saveLaserCorrSeedCrystal_;
   /*    edm::ESHandle<EcalLaserDbService> laser; */
      
      
/* }; */
