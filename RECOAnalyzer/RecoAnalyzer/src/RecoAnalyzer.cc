// -*- C++ -*-
//
// Package:    RecoAnalyzer
// Class:      RecoAnalyzer
// 
/**\class RecoAnalyzer RecoAnalyzer.cc RECOAnalyzer/RecoAnalyzer/src/RecoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yong Yang,32 3-A03,+41227679643,
//         Created:  Fri May 13 16:42:13 CEST 2011
// $Id$
//
//
 
#include "RECOAnalyzer/RecoAnalyzer/interface/RecoAnalyzer.h"

using namespace std;
using namespace edm;
using namespace reco;



class PtSorter {
public:
  template <class T> bool operator() ( const T& a, const T& b ) {
    return ( a.pt() > b.pt() );
  }
};

class EtSorter {
public:
  template <class T> bool operator() ( const T& a, const T& b ) {
    return ( a.et() > b.et() );
  }
};



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecoAnalyzer::RecoAnalyzer(const edm::ParameterSet& iConfig)

{
  
  cout<<"RecoAnalyzer initialzied.. " <<endl;
  
  //now do what ever initialization is needed
    

  //HLT 
  const edm::InputTag dtrigTag("TriggerResults","","HLT");
  trigResultTag_    = iConfig.getUntrackedParameter<edm::InputTag> ("trigResultTag",dtrigTag); 
  const edm::InputTag dtrigSumTag("hltTriggerSummaryAOD","","HLT");
  trigSummaryTag_ = iConfig.getUntrackedParameter<edm::InputTag>("trigSummaryTag",dtrigSumTag);
  
  //vertex
  m_vertexSrc =  iConfig.getUntrackedParameter<edm::InputTag>("vertex",edm::InputTag("offlinePrimaryVerticesWithBS"));
  m_vertexSrc2 =  iConfig.getUntrackedParameter<edm::InputTag>("vertex2",edm::InputTag("offlinePrimaryVertices"));
  
  saveTrk_ = iConfig.getUntrackedParameter<bool>("saveAllTracks",true);
  saveLaserCorrSeedCrystal_ = iConfig.getUntrackedParameter<bool>("saveLaserCorrSeedCrystal",false);
  saveAllRecHitsSC_ = iConfig.getUntrackedParameter<bool>("saveAllRecHitsSC",false);
  
  m_beamSpot = iConfig.getUntrackedParameter<edm::InputTag> ("beamSpotInputTag",edm::InputTag("offlineBeamSpot")); 
  m_tracksSrc  = iConfig.getUntrackedParameter<edm::InputTag>("tracks",edm::InputTag("generalTracks"));
  m_muonSrc    = iConfig.getUntrackedParameter<edm::InputTag>("muons",edm::InputTag("muons"));
  
  muIsoTrkMap_ = iConfig.getUntrackedParameter<edm::InputTag>("muIsoTrkMap", edm::InputTag("muIsoDepositTk")); 
  muIsoEcalMap_ = iConfig.getUntrackedParameter<edm::InputTag>("muIsoEcalMap", edm::InputTag("muIsoDepositCalByAssociatorTowers:ecal")); 
  muIsoHcalMap_ = iConfig.getUntrackedParameter<edm::InputTag>("muIsoHcalMap", edm::InputTag("muIsoDepositCalByAssociatorTowers:hcal")); 
  muIsoHoMap_ = iConfig.getUntrackedParameter<edm::InputTag>("muIsoHoMap", edm::InputTag("muIsoDepositCalByAssociatorTowers:ho")); 


  m_barrelSC = iConfig.getUntrackedParameter<edm::InputTag>("barrelSC",edm::InputTag("correctedHybridSuperClusters"));
  m_endcapSC = iConfig.getUntrackedParameter<edm::InputTag>("endcapSC",edm::InputTag("correctedMulti5x5SuperClustersWithPreshower"));
  
  m_ecalRecHitBarrel    = iConfig.getUntrackedParameter<edm::InputTag>("ecalRecHitBarrel",edm::InputTag("reducedEcalRecHitsEB"));
  m_ecalRecHitEndcap    = iConfig.getUntrackedParameter<edm::InputTag>("ecalRecHitEndcap",edm::InputTag("reducedEcalRecHitsEE"));


  m_photonSrc = iConfig.getUntrackedParameter<edm::InputTag>("photon",edm::InputTag("photons"));
  m_electronSrc = iConfig.getUntrackedParameter<edm::InputTag>("electron",edm::InputTag("gsfElectrons"));
  
  m_PileupSrc = iConfig.getUntrackedParameter<edm::InputTag>("pileup",edm::InputTag("addPileupInfo"));
  m_rhoCollection =  iConfig.getUntrackedParameter<edm::InputTag>("rhoCorrection",edm::InputTag("kt6PFJets","rho"));
  m_rhoCollection2 =  iConfig.getUntrackedParameter<edm::InputTag>("rhoCorrection2",edm::InputTag("kt6PFJets","rho"));
  
  m_allConversionsColl =  iConfig.getUntrackedParameter<edm::InputTag>("ConvertedPhotonColl",edm::InputTag("allConversions"));
  
  
  debug_ = iConfig.getParameter<int> ("debugLevel");
  
  hlt_pathName = new std::vector<std::string>; hlt_pathName->clear();
  hlt_candpath = new std::vector<std::vector<unsigned short> >; hlt_candpath->clear();
  hltfilter_candpath = new std::vector<std::vector<unsigned short> >; hltfilter_candpath->clear();
  hlt_bitFired = new std::vector<unsigned short>; hlt_bitFired->clear();
  hlt_label_names = new std::vector<std::string>; hlt_label_names->clear();
  
  m_l1GtRecordInputTag = iConfig.getUntrackedParameter<edm::InputTag>("L1GtRecordInputTag",edm::InputTag("gtDigis"));
  m_l1IsolatedTag = iConfig.getUntrackedParameter< edm::InputTag > ("l1IsolatedTag",edm::InputTag("l1extraParticles","Isolated"));
  m_l1NonIsolatedTag = iConfig.getUntrackedParameter< edm::InputTag > ("l1NonIsolatedTag",edm::InputTag("l1extraParticles","NonIsolated"));
  m_l1MuTag = iConfig.getUntrackedParameter< edm::InputTag > ("l1MuTag",edm::InputTag("l1extraParticles"));
  
  
  
  ///jetCorrectionService  = iConfig.getUntrackedParameter<std::string>  ("JetCorrectionService","ak5CaloL2L3"); 

  outputFile_   = iConfig.getParameter<std::string>("outputFile");
  rootFile_ = new TFile(outputFile_.c_str(),"RECREATE"); // open output file to store root-trees. 


  nEventsProcessed  =0;
  
  
}


RecoAnalyzer::~RecoAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
  
  runNumber      = iEvent.id().run();
  evtNumber    = iEvent.id().event();
  lumiBlock = iEvent.luminosityBlock();
  bunchX          = iEvent.bunchCrossing();
  orbitNumber       = iEvent.orbitNumber();
  const edm::Timestamp jtime = iEvent.time();
  evtTime = jtime.value() >> 32;
  isRealData = iEvent.isRealData();


  
  //if(saveLaserCorrSeedCrystal_ ){
  ///edm::ESHandle<EcalLaserDbService> laser;
  iSetup.get<EcalLaserDbRecord>().get(laser);
  

  if( !isRealData && nEventsProcessed ==0 ){
    Analysis->Branch("photongenphtmatch",photongenphtmatch,"photongenphtmatch[nPhoton][3]/F");
    Analysis->Branch("photongenelematch",photongenelematch,"photongenelematch[nPhoton][3]/F");
    Analysis->Branch("photongenphtconv",photongenphtconv,"photongenphtconv[nPhoton][4]/F");
    Analysis->Branch("photonPartonmatch",photonPartonmatch,"photonPartonmatch[nPhoton][6]/F");
    Analysis->Branch("photonPartonmatchp",photonPartonmatchp,"photonPartonmatchp[nPhoton][3]/F");
    
    Analysis->Branch("electrongenelematch",electrongenelematch,"electrongenelematch[nElectron][4]/F");
    Analysis->Branch("electrongentrkmatch",electrongentrkmatch,"electrongentrkmatch[nElectron][4]/F");
    Analysis->Branch("electronsimelematch",electronsimelematch,"electronsimelematch[nElectron][7]/F");
    Analysis->Branch("electronPartonmatch",electronPartonmatch,"electronPartonmatch[nElectron][6]/F");

    
    Analysis->Branch("muonGenmumatch",muonGenmumatch,"muonGenmumatch[nMuon][4]/F");
    Analysis->Branch("muonGentrkmatch",muonGentrkmatch,"muonGentrkmatch[nMuon][4]/F");
    Analysis->Branch("muonSimmumatch",muonSimmumatch,"muonSimmumatch[nMuon][7]/F");
    Analysis->Branch("muonPartonmatch",muonPartonmatch,"muonPartonmatch[nMuon][6]/F");


    Analysis->Branch("ak5CaloJetPartonmatch",ak5CaloJetPartonmatch,"ak5CaloJetPartonmatch[nak5CaloJet][6]/F");
    Analysis->Branch("ak7CaloJetPartonmatch",ak7CaloJetPartonmatch,"ak7CaloJetPartonmatch[nak7CaloJet][6]/F");
    Analysis->Branch("ak5PFJetPartonmatch",ak5PFJetPartonmatch,"ak5PFJetPartonmatch[nak5PFJet][6]/F");
    Analysis->Branch("ak7PFJetPartonmatch",ak7PFJetPartonmatch,"ak7PFJetPartonmatch[nak7PFJet][6]/F");
    
    cout<<"additional branch for MC matchinng defined " <<endl; 

    
  ///PDF
  Analysis->Branch("pdfidfirst",&pdfidfirst,"pdfidfirst/I");
  Analysis->Branch("pdfidsecond",&pdfidsecond,"pdfidsecond/I");
  Analysis->Branch("pdfxfirst",&pdfxfirst,"pdfxfirst/F");
  Analysis->Branch("pdfxsecond",&pdfxsecond,"pdfxsecond/F");
  Analysis->Branch("pdfscalePDF",&pdfscalePDF,"pdfscalePDF/F");
  Analysis->Branch("xPDFfirst",&xPDFfirst,"xPDFfirst/F");
  Analysis->Branch("xPDFsecond",&xPDFsecond,"xPDFsecond/F");
  
  ///signal procID 
  Analysis->Branch("signalProcessID",&signalProcessID,"signalProcessID/I");
  Analysis->Branch("qScale",&qScale,"qScale/F");
  
  
  ///gen electron, muon,photon
  Analysis->Branch("nGenPht",&nGenPht,"nGenPht/I");
  Analysis->Branch("etaGenPht",etaGenPht,"etaGenPht[nGenPht]/F");
  Analysis->Branch("phiGenPht",phiGenPht,"phiGenPht[nGenPht]/F");
  Analysis->Branch("ptGenPht",ptGenPht,"ptGenPht[nGenPht]/F");
  Analysis->Branch("vxGenPht",vxGenPht,"vxGenPht[nGenPht]/F");
  Analysis->Branch("vyGenPht",vyGenPht,"vyGenPht[nGenPht]/F");
  Analysis->Branch("vzGenPht",vzGenPht,"vzGenPht[nGenPht]/F");
  Analysis->Branch("pidmomGenPht",pidmomGenPht,"pidmomGenPht[nGenPht]/I");
  Analysis->Branch("pidmom2GenPht",pidmom2GenPht,"pidmom2GenPht[nGenPht]/I");
  Analysis->Branch("pidmom3GenPht",pidmom3GenPht,"pidmom3GenPht[nGenPht]/I");
  Analysis->Branch("statusGenPht",statusGenPht,"statusGenPht[nGenPht]/I");
  
  
  Analysis->Branch("nGenEle",&nGenEle,"nGenEle/I");
  Analysis->Branch("etaGenEle",etaGenEle,"etaGenEle[nGenEle]/F");
  Analysis->Branch("phiGenEle",phiGenEle,"phiGenEle[nGenEle]/F");
  Analysis->Branch("ptGenEle",ptGenEle,"ptGenEle[nGenEle]/F");
  Analysis->Branch("pidmomGenEle",pidmomGenEle,"pidmomGenEle[nGenEle]/I");
  Analysis->Branch("pidmom2GenEle",pidmom2GenEle,"pidmom2GenEle[nGenEle]/I");
  Analysis->Branch("pidmom3GenEle",pidmom3GenEle,"pidmom3GenEle[nGenEle]/I");
  Analysis->Branch("statusGenEle",statusGenEle,"statusGenEle[nGenEle]/I");
  Analysis->Branch("vxGenEle",vxGenEle,"vxGenEle[nGenEle]/F");
  Analysis->Branch("vyGenEle",vyGenEle,"vyGenEle[nGenEle]/F");
  Analysis->Branch("vzGenEle",vzGenEle,"vzGenEle[nGenEle]/F");
  
  Analysis->Branch("nGenMu",&nGenMu,"nGenMu/I");
  Analysis->Branch("etaGenMu",etaGenMu,"etaGenMu[nGenMu]/F");
  Analysis->Branch("phiGenMu",phiGenMu,"phiGenMu[nGenMu]/F");
  Analysis->Branch("ptGenMu",ptGenMu,"ptGenMu[nGenMu]/F");
  Analysis->Branch("pidmomGenMu",pidmomGenMu,"pidmomGenMu[nGenMu]/I");
  Analysis->Branch("pidmom2GenMu",pidmom2GenMu,"pidmom2GenMu[nGenMu]/I");
  Analysis->Branch("pidmom3GenMu",pidmom3GenMu,"pidmom3GenMu[nGenMu]/I");
  Analysis->Branch("statusGenMu",statusGenMu,"statusGenMu[nGenMu]/I");
  Analysis->Branch("vxGenMu",vxGenMu,"vxGenMu[nGenMu]/F");
  Analysis->Branch("vyGenMu",vyGenMu,"vyGenMu[nGenMu]/F");
  Analysis->Branch("vzGenMu",vzGenMu,"vzGenMu[nGenMu]/F");

    
    
    Analysis->Branch("genhiggsm",&genhiggsm,"genhiggsm/F");
    Analysis->Branch("genhiggspt",&genhiggspt,"genhiggspt/F");
    Analysis->Branch("genhiggseta",&genhiggseta,"genhiggseta/F");
    Analysis->Branch("genhiggsphi",&genhiggsphi,"genhiggsphi/F");
    Analysis->Branch("genhiggsvx",&genhiggsvx,"genhiggsvx/F");
    Analysis->Branch("genhiggsvy",&genhiggsvy,"genhiggsvy/F");
    Analysis->Branch("genhiggsvz",&genhiggsvz,"genhiggsvz/F");
    
    ///Analysis->Branch("genhiggsstatus",&genhiggsstatus,"genhiggsstatus/I");
    
  }

  
  nEventsProcessed ++; 
  //if( evtNumber != 30470 && evtNumber != 47710) return; 
  

  /////////////////////////// ================== HLT ===================== ////////////////////
  
  edm::Handle<trigger::TriggerEvent> triggerObj;
  edm::Handle<trigger::TriggerEventWithRefs> triggerObjWithRef;
  
  
  iEvent.getByLabel(trigSummaryTag_, triggerObj);
  hlt_bitFired->clear();
  bool changed = false;
  configProvider.init(iEvent.getRun(),iSetup,trigResultTag_.process(),changed);
  edm::Handle<edm::TriggerResults> h_triggerResults_HLT1;
  iEvent.getByLabel(trigResultTag_, h_triggerResults_HLT1);
  hlt_pathName->clear();
  if (h_triggerResults_HLT1.isValid()) {
    if(debug_ > 9) {
      std::cout << "Fill names HLT1" << std::endl;
    }
    for (size_t i = 0; i < configProvider.size(); ++i){
      hlt_pathName->push_back(configProvider.triggerName(i));
    }
    
    // Trigger Results
    if(debug_ > 99){
      std::cout << "### Trigger Results 1 :" << trigResultTag_.process() << std::endl;
    }
    for (size_t i = 0; i < configProvider.size(); ++i) {
      if(debug_ > 99) {
	std::cout << i << "\t" << configProvider.triggerName(i) << " " << (h_triggerResults_HLT1->accept(i) ? "passed" : "failed") << std::endl;
      }
      
      if(h_triggerResults_HLT1->accept(i)){
        hlt_bitFired->push_back((unsigned short)(i));
      }
    }
  }
  if(!triggerObj.isValid()) 
    throw(cms::Exception("Release Validation Error") << "RAW-type HLT results not found" );
  
  // This can be improved doing it only when necessary
  theHLTLabels.clear();
  hlt_label_names->clear();
  for(int i=0; i<triggerObj->sizeFilters(); ++i) {
    if(debug_ > 9) {
      std::cout << "trigojbec fitlerTag "<<  triggerObj->filterTag(i) << " filterTagLabel "<< triggerObj->filterTag(i).label() <<endl; 
    }
    theHLTLabels.push_back( triggerObj->filterTag(i));
    
    ////to save the lable names of all fitlers
    hlt_label_names->push_back( triggerObj->filterTag(i).label());
    
  }
  
  hlt_pt->clear();
  hlt_eta->clear();
  hlt_phi->clear();
  hlt_en->clear();
  
  ///  hlt_p4->Clear();
  hlt_candpath->clear();
  hltfilter_candpath->clear();
  
  
  trigger::TriggerObjectCollection triggerObjs = triggerObj->getObjects();
  if(debug_ > 99) {
    std::cout << "Trigger Objects found " << triggerObjs.size() << std::endl;
  }
  
  ///loop all all trigger objects and check if it fires the HLT, if yes, save the list of trigger bits 
  for (unsigned int iCand=0; iCand<triggerObjs.size(); ++iCand ) { 

  //   if (hlt_n >= MAX_HLT) {
//       std::cout << "GlobeHLT: WARNING TOO MANY HLT CANDIDATES:  " << hlt_n << " found (allowed " << MAX_HLT << ")" << std::endl;
//       break;
//     }
    
    std::vector<unsigned short> temp; /// to check all HLT_Path which has been fired by this trigobject ( icand) 
        
    std::vector<unsigned short> temp1;  ///to check all HLT_Filter which has been fired by this trigojbec ( icand) 
    
    trigger::TriggerObject object = triggerObjs[iCand];

    for(unsigned int n=0; n<theHLTLabels.size(); n++) { ///loop over all filters 

      trigger::size_type index = triggerObj->filterIndex(theHLTLabels[n]);
      
      // Check HLT( atcually filter) 
      bool firedHLT = false;
      if (!(index >= triggerObj->sizeFilters())) {
        const trigger::Keys & k = triggerObj->filterKeys(index);
        for (trigger::Keys::const_iterator ki = k.begin(); ki !=k.end(); ++ki ) {
          if(*ki == iCand)
            firedHLT = true;
        }
      }
      
      
      if(firedHLT) {
	temp1.push_back(n);
	
        if(debug_ > 99) {
	  std::cout << n << "\t" << theHLTLabels[n].label() << "\t fired" << std::endl;
	}
      
        for (unsigned int i=0; i<configProvider.size(); i++) {  ///check all HLT_Path 
          
          unsigned int nModules = configProvider.moduleLabels(i).size();
          unsigned int moduleIndex = configProvider.moduleIndex(i, theHLTLabels[n].label());
          //std::cout << nModules << " " << moduleIndex << std::endl;
          //for(unsigned int y=0; y<nModules; y++) {
          //  std::cout << y << " " << configProvider.moduleLabel(i, y) << " " << theHLTLabels[n].label() << " " << configProvider.triggerName(i) << std::endl;
          //  if (configProvider.moduleLabel(i, y) == theHLTLabels[n].label()) 
          //    std::cout << "MATCH !!!!"  << " " << y << " " << nModules << std::endl;
          // }
          //std::cout << "------------" << std::endl;
          if ((nModules - moduleIndex) == 2) {
            //std::cout <<  theHLTLabels[n].label() << " " << std::endl;;
            //std::cout << configProvider.moduleLabel(i, moduleIndex) << " " << nModules << " " << moduleIndex << std::endl;
            //std::cout <<  configProvider.triggerName(i) << std::endl;
            //std::cout << std::endl;
            temp.push_back(i);
          }
        }
	
      }
      
    } //end of loop over all hlt filters
    
    
    if(temp1.size() !=0 ){ ///if this object pass more than one filter
      hltfilter_candpath->push_back(temp1);  /// what stored is <0,2 > , <2,3> ....  each number refer to the index inside the vector of hlt_label_names-
      
      if(temp.size() ==0){
	temp.push_back(-1); ///dummy value in case of 0 
      }
      hlt_candpath->push_back(temp);  ///what stored is  < 0,2 > , <-1>,  , ... each number refer to the index inside the vector of hlt_pathName
      hlt_pt->push_back(object.pt());
      hlt_eta->push_back(object.eta());
      hlt_phi->push_back(object.phi());
      hlt_en->push_back(object.energy());
    }
    
    //  // Skip if no triggers were fired
    //     if(temp.size() != 0) {
    //       hlt_candpath->push_back(temp);
    
    //       // Set HLT candidate p4
    //       //TLorentzVector lv(object.px(), object.py(), object.pz(), 0);
    //       //new ((*hlt_p4)[hlt_n]) TLorentzVector();
    //       //((TLorentzVector *)hlt_p4->At(hlt_n))->SetXYZT(object.px(), object.py(), object.pz(), object.energy());
    
    
    //       hlt_n++;          
    //     } // Store candidate which fired at least 1 HLT
    
  } // TriggerCandidate's Loop
  
  if(debug_ > 99) std::cout << "Trigger Objects stored " << hlt_pt->size() << std::endl;
  
  
  /////////////////////////// ================== END OF HLT ===================== ////////////////////
  
  
  //PILE-UP
  if( !isRealData){
    pileupBunchX->clear();
    pileupNInteraction->clear();
    
    Handle<std::vector< PileupSummaryInfo > >  PupInfo;
    iEvent.getByLabel(m_PileupSrc, PupInfo);
    std::vector<PileupSummaryInfo>::const_iterator PVI;
    for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
      if(debug_> 1) std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
      pileupBunchX->push_back(PVI->getBunchCrossing());
      pileupNInteraction->push_back(PVI->getPU_NumInteractions());
    }
  }
    
  
  
  edm::Handle<double> rhoHandle;
  rho = 0; 
  try{
    iEvent.getByLabel(m_rhoCollection, rhoHandle);
    rho = *(rhoHandle.product()); ///crash here ( re-run is needed) 
    //cout<<*(rhoHandle.product()) <<endl; 
    
  }catch(std::exception& ex ){
    cout<<"rhoCollection not working.."<<endl;
  }
  
  
  edm::Handle<double> rhoHandle2;
  rhoEtaMax44 = 0; 
  try{
    iEvent.getByLabel(m_rhoCollection2, rhoHandle2);
    rhoEtaMax44 = *(rhoHandle2.product()); ///crash here ( re-run is needed) 
    //cout<<*(rhoHandle.product()) <<endl; 
    
  }catch(std::exception& ex ){
    cout<<"rhoCollection2 not working.."<<endl;
  }
  
  /////////////// ============== Generator information ================= ////////////////
  partonList.clear();
  if( !isRealData){
    
    edm::Handle<GenEventInfoProduct> geninfos;
    try{
      iEvent.getByLabel("generator",geninfos);
    
      pdfidfirst= geninfos->pdf()->id.first; 
      pdfidsecond= geninfos->pdf()->id.second;
      pdfxfirst= geninfos->pdf()->x.first; 
      pdfxsecond= geninfos->pdf()->x.second; 
      xPDFfirst= geninfos->pdf()->xPDF.first;
      xPDFsecond= geninfos->pdf()->xPDF.second; 
      pdfscalePDF= geninfos->pdf()->scalePDF; 
      signalProcessID= int( geninfos->signalProcessID());
      //ptHAT = geninfos->binningValues()[0]; //this is ptHAT, same as geninfos->qScale() 
      qScale= geninfos->qScale();
      
      if(debug_>=2) cout<<"ptHat: "<< qScale <<" "<< signalProcessID <<" "<<endl; 
      
    }catch(std::exception& ex ){
      cout<<"geninfos not working.."<<endl;
    }    
    
    genhiggsm = 0; 
    genhiggspt = 0; 
    genhiggseta = 0; 
    genhiggsphi = 0;
    genhiggsstatus = 0; 
    genhiggsvx =0;
    genhiggsvy =0;
    genhiggsvz =0;
    
    Handle<GenParticleCollection> genPart;
    try {
      iEvent.getByLabel( "genParticles", genPart);
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
	
	pidFDauMCpart[nMCpart] = 0; 
	if(nDauMCpart[nMCpart]>0){ ///for gluon, this is mostly 92
	  pidFDauMCpart[nMCpart] = p.daughter(0)->pdgId(); 
	}
	
	if( p.pdgId() == 25  && p.status() ==2 ){ //two status 3 and 2 and the one with 2 is saved
	  genhiggsm = p.mass();
	  genhiggspt = p.pt();
	  genhiggseta = p.eta();
	  genhiggsphi = p.phi();
	  genhiggsstatus = p.status();
	  genhiggsvx = p.vx();
	  genhiggsvy = p.vy();
	  genhiggsvz = p.vz();
	  if(debug_>0) cout<<"mcfoundhiggs " << runNumber <<" "<<  evtNumber <<" "<<  i <<" "<< p.status()<<" "<< p.pdgId()<<" "<< p.pt()<<endl; 
	  
	}
	
	// 	if( evtNumber == 30470 || evtNumber == 47710){
	// 	  if( p.pdgId() == 22 || p.pdgId() == 25){
	// 	    cout<<"chekcrunhiggs " << runNumber <<" "<<  evtNumber <<" "<<  i <<" "<< p.pdgId()<<" "<< p.pt()<<endl; 
	// 	  }
	// 	}
	

	
	if( nMCpart >= 6){// Don't take into account first 6 particles in generator list, this list is used for quark/gluon matching
	  int pid = abs(p.pdgId()); 
	  int status = p.status();
	  if(  (( pid >=1 && pid <=5) || pid==21) && status ==3  ){
	    partonList.push_back(nMCpart);
	  }
	  else if( pidFDauMCpart[nMCpart] == 91 || pidFDauMCpart[nMCpart] ==92 ){
	    partonList.push_back(nMCpart);
	  }
	}	
	
	nMCpart++;
      }
    }catch(std::exception& ex ){
      cout<<"genParticles not working.."<<endl;
    }   
    
    
    //find mother index
    for(int i=0; i< int(genPart->size())&& i <MAXMC; i++){
      const Candidate & p = (*genPart)[i];
      //int pid = p.pdgId();
      const Candidate * mom = p.mother();
      if( mom !=NULL){
	pidmomMCpart[i] = mom->pdgId();
	barcodemomMCpart[i] = indexofParticle(mom->px(),mom->pz(),mom->status());
      }
    }
    
    
//     if( evtNumber == 30470 || evtNumber == 47710){
//       for(int j=0; j< nMCpart; j++){
// 	const GenParticle & p = (*genPart)[j];
// 	int indmom = getMotherIndex(j);
// 	if(indmom>=0){
// 	  if( pidMCpart[indmom] ==25){
// 	    cout<<"foundhiggsdecay " << j<<" "<< p.pdgId()<<" "<< p.pt()<<endl; 
// 	  }
// 	}
//       }
    
//     }
    
    //Save all gen-level electron muon, photon (with pt > 1) 
    nGenPht = 0; 
    nGenMu = 0; 
    nGenEle = 0; 
    for(int i=0; i< int(genPart->size())&& i <MAXMC; i++){
      const Candidate & p = (*genPart)[i];
	  	  
      if( ( abs(p.pdgId()) == 11 ) && nGenMu < MAXGenSaved ){
	etaGenMu[nGenMu] = p.eta();
	phiGenMu[nGenMu] = p.phi();
	ptGenMu[nGenMu] = p.pt();
	pidmomGenMu[nGenMu] = 0; 
	pidmom2GenMu[nGenMu] = 0;
	statusGenMu[nGenMu] = p.status();
	pidGenMu[nGenMu] = p.pdgId();
	vxGenMu[nGenMu] = p.vx();
	vyGenMu[nGenMu] = p.vy();
	vzGenMu[nGenMu] = p.vz();
	
	const Candidate * mom = p.mother();
	if( mom !=NULL){
	  pidmomGenMu[nGenMu] = mom->pdgId();
	  const Candidate * mom2 = mom->mother();
	  if( mom2 !=NULL){
	    pidmom2GenMu[nGenMu] = mom2->pdgId();
	  }
	}
	int indmom = getMotherIndex(i);
	if(indmom>=0){
	  pidmom3GenMu[nGenMu] = pidMCpart[indmom];
	}
	nGenMu ++; 
      }
      
      if( ( abs(p.pdgId()) == 13 ) && nGenEle < MAXGenSaved ){
	etaGenEle[nGenEle] = p.eta();
	phiGenEle[nGenEle] = p.phi();
	ptGenEle[nGenEle] = p.pt();
	pidmomGenEle[nGenEle] = 0; 
	pidmom2GenEle[nGenEle] = 0;
	statusGenEle[nGenEle] = p.status();
	pidGenEle[nGenEle] = p.pdgId();
	vxGenEle[nGenEle] = p.vx();
	vyGenEle[nGenEle] = p.vy();
	vzGenEle[nGenEle] = p.vz();
	const Candidate * mom = p.mother();
	if( mom !=NULL){
	  pidmomGenEle[nGenEle] = mom->pdgId();
	  const Candidate * mom2 = mom->mother();
	  if( mom2 !=NULL){
	    pidmom2GenEle[nGenEle] = mom2->pdgId();
	  }
	}
	int indmom = getMotherIndex(i);
	if(indmom>=0){
	  pidmom3GenEle[nGenEle] = pidMCpart[indmom];   ///this is for cross check, trace back all untill find the mother ( not id itself). 
	}
	nGenEle ++; 
      }
      
      
      if(p.pdgId() == 22 ){
	
	bool tosave = false; 
	if( p.pt()>1) tosave = true; 



	const Candidate * mom = p.mother();
	if( mom !=NULL){
	  if( mom->pdgId() <100) tosave = true; 
	}
	if( tosave  && nGenPht < MAXGenSaved ) {
	  etaGenPht[nGenPht] = p.eta();
	  phiGenPht[nGenPht] = p.phi();
	  ptGenPht[nGenPht] = p.pt();
	  pidmomGenPht[nGenPht] = 0; 
	  pidmom2GenPht[nGenPht] = 0;
	  statusGenPht[nGenPht] = p.status();
	  vxGenPht[nGenPht] = p.vx();
	  vyGenPht[nGenPht] = p.vy();
	  vzGenPht[nGenPht] = p.vz();
	
	  const Candidate * mom = p.mother();
	  if( mom !=NULL){
	    pidmomGenPht[nGenPht] = mom->pdgId();
	    const Candidate * mom2 = mom->mother();
	    if( mom2 !=NULL){
	      pidmom2GenPht[nGenPht] = mom2->pdgId();
	    }
	  }
	  int indmom = getMotherIndex(i);
	  if(indmom>=0){
	    pidmom3GenPht[nGenPht] = pidMCpart[indmom];
	  }
	  nGenPht ++; 
	  
	}
	
      }
    }
    
    
    bool runsimTrack = false; //not in AOD

    if( runsimTrack ){
      // Simulated tracks (i.e. GEANT particles).
      Handle<edm::SimTrackContainer> psimtracks; 
      // Get the associated vertices
      Handle<SimVertexContainer> simvertices;
      try{
	iEvent.getByLabel("g4SimHits", simvertices);
      }catch( std::exception& ex){
	cout<<" SimVertexContainer not working.."<<endl;
      }
      
      try{
	iEvent.getByLabel("g4SimHits", psimtracks);
	
	// Need to check that SimTrackContainer is sorted; otherwise, copy and sort :-(
	std::auto_ptr<SimTrackContainer> simtracksTmp;
	const SimTrackContainer * simtracksSorted = &* psimtracks;
	
	if (!__gnu_cxx::is_sorted(psimtracks->begin(), psimtracks->end(), LessById())) {
	  simtracksTmp.reset(new SimTrackContainer(*psimtracks));
	  std::sort(simtracksTmp->begin(), simtracksTmp->end(), LessById());
	  simtracksSorted = &* simtracksTmp;
	}
	
      
	int nMCpartGen = nMCpart; 
	
	
	map<unsigned int,int> map_simtrackid; 
    
	// loop through all simParticles
	for (SimTrackContainer::const_iterator iM = psimtracks->begin(); 
	     iM != psimtracks->end(); ++iM) {
      
	  // Skip PYTHIA tracks.
	  if (iM->genpartIndex() != -1) continue; 
	
	
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
		if(abs(iM->type())==11 && pidMCpart[it->genpartIndex() -1] == 22) {
		  convPhtMCpart[it->genpartIndex() -1][0] = 1; ///this gen photon converts
		  if( vtxXMCpart[nMCpartGen + nSimTrack] != 0){
		    if(  convPhtMCpart[it->genpartIndex() -1][1] < 0) convPhtMCpart[it->genpartIndex() -1][1] = nMCpartGen + nSimTrack; ///the first electon of this gen photon converts
		    else if(  convPhtMCpart[it->genpartIndex() -1][2] < 0) convPhtMCpart[it->genpartIndex() -1][2] = nMCpartGen + nSimTrack; ///the second electon of this gen photon converts
		  }
		}
	      
		vector<int>::iterator iit = find(indpht.begin(),indpht.end(),it->genpartIndex() -1);
		if( iit == indpht.end()){
		  indpht.push_back(it->genpartIndex() -1); 
		}
		if( iM->type() == 11) inde1.push_back(nMCpartGen+nSimTrack); 
		else inde2.push_back(nMCpartGen+nSimTrack); 
	      }else{
		map<unsigned int,int>::iterator intt = map_simtrackid.find(idx);
		if( intt != map_simtrackid.end()){
		  barcodemomMCpart[nMCpartGen+nSimTrack]  = intt->second;
		}else{
		  barcodemomMCpart[nMCpartGen+nSimTrack] = -3; ////check this if any
		}
	      
	      }
	    } ///found parent in simtrack
	  
	  }
	
	  nSimTrack ++; 
	
	
	} ///end of 2nd loop 
      
      

      }catch( std::exception& ex){
	cout<<"SimTrackContainer  not working.."<<endl;
      }
    
    }
    

  }
  /////////// =================END OF Generator information ================ ////////////////////
    







  ////////////////////// L1 ojbects ////////////////

  nL1EMIso=0;
  nL1EMnonIso = 0;
  nL1Muon = 0;
  
  edm::Handle< l1extra::L1EmParticleCollection > emIsolColl ;
  try{
    iEvent.getByLabel(m_l1IsolatedTag, emIsolColl ) ;
    for( l1extra::L1EmParticleCollection::const_iterator emItr = emIsolColl->begin(); emItr != emIsolColl->end() && nL1EMIso < MAXL1OBJ ;++emItr){
      L1EMIso_e[nL1EMIso] = emItr->energy();
      L1EMIso_et[nL1EMIso] = emItr->et();
      L1EMIso_eta[nL1EMIso] = emItr->eta();
      L1EMIso_phi[nL1EMIso] = emItr->phi();
      nL1EMIso++;
    }
    
  }catch( std::exception& ex ) {
    cout<<"L1 em iso not working.."<<endl;
  }
  edm::Handle< l1extra::L1EmParticleCollection > emNonIsolColl ;
  
  try{
    iEvent.getByLabel(m_l1NonIsolatedTag, emNonIsolColl ) ;
    for( l1extra::L1EmParticleCollection::const_iterator emItr = emNonIsolColl->begin(); emItr != emNonIsolColl->end() && nL1EMnonIso <MAXL1OBJ  ;++emItr){
      L1EMnonIso_e[nL1EMnonIso] = emItr->energy();
      L1EMnonIso_et[nL1EMnonIso] = emItr->et();
      L1EMnonIso_eta[nL1EMnonIso] = emItr->eta();
      L1EMnonIso_phi[nL1EMnonIso] = emItr->phi();
      nL1EMnonIso++;
    }
  }catch( std::exception& ex ) {
    cout<<"L1 em noniso not working.."<<endl;
  }
  
  //Get the L1 Muon Collection
  edm::Handle< l1extra::L1MuonParticleCollection > muColl ;
  try{
    iEvent.getByLabel(m_l1MuTag, muColl ) ;
    for( l1extra::L1MuonParticleCollection::const_iterator emItr = muColl->begin(); emItr != muColl->end() && nL1Muon < MAXL1OBJ;++emItr){
      L1Muon_e[nL1Muon] = emItr->energy();
      L1Muon_pt[nL1Muon] = emItr->pt();
      L1Muon_eta[nL1Muon] = emItr->eta();
      L1Muon_phi[nL1Muon] = emItr->phi();
      nL1Muon++;
    }
  }catch( std::exception& ex ) {
    cout<<"L1 muon  not working.."<<endl;
  }
  
  


  //////========== BeamSpot ==============////////
  reco::TrackBase::Point beamPoint(0,0, 0);
  reco::BeamSpot beamSpot;
  edm::Handle<reco::BeamSpot> recoBeamSpotHandle;
  try{
    iEvent.getByLabel(m_beamSpot,recoBeamSpotHandle);
    beamSpot = *recoBeamSpotHandle;
    
    beamPoint = beamSpot.position();  //used in below as well 
    beamSpotX = beamSpot.position().X();
    beamSpotY = beamSpot.position().Y();
    beamSpotZ = beamSpot.position().Z();
    
    if(debug_ >=2) cout<<"beamPoint: "<< beamSpot.position().X() <<" "<< beamSpot.position().Y() <<" "<< beamSpot.position().Z() <<endl; 
    
  }catch(std::exception& ex ){
    cout<<"beamSpot not working.."<<endl;
  }
  //////========== End of BeamSpot ==============////////
   //////========== tracks ==============////////
  nTrack = 0; 
  int nHighPurityTrack = 0; 
  trackhighPurityFraction = -1; 
  reco::TrackCollection ptSortedtracks; 
  
  edm::Handle<TrackCollection> tracks;
  try{
    iEvent.getByLabel(m_tracksSrc, tracks);
    
    ptSortedtracks = *tracks; 
    std::stable_sort( ptSortedtracks.begin(), ptSortedtracks.end(), PtSorter());
    TrackCollection::const_iterator ittrk = ptSortedtracks.begin();
    nTrack = 0; 
    for(;  ittrk != ptSortedtracks.end() ; ++ittrk){
      if(ittrk->quality(reco::Track::highPurity)){
	nHighPurityTrack ++; 
      }

      if( nTrack < nTrackMAX){
	int qualityFlag = 0;
	if(ittrk->quality(reco::Track::undefQuality))
	  qualityFlag = qualityFlag | 1 << 0; 
	if(ittrk->quality(reco::Track::loose))
	  qualityFlag = qualityFlag | 1 << 1; 
	if(ittrk->quality(reco::Track::tight))
	  qualityFlag = qualityFlag | 1 << 2; 
	if(ittrk->quality(reco::Track::highPurity)){
	  qualityFlag = qualityFlag | 1 << 3; 
	}
	if(ittrk->quality(reco::Track::confirmed))
	  qualityFlag = qualityFlag | 1 << 4; 
	if(ittrk->quality(reco::Track::goodIterative))
	  qualityFlag = qualityFlag | 1 << 5; 
	if(ittrk->quality(reco::Track::qualitySize))
	  qualityFlag = qualityFlag | 1 << 6; 
      
	trackd0[nTrack] = ittrk->d0();
	trackd0Error[nTrack] = ittrk->d0Error();
	trackdz[nTrack] = ittrk->dz();
	trackdzError[nTrack] = ittrk->dzError();
	trackpt[nTrack]= ittrk->pt();
	trackcharge[nTrack]= ittrk->charge();
	trackpx[nTrack]= ittrk->px();
	trackpy[nTrack]= ittrk->py();
	trackpz[nTrack]= ittrk->pz();
	trackptError[nTrack]= ittrk->ptError();
	tracketa[nTrack]= ittrk->eta();
	trackphi[nTrack]= ittrk->phi();
	trackvx[nTrack]= ittrk->vx();
	trackvy[nTrack]= ittrk->vy();
	trackvz[nTrack]= ittrk->vz();
	tracknormalizedChi2[nTrack] = ittrk->normalizedChi2();
	tracknumberOfValidHits[nTrack] = ittrk->numberOfValidHits();
	trackndof[nTrack] = ittrk->ndof(); 
	///tracknumberOfValidHits[nTrack] = ittrk->numberOfValidHits(); 
	tracknumberOfValidPixelHits[nTrack] = ittrk->hitPattern().numberOfValidPixelHits(); 
	tracknumberOfValidStripHits[nTrack] = ittrk->hitPattern().numberOfValidStripHits(); 
	tracknumberOfValidTrackerHits[nTrack] = ittrk->hitPattern().numberOfValidTrackerHits();
	trackalgo[nTrack] = int(ittrk->algo()) ;
	trackqualityFlagTracks[nTrack] = qualityFlag; 
	
	
	nTrack ++; 

      }
      
      
    }
    if( int(tracks->size()) >0){
      trackhighPurityFraction = nHighPurityTrack *1.0 / tracks->size(); 
    }
    
    if(debug_ >=2) cout<<" trackhighPurityFraction "<< nHighPurityTrack <<" "<< nTrack<<" "<< ptSortedtracks.size()<<" "<< trackhighPurityFraction <<endl; 
    
    
  }catch(std::exception& ex ){
    cout<<"track not working.."<<endl;
  }
  //////========== End of tracks ==============////////
  
  //////========== vertex ==============////////
  
  try{
    
    nVertex = 0; 
    
    //vtx_tkind->clear();
    vertex_trkind->clear();
    vertex_trkWeight->clear();
        
    Handle<VertexCollection> vertices;
    iEvent.getByLabel(m_vertexSrc, vertices);
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
      
      
      double vertex_sumTrks = 0.0;
      double vertex_sumTrks2 = 0.0;
      std::vector<short> temp;
      std::vector<float> temp_float;
      //for (reco::Vertex::trackRef_iterator vertex_curTrack = (*vertices)[j].tracks_begin(); vertex_curTrack!=(*vertices)[j].tracks_end(); vertex_curTrack++) {
      for (std::vector<reco::TrackBaseRef>::const_iterator vertex_curTrack = (*vertices)[j].tracks_begin(); vertex_curTrack!=(*vertices)[j].tracks_end(); vertex_curTrack++) {
	vertex_sumTrks += (*vertex_curTrack)->pt();
	vertex_sumTrks2 += (*vertex_curTrack)->pt() * (*vertex_curTrack)->pt(); 
	
	int index = findIndexofTRK( (*vertex_curTrack)->pt(), (*vertex_curTrack)->eta(), (*vertex_curTrack)->phi());
	bool ismatched = index >=0 ; 
	if( ismatched ){
	  temp.push_back(index);
	  temp_float.push_back((*vertices)[j].trackWeight(*vertex_curTrack) );
	}
	// 	bool ismatched2 = false; 
// 	//for(reco::TrackCollection::size_type t = 0; t<ptSortedtracks.size() ; ++t) { ///doesn't work for ptSortedtrack 
// 	for(reco::TrackCollection::size_type t = 0; t<tracks->size(); ++t) {
//           reco::TrackRef track(tracks, t);
// 	  //reco::TrackRef track(&ptSortedtracks, t);
// 	  if (&(**vertex_curTrack) == &(*track)) {
// 	    ismatched2 = true; 
// 	    ///temp.push_back(index);
// 	    ///temp_float.push_back((*vertices)[j].trackWeight(track));
// 	    break; 
// 	  }
// 	}
// 	if( ismatched != ismatched2){
// 	  cout<<"check trk match " <<endl; 
// 	  exit(1);
// 	}
		
	if( ! ismatched ) {
          temp.push_back(-9999);
          temp_float.push_back(-9999);
        }
      }
      
      vertex_trkind->push_back(temp);
      vertex_trkWeight->push_back(temp_float);
      
      vertexsumtrackspt[nVertex] = vertex_sumTrks; 
      vertexsumtracksptSquare[nVertex] = vertex_sumTrks2; 
            
      ///sumPt of refitted Tracks   NA 
      if((*vertices)[j].hasRefittedTracks()){
	vector<Track> vTrack = (*vertices)[j].refittedTracks();
	double sumpt = 0; 
	double sumpt2 = 0; 
	for( vector<Track>::const_iterator itrk= vTrack.begin(); itrk != vTrack.end(); itrk++){
	  sumpt += itrk->pt();
	  sumpt2 += itrk->pt() * itrk->pt();
	}
	vertexsumRefittedtrackspt[nVertex] = sumpt; 
	vertexsumRefittedtracksptSquare[nVertex] = sumpt2; 
      
      }else{
	vertexsumRefittedtrackspt[nVertex] = -1;
	vertexsumRefittedtracksptSquare[nVertex] = -1; 
      }

      if(debug_ >=2) cout<<"vertex "<< nVertex <<" "<< vertexx[nVertex]<<" "<< vertexy[nVertex]<<" "<<vertexz[nVertex]<<endl; 
      
      
      nVertex++ ; 
    }
    
  }catch(std::exception& ex ){
    cout<<"vertex not working.."<<endl;
  }
  //////========== End of veretx ==============////////
  

  phyDeclared = filter(iEvent, iSetup);  //save this into root tree
  
  
  
  try{
    
    nVertexNoBS = 0; 
        
    Handle<VertexCollection> vertices;
    iEvent.getByLabel(m_vertexSrc2, vertices);
    for (int j = 0; j < int(vertices->size()) && j < MAXVX; j++){
      vertexNoBSx[nVertexNoBS] = (*vertices)[j].x();
      vertexNoBSy[nVertexNoBS] = (*vertices)[j].y();
      vertexNoBSz[nVertexNoBS] = (*vertices)[j].z();
      vertexNoBSchi2[nVertexNoBS] = (*vertices)[j].chi2();
      vertexNoBSndof[nVertexNoBS] = (*vertices)[j].ndof(); 
      vertexNoBSnormalizedChi2[nVertexNoBS] = (*vertices)[j].normalizedChi2(); 
      vertexNoBStrackSize[nVertexNoBS] = (*vertices)[j].tracksSize(); 
      vertexNoBSisFake[nVertexNoBS] = (*vertices)[j].isFake();
      vertexNoBSisValid[nVertexNoBS] = (*vertices)[j].isValid();
      
      nVertexNoBS++ ; 
      
    }
  }
  catch(std::exception& ex ){
    cout<<"vertex nobs not working.."<<endl;
  }
  //////========== End of veretx ==============////////
  


  
  
  Handle<reco::SuperClusterCollection> scBarrelHandle;
  nSC = 0; 
  try{
    iEvent.getByLabel(m_barrelSC, scBarrelHandle);
    reco::SuperClusterCollection::const_iterator it_sc=  scBarrelHandle->begin();
    for(;it_sc != scBarrelHandle->end() && nSC<MAXSC; it_sc++,nSC++) {
      etaSC[nSC] = it_sc->position().eta(); 
      phiSC[nSC] = it_sc->position().phi(); 
      eSC[nSC] = it_sc->energy();
      flagSC[nSC] = 0; 
    }
  }catch ( std::exception& ex ) {
    cout<<"barrelSC  not working.."<<endl;
  }
  Handle<reco::SuperClusterCollection> scEndcapHandle;
  try{
    iEvent.getByLabel(m_endcapSC, scEndcapHandle);
    reco::SuperClusterCollection::const_iterator it_sc=  scEndcapHandle->begin();
    for(;it_sc != scEndcapHandle->end() && nSC<MAXSC; it_sc++,nSC++) {
      etaSC[nSC] = it_sc->position().eta(); 
      phiSC[nSC] = it_sc->position().phi(); 
      eSC[nSC] = it_sc->energy();
      
    }
  }catch ( std::exception& ex ) {
    cout<<"endcapSC  not working.."<<endl;
  }
  
  
  
  edm::Handle<EcalRecHitCollection> Brechit;//EB                                                                                                     
  edm::Handle<EcalRecHitCollection> Erechit;//endcap                                                                                                     
  iEvent.getByLabel(m_ecalRecHitBarrel,Brechit);
  iEvent.getByLabel(m_ecalRecHitEndcap,Erechit);
  EcalClusterLazyTools lazyTool( iEvent, iSetup, m_ecalRecHitBarrel,m_ecalRecHitEndcap);
  const EcalRecHitCollection* barrelRecHits= Brechit.product();  ///used in reco::electrons for swiss cross                                              
  const EcalRecHitCollection* endcapRecHits= Erechit.product();
  
  float res[100];

  //Get the channel status from the db
  edm::ESHandle<EcalChannelStatus> chStatus;
  iSetup.get<EcalChannelStatusRcd>().get(chStatus);
  
  edm::ESHandle<EcalSeverityLevelAlgo> sevlv;
  iSetup.get<EcalSeverityLevelAlgoRcd>().get(sevlv);
  
  

  edm::Handle<reco::GsfElectronCollection> electrons;
  
  vector<int> indPTSortedEle; 
  
  
  try{
    iEvent.getByLabel(m_electronSrc, electrons);
    reco::GsfElectronCollection etSortedele = *electrons; 
    std::stable_sort( etSortedele.begin(), etSortedele.end(), EtSorter());
    
    
//     GsfElectronCollection::const_iterator itele = etSortedele.begin(); 
//     GsfElectronCollection::const_iterator itele2 = etSortedele.begin();  
//     ///flag duplicated ele, ( same superCluster) First
//     vector<int> indDup; 
//     itele = etSortedele.begin(); 
//     for(int j=0; itele != etSortedele.end() ; ++itele,j++){
//       /////  already flag as dups.
//       vector<int>::iterator itd = find(indDup.begin(),indDup.end(),j);
//       if(itd != indDup.end()) continue; 
      
//       reco::SuperClusterRef sclRef=itele->superCluster();
//       float eop = itele->eSuperClusterOverP();
//       int k=0; 
//       for(itele2 = etSortedele.begin(); itele2 != etSortedele.end() ; ++itele2,k++){
// 	reco::SuperClusterRef sclRef2=itele2->superCluster();
// 	if(k == j) continue; 
// 	float eop2 = itele2->eSuperClusterOverP();
// 	if(sclRef == sclRef2){
// 	  int tmp = j; 
// 	  if(fabs(eop-1) < fabs(eop2-1)) tmp = k; ////j is better 
// 	  else if(fabs(eop-1) > fabs(eop2-1)) tmp = j; 
// 	  else tmp = j<k? k:j; 
// 	  vector<int>::iterator itd = find(indDup.begin(),indDup.end(),tmp);
// 	  if(itd ==indDup.end()) indDup.push_back(tmp); 
	  
// 	}
//       }
//     }
    
//     if(int(indDup.size()) >=1){ ///print out if find the duplicated electron
//       cout<<"founddupele:  " <<" "<< runNumber <<" "<< evtNumber <<" "<< lumiBlock <<endl; 
//       for(int j=0; j< int(indDup.size()) ; j++){
// 	cout<< indDup[j]<<" "; 
//       }
//       cout<<endl; 
//       itele = etSortedele.begin(); 
//       for( int j=0; itele != etSortedele.end(); ++itele,j++){
// 	cout<<"j: "<< j<<" "<< itele->superCluster()->eta() <<" "<< itele->superCluster()->phi() <<" "<< itele->eSuperClusterOverP()<<" "<< itele->ecalDrivenSeed() <<" "<< itele->trackerDrivenSeed() <<endl; 
//       }
      
//     }
    
 //    ///find index of original electrons
//     vector<int> indPTSorted; 
    
//     itele = etSortedele.begin(); 
//     for( int j=0; itele != etSortedele.end(); ++itele,j++){
      
//       vector<int>::iterator itd = find(indDup.begin(),indDup.end(),j);
//       if(itd != indDup.end()) continue; 
     
//       itele2 = electrons->begin(); 
      
//      int indfound = -1; 
//      float drMin = 0.1; 
//      for( int k=0; itele2 != electrons->end(); ++itele2,k++){
//        vector<int>::iterator itd = find(indPTSorted.begin(),indPTSorted.end(),k);
//        if(itd != indPTSorted.end()) continue; 
//        float dr = GetDeltaR(itele->eta(),itele2->eta(),itele->phi(),itele2->phi());
//        if(dr<drMin){
// 	 drMin = dr; 
// 	 indfound = k; 
//        }
//      }
//      if(indfound>-1) indPTSorted.push_back(indfound);
//      else{
//        cout<<"wrong_indPTsrotedelectron_notfound.."<<itele->pt()<<" "<<itele->eta()<<endl;
//      }
//     }
    
    //NOW SAVE all  ( clean duplicates offline if needed) 
    electronscbclusterenergy->clear();

    if( saveAllRecHitsSC_){
      electronscrechitenergy->clear();
      electronscrechitfraction->clear();
      electronscrechitieta->clear();
      electronscrechitiphi->clear();
      electronscrechitlaserCorr->clear();
    }
    
    nElectron = 0; 
    for(GsfElectronCollection::const_iterator jj = etSortedele.begin(); jj< etSortedele.end() && nElectron < nElectronMAX ; jj++){
     electronpt[nElectron] = jj->pt();
     electroneta[nElectron] = jj->eta();
     electronphi[nElectron] = jj->phi();
     electroncharge[nElectron]                        = jj->charge();
     
     
     
     if(debug_>1) cout<<"electron : " << jj->pt() <<endl; 
     
     
     electronvertexx[nElectron] = jj->vertex().x();
     electronvertexy[nElectron] = jj->vertex().y();
     electronvertexz[nElectron] = jj->vertex().z();

     reco::SuperClusterRef scRef = jj->superCluster();
     electronscrawEnergy[nElectron] = scRef->rawEnergy();
     electronsceta[nElectron] = scRef->eta();
     electronscphi[nElectron] = scRef->phi();
     electronscenergy[nElectron] = scRef->energy();
     electronisEcalEnergyCorrected[nElectron] = jj->isEcalEnergyCorrected();
     electronecalEnergy[nElectron] = jj->ecalEnergy();
          
     electronscpreshowerEnergy[nElectron] = scRef->preshowerEnergy();
     electronscphiWidth[nElectron] = scRef->phiWidth();
     electronscetaWidth[nElectron] = scRef->etaWidth();
     electronscclusterSize[nElectron] = scRef->clustersSize();
     std::vector<std::pair<DetId, float> > v_id = scRef->hitsAndFractions();
     electronscnhits[nElectron] = int(v_id.size());
     vector<float> bcenergy;
     for(reco::CaloCluster_iterator cIt = scRef->clustersBegin(); cIt != scRef->clustersEnd(); ++cIt) {
       const reco::CaloClusterPtr cc = *cIt;
       bcenergy.push_back( cc->energy());
     }
     electronscbclusterenergy->push_back(bcenergy);
     
     electroncaloPositionx[nElectron] = scRef->position().x();
     electroncaloPositiony[nElectron] = scRef->position().y();
     electroncaloPositionz[nElectron] = scRef->position().z();
     
     
     electronfiduficalFlag[nElectron] = 0; 
     if( jj->isEB()){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 0; 
     }
     if( jj->isEE()){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 1; 
     }
     if( jj->isGap()){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 2; 
     }
     if( jj->isEBEEGap()){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 3; 
     }
     if( jj->isEBEtaGap() ){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 4; 
     }
     if( jj->isEBPhiGap() ){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 5; 
     }
     if( jj->isEEGap() ){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 6; 
     }
     if( jj->isEEDeeGap() ){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 7; 
     }
     if( jj->isEERingGap() ){
       electronfiduficalFlag[nElectron] = electronfiduficalFlag[nElectron]  | 1 << 8; 
     }
          
     electrontrackerDrivenSeed[nElectron]             = jj->trackerDrivenSeed();
     electronecalDrivenSeed[nElectron]                = jj->ecalDrivenSeed();
     
     electronfbrem[nElectron]                         = jj->fbrem();
     electronnumberOfBrems[nElectron] = jj->numberOfBrems();
     //isolation 
     electrondr03TkSumPt[nElectron]                   = jj->dr03TkSumPt();
     electrondr03EcalRecHitSumEt[nElectron]           = jj->dr03EcalRecHitSumEt();
     electrondr03HcalDepth1TowerSumEt[nElectron]      = jj->dr03HcalDepth1TowerSumEt();
     electrondr03HcalDepth2TowerSumEt[nElectron]      = jj->dr03HcalDepth2TowerSumEt();
     electrondr03HcalTowerSumEt[nElectron]            = jj->dr03HcalTowerSumEt();

     electrondr04EcalRecHitSumEt[nElectron]           = jj->dr04EcalRecHitSumEt();
     electrondr04HcalDepth1TowerSumEt[nElectron]      = jj->dr04HcalDepth1TowerSumEt();
     electrondr04HcalDepth2TowerSumEt[nElectron]      = jj->dr04HcalDepth2TowerSumEt();
     electrondr04HcalTowerSumEt[nElectron]            = jj->dr04HcalTowerSumEt();
     electrondr04TkSumPt[nElectron]                   = jj->dr04TkSumPt();
     
     electronhcalDepth1OverEcal[nElectron]= jj->hcalDepth1OverEcal();
     electronhcalDepth2OverEcal[nElectron] = jj->hcalDepth2OverEcal();
     electronhcalOverEcal[nElectron] = jj->hcalOverEcal();
     
     //showershape
     electrone1x5[nElectron]                          = jj->e1x5();
     electrone2x5Max[nElectron]                       = jj->e2x5Max();
     electrone5x5[nElectron]                          = jj->e5x5();
     electronsigmaIetaIeta[nElectron]                 = jj->sigmaIetaIeta();
     electrone3x3[nElectron] = lazyTool.e3x3(*(scRef->seed())); 
     electroneMax[nElectron] = lazyTool.eMax( *(scRef->seed())); 
     
     
      ///more variables from lazyTools
     electroneLeft[nElectron] = lazyTool.eLeft( *(scRef->seed()));
     electroneRight[nElectron] = lazyTool.eRight( *(scRef->seed()));
     electroneBottom[nElectron] = lazyTool.eBottom( *(scRef->seed()));
     electroneTop[nElectron] = lazyTool.eTop( *(scRef->seed()));
     
      electrone1x3[nElectron] = lazyTool.e1x3( *(scRef->seed()) );
      electrone3x1[nElectron] = lazyTool.e3x1( *(scRef->seed()) );
      electrone2x2[nElectron] = lazyTool.e2x2( *(scRef->seed()) );
      electrone3x2[nElectron] = lazyTool.e3x2( *(scRef->seed()) );
      electrone4x4[nElectron] = lazyTool.e4x4( *(scRef->seed()) );
      electrone2x5Right[nElectron] = lazyTool.e2x5Right( *(scRef->seed()) );
      electrone2x5Left[nElectron] = lazyTool.e2x5Left( *(scRef->seed()) );
      electrone2x5Top[nElectron] = lazyTool.e2x5Top( *(scRef->seed()) );
      electrone2x5Bottom[nElectron] = lazyTool.e2x5Bottom( *(scRef->seed()) );
      ////electrone2x5Max[nElectron] = lazyTool.e2x5Max( *(scRef->seed()) );


     
     ///ID
     electroneSuperClusterOverP[nElectron] = jj->eSuperClusterOverP();
     electroneSeedClusterOverP[nElectron]             = jj->eSeedClusterOverP();
     electroneSeedClusterOverPout[nElectron]          = jj->eSeedClusterOverPout();
     electroneEleClusterOverPout[nElectron]           = jj->eEleClusterOverPout();
     electrondeltaEtaSuperClusterTrackAtVtx[nElectron]= jj->deltaEtaSuperClusterTrackAtVtx();
     electrondeltaEtaSeedClusterTrackAtCalo[nElectron]= jj->deltaEtaSeedClusterTrackAtCalo();
     electrondeltaEtaEleClusterTrackAtCalo[nElectron] = jj->deltaEtaEleClusterTrackAtCalo();
     electrondeltaPhiSuperClusterTrackAtVtx[nElectron]= jj->deltaPhiSuperClusterTrackAtVtx();
     electrondeltaPhiSeedClusterTrackAtCalo[nElectron]= jj->deltaPhiSeedClusterTrackAtCalo();
     electrondeltaPhiEleClusterTrackAtCalo[nElectron] = jj->deltaPhiEleClusterTrackAtCalo();
     electronclassification[nElectron] = jj->classification();
     
     electronmva[nElectron]                           = jj->mva();  
     electronnumberOfTracks[nElectron]                = jj->numberOfTracks();  ///0 ??
     
          ///conversion
     electronconvDist[nElectron] = jj->convDist();
     electronconvDcot[nElectron] = jj->convDcot();
     electronconvRadius[nElectron] = jj->convRadius();
     electronconvFlags[nElectron] =  jj->convFlags();
     
     ///gsfTrack
     reco::GsfTrackRef trackref = jj->gsfTrack();
     if(trackref.isNonnull()){
       const reco::HitPattern& p_inner = trackref->trackerExpectedHitsInner();
       electronExpectedHitsInnernumberOfHits[nElectron] = p_inner.numberOfHits();
       electronExpectedHitsOuternumberOfHits[nElectron] = trackref->trackerExpectedHitsOuter().numberOfHits();
       
       electrongsfTrackvx[nElectron] = trackref->vx();
       electrongsfTrackvy[nElectron] = trackref->vy();
       electrongsfTrackvz[nElectron] = trackref->vz();
       electrongsfTracknormalizedChi2[nElectron] =  trackref->normalizedChi2();
       electrongsfTrackdxybeamSpot[nElectron] =  trackref->dxy(beamPoint);
       electrongsfTrackdzbeamSpot[nElectron] =  trackref->dz(beamPoint);
       electrongsfTracknumberOfValidHits[nElectron] =  trackref->numberOfValidHits(); 
       electrongsfTracknumberOfLostHits[nElectron]           = trackref->numberOfLostHits();
       electrongsfTracknumberOfValidPixelHits[nElectron] = trackref->hitPattern().numberOfValidPixelHits();
       electrongsfTracknumberOfValidTrackerHits[nElectron] = trackref->hitPattern().numberOfValidTrackerHits();
       
       electrongsfTrackpt[nElectron] = trackref->pt();
       electrongsfTracketa[nElectron] = trackref->eta();
       electrongsfTrackphi[nElectron] = trackref->phi();
       
     }else{
       //electronExpectednumberOfHits[nElectron] = -99; 
       electronExpectedHitsInnernumberOfHits[nElectron] = -99; 
       electronExpectedHitsOuternumberOfHits[nElectron] = -99; 
       
       electrongsfTrackvx[nElectron] = -99; 
       electrongsfTrackvy[nElectron] = -99; 
       electrongsfTrackvz[nElectron] = -99; 
       electrongsfTracknormalizedChi2[nElectron] =  -99; 
       electrongsfTrackdxybeamSpot[nElectron] =  -99;
       electrongsfTrackdzbeamSpot[nElectron] =  -99;
       electrongsfTracknumberOfValidHits[nElectron] = -99; 
       electrongsfTracknumberOfLostHits[nElectron] = -99; 
       electrongsfTracknumberOfValidPixelHits[nElectron] = -99; 
       electrongsfTracknumberOfValidTrackerHits[nElectron] = -99; 
       electrongsfTrackpt[nElectron] = -99;
       electrongsfTracketa[nElectron] = -99; 
       electrongsfTrackphi[nElectron] = -99; 
     }
     
     vector<float> vlat = lazyTool.lat(*(scRef->seed()) );
     for(int k=0 ; k<3; k++){
       electronlat[nElectron][k] = vlat[k];
     }
     
     vector<float> covariances = lazyTool.covariances(*(scRef->seed()));
     electronCovEtaEta[nElectron] = covariances[0];
     electronCovEtaPhi[nElectron] = covariances[1];
     electronCovPhiPhi[nElectron] = covariances[2];
     
     vector<float> lcovariances = lazyTool.localCovariances(*(scRef->seed()));
     electronCovIEtaIEta[nElectron] = lcovariances[0];
     electronCovIEtaIPhi[nElectron] = lcovariances[1];
     electronCovIPhiIPhi[nElectron] = lcovariances[2];
     
     
     ///idmax's   
     
     DetId idmax = lazyTool.getMaximum(*(scRef->seed())).first; 
     const EcalRecHitCollection & rechits  = idmax.subdetId() == EcalBarrel ? *Brechit : *Erechit; 
     EcalRecHitCollection::const_iterator it = rechits.find( idmax );
     
     ///for spike's
     electronswissCross[nElectron] = EcalTools::swissCross(idmax,rechits,1,false); //compute only if Et(seed) > 1GeV, ieta==85 as well

     if(idmax.subdetId() == EcalBarrel){
       EBDetId det = EBDetId(idmax);
       electronieta[nElectron] = det.ieta();
       electroniphi[nElectron] = det.iphi();
     }else{
       EEDetId det = EEDetId(idmax);
       electronieta[nElectron] = det.ix();
       electroniphi[nElectron] = det.iy();
     }
     electronE2overE9[nElectron] = mye2overe9(idmax,rechits);
     
     if( it != rechits.end() ) { 
       electronseedtime[nElectron] = it->time(); 
       electronseedoutOfTimeChi2[nElectron]= it->outOfTimeChi2();
       electronseedchi2[nElectron]= it->chi2();
       electronseedrecoFlag[nElectron]= it->recoFlag();
       electronseedseverityLevel[nElectron]= sevlv->severityLevel(idmax, rechits);
       if( saveLaserCorrSeedCrystal_) electronseedlaserCorr[nElectron] = laser->getLaserCorrection(idmax,iEvent.time());
     }else{
       electronseedtime[nElectron] = -99; 
       electronseedoutOfTimeChi2[nElectron]= -99; 
       electronseedchi2[nElectron]= -99; 
       electronseedrecoFlag[nElectron]= -99; 
       electronseedseverityLevel[nElectron]= -99;
       electronseedlaserCorr[nElectron] = -99;
     }
     
     if( saveAllRecHitsSC_)       {
       const std::vector<std::pair<DetId,float> >& hits = scRef->hitsAndFractions();
       vector<float> rhenergy;
       vector<float> rhfraction;
       vector<short> rhieta;
       vector<short> rhiphi;
       vector<float> rhlasercorr; 
       for(std::vector<std::pair<DetId,float> >::const_iterator rh = hits.begin(); rh!=hits.end(); 	   ++rh){
	 float rhLaserCorrection = -1.;
	 if ((*rh).first.subdetId()== EcalBarrel) {
	   EBRecHitCollection::const_iterator itrechit = barrelRecHits->find((*rh).first);
	   if( itrechit != barrelRecHits->end()){
	     EBDetId barrelId (itrechit->id ());
	     rhLaserCorrection = laser->getLaserCorrection(barrelId, iEvent.time());
	     rhenergy.push_back(itrechit->energy());

	     ///only matters for tracker driven ele
	     if( !jj->ecalDrivenSeed() ) rhfraction.push_back((*rh).second);
	     
	     rhlasercorr.push_back(rhLaserCorrection);
	     rhieta.push_back(barrelId.ieta());
	     rhiphi.push_back(barrelId.iphi());
	   }
	 }
	 if ((*rh).first.subdetId()== EcalEndcap) {
	   EERecHitCollection::const_iterator itrechit = endcapRecHits->find((*rh).first);
	   if( itrechit != endcapRecHits->end()){
	     EEDetId endcapId (itrechit->id ());
	     rhLaserCorrection = laser->getLaserCorrection(endcapId, iEvent.time());
	     rhenergy.push_back(itrechit->energy());
	     if( !jj->ecalDrivenSeed() ) rhfraction.push_back((*rh).second);
	     rhlasercorr.push_back(rhLaserCorrection);
	     rhieta.push_back(endcapId.ix());
	     rhiphi.push_back(endcapId.iy());
	   }
	 }
       }
       if( jj->ecalDrivenSeed() ) {
	 rhfraction.push_back(1); //dummy value not to be used
       }
       electronscrechitenergy->push_back(rhenergy);
       electronscrechitfraction->push_back(rhfraction);
       electronscrechitlaserCorr->push_back(rhlasercorr);
       electronscrechitieta->push_back(rhieta);
       electronscrechitiphi->push_back(rhiphi);
       
     }
     


     // MIG                                                                                     
     if( saveLaserCorrSeedCrystal_)       {
       const std::vector<std::pair<DetId,float> >& hits = scRef->hitsAndFractions();
       float sumRecHitE = 0.;
       float sumLaserCorrectionRecHitE = 0.;
       for(std::vector<std::pair<DetId,float> >::const_iterator rh = hits.begin(); rh!=hits.end(); 	   ++rh){
	 float rhLaserCorrection = -1.;
	 if ((*rh).first.subdetId()== EcalBarrel) {
	   EBRecHitCollection::const_iterator itrechit = barrelRecHits->find((*rh).first);
	   if( itrechit != barrelRecHits->end()){
	     EBDetId barrelId (itrechit->id ());
	     rhLaserCorrection = laser->getLaserCorrection(barrelId, iEvent.time());
	     sumRecHitE += itrechit->energy();
	     sumLaserCorrectionRecHitE += itrechit->energy() * rhLaserCorrection;
	   }
	 }
	 if ((*rh).first.subdetId()== EcalEndcap) {
	   EERecHitCollection::const_iterator itrechit = endcapRecHits->find((*rh).first);
	   if( itrechit != endcapRecHits->end()){
	     EEDetId endcapId (itrechit->id ());
	     rhLaserCorrection = laser->getLaserCorrection(endcapId, iEvent.time());
	     sumRecHitE += itrechit->energy();
	     sumLaserCorrectionRecHitE += itrechit->energy() * rhLaserCorrection;
	   }
	 }
       }
       float laseraver= sumRecHitE >0 ? sumLaserCorrectionRecHitE/sumRecHitE : -1;
       electronaverlaserCorr[nElectron] = laseraver;
     }
     
     
     electronscindex[nElectron] = findIndexofSC(scRef->energy(),scRef->eta(),scRef->phi());

     
     if( !isRealData){
       MatchToGenElectron(jj->eta(),jj->phi(),res);
       for(int n=0; n<4; n++) electrongenelematch[nElectron][n] = res[n];
       for(int n=0; n<4; n++) electrongentrkmatch[nElectron][n] = res[n+4];
       for(int n=0; n<7; n++) electronsimelematch[nElectron][n] = res[n+8];
       
       ///parton matching
       float drCone[3] = {0.3,0.4,0.5}; 
       for(int n=0; n<6; n++){
	 electronPartonmatch[nElectron][n] = 0; 
       }
       for(int n=0; n<3; n++){
	 int tmp1 = partonMatchingAlgo(jj->eta(),jj->phi(),drCone[n]);
	 if(tmp1>0){
	   electronPartonmatch[nElectron][n*2] = pidMCpart[tmp1];
	   electronPartonmatch[nElectron][n*2+1] = statusMCpart[tmp1];
	 }
       }
     }
     
     nElectron ++; 
    }
    
    
  }catch( std::exception& ex){
    cout<<"electron  not working.."<<endl;
  }
  


 //////////////==========START OF Conversion =============///////////
  
  // get collections
  edm::Handle<reco::ConversionCollection> convH;
  iEvent.getByLabel(m_allConversionsColl, convH);
  
  
  nConv =0; 
  convnHitsBeforeVtx->clear();
  for( reco::ConversionCollection::const_iterator  iConv = convH->begin(); iConv != convH->end() && nConv < nConvMAX; iConv++) {
    
    //if (conv_n >= MAX_CONVERTEDPHOTONS) {
    // std::cout << "GlobeConversions: WARNING TOO MANY CONVERSIONS: " << convH->size() << " (allowed " << MAX_CONVERTEDPHOTONS << ")" << std::endl;
    // break;
    // }
    

    reco::Conversion localConv = reco::Conversion(*iConv);
    convrefittedPair4Momentumeta[nConv] = localConv.refittedPair4Momentum().eta();
    convrefittedPair4Momentumphi[nConv] = localConv.refittedPair4Momentum().phi();
    convrefittedPair4Momentumpt[nConv] = localConv.refittedPair4Momentum().pt();
    convrefittedPair4Momentumenergy[nConv] = localConv.refittedPair4Momentum().energy();
    
    convcaloClustersize[nConv] = localConv.caloCluster().size(); 
    if ( localConv.caloCluster().size() > 0 ) {
      convcaloCluster0x[nConv] = (*localConv.caloCluster()[0]).x();
      convcaloCluster0y[nConv] = (*localConv.caloCluster()[0]).y();
      convcaloCluster0z[nConv] = (*localConv.caloCluster()[0]).z();
    }else{
      convcaloCluster0x[nConv] = -99; 
      convcaloCluster0y[nConv] = -99; 
      convcaloCluster0z[nConv] = -99; 
    }
    
    convpairMomentumx[nConv] = localConv.pairMomentum().x();
    convpairMomentumy[nConv] = localConv.pairMomentum().y();
    convpairMomentumz[nConv] = localConv.pairMomentum().z();
    convrefittedPairMomentumx[nConv] = localConv.refittedPairMomentum().x(); 
    convrefittedPairMomentumy[nConv] = localConv.refittedPairMomentum().y(); 
    convrefittedPairMomentumz[nConv] = localConv.refittedPairMomentum().z(); 
    
    if(debug_ >1) cout<<" ConversionCollection " << nConv <<" "<< convcaloClustersize[nConv] <<endl; 
    

    convconversionVertexisValid[nConv]=localConv.conversionVertex().isValid();
    if(localConv.conversionVertex().isValid()){
      reco::Vertex vtx=localConv.conversionVertex();
      convconversionVertexx[nConv] = vtx.x();
      convconversionVertexy[nConv] = vtx.y();
      convconversionVertexz[nConv] = vtx.z();
      convconversionVertexchi2[nConv] = vtx.chi2();
      convconversionVertexChiSquaredProbability[nConv]=ChiSquaredProbability(vtx.chi2(), vtx.ndof());
      convconversionVertexxError[nConv]= vtx.xError();
      convconversionVertexyError[nConv]= vtx.yError();
      convconversionVertexzError[nConv]= vtx.zError();
      convconversionVertexnTracks[nConv]=localConv.nTracks();
      convconversionVertexMVAout[nConv]=localConv.MVAout();
    }else{
      convconversionVertexx[nConv] = -99; 
      convconversionVertexy[nConv] = -99;
      convconversionVertexz[nConv] = -99;
      convconversionVertexchi2[nConv] = -99;
      convconversionVertexChiSquaredProbability[nConv]= -99;
      convconversionVertexxError[nConv]= -99;
      convconversionVertexyError[nConv]= -99;
      convconversionVertexzError[nConv]= -99;
      convconversionVertexnTracks[nConv]= -99;
      convconversionVertexMVAout[nConv]=-99;
    }
    convnTracks[nConv]=localConv.nTracks();
    
    if( localConv.nTracks()) {
      const std::vector<edm::RefToBase<reco::Track> > tracks = localConv.tracks();
      for (unsigned int i=0; i<tracks.size(); i++) {
	if(i==0) {
	  conv_track1_dz[nConv]=tracks[i]->dz();
	  conv_track1_dzError[nConv]=tracks[i]->dzError();
	  conv_track1_charge[nConv]=tracks[i]->charge();
	  conv_track1_algo[nConv] = tracks[i]->algo();

	  conv_track1_pt[nConv] = tracks[i]->pt();
	  conv_track1_eta[nConv] = tracks[i]->eta();
	  conv_track1_phi[nConv] = tracks[i]->phi();
	  
	}
	else if(i==1) {
	  conv_track2_dz[nConv]=tracks[i]->dz();
	  conv_track2_dzError[nConv]=tracks[i]->dzError();
	  conv_track2_charge[nConv]=tracks[i]->charge();
	  conv_track2_algo[nConv] = tracks[i]->algo();

	  conv_track2_pt[nConv] = tracks[i]->pt();
	  conv_track2_eta[nConv] = tracks[i]->eta();
	  conv_track2_phi[nConv] = tracks[i]->phi();

	}
      }
    }else{
      conv_track1_dz[nConv]= -99; 
      conv_track1_dzError[nConv]=-99;
      conv_track1_charge[nConv]=-99;
      conv_track1_algo[nConv]=-99;

      conv_track1_pt[nConv] = -99;
      conv_track1_eta[nConv] = -99;
      conv_track1_phi[nConv] = -99;

      conv_track2_dz[nConv]= -99; 
      conv_track2_dzError[nConv]=-99;
      conv_track2_charge[nConv]=-99;
      conv_track2_algo[nConv]=-99;

      conv_track2_pt[nConv] = -99;
      conv_track2_eta[nConv] = -99;
      conv_track2_phi[nConv] = -99;
            
    }
    if(debug_ >1) cout<<" ConversionCollection " << nConv <<" "<< convcaloClustersize[nConv] <<" "<< convnTracks[nConv] <<endl; 
    
    
    convpairInvariantMass[nConv]=localConv.pairInvariantMass();
    convpairCotThetaSeparation[nConv]=localConv.pairCotThetaSeparation();
    // will work in 420 conv_eoverp[nConv]=localConv.EoverPrefittedTracks();
    convEoverPrefittedTracks[nConv]=localConv.EoverPrefittedTracks();
    convzOfPrimaryVertexFromTracks[nConv]=localConv.zOfPrimaryVertexFromTracks();
    convdistOfMinimumApproach[nConv]=localConv.distOfMinimumApproach();
    convdPhiTracksAtVtx[nConv]=localConv.dPhiTracksAtVtx();
    convdPhiTracksAtEcal[nConv]=localConv.dPhiTracksAtEcal();
    convdEtaTracksAtEcal[nConv]=localConv.dEtaTracksAtEcal();
    
    if(debug_ >1) cout<<" ConversionCollection " << nConv <<" "<< convcaloClustersize[nConv] <<" "<< convnTracks[nConv] <<" "<<convdPhiTracksAtVtx[nConv] <<endl; 
    
    
    std::vector<unsigned short> tmp;
    for (unsigned int i=0; i<localConv.nHitsBeforeVtx().size(); ++i) {
      tmp.push_back(static_cast<unsigned short>(localConv.nHitsBeforeVtx()[i]));
    }
    convnHitsBeforeVtx->push_back(tmp);
    
    convnSharedHits[nConv] = localConv.nSharedHits();
    
    if(localConv.tracks().size() > 0) {
      conv_track1_d0[nConv]=localConv.tracksSigned_d0()[0];
      conv_track1_pout[nConv]=sqrt(localConv.tracksPout()[0].Mag2());
      conv_track1_pin[nConv]=sqrt(localConv.tracksPin()[0].Mag2());
    }else{
      conv_track1_d0[nConv]= -99; 
      conv_track1_pout[nConv]=-99; 
      conv_track1_pin[nConv]=-99; 
    }
    
    if(localConv.tracks().size() > 1) {
      conv_track2_d0[nConv]=localConv.tracksSigned_d0()[1];
      conv_track2_pout[nConv]=sqrt(localConv.tracksPout()[1].Mag2());
      conv_track2_pin[nConv]=sqrt(localConv.tracksPin()[1].Mag2());
    }else{
      conv_track2_d0[nConv]= -99; 
      conv_track2_pout[nConv]=-99; 
      conv_track2_pin[nConv]=-99; 
    }
    
    nConv ++; 
  }
  
  
  //////////////==========END OF Conversion =============///////////
  



  //////////////==========START OF Photon =============///////////

  Handle<PhotonCollection> photons;
  vector<int> indPtPhtSorted; 
  try {
    iEvent.getByLabel(m_photonSrc, photons);
    
    reco::PhotonCollection etSortedPht = *photons; 
    std::stable_sort( etSortedPht.begin(), etSortedPht.end(), PtSorter());
    
  //   ///index of sorted in PT
//     vector<int> indPt; 
//     for(int j=0; j< int(photons->size()); j++){
      
//       float eta = etSortedPht[j].eta();
//       float phi = etSortedPht[j].phi();
      
//       float drMin = 0.01; 
//       int indMin = -1; 
//       for(int k=0; k< int(photons->size()); k++){
// 	vector<int>::iterator itd = find(indPt.begin(),indPt.end(),k);
// 	if(itd != indPt.end()) continue; 
// 	float dr = GetDeltaR(eta,(*photons)[k].eta(),phi,(*photons)[k].phi());
// 	if(dr<drMin){
// 	  drMin = dr; 
// 	  indMin = k; 
// 	}
//       }
//       if(indMin<0){cout<<"wrong..sortPtphoton.."<<endl;}
//       else indPt.push_back(indMin);
      
//     }
    
    photonscbclusterenergy->clear();
    
    nPhoton = 0; 
    for (reco::PhotonCollection::const_iterator jj = etSortedPht.begin(); jj != etSortedPht.end()&& nPhoton <nPhotonMAX; jj++){
      ///for (int jj = 0; jj < int(photons->size()) && nPhoton <nPhotonMAX; jj++){
      ///int j = indPt[jj];
      ///reco::PhotonRef eRef(photons,j);
      /////const reco::Candidate &c = dynamic_cast<const reco::Candidate&>((*photons)[j]);
      reco::SuperClusterRef scRef=jj->superCluster();
      photonpt[nPhoton] = jj->pt(); 
      photonenergy[nPhoton] = jj->energy(); 
      photoneta[nPhoton] = jj->eta(); 
      photonphi[nPhoton] = jj->phi(); 
      
      if(debug_>1) cout<<"photon: " << jj->pt() <<endl; 
      
      //vertex
      photonvertexx[nPhoton] = jj->vertex().x();
      photonvertexy[nPhoton] = jj->vertex().y();
      photonvertexz[nPhoton] = jj->vertex().z();
      
      photonhasPixelSeed[nPhoton] =  jj->hasPixelSeed();
      photonhasConversionTracks[nPhoton] = jj->hasConversionTracks(); 
      
      ///superCluster's 
      photonscrawEnergy[nPhoton] = scRef->rawEnergy();
      photonscpreshowerEnergy[nPhoton] = scRef->preshowerEnergy();
      photonscphiWidth[nPhoton] = scRef->phiWidth();
      photonscetaWidth[nPhoton] = scRef->etaWidth();
      photonscclusterSize[nPhoton] = scRef->clustersSize();
      vector<float> bcenergy; 
      for(reco::CaloCluster_iterator cIt = scRef->clustersBegin(); cIt != scRef->clustersEnd(); ++cIt)	{
	const reco::CaloClusterPtr cc = *cIt; 
	bcenergy.push_back( cc->energy());
      }
      photonscbclusterenergy->push_back(bcenergy);
      
      std::vector<std::pair<DetId, float> > v_id = scRef->hitsAndFractions();
      photonscnhits[nPhoton] = int(v_id.size());
      
      
      photonscenergy[nPhoton] = scRef->energy();
      photonsceta[nPhoton] = scRef->eta();
      photonscphi[nPhoton] = scRef->phi();
      photoncaloPositionx[nPhoton] = jj->caloPosition().x(); 
      photoncaloPositiony[nPhoton] = jj->caloPosition().y(); 
      photoncaloPositionz[nPhoton] = jj->caloPosition().z(); 
      
      
      //showershape
      photone3x3[nPhoton] = jj->e3x3();
      photone1x5[nPhoton] = jj->e1x5();
      photone2x5[nPhoton] = jj->e2x5();
      photone5x5[nPhoton] = jj->e5x5();
      photonmaxEnergyXtal[nPhoton] = jj->maxEnergyXtal();
      photonr9[nPhoton] = jj->r9();
      photonsigmaIetaIeta[nPhoton] = jj->sigmaIetaIeta();
      
      
      ///more variables from lazyTools
      photoneLeft[nPhoton] = lazyTool.eLeft( *(scRef->seed()));
      photoneRight[nPhoton] = lazyTool.eRight( *(scRef->seed()));
      photoneBottom[nPhoton] = lazyTool.eBottom( *(scRef->seed()));
      photoneTop[nPhoton] = lazyTool.eTop( *(scRef->seed()));
      photone1x3[nPhoton] = lazyTool.e1x3( *(scRef->seed()) );
      photone3x1[nPhoton] = lazyTool.e3x1( *(scRef->seed()) );
      photone2x2[nPhoton] = lazyTool.e2x2( *(scRef->seed()) );
      photone3x2[nPhoton] = lazyTool.e3x2( *(scRef->seed()) );
      photone4x4[nPhoton] = lazyTool.e4x4( *(scRef->seed()) );
      photone2x5Right[nPhoton] = lazyTool.e2x5Right( *(scRef->seed()) );
      photone2x5Left[nPhoton] = lazyTool.e2x5Left( *(scRef->seed()) );
      photone2x5Top[nPhoton] = lazyTool.e2x5Top( *(scRef->seed()) );
      photone2x5Bottom[nPhoton] = lazyTool.e2x5Bottom( *(scRef->seed()) );
      photone2x5Max[nPhoton] = lazyTool.e2x5Max( *(scRef->seed()) );
      
      ///vector<float>
      ///photonenergyBasketFractionEta[nPhoton] = lazyTool.energyBasketFractionEta(*(scRef->seed()) );
      ///photonenergyBasketFractionPhi[nPhoton] = lazyTool.energyBasketFractionPhi(*(scRef->seed()) );
      
      
      vector<float> vlat = lazyTool.lat(*(scRef->seed()) );
      for(int k=0 ; k<3; k++){
	photonlat[nPhoton][k] = vlat[k];
      }
      
      vector<float> covariances = lazyTool.covariances(*(scRef->seed()));
      photonCovEtaEta[nPhoton] = covariances[0];
      photonCovEtaPhi[nPhoton] = covariances[1];
      photonCovPhiPhi[nPhoton] = covariances[2];
      
      vector<float> lcovariances = lazyTool.localCovariances(*(scRef->seed()));
      photonCovIEtaIEta[nPhoton] = lcovariances[0];
      photonCovIEtaIPhi[nPhoton] = lcovariances[1];
      photonCovIPhiIPhi[nPhoton] = lcovariances[2];

      
      //       vector<float> sclcovariances = lazyTool.scLocalCovariances(*scRef);
      //       photonscCovIEtaIEta[nPhoton] = sclcovariances[0];
      //       photonscCovIEtaIPhi[nPhoton] = sclcovariances[1];
      //       photonscCovIPhiIPhi[nPhoton] = sclcovariances[2];
      
      
      photonzernike20[nPhoton] =       lazyTool.zernike20( *(scRef->seed()) );
      photonzernike42[nPhoton] =       lazyTool.zernike42( *(scRef->seed()) );
      
      
      
      //isolation
      photonhadronicOverEm[nPhoton]                = jj->hadronicOverEm();
      photonecalRecHitSumEtConeDR03[nPhoton]       = jj->ecalRecHitSumEtConeDR03();
      photonhcalDepth1TowerSumEtConeDR03[nPhoton]  = jj->hcalDepth1TowerSumEtConeDR03();
      photonhcalDepth2TowerSumEtConeDR03[nPhoton]  = jj->hcalDepth2TowerSumEtConeDR03();
      photonhcalTowerSumEtConeDR03[nPhoton]        = jj->hcalTowerSumEtConeDR03();
      photontrkSumPtHollowConeDR03[nPhoton]        = jj->trkSumPtHollowConeDR03();
      photontrkSumPtSolidConeDR03[nPhoton]         = jj->trkSumPtSolidConeDR03();
      photonnTrkHollowConeDR03[nPhoton]            = jj->nTrkHollowConeDR03();
      photonnTrkSolidConeDR03[nPhoton]             = jj->nTrkSolidConeDR03();
      
      photonecalRecHitSumEtConeDR04[nPhoton]       = jj->ecalRecHitSumEtConeDR04();
      photonhcalDepth1TowerSumEtConeDR04[nPhoton]  = jj->hcalDepth1TowerSumEtConeDR04();
      photonhcalDepth2TowerSumEtConeDR04[nPhoton]  = jj->hcalDepth2TowerSumEtConeDR04();
      photonhcalTowerSumEtConeDR04[nPhoton]        = jj->hcalTowerSumEtConeDR04();
      photontrkSumPtHollowConeDR04[nPhoton]        = jj->trkSumPtHollowConeDR04();
      photontrkSumPtSolidConeDR04[nPhoton]         = jj->trkSumPtSolidConeDR04();
      photonnTrkHollowConeDR04[nPhoton]            = jj->nTrkHollowConeDR04();
      photonnTrkSolidConeDR04[nPhoton]             = jj->nTrkSolidConeDR04();
            
      ///idmax's 
      DetId idmax = lazyTool.getMaximum(*(scRef->seed())).first; 
      const EcalRecHitCollection & rechits  = idmax.subdetId() == EcalBarrel ? *Brechit : *Erechit; 
      
      ///for spike's
      photonswissCross[nPhoton] = EcalTools::swissCross(idmax,rechits,1,false); //compute only if Et(seed) > 1GeV, ieta==85 as well
      photonieta[nPhoton] = -100;
      photoniphi[nPhoton] = -100;
      if(idmax.subdetId() == EcalBarrel){
	EBDetId det = EBDetId(idmax);
	photonieta[nPhoton] = det.ieta();
	photoniphi[nPhoton] = det.iphi();
      }else{
	EEDetId det = EEDetId(idmax);
	photonieta[nPhoton] = det.ix();
	photoniphi[nPhoton] = det.iy();
      }
      
      photonE2overE9[nPhoton] = mye2overe9(idmax,rechits);
            
      EcalRecHitCollection::const_iterator it = rechits.find( idmax );
      if( it != rechits.end() ) { 
	photonseedtime[nPhoton] = it->time(); 
	photonseedoutOfTimeChi2[nPhoton]= it->outOfTimeChi2();
	photonseedchi2[nPhoton]= it->chi2();
	photonseedrecoFlag[nPhoton]= it->recoFlag();
	photonseedseverityLevel[nPhoton]= sevlv->severityLevel(idmax, rechits);
      }else{
	photonseedtime[nPhoton] = -99; 
	photonseedoutOfTimeChi2[nPhoton]= -99; 
	photonseedchi2[nPhoton]= -99; 
	photonseedrecoFlag[nPhoton]= -99; 
	photonseedseverityLevel[nPhoton]= -99;
      }
      
      //difucialFlag 
      photonfiducialFlag[nPhoton] = 0; 
      if( jj->isEB()){
	photonfiducialFlag[nPhoton] = photonfiducialFlag[nPhoton]  | 1 << 0; 
      }
      if( jj->isEE()){
	photonfiducialFlag[nPhoton] = photonfiducialFlag[nPhoton]  | 1 << 1; 
      }
      if( jj->isEBEEGap()){
	photonfiducialFlag[nPhoton] = photonfiducialFlag[nPhoton]  | 1 << 2; 
      }
      if( jj->isEBEtaGap() ){
	photonfiducialFlag[nPhoton] = photonfiducialFlag[nPhoton]  | 1 << 3; 
      }
      if( jj->isEBPhiGap() ){
	photonfiducialFlag[nPhoton] = photonfiducialFlag[nPhoton]  | 1 << 4; 
      }
      if( jj->isEEGap() ){
	photonfiducialFlag[nPhoton] = photonfiducialFlag[nPhoton]  | 1 << 5; 
      }
      if( jj->isEEDeeGap() ){
	photonfiducialFlag[nPhoton] = photonfiducialFlag[nPhoton]  | 1 << 6; 
      }
      if( jj->isEERingGap() ){
	photonfiducialFlag[nPhoton] = photonfiducialFlag[nPhoton]  | 1 << 7; 
      }
      
      
      //conversion's
      reco::ConversionRefVector conversions = jj->conversions();
      photonconversionsize[nPhoton] = int(conversions.size());
      if( photonconversionsize[nPhoton] >=1){
	reco::ConversionRef aConv=conversions[0];
	reco::Vertex vtx=aConv->conversionVertex();
	photonconversionVertexx[nPhoton] = vtx.position().x();
	photonconversionVertexy[nPhoton] = vtx.position().y();
	photonconversionVertexz[nPhoton] = vtx.position().z();
	
	
	photonconversionpairInvariantMass[nPhoton]=aConv->pairInvariantMass();
	photonconversionpairCotThetaSeparation[nPhoton]=aConv->pairCotThetaSeparation();
	photonconversionEoverPrefittedTracks[nPhoton]=aConv->EoverPrefittedTracks();
	photonconversionzOfPrimaryVertexFromTracks[nPhoton]=aConv->zOfPrimaryVertexFromTracks();
	photonconversiondistOfMinimumApproach[nPhoton]=aConv->distOfMinimumApproach();
	photonconversiondPhiTracksAtVtx[nPhoton]=aConv->dPhiTracksAtVtx();
	photonconversiondPhiTracksAtEcal[nPhoton]=aConv->dPhiTracksAtEcal();
	photonconversiondEtaTracksAtEcal[nPhoton]=aConv->dEtaTracksAtEcal();
	photonconversionnTracks[nPhoton]=aConv->nTracks();
	photonconversionMVAout[nPhoton]=aConv->MVAout();
	
	photonconversionrefittedPairMomentumx[nPhoton] = aConv->refittedPairMomentum().x();
	photonconversionrefittedPairMomentumy[nPhoton] = aConv->refittedPairMomentum().y();
	photonconversionrefittedPairMomentumz[nPhoton] = aConv->refittedPairMomentum().z();
	
	
	photonconversionVertexisValid[nPhoton]=vtx.isValid();
	if( vtx.isValid() ){
	  photonconversionVertexchi2[nPhoton]=vtx.chi2();
	  photonconversionChiSquaredProbability[nPhoton]=ChiSquaredProbability(vtx.chi2(), vtx.ndof());
	  
	  if( aConv->nTracks() > 0){
	    const std::vector<edm::RefToBase<reco::Track> > tracks = aConv->tracks();
	    for (unsigned int i=0; i<tracks.size(); i++) {
	      if(i==0) {
		photonconversion_track1_dz[nPhoton]=tracks[i]->dz();
		photonconversion_track1_dzError[nPhoton]=tracks[i]->dzError();
		photonconversion_track1_charge[nPhoton]=tracks[i]->charge();
		photonconversion_track1_d0[nPhoton]=aConv->tracksSigned_d0()[0];
		photonconversion_track1_tracksPout[nPhoton]=sqrt(aConv->tracksPout()[0].Mag2());
		photonconversion_track1_tracksPin[nPhoton]=sqrt(aConv->tracksPin()[0].Mag2());
		photonconversion_track1_algo[nPhoton] = tracks[i]->algo();
	      }
	      else if(i==1) {
		photonconversion_track2_dz[nPhoton]=tracks[i]->dz();
		photonconversion_track2_dzError[nPhoton]=tracks[i]->dzError();
		photonconversion_track2_charge[nPhoton]=tracks[i]->charge();
		photonconversion_track2_d0[nPhoton]=aConv->tracksSigned_d0()[1];
		photonconversion_track2_tracksPout[nPhoton]=sqrt(aConv->tracksPout()[1].Mag2());
		photonconversion_track2_tracksPin[nPhoton]=sqrt(aConv->tracksPin()[1].Mag2());
		photonconversion_track2_algo[nPhoton] = tracks[i]->algo();
	      }
	    }
	  }else{
	    photonconversion_track1_dz[nPhoton]= -98; 
	    photonconversion_track1_dzError[nPhoton]=-98;
	    photonconversion_track1_charge[nPhoton]=-98;
	    photonconversion_track1_d0[nPhoton]=-98;
	    photonconversion_track1_tracksPout[nPhoton]=-98;
	    photonconversion_track1_tracksPin[nPhoton]=-98;
	    photonconversion_track1_algo[nPhoton] = -98;
	    photonconversion_track2_dz[nPhoton]= -98; 
	    photonconversion_track2_dzError[nPhoton]=-98;
	    photonconversion_track2_charge[nPhoton]=-98;
	    photonconversion_track2_d0[nPhoton]=-98;
	    photonconversion_track2_tracksPout[nPhoton]=-98;
	    photonconversion_track2_tracksPin[nPhoton]=-98;
	    photonconversion_track2_algo[nPhoton] = -98;
	  }
	}
	
      }else{
	photonconversionVertexx[nPhoton] = -99; 
	photonconversionVertexy[nPhoton] = -99;
	photonconversionVertexz[nPhoton] = -99; 

	photonconversionrefittedPairMomentumx[nPhoton] = -99;
	photonconversionrefittedPairMomentumy[nPhoton] = -99;
	photonconversionrefittedPairMomentumz[nPhoton] = -99;
	
	photonconversionpairInvariantMass[nPhoton]= -99; 
	photonconversionpairCotThetaSeparation[nPhoton]=-99; 
	photonconversionEoverPrefittedTracks[nPhoton]=-99; 
	photonconversionzOfPrimaryVertexFromTracks[nPhoton]=-99;
	photonconversiondistOfMinimumApproach[nPhoton]=-99;
	photonconversiondPhiTracksAtVtx[nPhoton]=-99;
	photonconversiondPhiTracksAtEcal[nPhoton]=-99;
	photonconversiondEtaTracksAtEcal[nPhoton]=-99;
	photonconversionnTracks[nPhoton]=-99;
	photonconversionMVAout[nPhoton]=-99;
	photonconversionVertexisValid[nPhoton]=-99;
	photonconversionVertexchi2[nPhoton]=-99;
	photonconversionChiSquaredProbability[nPhoton]=-99;
	
	photonconversion_track1_dz[nPhoton]= -99; 
	photonconversion_track1_dzError[nPhoton]=-99;
	photonconversion_track1_charge[nPhoton]=-99;
	photonconversion_track1_d0[nPhoton]=-99;
	photonconversion_track1_tracksPout[nPhoton]=-99;
	photonconversion_track1_tracksPin[nPhoton]=-99;
	photonconversion_track1_algo[nPhoton]=-99;
	
	photonconversion_track2_dz[nPhoton]= -99; 
	photonconversion_track2_dzError[nPhoton]=-99;
	photonconversion_track2_charge[nPhoton]=-99;
	photonconversion_track2_d0[nPhoton]=-99;
	photonconversion_track2_tracksPout[nPhoton]=-99;
	photonconversion_track2_tracksPin[nPhoton]=-99;
	photonconversion_track2_algo[nPhoton]=-99;
      }
      
      
      photonscindex[nPhoton] = findIndexofSC(scRef->energy(),scRef->eta(),scRef->phi());
      
      photonhasMatchedPromptElectron[nPhoton] =  ConversionTools::hasMatchedPromptElectron(scRef, electrons, convH, beamPoint);
      
      if( !isRealData){
	
	MatchToGenPhoton( jj->eta(),jj->phi(),res);
	for(int n=0;n<3; n++){
	  photongenphtmatch[nPhoton][n] = res[n];  //matched to a gen-photon 
	}
	for(int n=0; n<3; n++){
	  photongenelematch[nPhoton][n] = res[n+7];  //matched to a real electron
	}
	for(int n=0; n<4; n++){
	  photongenphtconv[nPhoton][n] = res[n+10];
	}
	///partonMatcher
	float drCone[3] = {0.3,0.4,0.5}; 
	for(int n=0; n<6; n++){
	  photonPartonmatch[nPhoton][n] = 0; 
	  ///store the p of the matched parton
	  if(n<3) photonPartonmatchp[nPhoton][n] = 0; 
	  
	}
	for(int n=0; n<3; n++){
	  int tmp1 = partonMatchingAlgo(jj->eta(),jj->phi(),drCone[n]);
	  if(tmp1>0){
	    photonPartonmatch[nPhoton][n*2] = pidMCpart[tmp1];
	    photonPartonmatch[nPhoton][n*2+1] = statusMCpart[tmp1];
	    photonPartonmatchp[nPhoton][n] = ptMCpart[tmp1]/sin(2*atan(exp(-etaMCpart[tmp1])));
	  }
	}
      }
      
      if(debug_>1) cout<<"photon sc index  " <<  photonscindex[nPhoton] <<endl; 
            
      
      nPhoton ++; 
    }

  }catch( std::exception& ex){
    cout<<"photon  not working.."<<endl;
  }
  
  
  //////////////==========End OF Photon =============///////////
  
  
 
  
  
  
  
  //////////////==========START OF Electron =============///////////
  


  //////////////==========START OF MUON =============///////////
  ///Muon iso map
  edm::Handle<reco::IsoDepositMap> mtkMapH;
  iEvent.getByLabel(muIsoTrkMap_, mtkMapH);
  
  edm::Handle<reco::IsoDepositMap> mecalMapH;
  iEvent.getByLabel(muIsoEcalMap_, mecalMapH);
  
  
  edm::Handle<reco::IsoDepositMap> mhcalMapH;
  iEvent.getByLabel(muIsoHcalMap_, mhcalMapH);
  
  edm::Handle<reco::IsoDepositMap> mhoMapH;
  iEvent.getByLabel(muIsoHoMap_, mhoMapH);
  

  
  try{
    edm::Handle<edm::View<reco::Muon> > muonsHandle; 
        
    Handle<reco::MuonCollection> tempmuons;
    iEvent.getByLabel(m_muonSrc, tempmuons);
    reco::MuonCollection ptSortedmuons = *tempmuons;
    std::stable_sort( ptSortedmuons.begin(), ptSortedmuons.end(), PtSorter());
    iEvent.getByLabel(m_muonSrc, muonsHandle);
    const edm::View<reco::Muon>& muons = *muonsHandle;
    vector<int> indPTSorted;
    reco::MuonCollection::const_iterator itmuon =  ptSortedmuons.begin();
    for (; itmuon != ptSortedmuons.end(); ++itmuon){
      int jfound = -1;
      float drMin = 0.01;
      for(int j=0; j<int(muons.size());j++){
        vector<int>::iterator itd = find(indPTSorted.begin(),indPTSorted.end(),j);
        if(itd != indPTSorted.end()) continue;
        float dr = GetDeltaR(itmuon->eta(),muons[j].eta(),itmuon->phi(),muons[j].phi());
        if(dr<drMin){
          drMin = dr;
          jfound = j;
        }
      }
      if(jfound >-1) indPTSorted.push_back(jfound);
      else{
        cout<<"wrong_indPTsorted_muon not found..."<<itmuon->pt()<<" "<<itmuon->eta()<<endl;
      }
    }
    
    nMuon = 0; 

    for(MuonCollection::const_iterator jj = ptSortedmuons.begin(); jj< ptSortedmuons.end() && nMuon < nMuonMAX ; jj++){
      ///      for (int j = 0; j< int(muons.size()) && nMuon < nMuonMAX;j++){
      
      int jsort = indPTSorted[nMuon];
      
      muonpt[nMuon] =  jj->pt();
      muoneta[nMuon] =  jj->eta();
      muoncharge[nMuon] =  jj->charge();
      muonphi[nMuon] =  jj->phi();
      
      if(debug_ >=2)  cout<<" muon "<< nMuon <<" "<< muonpt[nMuon] << " "<< muons[jsort].pt()<<endl; 
      
      ///recoAlgorithm
      muonRecoAlgo[nMuon] = 0; 
      if( jj->isGlobalMuon()){
	muonRecoAlgo[nMuon] =  muonRecoAlgo[nMuon] | 1 << 0; 
      }
      if( jj->isTrackerMuon()){
	muonRecoAlgo[nMuon] =  muonRecoAlgo[nMuon] | 1 << 1; 
      }
      if( jj->isStandAloneMuon()){
	muonRecoAlgo[nMuon] =  muonRecoAlgo[nMuon] | 1 << 2; 
      }
      if( jj->isCaloMuon()){
	muonRecoAlgo[nMuon] =  muonRecoAlgo[nMuon] | 1 << 3; 
      }
      
      ///isGoodMuon()
      muonisGoodMuon[nMuon] = 0; 
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::TMLastStationLoose)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<0;  
      }
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::AllArbitrated)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<1;  
      }
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::GlobalMuonPromptTight)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<2;  
      }
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::TMLastStationLoose)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<3;  
      }
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::TMLastStationTight)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<4;  
      }
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::TM2DCompatibilityLoose)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<5;  
      }
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::TM2DCompatibilityTight)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<6;  
      }
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::TMOneStationLoose)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<7;  
      }
      if( muon::isGoodMuon(ptSortedmuons[nMuon],muon::TMOneStationTight)){
	muonisGoodMuon[nMuon] = muonisGoodMuon[nMuon] | 1 <<8;  
      }
      
      
      ///muon vertex 
      muonvx[nMuon] = jj->vx(); 
      muonvy[nMuon] = jj->vy(); 
      muonvz[nMuon] = jj->vz(); 
      
      ///number of matches and chambles
      muonnumberOfMatches[nMuon] = jj->numberOfMatches(); 
      muonnumberOfChambers[nMuon] = jj->numberOfChambers();   
      
      ///global track quantities
      if(! jj->globalTrack().isNull() ){
	
	muonglobalTrackvx[nMuon] = jj->globalTrack()->vx();
	muonglobalTrackvy[nMuon] = jj->globalTrack()->vy();
	muonglobalTrackvz[nMuon] = jj->globalTrack()->vz();
	muonglobalTrackdzbeamSpot[nMuon] = jj->globalTrack()->dz(beamPoint);
	muonglobalTrackdxybeamSpot[nMuon] = jj->globalTrack()->dxy(beamPoint);
	muonglobalTracknumberOfValidPixelHits[nMuon] = jj->globalTrack()->hitPattern().numberOfValidPixelHits();
	muonglobalTracknumberOfValidTrackerHits[nMuon] = jj->globalTrack()->hitPattern().numberOfValidTrackerHits();
	muonglobalTracknormalizedChi2[nMuon] = jj->globalTrack()->normalizedChi2(); 
	muonglobalTracknumberOfValidMuonHits[nMuon] = jj->globalTrack()->hitPattern().numberOfValidMuonHits(); 
		
      }else{
	muonglobalTrackvx[nMuon] = -99; 
	muonglobalTrackvy[nMuon] = -99; 
	muonglobalTrackvz[nMuon] = -99; 
	muonglobalTrackdzbeamSpot[nMuon] = -99; 
	muonglobalTrackdxybeamSpot[nMuon] = -99; 
	muonglobalTracknumberOfValidPixelHits[nMuon] = -99; 
	muonglobalTracknumberOfValidTrackerHits[nMuon] = -99; 
	muonglobalTracknormalizedChi2[nMuon] = -99; 
	muonglobalTracknumberOfValidMuonHits[nMuon] = -99; 
      }
      

      
      ///STA Muon
      if(! jj->outerTrack().isNull() ){
	muonouterTracknumberOfValidMuonHits[nMuon] = jj->outerTrack()->hitPattern().numberOfValidMuonHits();
	muonouterTrackpt[nMuon] = jj->outerTrack()->pt();
	muonouterTracketa[nMuon] = jj->outerTrack()->eta();
	muonouterTrackphi[nMuon] = jj->outerTrack()->phi();
      }else{
	muonouterTracknumberOfValidMuonHits[nMuon] = -99;
	muonouterTrackpt[nMuon] = -99;
	muonouterTracketa[nMuon] = -99;
	muonouterTrackphi[nMuon] = -99;
      }
      
      //trackerMuon
      if(! jj->innerTrack().isNull() ){
	muoninnerTracknumberOfValidPixelHits[nMuon] =  jj->innerTrack()->hitPattern().numberOfValidPixelHits();
	muoninnerTracknumberOfValidTrackerHits[nMuon] =  jj->innerTrack()->hitPattern().numberOfValidTrackerHits();
	muoninnerTrackpt[nMuon] = jj->innerTrack()->pt();
	muoninnerTracketa[nMuon] = jj->innerTrack()->eta();
	muoninnerTrackphi[nMuon] = jj->innerTrack()->phi();
      }else{
	muoninnerTracknumberOfValidPixelHits[nMuon] =  -99;
	muoninnerTracknumberOfValidTrackerHits[nMuon] =  -99;
	muoninnerTrackpt[nMuon] = -99;
	muoninnerTracketa[nMuon] = -99;
	muoninnerTrackphi[nMuon] = -99;
      }
      
      
      ////isolation variables dr03 and dr05, and dr04 ( from isomap) 
      muonisolationR03sumPt[nMuon] = jj->isolationR03().sumPt;
      muonisolationR03emEt[nMuon] = jj->isolationR03().emEt;
      muonisolationR03hadEt[nMuon] = jj->isolationR03().hadEt;
      muonisolationR03hoEt[nMuon] = jj->isolationR03().hoEt;
      
      muonisolationR05sumPt[nMuon] = jj->isolationR05().sumPt;
      muonisolationR05emEt[nMuon] = jj->isolationR05().emEt;
      muonisolationR05hadEt[nMuon] = jj->isolationR05().hadEt;
      muonisolationR05hoEt[nMuon] = jj->isolationR05().hoEt;
      
      
      const reco::IsoDeposit tkDep((*mtkMapH)[muons.refAt(jsort)]);
      const reco::IsoDeposit ecalDep((*mecalMapH)[muons.refAt(jsort)]);
      const reco::IsoDeposit hcalDep((*mhcalMapH)[muons.refAt(jsort)]);
      const reco::IsoDeposit hoDep((*mhoMapH)[muons.refAt(jsort)]);
      muonisolationR04sumPt[nMuon] =  tkDep.depositAndCountWithin(0.4).first; 
      muonisolationR04emEt[nMuon] = ecalDep.depositAndCountWithin(0.4).first; 
      muonisolationR04hadEt[nMuon]  = hcalDep.depositAndCountWithin(0.4).first;
      muonisolationR04hoEt[nMuon] = hoDep.depositAndCountWithin(0.4).first;
      
      
      if(debug_ >=2) cout<<" muon "<< nMuon <<" "<< muonpt[nMuon] <<endl; 
      
      
      //generator-leve matching 
      if( !isRealData){
	///gen matching
	MatchToGenMuon(jj->eta(),jj->phi(),res);
	for(int n=0;n<4; n++){
	  muonGenmumatch[nMuon][n] = res[n];
	}
	for(int n=4;n<8; n++){
	  muonGentrkmatch[nMuon][n-4] = res[n];
	}
	for(int n=8;n<15; n++){
	  muonSimmumatch[nMuon][n-8] = res[n];
	}
	
	///parton matching
	float drCone[3] = {0.3,0.4,0.5}; 
	for(int n=0; n<6; n++){
	  muonPartonmatch[nMuon][n] = 0; 
	}
	for(int n=0; n<3; n++){
	  int tmp1 = partonMatchingAlgo(jj->eta(),jj->phi(),drCone[n]);
	  if(tmp1>0){
	    muonPartonmatch[nMuon][n*2] = pidMCpart[tmp1];
	    muonPartonmatch[nMuon][n*2+1] = statusMCpart[tmp1];
	  }
	}
      }
      
      nMuon ++; 
    }
    
  }catch( std::exception& ex){
    cout<<"muons  not working.."<<endl;
  }
  //////////////==========END OF MUON =============///////////
  
  
  //////////////==========START OF MET =============///////////
  
  edm::Handle< edm::View<reco::CaloMET> > caloMEThandle;
  try{
    iEvent.getByLabel("met", caloMEThandle);
    
    caloMETet                 = (caloMEThandle->front() ).et();
    caloMETsumEt                 = (caloMEThandle->front() ).sumEt();
    caloMETphi            = (caloMEThandle->front() ).phi();
    
  }catch ( std::exception& ex ) {
    
    caloMETet = -99; 
    caloMETsumEt  = -99; 
    caloMETphi  = -99; 
    caloMETsig = -99; 
    cout<<" caloMET not working.."<<endl;
  }
    
  // MET object that corrects the basic calorimeter MET for muons
  edm::Handle< edm::View<reco::CaloMET> > muCorrMEThandle;
  try{
    iEvent.getByLabel("corMetGlobalMuons", muCorrMEThandle);
    muCorrMETet                 = (muCorrMEThandle->front() ).et();
    muCorrMETsumEt                 = (muCorrMEThandle->front() ).sumEt();
    muCorrMETphi            = (muCorrMEThandle->front() ).phi();
  }catch ( std::exception& ex ) {
    muCorrMETet = -99; 
    muCorrMETsumEt  = -99; 
    muCorrMETphi  = -99; 
    cout<<" muCorrMET not working.."<<endl;
  }
  
  
  // MET object that corrects the basic calorimeter MET for muons and tracks
  edm::Handle< edm::View<reco::MET> > tcMEThandle;
  try{
    iEvent.getByLabel("tcMet", tcMEThandle);
    tcMETet                 = (tcMEThandle->front() ).et();
    tcMETsumEt                 = (tcMEThandle->front() ).sumEt();
    tcMETphi            = (tcMEThandle->front() ).phi();
  }catch ( std::exception& ex ) {
    tcMETet = -99; 
    tcMETsumEt  = -99; 
    tcMETphi  = -99; 
    cout<<" tcMET not working.."<<endl;
  }

  // MET object built as the (negative) vector sum of all particles (PFCandidates) reconstructed in the event
  edm::Handle< edm::View<reco::PFMET> > pfMEThandle;
  try{
    iEvent.getByLabel("pfMet", pfMEThandle);
    pfMETet                 = (pfMEThandle->front() ).et();
    pfMETsumEt                 = (pfMEThandle->front() ).sumEt();
    pfMETphi            = (pfMEThandle->front() ).phi();
  }catch ( std::exception& ex ) {

    pfMETet = -99; 
    pfMETsumEt  = -99; 
    pfMETphi  = -99; 
    cout<<" pfMET not working.."<<endl;
  }
  //////////////==========END OF MET =============///////////
  

  ////////////===========Jet==============////////////

  nak5CaloJet = 0; 
  try{
    //const JetCorrector* corrector = JetCorrector::getJetCorrector ("ak5CaloL1FastL2L3",iSetup);
    Handle<CaloJetCollection> jets;
    iEvent.getByLabel("ak5CaloJets", jets);
    reco::CaloJetCollection etSortedjets = *jets; 
    std::stable_sort( etSortedjets.begin(), etSortedjets.end(), EtSorter());
    CaloJetCollection::const_iterator itjet = etSortedjets.begin();
    for(; itjet != etSortedjets.end(); ++itjet ){
      if(itjet->et()< 15 ) continue;
      if(nak5CaloJet >= nCJetMAX) break; 
      ak5CaloJetet[nak5CaloJet] = itjet->et();
      ak5CaloJeteta[nak5CaloJet] = itjet->eta();
      ak5CaloJetphi[nak5CaloJet] = itjet->phi();
      ak5CaloJetenergyFractionHadronic[nak5CaloJet] = itjet->energyFractionHadronic();
      ak5CaloJetn90[nak5CaloJet] = itjet->n90();
      ak5CaloJetn60[nak5CaloJet] = itjet->n60();
      ak5CaloJettowersArea[nak5CaloJet] = itjet->towersArea();
      ak5CaloJetCalosize[nak5CaloJet] = itjet->getCaloConstituents().size(); 
      ak5CaloJetmaxEInEmTowers[nak5CaloJet] = itjet->maxEInEmTowers();
      ak5CaloJetmaxEInHadTowers[nak5CaloJet] = itjet->maxEInHadTowers();
      

      jetIDHelper.calculate(iEvent, *itjet);
      ///more vairables from jetID
      ak5CaloJetfHPD[nak5CaloJet] = jetIDHelper.fHPD();
      ak5CaloJetfRBX[nak5CaloJet] = jetIDHelper.fRBX();
      ak5CaloJetn90Hits[nak5CaloJet] = jetIDHelper.n90Hits();
      ak5CaloJetrestrictedEMF[nak5CaloJet] = jetIDHelper.restrictedEMF() ; 
      // cout<<"ak5CaloJetet jet id " << ak5CaloJetet[nak5CaloJet] << " "<< ak5CaloJetfHPD[nak5CaloJet] <<endl; 
      //  ak5CaloJetcorrection[nak5CaloJet] = corrector->correction(itjet->p4());
      // cout<<"ak5CaloJetcorr " << ak5CaloJetcorrection[nak5CaloJet] <<endl;
      
      if( !isRealData){
	///parton matching
	float drCone[3] = {0.3,0.4,0.5}; 
	for(int n=0; n<6; n++){
	  ak5CaloJetPartonmatch[nak5CaloJet][n] = 0; 
	}
	for(int n=0; n<3; n++){
	  int tmp1 = partonMatchingAlgo(itjet->eta(),itjet->phi(),drCone[n]);
	  if(tmp1>0){
	    ak5CaloJetPartonmatch[nak5CaloJet][n*2] = pidMCpart[tmp1];
	    ak5CaloJetPartonmatch[nak5CaloJet][n*2+1] = statusMCpart[tmp1];
	  }
	}
      }
      
      nak5CaloJet++; 
      
    } 
  }catch ( std::exception& ex ) {
    cout<<" ak5CaloJets  not working.."<<endl;
  }
  
  nak7CaloJet = 0; 
  try{
    //const JetCorrector* corrector = JetCorrector::getJetCorrector ("ak7CaloL1FastL2L3",iSetup);
    Handle<CaloJetCollection> jets;
    iEvent.getByLabel("ak7CaloJets", jets);
    reco::CaloJetCollection etSortedjets = *jets; 
    std::stable_sort( etSortedjets.begin(), etSortedjets.end(), EtSorter());
    CaloJetCollection::const_iterator itjet = etSortedjets.begin();
    for(; itjet != etSortedjets.end(); ++itjet ){
      if(itjet->et()< 15 ) continue;
      if(nak7CaloJet >= nCJetMAX) break; 
      ak7CaloJetet[nak7CaloJet] = itjet->et();
      ak7CaloJeteta[nak7CaloJet] = itjet->eta();
      ak7CaloJetphi[nak7CaloJet] = itjet->phi();
      ak7CaloJetenergyFractionHadronic[nak7CaloJet] = itjet->energyFractionHadronic();
      ak7CaloJetn90[nak7CaloJet] = itjet->n90();
      ak7CaloJetn60[nak7CaloJet] = itjet->n60();
      ak7CaloJettowersArea[nak7CaloJet] = itjet->towersArea();
      ak7CaloJetCalosize[nak7CaloJet] = itjet->getCaloConstituents().size(); 
      ak7CaloJetmaxEInEmTowers[nak7CaloJet] = itjet->maxEInEmTowers();
      ak7CaloJetmaxEInHadTowers[nak7CaloJet] = itjet->maxEInHadTowers();
      

      jetIDHelper.calculate(iEvent, *itjet);
      ///more vairables from jetID
      ak7CaloJetfHPD[nak7CaloJet] = jetIDHelper.fHPD();
      ak7CaloJetfRBX[nak7CaloJet] = jetIDHelper.fRBX();
      ak7CaloJetn90Hits[nak7CaloJet] = jetIDHelper.n90Hits();
      ak7CaloJetrestrictedEMF[nak7CaloJet] = jetIDHelper.restrictedEMF() ; 
      // cout<<"ak7CaloJetet jet id " << ak7CaloJetet[nak7CaloJet] << " "<< ak7CaloJetfHPD[nak7CaloJet] <<endl; 
      //  ak7CaloJetcorrection[nak7CaloJet] = corrector->correction(itjet->p4());
      // cout<<"ak7CaloJetcorr " << ak7CaloJetcorrection[nak7CaloJet] <<endl;
      
      if( !isRealData){
	///parton matching
	float drCone[3] = {0.3,0.4,0.5}; 
	for(int n=0; n<6; n++){
	  ak7CaloJetPartonmatch[nak7CaloJet][n] = 0; 
	}
	for(int n=0; n<3; n++){
	  int tmp1 = partonMatchingAlgo(itjet->eta(),itjet->phi(),drCone[n]);
	  if(tmp1>0){
	    ak7CaloJetPartonmatch[nak7CaloJet][n*2] = pidMCpart[tmp1];
	    ak7CaloJetPartonmatch[nak7CaloJet][n*2+1] = statusMCpart[tmp1];
	  }
	}
      }
      
      nak7CaloJet++; 
      
    } 
  }catch ( std::exception& ex ) {
    cout<<" ak7CaloJets  not working.."<<endl;
  }
  
  
  nak5PFJet = 0; 
  try{
    // const JetCorrector* corrector = JetCorrector::getJetCorrector ("ak5PFL1FastL2L3",iSetup);
    
    Handle<PFJetCollection> jets;
    iEvent.getByLabel("ak5PFJets", jets);
    reco::PFJetCollection etSortedjets = *jets; 
    std::stable_sort( etSortedjets.begin(), etSortedjets.end(), EtSorter());
    PFJetCollection::const_iterator itjet = etSortedjets.begin();
    for(; itjet != etSortedjets.end(); ++itjet ){
      if(itjet->et()< 15 ) continue;
      if(nak5PFJet >= nCJetMAX) break; 
      ak5PFJetet[nak5PFJet] = itjet->et();
      ak5PFJeteta[nak5PFJet] = itjet->eta();
      ak5PFJetphi[nak5PFJet] = itjet->phi();
      ak5PFJetchargedHadronEnergyFraction[nak5PFJet] = itjet->chargedHadronEnergyFraction();
      ak5PFJetchargedMuEnergyFraction[nak5PFJet] = itjet->chargedMuEnergyFraction();
      ak5PFJetneutralEmEnergyFraction[nak5PFJet] = itjet->neutralEmEnergyFraction();
      ak5PFJetchargedMultiplicity[nak5PFJet] = itjet->chargedMultiplicity();
      ak5PFJetneutralMultiplicity[nak5PFJet] = itjet->neutralMultiplicity();
      ak5PFJetPFsize[nak5PFJet] = itjet->getPFConstituents().size(); 
      

      ////cout<<" ak5PFJetet " << ak5PFJetet[nak5PFJet] <<endl; 

     //  jetIDHelper.calculate(iEvent, *itjet);
//       ///more vairables from jetID
//       ak5PFJetfHPD[nak5PFJet] = jetIDHelper.fHPD();
//       ak5PFJetfRBX[nak5PFJet] = jetIDHelper.fRBX();
//       ak5PFJetn90Hits[nak5PFJet] = jetIDHelper.n90Hits();
//       ak5PFJetrestrictedEMF[nak5PFJet] = jetIDHelper.restrictedEMF() ; 
//       cout<<"ak5PFJetet jet id " << ak5PFJetet[nak5PFJet] << " "<< ak5PFJetfHPD[nak5PFJet] <<endl; 
      //ak5PFJetcorrection[nak5PFJet] = corrector->correction(itjet->p4());
      //cout<<"ak5PFJetcorr " << ak5PFJetcorrection[nak5PFJet] <<endl;
      
      if( !isRealData){
	///parton matching
	float drCone[3] = {0.3,0.4,0.5}; 
	for(int n=0; n<6; n++){
	  ak5PFJetPartonmatch[nak5PFJet][n] = 0; 
	}
	for(int n=0; n<3; n++){
	  int tmp1 = partonMatchingAlgo(itjet->eta(),itjet->phi(),drCone[n]);
	  if(tmp1>0){
	    ak5PFJetPartonmatch[nak5PFJet][n*2] = pidMCpart[tmp1];
	    ak5PFJetPartonmatch[nak5PFJet][n*2+1] = statusMCpart[tmp1];
	  }
	}
      }
      
      nak5PFJet++; 
      
    } 
  }catch ( std::exception& ex ) {
    cout<<" ak5PFJets  not working.."<<endl;
  }
  
  
  nak7PFJet = 0; 
  try{
    // const JetCorrector* corrector = JetCorrector::getJetCorrector ("ak7PFL1FastL2L3",iSetup);
    
    Handle<PFJetCollection> jets;
    iEvent.getByLabel("ak7PFJets", jets);
    reco::PFJetCollection etSortedjets = *jets; 
    std::stable_sort( etSortedjets.begin(), etSortedjets.end(), EtSorter());
    PFJetCollection::const_iterator itjet = etSortedjets.begin();
    for(; itjet != etSortedjets.end(); ++itjet ){
      if(itjet->et()< 15 ) continue;
      if(nak7PFJet >= nCJetMAX) break; 
      ak7PFJetet[nak7PFJet] = itjet->et();
      ak7PFJeteta[nak7PFJet] = itjet->eta();
      ak7PFJetphi[nak7PFJet] = itjet->phi();
      ak7PFJetchargedHadronEnergyFraction[nak7PFJet] = itjet->chargedHadronEnergyFraction();
      ak7PFJetchargedMuEnergyFraction[nak7PFJet] = itjet->chargedMuEnergyFraction();
      ak7PFJetneutralEmEnergyFraction[nak7PFJet] = itjet->neutralEmEnergyFraction();
      ak7PFJetchargedMultiplicity[nak7PFJet] = itjet->chargedMultiplicity();
      ak7PFJetneutralMultiplicity[nak7PFJet] = itjet->neutralMultiplicity();
      ak7PFJetPFsize[nak7PFJet] = itjet->getPFConstituents().size(); 
      

      //// cout<<" ak7PFJetet " << ak7PFJetet[nak7PFJet] <<endl; 

     //  jetIDHelper.calculate(iEvent, *itjet);
//       ///more vairables from jetID
//       ak7PFJetfHPD[nak7PFJet] = jetIDHelper.fHPD();
//       ak7PFJetfRBX[nak7PFJet] = jetIDHelper.fRBX();
//       ak7PFJetn90Hits[nak7PFJet] = jetIDHelper.n90Hits();
//       ak7PFJetrestrictedEMF[nak7PFJet] = jetIDHelper.restrictedEMF() ; 
//       cout<<"ak7PFJetet jet id " << ak7PFJetet[nak7PFJet] << " "<< ak7PFJetfHPD[nak7PFJet] <<endl; 
      //ak7PFJetcorrection[nak7PFJet] = corrector->correction(itjet->p4());
      //cout<<"ak7PFJetcorr " << ak7PFJetcorrection[nak7PFJet] <<endl;
      
      if( !isRealData){
	///parton matching
	float drCone[3] = {0.3,0.4,0.5}; 
	for(int n=0; n<6; n++){
	  ak7PFJetPartonmatch[nak7PFJet][n] = 0; 
	}
	for(int n=0; n<3; n++){
	  int tmp1 = partonMatchingAlgo(itjet->eta(),itjet->phi(),drCone[n]);
	  if(tmp1>0){
	    ak7PFJetPartonmatch[nak7PFJet][n*2] = pidMCpart[tmp1];
	    ak7PFJetPartonmatch[nak7PFJet][n*2+1] = statusMCpart[tmp1];
	  }
	}
      }
      
      nak7PFJet++; 
      
    } 
  }catch ( std::exception& ex ) {
    cout<<" ak7PFJets  not working.."<<endl;
  }
  
  ///Fill the branches
  Analysis->Fill();
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
RecoAnalyzer::beginJob()
{


  cout<<" RecoAnalyzer::beginJob( " <<endl; 
  

  
  Analysis = new TTree("Analysis","Reco Analysis");
  
  
  //Event info
  Analysis->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
  Analysis->Branch("runNumber",&runNumber,"runNumber/I");
  Analysis->Branch("evtNumber",&evtNumber,"evtNumber/I");
  Analysis->Branch("bunchX",&bunchX,"bunchX/I");
  Analysis->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
  Analysis->Branch("evtTime",&evtTime,"evtTime/I");
  Analysis->Branch("isRealData",&isRealData,"isRealData/I");
  
  
  hlt_bitFired = new std::vector<unsigned short>; hlt_bitFired->clear();
  hlt_pathName = new std::vector<std::string>; hlt_pathName->clear();
  
  hlt_en = new std::vector<float>; hlt_en->clear();
  hlt_pt = new std::vector<float>; hlt_pt->clear();
  hlt_eta = new std::vector<float>; hlt_eta->clear();
  hlt_phi = new std::vector<float>; hlt_phi->clear();
    
  //HLT
  // Event Trigger
  ///hlt_p4  = new TClonesArray("TLorentzVector",   MAXHLTbits);
  
  Analysis->Branch("hlt_bitFired", "std::vector<unsigned short>", &hlt_bitFired);  ///vector of those bit which fired the trigger 
  Analysis->Branch("hlt_pathName", "std::vector<std::string>", &hlt_pathName);/// vector of string for all trigger path
  
  
  // Trigger Candidates
  //Analysis->Branch("hlt_n", &hlt_n, "hlt_n/I");
  //Analysis->Branch("hlt_en", hlt_en, "hlt_en[hlt_n]/F");
  //Analysis->Branch("hlt_pt", hlt_pt, "hlt_pt[hlt_n]/F");
  //Analysis->Branch("hlt_eta", hlt_eta, "hlt_eta[hlt_n]/F");
  //Analysis->Branch("hlt_phi", hlt_phi, "hlt_phi[hlt_n]/F");
  
  Analysis->Branch("hlt_en", "std::vector<float>", &hlt_en);
  Analysis->Branch("hlt_pt", "std::vector<float>", &hlt_pt);
  Analysis->Branch("hlt_eta", "std::vector<float>", &hlt_eta);
  Analysis->Branch("hlt_phi", "std::vector<float>", &hlt_phi);
  
  
  hlt_candpath = new std::vector<std::vector<unsigned short> >; hlt_candpath->clear();
  hltfilter_candpath = new std::vector<std::vector<unsigned short> >; hltfilter_candpath->clear();
  ///Analysis->Branch("hlt_p4", "TClonesArray", &hlt_p4, 32000, 0); ///p4 of all trigger objects in this event
  Analysis->Branch("hlt_candpath", "std::vector<std::vector<unsigned short> >", &hlt_candpath); // vector of vector ( each object's trigger path)
  Analysis->Branch("hltfilter_candpath", "std::vector<std::vector<unsigned short> >", &hltfilter_candpath); // vector of vector ( each object's trigger path)
  
 //  Analysis->Branch("hlt_candpath_1"    , &hlt_candpath_1    , "hlt_candpath_1[hlt_n]/I");
//   Analysis->Branch("hlt_candpath_2"    , &hlt_candpath_2    , "hlt_candpath_2[hlt_n]/I");
//   Analysis->Branch("hlt_candpath_3"    , &hlt_candpath_3    , "hlt_candpath_3[hlt_n]/I");
//   Analysis->Branch("hlt_candpath_4"    , &hlt_candpath_4    , "hlt_candpath_4[hlt_n]/I");
//   Analysis->Branch("hlt_candpath_5"    , &hlt_candpath_5    , "hlt_candpath_5[hlt_n]/I");
//   Analysis->Branch("hlt_candpath_6"    , &hlt_candpath_6    , "hlt_candpath_6[hlt_n]/I");
  
//   Analysis->Branch("hlt_candpath"    , &hlt_candpath    , "hlt_candpath[4][hlt_n]/i");
  
  hlt_label_names = new std::vector<std::string>; hlt_label_names ->clear();
  Analysis->Branch("hlt_label_names", "std::vector<std::string>", &hlt_label_names);   ////  vector of string for all trig object's fiter names
  
  

  pileupBunchX = new std::vector<short>; pileupBunchX->clear();
  pileupNInteraction = new std::vector<short>; pileupNInteraction->clear();
  
  
  //nPileUP
  Analysis->Branch("pileupBunchX","std::vector<short>", &pileupBunchX);
  Analysis->Branch("pileupNInteraction","std::vector<short>", &pileupNInteraction);
  
  
  Analysis->Branch("rho", &rho,"rho/F");
  Analysis->Branch("rhoEtaMax44", &rhoEtaMax44,"rhoEtaMax44/F");
  
  
  /////Photon
  Analysis->Branch("nPhoton",&nPhoton, "nPhoton/I");

  Analysis->Branch("photonpt",photonpt,"photonpt[nPhoton]/F");
Analysis->Branch("photonenergy",photonenergy,"photonenergy[nPhoton]/F");
Analysis->Branch("photoneta",photoneta,"photoneta[nPhoton]/F");
Analysis->Branch("photonphi",photonphi,"photonphi[nPhoton]/F");
Analysis->Branch("photonvertexx",photonvertexx,"photonvertexx[nPhoton]/F");
Analysis->Branch("photonvertexy",photonvertexy,"photonvertexy[nPhoton]/F");
Analysis->Branch("photonvertexz",photonvertexz,"photonvertexz[nPhoton]/F");
Analysis->Branch("photonhasPixelSeed",photonhasPixelSeed,"photonhasPixelSeed[nPhoton]/I");
Analysis->Branch("photonhasConversionTracks",photonhasConversionTracks,"photonhasConversionTracks[nPhoton]/I");

 

 Analysis->Branch("photonscrawEnergy",photonscrawEnergy,"photonscrawEnergy[nPhoton]/F");
 Analysis->Branch("photonscpreshowerEnergy",photonscpreshowerEnergy,"photonscpreshowerEnergy[nPhoton]/F");
 Analysis->Branch("photonscphiWidth",photonscphiWidth,"photonscphiWidth[nPhoton]/F");
 Analysis->Branch("photonscetaWidth",photonscetaWidth,"photonscetaWidth[nPhoton]/F");
  
 Analysis->Branch("photonscenergy",photonscenergy,"photonscenergy[nPhoton]/F");
  

Analysis->Branch("photonsceta",photonsceta,"photonsceta[nPhoton]/F");
Analysis->Branch("photonscphi",photonscphi,"photonscphi[nPhoton]/F");
Analysis->Branch("photone3x3",photone3x3,"photone3x3[nPhoton]/F");
Analysis->Branch("photone1x5",photone1x5,"photone1x5[nPhoton]/F");
Analysis->Branch("photone2x5",photone2x5,"photone2x5[nPhoton]/F");
Analysis->Branch("photone5x5",photone5x5,"photone5x5[nPhoton]/F");
Analysis->Branch("photonmaxEnergyXtal",photonmaxEnergyXtal,"photonmaxEnergyXtal[nPhoton]/F");
Analysis->Branch("photonr9",photonr9,"photonr9[nPhoton]/F");
Analysis->Branch("photonsigmaIetaIeta",photonsigmaIetaIeta,"photonsigmaIetaIeta[nPhoton]/F");
Analysis->Branch("photonhadronicOverEm",photonhadronicOverEm,"photonhadronicOverEm[nPhoton]/F");
Analysis->Branch("photonecalRecHitSumEtConeDR03",photonecalRecHitSumEtConeDR03,"photonecalRecHitSumEtConeDR03[nPhoton]/F");
Analysis->Branch("photonhcalDepth1TowerSumEtConeDR03",photonhcalDepth1TowerSumEtConeDR03,"photonhcalDepth1TowerSumEtConeDR03[nPhoton]/F");
Analysis->Branch("photonhcalDepth2TowerSumEtConeDR03",photonhcalDepth2TowerSumEtConeDR03,"photonhcalDepth2TowerSumEtConeDR03[nPhoton]/F");
Analysis->Branch("photonhcalTowerSumEtConeDR03",photonhcalTowerSumEtConeDR03,"photonhcalTowerSumEtConeDR03[nPhoton]/F");
Analysis->Branch("photontrkSumPtHollowConeDR03",photontrkSumPtHollowConeDR03,"photontrkSumPtHollowConeDR03[nPhoton]/F");
Analysis->Branch("photontrkSumPtSolidConeDR03",photontrkSumPtSolidConeDR03,"photontrkSumPtSolidConeDR03[nPhoton]/F");
Analysis->Branch("photonnTrkHollowConeDR03",photonnTrkHollowConeDR03,"photonnTrkHollowConeDR03[nPhoton]/I");
Analysis->Branch("photonnTrkSolidConeDR03",photonnTrkSolidConeDR03,"photonnTrkSolidConeDR03[nPhoton]/I");
Analysis->Branch("photonecalRecHitSumEtConeDR04",photonecalRecHitSumEtConeDR04,"photonecalRecHitSumEtConeDR04[nPhoton]/F");
Analysis->Branch("photonhcalDepth1TowerSumEtConeDR04",photonhcalDepth1TowerSumEtConeDR04,"photonhcalDepth1TowerSumEtConeDR04[nPhoton]/F");
Analysis->Branch("photonhcalDepth2TowerSumEtConeDR04",photonhcalDepth2TowerSumEtConeDR04,"photonhcalDepth2TowerSumEtConeDR04[nPhoton]/F");
Analysis->Branch("photonhcalTowerSumEtConeDR04",photonhcalTowerSumEtConeDR04,"photonhcalTowerSumEtConeDR04[nPhoton]/F");
Analysis->Branch("photontrkSumPtHollowConeDR04",photontrkSumPtHollowConeDR04,"photontrkSumPtHollowConeDR04[nPhoton]/F");
Analysis->Branch("photontrkSumPtSolidConeDR04",photontrkSumPtSolidConeDR04,"photontrkSumPtSolidConeDR04[nPhoton]/F");
Analysis->Branch("photonnTrkHollowConeDR04",photonnTrkHollowConeDR04,"photonnTrkHollowConeDR04[nPhoton]/I");
Analysis->Branch("photonnTrkSolidConeDR04",photonnTrkSolidConeDR04,"photonnTrkSolidConeDR04[nPhoton]/I");
Analysis->Branch("photonseedtime",photonseedtime,"photonseedtime[nPhoton]/F");
Analysis->Branch("photonseedoutOfTimeChi2",photonseedoutOfTimeChi2,"photonseedoutOfTimeChi2[nPhoton]/F");
Analysis->Branch("photonseedchi2",photonseedchi2,"photonseedchi2[nPhoton]/F");
Analysis->Branch("photonseedrecoFlag",photonseedrecoFlag,"photonseedrecoFlag[nPhoton]/I");
 Analysis->Branch("photonseedseverityLevel",photonseedseverityLevel,"photonseedseverityLevel[nPhoton]/I");
 Analysis->Branch("photonfiducialFlag",photonfiducialFlag,"photonfiducialFlag[nPhoton]/I");
 Analysis->Branch("photonconversionsize",photonconversionsize,"photonconversionsize[nPhoton]/I");
 Analysis->Branch("photonscindex",photonscindex,"photonscindex[nPhoton]/I");
 
 Analysis->Branch("photonhasMatchedPromptElectron",photonhasMatchedPromptElectron,"photonhasMatchedPromptElectron[nPhoton]/I");
 
 
//  Analysis->Branch("photongenphtmatch",photongenphtmatch,"photongenphtmatch[nPhoton][3]/F");
//  Analysis->Branch("photongenelematch",photongenelematch,"photongenelematch[nPhoton][3]/F");
//  Analysis->Branch("photongenphtconv",photongenphtconv,"photongenphtconv[nPhoton][4]/F");
//  Analysis->Branch("photonPartonmatch",photonPartonmatch,"photonPartonmatch[nPhoton][6]/F");
//  Analysis->Branch("photonPartonmatchp",photonPartonmatchp,"photonPartonmatchp[nPhoton][3]/F");
 
 
 Analysis->Branch("photonieta",photonieta,"photonieta[nPhoton]/I");
 Analysis->Branch("photoniphi",photoniphi,"photoniphi[nPhoton]/I");

 Analysis->Branch("photonswissCross",photonswissCross,"photonswissCross[nPhoton]/F");
 Analysis->Branch("photonE2overE9",photonE2overE9,"photonE2overE9[nPhoton]/F");
 

 Analysis->Branch("photoncaloPositionx",photoncaloPositionx,"photoncaloPositionx[nPhoton]/F");
 Analysis->Branch("photoncaloPositiony",photoncaloPositiony,"photoncaloPositiony[nPhoton]/F");
 Analysis->Branch("photoncaloPositionz",photoncaloPositionz,"photoncaloPositionz[nPhoton]/F");
 
 

 ///Photon's conversion
 Analysis->Branch("photonconversionsize",photonconversionsize,"photonconversionsize[nPhoton]/I");
Analysis->Branch("photonconversionVertexx",photonconversionVertexx,"photonconversionVertexx[nPhoton]/F");
Analysis->Branch("photonconversionVertexy",photonconversionVertexy,"photonconversionVertexy[nPhoton]/F");
Analysis->Branch("photonconversionVertexz",photonconversionVertexz,"photonconversionVertexz[nPhoton]/F");
 
 Analysis->Branch("photonconversionrefittedPairMomentumx",photonconversionrefittedPairMomentumx,"photonconversionrefittedPairMomentumx[nPhoton]/F");
 Analysis->Branch("photonconversionrefittedPairMomentumy",photonconversionrefittedPairMomentumy,"photonconversionrefittedPairMomentumy[nPhoton]/F");
 Analysis->Branch("photonconversionrefittedPairMomentumz",photonconversionrefittedPairMomentumz,"photonconversionrefittedPairMomentumz[nPhoton]/F");

Analysis->Branch("photonconversionpairInvariantMass",photonconversionpairInvariantMass,"photonconversionpairInvariantMass[nPhoton]/F");
Analysis->Branch("photonconversionpairCotThetaSeparation",photonconversionpairCotThetaSeparation,"photonconversionpairCotThetaSeparation[nPhoton]/F");
Analysis->Branch("photonconversionEoverPrefittedTracks",photonconversionEoverPrefittedTracks,"photonconversionEoverPrefittedTracks[nPhoton]/F");
Analysis->Branch("photonconversionzOfPrimaryVertexFromTracks",photonconversionzOfPrimaryVertexFromTracks,"photonconversionzOfPrimaryVertexFromTracks[nPhoton]/F");
Analysis->Branch("photonconversiondistOfMinimumApproach",photonconversiondistOfMinimumApproach,"photonconversiondistOfMinimumApproach[nPhoton]/F");
Analysis->Branch("photonconversiondPhiTracksAtVtx",photonconversiondPhiTracksAtVtx,"photonconversiondPhiTracksAtVtx[nPhoton]/F");
Analysis->Branch("photonconversiondPhiTracksAtEcal",photonconversiondPhiTracksAtEcal,"photonconversiondPhiTracksAtEcal[nPhoton]/F");
Analysis->Branch("photonconversiondEtaTracksAtEcal",photonconversiondEtaTracksAtEcal,"photonconversiondEtaTracksAtEcal[nPhoton]/F");
Analysis->Branch("photonconversionnTracks",photonconversionnTracks,"photonconversionnTracks[nPhoton]/I");
Analysis->Branch("photonconversionMVAout",photonconversionMVAout,"photonconversionMVAout[nPhoton]/F");
Analysis->Branch("photonconversionVertexisValid",photonconversionVertexisValid,"photonconversionVertexisValid[nPhoton]/I");
Analysis->Branch("photonconversionVertexchi2",photonconversionVertexchi2,"photonconversionVertexchi2[nPhoton]/F");
Analysis->Branch("photonconversionChiSquaredProbability",photonconversionChiSquaredProbability,"photonconversionChiSquaredProbability[nPhoton]/F");
Analysis->Branch("photonconversion_track1_dz",photonconversion_track1_dz,"photonconversion_track1_dz[nPhoton]/F");
Analysis->Branch("photonconversion_track1_dzError",photonconversion_track1_dzError,"photonconversion_track1_dzError[nPhoton]/F");
Analysis->Branch("photonconversion_track1_charge",photonconversion_track1_charge,"photonconversion_track1_charge[nPhoton]/I");
Analysis->Branch("photonconversion_track1_d0",photonconversion_track1_d0,"photonconversion_track1_d0[nPhoton]/F");
Analysis->Branch("photonconversion_track1_tracksPout",photonconversion_track1_tracksPout,"photonconversion_track1_tracksPout[nPhoton]/F");
Analysis->Branch("photonconversion_track1_tracksPin",photonconversion_track1_tracksPin,"photonconversion_track1_tracksPin[nPhoton]/F");
 Analysis->Branch("photonconversion_track1_algo",photonconversion_track1_algo,"photonconversion_track1_algo[nPhoton]/I");


Analysis->Branch("photonconversion_track2_dz",photonconversion_track2_dz,"photonconversion_track2_dz[nPhoton]/F");
Analysis->Branch("photonconversion_track2_dzError",photonconversion_track2_dzError,"photonconversion_track2_dzError[nPhoton]/F");
Analysis->Branch("photonconversion_track2_charge",photonconversion_track2_charge,"photonconversion_track2_charge[nPhoton]/I");
Analysis->Branch("photonconversion_track2_d0",photonconversion_track2_d0,"photonconversion_track2_d0[nPhoton]/F");
Analysis->Branch("photonconversion_track2_tracksPout",photonconversion_track2_tracksPout,"photonconversion_track2_tracksPout[nPhoton]/F");
Analysis->Branch("photonconversion_track2_tracksPin",photonconversion_track2_tracksPin,"photonconversion_track2_tracksPin[nPhoton]/F");
 Analysis->Branch("photonconversion_track2_algo",photonconversion_track2_algo,"photonconversion_track2_algo[nPhoton]/I");


//more added
 Analysis->Branch("photonscnhits",photonscnhits,"photonscnhits[nPhoton]/I");
 Analysis->Branch("photonscclusterSize",photonscclusterSize,"photonscclusterSize[nPhoton]/I");

 photonscbclusterenergy = new std::vector<std::vector<float> >; photonscbclusterenergy->clear();
 Analysis->Branch("photonscbclusterenergy", "std::vector<std::vector<float> >", &photonscbclusterenergy);
 

 Analysis->Branch("photoneLeft",photoneLeft,"photoneLeft[nPhoton]/F");
 Analysis->Branch("photoneRight",photoneRight,"photoneRight[nPhoton]/F");
 Analysis->Branch("photoneBottom",photoneBottom,"photoneBottom[nPhoton]/F");
 Analysis->Branch("photoneTop",photoneTop,"photoneTop[nPhoton]/F");

Analysis->Branch("photone1x3",photone1x3,"photone1x3[nPhoton]/F");
Analysis->Branch("photone3x1",photone3x1,"photone3x1[nPhoton]/F");
Analysis->Branch("photone2x2",photone2x2,"photone2x2[nPhoton]/F");
Analysis->Branch("photone3x2",photone3x2,"photone3x2[nPhoton]/F");
Analysis->Branch("photone4x4",photone4x4,"photone4x4[nPhoton]/F");
Analysis->Branch("photone2x5Right",photone2x5Right,"photone2x5Right[nPhoton]/F");
Analysis->Branch("photone2x5Left",photone2x5Left,"photone2x5Left[nPhoton]/F");
Analysis->Branch("photone2x5Top",photone2x5Top,"photone2x5Top[nPhoton]/F");
Analysis->Branch("photone2x5Bottom",photone2x5Bottom,"photone2x5Bottom[nPhoton]/F");
Analysis->Branch("photone2x5Max",photone2x5Max,"photone2x5Max[nPhoton]/F");
////Analysis->Branch("photonenergyBasketFractionEta",photonenergyBasketFractionEta,"photonenergyBasketFractionEta[nPhoton]/F");
////Analysis->Branch("photonenergyBasketFractionPhi",photonenergyBasketFractionPhi,"photonenergyBasketFractionPhi[nPhoton]/F");
Analysis->Branch("photonlat",photonlat,"photonlat[nPhoton][3]/F");
Analysis->Branch("photonCovEtaEta",photonCovEtaEta,"photonCovEtaEta[nPhoton]/F");
Analysis->Branch("photonCovEtaPhi",photonCovEtaPhi,"photonCovEtaPhi[nPhoton]/F");
Analysis->Branch("photonCovPhiPhi",photonCovPhiPhi,"photonCovPhiPhi[nPhoton]/F");
////Analysis->Branch("photonCovIEtaIEta",photonCovIEtaIEta,"photonCovIEtaIEta[nPhoton]/F"); //same as sigietaieta saved

Analysis->Branch("photonCovIEtaIPhi",photonCovIEtaIPhi,"photonCovIEtaIPhi[nPhoton]/F");
Analysis->Branch("photonCovIPhiIPhi",photonCovIPhiIPhi,"photonCovIPhiIPhi[nPhoton]/F");

// Analysis->Branch("photonscCovIEtaIEta",photonscCovIEtaIEta,"photonscCovIEtaIEta[nPhoton]/F");
// Analysis->Branch("photonscCovIEtaIPhi",photonscCovIEtaIPhi,"photonscCovIEtaIPhi[nPhoton]/F");
// Analysis->Branch("photonscCovIPhiIPhi",photonscCovIPhiIPhi,"photonscCovIPhiIPhi[nPhoton]/F");

Analysis->Branch("photonzernike20",photonzernike20,"photonzernike20[nPhoton]/F");
Analysis->Branch("photonzernike42",photonzernike42,"photonzernike42[nPhoton]/F");



  
  ///electron
  Analysis->Branch("nElectron",&nElectron, "nElectron/I");
  Analysis->Branch("electronpt",electronpt,"electronpt[nElectron]/F");
  Analysis->Branch("electroneta",electroneta,"electroneta[nElectron]/F");
  Analysis->Branch("electronphi",electronphi,"electronphi[nElectron]/F");
  Analysis->Branch("electroncharge",electroncharge,"electroncharge[nElectron]/I");
  Analysis->Branch("electronvertexx",electronvertexx,"electronvertexx[nElectron]/F");
  Analysis->Branch("electronvertexy",electronvertexy,"electronvertexy[nElectron]/F");
  Analysis->Branch("electronvertexz",electronvertexz,"electronvertexz[nElectron]/F");
  Analysis->Branch("electronscrawEnergy",electronscrawEnergy,"electronscrawEnergy[nElectron]/F");
  Analysis->Branch("electronsceta",electronsceta,"electronsceta[nElectron]/F");
  Analysis->Branch("electronscphi",electronscphi,"electronscphi[nElectron]/F");
  Analysis->Branch("electronscenergy",electronscenergy,"electronscenergy[nElectron]/F");
  
  Analysis->Branch("electronisEcalEnergyCorrected",electronisEcalEnergyCorrected,"electronisEcalEnergyCorrected[nElectron]/I");
  Analysis->Branch("electronecalEnergy",electronecalEnergy,"electronecalEnergy[nElectron]/F");
  Analysis->Branch("electronscpreshowerEnergy",electronscpreshowerEnergy,"electronscpreshowerEnergy[nElectron]/F");
  Analysis->Branch("electronscphiWidth",electronscphiWidth,"electronscphiWidth[nElectron]/F");
  Analysis->Branch("electronscetaWidth",electronscetaWidth,"electronscetaWidth[nElectron]/F");
  Analysis->Branch("electronscnhits",electronscnhits,"electronscnhits[nElectron]/I");
  Analysis->Branch("electronscclusterSize",electronscclusterSize,"electronscclusterSize[nElectron]/I");
  
  electronscbclusterenergy = new std::vector<std::vector<float> >; electronscbclusterenergy->clear();
  Analysis->Branch("electronscbclusterenergy", "std::vector<std::vector<float> >", &electronscbclusterenergy);

  if( saveAllRecHitsSC_){
    electronscrechitenergy = new std::vector<std::vector<float> >; electronscrechitenergy->clear();
    Analysis->Branch("electronscrechitenergy", "std::vector<std::vector<float> >", &electronscrechitenergy);
    electronscrechitfraction = new std::vector<std::vector<float> >; electronscrechitfraction->clear();
    Analysis->Branch("electronscrechitfraction", "std::vector<std::vector<float> >", &electronscrechitfraction);
    

    electronscrechitieta = new std::vector<std::vector<short> >; electronscrechitieta->clear();
    Analysis->Branch("electronscrechitieta", "std::vector<std::vector<short> >", &electronscrechitieta);
    electronscrechitiphi = new std::vector<std::vector<short> >; electronscrechitiphi->clear();
    Analysis->Branch("electronscrechitiphi", "std::vector<std::vector<short> >", &electronscrechitiphi);
    electronscrechitlaserCorr = new std::vector<std::vector<float> >; electronscrechitlaserCorr->clear();
    Analysis->Branch("electronscrechitlaserCorr", "std::vector<std::vector<float> >", &electronscrechitlaserCorr);
  }
  
  Analysis->Branch("electroncaloPositionx",electroncaloPositionx,"electroncaloPositionx[nElectron]/F");
  Analysis->Branch("electroncaloPositiony",electroncaloPositiony,"electroncaloPositiony[nElectron]/F");
  Analysis->Branch("electroncaloPositionz",electroncaloPositionz,"electroncaloPositionz[nElectron]/F");
 
  
  Analysis->Branch("electronfiduficalFlag",&electronfiduficalFlag,"electronfiduficalFlag[nElectron]/I");
  Analysis->Branch("electrontrackerDrivenSeed",electrontrackerDrivenSeed,"electrontrackerDrivenSeed[nElectron]/I");
  Analysis->Branch("electronecalDrivenSeed",electronecalDrivenSeed,"electronecalDrivenSeed[nElectron]/I");
  Analysis->Branch("electronfbrem",electronfbrem,"electronfbrem[nElectron]/F");
  Analysis->Branch("electronnumberOfBrems",electronnumberOfBrems,"electronnumberOfBrems[nElectron]/I");
  Analysis->Branch("electrondr03TkSumPt",electrondr03TkSumPt,"electrondr03TkSumPt[nElectron]/F");
  Analysis->Branch("electrondr03EcalRecHitSumEt",electrondr03EcalRecHitSumEt,"electrondr03EcalRecHitSumEt[nElectron]/F");
  Analysis->Branch("electrondr03HcalDepth1TowerSumEt",electrondr03HcalDepth1TowerSumEt,"electrondr03HcalDepth1TowerSumEt[nElectron]/F");
  Analysis->Branch("electrondr03HcalDepth2TowerSumEt",electrondr03HcalDepth2TowerSumEt,"electrondr03HcalDepth2TowerSumEt[nElectron]/F");
  Analysis->Branch("electrondr03HcalTowerSumEt",electrondr03HcalTowerSumEt,"electrondr03HcalTowerSumEt[nElectron]/F");
  Analysis->Branch("electrondr04EcalRecHitSumEt",electrondr04EcalRecHitSumEt,"electrondr04EcalRecHitSumEt[nElectron]/F");
  Analysis->Branch("electrondr04HcalDepth1TowerSumEt",electrondr04HcalDepth1TowerSumEt,"electrondr04HcalDepth1TowerSumEt[nElectron]/F");
  Analysis->Branch("electrondr04HcalDepth2TowerSumEt",electrondr04HcalDepth2TowerSumEt,"electrondr04HcalDepth2TowerSumEt[nElectron]/F");
  Analysis->Branch("electrondr04HcalTowerSumEt",electrondr04HcalTowerSumEt,"electrondr04HcalTowerSumEt[nElectron]/F");
  Analysis->Branch("electrondr04TkSumPt",electrondr04TkSumPt,"electrondr04TkSumPt[nElectron]/F");
  Analysis->Branch("electronhcalDepth1OverEcal",electronhcalDepth1OverEcal,"electronhcalDepth1OverEcal[nElectron]/F");
  Analysis->Branch("electronhcalDepth2OverEcal",electronhcalDepth2OverEcal,"electronhcalDepth2OverEcal[nElectron]/F");
  Analysis->Branch("electronhcalOverEcal",electronhcalOverEcal,"electronhcalOverEcal[nElectron]/F");
  Analysis->Branch("electrone1x5",electrone1x5,"electrone1x5[nElectron]/F");
  Analysis->Branch("electrone2x5Max",electrone2x5Max,"electrone2x5Max[nElectron]/F");
  Analysis->Branch("electrone5x5",electrone5x5,"electrone5x5[nElectron]/F");
  Analysis->Branch("electronsigmaIetaIeta",electronsigmaIetaIeta,"electronsigmaIetaIeta[nElectron]/F");
  Analysis->Branch("electrone3x3",electrone3x3,"electrone3x3[nElectron]/F");
  Analysis->Branch("electroneMax",electroneMax,"electroneMax[nElectron]/F");
  
  Analysis->Branch("electroneLeft",electroneLeft,"electroneLeft[nElectron]/F");
  Analysis->Branch("electroneRight",electroneRight,"electroneRight[nElectron]/F");
  Analysis->Branch("electroneBottom",electroneBottom,"electroneBottom[nElectron]/F");
  Analysis->Branch("electroneTop",electroneTop,"electroneTop[nElectron]/F");

Analysis->Branch("electrone1x3",electrone1x3,"electrone1x3[nElectron]/F");
Analysis->Branch("electrone3x1",electrone3x1,"electrone3x1[nElectron]/F");
Analysis->Branch("electrone2x2",electrone2x2,"electrone2x2[nElectron]/F");
Analysis->Branch("electrone3x2",electrone3x2,"electrone3x2[nElectron]/F");
Analysis->Branch("electrone4x4",electrone4x4,"electrone4x4[nElectron]/F");
Analysis->Branch("electrone2x5Right",electrone2x5Right,"electrone2x5Right[nElectron]/F");
Analysis->Branch("electrone2x5Left",electrone2x5Left,"electrone2x5Left[nElectron]/F");
Analysis->Branch("electrone2x5Top",electrone2x5Top,"electrone2x5Top[nElectron]/F");
Analysis->Branch("electrone2x5Bottom",electrone2x5Bottom,"electrone2x5Bottom[nElectron]/F");
////Analysis->Branch("electrone2x5Max",electrone2x5Max,"electrone2x5Max[nElectron]/F");
Analysis->Branch("electronlat",electronlat,"electronlat[nElectron][3]/F");
Analysis->Branch("electronCovEtaEta",electronCovEtaEta,"electronCovEtaEta[nElectron]/F");
Analysis->Branch("electronCovEtaPhi",electronCovEtaPhi,"electronCovEtaPhi[nElectron]/F");
Analysis->Branch("electronCovPhiPhi",electronCovPhiPhi,"electronCovPhiPhi[nElectron]/F");
Analysis->Branch("electronCovIEtaIPhi",electronCovIEtaIPhi,"electronCovIEtaIPhi[nElectron]/F");
Analysis->Branch("electronCovIPhiIPhi",electronCovIPhiIPhi,"electronCovIPhiIPhi[nElectron]/F");
Analysis->Branch("electronzernike20",electronzernike20,"electronzernike20[nElectron]/F");
Analysis->Branch("electronzernike42",electronzernike42,"electronzernike42[nElectron]/F");




  Analysis->Branch("electroneSuperClusterOverP",electroneSuperClusterOverP,"electroneSuperClusterOverP[nElectron]/F");
  Analysis->Branch("electroneSeedClusterOverP",electroneSeedClusterOverP,"electroneSeedClusterOverP[nElectron]/F");
  Analysis->Branch("electroneSeedClusterOverPout",electroneSeedClusterOverPout,"electroneSeedClusterOverPout[nElectron]/F");
  Analysis->Branch("electroneEleClusterOverPout",electroneEleClusterOverPout,"electroneEleClusterOverPout[nElectron]/F");
  Analysis->Branch("electrondeltaEtaSuperClusterTrackAtVtx",electrondeltaEtaSuperClusterTrackAtVtx,"electrondeltaEtaSuperClusterTrackAtVtx[nElectron]/F");
  Analysis->Branch("electrondeltaEtaSeedClusterTrackAtCalo",electrondeltaEtaSeedClusterTrackAtCalo,"electrondeltaEtaSeedClusterTrackAtCalo[nElectron]/F");
  Analysis->Branch("electrondeltaEtaEleClusterTrackAtCalo",electrondeltaEtaEleClusterTrackAtCalo,"electrondeltaEtaEleClusterTrackAtCalo[nElectron]/F");
  Analysis->Branch("electrondeltaPhiSuperClusterTrackAtVtx",electrondeltaPhiSuperClusterTrackAtVtx,"electrondeltaPhiSuperClusterTrackAtVtx[nElectron]/F");
  Analysis->Branch("electrondeltaPhiSeedClusterTrackAtCalo",electrondeltaPhiSeedClusterTrackAtCalo,"electrondeltaPhiSeedClusterTrackAtCalo[nElectron]/F");
  Analysis->Branch("electrondeltaPhiEleClusterTrackAtCalo",electrondeltaPhiEleClusterTrackAtCalo,"electrondeltaPhiEleClusterTrackAtCalo[nElectron]/F");
  Analysis->Branch("electronclassification",electronclassification,"electronclassification[nElectron]/I");
  Analysis->Branch("electronmva",electronmva,"electronmva[nElectron]/F");
  Analysis->Branch("electronnumberOfTracks",electronnumberOfTracks,"electronnumberOfTracks[nElectron]/I");
  Analysis->Branch("electronconvDist",electronconvDist,"electronconvDist[nElectron]/F");
  Analysis->Branch("electronconvDcot",electronconvDcot,"electronconvDcot[nElectron]/F");
  Analysis->Branch("electronconvRadius",electronconvRadius,"electronconvRadius[nElectron]/F");
  Analysis->Branch("electronconvFlags",electronconvFlags,"electronconvFlags[nElectron]/I");
  //Analysis->Branch("electronExpectednumberOfHits",electronExpectednumberOfHits,"electronExpectednumberOfHits[nElectron]/I");
  Analysis->Branch("electronExpectedHitsInnernumberOfHits",electronExpectedHitsInnernumberOfHits,"electronExpectedHitsInnernumberOfHits[nElectron]/I");
  Analysis->Branch("electronExpectedHitsOuternumberOfHits",electronExpectedHitsOuternumberOfHits,"electronExpectedHitsOuternumberOfHits[nElectron]/I");

  Analysis->Branch("electrongsfTrackvx",electrongsfTrackvx,"electrongsfTrackvx[nElectron]/F");
  Analysis->Branch("electrongsfTrackvy",electrongsfTrackvy,"electrongsfTrackvy[nElectron]/F");
  Analysis->Branch("electrongsfTrackvz",electrongsfTrackvz,"electrongsfTrackvz[nElectron]/F");
  Analysis->Branch("electrongsfTracknormalizedChi2",electrongsfTracknormalizedChi2,"electrongsfTracknormalizedChi2[nElectron]/F");
  Analysis->Branch("electrongsfTrackdxybeamSpot",electrongsfTrackdxybeamSpot,"electrongsfTrackdxybeamSpot[nElectron]/F");
  Analysis->Branch("electrongsfTrackdzbeamSpot",electrongsfTrackdzbeamSpot,"electrongsfTrackdzbeamSpot[nElectron]/F");
  Analysis->Branch("electrongsfTracknumberOfValidHits",electrongsfTracknumberOfValidHits,"electrongsfTracknumberOfValidHits[nElectron]/I");
  Analysis->Branch("electrongsfTracknumberOfLostHits",electrongsfTracknumberOfLostHits,"electrongsfTracknumberOfLostHits[nElectron]/I");
  Analysis->Branch("electrongsfTracknumberOfValidPixelHits",electrongsfTracknumberOfValidPixelHits,"electrongsfTracknumberOfValidPixelHits[nElectron]/I");
  Analysis->Branch("electrongsfTracknumberOfValidTrackerHits",electrongsfTracknumberOfValidTrackerHits,"electrongsfTracknumberOfValidTrackerHits[nElectron]/I");
  Analysis->Branch("electrongsfTrackpt",electrongsfTrackpt,"electrongsfTrackpt[nElectron]/F");
  Analysis->Branch("electrongsfTracketa",electrongsfTracketa,"electrongsfTracketa[nElectron]/F");
  Analysis->Branch("electrongsfTrackphi",electrongsfTrackphi,"electrongsfTrackphi[nElectron]/F");

  

  if(saveLaserCorrSeedCrystal_){
    Analysis->Branch("electronseedlaserCorr",electronseedlaserCorr,"electronseedlaserCorr[nElectron]/F");
    Analysis->Branch("electronaverlaserCorr",electronaverlaserCorr,"electronaverlaserCorr[nElectron]/F");
  }
  
  Analysis->Branch("electronseedtime",electronseedtime,"electronseedtime[nElectron]/F");
  Analysis->Branch("electronseedoutOfTimeChi2",electronseedoutOfTimeChi2,"electronseedoutOfTimeChi2[nElectron]/F");
  Analysis->Branch("electronseedchi2",electronseedchi2,"electronseedchi2[nElectron]/F");
  Analysis->Branch("electronseedrecoFlag",electronseedrecoFlag,"electronseedrecoFlag[nElectron]/I");
  Analysis->Branch("electronseedseverityLevel",electronseedseverityLevel,"electronseedseverityLevel[nElectron]/I");
  Analysis->Branch("electronscindex",electronscindex,"electronscindex[nElectron]/I");

  Analysis->Branch("electronieta",electronieta,"electronieta[nElectron]/I");
  Analysis->Branch("electroniphi",electroniphi,"electroniphi[nElectron]/I");
  Analysis->Branch("electronswissCross",electronswissCross,"electronswissCross[nElectron]/F");
  Analysis->Branch("electronE2overE9",electronE2overE9,"electronE2overE9[nElectron]/F");
  
  
  ///photonp4 = new TClonesArray("TLorentzVector", nMaxPhoton);
  ///Analysis->Branch("photonp4", "TClonesArray", &photonp4, 32000, 0);
    
  ///Muon
  Analysis->Branch("nMuon",&nMuon,"nMuon/I");
  ///muonp4 = new TClonesArray("TLorentzVector", nMuonMAX);  this take lots of space.
  ////Analysis->Branch("muonp4", "TClonesArray", &muonp4, 32000, 0);
  
  Analysis->Branch("muonpt",muonpt,"muonpt[nMuon]/F");
  Analysis->Branch("muoncharge",muoncharge,"muoncharge[nMuon]/I");
  Analysis->Branch("muoneta",muoneta,"muoneta[nMuon]/F");
  Analysis->Branch("muonphi",muonphi,"muonphi[nMuon]/F");
  
  Analysis->Branch("muonRecoAlgo",muonRecoAlgo,"muonRecoAlgo[nMuon]/I");
  Analysis->Branch("muonisGoodMuon",muonisGoodMuon,"muonisGoodMuon[nMuon]/I");
  Analysis->Branch("muonvx",muonvx,"muonvx[nMuon]/F");
  Analysis->Branch("muonvy",muonvy,"muonvy[nMuon]/F");
  Analysis->Branch("muonvz",muonvz,"muonvz[nMuon]/F");
  Analysis->Branch("muonnumberOfMatches",muonnumberOfMatches,"muonnumberOfMatches[nMuon]/I");
  Analysis->Branch("muonnumberOfChambers",muonnumberOfChambers,"muonnumberOfChambers[nMuon]/I");
  Analysis->Branch("muonglobalTrackvx",muonglobalTrackvx,"muonglobalTrackvx[nMuon]/F");
  Analysis->Branch("muonglobalTrackvy",muonglobalTrackvy,"muonglobalTrackvy[nMuon]/F");
  Analysis->Branch("muonglobalTrackvz",muonglobalTrackvz,"muonglobalTrackvz[nMuon]/F");
  Analysis->Branch("muonglobalTrackdzbeamSpot",muonglobalTrackdzbeamSpot,"muonglobalTrackdzbeamSpot[nMuon]/F");
  Analysis->Branch("muonglobalTrackdxybeamSpot",muonglobalTrackdxybeamSpot,"muonglobalTrackdxybeamSpot[nMuon]/F");
  Analysis->Branch("muonglobalTracknumberOfValidPixelHits",muonglobalTracknumberOfValidPixelHits,"muonglobalTracknumberOfValidPixelHits[nMuon]/I");
  Analysis->Branch("muonglobalTracknumberOfValidTrackerHits",muonglobalTracknumberOfValidTrackerHits,"muonglobalTracknumberOfValidTrackerHits[nMuon]/I");
  Analysis->Branch("muonglobalTracknormalizedChi2",muonglobalTracknormalizedChi2,"muonglobalTracknormalizedChi2[nMuon]/F");
  Analysis->Branch("muonglobalTracknumberOfValidMuonHits",muonglobalTracknumberOfValidMuonHits,"muonglobalTracknumberOfValidMuonHits[nMuon]/I");
  Analysis->Branch("muonisolationR03sumPt",muonisolationR03sumPt,"muonisolationR03sumPt[nMuon]/F");
  Analysis->Branch("muonisolationR03emEt",muonisolationR03emEt,"muonisolationR03emEt[nMuon]/F");
  Analysis->Branch("muonisolationR03hadEt",muonisolationR03hadEt,"muonisolationR03hadEt[nMuon]/F");
  Analysis->Branch("muonisolationR03hoEt",muonisolationR03hoEt,"muonisolationR03hoEt[nMuon]/F");
  Analysis->Branch("muonisolationR04sumPt",muonisolationR04sumPt,"muonisolationR04sumPt[nMuon]/F");
  Analysis->Branch("muonisolationR04emEt",muonisolationR04emEt,"muonisolationR04emEt[nMuon]/F");
  Analysis->Branch("muonisolationR04hadEt",muonisolationR04hadEt,"muonisolationR04hadEt[nMuon]/F");
  Analysis->Branch("muonisolationR04hoEt",muonisolationR04hoEt,"muonisolationR04hoEt[nMuon]/F");
  Analysis->Branch("muonisolationR05sumPt",muonisolationR05sumPt,"muonisolationR05sumPt[nMuon]/F");
  Analysis->Branch("muonisolationR05emEt",muonisolationR05emEt,"muonisolationR05emEt[nMuon]/F");
  Analysis->Branch("muonisolationR05hadEt",muonisolationR05hadEt,"muonisolationR05hadEt[nMuon]/F");
  Analysis->Branch("muonisolationR05hoEt",muonisolationR05hoEt,"muonisolationR05hoEt[nMuon]/F");
  
  Analysis->Branch("muoninnerTracknumberOfValidPixelHits",muoninnerTracknumberOfValidPixelHits,"muoninnerTracknumberOfValidPixelHits[nMuon]/I");
  Analysis->Branch("muoninnerTracknumberOfValidTrackerHits",muoninnerTracknumberOfValidTrackerHits,"muoninnerTracknumberOfValidTrackerHits[nMuon]/I");
  Analysis->Branch("muoninnerTrackpt",muoninnerTrackpt,"muoninnerTrackpt[nMuon]/F");
Analysis->Branch("muoninnerTracketa",muoninnerTracketa,"muoninnerTracketa[nMuon]/F");
  Analysis->Branch("muoninnerTrackphi",muoninnerTrackphi,"muoninnerTrackphi[nMuon]/F");
  
  Analysis->Branch("muonouterTracknumberOfValidMuonHits",muonouterTracknumberOfValidMuonHits,"muonouterTracknumberOfValidMuonHits[nMuon]/I");
  Analysis->Branch("muonouterTrackpt",muonouterTrackpt,"muonouterTrackpt[nMuon]/F");
  Analysis->Branch("muonouterTracketa",muonouterTracketa,"muonouterTracketa[nMuon]/F");
  Analysis->Branch("muonouterTrackphi",muonouterTrackphi,"muonouterTrackphi[nMuon]/F");
  
  
  ///vertex
  Analysis->Branch("nVertex",&nVertex,"nVertex/I");
  Analysis->Branch("vertexx",vertexx,"vertexx[nVertex]/F");
  Analysis->Branch("vertexy",vertexy,"vertexy[nVertex]/F");
  Analysis->Branch("vertexz",vertexz,"vertexz[nVertex]/F");
  Analysis->Branch("vertexchi2",vertexchi2,"vertexchi2[nVertex]/F");
  Analysis->Branch("vertexndof",vertexndof,"vertexndof[nVertex]/F");
  Analysis->Branch("vertexnormalizedChi2",vertexnormalizedChi2,"vertexnormalizedChi2[nVertex]/F");
  Analysis->Branch("vertextrackSize",vertextrackSize,"vertextrackSize[nVertex]/I");
  Analysis->Branch("vertexisFake",vertexisFake,"vertexisFake[nVertex]/I");
  Analysis->Branch("vertexisValid",vertexisValid,"vertexisValid[nVertex]/I");
  Analysis->Branch("vertexsumtrackspt",vertexsumtrackspt,"vertexsumtrackspt[nVertex]/F");
  Analysis->Branch("vertexsumtracksptSquare",vertexsumtracksptSquare,"vertexsumtracksptSquare[nVertex]/F");
  Analysis->Branch("vertexsumRefittedtrackspt",vertexsumRefittedtrackspt,"vertexsumRefittedtrackspt[nVertex]/F");
  Analysis->Branch("vertexsumRefittedtracksptSquare",vertexsumRefittedtracksptSquare,"vertexsumRefittedtracksptSquare[nVertex]/F");
  
  vertex_trkind = new std::vector<std::vector<short> >; vertex_trkind->clear();
  Analysis->Branch("vertex_trkind", "std::vector<std::vector<short> >", &vertex_trkind);
  vertex_trkWeight = new std::vector<std::vector<float> >; vertex_trkWeight->clear();
  Analysis->Branch("vertex_trkWeight", "std::vector<std::vector<float> >", &vertex_trkWeight);
  

  ///vetex nobs
  Analysis->Branch("nVertexNoBS",&nVertexNoBS,"nVertexNoBS/I");
  Analysis->Branch("vertexNoBSx",vertexNoBSx,"vertexNoBSx[nVertexNoBS]/F");
  Analysis->Branch("vertexNoBSy",vertexNoBSy,"vertexNoBSy[nVertexNoBS]/F");
  Analysis->Branch("vertexNoBSz",vertexNoBSz,"vertexNoBSz[nVertexNoBS]/F");
  Analysis->Branch("vertexNoBSchi2",vertexNoBSchi2,"vertexNoBSchi2[nVertexNoBS]/F");
  Analysis->Branch("vertexNoBSndof",vertexNoBSndof,"vertexNoBSndof[nVertexNoBS]/F");
  Analysis->Branch("vertexNoBSnormalizedChi2",vertexNoBSnormalizedChi2,"vertexNoBSnormalizedChi2[nVertexNoBS]/F");
  Analysis->Branch("vertexNoBStrackSize",vertexNoBStrackSize,"vertexNoBStrackSize[nVertexNoBS]/I");
  Analysis->Branch("vertexNoBSisFake",vertexNoBSisFake,"vertexNoBSisFake[nVertexNoBS]/I");
  Analysis->Branch("vertexNoBSisValid",vertexNoBSisValid,"vertexNoBSisValid[nVertexNoBS]/I");
  
  
  //physics declared
  Analysis->Branch("phyDeclared",&phyDeclared,"phyDeclared/I");

  Analysis->Branch("nTrack",&nTrack,"nTrack/I");
  Analysis->Branch("trackhighPurityFraction",&trackhighPurityFraction,"trackhighPurityFraction/F");
  if(saveTrk_){
    //track
  Analysis->Branch("trackd0",trackd0,"trackd0[nTrack]/F");
  Analysis->Branch("trackd0Error",trackd0Error,"trackd0Error[nTrack]/F");
  Analysis->Branch("trackdz",trackdz,"trackdz[nTrack]/F");
  Analysis->Branch("trackdzError",trackdzError,"trackdzError[nTrack]/F");
  Analysis->Branch("trackpt",trackpt,"trackpt[nTrack]/F");
  Analysis->Branch("trackcharge",trackcharge,"trackcharge[nTrack]/I");
  Analysis->Branch("trackpx",trackpx,"trackpx[nTrack]/F");
  Analysis->Branch("trackpy",trackpy,"trackpy[nTrack]/F");
  Analysis->Branch("trackpz",trackpz,"trackpz[nTrack]/F");
  Analysis->Branch("trackptError",trackptError,"trackptError[nTrack]/F");
  Analysis->Branch("tracketa",tracketa,"tracketa[nTrack]/F");
  Analysis->Branch("trackphi",trackphi,"trackphi[nTrack]/F");
  Analysis->Branch("trackvx",trackvx,"trackvx[nTrack]/F");
  Analysis->Branch("trackvy",trackvy,"trackvy[nTrack]/F");
  Analysis->Branch("trackvz",trackvz,"trackvz[nTrack]/F");
  Analysis->Branch("tracknormalizedChi2",tracknormalizedChi2,"tracknormalizedChi2[nTrack]/F");
  Analysis->Branch("tracknumberOfValidHits",tracknumberOfValidHits,"tracknumberOfValidHits[nTrack]/I");
  Analysis->Branch("trackndof",trackndof,"trackndof[nTrack]/I");
  Analysis->Branch("tracknumberOfValidPixelHits",tracknumberOfValidPixelHits,"tracknumberOfValidPixelHits[nTrack]/I");
  Analysis->Branch("tracknumberOfValidStripHits",tracknumberOfValidStripHits,"tracknumberOfValidStripHits[nTrack]/I");
  Analysis->Branch("tracknumberOfValidTrackerHits",tracknumberOfValidTrackerHits,"tracknumberOfValidTrackerHits[nTrack]/I");
  Analysis->Branch("trackalgo",trackalgo,"trackalgo[nTrack]/I");
  Analysis->Branch("trackqualityFlagTracks",trackqualityFlagTracks,"trackqualityFlagTracks[nTrack]/I");
  }
  
  
  ///All Conversions
  Analysis->Branch("nConv",&nConv,"nConv/I");
  Analysis->Branch("convrefittedPair4Momentumeta",convrefittedPair4Momentumeta,"convrefittedPair4Momentumeta[nConv]/F");
  Analysis->Branch("convrefittedPair4Momentumphi",convrefittedPair4Momentumphi,"convrefittedPair4Momentumphi[nConv]/F");
  Analysis->Branch("convrefittedPair4Momentumpt",convrefittedPair4Momentumpt,"convrefittedPair4Momentumpt[nConv]/F");
  Analysis->Branch("convrefittedPair4Momentumenergy",convrefittedPair4Momentumenergy,"convrefittedPair4Momentumenergy[nConv]/F");
  Analysis->Branch("convnTracks",convnTracks,"convnTracks[nConv]/I");
  Analysis->Branch("convcaloClustersize",convcaloClustersize,"convcaloClustersize[nConv]/I");
  
  Analysis->Branch("convcaloCluster0x",convcaloCluster0x,"convcaloCluster0x[nConv]/F");
  Analysis->Branch("convcaloCluster0y",convcaloCluster0y,"convcaloCluster0y[nConv]/F");
  Analysis->Branch("convcaloCluster0z",convcaloCluster0z,"convcaloCluster0z[nConv]/F");
  
  
  Analysis->Branch("convconversionVertexisValid",convconversionVertexisValid,"convconversionVertexisValid[nConv]/I");
  Analysis->Branch("convpairMomentumx",convpairMomentumx,"convpairMomentumx[nConv]/F");
  Analysis->Branch("convpairMomentumy",convpairMomentumy,"convpairMomentumy[nConv]/F");
  Analysis->Branch("convpairMomentumz",convpairMomentumz,"convpairMomentumz[nConv]/F");
  Analysis->Branch("convrefittedPairMomentumx",convrefittedPairMomentumx,"convrefittedPairMomentumx[nConv]/F");
  Analysis->Branch("convrefittedPairMomentumy",convrefittedPairMomentumy,"convrefittedPairMomentumy[nConv]/F");
  Analysis->Branch("convrefittedPairMomentumz",convrefittedPairMomentumz,"convrefittedPairMomentumz[nConv]/F");
  Analysis->Branch("convconversionVertexx",convconversionVertexx,"convconversionVertexx[nConv]/F");
  Analysis->Branch("convconversionVertexy",convconversionVertexy,"convconversionVertexy[nConv]/F");
  Analysis->Branch("convconversionVertexz",convconversionVertexz,"convconversionVertexz[nConv]/F");
  Analysis->Branch("convconversionVertexchi2",convconversionVertexchi2,"convconversionVertexchi2[nConv]/F");
  Analysis->Branch("convconversionVertexChiSquaredProbability",convconversionVertexChiSquaredProbability,"convconversionVertexChiSquaredProbability[nConv]/F");
  Analysis->Branch("convconversionVertexxError",convconversionVertexxError,"convconversionVertexxError[nConv]/F");
  Analysis->Branch("convconversionVertexyError",convconversionVertexyError,"convconversionVertexyError[nConv]/F");
  Analysis->Branch("convconversionVertexzError",convconversionVertexzError,"convconversionVertexzError[nConv]/F");
  Analysis->Branch("convconversionVertexnTracks",convconversionVertexnTracks,"convconversionVertexnTracks[nConv]/F");
  Analysis->Branch("convconversionVertexMVAout",convconversionVertexMVAout,"convconversionVertexMVAout[nConv]/F");
  Analysis->Branch("conv_track1_dz",conv_track1_dz,"conv_track1_dz[nConv]/F");
  Analysis->Branch("conv_track1_dzError",conv_track1_dzError,"conv_track1_dzError[nConv]/F");
  Analysis->Branch("conv_track1_charge",conv_track1_charge,"conv_track1_charge[nConv]/I");
  Analysis->Branch("conv_track2_dz",conv_track2_dz,"conv_track2_dz[nConv]/F");
  Analysis->Branch("conv_track2_dzError",conv_track2_dzError,"conv_track2_dzError[nConv]/F");
  Analysis->Branch("conv_track2_charge",conv_track2_charge,"conv_track2_charge[nConv]/I");
  Analysis->Branch("convpairInvariantMass",convpairInvariantMass,"convpairInvariantMass[nConv]/F");
  Analysis->Branch("convpairCotThetaSeparation",convpairCotThetaSeparation,"convpairCotThetaSeparation[nConv]/F");
  Analysis->Branch("convEoverPrefittedTracks",convEoverPrefittedTracks,"convEoverPrefittedTracks[nConv]/F");
  Analysis->Branch("convzOfPrimaryVertexFromTracks",convzOfPrimaryVertexFromTracks,"convzOfPrimaryVertexFromTracks[nConv]/F");
  Analysis->Branch("convdistOfMinimumApproach",convdistOfMinimumApproach,"convdistOfMinimumApproach[nConv]/F");
  Analysis->Branch("convdPhiTracksAtVtx",convdPhiTracksAtVtx,"convdPhiTracksAtVtx[nConv]/F");
  Analysis->Branch("convdPhiTracksAtEcal",convdPhiTracksAtEcal,"convdPhiTracksAtEcal[nConv]/F");
  Analysis->Branch("convdEtaTracksAtEcal",convdEtaTracksAtEcal,"convdEtaTracksAtEcal[nConv]/F");
  
  Analysis->Branch("convnSharedHits",convnSharedHits,"convnSharedHits[nConv]/I");
  Analysis->Branch("conv_track1_d0",conv_track1_d0,"conv_track1_d0[nConv]/F");
  Analysis->Branch("conv_track1_pout",conv_track1_pout,"conv_track1_pout[nConv]/F");
  Analysis->Branch("conv_track1_pin",conv_track1_pin,"conv_track1_pin[nConv]/F");
  Analysis->Branch("conv_track1_algo",conv_track1_algo,"conv_track1_algo[nConv]/I");
  Analysis->Branch("conv_track2_d0",conv_track2_d0,"conv_track2_d0[nConv]/F");
  Analysis->Branch("conv_track2_pout",conv_track2_pout,"conv_track2_pout[nConv]/F");
  Analysis->Branch("conv_track2_pin",conv_track2_pin,"conv_track2_pin[nConv]/F");
  Analysis->Branch("conv_track2_algo",conv_track2_algo,"conv_track2_algo[nConv]/I");

  Analysis->Branch("conv_track1_pt",conv_track1_pt,"conv_track1_pt[nConv]/F");
  Analysis->Branch("conv_track1_eta",conv_track1_eta,"conv_track1_eta[nConv]/F");
  Analysis->Branch("conv_track1_phi",conv_track1_phi,"conv_track1_phi[nConv]/F");
  Analysis->Branch("conv_track2_pt",conv_track2_pt,"conv_track2_pt[nConv]/F");
  Analysis->Branch("conv_track2_eta",conv_track2_eta,"conv_track2_eta[nConv]/F");
  Analysis->Branch("conv_track2_phi",conv_track2_phi,"conv_track2_phi[nConv]/F");
  

  
  convnHitsBeforeVtx = new std::vector<std::vector<unsigned short> >; convnHitsBeforeVtx->clear();
  Analysis->Branch("convnHitsBeforeVtx", "std::vector<std::vector<unsigned short> >", &convnHitsBeforeVtx);
    
  ///beamspot
  Analysis->Branch("beamSpotX",&beamSpotX,"beamSpotX/F");
  Analysis->Branch("beamSpotY",&beamSpotY,"beamSpotY/F");
  Analysis->Branch("beamSpotZ",&beamSpotZ,"beamSpotZ/F");
  
  
  
  ///jet
  Analysis->Branch("nak5CaloJet",&nak5CaloJet,"nak5CaloJet/I");
  Analysis->Branch("ak5CaloJetet",ak5CaloJetet,"ak5CaloJetet[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJeteta",ak5CaloJeteta,"ak5CaloJeteta[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetphi",ak5CaloJetphi,"ak5CaloJetphi[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetenergyFractionHadronic",ak5CaloJetenergyFractionHadronic,"ak5CaloJetenergyFractionHadronic[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetn90",ak5CaloJetn90,"ak5CaloJetn90[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetn60",ak5CaloJetn60,"ak5CaloJetn60[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJettowersArea",ak5CaloJettowersArea,"ak5CaloJettowersArea[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetCalosize",ak5CaloJetCalosize,"ak5CaloJetCalosize[nak5CaloJet]/I");
  Analysis->Branch("ak5CaloJetfHPD",ak5CaloJetfHPD,"ak5CaloJetfHPD[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetfRBX",ak5CaloJetfRBX,"ak5CaloJetfRBX[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetn90Hits",ak5CaloJetn90Hits,"ak5CaloJetn90Hits[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetrestrictedEMF",ak5CaloJetrestrictedEMF,"ak5CaloJetrestrictedEMF[nak5CaloJet]/F");
  Analysis->Branch("ak5CaloJetcorrection",ak5CaloJetcorrection,"ak5CaloJetcorrection[nak5CaloJet]/F");
  
  
    
  ///jet
  Analysis->Branch("nak7CaloJet",&nak7CaloJet,"nak7CaloJet/I");
  Analysis->Branch("ak7CaloJetet",ak7CaloJetet,"ak7CaloJetet[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJeteta",ak7CaloJeteta,"ak7CaloJeteta[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetphi",ak7CaloJetphi,"ak7CaloJetphi[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetenergyFractionHadronic",ak7CaloJetenergyFractionHadronic,"ak7CaloJetenergyFractionHadronic[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetn90",ak7CaloJetn90,"ak7CaloJetn90[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetn60",ak7CaloJetn60,"ak7CaloJetn60[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJettowersArea",ak7CaloJettowersArea,"ak7CaloJettowersArea[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetCalosize",ak7CaloJetCalosize,"ak7CaloJetCalosize[nak7CaloJet]/I");
  Analysis->Branch("ak7CaloJetfHPD",ak7CaloJetfHPD,"ak7CaloJetfHPD[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetfRBX",ak7CaloJetfRBX,"ak7CaloJetfRBX[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetn90Hits",ak7CaloJetn90Hits,"ak7CaloJetn90Hits[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetrestrictedEMF",ak7CaloJetrestrictedEMF,"ak7CaloJetrestrictedEMF[nak7CaloJet]/F");
  Analysis->Branch("ak7CaloJetcorrection",ak7CaloJetcorrection,"ak7CaloJetcorrection[nak7CaloJet]/F");

  
  
  
  Analysis->Branch("nak5PFJet",&nak5PFJet,"nak5PFJet/I");
  Analysis->Branch("ak5PFJetet",ak5PFJetet,"ak5PFJetet[nak5PFJet]/F");
  Analysis->Branch("ak5PFJeteta",ak5PFJeteta,"ak5PFJeteta[nak5PFJet]/F");
  Analysis->Branch("ak5PFJetphi",ak5PFJetphi,"ak5PFJetphi[nak5PFJet]/F");
  Analysis->Branch("ak5PFJetchargedHadronEnergyFraction",ak5PFJetchargedHadronEnergyFraction,"ak5PFJetchargedHadronEnergyFraction[nak5PFJet]/F");
  Analysis->Branch("ak5PFJetchargedMuEnergyFraction",ak5PFJetchargedMuEnergyFraction,"ak5PFJetchargedMuEnergyFraction[nak5PFJet]/F");
  Analysis->Branch("ak5PFJetneutralEmEnergyFraction",ak5PFJetneutralEmEnergyFraction,"ak5PFJetneutralEmEnergyFraction[nak5PFJet]/F");
  Analysis->Branch("ak5PFJetchargedMultiplicity",ak5PFJetchargedMultiplicity,"ak5PFJetchargedMultiplicity[nak5PFJet]/I");
  Analysis->Branch("ak5PFJetneutralMultiplicity",ak5PFJetneutralMultiplicity,"ak5PFJetneutralMultiplicity[nak5PFJet]/I");
  Analysis->Branch("ak5PFJetPFsize",ak5PFJetPFsize,"ak5PFJetPFsize[nak5PFJet]/I");

  
//   Analysis->Branch("ak5PFJetfHPD",ak5PFJetfHPD,"ak5PFJetfHPD[nak5PFJet]/F");
//   Analysis->Branch("ak5PFJetfRBX",ak5PFJetfRBX,"ak5PFJetfRBX[nak5PFJet]/F");
//   Analysis->Branch("ak5PFJetn90Hits",ak5PFJetn90Hits,"ak5PFJetn90Hits[nak5PFJet]/F");
//   Analysis->Branch("ak5PFJetrestrictedEMF",ak5PFJetrestrictedEMF,"ak5PFJetrestrictedEMF[nak5PFJet]/F");
  Analysis->Branch("ak5PFJetcorrection",ak5PFJetcorrection,"ak5PFJetcorrection[nak5PFJet]/F");

    
  
  
  Analysis->Branch("nak7PFJet",&nak7PFJet,"nak7PFJet/I");
  Analysis->Branch("ak7PFJetet",ak7PFJetet,"ak7PFJetet[nak7PFJet]/F");
  Analysis->Branch("ak7PFJeteta",ak7PFJeteta,"ak7PFJeteta[nak7PFJet]/F");
  Analysis->Branch("ak7PFJetphi",ak7PFJetphi,"ak7PFJetphi[nak7PFJet]/F");
  Analysis->Branch("ak7PFJetchargedHadronEnergyFraction",ak7PFJetchargedHadronEnergyFraction,"ak7PFJetchargedHadronEnergyFraction[nak7PFJet]/F");
  Analysis->Branch("ak7PFJetchargedMuEnergyFraction",ak7PFJetchargedMuEnergyFraction,"ak7PFJetchargedMuEnergyFraction[nak7PFJet]/F");
  Analysis->Branch("ak7PFJetneutralEmEnergyFraction",ak7PFJetneutralEmEnergyFraction,"ak7PFJetneutralEmEnergyFraction[nak7PFJet]/F");
  Analysis->Branch("ak7PFJetchargedMultiplicity",ak7PFJetchargedMultiplicity,"ak7PFJetchargedMultiplicity[nak7PFJet]/I");
  Analysis->Branch("ak7PFJetneutralMultiplicity",ak7PFJetneutralMultiplicity,"ak7PFJetneutralMultiplicity[nak7PFJet]/I");
  Analysis->Branch("ak7PFJetPFsize",ak7PFJetPFsize,"ak7PFJetPFsize[nak7PFJet]/I");
  
//   Analysis->Branch("ak7PFJetfHPD",ak7PFJetfHPD,"ak7PFJetfHPD[nak7PFJet]/F");
//   Analysis->Branch("ak7PFJetfRBX",ak7PFJetfRBX,"ak7PFJetfRBX[nak7PFJet]/F");
//   Analysis->Branch("ak7PFJetn90Hits",ak7PFJetn90Hits,"ak7PFJetn90Hits[nak7PFJet]/F");
//   Analysis->Branch("ak7PFJetrestrictedEMF",ak7PFJetrestrictedEMF,"ak7PFJetrestrictedEMF[nak7PFJet]/F");
  Analysis->Branch("ak7PFJetcorrection",ak7PFJetcorrection,"ak7PFJetcorrection[nak7PFJet]/F");
  /////  Analysis->Branch("ak7PFJetPartonmatch",ak7PFJetPartonmatch,"ak7PFJetPartonmatch[nak7PFJet][6]/F");
  
  

  Analysis->Branch("nL1EMIso",&nL1EMIso,"nL1EMIso/I");
  Analysis->Branch("L1EMIso_e",L1EMIso_e,"L1EMIso_e[nL1EMIso]/F");
  Analysis->Branch("L1EMIso_et",L1EMIso_et,"L1EMIso_et[nL1EMIso]/F");
  Analysis->Branch("L1EMIso_eta",L1EMIso_eta,"L1EMIso_eta[nL1EMIso]/F");
  Analysis->Branch("L1EMIso_phi",L1EMIso_phi,"L1EMIso_phi[nL1EMIso]/F");
  // L1 EM non-isolated objects
  Analysis->Branch("nL1EMnonIso",&nL1EMnonIso,"nL1EMnonIso/I");
  Analysis->Branch("L1EMnonIso_e",L1EMnonIso_e,"L1EMnonIso_e[nL1EMnonIso]/F");
  Analysis->Branch("L1EMnonIso_et",L1EMnonIso_et,"L1EMnonIso_et[nL1EMnonIso]/F");
  Analysis->Branch("L1EMnonIso_eta",L1EMnonIso_eta,"L1EMnonIso_eta[nL1EMnonIso]/F");
  Analysis->Branch("L1EMnonIso_phi",L1EMnonIso_phi,"L1EMnonIso_phi[nL1EMnonIso]/F");
  
  // L1 Muon objects
  Analysis->Branch("nL1Muon",&nL1Muon,"nL1Muon/I");
  Analysis->Branch("L1Muon_e",&L1Muon_e,"L1Muon_e[nL1Muon]/F");
  Analysis->Branch("L1Muon_pt",&L1Muon_pt,"L1Muon_pt[nL1Muon]/F");
  Analysis->Branch("L1Muon_eta",&L1Muon_eta,"L1Muon_eta[nL1Muon]/F");
  Analysis->Branch("L1Muon_phi",&L1Muon_phi,"L1Muon_phi[nL1Muon]/F");
  

  
  //MET
  Analysis->Branch("caloMETet",&caloMETet,"caloMETet/F");
  Analysis->Branch("caloMETsumEt",&caloMETsumEt,"caloMETsumEt/F");
  Analysis->Branch("caloMETphi",&caloMETphi,"caloMETphi/F");
  Analysis->Branch("muCorrMETsumEt",&muCorrMETsumEt,"muCorrMETsumEt/F");
  Analysis->Branch("muCorrMETet",&muCorrMETet,"muCorrMETet/F");
  Analysis->Branch("muCorrMETphi",&muCorrMETphi,"muCorrMETphi/F");
  Analysis->Branch("tcMETsumEt",&tcMETsumEt,"tcMETsumEt/F");
  Analysis->Branch("tcMETet",&tcMETet,"tcMETet/F");
  Analysis->Branch("tcMETphi",&tcMETphi,"tcMETphi/F");
  Analysis->Branch("pfMETsumEt",&pfMETsumEt,"pfMETsumEt/F");
  Analysis->Branch("pfMETet",&pfMETet,"pfMETet/F");
  Analysis->Branch("pfMETphi",&pfMETphi,"pfMETphi/F");
  
  

  cout<<" RecoAnalyzer tree branches defined.. " <<endl; 
  
  

}



////test the particle-to-parton Matching
///retrun index in MCpart array
int RecoAnalyzer::partonMatchingAlgo(float eta, float phi, float dr){
  
  
  int tempParticle = -1;
  int tempPartonHighestPt = -1;
  int tempNearest = -1;
  float maxPt = 0;
  float minDr = 1000;
  ///bool foundPriority = false;
  
  float coneSizeToAssociate = dr; 
  
  if(  nMCpart <1) return -99; 
  if( int(partonList.size()) <1) return -99; 
  
    
  for (int m = 0; m<int(partonList.size()); m++){
    
    int j = partonList[m];
    
    int pid = abs(pidMCpart[j]);
    
    /// if( pid >= 6 && pid != 21 ) continue; 
    float dist = GetDeltaR(eta,etaMCpart[j],phi,phiMCpart[j]);
    
    if( pidFDauMCpart[j] == 91 || pidFDauMCpart[j] == 92){
      
      
      if( dist < coneSizeToAssociate ){
	if( dist < minDr ) {
	  minDr = dist;
	  tempNearest = j;
	}
	
	if( tempParticle == -1 && pid == 4  ) tempParticle = j;
	if(  pid == 5    ) tempParticle = j;
	if( ptMCpart[j] > maxPt ) {
	  maxPt = ptMCpart[j];
	  tempPartonHighestPt = j;
	}
      }
    }
  }
  
  if(tempParticle ==-1) tempParticle = tempPartonHighestPt; 
  
  if( tempParticle <0) return 0; 
  else return tempParticle; 
    
  
}


int RecoAnalyzer::findIndexofTRK(float pt, float eta, float phi){
  
  float drmin = 0.00001; 
  int flag = -1; 
  for(int j=0; j< nTrack; j++){
    if( fabs( pt - trackpt[j]) > 0.00001) continue; 
    float dr = GetDeltaR(eta,tracketa[j],phi,trackphi[j]);
    if( dr <drmin){
      drmin = drmin; 
      flag = j; 
    }
  }
  return flag;
  
}


int RecoAnalyzer::findIndexofSC(float en,float eta,float phi){
  int flag =-1; 
  int nfound =0;
  float err = 0.0001; 
  for( int j=0; j<nSC; j++){
    if( fabs(eSC[j]-en)<err && fabs(etaSC[j]-eta)<err && fabs(phiSC[j]-phi)<err){
      flag = j; 
      nfound++; 
    }
  }
  return flag; 
}


int RecoAnalyzer::getMotherIndex(int j){
  if(j<0) {
    cout<<"no mother. input -1!!!!!!!!"<<endl;
    return -1; 
  }
  
  int indmom = barcodemomMCpart[j];
  
  while( indmom >=0 && pidMCpart[indmom]==pidMCpart[j] ){ //if it is the same pid 
    //int st = statusMCpart[indmom];
    indmom = barcodemomMCpart[indmom];
  }
  
  return indmom; 
  
  
}


void RecoAnalyzer::MatchToGenElectron(float eta, float phi, float res[]){
  
  for(int j=0; j<20; j++) res[j] = 1; 
  float drMin = 0.3; 
  int indMin = -1; 

  float drMin1 = 0.3; 
  int indMin1 = -1; 

  
  float drMin2 = 0.3; 
  int indMin2 = -1;

  
  if( nMCpart <1) return; 

  for(int j=1; j<nMCpart; j++){

    float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
    
    if( ptMCpart[j] > 3 && statusMCpart[j] ==1 && abs(pidMCpart[j])==11){ ///gen real electron
      if(dr <drMin){
	drMin = dr; 
	indMin = j; 
      }
    }
    
    if( ptMCpart[j] > 3 && statusMCpart[j] == 1 &&  abs(chargeMCpart[j]) == 1  && abs(pidMCpart[j]) != 11){ //gen charged, not electron  
      if(dr <drMin1){
	drMin1 = dr; 
	indMin1 = j; 
      }
    }
    if(statusMCpart[j] == -99 && barcodemomMCpart[j]>0 && ptMCpart[j] >3 && abs(pidMCpart[j]) == 11){ ///from photon conversion
      if( statusMCpart[barcodemomMCpart[j]] ==1 && pidMCpart[barcodemomMCpart[j]] == 22){
	if( dr < drMin2){
	  drMin2 = dr; 
	  indMin2 = j; 
	}
      }
      
    } 
  }
  

  if( indMin>0){ ///matched to a gen-electron
    res[0] = drMin; 
    res[1] = ptMCpart[indMin];
    res[2] = pidMCpart[indMin];
    int indm = getMotherIndex(indMin); 
    if(indm >=0){
      res[3] = pidMCpart[indm];
    }else{
      res[3] = 0; 
    }
  }
  
  if( indMin1>0){///matched to a gen-charged
    res[4] = drMin1; 
    res[5] = ptMCpart[indMin1];
    res[6] = pidMCpart[indMin1];
    int indm = getMotherIndex(indMin1); 
    if(indm >=0){
      res[7] = pidMCpart[indm];
    }else{
      res[7] = 0; 
    }
  }
  if( indMin2>0){ ///matched to electron from photon conv.
    res[8] = drMin2; 
    res[9] = ptMCpart[indMin2];
    res[10] = pidMCpart[indMin2];
    res[11] = vtxXMCpart[indMin2];
    res[12] = vtxYMCpart[indMin2];
    res[13] = vtxZMCpart[indMin2];
    int indmm = getMotherIndex(barcodemomMCpart[indMin2]);
    if(indmm >=0){
      res[14] = pidMCpart[indmm]; ///photon's mother
    }else{
      res[14] = 0; 
    }
  }
    
}

void RecoAnalyzer::MatchToGenPhoton(float eta, float phi, float res[]){
  
  
  float drMin = 0.2; 
  int indMin = -1; 
  for(int j=0; j<20; j++){
    res[j] = 1; 
  }
  
  if( nMCpart <1) return; 
  
  for(int j=0; j<nMCpart; j++){
    if(statusMCpart[j] !=1) continue;
    if( ptMCpart[j]< 3 ) continue;  ///don't matching to very small pt photon's
    if( abs(pidMCpart[j]) != 22 ) continue; 
    float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);

    if(dr <drMin){
      drMin = dr; 
      indMin = j; 
    }
    
  }
  
  
  if(indMin>=0){
    res[0] = drMin; 
    ///res[1] = ptMCpart[indMin]/sin(2*atan(exp(-etaMCpart[indMin])));
    //res[2] = etaMCpart[indMin];
    //res[3] = phiMCpart[indMin];
    res[1] = ptMCpart[indMin];
    ///mother's index and pid
    ///for qcd events,45->67.
    int indmom = getMotherIndex(indMin);
    //res[5] = indmom;
    if( indmom>=0){
      res[2] = pidMCpart[indmom];
    }else{
      res[2] = 0;
    }
    ///NOTE pythia 8 have photon's mother id with 22 , so change from 1,6, 21 13, 11, to >=1 && <30 
    ///if it's pi0/eta check decay angle of two photon
    ///if it's not from ISR or FSR.
    
    
  }
  
  float drMin1 = 0.2; 
  int indMin1 = -1; 
  
  for(int j=0; j<nMCpart; j++){
    
    if(statusMCpart[j] !=1) continue;
    if( abs(pidMCpart[j]) !=11) continue;  ///matched to electron
    
    if( ptMCpart[j]<10) continue; 
    
    float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
    
    if(dr <drMin1){
      drMin1 = dr; 
      indMin1 = j; 
    }
  }
  
  
  if( indMin1>0){ //matched to a real electron
    res[7] = drMin1; 
    res[8] = ptMCpart[indMin1];
    int indm = getMotherIndex(indMin1); 
    if(indm >=0){
      res[9] = pidMCpart[indm];
    }else{
      res[9] = 0; 
    }
  }
  
  
  if(indMin>=0){
    if( convPhtMCpart[indMin][0] > 0){
      //conv vtx
      res[10] = 1; 
      res[11] = vtxXMCpart[convPhtMCpart[indMin][1]]; 
      res[12] = vtxYMCpart[convPhtMCpart[indMin][1]]; 
      res[13] = vtxZMCpart[convPhtMCpart[indMin][1]]; 
    }
    else {
      res[10] = 0; 
      res[11] = 0; 
      res[12] = 0; 
      res[13] = 0; 
    }
  }else{
    res[10] = 0; 
    res[11] = 0; 
    res[12] = 0; 
    res[13] = 0; 
  }
  
  
  
}



int RecoAnalyzer::indexofParticle(float px, float pz, int status){
  
  float err = 0.00001; 
  
  
  int indMin = -1; 
  float errMin = 0.002; 
  for( int j=0; j<nMCpart; j++){
    if(fabs(px-pxMCpart[j])<err && fabs(pz-pzMCpart[j])<err
       && status == statusMCpart[j]){
      float err1 = fabs(px-pxMCpart[j]) + fabs(pz-pzMCpart[j]);
      if( err1 < errMin){
	errMin = err1; 
	indMin = j; 
      }
    }
  }
  
  return indMin; 
  
  //   if(int(nn.size()) !=1){ ///diquarks or gluon has two copy sometimesss!!! 
  //     cout<<"wrong! indexofParticle: n: "<<int(nn.size())<<endl;
  //     for(int j=0; j<int(nn.size()); j++){
  //       cout<<"  dup j: "<<j<<" ind "<<nn[j]<<endl;
  //     }
  //   }
  //return nn[0]; 
  
}


double RecoAnalyzer::DeltaPhi(double v1, double v2)
{ // Computes the correctly normalized phi difference
  // v1, v2 = phi of object 1 and 2
  double diff = v1 - v2;
  
  while (diff >acos(-1)) diff -= 2*acos(-1);
  while (diff <= -acos(-1)) diff += 2*acos(-1);
  
  return diff; 
}


double RecoAnalyzer::GetDeltaR(double eta1, double eta2, double phi1, double phi2){ 
  // Computes the DeltaR of two objects from their eta and phi values
  
  return sqrt( (eta1-eta2)*(eta1-eta2) + DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
  
}


void RecoAnalyzer::MatchToGenMuon(float eta, float phi, float res[]){
  
  float drMin = 0.3; 
  int indMin = -1; 
  
  float drMin1 = 0.3; 
  int indMin1 = -1; 
  
  float drMin2 = 0.3; 
  int indMin2 = -1; 
  
  for(int j=0; j<20; j++){
    res[j] = 1; 
  }

  if( nMCpart <1) return; 
    

  for(int j=6; j<nMCpart; j++){
    
    if( ptMCpart[j] < 3 ) continue; 
        
    ///gen real muon 
    if( abs(pidMCpart[j]) == 13 && statusMCpart[j] ==1 ){
      float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
      if(dr <drMin){
	drMin = dr; 
	indMin = j; 
      }
    }
    
    ///gen real charged non-muon , non-eletron particles
    if(statusMCpart[j] ==1 && abs(chargeMCpart[j]) ==1 && abs(pidMCpart[j]) != 13 && abs(pidMCpart[j]) !=11){
      float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
     if(dr <drMin1){
	drMin1 = dr; 
	indMin1 = j; 
     }
    }
    ////sim real moun from gen
    if(statusMCpart[j]== -99 && abs(pidMCpart[j]) == 13 && barcodemomMCpart[j] >0 ){
      float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
      if( dr< drMin2){
	drMin2 = dr; 
	indMin2 = j; 
      }
    }
    
  }
  
  
  if(indMin>0){ ///if found a real muon in gen-level
    res[0] = drMin; 
    res[1] = ptMCpart[indMin];
    res[2] = pidMCpart[indMin];

    int indm = getMotherIndex(indMin); 
    if( indm >=0) {
      res[3] = pidMCpart[indm];
    }else{
      res[3] = 0; 
    }
    
  }
  
  if(indMin1>0){ ///if found a charged particle in gen-level
    res[4] = drMin1; 
    res[5] = ptMCpart[indMin1];
    res[6] = pidMCpart[indMin1];
    int indm = getMotherIndex(indMin1);
    if( indm>=0){
      res[7] = pidMCpart[indm];
    }else{
      res[7] = 0; 
    }
    
  }
  
  if( indMin2>0){ ///if found a muon in sim-level
    res[8] = drMin2; 
    res[9] = ptMCpart[indMin2];
    res[10] = pidMCpart[indMin2];
    res[11] = vtxXMCpart[indMin2];
    res[12] = vtxYMCpart[indMin2];
    res[13] = vtxZMCpart[indMin2];
    int indm = getMotherIndex(indMin2);
    if( indm>=0){
      ///res[14] = pidMCpart[barcodemomMCpart[indMin2]]; /// muon's mother, mostly should be pi0/kion decay with staus ==1 
      res[14] =  pidMCpart[indm];
    }else{
      res[14] = 0; 
    }
  }
  
  
}


double RecoAnalyzer::etaTransformation(  float EtaParticle , float Zvertex)  {
  
  //---Definitions
  const float pi = 3.1415927;

  //---Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479; 
   
  //---ETA correction

  float Theta = 0.0  ; 
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+pi ;
  double ETA = - log(tan(0.5*Theta));
         
  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - Zvertex ;
      float RR = Zlen/sinh(EtaParticle); 
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+pi ;
      ETA = - log(tan(0.5*Theta));		      
    } 
  //---Return the result
  return ETA;
  //---end
}



//////////from HiggsAnalysis/HiggsTo2photons/h2gglobe/GeneralFunctions_cc.h
////retrun the index to assced conv..[]
int  RecoAnalyzer::matchPhotonToConversion( int lpho) {

  int result=-99;
  double conv_eta=-999.;
  double conv_phi=-999.;
  
  //float sc_eta  = ((TVector3 *) pho_calopos->At(lpho))->Eta();
  //  float sc_phi  = ((TVector3 *) pho_calopos->At(lpho))->Phi();
  float sc_eta = photonsceta[lpho];
  float sc_phi = photonscphi[lpho];
  

  //TLorentzVector * p4 = (TLorentzVector *) pho_p4->At(lpho);
  // float et = p4->Energy() / cosh(sc_eta);
  //float et = photonenergy[lpho] / cosh(sc_eta);
  
  //  cout << " photon index " << lpho << " eta " <<sc_eta<< " phi " << sc_phi << " et " << et << endl; 
  
  float detaMin=999.;
  float dphiMin=999.;   
  ///float dRMin = 999.;

  float mconv_pt=-999999;
  int iMatch=-1;     
  
  if(debug_ > 10 )  cout << "   LoopAll::matchPhotonToConversion conv_n " << nConv << endl; 
  for(int iconv=0; iconv<nConv; iconv++) {
    
    //vector<TVector3> refittedPairMomentum

    ///refittedPairMomentum.push_back(*((TVector3*) conv_refitted_momentum->At(iconv)));
    
    TVector3 refittedPairMomentum(convrefittedPairMomentumx[iconv],convrefittedPairMomentumy[iconv],convrefittedPairMomentumz[iconv]);
    float conv_pt =  refittedPairMomentum.Pt();
    if (conv_pt < 1 ) continue;    
    if ( convconversionVertexisValid[iconv] == 0 ) continue;
    
    conv_phi  = refittedPairMomentum.Phi();
    float eta  = refittedPairMomentum.Eta();
    conv_eta = etaTransformation(eta, photonconversionzOfPrimaryVertexFromTracks[iconv] ); //transfer eta to eta at ECAL 

    //float dPhi =conv_phi - sc_phi;       
    //double delta_phi = phiNorm (dPhi);
    double delta_phi = DeltaPhi(conv_phi, sc_phi);
    double delta_eta = conv_eta - sc_eta;
    
    //cout << " delta_eta " << delta_eta << " delta_phi " << delta_phi << endl;
    //delta_phi=pow(delta_phi,2);
    //delta_eta=pow(delta_eta,2);
    //float dR = sqrt( delta_phi+delta_eta); 
    
    if ( fabs(delta_eta) < detaMin && fabs(delta_phi) < dphiMin ) {  ///why not using a deltaR matching  ??
      //    if ( dR < dRMin ) {
      detaMin=  fabs(delta_eta);
      dphiMin=  fabs(delta_phi);
      //dRMin=dR;
      iMatch=iconv;
      mconv_pt = conv_pt;
    }
    
  }
  //  cout << " minimized conversion index " << iMatch << " eta " <<conv_eta<< " phi " << conv_phi <<endl; 
  if ( detaMin < 0.1 && dphiMin < 0.1 ) {
    //  if ( dRMin< 0.1 ) {
    if(debug_)    cout << " matched conversion index " << iMatch << " eta " <<conv_eta<< " phi " << conv_phi << " pt " << mconv_pt << endl; 	
    result = iMatch;
  } else {
    result = -1;
  }
  
  return result;
    
}



///from  HggVertexFromConversions::vtxZ(const PhotonInfo & pho)
double RecoAnalyzer::HggVertexFromConversionsvtxZ(int j)
{
  
  bool useAllConvs= true; 
  
  
  double phtConX = 0; 
  double phtConY = 0; 
  double phtConZ = 0; 
  if( useAllConvs && photonmatchToallConv[j] >=0){
    phtConX = convconversionVertexx[photonmatchToallConv[j]];
    phtConY = convconversionVertexy[photonmatchToallConv[j]];
    phtConZ = convconversionVertexz[photonmatchToallConv[j]];
  }else{
    phtConX = photonconversionVertexx[j];
    phtConY = photonconversionVertexy[j];
    phtConZ = photonconversionVertexz[j];
  }
  
  // get the z from conversions
  double deltaX1 =  photoncaloPositionx[j] - phtConX;
  double deltaY1 =  photoncaloPositiony[j] - phtConY;
  double deltaZ1 =  photoncaloPositionz[j] - phtConZ;
  double R1 = sqrt(deltaX1*deltaX1+deltaY1*deltaY1);
  double tantheta = R1/deltaZ1;
  
  double deltaX2 = phtConX-beamSpotX; 
  double deltaY2 = phtConY-beamSpotY; 
  double R2 = sqrt(deltaX2*deltaX2+deltaY2*deltaY2);
  double deltaZ2 = R2/tantheta;
  double higgsZ =  photoncaloPositionz[j]-deltaZ1-deltaZ2;
  return higgsZ;
  
}



///from  HggVertexFromConversions::vtxdZ(const PhotonInfo & pho)
double RecoAnalyzer::HggVertexFromConversionsvtxdZ(int j){
  
  ///// attribute the error depending on the tracker region
  double dz=-99999;
  
  bool useAllConvs= true; 
  
  double phtConX = 0; 
  double phtConY = 0; 
  double phtConZ = 0; 
  if( useAllConvs && photonmatchToallConv[j] >=0){
    phtConX = convconversionVertexx[photonmatchToallConv[j]];
    phtConY = convconversionVertexy[photonmatchToallConv[j]];
    phtConZ = convconversionVertexz[photonmatchToallConv[j]];
  }else{
    phtConX = photonconversionVertexx[j];
    phtConY = photonconversionVertexy[j];
    phtConZ = photonconversionVertexz[j];
  }
  
  float sigmaPix_ =0.06;
  float sigmaTib_ =0.67;
  float sigmaTob_ =2.04;
  float sigmaFwd1_ =0.18;
  float sigmaFwd2_ =0.61;
  float sigmaFwd3_ =0.99;
  
  
  ///if ( pho.iDet() ==1 ) { // barrel
  if ( fabs(photonsceta[j])<1.48 ) { // barrel

    double perp = sqrt( phtConX*phtConX + phtConY*phtConY);
    if ( perp <=15 ) {
      dz=sigmaPix_;
    } else if ( perp > 15 && perp <=60 ) {
      dz=sigmaTib_;
    } else {
      dz=sigmaTob_;
    }

  } else { // endcap
    
    if ( fabs(phtConZ ) <=50 ) {
      dz=sigmaFwd1_;
    } else if (phtConZ > 50 &&  phtConZ <= 100 ) {
      dz=sigmaFwd2_;
    } else {
      dz=sigmaFwd3_;
    }
  }
  
  return dz;
  
}


TLorentzVector RecoAnalyzer::photonp4wrtvertex(int indpht, int indvtx){
  
  TVector3 vPos(vertexx[indvtx],vertexy[indvtx],vertexz[indvtx]);
  TVector3 caloPosition(photoncaloPositionx[indpht],photoncaloPositiony[indpht],photoncaloPositionz[indpht]);
  TVector3 direction = caloPosition - vPos;
  TVector3 p = direction.Unit() * photonenergy[indpht];
  TLorentzVector p4(p.x(),p.y(),p.z(),photonenergy[indpht]);
  return p4;
  
}



////from HiggsAnalysis/HiggsTo2photons/h2gglobe/GeneralFunctions_cc.h
//// LoopAll::vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, 
////					  PhotonInfo & pho1, PhotonInfo & pho2, std::vector<std::string> & vtxVarNames)

//Input index of the two photon
void RecoAnalyzer::higgsVertexAnalysis(int indpht1, int indpht2){
  

  ///first build quantities /// inside HggVertexAnalyzer::analyze()
  //analyze(const VertexInfoAdapter & e, const PhotonInfo & p1, const PhotonInfo & p2)
  
  std::vector<unsigned short> vtxTracksBuf;
  std::vector<int> vtxTracksSizeBuf;
  //if( ! e.hasVtxTracks() ) {
  if( false ) {
 //    vtxTracksBuf.resize(nvtx*e.ntracks());
//     vtxTracksSizeBuf.resize(nVertex,0);
//     for(int it=0; it<e.ntracks(); ++it) {
//       int vid = e.tkVtxId(it);
      
//       if( it < 5) cout<<" tkVtxId " << it <<" "<<  vid <<endl; 
      
//       int & ntks = vtxTracksSizeBuf[vid];
//       vtxTracksBuf[ vid*e.ntracks() + ntks ] = it;
//       ++ntks;
//     }
  }
  
  vector<int> preselection_; 
  preselection_.clear();

  bool usehighPurityOnly = false; 
    
  float maxD0Signif = 999; 
  float maxDzSignif = 999; 

  bool rescaleTkPtByError = false; 
  float trackCountThr = 0; 
  
  bool removeTracksInCone = true; 
  float coneSize = 0.05; 
  
  std::vector<float> ptbal_;
  std::vector<float> thrust_;
  std::vector<float> sumpt_;
  std::vector<float> sumpt2_;
  std::vector<float> sumawy_;
  std::vector<float> sumtwd_;
  std::vector<float> sumtrv_;
  std::vector<float> sumweight_;
  std::vector<float> ptmax_;
  std::vector<float> nchthr_;
  std::vector<float> nch_;
  std::vector<std::vector<float> > tksPt_;
  std::vector<TMatrixDSym> sphers_;
  std::vector<float> sumpr_;
  std::vector<float> spher_;
  std::vector<float> tspher_;
  std::vector<float> aplan_;
  std::vector<float> threejetC_;
  std::vector<float> fourjetD_;
  
  std::vector<TVector3> vtxP_;
  std::vector<TVector2> vtxPt_;
  std::vector<float> ptvtx_;
  std::vector<float> diphopt_; 
  std::vector<float> diPhotonPz_;
  std::vector<float> acosA_;
  std::vector<float> ptasym_;
  
  std::vector<float> ptmax3_;
  std::vector<float> ptratio_;
  std::vector<float> pzasym_;
  std::vector<float> awytwdasym_;
  std::vector<TVector2> diPhotonPt_;
  
  
  ptbal_.clear(); ptbal_.resize(nVertex,0.);
  thrust_.clear(); thrust_.resize(nVertex,0.);
  sumpt_.clear(); sumpt_.resize(nVertex,0.);
  sumpt2_.clear(); sumpt2_.resize(nVertex,0.);
  sumawy_.clear(); sumawy_.resize(nVertex,0.);
  sumtwd_.clear(); sumtwd_.resize(nVertex,0.);
  sumtrv_.clear(); sumtrv_.resize(nVertex,0.);
  sumweight_.clear(); sumweight_.resize(nVertex,0.);
  ptmax_.clear(); ptmax_.resize(nVertex,0.);
  nchthr_.clear(); nchthr_.resize(nVertex,0.);
  nch_.clear(); nch_.resize(nVertex,0.);
  vtxP_.clear(); vtxP_.resize(nVertex,0.);
  tksPt_.clear(); tksPt_.resize(nVertex, vector<float>(1));
  sphers_.clear(); sphers_.resize(nVertex,TMatrixDSym(3));
  sumpr_.clear(); sumpr_.resize(nVertex,0.);
  spher_.clear(); spher_.resize(nVertex,0.);
  tspher_.clear(); tspher_.resize(nVertex,0.);
  aplan_.clear(); aplan_.resize(nVertex,0.);
  threejetC_.clear(); threejetC_.resize(nVertex,0.);
  fourjetD_.clear(); fourjetD_.resize(nVertex,0.);
	
  diphopt_.clear(); diphopt_.resize(nVertex);
  diPhotonPt_.clear(); diPhotonPt_.resize(nVertex);
  vtxPt_.clear(); vtxPt_.resize(nVertex);
  ptvtx_.clear(); ptvtx_.resize(nVertex);
  diPhotonPz_.clear(); diPhotonPz_.resize(nVertex);
	
  acosA_.clear(); acosA_.resize(nVertex);
  ptasym_.clear(); ptasym_.resize(nVertex);
	
  ptmax3_.clear(); ptmax3_.resize(nVertex);
  thrust_.clear(); thrust_.resize(nVertex);
	
  ptratio_.clear(); ptratio_.resize(nVertex);
  pzasym_.clear(); pzasym_.resize(nVertex);
	
  awytwdasym_.clear(); awytwdasym_.resize(nVertex);
  


  static const float spherPwr_ = 1.5; 
  
  

  std::vector<TLorentzVector> diPhoton_;
  diPhoton_.clear(); diPhoton_.resize(nVertex);
  for(int i=0; i<nVertex; ++i) {
    diPhoton_[i] = photonp4wrtvertex(indpht1,i)
      + photonp4wrtvertex(indpht2,i); 
  }
  
  
  // filling loop over vertexes
  for(int n=0; n<nVertex; n++) {
    
    //foreach vertex loopover all the it's associated  track
    for(int t =0 ; t < vertextrackSize[n]; t++){
      
      float tkWeight = vertex_trkWeight->at(n)[t];
      int indtrk = vertex_trkind->at(n)[t];

      if( indtrk <0 || indtrk >= nTrack){ // 
	continue; 
      }
      
      if( usehighPurityOnly ){
	int qua = trackqualityFlagTracks[indtrk]; 
	
	int highPurityFlag = 3; 
	bool isHP = false; 
	if( ( qua & 1 << highPurityFlag) > 0){
	  isHP = true; 
	}
	if( !isHP) continue; 
      }
      
      if( trackd0[indtrk]/ trackd0Error[indtrk] > maxD0Signif ) continue; 
      if( trackdz[indtrk]/ trackdzError[indtrk] > maxDzSignif ) continue; 
      
      float tkPt = trackpt[indtrk];
      float modpt = tkPt > trackptError[indtrk] ? tkPt - trackptError[indtrk]  : 0.;
      
      TVector3 tkPVec(trackpx[indtrk],trackpy[indtrk],trackpz[indtrk]);
      
      TVector2 tkPtVec(trackpx[indtrk],trackpy[indtrk]);
            
      
      if(rescaleTkPtByError){
	if( modpt ==0) continue; 
	float ptcorr = modpt / tkPt; 
	tkPtVec *= ptcorr;
	tkPt = modpt;
      }
      sumpt2_[n] += tkPtVec.Mod2();
      sumpt_[n] += tkPt;
      if(tkPt > trackCountThr) nchthr_[n] += 1;
      nch_[n] += 1;
      
      
      if ( removeTracksInCone){
	float dr1 = tkPVec.DeltaR( photonp4wrtvertex(indpht1,n).Vect());
	float dr2 = tkPVec.DeltaR( photonp4wrtvertex(indpht2,n).Vect());
	if( dr1 < coneSize || dr2 < coneSize){
	  continue; 
	}
      }
      ptbal_[n] -= tkPtVec * diPhoton_[n].Vect().XYvector().Unit();
      float cosTk = tkPVec.Unit() * diPhoton_[n].Vect().Unit();
      float val = tkPtVec.Mod();
      if ( cosTk < -0.5 )	{
	sumawy_[n] += val;
      } else if ( cosTk > 0.5 ){
	sumtwd_[n] += val;
      } else {
	sumtrv_[n] += val;
      }
      sumweight_[n] += tkWeight;
      vtxP_[n] += tkPVec;
      tksPt_[n].push_back(tkPt);
      
      Float_t p[3] = {0.,0.,0.};
      tkPVec.GetXYZ(p);
      for(int j=3; j--;){
	for(int k=j+1; k--;){
	  (sphers_[n])[j][k] += pow(tkPVec.Mag(),spherPwr_-2.) * p[j]*p[k];
	}
      }
      sumpr_[n] += pow(tkPVec.Mag(),spherPwr_);
      
      
      //if( ( usehighPurityOnly && !e.tkIsHighPurity(tid)  )
      ///|| fabs(e.tkd0(tid,n)/e.tkd0Err(tid,n)) > params_.maxD0Signif 
	    //|| fabs(e.tkdz(tid,n)/e.tkdzErr(tid,n)) > params_.maxDzSignif ) {
	//// std::cerr << "skipping [line" << __LINE__ << "] ";
	//// std::cerr << (params_.highPurityOnly && !e.tkIsHighPurity(tid)  )
	//// 	  << " " 
	//// 	  << (fabs(e.tkd0(tid,n)/e.tkd0Err(tid,n)) > params_.maxD0Signif )
	//// 	  << " " 
	//// 	  << (fabs(e.tkdz(tid,n)/e.tkdzErr(tid,n)) > params_.maxDzSignif) << std::endl;
      
      
    }
    
    sphers_[n] *= 1./sumpr_[n];
    
    TVectorD eigVals(3);
    eigVals = TMatrixDSymEigen(sphers_[n]).GetEigenValues();
    
    spher_[n] = 1.5 * (eigVals[1]+eigVals[2]);
    tspher_[n] = 2. * eigVals[1] / (eigVals[0]+eigVals[1]);
    aplan_[n] = 1.5 * eigVals[2];
    
    threejetC_[n] = 3. * (eigVals[0]*eigVals[1] + eigVals[0]*eigVals[2] + eigVals[1]*eigVals[2]);
    fourjetD_[n] = 27. * eigVals[0]*eigVals[1]*eigVals[2];
    
    
    diPhotonPt_[n] = diPhoton_[n].Vect().XYvector();
    diphopt_[n]    = diPhotonPt_[n].Mod();
    vtxPt_[n]      = vtxP_[n].XYvector();
    ptvtx_[n]      = vtxPt_[n].Mod();
    diPhotonPz_[n]   = diPhoton_[n].Vect().Pz();
    
    sort(tksPt_[n].begin(), tksPt_[n].end(), greater<float>());
    
    acosA_[n] =  	acos(vtxPt_[n].Unit() * diPhotonPt_[n].Unit());
    ptasym_[n] = 	(vtxPt_[n].Mod() - diPhotonPt_[n].Mod())/(vtxPt_[n].Mod() + diPhotonPt_[n].Mod());
    
    ptmax_ [n] = 	tksPt_[n][0];
    ptmax3_[n] = 	accumulate(tksPt_[n].begin(),tksPt_[n].begin() + min(tksPt_[n].size(),(size_t)3), 0.0) ;
    thrust_[n] = 	ptbal_[n]/sumpt_[n];
    
    ptratio_[n] = 	vtxPt_[n].Mod()/diPhotonPt_[n].Mod();
    pzasym_[n] = 	fabs( (vtxP_[n].Pz() - diPhotonPz_[n])/(vtxP_[n].Pz() + diPhotonPz_[n]) );
    
    awytwdasym_[n] = (sumawy_[n]-sumtwd_[n])/(sumawy_[n]+sumtwd_[n]);
    
    
  } //end of loop over each rec vertex
    
  
  bool useAllConvs= true; 
  
  // preselect vertices : all vertices
  std::vector<int> preselAll;
  for(int i=0; i< nVertex; i++) {
    preselAll.push_back(i); 
  }
  
  float zconv = 0; 
  float dzconv = 0;
  std::vector<int> preselConv;
  
  //// VertexAnalysis/src/PhotonInfo.cc:  if (  nTracks_ == 2  &&  convVtxValid_ &&   convVtxChi2Prob_ > 0.0005 )  isAConversion=true;
  //bool  PhotonInfo::isAConversion() {
  // bool isAConversion=false;
  // if (  nTracks_ == 2  &&  convVtxValid_ &&   convVtxChi2Prob_ > 0.0005 )  isAConversion=true;
  // return isAConversion;

  //}

  int iConv1 = useAllConvs ? matchPhotonToConversion(indpht1) : -1;
  int iConv2 = useAllConvs ? matchPhotonToConversion(indpht2) : -1;
  
  //saved iConv1 and iConv2 as well
  photonmatchToallConv[indpht1] = iConv1; 
  photonmatchToallConv[indpht2] = iConv2; 
  
  
  bool ispht1AConversion; 
  if( iConv1 >=0 ){
    ispht1AConversion = convnTracks[iConv1] == 2 && convconversionVertexisValid[iConv1] && convconversionVertexChiSquaredProbability[iConv1] > 0.0005; 
  }else{
    ispht1AConversion = photonconversionnTracks[indpht1] ==2 && photonconversionVertexisValid[indpht1] ==1  && photonconversionChiSquaredProbability[indpht1] > 0.0005; 
  }
  bool ispht2AConversion; 
  
  if( iConv2 >=0 ){
    ispht2AConversion = convnTracks[iConv2] == 2 && convconversionVertexisValid[iConv2] && convconversionVertexChiSquaredProbability[iConv2] > 0.0005; 
  } else{
    ispht2AConversion= photonconversionnTracks[indpht2] ==2 && photonconversionVertexisValid[indpht2] ==1  && photonconversionChiSquaredProbability[indpht2] > 0.0005; 
  }
  if( ispht1AConversion || ispht2AConversion){
    
    if( ispht1AConversion && ! ispht2AConversion){
      //zconv  = vtxAnaFromConv.vtxZ(pho1);
      //dzconv = vtxAnaFromConv.vtxdZ(pho1);
      
      zconv = HggVertexFromConversionsvtxZ(indpht1);
      dzconv = HggVertexFromConversionsvtxdZ(indpht1);
    }
    if( ! ispht1AConversion && ispht2AConversion){
      zconv = HggVertexFromConversionsvtxZ(indpht2);
      dzconv = HggVertexFromConversionsvtxdZ(indpht2);
    }
    if(  ispht1AConversion && ispht2AConversion){
      float z1 = HggVertexFromConversionsvtxZ(indpht1);
      float dz1 = HggVertexFromConversionsvtxdZ(indpht1);
      float z2 = HggVertexFromConversionsvtxZ(indpht2);
      float dz2 = HggVertexFromConversionsvtxdZ(indpht2);
      
      zconv  = sqrt ( 1./(1./dz1/dz1 + 1./dz2/dz2 )*(z1/dz1/dz1 + z2/dz2/dz2) ) ;  // weighted average
      dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
      
    }
    // preselect vertices : only vertices in a window zconv +/- dzconv
    for(int i=0; i < nVertex ; i++) {
      //TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(i);
      //if ( fabs(zconv - vtxpos->Z() ) < dzconv ) 
      if( fabs( zconv - vertexz[i] ) < dzconv){
	preselConv.push_back(i); 
      }
    }
    
  }
  
  ///now preslected vertex 
  vector<int> preselection; 
  if ( preselConv.size()==0 ) 
    preselection = preselAll; 
  else 
    preselection = preselConv; 
  
  
  //std::vector<int> rankprod = vtxAna.rankprod(vtxVarNames);
  //return rankprod;
  
  
  
  
}









float RecoAnalyzer::recHitE( const  DetId id,  const EcalRecHitCollection &recHits )
{
  if ( id == DetId(0) ) {
    return 0;
  } else {
    EcalRecHitCollection::const_iterator it = recHits.find( id );
    if ( it != recHits.end() ) return (*it).energy();
  }
  return 0;
}


float RecoAnalyzer::recHitE( const DetId id, const EcalRecHitCollection & recHits,int di, int dj )
{
  // in the barrel:   di = dEta   dj = dPhi
  // in the endcap:   di = dX     dj = dY
  
  DetId nid;
  if( id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy( id, di, dj );
  else if( id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy( id, di, dj );
  
  return ( nid == DetId(0) ? 0 : recHitE( nid, recHits ) );
}




// ------------ method called to produce the data  ------------
bool RecoAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool accept = false; 
  
  if (iEvent.isRealData()) {
    // for real data, access the "physics enabled" bit in the L1 GT data
    edm::Handle<L1GlobalTriggerReadoutRecord> h_gtDigis;
    if (not iEvent.getByLabel(m_l1GtRecordInputTag, h_gtDigis)) {
      ////edm::LogWarning(h_gtDigis.whyFailed()->category()) << h_gtDigis.whyFailed()->what();
      // not enough informations to make a decision - reject the event
      return false;
    } else {
      L1GtFdlWord fdlWord = h_gtDigis->gtFdlWord();
      if (fdlWord.physicsDeclared() == 1) 
        accept = true;
    }
  } else {
    // for MC, assume the "physics enabled" bit to be always set
    accept = true;
  }
  
  return accept;
  
  
}










float RecoAnalyzer::recHitApproxEt(  const DetId id,  const EcalRecHitCollection &recHits )
{
  // for the time being works only for the barrel
  if ( id.subdetId() == EcalBarrel ) {
    return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
  }
  return 0;
}

double RecoAnalyzer::mye2overe9( const DetId id, const EcalRecHitCollection & recHits)
{
  ///////////start calculating e2/e9                                                                                                                                   
  ////http://cmslxr.fnal.gov/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc#240                                                                    
  // compute e2overe9                                                                                                                                                  
  //                                                                                                                                                                   
  //   | | | |                                                                                                                                                         
  //   +-+-+-+                                                                                                                                                         
  //   | |1|2|                                                                                                                                                         
  //   +-+-+-+                                                                                                                                                         
  //   | | | |                                                                                                                                                         
  //                                                                                                                                                                   
  //   1 - input hit,  2 - highest energy hit in a 3x3 around 1  
  //                                                                                                                                                                     
  //   rechit 1 must have E_t > recHitEtThreshold                                                                                                                        
  //   rechit 2 must have E_t > recHitEtThreshold2                                                                                                                       
  //                                                                                                                                                                     
  //   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2                                                                        
  //                                                                                                                                                                     
  //   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2                                                              
  //   *provided* that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0                                                             

  
  float recHitEtThreshold = 10.0;
  float recHitEtThreshold2 = 1.0;
  bool avoidIeta85=false;
  bool KillSecondHit=true;
  

  if ( id.subdetId() == EcalBarrel ) {

    EBDetId ebId( id );
    // avoid recHits at |eta|=85 where one side of the neighbours is missing                                                                                             
    if ( abs(ebId.ieta())==85 && avoidIeta85) return 0;

    // select recHits with Et above recHitEtThreshold                                                                                                                    
    float e1 = recHitE( id, recHits );
    float ete1=recHitApproxEt( id, recHits );


    // check that rechit E_t is above threshold                                                                                                                          
    if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) return 0;

    if (ete1 < recHitEtThreshold && !KillSecondHit ) return 0;

    float e2=-1;
    float ete2=0;
    float s9 = 0;

    // coordinates of 2nd hit relative to central hit 
    int e2eta=0;
    int e2phi=0;

    // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1                                                                                                                         

    for ( int deta = -1; deta <= +1; ++deta ) {
      for ( int dphi = -1; dphi <= +1; ++dphi ) {

	// compute 3x3 energy                                                                                                                                            
	float etmp=recHitE( id, recHits, deta, dphi );
	s9 += etmp;

	EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
	float eapproxet=recHitApproxEt( idtmp, recHits );
	
	// remember 2nd highest energy deposit (above threshold) in 3x3 array                                                                                            
	if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {
	  e2=etmp;
	  ete2=eapproxet;
	  e2eta=deta;
	  e2phi=dphi;

	}

      }
    }

    if ( e1 == 0 )  return 0;

    // return 0 if 2nd hit is below threshold                                                                                                                            
    if ( e2 == -1 ) return 0;

    // compute e2/e9 centered around 1st hit                                                                                                                             
    float e2nd=e1+e2;
    float e2e9=0;
    if (s9!=0) e2e9=e2nd/s9;

    // if central hit has higher energy than 2nd hit                                                                                                                     
    //  return e2/e9 if 1st hit is above E_t threshold                                                                                                                   

    if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;

    // if second hit has higher energy than 1st hit                                                                                                                      
    if ( e2 > e1 ) {


      // return 0 if user does not want to flag 2nd hit, or                                                                                                              
      // hits are below E_t thresholds - note here we                                                                                                                    
      // now assume the 2nd hit to be the leading hit.                                                                                                                   

      if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {
	return 0;

      }
            
      else {

	// LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2                                                                                                                     
	float s92nd=0;

	float e2nd_prime=0;
	int e2prime_eta=0;
	int e2prime_phi=0;

	EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);

	for ( int deta = -1; deta <= +1; ++deta ) {
	  for ( int dphi = -1; dphi <= +1; ++dphi ) {
	    // compute 3x3 energy                                                                                                                                        
	    float etmp=recHitE( secondid, recHits, deta, dphi );
	    s92nd += etmp;

	    if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
	      e2nd_prime=etmp;
	      e2prime_eta=deta;
	      e2prime_phi=dphi;
	    }

	  }
	}

	// if highest energy hit around E2 is not the same as the input hit, return 0;                                                                                   
	if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi))
	  {
	    return 0;
	  }


	// compute E2/E9 around second hit                                                                                                                               
	float e2e9_2=0;
	if (s92nd!=0) e2e9_2=e2nd/s92nd;

	//   return the value of E2/E9 calculated around 2nd hit                                                                                                         
	return e2e9_2;


      }

    }

  } else if ( id.subdetId() == EcalEndcap ) {
    // only used for EB at the moment                                                                                                                                    
    return 0;
  }
  return 0;
}//double Analyser::mye2overe9( const DetId id, const EcalRecHitCollection & recHits)                                                                                  



// ------------ method called once each job just after ending the event loop  ------------
void 
RecoAnalyzer::endJob() 
{

  rootFile_->cd();
  Analysis->Write();
  
  rootFile_->Write() ;
  rootFile_->Close() ;

}

// ------------ method called when starting to processes a run  ------------
void 
RecoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
RecoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
RecoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
RecoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoAnalyzer);
