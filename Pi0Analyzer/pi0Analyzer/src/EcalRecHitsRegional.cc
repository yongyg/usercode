#include "Pi0Analyzer/pi0Analyzer/interface/EcalRecHitsRegional.h"
#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "DataFormats/HLTReco/interface/TriggerFilterObjectWithRefs.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "FWCore/Framework/interface/ESHandle.h"


using namespace edm;
using namespace std;

EcalRecHitsRegional::EcalRecHitsRegional(const edm::ParameterSet& iConfig)
{
  barrelHits_ = iConfig.getParameter< edm::InputTag > ("barrelHitCollection");
  endcapHits_ = iConfig.getParameter< edm::InputTag > ("endcapHitCollection");

  selectedBarrelHits_ =
    iConfig.getParameter< std::string > ("selectedBarrelHitCollection");
  selectedEndcapHits_ =
    iConfig.getParameter< std::string > ("selectedEndcapHitCollection");


  RegionalMatch_ = iConfig.getParameter<bool>("RegionalMatch");
  ptMinEMObj_ = iConfig.getParameter<double>("ptMinEMObj");
  EMregionEtaMargin_ = iConfig.getParameter<double>("EMregionEtaMargin"); ///0.25, //              
  EMregionPhiMargin_ = iConfig.getParameter<double>("EMregionPhiMargin"); ///0.4                    
  l1IsolatedTag_ = iConfig.getParameter< edm::InputTag > ("l1IsolatedTag");
  l1NonIsolatedTag_ = iConfig.getParameter< edm::InputTag > ("l1NonIsolatedTag");
  debug_ = iConfig.getParameter<int> ("debugLevel");

  TheMapping = new EcalElectronicsMapping();
  first_ = true;

  //register your products                                                                             
  produces< EBRecHitCollection >(selectedBarrelHits_);
  produces< EERecHitCollection >(selectedEndcapHits_);
  

}


EcalRecHitsRegional::~EcalRecHitsRegional()
{
}

// ------------ method called to produce the data  ------------
void
EcalRecHitsRegional::produce(edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  if (first_) {
    edm::ESHandle< EcalElectronicsMapping > ecalmapping;
    iSetup.get< EcalMappingRcd >().get(ecalmapping);
    const EcalElectronicsMapping* TheMapping_ = ecalmapping.product();
    *TheMapping = *TheMapping_;
    first_ = false;
    
  }
  
  
  vector<int>::iterator it; 
  
  
  FEDListUsed.clear();
  if( RegionalMatch_ ){

    // Get the L1 CaloGeometry
    edm::ESHandle<L1CaloGeometry> l1CaloGeom ;
    iSetup.get<L1CaloGeometryRecord>().get(l1CaloGeom) ;
    edm::Handle< l1extra::L1EmParticleCollection > l1EGIso;
    edm::Handle< l1extra::L1EmParticleCollection > l1EGNonIso ;
    ///Get L1 EG ojbects
    try{
      iEvent.getByLabel(l1IsolatedTag_, l1EGIso ) ;
      for( l1extra::L1EmParticleCollection::const_iterator emItr = l1EGIso->begin();
           emItr != l1EGIso->end() ;++emItr ){
      
        if(  emItr -> pt() < ptMinEMObj_ ) continue; 
      

        if( debug_>=2) cout<<"L1EGISO: "<<  emItr -> pt() <<" "<< emItr ->eta()<<" "<<emItr -> phi()<<endl;

        int etaIndex = emItr->gctEmCand()->etaIndex() ;
        int phiIndex = emItr->gctEmCand()->phiIndex() ;
        double etaLow  = l1CaloGeom->etaBinLowEdge( etaIndex ) ;
        double etaHigh = l1CaloGeom->etaBinHighEdge( etaIndex ) ;
        double phiLow  = l1CaloGeom->emJetPhiBinLowEdge( phiIndex ) ;
        double phiHigh = l1CaloGeom->emJetPhiBinHighEdge( phiIndex ) ;
    
        std::vector<int> feds = ListOfFEDS(etaLow, etaHigh, phiLow, phiHigh, EMregionEtaMargin_, EMregionPhiMargin_);
        for (int n=0; n < (int)feds.size(); n++) {
          int fed = feds[n];
          it = find(FEDListUsed.begin(),FEDListUsed.end(),fed);
          if( it == FEDListUsed.end()){
            FEDListUsed.push_back(fed);
          }
        }
      }
      

    }catch(std::exception& ex ){
      cout<<"l1EGIso not working.."<<endl;
    }
      
    try{
      iEvent.getByLabel(l1NonIsolatedTag_, l1EGNonIso ) ;
       
      for( l1extra::L1EmParticleCollection::const_iterator emItr = l1EGNonIso->begin();
           emItr != l1EGNonIso->end() ;++emItr ){
    
        if( emItr -> pt()< ptMinEMObj_ ) continue; 
      
        if( debug_>=2) cout<<"L1EGnonIso: "<<  emItr -> pt() <<" "<< emItr ->eta()<<" "<<emItr -> phi()<<endl;
      
        int etaIndex = emItr->gctEmCand()->etaIndex() ;
        int phiIndex = emItr->gctEmCand()->phiIndex() ;
        double etaLow  = l1CaloGeom->etaBinLowEdge( etaIndex ) ;
        double etaHigh = l1CaloGeom->etaBinHighEdge( etaIndex ) ;
        double phiLow  = l1CaloGeom->emJetPhiBinLowEdge( phiIndex ) ;
        double phiHigh = l1CaloGeom->emJetPhiBinHighEdge( phiIndex ) ;
    
        std::vector<int> feds = ListOfFEDS(etaLow, etaHigh, phiLow, phiHigh, EMregionEtaMargin_, EMregionPhiMargin_);
        for (int n=0; n < (int)feds.size(); n++) {
          int fed = feds[n];
          it = find(FEDListUsed.begin(),FEDListUsed.end(),fed);
          if( it == FEDListUsed.end()){
            FEDListUsed.push_back(fed);
          }
        }
      }
      
    }catch(std::exception& ex ){
      cout<<"l1EGnonIso not working.."<<endl;
    }
    
    
  
  }   //// end of getting FED List
  ///separate into barrel and endcap to speed up when checking
  FEDListUsedBarrel.clear();
  FEDListUsedEndcap.clear();
  for(  int j=0; j< int(FEDListUsed.size());j++){
    int fed = FEDListUsed[j];
    if( fed >= 10 && fed <= 45){
      FEDListUsedBarrel.push_back(fed);
    }else FEDListUsedEndcap.push_back(fed);
  }
  


//   edm::ESHandle<EcalChannelStatus> csHandle;
//   if (! useRecoFlag_) iSetup.get<EcalChannelStatusRcd>().get(csHandle);
//   const EcalChannelStatus& channelStatus = *csHandle; 
  


  Handle<EBRecHitCollection> barrelRecHitsHandle;
  Handle<EERecHitCollection> endcapRecHitsHandle;

  
  iEvent.getByLabel(barrelHits_,barrelRecHitsHandle);
  iEvent.getByLabel(endcapHits_,endcapRecHitsHandle);
 
  //Create empty output collections
  std::auto_ptr< EBRecHitCollection > selectedEBRecHitCollection( new EBRecHitCollection );
  std::auto_ptr< EERecHitCollection > selectedEERecHitCollection( new EERecHitCollection );


  //Select interesting EcalRecHits (barrel)
  EBRecHitCollection::const_iterator itb;
  for (itb=barrelRecHitsHandle->begin(); itb!=barrelRecHitsHandle->end(); itb++) {
    

      EBDetId det = itb->id();
      if (RegionalMatch_){
        int fed = TheMapping->DCCid(det);
        it = find(FEDListUsedBarrel.begin(),FEDListUsedBarrel.end(),fed);
        if(it == FEDListUsedBarrel.end()) continue; 
      }
   
      selectedEBRecHitCollection->push_back(*itb);
      
  }
  
  //Select interesting EcalRecHits (endcaps)
  EERecHitCollection::const_iterator ite;
  for (ite=endcapRecHitsHandle->begin(); ite!=endcapRecHitsHandle->end(); ite++) {

    EEDetId det = ite->id();
    
      if (RegionalMatch_){
        EcalElectronicsId elid = TheMapping->getElectronicsId(det);
        int fed = elid.dccId();
        it = find(FEDListUsedEndcap.begin(),FEDListUsedEndcap.end(),fed);
        if(it == FEDListUsedEndcap.end()) continue; 
      }
      
      selectedEERecHitCollection->push_back(*ite);
      
  }
  
  if( debug_>=1) cout<<" selectedEB/EErechits: "<< selectedEBRecHitCollection->size()<<" "<<selectedEERecHitCollection->size()<<endl; 
    
  
  //Put selected information in the event
  iEvent.put( selectedEBRecHitCollection, selectedBarrelHits_);
  iEvent.put( selectedEERecHitCollection, selectedEndcapHits_);
  
  
}



//////FED list this is obsolete 
std::vector<int> EcalRecHitsRegional::ListOfFEDS(double etaLow, double etaHigh, double phiLow, 
                                                 double phiHigh, double etamargin, double phimargin)
{
  
        std::vector<int> FEDs;

        if (phimargin > Geom::pi()) phimargin =  Geom::pi() ;


        if (debug_>=2) std::cout << " etaLow etaHigh phiLow phiHigh " << etaLow << " " << 
                        etaHigh << " " << phiLow << " " << phiHigh << std::endl;

        etaLow -= etamargin;
        etaHigh += etamargin;
        double phiMinus = phiLow - phimargin;
        double phiPlus = phiHigh + phimargin;

        bool all = false;
        double dd = fabs(phiPlus-phiMinus);
        if (debug_>=2) std::cout << " dd = " << dd << std::endl;
        if (dd > 2.*Geom::pi() ) all = true;

        while (phiPlus > Geom::pi()) { phiPlus -= 2.*Geom::pi() ; }
        while (phiMinus < 0) { phiMinus += 2.*Geom::pi() ; }
        if ( phiMinus > Geom::pi()) phiMinus -= 2.*Geom::pi() ;

        double dphi = phiPlus - phiMinus;
        if (dphi < 0) dphi += 2.*Geom::pi() ;
        if (debug_>=2) std::cout << "dphi = " << dphi << std::endl;
        if (dphi > Geom::pi()) {
                int fed_low1 = TheMapping -> GetFED(etaLow,phiMinus*180./Geom::pi());
                int fed_low2 = TheMapping -> GetFED(etaLow,phiPlus*180./Geom::pi());
                if (debug_>=2) std::cout << "fed_low1 fed_low2 " << fed_low1 << " " << fed_low2 << std::endl;
                if (fed_low1 == fed_low2) all = true;
                int fed_hi1 = TheMapping -> GetFED(etaHigh,phiMinus*180./Geom::pi());
                int fed_hi2 = TheMapping -> GetFED(etaHigh,phiPlus*180./Geom::pi());
                if (debug_>=2) std::cout << "fed_hi1 fed_hi2 " << fed_hi1 << " " << fed_hi2 << std::endl;
                if (fed_hi1 == fed_hi2) all = true;
        }

        if (all) {
                if (debug_>=2) std::cout << " unpack everything in phi ! " << std::endl;
                phiMinus = -20 * Geom::pi() / 180.;  // -20 deg
                phiPlus = -40 * Geom::pi() / 180.;  // -20 deg
        }

        if (debug_>=2) std::cout << " with margins : " << etaLow << " " << etaHigh << " " << 
                        phiMinus << " " << phiPlus << std::endl;


        const EcalEtaPhiRegion ecalregion(etaLow,etaHigh,phiMinus,phiPlus);

        FEDs = TheMapping -> GetListofFEDs(ecalregion);

/*
        if (debug_) {
        int nn = (int)FEDs.size();
           for (int ii=0; ii < nn; ii++) {
           std::cout << "unpack fed " << FEDs[ii] << std::endl;
           }
           }
*/

        return FEDs;

}

//define this as a plug-in
///DEFINE_FWK_MODULE(EcalRecHitsRegional); 
