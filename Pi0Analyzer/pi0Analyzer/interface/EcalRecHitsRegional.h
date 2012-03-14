#ifndef _EcalRecHitsRegional_H
#define _EcalRecHitsRegional_H


#include <memory>
#include <time.h>

#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDProducer.h"
#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/EventSetup.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"

#include "DataFormats/EcalRecHit/interface/EcalRecHitCollections.h"
#include "RecoEcal/EgammaCoreTools/interface/PositionCalc.h"
#include "Geometry/CaloTopology/interface/CaloSubdetectorTopology.h"

//
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"
#include "DataFormats/Math/interface/Point3D.h"
#include "DataFormats/CaloRecHit/interface/CaloID.h"

//Ecal status
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalChannelStatus.h"



#include "DataFormats/L1Trigger/interface/L1EmParticle.h"
#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"
#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"


#include "TVector3.h"
#include <vector>




class EcalRecHitsRegional : public edm::EDProducer 
{
  public:

      EcalRecHitsRegional(const edm::ParameterSet& ps);

      ~EcalRecHitsRegional();
      
     
      
      virtual void produce(edm::Event&, const edm::EventSetup&);
      
 private:
      
      std::vector<int> ListOfFEDS(double etaLow, double etaHigh, double phiLow,
                                  double phiHigh, double etamargin, double phimargin);
      

      // ----------member data ---------------------------
      
      
      edm::InputTag barrelHits_;
      edm::InputTag endcapHits_;
      std::string selectedBarrelHits_;
      std::string selectedEndcapHits_;

      EcalElectronicsMapping* TheMapping;
 
      edm::InputTag l1IsolatedTag_;
      edm::InputTag l1NonIsolatedTag_;
      std::vector<int> FEDListUsed; ///by regional objects.  ( em, jet, etc)
  
      std::vector<int> FEDListUsedBarrel; 
      std::vector<int> FEDListUsedEndcap; 
  
      bool RegionalMatch_;
  
      double ptMinEMObj_ ; 
  
      double EMregionEtaMargin_;
      double EMregionPhiMargin_;
  
      bool first_; 
  
      int debug_;
   
      
};


#endif
