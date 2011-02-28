#ifndef _EcalRecHitsRegional_H
#define _EcalRecHitsRegional_H

// -*- C++ -*-
//
// Package:    EcalRecHitsRegional
// Class:      EcalRecHitsRegional
// 
// \class EcalRecHitsRegional EcalRecHitsRegional.cc 

/* Description: 
This is the filter to reduce full-unpacked EcalRechits (What one gets normally from RECO) to the "regionally unpacked" EcalRecHits. 
There is not actually unpacking going on in this filter, it just selects the EcalRecHits which should have been unpacked with some 
defined threshold in the RawToRecHit process and save those EcalRecHits. 

The main puropose is that one can run from RECO, instead of from RAW.  

Right now Feb.28.11, only works for L1-egamma seeds regional unpacking. 

*/

//
// Original Author:  Yong Yang 
//
//


// system include files
#include <memory>

// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "HLTrigger/HLTcore/interface/HLTFilter.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "Geometry/EcalMapping/interface/EcalElectronicsMapping.h"



#include "DataFormats/L1Trigger/interface/L1EmParticle.h"

#include "DataFormats/L1Trigger/interface/L1EmParticleFwd.h"


#include "Geometry/EcalMapping/interface/EcalMappingRcd.h"
#include "CondFormats/L1TObjects/interface/L1CaloGeometry.h"
#include "CondFormats/DataRecord/interface/L1CaloGeometryRecord.h"


//
// class decleration
//

class EcalRecHitsRegional : public HLTFilter {
   public:
      explicit EcalRecHitsRegional(const edm::ParameterSet&);
      ~EcalRecHitsRegional();


      virtual bool filter(edm::Event &, const edm::EventSetup&);
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
