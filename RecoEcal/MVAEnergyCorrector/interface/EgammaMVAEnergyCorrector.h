// $Id $
//
// EgammaMVAEnergyCorrector
//
// Helper Class for applying regression-based energy corrections
//
// Authors: Y.Yang
//--------------------------------------------------------------------------------------------------

#ifndef EgammaMVAEnergyCorrector_H
#define EgammaMVAEnergyCorrector_H

#include "PhotonFix.h"

// system include files
#include <memory>
// user include files
#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "DataFormats/EgammaCandidates/interface/PhotonFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "RecoEcal/EgammaCoreTools/interface/EcalClusterLazyTools.h"
#include "DataFormats/EgammaCandidates/interface/Photon.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EcalDetId/interface/EcalSubdetector.h"

#include "TROOT.h"
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
using namespace TMVA;
using namespace std;
using namespace edm;
using namespace reco;


//class EcalClusterLazyTools;

//
// class declaration
//
class EgammaMVAEnergyCorrector{
 public:
  EgammaMVAEnergyCorrector();
  ~EgammaMVAEnergyCorrector();
  
  void Initialize(const edm::EventSetup &iSetup, vector<std::string> regweights);
  Bool_t IsInitialized() const { return fIsInitialized; }
  
  std::pair<float,float> CorrectedEnergyWithError(const reco::Photon &p,EcalClusterLazyTools &clustertools, int nvtx, float rho);
  std::pair<float,float> CorrectedEnergyWithError(const reco::GsfElectron &e, EcalClusterLazyTools &clustertools,int nvtx, float rho);
  
 protected:
  

  
  TMVA::Reader *reader_eb;
  TMVA::Reader *reader_ebErr;
  TMVA::Reader *reader_ee;
  TMVA::Reader *reader_eeErr;
  Bool_t fIsInitialized;
  Float_t *fVals;
  
  TString methodName; 

};

#endif
