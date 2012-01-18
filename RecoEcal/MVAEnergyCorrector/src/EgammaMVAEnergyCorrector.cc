// -*- C++ -*-
//
// Package:    EgammaMVAEnergyCorrector
// Class:      EgammaMVAEnergyCorrector
// 
/**\class EgammaMVAEnergyCorrector EgammaMVAEnergyCorrector.cc RecoEcal/MVAEnergyCorrector/src/EgammaMVAEnergyCorrector.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yong Yang
//         Created:  Wed Nov 23 09:31:00 CST 2011
// $Id: EgammaMVAEnergyCorrector.cc,v 1.4 2012/01/18 14:11:30 yangyong Exp $
//
//

#include "RecoEcal/MVAEnergyCorrector/interface/EgammaMVAEnergyCorrector.h"

void EgammaMVAEnergyCorrector::Initialize(const edm::EventSetup &iSetup, vector<std::string> regweights){
  
  
  fIsInitialized = kTRUE;
  PhotonFix::initialiseGeometry(iSetup);
  
  if (fVals) delete [] fVals;
  if (reader_eb) delete reader_eb;
  if (reader_ebErr) delete reader_ebErr;
  if (reader_ee) delete reader_ee;
  if (reader_eeErr) delete reader_eeErr;
  
  methodName = TString("BDTG") + TString(" method");
  
  reader_eb = new TMVA::Reader( "!Color:!Silent" );
  reader_ebErr = new TMVA::Reader( "!Color:!Silent" );
  reader_ee = new TMVA::Reader( "!Color:!Silent" );
  reader_eeErr = new TMVA::Reader( "!Color:!Silent" );
  
  fVals = new Float_t[26];
  vector<TString> ebvarname; 
  ebvarname.push_back("escraw");
  ebvarname.push_back("etasc");
  ebvarname.push_back("phisc");
  ebvarname.push_back("e5x5/escraw");    
  ebvarname.push_back("r9");    
  ebvarname.push_back("e2x5left/escraw");
  ebvarname.push_back("e2x5right/escraw");
  ebvarname.push_back("e2x5top/escraw");
  ebvarname.push_back("e2x5bottom/escraw");
  ebvarname.push_back("emax/escraw");
  ebvarname.push_back("eleft/escraw");
  ebvarname.push_back("eright/escraw");
  ebvarname.push_back("etop/escraw");
  ebvarname.push_back("ebottom/escraw");
  ebvarname.push_back("sigietaieta");
  ebvarname.push_back("scetawidth");
  ebvarname.push_back("scphiwidth");
  ebvarname.push_back("etaSGap");
  ebvarname.push_back("phiSGap");
  ebvarname.push_back("etaMGap");
  ebvarname.push_back("phiMGap");
  ebvarname.push_back("etaCGap");
  ebvarname.push_back("phiCGap");
  ebvarname.push_back("nVertex");
  ebvarname.push_back("rho");

  vector<TString> eevarname; 
  for(unsigned int j =0; j< ebvarname.size(); j++){
    if( j==3){
      eevarname.push_back("eps/escraw");
    }
    if( j >= 17 && j <= 22) continue; 
    eevarname.push_back(ebvarname[j]);
  }
  for(unsigned int j =0; j< ebvarname.size(); j++){
    reader_eb->AddVariable(ebvarname[j],&fVals[j]);
    reader_ebErr->AddVariable(ebvarname[j],&fVals[j]);
  }
  reader_ebErr->AddVariable("escregcorr/escraw",&fVals[25]);

  for(unsigned int j =0; j< eevarname.size(); j++){
    reader_ee->AddVariable(eevarname[j],&fVals[j]);
    reader_eeErr->AddVariable(eevarname[j],&fVals[j]);
  }
  reader_eeErr->AddVariable("escregcorr/(escraw+eps)",&fVals[20]);
  reader_eb->BookMVA(methodName,regweights[0].c_str());
  reader_ebErr->BookMVA(methodName,regweights[1].c_str());
  reader_ee->BookMVA(methodName,regweights[2].c_str());
  reader_eeErr->BookMVA(methodName,regweights[3].c_str());
  
}



//--------------------------------------------------------------------------------------------------
EgammaMVAEnergyCorrector::EgammaMVAEnergyCorrector() :
  reader_eb(0),
  reader_ebErr(0),
  reader_ee(0),
  reader_eeErr(0),
  fIsInitialized(kFALSE),
  fVals(0)
{
  // Constructor.
}



EgammaMVAEnergyCorrector::~EgammaMVAEnergyCorrector()
{
  if (fVals) delete [] fVals;
  if (reader_eb) delete reader_eb;
  if (reader_ebErr) delete reader_ebErr;
  if (reader_ee) delete reader_ee;
  if (reader_eeErr) delete reader_eeErr;
}


std::pair<float,float> EgammaMVAEnergyCorrector::CorrectedEnergyWithError(const Photon &p, EcalClusterLazyTools &clustertools, int nvtx, float rho) {
  
  const SuperCluster &s = *p.superCluster();
  const BasicCluster &b = *s.seed();
  PhotonFix phfix(s.eta(),s.phi()); 
  Bool_t isbarrel =  b.hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
  
  if (isbarrel) {
    fVals[0]  = s.rawEnergy();
    fVals[1]  = s.eta();
    fVals[2]  = s.phi();
    fVals[3]  = p.e5x5() /s.rawEnergy();
    fVals[4]  = p.r9();
    fVals[5]  = clustertools.e2x5Left(b);
    fVals[6]  = clustertools.e2x5Right(b);
    fVals[7]  = clustertools.e2x5Top(b);
    fVals[8]  = clustertools.e2x5Bottom(b);
    fVals[9]  = p.maxEnergyXtal();
    fVals[10]  = clustertools.eLeft(b);
    fVals[11]  = clustertools.eRight(b);
    fVals[12]  = clustertools.eTop(b);
    fVals[13]  = clustertools.eBottom(b);
    for(int j=5; j<=13; j++){
      fVals[j] /= s.rawEnergy();
    }
    
    fVals[14]  = p.sigmaIetaIeta();
    fVals[15]  = s.etaWidth();
    fVals[16]  = s.phiWidth();
    fVals[17]  = phfix.etaS();
    fVals[18]  = phfix.phiS();
    fVals[19]  = phfix.etaM();
    fVals[20]  = phfix.phiM();    
    fVals[21]  = phfix.etaC();
    fVals[22]  = phfix.phiC();
    fVals[23]  = nvtx; 
    fVals[24]  = rho; 
    
  }else{
    
    fVals[0]  = s.rawEnergy();
    fVals[1]  = s.eta();
    fVals[2]  = s.phi();
    fVals[3]  = s.preshowerEnergy()/s.rawEnergy();
    fVals[4]  = p.e5x5() / s.rawEnergy(); 
    fVals[5]  = p.r9();
    fVals[6]  = clustertools.e2x5Left(b);
    fVals[7]  = clustertools.e2x5Right(b);
    fVals[8]  = clustertools.e2x5Top(b);
    fVals[9]  = clustertools.e2x5Bottom(b);
    fVals[10]  = p.maxEnergyXtal();
    fVals[11]  = clustertools.eLeft(b);
    fVals[12]  = clustertools.eRight(b);
    fVals[13]  = clustertools.eTop(b);
    fVals[14]  = clustertools.eBottom(b);
    for(int j=6; j<=14; j++){
      fVals[j] /= s.rawEnergy();
    }
    
    fVals[15]  = p.sigmaIetaIeta();
    fVals[16]  = s.etaWidth();
    fVals[17]  = s.phiWidth();
    fVals[18]  = nvtx; 
    fVals[19]  = rho; 
    
  }
  
  //const Float_t varscale = 1.253;
  const Float_t varscale = 1; //this number needs to be re-studied.right now put 1 
  
  float corr; 
  float corrErr; 
  float eraw ; 
  float ecorr; 
  float ecorrErr; 
  if (isbarrel) {
    corr  = (reader_eb->EvaluateRegression(methodName))[0];
    fVals[25] = corr; 
    corrErr  = (reader_ebErr->EvaluateRegression(methodName))[0];
    eraw = s.rawEnergy() ; 
  }else{
    corr  = (reader_ee->EvaluateRegression(methodName))[0];
    fVals[20] = corr;
    corrErr  = (reader_eeErr->EvaluateRegression(methodName))[0];
    eraw = s.rawEnergy()+s.preshowerEnergy();
  }
  ecorr = eraw * corr; 
  ecorrErr = ecorr * corrErr * varscale; ///error is trained to |escregcorr-etrue|/escregcorr
  
  return std::pair<float,float>(ecorr,ecorrErr);
  
}


std::pair<float,float> EgammaMVAEnergyCorrector::CorrectedEnergyWithError(const GsfElectron &e, EcalClusterLazyTools &clustertools, int nvtx, float rho) {
  
  const SuperCluster &s = *e.superCluster();
  ///for tracker-driven only, no correction
  if (!e.ecalDrivenSeed()) return std::pair<double,double>(s.energy(),0);
  
  const BasicCluster &b = *s.seed();
  PhotonFix phfix(s.eta(),s.phi()); 
  Bool_t isbarrel =  b.hitsAndFractions().at(0).first.subdetId()==EcalBarrel;
  
  if (isbarrel) {
    fVals[0]  = s.rawEnergy();
    fVals[1]  = s.eta();
    fVals[2]  = s.phi();
    fVals[3]  = e.e5x5();
    fVals[4]  = clustertools.e3x3(b); 
    fVals[5]  = clustertools.e2x5Left(b);
    fVals[6]  = clustertools.e2x5Right(b);
    fVals[7]  = clustertools.e2x5Top(b);
    fVals[8]  = clustertools.e2x5Bottom(b);
    fVals[9]  = clustertools.eMax(b);
    fVals[10]  = clustertools.eLeft(b);
    fVals[11]  = clustertools.eRight(b);
    fVals[12]  = clustertools.eTop(b);
    fVals[13]  = clustertools.eBottom(b);
    for(int j=3; j<=13; j++){
      fVals[j] /= s.rawEnergy();
    }

    fVals[14]  = e.sigmaIetaIeta();
    fVals[15]  = s.etaWidth();
    fVals[16]  = s.phiWidth();
    fVals[17]  = phfix.etaS();
    fVals[18]  = phfix.phiS();
    fVals[19]  = phfix.etaM();
    fVals[20]  = phfix.phiM();    
    fVals[21]  = phfix.etaC();
    fVals[22]  = phfix.phiC();
    fVals[23]  = nvtx; 
    fVals[24]  = rho; 
    
  }else{
    
    fVals[0]  = s.rawEnergy();
    fVals[1]  = s.eta();
    fVals[2]  = s.phi();
    fVals[3]  = s.preshowerEnergy(); 
    fVals[4]  = e.e5x5();
    fVals[5]  = clustertools.e3x3(b);
    fVals[6]  = clustertools.e2x5Left(b);
    fVals[7]  = clustertools.e2x5Right(b);
    fVals[8]  = clustertools.e2x5Top(b);
    fVals[9]  = clustertools.e2x5Bottom(b);
    fVals[10]  = clustertools.eMax(b);
    fVals[11]  = clustertools.eLeft(b);
    fVals[12]  = clustertools.eRight(b);
    fVals[13]  = clustertools.eTop(b);
    fVals[14]  = clustertools.eBottom(b);
    for(int j=3; j<=14; j++){
      fVals[j] /= s.rawEnergy();
    }
    fVals[15]  = e.sigmaIetaIeta();
    fVals[16]  = s.etaWidth();
    fVals[17]  = s.phiWidth();
    fVals[18]  = nvtx; 
    fVals[19]  = rho; 
    
  }
  
  //const Float_t varscale = 1.253;
  const Float_t varscale = 1; //this number needs to be re-studied.right now put 1 
  float corr; 
  float corrErr; 
  float eraw ; 
  float ecorr; 
  float ecorrErr; 
  if (isbarrel) {
    corr  = (reader_eb->EvaluateRegression(methodName))[0];
    fVals[25] = corr;
    corrErr  = (reader_ebErr->EvaluateRegression(methodName))[0];
    eraw = s.rawEnergy() ; 
  }else{
    corr  = (reader_ee->EvaluateRegression(methodName))[0];
    fVals[20] = corr;
    corrErr  = (reader_eeErr->EvaluateRegression(methodName))[0];
    eraw = s.rawEnergy()+s.preshowerEnergy();
  }
  ecorr = eraw * corr; 
  ecorrErr = ecorr * corrErr * varscale; ///error is trained to |escregcorr-etrue|/escregcorr
  return std::pair<float,float>(ecorr,ecorrErr);
  
}



