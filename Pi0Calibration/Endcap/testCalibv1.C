#include <iomanip>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TString.h"
#include "TFile.h"
#include "TMath.h"
#include "TBranch.h"
#include "TChain.h"
#include <iostream>
#include <algorithm>
#include <map>
#include "TLorentzVector.h"
#include "TVector3.h"
#include "TApplication.h"
#include "TSystem.h"
#include "TF1.h"
#include "TGraphErrors.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TRandom.h"
#include "TRandom3.h"
#include "TMultiLayerPerceptron.h"
#include "TMLPAnalyzer.h"
#include "TMatrixD.h"
#include "TMatrix.h"
#include "TMatrixDSymEigen.h"
#include "TMatrixDSym.h"
#include "TMatrixTSym.h"
#include <vector>
#include "TArrow.h"
#include "TLegend.h"
#include "TROOT.h"
#include <stdlib.h>
#include "TStyle.h"
#include "TLatex.h"
using namespace std;

///for L3 method, 
double WTSUM[2][101][101];  /// wt = pow( en_crystal/ en3x3,2);
double CORSUM[2][101][101];  // wt * (mean/peak)*(mean*peak)

float C0[2][101][101];


//Common files                                                                                         
#include "../Common/variables.h"
#include "../Common/utils.cc"
#include "../Common/PositionCalc.cc"
#include "../Common/getGoodLS.cc"

#include "foldEndcap.cc"

#include "setBranchAddress.cc"
#include "datachain.cc"
#include "calibUtils.cc"
#include "getpreCalibConst.cc"
#include "correction.cc"


///int nCounted[2][101][101][100]; 

double binwidth = 0.005; 
double xHighLimit = 0.3; 
double xLowLimit = 0; 
int nbinMax = int(xHighLimit/binwidth + 0.1);



///for endcap 
void testCalibv1(int test_dataflag,int test_pizEta, int test_calibStep, int test_calibIteration, int test_evtRange){
  workingDirectory = "MyJobWorkingDirectory";
  
  doBarrel =2; 

  nMaxRingIC = 20; 
  if(pizEta==2){
    nMaxRingIC = 30;
  }
  
  

  readInterCalibEndcap_GR09_V8();
  
  
  dataflag = test_dataflag; 
  pizEta  = test_pizEta; 
  evtRange = test_evtRange; 
  stepc = test_calibStep;
  iter = test_calibIteration;

  if(pizEta!=1 && pizEta!=2){
    cout<<"pizEta 1 or 2  " << pizEta <<endl;
    return;
  }
  
  if( pizEta==2){ //this is for the output text 
    binwidth = 0.01; 
    xHighLimit = 0.9; 
    xLowLimit = 0.2; 
    nbinMax = int( (xHighLimit-xLowLimit) /binwidth + 0.1);
  }
  
  cout<<"binwidth: "<< binwidth <<" "<< xLowLimit <<" "<< xHighLimit <<" "<<nbinMax <<endl; 
  


  cout<<"initilization counts.." <<endl;
  for(int iz = 0; iz<2; iz++){
    for(int j=0; j<101; j++){
      for(int k=0; k<101; k++){
	for(int n=0; n< nbinMax; n++){
	  nCountedEE[iz][j][k][n] = 0; 
	}
      }
    }
  }
  
  for(int j=0; j<2; j++){
    for(int k=0; k< kEndcEtaRings; k++){
      corrfactorEtaRings[j][k] = 1; 
    } 
    for(int k=0; k< 101; k++){
      for(int n=0; n< 101; n++){
	corrfactoriZiXiY[j][k][n] = 1; 
      }
    }
  }


  ///x/y/z of each crystal to re-calculate X/Y/Z/ of cluster                                         
  get_xyzEBrechits(); 
  setEtaRingBoundaryEndcap();  //for eta-ring 
  
  if(stepc==1 && iter>=2 ){ ///this is to dervie eta correction.
    getcorrFactorEtaOfEachStep();
  }

  if(stepc==2){ //this is to get the final correction of eta-rings
    getcorrFactorEta(); 
    if( iter>=2){
      getcorrFactorIzIxIyOfEachStep();
    }
  }
  
  for(int iz = 0; iz<2; iz++){
    for(int j=0; j<101; j++){
      for(int k=0; k<101; k++){
	C0[iz][j][k] = 1; 
	if( validRecHitEndCap[iz][j][k] <1) continue;
	if( stepc ==2 && iter>=2){
	  C0[iz][j][k] =  corrfactoriZiXiY[iz][j][k];
	}
      }
    }
  }
    
  if(pizEta==1){
    fChain = new TChain("pizSele");
  }else{
    fChain = new TChain("etaSele");
  }

  addDataChain();
  
  setBranchAddress();
  

  totalEntries = fChain->GetEntries(); 
  
  TString filename; 
  
  /// in the endcap, eta depdendent cuts. 
  float ptminCut_eta[4]; 
  float ptpairCut_eta[4]; 
  float s4s9minCut_eta[4]; 
  float s9s25minCut_eta[4]; 
  float isoCut_eta[4]; 
  float ptminCutHigh_eta[4] = {1E6,1E6,1E6,1E6};
  float ptpairCutHigh_eta[4] = {1E6,1E6,1E6,1E6};
  
  
  if( pizEta ==2){
    ptminCut_eta[0] = 1.0; 
    ptminCut_eta[1] = 1.0; 
    ptminCut_eta[2] = 0.7; 
    ptminCut_eta[3] = 0.7; 
    
    ptpairCut_eta[0] = 3.0; 
    ptpairCut_eta[1] = 3.0; 
    ptpairCut_eta[2] = 3.0; 
    ptpairCut_eta[3] = 3.0; 
    
    s4s9minCut_eta[0] = 0.9;
    s4s9minCut_eta[1] = 0.9;
    s4s9minCut_eta[2] = 0.9;
    s4s9minCut_eta[3] = 0.9;
    
    s9s25minCut_eta[0] = 0.85; 
    s9s25minCut_eta[1] = 0.85; 
    s9s25minCut_eta[2] = 0.85; 
    s9s25minCut_eta[3] = 0.85; 
    
    isoCut_eta[0] = 0.5; 
    isoCut_eta[1] = 0.5; 
    isoCut_eta[2] = 0.5; 
    isoCut_eta[3] = 0.5; 
  }
  
  if( pizEta ==1){
    ///right now. |eta|<2.5 pt >0.6 and ptpair >2.5 are online
    /// |eta|>2.5 , pt>0.5 and 1 < ptpair < 2.5 
    /// try upper cut for 2<|eta| <2.5 as well
    ptminCut_eta[0] = 0.6; 
    ptminCut_eta[1] = 0.6; 
    ptminCut_eta[2] = 0.5; 
    ptminCut_eta[3] = 0.5; 
    
    ptpairCut_eta[0] = 2.5; 
    ptpairCut_eta[1] = 2.5; 
    ptpairCut_eta[2] = 1.; 
    ptpairCut_eta[3] = 1.; 
      
    ptpairCutHigh_eta[0] = 999; 
    ptpairCutHigh_eta[1] = 4; 
    ptpairCutHigh_eta[2] = 2.5; 
    ptpairCutHigh_eta[3] = 2.5; 
    
    s4s9minCut_eta[0] = 0.9; 
    s4s9minCut_eta[1] = 0.9;
    s4s9minCut_eta[2] = 0.9;
    s4s9minCut_eta[3] = 0.9;
    
    s9s25minCut_eta[0] = -1; 
    s9s25minCut_eta[1] = -1; 
    s9s25minCut_eta[2] = -1; 
    s9s25minCut_eta[3] = -1; 
    
    isoCut_eta[0] = 0.5; 
    isoCut_eta[1] = 0.5; 
    isoCut_eta[2] = 0.5; 
    isoCut_eta[3] = 0.5; 
    
  }
  
  filename = TString(Form("testCalibv1.dflag%d.pe%d.step%d.iter%d.r%d.root",dataflag,pizEta, stepc,iter,evtRange));
  cout<<filename<<endl;
  TFile *fnew = new TFile(filename,"recreate");
  filename = TString(Form("testCalibv1.dflag%d.pe%d.step%d.iter%d.r%d.txt",dataflag,pizEta, stepc,iter ,evtRange));
  txtout.open(filename,ios::out);
  
  TH1F *hh_mpair_etaRing[2][40];
  for(int k=0; k<2; k++){
    for(int j=0; j< kEndcEtaRings; j++){
      filename = TString(Form("hh_mpair_etaRing_%d_%d",k,j));
      hh_mpair_etaRing[k][j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
    }
  }
  

  float res[10];
  
  
  cout<<"totalEntries: "<< dataflag<<" "<< totalEntries <<endl; 

  
  int start_entry  = 0; 
  int end_entry = totalEntries  -1; 
  
  
  cout<<"start/end entry: "<< test_evtRange <<" "<< start_entry <<" "<< end_entry <<endl; 
  cout<<"nbinMax  " << nbinMax <<endl; 

  //sigma_side[0] = 0.022; 
  //sigma_side[1] = 0.022; 
  
  ///mean_side is calculated inside corr_endcap.cc
  meanMass = 0.5*( mean_side[0] + mean_side[1]);
  
  cout<<"mean_sigma_side: "<<mean_side[0]<<" "<< mean_side[1] <<" "<< sigma_side[0]<<" "<< sigma_side[1]<< " mean " << meanMass <<endl; 
  
  int nEventsCount = 0;


  ///every time when run the code please check if you have the updated file 
  vector<string> certfiles;
  string certfile = string(workingDirectory) + string("/ecalGoodesGoodLumiBlocks.txt");
  certfiles.push_back(certfile);
  getLSrangeofEachRuns(certfiles);
  int curLumiBlock = -1;
  int curRun = -1;
  bool goodCurLumiBlock = false;


  for(entry = start_entry; entry <= end_entry; entry++){
    fChain->GetEntry(entry);
    nEventsCount ++; 
  
    
    vector<int>::iterator it = find(goodRunList.begin(),goodRunList.end(),runNumber);
    if( it == goodRunList.end()){
      continue;
    }
    if( curLumiBlock != lumiBlock || curRun != runNumber){ /// a new lumiBlock  or starting of a new Run                  
      curLumiBlock = lumiBlock;
      curRun =  runNumber;
      goodCurLumiBlock = checkLumiBlockofRun();  //check this lumiBlock                                                   
    }
    if( ! goodCurLumiBlock) continue;
    




    if(entry % 100000==0) cout<<"entry: "<<entry<<" "<<mpair<<endl; 
    int ieta1 = ietaXtalClus1[0];
    int ieta2 = ietaXtalClus2[0];
    int iphi1 = iphiXtalClus1[0];
    int iphi2 = iphiXtalClus2[0];
    
    int izz1 = izXtalClus1 <0 ? 0: 1;
    int izz2 = izXtalClus2 <0 ? 0: 1;
    
    int ietaRegion = -1; 
    if(fabs(etapair) < 2) ietaRegion = 0; 
    else if( fabs(etapair) < 2.5) ietaRegion = 1; 
    else if( fabs(etapair) < 2.8) ietaRegion = 2; 
    else ietaRegion = 3; 
    
    ////selections
    if( s4s9min <  s4s9minCut_eta[ietaRegion] ) continue; 
    if( ! ( ptmin > ptminCut_eta[ietaRegion] && ptmin <  ptminCutHigh_eta[ietaRegion] ) ) continue; 
    if( ! ( ptpair >  ptpairCut_eta[ietaRegion] && ptpair < ptpairCutHigh_eta[ietaRegion] )) continue; 
    if( isolation > isoCut_eta[ietaRegion] ) continue; 
    if( pizEta==2 && s9s25min <  s9s25minCut_eta[ietaRegion]) continue; 
    
    ///C0 is the updated calconst after each step 
    calcPairClusterUpdated(C0,res);
    if( res[0] <=0) continue; 
    float mpair_new = res[0];
    float ensum1 = res[1];
    float ensum2 = res[2];
    
    int bin = int( (mpair_new - xLowLimit )/binwidth); 
    if( mpair_new >= xLowLimit && bin>=0 && bin < nbinMax ){
      nCountedEE[izz1][ieta1][iphi1][bin] ++; 
      nCountedEE[izz2][ieta2][iphi2][bin] ++; 
    }
    
    int ietaring1 = iRingEndCap(izXtalClus1,ietaXtalClus1[0],iphiXtalClus1[0]);
    int ietaring2 = iRingEndCap(izXtalClus2,ietaXtalClus2[0],iphiXtalClus2[0]);
    hh_mpair_etaRing[izz1][ietaring1]->Fill(mpair_new);
    hh_mpair_etaRing[izz2][ietaring2]->Fill(mpair_new);
    
    
    
    int izside = izXtalClus1 <0 ? 0: 1;
    
    if( mpair_new > mean_side[izside] -2.*sigma_side[izside] && mpair_new < mean_side[izside] + 2.*sigma_side[izside] ){
      
      for(Int_t ixtal=0; ixtal  < nxtClus1; ixtal++){
	int ieta = ietaXtalClus1[ixtal];
	int iphi = iphiXtalClus1[ixtal];
	int izz = izXtalClus1 <0 ? 0: 1; 
	float en = eXtalClus1[ixtal]*C0[izz][ieta][iphi];
	int ietaRing = iRingEndCap(izXtalClus1,ietaXtalClus1[ixtal],iphiXtalClus1[ixtal]); 
	if( stepc ==2 ){
	  en *= corrfactorEtaRings[izz][ietaRing];
	}
	double wt = en / ensum1;
	if( wt >0.01){
	  WTSUM[izz][ieta][iphi] += pow(wt,2); 
	  CORSUM[izz][ieta][iphi] += pow(wt,2) * pow( meanMass/ mpair_new,2);
	}
      }
      for(Int_t ixtal=0; ixtal  < nxtClus2; ixtal++){
	int ieta = ietaXtalClus2[ixtal];
	int iphi = iphiXtalClus2[ixtal];
	int izz = izXtalClus2 <0 ? 0: 1; 
	float en = eXtalClus2[ixtal]*C0[izz][ieta][iphi];
	int ietaRing = iRingEndCap(izXtalClus2,ietaXtalClus2[ixtal],iphiXtalClus2[ixtal]); 
	if( stepc == 2 ){
	  en *= corrfactorEtaRings[izz][ietaRing];
	}
	double wt = en / ensum2;
	if( wt >0.01){
	  WTSUM[izz][ieta][iphi] += pow(wt,2); 
	  CORSUM[izz][ieta][iphi] += pow(wt,2) * pow( meanMass/ mpair_new,2);
	}
      }
      
    }
    
  }
  

  for(int iz = 0; iz<2; iz++){
    for(int j= 0; j<101; j++){
      for(int k= 0; k<101; k++){
	if( validRecHitEndCap[iz][j][k] <1) continue; 
	txtout<<"j "<<iz<<" "<<j<<" "<<k<<" "; 
	
	if(stepc==1){
	  for(int n =1; n<= nbinMax; n++){
	    txtout<<nCountedEE[iz][j][k][n-1]<<" "; 
	  }
	  txtout<<endl; 
	}else if( stepc==2){
	  txtout<<WTSUM[iz][j][k]<<" "<<CORSUM[iz][j][k]<<endl; 
	}else{
	  cout<<"stepc NA!"<<endl; 
	  return; 
	}
      }
    }
  }
  
  
  txtout<<"totalEntries: "<< dataflag<<" "<< totalEntries << " "<< nEventsCount <<endl; 
  txtout<<"start/end entry: "<< test_evtRange <<" "<< start_entry <<" "<< end_entry <<endl; 
  
  fnew->Write();
  fnew->Close();
  
  
}
