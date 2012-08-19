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

static const int ntt_tot = 2448; 
static const int nsm_tot = 36; 

TH1F *hhEnCorr[100]; 
TH1F *hhEnCorr_ietaComb[100]; 

int debug = 1; 

///for L3 method, 
double WTSUM[170][360];  /// wt = pow( en_crystal/ en3x3,2);
double CORSUM[170][360];  // wt * (mean/peak)*(mean*peak)

float C0[170][360];

//Common files
#include "../Common/variables.h"
#include "../Common/utils.cc"
#include "../Common/PositionCalc.cc"
#include "../Common/getGoodLS.cc"

#include "setBranchAddress.cc"


#include "datachain.cc"
#include "calibUtils.cc"
#include "correction.cc"
#include "getCrystaldeadflag.cc"


double binwidth = 0.005; 
double xHighLimit = 0.3; 
double xLowLimit = 0; 
int nbinMax = int(xHighLimit/binwidth + 0.1);

void testCalibv1(int test_dataflag,int test_pizEta, int test_calibStep, int test_calibIteration, int test_evtRange){
  workingDirectory = "MyJobWorkingDirectory";
  cout.precision(10);
  
  doBarrel =1; 
  pizEta = test_pizEta;
  
  
  if(pizEta!=1 && pizEta!=2){
    cout<<"pizEta 1 or 2  " << pizEta <<endl; 
    return; 
  }
    
  
  if( pizEta==2){ //this is for the output text 
    binwidth = 0.01; 
    xHighLimit = 0.8; 
    xLowLimit = 0.3; 
  }
  
  nbinMax = int( (xHighLimit-xLowLimit) /binwidth + 0.0001);
  cout<<"binwidth: "<< binwidth <<" "<< xLowLimit <<" "<< xHighLimit <<" "<<nbinMax <<endl; 
  
  
  cout<<"initilization counts.." <<endl;
  for(int j=0; j<170; j++){
    for(int k=0; k<360; k++){
      for(int n=0; n< nbinMax+1; n++){
        nCounted[j][k][n] = 0;
      }
    }
  }
  initialize_barrel();
  dataflag = test_dataflag; 
  evtRange = test_evtRange; 
  stepc = test_calibStep; 
  iter = test_calibIteration; 
  if(stepc <1 || stepc >4){
    cout<<"calibration step 1 to 4  " << stepc <<endl; 
    return; 
  }
  if(stepc==1 && iter>=2){
    cout<<"step one run 1 iteraction only"<<endl; 
    return; 
  }
  
  if(stepc==4){ ///last step apply dead crystal correction
    getCorrFactorDead_Barrel_CorrectionOnClusterEnergy();
  }
  
  if(stepc>=2){ ///after 1st step, get the SM scale correction 
    getSMcorrFactor(); 
  }
  
  if(stepc==2 && iter>=2){
    getcorrFactorEtaOfEachStep(); ///eta-correction step 
  }
  if(stepc==3 &&  iter>=1){
    getcorrFactorEta();   ///now get the eta-correction 
    if( iter>=2){ 
      getcorrFactorPhiOfEachStep(); //phi-correction step 
    }
  }
  
  if( stepc==4){
    getcorrFactorEta();
    getcorrFactorPhi();

    if(iter>=2){ //the final IC step 
      getcorrFactorIetaIphiOfEachStep();
      
      //now update C0[][] 
      for(int j=0; j<170; j++){
	for(int k=0; k< 360; k++){
	  C0[j][k] = corrfactoriEtaiPhi[j][k]; 
	}
      }
    }
  }

  ///x/y/z of each crystal to re-calculate X/Y/Z/ of cluster
  get_xyzEBrechits(); 
  
  if(pizEta==1){
    fChain = new TChain("pizSelb");
  }else{
    fChain = new TChain("etaSelb");
  }
  
  addDataChain();

  ///stepc==3 iter==11 to derive the dead crystal's correction
  // makeDeadCrystalFlag.C has to be runn
  if( (stepc==3 && iter==11)||stepc==4){ 
    getCrystaldeadflagBarrel(); 
  }
  
  setBranchAddress();
  
  totalEntries = fChain->GetEntries();
  cout<<"dataflag "<< dataflag <<" evtRange  "<< evtRange <<" totalEvents "<< totalEntries <<" Nfiles "<< fChain->GetNtrees()<<endl; 

  TString filename; 
  filename = TString(Form("testCalibv1.dflag%d.pe%d.step%d.iter%d.r%d.root",dataflag,pizEta, stepc,iter,evtRange));
  cout<<filename<<endl; 
  TFile *fnew = new TFile(filename,"recreate");
  filename = TString(Form("testCalibv1.dflag%d.pe%d.step%d.iter%d.r%d.txt",dataflag,pizEta, stepc,iter,evtRange));
  txtout.open(filename,ios::out);
  
  if(stepc<=3 ){
    txtout.precision(10);
  }else{
    txtout.precision(16);
  }
  
  float res[10];
  
  ///mass peak on SM 
  TH1F *hh_mpair_sm[36]; /// seed in the same SM
  for(int j=0; j<36; j++){
    filename = TString(Form("hh_mpair_sm_%d",j));
    hh_mpair_sm[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
  }
  
  int start_entry  = 0; 
  int end_entry = totalEntries  -1; 
  cout<<"start/end entry: "<< test_evtRange <<" "<< start_entry <<" "<< end_entry <<endl; 

  if(stepc==4){
    cout<<"L3 window pizEta: "<< pizEta <<" mean "<< meanMass <<" sigma "<< sigmaMass <<endl; 
  }

  int nEventsCount = 0; 

  ///end_entry = 100000;


  ///every time when run the code please check if you have the updated file 
  vector<string> certfiles;
  string certfile = string(workingDirectory) + string("/ecalGoodLumiBlocks.txt");
  certfiles.push_back(certfile);
  getLSrangeofEachRuns(certfiles);
  int curLumiBlock = -1;
  int curRun = -1;
  bool goodCurLumiBlock = false;
  


  for(entry = start_entry; entry <= end_entry; entry++){

    if(entry%100000==0) cout<<"entry " << entry <<endl; 
    
    
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
    



    int ieta1 = ietaXtalClus1[0];
    int iphi1 = iphiXtalClus1[0];
    convxtalid(iphi1,ieta1);
    
    int ieta2 = ietaXtalClus2[0];
    int iphi2 = iphiXtalClus2[0];
    convxtalid(iphi2,ieta2);
    
    if( (stepc==3 && iter==11) || stepc==4 ){ 
      ////check dead flag 
      if( flag_ietaiphi[ieta1+85][iphi1] >0) continue; 
      if( flag_ietaiphi[ieta2+85][iphi2] >0) continue; 
    }
    
    int ism1 = convertIetaIphiToSMNumber(ietaXtalClus1[0],iphiXtalClus1[0]);
    int ism2 = convertIetaIphiToSMNumber(ietaXtalClus2[0],iphiXtalClus2[0]);
    
    ///C0 is the updated calconst after each step 
    calcPairClusterUpdated(C0,res);
    //    cout<<"res[0] " << res[0] <<" "<<mpair <<endl; 

    if( res[0] <=0){
      continue; ///this barely happen. 
    }
    
    
    float mpair_new = res[0];
    float ensum1 = res[1];
    float ensum2 = res[2];
    
    if( ism1 == ism2){ ///to derive the SM correction 
      hh_mpair_sm[ism1-1]->Fill(mpair_new);
    }
    
    int bin = int( (mpair_new - xLowLimit )/binwidth); 
    if( mpair_new >= xLowLimit && bin>=0 && bin < nbinMax ){
      nCounted[ieta1+85][iphi1][bin] ++; 
      nCounted[ieta2+85][iphi2][bin] ++; 
    }
    
    if(stepc==4 && mpair_new > meanMass -2.*sigmaMass && mpair_new < meanMass + 2.*sigmaMass ){
      
      for(Int_t ixtal=0; ixtal  < nxtClus1; ixtal++){
	int ieta = ietaXtalClus1[ixtal];
	int iphi = iphiXtalClus1[ixtal];
	convxtalid(iphi,ieta); /// now iphi[0,359]; ieta[-85,84];
	if( (stepc==3 && iter==11) || stepc==4 ){
	  if( flag_ietaiphi[ieta+85][iphi] > 0) continue; 
	}
	float en = eXtalClus1[ixtal]*C0[ieta+85][iphi];
	
	if(stepc >=2){
	  en *= corrfactorEta[ieta+85];
	  en *= corrfactorPhi[iphi];
	  int ism = convertIetaIphiToSMNumber(ietaXtalClus1[ixtal],iphiXtalClus1[ixtal]);
	  en *= corrfactorSM[ism-1]; 
	}
	
	
	double wt = en / ensum1;
	if( wt >0.01){
	  WTSUM[ieta+85][iphi] += pow(wt,2); 
	  CORSUM[ieta+85][iphi] += pow(wt,2) * pow( meanMass/ mpair_new,2);
	}
	
      }
      
      for(Int_t ixtal=0; ixtal  < nxtClus2; ixtal++){
	int ieta = ietaXtalClus2[ixtal];
	int iphi = iphiXtalClus2[ixtal];
	convxtalid(iphi,ieta); /// now iphi[0,359]; ieta[-85,84];
	if( (stepc==3 && iter==11) || stepc==4 ){
	  if( flag_ietaiphi[ieta+85][iphi] > 0) continue;
	}
	
	float en = eXtalClus2[ixtal]*C0[ieta+85][iphi];
	if(stepc >=2){
	  en *= corrfactorEta[ieta+85];
	  en *= corrfactorPhi[iphi];
	  int ism = convertIetaIphiToSMNumber(ietaXtalClus2[ixtal],iphiXtalClus2[ixtal]);
	  en *= corrfactorSM[ism-1]; 
	}
	double wt = en / ensum2;
	if( wt >0.01){
	  WTSUM[ieta+85][iphi] += pow(wt,2);
  	  CORSUM[ieta+85][iphi] += pow(wt,2) * pow( meanMass/ mpair_new,2);
	}
      }
    }//inside mass peak
    
  }
  
  ///now print the counts in each bin of the mass distribution into the text file, 
  for(int j=0; j<170; j++){
    for(int k=0; k< 360; k++){
      txtout<<"j "<<j<<" "<<k<<" "; 
      if(stepc<=3){
	for(int n =1; n<= nbinMax; n++){
	  txtout<<nCounted[j][k][n-1]<<" "; 
	}
	txtout<<endl; 
      }
      else{
	txtout<<WTSUM[j][k]<<" "<<CORSUM[j][k]<<endl;
      }
    }
  }
  
  
  txtout<<"totalEntries: "<< dataflag<<" "<< totalEntries <<" "<< nEventsCount <<endl; 
  txtout<<"start/end entry: "<< test_evtRange <<" "<< start_entry <<" "<< end_entry <<endl; 
  
  fChain->Delete();
  

  fnew->Write();
  fnew->Close();
  
  
}
