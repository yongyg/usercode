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

//Common files                                                                                         
#include "../Common/variables.h"
#include "../Common/utils.cc"
#include "../Common/fitpeak.cc"

#include "foldEndcap.cc"
#include "getpreCalibConst.cc"
#include "correction.cc"


double binwidth = 0.005; 
double xHighLimit = 0.3; 
double xLowLimit = 0;
int nbinMax = int(xHighLimit/binwidth);


TH1F *hh_mpair_etaRing[2][40];
TH1F *hh_mpair1_etaRing[2][40];

///for L3 method,
double WTSUM[2][101][101];  /// wt = pow( en_crystal/ en3x3,2);
double CORSUM[2][101][101];  // wt * (mean/peak)*(mean*peak)


void deriveCalibConst(int test_dataflag,int test_pizEta, int test_calibStep, int test_calibIteration, int nEventRange){
  workingDirectory = "MyJobWorkingDirectory";
  cout.precision(10);

  
  get_xyzEBrechits();
  setEtaRingBoundaryEndcap();
  readInterCalibEndcap_GR09_V8();

  dataflag = test_dataflag; 
  pizEta  = test_pizEta; 
  stepc = test_calibStep;
  iter = test_calibIteration;
    
  
  if( pizEta==2){ //this is for the output text 
    binwidth = 0.01; 
    xHighLimit = 0.9; 
    xLowLimit = 0.2; 
    
  }
  
  nbinMax = int( (xHighLimit-xLowLimit) /binwidth + 0.0001);
  
  cout<<"binwidth: "<< binwidth <<" "<< xLowLimit <<" "<< xHighLimit <<" "<<nbinMax <<endl; 
  
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
  for(int iz = 0; iz<2; iz++){
    for(int j=0; j<101; j++){
      for(int k=0; k<101; k++){
	WTSUM[iz][j][k] = 0;
	CORSUM[iz][j][k] = 0;
      }
    }
  }
  
  if(stepc == 1 && iter >=2){ ///this is to dervie eta correction. 
    getcorrFactorEtaOfEachStep();
  }
  if( stepc ==2  ){
    getcorrFactorEta(); 
    if(iter>=2){
      getcorrFactorIzIxIyOfEachStep();
    }
  }
  
  
  TString filename; 
  TString checkfilename; 
  
  string tmps; 
  int ieta,iphi;
  int binCont; 
  int izz; 

  double wtsum;
  double corsum; 
  TString filenamepart = TString(Form("dflag%d.pe%d.step%d.iter%d",dataflag,pizEta,stepc,iter));
  for(int n=1; n<= nEventRange; n++){
    filename = TString("calibres/testCalibv1.") + filenamepart + TString(Form(".r%d.txt",n));
    checkfilename = TString("ls calibres/testCalibv1.") + filenamepart + TString(Form(".r%d.txt",n));
    if(gSystem->Exec(checkfilename)!=0){
      cout<<"file not found.."<<endl; 
      exit(1); 
    }
    ifstream txtin(filename,ios::in);
    if (txtin.fail()) {
      cout<<"error "<< filename <<" can not be opened."<<endl; 
      exit(1);
    }
    
    while(txtin.good()){
      txtin>>tmps >> izz >> ieta>>iphi;
      if(stepc==1){
	for(int k=1; k<= nbinMax; k++){
	  txtin>>binCont; 
	  nCountedEE[izz][ieta][iphi][k-1] += binCont; 
	}
      }
      if(stepc==2){
	txtin >>  wtsum >> corsum; 
	WTSUM[izz][ieta][iphi] += wtsum; 
	CORSUM[izz][ieta][iphi] += corsum; 
      }
      if( izz ==1 && ieta == 100 && iphi ==60) break; //the last line 
    }
    txtin.close();
  }
  
  filename = TString("calibres/deriveCalibConst.") + filenamepart + TString(".root");
  cout<<filename<<endl; 
  TFile *fnew = new TFile(filename,"recreate");
  filename = TString("calibres/deriveCalibConst.") + filenamepart + TString(".txt");
  txtout.open(filename,ios::out);
  txtout.precision(10);
  
    ///for ieta/iphi correction pi0 
  double xlowFit = 0.04; 
  double xhighFit = 0.3; 
  int nPowFit = 4; 
  int pizFit = 0; 

  ///for alca [0.04,0.06] new cuts [0.04,0.23];
  if(pizEta==1 ){
    xlowFit = 0.06;
    xhighFit = 0.25;
    nPowFit = 2;
  }
  if(pizEta==2){
    xlowFit = 0.3; 
    xhighFit = 0.8; 
    nPowFit =3; 
    pizFit =1; 
  }
  
  
  cout<<"fit ragne: "<< xlowFit <<" "<< xhighFit <<" "<< nPowFit <<endl; 
  
  for(int k=0; k<2; k++){
    for(int j=0; j< kEndcEtaRings; j++){
      filename = TString(Form("hh_mpair_etaRing_%d_%d",k,j));
      hh_mpair_etaRing[k][j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit); 
    }
  }
  

  
  TH1F *hh_res_ieta[2][4];
  for(int k=0; k<2; k++){
    for(int j=0; j<4;j++){
      filename = TString(Form("hh_res_ieta_%d_%d",k,j));
      hh_res_ieta[k][j] = new TH1F(filename,filename,100,0,100);
    }
  }
  
  
  TH1F *hh_corr_etaRing[2]; 
  for(int k=0; k<2; k++){
    filename = TString(Form("hh_corr_etaRing_%d",k));
    hh_corr_etaRing[k] = new TH1F(filename,filename,kEndcEtaRings,0,kEndcEtaRings); 
  }

  TH1F *hhmpair = new TH1F("hhmpair","hhmpair",nbinMax,xLowLimit,xHighLimit);
  
  
  double res[10];

  nMaxRingIC = 20; 
  if(pizEta==2){
    nMaxRingIC = 30;
  }
  
  if(stepc==1){
    
    for(int iz=0; iz<2; iz++){
      for(int j=1; j< 101; j++){
	for(int k=1; k< 101; k++){
	  if( validRecHitEndCap[iz][j][k] <1 ) continue; 
	  for(int b=1; b<= nbinMax; b++){
	    hhmpair->SetBinContent(b,nCountedEE[iz][j][k][b-1]);
	  }
	  int ietaring = iRingEndCap(2*iz-1,j,k); 
	  hh_mpair_etaRing[iz][ietaring]->Add(hhmpair);
	}
      }
    }
    
    
    float mean_etaRing[2] = {0};

    for(int iz=0; iz<2; iz++){
      for(int j=0; j< kEndcEtaRings; j++){
	
	if( hh_mpair_etaRing[iz][j]->Integral()<1000) continue; 
	
	if(pizEta==1 && j<= nMaxRingIC ){ ////above 20 no calibration from pi0s 
	  pi0_mfitpeak(hh_mpair_etaRing[iz][j],pizFit,xlowFit,xhighFit,nPowFit,res,"","ieta");
	}else if(pizEta==2 && j<= nMaxRingIC){
	  pi0_mfitpeak(hh_mpair_etaRing[iz][j],pizFit,xlowFit,xhighFit,nPowFit,res,"","ieta");
	}

	for(int k=0; k<4;k++){
	  hh_res_ieta[iz][k]->SetBinContent(j+1,res[2*k]);
	  hh_res_ieta[iz][k]->SetBinError(j+1,res[2*k+1]);
	}
	if( j<= nMaxRingIC ){
	  mean_etaRing[iz] += res[0];
	}
	
      }      
    }
    
    for(int iz=0; iz<2; iz++){
      mean_etaRing[iz] /= (nMaxRingIC+1);
    }
    
    
    cout<<"mean_etaRing: "<< mean_etaRing[0] <<" "<< mean_etaRing[1] <<endl; 
    //in the end, print out the new calibration constants from this step 
    if( stepc==1 ){ /// derive eta correction from this step
      for(int iz=0; iz<2; iz++){
	for(int j=0; j< kEndcEtaRings ; j++){
	  if( j <= nMaxRingIC){
	    float tmp = corrfactorEtaRings[iz][j] * mean_etaRing[iz] / hh_res_ieta[iz][0]->GetBinContent(j+1);
	    hh_corr_etaRing[iz]->SetBinContent(j+1,tmp);
	  }else{
	    hh_corr_etaRing[iz]->SetBinContent(j+1,hh_corr_etaRing[iz]->GetBinContent(nMaxRingIC+1));
	  }
	}
      }
    }
  }
  
  if(stepc==2){


    for(int n=1; n<= nEventRange; n++){
      filename = TString("calibres/testCalibv1.") + filenamepart + TString(Form(".r%d.root",n));
      checkfilename = TString("ls calibres/testCalibv1.") + filenamepart + TString(Form(".r%d.root",n));
      if( gSystem->Exec(checkfilename)!=0){
	cout<<"not found."<<endl;
      exit(1);
      }
      TFile *fff = new TFile(filename,"read");
      for(int j=0; j<2; j++){
	for(int k=0; k<=nMaxRingIC; k++){
	  filename = TString(Form("hh_mpair_etaRing_%d_%d",j,k));
	  TH1F *hhtmp = (TH1F*)fff->Get(filename);
	  if(hhtmp!=NULL){
	    hh_mpair_etaRing[j][k]->Add(hhtmp);
	  }
	}
      }
      fff->Close();
    }
    for(int iz=0; iz<2; iz++){
      for(int j=0; j< kEndcEtaRings; j++){
	if( hh_mpair_etaRing[iz][j]->Integral()<1000) continue; 
	if(pizEta==1 && j<= nMaxRingIC ){ ////above 20 no calibration from pi0s 
	  pi0_mfitpeak(hh_mpair_etaRing[iz][j],pizFit,xlowFit,xhighFit,nPowFit,res,"","ieta");
	}else if(pizEta==2 && j<= nMaxRingIC){
	  pi0_mfitpeak(hh_mpair_etaRing[iz][j],pizFit,xlowFit,xhighFit,nPowFit,res,"","ieta");
	}
	for(int k=0; k<4;k++){
	  hh_res_ieta[iz][k]->SetBinContent(j+1,res[2*k]);
	  hh_res_ieta[iz][k]->SetBinError(j+1,res[2*k+1]);
	}
      }      
    }
    

    //printout the IC     
    for(int iz=0; iz<2; iz++){
      for(int j=0; j<101; j++){
	for(int k=0; k< 101; k++){
	  if( validRecHitEndCap[iz][j][k] <1) continue;
	  double tmp = -1; 
	  
	  int ietaRing = iRingEndCap(2*iz-1,j,k);
	  if( WTSUM[iz][j][k] > 0 && ietaRing <= nMaxRingIC ){
	    tmp = corrfactoriZiXiY[iz][j][k] * CORSUM[iz][j][k]/ WTSUM[iz][j][k];
	  }
	  if( j==0 && k<=1){
	    cout<<"j/k "<< j<<" "<<k<<" "<< corrfactoriZiXiY[iz][j][k] <<" " << CORSUM[iz][j][k] <<" " <<  WTSUM[iz][j][k] <<" "<< tmp <<endl; 
	  }
	  txtout<<iz<<" "<<j<<" "<<k<<" "<<tmp <<endl; 
	}
      }
    }
  }

  
  fnew ->Write();
  fnew ->Close();
  
  cout<<"end.." <<endl; 
  
}
