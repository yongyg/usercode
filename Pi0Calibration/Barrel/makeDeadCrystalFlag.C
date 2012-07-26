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

#include "../Common/variables.h"
#include "../Common/utils.cc"


void checkHowManyDeadCrystalNext(int j, int k, int flagmap_ietaiphi[170][360],int &ncorner, int &nside, int &nflagdeadNextTo){
  
  ///input i, j after convxvalid 
  
  if( !( j>=0 && j<= 169 && k>=0 && k<= 359)){
    cout<<"input checkHowManyDeadCrystalNext j/k: "<< j<<" "<<k<<endl; 
    exit(1);
  }
  

  int nbad = 0; 
  nside = 0;
  ncorner = 0;
  
  
  
  

  for(int nieta = -1; nieta <=1; nieta++){
    int jj= j+ nieta; 
    if( jj <0 || jj >=170) continue; 
    
    for(int niphi = -1; niphi <=1; niphi++){
      
      if(nieta==0 && niphi==0) continue; //itself
      

      int kk = k + niphi; 
      if( kk == -1) kk = 359; 
      if( kk == 360) kk = 0; 
      
      
      if(flagmap_ietaiphi[jj][kk]>0 ){
	nbad ++;
	int deta = diff_neta(j,jj);
	int dphi = diff_nphi(k,kk);
	
	if( deta ==1 && dphi ==1) ncorner ++;
	else nside ++;
	
      }
    }
    
  }
  
  nflagdeadNextTo = 0; 

  if(nbad==0){
    nflagdeadNextTo = 0; 
  }
  else if(nbad==1){
    if(ncorner ==1) {
      nflagdeadNextTo = 1; 
    }else{
      nflagdeadNextTo = 2; 
    }
  }else if(nbad==2){
    if(ncorner ==2) {
      nflagdeadNextTo = 3; 
    }else if( nside ==2) { 
      nflagdeadNextTo = 4; 
    }else { 
      nflagdeadNextTo = 5; 
    }
  }else if(nbad==3){
    if(ncorner ==3) {
      nflagdeadNextTo = 6; 
    }else if( nside ==3) { 
      nflagdeadNextTo = 7; 
    }else if( ncorner ==1 && nside ==2) {
      nflagdeadNextTo = 8; 
    }else if( ncorner ==2 && nside ==1) {  
      nflagdeadNextTo = 9; 
    }
  }else if( nbad==4){
    if(ncorner ==4 ) {
      nflagdeadNextTo = 10; 
    }else if( nside ==4) {
      nflagdeadNextTo = 11; 
    }else if( ncorner ==2 && nside ==2) {  
      nflagdeadNextTo = 12; 
    }else if( ncorner ==1 && nside ==3) {
      nflagdeadNextTo = 13; 
    }else if( ncorner ==3 && nside ==1) {  
      nflagdeadNextTo = 14;
    }
  }else{
    cout<<"nbad? "<< nbad <<" "<< j<<" "<<k<<endl; 
    exit(1); 
  }
    
  
}


/// code to make deadflag map for each crystal in the Barrel 
/// simply to check the count of n.b. of pi0/eta candidates seeded in each crystal. 
/// should work OK when use all statistics for a given calibration period. 

void makeDeadCrystalFlag(char *inputfile, int test_dataflag){
  dataflag = test_dataflag;
    
  TFile *ff = new TFile(inputfile,"read");
  TString output = TString(Form("crystal_deadflag_eb_dflag%d.txt",dataflag));
  cout<<output<<endl; 
  txtout.open(output,ios::out);
  
  
  TH2F *hhtmp = (TH2F*)ff->Get("hh_res_count_0");
  int flagmap_ietaiphi[170][360]; 
  
  for(int j=0; j<170; j++){
    int bx = j-85;
    if( bx>=0) bx +=1;
    for(int k=0; k< 360; k++){
      
      int by = k;
      if( k== 0) by = 360;
      
      int count = hhtmp->GetBinContent(bx+85+1,by);
      
      if( count >0 && count < 50) {
	cout<<"pls check small count!!! " << bx <<" "<<by <<" "<< count <<endl; 
      }
      
      if( count >10){
	flagmap_ietaiphi[j][k] = 0; 
      }else{
	flagmap_ietaiphi[j][k] = 1; 
      }
      
      
    }
  }
    
  int ncorner; 
  int nside; 
  int nflagdeadNextTo;
  

  for(int ieta=-85; ieta<=85; ieta++){
    if( ieta == 0) continue; 
    for(int iphi=1; iphi<=360; iphi++){
      int j = ieta; 
      int k = iphi; 
      convxtalid(k,j);
      j += 85; 
      
      if(flagmap_ietaiphi[j][k] <1){
	checkHowManyDeadCrystalNext(j,k, flagmap_ietaiphi, ncorner, nside, nflagdeadNextTo); 
	txtout<<j<<" "<<k<<" "<< flagmap_ietaiphi[j][k] <<" "<< ncorner <<" "<< nside << " "<< nflagdeadNextTo <<endl;
      }else{
	txtout<<j<<" "<<k<<" "<< flagmap_ietaiphi[j][k] <<" "<< -1 <<" "<< -1 <<" "<< -1 <<endl;
      }
    }
  }
  
  
  
}
