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



void split_job(char *inputList= "list.txt", int nEventEach = 1E7, int startEvtRange = 0 , int pizEta = 1, 
	       int barrelEndcap=1, char *name = "pizeb2012CGRPV39"){
  
  TString filename; 
    
  ifstream inputcc(inputList,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << inputList<<endl; 
    exit(1);
  }
  
  
  string tempfile; 
  
  filename = TString("dataChain/datachain_") + TString(name) + ".cc";
  ofstream txtout(filename,ios::out);
  bool firstChain = true; 
  
  int nchain = startEvtRange; 
  

  TChain *fChain ; 


  txtout<<"void datachain_"<<name<<"(){"<<endl; 

  bool closed = false; 
  
  while (inputcc.good()){
    
    inputcc >> tempfile; 

    if( inputcc.eof()) break; 
    
    int nrun = tempfile.find("run");

    if(nrun <=0){
      cout<<"nruN: "<< nrun <<" "<< tempfile.c_str()<<endl; 
      exit(1);
    }
    
    string run = tempfile.substr(nrun+3,6); 
    int runNumber = atoi(run.c_str());
    //here the runumber you might want to change 
    if( runNumber< 198346){
      filename = TString("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/11/store/user/wangj/data/pizdata/") + TString(run) + TString("/") +  TString(tempfile.c_str());
    }else{
      filename = TString("dcap://cmsgridftp.fnal.gov:24125/pnfs/fnal.gov/usr/cms/WAX/resilient/cmorgoth/data/pizdata/") + TString(run) + TString("/") + TString(tempfile.c_str());
    }
    
    cout<<filename<<endl; 
    
    if( firstChain){
      
      if(barrelEndcap==1){
	if(pizEta==1){
	  fChain = new TChain("pizSelb");
	}else{
	  fChain = new TChain("etaSelb");
	}
      }else{
	if(pizEta==1){
          fChain = new TChain("pizSele");
        }else{
          fChain = new TChain("etaSele");
        }
      }
      
      
      txtout<<"if( evtRange == "<<nchain+1<<") {"<<endl; 
      firstChain = false; 

      
      closed = false; 
    }
    
    
    txtout<<"fChain->Add(\""<<filename<<"\");"<<endl;
    
    fChain->Add(filename); 
    
    int nevents = fChain->GetEntries();
    
    if( nevents > nEventEach) {
      closed = true; 
      txtout<<"};"<<endl; 
      
      firstChain = true; 
      
      fChain->Delete();
      
      cout<<"nchain: "<< nchain+1 << " "<< nevents <<endl; 
      nchain ++; 
    }
    
  }

  if(!closed){
    txtout<<"};"<<endl;
  }
  txtout<<"}" <<endl; 
  
  cout<<"total chains: "<< nchain <<endl; 
  
}
