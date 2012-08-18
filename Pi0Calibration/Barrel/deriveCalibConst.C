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
#include "getCrystaldeadflag.cc"
#include "correction.cc"


double binwidth = 0.005; 
double xHighLimit = 0.3; 
double xLowLimit = 0;
int nbinMax = int(xHighLimit/binwidth);

///for L3 method, 
double WTSUM[170][360];    /// wt = pow( en_crystal/ en3x3,2);
double CORSUM[170][360];   // wt * (mean/peak)*(mean*peak)

void deriveCalibConst(int test_dataflag,int test_pizEta, int test_calibStep, int test_calibIteration, int nEventRange){
  workingDirectory = "MyJobWorkingDirectory";
  cout.precision(10);

  pizEta = test_pizEta;
  dataflag = test_dataflag; 
  stepc = test_calibStep; 
  iter = test_calibIteration; 
  if( pizEta ==2){
    binwidth = 0.01; 
    xHighLimit = 0.8; 
    xLowLimit = 0.3; 
  }
  
  nbinMax = int( (xHighLimit-xLowLimit) /binwidth + 0.0001);
  cout<<"binwidth: "<< binwidth <<" "<< xLowLimit <<" "<< xHighLimit <<" "<<nbinMax <<endl; 
  
  if( (stepc==3 && iter==11)|| stepc==4){
    getCrystaldeadflagBarrel(); 
  }
  
  

  for(int j=0;j<38; j++){
    for(int k=0;k<85; k++){
      corrfactorIetaSM[j][k] = 1; 
    }
  }

  for(int j=0; j< 170; j++){
    corrfactorEta[j] = 1; 
    corrfactorEtatb[j] = 1; 
    corrfactorEtaco[j] = 1; 
  }
  for(int j=0; j< 360; j++){
    corrfactorPhi[j] = 1; 
  }
    
  for(int j=0; j<36; j++){
    corrfactorSM[j] = 1; 
  }
    
  for(int j=0; j< 20; j++){
    corrfactorDead[j] = 1;
  }
  
  
  for(int j=0; j<170; j++){
    for(int k=0; k<360; k++){
      flagiEtaiPhi[j][k] = 0; 
      corrfactoriEtaiPhi[j][k] = 1;  //for step1 only actually 
      WTSUM[j][k] = 0;
      CORSUM[j][k] = 0;
    }
  }

  if(stepc==2 && iter>=2){
    getcorrFactorEtaOfEachStep(); ///eta-correction step 
  }

  if(stepc==3 &&  iter>=2){
    getcorrFactorPhiOfEachStep(); //phi-correction step 
  }
  
  if(stepc==4 &&  iter>=2){ ///ic- step 
    getcorrFactorIetaIphiOfEachStep(); 
  }
  
  
  TString filename; 
  TString checkfilename; 
  
  string tmps; 
  int ieta,iphi;
  int binCont; 
  
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
    
    int nXtal = 1; 
    while(txtin.good()){
      
      txtin>>tmps >> ieta>>iphi;
      if(stepc<=3){
	for(int k=1; k<= nbinMax; k++){
	  txtin>>binCont; 
	  nCounted[ieta][iphi][k-1] += binCont; 
	}
      }else{
	txtin >> wtsum >> corsum; 
	WTSUM[ieta][iphi] += wtsum; 
	CORSUM[ieta][iphi] += corsum; 
      }
      
      nXtal++; 
      if( nXtal > 61200) break; 
    }
    
    txtin.close();
  }
  
  filename = TString("calibres/deriveCalibConst.") + filenamepart + TString(".root");
  TFile *fnew = new TFile(filename,"recreate");
  filename = TString("calibres/deriveCalibConst.") + filenamepart + TString(".txt");
  txtout.open(filename,ios::out);
  txtout.precision(10);
  cout<<filename<<endl; 
  double xlowFit = 0.04; 
  double xhighFit = 0.3; 
  int nPowFit = 4; 
  int pizFit = 0; 

  ///for alca [0.04,0.06] new cuts [0.04,0.23];
  if(pizEta==1 ){
    xlowFit = 0.06;
    xhighFit = 0.22;
    nPowFit = 3;
  }
  
  if(pizEta==2){
    xlowFit = 0.35; 
    xhighFit = 0.75; 
    nPowFit =3; 
    pizFit =1; 
  }
  
  cout<<"fit ragne: "<< xlowFit <<" "<< xhighFit <<" "<< nPowFit <<endl; 
  
  
  TH1F *hhmpair = new TH1F("hhmpair","hhmpair",nbinMax,xLowLimit,xHighLimit);
  
  double res[10];
  
  //simply count n.b. of pi0 inside each crystal. 
  TH2F *hh_res_count[2]; 
  for(int j=0; j<2; j++){
    filename = TString(Form("hh_res_count_%d",j));
    hh_res_count[j] = new TH2F(filename,filename,171,-85,86,360,1,361);
  }
    
  TH1F *hh_mpair_ieta[170]; 
  TH1F *hh_mpair_ietatb[170]; 
  TH1F *hh_mpair_ietaco[170]; 
  
  for(int j=0; j<170; j++){
    filename = TString(Form("hh_mpair_ieta_%d",j));
    hh_mpair_ieta[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
    
    filename = TString(Form("hh_mpair_ietatb_%d",j));
    hh_mpair_ietatb[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
    filename = TString(Form("hh_mpair_ietaco_%d",j));
    hh_mpair_ietaco[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
  }

  TH1F *hh_res_sm[4]; 
  
  TH1F *hh_res_ietaSM[38][4];
  TH1F *hh_corr_ietaSM[38];

  TH1F *hh_corr_ieta = new TH1F("hh_corr_ieta","hh_corr_ieta",171,-85,86);
  TH1F *hh_corr_ietatb = new TH1F("hh_corr_ietatb","hh_corr_ietatb",171,-85,86);
  TH1F *hh_corr_ietaco = new TH1F("hh_corr_ietaco","hh_corr_ietaco",171,-85,86);
  
  
  for(int j=0; j<38; j++){
    filename = TString(Form("hh_corr_ietaSM_%d",j));
    hh_corr_ietaSM[j]=new TH1F(filename,filename,85,1,86);
    for(int k=0;k<4; k++){
      filename = TString(Form("hh_res_ietaSM_%d_%d",j,k));
      hh_res_ietaSM[j][k] = new TH1F(filename,filename,85,1,86);
    }
  }

  
  float xbinLow[35];
  for(int j=0; j< 35; j++){
    float xbin = -85 + 5*j;
    xbinLow[j] = xbin;
  }
  float xbinLow1[35];
  for(int j=0; j<= 17 ; j++){
    float xbin =  5 * j;
    xbinLow1[j] = xbin+1;
  }
  TH1F *hh_res_ietaTT[4];
  TH1F *hh_res_ietaTTAbs[4];
  for(int j=0; j<4; j++){
    filename = TString(Form("hh_res_ietaTT_%d",j));
    hh_res_ietaTT[j] = new TH1F(filename,filename,34,xbinLow);
    filename = TString(Form("hh_res_ietaTTAbs_%d",j));
    hh_res_ietaTTAbs[j] = new TH1F(filename,filename,17,xbinLow1);
  }


  TH1F *hh_res_ieta[4];
  TH1F *hh_res_ietatb[4];
  TH1F *hh_res_ietaco[4];
  
  for(int j=0; j<4; j++){
    filename = TString(Form("hh_res_ieta_%d",j));
    hh_res_ieta[j] = new TH1F(filename,filename,171,-85,86);
    filename = TString(Form("hh_res_ietatb_%d",j));
    hh_res_ietatb[j] = new TH1F(filename,filename,171,-85,86);
    filename = TString(Form("hh_res_ietaco_%d",j));
    hh_res_ietaco[j] = new TH1F(filename,filename,171,-85,86);
  }
  TH1F *hh_corr_iphi = new TH1F("hh_corr_iphi","hh_corr_iphi",360,1,361);
  TH1F *hh_corr_iphismtb = new TH1F("hh_corr_iphismtb","hh_corr_iphismtb",20,1,21);
  TH1F *hh_res_iphismtb[4];
  for(int j=0; j<4; j++){
    filename = TString(Form("hh_res_iphismtb_%d",j));
    hh_res_iphismtb[j] = new TH1F(filename,filename,20,1,21);
  }

  TH1F *hh_res_iphi[4];
  for(int j=0; j<4; j++){
    filename = TString(Form("hh_res_iphi_%d",j));
    hh_res_iphi[j] = new TH1F(filename,filename,360,1,361);
  }
  
  
  //// dead crystal's correction
  TH1F *hh_mpair_deadv1[15];
  for(int j=0; j< 15; j++){
    filename = TString(Form("hh_mpair_deadv1_%d",j));
    hh_mpair_deadv1[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
  }
  TH1F *hh_corr_deadv1 = new TH1F("hh_corr_deadv1","hh_corr_deadv1",15,0,15); 
  
  ///mass peak on SM 
  TH1F *hh_mpair_sm[36]; 
  for(int j=0; j<36; j++){
    filename = TString(Form("hh_mpair_sm_%d",j));
    hh_mpair_sm[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
  }
  
  ///on each TT or abs(TT)                                                                                                
  TH1F *hh_mpair_ietaTT[35];
  TH1F *hh_mpair_ietaTTAbs[17];
  for(int j=0; j< 35; j++){
    filename = TString(Form("hh_mpair_ietaTT_%d",j));
    hh_mpair_ietaTT[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
  }
  for(int j=0; j< 17; j++){
    filename = TString(Form("hh_mpair_ietaTTAbs_%d",j));
    hh_mpair_ietaTTAbs[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
  }
  
  
  ////ieta of each SM
  TH1F *hh_mpair_ietaSM[38][85]; //// 36 for tb, 37 for others
  for(int j=0;j<38; j++){
    for(int k=0; k<85; k++){
      filename = TString(Form("hh_mpair_ietaSM_%d_%d",j,k));
      hh_mpair_ietaSM[j][k] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
    }
  }
  ///combined testbeam SM mpair vs iphi 20 only
  TH1F *hh_mpair_iphi_smtb[20];
  for(int j=0; j< 20; j++){
    filename = TString(Form("hh_mpair_iphi_smtb_%d",j));
    hh_mpair_iphi_smtb[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
  }
  
  
  TH1F *hh_mpair_iphi[360]; 
  for(int j=0; j<360; j++){
    filename = TString(Form("hh_mpair_iphi_%d",j));
    hh_mpair_iphi[j] = new TH1F(filename,filename,nbinMax,xLowLimit,xHighLimit);
  }
  TH1F *hh_corr_sm = new TH1F("hh_corr_sm","correction to SM",36,1,37);
  for(int j=0; j<4; j++){
    filename = TString(Form("hh_res_sm_%d",j));
    hh_res_sm[j] = new TH1F(filename,filename,36,1,37);
  }
  
  double meancorr_sm = 0; 
  cout<<"adding hhmpair vs ieta "<<endl; 
  for(int n=1; n<= nEventRange; n++){
    filename = TString("calibres/testCalibv1.") + filenamepart + TString(Form(".r%d.root",n));
    checkfilename = TString("ls calibres/testCalibv1.") + filenamepart + TString(Form(".r%d.root",n));
    if( gSystem->Exec(checkfilename)!=0){
      cout<<"not found."<<endl; 
      exit(1); 
    }
    TFile *fff = new TFile(filename,"read");
    for(int j=0; j<36; j++){
      filename = TString(Form("hh_mpair_sm_%d",j));
      TH1F *hhtmp = (TH1F*)fff->Get(filename);
      hh_mpair_sm[j]->Add(hhtmp);
    }
    fff->Close();
  }
  
  cout<<"fitting mpair vs ism "<<endl; 
  
  for(int j=0; j<36; j++){
      
    int bx = j+1; 
    if( hh_mpair_sm[j]->Integral()<100){
      cout<<"check mpair_sm events too small!!!: "<< j<<" "<<hh_mpair_sm[j]->Integral() <<endl;
      exit(1);
    }
    pi0_mfitpeak(hh_mpair_sm[j],pizFit,xlowFit,xhighFit,nPowFit,res,"","ism");
    meancorr_sm += res[0];
    for(int n=0; n<4; n++){
      hh_res_sm[n]->SetBinContent(bx,res[2*n]);
      hh_res_sm[n]->SetBinError(bx,res[2*n+1]);
    }
  }
  meancorr_sm /= 36; 
  cout<<" Average peak poistion on all SMs:  "<< meancorr_sm<<endl; 
  for(int j=0; j<36; j++){
    int bx = j+1; 
    double tmp = corrfactorSM[j] *  meancorr_sm / hh_res_sm[0]->GetBinContent(bx) ;
    hh_corr_sm->SetBinContent(bx,tmp);
  }
  double mean_allIeta = 0; 
  double mean_allIetatb = 0; 
  double mean_allIetaco = 0; 
  double mean_allIphi = 0; 
  
  if(stepc<=3){
    for(int j=0; j<170; j++){
      int bx = j-85; 	
      if( bx>=0) bx +=1; 
      for(int k=0; k< 360; k++){
	for(int b=1; b<= nbinMax; b++){
	  hhmpair->SetBinContent(b,nCounted[j][k][b-1]);
	}
	hh_mpair_ieta[j]->Add(hhmpair); 
	hh_mpair_iphi[k]->Add(hhmpair);
	
	int itteta = j/5;
        hh_mpair_ietaTT[itteta]->Add(hhmpair);
        int ietaAbs = (abs(bx)-1)/5;
        hh_mpair_ietaTTAbs[ietaAbs]->Add(hhmpair);
	

	int by = k;
        if( k== 0) by = 360;
	
	int iSM = (by-1)/20+1;
	if( bx<0) iSM += 18;
	bool TBSM = iSM ==1 || iSM ==2 || iSM == 11 || iSM == 24  || iSM == 21  || iSM == 15; 
	if( TBSM){
	  hh_mpair_iphi_smtb[(by-1)%20]->Add(hhmpair);
	}
	hh_mpair_ietaSM[iSM-1][abs(bx)-1]->Add(hhmpair);
	if (TBSM){
	  hh_mpair_ietaSM[36][abs(bx)-1]->Add(hhmpair);
	  hh_mpair_ietatb[j]->Add(hhmpair);
	}else{
	  hh_mpair_ietaco[j]->Add(hhmpair);
	  hh_mpair_ietaSM[37][abs(bx)-1]->Add(hhmpair);
	}
	
	if(stepc==3 && iter ==11){ ///Only stepc3 and iter==11 to derive the dead crystal's correction
	  ///check crystals near dead crystals
	  if(flag_ietaiphi[j][k] < 1) { //itself is good
	    int ndeadflag = ndeadflag_ietaiphi[j][k]; 
	    if( ndeadflag <0){
	      cout<<"error.. ndeadflag: "<< ndeadflag<<" "<<j<<" "<<k<<endl; 
	      return; 
	    }
	    if(ndeadflag>=15){
	      cout<<"error pls check deadcrystal map ndeadflag > than defined!!" << ndeadflag <<endl; 
	      return; 
	    }
	    hh_mpair_deadv1[ndeadflag]->Add(hhmpair);
	  }
	}
	
	hh_res_count[0]->SetBinContent(bx+85+1,by,hhmpair->Integral()); /// all 
	hh_res_count[1]->SetBinContent(bx+85+1,by,hhmpair->Integral(int(0.1/binwidth),int(0.16/binwidth))); 
      }
      
      pi0_mfitpeak(hh_mpair_ieta[j],pizFit,xlowFit,xhighFit,nPowFit,res,"","ieta");
      mean_allIeta += res[0]; 
      for(int n=0; n<4; n++){
	hh_res_ieta[n]->SetBinContent(bx+85+1,res[2*n]);
	hh_res_ieta[n]->SetBinError(bx+85+1,res[2*n+1]);
      }
        
      pi0_mfitpeak(hh_mpair_ietatb[j],pizFit,xlowFit,xhighFit,nPowFit,res,"","ietatb");
      mean_allIetatb += res[0]; 
      for(int n=0;n<4;n++){
	hh_res_ietatb[n]->SetBinContent(bx+85+1,res[2*n]);
	hh_res_ietatb[n]->SetBinError(bx+85+1,res[2*n+1]);
      }
      
      pi0_mfitpeak(hh_mpair_ietaco[j],pizFit,xlowFit,xhighFit,nPowFit,res,"","ietaco");
      mean_allIetaco += res[0]; 
      for(int n=0; n<4;n++){
	hh_res_ietaco[n]->SetBinContent(bx+85+1,res[2*n]);
	hh_res_ietaco[n]->SetBinError(bx+85+1,res[2*n+1]);
      }
    } 
    
    mean_allIeta /= 170; 
    mean_allIetatb /= 170; 
    mean_allIetaco /= 170; 
    
    cout<<"Average peak positions of all Ieta: "<< mean_allIeta << " TBSM "<< mean_allIetatb <<" OtherSM "<< mean_allIetaco <<endl;
    
    
    for(int k=0; k< 17 ; k++){
      int bx = k+1;
      if( hh_mpair_ietaTTAbs[k]->Integral()<1000) continue;
      pi0_mfitpeak(hh_mpair_ietaTTAbs[k],pizFit,xlowFit,xhighFit,nPowFit,res,"","ietaTT");
      for(int n=0; n<4; n++){
        hh_res_ietaTTAbs[n]->SetBinContent(bx,res[2*n]);
        hh_res_ietaTTAbs[n]->SetBinError(bx,res[2*n+1]);
      }
    }
    for(int k=0; k< 34 ; k++){
      int bx = k+1;
      if( hh_mpair_ietaTT[k]->Integral()<1000) continue;
      pi0_mfitpeak(hh_mpair_ietaTT[k],pizFit,xlowFit,xhighFit,nPowFit,res,"","ietaTT");
      for(int n=0; n<4; n++){
        hh_res_ietaTT[n]->SetBinContent(bx,res[2*n]);
        hh_res_ietaTT[n]->SetBinError(bx,res[2*n+1]);
      }
    }
    


    double meanIetaSM[38]={0};

    cout<<"fitting peak of each ieta in each SM" <<endl; 
    ///peak vs ieta inside each SM. ( j==36--> SMTB,  j=37 --> others) 
    for(int j=0; j<38; j++){
      for(int k=0; k<85; k++){
	if( hh_mpair_ietaSM[j][k]->Integral()<1000) continue;
	pi0_mfitpeak(hh_mpair_ietaSM[j][k],pizFit,xlowFit,xhighFit,nPowFit,res,"","ietaSM");
	int by = k+1; 
	meanIetaSM[j] += res[0] / 85;
	for(int n=0; n<4; n++){
	  hh_res_ietaSM[j][n]->SetBinContent(by,res[2*n]);
	  hh_res_ietaSM[j][n]->SetBinError(by,res[2*n+1]);
	}
      }
    }

    ///mass peak vs iphi (SM TB ) 
    double mean_allIphismtb = 0; 
    for(int k=0; k<20; k++){
      if( hh_mpair_iphi_smtb[k]->Integral()<1000) continue; 
      pi0_mfitpeak(hh_mpair_iphi_smtb[k],pizFit,xlowFit,xhighFit,nPowFit,res,"","iphi");
      mean_allIphismtb += res[0] /20; 
      int by = k+1; 
      for(int n=0; n<4; n++){
	hh_res_iphismtb[n]->SetBinContent(by,res[2*n]);
	hh_res_iphismtb[n]->SetBinError(by,res[2*n+1]);
      }
    }
    cout<<"Average peak positions of all Iphi of SMTB "<< mean_allIphismtb <<endl;
    
    ///mas peak vs iphi  all SM .
    for(int k=0; k< 360; k++){
      int by = k; 
      if( k== 0) by = 360; 
      if(hh_mpair_iphi[k]->Integral()<1000) continue;
      pi0_mfitpeak(hh_mpair_iphi[k],pizFit,xlowFit,xhighFit,nPowFit,res,"","iphi");
      mean_allIphi += res[0];
      for(int n=0; n<4; n++){
	hh_res_iphi[n]->SetBinContent(by,res[2*n]);
	hh_res_iphi[n]->SetBinError(by,res[2*n+1]);
      }
    }
    mean_allIphi /= 360; 
    cout<<"Average peak positions of all Iphi: "<< mean_allIphi <<endl;
    ///correction vs ieta , ieta_SMTB, ieta_SMothers
    for(int j=0; j< 170; j++){
      int bx = j-85; 	
      if( bx>=0) bx +=1; 
      int binx = bx+ 85 + 1; 
      double tmp = corrfactorEta[j] * pow( mean_allIeta/ hh_res_ieta[0]->GetBinContent(binx),1);
      double tmpErr = tmp * hh_res_ieta[0]->GetBinError(binx)/ hh_res_ieta[0]->GetBinContent(binx);
      hh_corr_ieta->SetBinContent(binx,tmp);
      hh_corr_ieta->SetBinError(binx,tmpErr);
      tmp = corrfactorEta[j] * pow( mean_allIetatb/ hh_res_ietatb[0]->GetBinContent(binx),1);
      tmpErr = tmp * hh_res_ietatb[0]->GetBinError(binx)/ hh_res_ietatb[0]->GetBinContent(binx);
      hh_corr_ietatb->SetBinContent(binx,tmp);
      hh_corr_ietatb->SetBinError(binx,tmpErr);
      tmp = corrfactorEta[j] * pow( mean_allIetaco/ hh_res_ietaco[0]->GetBinContent(binx),1);
      tmpErr = tmp * hh_res_ietaco[0]->GetBinError(binx)/ hh_res_ietaco[0]->GetBinContent(binx);
      hh_corr_ietaco->SetBinContent(binx,tmp);
      hh_corr_ietaco->SetBinError(binx,tmpErr);
    }
    
    /// ieta of each SM 
    for(int j=0; j<38; j++){
      if( j==36) cout<<"Average peak positions of all Ieta in SMTB "<<j+1<<":  "<<meanIetaSM[j]<<endl;
      if( j==37) cout<<"Average peak positions of all Ieta in other SMs "<<j+1<<":  "<<meanIetaSM[j]<<endl;
      
      for(int k=0; k<85; k++){
	double tmp = corrfactorEta[k+85] * pow( meanIetaSM[j]/hh_res_ietaSM[j][0]->GetBinContent(k+1),1);
	hh_corr_ietaSM[j]->SetBinContent(k+1,tmp);
	double tmpErr = tmp * hh_res_ietaSM[j][0]->GetBinError(k+1)/hh_res_ietaSM[j][0]->GetBinContent(k+1);
	hh_corr_ietaSM[j]->SetBinError(k+1,tmpErr);
      }
    }
    
    ////iphi of SMTB 
    for(int k=1; k<=20; k++){
      double tmp = corrfactorPhi[k] * mean_allIphismtb / hh_res_iphismtb[0]->GetBinContent(k); 
      double tmpErr = tmp *  hh_res_iphismtb[0]->GetBinError(k)/ hh_res_iphismtb[0]->GetBinContent(k); 
      ///if( k==1 || k==20) cout<<"corrfactorPhi["<<k<<"]"<< corrfactorPhi[k] <<" new "<< tmp <<endl; 
      hh_corr_iphismtb ->SetBinContent(k,tmp);
      hh_corr_iphismtb ->SetBinError(k,tmpErr);
    }
    
    ///iphi of all SMs
    for(int k=0; k< 360; k++){
      int by = k; 
      if( k== 0) by = 360; 
	double tmp = corrfactorPhi[k] * pow( mean_allIphi/ hh_res_iphi[0]->GetBinContent(by),1);
	///txtout<<k<<" "<< tmp <<endl; 
	double tmpErr = tmp * hh_res_iphi[0]->GetBinError(by)/hh_res_iphi[0]->GetBinContent(by);
	hh_corr_iphi->SetBinContent(by,tmp);
	hh_corr_iphi->SetBinError(by,tmpErr);
    }

    if(hh_mpair_deadv1[0]->Integral()>10000){
      cout<<"check mass with dead crystals.."<<endl; 
      pi0_mfitpeak(hh_mpair_deadv1[0],pizFit,xlowFit,xhighFit,nPowFit,res,"","Dead0");
      double mpeakdead0 = res[0]; 
      cout<<"mpeak dead0 "<< mpeakdead0 <<endl; 
      for(int b=1; b<= hh_corr_deadv1->GetNbinsX(); b++){
	if( hh_mpair_deadv1[b-1]->Integral() > 1 ) {
	  if( (pizEta==1 &&hh_mpair_deadv1[b-1]->Integral() < 500 ) || (pizEta==2 &&hh_mpair_deadv1[b-1]->Integral() < 1000 ) ) {
	    cout<<"warning. too small statistic hh_mpair_deadv1 "<< b <<" "<< hh_mpair_deadv1[b-1]->Integral() <<endl;  
	    hh_corr_deadv1->SetBinContent(b,1);
	  }else{
	    pi0_mfitpeak(hh_mpair_deadv1[b-1],pizFit,xlowFit,xhighFit,nPowFit,res,"","Dead");
	    double tmp1 = corrfactorDead[b-1] * mpeakdead0/res[0]; 
	    hh_corr_deadv1->SetBinContent(b, tmp1); 
	    cout<<"mpeak dead"<<b<<" "<< res[0] <<" correction "<<tmp1<<endl;  
	  }
	}else{
	  hh_corr_deadv1->SetBinContent(b,1); 
	}
	if( pizEta==2 && b == 8 ){ ///eta no visible peak
	  hh_corr_deadv1->SetBinContent(b,1);
	}
	
      }
    }
    
  } ////for calibration step <=3  using Fit 
  
  
  if(stepc==4){ ///L3 method for last step 
    
    for(int j=0; j<170; j++){
      for(int k=0; k< 360; k++){
	double tmp = -1; 
	if( WTSUM[j][k] > 0){
	  tmp = corrfactoriEtaiPhi[j][k] * CORSUM[j][k]/ WTSUM[j][k];
	}
	if( j==0 && k<=1){
	  cout<<"j/k "<< j<<" "<<k<<" "<< corrfactoriEtaiPhi[j][k] <<" " << CORSUM[j][k] <<" " <<  WTSUM[j][k] <<" "<< tmp <<endl; 
	}
	txtout<<j<<" "<<k<<" "<<tmp <<endl; 
      }
    }
  }
  
  
  fnew ->Write();
  fnew ->Close();
  
  cout<<"end.." <<endl; 
  
}
