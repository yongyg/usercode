#include "rootheader.h"
#include "testSelection_cluster.h"

#include "common_functions.cc"
#include "foldEndcap.cc"

TRandom3 *rgen_;

#include "usefullcode.cc"

float interCalib_preCalib[170][360];
float interCalibEndcap_preCalib[2][101][101];

void copyConstant(float c1[170][360], float c2[170][360]){
  
  for(int j=0; j<170; j++){
    for(int k=0; k<360; k++){
      c2[j][k] = c1[j][k];
    }
  }
  
}

void isAtEcalBarrelModuleCracks(int ieta, int iphi, bool &etaCracks, bool &phiCracks){
  int absieta = abs(ieta);
  
  etaCracks = false; 
  phiCracks = false;
  
  if( absieta == 25 || absieta == 26 || absieta == 45 || absieta == 46 || absieta == 65 || absieta == 66) etaCracks = true; 
  if( iphi %20 == 1 || iphi %20 == 0) phiCracks = true; 
  
  
}

float ic[50][170][360]; ///at most 20 set

float icv1[50][170][360];

int ndeadflagietaiphi_ic[50][170][360];

float icwt_period[50][170][360]; //pi0&eta combined for each period



bool isTestBeamSM(int iSM){
    


  //if( iSM == 1 || iSM == 2  || iSM == 11 || iSM == 15  || iSM == 21 || iSM == 24) return true; 
  
  if( iSM == 1  || iSM == 11 || iSM == 15  || iSM == 24) return true; 
  
  
  else return false; 
    
}


#include "effSigma.C"

#include "getCrystaldeadflag.cc"

#include "getCalibConstants.cc"

#include "gausfit.cc"

float corrGR10_P_V10_over_GR09_R_V6A[170][360];
float corrGR_R_311_V1A_over_GR09_R_V6A[170][360];
float corrGR_R_311_V1A_over_GR10_P_V10[170][360];


float CphiCorrin[170][360];
float CBSCorrin[170][360];

float fitWind = 3; 
float sigmaTB = 0.55; 
//float sigmaTB = 0.44; 
//float sigmaTB = 0.95;

void fitHistogram(TH1F *h1, double res[]){
  
  float mean = h1->GetMean();
  float meanErr = h1->GetMeanError();
  float rmsEff = 100 * h1->GetRMS();
  float rmsEffErr = 100 * h1->GetRMSError();
  float rmsGaus = 100 * h1->GetRMS();
  float rmsGausErr = 100 * h1->GetRMSError();
  if( h1->Integral()>50){
    double resd[10];
    effSigma(h1,resd);
    rmsEff = resd[2]*100;
    float resf[20];
    
    fitgauswind2(h1,fitWind,fitWind,resf);
    //fitgauswindRefit(h1,resf);
    rmsGaus = resf[6]*100;
    rmsGausErr = resf[3] * 100;
    rmsEffErr = resf[3] * 100;
    mean = resf[0];
    meanErr = resf[1];
  }
  res[0] = mean; 
  res[1] = meanErr; 
  res[2] = rmsEff;
  res[3] = rmsEffErr; 
  res[4] = rmsGaus; 
  res[5] = rmsGausErr; 
}



float cpizv1[170][360];
float cpizv2[170][360];
float cpizv1in[170][360];
float cpizv2in[170][360];

float cetav1[170][360];
float cetav2[170][360];
float cetav1in[170][360];
float cetav2in[170][360];


float cpizoutput[170][360];
float cetaoutput[170][360];
float cpizetaoutput[170][360];

float cpizetaphibsoutput[170][360];
float cpizetaphibsoutputv1[170][360];

float cpizetaphibsprecaliboutput[170][360];
float cpizetaphibsprecaliboutputFinal[170][360];
float cpizetaphibsprecaliboutputv1[170][360];


float cphibsoutput[170][360];
float cphibsoutputv1[170][360];


float cpizoutputv1[170][360];
float cetaoutputv1[170][360];
float cpizetaoutputv1[170][360];


void rescaleConstTo(float Cin[170][360], float Corr[170][360]){
  
  float mean = 0; 
  int nmean = 0; 
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if( Cin[j][k]> 0){
	Cin[j][k] /= Corr[j][k];
	nmean ++; 
	mean += Cin[j][k];
      }
    }
  }
  mean /= nmean; 
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if( Cin[j][k]> 0){
	Cin[j][k] /= mean; 
      }
    }
  }
}

///Cprecaib * c1 == c_new * c2 
///c2 = c1 / ( c_new / cprecalib)

void rescaleConstTov1(float Cin[170][360], float Corr[170][360]){
  
  float mean = 0; 
  int nmean = 0; 
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if( Cin[j][k]> 0){
	Cin[j][k] *= Corr[j][k];
	nmean ++; 
	mean += Cin[j][k];
      }
    }
  }
  mean /= nmean; 
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if( Cin[j][k]> 0){
	Cin[j][k] /= mean; 
      }
    }
  }
}


///code to combined different pi0/eta EB IC for each laserTag


void combineCalibConstantv2(){
  
  

  readInterCalibConstEBSimple("interCalibEB_GR_P_V39.txt",interCalib_preCalib); 
  

  
  map<int, string> smScaleFiles; 


  
  
  //SM-scale files Pi0 
  
  smScaleFiles[0] = "calibres/deriveCalibConst.dflag2.pe1.step1.iter1.root";
  smScaleFiles[1] = "calibres/deriveCalibConst.dflag3.pe2.step1.iter1.root";
 
  // number of IC periods
  int nIC = 1; 
  


  map<int,string> icFiles; 
  
  //IC Pi0 
  icFiles[0] = "calibres/deriveCalibConst.dflag2.pe1.step4.iter30.txt";
  //IC Eta 
  icFiles[1] = "calibres/deriveCalibConst.dflag3.pe2.step4.iter30.txt";
  
  //crystal dead flag
  getCrystaldeadflagBarrel_v1("crystal_deadflag_eb_dflag2.txt",ndeadflagietaiphi_ic[0]); 
  getCrystaldeadflagBarrel_v1("crystal_deadflag_eb_dflag3.txt",ndeadflagietaiphi_ic[1]); 


    
  vector<string> inputfileStat; 
  inputfileStat.push_back("calibres/deriveCalibConst.dflag2.pe1.step3.iter11.root");
  inputfileStat.push_back("calibres/deriveCalibConst.dflag3.pe2.step3.iter11.root");



  ofstream txtoutTocheck("combinedCalibConstantv2.txt");
    
  
  int nConstantSet = int(icFiles.size());
  
  cout<<" nConstantSet " << nConstantSet <<endl; 
  
  
  map<int,TH1F*> hh_smscales; 
  float cc[170][360];


  int nMaxICSet = 50; 
  if( nConstantSet > nMaxICSet){
    cout<<"more than  " << nMaxICSet<<" "<< nConstantSet <<endl; 
    return; 
  }
  
  
  
  float icwt[170][360];
  float icwtv1[170][360];
  
  

  //   float icwt1[170][360];
  //   float icwt2[170][360];
  
  

  for(int n=0; n< nConstantSet; n++){
    string filename = smScaleFiles[n];
    TFile *ff = new TFile(filename.c_str(),"read");
    TH1F *hhtmp = (TH1F*)ff->Get("hh_corr_sm");
    hh_smscales[n] = hhtmp;
    
    filename = icFiles[n];

    readInterCalibConstEBSimplev1(filename.c_str(),cc);
  
    for(int j=0; j< 170; j++){
      for(int k=0; k< 360; k++){
	ic[n][j][k]  = cc[j][k];
      }
    }
    
    
    NormIetaAbsToUnitTestBeamSMOnly(ic[n]);
       
    if(n==0) cout<<"check ic[0][1] 1 " << ic[n][0][1]<<endl; 

    SetSMScale(ic[n], hh_smscales[n]);

    if(n==0) cout<<"check ic[0][1] 2 " << ic[n][0][1]<<endl; 
        
    NormCrystalDeadFlag_v1(ic[n],ndeadflagietaiphi_ic[n]);

    if(n==0) cout<<"check ic[0][1] 3 " << ic[n][0][1]<<endl; 
    
    
    copyConstant(ic[n],icv1[n]);
    
    /////for estimating precision
    NormSMScaleToUnit(icv1[n]);
  }

  //return; 
    

  TFile *fnew = new TFile("combineCalibConstantv2.root","recreate");
  
  
  TH1F *hh_csm[51][36][3]; //for each SM, all/centra/outer
  
  
  //|ieta|<=45 , all, removing eta phi boundaries, at eta moduaries, at phi boduaris 
  TH1F *hh_csmtb_ietaMod12[50][5];
  for(int n=0; n< nConstantSet; n++){
    for(int k=0; k<5; k++){
      TString histname = TString(Form("hh_csmtb_ietaMod12_%d_%d",n,k));
      hh_csmtb_ietaMod12[n][k] = new TH1F(histname,histname,500,0,2);
    }
  }
  

  
  TH1F *hh_c_ieta[50][170];
  TH1F *hh_c_ietaAbs[50][85];
  
  
  
    
  TH1F *hh_csmtb_ietaTT[50][34];
  TH1F *hh_csmtb_ietaTTAbs[50][17];
  TH1F *hh_csmco_ietaTTAbs[50][17];

  TH1F *hh_csmco_ietaTT[50][34];

  TH1F *hh_csmall_ietaTTAbs[50][17];
  
  
  
  TH1F *hh_wtavg_csmtb_ietaTT[34];
  TH1F *hh_wtavg_csmtb_ietaTTAbs[17];
  TH1F *hh_wtavg_csmco_ietaTT[34];
  TH1F *hh_wtavg_csmco_ietaTTAbs[17];

  TH1F *hh_wtavg_csmall_ietaTTAbs[17];

  
  TH1F *hh_res_csmtbietaTT[50][3];
  TH1F *hh_res_csmtbietaTTAbs[50][3];
  
  TH1F *hh_res_csmcoietaTT[50][3];
  TH1F *hh_res_csmcoietaTTAbs[50][3];
  
  TH1F *hh_res_csmallietaTTAbs[50][3];
  

  
  TH1F *hh_res_wtavg_csmtbietaTT[3]; //mean/rmsEff/rmsGaus

  TH1F *hh_res_wtavg_csmtbietaTTAbs[3]; //mean/rmsEff/rmsGaus
  
  TH1F *hh_res_wtavg_csmallietaTTAbs[3]; //mean/rmsEff/rmsGaus
  
  
  TH1F *hh_res_wtavg_csmcoietaTT[3]; //mean/rmsEff/rmsGaus
  TH1F *hh_res_wtavg_csmcoietaTTAbs[3]; //mean/rmsEff/rmsGaus
  
  TH1F *hh_diff_csmtb_ietaTTAbs[50][50][17];
  TH1F *hh_diff_csmtb_ietaTT[50][50][34];
  TH1F *hh_diff_csmco_ietaTTAbs[50][50][17];
  TH1F *hh_diff_csmco_ietaTT[50][50][34];
  TH1F *hh_diff_csmall_ietaTT[50][50][34];
  TH1F *hh_diff_csmall_ietaTTAbs[50][50][17];
  
  TH1F *hh_diff_csm[50][50][36][3]; //by each SM, all/central/outer
    
  
  TH1F *hh_res_csm[51][3][3];
  
  TH2F *hh2_c[51];
  TH2F *hh2_c_period[51];


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
  

  TH1F *hh_c_ietaAbs_period[50][85];
  TH1F *hh_c_ietaTTAbs_period[50][17];
  TH1F *hh_c_smtbietaTTAbs_period[50][17];
  
  
  for(int n=0; n< nIC+1; n++){
    for(int k=0; k< 85; k++){
      TString histname = TString(Form("hh_c_ietaAbs_period_%d_%d",n,k));
      hh_c_ietaAbs_period[n][k] = new TH1F(histname,histname,500,0,2);
    }
    for(int k=0; k< 17; k++){
      TString histname = TString(Form("hh_c_ietaTTAbs_period_%d_%d",n,k));
      hh_c_ietaTTAbs_period[n][k] = new TH1F(histname,histname,500,0,2);
      histname = TString(Form("hh_c_smtbietaTTAbs_period_%d_%d",n,k));
      hh_c_smtbietaTTAbs_period[n][k] = new TH1F(histname,histname,500,0,2);
    }
  }

    
  
  TH1F *hh_res_cietaAbs_period[50][3];
  TH1F *hh_res_cietaTTAbs_period[50][3];
  TH1F *hh_res_csmtbietaTTAbs_period[50][3];
  for(int n=0; n< nIC+1; n++){
    for(int k=0;k<3; k++){
      TString histname = TString(Form("hh_res_cietaAbs_period_%d_%d",n,k));
      hh_res_cietaAbs_period[n][k] = new TH1F(histname,histname,85,1,86);
      histname = TString(Form("hh_res_cietaTTAbs_period_%d_%d",n,k));
      hh_res_cietaTTAbs_period[n][k] = new TH1F(histname,histname,17,xbinLow1);
      histname = TString(Form("hh_res_csmtbietaTTAbs_period_%d_%d",n,k));
      hh_res_csmtbietaTTAbs_period[n][k] = new TH1F(histname,histname,17,xbinLow1);
    }
  }
  

  
  for(int n=0; n< nIC; n++){
    TString histname = TString(Form("hh2_c_period_%d",n));
    hh2_c_period[n] = new TH2F(histname,histname,171,-85,86,360,1,361);
   
  }

  
  for(int n=0; n< nConstantSet+1; n++){
    TString histname = TString(Form("hh2_c_%d",n));
    hh2_c[n] = new TH2F(histname,histname,171,-85,86,360,1,361);
   
  }
  
  for(int n=0; n< nConstantSet+1; n++){
    for(int j=0;j<36;j++){
      for(int k=0;k<3; k++){
	TString histname = TString(Form("hh_csm_%d_%d_%d",n,j,k));
	hh_csm[n][j][k] = new TH1F(histname,histname,500,0,2);
      }
    }
    
    for(int j=0;j<3; j++){
      for(int k=0;k<3; k++){
	TString histname = TString(Form("hh_res_csm_%d_%d_%d",n,j,k));
	hh_res_csm[n][j][k] = new TH1F(histname,histname,36,1,37);
      }
    }
  }
  


  for(int n=0; n< nConstantSet; n++){
    

    for(int j=0; j<17; j++){
      TString histname = TString(Form("hh_csmtb_ietaTTAbs_%d_%d",n,j));
      hh_csmtb_ietaTTAbs[n][j] = new TH1F(histname,histname,500,0,2);
      histname = TString(Form("hh_csmco_ietaTTAbs_%d_%d",n,j));
      hh_csmco_ietaTTAbs[n][j] = new TH1F(histname,histname,500,0,2);

      histname = TString(Form("hh_csmall_ietaTTAbs_%d_%d",n,j));
      hh_csmall_ietaTTAbs[n][j] = new TH1F(histname,histname,500,0,2);
      
    }
    
    
    for(int j=0; j<34; j++){
      TString histname = TString(Form("hh_csmtb_ietaTT_%d_%d",n,j));
      hh_csmtb_ietaTT[n][j] = new TH1F(histname,histname,500,0,2);
      histname = TString(Form("hh_csmco_ietaTT_%d_%d",n,j));
      hh_csmco_ietaTT[n][j] = new TH1F(histname,histname,500,0,2);
    }
  }
  
  for(int j=0; j<34; j++){
    TString histname = TString(Form("hh_wtavg_csmtb_ietaTT_%d",j));
    hh_wtavg_csmtb_ietaTT[j] = new TH1F(histname,histname,500,0,2);
    histname = TString(Form("hh_wtavg_csmco_ietaTT_%d",j));
    hh_wtavg_csmco_ietaTT[j] = new TH1F(histname,histname,500,0,2);
  }
  for(int j=0; j<17; j++){
    TString histname = TString(Form("hh_wtavg_csmtb_ietaTTAbs_%d",j));
    hh_wtavg_csmtb_ietaTTAbs[j] = new TH1F(histname,histname,500,0,2);
    histname = TString(Form("hh_wtavg_csmco_ietaTTAbs_%d",j));
    hh_wtavg_csmco_ietaTTAbs[j] = new TH1F(histname,histname,500,0,2);
    
    histname = TString(Form("hh_wtavg_csmall_ietaTTAbs_%d",j));
    hh_wtavg_csmall_ietaTTAbs[j] = new TH1F(histname,histname,500,0,2);
    
  }
  
 
  TH1F *hh_c_deadflag[50][20]; 
  
  TH1F *hh_c_deadflag_period[50][20];
  TH1F *hh_cc_deadflag_period[50][20];
  TH1F *hh_cc1_deadflag_period[50][20];
  
  for(int j=0; j< nConstantSet; j++){
    for(int k=0; k<20; k++){
      TString histname = TString(Form("hh_c_deadflag_%d_%d",j,k));
      hh_c_deadflag[j][k] = new TH1F(histname,histname,500,0,2);
      
      histname = TString(Form("hh_c_deadflag_period_%d_%d",j,k));
      hh_c_deadflag_period[j][k] = new TH1F(histname,histname,500,0,2);

      histname = TString(Form("hh_cc_deadflag_period_%d_%d",j,k));
      hh_cc_deadflag_period[j][k] = new TH1F(histname,histname,500,0,2);
      histname = TString(Form("hh_cc1_deadflag_period_%d_%d",j,k));
      hh_cc1_deadflag_period[j][k] = new TH1F(histname,histname,500,0,2);
      
    }
  }
  

  for(int j=0; j< nConstantSet; j++){
    for(int k=0; k<3; k++){

      TString histname = TString(Form("hh_res_csmtbietaTT_%d_%d",j,k));
      hh_res_csmtbietaTT[j][k] = new TH1F(histname,histname,34,xbinLow);
      histname = TString(Form("hh_res_csmtbietaTTAbs_%d_%d",j,k));
      hh_res_csmtbietaTTAbs[j][k] = new TH1F(histname,histname,17,xbinLow1);

      
      histname = TString(Form("hh_res_csmallietaTTAbs_%d_%d",j,k));
      hh_res_csmallietaTTAbs[j][k] = new TH1F(histname,histname,17,xbinLow1);
      
      
      histname = TString(Form("hh_res_csmcoietaTT_%d_%d",j,k));
      hh_res_csmcoietaTT[j][k] = new TH1F(histname,histname,34,xbinLow);
      histname = TString(Form("hh_res_csmcoietaTTAbs_%d_%d",j,k));
      hh_res_csmcoietaTTAbs[j][k] = new TH1F(histname,histname,17,xbinLow1);
      
    }
  }
  
  for(int k=0; k<3; k++){
    TString histname = TString(Form("hh_res_wtavg_csmtbietaTT_%d",k));
    hh_res_wtavg_csmtbietaTT[k] = new TH1F(histname,histname,34,xbinLow);
    histname = TString(Form("hh_res_wtavg_csmtbietaTTAbs_%d",k));
    hh_res_wtavg_csmtbietaTTAbs[k] = new TH1F(histname,histname,17,xbinLow1);
    
    histname = TString(Form("hh_res_wtavg_csmallietaTTAbs_%d",k));
    hh_res_wtavg_csmallietaTTAbs[k] = new TH1F(histname,histname,17,xbinLow1);

    histname = TString(Form("hh_res_wtavg_csmcoietaTT_%d",k));
    hh_res_wtavg_csmcoietaTT[k] = new TH1F(histname,histname,34,xbinLow);
    histname = TString(Form("hh_res_wtavg_csmcoietaTTAbs_%d",k));
    hh_res_wtavg_csmcoietaTTAbs[k] = new TH1F(histname,histname,17,xbinLow1);
  }
  TH2F *hh2_diff[50][50];
  
  
  TH1F *hh_res_diff_csmtbietaTTAbs[50][50][3];
  TH1F *hh_res_diff_csmcoietaTTAbs[50][50][3];
  TH1F *hh_res_diff_csmtbietaTT[50][50][3];
  TH1F *hh_res_diff_csmcoietaTT[50][50][3];
  TH1F *hh_res_diff_csmallietaTT[50][50][3];
  TH1F *hh_res_diff_csmallietaTTAbs[50][50][3];
  

  for(int n=0; n< nConstantSet; n++){
    for(int k=n+1; k< nConstantSet; k++){
      TString histname = TString (Form("hh2_diff_%dand%d",n,k));
      hh2_diff[n][k] = new TH2F(histname,histname,171,-85,86,360,1,361);
            
      for(int j=0; j<36; j++){
	for(int m=0; m<3; m++){
	  histname = TString ( Form("hh_diff_csm_%dand%d_%d_%d",n,k,j,m));
	  hh_diff_csm[n][k][j][m] = new TH1F(histname,histname,50,-0.1,0.1);
	}
      }
      
      for(int j=0; j<17; j++){
	histname = TString (Form("hh_diff_csmtb_ietaTTAbs_%dand%d_%d",n,k,j));
	hh_diff_csmtb_ietaTTAbs[n][k][j] = new TH1F(histname,histname,500,-0.1,0.1);
	histname = TString (Form("hh_diff_csmco_ietaTTAbs_%dand%d_%d",n,k,j));
	hh_diff_csmco_ietaTTAbs[n][k][j] = new TH1F(histname,histname,500,-0.1,0.1);
	histname = TString (Form("hh_diff_csmall_ietaTTAbs_%dand%d_%d",n,k,j));
	hh_diff_csmall_ietaTTAbs[n][k][j] = new TH1F(histname,histname,500,-0.1,0.1);

      }
      for(int j=0; j<34; j++){
	histname = TString (Form("hh_diff_csmtb_ietaTT_%dand%d_%d",n,k,j));
	hh_diff_csmtb_ietaTT[n][k][j] = new TH1F(histname,histname,500,-0.1,0.1);
	histname = TString (Form("hh_diff_csmco_ietaTT_%dand%d_%d",n,k,j));
	hh_diff_csmco_ietaTT[n][k][j] = new TH1F(histname,histname,500,-0.1,0.1);

	histname = TString (Form("hh_diff_csmall_ietaTT_%dand%d_%d",n,k,j));
	hh_diff_csmall_ietaTT[n][k][j] = new TH1F(histname,histname,500,-0.1,0.1);
	
      }
      for(int j=0;j<3; j++){
	histname = TString(Form("hh_res_diff_csmtbietaTTAbs_%dand%d_%d",n,k,j));
	hh_res_diff_csmtbietaTTAbs[n][k][j] = new TH1F(histname,histname,17,xbinLow1);
	histname = TString(Form("hh_res_diff_csmcoietaTTAbs_%dand%d_%d",n,k,j));
	hh_res_diff_csmcoietaTTAbs[n][k][j] = new TH1F(histname,histname,17,xbinLow1);
	histname = TString(Form("hh_res_diff_csmtbietaTT_%dand%d_%d",n,k,j));
	hh_res_diff_csmtbietaTT[n][k][j] = new TH1F(histname,histname,34,xbinLow);
	histname = TString(Form("hh_res_diff_csmcoietaTT_%dand%d_%d",n,k,j));
	hh_res_diff_csmcoietaTT[n][k][j] = new TH1F(histname,histname,34,xbinLow);

	histname = TString(Form("hh_res_diff_csmallietaTT_%dand%d_%d",n,k,j));
	hh_res_diff_csmallietaTT[n][k][j] = new TH1F(histname,histname,34,xbinLow);
	histname = TString(Form("hh_res_diff_csmallietaTTAbs_%dand%d_%d",n,k,j));
	hh_res_diff_csmallietaTTAbs[n][k][j] = new TH1F(histname,histname,17,xbinLow1);
      }
      
    }
  }


  TH2F *hh2_largeICdiff[50][50];
  for(int n=0; n< nIC; n++){
    for(int k=n+1; k< nIC; k++){
      TString histname = TString (Form("hh2_largeICdiff_%dand%d",n,k));
      hh2_largeICdiff[n][k] =new TH2F(histname,histname,171,-85,86,360,1,361);
    }
  }
  TH2F *hh2_largeICdiff_all =new TH2F("hh2_largeICdiff_all","hh2_largeICdiff_all",171,-85,86,360,1,361);

  
  
  TH1F *hh_res_ieta[50][4];////[4] means the stat error.
  
  for(int n=0; n< int(inputfileStat.size()); n++){
    for(int k=0;k<4;k++){
      string filename = string(Form("hh_res_ieta_%d_%d",n,k));
      hh_res_ieta[n][k] =new TH1F(filename.c_str(),filename.c_str(),171,-85,86);
    }
  }
  TH1F *hh_statErr_ietaAbs[50];
  for(int n=0; n< int(inputfileStat.size()); n++){
    string filename = string(Form("hh_statErr_ietaAbs_%d",n));
    hh_statErr_ietaAbs[n] =new TH1F(filename.c_str(),filename.c_str(),85,1,86);
  }
  TH1F *hh_statErr_ietaTTAbs[50];
  for(int n=0; n< int(inputfileStat.size()); n++){
    string filename = string(Form("hh_statErr_ietaTTAbs_%d",n));
    hh_statErr_ietaTTAbs[n] =new TH1F(filename.c_str(),filename.c_str(),17,xbinLow1);
  }

  TH1F *hh_statErr_ietaAbs_period[50]; //pi0&eta combined for each period
  for(int j=0; j< nIC; j++){
    string filename = string(Form("hh_statErr_ietaAbs_period_%d",j));
    hh_statErr_ietaAbs_period[j]=new TH1F(filename.c_str(),filename.c_str(),85,1,86);
  }
  
  
  ///using MC-based forumula + 0.5/1 % sys.
  // 7.4/resolution * 17 /sqrt(N) * sqrt( 1+ 1.8/sob) + 0.5/1% 
  
 for(int j=0; j< int(inputfileStat.size());j++){
    string filename = inputfileStat[j];
    TFile *f1 = new TFile(filename.c_str(),"read");
    for(int k=0;k<4;k++){
      string histname = string (Form("hh_res_ieta_%d",k));
      TH1F *hhtmp = (TH1F*)f1->Get(histname.c_str());
      if(hhtmp==0){
	cout<<"empty hh_res_ieta_ ! "<<endl; 
	return; 
      }
      
      hh_res_ieta[j][k]->Add(hhtmp);
    }

    if( f1->Get("hh_res_ietaTTAbs_0") !=0){
      TH1F *htmp1 = (TH1F*)f1->Get("hh_res_ietaTTAbs_0");
      TH1F *htmp2 = (TH1F*)f1->Get("hh_res_ietaTTAbs_1");
      TH1F *htmp3 = (TH1F*)f1->Get("hh_res_ietaTTAbs_2");
      TH1F *htmp4 = (TH1F*)f1->Get("hh_res_ietaTTAbs_3");
      for(int b=1; b<= 17; b++){
	float m0 = htmp1->GetBinContent(b);
	float sig = htmp3->GetBinContent(b);
	float sb = htmp4->GetBinContent(b);
	float n = htmp2->GetBinContent(b)/ (360*5*2*2);
	float statErr = (sig/m0*100)/7.4 * 17/sqrt(n) * sqrt( 1+ 1.8/sb); 
	hh_statErr_ietaTTAbs[j]->SetBinContent(b,statErr);
      }
    }
    
  }
 
 
  
  
  for(int j=0; j< int(inputfileStat.size());j++){

    for(int n=1; n<=85 ; n++){
      float npiz = 0.5*(hh_res_ieta[j][1]->GetBinContent(87+n-1) + hh_res_ieta[j][1]->GetBinContent(85-(n-1)) ) /(360*2) ; 
      float sob = 0.5*(hh_res_ieta[j][3]->GetBinContent(87+n-1) + hh_res_ieta[j][3]->GetBinContent(85-(n-1)) );
      float reso = 0.5*(hh_res_ieta[j][2]->GetBinContent(87+n-1) + hh_res_ieta[j][2]->GetBinContent(85-(n-1)) )/hh_res_ieta[j][0]->GetBinContent(87+n-1);

      if(npiz<=0){
	cout<<"empty histogram!!" <<endl; 
	return; 
      }
      
      float statErr = (reso*100)/7.4 * 17/sqrt(npiz) * sqrt( 1+ 1.8/sob); 

      cout<<"statErr " << reso <<" "<< npiz <<" "<< sob <<" "<< statErr <<endl; 

      hh_statErr_ietaAbs[j]->SetBinContent(n,statErr);
    }
  }
  
  //return; 

  
  ///the combined IC from all ICs
  ofstream txtout("interCalibConstants.combinedPi0EtaAllPeriod.EcalBarrel.txt",ios::out);
  

  ofstream txtout_period[50]; //pi0eta combined for each period
  for(int j=0; j< nIC; j++){
    string filename = string(Form("interCalibConstants.combinedPi0EtaPeriod%d.EcalBarrel.txt",j));
    txtout_period[j].open(filename.c_str(),ios::out);
  }
  
  ofstream txtout_periodv1[50]; //pi0eta combined for each period / averaged all 
  for(int j=0; j< nIC; j++){
    string filename = string(Form("interCalibConstants.combinedPi0EtaPeriod%d.EcalBarrel.txtv1",j));
    txtout_periodv1[j].open(filename.c_str(),ios::out);
  }
  

  cout<<"fill" <<endl; 

  
  

  for(int n=0; n< nConstantSet; n++){
    
    for(int j=0; j<170; j++){
      
      int ieta = j-85; 
      if( ieta >=0) ieta += 1; 
      
      for(int k=0; k<360; k++){
	int iphi = k; 
	if( k==0) iphi = 360; 
	int ietaTTAbs = (abs(ieta)-1)/5; 
	int iSM = (iphi-1)/20+1;
	if( ieta<0) iSM += 18;
	int smTB = isTestBeamSM(iSM);
	int beta = ieta+85+1; 
	int bphi = iphi; 
	float c = icv1[n][j][k]; ///sm normalized to unit
	if( c>0){
	  
	  int ndflag = -1; 
	  ndflag = ndeadflagietaiphi_ic[n][j][k];
	  
	  if( ndflag <0){
	    cout<<"wrong dflag? " << ndflag <<" "<<n<<" "<<j<<" "<<k<<endl; 
	    return; 
	  }
	  hh_c_deadflag[n][ndflag]->Fill(c);
	  	  
	  
	  hh2_c[n]->SetBinContent(beta,bphi,ic[n][j][k]);
	  
	  int ietaflag[3] = {abs(ieta)<25,abs(ieta)>=25 && abs(ieta)<60,abs(ieta)>=60};
	  for(int jj=0;jj<3; jj++){
	    if( ietaflag[jj]==0) continue; 
	    hh_csm[n][iSM-1][jj]->Fill(c);
	  }
	  
	  if( smTB){
	    hh_csmtb_ietaTTAbs[n][ietaTTAbs]->Fill(c);
	    hh_csmtb_ietaTT[n][j/5]->Fill(c);
	    
	    bool etaCracks; 
	    bool phiCracks; 
	    isAtEcalBarrelModuleCracks(ieta,iphi,etaCracks,phiCracks);
	    int crackFlag[5] = {1,!etaCracks && !phiCracks, etaCracks, phiCracks, etaCracks || phiCracks};
	    for(int t =0; t<5; t++){
	      if(crackFlag[t]==0) continue; 
	      hh_csmtb_ietaMod12[n][t]->Fill(c);
	    }
	    
	  }else{
	    hh_csmco_ietaTTAbs[n][ietaTTAbs]->Fill(c);
	    hh_csmco_ietaTT[n][j/5]->Fill(c);
	  }
	  
	  hh_csmall_ietaTTAbs[n][ietaTTAbs]->Fill(c);

	}else{
	  hh2_c[n]->SetBinContent(beta,bphi,-1);
	}
	
      }

    }
  }
  
  
  for(int j=0; j<170; j++){
      
    int ieta = j-85; 
    if( ieta >=0) ieta += 1; 
      
    for(int k=0; k<360; k++){
      int iphi = k; 
      if( k==0) iphi = 360; 
      int ietaTTAbs = (abs(ieta)-1)/5; 
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      int smTB = isTestBeamSM(iSM);
      int beta = ieta+85+1; 
      int bphi = iphi; 
      
      for(int n1 =0; n1 < nConstantSet; n1++){
	for(int n2 =n1+1; n2 < nConstantSet; n2++){
	  if(ic[n1][j][k] >0 && ic[n2][j][k] > 0){

	    float diff = ic[n1][j][k] - ic[n2][j][k]; 
	    float average = 0.5*(ic[n1][j][k] + ic[n2][j][k]);
	    diff /= average; 
	    	    
	    hh2_diff[n1][n2] ->SetBinContent(beta,bphi, diff);
	    
	    if(smTB){
	      hh_diff_csmtb_ietaTTAbs[n1][n2][ietaTTAbs]->Fill(diff);
	      hh_diff_csmtb_ietaTT[n1][n2][j/5]->Fill(diff);
	    }else{
	      hh_diff_csmco_ietaTTAbs[n1][n2][ietaTTAbs]->Fill(diff);
	      hh_diff_csmco_ietaTT[n1][n2][j/5]->Fill(diff);
	    }
	    hh_diff_csmall_ietaTT[n1][n2][j/5]->Fill(diff);
	    hh_diff_csmall_ietaTTAbs[n1][n2][ietaTTAbs]->Fill(diff);
	    	    
	  }else{
	    hh2_diff[n1][n2] ->SetBinContent(beta,bphi,-1);
	  }
	}
      }
    }
  }
    
  
   
  cout<<"fitting "<<endl; 
  double resfit[10];
  
  
  
  for(int n1 =0; n1 < nConstantSet; n1++){
     for(int n2 =n1+1; n2 < nConstantSet; n2++){
       
       for(int j=0;j<34; j++){
	 fitHistogram(hh_diff_csmtb_ietaTT[n1][n2][j],resfit);
	 for(int n=0;n<3; n++){
	   hh_res_diff_csmtbietaTT[n1][n2][n]->SetBinContent(j+1,resfit[2*n]);
	   hh_res_diff_csmtbietaTT[n1][n2][n]->SetBinError(j+1,resfit[2*n+1]);
	 }
	 fitHistogram(hh_diff_csmco_ietaTT[n1][n2][j],resfit);
	 for(int n=0;n<3; n++){
	   hh_res_diff_csmcoietaTT[n1][n2][n]->SetBinContent(j+1,resfit[2*n]);
	   hh_res_diff_csmcoietaTT[n1][n2][n]->SetBinError(j+1,resfit[2*n+1]);
	 }
	 fitHistogram(hh_diff_csmall_ietaTT[n1][n2][j],resfit);
	 for(int n=0;n<3; n++){
	   hh_res_diff_csmallietaTT[n1][n2][n]->SetBinContent(j+1,resfit[2*n]);
	   hh_res_diff_csmallietaTT[n1][n2][n]->SetBinError(j+1,resfit[2*n+1]);
	 }

       }
       for(int j=0;j<17; j++){
	 fitHistogram(hh_diff_csmtb_ietaTTAbs[n1][n2][j],resfit);
	 for(int n=0;n<3; n++){
	   hh_res_diff_csmtbietaTTAbs[n1][n2][n]->SetBinContent(j+1,resfit[2*n]);
	   hh_res_diff_csmtbietaTTAbs[n1][n2][n]->SetBinError(j+1,resfit[2*n+1]);
	 }
	 fitHistogram(hh_diff_csmco_ietaTTAbs[n1][n2][j],resfit);
	 for(int n=0;n<3; n++){
	   hh_res_diff_csmcoietaTTAbs[n1][n2][n]->SetBinContent(j+1,resfit[2*n]);
	   hh_res_diff_csmcoietaTTAbs[n1][n2][n]->SetBinError(j+1,resfit[2*n+1]);
	 }
	 fitHistogram(hh_diff_csmall_ietaTTAbs[n1][n2][j],resfit);
	 for(int n=0;n<3; n++){
	   hh_res_diff_csmallietaTTAbs[n1][n2][n]->SetBinContent(j+1,resfit[2*n]);
	   hh_res_diff_csmallietaTTAbs[n1][n2][n]->SetBinError(j+1,resfit[2*n+1]);
	 }
	 
       }
       
     }
  }
  
  
  for(int k=0;k<nConstantSet; k++){
    
    
    for(int j=0;j<34; j++){
      fitHistogram(hh_csmtb_ietaTT[k][j],resfit);
      for(int n=0;n<3; n++){
	hh_res_csmtbietaTT[k][n]->SetBinContent(j+1,resfit[2*n]);
	hh_res_csmtbietaTT[k][n]->SetBinError(j+1,resfit[2*n+1]);
      }
    }
    for(int j=0;j<17; j++){
      fitHistogram(hh_csmtb_ietaTTAbs[k][j],resfit);
      for(int n=0;n<3; n++){
	hh_res_csmtbietaTTAbs[k][n]->SetBinContent(j+1,resfit[2*n]);
	hh_res_csmtbietaTTAbs[k][n]->SetBinError(j+1,resfit[2*n+1]);
      }
    }
    
    for(int j=0;j<17; j++){
      fitHistogram(hh_csmall_ietaTTAbs[k][j],resfit);
      for(int n=0;n<3; n++){
	hh_res_csmallietaTTAbs[k][n]->SetBinContent(j+1,resfit[2*n]);
	hh_res_csmallietaTTAbs[k][n]->SetBinError(j+1,resfit[2*n+1]);
      }
    }
    
    

    for(int j=0;j<34; j++){
      fitHistogram(hh_csmco_ietaTT[k][j],resfit);
      for(int n=0;n<3; n++){
	hh_res_csmcoietaTT[k][n]->SetBinContent(j+1,resfit[2*n]);
	hh_res_csmcoietaTT[k][n]->SetBinError(j+1,resfit[2*n+1]);
      }
    }
    for(int j=0;j<17; j++){
      fitHistogram(hh_csmco_ietaTTAbs[k][j],resfit);
      for(int n=0;n<3; n++){
	hh_res_csmcoietaTTAbs[k][n]->SetBinContent(j+1,resfit[2*n]);
	hh_res_csmcoietaTTAbs[k][n]->SetBinError(j+1,resfit[2*n+1]);
      }
    }
    
  }
  
  cout<<"combining " <<endl; 
  
  
  float wtSumC_period[50] = {0};
  float wtSumS_period[50] = {0};
  


  for(int j=0; j<170; j++){
    int ieta = j-85; 
    if( ieta >=0) ieta += 1; 
    for(int k=0; k<360; k++){
      int iphi = k; 
      if( k==0) iphi = 360; 
      int ietaTTAbs = (abs(ieta)-1)/5; 
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      int smTB = isTestBeamSM(iSM);
      int beta = ieta+85+1; 
      int bphi = iphi; 
      

      float wtSumC = 0; 
      float wtSumS = 0; 

      for(int n=0; n< nIC; n++){ ///for each period
	wtSumC_period[n] = 0; 
	wtSumS_period[n] = 0; 
      }
      
      
      for(int n=0; n< nConstantSet; n++){
	float c = ic[n][j][k];


	
	
	//float sigma = hh_res_csmtbietaTTAbs[n][2]->GetBinContent(ietaTTAbs+1);
	//if( sigma > sigmaTB){
	//sigma = sqrt( sigma * sigma - sigmaTB * sigmaTB);
	//}
	float statErr = 0; 
	float sysErr = 0.5; 
	if(abs(ieta)>=60) sysErr = 1; 
	
	// 	if(n>=0 && n<=1){
	// 	  statErr = hh_res_diff_csmallietaTTAbs[0][1][2]->GetBinContent(ietaTTAbs+1)/sqrt(2);
	// 	}else if( n>=2 && n<=3){
	// 	  statErr = hh_res_diff_csmallietaTTAbs[2][3][2]->GetBinContent(ietaTTAbs+1)/sqrt(2);
	// 	}
	

	//now use MC-predicted precision
	statErr = hh_statErr_ietaAbs[n]->GetBinContent(abs(ieta));
	
	float sigma = sqrt( statErr * statErr + sysErr * sysErr);
	
	//if( ieta==1) cout<<"sigma " << ieta<<" n"<< n<<" "<< sigma <<" "<< statErr <<" "<< sysErr <<endl; 
	
	if( c > 0){
	    
	  float tmp1 = c/ ( sigma * sigma);
	  float tmp2 = 1/(sigma * sigma);
	  int nperiod = n% nIC; 
	  
	  wtSumC_period[nperiod] += tmp1; 
	  wtSumS_period[nperiod] += tmp2; 
	  
	  if( c> 0){ //all combined
	    wtSumC += tmp1; 
	    wtSumS += tmp2; 
	  }
	}
      }
      
      
      if( wtSumC > 0){  //combined all
	icwt[j][k] = wtSumC / wtSumS; 
	hh2_c[nConstantSet]->SetBinContent(beta,bphi,icwt[j][k]);
      }else{
	icwt[j][k] = -1;
	hh2_c[nConstantSet]->SetBinContent(beta,bphi,-1);
      }
      
      for(int n=0; n< nIC; n++){ ///for each period
	if( wtSumC_period[n] > 0){
	  icwt_period[n][j][k] = wtSumC_period[n] / wtSumS_period[n]; 
	  hh2_c_period[n]->SetBinContent(beta,bphi,icwt_period[n][j][k]);
	}else{
	  icwt_period[n][j][k] = -1;
	  hh2_c_period[n]->SetBinContent(beta,bphi,-1);
	}
      }
      
    }
  }
  
  

  float statErr_allCombined[85] = {0} ; 
  
  for(int b=1; b<= hh_statErr_ietaAbs[0]->GetNbinsX(); b++){
    float sumStatErr2 = 0; 
    for(int n=0; n< nConstantSet; n++){
      float statErr = hh_statErr_ietaAbs[n]->GetBinContent(b);
      sumStatErr2 += 1./ ( statErr * statErr );
    }
    float statErr = sqrt( 1./ sumStatErr2 ); 
    statErr_allCombined[b-1] = statErr; 
    
    cout<<" statErr_allCombined " << b <<" "<< statErr_allCombined[b-1] << " % "<<endl; 

  }
  

  ///Stat error for each period
  for(int j=0; j< nIC; j++){
    for(int b=1; b<= hh_statErr_ietaAbs[0]->GetNbinsX();b++){
      float statErrPi0 =  hh_statErr_ietaAbs[j]->GetBinContent(b);
      float statErrEta =  hh_statErr_ietaAbs[j+nIC]->GetBinContent(b);
      float statErr = 1/(statErrPi0*statErrPi0) + 1./(statErrEta*statErrEta);
      statErr = sqrt(1/statErr);
      hh_statErr_ietaAbs_period[j]->SetBinContent(b,statErr);
    }
  }
  

  txtoutTocheck <<"checkme " <<endl; 
  
  for(int j=0; j<170; j++){
    
    int ieta = j-85; 
    if( ieta >=0) ieta += 1; 
    
   
    int absieta = abs(ieta);
    
    for(int k=0; k<360; k++){
      int iphi = k; 
      if( k==0) iphi = 360; 
      int beta = ieta+85+1; 
      int bphi = iphi; 
      int ietaTTAbs = (abs(ieta)-1)/5; 
      
      bool largeIC = false; 
      
      for(int n=0; n< nIC; n++){ ///for each period
	if( icwt_period[n][j][k] > 0 && (icwt_period[n][j][k] >1.2 || icwt_period[n][j][k] < 0.8)){
	  largeIC = true; 
	}
      }
      if(largeIC){
	txtoutTocheck<<"largeIC "<< j<<" "<<k <<" "<<ieta <<" "<<iphi<<" " ; 
	for(int n1 =0; n1 < nIC ; n1++){
	  txtoutTocheck<<icwt_period[n1][j][k]<<" "; 
	}
	txtoutTocheck<<endl; 
      }
      
      bool checkICdiff = false; 
      for(int n1=0; n1< nIC; n1++){ ///for each period
	for(int n2=n1+1;n2< nIC; n2++){ ///for each period
	  if( icwt_period[n1][j][k]>0 && icwt_period[n2][j][k]>0){
	    float reldiff = fabs( icwt_period[n1][j][k] - icwt_period[n2][j][k] )/ ( 0.5* (icwt_period[n1][j][k] + icwt_period[n2][j][k])); 

	    float statErr1 = hh_statErr_ietaAbs_period[n1]->GetBinContent(absieta)/100;
	    float statErr2 = hh_statErr_ietaAbs_period[n2]->GetBinContent(absieta)/100;
	    
	    ////assuming same sys.
	    float diff_statErr = sqrt( statErr1*statErr1 + statErr2*statErr2);
	    
	    checkICdiff = false; 
	   //  if( abs(ieta)<60 && reldiff >0.02){
// 	      checkICdiff = true; 
// 	    }
// 	    if( abs(ieta)>=60 && reldiff >0.1){
// 	      checkICdiff = true; 
//  	    }
	    checkICdiff = reldiff > 3 * diff_statErr ; 
	    if(checkICdiff){
	      hh2_largeICdiff[n1][n2]->SetBinContent(beta,bphi,1);
	    }else{
	      hh2_largeICdiff[n1][n2]->SetBinContent(beta,bphi,-1);
	    }
	    
	  }else{
	    hh2_largeICdiff[n1][n2]->SetBinContent(beta,bphi,-1);
	  }
	}
      }
      if(checkICdiff){
	hh2_largeICdiff_all->SetBinContent(beta,bphi,1);
	txtoutTocheck<<"checkICdiff "<< j<<" "<<k <<" "<<ieta <<" "<<iphi<<" " ; 
	for(int n1 =0; n1 < nIC ; n1++){
	  txtoutTocheck<<icwt_period[n1][j][k]<<" "; 
	}
	txtoutTocheck<<endl; 
      }else{
	hh2_largeICdiff_all->SetBinContent(beta,bphi,-1);
      }
      
    }
  }
  


  
  scaleMeanToUnit(icwt);
  ///for estimating precision of combined IC
  copyConstant(icwt,icwtv1);
  NormSMScaleToUnit(icwtv1);
    
  
  for(int n=0; n< nIC; n++){ ///for each period
    scaleMeanToUnit(icwt_period[n]);
  }
    
  
  for(int j=0; j<170; j++){
    int ieta = j-85; 
    if( ieta >=0) ieta += 1; 
    for(int k=0; k<360; k++){
      int iphi = k; 
      if( k==0) iphi = 360; 
      int ietaTTAbs = (abs(ieta)-1)/5; 
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      int smTB = isTestBeamSM(iSM);
      
      
      float sysErr = 0.5; 
      if(abs(ieta)>=60) sysErr = 1; 
      
      float icErr = sysErr; 
      float statErr = statErr_allCombined[abs(ieta)-1];
      icErr = sqrt(sysErr *sysErr + statErr * statErr);
      icErr /=100;
      
      
      if( icwt[j][k] > 0){
	txtout<<ieta<<" "<<iphi<<" "<< icwt[j][k]*interCalib_preCalib[j][k]<<" "<< icErr* icwt[j][k]*interCalib_preCalib[j][k]<<endl; 
	//txtout<<ieta<<" "<<iphi<<" "<< icwt[j][k]*interCalib_preCalib[j][k]<<endl; 
      }else{	
	txtout<<ieta<<" "<<iphi<<" "<< -1 <<" "<< 999 <<endl; 
	//txtout<<ieta<<" "<<iphi<<" "<< interCalib_preCalib[j][k]<<endl; 

      }
      
      float c = icwtv1[j][k]; ///sm normalized to unit
      
      
      int absieta = abs(ieta);
      
      for(int n=0; n< nIC+1; n++){ ///for each period and all combined

	if(n< nIC){
	  if( icwt_period[n][j][k] > 0){
	    hh_c_ietaAbs_period[n][absieta-1]->Fill(icwt_period[n][j][k]);
	    hh_c_ietaTTAbs_period[n][ietaTTAbs]->Fill(icwt_period[n][j][k]);

	    if( smTB){
	      hh_c_smtbietaTTAbs_period[n][ietaTTAbs]->Fill(icwt_period[n][j][k]);
	    }
	    int deadflag = ndeadflagietaiphi_ic[n][j][k]; 
	    if( deadflag<0){
	      cout<<"wrong deadlfag !!! n " << n <<" "<<endl; 
	      exit(1);
	    }
	    hh_c_deadflag_period[n][deadflag]->Fill(icwt_period[n][j][k]);

	    if(deadflag>0){
	      hh_c_deadflag_period[n][19]->Fill(icwt_period[n][j][k]);
	      hh_cc_deadflag_period[n][19]->Fill( interCalib_preCalib[j][k] * icwt_period[n][j][k]); 
	    }
	    
	  }
	}else{
	  if( icwt[j][k] > 0){
	    hh_c_ietaAbs_period[n][absieta-1]->Fill(icwt[j][k]);
	    hh_c_ietaTTAbs_period[n][ietaTTAbs]->Fill(icwt[j][k]);

	    int deadflag = ndeadflagietaiphi_ic[0][j][k]; 
	    if( deadflag<0){
	      cout<<"wrong deadlfag !!! comb " << n <<" "<<endl; 
	      exit(1);
	    }
	    hh_c_deadflag_period[n][deadflag]->Fill(icwt[j][k]);
	    if(deadflag>0){
	      hh_c_deadflag_period[n][19]->Fill(icwt[j][k]);
	    }
	    
	    if( smTB){
	      hh_c_smtbietaTTAbs_period[n][ietaTTAbs]->Fill(icwt[j][k]);
	    }

	  }
	}
      }
            

      for(int n=0; n< nIC; n++){ ///for each period
	statErr = hh_statErr_ietaAbs_period[n]->GetBinContent(absieta); 
	icErr = sqrt( statErr * statErr + sysErr * sysErr);
	icErr /= 100;
	
	if( icwt_period[n][j][k] > 0){
	  txtout_period[n]<<ieta<<" "<<iphi<<" "<< icwt_period[n][j][k]*interCalib_preCalib[j][k]<<" "<< icErr*icwt_period[n][j][k]*interCalib_preCalib[j][k] <<endl; 
	  //txtout_period[n]<<ieta<<" "<<iphi<<" "<< icwt_period[n][j][k]*interCalib_preCalib[j][k]<<endl; 
	}else{	
	  txtout_period[n]<<ieta<<" "<<iphi<<" "<< -1<<" "<< 999 <<endl; 
	  //txtout_period[n]<<ieta<<" "<<iphi<<" "<< interCalib_preCalib[j][k]<<" "<< 999 <<endl; 
	  
	}
      }
      
      
     
      if( c>0){

	int ietaflag[3] = {abs(ieta)<25,abs(ieta)>=25 && abs(ieta)<60,abs(ieta)>=60};
	for(int jj=0;jj<3; jj++){
	  if( ietaflag[jj]==0) continue; 
	  hh_csm[nConstantSet][iSM-1][jj]->Fill(c);
	}
	
	if( smTB){
	  hh_wtavg_csmtb_ietaTT[j/5]->Fill(c);
	  hh_wtavg_csmtb_ietaTTAbs[ietaTTAbs]->Fill(c);
	}else{
	  hh_wtavg_csmco_ietaTT[j/5]->Fill(c);
	  hh_wtavg_csmco_ietaTTAbs[ietaTTAbs]->Fill(c);
	}

	///all sm
	hh_wtavg_csmall_ietaTTAbs[ietaTTAbs]->Fill(c);
	
      }
      
    }
  }
  
  cout <<" fitting combined  " <<endl; 
  
  
  for(int n=0; n< nConstantSet+1; n++){
    for(int j=0; j<36; j++){
      
      for(int k=0;k<3;k++){
	fitHistogram(hh_csm[n][j][k],resfit);
	for(int jj=0;jj<3; jj++){
	  hh_res_csm[n][k][jj]->SetBinContent(j+1,resfit[2*jj]);
	  hh_res_csm[n][k][jj]->SetBinError(j+1,resfit[2*jj+1]);
	}
      }
    }
  }
    

  
  for(int n1=0; n1< nIC+1; n1++){ ///for each period and all combined
    for(int j=0; j< 85; j++){
      fitHistogram(hh_c_ietaAbs_period[n1][j],resfit);
      for(int n=0;n<3; n++){
	hh_res_cietaAbs_period[n1][n]->SetBinContent(j+1,resfit[2*n]);
	hh_res_cietaAbs_period[n1][n]->SetBinError(j+1,resfit[2*n+1]);
      }
    }
    for(int j=0; j< 17; j++){
      fitHistogram(hh_c_ietaTTAbs_period[n1][j],resfit);
      for(int n=0;n<3; n++){
	hh_res_cietaTTAbs_period[n1][n]->SetBinContent(j+1,resfit[2*n]);
	hh_res_cietaTTAbs_period[n1][n]->SetBinError(j+1,resfit[2*n+1]);
      }
      fitHistogram(hh_c_smtbietaTTAbs_period[n1][j],resfit);
      for(int n=0;n<3; n++){
	hh_res_csmtbietaTTAbs_period[n1][n]->SetBinContent(j+1,resfit[2*n]);
	hh_res_csmtbietaTTAbs_period[n1][n]->SetBinError(j+1,resfit[2*n+1]);
      }
    }
  }
  

  for(int j=0;j<34; j++){
    fitHistogram(hh_wtavg_csmtb_ietaTT[j],resfit);
    for(int n=0;n<3; n++){
      hh_res_wtavg_csmtbietaTT[n]->SetBinContent(j+1,resfit[2*n]);
      hh_res_wtavg_csmtbietaTT[n]->SetBinError(j+1,resfit[2*n+1]);
    }
    fitHistogram(hh_wtavg_csmco_ietaTT[j],resfit);
    for(int n=0;n<3; n++){
      hh_res_wtavg_csmcoietaTT[n]->SetBinContent(j+1,resfit[2*n]);
      hh_res_wtavg_csmcoietaTT[n]->SetBinError(j+1,resfit[2*n+1]);
    }
  }
  for(int j=0;j<17; j++){
    fitHistogram(hh_wtavg_csmtb_ietaTTAbs[j],resfit);
    for(int n=0;n<3; n++){
      hh_res_wtavg_csmtbietaTTAbs[n]->SetBinContent(j+1,resfit[2*n]);
      hh_res_wtavg_csmtbietaTTAbs[n]->SetBinError(j+1,resfit[2*n+1]);
    }


    fitHistogram(hh_wtavg_csmall_ietaTTAbs[j],resfit);
    for(int n=0;n<3; n++){
      hh_res_wtavg_csmallietaTTAbs[n]->SetBinContent(j+1,resfit[2*n]);
      hh_res_wtavg_csmallietaTTAbs[n]->SetBinError(j+1,resfit[2*n+1]);
    }

    fitHistogram(hh_wtavg_csmco_ietaTTAbs[j],resfit);
    for(int n=0;n<3; n++){
      hh_res_wtavg_csmcoietaTTAbs[n]->SetBinContent(j+1,resfit[2*n]);
      hh_res_wtavg_csmcoietaTTAbs[n]->SetBinError(j+1,resfit[2*n+1]);
    }
  }
  
  
  fnew->Write();
  fnew->Close();
  
}
