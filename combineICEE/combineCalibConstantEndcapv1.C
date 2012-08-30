#include "rootheader.h"
#include "testSelection_cluster.h"

#include "common_functions.cc"
#include "foldEndcap.cc"

TRandom3 *rgen_;

#include "usefullcode.cc"



float ccalibpretag[2][101][101];

float interCalib_preCalib[170][360];
float interCalibEndcap_preCalib[2][101][101];



bool isTestBeamSM(int iSM){
    


  //if( iSM == 1 || iSM == 2  || iSM == 11 || iSM == 15  || iSM == 21 || iSM == 24) return true; 
  
  if( iSM == 1  || iSM == 11 || iSM == 15  || iSM == 24) return true; 
  
  
  else return false; 
    
}


#include "effSigma.C"

#include "getCrystaldeadflag.cc"

#include "getCalibConstants.cc"

#include "gausfit.cc"


float fitWind = 3; 

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


float ic[50][2][101][101]; //at most 10 sets

int ndeadflag_ic[50][2][101][101];

float icwt_period[25][2][101][101];  ///pi0&eta combined of each period


void combinCalibConstantEndcapv1(){


  cout.setf(ios::fixed,ios::floatfield);
  cout.precision(8);
  

  for(int j=0; j<2; j++){
    for(int x =0; x<101; x++){
      for(int y=0; y<101; y++){
	validRecHitEndCap[j][x][y] = 0;
      }
    }
  }
  
  
  readInterCalibEndcap_GR09_V8();

  getInterCalibEndcapv1("interCalibEE_GR_P_V39.txt",ccalibpretag);
  
  
  get_xyzEBrechits();
  setEtaRingBoundaryEndcap();
  
  map<int,string> icFiles; 
  
  
  icFiles[0] = "calibres/deriveCalibConst.dflag2.pe1.step2.iter50.txt";
  icFiles[1] = "calibres/deriveCalibConst.dflag3.pe2.step2.iter50.txt";
  
  int nIC = 1; 
  
  getCrystaldeadflagEndcap_v1("deadflag/crystal_deadflag_ee_dflag2.txt",ndeadflag_ic[0]);
  getCrystaldeadflagEndcap_v1("deadflag/crystal_deadflag_ee_dflag3.txt",ndeadflag_ic[1]);
  
  
  int nConstantSet = int(icFiles.size());
  
  if( nConstantSet >50){
    cout<<"at most 50 IC " << nConstantSet <<endl; 
    return; 
  }
  
  
  
  float cc[2][101][101];
    
  TFile *fnew = new TFile("combineCalibConstantEndcapv1.root","recreate");
  
    

  TH1F *hh_cc_ietaring[50][3][40]; //file,ee-/ee+/both, ring, 
  TH1F *hh_res_cc_ietaring[50][3][3];
  
  //wted average cc
  TH1F *hh_ccwtavg_ietaring[3][40];
  TH1F *hh_res_ccwtavg_ietaring[3][3];

  TH1F *hh_ccwtavg_ietaring_period[25][3][40];
  TH1F *hh_res_ccwtavg_ietaring_period[25][3][3]; //of each period
  
  for(int n=0; n< nIC; n++){
    for(int j=0;j<3; j++){
      for(int k=0;k<40;k++){
	TString filename = TString (Form("hh_ccwtavg_ietaring_period_%d_%d_%d",n,j,k));
	hh_ccwtavg_ietaring_period[n][j][k] = new TH1F(filename,filename,200,0,4);
      }
    }
    for(int j=0;j<3; j++){
      for(int k=0;k<3;k++){
	TString filename = TString (Form("hh_res_ccwtavg_ietaring_period_%d_%d_%d",n,j,k));
	hh_res_ccwtavg_ietaring_period[n][j][k] = new TH1F(filename,filename,40,0,40);
      }
    }
  }
      
  
  TH2F *hh2_cc_xy[51][2];
  for(int j=0; j< nConstantSet+1; j++){
    for(int k=0;k<2; k++){
      TString filename = TString (Form("hh2_cc_xy_%d_%d",j,k));
      hh2_cc_xy[j][k] = new TH2F(filename,filename,100,1,101,100,1,101);
    }
  }
  
  TH2F *hh_corr_cc[50][10];
  for(int n=0; n<nConstantSet; n++){
    for(int k=n+1; k<nConstantSet; k++){
      TString filename = TString (Form("hh_corr_cc_%dand%d",n,k));
      hh_corr_cc[n][k] = new TH2F(filename,filename,400,0,2,400,0,2);
    }
  }
  
  
  for(int n=0; n<nConstantSet; n++){
    for(int j=0;j<3; j++){
      for(int k=0;k<40;k++){
	TString filename = TString (Form("hh_cc_ietaring_%d_%d_%d",n,j,k));
	hh_cc_ietaring[n][j][k] = new TH1F(filename,filename,200,0,4);
      }
    }
    for(int j=0;j<3; j++){
      for(int k=0;k<3;k++){
	TString filename = TString (Form("hh_res_cc_ietaring_%d_%d_%d",n,j,k));
	hh_res_cc_ietaring[n][j][k] = new TH1F(filename,filename,40,0,40);
      }
    }
  }
    
  for(int j=0;j<3; j++){
    for(int k=0;k<40;k++){
      TString filename = TString (Form("hh_ccwtavg_ietaring_%d_%d",j,k));
      hh_ccwtavg_ietaring[j][k] = new TH1F(filename,filename,200,0,4);
    }
  }
  for(int j=0;j<3; j++){
    for(int k=0;k<3;k++){
      TString filename = TString (Form("hh_res_ccwtavg_ietaring_%d_%d",j,k));
      hh_res_ccwtavg_ietaring[j][k] = new TH1F(filename,filename,40,0,40);
    }
  }

  
  //(c1 - c2) / ( average) 
  TH1F *hh_diff_cc_ietaring[50][50][3][40];
  TH1F *hh_res_diff_cc_ietaring[50][50][3][3];
  TH2F *hh2_diff_cc[50][50][2];
  
  
  for(int n=0; n<nConstantSet; n++){
    for(int k=n+1; k<nConstantSet; k++){
      for(int j=0;j<2;j++){
	TString filename = TString (Form("hh2_diff_cc_%dand%d_%d",n,k,j));
	hh2_diff_cc[n][k][j] = new TH2F(filename,filename,100,1,101,100,1,101);
      }
    }
  }
  
  
  for(int n=0; n<nConstantSet; n++){
    for(int k=n+1; k<nConstantSet; k++){

      for(int m=0; m<3; m++){
	for(int j=0;j<40;j++){
	  TString filename = TString (Form("hh_diff_cc_ietaring_%dand%d_%d_%d",n,k,m,j));
	  hh_diff_cc_ietaring[n][k][m][j] = new TH1F(filename,filename,100,-0.5,0.5);
	}
      }
      
      for(int m=0; m<3; m++){
	for(int j=0;j<3;j++){
	  TString filename = TString (Form("hh_res_diff_cc_ietaring_%dand%d_%d_%d",n,k,m,j));
	  hh_res_diff_cc_ietaring[n][k][m][j] = new TH1F(filename,filename,40,0,40);
	}
      }
    }
  }
  
  
  
  TH1F *hh_c_deadflag[51][30]; 
  TH1F *hh_c_deadflag_period[51][30]; 
  for(int j=0; j< nConstantSet+1; j++){
    for(int k=0; k<30; k++){
      TString histname = TString(Form("hh_c_deadflag_%d_%d",j,k));
      hh_c_deadflag[j][k] = new TH1F(histname,histname,500,0,2);

      histname = TString(Form("hh_c_deadflag_period_%d_%d",j,k));
      hh_c_deadflag_period[j][k] = new TH1F(histname,histname,500,0,2);
      
    }
  }
  
  
  

  ofstream txtoutcheck("combineCalibConstantEndcapv1.txt",ios::out);
    
  
  for(int n=0; n< nConstantSet; n++){
    string  filename = icFiles[n];


    getInterCalibEndcap(filename.c_str(),cc);


    
    NormCrystalDeadFlagEndcap_v1(cc,ndeadflag_ic[n]);
    
    

    scaleMeanToUnitEndcap(cc);
    for(int iz=0; iz<2; iz++){
      for(int j=0; j<101; j++){
	for(int k=0; k< 101; k++){
	  ic[n][iz][j][k] = cc[iz][j][k];

	  if( ndeadflag_ic[n][iz][j][k] < 0  && ic[n][iz][j][k] >0){
	    cout<<"warning (can be ignored) dead crystal  " <<ic[n][iz][j][k] <<endl; 
	    ic[n][iz][j][k] = -1; 
	  }
	  
	}
      }
    }
    
    
    
  }
  

  cout<<" nConstantSet " << nConstantSet <<endl; 
  

  for(int n=0; n< nConstantSet; n++){
    
     for(int iz=0; iz<2; iz++){
      for(int j=0; j<101; j++){
	for(int k=0; k< 101; k++){
	  
	  if( validRecHitEndCap[iz][j][k] <1) {
	    continue; 
	  }
	  
	  int iring = iRingEndCap(2*iz-1,j,k); ///input -1/1,  
	  float c = ic[n][iz][j][k];
	  if(c >0){
	    hh_cc_ietaring[n][iz][iring]->Fill(c);
	    hh_cc_ietaring[n][2][iring]->Fill(c);
	    hh2_cc_xy[n][iz]->SetBinContent(j,k,c);
	  }
	}
      }
     }
  }
  
  
  for(int n=0; n< nConstantSet; n++){
    for(int m=n+1; m < nConstantSet; m++){
      
      for(int iz=0; iz<2; iz++){
	for(int j=0; j<101; j++){
	  for(int k=0; k< 101; k++){
	    
	    if( validRecHitEndCap[iz][j][k] <1) {
	      hh2_diff_cc[n][m][iz]->SetBinContent(j,k,-1);
	      continue; 
	    }
	    int iring = iRingEndCap(2*iz-1,j,k); ///input -1/1,  
	    float c1 = ic[n][iz][j][k];
	    float c2 = ic[m][iz][j][k];
	    if(c1 >0 && c2 > 0 ){
	      hh_corr_cc[n][m]->Fill(c1,c2);
	      hh2_diff_cc[n][m][iz]->SetBinContent(j,k,c1-c2);
	      hh_diff_cc_ietaring[n][m][iz][iring]->Fill( (c1 -c2)/(0.5*(c1+c2)));
	      hh_diff_cc_ietaring[n][m][2][iring]->Fill( (c1 -c2)/(0.5*(c1+c2)));
	      
	    }else{
	      hh2_diff_cc[n][m][iz]->SetBinContent(j,k,-1);
	    }
	  }
	}
      }
    }
  }
  
  

  cout<<"fit " <<endl; 
  
  
   
  double resfit[20];

  

  
  for(int n=0; n< nConstantSet; n++){
    for(int m=n+1; m < nConstantSet; m++){
      for(int k=0; k< kEndcEtaRings; k++){
	
	for(int iz=0; iz<3; iz++){
	  fitHistogram(hh_diff_cc_ietaring[n][m][iz][k],resfit);
	  
	  for(int j=0; j<3; j++){
	    hh_res_diff_cc_ietaring[n][m][iz][j]->SetBinContent(k+1,resfit[2*j]);
	    hh_res_diff_cc_ietaring[n][m][iz][j]->SetBinError(k+1,resfit[2*j+1]);
	  }
	}
      }
    }
  }
  
  
  for(int n=0; n< nConstantSet; n++){
    for(int j=0; j<3; j++){
      for(int k=0; k< kEndcEtaRings; k++){
	
	fitHistogram(hh_cc_ietaring[n][j][k],resfit);
	for(int m=0;m<3; m++){
	  hh_res_cc_ietaring[n][j][m]->SetBinContent(k+1,resfit[2*m]);
	  hh_res_cc_ietaring[n][j][m]->SetBinError(k+1,resfit[2*m+1]);
	}
	
      }
    }
  }
  
  

  float sigmaSys = 2; ///precision of precalib+LC 
  
  
  cout<<"combine " <<endl; 
  float icwt[2][101][101];
  
  
  float wtSumC_period[25] = {0};
  float wtSumS_period[25] = {0};
    

  for(int iz=0; iz<2; iz++){
    for(int j=0; j<101; j++){
      for(int k=0; k< 101; k++){
	  
	if( validRecHitEndCap[iz][j][k] <1) {
	  continue; 
	}
	int iring = iRingEndCap(2*iz-1,j,k); ///input -1/1,  
	
	float wtSumC = 0; 
	float wtSumS = 0; 
	for(int n=0; n< nIC; n++){ ///for each period
	  wtSumC_period[n] = 0; 
	  wtSumS_period[n] = 0; 
	}
	
	for(int n=0; n< nConstantSet; n++){
	  float c = ic[n][iz][j][k]; 


	  float sigma = hh_res_cc_ietaring[n][2][2]->GetBinContent(iring+1);
	  sigma = sqrt( sigma * sigma + sigmaSys * sigmaSys);
	  
	  if( c > 0){
	    float tmp1 = c/ ( sigma * sigma);
	    float tmp2 = 1/(sigma * sigma);
	    
	    wtSumC += tmp1; 
	    wtSumS += tmp2; 
	    
	    int nperiod = n% nIC; 
	    wtSumC_period[nperiod] += tmp1; 
	    wtSumS_period[nperiod] += tmp2; 
	    
	    int deadflag = ndeadflag_ic[n][iz][j][k];
	    if(deadflag>=0){
	      hh_c_deadflag[n][deadflag]->Fill( c);
	    }

	  }
	}
	
	if( wtSumC > 0){
	  icwt[iz][j][k] = wtSumC / wtSumS; 


	  hh_ccwtavg_ietaring[iz][iring] ->Fill( icwt[iz][j][k] );
	  hh_ccwtavg_ietaring[2][iring] ->Fill( icwt[iz][j][k] );
	  
	}else{
	  icwt[iz][j][k] = -1;
	}
	for(int n=0; n< nIC; n++){ ///for each period
	  if( wtSumC_period[n] > 0){
	    icwt_period[n][iz][j][k] = wtSumC_period[n] / wtSumS_period[n]; 
	    
	    hh_ccwtavg_ietaring_period[n][iz][iring]->Fill( icwt_period[n][iz][j][k] );
	    hh_ccwtavg_ietaring_period[n][2][iring]->Fill( icwt_period[n][iz][j][k] );
	    
	  }else{
	    icwt_period[n][iz][j][k] = -1; 
	  }
	}
	
      }

    }
  }
  
  
  for(int j=0; j<3; j++){
    for(int k=0; k< kEndcEtaRings; k++){
      fitHistogram(hh_ccwtavg_ietaring[j][k],resfit);
      for(int m=0;m<3; m++){
	hh_res_ccwtavg_ietaring[j][m]->SetBinContent(k+1,resfit[2*m]);
	hh_res_ccwtavg_ietaring[j][m]->SetBinError(k+1,resfit[2*m+1]);
      }
    }
  }
    
  for(int n=0; n< nIC; n++){ ///for each period
    for(int j=0; j<3; j++){
      for(int k=0; k< kEndcEtaRings; k++){
	fitHistogram(hh_ccwtavg_ietaring_period[n][j][k],resfit);
	for(int m=0;m<3; m++){
	  hh_res_ccwtavg_ietaring_period[n][j][m]->SetBinContent(k+1,resfit[2*m]);
	  hh_res_ccwtavg_ietaring_period[n][j][m]->SetBinError(k+1,resfit[2*m+1]);
	}
      }
    }
  }
  
  
  scaleMeanToUnitEndcap(icwt);
    
  ofstream txtout("interCalibConstants.combined.EcalEndcap.txt",ios::out);
  
  
  ofstream txtout_period[50]; //pi0eta combined for each period
  for(int j=0; j< nIC; j++){
    string filename = string(Form("interCalibConstants.combinedPi0EtaPeriod%d.EcalEndcap.txt",j));
    txtout_period[j].open(filename.c_str(),ios::out);
  }
  
  
  cout<<"print out final " <<endl; 
  for(int iz=0; iz<2; iz++){
    for(int j=0; j<101; j++){
      for(int k=0; k< 101; k++){
        if( validRecHitEndCap[iz][j][k] <1) continue;

	
	int iring = iRingEndCap(2*iz-1,j,k); ///input -1/1,  
	
	for(int n=0; n< nIC; n++){ ///for each period
	  
	  float sigmaC = hh_res_ccwtavg_ietaring_period[n][2][2]->GetBinContent(iring+1);
	  
	  float cErr = sqrt( sigmaC *sigmaC + sigmaSys * sigmaSys);
	  cErr /=100;
	  
	  float c = icwt_period[n][iz][j][k];
	  if( c >0){
	    
	    int deadflag1 = ndeadflag_ic[n%nIC][iz][j][k];
	    int deadflag2 = ndeadflag_ic[nIC+n%nIC][iz][j][k];
	    if(deadflag1<0 && deadflag2<0){
	      cout<<"wrong deadflag ! " << n << " "<<iz <<" "<< j<<" "<<k <<endl; 
	      return; 
	    }
	    if(deadflag1>0 ){
	      hh_c_deadflag_period[n][deadflag1]->Fill( c);
	      hh_c_deadflag_period[n][19]->Fill( c);
	    }
	    
	  }
	  
	  if( c > 0){
	    txtout_period[n]<<j<<" "<<k<<" "<< 2*iz-1<<" "<< c*ccalibpretag[iz][j][k]<<" "<< cErr * c*ccalibpretag[iz][j][k] <<endl; 
	    
	  }else{	
	    txtout_period[n]<<j<<" "<<k<<" "<< 2*iz-1<<" "<<-1 <<" "<< 999 <<endl; 
	  }
	}
	
	
	float c = icwt[iz][j][k];
	if( c>0){
	  hh2_cc_xy[nConstantSet][iz]->SetBinContent(j,k,c);
	  
	  float sigmaC = hh_res_ccwtavg_ietaring[2][2]->GetBinContent(iring+1);
	  
	  float cErr = sqrt( sigmaC *sigmaC + sigmaSys * sigmaSys);
	  cErr /=100;
	  txtout<<j<<" "<<k<<" "<< 2*iz-1<<" "<< c * ccalibpretag[iz][j][k] <<" "<< cErr * c * ccalibpretag[iz][j][k] <<endl; 
	}else{
	  txtout<<j<<" "<<k<<" "<< 2*iz-1<<" "<<-1 <<" "<< 999 <<endl; 
	}
      }
    }
  }
  
  
  fnew->Write();
  fnew->Close();
  
  
  
}
