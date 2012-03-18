#include "rootheader.h"
#include "roofitheader.h"

#include "Pi0Analyzer.h"

TChain *fChain; 

#include "setBranchAddress.cc"
#include "getGoodLS.cc"

#include "physUtils.cc"
#include "selections.cc"

#include "loadl1seeds.cc"

#include "roodatasetth1.cc"

void print_L1rate(){
  
  int run = 178160;

//  int run = 177730;

  map<int,double> instantLum_run; 
  instantLum_run[178160] = 31.733/(456-51+1+534-458+1) *1./23.31 *1E36;
  instantLum_run[177730] = 105.532/(2237)  *1./23.31 *1E36;
  
  if( instantLum_run[run] == 0){
    cout<<"no instatenous lumi is given. " << run <<endl; 
    return; 
  }
  
  int nZeroBiasSelected = 110366;
  
  if( run == 177730){
    nZeroBiasSelected = 515727; 
  }
  
  TString filename = TString (Form("testSelection.run%d.root",run));
  TFile *f1 = new TFile(filename,"read");

  
  TH1F *hh_nL1firedAlca = (TH1F*)f1->Get("hh_nL1firedAlca");
  
  double lum = instantLum_run[run];
  double expectedLum[10] = {5E33,7E33};

  loadL1SeedsAndPrescale_5e33_7e33();


  float sumL1[2] = {0};
  
  for(int b=1; b<=100; b++){
    string l1 = hh_nL1firedAlca->GetXaxis()->GetBinLabel(b);
    if( l1 == ""){
      break; 
    }
    float nl1 = hh_nL1firedAlca->GetBinContent(b);
    
    float l1rate[2];
    float l1rateErr[2];
    for(int j=0; j<2; j++){

      float eff = nl1 / nZeroBiasSelected;
      float effErr = sqrt( eff *(1-eff)/nZeroBiasSelected);
      float projectedRate = expectedLum[j]/lum * eff * 11 * 1317; 
      float projectedRateErr = expectedLum[j]/lum * effErr * 11 * 1317; 
      l1rate[j] = projectedRate;
      l1rateErr[j] = projectedRateErr;
    }
    int ps1 = l1menu_prescale1[l1];
    int ps2 = l1menu_prescale2[l1];
    

    if(ps1 >=1 && ps2 >=1){
      TString output = TString (Form("%d & %2.2f & %d & %2.2f \\\\ \\hline ",ps1, l1rate[0]/ps1, ps2, l1rate[1]/ps2));
      cout<<l1.c_str()<<" & " << output<<endl; 
      
      sumL1[0] +=  l1rate[0]/ps1;
      sumL1[1] +=  l1rate[1]/ps2;
      
      
    }else{
      TString output = TString (Form(" -- & -- & -- & -- \\\\ \\hline "));
      cout<<l1.c_str()<<" & " << output<<endl; 
    }
    

  }
  //nL1Selected1 107 nL1AlcaSelected 34
  //nL1Selected2 73 nL1AlcaSelected 25
  
  cout<<"sumL1 " <<  sumL1[0] <<" "<<  sumL1[1]<<endl; 

  float nL1allps[2] = {107,73};
 
  if(run == 177730){
    nL1allps[0] = 293;
    nL1allps[1] = 211;
  }

  for(int j=0;j<2;j++){
    float n = nL1allps[j];
    float eff = n/ nZeroBiasSelected;
    float l1rate = expectedLum[j]/lum * eff * 11 * 1317; 
    cout<<"l1rate "<< l1rate <<endl; 
  }
  
}
