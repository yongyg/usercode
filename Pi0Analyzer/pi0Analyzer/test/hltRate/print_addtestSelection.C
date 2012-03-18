#include "rootheader.h"
#include "roofitheader.h"
#include "/afs/cern.ch/user/y/yangyong/macros/pi0_mfitpeak_npol_v1.C"

void print_addtestSelection(){
  
  TFile *fnew = new TFile("addtestSelection.root","read");
  
  ofstream txtout("print_addtestSelection.txt",ios::out);
  
  
  TH1F *hh_nL1fired = (TH1F*)fnew->Get("hh_nL1fired");
  TH1F *hh_nL1fired_passpizeb = (TH1F*)fnew->Get("hh_nL1fired_passpizeb");
  TH1F *hh_nL1fired_passpizee = (TH1F*)fnew->Get("hh_nL1fired_passpizee");
  TH1F *hh_nL1fired_passetaeb = (TH1F*)fnew->Get("hh_nL1fired_passetaeb");
  TH1F *hh_nL1fired_passetaee = (TH1F*)fnew->Get("hh_nL1fired_passetaee");
  
  TH1F *hh_nL1fired_effpizeb = new TH1F("hh_nL1fired_effpizeb","",200,0,200);
  hh_nL1fired_effpizeb->Sumw2();
  TH1F *hh_nL1fired_effpizee = new TH1F("hh_nL1fired_effpizee","",200,0,200);
  hh_nL1fired_effpizee->Sumw2();
  TH1F *hh_nL1fired_effetaeb = new TH1F("hh_nL1fired_effetaeb","",200,0,200);
  hh_nL1fired_effetaeb->Sumw2();
  TH1F *hh_nL1fired_effetaee = new TH1F("hh_nL1fired_effetaee","",200,0,200);
  hh_nL1fired_effetaee->Sumw2();
  
  cout<<hh_nL1fired->Integral()<<endl;
  hh_nL1fired_effpizeb->Divide(hh_nL1fired_passpizeb,hh_nL1fired,1,1,"B");
  hh_nL1fired_effpizee->Divide(hh_nL1fired_passpizee,hh_nL1fired,1,1,"B");
  cout<<hh_nL1fired->Integral()<<endl;
  
  hh_nL1fired_effetaeb->Divide(hh_nL1fired_passetaeb,hh_nL1fired,1,1,"B");
  hh_nL1fired_effetaee->Divide(hh_nL1fired_passetaee,hh_nL1fired,1,1,"B");
  
  //return; 
  
  vector<string> l1alca;
  
  ifstream txtin("L1alca.txt",ios::in);
  string l1;
  while(txtin.good()){
    txtin>> l1;
    l1alca.push_back(l1);
    if(txtin.eof()) break; 
  }
  
  float res[10];
  
  map<string,int> l1alca_binnb; 
  for(int j=0; j<int(l1alca.size());j++){
    string l1 = l1alca[j];
    for(int b=1; b<=200; b++){
      string label = hh_nL1fired->GetXaxis()->GetBinLabel(b);
      if( label == l1){
	l1alca_binnb[l1] = b; 
	break; 
      }
    }
  }
  for(int j=0; j<int(l1alca.size());j++){
    string l1 = l1alca[j];
    int b = l1alca_binnb[l1];
    
    //  for(int b=1; b<=200; b++){
    TString label = hh_nL1fired->GetXaxis()->GetBinLabel(b);
    if(label=="") break; 
    
    float eff1 = hh_nL1fired_effpizeb->GetBinContent(b) * 100; 
    float eff1Err = hh_nL1fired_effpizeb->GetBinError(b) * 100; 
      
    float eff2 = hh_nL1fired_effpizee->GetBinContent(b) * 100; 
    float eff2Err = hh_nL1fired_effpizee->GetBinError(b) * 100; 
    
    float eff3 = hh_nL1fired_effetaeb->GetBinContent(b) * 100; 
    float eff3Err = hh_nL1fired_effetaeb->GetBinError(b) * 100; 
      
    float eff4 = hh_nL1fired_effetaee->GetBinContent(b) * 100; 
    float eff4Err = hh_nL1fired_effetaee->GetBinError(b) * 100; 
    

    TString output = TString(Form("$%2.1f \\pm %2.1f$ &  $%2.1f \\pm %2.1f$ & $%2.1f \\pm %2.1f$ & $%2.1f \\pm %2.1f$ \\\\ \\hline ", eff1, eff1Err,eff2, eff2Err, eff3, eff3Err,eff4, eff4Err));    
    txtout<<label<<" & "<<output<<endl; 
    
    
  }
  

  
  for(int j=0; j<int(l1alca.size());j++){
    string l1 = l1alca[j];
    
    TString hname = TString("th1f_mpair_pizeb_") +TString(l1);
    TH1F *hhtmp = (TH1F*)fnew->Get(hname);
    
    if( hhtmp->Integral(0.06/0.005,0.23/0.005)<500){
      hhtmp->Rebin(2);
    }


    if( l1 == "L1_Mu5_DoubleEG6" || l1 == "L1_Mu12_EG7" || l1 == "L1_DoubleEG6_HTT125" || l1 == "L1_Mu5_DoubleEG5" || l1 == "L1_QuadJetC32" || l1 == "L1_QuadJetC36" || l1 == "L1_QuadJetC40"){
      pi0_mfitpeak(hhtmp,0,0.06,0.23,1,res,1,"plots",hname,-1,-1,"");

    }else{
      pi0_mfitpeak(hhtmp,0,0.06,0.23,3,res,1,"plots",hname,-1,-1,"");
    }
    float S = res[0];
    float Serr = res[1];
    float nL1 = hh_nL1fired->GetBinContent(l1alca_binnb[l1]);
    float effS = S/nL1;
    float effSerr = sqrt(effS*(1-effS)/nL1);
    
    
//     hname = TString("th1f_mpair_etaeb_") +TString(l1);
//     TH1F *hhtmp1 = (TH1F*)fnew->Get(hname);
//     hhtmp1->Rebin(2);
     float S2 = 0; 
     float S2err = 0; 
//     if( hhtmp1->Integral(0.35/0.01,0.75/0.01) >200){
//       pi0_mfitpeak(hhtmp1,1,0.35,0.75,3,res,1,"plots",hname,-1,-1,"");
//     }
    S2 = res[0];
    S2err = res[1];
    float effS2 = S2/nL1;
    float effS2err = sqrt(effS2*(1-effS2)/nL1);
    TString output = TString(Form("$%2.1f \\pm %2.1f$\\\\ \\hline ", effS*100, effSerr*100));
    txtout<<l1.c_str()<<" & "<<output<<endl;
    
  }
  

}
