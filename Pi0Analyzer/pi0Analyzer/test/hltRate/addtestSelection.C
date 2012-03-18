#include "rootheader.h"
#include "roofitheader.h"

#include "roodatasetth1.cc"

#include "pi0_mfitpeak_npol_v1.C"


void addtestSelection(){
  
  TFile *fnew = new TFile("addtestSelection.root","recreate");
  
  ///all L1 
  TH1F *hh_nL1fired = new TH1F("hh_nL1fired","L1 bit fired",200,0,200);
  TH1F *hh_nL1fired_passpizeb = new TH1F("hh_nL1fired_passpizeb","L1 bit fired && pizeb",200,0,200);
  TH1F *hh_nL1fired_passpizee = new TH1F("hh_nL1fired_passpizee","L1 bit fired && pizee",200,0,200);
  TH1F *hh_nL1fired_passetaeb = new TH1F("hh_nL1fired_passetaeb","L1 bit fired && etaeb",200,0,200);
  TH1F *hh_nL1fired_passetaee = new TH1F("hh_nL1fired_passetaee","L1 bit fired && etaee",200,0,200);
  
  
  //only alcal1 no ps, ps1,ps2
  TH1F *hh_nalcaL1Selected = new TH1F("hh_nalcaL1Selected","hh_nalcaL1Selected",10,0,10);
  TH1F *hh_nalcaL1Selected_passebpiz = new TH1F("hh_nalcaL1Selected_passebpiz","hh_nalcaL1Selected_passebpiz",10,0,10);
  TH1F *hh_nalcaL1Selected_passeepiz = new TH1F("hh_nalcaL1Selected_passeepiz","hh_nalcaL1Selected_passeepiz",10,0,10);
  TH1F *hh_nalcaL1Selected_passebeta = new TH1F("hh_nalcaL1Selected_passebeta","hh_nalcaL1Selected_passebeta",10,0,10);
  TH1F *hh_nalcaL1Selected_passeeeta = new TH1F("hh_nalcaL1Selected_passeeeta","hh_nalcaL1Selected_passeeeta",10,0,10);
  TH1F *hh_nalcaL1Selected_passall = new TH1F("hh_nalcaL1Selected_passall","hh_nalcaL1Selected_passall",10,0,10);
  
  

  ifstream txtin("l1algoName",ios::in);
  string l1; 
  vector<string> l1algo; 
  while(txtin.good()){
    txtin>> l1; 
    if(txtin.eof()) break; 
    l1algo.push_back(l1);
  }
  
  
  for(int j=0;j<3;j++){
    string histname = string(Form("mpair_ebpiz_l1alcaps_%d",j));
    makeTH1F(histname,200,0,1);
    histname = string(Form("mpair_eepiz_l1alcaps_%d",j));
    makeTH1F(histname,200,0,1);
    histname = string(Form("mpair_ebeta_l1alcaps_%d",j));
    makeTH1F(histname,200,0,1);
    histname = string(Form("mpair_eeeta_l1alcaps_%d",j));
    makeTH1F(histname,200,0,1);
  }

  for(int j=0; j<int(l1algo.size());j++){
    l1 = l1algo[j];
    string histname = "mpair_pizeb_" +l1;	
    makeTH1F(histname,200,0,1);
    histname = "mpair_pizee_" +l1;
    makeTH1F(histname,200,0,1);
    histname = "mpair_etaeb_" +l1;
    makeTH1F(histname,200,0,1);
    histname = "mpair_etaee_" +l1;
    makeTH1F(histname,200,0,1);
  }
  
  for(int r=1; r<=162; r++){
    TString filename = TString(Form("res/testSelection.MinimumBiasRun2011B-v1RAW.r%d.root",r));


    TString checkfilename = TString(Form("ls res/testSelection.MinimumBiasRun2011B-v1RAW.r%d.root",r));
    if( gSystem->Exec(checkfilename) !=0){
      continue; 
    }
    
    
    TFile *ff = new TFile(filename,"read");
    
    
    TH1F *hhtmp1 = (TH1F*)ff->Get("hh_nL1fired");
    hh_nL1fired->Add(hhtmp1);
    TH1F *hhtmp2 = (TH1F*)ff->Get("hh_nL1fired_passpizeb");
    hh_nL1fired_passpizeb->Add(hhtmp2);
    TH1F *hhtmp3 = (TH1F*)ff->Get("hh_nL1fired_passpizee");
    hh_nL1fired_passpizee->Add(hhtmp3);
    TH1F *hhtmp4 = (TH1F*)ff->Get("hh_nL1fired_passetaeb");
    hh_nL1fired_passetaeb->Add(hhtmp4);
    TH1F *hhtmp5 = (TH1F*)ff->Get("hh_nL1fired_passetaee");
    hh_nL1fired_passetaee->Add(hhtmp5);
    
    if(r==1){
      for(int b=1; b<= 200; b++){
	TString label = hhtmp1->GetXaxis()->GetBinLabel(b);
	hh_nL1fired->GetXaxis()->SetBinLabel(b,label);
	hh_nL1fired_passpizeb->GetXaxis()->SetBinLabel(b,label);
	hh_nL1fired_passpizee->GetXaxis()->SetBinLabel(b,label);
	hh_nL1fired_passetaeb->GetXaxis()->SetBinLabel(b,label);
	hh_nL1fired_passetaee->GetXaxis()->SetBinLabel(b,label);
      }
    }
    
    TH1F *hhtmp6 = (TH1F*)ff->Get("hh_nalcaL1Selected");
    hh_nalcaL1Selected->Add(hhtmp6);
    TH1F *hhtmp7 = (TH1F*)ff->Get("hh_nalcaL1Selected_passebpiz");
    hh_nalcaL1Selected_passebpiz->Add(hhtmp7);
    TH1F *hhtmp8 = (TH1F*)ff->Get("hh_nalcaL1Selected_passeepiz");
    hh_nalcaL1Selected_passeepiz->Add(hhtmp8);
    TH1F *hhtmp9 = (TH1F*)ff->Get("hh_nalcaL1Selected_passebeta");
    hh_nalcaL1Selected_passebeta->Add(hhtmp9);
    TH1F *hhtmp10 = (TH1F*)ff->Get("hh_nalcaL1Selected_passeeeta");
    hh_nalcaL1Selected_passeeeta->Add(hhtmp10);
    TH1F *hhtmp11 = (TH1F*)ff->Get("hh_nalcaL1Selected_passall");
    hh_nalcaL1Selected_passall->Add(hhtmp11);
    
    for(int j=0;j<3;j++){
      string histname = string(Form("mpair_ebpiz_l1alcaps_%d",j));
      TString hname = TString("th1f_") + histname; 
      TH1F *hh1 = (TH1F*)ff->Get(hname);
      th1f_map[histname]->Add(hh1);
      

      histname = string(Form("mpair_eepiz_l1alcaps_%d",j));
      hname = TString("th1f_") + histname; 
      TH1F *hh2 = (TH1F*)ff->Get(hname);
      th1f_map[histname]->Add(hh2);
      
      
      histname = string(Form("mpair_ebeta_l1alcaps_%d",j));
      hname = TString("th1f_") + histname; 
      TH1F *hh3 = (TH1F*)ff->Get(hname);
      th1f_map[histname]->Add(hh3);
      
      histname = string(Form("mpair_eeeta_l1alcaps_%d",j));
      hname = TString("th1f_") + histname; 
      TH1F *hh4 = (TH1F*)ff->Get(hname);
      th1f_map[histname]->Add(hh4);
      

    }
    
    
    
    for(int j=0; j<int(l1algo.size());j++){
      l1 = l1algo[j];
      string histname = "mpair_pizeb_" +l1;	
      TString hname = TString("th1f_mpair_pizeb_") +TString(l1);
      TH1F *hh1 = (TH1F*)ff->Get(hname);
      if( !hh1){
	cout<<"non hist " << r<<" "<< hname<<endl; 
	return; 
      }
      
      th1f_map[histname]->Add(hh1);
      
      histname = "mpair_pizee_" +l1;	
      hname = TString("th1f_mpair_pizee_") +TString(l1);
      TH1F *hh2 = (TH1F*)ff->Get(hname);
      th1f_map[histname]->Add(hh2);

      histname = "mpair_etaeb_" +l1;	
      hname = TString("th1f_mpair_etaeb_") +TString(l1);
      TH1F *hh3 = (TH1F*)ff->Get(hname);
      th1f_map[histname]->Add(hh3);
      
      histname = "mpair_etaee_" +l1;	
      hname = TString("th1f_mpair_etaee_") +TString(l1);
      TH1F *hh4 = (TH1F*)ff->Get(hname);
      th1f_map[histname]->Add(hh4);
      
    }
    
  }
  

  fnew->Write();
  fnew->Close();
  
  
//    "L1_SingleEG5",
//     "L1_SingleEG7",
//     "L1_SingleEG12",
//     "L1_SingleEG20",
//     "L1_SingleEG22",
//     "L1_SingleEG24",
//     "L1_SingleEG30",
//     "L1_DoubleEG_13_7",
//     "L1_TripleEG7",
//     "L1_TripleEG_12_7_5",
//     "L1_DoubleEG5",
//     "L1_TripleJet_64_44_24_VBF",
//     "L1_TripleJet_64_48_28_VBF",
//     "L1_TripleJetC_52_28_28",
//     "L1_QuadJetC32",
//     "L1_QuadJetC36",
//     "L1_QuadJetC40",
//     "L1_DoubleEG6_HTT100",
//     "L1_DoubleEG6_HTT125",
//     "L1_EG8_DoubleJetC20",
//     "L1_Mu12_EG7",
//     "L1_MuOpen_EG12",
//     "L1_DoubleMu3p5_EG5",
//     "L1_DoubleMu5_EG5",
//     "L1_Mu12_EG7",
//     "L1_Mu5_DoubleEG5",
//     "L1_Mu5_DoubleEG6",
//     "L1_MuOpen_EG5"
  

}
