#include "rootheader.h"
#include "roofitheader.h"

#include "Pi0Analyzer.h"

TChain *fChain; 

#include "setBranchAddress.cc"

int evtRange;

//#include "datachain.cc"
#include "datachain.cc"


#include "getGoodLS.cc"

#include "physUtils.cc"
#include "selections.cc"

#include "loadl1seeds.cc"

#include "roodatasetth1.cc"



void testSelection(char *dataset,int test_evtRange){
  evtRange = test_evtRange; 
  
  
  l1bitFired_prescale1 = new vector<unsigned short> ; 
  l1bitFired_prescale2 = new vector<unsigned short> ; 
  l1bitFired_prescale3 = new vector<unsigned short> ; 

  map<int,double> instantLum_run; 

  instantLum_run[178160] = 31.733/(483) *1./23.31 *1E36;
  

//   if( instantLum_run[run] == 0){
//     cout<<"warning! no instatenous lumi is given. " << run <<endl; 
//     instantLum_run[run] = 3E33; //some dummy value
//     //return; 
//   }
  
  
  fChain = new TChain("clusters");
  
  TString filename; 

  ///filename = TString("/mnt/hadoop/user/yangyong/data/pizdata_rerun2012L1menu_2010b/") + TString(Form("%d",run)) + TString(Form("/clusters_%s_*.root",dataset));
  //cout<<filename<<endl; 
  //fChain->Add(filename);
  
  //  datachain(run);

  datachain();
  
  
  vector<string> certfiles;
  certfiles.push_back("/afs/cern.ch/cms/cit/yongy/data/JSON/Cert_160404-180252_7TeV_All2011_Nov30ReReco_v1.txtv1");
  getLSrangeofEachRuns(certfiles);
  
  
  setBranchAddress();
  
  
  
  filename = TString(Form("testSelection.%s.r%d.root",dataset,evtRange));
  TFile *fnew = new TFile(filename,"recreate");
  
  filename = TString(Form("testSelection.%s.r%d.txt",dataset,evtRange));
  ofstream txtout(filename,ios::out);
     
//   TH1F *hh_mpair_ebpiz = new TH1F("hh_mpair_ebpiz", "hh_mpair_ebpiz", 200,0,1);
//   TH1F *hh_mpair_ebeta = new TH1F("hh_mpair_ebeta", "hh_mpair_ebeta", 200,0,1);
//   TH1F *hh_mpair_eepiz = new TH1F("hh_mpair_eepiz", "hh_mpair_eepiz", 200,0,1);
//   TH1F *hh_mpair_eeeta = new TH1F("hh_mpair_eeeta", "hh_mpair_eeeta", 200,0,1);

  

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
  
  TH1F *hh_nalcaL1Selected = new TH1F("hh_nalcaL1Selected","hh_nalcaL1Selected",10,0,10);
  TH1F *hh_nalcaL1Selected_passebpiz = new TH1F("hh_nalcaL1Selected_passebpiz","hh_nalcaL1Selected_passebpiz",10,0,10);
  TH1F *hh_nalcaL1Selected_passeepiz = new TH1F("hh_nalcaL1Selected_passeepiz","hh_nalcaL1Selected_passeepiz",10,0,10);
  TH1F *hh_nalcaL1Selected_passebeta = new TH1F("hh_nalcaL1Selected_passebeta","hh_nalcaL1Selected_passebeta",10,0,10);
  TH1F *hh_nalcaL1Selected_passeeeta = new TH1F("hh_nalcaL1Selected_passeeeta","hh_nalcaL1Selected_passeeeta",10,0,10);
  TH1F *hh_nalcaL1Selected_passall = new TH1F("hh_nalcaL1Selected_passall","hh_nalcaL1Selected_passall",10,0,10);
  
  
  
  TH1F *hh_nL1fired = new TH1F("hh_nL1fired","L1 bit fired",200,0,200);
  TH1F *hh_nL1firedps1 = new TH1F("hh_nL1firedps1","L1 bit fired",200,0,200);
  TH1F *hh_nL1firedps2 = new TH1F("hh_nL1firedps2","L1 bit fired",200,0,200);
  
  TH1F *hh_nL1fired_passpizeb = new TH1F("hh_nL1fired_passpizeb","L1 bit fired && pizeb",200,0,200);
  TH1F *hh_nL1fired_passpizee = new TH1F("hh_nL1fired_passpizee","L1 bit fired && pizee",200,0,200);
  TH1F *hh_nL1fired_passetaeb = new TH1F("hh_nL1fired_passetaeb","L1 bit fired && etaeb",200,0,200);
  TH1F *hh_nL1fired_passetaee = new TH1F("hh_nL1fired_passetaee","L1 bit fired && etaee",200,0,200);
    
  TH1F *hh_nL1firedAlca = new TH1F("hh_nL1firedAlca","L1 bit fired",200,0,200);
  TH1F *hh_nL1firedAlca_passpizeb = new TH1F("hh_nL1firedAlca_passpizeb","L1 bit fired && pizeb",200,0,200);
  TH1F *hh_nL1firedAlca_passpizee = new TH1F("hh_nL1firedAlca_passpizee","L1 bit fired && pizee",200,0,200);
  TH1F *hh_nL1firedAlca_passetaeb = new TH1F("hh_nL1firedAlca_passetaeb","L1 bit fired && etaeb",200,0,200);
  TH1F *hh_nL1firedAlca_passetaee = new TH1F("hh_nL1firedAlca_passetaee","L1 bit fired && etaee",200,0,200);
  
  
  int psAlcaPi0EB = 1; 
  int psAlcaEtaEB = 1; 
  int psAlcaPi0EE = 1;
  int psAlcaEtaEE = 1;
  
  
  
  string alcaPi0UsedL1[200] = {
    "L1_SingleEG5", 
    "L1_SingleEG7",
    "L1_SingleEG12",
    "L1_SingleEG20",
    "L1_SingleEG22",
    "L1_SingleEG24",
    "L1_SingleEG30",
    "L1_DoubleEG_13_7",
    "L1_TripleEG7",
    "L1_TripleEG_12_7_5",
    "L1_DoubleEG5",
    "L1_TripleJet_64_44_24_VBF",
    "L1_TripleJet_64_48_28_VBF",
    "L1_TripleJetC_52_28_28",
    "L1_QuadJetC32",
    "L1_QuadJetC36",
    "L1_QuadJetC40",
    "L1_DoubleEG6_HTT100",
    "L1_DoubleEG6_HTT125",
    "L1_EG8_DoubleJetC20",
    "L1_Mu12_EG7",
    "L1_MuOpen_EG12",
    "L1_DoubleMu3p5_EG5",
    "L1_DoubleMu5_EG5",
    "L1_Mu12_EG7",
    "L1_Mu5_DoubleEG5",
    "L1_Mu5_DoubleEG6",
    "L1_MuOpen_EG5"
  };
  
  int nL1alca = 28;
  int isL1AlcaFired[200];
  int isL1Fired[200];
  
  int npizSel = 0; 
  int npizEB = 0; 
  int npizEE = 0; 
  
  int netaSel = 0; 
  int netaEB = 0; 
  int netaEE = 0; 

  
  ///for all l1
  int nL1bitsfired[200] = {0};
  int nL1bitsfiredps1[200] = {0};
  int nL1bitsfiredps2[200] = {0};
  

  ///only for alca
  int nL1bitsfiredAlca[200] = {0};
  int nL1bitsfiredAlcaps1[200] = {0};
  int nL1bitsfiredAlcaps2[200] = {0};
  
  int nL1bitsfiredAlca_passSelpizEB[200] = {0};
  int nL1bitsfiredAlca_passSelpizEE[200] = {0};
  int nL1bitsfiredAlca_passSeletaEB[200] = {0};
  int nL1bitsfiredAlca_passSeletaEE[200] = {0};
  
  
  int nL1bitsfired_passSelpizEB[200] = {0};
  int nL1bitsfired_passSelpizEE[200] = {0};
  int nL1bitsfired_passSeletaEB[200] = {0};
  int nL1bitsfired_passSeletaEE[200] = {0};
  
  
  
  
  bool goodCurLumiBlock = false; 
  int curLumiBlock = -1; 
  int curRun = -1; 
  vector<int> runNotUsed; 
  
  int nJsonSelected = 0; 
  int nZeroBiasSelected = 0; 

  int nL1Selected[10] = {0};
  int nL1AlcaSelected[10] = {0};
  int nL1AlcaSelected_passpizeb[10] = {0};
  int nL1AlcaSelected_passpizee[10] = {0};
  int nL1AlcaSelected_passetaeb[10] = {0};
  int nL1AlcaSelected_passetaee[10] = {0};
  
  
  loadL1SeedsAndPrescale_5e33_7e33();
  
  int totalEntries = fChain->GetEntries();
  
  cout<<" totalEntries " << totalEntries <<endl; 
  
  //totalEntries = 10000;
    
  for(int entry=0; entry< totalEntries; entry++){
    fChain->GetEntry(entry);

    
    ///define th1f
    if(entry==0){
      for(int j=0; j< int( l1algoName->size()); j++ ){
	string l1 = l1algoName->at(j);
	string histname = "mpair_pizeb_" +l1;	
	makeTH1F(histname,200,0,1);
	histname = "mpair_pizee_" +l1;
	makeTH1F(histname,200,0,1);
	histname = "mpair_etaeb_" +l1;
	makeTH1F(histname,200,0,1);
	histname = "mpair_etaee_" +l1;
	makeTH1F(histname,200,0,1);
      }
    }
    
    
    if(entry==0){
      checkMissingL1Prescale();
      cout<<" l1bitFired " << l1bitFired->size() <<" "<< l1algoName->size()<<endl; 
      for(int j=0; j< int( l1algoName->size()); j++ ){
	cout<<"l1algoName " << j <<" "<< l1algoName->at(j).c_str()<<endl; 
      }
      cout<<" hlt_bitFired " << hlt_bitFired->size()<<" "<< hlt_pathName->size() <<endl; 
    }
    
    if(entry % 1000 ==0) cout<<" entry " << entry <<endl; 
    
    
    ///JSON file good run and LS
    if( isRealData){
      vector<int>::iterator it = find(goodRunList.begin(),goodRunList.end(),runNumber); 
      if( it == goodRunList.end()){
	vector<int>::iterator it2 = find(runNotUsed.begin(),runNotUsed.end(),runNumber); 
	if( it2 == runNotUsed.end()){
	  runNotUsed.push_back(runNumber); 
	}
	continue; 
      }
      if( curLumiBlock != lumiBlock || curRun != runNumber){ /// a new lumiBlock  or starting of a new Run
	curLumiBlock = lumiBlock; 
	curRun =  runNumber; 
	goodCurLumiBlock = checkLumiBlockofRun();  //check this lumiBlock
      }
      if( ! goodCurLumiBlock) continue; 
      
    }
    nJsonSelected ++; 
    
    ////ZeroBias is fired at HLT
    
    if( ! isTriggerPathFiredv1("HLT_ZeroBias")) continue; 
    nZeroBiasSelected ++; 
    
    
    bool alcaL1Passed[3]; //nops ps1 ps2
    for(int j=0; j<3; j++){
      alcaL1Passed[j] = false;
    } 
    
    
    
    for(unsigned int j=0; j< l1algoName->size() ;j++){
      string l1 = l1algoName->at(j);
      isL1Fired[j] =  isL1PathFired(l1);
      if( isL1Fired[j]){
	nL1bitsfired[j] ++; 
      }
    }
    
    bool passAlcaPi0L1 = false; 
    for(int n=0; n< nL1alca; n++){
      string l1name = alcaPi0UsedL1[n];
      isL1AlcaFired[n] = isL1PathFired(l1name);
      if(isL1AlcaFired[n]){
	passAlcaPi0L1 = true; 
	alcaL1Passed[0] = true; 
	nL1bitsfiredAlca[n]++;
      }
    }
    prescale_L1seeds();
    
    
    for(unsigned int j=0; j< l1algoName->size() ;j++){
      string l1 = l1algoName->at(j);
      if( isL1PathFired_prescale(l1,1) ){
	nL1bitsfiredps1[j] ++; 
      }
      if( isL1PathFired_prescale(l1,2) ){
	nL1bitsfiredps2[j] ++; 
      }
    }
    
    
    for(int n=0; n< nL1alca; n++){
      string l1name = alcaPi0UsedL1[n];
      if(isL1PathFired_prescale(l1name,1) ){
	nL1bitsfiredAlcaps1[n]++;
	alcaL1Passed[1] = true; 
      }
      if(isL1PathFired_prescale(l1name,2) ){
	nL1bitsfiredAlcaps2[n]++;
	alcaL1Passed[2] = true; 
      }
    }
        
    for(int j=0; j<5; j++){
      if(alcaL1Passed[j])  nL1Selected[j] += 1; 
    }
    
    
    bool toomanyEB = nSeedsEB > 200 || n3x3ClusEB > 30 ; 
    bool toomanyEE = nSeedsEE > 200 || n3x3ClusEE > 30 ; 
    if( toomanyEB ){
      n3x3ClusEB = 0; 
    }
    if( toomanyEE ){
      n3x3ClusEE = 0; 
    }
    
    
    vector<float> mpizseleb = selection_EB_piz();
    bool passMassPi0EB = false; 
    for(int j=0; j< int(mpizseleb.size()); j++){
      mpair = mpizseleb[j];
      if( mpair >0.04 && mpair < 0.23){
	passMassPi0EB = true; 
      }
      for(int n=0; n<3; n++){
	if( alcaL1Passed[n]){
	  string histname = string(Form("mpair_ebpiz_l1alcaps_%d",n));
	  fillTH1F(histname,mpair,1);
	}
      }
      for(unsigned int n=0; n< l1algoName->size() ;n++){
	string l1 = l1algoName->at(n);
	string histname = "mpair_pizeb_" +l1;
	if( isL1Fired[n]){
	  fillTH1F(histname,mpair,1);
	}
      }
      
    }
    
    bool passMassPi0EE = false; 
    vector<float> mpizselee = selection_EE_piz();
    for(int j=0; j< int(mpizselee.size()); j++){
      mpair = mpizselee[j];
      if( mpair >0.05 && mpair < 0.3){
	passMassPi0EE = true; 
      }
      
      for(int n=0; n<3; n++){
	if( alcaL1Passed[n]){
	  string histname = string(Form("mpair_eepiz_l1alcaps_%d",n));
	  fillTH1F(histname,mpair,1);
	}
      }
      for(unsigned int n=0; n< l1algoName->size() ;n++){
	string l1 = l1algoName->at(n);
	string histname = "mpair_pizee_" +l1;
	if( isL1Fired[n]){
	  fillTH1F(histname,mpair,1);
	}
      }
      
    }
    
    bool passPi0EB =  passMassPi0EB;
    bool passPi0EE =  passMassPi0EE;
    if( passPi0EB){
      npizEB ++; 
    }
    if( passPi0EE){
      npizEE ++; 
    }
    bool passMassEtaEB = false; 
    vector<float> metaseleb = selection_EB_eta();
    for(int j=0; j< int(metaseleb.size()); j++){
      mpair = metaseleb[j];
      if( mpair >0.3 && mpair < 0.8){
	passMassEtaEB = true; 
      }
      for(int n=0; n<3; n++){
	if( alcaL1Passed[n]){
	  string histname = string(Form("mpair_ebeta_l1alcaps_%d",n));
	  fillTH1F(histname,mpair,1);
	}
      }
      for(unsigned int n=0; n< l1algoName->size() ;n++){
	string l1 = l1algoName->at(n);
	string histname = "mpair_etaeb_" +l1;
	if( isL1Fired[n]){
	  fillTH1F(histname,mpair,1);
	}
      }
      
    }
    
    bool passMassEtaEE = false; 
    vector<float> metaselee = selection_EE_eta();
    for(int j=0; j< int(metaselee.size()); j++){
      mpair = metaselee[j];
      
      for(int n=0; n<3; n++){
	if( alcaL1Passed[n]){
	  string histname = string(Form("mpair_eeeta_l1alcaps_%d",n));
	  fillTH1F(histname,mpair,1);
	}
      }
      for(unsigned int n=0; n< l1algoName->size() ;n++){
	string l1 = l1algoName->at(n);
	string histname = "mpair_etaee_" +l1;
	if( isL1Fired[n]){
	  fillTH1F(histname,mpair,1);
	}
      }
      if( mpair >0.2 && mpair < 0.9){
	passMassEtaEE = true; 
      }
    }
    
    bool passEtaEB = passMassEtaEB;
    bool passEtaEE = passMassEtaEE;
    
    if( passEtaEB){
      netaEB ++; 
    }
    if( passEtaEE){
      netaEE ++; 
    }
    
    ///for only alca 
    for(int n=0; n< nL1alca; n++){
      if(isL1AlcaFired[n]==1){
	if(passPi0EB) nL1bitsfiredAlca_passSelpizEB[n]++;
	if(passPi0EE) nL1bitsfiredAlca_passSelpizEE[n]++;
	if(passEtaEB) nL1bitsfiredAlca_passSeletaEB[n]++;
	if(passEtaEE) nL1bitsfiredAlca_passSeletaEE[n]++;
      }
    }
    
    ////for all L1 bit
    for(unsigned int j=0; j< l1algoName->size() ;j++){
      if(isL1Fired[j]==1){
	if( passPi0EB) nL1bitsfired_passSelpizEB[j]++;
	if( passPi0EE) nL1bitsfired_passSelpizEE[j]++;
	if( passEtaEB) nL1bitsfired_passSeletaEB[j]++;
	if( passEtaEE) nL1bitsfired_passSeletaEE[j]++;
      }
    }
    
    
    if( passPi0EB){
      if( npizEB % psAlcaPi0EB == 0){
	passPi0EB = true; 
      }else{
	passPi0EB = false; 
      }
    }
    if( passPi0EE){
      if( npizEE % psAlcaPi0EE == 0){
	passPi0EE = true; 
      }else{
	passPi0EE = false; 
      }
    }
    if( passEtaEB){
      if( netaEB % psAlcaEtaEB == 0){
	passEtaEB = true; 
      }else{
	passEtaEB = false; 
      }
    }
    if( passEtaEE){
      if( netaEE % psAlcaEtaEE == 0){
	passEtaEE = true; 
      }else{
	passEtaEE = false; 
      }
    }
        
    bool passAlca = passPi0EB || passEtaEB || passPi0EE || passEtaEE;
    
    for(int j=0; j<3; j++){
      if(alcaL1Passed[j]) {
	if( passAlca){
	  nL1AlcaSelected[j] ++; 
	}
	if(passPi0EB){
	  nL1AlcaSelected_passpizeb[j] ++; 
	}
	if(passPi0EE){
	  nL1AlcaSelected_passpizee[j] ++; 
	}
	if(passEtaEB){
	  nL1AlcaSelected_passetaeb[j] ++; 
	}
	if(passEtaEE){
	  nL1AlcaSelected_passetaee[j] ++; 
	}
	
      }
    }
    
  }
  
  ///for only alca l1
  for(int n=0; n< nL1alca; n++){
    string l1 = alcaPi0UsedL1[n];
    int b = n+1;
    hh_nL1firedAlca->SetBinContent(b,nL1bitsfiredAlca[n]);
    hh_nL1firedAlca_passpizeb->SetBinContent(b,nL1bitsfiredAlca_passSelpizEB[b-1]);
    hh_nL1firedAlca_passpizee->SetBinContent(b,nL1bitsfiredAlca_passSelpizEE[b-1]);
    hh_nL1firedAlca_passetaeb->SetBinContent(b,nL1bitsfiredAlca_passSeletaEB[b-1]);
    hh_nL1firedAlca_passetaee->SetBinContent(b,nL1bitsfiredAlca_passSeletaEE[b-1]);
    hh_nL1firedAlca->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1firedAlca_passpizeb->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1firedAlca_passpizee->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1firedAlca_passetaeb->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1firedAlca_passetaee->GetXaxis()->SetBinLabel(b,l1.c_str());
  }
  
  
  ///for all L1s
  for(int b=1; b<= int(l1algoName->size()) && b<=200; b++){
    string l1 = l1algoName->at(b-1);
    hh_nL1fired->SetBinContent(b,nL1bitsfired[b-1]);
    hh_nL1firedps1->SetBinContent(b,nL1bitsfiredps1[b-1]);
    hh_nL1firedps2->SetBinContent(b,nL1bitsfiredps2[b-1]);

    hh_nL1fired_passpizeb->SetBinContent(b,nL1bitsfired_passSelpizEB[b-1]);
    hh_nL1fired_passpizee->SetBinContent(b,nL1bitsfired_passSelpizEE[b-1]);
    hh_nL1fired_passetaeb->SetBinContent(b,nL1bitsfired_passSeletaEB[b-1]);
    hh_nL1fired_passetaee->SetBinContent(b,nL1bitsfired_passSeletaEE[b-1]);
    
    hh_nL1fired->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1firedps1->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1firedps2->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1fired_passpizeb->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1fired_passpizee->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1fired_passetaeb->GetXaxis()->SetBinLabel(b,l1.c_str());
    hh_nL1fired_passetaee->GetXaxis()->SetBinLabel(b,l1.c_str());
    
  }
  
  
  txtout<<" run " <<" totalEvents "<< totalEntries << " goodLum "<< nJsonSelected<<" ZeroBias "<< nZeroBiasSelected <<endl; 
  
  for(int n=0; n< nL1alca; n++){
    string l1name = alcaPi0UsedL1[n];
    txtout<<l1name.c_str()<<" " << nL1bitsfiredAlca[n]<<" "<< nL1bitsfiredAlcaps1[n]<<" "<< nL1bitsfiredAlcaps2[n]<<endl; 
    if( nL1bitsfiredAlca[n]>0){
      
      float effperL1_pizeb = 1.0* nL1bitsfiredAlca_passSelpizEB[n]/ nL1bitsfiredAlca[n];
      float effperL1_pizebErr = sqrt( effperL1_pizeb *(1-effperL1_pizeb)/nL1bitsfiredAlca[n]);
      float effperL1_pizee = 1.0* nL1bitsfiredAlca_passSelpizEE[n]/ nL1bitsfiredAlca[n];
      float effperL1_pizeeErr = sqrt( effperL1_pizee *(1-effperL1_pizee)/nL1bitsfiredAlca[n]);
      float effperL1_etaeb = 1.0 * nL1bitsfiredAlca_passSeletaEB[n]/ nL1bitsfiredAlca[n];
      float effperL1_etaebErr = sqrt( effperL1_pizeb *(1-effperL1_etaeb)/nL1bitsfiredAlca[n]);
      float effperL1_etaee = 1.0 * nL1bitsfiredAlca_passSeletaEE[n]/ nL1bitsfiredAlca[n];
      float effperL1_etaeeErr = sqrt( effperL1_pizeb *(1-effperL1_etaee)/nL1bitsfiredAlca[n]);
      
      txtout<<"efficiency_per_l1 " << l1name.c_str() <<" pizeb "<< effperL1_pizeb<<"+/-"<<effperL1_pizebErr<<" pizee "<< effperL1_pizee<<"+/-"<<effperL1_pizeeErr <<" etaeb " << effperL1_etaeb<<"+/-"<<effperL1_etaebErr<<" etaee " << effperL1_etaee<<"+/-"<<effperL1_etaeeErr<<endl; 
    }
  }
  
  txtout<<"npizSeleb/ee "<<" "<< npizEB <<" "<< npizEE <<" netaSeleb/ee " <<" "<< netaEB <<" "<< netaEE <<endl; 
  
  for(int j=0; j<3; j++){
    hh_nalcaL1Selected->SetBinContent(j+1,nL1Selected[j]);
    hh_nalcaL1Selected_passebpiz->SetBinContent(j+1,nL1AlcaSelected_passpizeb[j]);
    hh_nalcaL1Selected_passeepiz->SetBinContent(j+1,nL1AlcaSelected_passpizee[j]);
    hh_nalcaL1Selected_passebeta->SetBinContent(j+1,nL1AlcaSelected_passetaeb[j]);
    hh_nalcaL1Selected_passeeeta->SetBinContent(j+1,nL1AlcaSelected_passetaee[j]);
    hh_nalcaL1Selected_passall->SetBinContent(j+1,nL1AlcaSelected[j]);
    
    txtout<<"nL1Selected" << j<<" "<< nL1Selected[j]<<" nL1AlcaSelected "<<  nL1AlcaSelected[j]<<endl; 
  }
  
  
//   double lum = instantLum_run[run];
//   txtout<<"run " << run <<" lum " << lum <<endl; 
//   double expectedLum[10] = {5E33,7E33};
//   for(int j=1; j<3; j++){
//     float projectedRate = expectedLum[j-1]/lum * nL1AlcaSelected[j]/ nZeroBiasSelected * 11 * 1317; 
//     txtout<<" projectedRate at " << expectedLum[j-1]<<" is "<< projectedRate<<" kHz " <<endl; 
//   }
  
  
  fnew->Write();
  fnew->Close();
  
  
}
