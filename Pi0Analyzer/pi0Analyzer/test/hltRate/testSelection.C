#include "rootheader.h"
#include "Pi0Analyzer.h"

TChain *fChain; 

#include "setBranchAddress.cc"
#include "getGoodLS.cc"

#include "physUtils.cc"
#include "selections.cc"

#include "loadl1seeds.cc"

void testSelection(int run){
  
  l1bitFired_prescale1 = new vector<unsigned short> ; 
  l1bitFired_prescale2 = new vector<unsigned short> ; 
  l1bitFired_prescale3 = new vector<unsigned short> ; 
  
  
  fChain = new TChain("clusters");

  //fChain->Add("/uscms/home/marat/lpcegm/trypiz520/CMSSW_5_2_0/src/test/178160/crab_0_120315_013419/res/pizeta_6_1_l4k.root");
  //fChain->Add("/uscms/home/marat/lpcegm/trypiz520/CMSSW_5_2_0/src/test/178160/crab_0_120315_013419/res/pizeta_*.root");
  //fChain->Add("test.root");
  if( run == 178160){
    fChain->Add("crab_jobs/178160/crab_0_120315_124856/res/*root");
  }
  
  
  vector<string> certfiles;
  certfiles.push_back("Cert_160404-180252_7TeV_All2011_Nov30ReReco_v1.txtv1");
  getLSrangeofEachRuns(certfiles);
    
  
  setBranchAddress();
  

  
  TString filename = TString(Form("testSelection.run%d.root",run));
  TFile *fnew = new TFile(filename,"recreate");
  
  filename = TString(Form("testSelection.run%d.txt",run));
  ofstream txtout(filename,ios::out);
     
  TH1F *hh_mpair_ebpiz = new TH1F("hh_mpair_ebpiz", "hh_mpair_ebpiz", 200,0,1);
  TH1F *hh_mpair_ebeta = new TH1F("hh_mpair_ebeta", "hh_mpair_ebeta", 200,0,1);
  TH1F *hh_mpair_eepiz = new TH1F("hh_mpair_eepiz", "hh_mpair_eepiz", 200,0,1);
  TH1F *hh_mpair_eeeta = new TH1F("hh_mpair_eeeta", "hh_mpair_eeeta", 200,0,1);
    
  TH1F *hh_mpair_ebpiz_L1alca[200] ;
  TH1F *hh_mpair_ebeta_L1alca[200] ;
  for(int j=0; j< 200; j++){
    TString histname = TString (Form("hh_mpair_ebpiz_L1alca%d",j));
    hh_mpair_ebpiz_L1alca[j] = new TH1F(histname,histname,200,0,1);
    histname = TString (Form("hh_mpair_ebeta_L1alca_%d",j));
    hh_mpair_ebeta_L1alca[j] = new TH1F(histname,histname,200,0,1);
  }
  
  TH1F *hh_mpair_eepiz_L1alca[200] ;
  TH1F *hh_mpair_eeeta_L1alca[200] ;
  for(int j=0; j< 200; j++){
    TString histname = TString (Form("hh_mpair_eepiz_L1alca%d",j));
    hh_mpair_eepiz_L1alca[j] = new TH1F(histname,histname,200,0,1);
    histname = TString (Form("hh_mpair_eeeta_L1alca_%d",j));
    hh_mpair_eeeta_L1alca[j] = new TH1F(histname,histname,200,0,1);
  }

  TH1F *hh_nL1fired_alca = new TH1F("hh_nL1fired_alca","L1 bit fired",200,0,200);
  TH1F *hh_nL1fired_alcaps1 = new TH1F("hh_nL1fired_alcaps1","L1 bit fired",200,0,200);
  TH1F *hh_nL1fired_alcaps2 = new TH1F("hh_nL1fired_alcaps2","L1 bit fired",200,0,200);

   
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

  int npizSel = 0; 
  int npizEB = 0; 
  int npizEE = 0; 
  
  int netaSel = 0; 
  int netaEB = 0; 
  int netaEE = 0; 

  
  int nL1bitsfiredAlca[200] = {0};
  int nL1bitsfiredAlcaps1[200] = {0};
  int nL1bitsfiredAlcaps2[200] = {0};

  int nL1bitsfiredAlca_passSelpizEB[200] = {0};
  int nL1bitsfiredAlca_passSelpizEE[200] = {0};
  int nL1bitsfiredAlca_passSeletaEB[200] = {0};
  int nL1bitsfiredAlca_passSeletaEE[200] = {0};
  
  


  
  bool goodCurLumiBlock = false; 
  int curLumiBlock = -1; 
  int curRun = -1; 
  vector<int> runNotUsed; 
  
  int nJsonSelected = 0; 
  int nZeroBiasSelected = 0; 

  int nL1Selected[10] = {0};
  int nL1AlcaSelected[10] = {0};
  
  
  loadL1SeedsAndPrescale_5e33_7e33();
  
  int totalEntries = fChain->GetEntries();
  
  cout<<" totalEntries " << totalEntries <<endl; 
  
  //totalEntries = 10000;
    
  for(int entry=0; entry< totalEntries; entry++){
    fChain->GetEntry(entry);
    
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
    
    
    vector<float> mpizseleb = selection_EB_piz();
    bool passMassPi0EB = false; 
    for(int j=0; j< int(mpizseleb.size()); j++){
      mpair = mpizseleb[j];
      hh_mpair_ebpiz->Fill(mpair);
      if( mpair >0.04 && mpair < 0.23){
	passMassPi0EB = true; 
      }
      for(int n=0; n< nL1alca; n++){
	if( isL1AlcaFired[n]){
	  hh_mpair_ebpiz_L1alca[n] ->Fill(mpair);
	}
      }
    }

    bool passMassPi0EE = false; 
    vector<float> mpizselee = selection_EE_piz();
    for(int j=0; j< int(mpizselee.size()); j++){
      mpair = mpizselee[j];
      hh_mpair_eepiz->Fill(mpair);
      if( mpair >0.05 && mpair < 0.3){
	passMassPi0EE = true; 
      }
       
      for(int n=0; n< nL1alca; n++){
	if( isL1AlcaFired[n]){
	  hh_mpair_eepiz_L1alca[n] ->Fill(mpair);
	}
      }
    }
    bool passPi0EB = ( !toomanyEB &&  passMassPi0EB);
    bool passPi0EE = ( !toomanyEE &&  passMassPi0EE);
    if( passPi0EB){
      npizEB ++; 

      for(int n=0; n< nL1alca; n++){
	if(isL1AlcaFired[n]==1){
	  nL1bitsfiredAlca_passSelpizEB[n]++;
	}
      }
      
    }
    if( passPi0EE){
      npizEE ++; 
      for(int n=0; n< nL1alca; n++){
	if(isL1AlcaFired[n]==1){
	  nL1bitsfiredAlca_passSelpizEE[n]++;
	}
      }
    }
    bool passMassEtaEB = false; 
    vector<float> metaseleb = selection_EB_eta();
    for(int j=0; j< int(metaseleb.size()); j++){
      mpair = metaseleb[j];
      hh_mpair_ebeta->Fill(mpair);
      if( mpair >0.3 && mpair < 0.8){
	passMassEtaEB = true; 
      }
      for(int n=0; n< nL1alca; n++){
	if( isL1AlcaFired[n]){
	  hh_mpair_ebeta_L1alca[n] ->Fill(mpair);
	}
      }
    }
    
    bool passMassEtaEE = false; 
    vector<float> metaselee = selection_EE_eta();
    for(int j=0; j< int(metaselee.size()); j++){
      mpair = metaselee[j];
      hh_mpair_eeeta->Fill(mpair);

      if( mpair >0.2 && mpair < 0.9){
	passMassEtaEE = true; 
      }
      for(int n=0; n< nL1alca; n++){
	if( isL1AlcaFired[n]){
	  hh_mpair_eeeta_L1alca[n] ->Fill(mpair);
	}
      }
    }
    bool passEtaEB = ( !toomanyEB && passMassEtaEB);
    bool passEtaEE = ( !toomanyEE && passMassEtaEE) ;
    
    if( passEtaEB){
      netaEB ++; 
      for(int n=0; n< nL1alca; n++){
	if( isL1AlcaFired[n]){
	  nL1bitsfiredAlca_passSeletaEB[n]++;
	}
      }
    }
    if( passEtaEE){
      netaEE ++; 
      for(int n=0; n< nL1alca; n++){
	if( isL1AlcaFired[n]){
	  nL1bitsfiredAlca_passSeletaEE[n]++;
	}
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
      }
    }
    
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
    txtout<<"nL1Selected" << j<<" "<< nL1Selected[j]<<" nL1AlcaSelected "<<  nL1AlcaSelected[j]<<endl; 
  }
  
  double lum = 2.82E33; 
  double expectedLum[10] = {5E33,7E33};
  for(int j=1; j<3; j++){
    float projectedRate = expectedLum[j-1]/lum * nL1AlcaSelected[j]/ nZeroBiasSelected * 11 * 1317; 
    txtout<<" projectedRate at " << expectedLum[j-1]<<" is "<< projectedRate<<" kHz " <<endl; 
  }
  
   
  fnew->Write();
  fnew->Close();
  
  
}
