#include "rootheader.h"
#include "roofitheader.h"
#include "RecoAnalyzer.h"


#ifdef __MAKECINT__
#pragma link C++ class vector<vector<float> >+;
#pragma link C++ class vector<vector<unsigned short> >+;
#endif



#include "TMVAGui.C"

#if not defined(__CINT__) || defined(__MAKECINT__)
#include "TMVA/Tools.h"
#include "TMVA/Factory.h"
#include "TMVA/Reader.h"
#endif
using namespace TMVA;



#include "testSelection.h"

#include "ecalGap.cc"
#include "egammMVACorrection.cc"

int debug_ = -2; 


#include "setBranchAddress_ele.cc"

#include "physUtils.cc"
#include "usefullcode.cc"

#include "getGoodLS.cc"

#include "electronSelection.cc"

///re-weighting code 
#include "LumiReWeighting.cc" 


void testSelectionZee(char *test_datasetname, int test_evtRange){
  
  
  rgen_ = new TRandom3(0);
  cout<< rgen_->Gaus(0,0.01) <<endl; 
  
  
  dataversion = "v3";
  datasetname = test_datasetname;
  evtRange = test_evtRange; 
  
  
  egammaMVACorrection_LoadWeights();
    
  fChain  =new TChain("Analysis");
  TString filename ; 

  /// those files are preselected with two electron Et > 25GeV
  filename = TString( Form("/castor/cern.ch/user/y/yangyong/data/Run2011A/HiggsAnalysis/dielectronSkimmed/testAnalysisZee.%s.%s.allconv1.vtxmethod2.presel5.vtxtmva0.r%d.root",dataversion.c_str(),datasetname.c_str(),evtRange));
  cout<<filename<<endl;
  fChain->Add(filename);
  
  
  TChain *fevtInfo = new TChain("evtInfo");
  fevtInfo->Add(filename);
  fevtInfo->SetBranchAddress("totalNumberofEvents", &totalNumberofEvents);
  fevtInfo->SetBranchAddress("preselectedNumberofEvents", &preselectedNumberofEvents);
  fevtInfo->SetBranchAddress("datasettag", &datasettag);
  
  fevtInfo->GetEntry(0);
  string dname = datasettag->at(0);
  string dv = datasettag->at(1);
  int ntotal = totalNumberofEvents; 
  int nskimmed = preselectedNumberofEvents; 
  for(int n=1; n< fevtInfo->GetEntries(); n++){
    fevtInfo->GetEntry(n);
    ntotal += totalNumberofEvents; 
    nskimmed += preselectedNumberofEvents; 
  }
  
  
  setBranchAddress_ele();
  
  totalEntries = fChain->GetEntries();
  
  cout<<"totalEntries " << totalEntries <<endl; 
  
  ///get json file good lumis for each run 
  vector<string> certFile; 
  certFile.push_back("JsonFile/Cert_160404-163869_7TeV_May10ReReco_Collisions11_CMSSWConfig_v2.txtv1");
  certFile.push_back("JsonFile/Cert_160404-167913_7TeV_PromptReco_Collisions11_CMSSWConfig.txt.afterMay10");
  certFile.push_back("JsonFile/Cert_170249-172619_7TeV_ReReco5Aug_Collisions11_CMSSWConfig.txtv1");
  certFile.push_back("JsonFile/Cert_160404-173692_7TeV_PromptReco_Collisions11_CMSSWConfig.txtv1.run172620up");
  getLSrangeofEachRuns(certFile);

  
  TH1D *hh_pileupmc = new TH1D("hh_pileupmc","hh_npileup mc",51,-0.5,50.5);
  vector<double> weights_pileup;
  if( datasetname.find("SIM") != string::npos){
    /////for PU reweighting (added root files together)
    filename = TString( Form("PileUPDistrubtion/testAnalysis.%s.%s.puHist.root",dataversion.c_str(),datasetname.c_str()));
    TFile *fpumc = new TFile(filename,"read");
    TH1D *hh_npileupbx0 = (TH1D*)fpumc->Get("hh_npileupbx0");
    TH1D *hh_npileupbx1 = (TH1D*)fpumc->Get("hh_npileupbx1");
    TH1D *hh_npileupbxm1 = (TH1D*)fpumc->Get("hh_npileupbxm1");
    TH1D *hh_npileupbxav = (TH1D*)fpumc->Get("hh_npileupbxav"); //right now using averaged
    TFile *fpudata = new TFile("PileUPDistrubtion/Pileup_2011_EPS_8_jul.root","read");
    
    TH1D *hh_npileupdata = (TH1D*)fpudata->Get("pileup");
    weights_pileup = LumiReWeighting_getWeights(hh_npileupdata, hh_npileupbxav);
    fpumc->Close();
    fpudata->Close();
  }
  
  
  filename = TString (Form("selres/testSelectionZee.%s.%s.r%d.root",dataversion.c_str(), datasetname.c_str(),evtRange));
  cout<<filename<<endl; 
  TFile *fnew = new TFile(filename,"recreate");
  
  TH1I *hh_nevents = new TH1I("hh_nevents","total and skimm events",2,0,2);
  hh_nevents->SetBinContent(1,ntotal);
  hh_nevents->SetBinContent(2,nskimmed);
  
  
  vector<string> mpair_var;
  mpair_var.push_back("mpair_ebeb");
  mpair_var.push_back("mpair_ebeb_highr9");
  mpair_var.push_back("mpair_ebeb_lowr9");
  mpair_var.push_back("mpair_eeee");
  mpair_var.push_back("mpair_eeee_highr9");
  mpair_var.push_back("mpair_eeee_lowr9");
  mpair_var.push_back("mpair_eeee_highr9_loweps"); //eps_escraw < 0.05
  mpair_var.push_back("mpair_eeee_highr9_higheps");//eps_escraw > 0.05
  mpair_var.push_back("mpair_eeee_lowr9_loweps"); //eps_escraw < 0.05
  mpair_var.push_back("mpair_eeee_lowr9_higheps");//eps_escraw > 0.05
  
  
  
  RooRealVar *rv_mass = new RooRealVar("rv_mass","mass",100,0,1000);
  RooRealVar *rv_weight = new RooRealVar("rv_weight","weight",1.0,0,1E6);
  map<string,RooDataSet*> rhs_map_mpair;
  map<string,RooDataSet*> rhs_map_mpair_corr;
  for(int j=0; j< int( mpair_var.size()); j++){
    string mpairname = mpair_var[j];
    TString rname = TString( Form("rhs_%s",mpairname.c_str()) );
    rhs_map_mpair[mpairname] = new RooDataSet(rname,rname,RooArgList(*rv_mass,*rv_weight),rv_weight->GetName());
    mpairname = string (Form("%s_corr",mpair_var[j].c_str()));
    rname = TString( Form("rhs_%s",mpairname.c_str()) );
    rhs_map_mpair_corr[mpairname] = new RooDataSet(rname,rname,RooArgList(*rv_mass,*rv_weight),rv_weight->GetName());
  }
  map<string,TH1F*> hh_map_mpair; 
  map<string,TH1F*> hh_map_mpaircorr; 
  for(int j=0; j< int( mpair_var.size()); j++){
    string mpairname = mpair_var[j];
    filename = TString(Form("hh_%s",mpairname.c_str()));
    TH1F *hhtmp = new TH1F(filename,filename,2000,0,200);
    hh_map_mpair[mpairname] = hhtmp;
    hh_map_mpair[mpairname] ->Sumw2();
    
    mpairname = string (Form("%s_corr",mpair_var[j].c_str()));
    filename = TString(Form("hh_%s",mpairname.c_str())); 
    TH1F *hhtmp2 = new TH1F(filename,filename,2000,0,200);
    hh_map_mpaircorr[mpairname] = hhtmp2;
    hh_map_mpaircorr[mpairname] ->Sumw2();
  }
  
  
  
  bool goodCurLumiBlock = false; 
  int curLumiBlock = -1; 
  int curRun = -1; 
  vector<int> runNotUsed; 
  
  
  
  float en[2];
  float encorr[2];
  float pt[2];
  float eta[2];
  float phi[2];
  float res[20];
  float scr9[2];
  float scrps[2];
  
  float mpaircorr; 
  
  
  int nsel = 0; 
  
  bool firstEvent = true; 

  ////loop over all events ( di-photon skimmed)
  for(entry = 0; entry < totalEntries; entry ++){
    fChain->GetEntry(entry);
    if(entry % 1000 ==0 ) cout<<" entry " <<entry<<endl; 
    
    if( firstEvent){
      loadecalGapCoordinates();
      firstEvent = false; 
    }

    ///JSON file good run and LS
    if( isRealData ){
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
    
    for(int j=0; j< nElectron; j++){
      electronecalEnergy[j] = electronscenergy[j];  //right now ecalEnerg = scEnergy  
      float corr = getElectronCorrection(j); ///correction from regression 
      electronecalEnergyCorr[j] = (electronscrawEnergy[j] + electronscpreshowerEnergy[j]) * corr;
    }
    
    
    vector<int> indsel; 
    for(int j=0; j< nElectron; j++){
      if( electronecalDrivenSeed[j] ==0) continue; 
      float et = electronscenergy[j] * sin(2*atan(exp(-electroneta[j])));
      if( et < 30 ) continue;  ///ET > 30 
      if(! passElectronID(j,3)) continue; ///right now WP80 
      indsel.push_back(j);
    }
    if( int(indsel.size()) <2) continue; 
    
    nsel ++; 
    
    
    for(int j=0; j<2;j++){
      int k = indsel[j];
      en[j] = electronecalEnergy[k];
      encorr[j] = electronecalEnergyCorr[k];
      eta[j] = electrongsfTracketa[k];
      phi[j] = electrongsfTrackphi[k];
      scr9[j] = electrone3x3[k]/electronscrawEnergy[k];
      scrps[j] = electronscpreshowerEnergy[k] / electronscrawEnergy[k];
    }
    
    calcPairObjects(11,11,en,eta,phi,res);
    mpair = res[0];
    calcPairObjects(11,11,encorr,eta,phi,res);
    mpaircorr = res[0];
    
    
    double evtweight =1; 


    ///get pile up weghit 
    if( ! isRealData){
      int npubxav = 0; 
      for(unsigned int j =0; j < pileupBunchX->size(); j++){
	int BX = pileupBunchX->at(j);
	int npv = pileupNInteraction->at(j);
	if( BX >=-1 && BX <= 1){
	  npubxav += npv; 
	}
      }
      npvbxav = npubxav*1.0/3;
      int bin_weight = hh_pileupmc->FindFixBin(npvbxav) -1;
      double puwt = LumiReWeighting_weightINT_v2(bin_weight, weights_pileup);
      
      evtweight = puwt; 
      
    }
    
    
    int ind1 = indsel[0];
    int ind2 = indsel[1];
    int bothEB =  fabs(electronsceta[ind1]) < 1.482 && fabs(electronsceta[ind2]) < 1.482;
    int bothEE =  fabs(electronsceta[ind1]) > 1.482 && fabs(electronsceta[ind2]) > 1.482;
    int bothHighR9 = scr9[0] > 0.94 && scr9[1] > 0.94; 
    int bothLowR9 = scr9[0] < 0.94 && scr9[1] < 0.94; 
    int bothLowEPS = scrps[0] < 0.05 && scrps[1] < 0.05; 
    int bothHighEPS = scrps[0] > 0.05 && scrps[1] > 0.05; 
    
    
    
    bool fillflag[20] = {bothEB,
			 bothEB && bothHighR9,
			 bothEB && bothLowR9, 
			 bothEE, 
			 bothEE && bothHighR9, 
			 bothEE && bothLowR9,
			 bothEE && bothHighR9 && bothLowEPS, 
			 bothEE && bothHighR9 && bothHighEPS,
			 bothEE && bothLowR9 && bothLowEPS, 
			 bothEE && bothLowR9 && bothHighEPS};
    
    for(int j=0; j< int( mpair_var.size()); j++){
      if( fillflag[j]==false) continue;
      string mpairname = mpair_var[j];
      rv_mass->setVal(mpair);
      hh_map_mpair[mpairname]->Fill(mpair,evtweight);
      rhs_map_mpair[mpairname]->add(*rv_mass,evtweight);
      rv_mass->setVal(mpaircorr);
      mpairname = TString( Form("%s_corr",mpair_var[j].c_str()));
      hh_map_mpaircorr[mpairname]->Fill(mpaircorr,evtweight);
      rhs_map_mpair_corr[mpairname]->add(*rv_mass,evtweight);
    }
    
  }

  cout <<" nsel" << nsel <<endl; 
  
  
  fnew->cd();
  //RooContainer_Save();
  RooWorkspace *w = new RooWorkspace("zeeShape","workspace") ;
  for (std::map<string,RooDataSet*>::iterator it_data = rhs_map_mpair.begin()
	 ;it_data != rhs_map_mpair.end();it_data++)	{
    w->import(*it_data->second);
  }
  for (std::map<string,RooDataSet*>::iterator it_data = rhs_map_mpair_corr.begin()
	 ;it_data != rhs_map_mpair_corr.end();it_data++)	{
    w->import(*it_data->second);
  }
  w->Write();
  
  fnew->Write();
  fnew->Close();
}
