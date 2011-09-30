#include "rootheader.h"
#include "roofitheader.h"


#include "usefullcode.cc"

using namespace RooFit ;

#include "usefullcoderoofit.cc"


///combined the workspace produced by testSelectionZee.C 
void makeZeeWsMCShape(char *test_datasetname, int total_evtRange){
  
  
  TString filename ;
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
  
  
  map<string,RooDataSet*> rhs_map_mpair;
  map<string,RooDataSet*> rhs_map_mpair_corr;

  RooRealVar *rv_mass = new RooRealVar("rv_mass","mass",100,60,120);
  RooRealVar *rv_weight = new RooRealVar("rv_weight","weight",1.0,0,1E6);
  
  for(int j=0; j< int( mpair_var.size()); j++){
    string mpairname = mpair_var[j];
    TString rname = TString( Form("rhs_%s",mpairname.c_str()) );
    rhs_map_mpair[mpairname] = new RooDataSet(rname,rname,RooArgList(*rv_mass,*rv_weight),rv_weight->GetName());
    
    mpairname = string (Form("%s_corr",mpair_var[j].c_str()));
    rname = TString( Form("rhs_%s",mpairname.c_str()) );
    rhs_map_mpair_corr[mpairname] = new RooDataSet(rname,rname,RooArgList(*rv_mass,*rv_weight),rv_weight->GetName());
    
  }
  
  float mMin = 0; 
  float mMax = 10000; 
    
  for(int n=1; n<=total_evtRange; n++){
    filename = TString (Form("selres/testSelectionZee.v3.%s.r%d.root",test_datasetname,n));
    cout<<filename<<endl; 
    TFile *f = new TFile(filename,"read");
    RooWorkspace* w = (RooWorkspace*) f->Get("zeeShape") ;
    
    for(int j=0; j< int( mpair_var.size()); j++){
      string mpairname = mpair_var[j];
      TString rname = TString( Form("rhs_%s",mpairname.c_str()) );
      RooDataSet *rd = (RooDataSet*)w->data(rname);
      appendcwd(rhs_map_mpair[mpairname],rd,rv_mass,rv_weight,mMin, mMax,1);
      
      rname = TString( Form("rhs_%s",mpairname.c_str()) ) + TString("_corr");
      rd = (RooDataSet*)w->data(rname);
      mpairname = mpair_var[j] + string("_corr");
      appendcwd(rhs_map_mpair_corr[mpairname],rd,rv_mass,rv_weight,mMin, mMax,1);

    }
  }
  
  filename = TString(Form("res/ZeeShape.v3.%s.r1to%d.root",test_datasetname, total_evtRange));
  cout<<filename<<endl; 
  TFile *fnew = new TFile(filename,"recreate");
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
