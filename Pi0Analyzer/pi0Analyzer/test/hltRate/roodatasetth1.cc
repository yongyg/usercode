

map<string,RooDataSet*> rds_map;
map<string,TH1F*> th1f_map;
map<string,TH2F*> th2f_map;
RooRealVar *rv_mass; 
RooRealVar *rv_weight; 




void makeRootDataSet(string sname){
  

  std::map<string,RooDataSet*>::iterator it_data = rds_map.find(sname);
  if( it_data != rds_map.end()){
    cout<<" makeRootDataSet already rds defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  
  //string sname = string (Form("cc%d_sm%d_ietabin%d",n,j,k));
  TString tname = TString("rds_") +  TString(sname);
  rds_map[sname] = new RooDataSet(tname,tname,RooArgList(*rv_mass,*rv_weight),rv_weight->GetName());
  
}

void makeRootDataSetAndTH1F(string sname, int nbins, float xlow, float xhigh){


  std::map<string,RooDataSet*>::iterator it_data = rds_map.find(sname);
  if( it_data != rds_map.end()){
    cout<<" makeRootDataSetAndTH1F already rds defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  
  //string sname = string (Form("cc%d_sm%d_ietabin%d",n,j,k));
  TString tname = TString("rds_") +  TString(sname);
  rds_map[sname] = new RooDataSet(tname,tname,RooArgList(*rv_mass,*rv_weight),rv_weight->GetName());
  tname = TString("th1f_") +  TString(sname);
  
  std::map<string,TH1F*>::iterator it_data2 = th1f_map.find(sname);
  if( it_data2 != th1f_map.end()){
    cout<<" makeRootDataSetAndTH1F already th1f defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  
  th1f_map[sname] = new TH1F(tname,tname,nbins, xlow, xhigh);
  th1f_map[sname] ->Sumw2();
  
}


void makeTH2F(string sname, int nbinx, float xlow, float xhigh, int nbiny, float ylow, float yhigh){
  //string sname = string (Form("cc%d_sm%d_ietabin%d",n,j,k));

  std::map<string,TH2F*>::iterator it_data2 = th2f_map.find(sname);
  if( it_data2 != th2f_map.end()){
    cout<<" makeTH2F already th2f defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  TString tname = TString("th2f_") +  TString(sname);
  th2f_map[sname] = new TH2F(tname,tname,nbinx, xlow, xhigh,nbiny,ylow,yhigh);
  
}


void makeTH1F(string sname, int nbinx, float xlow, float xhigh){
  
  std::map<string,TH1F*>::iterator it_data2 = th1f_map.find(sname);
  if( it_data2 != th1f_map.end()){
    cout<<" makeTH1F already th1f defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  
  TString tname = TString("th1f_") +  TString(sname);
  th1f_map[sname] = new TH1F(tname,tname,nbinx, xlow, xhigh);
  th1f_map[sname] ->Sumw2();
}

void makeTH1F(string sname, int nbinx, float xlow[]){
  
  std::map<string,TH1F*>::iterator it_data2 = th1f_map.find(sname);
  if( it_data2 != th1f_map.end()){
    cout<<" makeTH1F already th1f defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  
  TString tname = TString("th1f_") +  TString(sname);
  th1f_map[sname] = new TH1F(tname,tname,nbinx, xlow);
  th1f_map[sname] ->Sumw2();
}


void fillTH1F(string sname, float val, float wt){
  std::map<string,TH1F*>::iterator it_data2 = th1f_map.find(sname);
  if( it_data2 == th1f_map.end()){
    cout<<" fillTH1F not defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  th1f_map[sname]->Fill(val,wt);
  
}

void fillTH2F(string sname, float valx, float valy, float wt){
  std::map<string,TH2F*>::iterator it_data2 = th2f_map.find(sname);
  if( it_data2 == th2f_map.end()){
    cout<<" fillTH2F not defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  th2f_map[sname]->Fill(valx,valy,wt);
  
}

void setBinContentTH2F(string sname, int binx, int biny, float val){
  std::map<string,TH2F*>::iterator it_data2 = th2f_map.find(sname);
  if( it_data2 == th2f_map.end()){
    cout<<" seBinContentTH2F not defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  th2f_map[sname]->SetBinContent(binx, biny, val);
  
}


void setBinContentTH1F(string sname, int binx, float val,float err){
  std::map<string,TH1F*>::iterator it_data2 = th1f_map.find(sname);
  if( it_data2 == th1f_map.end()){
    cout<<" seBinContentTH1F not defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  th1f_map[sname]->SetBinContent(binx, val);
  th1f_map[sname]->SetBinError(binx, err);
  
}


void fillRootDataSet(string sname,float val, float wt){

  std::map<string,RooDataSet*>::iterator it_data = rds_map.find(sname);
  if( it_data == rds_map.end()){
    cout<<" fillRootDataSet rds not defined! " << sname.c_str()<<endl; 
    exit(1);
  }
  rv_mass->setVal(val);
  rds_map[sname]->add(*rv_mass,wt);
}

void fillRootDataSetAndTH1F(string sname,float val, float wt){


  std::map<string,RooDataSet*>::iterator it_data = rds_map.find(sname);
  if( it_data == rds_map.end()){
    cout<<" fillRootDataSetAndTH1F rds not defined! " << sname.c_str()<<endl; 
    exit(1);
  }
    
  rv_mass->setVal(val);
  ///string sname = string (Form("sumcc_smtbTTAbs_%d",ietaTTAbs+1));
  rds_map[sname]->add(*rv_mass,wt);
  th1f_map[sname]->Fill(val);
}
