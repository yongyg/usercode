


/////simply copy and scale weight 
RooDataSet *cwdset(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,Double_t weightscale) {
  
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();


    if( val > 700){
      cout<<"val " << val <<endl; 
    }
    
    if( val >vmin && val < vmax){
      mvar->setVal(val);
      outdata->add(*mvar,weightscale*indata->weight());
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
    
}


TH1F *convertRooDataSetToTH1F(RooDataSet *indata,RooRealVar *mvar, TString name,  double vmin, double vmax,Double_t weightscale){
  TH1F *outdata = new TH1F(name,name,60,vmin,vmax);
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    if( val >vmin && val < vmax){
      //mvar->setVal(val);
      //outdata->add(*mvar,weightscale*indata->weight());
      outdata->Fill(val,weightscale*indata->weight());
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
    
}


/////simply copy and scale weight 
RooDataSet *cwdsetv1(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,Double_t weightscale,  double scaleMass) {
  
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    val *= scaleMass;
    if( val >vmin && val < vmax){
      mvar->setVal(val);
      outdata->add(*mvar,weightscale*indata->weight());
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
  
}

/////simply copy , choose weight ( for data , choose runNumber)
RooDataSet *cwdsetv2(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, TString name, double vmin, double vmax,double wmin, double wmax,Double_t weightscale, double scaleMass) {
  
  RooDataSet *outdata = new RooDataSet(name,"",RooArgList(*mvar,*wvar),wvar->GetName());
  for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
    const RooArgSet *ent = indata->get(ient);
    double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
    val *= scaleMass;
    if(indata->weight() >= wmin && indata->weight() <= wmax &&  val >vmin && val < vmax){
      mvar->setVal(val);
      outdata->add(*mvar,weightscale*1);
      //cout<<"val " << val <<" "<< indata ->weight()<<endl; 
    }
  }
  
  return outdata;
  
}


/// addpend and scalw weight 
void appendcwd(RooDataSet *outdata, RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar,  double vmin, double vmax,Double_t weightscale) {
  
  
  
    for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
      const RooArgSet *ent = indata->get(ient);
      double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
      if(val >vmin && val < vmax){
	mvar->setVal(val);

	//cout<<"check " << val <<" " << indata->weight()<<endl;  
	
	outdata->add(*mvar,weightscale*indata->weight());
      }
      //      if( ient < 10) cout<<"checke " << mvar->getVal()<<" "<< indata->weight() <<endl; 
      
    }
    ///return outdata;
}




/// addpend and scalw weight 
void printRooDataSetToFile(RooDataSet *indata, RooRealVar *mvar, RooRealVar *wvar, char *filename) {
  
  ofstream testme(filename,ios::out);
  
    for (Int_t ient=0; ient<indata->numEntries(); ++ient) {
      const RooArgSet *ent = indata->get(ient);
      double val = static_cast<RooAbsReal*>(ent->find(mvar->GetName()))->getVal();
      testme<<  val<<" "<< indata->weight() <<endl; 
      
    }
    ///return outdata;
}



