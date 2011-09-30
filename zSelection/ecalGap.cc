

// EB gap positions
double _barrelCGap[169][360][2];
double _barrelSGap[33][180][2];
double _barrelTGap[33][180][2];
double _barrelMGap[7][18][2];

// EE crystal existence and gap positions
bool   _endcapCrystal[100][100];
double _endcapCGap[2][7080][2];
double _endcapSGap[2][264][2];
double _endcapMGap[2][1][2];

double _endcapCGapxy[2][7080][2];
double _endcapSGapxy[2][264][2];
double _endcapMGapxy[2][1][2];


// Actual data for each instantiated object
unsigned _be,_hl;
double _e,_eta,_phi,_r9;
double _aC,_aS,_aM,_bC,_bS,_bM;
double _aT,_bT;
double _xyaC,_xyaS,_xyaM,_xybC,_xybS,_xybM;


void loadecalGapCoordinates(){
  
  TChain *feg = new TChain("ecalGap");

  if(isRealData){
    feg->Add("GapGeometryFile/ecalGap.cmssw425.GeometryDB.FT_R_42_V17A.root"); ///data 
  }else{
    feg->Add("GapGeometryFile/ecalGap.cmssw425.GeometryDB.START42_V12_new.root"); //MC 
  }
  

  feg->SetBranchAddress("_barrelCGap", _barrelCGap);
  feg->SetBranchAddress("_barrelSGap", _barrelSGap);
  feg->SetBranchAddress("_barrelTGap", _barrelTGap);
  feg->SetBranchAddress("_barrelMGap", _barrelMGap);
  feg->SetBranchAddress("_endcapCGap", _endcapCGap);
  feg->SetBranchAddress("_endcapSGap", _endcapSGap);
  feg->SetBranchAddress("_endcapMGap", _endcapMGap);
  feg->SetBranchAddress("_endcapCGapxy", _endcapCGapxy);
  feg->SetBranchAddress("_endcapSGapxy", _endcapSGapxy);
  feg->SetBranchAddress("_endcapMGapxy", _endcapMGapxy);

  
  feg->GetEntry(0);
    
  cout<<"ecalGap Positions loaded " << isRealData <<" "<< endl; 
  
}

const double _onePi = acos(-1);
const double _twoPi = 2.0 * acos(-1);


double dPhi(double f0, double f1) {
  double df(f0-f1);
  if(df> _onePi) df-=_twoPi;
  if(df<-_onePi) df+=_twoPi;
  return df;
}

double aPhi(double f0, double f1) {
  double af(0.5*(f0+f1));
  if(fabs(dPhi(af,f0))>0.5*_onePi) {
    if(af>=0.0) af-=_onePi;
    else        af+=_onePi;
  }
  
  assert(fabs(dPhi(af,f0))<0.5*_onePi);
  assert(fabs(dPhi(af,f1))<0.5*_onePi);
  
  return af;
}


double xZ(double eta, double phi){
  assert(_be==1);
  return asinh(cos(phi)/sinh(eta));
}

double yZ(double eta, double phi) {
  assert(_be==1);
  return asinh(sin(phi)/sinh(eta));
}



void getGapCoordinates(double eta, double phi){
  // Check constants have been set up
  //assert(_initialised);
  
  // Determine if EB or EE
  _be=(fabs(eta)<1.482?0:1);
  
  //  // Determine if high or low R9
  //   if(_be==0) _hl=(_r9>=0.94?0:1);
  //   else       _hl=(_r9>=0.95?0:1);
  
  // Coordinates relative to cracks
  double r2Min;
  if(_be==0) {
    
    r2Min=1.0e6;
    for(unsigned i(0);i<169;i++) {
      for(unsigned j(0);j<360;j++) {
	double de(eta-_barrelCGap[i][j][0]);
	double df(dPhi(phi,_barrelCGap[i][j][1]));
	double r2(de*de+df*df);
	
	if(r2<r2Min) {
	  r2Min=r2;
	  if(i>=84) {
	    _aC= de;
	    _bC=-df;
	  } else {
	    _aC=-de;
	    _bC= df;
	  }
	}
      }
    }
    
    r2Min=1.0e6;
    for(unsigned i(0);i<33;i++) {
      for(unsigned j(0);j<180;j++) {
	double de(eta-_barrelSGap[i][j][0]);
	double df(dPhi(phi,_barrelSGap[i][j][1]));
	double r2(de*de+df*df);
	
	if(r2<r2Min) {
	  r2Min=r2;
	  if(i>=16) {
	    _aS= de;
	    _bS=-df;
	    
	  } else {
	    _aS=-de;
	    _bS= df;
	  }
	}
      }
    }
    
    r2Min=1.0e6;
    for(unsigned i(0);i<7;i++) {
      for(unsigned j(0);j<18;j++) {
	double de(eta-_barrelMGap[i][j][0]);
	double df(dPhi(phi,_barrelMGap[i][j][1]));
	double r2(de*de+df*df);
	
	if(r2<r2Min) {
	  r2Min=r2;
	  if(i>=3) {
	    _aM= de;
	    _bM=-df;
	  } else {
	    _aM=-de;
	    _bM= df;
	  }
	}
      }
    }
    
    
    r2Min=1.0e6;
    for(unsigned i(0);i<33;i++) {
      for(unsigned j(0);j<72;j++) {
	double de(eta-_barrelTGap[i][j][0]);
	double df(dPhi(phi,_barrelTGap[i][j][1]));
	double r2(de*de+df*df);
	
	if(r2<r2Min) {
	  r2Min=r2;
	  if(i>=16) {
	    _aT= de;
	    _bT=-df;
	  } else {
	    _aT=-de;
	    _bT= df;
	  }
	}
      }
    }
    
  } else {
    unsigned iz(eta>=0.0?0:1);
    double r[2]={xZ(eta,phi),yZ(eta,phi)};
    
    r2Min=1.0e6;
    for(unsigned i(0);i<7080;i++) {
      double dx(r[0]-_endcapCGap[iz][i][0]);
      double dy(r[1]-_endcapCGap[iz][i][1]);
      double r2(dx*dx+dy*dy);

      if(r2<r2Min) {
	r2Min=r2;
	if(r[0]>0.0) _aC= dx;
	else         _aC=-dx;
	if(r[1]>0.0) _bC= dy;
	else         _bC=-dy;
      }
    }
    
    r2Min=1.0e6;
    for(unsigned i(0);i<264;i++) {
      double dx(r[0]-_endcapSGap[iz][i][0]);
      double dy(r[1]-_endcapSGap[iz][i][1]);
      double r2(dx*dx+dy*dy);

      if(r2<r2Min) {
	r2Min=r2;
	if(r[0]>0.0) _aS= dx;
	else         _aS=-dx;
	if(r[1]>0.0) _bS= dy;
	else         _bS=-dy;
      }
    }
    
    r2Min=1.0e6;
    for(unsigned i(0);i<1;i++) {
      double dx(r[0]-_endcapMGap[iz][i][0]);
      double dy(r[1]-_endcapMGap[iz][i][1]);
      double r2(dx*dx+dy*dy);

      if(r2<r2Min) {
	r2Min=r2;
	if(iz==0) {_aM= dx;_bM= dy;}
	else      {_aM=-dx;_bM=-dy;}
      }
    }
  }
}



void getGapCoordinatesXY(double x, double y,double eta){
  // Check constants have been set up
  //assert(_initialised);
  
  // Determine if EB or EE
  ///_be=(fabs(eta)<1.482?0:1);
  
  if( fabs(eta) < 1.482) return;  //for endcap only
  
  unsigned iz(eta>=0.0?0:1);
  double r[2]={x,y};
  
  double r2Min=1.0e6;
  for(unsigned i(0);i<7080;i++) {
    double dx(r[0]-_endcapCGapxy[iz][i][0]);
    double dy(r[1]-_endcapCGapxy[iz][i][1]);
    double r2(dx*dx+dy*dy);

    if(r2<r2Min) {
      r2Min=r2;
      if(r[0]>0.0) _xyaC= dx;
      else         _xyaC=-dx;
      if(r[1]>0.0) _xybC= dy;
      else         _xybC=-dy;
    }
  }
    
  r2Min=1.0e6;
  for(unsigned i(0);i<264;i++) {
    double dx(r[0]-_endcapSGapxy[iz][i][0]);
    double dy(r[1]-_endcapSGapxy[iz][i][1]);
    double r2(dx*dx+dy*dy);

    if(r2<r2Min) {
      r2Min=r2;
      if(r[0]>0.0) _xyaS= dx;
      else         _xyaS=-dx;
      if(r[1]>0.0) _xybS= dy;
      else         _xybS=-dy;
    }
  }
    
  r2Min=1.0e6;
  for(unsigned i(0);i<1;i++) {
    double dx(r[0]-_endcapMGapxy[iz][i][0]);
    double dy(r[1]-_endcapMGapxy[iz][i][1]);
    double r2(dx*dx+dy*dy);

    if(r2<r2Min) {
      r2Min=r2;
      if(iz==0) {_xyaM= dx;_xybM= dy;}
      else      {_xyaM=-dx;_xybM=-dy;}
    }
  }
  
}

