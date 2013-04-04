

int convertIetaIphiToSMNumber(int ieta, int iphi){
  
  if(  ! (abs(ieta)>=1 && abs(ieta)<=85 && iphi>=1 && iphi<=360)) {
    cout<<" convertSMNumber " << ieta <<" "<<iphi<<endl; 
    exit(1);
  }
  int ism = (iphi-1)/20+1;
  if( ieta<0) ism += 18;
  return ism; 
  
}

// Declaration of leaf types
double         xEBAll[170][360];
double         yEBAll[170][360];
double         zEBAll[170][360];

double         xctEBAll[170][360];
double         yctEBAll[170][360];
double         zctEBAll[170][360];



double         etaEBAll[170][360];
double         phiEBAll[170][360];
double         coxEBAll[170][360][8];
double         coyEBAll[170][360][8];
double         cozEBAll[170][360][8];
double         xEEAll[2][101][101];
double         yEEAll[2][101][101];
double         zEEAll[2][101][101];

//in the 2d plane
double         coxnewEBAll[170][360][8];
double         coynewEBAll[170][360][8];


double         xctEEAll[2][101][101];
double         yctEEAll[2][101][101];
double         zctEEAll[2][101][101];


double         etaEEAll[2][101][101];
double         phiEEAll[2][101][101];

double         coxEEAll[2][101][101][8];
double         coyEEAll[2][101][101][8];
double         cozEEAll[2][101][101][8];

double         coxnewEEAll[2][101][101][8];
double         coynewEEAll[2][101][101][8];


//information on the crystal's x/y/z needed to compute cluster position X/Y/Z after
//each inter-calibration step. 
void get_xyzEBrechits(){
  TChain *ch = new TChain("Analysis");
  ch->Add("xyzECAL.root");
  
  ch->SetBranchAddress("xEBAll",xEBAll);
  ch->SetBranchAddress("yEBAll",yEBAll);
  ch->SetBranchAddress("zEBAll",zEBAll);

  ch->SetBranchAddress("coxEBAll",coxEBAll);
  ch->SetBranchAddress("coyEBAll",coyEBAll);
  ch->SetBranchAddress("cozEBAll",cozEBAll);
  

  ch->SetBranchAddress("etaEBAll",etaEBAll);
  ch->SetBranchAddress("phiEBAll",phiEBAll);
  
//   ch->SetBranchAddress("dxEBAll",dxEBAll);
//   ch->SetBranchAddress("dyEBAll",dyEBAll);
//   ch->SetBranchAddress("dzEBAll",dzEBAll);
    
  ch->SetBranchAddress("xEEAll",xEEAll);
  ch->SetBranchAddress("yEEAll",yEEAll);
  ch->SetBranchAddress("zEEAll",zEEAll);


  ch->SetBranchAddress("coxEEAll",coxEEAll);
  ch->SetBranchAddress("coyEEAll",coyEEAll);
  ch->SetBranchAddress("cozEEAll",cozEEAll);

  
//   ch->SetBranchAddress("dxEEAll",dxEEAll);
//   ch->SetBranchAddress("dyEEAll",dyEEAll);
//   ch->SetBranchAddress("dzEEAll",dzEEAll);
  
  
  
  ch->SetBranchAddress("etaEEAll",etaEEAll);
  ch->SetBranchAddress("phiEEAll",phiEEAll);
  
  
  

  
  ch->GetEntry(0);

  ch->Delete();
  
  cout<<"got xyzEcal.."<<endl; 
  
  
    
}

/////////////////////////////////////////////////////////////////////////////////////////////////

void convxtalid(Int_t &nphi,Int_t &neta)
{
  // Changed to what Yong's convention; output will give just two indices
  // phi is unchanged; only eta now runs from
  //
  // 03/01/2008 changed to the new definition in CMSSW. The output is still the same...
  // Barrel only
  // Output nphi 0...359; neta 0...84; nside=+1 (for eta>0), or 0 (for eta<0).
  // neta will be [-85,-1] , or [0,84], the minus sign indicates the z<0 side.

     if(neta > 0) neta -= 1;
     if(nphi > 359) nphi=nphi-360;

  // final check
   if(nphi >359 || nphi <0 || neta< -85 || neta > 84)
     {
    cout <<" output not in range: "<<  nphi <<  " " << neta <<  " " <<endl;
    exit(1);
     }
} //end of convxtalid


// Calculate the distance in xtals taking into account possibly different sides
// change to coincide with yongs definition
Int_t diff_neta(Int_t neta1, Int_t neta2){
    Int_t mdiff;
    mdiff=abs(neta1-neta2);
    return mdiff;
}
 
// Calculate the absolute distance in xtals taking into account the periodicity of the Barrel
Int_t diff_nphi(Int_t nphi1,Int_t nphi2) {
  Int_t mdiff;
   mdiff=abs(nphi1-nphi2);
   if (mdiff > (360-abs(nphi1-nphi2))) mdiff=(360-abs(nphi1-nphi2));
   return mdiff;
}

// Calculate the distance in xtals taking into account possibly different sides
// Then the distance would be from the 1st to the 2nd argument
// _s means that it gives the sign; the only difference from the above !
// also changed to coincide with Yong's definition
Int_t diff_neta_s(Int_t neta1, Int_t neta2){
  Int_t mdiff;
  mdiff=(neta1-neta2);
  return mdiff;
}
 
// Calculate the distance in xtals taking into account the periodicity of the Barrel
Int_t diff_nphi_s(Int_t nphi1,Int_t nphi2) {
   Int_t mdiff;
   if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) {
     mdiff=nphi1-nphi2;
   }
   else {
   mdiff=360-abs(nphi1-nphi2);
   if(nphi1>nphi2) mdiff=-mdiff;
   }
   return mdiff;
}


///to access xEBAll[ieta][iphi]
/// input ieta -85,0,84
int getIndetaxyzEBAll(int ieta){
  return ieta+85; 
}

////something not consistent with 167,152?


///input 0, 359 after convxtalid
int getIndphixyzEBAll(int iphi){
  
  iphi = iphi-1; 
  if(iphi<0) return 359; 
  else return iphi;
  
}

double DeltaPhi(double phi1, double phi2){
  
  double diff = fabs(phi2 - phi1);
  
  while (diff >acos(-1)) diff -= 2*acos(-1);
  while (diff <= -acos(-1)) diff += 2*acos(-1);
  
  return diff; 
  
}


double GetDeltaR(double eta1, double phi1,double eta2, double phi2){
  
  return sqrt( (eta1-eta2)*(eta1-eta2) 
	       + DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
  
}




Float_t getcosd(Float_t eta1, Float_t phi1, Float_t eta2, Float_t phi2) {
  Float_t theta1 = 2*atan(exp(-eta1));
  Float_t theta2 = 2*atan(exp(-eta2));
  Float_t cosd;
  Float_t dphi = DeltaPhi(phi1,phi2);
  cosd = cos(theta1)*cos(theta2)+sin(theta1)*sin(theta2)*cos(dphi);  //opening angle 
  return cosd;
}
void separation(Float_t sceta1, Float_t scphi1, Float_t sceta2, Float_t scphi2, Float_t &dr)
{
  float dphi=fabs(scphi1-scphi2);
  if(dphi > (2*acos(-1)-fabs(scphi1-scphi2))) dphi = ( 2*acos(-1)-fabs(scphi1-scphi2));
  dr=sqrt((sceta1- sceta2)*(sceta1- sceta2)+dphi*dphi);
}



void phinorm(Float_t & PHI)
{
  while (PHI<0)  PHI= PHI + 2*acos(-1);
  if(PHI>2*acos(-1))  PHI= PHI - 2*acos(-1);
  
}


////change to [-pi,pi];
float phinorm2(float phi){
  while( phi > acos(-1) ) phi -= 2*acos(-1);
  while(phi< -acos(-1)) phi += 2*acos(-1);

  return phi;


}



void convertindTTindCrystal(int ietaT,int iphiT, int &ieta, int &iphi){

  if( iphiT >= 71) iphi = (iphiT-71)*5+3;
  else if( iphiT >=1 && iphiT <=70){
    iphi = (iphiT+1)*5+3;
  }else{
    cout<<"fatal error ..iphiT "<<iphiT<<endl;
    exit(1);
  }

  ieta = abs(ietaT)*5-2;

  if(ietaT<0) ieta *= -1;
  

  if( ieta >85 || ieta <-85 || ieta ==0){
    cout<<"convertindTTindCrystal...  "<<ieta<<" "<<ietaT<<endl;
    exit(1);
  }
  
  if( iphi <0 || iphi >360){
    cout<<"convertindTTindCrystal...  "<<iphi<<" "<<iphiT<<endl;
    exit(1);
  }
  
}


double ecalEta(double EtaParticle ,double Zvertex, double RhoVertex){
  
  
  //  const Double_t PI    = 3.1415927;
  double PI    = acos(-1);
  
  //---Definitions for ECAL
  double R_ECAL           = 136.5;
  double Z_Endcap         = 328.0;
  double etaBarrelEndcap  = 1.479; 

  if (EtaParticle!= 0.)
    {
      double Theta = 0.0  ;
      double ZEcal = (R_ECAL-RhoVertex)*sinh(EtaParticle)+Zvertex;
      
      if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
      if(Theta<0.0) Theta = Theta+PI;

      double ETA = - log(tan(0.5*Theta));
      
      if( fabs(ETA) > etaBarrelEndcap )
	{
	  double Zend = Z_Endcap ;
	  if(EtaParticle<0.0 )  Zend = -Zend ;
	  double Zlen = Zend - Zvertex ;
	  double RR = Zlen/sinh(EtaParticle);
	  Theta = atan((RR+RhoVertex)/Zend);
	  if(Theta<0.0) Theta = Theta+PI;
	  ETA = - log(tan(0.5*Theta));
	}
      return ETA;
    }
  else
    {
      return EtaParticle;
    }
}



double ecalPhi(double phi,double x0,double y0){
  
  //double R_ECAL = 136.5; ///cm 
  double r = 136.5; 

  double r0 = sqrt(x0*x0 + y0*y0);
  
  if(r0<1E-5) return phi; 
  
  if( r0 >= r){
    cout<<"warning. ecalPhi vtx outside ecal return input" << r0 <<" "<< r <<endl;
    return phi; 
  }
  
  double theta0 ;
  if(fabs(y0)>0) theta0= y0/fabs(y0) * acos(x0/r0);
  else theta0 = acos(x0/r0);
  
  ///  cout<<theta0<<" "<<phi<<endl;
  
  double theta = phi + asin( r0/r *sin(theta0-phi));

  //phinorm2(theta);
  double PI    = acos(-1);
  while ( theta < -PI) theta += PI; 
  while ( theta > PI) theta -= PI; 
  
  return theta; 
  
  
}



int simpleEcalID(int ieta, int iphi,int iz=0){
  if(iz==0) return ieta*360+iphi;
  else{
    return iz*(100000 + ieta * 100 + iphi);
  }
  
}


//return 5x5 crystals around maximum
vector<int> get5x5CrystalSim(float geta,float gphi,float vx,float vy,float vz){
  ////EB;-85,85  , 360 ;  ieta*360+iphi
  /// EE: 100*100;  iz*(100000 + ix * iy)
  
  float etag = ecalEta(geta,vz,sqrt(vx*vx+vy+vy));
  float phig = ecalPhi(gphi,vx,vy);
  
  if(debug) cout<<" eta/phi " << geta<<" "<<gphi<<" ecal "<<etag <<" "<<phig<<" "<<vz<<" e "<< ptGenPht[0]/sin(2*atan(exp(-geta)))<<endl;
  
  vector<int> vj; 
  if( fabs(etag)<1.5){

    float drmin = 0.1; 
    float emax = 0; 
    int jmax = -1; 
    for(int j=0; j< nsimEB;j++){
      
      int ieta1 = ietasimEB[j];
      int iphi1 = iphisimEB[j];
      convxtalid(iphi1,ieta1);
      int ieta = ieta1+85; 
      int iphi = getIndphixyzEBAll(iphi1);
      TVector3 vv(xEBAll[ieta][iphi],yEBAll[ieta][iphi],zEBAll[ieta][iphi]);
      float dr = GetDeltaR(etag,phig,vv.Eta(),vv.Phi());
      if(dr>drmin) continue; 
      float e = esumsimEB[j];
      if(emax <e) {
	emax = e; 
	jmax = j; 
      }
    }

    if(jmax<0) return vj; 

    vj.push_back(jmax);
    int ietam = ietasimEB[jmax];
    int iphim = iphisimEB[jmax];
    convxtalid(iphim,ietam);
    
    float esum = esumsimEB[jmax];
    
    for(int j=0; j< nsimEB;j++){
      if(j==jmax) continue; 
      int ieta1 = ietasimEB[j];
      int iphi1 = iphisimEB[j];
      convxtalid(iphi1,ieta1);

      if(diff_neta(ietam,ieta1)<=2 && diff_nphi(iphim,iphi1)<=2){
	vj.push_back(j);
	esum += esumsimEB[j];
      }
    }
    if(debug) cout<<"esum " << esum<<endl;
    
  }else{
    
    float drmin = 0.1; 
    float emax = 0; 
    int jmax = -1; 
    for(int j=0; j< nsimEE;j++){
      
      int ieta = ixsimEE[j];
      int iphi = iysimEE[j];
      int indz = izsimEE[j]>0; 
      TVector3 vv(xEEAll[indz][ieta][iphi],yEEAll[indz][ieta][iphi],zEEAll[indz][ieta][iphi]);
      float dr = GetDeltaR(etag,phig,vv.Eta(),vv.Phi());
      if(dr>drmin) continue; 
      float e = esumsimEE[j];
      if(emax <e) {
	emax = e; 
	jmax = j; 
      }
    }
    
    if(jmax<0) return vj; 

    vj.push_back(jmax);
    int ietam = ixsimEE[jmax];
    int iphim = iysimEE[jmax];

    float esum = esumsimEE[jmax];
    

    for(int j=0; j< nsimEE;j++){
      if(j==jmax) continue; 
      int ieta1 = ixsimEE[j];
      int iphi1 = iysimEE[j];
      if(diff_neta(ietam,ieta1)<=2 && diff_nphi(iphim,iphi1)<=2){
	vj.push_back(j);
	esum += esumsimEE[j];
      }
    }
    if(debug) cout<<"esum " << esum <<endl;

  }
  
  return vj;

}

//return 5x5 crystals around maximum
vector<int> get5x5CrystalStep(float geta,float gphi,float vx,float vy,float vz){
  ////EB;-85,85  , 360 ;  ieta*360+iphi
  /// EE: 100*100;  iz*(100000 + ix * iy)
  
  float etag = ecalEta(geta,vz,sqrt(vx*vx+vy+vy));
  float phig = ecalPhi(gphi,vx,vy);
  
  if(debug) cout<<" eta/phi " << geta<<" "<<gphi<<" ecal "<<etag <<" "<<phig<<" "<<vz<<" e "<< ptGenPht[0]/sin(2*atan(exp(-geta)))<<endl;
  
  vector<int> vj; 
  if( fabs(etag)<1.5){

    float drmin = 0.1; 
    float emax = 0; 
    int jmax = -1; 
    if(debug) cout<<"ng4EB " << ng4EB<<endl; 

    for(int j=0; j< ng4EB;j++){
      
      int ieta1 = ietag4EB[j];
      int iphi1 = iphig4EB[j];
      convxtalid(iphi1,ieta1);
      int ieta = ieta1+85; 
      int iphi = getIndphixyzEBAll(iphi1);
      TVector3 vv(xEBAll[ieta][iphi],yEBAll[ieta][iphi],zEBAll[ieta][iphi]);
      float dr = GetDeltaR(etag,phig,vv.Eta(),vv.Phi());
      if(dr>drmin) continue; 
      float e = esumg4EB[j];
      if(emax <e) {
	emax = e; 
	jmax = j; 
      }
    }

    if(jmax<0) return vj; 

    if(debug) cout<<"jmax "<< jmax <<endl; 

    vj.push_back(jmax);
    int ietam = ietag4EB[jmax];
    int iphim = iphig4EB[jmax];
    convxtalid(iphim,ietam);
    float esum = esumg4EB[jmax];

//     vector<float> ev = eg4EB->at(jmax);
//     float testesum = 0; 
//     for(int j=0; j<int(ev.size()); j++){
//       testesum += ev[j];
//     }

    if(debug) cout<<"emaxeb; " << esum <<endl; 

    for(int j=0; j< ng4EB;j++){
      if(j==jmax) continue; 
      int ieta1 = ietag4EB[j];
      int iphi1 = iphig4EB[j];
      convxtalid(iphi1,ieta1);

      if(diff_neta(ietam,ieta1)<=2 && diff_nphi(iphim,iphi1)<=2){
	vj.push_back(j);
	esum += esumg4EB[j];
      }
    }
    if(debug) cout<<"esum " << esum<<endl;
    
  }else{
    
    float drmin = 0.1; 
    float emax = 0; 
    int jmax = -1; 
    for(int j=0; j< ng4EE;j++){
      
      int ieta = ixg4EE[j];
      int iphi = iyg4EE[j];
      int indz = izg4EE[j]>0; 
      TVector3 vv(xEEAll[indz][ieta][iphi],yEEAll[indz][ieta][iphi],zEEAll[indz][ieta][iphi]);
      float dr = GetDeltaR(etag,phig,vv.Eta(),vv.Phi());
      if(dr>drmin) continue; 
      float e = esumg4EE[j];
      if(emax <e) {
	emax = e; 
	jmax = j; 
      }
    }
    
    if(jmax<0) return vj; 

    vj.push_back(jmax);
    int ietam = ixg4EE[jmax];
    int iphim = iyg4EE[jmax];

    float esum = esumg4EE[jmax];
    if(debug) cout<<"emaxee; " << esum <<endl;

    for(int j=0; j< ng4EE;j++){
      if(j==jmax) continue; 
      int ieta1 = ixg4EE[j];
      int iphi1 = iyg4EE[j];
      if(diff_neta(ietam,ieta1)<=2 && diff_nphi(iphim,iphi1)<=2){
	vj.push_back(j);
	esum += esumg4EE[j];
      }
    }
    if(debug)cout<<"esum " << esum <<endl;

  }
  
  return vj;

}
