

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
  
  //double diff = fabs(phi2 - phi1);
  double diff = phi1 - phi2;
  
  while (diff >acos(-1)) diff -= 2*acos(-1);
  while (diff <= -acos(-1)) diff += 2*acos(-1);
  
  return diff; 
  
}


double GetDeltaR(double eta1, double eta2, double phi1, double phi2){
  
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



void calcPairObjects(int pid1, int pid2,float en[],float eta[],float phi[],float res[]){
  
  TLorentzVector vpht[2];
  
  TLorentzVector vpair; 
  
  float mass[2]={0,0};
  
  if(pid1==11) mass[0] = 0.511*0.001; 
  else if(pid1==13) mass[0] = 0.105658;
  else if(pid1==22) mass[0] = 0;
  else{
    cout<<"calcPairMuPht, mass.pid1 "<<endl;
    exit(1);
  }
  
  if(pid2==11) mass[1] = 0.511*0.001; 
  else if(pid2==13) mass[1] = 0.105658;
  else if(pid2 ==22) mass[1] = 0; 
  else {
    cout<<"calcPairMuPht, mass. pid2"<<endl;
    exit(1);
  }
  

  if( en[0] < mass[0] || en[1] < mass[1]) {
    cout<<"warning calcPairObjects E<m? "<< en[0]<<" "<<en[1] <<" "<< pid1<<" "<<pid2<<endl; 
    exit(1); 
  }
  
  for( int j= 0; j<2; j++){
    
    float p = sqrt( en[j] * en[j] - mass[j]*mass[j]);
    float px = p / cosh(eta[j]) * cos(phi[j]);
    float py = p / cosh(eta[j]) * sin(phi[j]);
    float pz = p * tanh(eta[j]);
    vpht[j].SetXYZT(px,py,pz,en[j]);
  }
  
  vpair = vpht[0] + vpht[1];

  res[0] = vpair.M();
  res[1] = vpair.Eta();
  res[2] = vpair.Phi();
  res[3] = vpair.Pt();
  
  
  
  
}



void calcPairPhoton(float en[],float eta[],float phi[],float res[]){
  
  TLorentzVector vpht[2];
  
  TLorentzVector vpair; 
  
  
  for( int j= 0; j<2; j++){
    float e = en[j];
    float px = e /cosh(eta[j]) * cos(phi[j]);
    float py = e /cosh(eta[j]) * sin(phi[j]);
    float pz = e * tanh(eta[j]);
    vpht[j].SetXYZT(px,py,pz,e);
  }
  vpair = vpht[0] + vpht[1];
  res[0] = vpair.M();
  res[1] = vpair.Eta();
  res[2] = vpair.Phi();
  res[3] = vpair.Pt();
  
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




///input phi : [-pi,pi];
////this is the not exactly the same as; 
////RecoEcal/EgammaCoreTools/src/LogPositionCalc.cc
///need x,y,z information of each crystal to do that. 

/////center iphi = 190 or 191 phi change signs. 3.133, -3.133

void simpleLogWeightedEtaPhi(int nxt, float esum, float energy[],float eta[],float phi[],int phimax, float res[]){
  
  
  
  float etasum = 0; 
  float phisum = 0; 
  float wtsum = 0; 
  
  if(phimax==190 || phimax==191){
    for( int j=0; j<nxt; j++){
      phinorm(phi[j]);
    }
  }
    
  
  for( int j=0; j<nxt; j++){
    float mw=4.2+log(fabs(energy[j])/esum);
    if(mw < 0.) mw=0.;
    wtsum += mw; 
    etasum += mw * eta[j];
    phisum += mw * phi[j];
    
  }
  
  
  etasum /= wtsum; 
  phisum /= wtsum; 
  
  ///change to [-pi,pi]
  phinorm2(phisum);
  
  
  res[0] = etasum; 
  res[1] = phisum; 
  
  ////cluster shape

  

  
}

// // ///cluster shape SigmaEtaEta, SigmaEtaPhi,SigmaPhiPhi
// // ///


void Calculate_ClusterCovariance(int nxt, float esum,float ceta,float cphi,float en[],float eta[],float phi[],float res[]){
  
  double numeratorEtaEta = 0;
  double numeratorEtaPhi = 0;
  double numeratorPhiPhi = 0;
  double denominator     = 0;
  
  for( int j=0; j<nxt; j++){
    
    float dPhi = phi[j] - cphi; 

    if(dPhi > acos(-1)) dPhi = 2*acos(-1) - dPhi; 
    if(dPhi <-acos(-1)) dPhi = 2*acos(-1) + dPhi; 
    
    float dEta = eta[j] - ceta; 
    
    
    float w=4.2+log(fabs(en[j])/esum);
    if(w < 0.) w=0.;
    
    
    denominator += w;
    numeratorEtaEta += w * dEta * dEta;
    numeratorEtaPhi += w * dEta * dPhi;
    numeratorPhiPhi += w * dPhi * dPhi;
    
  }
  
     
  res[0] = numeratorEtaEta / denominator;
  res[1] = numeratorEtaPhi / denominator;
  res[2] = numeratorPhiPhi / denominator;
  
    


}
   



float Calculate_LAT(int nxt, float xclus, float yclus, float zclus,float en[],float x[],float y[],float z[]){
  
  if( nxt <3) return 10; 
  
  TVector3 clVect(xclus,yclus,zclus);
  
  TVector3 clDir = clVect; 
  clDir *= 1.0/clDir.Mag();
  
  float redmoment = 0; 
  
  float e12 = en[0] + en[1]; 
  TVector3 gblPos; 
  for( int j=2; j< nxt; j++){
    gblPos.SetXYZ(x[j],y[j],z[j]);
    TVector3 diff = gblPos - clVect;

    TVector3 DigiVect =  diff - diff.Dot(clDir)*clDir;
    float r = DigiVect.Mag();
    redmoment += r*r*fabs(en[j]);
  }
  
  float lat = redmoment/(redmoment + 2.19*2.19*e12);

  return lat; 
  

}


///change [-85,84] to bin 1, to bin 170
///change [0,359] to bin 1, 260; 

void getBinEtaPhi(int eta, int phi, float res[]){
  
  if( eta <0) {
    res[0] = eta + 85 +1+0.001; 
    
    res[2] = eta; 
    
  }
  else {
    res[0] = eta + 85 +2+0.001; 
    res[2] = eta +1; 
  }
  
  if( phi==0) {
    res[1] = 360; 
    res[3] = 360; 
  }
  else {
    res[1] = phi+0.001; 
    res[3] = phi; 
  }
  
  
  

}



////check crystal at boader of em obj boader

int IsatBoaderEMObjPhi(int iphi){
  
  // if( iphi ==0 || iphi ==1) return 1; 

  if( iphi%20 ==0 || (iphi-1)%20 ==0) return 1; 
  
  return 0; 

}

int IsatBoaderEMObjEta(int ieta){
  
  if( ieta ==0 || ieta == -1) return 1; 
  if( ieta ==19 || ieta==20) return 1; 
  if( ieta ==39 || ieta==40) return 1; 
  if( ieta == 59 || ieta ==60) return 1; 
  
  
  if( ieta==-20 || ieta==-21) return 1; 
  if( ieta==-40 || ieta==-41) return 1;
  if( ieta==-60 || ieta==-61) return 1;
    
  
  return 0; 

}



//convert ietaTT , iphiTT to ieta,iphi of central crystal

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


////invariant mass of two obejcts given, pid each. 
void calcPairPtEtaPhi(float pt[2],float eta[2],float phi[2],int pid1,int pid2,float res[]){
  //void calcPairMuPht(double pt[2],double eta[2],double phi[2],int pid1,int pid2,double res[]){
  
  TLorentzVector vpht[2];
  
  TLorentzVector vpair; 
  
  double mass[2]={0,0};
  
  if(pid1==11) mass[0] = 0.511*0.001; 
  else if(pid1==13) mass[0] = 0.105658;
  else if(pid1==22) mass[0] = 0;
  else{
    cout<<"calcPairMuPht, mass.pid1 "<<endl;
    exit(1);
  }
  
  if(pid2==11) mass[1] = 0.511*0.001; 
  else if(pid2==13) mass[1] = 0.105658;
  else if(pid2 ==22) mass[1] = 0; 
  else {
    cout<<"calcPairMuPht, mass. pid2"<<endl;
    exit(1);
  }
  
  
  for( int j= 0; j<2; j++){
    
    float p = pt[j] * cosh(eta[j]);
    float e = sqrt(p*p+mass[j]*mass[j]);
    float px = p /cosh(eta[j]) * cos(phi[j]);
    float py = p /cosh(eta[j]) * sin(phi[j]);
    float pz = p * tanh(eta[j]);
    vpht[j].SetXYZT(px,py,pz,e);
  }
  
  vpair = vpht[0] + vpht[1];
  
  res[0] = vpair.M();
  res[1] = vpair.Eta();
  res[2] = vpair.Phi();
  res[3] = vpair.Pt();
  res[4] = vpht[0].DeltaR(vpht[1]);
  
  
  ///cout<<"dot: "<< vpht[0].E() <<" "<< vpht[1].E()<< " "<< vpht[0].X() <<" "<< vpht[1].X()<<" "<< vpht[0].Y() <<" "<< vpht[1].Y() <<" "<< vpht[0].Z() <<" "<< vpht[1].Z()<<endl; 
  
  
  // res[5] = vpht[0].Dot(vpht[1])/vpht[0].P()/vpht[1].P();
  
  
  
}




void calcPairEtPtEtaPhi(float et[2],float pt[2],float eta[2],float phi[2],float res[]){
  
  TLorentzVector vpht[2];
  TLorentzVector vpair; 
  
  for( int j= 0; j<2; j++){
    float e  = et[j] *cosh(eta[j]);
    float p = pt[j] * cosh(eta[j]);
    float px = pt[j] * cos(phi[j]); 
    float py = pt[j] * sin(phi[j]); 
    float pz = p * tanh(eta[j]);
    vpht[j].SetXYZT(px,py,pz,e);
  }
  vpair = vpht[0] + vpht[1];
  res[0] = vpair.M();
  res[1] = vpair.Eta();
  res[2] = vpair.Phi();
  res[3] = vpair.Pt();
  res[4] = vpht[0].DeltaR(vpht[1]);
  
}



////invariant mass of two obejcts given, pid each. 
//void calcMass3Objects(float pt[3],float eta[3],float phi[3],int pid1,int pid2,int pid3,float res[]){
void calcMass3Objects(float pt[3],float eta[3],float phi[3],int pid1,int pid2,int pid3,float res[]){
  
  TLorentzVector vpht[3];
  
  TLorentzVector vpair; 
  
  float mass[3]={0,0,0};
  
  if(pid1==11) mass[0] = 0.511*0.001; 
  else if(pid1==13) mass[0] = 0.105658;
  else if(pid1==22) mass[0] = 0;
  else{
    cout<<"calcPairMuPht, mass.pid1 "<<endl;
    exit(1);
  }
  
  if(pid2==11) mass[1] = 0.511*0.001; 
  else if(pid2==13) mass[1] = 0.105658;
  else if(pid2 ==22) mass[1] = 0; 
  else {
    cout<<"calcPairMuPht, mass. pid2"<<endl;
    exit(1);
  }
  
  if(pid3==11) mass[2] = 0.511*0.001; 
  else if(pid3==13) mass[2] = 0.105658;
  else if(pid3 ==22) mass[2] = 0; 
  else {
    cout<<"calcPairMuPht, mass. pid3"<<endl;
    exit(1);
  }
    
  
  for( int j= 0; j<3; j++){
    
    float p = pt[j] * cosh(eta[j]);
    float e = sqrt(p*p+mass[j]*mass[j]);
    float px = p /cosh(eta[j]) * cos(phi[j]);
    float py = p /cosh(eta[j]) * sin(phi[j]);
    float pz = p * tanh(eta[j]);
    
    vpht[j].SetXYZT(px,py,pz,e);
    
  }
  
  vpair = vpht[0] + vpht[1] + vpht[2];
  
  res[0] = vpair.M();
  res[1] = vpair.Eta();
  res[2] = vpair.Phi();
  res[3] = vpair.Pt();
  res[4] = vpht[0].DeltaR(vpht[1]);
  
  
  

  
}






void printRatioTwoHist(TH1F *h1, TH1F *h2){
  
  int nbins = h1->GetNbinsX();
  
  float xmin = h1->GetXaxis()->GetXmin() ;
  float xmax = h1->GetXaxis()->GetXmax() ;
  
  TH1F *hh = new TH1F("hh","eff",nbins,xmin,xmax);
  
    
  for(int j=1; j<= nbins; j++){
    
    float y1 = h1->GetBinContent(j);
    float y2 = h2->GetBinContent(j);

    float eff = 0; 
    float effErr = 0; 
    if(y1>0 && y2>0){
      eff = y2/y1; 
      effErr = sqrt(eff*(1-eff)/y1);
      
      hh->SetBinContent(j,eff);
      hh->SetBinError(j,effErr);
    }
  }
  
  // TCanvas *can0 = new TCanvas("can0","c000",200,10,550,500);
  hh->Draw();
  
  
}

double etaTransformation(  float EtaParticle , float Zvertex)  {
  
  //---Definitions
  const float pi = 3.1415927;

  //---Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479; 
   
  //---ETA correction

  float Theta = 0.0  ; 
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+pi ;
  double ETA = - log(tan(0.5*Theta));
         
  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - Zvertex ;
      float RR = Zlen/sinh(EtaParticle); 
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+pi ;
      ETA = - log(tan(0.5*Theta));		      
    } 
  //---Return the result
  return ETA;
  //---end
}

////transform eta ( z, pho), to eta at ecal ( w.r.t 0,0,0,)
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



TVector3 ecalSurfaceXYZPositionGivenEtaPhi(float eta, float phi){

  //---Definitions for ECAL
  double R_ECAL           = 136.5;
  double Z_Endcap         = 328.0;

  if( fabs(eta)<1.48){
    float z = R_ECAL * sinh(eta);
    float x = R_ECAL * cos(phi);
    float y = R_ECAL * sin(phi);
    return TVector3(x,y,z);
  }else{
    float z = Z_Endcap;
    float x = R_ECAL / sinh(eta) * cos(phi);
    float y = R_ECAL / sinh(eta) * sin(phi);
    return TVector3(x,y,z);
  }
  
  
}
