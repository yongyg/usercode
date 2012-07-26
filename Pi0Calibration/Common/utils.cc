int convertIetaIphiToSMNumber(int ieta, int iphi){
  
  if(  ! (abs(ieta)>=1 && abs(ieta)<=85 && iphi>=1 && iphi<=360)) {
    cout<<" convertSMNumber " << ieta <<" "<<iphi<<endl; 
    exit(1);
  }
  int ism = (iphi-1)/20+1;
  if( ieta<0) ism += 18;
  return ism; 
  
}



//information on the crystal's x/y/z needed to compute cluster position X/Y/Z after
//each inter-calibration step. 
void get_xyzEBrechits(){
  TChain *ch = new TChain("Analysis");
  ch->Add("/uscms_data/d2/yongy/tmp1/CMSSW_5_2_4/src/Pi0Calibration/Common/xyzECAL.root"); 
  
  ch->SetBranchAddress("xEBAll",xEBAll);
  ch->SetBranchAddress("yEBAll",yEBAll);
  ch->SetBranchAddress("zEBAll",zEBAll);
  ch->SetBranchAddress("etaEBAll",etaEBAll);
  ch->SetBranchAddress("phiEBAll",phiEBAll);
  
  ch->SetBranchAddress("dxEBAll",dxEBAll);
  ch->SetBranchAddress("dyEBAll",dyEBAll);
  ch->SetBranchAddress("dzEBAll",dzEBAll);
    
  ch->SetBranchAddress("xEEAll",xEEAll);
  ch->SetBranchAddress("yEEAll",yEEAll);
  ch->SetBranchAddress("zEEAll",zEEAll);
  
  ch->SetBranchAddress("dxEEAll",dxEEAll);
  ch->SetBranchAddress("dyEEAll",dyEEAll);
  ch->SetBranchAddress("dzEEAll",dzEEAll);
  
  
  
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
