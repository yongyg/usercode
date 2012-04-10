
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



void sortPtvector(vector<float> pt, vector<int> &ind){
  
  ind.clear();

  int n = int( pt.size());

  vector<float> pt0; 
  for(int j=0; j< n ; j++){
    pt0.push_back(pt[j]);
  }
  
  sort(pt.begin(),pt.end());
  
  for(int j = n -1; j>=0; j--){
    
    float dptmin = 0.1; 
    int jmin = -1; 
    for(int k=0; k< n; k++){
      
      float dpt = fabs( pt[j] - pt0[k]); 
      vector<int>::iterator it = find( ind.begin(),ind.end(),k);
      if( it != ind.end()) continue; 
      
      if( dpt < dptmin){
	dptmin = dpt; 
	jmin = k; 
      }

    }
    ind.push_back(jmin);
    
  }
  
  if(int(ind.size()) != n){
    cout<<"erorr_sortPtvector.."<<endl; 
    exit(1);
  }
  
}


float  calculateS4(int nn, float en[],int ieta[], int iphi[]){
  ///shower-shape
  Float_t S4_0 =0.; 	 Float_t S4_1 =0.; 
  Float_t S4_2 =0.; 	 Float_t S4_3 =0.; 
  
  
  for(int n=0; n< nn; n++){
    float energy = en[n]; 
    int deta = diff_neta_s(ieta[0],ieta[n]);
    int dphi = diff_nphi_s(iphi[0],iphi[n]);
    if( abs(deta) <=1 && abs(dphi) <= 1){
      if( dphi <= 0 && deta <= 0) S4_0+=energy;
      if( dphi >= 0 && deta <= 0) S4_1+=energy;
      if( dphi <= 0 && deta >= 0) S4_2+=energy;
      if( dphi >= 0 && deta >= 0) S4_3+=energy;
    }
  }
  Float_t S4max=S4_0;
  if(S4_1 > S4max) S4max=S4_1;
  if(S4_2 > S4max) S4max=S4_2;
  if(S4_3 > S4max) S4max=S4_3;
  
  return S4max; 
  
}

void calculateS4_v1(int nn, float en[],int ieta[], int iphi[],float res[]){
  ///shower-shape
  Float_t S4_0 =0.; 	 Float_t S4_1 =0.; 
  Float_t S4_2 =0.; 	 Float_t S4_3 =0.; 

  float S9 = 0; 
    
  for(int n=0; n< nn; n++){
    float energy = en[n]; 
    int deta = diff_neta_s(ieta[0],ieta[n]);
    int dphi = diff_nphi_s(iphi[0],iphi[n]);
    if( abs(deta) <=1 && abs(dphi) <= 1){

      S9 += energy; 
      
      if( dphi <= 0 && deta <= 0) S4_0+=energy;
      if( dphi >= 0 && deta <= 0) S4_1+=energy;
      if( dphi <= 0 && deta >= 0) S4_2+=energy;
      if( dphi >= 0 && deta >= 0) S4_3+=energy;
    }
  }
  Float_t S4max=S4_0;
  if(S4_1 > S4max) S4max=S4_1;
  if(S4_2 > S4max) S4max=S4_2;
  if(S4_3 > S4max) S4max=S4_3;
  
  ////  return S4max; 
  
  res[0] = S4max; 
  res[1] = S9; 
  

}


///input ieta iphi originall 

void calculateS4_v2(int nn, float en[],int ieta[], int iphi[],float res[]){
  ///shower-shape
  Float_t S4_0 =0.; 	 Float_t S4_1 =0.; 
  Float_t S4_2 =0.; 	 Float_t S4_3 =0.; 

  float S9 = 0; 

  int ieta_seed = ieta[0];
  int iphi_seed = iphi[0];
  
  convxtalid(iphi_seed,ieta_seed);
  

    
  for(int n=0; n< nn; n++){
    float energy = en[n]; 

    int eta = ieta[n];
    int phi = iphi[n];
    convxtalid(phi,eta);
    
    int deta = diff_neta_s(ieta_seed,eta);
    int dphi = diff_nphi_s(iphi_seed,phi);
    
    if( abs(deta) <=1 && abs(dphi) <= 1){

      S9 += energy; 
      
      if( dphi <= 0 && deta <= 0) S4_0+=energy;
      if( dphi >= 0 && deta <= 0) S4_1+=energy;
      if( dphi <= 0 && deta >= 0) S4_2+=energy;
      if( dphi >= 0 && deta >= 0) S4_3+=energy;
    }
  }
  Float_t S4max=S4_0;
  if(S4_1 > S4max) S4max=S4_1;
  if(S4_2 > S4max) S4max=S4_2;
  if(S4_3 > S4max) S4max=S4_3;
  
  ////  return S4max; 
  
  res[0] = S4max; 
  res[1] = S9; 
  

}
