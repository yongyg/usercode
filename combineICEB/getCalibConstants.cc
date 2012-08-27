

void scaleMeanToUnit(float C[170][360]){

  
  double mean = 0; 
  int n = 0; 
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      
      if( C[j][k] >0){
	mean += C[j][k];
	n ++; 
      }	
    }
  }
  mean /= n;
  cout<<"mean " << mean<<" "<<n<<endl; 

  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if(C[j][k] >0){
	C[j][k] /= mean; 
      }
    }
  }
  
}







void NormIetaToUnitTestBeamSMOnly(float C[170][360]){
  
  float mean_ieta[170] = {0}; 
  int nmean_ieta[170] = {0};

  float mean_ietaTT[170] = {0}; 
  int nmean_ietaTT[170] = {0};
  

  for(int j=0; j< 170; j++){
    
    int ieta = j-85; 
    if( ieta >=0) ieta += 1; 
    
    for(int k=0; k< 360; k++){
            
      int iphi = k; 
      if( k==0) iphi = 360; 
      
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      int smTB = isTestBeamSM(iSM);
      if( smTB ==0) continue; 
      
      //if( C[j][k] >0 ){
      if( C[j][k] >0.5 && C[j][k] < 1.5 ){
	mean_ieta[j] += C[j][k];
	nmean_ieta[j] ++; 

	mean_ietaTT[j/5] += C[j][k];
	nmean_ietaTT[j/5] ++; 
	
      }
      
    }
  }
  
  
  for(int j=0; j< 170; j++){
    mean_ieta[j] /= nmean_ieta[j];
    //cout<<" mean_ieta[" << j<<"] "<<  mean_ieta[j] <<endl; 
  }
   
//   for(int j=0; j< 34; j++){
//     mean_ietaTT[j] /= nmean_ietaTT[j];
//     cout<<" mean_ietaTT[" << j<<"] "<<  mean_ietaTT[j] <<endl; 
//   }
  
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      
      if(  C[j][k] >0 ){
	C[j][k] /=  mean_ieta[j] ; 
      }
    }
  }
  
}



void NormIetaAbsToUnitTestBeamSMOnly(float C[170][360]){
  
  
  float mean_ieta[170] = {0}; 
  int nmean_ieta[170] = {0};

  float mean_ietaTT[170] = {0}; 
  int nmean_ietaTT[170] = {0};
  

  for(int j=0; j< 170; j++){
    
    int ieta = j-85; 
    if( ieta >=0) ieta += 1; 
    
    for(int k=0; k< 360; k++){
            
      int iphi = k; 
      if( k==0) iphi = 360; 
      
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      int smTB = isTestBeamSM(iSM);
      if( smTB ==0) continue; 
      
      int ietaAbs = abs(ieta);
      
      //if( C[j][k] >0 ){
      if( C[j][k] >0.5 && C[j][k] < 1.5 ){
	mean_ieta[ietaAbs-1] += C[j][k];
	nmean_ieta[ietaAbs-1] ++; 
	
	mean_ietaTT[j/5] += C[j][k];
	nmean_ietaTT[j/5] ++; 
	
      }
    }
  }
  
  
  for(int j=0; j< 85; j++){
    mean_ieta[j] /= nmean_ieta[j];
    //cout<<" mean_ietaAbs[" << j<<"] "<<  mean_ieta[j] <<endl; 
  }
  
//   for(int j=0; j< 34; j++){
//     mean_ietaTT[j] /= nmean_ietaTT[j];
//     cout<<" mean_ietaTT[" << j<<"] "<<  mean_ietaTT[j] <<endl; 
//   }
  
  for(int j=0; j< 170; j++){
    
    int ieta = j-85; 
    if( ieta >=0) ieta += 1; 
    
    for(int k=0; k< 360; k++){
      
      int ietaAbs = abs(ieta);
      
      if(  C[j][k] >0 ){
		
	C[j][k] /=  mean_ieta[ietaAbs-1] ; 
	
      }
    }
  }
  
}



void NormCrystalDeadFlag_v1(float C[170][360], int ndeadflagietaiphi[170][360]){
  
  float mean_deadflag[20] = {0};
  float nmean_deadflag[20] = {0};
		       
  int ndeadCrystals = 0; 
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if( C[j][k] <0){
	ndeadCrystals++; 
	continue; 
      }
      int ndeadflag = ndeadflagietaiphi[j][k];
      if( ndeadflag >=0){
	mean_deadflag[ndeadflag] += C[j][k];
	nmean_deadflag[ndeadflag] ++; 
      }
    }
  }
  cout<<" NormCrystalDeadFlag ndeadCrystals: " << ndeadCrystals <<endl; 
  for(int j=0; j<20; j++){
    if( nmean_deadflag[j] >0){
      mean_deadflag[j] /= nmean_deadflag[j];
      cout<<" mean_deadflag "<<j<<"  "<< mean_deadflag[j]<< " "<< nmean_deadflag[j] <<endl;  
    }
  }
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if( C[j][k] <0) continue;
      int ndeadflag = ndeadflagietaiphi[j][k];
      if( ndeadflag >=0){
	C[j][k] /= mean_deadflag[ndeadflag];
      }
    }
  }
  
}

void NormCrystalDeadFlag(float C[170][360]){
  
  float mean_deadflag[20] = {0};
  float nmean_deadflag[20] = {0};
		       
  int ndeadCrystals = 0; 
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if( C[j][k] <0){
	ndeadCrystals++; 
	continue; 
      }
      int ndeadflag = ndeadflag_ietaiphi[j][k];
      if( ndeadflag >=0){
	mean_deadflag[ndeadflag] += C[j][k];
	nmean_deadflag[ndeadflag] ++; 
      }
    }
  }
  cout<<" NormCrystalDeadFlag ndeadCrystals: " << ndeadCrystals <<endl; 
  for(int j=0; j<20; j++){
    if( nmean_deadflag[j] >0){
      mean_deadflag[j] /= nmean_deadflag[j];
      cout<<" mean_deadflag "<<j<<"  "<< mean_deadflag[j]<< " "<< nmean_deadflag[j] <<endl;  
    }
  }
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      if( C[j][k] <0) continue;
      int ndeadflag = ndeadflag_ietaiphi[j][k];
      if( ndeadflag >=0){
	C[j][k] /= mean_deadflag[ndeadflag];
      }
    }
  }
  
}



void SetSMScale(float C[170][360],TH1F *hhtmp){
  float mean_ism[36] = {0};
  float nmean_ism[36]= {0};

  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){

      int ieta = j-85;
      if( j >=0) ieta += 1;
      int iphi = k;
      if( k==0) iphi = 360;
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      
      if(  C[j][k] >0 ){
        mean_ism[iSM-1 ]  += C[j][k];
        nmean_ism[iSM-1] ++;
      }
    }
  }

  for(int j=0; j<36; j++){
    mean_ism[j] /=  nmean_ism[j];
    ///    cout<<"sm mean: " << j+1 <<" "<< mean_ism[j] <<" "<< nmean_ism[j]<<endl;
  }
  
  
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){

      int ieta = j-85;
      if( j >=0) ieta += 1;
      int iphi = k;
      if( k==0) iphi = 360;
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      
      if(  C[j][k] >0 ){
        C[j][k] *= hhtmp->GetBinContent(iSM);
      }
    }
  }
  
  
}







void NormSMScaleToUnit(float C[170][360]){
  float mean_ism[36] = {0}; 
  float nmean_ism[36]= {0};
  
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){
      
      int ieta = j-85; 
      if( j >=0) ieta += 1; 
      int iphi = k; 
      if( k==0) iphi = 360; 
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      
      
      if(  C[j][k] >0 ){
	mean_ism[iSM-1 ]  += C[j][k]; 
	nmean_ism[iSM-1] ++; 
      }
    }
  }
  
  for(int j=0; j<36; j++){
    mean_ism[j] /=  nmean_ism[j];
    //cout<<"sm mean: " << j+1 <<" "<< mean_ism[j] <<" "<< nmean_ism[j]<<endl; 
  }
  
  for(int j=0; j< 170; j++){
    for(int k=0; k< 360; k++){

      int ieta = j-85; 
      if( j >=0) ieta += 1; 
      int iphi = k; 
      if( k==0) iphi = 360; 
      int iSM = (iphi-1)/20+1;
      if( ieta<0) iSM += 18;
      
      if(  C[j][k] >0 ){
	C[j][k] /=  mean_ism[iSM-1];
      }
      
    }
  }
  
}



void readInterCalibConstEBSimple(const char *input,float C[170][360]){
  
  ifstream txtin1(input,ios::in);
  if (txtin1.fail()){
    cout<<"error open file barrel.. " << input<<endl;
    exit(1);
  }
  
  
  int eta;
  int phi;
  
  float cc; 

  int n =0; 
  while(txtin1.good()){
    
    
    txtin1 >> eta >> phi >> cc; 
    convxtalid(phi,eta);
    eta += 85;
    C[eta][phi] = cc; 
    
    n++; 
    if( n>= 61200) break; 
    
  }
  
}





void readInterCalibConstEBSimplev1(const char *input,float C[170][360]){
  
  ifstream txtin1(input,ios::in);
  if (txtin1.fail()){
    cout<<"error open file barrel.. " << input<<endl;
    exit(1);
  }
  
  
  int eta;
  int phi;
  
  float cc; 
  double mean = 0; 
  int ngood = 0; 
  int n =0; 
  while(txtin1.good()){
    
    txtin1 >> eta >> phi >> cc; 
    
    if( !(eta>=0 && eta <= 169 && phi >=0 && phi<=359)) {
      cout<<"wrong eta/phi " << eta<<" "<<phi <<endl; 
      exit(1);
    }

    C[eta][phi] = cc; 
    if( cc >0){
      mean += cc; 
      ngood ++; 
    }
    n++; 
    if( n>= 61200) break; 
  }
  mean /= ngood; 
  cout<<"read ngood "<< ngood <<" mean: "<< mean<<endl; 
  
  for(int j=0; j<170; j++){
    for(int k=0; k<360;k++){
      if(C[j][k] >0) {
	C[j][k] /= mean;
      }
    }
  }
  
}




