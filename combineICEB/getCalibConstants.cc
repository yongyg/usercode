

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




// void scaleMeanToUnitEndcap(float C[2][101][101]){

  
//   double mean = 0; 
//   int n = 0; 
//   for(int j=0; j< 170; j++){
//     for(int k=0; k< 360; k++){
      
//       if( C[j][k] >0){
// 	mean += C[j][k];
// 	n ++; 
//       }	
//     }
//   }
//   mean /= n;
//   cout<<"mean " << mean<<" "<<n<<endl; 

//   for(int j=0; j< 170; j++){
//     for(int k=0; k< 360; k++){
//       if(C[j][k] >0){
// 	C[j][k] /= mean; 
//       }
//     }
//   }
  
// }




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

    int iphi = phi;
    int ieta = eta;

    // cout<< eta <<" "<< phi <<" "<< cc <<endl; 
    
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

    int iphi = phi;
    int ieta = eta;
    
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



void readInterCalibBarre_Phiv1new(){
  
  char *file_input = new char[500];
  sprintf(file_input,"phisym_EcalIntercalib-160M-corrected.txt");
  
  sprintf(file_input,"EcalIntercalibConstants_Run2010A_corr.xml.txt");
  
  cout<<"READING fro inter-calib file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  float c1; 
  for(int ieta=-85; ieta <= 85; ieta++){
    if( ieta==0) continue; 
    for(int iphi=1; iphi <= 360; iphi++){
      inputcc >> c1; 
      
      int eta = ieta; 
      int phi = iphi; 
      convxtalid(phi,eta);
      
      if( c1 <=0){
	c1 =1; 
      }
      
      if( c1 == 1 ){
	c1 = - 1; 
      }
      if( c1 <0.6 || c1 >2.5){
	c1 = -1; 
      }
      
      
      CphiCorr[eta+85][phi] = c1; 
    }
  }
    
}


void readInterCalibBarre_Phi2011a(){
  
  char *file_input = new char[500];

  
  sprintf(file_input,"EcalIntercalib-2011A-corrected.xmlv1");
  
  cout<<"READING fro inter-calib file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  float c1; 
  for(int ieta=-85; ieta <= 85; ieta++){
    if( ieta==0) continue; 
    for(int iphi=1; iphi <= 360; iphi++){
      inputcc >> c1; 
      
      int eta = ieta; 
      int phi = iphi; 
      convxtalid(phi,eta);
      
      if( c1 <=0){
	c1 =1; 
      }
      
      if( c1 == 1 ){
	c1 = - 1; 
      }
      if( c1 <0.6 || c1 >2.5){
	c1 = -1; 
      }
      
      
      CphiCorr[eta+85][phi] = c1; 
    }
  }
  


}




void readInterCalibBarre_Phiv1(){
    
  char *file_input = new char[500];
  sprintf(file_input,"phisym_EcalIntercalib-160M-corrected.txt");
  
  
  cout<<"READING fro inter-calib file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  int eta; 
  int phi; 
  float c = 1; 
  float c1 = 1; 
  
  int n = 0; 

  eta  = -85; 
  phi = 0; 
  
  float mean = 0; 
  
  int cureta = -85;
  
  int ngood = 0; 
  
  
  
  while (inputcc.good()){

    inputcc >> c1; 

    if(c1==1) c1 = -1; 
    
    
    if(  n% 360 ==0 && n>0 ) {
      eta++; 
      if(eta==0) eta = 1;  
    }
    
    if( eta == cureta){
      phi ++; 
    }else{
      cureta = eta; 
      phi = 1; 
      
    }
    
    
    n++; 
    
    
    int iSM = (phi-1)/20+1;
    if( eta<0) iSM += 18;
        
    
    int ieta = eta; 
    int iphi = phi; 
    
    
    convxtalid(iphi,ieta);
    
    ///txtout<<eta<<" "<<phi<<endl; 
    
    CphiCorr[ieta+85][iphi] = c1;  
    
    
    if(c1 >0 ){
      mean += c1; 
      ngood ++; 
    }
    
    
    if(n >= 360*170) break; 
    
  }
  
  
  mean /=  ngood; 
  cout<<n<<" phi  constants read. "<< mean <<endl; 
  
  
//   for(int j=0; j< 170; j++){
//     mean_ieta[j] /= nmean_ieta[j];
//     txtout<<"mean_ieta_phisym "<< j<<" "<< mean_ieta[j] <<" "<< nmean_ieta[j] <<endl; 
//   }
  
  
//   for(int j=0; j< 170; j++){
//     for(int k=0; k< 360; k++){
      
//       if(  CphiCorr[j][k] >0 ){
// 	CphiCorr[j][k] /=  mean_ieta[j] ; 
	
//       }
//       Cphi[j][k] = 	CphiCorr[j][k] *  iCalConst_data[j][k];
//     }
//   }
  
//   NormSMScaleToUnit(CphiCorr);
  
}


///read initi calibration constants for each crystal
void readInterCalibBarre_BSv1(){
  
  
  char *file_input = new char[500];
  sprintf(file_input,"IC_barrel_splash09_AllRuns.txt");
  
  
  cout<<"READING fro inter-calib bsfile "<<file_input<<endl;
  txtout<<"READING fro inter-calib bsfile "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  int eta; 
  int phi; 
  float c = 1; 
  float c1 = 1; 
  
  int n = 0; 
  int ngood = 0; 
  

  float mean = 0; 
   
  
  
  
  float mean_ieta[170] = {0}; 
  int nmean_ieta[170] = {0};
  

  while (inputcc.good()){

    inputcc >> eta >> phi >>c1; 
    
    //  cout<<phi<<" "<<eta<<" "<<c<<endl;
    
    convxtalid(phi,eta);
    
    eta += 85; 
    
    if( c1 == 1){
      c1 = -1; 
    }
    
    if( c1 <0.6 || c1 >2.5){
      c1 = -1; 
    }
    
    
    ///iCalConstPi0[eta][phi] = c1 ; 
    
    CBSCorr[eta][phi] = c1;  
    
    if(c1 >0 ){
      mean += c1; 
      ngood ++;
    }

    if(c1 >0 ){
      mean_ieta[eta] += c1; 
      nmean_ieta[eta] ++; 
    }
    
    n++; 
    
    /// hh_imcEB->Fill(c);
    if(n >= 360*170) break; 
    
  }
  
  
  
  mean /= ngood; 
  
  
    
  cout<<n<<" bs  constants read. "<< mean << endl;
  
    
  
//   for(int j=0; j< 170; j++){
//     for(int k=0; k< 360; k++){
      
//       if(  CBSCorr[j][k] >0 ){
// 	CBSCorr[j][k] /=  mean_ieta[j];
//       }
      
//       CBS[j][k] = CBSCorr[j][k] * iCalConst_data[j][k];
      
      
//     }
//   }
//   NormSMScaleToUnit(CBSCorr);
  
  
  
}





///read initi calibration constants for each crystal
void readInterCalibBarre_GR09_V8(){


  
  char *file_input = new char[500];
  //file_input ="/data/yangyong/data/p
  file_input = "/afs/cern.ch/cms/cit/yongy/interCalib_GR09_R_V6A_barrel.txt";
  
  ///sprintf(file_input,"/uscms_data/d2/yongy/calib/backup/interCalib_STARTUP3X_V8P_barrel_mc.txt");  ///now the same used in data
  
  //sprintf(file_input,"/uscms/home/yongy/work/interCalib_STARTUP3X_V8P_barrel.txt");
  
  //  sprintf(file_input,"/uscms/home/yongy/work/interCalib_GR_R_36X_V11A_Barrel.txt");

  cout<<"READING fro inter-calib file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  int eta; 
  int phi; 
  float c = 1; 
  float c1 = 1; 
  
  int n = 0; 
  
  double mean = 0; 
  

  while (inputcc.good()){

    inputcc >> eta >> phi >>c1; 
    
    //  cout<<phi<<" "<<eta<<" "<<c<<endl;
    
    convxtalid(phi,eta);
    
    eta += 85; 
    
    ///iCalConst[eta][phi] = c; 
    
    //iCalConst_data[eta][phi] = c1 * 0.038012/ 0.038;  ///Now this is exactly the same as the inter-calibration used in data. 
    
    interCalib_preCalib[eta][phi] = c1 ; /// now OK iN MC startupv8p. 
    
    
    n++; 
    //  cout<<phi<<" "<<eta<<" "<<c<<endl;
    mean += c1; 
    
    /// hh_imcEB->Fill(c);
    if(n >= 360*170) break; 
    
  }
  
 
  mean /= 61200; 
    
  cout<<n<<" constants read. "<< mean <<endl; 
  
  
}








///read initi calibration constants for each crystal
void readInterCalibEndcap_GR09_V8(){
  char *file_input = new char[500];
  //file_input ="/data/yangyong/data/p
  ///file_input = "/afs/cern.ch/cms/cit/yongy/interCalib_GR09_R_V6A_barrel.txt";
  
  file_input = "/afs/cern.ch/cms/cit/yongy/interCalib_GR09_H_V6OFF_endcap.txt";
    
  cout<<"READING fro inter-calib file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  int eta; 
  int phi; 
  int iz; 
  float c = 1; 
  int n = 0; 
  
  
  while (inputcc.good()){
    
    inputcc >> iz >> eta >> phi  >> c; 
    
    int izz = iz < 0? 0:1; 

    interCalibEndcap_preCalib[izz][eta][phi] = c; 
    validRecHitEndCap[izz][eta][phi] = 1; 
    
    n ++; 
    
  }
  
  cout<<n<<" constants read. "<<endl;
  
}


void getInterCalibEndcapv1(const char *inputfile,float C[2][101][101]){

  
    
  ifstream inputcc(inputfile);
  if (inputcc.fail()) {
    cout<<"error "<< inputfile<<" can not be opened."<<endl; 
    exit(1);
  }
  int eta; 
  int phi; 
  int iz; 
  float c = 1; 
  int n = 0; 
  
  double mean = 0; 
  
  while (inputcc.good()){

    inputcc >> iz >> eta >> phi  >> c; 
    
    int izz = iz < 0? 0:1; 

    C[izz][eta][phi] = c; 
    mean += c; 
    
    n ++; 

    if( n>= 14648){
      break; 
    }
    
  }
  
  cout<<n<<" constants read. "<< mean /n <<endl;
  
}



void getInterCalibEndcap(const char *inputfile, float C[2][101][101]){
  
  
  ifstream inputcc(inputfile);
  if (inputcc.fail()) {
    cout<<"error "<< inputfile<<" can not be opened."<<endl; 
    exit(1);
  }
  
  int ix,iy,iz; 
  float cc; 

  float ccMeaniring[50] ;
  int nccMeaniring[50];
  

//  for(int j=0; j<2;j++){
  for(int k=0; k< kEndcEtaRings; k++){
    ccMeaniring[k] = 0; 
    nccMeaniring[k] = 0; 
  }
  //  }
  
  int nread = 0; 
  while( inputcc.good()){
    inputcc >> iz >> ix >> iy >> cc; 
    C[iz][ix][iy] = cc; 
    int iring = iRingEndCap(2*iz-1,ix,iy); ///input -1/1, 
    if( cc > 0){
      ccMeaniring[iring] += cc; 
      nccMeaniring[iring] += 1; 
    }
    
    nread ++; 
    if( iz==1 && ix == 100 && iy == 60) break; 
  }
  
  cout<<"nb of crystals endcap " << nread <<endl; 
  
  
  //for(int j=0; j<2;j++){
  for(int k=0; k< kEndcEtaRings; k++){
    ccMeaniring[k] /= nccMeaniring[k];
    //  }
  }
  
  double meanall = 0; 
  int nall = 0; 
  for(int iz=0; iz<2; iz++){
    for(int j=0; j<101; j++){
      for(int k=0; k< 101; k++){
        if( validRecHitEndCap[iz][j][k] <1) continue;
        int iring = iRingEndCap(2*iz-1,j,k); ///input -1/1,                      
        if( C[iz][j][k] >0) {
          C[iz][j][k] /= ccMeaniring[iring];
	  meanall += C[iz][j][k];
	  nall ++; 
	}
	
      }
    }
  }
  
  meanall /= nall; 
  cout<<" corrfactoriZiXiY read " << meanall <<" "<< nall <<endl;
  
  
}



void NormCrystalDeadFlagEndcap_v1(float C[2][101][101],int ndeadflag[2][101][101]){
  
  float mean_deadflag[20] = {0};
  float nmean_deadflag[20] = {0};
		       
  int ndeadCrystals = 0;

  
  for(int iz =0; iz < 2; iz++){
    for(int j=0; j< 101; j++){
      for(int k=0; k<101; k++){
	
	if( validRecHitEndCap[iz][j][k] ==0) continue; 
	if( C[iz][j][k] <0){
	  continue; 
	}
	int ndead = ndeadflag[iz][j][k];
	if( ndead >=0){
	  mean_deadflag[ndead] += C[iz][j][k];
	  nmean_deadflag[ndead] ++; 
	}
	
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
  
  
  for(int iz =0; iz < 2; iz++){
    for(int j=0; j< 101; j++){
      for(int k=0; k<101; k++){
	
	if( validRecHitEndCap[iz][j][k] ==0) continue; 
	if( C[iz][j][k] <0){
	  continue; 
	}
	
	int ndead = ndeadflag[iz][j][k];
	if( ndead >=0){
	  C[iz][j][k] /= mean_deadflag[ndead];
	}
	
      }
    }
  }
  
  
  
}






void NormCrystalDeadFlagEndcap(float C[2][101][101]){
  
  float mean_deadflag[20] = {0};
  float nmean_deadflag[20] = {0};
		       
  int ndeadCrystals = 0;

  
  for(int iz =0; iz < 2; iz++){
    for(int j=0; j< 101; j++){
      for(int k=0; k<101; k++){
	
	if( validRecHitEndCap[iz][j][k] ==0) continue; 
	if( C[iz][j][k] <0){
	  continue; 
	}
	int ndeadflag = ndeadflag_endcap[iz][j][k];
	if( ndeadflag >=0){
	  mean_deadflag[ndeadflag] += C[iz][j][k];
	  nmean_deadflag[ndeadflag] ++; 
	}
	
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
  
  
  for(int iz =0; iz < 2; iz++){
    for(int j=0; j< 101; j++){
      for(int k=0; k<101; k++){
	
	if( validRecHitEndCap[iz][j][k] ==0) continue; 
	if( C[iz][j][k] <0){
	  continue; 
	}
	
	int ndeadflag = ndeadflag_endcap[iz][j][k];
	if( ndeadflag >=0){
	  C[iz][j][k] /= mean_deadflag[ndeadflag];
	}
	
      }
    }
  }
  
  
  
}



void scaleMeanToUnitEndcap(float C[2][101][101]){
  
  
  double mean = 0; 
  int n = 0; 
  
  
  for(int iz =0; iz < 2; iz++){
    for(int j=0; j< 101; j++){
      for(int k=0; k<101; k++){
	if( validRecHitEndCap[iz][j][k] ==0) continue; 
	if( C[iz][j][k] <0){
	  continue; 
	}
	mean += C[iz][j][k]; 
	n ++; 
      }
    }
  }
    
  mean /= n;
  cout<<"mean " << mean<<" "<<n<<endl; 
   
  
  for(int iz =0; iz < 2; iz++){
    for(int j=0; j< 101; j++){
      for(int k=0; k<101; k++){
	if( validRecHitEndCap[iz][j][k] ==0) continue; 
	if( C[iz][j][k] <0){
	  continue; 
	}
	
	C[iz][j][k] /= mean; 
      }
    }
  }
  
}
