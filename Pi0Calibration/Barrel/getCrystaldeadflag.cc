
void getCrystaldeadflagBarrel(){
  
  ///TString file_input = "/uscms/home/yongy/work/crystal_deadflag_eb_dflag64.txt";
  TString file_input = TString(workingDirectory) + TString(Form("/crystal_deadflag_eb_dflag%d.txt",dataflag));
  
  cout<<"READING for dead cystal barrel file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  int ieta; 
  int iphi; 
  int flag; 
  
  int n = 0; 
  int ndead = 0; 

  int ncorner;
  int nside; 

  int ndeadflag; 

  //format
  /// ieta, iphi ( after convxtalid, ieta+85 ) , flag ( 0/1) ncorner, nside, ndeadflag
    
  while (inputcc.good()){

    inputcc >> ieta >> iphi >> flag >> ncorner >> nside  >> ndeadflag; 
    flag_ietaiphi[ieta][iphi] = flag; 
    ndeadcorner_ietaiphi[ieta][iphi] = ncorner; 
    ndeadside_ietaiphi[ieta][iphi] = nside;
    
    ndeadflag_ietaiphi[ieta][iphi] = ndeadflag; 
    
    
    n ++; 
    if( flag ==1) ndead ++; 
    
    if( n>= 61200) break; 
  }
  
  cout<<"n dead xtals "<< ndead <<endl; 
  
  
}



void getCorrFactorDead(){
  
  for(int j=0; j< 20; j++){
    corrfactorDead[j] = 1; 
  }
  
}

void getCorrFactorDead_Barrel_CorrectionOnClusterEnergy(){
  
  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step3.iter11.root",dataflag,pizEta));
  cout<<" Reading dead crystal correction file " << filename<<endl; 
  TFile *ff =new TFile(filename,"read"); 

  if( ff->Get("hh_corr_deadv1") == NULL){
    cout<<"no hh_corr_deadv1 " <<endl; 
    exit(1);
  }
  
  TH1F *hh_corr_dead = (TH1F*)ff->Get("hh_corr_deadv1");
  
  for(int b=1; b<= hh_corr_dead->GetNbinsX(); b++){
    corrfactorDead[b-1] = hh_corr_dead->GetBinContent(b); 
    if( b==1 || corrfactorDead[b-1]!=1) cout<<"corr_dead"<<b-1<<" "<< corrfactorDead[b-1]<<endl; 
  }
  cout<<"correction of dead crystal read.."<<endl; 
    
  ///read the final mean and width for L3 method
  TH1F *hhsigma = (TH1F*)ff->Get("hh_res_ieta_2");
  int b1; 
  int b2; 
  int b3; 
  int b4;
  int nbinforsigma;
  if(doPizEta==2){
    ///for eta, use the mean value of ieta30 to40
    b1 = hhsigma->GetXaxis()->FindBin(30);
    b2 = hhsigma->GetXaxis()->FindBin(40);
    b3 = hhsigma->GetXaxis()->FindBin(-40);
    b4 = hhsigma->GetXaxis()->FindBin(-30);
    nbinforsigma = 2*11; 
  }else{
    b1 = hhsigma->GetXaxis()->FindBin(10);
    b2 = hhsigma->GetXaxis()->FindBin(20);
    b3 = hhsigma->GetXaxis()->FindBin(-20);
    b4 = hhsigma->GetXaxis()->FindBin(-10);
    nbinforsigma = 2*11;
  }
  sigmaMass = 0;
  for(int b=b1; b<= b2; b++){
    sigmaMass += hhsigma->GetBinContent(b);
  }
  for(int b=b3; b<= b4; b++){
    sigmaMass += hhsigma->GetBinContent(b);
  }
  sigmaMass /= nbinforsigma;

  float mean1 = 0;
  int nbinformean = 0;
  float wtMeanS = 0; 
  float wtMeanD = 0; 
  
  TH1F *hhpeak =  (TH1F*)ff->Get("hh_res_ietaSM_36_0"); 
  if(hhpeak == NULL){
    cout<<"error getCorrFactorDead_Barrel_CorrectionOnClusterEnergy null hist "  <<endl; 
    exit(1);
  }
  for(int b=1; b<= hhpeak->GetNbinsX();b++){
    if( hhpeak->GetBinContent(b)>0.05 && hhpeak->GetBinContent(b) < 1){ //mainly for check 0 
      mean1 += hhpeak->GetBinContent(b);
      float mm = hhpeak->GetBinContent(b); 
      float ms = hhpeak->GetBinError(b);
      if( ms / mm < 0.05){ //in case bad fit 
	wtMeanS +=  mm/ (ms* ms); 
	wtMeanD += 1/ (ms*ms);
      }
      nbinformean ++; 
    }
  }
  mean1 /= nbinformean; 
  meanMass = wtMeanS / wtMeanD; 
  
  
}
