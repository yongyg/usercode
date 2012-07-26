
void getSMcorrFactor(){
  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step1.iter1.root",dataflag,pizEta));
  cout<<"READING from a SM correction file "<<filename<<endl; 
  TFile *fftmp = new TFile(filename,"read"); 
  TH1F *hhtmp =(TH1F*)fftmp->Get("hh_corr_sm");
  if(hhtmp == NULL){
    cout<<" getSMcorrFactor null hist!!! " <<endl; 
    exit(1);
  }
  double meancorr = 0; 
  int nbinsx = hhtmp->GetNbinsX();
  if(nbinsx !=36){
    cout<<"warning! nbinsx sm correction: "<< nbinsx<<endl; 
  }
  for(int b =1; b<= nbinsx; b++){
    meancorr += hhtmp->GetBinContent(b); 
  }
  
  meancorr /= nbinsx; 
  cout<<"meancorr sm correction: "<< meancorr <<endl; 
  
  for(int b =1; b<= nbinsx; b++){
    float tmp = hhtmp->GetBinContent(b)/ meancorr; 
    corrfactorSM[b-1] = tmp; 
    if( b==30 || b==31) cout<<"correction factor SM: "<<b <<"  "<<  tmp <<endl; 
  }
    
}





void getcorrFactorEtaOfEachStep(){
  cout<<"now getting eta correction file." <<endl; 
  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step2.iter%d.root",dataflag,pizEta,iter-1));  

  cout<<"READING from a eta ring correction file "<<filename<<endl; 
  
  TFile *ffceta = new TFile(filename,"read");
  TString histname = "hh_corr_ietaSM_36";
  TH1F *hh_corr_ieta  = (TH1F*)ffceta->Get(histname);
  if( hh_corr_ieta == NULL){
    cout<<"getcorrFactorEtaOfEachStep no histogram ! " <<histname<<endl; 
    exit(1);
  }
  for(int b=1; b<= 85; b++){
    corrfactorEta[b-1] = hh_corr_ieta->GetBinContent(85-b+1);
    corrfactorEta[b-1+85] = hh_corr_ieta->GetBinContent(b);
  }
  
  
  double meancorr = 0;
  for(int j=0; j< 170;j++){
    meancorr += corrfactorEta[j];
    if( ! ( corrfactorEta[j]> 0.5 && corrfactorEta[j] < 1.5 ) ) {
      cout<<"warning pretty large corrfactorEta  " <<j <<" "<< corrfactorEta[j]<<endl; 
      exit(1);
    }
    if( j==0 || j== 84 || j== 85 || j==169 ) cout<<" corrfactorEta["<<j<<"]" << corrfactorEta[j]<<endl;
  }
  meancorr /= 170; 
  cout<<"  meancorr_ieta " << meancorr <<endl; 
  for(int j=0; j< 170;j++){
    corrfactorEta[j] /= meancorr; 
  }
  
}




//this file reads the eta/phi correction factors
/// corrfactorEta[170]; and corrfactorPhi[360];

///this is after all steps done. 

void getcorrFactorEta(){
  
  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step2.iter%d.root",dataflag,pizEta,10));
  cout<<"READING from a final eta ring correction file "<<filename<<endl; 
  
  TFile *ffceta = new TFile(filename,"read");
  TH1F *hh_corr  = (TH1F*)ffceta->Get("hh_corr_ietaSM_36"); ///Abs(ieta);
  if(hh_corr == NULL){
    cout<<" getcorrFactorEta NULL histogram  " <<endl; 
    exit(1);
  }
  for(int b=1; b<= 85; b++){
    corrfactorEta[b-1] = hh_corr->GetBinContent(85-b+1);
    corrfactorEta[b-1+85] = hh_corr->GetBinContent(b);
  }
  
  double meancorr = 0;
  for(int j=0; j< 170;j++){
    meancorr += corrfactorEta[j];
    if( ! ( corrfactorEta[j]> 0.5 && corrfactorEta[j] < 1.5 ) ) {
      cout<<"warning large corrfactorEta " <<j <<" "<< corrfactorEta[j]<<endl; 
      exit(1);
    }
    if( j==0 || j== 84 || j== 85 || j==169 ) cout<<" corrfactorEta["<<j<<"]" << corrfactorEta[j]<<endl;
  }
  meancorr /= 170; 
  
  cout<<"  meancorr_ieta " << meancorr <<endl; 
  for(int j=0; j< 170;j++){
    corrfactorEta[j] /= meancorr; 
  }
  
}


///to derive phi correction 
void getcorrFactorPhiOfEachStep(){
  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step3.iter%d.root",dataflag,pizEta,iter-1));
  
  cout<<"READING from a phi ring correction file "<<filename<<endl; 
  TFile *fpcorr = new TFile(filename,"read");
  if( fpcorr->Get("hh_corr_iphismtb") == NULL){
    cout<<"getcorrFactorPhiOfEachStep no histogram !" <<endl; 
    exit(1);
  }
  TH1F *hh_corr_iphismtb = (TH1F*)fpcorr->Get("hh_corr_iphismtb");
  for(int j=1; j< 360; j++){
    int by = (j-1)%20 +1; 
    corrfactorPhi[j] = hh_corr_iphismtb->GetBinContent(by);
  }
  corrfactorPhi[0] = hh_corr_iphismtb->GetBinContent(20);
  double meancorr = 0;
  for(int j=0; j< 360; j++){
    if( j==0 || j==1 ||  j== 20 || j== 21 || j==40 || j== 359 ) cout<<" corrfactorPhi[" <<j<<"]: " <<corrfactorPhi[j]<<endl; 
    meancorr += corrfactorPhi[j];
  }
  meancorr /= 360; 
  cout<<" corrfactorPhi meancorr " << meancorr <<endl; 
  
  for(int j=0; j< 360; j++){
    corrfactorPhi[j] /= meancorr;     
  }
  
  
}


///to get phi correction after all steps 
void getcorrFactorPhi(){
  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step3.iter%d.root",dataflag,pizEta,10));
  cout<<"READING from a final phi ring correction file "<<filename<<endl; 
  
  TFile *fpcorr = new TFile(filename,"read");
  if(fpcorr->Get("hh_corr_iphismtb") == NULL){
    cout<<"correction of iphi histograms not found.." <<endl; 
    exit(1);
  }
  
  TH1F *hh_corr_iphismtb = (TH1F*)fpcorr->Get("hh_corr_iphismtb");
  for(int j=1; j< 360; j++){
    int by = (j-1)%20 +1; 
    corrfactorPhi[j] = hh_corr_iphismtb->GetBinContent(by);
  }
  corrfactorPhi[0] = hh_corr_iphismtb->GetBinContent(20);
  double meancorr = 0;
  for(int j=0; j< 360; j++){
    if( j==0 || j==1 ||  j== 20 || j== 21 || j==40 || j== 359 ) cout<<" corrfactorPhi[" <<j<<"]: " <<corrfactorPhi[j]<<endl; 
    meancorr += corrfactorPhi[j];
  }
  meancorr /= 360; 
  cout<<" corrfactorPhi meancorr " << meancorr <<endl; 
  
  for(int j=0; j< 360; j++){
    corrfactorPhi[j] /= meancorr;     
  }
  

}
  




/// this is the derived correction to pre calbibrations constants 
void getcorrFactorIetaIphiOfEachStep(){
  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step4.iter%d.txt",dataflag,pizEta,iter-1));
  cout<<"READING from a ieta/iphi correction file "<<filename<<endl; 
  ifstream inputcc(filename);
  if (inputcc.fail()) {
    cout<<"error "<< filename <<" can not be opened."<<endl; 
    exit(1);
  }
  
  int ieta; 
  int iphi; 
  double corr; 
  double meancorr = 0; 
  int nread = 0; 
  int ngood = 0; 
  while( inputcc.good()){
    inputcc >> ieta >> iphi >> corr ; 
    if( (iphi <0 || iphi >=360) || (ieta<0 || ieta >=170 )){
      cout<<"error input corrPhi ieta[0,169), iphi [0,360): "<< ieta <<" "<< iphi <<endl; 
      exit(1);
    }
    corrfactoriEtaiPhi[ieta][iphi]  = corr; 
    
    if( ieta == 0 && iphi == 0 ) cout<<ieta<<" "<< iphi <<" "<< corr <<endl; 
    
    if( corr < 0 ){
      flagiEtaiPhi[ieta][iphi]= -1; 
    }else {
      flagiEtaiPhi[ieta][iphi]= 0; 
      meancorr += corrfactoriEtaiPhi[ieta][iphi];
      ngood ++; 
    }
    
    nread ++; 
    if( nread >= 61200){
      break; 
    }
  }
  
  meancorr /= ngood; 
  
  cout<<"meancorr ieta/iphi correction: "<< meancorr <<" "<< ngood <<endl; 
  
  for(int j=0; j<170; j++){
    for(int k=0; k<360; k++){
      
      if( corrfactoriEtaiPhi[j][k]>0){
	corrfactoriEtaiPhi[j][k] /= meancorr; 
      }else{
	corrfactoriEtaiPhi[j][k] = 1; 
      }
    }
  } 
  
}
