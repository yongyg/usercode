

void getcorrFactorEtaOfEachStep(){
  

  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step1.iter%d.root",dataflag,pizEta,iter-1));
  cout<<"READING from a eta ring correction file "<<filename<<endl; 
  TFile *ff1 = new TFile(filename,"read"); 
  
  float mean_etaRing[2] = {0};
  
  for(int iz=0; iz<2; iz++){
    filename = TString(Form("hh_corr_etaRing_%d",iz));
    TH1F *hhtmp = (TH1F*)ff1->Get(filename); 
    if(hhtmp==0){
      cout<<filename<<" NA! " <<endl; 
      exit(1);
    }
    
    for(int j=0; j< kEndcEtaRings; j++){
      corrfactorEtaRings[iz][j] = hhtmp->GetBinContent(j+1); 
      if(j<= nMaxRingIC) mean_etaRing[iz] += corrfactorEtaRings[iz][j];
      if(j<2) cout<<"etacorrRing: "<< iz <<" "<<j<<" "<< corrfactorEtaRings[iz][j]<<endl;  
    }
  }
  
  for(int iz=0; iz<2; iz++){
    mean_etaRing[iz] /= (nMaxRingIC+1);
  }
  
  
  cout<<"mean_etaRing_corr: "<< mean_etaRing[0] <<" "<< mean_etaRing[1] <<endl;
}



void getcorrFactorEta(){
  

  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step1.iter10.root",dataflag,pizEta));
  cout<<"READING from a eta ring correction file "<<filename<<endl;
  TFile *ff1 = new TFile(filename,"read");
  float mean_etaRing[2] = {0};
  
  ///use ietaring from 10 - 30
  float peakee[2]={0};
  float widthee[2]={0};
  
  for(int iz=0; iz<2; iz++){
    filename = TString(Form("hh_corr_etaRing_%d",iz));
    TH1F *hhtmp = (TH1F*)ff1->Get(filename); 
    
    if(hhtmp==0 ){
      cout<<filename<<" NA! " <<endl;
      exit(1);
    }
    
    filename = TString(Form("hh_res_ieta_%d_2",iz));
    TH1F *hhtmp1 = (TH1F*)ff1->Get(filename);
    if(hhtmp1==0 ){
      cout<<filename<<" NA! " <<endl;
      exit(1);
    }
    
    filename = TString(Form("hh_res_ieta_%d_0",iz));
    TH1F *hhtmp2 = (TH1F*)ff1->Get(filename);
    if(hhtmp2==0 ){
      cout<<filename<<" NA! " <<endl;
      exit(1);
    }
    
    for(int j=0; j< kEndcEtaRings; j++){
      corrfactorEtaRings[iz][j] = hhtmp->GetBinContent(j+1); 
      if( j>=5 && j<=15){ ////good rings.
	peakee[iz] += hhtmp2->GetBinContent(j+1);
	mean_etaRing[iz] += corrfactorEtaRings[iz][j];
	widthee[iz] += hhtmp1->GetBinContent(j+1);
      }
      if(j<2) cout<<"etacorrRing: "<< iz <<" "<<j<<" "<< corrfactorEtaRings[iz][j]<< " "<< hhtmp1->GetBinContent(j+1) <<endl; 
      
    }
  }
  
  for(int iz=0; iz<2; iz++){
    mean_etaRing[iz] /= 11; 
    mean_side[iz] = peakee[iz]/ 11; 
    sigma_side[iz] = widthee[iz]/11; 
  }
  
  cout<<"mean_etaRing_corr: "<< mean_side[0] <<" "<< mean_side[1] <<"  width " << sigma_side[0]<<" "<< sigma_side[1]<<endl; 
  
}




void getcorrFactorIzIxIyOfEachStep(){

  TString filename = TString(workingDirectory) + TString(Form("/calibres/deriveCalibConst.dflag%d.pe%d.step2.iter%d.txt",dataflag,pizEta,iter-1));
  cout<<"READING from a ix/iy correction file "<<filename<<endl; 
  ifstream inputcc(filename);
  if (inputcc.fail()) {
    cout<<"error "<< filename <<" can not be opened."<<endl; 
    exit(1);
  }
  int ix,iy,iz; 
  float cc; 
  float ccMeaniring[50] ={0};
  int nccMeaniring[50]={0};
  
  while( inputcc.good()){
    inputcc >> iz >> ix >> iy >> cc; 
    corrfactoriZiXiY[iz][ix][iy] = cc; 
    
    int iring = iRingEndCap(2*iz-1,ix,iy); ///input -1/1, 
    if( cc > 0){
      ccMeaniring[iring] += cc; 
      nccMeaniring[iring] += 1; 
    }
    if( iz==1 && ix == 100 && iy == 60) break; 
  }
  for(int k=0; k< kEndcEtaRings; k++){
    ccMeaniring[k] /= nccMeaniring[k];
  }
  double meanall = 0; 
  int nall = 0; 
  for( iz=0; iz<2; iz++){
    for(int j=0; j<101; j++){
      for(int k=0; k< 101; k++){
	if( validRecHitEndCap[iz][j][k] <1) continue;
	int iring = iRingEndCap(2*iz-1,j,k); ///input -1/1,                      
	if( corrfactoriZiXiY[iz][j][k] >0) {
	  corrfactoriZiXiY[iz][j][k] /= ccMeaniring[iring];
	}else{
	  corrfactoriZiXiY[iz][j][k] = 1; 
	}
	meanall += corrfactoriZiXiY[iz][j][k]; 
	nall ++; 
      }
    }
  }
  
  meanall /= nall; 
  cout<<" corrfactoriZiXiY read " << meanall <<" nall "<< nall <<endl;
  
}
