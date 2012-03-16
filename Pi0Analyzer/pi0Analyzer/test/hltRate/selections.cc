int ntrigl1Warining =0;
int ntrigWarining =0;

int testL1prescale;

vector<unsigned short> *l1bitFired_prescale1;
vector<unsigned short> *l1bitFired_prescale2;
vector<unsigned short> *l1bitFired_prescale3;


bool isL1PathFired(string pathname){

  vector<unsigned int> nfound ;
  for(unsigned int j=0; j< l1algoName->size() ;j++){
    string l1 = l1algoName->at(j);
    if( l1 == pathname){
      nfound.push_back(j);
    }
  }
  if( nfound.size() ==0 ){

    if( ntrigl1Warining < 10){
      cout<<"not found any L1 path with name return false ! " << pathname.c_str()<<" "<<runNumber <<" "<<lumiBlock<<endl; 
      ntrigl1Warining ++;
    }
    return false;
  }
  if(  nfound.size() > 1){
    cout<<"found more than one path with the same name pls check! " << pathname.c_str()<<endl;
    for(int n=0; n< int(nfound.size()); n++){
      cout<< l1algoName->at(nfound[n]).c_str()<<endl;
    }
    exit(1);
  }
  
  vector<unsigned short>::const_iterator b = find(l1bitFired->begin(),l1bitFired->end(),nfound[0]);
  if( b != l1bitFired->end()){
    return true;
  }else return false;
    
}


bool isL1PathFired_prescale(string pathname, int prescale){

  vector<unsigned int> nfound ;
  for(unsigned int j=0; j< l1algoName->size() ;j++){
    string l1 = l1algoName->at(j);
    if( l1 == pathname){
      nfound.push_back(j);
    }
  }
  
  if( nfound.size() ==0 ){
    
    if( ntrigl1Warining < 10){
      cout<<"not found any L1 path with name prescale return false ! " << pathname.c_str()<<" "<<runNumber <<" "<<lumiBlock<<endl; 
      ntrigl1Warining ++;
    }
    return false;
    
  }
  if(  nfound.size() > 1){
    cout<<"found more than one path with the same name pls check! " << pathname.c_str()<<endl;
    for(int n=0; n< int(nfound.size()); n++){
      cout<< l1algoName->at(nfound[n]).c_str()<<endl;
    }
    exit(1);
  }
  

  if( prescale ==1){
    vector<unsigned short>::const_iterator b = find(l1bitFired_prescale1->begin(),l1bitFired_prescale1->end(),nfound[0]);
    if( b != l1bitFired_prescale1->end()){
      return true;
    }else return false;
  }else if(  prescale ==2){
    vector<unsigned short>::const_iterator b = find(l1bitFired_prescale2->begin(),l1bitFired_prescale2->end(),nfound[0]);
    if( b != l1bitFired_prescale2->end()){
      return true;
    }else return false;
  }else{
    cout<<"ps NA " << prescale <<endl; 
    exit(1);
  }
  
}


////check if one particular trigger path is fired , input the HLT name 
bool isTriggerPathFiredv1(string pathname){
  
  
  vector<unsigned int> nfound ; 
  for(unsigned int j=0; j< hlt_pathName->size() ;j++){
    string hlt = hlt_pathName->at(j);
    if( hlt.find(pathname)!=string::npos){
      nfound.push_back(j);
    }
  }
  
  if( nfound.size() ==0 ){

    if( ntrigWarining < 10){
      cout<<"not found any HLT path with name return true ! " << pathname.c_str()<<endl; 
      ntrigWarining ++; 
    }
    return true;
    
  }
  if(  nfound.size() > 1){
    cout<<"found more than one path with the same name pls check! " << pathname.c_str()<<endl; 
    exit(1);
  }
  
  vector<unsigned short>::const_iterator b = find(hlt_bitFired->begin(),hlt_bitFired->end(),nfound[0]);
  if( b != hlt_bitFired->end()){
    return true; 
  }else return false; 
  
  
}




TLorentzVector photonp4barrel(int j){
  
  
  TVector3 caloPosition(x3x3ClusEB[j],y3x3ClusEB[j],z3x3ClusEB[j]);
  TVector3 direction = caloPosition; 
  TVector3 p = direction.Unit() * e3x3ClusEB[j];
  TLorentzVector p4(p.x(),p.y(),p.z(), e3x3ClusEB[j]);
  return p4;
  
}

TLorentzVector photonp4endcap(int j){
  
  TVector3 caloPosition(x3x3ClusEE[j],y3x3ClusEE[j],z3x3ClusEE[j]);
  TVector3 direction = caloPosition; 
  TVector3 p = direction.Unit() * e3x3ClusEE[j];
  TLorentzVector p4(p.x(),p.y(),p.z(), e3x3ClusEE[j]);
  return p4;
  
}

float isolationEB(float eta, float phi, int j1, int j2, int pizEta){
  
  
  float sum = 0; 
  
  float detaCut = 0.05; 
  float drCut = 0.2; 
  
  if( pizEta ==2){
    detaCut = 0.1; 
    drCut = 0.3;
  }
  
  for(int j=0; j< n3x3ClusEB; j++){
    if( j== j1 || j== j2) continue; 
    TLorentzVector phtp4 = photonp4barrel(j);
    if( phtp4.Pt() < 0.5 ) continue; 
    
    float deta = fabs( eta - phtp4.Eta()); 
    float dr = GetDeltaR(eta,phtp4.Eta(), phi, phtp4.Phi());
    
    if(deta< detaCut  && dr < drCut ){
      sum += phtp4.Pt(); 
    }
  }
  
  return sum; 
  
}




float isolationEE(float eta, float phi, int j1, int j2, int pizEta){
  
  
  float sum = 0; 
  
  float detaCut = 0.05; 
  float drCut = 0.2; 
  
  if( pizEta ==2){
    detaCut = 0.1; 
    drCut = 0.3;
  }
  
  for(int j=0; j< n3x3ClusEE; j++){
    if( j== j1 || j== j2) continue; 
    TLorentzVector phtp4 = photonp4endcap(j);
    if( phtp4.Pt() < 0.5 ) continue; 
    
    float deta = fabs( eta - phtp4.Eta()); 
    float dr = GetDeltaR(eta,phtp4.Eta(), phi, phtp4.Phi());
    
    if(deta< detaCut  && dr < drCut ){
      sum += phtp4.Pt(); 
    }
  }
  
  return sum; 
  
}




vector<float> selection_EE_piz(){

  //   float ptminCut = 1.3; 
  //   float ptpairCut = 2.6;
  
  float s4s9minCut = 0.90; 
  float isoCut = 0.5; 
  float selePtGammaEndCap_region1_ = 0.6; 
  float selePtPairEndCap_region1_ = 2.5; 
  float selePtGammaEndCap_region2_ = 0.6; 
  float selePtPairEndCap_region2_ = 2.5; 

  float selePtGammaEndCap_region3_ = 0.5; 

  //float selePtPairEndCap_region3_ = 1.0; 
  float selePtPairEndCap_region3_ = 99.0; 

  float selePtPairMaxEndCap_region3_ = 2.5; 
  
  
  vector<float> mselected; 
  
  float region1_EndCap_ = 2.0; 
  float region2_EndCap_ = 2.5; 
  
  
  for(int j=0; j< n3x3ClusEE; j++){

    if( s4s93x3ClusEE[j]  < s4s9minCut) continue; 
    TLorentzVector p4_pht1 = photonp4endcap(j);

    for(int k=j+1; k< n3x3ClusEE;  k++){
      
      if( s4s93x3ClusEE[k]  < s4s9minCut) continue; 
      TLorentzVector p4_pht2 = photonp4endcap(k);
      
      TLorentzVector p4_pair = p4_pht1 + p4_pht2 ; 
      ptpair = p4_pair.Pt();
      
      ptmin = p4_pht1.Pt() <  p4_pht2.Pt()  ?  p4_pht1.Pt() :  p4_pht2.Pt(); 
            
      etapair = fabs(p4_pair.Eta());
      if(etapair <= region1_EndCap_){
	if(ptmin < selePtGammaEndCap_region1_ || ptpair < selePtPairEndCap_region1_) continue; 
      }else if( etapair <= region2_EndCap_){
	if(ptmin < selePtGammaEndCap_region2_ || ptpair < selePtPairEndCap_region2_) continue;
      }else{
	  if(ptmin < selePtGammaEndCap_region3_ || ptpair < selePtPairEndCap_region3_) continue;
	  if(ptpair > selePtPairMaxEndCap_region3_ ) continue; 
      }
      
      float iso = isolationEE(p4_pair.Eta(), p4_pair.Phi(), j,k, 1);
      if( iso / ptpair > isoCut ) continue; 
      
       mpair = p4_pair.M();
      
      //cout<<"check mpair " << mpair <<" "<< s4s93x3ClusEE[j]<<" "<< s4s93x3ClusEE[k]<<" "<< iso/ptpair<<endl; 

      mselected.push_back(mpair);
      
      
    }
  }
  
  return mselected; 
  
  
}




vector<float> selection_EE_eta(){
  
  
  float s4s9minCut = 0.90; 
  float s9s25minCut =0.85;
  float isoCut = 0.5; 
  float selePtGammaEndCap_region1_ = 1.0; 
  float selePtPairEndCap_region1_ = 3.0; 
  float selePtGammaEndCap_region2_ = 1.0; 
  float selePtPairEndCap_region2_ = 3.0; 

  float selePtGammaEndCap_region3_ = 0.7; 
  float selePtPairEndCap_region3_ = 3.0; 
  float selePtPairMaxEndCap_region3_ = 99999;
    
  vector<float> mselected; 
  
  float region1_EndCap_ = 2.0; 
  float region2_EndCap_ = 2.5; 
  
  
  for(int j=0; j< n3x3ClusEE; j++){

    if( s4s93x3ClusEE[j]  < s4s9minCut) continue; 
    if( s9s253x3ClusEE[j]  < s9s25minCut) continue; 
    
    TLorentzVector p4_pht1 = photonp4endcap(j);

    for(int k=j+1; k< n3x3ClusEE;  k++){
      
      if( s4s93x3ClusEE[k]  < s4s9minCut) continue; 
      if( s9s253x3ClusEE[k]  < s9s25minCut) continue; 
      
      TLorentzVector p4_pht2 = photonp4endcap(k);
      
      TLorentzVector p4_pair = p4_pht1 + p4_pht2 ; 
      ptpair = p4_pair.Pt();
      
      etapair = fabs(p4_pair.Eta());

      ptmin = p4_pht1.Pt() <  p4_pht2.Pt()  ?  p4_pht1.Pt() :  p4_pht2.Pt(); 
      

      if(etapair <= region1_EndCap_){
	if(ptmin < selePtGammaEndCap_region1_ || ptpair < selePtPairEndCap_region1_) continue; 
      }else if( etapair <= region2_EndCap_){
	if(ptmin < selePtGammaEndCap_region2_ || ptpair < selePtPairEndCap_region2_) continue;
      }else{
	  if(ptmin < selePtGammaEndCap_region3_ || ptpair < selePtPairEndCap_region3_) continue;
	  if(ptpair > selePtPairMaxEndCap_region3_ ) continue; 
      }
      
      float iso = isolationEE(p4_pair.Eta(), p4_pair.Phi(), j,k, 2);
      if( iso / ptpair > isoCut ) continue; 
      
      mpair = p4_pair.M();
      
      mselected.push_back(mpair);
      
    }
  }
  
  return mselected; 
  
  
}




vector<float> selection_EB_piz(){
  
  float ptminCut = 1.3; 
  float ptpairCut = 2.6;
  float s4s9minCut = 0.83; 
  float isoCut = 0.5; 
  
  vector<float> mselected; 
  
  
  for(int j=0; j< n3x3ClusEB; j++){

    if( s4s93x3ClusEB[j]  < s4s9minCut) continue; 

    TLorentzVector p4_pht1 = photonp4barrel(j);
    if( p4_pht1.Pt() < ptminCut ) continue; 
    

    for(int k=j+1; k< n3x3ClusEB;  k++){
      
      if( s4s93x3ClusEB[k]  < s4s9minCut) continue; 
      TLorentzVector p4_pht2 = photonp4barrel(k);
      if( p4_pht2.Pt() < ptminCut ) continue; 
      

      TLorentzVector p4_pair = p4_pht1 + p4_pht2 ; 
      ptpair = p4_pair.Pt();
      if( ptpair < ptpairCut ) continue; 
      
      float iso = isolationEB(p4_pair.Eta(), p4_pair.Phi(), j,k, 1);
      if( iso / ptpair > isoCut ) continue; 
      
      mpair = p4_pair.M();
  
      //cout<<"mpairebsel " << mpair << " "<< iso/ptpair <<" "<< s4s93x3ClusEB[j]<<" "<< s4s93x3ClusEB[k]<<endl; 
    
      mselected.push_back(mpair);
      
      
    }
  }
  
  return mselected; 
  
  
}



vector<float> selection_EB_eta(){
  
  float ptminCut = 1.2; 
  float ptpairCut = 4.0;
  float s4s9minCut = 0.87; 
  float s9s25minCut = 0.8; 
  float isoCut = 0.5; 
  
  
  vector<int> indClusPi0Candidates;  ///those clusters identified as pi0s                            

  for(int j=0; j< n3x3ClusEB; j++){
    TLorentzVector p4_pht1 = photonp4barrel(j);
    for(int k=j+1; k< n3x3ClusEB;  k++){
      TLorentzVector p4_pht2 = photonp4barrel(k);
      TLorentzVector p4_pair = p4_pht1 + p4_pht2 ; 
       mpair = p4_pair.M();
      
      if( mpair >0.084 && mpair < 0.156){
	int indtmp[2] = {j,k};
	 for(int n=0; n<2; n++){
	   std::vector<int>::iterator it = find(indClusPi0Candidates.begin(),indClusPi0Candidates.end(),indtmp[n]);
	   if( it == indClusPi0Candidates.end()) {
	     indClusPi0Candidates.push_back(indtmp[n]);
	   }
	 }
      }
      
    }
  }
  
  
  vector<float> mselected; 
  
  for(int j=0; j< n3x3ClusEB; j++){

    
    if( s4s93x3ClusEB[j]  < s4s9minCut) continue; 
    if( s9s253x3ClusEB[j]  < s9s25minCut) continue; 
    TLorentzVector p4_pht1 = photonp4barrel(j);
    
    
    
    std::vector<int>::iterator it = find(indClusPi0Candidates.begin(),indClusPi0Candidates.end(),j);
    if( it != indClusPi0Candidates.end()) continue;
    
    if( p4_pht1.Pt() < ptminCut ) continue; 

    for(int k=j+1; k< n3x3ClusEB;  k++){

      if( s4s93x3ClusEB[k]  < s4s9minCut) continue; 
      if( s9s253x3ClusEB[k]  < s9s25minCut) continue; 
      
      
      it = find(indClusPi0Candidates.begin(),indClusPi0Candidates.end(),k);
      if( it != indClusPi0Candidates.end()) continue;

      TLorentzVector p4_pht2 = photonp4barrel(k);
      if( p4_pht2.Pt() < ptminCut ) continue; 
      
      
      TLorentzVector p4_pair = p4_pht1 + p4_pht2 ; 
      ptpair = p4_pair.Pt();
      if( ptpair < ptpairCut ) continue; 
      
      float iso = isolationEB(p4_pair.Eta(), p4_pair.Phi(), j,k, 2);
      if( iso / ptpair > isoCut ) continue; 
      
      mpair = p4_pair.M();
      
      mselected.push_back(mpair);
      
      
    }
  }
  
  return mselected; 
  
  
}
