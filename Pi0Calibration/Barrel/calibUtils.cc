
void initialize_barrel(){
  for(int j=0; j<170; j++){
    for(int k=0; k<360; k++){
      corrfactoriEtaiPhi[j][k] = 1;
    }
  }
  for(int j=0; j<36; j++){
    corrfactorSM[j] = 1;
  }
  for(int j=0; j< 170; j++){
    corrfactorEta[j] = 1;
  }
  for(int j=0; j< 360; j++){
    corrfactorPhi[j] = 1;
  }
  for(int j=0; j<20; j++){
    corrfactorDead[j] = 1;
  }
  for(int j=0; j<170; j++){
    for(int k=0; k<360; k++){
      C0[j][k] = 1;  ///first step set to be 1 

      WTSUM[j][k] = 0;
      CORSUM[j][k] = 0;
    }
  }
}



//Input two clusters, compute all variables from that
///with new calibration calconst after each step

void calcPairClusterUpdated(float CalInput[170][360],float res[]){
  
  Float_t msum9[2], meta[2], mphi[2];
  
  
  int ietaxt[9]; //after convxtlid
  int iphixt[9];
  
  float ext[9];
  
  float posi[3];
  
  
  msum9[0] = 0; 
  int ixtal = 0; 

  for(Int_t jxtal=0; jxtal  < nxtClus1; jxtal++){
    int ieta = ietaXtalClus1[jxtal];
    int iphi = iphiXtalClus1[jxtal];
    convxtalid(iphi,ieta); /// now iphi[0,359]; ieta[-85,84];

    if( (stepc==3 && iter==11) || stepc==4 ){
      if(flag_ietaiphi[ieta+85][iphi] > 0) continue; 
    }
    
    ietaxt[ixtal] = ieta;
    iphixt[ixtal] = iphi;

    ext[ixtal] = eXtalClus1[jxtal] * CalInput[ieta+85][iphi];
    int ism = convertIetaIphiToSMNumber(ietaXtalClus1[jxtal],iphiXtalClus1[jxtal]);
    if(stepc >=2){
      ext[ixtal] *= corrfactorEta[ieta+85];
      ext[ixtal] *= corrfactorPhi[iphi];
      ext[ixtal] *= corrfactorSM[ism-1]; 
    }

    msum9[0] += ext[ixtal];
    ixtal ++; 
    
  } 
  
  if(msum9[0]<=0  ) {
    res[0] = -1; /// return now.
    return; 
  }
  
  calcClusterLocationBarrel(msum9[0],ext,ietaxt,iphixt,ixtal,posi);
  TVector3 vtmp(posi[0]- vBeamSpot[0] ,posi[1]- vBeamSpot[1],posi[2]- vBeamSpot[2]);
  meta[0] = vtmp.Eta();
  mphi[0] = vtmp.Phi();

  
  msum9[1] = 0; 
  ixtal = 0; 
  for(Int_t jxtal=0; jxtal  < nxtClus2; jxtal++){
    int ieta = ietaXtalClus2[jxtal];
    int iphi = iphiXtalClus2[jxtal];
    
    convxtalid(iphi,ieta); /// now iphi[0,359]; ieta[-85,84];
    if( (stepc==3 && iter==11) || stepc==4 ){
      if(flag_ietaiphi[ieta+85][iphi] > 0) continue;
    }
    
    ietaxt[ixtal] = ieta;
    iphixt[ixtal] = iphi;
    ext[ixtal] = eXtalClus2[jxtal] * CalInput[ieta+85][iphi];
    int ism = convertIetaIphiToSMNumber(ietaXtalClus2[jxtal],iphiXtalClus2[jxtal]);
    if(stepc >=2){
      ext[ixtal] *= corrfactorEta[ieta+85];
      ext[ixtal] *= corrfactorPhi[iphi];
      ext[ixtal] *= corrfactorSM[ism-1]; 
    }
    
    msum9[1] += ext[ixtal];
    ixtal ++; 
    
  } 
  
  if(msum9[1]<=0 ) {
    res[0] = -1; /// return now.
    return; 
  }
  
  calcClusterLocationBarrel(msum9[1],ext,ietaxt,iphixt,ixtal,posi);
  TVector3 vtmp2(posi[0]- vBeamSpot[0] ,posi[1]- vBeamSpot[1],posi[2]- vBeamSpot[2]);
  meta[1] = vtmp2.Eta();
  mphi[1] = vtmp2.Phi();
  
  
  ///Now compute invariant mass
  float cosd = getcosd(meta[0],mphi[0],meta[1],mphi[1]);
  res[0] = sqrt(2 * msum9[0] * msum9[1] * (1-cosd)); 
  
  if( stepc==4){
    ///correction due to dead crystals
    int seed1eta = ietaXtalClus1[0];
    int seed1phi = iphiXtalClus1[0];
    convxtalid(seed1phi,seed1eta);
    int seed2eta = ietaXtalClus2[0];
    int seed2phi = iphiXtalClus2[0];
    convxtalid(seed2phi,seed2eta);
    
    int ndeadflag1 = ndeadflag_ietaiphi[seed1eta+85][seed1phi];
    int ndeadflag2 = ndeadflag_ietaiphi[seed2eta+85][seed2phi];
    
    if(flag_ietaiphi[seed1eta+85][seed1phi] < 1 && flag_ietaiphi[seed2eta+85][seed2phi] <1){
      if(ndeadflag1<0 || ndeadflag2<0){
	cout<<"fatal error deadflag! "<< ndeadflag1 <<" "<<ietaXtalClus1[0] <<" "<< iphiXtalClus1[0] <<" "<< ndeadflag2 <<" "<<ietaXtalClus2[0] <<" "<< iphiXtalClus2[0] <<endl; 
	exit(1); 
      }
      res[0] *= corrfactorDead[ndeadflag1] * corrfactorDead[ndeadflag2] ; 
    }
  }
  
  res[1] = msum9[0];
  res[2] = msum9[1];
  
}


