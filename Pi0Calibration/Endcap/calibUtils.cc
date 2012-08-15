

//Input two clusters, compute all variables from that
///with new calibration calconst after each step
void calcPairClusterUpdated(float CalInput[2][101][101],float res[]){
  
  Float_t msum9[2], meta[2], mphi[2];
  
  
  int ietaxt[9]; //after convxtlid
  int iphixt[9];
  int izxt[9];
  
  float ext[9];
  
  float posi[3];
  
  
  msum9[0] = 0; 
  
  for(Int_t ixtal=0; ixtal  < nxtClus1; ixtal++){
    int ieta = ietaXtalClus1[ixtal];
    int iphi = iphiXtalClus1[ixtal];
    ietaxt[ixtal] = ieta;
    iphixt[ixtal] = iphi;
    izxt[ixtal] = izXtalClus1;
    int izz = izXtalClus1 <0 ? 0: 1; 
    ext[ixtal] = eXtalClus1[ixtal]*CalInput[izz][ieta][iphi];
    int ietaRing = iRingEndCap(izXtalClus1,ietaXtalClus1[ixtal],iphiXtalClus1[ixtal]); 
    ext[ixtal] *= corrfactorEtaRings[izz][ietaRing];
    msum9[0] += ext[ixtal];
    
  } 
  
  if(msum9[0]<=0  ) {
    res[0] = -1; /// return now.
    return; 
  }
  
  calcClusterLocationEndcap(msum9[0],ext,ietaxt,iphixt,izxt,nxtClus1,posi);
  
  msum9[0] += (infoESX[0][0] + infoESY[0][0]); ///add ES energy
  
  TVector3 vtmp(posi[0]- vBeamSpot[0] ,posi[1]- vBeamSpot[1],posi[2]- vBeamSpot[2]);
  meta[0] = vtmp.Eta();
  mphi[0] = vtmp.Phi();
  
  msum9[1] = 0; 
  
  for(Int_t ixtal=0; ixtal  < nxtClus2; ixtal++){
    int ieta = ietaXtalClus2[ixtal];
    int iphi = iphiXtalClus2[ixtal];
    ietaxt[ixtal] = ieta;
    iphixt[ixtal] = iphi;
    izxt[ixtal] = izXtalClus2;
    int izz = izXtalClus2 <0 ? 0: 1; 
    ext[ixtal] = eXtalClus2[ixtal]*CalInput[izz][ieta][iphi];
    int ietaRing = iRingEndCap(izXtalClus2,ietaXtalClus2[ixtal],iphiXtalClus2[ixtal]); 
    ext[ixtal] *= corrfactorEtaRings[izz][ietaRing];
    msum9[1] += ext[ixtal];
  } 
  if(msum9[1]<=0 ) {
    res[0] = -1; /// return now.
    return; 
  }

  calcClusterLocationEndcap(msum9[1],ext,ietaxt,iphixt,izxt,nxtClus2,posi);
  msum9[1] += (infoESX[1][0] + infoESY[1][0]);
  
  TVector3 vtmp2(posi[0]- vBeamSpot[0] ,posi[1]- vBeamSpot[1],posi[2]- vBeamSpot[2]);
  meta[1] = vtmp2.Eta();
  mphi[1] = vtmp2.Phi();

  ///Now compute invariant mass
  float cosd = getcosd(meta[0],mphi[0],meta[1],mphi[1]);
  res[0] = sqrt(2 * msum9[0] * msum9[1] * (1-cosd)); 
  res[1] = msum9[0];
  res[2] = msum9[1];
  
}

