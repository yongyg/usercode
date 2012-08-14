
///input iz, ix, iy, 
///iz = -1 /1; ix, iy >=1 && <=100; 
int iRingEndCap(int iz, int ix,int iy){

  if(abs(iz) !=1 || !(ix>=1&&ix<=100) || !(iy>=1&&iy<=100)){
    cout<<"iRingEndcap.."<<iz<<" "<<ix <<" "<<iy<<endl; 
    exit(1);
  }
  
  
  int izz = iz <0? 0: 1; 

  int k = -1; 
  for (int ring=0; ring<kEndcEtaRings; ring++) {
    
    if( fabs(etaEEAll[izz][ix][iy])>etaBoundary_eezside[izz][ring] &&  fabs(etaEEAll[izz][ix][iy]) < etaBoundary_eezside[izz][ring+1]){
      k = ring; 
      break; 
    }
  }
  
  if(k<0){
    cout<<"error.iringendcap;"<<iz<<" "<<ix<<" "<<iy<<" "<<etaEEAll[izz][ix][iy]<<endl;
    exit(1);
  }
  
  return k; 
  
  
}



///eta lower and upper value of each eta ring in endcap
void setEtaRingBoundaryEndcap(){


  
  ///
  for(int k=0; k<kEndcEtaRings; k++){

    ///eta_ring_ee[k] = etaEEAll[1][k+1][51];
    //updated to take average of both sized
    eta_ring_ee[k] = 0.5 * ( fabs(etaEEAll[0][k+1][51]) + fabs(etaEEAll[1][k+1][51]));
    
    ///in data better use both size
    eta_ring_eezside[0][k] = fabs(etaEEAll[0][k+1][51]);
    eta_ring_eezside[1][k] = etaEEAll[1][k+1][51];
    
    if(fabs(eta_ring_ee[k])<1.482){
      cout<<"eroor..etaring..endcap? "<<eta_ring_ee[k] <<endl;
      exit(1);
    }
    //  cout<<"k: "<<eta_ring_ee[k]<<endl;
  }
  
  etaBoundary_ee[0] = 1.482;
  etaBoundary_ee[39]=4.;
  
  etaBoundary_eezside[0][0] = 1.482;
  etaBoundary_eezside[0][39]=4.;
  etaBoundary_eezside[1][0] = 1.482;
  etaBoundary_eezside[1][39]=4.;
  

  for (int ring=1; ring<kEndcEtaRings; ring++) {
    etaBoundary_ee[ring] = (eta_ring_ee[ring]+eta_ring_ee[ring-1])/2.;
    etaBoundary_eezside[0][ring] = (eta_ring_eezside[0][ring]+eta_ring_eezside[0][ring-1])/2.;
    etaBoundary_eezside[1][ring] = (eta_ring_eezside[1][ring]+eta_ring_eezside[1][ring-1])/2.;
  }
  
  ///# of xtal each ring
  for(int iz =0; iz<2 ;iz++){
    
    for(int ix=1; ix<=100; ix++){
      
      for(int iy=1; iy<=100; iy++){
	
	if(fabs(etaEEAll[iz][ix][iy])>1.482){
	  
	  int iring = iRingEndCap(2*iz-1,ix,iy);
	  nxtal_ring_ee[iring] ++; 
	}
	
      }
    }
    
  }
  
  ///should be 14648; 
  int ntotal = 0; 
  
  for(int k=0; k<kEndcEtaRings; k++){
    ntotal += nxtal_ring_ee[k]; 
  }
  if( ntotal != 14648 ){
    cout<<"checkntotalendcap " << ntotal <<endl; 
    exit(1);
  }
  
}
