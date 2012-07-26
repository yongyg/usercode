

////to make the cmssw official-like cluster position


void calcClusterLocationBarrel(float eTot,float energy[],int ietaHit[],int iphiHit[],int nhit,float posi[]){
  
  float param_T0_barl_ = 7.4; 
  // float param_T0_endc_ = 3.1; 
  //  float param_t0_endcPresh = 1.2; 
  float param_W0_ = 4.2; 
  float param_X0_ = 0.89; 
  
  ///no preshower..
  ///bool preshower = false;
  
  if(eTot<=0 || nhit <1){
    cout<<"fatal error calcClusterLocation ..negative energy "<<eTot<<" "<<nhit<<endl;
    exit(1);
  }
  
  
  int nmaxhit = 0; 
  float eMax = 0; 

  for(int j=0; j<nhit; j++){
    if(eMax < energy[j]){
      eMax = energy[j];
      nmaxhit = j; 
    }
  }
  

  int ieta = ietaHit[nmaxhit] + 85; 
  int iphi = iphiHit[nmaxhit]; 
  iphi = getIndphixyzEBAll(iphi);
  
  TVector3 vmaxHit(xEBAll[ieta][iphi],yEBAll[ieta][iphi],zEBAll[ieta][iphi]);
    
  float T0 = param_T0_barl_;
  
  if(fabs(vmaxHit.Eta())>1.482){
    cout<<"error..barrel..?"<<" "<< vmaxHit.Eta()<<" "<<ieta <<" "<<iphi <<" "<<etaEBAll[ieta][iphi]<<endl; 
    exit(1);
  }
  
  float depth =  param_X0_ * (T0 + log(eTot));
  
  TVector3 dvmaxHit(depth*dxEBAll[ieta][iphi],depth*dyEBAll[ieta][iphi],depth*dzEBAll[ieta][iphi]);
  
  TVector3 center_pos = vmaxHit + dvmaxHit;
  
  
  
  // Loop over hits and get weights
  double weight = 0;
  double total_weight = 0;
  
  double center_phi = center_pos.Phi();
  double center_theta = center_pos.Theta();
  
  double delta_theta = 0;
  double delta_phi = 0;

  float w_min = 0.;
  
  double dphi = 0;
  
  const double k_PI = acos(-1);
  
  
  for(int j=0; j<nhit; j++){
    
    double e_j = energy[j];
    
    if(e_j<0) continue; 
    
    
    weight = param_W0_ +log ( e_j/eTot); 
    if(weight<w_min) weight = w_min;  
    
    total_weight += weight;

    ieta = ietaHit[j] + 85; 
    iphi = iphiHit[j]; 
    iphi = getIndphixyzEBAll(iphi);
    TVector3 vHit(xEBAll[ieta][iphi],yEBAll[ieta][iphi],zEBAll[ieta][iphi]);
    TVector3 dvHit(depth*dxEBAll[ieta][iphi],depth*dyEBAll[ieta][iphi],depth*dzEBAll[ieta][iphi]);
    TVector3 jth_pos = vHit + dvHit;
    
    
    delta_theta += weight * (jth_pos.Theta() - center_theta);
    dphi = (jth_pos.Phi() - center_phi);
    
    // Check the 2*pi problem for delta_phi
    if (dphi > k_PI)
      dphi -= 2.*k_PI;
    if (dphi < -k_PI)
      dphi += 2.*k_PI;
    
    delta_phi += dphi*weight;    
  }
  
  delta_theta /= total_weight;
  delta_phi /= total_weight;
  
  double cluster_theta = center_theta + delta_theta;
  double cluster_phi = center_phi + delta_phi;

  // Check the 2*pi problem for cluster_phi
  if (cluster_phi > k_PI)
    cluster_phi -= 2.*k_PI;
  if (cluster_phi < -k_PI)
    cluster_phi += 2.*k_PI;
  
  double radius = sqrt(center_pos.x()*center_pos.x()
		       + center_pos.y()*center_pos.y()
		       + center_pos.z()*center_pos.z());
  
  float xpos = radius*cos(cluster_phi)*sin(cluster_theta);
  float ypos = radius*sin(cluster_phi)*sin(cluster_theta);
  float zpos = radius*cos(cluster_theta);
  
  
  posi[0] = xpos; 
  posi[1] = ypos; 
  posi[2] = zpos; 
  
  


}




void calcClusterLocationEndcap(float eTot,float energy[],int ixHit[],int iyHit[],int izHit[],int nhit,float posi[]){
  
  //float param_T0_barl_ = 7.4; 
  float param_T0_endc_ = 3.1; 
  float param_t0_endcPresh = 1.2; 
  float param_W0_ = 4.2; 
  float param_X0_ = 0.89; 
  
  ///no preshower..
  bool preshower = false;
  
  if(eTot<=0 || nhit <1){
    cout<<"fatal error calcClusterLocationEndcap ..negative energy "<<eTot<<" "<<nhit<<endl;
    exit(1);
  }
  
  
  int nmaxhit = 0; 
  float eMax = 0; 

  for(int j=0; j<nhit; j++){
    if(eMax < energy[j]){
      eMax = energy[j];
      nmaxhit = j; 
    }
  }
  

  int ix = ixHit[nmaxhit];
  int iy = iyHit[nmaxhit]; 
  int iz = izHit[nmaxhit]; 
  int indz = iz > 0; 
  
  TVector3 vmaxHit(xEEAll[indz][ix][iy],yEEAll[indz][ix][iy],zEEAll[indz][ix][iy]);
  
  
  float T0 = param_T0_endc_;
  if(preshower){
    T0 = param_t0_endcPresh;
  }  
    
  if(fabs(vmaxHit.Eta())<1.482){
    cout<<"error..encap..?"<<vmaxHit.Eta()<<endl;
    exit(1);
  }
  
  float depth =  param_X0_ * (T0 + log(eTot));
  
  TVector3 dvmaxHit(depth*dxEEAll[indz][ix][iy],depth*dyEEAll[indz][ix][iy],depth*dzEEAll[indz][ix][iy]);
  
  TVector3 center_pos = vmaxHit + dvmaxHit;
  
  
  
  
  // Loop over hits and get weights
  double weight = 0;
  double total_weight = 0;

  double center_phi = center_pos.Phi();
  double center_theta = center_pos.Theta();
  
  double delta_theta = 0;
  double delta_phi = 0;

  float w_min = 0.;
  
  double dphi = 0;
  
  const double k_PI = acos(-1);
    

  for(int j=0; j<nhit; j++){
    
    double e_j = energy[j];
    
    if(e_j<0) continue; 
    
    
    weight = param_W0_ +log ( e_j/eTot); 
    if(weight<w_min) weight = w_min;  
    
    total_weight += weight;

    //    ieta = ietaHit[j] + 85; 
    /// iphi = iphiHit[j]; 

    ix = ixHit[j];
    iy = iyHit[j];
    iz = izHit[nmaxhit]; 
    indz = iz > 0;    
    
    TVector3 vHit(xEEAll[indz][ix][iy],yEEAll[indz][ix][iy],zEEAll[indz][ix][iy]);
    TVector3 dvHit(depth*dxEEAll[indz][ix][iy],depth*dyEEAll[indz][ix][iy],depth*dzEEAll[indz][ix][iy]);
    

    TVector3 jth_pos = vHit + dvHit;
    
    
    delta_theta += weight * (jth_pos.Theta() - center_theta);
    dphi = (jth_pos.Phi() - center_phi);
    
    // Check the 2*pi problem for delta_phi
    if (dphi > k_PI)
      dphi -= 2.*k_PI;
    if (dphi < -k_PI)
      dphi += 2.*k_PI;
    
    delta_phi += dphi*weight;    
  }
  
  delta_theta /= total_weight;
  delta_phi /= total_weight;
  
  double cluster_theta = center_theta + delta_theta;
  double cluster_phi = center_phi + delta_phi;

  // Check the 2*pi problem for cluster_phi
  if (cluster_phi > k_PI)
    cluster_phi -= 2.*k_PI;
  if (cluster_phi < -k_PI)
    cluster_phi += 2.*k_PI;
  
  double radius = sqrt(center_pos.x()*center_pos.x()
		       + center_pos.y()*center_pos.y()
		       + center_pos.z()*center_pos.z());
  
  float xpos = radius*cos(cluster_phi)*sin(cluster_theta);
  float ypos = radius*sin(cluster_phi)*sin(cluster_theta);
  float zpos = radius*cos(cluster_theta);
  
  
  posi[0] = xpos; 
  posi[1] = ypos; 
  posi[2] = zpos; 
  
    
}
