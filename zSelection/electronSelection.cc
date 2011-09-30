
bool passElectronID(int j, int cut){
      
  float sigiee = electronsigmaIetaIeta[j];
  float deta = electrondeltaEtaSuperClusterTrackAtVtx[j];
  float dphi = electrondeltaPhiSuperClusterTrackAtVtx[j];
  hoe = electronhcalOverEcal[j];
    
  float trkiso = electrondr03TkSumPt[j];
  float ecaliso = electrondr03EcalRecHitSumEt[j];
  float hcaliso = electrondr03HcalTowerSumEt[j];
  float hcalisodep1 = electrondr03HcalDepth1TowerSumEt[j];
  float hcalisodep2 = electrondr03HcalDepth2TowerSumEt[j];
  
  float pt = electronpt[j];
  
  int nmissinghits = electronExpectedHitsInnernumberOfHits[j];
  
  float dist = electronconvDist[j];
  float dcottheta = electronconvDcot[j];
  
  
  float e2x5oe5x5 = electrone2x5Max[j]/ electrone5x5[j];
  float e1x5oe5x5 = electrone1x5[j]/electrone5x5[j];
  
  float et = electronscenergy[j] * sin(2*atan(exp(-electroneta[j])));
  etasc = electronsceta[j];
  
  
  if(cut==0){
    return true; 
  }
  else if(cut ==1){
    
    //WP95
    if( fabs(etasc) < 1.48){
      
      if( sigiee < 0.01 && fabs(deta) < 0.007 && fabs(dphi) < 0.8  && hoe < 0.15 && trkiso/pt < 0.15 
	  && ecaliso/pt < 2 && hcaliso/pt < 0.12 && nmissinghits <=1 ) return true; 
      else return false; 
      
    }else{
      if( sigiee < 0.03   && fabs(deta) < 0.01 && fabs(dphi) < 0.7  && hoe < 0.07 && trkiso/pt < 0.08
	  && ecaliso/pt < 0.06 && hcaliso/pt < 0.05 && nmissinghits <=1 ) return true; 
      else return false; 
      
    }
    
  }
  else if(cut ==2){
    
    if( fabs(etasc) < 1.48){
      
      if( fabs(deta) < 0.007 && fabs(dphi) < 0.8  &&  nmissinghits <=1 ) return true; 
      else return false; 
      
    }else{
      if( fabs(deta) < 0.01 && fabs(dphi) < 0.7 && nmissinghits <=1 ) return true; 
      else return false; 
    }
    
  }

  else if(cut ==-3){ ///corrected for PU isolation, set rho = 0 

    rho = 0; 

    //WP80
    if( fabs(etasc) < 1.48){
      
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02)  &&   sigiee < 0.01 && fabs(deta) < 0.004 && fabs(dphi) < 0.06  && hoe < 0.04 && (trkiso- rho * 0 )/pt < 0.09 
	  && (ecaliso- rho * 0.096) /pt < 0.07 && (hcaliso - rho * 0.020) /pt < 0.1 && nmissinghits <=0 ) return true; 
      else return false; 
      
    }else{
      if(  !(fabs(dist)<0.02 && fabs(dcottheta)<0.02) && sigiee < 0.03 && fabs(deta) < 0.007 && fabs(dphi) < 0.03   && hoe < 0.15 && (trkiso - rho * 0) /pt < 0.04
	   && (ecaliso - rho * 0.044) /pt < 0.05 && (hcaliso - rho * 0.041) /pt < 0.025 && nmissinghits <=0 ) return true; 
      else return false; 
      
    }
    
  }


  else if(cut ==3){ ///corrected for PU isolation
    //WP80
    if( fabs(etasc) < 1.48){
      
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02)  &&   sigiee < 0.01 && fabs(deta) < 0.004 && fabs(dphi) < 0.06  && hoe < 0.04 && (trkiso- rho * 0 )/pt < 0.09 
	  && (ecaliso- rho * 0.096) /pt < 0.07 && (hcaliso - rho * 0.020) /pt < 0.1 && nmissinghits <=0 ) return true; 
      else return false; 
      
    }else{
      if(  !(fabs(dist)<0.02 && fabs(dcottheta)<0.02) && sigiee < 0.03 && fabs(deta) < 0.007 && fabs(dphi) < 0.03   && hoe < 0.15 && (trkiso - rho * 0) /pt < 0.04
	   && (ecaliso - rho * 0.044) /pt < 0.05 && (hcaliso - rho * 0.041) /pt < 0.025 && nmissinghits <=0 ) return true; 
      else return false; 
      
    }
    
  }
  
  
  


  










  else if(cut ==4){ //WP80 no isolation,no sigee
        
    //WP80
    if( fabs(etasc) < 1.48){
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02) && fabs(deta) < 0.004 && fabs(dphi) < 0.06  &&  nmissinghits <=1 ) return true; 
      else return false; 
      
    }else{
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02) && fabs(deta) < 0.007 && fabs(dphi) < 0.03 && nmissinghits <=1 ) return true; 
      else return false; 
    }
  }
  
  else if(cut ==5){
    //WP60
    if( fabs(etasc) < 1.48){
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02)  &&   sigiee < 0.01 && fabs(deta) < 0.025 && fabs(dphi) < 0.04  && hoe < 0.025 && trkiso/pt < 0.04 
	  && ecaliso/pt < 0.04 && hcaliso/pt < 0.03 && nmissinghits <=0 ) return true; 
      else return false; 
      
    }else{
      if(  !(fabs(dist)<0.02 && fabs(dcottheta)<0.02) && sigiee < 0.03 && fabs(deta) < 0.005 && fabs(dphi) < 0.02   && hoe < 0.025 && trkiso/pt < 0.025
	   && ecaliso/pt < 0.02 && hcaliso/pt < 0.02 && nmissinghits <=0 ) return true; 
      else return false; 
      
    }
    
  }
  else if(cut ==6){
    //WP60 no sigiee no isolation 
    if( fabs(etasc) < 1.48){
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02) && fabs(deta) < 0.004 && fabs(dphi) < 0.025  &&  nmissinghits <=1 ) return true; 
      else return false; 
      
    }else{
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02) && fabs(deta) < 0.005 && fabs(dphi) < 0.02 && nmissinghits <=1 ) return true; 
      else return false; 
    }
    
  }
  
  else if(cut ==7){ //WW from xiesi
    if( fabs(etasc) < 1.48){
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02 ) && sigiee < 0.01 && fabs(deta) < 0.004 && fabs(dphi) < 0.06  && hoe < 0.04 && ( trkiso + max(ecaliso-1,float(0)) + hcaliso)/pt < 0.1
	  ) return true; 
      else return false; 
    }else{
      if( !(fabs(dist)<0.02 && fabs(dcottheta)<0.02 ) && sigiee < 0.03   && fabs(deta) < 0.007 && fabs(dphi) < 0.3  && hoe < 0.025 && ( trkiso + ecaliso + hcaliso)/pt < 0.1 ) return true; 
      else return false; 
    }
  }
  
  
  ///HEEP ID
  
  else if(cut ==10){
    

    
    if( fabs(etasc) < 1.48){
      
      if(fabs(etasc)<1.442 && fabs(deta) < 0.005 && fabs(dphi)<0.09 && hoe < 0.05 && (e2x5oe5x5>0.94|| e1x5oe5x5>0.83) && (ecaliso + hcalisodep1) < 2 + 0.03 * et && trkiso < 7.5) return true; 
      else {

	//if( fabs(etasc)<1.442) {cout<<" check "<< etasc<<" "<< deta <<" "<< dphi <<" "<< hoe <<" "<< e2x5oe5x5ele <<" "<< ecaliso+hcalisodep1<<" "<< 2 + 0.03 * et<<" "<<trkiso <<" "<< endl; 
	//	}
	
	return false; 
      }
    }else{
      if(fabs(etasc)<2.5 && fabs(etasc)>1.560 && fabs(deta) < 0.007 && fabs(dphi)<0.09 && hoe < 0.05  && trkiso < 15 
	 && hcalisodep2 < 0.5  && sigiee < 0.03 
	 ){
	if( et <50){
	  if( ecaliso + hcalisodep1  < 2.5) return true; 
	  else return false; 
	}else{
	  if( ecaliso + hcalisodep1  < 2.5 + 0.03*(et-50) ) return true; 
	  else return false; 
	}
      }else return false; 
    }
  }
  
  
  
  else{
    cout<<"passElectronID NA.."  <<endl;
    exit(1);
  }
  
  
  
}
