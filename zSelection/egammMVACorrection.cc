
///photon/electron TMVA correction


void egammaMVACorrection_LoadWeights(){
  
  reader_eb_unconv = new TMVA::Reader( "!Color:!Silent" );    
  reader_eb_conv = new TMVA::Reader( "!Color:!Silent" );   
  reader_ee_unconv = new TMVA::Reader( "!Color:!Silent" );    
  reader_ee_conv = new TMVA::Reader( "!Color:!Silent" );  

  methodName = TString("BDTG") + TString(" method");
  
  
  
    
    
  ///EB trained on unconv photon 
  reader_eb_unconv->AddVariable("escraw", &escraw); 
  reader_eb_unconv->AddVariable("e5x5/escraw", &e5x5_escraw); 
    
  reader_eb_unconv->AddVariable("etasc", &etasc);
  reader_eb_unconv->AddVariable("phisc",&phisc);

  reader_eb_unconv->AddVariable("scbcfmax",&scbcfmax);
    
  reader_eb_unconv->AddVariable("scnbc",&scnbcf);
  reader_eb_unconv->AddVariable("scncrystal",&scncrystalf);
    
  reader_eb_unconv->AddVariable("e2x5left/escraw",&e2x5left_escraw);
  reader_eb_unconv->AddVariable("e2x5right/escraw", &e2x5right_escraw);
  reader_eb_unconv->AddVariable("e2x5top/escraw", &e2x5top_escraw);
  reader_eb_unconv->AddVariable("e2x5bottom/escraw", &e2x5bottom_escraw);
    
    
  reader_eb_unconv->AddVariable("emax/escraw",&emax_escraw);
  reader_eb_unconv->AddVariable("eleft/escraw",&eleft_escraw);
  reader_eb_unconv->AddVariable("eright/escraw", &eright_escraw);
  reader_eb_unconv->AddVariable("etop/escraw", &etop_escraw);
  reader_eb_unconv->AddVariable("ebottom/escraw", &ebottom_escraw);
    
      
  reader_eb_unconv->AddVariable("sigietaieta",&sigietaieta);
  reader_eb_unconv->AddVariable("sigiphiiphi",&sigiphiiphi);
  reader_eb_unconv->AddVariable("sigcovietaiphi",&sigcovietaiphi);
    
  reader_eb_unconv->AddVariable("r9",&r9);
  reader_eb_unconv->AddVariable("etaSGap",&etaSGap);
  reader_eb_unconv->AddVariable("phiSGap",&phiSGap);
  reader_eb_unconv->AddVariable("etaMGap",&etaMGap);
  reader_eb_unconv->AddVariable("phiMGap",&phiMGap);
  reader_eb_unconv->AddVariable("etaCGap", &etaCGap);
  reader_eb_unconv->AddVariable("phiCGap",&phiCGap);
  reader_eb_unconv->AddVariable("nVertex",&nVertexf);
  reader_eb_unconv->AddVariable("rho",&rho);
  reader_eb_unconv->AddVariable("hoe",&hoe);
    
    
  ///EE unconv
  reader_ee_unconv->AddVariable("escraw", &escraw); 
  reader_ee_unconv->AddVariable("e5x5/escraw", &e5x5_escraw); 
  reader_ee_unconv->AddVariable("eps/escraw",&eps_escraw);
  reader_ee_unconv->AddVariable("etasc", &etasc);
  reader_ee_unconv->AddVariable("phisc",&phisc);
  reader_ee_unconv->AddVariable("scbcfmax",&scbcfmax);
  reader_ee_unconv->AddVariable("scnbc",&scnbcf);
  reader_ee_unconv->AddVariable("scncrystal",&scncrystalf);
  reader_ee_unconv->AddVariable("e2x5left/escraw",&e2x5left_escraw);
  reader_ee_unconv->AddVariable("e2x5right/escraw", &e2x5right_escraw);
  reader_ee_unconv->AddVariable("e2x5top/escraw", &e2x5top_escraw);
  reader_ee_unconv->AddVariable("e2x5bottom/escraw", &e2x5bottom_escraw);
    
  reader_ee_unconv->AddVariable("emax/escraw",&emax_escraw);
  reader_ee_unconv->AddVariable("eleft/escraw",&eleft_escraw);
  reader_ee_unconv->AddVariable("eright/escraw", &eright_escraw);
  reader_ee_unconv->AddVariable("etop/escraw", &etop_escraw);
  reader_ee_unconv->AddVariable("ebottom/escraw", &ebottom_escraw);
      
  reader_ee_unconv->AddVariable("sigietaieta",&sigietaieta);
  reader_ee_unconv->AddVariable("sigiphiiphi",&sigiphiiphi);
  reader_ee_unconv->AddVariable("sigcovietaiphi",&sigcovietaiphi);
    
  reader_ee_unconv->AddVariable("r9",&r9);
    
    
  reader_ee_unconv->AddVariable("nVertex",&nVertexf);
  reader_ee_unconv->AddVariable("rho",&rho);
      
  TString weightfile; 
  
  weightfile = TString("weights/TMVARegression_BDTG.weights.unconvphtEB.xml");
  reader_eb_unconv->BookMVA( methodName, weightfile ); 
  weightfile = TString("weights/TMVARegression_BDTG.weights.unconvphtEE.xml");
  reader_ee_unconv->BookMVA( methodName, weightfile ); 
    
    
  
}






float getElectronCorrection(int j){
  
  
  if( j <0  || j > nElectron){
    cout<<" getElectronCorrection " << j <<endl; 
    exit(1);
  }
  escraw = electronscrawEnergy[j]; 
  e5x5 = electrone5x5[j]; 
  r9 = electrone3x3[j]/ escraw;
  
  
  e2x2 = electrone2x2[j];
  e3x3 = electrone3x3[j];
  e1x3 = electrone1x3[j];
  e3x1 = electrone3x1[j];
  e4x4 = electrone4x4[j];

  e2x5right = electrone2x5Right[j];
  e2x5left = electrone2x5Left[j];
  e2x5top = electrone2x5Top[j];
  e2x5bottom = electrone2x5Bottom[j];
  e2x5max = electrone2x5Max[j];
  eleft = electroneLeft[j];
  eright = electroneRight[j];
  etop = electroneTop[j];
  ebottom = electroneBottom[j];
  
  
  eps = electronscpreshowerEnergy[j];
  esc = electronscenergy[j];
  etasc = electronsceta[j];
  phisc = electronscphi[j];
  scetawidth = electronscetaWidth[j];
  scphiwidth = electronscphiWidth[j];
  sigietaieta = electronsigmaIetaIeta[j];
  sigiphiiphi = sqrt(electronCovIPhiIPhi[j]);
  sigcovietaiphi = electronCovIEtaIPhi[j];
  emax = electroneMax[j];
  scnbc = electronscclusterSize[j];
  scncrystal = electronscnhits[j];
  scbcfmax = 0; 
  sigiphi_sigieta = sigietaieta >0 ? sigiphiiphi/sigietaieta :100;
  scsigphi_sigeta = scetawidth >0 ? scphiwidth / scetawidth: 100;
  
  e2x5left_escraw = e2x5left / escraw; 
  e2x5right_escraw = e2x5right / escraw; 
  e2x5top_escraw = e2x5top / escraw; 
  e2x5bottom_escraw = e2x5bottom / escraw; 
  
  
  vector<float> vebc = electronscbclusterenergy->at(j);
  
  
  for(int k=0; k< int( vebc.size());k++){
    if( scbcfmax < vebc[k]){
      scbcfmax = vebc[k];
    }
  }
  scbcfmax /= escraw; 
  
  
  ///float 
  scnbcf = scnbc; 
  nVertexf = nVertex;
  scncrystalf = scncrystal;
    
  TVector3 caloPosition(electroncaloPositionx[j],electroncaloPositiony[j],electroncaloPositionz[j]);
  
  getGapCoordinates(caloPosition.Eta(),caloPosition.Phi());
  etaCGap = _aC;
  phiCGap = _bC;
  etaSGap = _aS;
  phiSGap = _bS;
  etaMGap = _aM;
  phiMGap = _bM;
    
  eps_escraw = eps / escraw; 
  e5x5_escraw = e5x5 / escraw; 
  emax_escraw = emax / escraw; 
  eright_escraw = eright / escraw; 
  eleft_escraw = eleft /escraw; 
  etop_escraw = etop / escraw; 
  ebottom_escraw = ebottom / escraw; 
  hoe = electronhcalOverEcal[j];
  
  
  tmva = 0; 
  
  if( fabs(etasc) < 1.48){
    tmva = (reader_eb_unconv->EvaluateRegression(methodName))[0];
  }else{
    tmva = (reader_ee_unconv->EvaluateRegression(methodName))[0];
  }
  
  return tmva; 
  
  
  
}
