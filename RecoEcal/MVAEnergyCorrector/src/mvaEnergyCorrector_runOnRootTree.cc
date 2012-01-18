
//those functions are examples to load MVA weights and run on top of root-trees 


TMVA::Reader *reader_eb; 
TMVA::Reader *reader_ebErr; 
TMVA::Reader *reader_ee; 
TMVA::Reader *reader_eeErr; 
TString methodName;

float escraw;
float etasc;
float phisc;
float e5x5; 
float r9;
float eps; 
float e2x5left; 
float e2x5right; 
float e2x5top; 
float e2x5bottom; 
float emax; 
float eleft; 
float eright; 
float etop; 
float ebottom; 
float e5x5_escraw;
float e2x5left_escraw;
float e2x5right_escraw;
float e2x5top_escraw;
float e2x5bottom_escraw;
float emax_escraw;
float eleft_escraw;
float eright_escraw;
float etop_escraw;
float ebottom_escraw;
float eps_escraw;

float sigietaieta;
float scetawidth;
float scphiwidth;
float etaSGap;
float phiSGap;
float etaMGap;
float phiMGap;
float etaCGap;
float phiCGap;
float nVertexf;

float tmva;
float tmvaErr;

///weightsPath is a vector of string, note the order , eb ebErr, ee, eeErr
void egammaMVACorrection_LoadWeights(vector<TString> weightsFilePath){
  
  reader_eb = new TMVA::Reader( "!Color:!Silent" );    
  reader_ebErr = new TMVA::Reader( "!Color:!Silent" );    
  reader_ee = new TMVA::Reader( "!Color:!Silent" );    
  reader_eeErr = new TMVA::Reader( "!Color:!Silent" );    
  
  
  methodName = TString("BDTG") + TString(" method");
  
  
  ///EB
  reader_eb->AddVariable("escraw", &escraw); 
  reader_eb->AddVariable("etasc", &etasc);
  reader_eb->AddVariable("phisc",&phisc);
  reader_eb->AddVariable("e5x5/escraw", &e5x5_escraw); 
  reader_eb->AddVariable("r9",&r9);
  reader_eb->AddVariable("e2x5left/escraw",&e2x5left_escraw);
  reader_eb->AddVariable("e2x5right/escraw", &e2x5right_escraw);
  reader_eb->AddVariable("e2x5top/escraw", &e2x5top_escraw);
  reader_eb->AddVariable("e2x5bottom/escraw", &e2x5bottom_escraw);
  reader_eb->AddVariable("emax/escraw",&emax_escraw);
  reader_eb->AddVariable("eleft/escraw",&eleft_escraw);
  reader_eb->AddVariable("eright/escraw", &eright_escraw);
  reader_eb->AddVariable("etop/escraw", &etop_escraw);
  reader_eb->AddVariable("ebottom/escraw", &ebottom_escraw);
  reader_eb->AddVariable("sigietaieta",&sigietaieta);
  reader_eb->AddVariable("scetawidth",&scetawidth);
  reader_eb->AddVariable("scphiwidth",&scphiwidth);
  reader_eb->AddVariable("etaSGap",&etaSGap);
  reader_eb->AddVariable("phiSGap",&phiSGap);
  reader_eb->AddVariable("etaMGap",&etaMGap);
  reader_eb->AddVariable("phiMGap",&phiMGap);
  reader_eb->AddVariable("etaCGap", &etaCGap);
  reader_eb->AddVariable("phiCGap",&phiCGap);
  reader_eb->AddVariable("nVertex",&nVertexf);
  reader_eb->AddVariable("rho",&rho);
  
  //EBErr
  reader_ebErr->AddVariable("escraw", &escraw); 
  reader_ebErr->AddVariable("etasc", &etasc);
  reader_ebErr->AddVariable("phisc",&phisc);
  reader_ebErr->AddVariable("e5x5/escraw", &e5x5_escraw); 
  reader_ebErr->AddVariable("r9",&r9);
  reader_ebErr->AddVariable("e2x5left/escraw",&e2x5left_escraw);
  reader_ebErr->AddVariable("e2x5right/escraw", &e2x5right_escraw);
  reader_ebErr->AddVariable("e2x5top/escraw", &e2x5top_escraw);
  reader_ebErr->AddVariable("e2x5bottom/escraw", &e2x5bottom_escraw);
  reader_ebErr->AddVariable("emax/escraw",&emax_escraw);
  reader_ebErr->AddVariable("eleft/escraw",&eleft_escraw);
  reader_ebErr->AddVariable("eright/escraw", &eright_escraw);
  reader_ebErr->AddVariable("etop/escraw", &etop_escraw);
  reader_ebErr->AddVariable("ebottom/escraw", &ebottom_escraw);
  reader_ebErr->AddVariable("sigietaieta",&sigietaieta);
  reader_ebErr->AddVariable("scetawidth",&scetawidth);
  reader_ebErr->AddVariable("scphiwidth",&scphiwidth);
  reader_ebErr->AddVariable("etaSGap",&etaSGap);
  reader_ebErr->AddVariable("phiSGap",&phiSGap);
  reader_ebErr->AddVariable("etaMGap",&etaMGap);
  reader_ebErr->AddVariable("phiMGap",&phiMGap);
  reader_ebErr->AddVariable("etaCGap", &etaCGap);
  reader_ebErr->AddVariable("phiCGap",&phiCGap);
  reader_ebErr->AddVariable("nVertex",&nVertexf);
  reader_ebErr->AddVariable("rho",&rho);
  reader_ebErr->AddVariable("escregcorr/escraw",&tmva);
  


  //EE 
  reader_ee->AddVariable("escraw", &escraw); 
  reader_ee->AddVariable("etasc", &etasc);
  reader_ee->AddVariable("phisc",&phisc);
  reader_ee->AddVariable("eps/escraw", &eps_escraw); 
  reader_ee->AddVariable("e5x5/escraw", &e5x5_escraw); 
  reader_ee->AddVariable("r9",&r9);
  reader_ee->AddVariable("e2x5left/escraw",&e2x5left_escraw);
  reader_ee->AddVariable("e2x5right/escraw", &e2x5right_escraw);
  reader_ee->AddVariable("e2x5top/escraw", &e2x5top_escraw);
  reader_ee->AddVariable("e2x5bottom/escraw", &e2x5bottom_escraw);
  reader_ee->AddVariable("emax/escraw",&emax_escraw);
  reader_ee->AddVariable("eleft/escraw",&eleft_escraw);
  reader_ee->AddVariable("eright/escraw", &eright_escraw);
  reader_ee->AddVariable("etop/escraw", &etop_escraw);
  reader_ee->AddVariable("ebottom/escraw", &ebottom_escraw);
  reader_ee->AddVariable("sigietaieta",&sigietaieta);
  reader_ee->AddVariable("scetawidth",&scetawidth);
  reader_ee->AddVariable("scphiwidth",&scphiwidth);
  reader_ee->AddVariable("nVertex",&nVertexf);
  reader_ee->AddVariable("rho",&rho);

  //EEerr
  reader_eeErr->AddVariable("escraw", &escraw); 
  reader_eeErr->AddVariable("etasc", &etasc);
  reader_eeErr->AddVariable("phisc",&phisc);
  reader_eeErr->AddVariable("eps/escraw", &eps_escraw); 
  reader_eeErr->AddVariable("e5x5/escraw", &e5x5_escraw); 
  reader_eeErr->AddVariable("r9",&r9);
  reader_eeErr->AddVariable("e2x5left/escraw",&e2x5left_escraw);
  reader_eeErr->AddVariable("e2x5right/escraw", &e2x5right_escraw);
  reader_eeErr->AddVariable("e2x5top/escraw", &e2x5top_escraw);
  reader_eeErr->AddVariable("e2x5bottom/escraw", &e2x5bottom_escraw);
  reader_eeErr->AddVariable("emax/escraw",&emax_escraw);
  reader_eeErr->AddVariable("eleft/escraw",&eleft_escraw);
  reader_eeErr->AddVariable("eright/escraw", &eright_escraw);
  reader_eeErr->AddVariable("etop/escraw", &etop_escraw);
  reader_eeErr->AddVariable("ebottom/escraw", &ebottom_escraw);
  reader_eeErr->AddVariable("sigietaieta",&sigietaieta);
  reader_eeErr->AddVariable("scetawidth",&scetawidth);
  reader_eeErr->AddVariable("scphiwidth",&scphiwidth);
  reader_eeErr->AddVariable("nVertex",&nVertexf);
  reader_eeErr->AddVariable("rho",&rho);
  reader_eeErr->AddVariable("escregcorr/(escraw+eps)",&tmva);
  
  
  //Book MVA 
  reader_eb->BookMVA( methodName, weightsFilePath[0]);
  reader_ebErr->BookMVA( methodName, weightsFilePath[1]);
  reader_ee->BookMVA( methodName, weightsFilePath[2]);
  reader_eeErr->BookMVA( methodName, weightsFilePath[3]);
    
}


std::pair<float,float> getPhotonCorrection(int j){
  
  if( j <0  || j > nPhoton){
    cout<<" getPhotonCorrection " << j <<endl; 
    exit(1);
  }
  
  ////=== those variables on the right of the equations are saved in the root-trees
  escraw = photonscrawEnergy[j]; 
  etasc = photonsceta[j];
  phisc = photonscphi[j];
  e5x5 = photone5x5[j]; 
  r9 = photonr9[j];
  e2x5right = photone2x5Right[j];
  e2x5left = photone2x5Left[j];
  e2x5top = photone2x5Top[j];
  e2x5bottom = photone2x5Bottom[j];
  e2x5max = photone2x5Max[j];
  emax = photonmaxEnergyXtal[j];
  eleft = photoneLeft[j];
  eright = photoneRight[j];
  etop = photoneTop[j];
  ebottom = photoneBottom[j];
  eps = photonscpreshowerEnergy[j];
  sigietaieta = photonsigmaIetaIeta[j];
  scetawidth = photonscetaWidth[j];
  scphiwidth = photonscphiWidth[j];

  nVertexf = nVertex;
  ////===
  
  ////rho should be already there ( fChain->GetEntry(entry) ) 
  

  e5x5_escraw = e5x5 / escraw; 
  e2x5left_escraw = e2x5left / escraw; 
  e2x5right_escraw = e2x5right / escraw; 
  e2x5top_escraw = e2x5top / escraw; 
  e2x5bottom_escraw = e2x5bottom / escraw; 
  emax_escraw = emax / escraw; 
  eright_escraw = eright / escraw; 
  eleft_escraw = eleft /escraw; 
  etop_escraw = etop / escraw; 
  ebottom_escraw = ebottom / escraw; 
  eps_escraw = eps / escraw; 
  
  getGapCoordinates(etasc,phisc);
  
  etaCGap = _aC;
  phiCGap = _bC;
  etaSGap = _aS;
  phiSGap = _bS;
  etaMGap = _aM;
  phiMGap = _bM;
  
  
  
  tmva = 0; 
  tmvaErr = 0;
    
  if( fabs(etasc) < 1.48){
    tmva = (reader_eb->EvaluateRegression(methodName))[0];
    tmvaErr = (reader_ebErr->EvaluateRegression(methodName))[0];
  }else{
    tmva = (reader_ee->EvaluateRegression(methodName))[0];
    tmvaErr = (reader_eeErr->EvaluateRegression(methodName))[0];
  }
    
  float varscale = 1; 
  float ecorr = (escraw+eps)*tmva; 
  float ecorrErr = corr*tmvaErr * varscale; 
  return std::pair<float,float>(ecorr,ecorrErr);
  
}



std::pair<float,float> getElectronCorrection(int j){
  
  if( j <0  || j > nElectron){
    cout<<" getElectronCorrection " << j <<endl; 
    exit(1);
  }
  
  ////=== those variables on the right of the equations are saved in the root-trees
  escraw = electronscrawEnergy[j]; 
  etasc = electronsceta[j];
  phisc = electronscphi[j];
  e5x5 = electrone5x5[j]; 
  r9 = electrone3x3[j]/ escraw;
  e2x5right = electrone2x5Right[j];
  e2x5left = electrone2x5Left[j];
  e2x5top = electrone2x5Top[j];
  e2x5bottom = electrone2x5Bottom[j];
  e2x5max = electrone2x5Max[j];
  emax = electroneMax[j];
  eleft = electroneLeft[j];
  eright = electroneRight[j];
  etop = electroneTop[j];
  ebottom = electroneBottom[j];
  eps = electronscpreshowerEnergy[j];
  scetawidth = electronscetaWidth[j];
  scphiwidth = electronscphiWidth[j];
  sigietaieta = electronsigmaIetaIeta[j];
  ////===
  
  ////rho should be already there ( fChain->GetEntry(entry) ) 
  

  getGapCoordinates(etasc,scphi);
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
  e2x5left_escraw = e2x5left / escraw; 
  e2x5right_escraw = e2x5right / escraw; 
  e2x5top_escraw = e2x5top / escraw; 
  e2x5bottom_escraw = e2x5bottom / escraw; 
  nVertexf = nVertex;
  
  tmva = 0; 
  tmvaErr = 0; 
  
  if( fabs(etasc) < 1.48){
    tmva = (reader_eb->EvaluateRegression(methodName))[0];
    tmvaErr = (reader_ebErr->EvaluateRegression(methodName))[0];
  }else{
    tmva = (reader_ee->EvaluateRegression(methodName))[0];
    tmvaErr = (reader_eeErr->EvaluateRegression(methodName))[0];
  }
  
  
  float varscale = 1; 
  float ecorr = (escraw+eps)*tmva; 
  float ecorrErr = corr*tmvaErr * varscale; 
  return std::pair<float,float>(ecorr,ecorrErr);
  
}
