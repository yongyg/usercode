void setTrainingCut(){
  
  
  
  ///some cuts to remove outliers in the distributions
  
  if(trainecal==1){
    mycut = TString("abs(etasc)<1.48 && sigietaieta >0.005 && sigietaieta< 0.013 && scetawidth >0.002 && scetawidth < 0.025 && scphiwidth >0.005 && scphiwidth < 0.12 && r9 >0.2 && r9 < 1.0001  && e5x5/escraw > 0.2 && e5x5/escraw <1.0001 && escraw >20 && escraw <500  && e2x5left/escraw >=0 && e2x5left/escraw <0.6  && e2x5right/escraw >=0 && e2x5right/escraw <0.6  && e2x5top/escraw >=0 && e2x5top/escraw <0.6 && e2x5bottom/escraw >=0 && e2x5bottom/escraw<0.6 && etop/bce >=0 && etop/bce<0.5 && ebottom/bce >=0 && ebottom/bce<0.5 && eleft/bce >=0 && eleft/bce<0.5 && eright/bce >=0 && eright/bce<0.5 && emax/bce >0.1 && emax/bce<0.9&& hoe>=0 && hoe<0.10 && nVertex >=1 && nVertex < 40 &&  escraw/etrue>0.6 && escraw/etrue<1.2 && evtNumber%2==0");
  }
  else{
    mycut =  TString("abs(etasc)>1.48 && eps/escraw>=0 && eps/escraw<0.30 && sigietaieta >0.015 && sigietaieta< 0.045  && scetawidth >0.005 && scetawidth < 0.06 && scphiwidth >0.005 && scphiwidth < 0.12  && r9 >0.2 && r9 < 1.0001  && e5x5/escraw > 0.2 && e5x5/escraw <1.0001 && escraw >40 && escraw <500  && e2x5left/escraw >=0 && e2x5left/escraw <0.6  && e2x5right/escraw >=0 && e2x5right/escraw <0.6  && e2x5top/escraw >=0 && e2x5top/escraw <0.6  && e2x5bottom/escraw >=0 && e2x5bottom/escraw<0.6  && etop/escraw >=0 && etop/escraw<0.5  && ebottom/escraw >=0 && ebottom/escraw<0.5  && eleft/escraw >=0 && eleft/escraw<0.5  && eright/escraw >=0 && eright/escraw<0.5  && emax/escraw >0 && emax/escraw<0.9 && hoe>=0 && hoe<0.15 && nVertex >=1 && nVertex < 40 && (escraw+eps)/etrue>0.6 && (escraw+eps)/etrue<1.2 && evtNumber%2==0");
  }
  
  
}
