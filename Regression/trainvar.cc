void setTrainingVar(){
      
  vec_trainVar.clear();
  
  if(trainecal ==1){
    
    vec_trainVar.push_back(make_pair("escraw","F"));    
    vec_trainVar.push_back(make_pair("etasc","F"));    
    vec_trainVar.push_back(make_pair("phisc","F"));    
    vec_trainVar.push_back(make_pair("e5x5/escraw","F"));    
    vec_trainVar.push_back(make_pair("r9","F"));    
    vec_trainVar.push_back(make_pair("sigietaieta","F"));
    vec_trainVar.push_back(make_pair("scetawidth","F"));
    vec_trainVar.push_back(make_pair("scphiwidth","F"));
    vec_trainVar.push_back(make_pair("e2x5left/escraw","F"));
    vec_trainVar.push_back(make_pair("e2x5right/escraw","F"));
    vec_trainVar.push_back(make_pair("e2x5top/escraw","F"));
    vec_trainVar.push_back(make_pair("e2x5bottom/escraw","F"));
    vec_trainVar.push_back(make_pair("emax/escraw","F"));
    vec_trainVar.push_back(make_pair("eleft/escraw","F"));
    vec_trainVar.push_back(make_pair("eright/escraw","F"));
    vec_trainVar.push_back(make_pair("etop/escraw","F"));
    vec_trainVar.push_back(make_pair("ebottom/escraw","F"));
    vec_trainVar.push_back(make_pair("etaSGap","F"));
    vec_trainVar.push_back(make_pair("phiSGap","F"));
    vec_trainVar.push_back(make_pair("etaMGap","F"));
    vec_trainVar.push_back(make_pair("phiMGap","F"));
    vec_trainVar.push_back(make_pair("etaCGap","F"));
    vec_trainVar.push_back(make_pair("phiCGap","F"));
    vec_trainVar.push_back(make_pair("(etop-ebottom)/(etop+ebottom)","F"));
    vec_trainVar.push_back(make_pair("(eleft-eright)/(eleft+eright)","F"));
    vec_trainVar.push_back(make_pair("seedieta","I"));
    vec_trainVar.push_back(make_pair("seediphi","I"));
    vec_trainVar.push_back(make_pair("(abs(seedieta)<=25)*(seedieta%25) + (abs(seedieta)>25)*((seedieta-25*abs(seedieta)/seedieta)%20)","I"));
    vec_trainVar.push_back(make_pair("seediphi%20","I"));
    vec_trainVar.push_back(make_pair("dbceta","F"));
    vec_trainVar.push_back(make_pair("dbcphi","F"));
    vec_trainVar.push_back(make_pair("bcseedetacry","F"));
    vec_trainVar.push_back(make_pair("bcseedphicry","F"));
    vec_trainVar.push_back(make_pair("hoe","F"));
    vec_trainVar.push_back(make_pair("nVertex","I"));
    
    
  }
  else {
    
    vec_trainVar.push_back(make_pair("escraw","F"));    
    vec_trainVar.push_back(make_pair("etasc","F"));    
    vec_trainVar.push_back(make_pair("phisc","F"));    
    vec_trainVar.push_back(make_pair("eps/escraw","F"));    
    vec_trainVar.push_back(make_pair("e5x5/escraw","F"));    
    vec_trainVar.push_back(make_pair("r9","F"));    
    vec_trainVar.push_back(make_pair("sigietaieta","F"));
    vec_trainVar.push_back(make_pair("scetawidth","F"));
    vec_trainVar.push_back(make_pair("scphiwidth","F"));
    vec_trainVar.push_back(make_pair("e2x5left/escraw","F"));
    vec_trainVar.push_back(make_pair("e2x5right/escraw","F"));
    vec_trainVar.push_back(make_pair("e2x5top/escraw","F"));
    vec_trainVar.push_back(make_pair("e2x5bottom/escraw","F"));
    vec_trainVar.push_back(make_pair("emax/escraw","F"));
    vec_trainVar.push_back(make_pair("eleft/escraw","F"));
    vec_trainVar.push_back(make_pair("eright/escraw","F"));
    vec_trainVar.push_back(make_pair("etop/escraw","F"));
    vec_trainVar.push_back(make_pair("ebottom/escraw","F"));
    vec_trainVar.push_back(make_pair("nVertex","I"));
    
    
  }
  
 
  
}
