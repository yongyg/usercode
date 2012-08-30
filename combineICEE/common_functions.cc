

void get_xyzEBrechits(){
  TChain *ch = new TChain("Analysis");
  ch->Add("xyzECAL.root");
  
  ch->SetBranchAddress("xEBAll",xEBAll);
  ch->SetBranchAddress("yEBAll",yEBAll);
  ch->SetBranchAddress("zEBAll",zEBAll);
  ch->SetBranchAddress("etaEBAll",etaEBAll);
  ch->SetBranchAddress("phiEBAll",phiEBAll);

  
  ch->SetBranchAddress("dxEBAll",dxEBAll);
  ch->SetBranchAddress("dyEBAll",dyEBAll);
  ch->SetBranchAddress("dzEBAll",dzEBAll);
  
  
  ch->SetBranchAddress("xEEAll",xEEAll);
  ch->SetBranchAddress("yEEAll",yEEAll);
  ch->SetBranchAddress("zEEAll",zEEAll);
  
  ch->SetBranchAddress("dxEEAll",dxEEAll);
  ch->SetBranchAddress("dyEEAll",dyEEAll);
  ch->SetBranchAddress("dzEEAll",dzEEAll);
  
  
  
  ch->SetBranchAddress("etaEEAll",etaEEAll);
  ch->SetBranchAddress("phiEEAll",phiEEAll);
  
  
  

  
  ch->GetEntry(0);

  ch->Delete();
  
  cout<<"got xyzEcal.."<<endl; 
  
  
    
}


/////////////////////////////////////////////////////////////////////////////////////////////////

void convxtalid(Int_t &nphi,Int_t &neta)
{
  // Changed to what Yong's convention; output will give just two indices
  // phi is unchanged; only eta now runs from
  //
  // 03/01/2008 changed to the new definition in CMSSW. The output is still the same...
  // Barrel only
  // Output nphi 0...359; neta 0...84; nside=+1 (for eta>0), or 0 (for eta<0).
  // neta will be [-85,-1] , or [0,84], the minus sign indicates the z<0 side.

     if(neta > 0) neta -= 1;
     if(nphi > 359) nphi=nphi-360;

     // final check
   if(nphi >359 || nphi <0 || neta< -85 || neta > 84)
     {
       cout <<" output not in range: "<<  nphi <<  " " << neta <<  " " <<endl;
    exit(1);
     }
} //end of convxtalid


// Calculate the distance in xtals taking into account possibly different sides
// change to coincide with yongs definition
Int_t diff_neta(Int_t neta1, Int_t neta2){
    Int_t mdiff;
    mdiff=abs(neta1-neta2);
    return mdiff;
}
 
// Calculate the absolute distance in xtals taking into account the periodicity of the Barrel
Int_t diff_nphi(Int_t nphi1,Int_t nphi2) {
  Int_t mdiff;
   mdiff=abs(nphi1-nphi2);
   if (mdiff > (360-abs(nphi1-nphi2))) mdiff=(360-abs(nphi1-nphi2));
   return mdiff;
}

// Calculate the distance in xtals taking into account possibly different sides
// Then the distance would be from the 1st to the 2nd argument
// _s means that it gives the sign; the only difference from the above !
// also changed to coincide with Yong's definition
Int_t diff_neta_s(Int_t neta1, Int_t neta2){
  Int_t mdiff;
  mdiff=(neta1-neta2);
  return mdiff;
}
 
// Calculate the distance in xtals taking into account the periodicity of the Barrel
Int_t diff_nphi_s(Int_t nphi1,Int_t nphi2) {
   Int_t mdiff;
   if(abs(nphi1-nphi2) < (360-abs(nphi1-nphi2))) {
     mdiff=nphi1-nphi2;
   }
   else {
   mdiff=360-abs(nphi1-nphi2);
   if(nphi1>nphi2) mdiff=-mdiff;
   }
   return mdiff;
}

