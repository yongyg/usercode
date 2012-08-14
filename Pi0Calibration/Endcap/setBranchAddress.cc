
void setBranchAddress(){


  cout<<" setBranchAddress for endcap " <<endl; 
fChain->SetBranchAddress("mpair", &mpair);

fChain->SetBranchAddress("lumiBlock", &lumiBlock);
fChain->SetBranchAddress("runNumber", &runNumber);
fChain->SetBranchAddress("evtNumber", &evtNumber);
fChain->SetBranchAddress("nxtClus1", &nxtClus1);
fChain->SetBranchAddress("eXtalClus1", eXtalClus1);
fChain->SetBranchAddress("ietaXtalClus1", ietaXtalClus1);
fChain->SetBranchAddress("iphiXtalClus1", iphiXtalClus1);
fChain->SetBranchAddress("izXtalClus1", &izXtalClus1);
 
fChain->SetBranchAddress("infoESX", infoESX);
fChain->SetBranchAddress("infoESY", infoESY);
 
//fChain->SetBranchAddress("tXtalClus1", tXtalClus1);

fChain->SetBranchAddress("nxtClus2", &nxtClus2);
fChain->SetBranchAddress("eXtalClus2", eXtalClus2);
fChain->SetBranchAddress("ietaXtalClus2", ietaXtalClus2);
fChain->SetBranchAddress("iphiXtalClus2", iphiXtalClus2);
fChain->SetBranchAddress("izXtalClus2", &izXtalClus2);

fChain->SetBranchAddress("ptmin", &ptmin);
fChain->SetBranchAddress("ptpair", &ptpair);
fChain->SetBranchAddress("etapair", &etapair);
fChain->SetBranchAddress("s4s9min", &s4s9min);

if(pizEta==2){
 fChain->SetBranchAddress("s9s25min", &s9s25min);
}

fChain->SetBranchAddress("isolation", &isolation);
fChain->SetBranchAddress("vBeamSpot",vBeamSpot);
 
}
