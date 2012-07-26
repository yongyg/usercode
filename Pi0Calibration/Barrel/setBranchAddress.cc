
void setBranchAddress(){

fChain->SetBranchAddress("mpair", &mpair);
fChain->SetBranchAddress("lumiBlock", &lumiBlock);
fChain->SetBranchAddress("runNumber", &runNumber);
fChain->SetBranchAddress("evtNumber", &evtNumber);
fChain->SetBranchAddress("nxtClus1", &nxtClus1);
fChain->SetBranchAddress("eXtalClus1", eXtalClus1);
fChain->SetBranchAddress("ietaXtalClus1", ietaXtalClus1);
fChain->SetBranchAddress("iphiXtalClus1", iphiXtalClus1);
//fChain->SetBranchAddress("tXtalClus1", tXtalClus1);
fChain->SetBranchAddress("nxtClus2", &nxtClus2);
fChain->SetBranchAddress("eXtalClus2", eXtalClus2);
fChain->SetBranchAddress("ietaXtalClus2", ietaXtalClus2);
fChain->SetBranchAddress("iphiXtalClus2", iphiXtalClus2);
fChain->SetBranchAddress("vBeamSpot",vBeamSpot);

fChain->SetBranchAddress("xClus1",&xClus1);
fChain->SetBranchAddress("yClus1",&yClus1);
fChain->SetBranchAddress("zClus1",&zClus1);

fChain->SetBranchAddress("xClus2",&xClus2);
fChain->SetBranchAddress("yClus2",&yClus2);
fChain->SetBranchAddress("zClus2",&zClus2);


}
