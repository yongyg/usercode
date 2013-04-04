
void setbranchaddress_g4step(){

fChain->SetBranchAddress("ng4EB", &ng4EB);
fChain->SetBranchAddress("ietag4EB", ietag4EB);
fChain->SetBranchAddress("iphig4EB", iphig4EB);
fChain->SetBranchAddress("esumg4EB", esumg4EB);
fChain->SetBranchAddress("tming4EB", tming4EB);
fChain->SetBranchAddress("xtming4EB", xtming4EB);
fChain->SetBranchAddress("ytming4EB", ytming4EB);
fChain->SetBranchAddress("ztming4EB", ztming4EB);


fChain->SetBranchAddress("ng4EE", &ng4EE);
fChain->SetBranchAddress("ixg4EE", ixg4EE);
fChain->SetBranchAddress("iyg4EE", iyg4EE);
fChain->SetBranchAddress("izg4EE", izg4EE);
fChain->SetBranchAddress("esumg4EE", esumg4EE);
fChain->SetBranchAddress("tming4EE", tming4EE);
fChain->SetBranchAddress("xtming4EE", xtming4EE);
fChain->SetBranchAddress("ytming4EE", ytming4EE);
fChain->SetBranchAddress("ztming4EE", ztming4EE);


fChain->SetBranchAddress("eg4EB", &eg4EB);
fChain->SetBranchAddress("tg4EB", &tg4EB);
fChain->SetBranchAddress("idg4EB", &idg4EB);
fChain->SetBranchAddress("pidg4EB", &pidg4EB);
fChain->SetBranchAddress("parentidg4EB", &parentidg4EB);

fChain->SetBranchAddress("postxg4EB", &postxg4EB);
fChain->SetBranchAddress("postyg4EB", &postyg4EB);
fChain->SetBranchAddress("postzg4EB", &postzg4EB);
fChain->SetBranchAddress("prexg4EB", &prexg4EB);
fChain->SetBranchAddress("preyg4EB", &preyg4EB);
fChain->SetBranchAddress("prezg4EB", &prezg4EB);
fChain->SetBranchAddress("pretg4EB", &pretg4EB);
fChain->SetBranchAddress("preeg4EB", &preeg4EB);
fChain->SetBranchAddress("enterg4EB", &enterg4EB);
fChain->SetBranchAddress("leaveg4EB", &leaveg4EB);


fChain->SetBranchAddress("eg4EE", &eg4EE);
fChain->SetBranchAddress("tg4EE", &tg4EE);
fChain->SetBranchAddress("idg4EE", &idg4EE);
fChain->SetBranchAddress("pidg4EE", &pidg4EE);
fChain->SetBranchAddress("parentidg4EE", &parentidg4EE);

fChain->SetBranchAddress("postxg4EE", &postxg4EE);
fChain->SetBranchAddress("postyg4EE", &postyg4EE);
fChain->SetBranchAddress("postzg4EE", &postzg4EE);
fChain->SetBranchAddress("prexg4EE", &prexg4EE);
fChain->SetBranchAddress("preyg4EE", &preyg4EE);
fChain->SetBranchAddress("prezg4EE", &prezg4EE);
fChain->SetBranchAddress("pretg4EE", &pretg4EE);
fChain->SetBranchAddress("preeg4EE", &preeg4EE);
fChain->SetBranchAddress("enterg4EE", &enterg4EE);
fChain->SetBranchAddress("leaveg4EE", &leaveg4EE);

}

void setbranchaddress_g4sim(){

  
  fChain1->SetBranchAddress("nsimEB", &nsimEB);
  fChain1->SetBranchAddress("ietasimEB", ietasimEB);
  fChain1->SetBranchAddress("iphisimEB", iphisimEB);
  fChain1->SetBranchAddress("esumsimEB", esumsimEB);
  fChain1->SetBranchAddress("tminsimEB", tminsimEB);
  
  fChain1->SetBranchAddress("esimEB", &esimEB);
  fChain1->SetBranchAddress("tsimEB", &tsimEB);
  fChain1->SetBranchAddress("bsimEB", &bsimEB);
  fChain1->SetBranchAddress("nsimEE", &nsimEE);
  fChain1->SetBranchAddress("ixsimEE", ixsimEE);
  fChain1->SetBranchAddress("iysimEE", iysimEE);
  fChain1->SetBranchAddress("izsimEE", izsimEE);
  fChain1->SetBranchAddress("esumsimEE", esumsimEE);
  fChain1->SetBranchAddress("tminsimEE", tminsimEE);

  fChain1->SetBranchAddress("esimEE", &esimEE);
  fChain1->SetBranchAddress("tsimEE", &tsimEE);
  fChain1->SetBranchAddress("bsimEE", &bsimEE);

  fChain1->SetBranchAddress("lumiBlock", &lumiBlock);
  fChain1->SetBranchAddress("runNumber", &runNumber);
  fChain1->SetBranchAddress("evtNumber", &evtNumber);
  fChain1->SetBranchAddress("bunchX", &bunchX);
  fChain1->SetBranchAddress("orbitNumber", &orbitNumber);
  fChain1->SetBranchAddress("evtTime", &evtTime);
  fChain1->SetBranchAddress("isRealData", &isRealData);

  
  fChain1->SetBranchAddress("nGenPht", &nGenPht);
  fChain1->SetBranchAddress("etaGenPht", etaGenPht);
  fChain1->SetBranchAddress("phiGenPht", phiGenPht);
  fChain1->SetBranchAddress("ptGenPht", ptGenPht);
  fChain1->SetBranchAddress("vxGenPht", vxGenPht);
  fChain1->SetBranchAddress("vyGenPht", vyGenPht);
  fChain1->SetBranchAddress("vzGenPht", vzGenPht);
  fChain1->SetBranchAddress("pidmomGenPht", pidmomGenPht);
  fChain1->SetBranchAddress("pidmom2GenPht", pidmom2GenPht);
  fChain1->SetBranchAddress("indmom2GenPht", indmom2GenPht);
  fChain1->SetBranchAddress("pidmom3GenPht", pidmom3GenPht);
  fChain1->SetBranchAddress("statusGenPht", statusGenPht);
  fChain1->SetBranchAddress("convGenPht", convGenPht);
  
  
}
