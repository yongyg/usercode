This is the analyzer to make trees with selection for pi0/eta->gg

4 main trees saved 
pizSelb
pizSele
etaSelb
etaSele

one small-size tree saved as well perl lumiBlock for beamSpot information
evtInfobs 

What saved is per diphoton pair satisfying selection cuts, not per event. 


while one can understand the meaning of the variables easily( by looking at the source code sr/RecoAnalyzer.cc) , below are a few most important ones, 


float mpair; /// invariant mass of two 3x3 clusters ( calculated assuming (0,0,0) vertex posistion) 
float ptpair;  /// pt pair of pi0/eta candidate

float infoESX[2][8];  //preshower information for each 3x3 cluster 0 for the first one, 1 for the 2nd one. 

infoESX[][0]; --> Energy of sum of preshower clusters in x plane 
infoESX[][1/2/3]; --> weighted position (x,y,z) of the preshower cluster in x plane

infoESX[][4]; --> Energy of maximum preshower cluster
infoESX[][5/6/7]; --> position (x,y,z) of the preshower cluster wiht maximum energy in x plane

float xClus1; /// x position of the 1st 3x3 cluster
float yClus1; /// y position of the 1st 3x3 cluster
float zClus1; /// z position of the 1st 3x3 cluster

float xClus2; /// x position of the 2nd 3x3 cluster
float yClus2; /// y position of the 2nd 3x3 cluster
float zClus2; /// z position of the 2nd 3x3 cluster

int nxtClus1; /// number of crystals for 1st cluster

int eXtalClus1[nxtClus1]; // energy of each crystal for 1st cluster
int ietaXtalClus1[nxtClus1]; // ieta ( for barrel) of each crystal for 1st cluster, ix ( for endcap) 
int iphiXtalClus1[nxtClus1]; // iphi ( for barrel) of each crystal for 1st cluster, iy ( for endcap) 
int izXtalClus1; /// iz , for endcap only

int nxtClus2 ; //  number of crystals for 2nd cluster. others similar meaning as above


two python files. 

testmc.py  --> run GEN-SIM-RECO 

testdata.py --> run AlCaRaw Data with laser tag 



