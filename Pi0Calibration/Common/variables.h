    
TChain *fChain; 
int dataflag; 
int totalEntries; 
int entry; 
int doBarrel; 
int doEndcap; 
int pizEta;


// Declaration of leaf types
Int_t           lumiBlock;
Int_t           runNumber;
Int_t           evtNumber;
Int_t           evtTime;
vector<unsigned short> *l1bitFired;   // List of L1 bit fired 
Float_t         mpair; /// invariant mass of the pi0 candidate 
Float_t         ptpair; // pt of the pi0 candidate 
Float_t         etapair; ///eta of the pi0 candidate 
Float_t         ptmin; // minimal pt of two photons 
Float_t         isolation; // isolation variable 
Float_t         vBeamSpot[3]; // beam spot x,y, z 
Float_t         s4s9min; /// minimal s4/s9 of two photons 
Float_t         s9s25min; /// minimal s9/s25 of two photons 

Float_t         xClus1; /// X position of 3x3 cluster 1 
Float_t         yClus1;
Float_t         zClus1;
Float_t         xClus2;  //X position of cluster 2 
Float_t         yClus2;
Float_t         zClus2;
Int_t           nxtClus1; ///Number of crystal in cluster 1  ( maximum 9 / 25 for eta ) 
Float_t         eXtalClus1[25];   //[nxtClus1]  // energy of each crystal 
Float_t         laserCorrXtalClus1[25];   //[nxtClus1] /// laserCorrection value 
Int_t           ietaXtalClus1[25];   //[nxtClus1]  /// ieta 
Int_t           iphiXtalClus1[25];   //[nxtClus1]  /// iphi 
Float_t         tXtalClus1[25];   //[nxtClus1]  // recHit time 
Int_t           nxtClus2;
Float_t         laserCorrXtalClus2[25];   //[nxtClus2]
Float_t         eXtalClus2[25];   //[nxtClus2]
Int_t           ietaXtalClus2[25];   //[nxtClus2]
Int_t           iphiXtalClus2[25];   //[nxtClus2]
Float_t         tXtalClus2[25];   //[nxtClus2]
Int_t izXtalClus1; 
Int_t izXtalClus2; 


float xEBAll[170][360];  ///eta -85, -1, 1, 86, [0,169].  phi, 1,360, [0,359];
float yEBAll[170][360]; 
float zEBAll[170][360]; 
float etaEBAll[170][360]; 
float phiEBAll[170][360]; 

float dxEBAll[170][360];  ///eta -85, -1, 1, 86, [0,169].  phi, 1,360, [0,359];
float dyEBAll[170][360]; 
float dzEBAll[170][360]; 

float xEEAll[2][101][101];
float yEEAll[2][101][101];
float zEEAll[2][101][101];
float etaEEAll[2][101][101];
float phiEEAll[2][101][101];


float dxEEAll[2][101][101];
float dyEEAll[2][101][101];
float dzEEAll[2][101][101];

int doPizEta; 
ofstream txtout; 


TTree *pizTree; 
TTree *etaTree; 

double corrfactorEtatb[170];
double corrfactorEtaco[170];
double corrfactorEta[170];
double corrfactorPhi[360];
double corrfactoriEtaiPhi[170][360];
int flagiEtaiPhi[170][360];
double corrfactorSM[36];
double corrfactorIetaSM[38][85];

int stepc;  
int iter; 
int evtRange;

const int kEndcEtaRings = 39;
float eta_ring_ee[kEndcEtaRings];
float etaBoundary_ee[kEndcEtaRings+1];
int nxtal_ring_ee[kEndcEtaRings];
double corrfactorEtaRings[2][kEndcEtaRings]; 
double corrfactoriZiXiY[2][101][101];
float etaBoundary_eezside[2][kEndcEtaRings+1];
float eta_ring_eezside[2][kEndcEtaRings];
float infoESX[2][8];
float infoESY[2][8];
double peakwidthEtaRings[2][kEndcEtaRings];
float sigma_sideRing[2][40];
float mean_side[2] = {0.1264,0.1244};
float sigma_side[2] = {0.02,0.02};
int validRecHitEndCap[2][101][101];


int flag_ietaiphi[170][360];
int ndead_ietaiphi[170][360];
int ndeadcorner_ietaiphi[170][360];
int ndeadside_ietaiphi[170][360];
int ndeadflag_ietaiphi[170][360];
double corrfactorDead[20] ;

/* ///some correction to crystals */
/* double preCorrectionIetaIphi[170][360]; */
/* double corrfactorIetaSM[38][85]; */
/* double corrfactorPhiSide[2][360]; */
/* bool is2010RunB; */
/* bool is2010RunA; */
/* bool is2011RunA; */

int nCounted[170][360][100];
int nCountedEE[2][101][101][100];

float sigmaMass;
float meanMass;

string workingDirectory;


///for endcaps
int nMaxRingIC;
