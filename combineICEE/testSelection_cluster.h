    
TChain *fChain; 

int dataflag; 

int totalEntries; 


int totalEntriesL1; 


int entry; 


int doBarrel; 
int doEndcap; 

int trigger; 

  
  float ptminCut; 
float ptpairCut; 
float s4s9Cut; 
float s9s25Cut; 
float isoCut; 


 // L1 objects
      static const int MAXL1OBJ = 20;
      //     static const int MAXL1BITS = 128;
      // L1 EM isolated objects
      int nL1EMIso;
      float L1EMIso_e[MAXL1OBJ];
      float L1EMIso_et[MAXL1OBJ];
      float L1EMIso_eta[MAXL1OBJ];
      float L1EMIso_phi[MAXL1OBJ];
      // L1 EM non-isolated objects
      int nL1EMnonIso;
      float L1EMnonIso_e[MAXL1OBJ];
      float L1EMnonIso_et[MAXL1OBJ];
      float L1EMnonIso_eta[MAXL1OBJ];
      float L1EMnonIso_phi[MAXL1OBJ];

      int bunchX; 
      int orbitNumber; 
      int evtTime; 
      
      int runNumber; 
      int evtNumber; 
      int lumiBlock; 
    static const int MAXL1bits = 500;
      static const int MAXHLTbits = 500;
      int nL1bits;
      int L1bits[MAXL1bits];

int nL1Alca; 
int L1Alca[MAXL1bits];

      
      int nL1bitsTech;
      int L1bitsTech[MAXL1bits];

int psL1; 
int L1Cut; 

int dpsL1DecisionWord8E29v3; 
int dpsL1DecisionWord8E29v3_1; 
int dpsL1DecisionWord8E29v3_2; 
int dpsL1DecisionWord8E29v3_3; 
int dpsL1DecisionWord8E29v3_4; 
int dpsL1DecisionWord8E29v3_5; 
int dpsL1DecisionWord8E29v4; 

int dpsL1bits8E29v3[MAXL1bits];

int dpsL1bits8E29v3_1[MAXL1bits];
int dpsL1bits8E29v3_2[MAXL1bits];
int dpsL1bits8E29v3_3[MAXL1bits];
int dpsL1bits8E29v3_4[MAXL1bits];
int dpsL1bits8E29v3_5[MAXL1bits];


int dpsL1bits8E29v4[MAXL1bits];

int iCountL1NoPrescale8E29v3[500]; ///2009.06.27

int iCountL1NoPrescale8E29v3_1[500]; ///2009.06.27 EG PS = 10 
int iCountL1NoPrescale8E29v3_2[500]; ///2009.06.27 EG PS = 20 
int iCountL1NoPrescale8E29v3_3[500]; ///2009.06.27 EG PS = 50 
int iCountL1NoPrescale8E29v3_4[500]; ///2009.06.27 EG PS = 2 
int iCountL1NoPrescale8E29v3_5[500]; ///2009.06.27 EG PS = 4 


int iCountL1NoPrescale8E29v4[500]; ///2009.06.27, prescale 100 for EG5, Jet6,10,20


      float xOffPVwithBS;
      float yOffPVwithBS;
      float zOffPVwithBS;


      static const int MAXGenPIZ = 1000; 

      int nGenpi0; 
      int nGeneta; 
      float etaGenpi0[MAXGenPIZ];
      float phiGenpi0[MAXGenPIZ];
      float etGenpi0[MAXGenPIZ];
      float eGenpi0[MAXGenPIZ];
      float mGenpi0[MAXGenPIZ];
      float dr2phGenpi0[MAXGenPIZ];
      float dr2phGeneta[MAXGenPIZ];
      float etaGeneta[MAXGenPIZ];
      float phiGeneta[MAXGenPIZ];
      float etGeneta[MAXGenPIZ];
      float eGeneta[MAXGenPIZ];
      float ePhtGenpi0[MAXGenPIZ][2];
      float etaPhtGenpi0[MAXGenPIZ][2];
      float phiPhtGenpi0[MAXGenPIZ][2];
      float ePhtGeneta[MAXGenPIZ][2];
      float etaPhtGeneta[MAXGenPIZ][2];
      float phiPhtGeneta[MAXGenPIZ][2];
      float vtxGenpi0[MAXGenPIZ][3];
      float vtxPhtGenpi0[MAXGenPIZ][3];
      float vtxGeneta[MAXGenPIZ][3];
      float vtxPhtGeneta[MAXGenPIZ][3];


      int pidMomGenpi0[MAXGenPIZ];
      int barcodeMomGenpi0[MAXGenPIZ];
      
      int pidMomGeneta[MAXGenPIZ];
      int barcodeMomGeneta[MAXGenPIZ];

      int isConvPhtGeneta[MAXGenPIZ][2];
      float convPht1Geneta[MAXGenPIZ][6]; ///e, eta, phi, e, eta, phi
      float convPht2Geneta[MAXGenPIZ][6]; ///e, eta, phi, e, eta, phi

      int isConvPhtGenpi0[MAXGenPIZ][2];
      float convPht1Genpi0[MAXGenPIZ][6]; ///e, eta, phi, e, eta, phi
      float convPht2Genpi0[MAXGenPIZ][6]; ///e, eta, phi, e, eta, phi
      float convVtxPhtGenpi0[MAXGenPIZ][6];
      float convVtxPhtGeneta[MAXGenPIZ][6];
      
   static const int MAX3x3ClusEB = 2000;
      int n3x3ClusEB;
      float e3x3ClusEB[MAX3x3ClusEB];
      float eta3x3ClusEB[MAX3x3ClusEB];
      float phi3x3ClusEB[MAX3x3ClusEB];

      ///tmp by simplelog
      float leta3x3ClusEB[MAX3x3ClusEB];
      float lphi3x3ClusEB[MAX3x3ClusEB];
      

      float x3x3ClusEB[MAX3x3ClusEB];
      float y3x3ClusEB[MAX3x3ClusEB];
      float z3x3ClusEB[MAX3x3ClusEB];

      int nXt3x3ClusEB[MAX3x3ClusEB];
      float eXt3x3ClusEB[MAX3x3ClusEB][9];
      int ietaXt3x3ClusEB[MAX3x3ClusEB][9]; ///index from -85,0,84; ,
      int iphiXt3x3ClusEB[MAX3x3ClusEB][9]; ///phi, 0,359
      float etaXt3x3ClusEB[MAX3x3ClusEB][9];
      float phiXt3x3ClusEB[MAX3x3ClusEB][9];
      
      
      float xXt3x3ClusEB[MAX3x3ClusEB][9];
      float yXt3x3ClusEB[MAX3x3ClusEB][9];
      float zXt3x3ClusEB[MAX3x3ClusEB][9];
      float S43x3ClusEB[MAX3x3ClusEB];
      float S63x3ClusEB[MAX3x3ClusEB];
      float S253x3ClusEB[MAX3x3ClusEB];
      
      float s4s93x3ClusEB[MAX3x3ClusEB];
      float s6s93x3ClusEB[MAX3x3ClusEB];
      float s9s253x3ClusEB[MAX3x3ClusEB];
      
      
    
      static const int MAX3x3ClusEE = 2000;
      int n3x3ClusEE;
      float e3x3ClusEE[MAX3x3ClusEE];
      float eta3x3ClusEE[MAX3x3ClusEE];
      float phi3x3ClusEE[MAX3x3ClusEE];

      ///tmp by simplelog
      float leta3x3ClusEE[MAX3x3ClusEE];
      float lphi3x3ClusEE[MAX3x3ClusEE];
      float tXt3x3ClusEB[MAX3x3ClusEB][9];  

      float x3x3ClusEE[MAX3x3ClusEE];
      float y3x3ClusEE[MAX3x3ClusEE];
      float z3x3ClusEE[MAX3x3ClusEE];

      int nXt3x3ClusEE[MAX3x3ClusEE];
      float eXt3x3ClusEE[MAX3x3ClusEE][9];
      int ixXt3x3ClusEE[MAX3x3ClusEE][9]; ///index from -85,0,84; ,
      int iyXt3x3ClusEE[MAX3x3ClusEE][9]; ///phi, 0,359
      int izXt3x3ClusEE[MAX3x3ClusEE][9]; ///phi, 0,359
      
      float etaXt3x3ClusEE[MAX3x3ClusEE][9];
      float phiXt3x3ClusEE[MAX3x3ClusEE][9];
     float tXt3x3ClusEE[MAX3x3ClusEE][9];   
      
      float xXt3x3ClusEE[MAX3x3ClusEE][9];
      float yXt3x3ClusEE[MAX3x3ClusEE][9];
      float zXt3x3ClusEE[MAX3x3ClusEE][9];
      float S43x3ClusEE[MAX3x3ClusEE];
      float S63x3ClusEE[MAX3x3ClusEE];
      float S253x3ClusEE[MAX3x3ClusEE];
      
      float s4s93x3ClusEE[MAX3x3ClusEE];
      float s6s93x3ClusEE[MAX3x3ClusEE];
      float s9s253x3ClusEE[MAX3x3ClusEE];
      

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
int nCrystalCut;
int doPizEta; 
bool useGenVtxforRecoGenMatch;
int dataOrMC;

vector<int> L1bitsUsed; 
vector<int> L1bitsTechUsed; 

vector<int> L1bitsNotUsed; 
vector<int> L1bitsTechNotUsed; 

ofstream txtout; 
int useRecVtxforReco ;
int vetoRunList;
int evtNotUse;
int applyBPTX;
int procID;
int vetoBeamSrape;

  int nPxlHits ; 
      float clusVtxQual; 
      
      int    nHfTowersP     ;
      int    nHfTowersN     ; 

      float sumHfEsubEpPlus; 
      float sumHfEsubEpMinus; 
      
      float sumHfEaddEpPlus; 
      float sumHfEaddEpMinus; 
      
      
      int phyDeclared; 
      
      float highPurityTrackFrac; 
 static const int MAXTRK = 5000;
      int nTracks;
      float pxTracks[MAXTRK];
      float pyTracks[MAXTRK];
      float pzTracks[MAXTRK];
      float pTracks[MAXTRK];
      float ptTracks[MAXTRK];
      float etaTracks[MAXTRK];
      float phiTracks[MAXTRK];
      int pidTracks[MAXTRK];
      int chargeTracks[MAXTRK]; 
      float vxTracks[MAXTRK];
      float vyTracks[MAXTRK];
      float vzTracks[MAXTRK];
      int nhitsTracks[MAXTRK];
      float nChi2Tracks[MAXTRK];
      int algoTracks[MAXTRK];
      int ndofTracks[MAXTRK];
      int nValidhitsTracks[MAXTRK];
      int nValidpixelhitsTracks[MAXTRK];
      int nValidstriphitsTracks[MAXTRK];
      int qualityFlagTracks[MAXTRK];
      
      float vBeamSpot[3];
int hfCoindence;
int vetoHighHitPxl;
int recalcs4s9;
int applyEnCorr;
int corrSeedEnScale;
float eSeed;
int vetoSpike;
int applyPhyDeclare;
int validVtx;
int NtrkOffPVwithBS;
float ptminCut_ee;
float s4s9Cut_ee;
float s9s25Cut_ee;
int passEG3;
int passEG2;
int passEG4;
int passEG5;




///make above in one
///first one mean ncut
///at most 10 cuts, 100 steps each.
float variable[10];
const int ncuts = 5; 
///min max cut 
float init_cutValue[10][2];
int nstepCut[10];
float cutValues[10][100];

const int nstepsMAX = 31; 
//at most 20 steps.
int nSigPassed[nstepsMAX][nstepsMAX][nstepsMAX][nstepsMAX][nstepsMAX];
int nBkgPassed[nstepsMAX][nstepsMAX][nstepsMAX][nstepsMAX][nstepsMAX];


int IsSignal; 
int IsBkgrd; 
int indexBestCuts1[10];
int indexBestCuts2[10];

float isolation; 


///for differetn fom
int indexBestCuts[10][10];

//SB>0.5,0.6...
int indexBestCuts_SB[10][10][10];


///for different fom and cuts
TH1F *hh_fom[10][10]; 


  Float_t         mpair;
   Float_t         ptpair;
   Int_t           pizMatched;
   Int_t           etaMatched;
   Int_t           pizMatchedv1;
   Int_t           etaMatchedv1;
   Float_t         ptmin;
   Float_t         s4s9min;
   Float_t         s9s25min;
   Float_t         iso;
   Float_t         isoeta;
   Float_t         isov1;
   Float_t         isoetav1;
   Float_t         isov2;
   Float_t         isoetav2;
   Int_t           passL1;
   Int_t           nxt9[2];
   Int_t           seediEta[2];
   Int_t           seediPhi[2];

int doEcut; 

int testAbsoluteIso; 
float ptminOpt;
int izXt3x3ClusEB[MAX3x3ClusEB][9];
Int_t seediZ[2];
Float_t etapair;
int isoflag;
int ptCutSetee;
int testpizwind;
float eSeed_ee;
int fullunpack;


TTree *pizTree; 


TTree *etaTree; 

      float mpair_3x3; 
int nxtClus1; 
      float eXtalClus1[25];
      int ietaXtalClus1[25];
      int iphiXtalClus1[25];
      float tXtalClus1[25];
      int fXtalClus1[25];
      int nxtClus2; 
      float eXtalClus2[25];
      int ietaXtalClus2[25];
      int iphiXtalClus2[25];
      float tXtalClus2[25];
      int fXtalClus2[25];
int savePizTree;


int passBSC40or41;
int passBSC34;
int passHF;
int passBeamScrape;
int passVtx;
int passBeamHalo;
int passBSC0;
int passHighHitPxl;
int rmOverlap;
float mpairv1;
float phipair; 
float geninfo[6];
int start_entry;
int end_entry;
int calibMethod;
int applyCorrEta;
int applyCorrPhi;
int applyCorrDead; 


double corrfactorEta[170];
double corrfactorEtatb[170];
double corrfactorEtaco[170];
double corrfactorPhi[360];
double corrfactoriEtaiPhi[170][360];

///pre-calib ( all 1)  or phi-sym or BS 
double preCalibiEtaiPhi[170][360];


int iterStep;

int cutset;
int evtRange;
int flagiEtaiPhi[170][360];
int usePreCalib;
int useVtx;
float randGausiEtaiPhi[170][360];

const int kEndcEtaRings = 39;
float eta_ring_ee[kEndcEtaRings];
float etaBoundary_ee[kEndcEtaRings+1];
int nxtal_ring_ee[kEndcEtaRings];
float eta_ring_eezside[2][kEndcEtaRings];
float etaBoundary_eezside[2][kEndcEtaRings+1];


///calibration constant fro endcap
double corrfactorEtaRings[2][kEndcEtaRings]; 
double corrfactoriZiXiY[2][101][101];

int izXtalClus1[25];
int izXtalClus2[25];

double corrfactorSM[36];
int applyCorrSM;
TChain *fChainAlca;

double corrfactorSM_final[36];



float Cpi0Corrv1[170][360];
float Cpi0v1[170][360];

float Ccombv1[170][360];
float CcombCorrv1[170][360];

float Cphi[170][360];
float CBS[170][360];

float CphiCorr[170][360];
float CBSCorr[170][360];
int flag_ietaiphi[170][360];
int ndead_ietaiphi[170][360];
int ndeadcorner_ietaiphi[170][360];
int ndeadside_ietaiphi[170][360];
int ndeadflag_ietaiphi[170][360];
double corrfactorDead[20] ;
int evtNotUsed;

///some correction to crystals
double preCorrectionIetaIphi[170][360];
double corrfactorIetaSM[38][85];
double corrfactorPhiSide[2][360];
int validRecHitEndCap[2][101][101];
int flag_endcap[2][101][101];
int ndeadcorner_endcap[2][101][101];
int ndeadside_endcap[2][101][101];
int ndeadflag_endcap[2][101][101];
