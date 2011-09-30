int dataflag; 
int evtRange; 
int dataOrMC; 
int totalEntries; 
int entry; 

TChain *fChain; 

TTree *newtree; 
string dataversion;
string datasetname;

float detasctrkvtx;
float dphisctrkvtx;

float e1x5;



float photonet[nPhotonMAX];
float photonetvtx0[nPhotonMAX];
float photoncaloeta[nPhotonMAX];
float photoncalophi[nPhotonMAX];
int photoninbarrel[nPhotonMAX];

float electronecalEnergy[nElectronMAX];
float electronecalEnergyCorr[nElectronMAX];


TTree *evtInfo;

int totalNumberofEvents ; 
int preselectedNumberofEvents; 
vector<string> *datasettag; 

float EnergyScaleOffPhotonCat[4];
float EnergySmearingPhotonCat[4];

TRandom3 *rgen_;

std::map< int,std::vector<TH1*> > kFactors_;


///ID efficieny correction

std::map<std::string,TGraphAsymmErrors*> IDCorrection_; 

/* std::map<std::string, RooRealVar> m_real_var_; */
/* std::map<std::string, double> m_var_min_; */
/* std::map<std::string, double> m_var_max_; */
/* std::map<std::string,RooRealVar*> m_data_var_ptr_; */
/* std::map<std::string,RooRealVar*> m_weight_var_ptr_; */
/* std::map<std::string,RooDataSet> data_; */
/* std::map<std::string,std::string> data_obs_names_; */


 std::vector<float> ptunbal_;
std::vector<float> angleptunbaldipht_;
std::vector<float> asym_ptunbal_;

std::vector<float> ncha_;
std::vector<float> nchapt1_;


 std::vector<float> ptbal_;
  std::vector<float> thrust_;
  std::vector<float> sumpt_;
  std::vector<float> sumpt2_;
  std::vector<float> logsumpt2_;
  std::vector<float> sumawy_;
  std::vector<float> sumtwd_;
  std::vector<float> sumtrv_;
  std::vector<float> sumweight_;
  std::vector<float> ptmax_;
  std::vector<float> nchthr_;
  std::vector<float> nch_;
  std::vector<std::vector<float> > tksPt_;
  std::vector<TMatrixDSym> sphers_;
  std::vector<float> sumpr_;
  std::vector<float> spher_;
  std::vector<float> tspher_;
  std::vector<float> aplan_;
  std::vector<float> threejetC_;
  std::vector<float> fourjetD_;
  
  std::vector<TVector3> vtxP_;
  std::vector<TVector2> vtxPt_;
  std::vector<float> ptvtx_;
  std::vector<float> diphopt_; 
  std::vector<float> diPhotonPz_;
  std::vector<float> acosA_;
  std::vector<float> ptasym_;
  
  std::vector<float> ptmax3_;
  std::vector<float> ptratio_;
  std::vector<float> pzasym_;
  std::vector<float> awytwdasym_;
  std::vector<TVector2> diPhotonPt_;

  vector<int> indpresel; 
  vector<int> indaccept; 

std::vector<int> preselConv;


std::vector<int> rankprodvtx;

std::vector<int> tmvavtx;

  float ptbal;
float tspher;
float spher;
float aplan;
float threejetC;
float fourjetD;
float diphopt;
float acosA;
float ptasym;
float ptmax;
float ptmax3;
float thrust;
float ptratio;
float pzasym;
float sumwt;
float vtxpt;

float sumpt2;
float sumpt;
float awytwdasym;
 int preselectedbyconv; 
 int isgenmatched; 
 float ptunbal; 
 float angleptunbaldipht;
 float asym_ptunbal;
float ncha; 
float nchapt1; 
 
 float nchapt1f;
 
 float logsumpt2; 

float ptdiphtfirstvertex;


TMVA::Reader *reader; 
TMVA::Reader *reader2; 

TMVA::Reader *reader_eb_unconv;
TMVA::Reader *reader_eb_conv;
TMVA::Reader *reader_ee_unconv;
TMVA::Reader *reader_ee_conv;


TString methodName;

float tmva; 
float tmva2; 



int isbestvtx;
int isworstvtx;
float distvtxgen;

int testTMVA;


std::vector<std::vector<float> >* photontrkiso_recvtx_030_002_0000_10_01;
float photontrkiso_badvtx_040_002_0000_10_01[nPhotonMAX];
float photondrtotrk[nPhotonMAX];
int indvertexSelected;


vector<short> *indvertexSelected_allpairpresel;
vector<short> *indvertexSelectedtmva_allpairpresel;

////vector<vector<float> >* photontrkisoworstdr04;
///the aboved updated, for each photon, there is one number 
float photontrkisoworstdr04[nPhotonMAX];
int photonindvtxtrkisoworstdr04[nPhotonMAX];
vector<vector<float> >* photontrkisoselvtxdr03;
vector<vector<float> >* photontrkisoselvtxtmvadr03;
//vector<vector<short> >* photonindvtxtrkisoworstdr04;


int ngenphoton;
float etagenphoton[nPhotonMAX];
float phigenphoton[nPhotonMAX];
float ptgenphoton[nPhotonMAX];
float vxgenphoton[nPhotonMAX];
float vygenphoton[nPhotonMAX];
float vzgenphoton[nPhotonMAX];
int pidmomgenphoton[nPhotonMAX];
int pidmom2genphoton[nPhotonMAX];
int pidmom3genphoton[nPhotonMAX];
int statusgenphoton[nPhotonMAX];


int ngenelectron;
float etagenelectron[nElectronMAX];
float phigenelectron[nElectronMAX];
float ptgenelectron[nElectronMAX];
float vxgenelectron[nElectronMAX];
float vygenelectron[nElectronMAX];
float vzgenelectron[nElectronMAX];
int pidmomgenelectron[nElectronMAX];
int pidmom2genelectron[nElectronMAX];
int pidmom3genelectron[nElectronMAX];
int statusgenelectron[nElectronMAX];


float photon_convp[nPhotonMAX];
float photon_convpt[nPhotonMAX];
float photon_convpttrk1[nPhotonMAX];
float photon_convpttrk2[nPhotonMAX];
float photon_deltaphiconvsc[nPhotonMAX];
float photon_deltaetaconvsc[nPhotonMAX];
float photon_convrho[nPhotonMAX];
float photon_convz[nPhotonMAX];




  float isosumoet; 
  float isosumoetbad; 
  float trkisooet; 
  //float sieie;
  float hoe; 
  float r9; 
  float drtotrk; 
  int haspromptele; 
  
   ///new tree for singlePhotonTree
  float etrue; 
float etatrue; 
float phitrue;
float etaecaltrue;
float phiecaltrue;

float vxtrue;
float vytrue;
float vztrue;

float vx;
float vy;
float vz;

float ecalisodr03;
float ecalisodr04;
float hcalisodr03;
float hcalisodr04;


  float escraw; 
  float e5x5; 
  
  float e2x2;
  float e3x3;
  float e1x3;
  float e3x1;
  float e4x4;
  float e2x5;
  float e2x5right;
  float e2x5left;
  float e2x5top;
  float e2x5bottom;
  float e2x5max;
    
  float eps_escraw;

  float e3x3_e5x5;
  float eps_e5x5;
  float e5x5_escraw; 
  float emax_escraw; 
  float eleft_escraw; 
  float eright_escraw; 
  float etop_escraw; 
  float ebottom_escraw; 

float e2x5left_escraw;
float e2x5right_escraw;
float e2x5top_escraw;
float e2x5bottom_escraw;


  float emax_e5x5; 
  float eleft_e5x5; 
  float eright_e5x5; 
  float etop_e5x5; 
  float ebottom_e5x5; 

float etaele; 
float phiele;

  float eps; 
  float esc; 
float eele;

  float epht;
  float etasc; 
  float phisc; 
  float etapht;
  float phipht;
  float scetawidth; 
  float scphiwidth; 
  float sigietaieta; 
  float sigiphiiphi;
  float sigcovietaiphi;
  float emax; 
  int scnbc; 
  int scncrystal;
  float scbcfmax;
  int seedieta; 
  int seediphi; 
  float eleft; 
  float eright; 
float e2max; 
  float etop;
  float ebottom;
  float convp; 
  float convpt; 
  float convpttrk1; 
  float convpttrk2; 
  float deltaphi_convsc;
  float deltaeta_convsc; 
  
  float convrho;
  float convz; 
  
  float escrawoverconvp;
  float convleadingptoverpt;

float pileupwt;
  
  float seedietaf; 
  float seediphif; 
  float scnbcf; 
  float scncrystalf; 
float nVertexf; 


  int truephtmatched; 
  float mpair; 
  

  float etaCGap; 
  float phiCGap;
  float etaSGap; 
  float phiSGap;
  float etaMGap; 
  float phiMGap;

float xCGap; 
float xMGap; 
float xSGap; 
float yCGap; 
float yMGap; 
float ySGap; 


int preselcut;

int testPhotonCorr;
int testVtxID;
int testPhotonID;
int trueelematched;
int npvbx0;
float npvbxav;



int nWarning;
int eleclass;
float elefbrem;
int elenbrem;
float sigiphi_sigieta;
float scsigphi_sigeta;
int trainTarget;
int pidphtmom; 
float trkisodr03vtxbestGenMatched;
float trkisodr04vtxWorst;

bool  doEscaleSyst, doEresolSyst, doPhotonIdEffSyst, doVtxEffSyst, doR9Syst, doTriggerEffSyst, doKFactorSyst;
bool doPhotonEffSys;


std::map<std::string,int>::iterator it_sys;

int nSystSteps; 
float systRange;

float photonenergy0[nPhotonMAX];

float et1cut;
float et2cut;
int nPhotonPreSelected;
float add_sigmaEoverE[4];
float addErr_sigmaEoverE[4];
