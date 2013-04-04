Int_t           ng4EB;
Int_t           ietag4EB[61200];   //[ng4EB]
Int_t           iphig4EB[61200];   //[ng4EB]
Float_t           esumg4EB[61200];   //[ng4EB]
Float_t           tming4EB[61200];   //[ng4EB]
Float_t           xtming4EB[61200];   //[ng4EB]
Float_t           ytming4EB[61200];   //[ng4EB]
Float_t           ztming4EB[61200];   //[ng4EB]

Int_t           ng4EE;
Int_t           ixg4EE[14648];   //[ng4EE]
Int_t           iyg4EE[14648];   //[ng4EE]
Int_t           izg4EE[14648];   //[ng4EE]
Float_t         esumg4EE[14648];   //[ng4EE]
Float_t         tming4EE[14648];   //[ng4EE]
Float_t         xtming4EE[14648];   //[ng4EE]
Float_t         ytming4EE[14648];   //[ng4EE]
Float_t         ztming4EE[14648];   //[ng4EE]



vector<vector<float> > *eg4EB;
vector<vector<float> > *tg4EB;
vector<vector<int> > *idg4EB;
vector<vector<int> > *pidg4EB;
vector<vector<int> > *parentidg4EB;

vector<vector<float> > *postxg4EB;
vector<vector<float> > *postyg4EB;
vector<vector<float> > *postzg4EB;

vector<vector<float> > *prexg4EB;
vector<vector<float> > *preyg4EB;
vector<vector<float> > *prezg4EB;


vector<vector<float> > *pretg4EB;


vector<vector<float> > *preeg4EB;
vector<vector<int> > *enterg4EB;
vector<vector<int> > *leaveg4EB;


vector<vector<float> > *eg4EE;
vector<vector<float> > *tg4EE;
vector<vector<int> > *idg4EE;
vector<vector<int> > *pidg4EE;
vector<vector<int> > *parentidg4EE;

vector<vector<float> > *postxg4EE;
vector<vector<float> > *postyg4EE;
vector<vector<float> > *postzg4EE;

vector<vector<float> > *prexg4EE;
vector<vector<float> > *preyg4EE;
vector<vector<float> > *prezg4EE;

vector<vector<float> > *pretg4EE;
vector<vector<float> > *preeg4EE;
vector<vector<int> > *enterg4EE;
vector<vector<int> > *leaveg4EE;


//Event info
int lumiBlock;
int runNumber;
int evtNumber;
int bunchX;
int orbitNumber;
int evtTime;
int isRealData;


std::map<int,std::map<int, std::vector<float> > > eEBsim;
std::map<int,std::map<int, std::vector<float> > > tEBsim;
std::map<int,std::map<int, std::vector<float> > > bEBsim;

std::map<int,std::map<int, std::vector<float> > > eEEmsim;
std::map<int,std::map<int, std::vector<float> > > tEEmsim;
std::map<int,std::map<int, std::vector<float> > > bEEmsim;

std::map<int,std::map<int, std::vector<float> > > eEEpsim;
std::map<int,std::map<int, std::vector<float> > > tEEpsim;
std::map<int,std::map<int, std::vector<float> > > bEEpsim;


int nsimEB;
int ietasimEB[61200];
int iphisimEB[61200];
float esumsimEB[61200];
float tminsimEB[61200];


std::vector<std::vector<float> > *esimEB;
std::vector<std::vector<float> > *tsimEB;
std::vector<std::vector<float> > *bsimEB;

int nsimEE;
int ixsimEE[14648];
int iysimEE[14648];
int izsimEE[14648];
float esumsimEE[14648];
float tminsimEE[14648];

std::vector<std::vector<float> > *esimEE;
std::vector<std::vector<float> > *tsimEE;
std::vector<std::vector<float> > *bsimEE;


static const int MAXGenSaved = 1000;
//gen-leve phton
int nGenPht;
float etaGenPht[MAXGenSaved];
float phiGenPht[MAXGenSaved];
float ptGenPht[MAXGenSaved];
int pidmomGenPht[MAXGenSaved];
int pidmom2GenPht[MAXGenSaved];
int indmom2GenPht[MAXGenSaved];
int pidmom3GenPht[MAXGenSaved];
int statusGenPht[MAXGenSaved];
float vxGenPht[MAXGenSaved];
float vyGenPht[MAXGenSaved];
float vzGenPht[MAXGenSaved];
int convGenPht[MAXGenSaved];
