///////////////////////////////////////////////////////////////////////////////
// File: CaloSD.cc
// Description: Sensitive Detector class for calorimeters
///////////////////////////////////////////////////////////////////////////////

#include "SimG4CMS/Calo/interface/CaloSD.h"
#include "SimDataFormats/SimHitMaker/interface/CaloSlaveSD.h"
#include "SimG4Core/Notification/interface/TrackInformation.h"
#include "SimG4Core/Application/interface/EventAction.h"

#include "G4EventManager.hh"
#include "G4SDManager.hh"
#include "G4Step.hh"
#include "G4Track.hh"
#include "G4VProcess.hh"
#include "G4GFlashSpot.hh"
#include "G4ParticleTable.hh"


#include "DataFormats/DetId/interface/DetId.h"
#include "DataFormats/EcalDetId/interface/EBDetId.h"
#include "DataFormats/EcalDetId/interface/EEDetId.h"
#include "Geometry/CaloGeometry/interface/CaloGenericDetId.h"
#include "Geometry/CaloGeometry/interface/CaloSubdetectorGeometry.h"
#include "Geometry/CaloGeometry/interface/CaloCellGeometry.h"



TFile *g4f;
TTree *g4t;

bool saveTree; 

G4String csdname;

using namespace std;

int ncounter1;
int ncounter2;
std::map<int,std::map<int, std::vector<float> > > estepg4EB;
std::map<int,std::map<int, std::vector<float> > > tstepg4EB;
std::map<int,std::map<int, std::vector<int> > > idstepg4EB;
std::map<int,std::map<int, std::vector<int> > > pidstepg4EB;
std::map<int,std::map<int, std::vector<int> > > parentidstepg4EB;

std::map<int,std::map<int, std::vector<int> > > enterstepg4EB;
std::map<int,std::map<int, std::vector<int> > > leavestepg4EB;

std::map<int,std::map<int, std::vector<float> > > postxstepg4EB;
std::map<int,std::map<int, std::vector<float> > > postystepg4EB;
std::map<int,std::map<int, std::vector<float> > > postzstepg4EB;

std::map<int,std::map<int, std::vector<float> > > pretstepg4EB;
std::map<int,std::map<int, std::vector<float> > > prexstepg4EB;
std::map<int,std::map<int, std::vector<float> > > preystepg4EB;
std::map<int,std::map<int, std::vector<float> > > prezstepg4EB;

std::map<int,std::map<int, std::vector<float> > > preestepg4EB;

//std::map<int,std::map<int, std::vector<float> > > prepxstepg4EB;
//std::map<int,std::map<int, std::vector<float> > > prepystepg4EB;
//std::map<int,std::map<int, std::vector<float> > > prepzstepg4EB;
std::map<int,std::map<int, std::vector<int> > > prechazstepg4EB;
std::map<int,std::map<int, std::vector<float> > > premzstepg4EB;


std::map<int,std::map<int, std::vector<float> > > estepg4EEm;
std::map<int,std::map<int, std::vector<float> > > tstepg4EEm;
std::map<int,std::map<int, std::vector<float> > > estepg4EEp;
std::map<int,std::map<int, std::vector<float> > > tstepg4EEp;
std::map<int,std::map<int, std::vector<int> > > idstepg4EEm;
std::map<int,std::map<int, std::vector<int> > > pidstepg4EEm;
std::map<int,std::map<int, std::vector<int> > > parentidstepg4EEm;
std::map<int,std::map<int, std::vector<int> > > enterstepg4EEm;
std::map<int,std::map<int, std::vector<int> > > leavestepg4EEm;


std::map<int,std::map<int, std::vector<float> > > postxstepg4EEm;
std::map<int,std::map<int, std::vector<float> > > postystepg4EEm;
std::map<int,std::map<int, std::vector<float> > > postzstepg4EEm;

std::map<int,std::map<int, std::vector<float> > > pretstepg4EEm;

std::map<int,std::map<int, std::vector<float> > > prexstepg4EEm;
std::map<int,std::map<int, std::vector<float> > > preystepg4EEm;
std::map<int,std::map<int, std::vector<float> > > prezstepg4EEm;


std::map<int,std::map<int, std::vector<float> > > preestepg4EEm;

std::map<int,std::map<int, std::vector<int> > > prechazstepg4EEm;
std::map<int,std::map<int, std::vector<float> > > premzstepg4EEm;



std::map<int,std::map<int, std::vector<int> > > idstepg4EEp;
std::map<int,std::map<int, std::vector<int> > > pidstepg4EEp;
std::map<int,std::map<int, std::vector<int> > > parentidstepg4EEp;

std::map<int,std::map<int, std::vector<float> > > postxstepg4EEp;
std::map<int,std::map<int, std::vector<float> > > postystepg4EEp;
std::map<int,std::map<int, std::vector<float> > > postzstepg4EEp;
std::map<int,std::map<int, std::vector<int> > > enterstepg4EEp;
std::map<int,std::map<int, std::vector<int> > > leavestepg4EEp;
std::map<int,std::map<int, std::vector<float> > > pretstepg4EEp;
std::map<int,std::map<int, std::vector<float> > > prexstepg4EEp;
std::map<int,std::map<int, std::vector<float> > > preystepg4EEp;
std::map<int,std::map<int, std::vector<float> > > prezstepg4EEp;

std::map<int,std::map<int, std::vector<float> > > preestepg4EEp;

std::map<int,std::map<int, std::vector<int> > > prechazstepg4EEp;
std::map<int,std::map<int, std::vector<float> > > premzstepg4EEp;




int ng4EB; 
int ietag4EB[61200];
int iphig4EB[61200];
float esumg4EB[61200];
float tming4EB[61200];
float xtming4EB[61200];
float ytming4EB[61200];
float ztming4EB[61200];

std::vector<std::vector<float> > *eg4EB; 
std::vector<std::vector<float> > *tg4EB; 
std::vector<std::vector<int> > *idg4EB; 
std::vector<std::vector<int> > *pidg4EB; 
std::vector<std::vector<int> > *parentidg4EB; 
std::vector<std::vector<int> > *enterg4EB; 
std::vector<std::vector<int> > *leaveg4EB; 
std::vector<std::vector<float> > *preeg4EB; 



std::vector<std::vector<float> > *postxg4EB; 
std::vector<std::vector<float> > *postyg4EB; 
std::vector<std::vector<float> > *postzg4EB; 

std::vector<std::vector<float> > *prexg4EB; 
std::vector<std::vector<float> > *preyg4EB; 
std::vector<std::vector<float> > *prezg4EB; 
std::vector<std::vector<float> > *pretg4EB; 
std::vector<std::vector<int> > *prechag4EB; 
std::vector<std::vector<float> > *premg4EB; 


int ng4EE;
int ixg4EE[14648];
int iyg4EE[14648];
int izg4EE[14648];
float esumg4EE[14648];
float tming4EE[14648];
float xtming4EE[14648];
float ytming4EE[14648];
float ztming4EE[14648];


std::vector<std::vector<float> > *eg4EE; 
std::vector<std::vector<float> > *tg4EE; 
std::vector<std::vector<int> > *idg4EE; 
std::vector<std::vector<int> > *pidg4EE; 
std::vector<std::vector<int> > *parentidg4EE; 

std::vector<std::vector<int> > *enterg4EE; 
std::vector<std::vector<int> > *leaveg4EE; 

std::vector<std::vector<float> > *preeg4EE; 


std::vector<std::vector<float> > *postxg4EE; 
std::vector<std::vector<float> > *postyg4EE; 
std::vector<std::vector<float> > *postzg4EE; 

std::vector<std::vector<float> > *prexg4EE; 
std::vector<std::vector<float> > *preyg4EE; 
std::vector<std::vector<float> > *prezg4EE; 
std::vector<std::vector<float> > *pretg4EE; 
std::vector<std::vector<int> > *prechag4EE; 
std::vector<std::vector<float> > *premg4EE; 


int debugCaloSD = 0; 

int treeConfig = 1; 
int saveOnlyEnterCur = 1; 
//int saveOnlyEnterCur = 0; 

//#define DebugLog

CaloSD::CaloSD(G4String name, const DDCompactView & cpv,
        SensitiveDetectorCatalog & clg, 
        edm::ParameterSet const & p, const SimTrackManager* manager,
        int tSlice, bool ignoreTkID) : 
  SensitiveCaloDetector(name, cpv, clg, p),
  G4VGFlashSensitiveDetector(), theTrack(0), preStepPoint(0), eminHit(0), 
  eminHitD(0), m_trackManager(manager), currentHit(0), runInit(false),
  timeSlice(tSlice), ignoreTrackID(ignoreTkID), hcID(-1), theHC(0), 
  meanResponse(0) {
  //Add Hcal Sentitive Detector Names
  
  collectionName.insert(name);

  csdname = name; 

  if(debugCaloSD) cout<<"CaloSD::CaloSD " << name <<" tslice "  << tSlice <<endl;

  if( "EcalHitsEB" == name){
    if(debugCaloSD) cout<<"creat newtree" <<endl;

    ncounter1 = 0; 
    ncounter2 = 0; 
    

    g4f = new TFile("g4simhitsEcal.root","RECREATE");
    g4t = new TTree("G4SIM","g4sim");
    
    eg4EB = new std::vector<std::vector<float> >; eg4EB->clear();
    tg4EB = new std::vector<std::vector<float> >; tg4EB->clear();
    idg4EB = new std::vector<std::vector<int> >; idg4EB->clear();
    pidg4EB = new std::vector<std::vector<int> >; pidg4EB->clear();
    parentidg4EB = new std::vector<std::vector<int> >; parentidg4EB->clear();
    enterg4EB = new std::vector<std::vector<int> >; enterg4EB->clear();
    leaveg4EB = new std::vector<std::vector<int> >; leaveg4EB->clear();

    preeg4EB = new std::vector<std::vector<float> >; preeg4EB->clear();
    prechag4EB = new std::vector<std::vector<int> >; prechag4EB->clear();
    premg4EB = new std::vector<std::vector<float> >; premg4EB->clear();
    
    

    postxg4EB = new std::vector<std::vector<float> >; postxg4EB->clear();
    postyg4EB = new std::vector<std::vector<float> >; postyg4EB->clear();
    postzg4EB = new std::vector<std::vector<float> >; postzg4EB->clear();    
    prexg4EB = new std::vector<std::vector<float> >; prexg4EB->clear();
    preyg4EB = new std::vector<std::vector<float> >; preyg4EB->clear();
    prezg4EB = new std::vector<std::vector<float> >; prezg4EB->clear();    
    pretg4EB = new std::vector<std::vector<float> >; pretg4EB->clear();    
    prexg4EE = new std::vector<std::vector<float> >; prexg4EE->clear();
    preyg4EE = new std::vector<std::vector<float> >; preyg4EE->clear();
    prezg4EE = new std::vector<std::vector<float> >; prezg4EE->clear();    
    pretg4EE = new std::vector<std::vector<float> >; pretg4EE->clear();    

    eg4EE = new std::vector<std::vector<float> >; eg4EE->clear();
    tg4EE = new std::vector<std::vector<float> >; tg4EE->clear();
    idg4EE = new std::vector<std::vector<int> >; idg4EE->clear();
    pidg4EE = new std::vector<std::vector<int> >; pidg4EE->clear();
    parentidg4EE = new std::vector<std::vector<int> >; parentidg4EE->clear();
    enterg4EE = new std::vector<std::vector<int> >; enterg4EE->clear();
    leaveg4EE = new std::vector<std::vector<int> >; leaveg4EE->clear();
    

    postxg4EE = new std::vector<std::vector<float> >; postxg4EE->clear();
    postyg4EE = new std::vector<std::vector<float> >; postyg4EE->clear();
    postzg4EE = new std::vector<std::vector<float> >; postzg4EE->clear();    

    preeg4EE = new std::vector<std::vector<float> >; preeg4EE->clear();


    
    g4t->Branch("ng4EB", &ng4EB, "ng4EB/I");
    g4t->Branch("ietag4EB", ietag4EB, "ietag4EB[ng4EB]/I");
    g4t->Branch("iphig4EB", iphig4EB, "iphig4EB[ng4EB]/I");
    g4t->Branch("esumg4EB", esumg4EB, "esumg4EB[ng4EB]/F");
    g4t->Branch("tming4EB", tming4EB, "tming4EB[ng4EB]/F");
    g4t->Branch("xtming4EB", xtming4EB, "xtming4EB[ng4EB]/F");
    g4t->Branch("ytming4EB", ytming4EB, "ytming4EB[ng4EB]/F");
    g4t->Branch("ztming4EB", ztming4EB, "ztming4EB[ng4EB]/F");
    
    
    g4t->Branch("ng4EE", &ng4EE, "ng4EE/I");
    g4t->Branch("ixg4EE", ixg4EE, "ixg4EE[ng4EE]/I");
    g4t->Branch("iyg4EE", iyg4EE, "iyg4EE[ng4EE]/I");
    g4t->Branch("izg4EE", izg4EE, "izg4EE[ng4EE]/I");
    g4t->Branch("esumg4EE", esumg4EE, "esumg4EE[ng4EE]/F");
    g4t->Branch("tming4EE", tming4EE, "tming4EE[ng4EE]/F");
    g4t->Branch("xtming4EE", xtming4EE, "xtming4EE[ng4EE]/F");
    g4t->Branch("ytming4EE", ytming4EE, "ytming4EE[ng4EE]/F");
    g4t->Branch("ztming4EE", ztming4EE, "ztming4EE[ng4EE]/F");
    
    
    if(treeConfig==1){
    g4t->Branch("eg4EB", "std::vector<std::vector<float> >", &eg4EB);
    g4t->Branch("tg4EB", "std::vector<std::vector<float> >", &tg4EB);
    g4t->Branch("idg4EB", "std::vector<std::vector<int> >", &idg4EB);
    g4t->Branch("pidg4EB", "std::vector<std::vector<int> >", &pidg4EB);
    g4t->Branch("parentidg4EB", "std::vector<std::vector<int> >", &parentidg4EB);
    g4t->Branch("enterg4EB", "std::vector<std::vector<int> >", &enterg4EB);
    g4t->Branch("leaveg4EB", "std::vector<std::vector<int> >", &leaveg4EB);
    
    g4t->Branch("postxg4EB", "std::vector<std::vector<float> >", &postxg4EB);
    g4t->Branch("postyg4EB", "std::vector<std::vector<float> >", &postyg4EB);
    g4t->Branch("postzg4EB", "std::vector<std::vector<float> >", &postzg4EB);
    
    g4t->Branch("prexg4EB", "std::vector<std::vector<float> >", &prexg4EB);
    g4t->Branch("preyg4EB", "std::vector<std::vector<float> >", &preyg4EB);
    g4t->Branch("prezg4EB", "std::vector<std::vector<float> >", &prezg4EB);
    g4t->Branch("pretg4EB", "std::vector<std::vector<float> >", &pretg4EB);

    g4t->Branch("preeg4EB", "std::vector<std::vector<float> >", &preeg4EB);


    
    g4t->Branch("eg4EE", "std::vector<std::vector<float> >", &eg4EE);
    g4t->Branch("tg4EE", "std::vector<std::vector<float> >", &tg4EE);
    g4t->Branch("idg4EE", "std::vector<std::vector<int> >", &idg4EE);
    g4t->Branch("pidg4EE", "std::vector<std::vector<int> >", &pidg4EE);
    g4t->Branch("parentidg4EE", "std::vector<std::vector<int> >", &parentidg4EE);

    g4t->Branch("enterg4EE", "std::vector<std::vector<int> >", &enterg4EE);
    g4t->Branch("leaveg4EE", "std::vector<std::vector<int> >", &leaveg4EE);
    

    g4t->Branch("postxg4EE", "std::vector<std::vector<float> >", &postxg4EE);
    g4t->Branch("postyg4EE", "std::vector<std::vector<float> >", &postyg4EE);
    g4t->Branch("postzg4EE", "std::vector<std::vector<float> >", &postzg4EE);
    g4t->Branch("prexg4EE", "std::vector<std::vector<float> >", &prexg4EE);
    g4t->Branch("preyg4EE", "std::vector<std::vector<float> >", &preyg4EE);
    g4t->Branch("prezg4EE", "std::vector<std::vector<float> >", &prezg4EE);
    g4t->Branch("pretg4EE", "std::vector<std::vector<float> >", &pretg4EE);
    
    g4t->Branch("preeg4EE", "std::vector<std::vector<float> >", &preeg4EE);


    
    }

    if(debugCaloSD) cout<<"tree branch set " <<endl;

    saveTree = true; 
  }






  //Parameters
  edm::ParameterSet m_CaloSD = p.getParameter<edm::ParameterSet>("CaloSD");
  energyCut    = m_CaloSD.getParameter<double>("EminTrack")*GeV;
  tmaxHit      = m_CaloSD.getParameter<double>("TmaxHit")*ns;
  std::vector<double> eminHits = m_CaloSD.getParameter<std::vector<double> >("EminHits");
  std::vector<double> tmaxHits = m_CaloSD.getParameter<std::vector<double> >("TmaxHits");
  std::vector<std::string> hcn = m_CaloSD.getParameter<std::vector<std::string> >("HCNames");
  std::vector<int>   useResMap = m_CaloSD.getParameter<std::vector<int> >("UseResponseTables");
  std::vector<double> eminHitX = m_CaloSD.getParameter<std::vector<double> >("EminHitsDepth");
  suppressHeavy= m_CaloSD.getParameter<bool>("SuppressHeavy");
  kmaxIon      = m_CaloSD.getParameter<double>("IonThreshold")*MeV;
  kmaxProton   = m_CaloSD.getParameter<double>("ProtonThreshold")*MeV;
  kmaxNeutron  = m_CaloSD.getParameter<double>("NeutronThreshold")*MeV;
  checkHits    = m_CaloSD.getUntrackedParameter<int>("CheckHits", 25);
  useMap       = m_CaloSD.getUntrackedParameter<bool>("UseMap", true);
  int verbn    = m_CaloSD.getUntrackedParameter<int>("Verbosity", 0);
  corrTOFBeam  = m_CaloSD.getParameter<bool>("CorrectTOFBeam");
  double beamZ = m_CaloSD.getParameter<double>("BeamPosition")*cm;
  correctT     = beamZ/c_light/nanosecond;

  SetVerboseLevel(verbn);
  for (unsigned int k=0; k<hcn.size(); ++k) {
    if (name == (G4String)(hcn[k])) {
      if (k < eminHits.size()) eminHit = eminHits[k]*MeV;
      if (k < eminHitX.size()) eminHitD= eminHitX[k]*MeV;
      if (k < tmaxHits.size()) tmaxHit = tmaxHits[k]*ns;
      if (k < useResMap.size() && useResMap[k] > 0) meanResponse = new CaloMeanResponse(p);
      break;
    }
  }
#ifdef DebugLog
  LogDebug("CaloSim") << "***************************************************" 
                      << "\n"
                      << "*                                                 *" 
                      << "\n"
                      << "* Constructing a CaloSD  with name " << GetName()
                      << "\n"
                      << "*                                                 *" 
                      << "\n"
                      << "***************************************************";
#endif
  slave      = new CaloSlaveSD(name);
  currentID  = CaloHitID(timeSlice, ignoreTrackID);
  previousID = CaloHitID(timeSlice, ignoreTrackID);
  
  primAncestor = 0;
  cleanIndex = 0;
  totalHits = 0;
  forceSave = false;

  //
  // Now attach the right detectors (LogicalVolumes) to me
  //
  std::vector<std::string> lvNames = clg.logicalNames(name);
  this->Register();
  for (std::vector<std::string>::iterator it=lvNames.begin(); it !=lvNames.end(); ++it) {
    this->AssignSD(*it);
#ifdef DebugLog
    LogDebug("CaloSim") << "CaloSD : Assigns SD to LV " << (*it);
#endif
  }

  edm::LogInfo("CaloSim") << "CaloSD: Minimum energy of track for saving it " 
                          << energyCut/GeV  << " GeV" << "\n"
                          << "        Use of HitID Map " << useMap << "\n"
                          << "        Check last " << checkHits 
                          << " before saving the hit\n" 
                          << "        Correct TOF globally by " << correctT
                          << " ns (Flag =" << corrTOFBeam << ")\n"
                          << "        Save hits recorded before " << tmaxHit
                          << " ns and if energy is above " << eminHit/MeV
                          << " MeV (for depth 0) or " << eminHitD/MeV
                          << " MeV (for nonzero depths); Time Slice Unit " 
                          << timeSlice << " Ignore TrackID Flag " << ignoreTrackID;
}

CaloSD::~CaloSD() { 

  if(debugCaloSD) cout<<"save g4t "<<endl;

  if(saveTree){
    g4f->cd();
    g4t->Write();
    g4f->Write() ;
    g4f->Close() ;
  }
  saveTree = false; 
  

  if (slave)           delete slave; 
  if (theHC)           delete theHC;
  if (meanResponse)    delete meanResponse;
}

bool CaloSD::ProcessHits(G4Step * aStep, G4TouchableHistory * ) {
  
  NaNTrap( aStep ) ;
  
  if (aStep == NULL) {
    return true;
  } else {
    if (getStepInfo(aStep)) {
      if (hitExists() == false && edepositEM+edepositHAD>0.) 
        currentHit = createNewHit();
    }
  }
  return true;
} 

bool CaloSD::ProcessHits(G4GFlashSpot* aSpot, G4TouchableHistory*) { 

  if (aSpot != NULL) {   
    theTrack = const_cast<G4Track *>(aSpot->GetOriginatorTrack()->GetPrimaryTrack());
    G4int particleCode = theTrack->GetDefinition()->GetPDGEncoding();
    
    if (particleCode == emPDG ||
        particleCode == epPDG ||
        particleCode == gammaPDG ) {
      edepositEM  = aSpot->GetEnergySpot()->GetEnergy();
      edepositHAD = 0.;
    } else {
      edepositEM  = 0.;
      edepositHAD = 0.;
    }
 
    if (edepositEM>0.) {
      G4Step *      fFakeStep          = new G4Step();
      preStepPoint                     = fFakeStep->GetPreStepPoint();
      G4StepPoint * fFakePostStepPoint = fFakeStep->GetPostStepPoint();
      preStepPoint->SetPosition(aSpot->GetPosition());
      fFakePostStepPoint->SetPosition(aSpot->GetPosition());
      
      G4TouchableHandle fTouchableHandle   = aSpot->GetTouchableHandle();
      preStepPoint->SetTouchableHandle(fTouchableHandle);
      fFakeStep->SetTotalEnergyDeposit(aSpot->GetEnergySpot()->GetEnergy());
      
      double       time   = 0;
      unsigned int unitID = setDetUnitId(fFakeStep);
      int          primaryID = getTrackID(theTrack);
      uint16_t     depth = getDepth(fFakeStep);

      if (unitID > 0) {
        currentID.setID(unitID, time, primaryID, depth);
#ifdef DebugLog
        LogDebug("CaloSim") << "CaloSD:: GetSpotInfo for"
                            << " Unit 0x" << std::hex << currentID.unitID() 
                            << std::dec << " Edeposit = " << edepositEM << " " 
                            << edepositHAD;
#endif
        // Update if in the same detector, time-slice and for same track   
        if (currentID == previousID) {
	  updateHit(currentHit);
        } else {
          posGlobal = aSpot->GetEnergySpot()->GetPosition();
          // Reset entry point for new primary
          if (currentID.trackID() != previousID.trackID()) {
            entrancePoint  = aSpot->GetPosition();
            entranceLocal  = aSpot->GetTouchableHandle()->GetHistory()->
                                      GetTopTransform().TransformPoint(entrancePoint);
            incidentEnergy = theTrack->GetKineticEnergy();
#ifdef DebugLog
            LogDebug("CaloSim") << "CaloSD: Incident energy " 
                                << incidentEnergy/GeV << " GeV and" 
                                << " entrance point " << entrancePoint 
                                << " (Global) " << entranceLocal << " (Local)";
#endif
          }

          if (checkHit() == false) currentHit = createNewHit();
        }
      }
      delete  fFakeStep;
    }
    return true;
  } 
  return false;
}                                   

double CaloSD::getEnergyDeposit(G4Step* aStep) {
  return aStep->GetTotalEnergyDeposit();
}

void CaloSD::Initialize(G4HCofThisEvent * HCE) { 
  totalHits = 0;
  
#ifdef DebugLog
  edm::LogInfo("CaloSim") << "CaloSD : Initialize called for " << GetName(); 
#endif
  
  //This initialization is performed at the beginning of an event
  //------------------------------------------------------------
  theHC = new CaloG4HitCollection(GetName(), collectionName[0]);
  
  if (hcID<0) hcID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  HCE->AddHitsCollection(hcID, theHC);
}

void CaloSD::EndOfEvent(G4HCofThisEvent* ) {
  // clean the hits for the last tracks
  

  if(debugCaloSD) cout<<"CaloSD::EndOfEvent " << csdname<<" "<<ncounter2 <<endl; 
  
  if(ncounter2%6==0){
    if(debugCaloSD) cout<<"filltree " <<estepg4EEm[85][58].size() <<endl; 
    
    ng4EB = 0; 
    for(int j=-85; j<= 85; j++){
      if(j==0) continue;
      for(int k = 1; k<= 360; k++){
	int nhit = estepg4EB[j][k].size();
	if(nhit<1) continue; 
	vector<float> temp; 
	vector<float> temp1; 
	vector<int> temp2; 
	vector<int> temp2a; 
	vector<int> temp2b; 

	vector<float> temp3; 
	vector<float> temp4; 
	vector<float> temp5; 
	vector<int> temp6;
	vector<int> temp7;

	vector<float> temp8;
        vector<float> temp9;
        vector<float> temp10;

	vector<float> temp11;
	
	vector<float> temp12;
        vector<float> temp13;
        vector<float> temp14;
	vector<float> temp15;
	

	
	double esum = 0; 
	double tmin = 1E9;
	int indtmin = 0;  
	for(int n=0; n<nhit; n++){

	  esum += estepg4EB[j][k][n];
	  
	  if(tmin>tstepg4EB[j][k][n]){
	    tmin = tstepg4EB[j][k][n]; 
	    indtmin = n; 
	  }
	  
	  if( saveOnlyEnterCur==1){
	    if(enterstepg4EB[j][k][n]==0) continue; 
	  }
	  
	  temp.push_back(estepg4EB[j][k][n]);
	  temp1.push_back(tstepg4EB[j][k][n]);
	  temp2.push_back(idstepg4EB[j][k][n]);
	  temp2a.push_back(pidstepg4EB[j][k][n]);
	  temp2b.push_back(parentidstepg4EB[j][k][n]);
	  temp3.push_back(postxstepg4EB[j][k][n]);
	  temp4.push_back(postystepg4EB[j][k][n]);
	  temp5.push_back(postzstepg4EB[j][k][n]);
	  
	  temp6.push_back(enterstepg4EB[j][k][n]);
	  temp7.push_back(leavestepg4EB[j][k][n]);
	  
	  temp8.push_back(prexstepg4EB[j][k][n]);
	  temp9.push_back(preystepg4EB[j][k][n]);
	  temp10.push_back(prezstepg4EB[j][k][n]);
	  temp11.push_back(pretstepg4EB[j][k][n]);

	  temp12.push_back(preestepg4EB[j][k][n]);

	  
	  
	}
	
	ietag4EB[ng4EB] = j; 
	iphig4EB[ng4EB] = k; 
	esumg4EB[ng4EB] = esum/1000; 
	tming4EB[ng4EB] = tmin;
	xtming4EB[ng4EB] = postxstepg4EB[j][k][indtmin];
	ytming4EB[ng4EB] = postystepg4EB[j][k][indtmin];
	ztming4EB[ng4EB] = postzstepg4EB[j][k][indtmin];
	

	eg4EB->push_back(temp);
	tg4EB->push_back(temp1);
	idg4EB->push_back(temp2);
	pidg4EB->push_back(temp2a);
	parentidg4EB->push_back(temp2b);
	postxg4EB->push_back(temp3);
	postyg4EB->push_back(temp4);
	postzg4EB->push_back(temp5);

	enterg4EB->push_back(temp6);
	leaveg4EB->push_back(temp7);

	prexg4EB->push_back(temp8);
	preyg4EB->push_back(temp9);
	prezg4EB->push_back(temp10);
	pretg4EB->push_back(temp11);

	preeg4EB->push_back(temp12);

	

	ng4EB++;
      }
    }

    ng4EE = 0; 
    for(int j=1;j<=100; j++){
      for(int k=1; k<=100; k++){
	int nhit = estepg4EEm[j][k].size();
	if(nhit<1) continue;
        vector<float> temp;
        vector<float> temp1;
	vector<int> temp2; 
	vector<int> temp2a; 
	vector<int> temp2b; 

	vector<float> temp3; 
	vector<float> temp4; 
	vector<float> temp5; 

	vector<int> temp6; 
	vector<int> temp7; 
	
	vector<float> temp8;
        vector<float> temp9;
        vector<float> temp10;

	vector<float> temp11;

	vector<float> temp12;
        vector<float> temp13;
        vector<float> temp14;
	vector<float> temp15;


	double esum = 0; 
	double tmin = 1E9; 
	int indtmin = 0; 
        for(int n=0; n<nhit; n++){

	  esum += estepg4EEm[j][k][n];
	  if( tmin > tstepg4EEm[j][k][n]){
	    tmin = tstepg4EEm[j][k][n];
	    indtmin = n; 
	  }

	  
	  if( saveOnlyEnterCur==1){
            if(enterstepg4EEm[j][k][n]==0) continue;
          }
	  
	  temp.push_back(estepg4EEm[j][k][n]);
          temp1.push_back(tstepg4EEm[j][k][n]);
	  temp2.push_back(idstepg4EEm[j][k][n]);
	  temp2a.push_back(pidstepg4EEm[j][k][n]);
	  temp2b.push_back(parentidstepg4EEm[j][k][n]);
	  temp3.push_back(postxstepg4EEm[j][k][n]);
	  temp4.push_back(postystepg4EEm[j][k][n]);
	  temp5.push_back(postzstepg4EEm[j][k][n]);

	  temp6.push_back(enterstepg4EEm[j][k][n]);
	  temp7.push_back(leavestepg4EEm[j][k][n]);

	  temp8.push_back(prexstepg4EEm[j][k][n]);
	  temp9.push_back(preystepg4EEm[j][k][n]);
	  temp10.push_back(prezstepg4EEm[j][k][n]);
	  temp11.push_back(pretstepg4EEm[j][k][n]);
	  
	  temp12.push_back(preestepg4EEm[j][k][n]);


	  
        }

	ixg4EE[ng4EE] = j;
	iyg4EE[ng4EE] = k;
	izg4EE[ng4EE] = -1;
	esumg4EE[ng4EE] = esum/1000; 
	tming4EE[ng4EE] = tmin;
	xtming4EE[ng4EE] = postxstepg4EEm[j][k][indtmin];
	ytming4EE[ng4EE] = postystepg4EEm[j][k][indtmin];
	ztming4EE[ng4EE] = postzstepg4EEm[j][k][indtmin];
	
	
	eg4EE->push_back(temp);
        tg4EE->push_back(temp1);
	idg4EE->push_back(temp2);
	pidg4EE->push_back(temp2a);
	parentidg4EE->push_back(temp2b);
	postxg4EE->push_back(temp3);
	postyg4EE->push_back(temp4);
	postzg4EE->push_back(temp5);

	enterg4EE->push_back(temp6);
	leaveg4EE->push_back(temp7);
	prexg4EE->push_back(temp8);
	preyg4EE->push_back(temp9);
	prezg4EE->push_back(temp10);
	pretg4EE->push_back(temp11);

	preeg4EE->push_back(temp12);
	

	
	ng4EE++;
      }
    }
    for(int j=1;j<=100; j++){
      for(int k=1; k<=100; k++){
	int nhit = estepg4EEp[j][k].size();
	if(nhit<1) continue;
        vector<float> temp;
        vector<float> temp1;
	vector<int> temp2; 
	vector<int> temp2a; 
	vector<int> temp2b; 

	vector<float> temp3; 
	vector<float> temp4; 
	vector<float> temp5; 

	vector<int> temp6;
        vector<int> temp7;

	vector<float> temp8;
        vector<float> temp9;
        vector<float> temp10;

	vector<float> temp11;

	vector<float> temp12;
        vector<float> temp13;
        vector<float> temp14;
	vector<float> temp15;


	double esum = 0; 
	double tmin = 1E9;
	int indtmin = 0; 
        for(int n=0; n<nhit; n++){
          
	  esum += estepg4EEp[j][k][n];
	  if(tmin > tstepg4EEp[j][k][n]){
	    tmin = tstepg4EEp[j][k][n];
	    indtmin = n; 
	  }
	  
	  if( saveOnlyEnterCur==1){
            if(enterstepg4EEp[j][k][n]==0) continue;
          }
	  
	  temp.push_back(estepg4EEp[j][k][n]);
          temp1.push_back(tstepg4EEp[j][k][n]);
	  temp2.push_back(idstepg4EEp[j][k][n]);
	  temp2a.push_back(pidstepg4EEp[j][k][n]);
	  temp2b.push_back(parentidstepg4EEp[j][k][n]);
	  temp3.push_back(postxstepg4EEp[j][k][n]);
	  temp4.push_back(postystepg4EEp[j][k][n]);
	  temp5.push_back(postzstepg4EEp[j][k][n]);
	  temp6.push_back(enterstepg4EEp[j][k][n]);
          temp7.push_back(leavestepg4EEp[j][k][n]);
	  temp8.push_back(prexstepg4EEp[j][k][n]);
	  temp9.push_back(preystepg4EEp[j][k][n]);
	  temp10.push_back(prezstepg4EEp[j][k][n]);
	  temp11.push_back(pretstepg4EEp[j][k][n]);
	  
	  
	  temp12.push_back(preestepg4EEp[j][k][n]);



        }
	ixg4EE[ng4EE] = j;
        iyg4EE[ng4EE] = k;
        izg4EE[ng4EE] = 1;
	esumg4EE[ng4EE] = esum/1000; 
	tming4EE[ng4EE] = tmin;
	xtming4EE[ng4EE] = postxstepg4EEp[j][k][indtmin];
	ytming4EE[ng4EE] = postystepg4EEp[j][k][indtmin];
	ztming4EE[ng4EE] = postzstepg4EEp[j][k][indtmin];
	
	eg4EE->push_back(temp);
        tg4EE->push_back(temp1);
	idg4EE->push_back(temp2);
	pidg4EE->push_back(temp2a);
	parentidg4EE->push_back(temp2b);
	postxg4EE->push_back(temp3);
	postyg4EE->push_back(temp4);
	postzg4EE->push_back(temp5);
	
	enterg4EE->push_back(temp6);
        leaveg4EE->push_back(temp7);
	prexg4EE->push_back(temp8);
	preyg4EE->push_back(temp9);
	prezg4EE->push_back(temp10);
	pretg4EE->push_back(temp11);
	
	preeg4EE->push_back(temp12);

	
	

        ng4EE++;
      }
    }

    if(debugCaloSD) cout<<"filling tree " << ng4EB <<" "<< eg4EB->size() <<" "<< ng4EE<<" "<<eg4EE->size()<<endl;
    g4t->Fill();
    
    
  }
  ncounter2 ++; 
  if(ncounter2 ==6){
    ncounter2 = 0; 
  }

  
  cleanHitCollection();
  
#ifdef DebugLog
  edm::LogInfo("CaloSim") << "CaloSD: EndofEvent entered with " << theHC->entries()
                          << " entries";
#endif
  //  TimeMe("CaloSD:sortAndMergeHits",false);
}

void CaloSD::clear() {} 

void CaloSD::DrawAll() {} 

void CaloSD::PrintAll() {
#ifdef DebugLog
  edm::LogInfo("CaloSim") << "CaloSD: Collection " << theHC->GetName();
#endif
  theHC->PrintAllHits();
} 

void CaloSD::fillHits(edm::PCaloHitContainer& c, std::string n) {
  if (slave->name() == n) c=slave->hits();
  slave->Clean();
}

bool CaloSD::getStepInfo(G4Step* aStep) {  

  preStepPoint = aStep->GetPreStepPoint(); 
  theTrack     = aStep->GetTrack();   
  
  G4int particleCode = theTrack->GetDefinition()->GetPDGEncoding();
  if (particleCode == emPDG ||
      particleCode == epPDG ||
      particleCode == gammaPDG ) {
    edepositEM  = getEnergyDeposit(aStep);
    edepositHAD = 0.;
  } else {
    edepositEM  = 0.;
    edepositHAD = getEnergyDeposit(aStep);
  }
  
  double       time  = (aStep->GetPostStepPoint()->GetGlobalTime())/nanosecond;
  unsigned int unitID= setDetUnitId(aStep);
  uint16_t     depth = getDepth(aStep);
  int          primaryID = getTrackID(theTrack);
  int parentID = theTrack->GetParentID();


  bool flag = (unitID > 0);

  double       timepre  = (aStep->GetPreStepPoint()->GetGlobalTime())/nanosecond;
  bool debug = true;
  
  if(flag && debug){

    G4ThreeVector pos1 = preStepPoint->GetPosition();
    G4ThreeVector pos2 = aStep->GetPostStepPoint()->GetPosition();

    const DetId detId(unitID);
    if(!detId.null() && detId.det() == DetId::Ecal){

      G4TouchableHandle touch1 =       preStepPoint->GetTouchableHandle();
      G4VPhysicalVolume* volume = touch1->GetVolume();
      G4String name = volume->GetName();
      int enterCurrent = 0; 
      G4String name2 = preStepPoint->GetPhysicalVolume()->GetLogicalVolume()->GetName(); 
      G4Material* material = preStepPoint  ->GetMaterial();
      G4String matname  = material->GetName();
      
      if (preStepPoint->GetStepStatus() == fGeomBoundary) {
	enterCurrent = 1; 
      }
      int leaveCurrent = 0; 
      if( aStep->GetPostStepPoint()->GetStepStatus() == fGeomBoundary){
	leaveCurrent = 1; 
      }
      
      
      G4VPhysicalVolume* mother = touch1->GetVolume(depth=1);
      G4String mothername = mother->GetName();
      
      if( EcalBarrel == detId.subdetId()){
	EBDetId eb(detId);
	int x = eb.ieta();
	int y = eb.iphi();
	estepg4EB[x][y].push_back(edepositEM);
        tstepg4EB[x][y].push_back(time);
	idstepg4EB[x][y].push_back(primaryID);
	pidstepg4EB[x][y].push_back(particleCode);
	parentidstepg4EB[x][y].push_back(parentID);
	enterstepg4EB[x][y].push_back(enterCurrent);
	leavestepg4EB[x][y].push_back(leaveCurrent);
	
	postxstepg4EB[x][y].push_back(pos2.x());
	postystepg4EB[x][y].push_back(pos2.y());
	postzstepg4EB[x][y].push_back(pos2.z());
	prexstepg4EB[x][y].push_back(pos1.x());
	preystepg4EB[x][y].push_back(pos1.y());
	prezstepg4EB[x][y].push_back(pos1.z());
	pretstepg4EB[x][y].push_back(timepre);
	
	preestepg4EB[x][y].push_back(aStep->GetPreStepPoint()->GetTotalEnergy());
	
	
	if(debugCaloSD && enterCurrent ){
	  if(x==8 && y==9){
	    cout<<"eb1xtal "<<enterCurrent<<" "<<aStep->GetPreStepPoint()->GetTotalEnergy()<<" "<< name <<" "<<name2<<" "
		<<matname<<" "<<mothername<<" ch/m "<< preStepPoint->GetCharge()<<" "<<preStepPoint->GetMass()<<" "<<particleCode <<" "<<parentID<<endl;
	  }
	  if(x==9 && y==9){
	    cout<<"eb2xtal "<<enterCurrent<<" "<<aStep->GetPreStepPoint()->GetTotalEnergy()<<" "<< name <<" "<<name2<<" "<<" "<<matname<<" "<<mothername<<endl;
	  }
	  
	}

	
      }else if(EcalEndcap == detId.subdetId() ) {
        EEDetId ee(detId);
	int x = ee.ix();
	int y = ee.iy();
	
	if( ee.zside() == -1 ){
	  estepg4EEm[x][y].push_back(edepositEM);
	  tstepg4EEm[x][y].push_back(time);
	  idstepg4EEm[x][y].push_back(primaryID);
	  pidstepg4EEm[x][y].push_back(particleCode);
	  enterstepg4EEm[x][y].push_back(enterCurrent);
	  leavestepg4EEm[x][y].push_back(leaveCurrent);
	  parentidstepg4EEm[x][y].push_back(parentID);
	  
	  postxstepg4EEm[x][y].push_back(pos2.x());
	  postystepg4EEm[x][y].push_back(pos2.y());
	  postzstepg4EEm[x][y].push_back(pos2.z());
	  prexstepg4EEm[x][y].push_back(pos1.x());
	  preystepg4EEm[x][y].push_back(pos1.y());
	  prezstepg4EEm[x][y].push_back(pos1.z());
	  pretstepg4EEm[x][y].push_back(timepre);
	  
	  preestepg4EEm[x][y].push_back(aStep->GetPreStepPoint()->GetTotalEnergy());

	
	  
	  if(debugCaloSD && ee.ix() == 85 && ee.iy() == 58 && ee.zside() == -1 ) {
	    std::cout<<"ee_step_xtal1: "<< edepositEM<<" "<<time<<" "<<primaryID<<" "<<pos1.x()<<" "<<pos1.y()<<" "<<pos1.z()<<" "<<pos2.x()<<" "<<pos2.y()<<" "<<pos2.z()<<" "<<name<<endl;
	    if(enterCurrent==1){ //theTrack->pos is the same as pos2 and theTrack->GetGlobaltime is the same as time
	      /// theTrack->GetTotalEnergy() == aStep->GetPostStepPoint()->GetTotalEnergy();
	      /// aStep->GetPreStepPoint()->GetTotalEnergy() is larger, the delta is ~ the deposted energy of this step
	      cout<<"xtal1_entercur "<< pos1.x()<<" "<<pos1.y()<<" "<<pos1.z()<<" "<< pos2.x()<<" "<<pos2.y()<<" "<<pos2.z()<<" "<<name<<" en "<<theTrack->GetTotalEnergy()<<" "<<theTrack->GetKineticEnergy()<<" trkt "<<theTrack->GetGlobalTime()<<" trkpos "<< theTrack->GetPosition().x()<<" "<<theTrack->GetPosition().y()<<" "<<theTrack->GetPosition().z()<<" etot "<< aStep->GetPostStepPoint()->GetTotalEnergy()<<" "<<aStep->GetPreStepPoint()->GetTotalEnergy()<<" "<<endl;
	    }
	    if(leaveCurrent==1){
	      cout<<"xtal1_leavecur "<< pos1.x()<<" "<<pos1.y()<<" "<<pos1.z()<<" "<< pos2.x()<<" "<<pos2.y()<<" "<<pos2.z()<<" "<<name<<endl;
	    }
	    
	  }


	  if(debugCaloSD && ee.ix() == 85 && ee.iy() == 59 && ee.zside() == -1 ) std::cout<<"ee_step_xtal2: "<< edepositEM<<" "<<time<<" "<<primaryID<<" "<<pos1.x()<<" "<<pos1.y()<<" "<<pos1.z()<<" "<<pos2.x()<<" "<<pos2.y()<<" "<<pos2.z()<<" "<<name<<endl;
	  
	}else{
	  estepg4EEp[x][y].push_back(edepositEM);
	  tstepg4EEp[x][y].push_back(time);
	  idstepg4EEp[x][y].push_back(primaryID);
	  pidstepg4EEp[x][y].push_back(particleCode);
	  enterstepg4EEp[x][y].push_back(enterCurrent);
	  leavestepg4EEp[x][y].push_back(leaveCurrent);
	  parentidstepg4EEp[x][y].push_back(parentID);

	  postxstepg4EEp[x][y].push_back(pos2.x());
	  postystepg4EEp[x][y].push_back(pos2.y());
	  postzstepg4EEp[x][y].push_back(pos2.z());
	  prexstepg4EEp[x][y].push_back(pos1.x());
	  preystepg4EEp[x][y].push_back(pos1.y());
	  prezstepg4EEp[x][y].push_back(pos1.z());
	  pretstepg4EEp[x][y].push_back(timepre);
	  preestepg4EEp[x][y].push_back(aStep->GetPreStepPoint()->GetTotalEnergy());
	

	}
	
	
      }
    }
  }
  


  if (flag) {
    currentID.setID(unitID, time, primaryID, depth);
#ifdef DebugLog
    G4TouchableHistory* touch =(G4TouchableHistory*)(theTrack->GetTouchable());
    LogDebug("CaloSim") << "CaloSD:: GetStepInfo for"
                        << " PV "     << touch->GetVolume(0)->GetName()
                        << " PVid = " << touch->GetReplicaNumber(0)
                        << " MVid = " << touch->GetReplicaNumber(1)
                        << " Unit   " << currentID.unitID() 
                        << " Edeposit = " << edepositEM << " " << edepositHAD;
  } else {
    G4TouchableHistory* touch =(G4TouchableHistory*)(theTrack->GetTouchable());
    LogDebug("CaloSim") << "CaloSD:: GetStepInfo for"
                        << " PV "     << touch->GetVolume(0)->GetName()
                        << " PVid = " << touch->GetReplicaNumber(0)
                        << " MVid = " << touch->GetReplicaNumber(1)
                        << " Unit   " << std::hex << unitID << std::dec 
                        << " Edeposit = " << edepositEM << " " << edepositHAD;
#endif
  }
  return flag;
}

G4ThreeVector CaloSD::setToLocal(G4ThreeVector global, const G4VTouchable* touch) {

  G4ThreeVector localPoint = touch->GetHistory()->GetTopTransform().TransformPoint(global);
  
  return localPoint;  
}

G4bool CaloSD::hitExists() {
#ifdef DebugLog
  if (currentID.trackID()<1)
    edm::LogWarning("CaloSim") << "***** CaloSD error: primaryID = " 
                               << currentID.trackID()
                               << " maybe detector name changed";
#endif  
  // Update if in the same detector, time-slice and for same track   
  if (currentID == previousID) {
    updateHit(currentHit);
    return true;
  }
  
  // Reset entry point for new primary
  posGlobal = preStepPoint->GetPosition();
  if (currentID.trackID() != previousID.trackID()) 
    resetForNewPrimary(preStepPoint->GetPosition(), preStepPoint->GetKineticEnergy());
  
  return checkHit();
}

G4bool CaloSD::checkHit() {  
  //look in the HitContainer whether a hit with the same ID already exists:
  bool       found = false;
  if (useMap) {
    std::map<CaloHitID,CaloG4Hit*>::const_iterator it = hitMap.find(currentID);
    if (it != hitMap.end()) {
      currentHit = it->second;
      found      = true;
    }
  } else {
    if (checkHits <= 0) return false;
    int  minhit= (theHC->entries()>checkHits ? theHC->entries()-checkHits : 0);
    int  maxhit= theHC->entries()-1;
    
    for (int j=maxhit; j>minhit&&!found; --j) {
      if ((*theHC)[j]->getID() == currentID) {
        currentHit = (*theHC)[j];
        found      = true;
      }
    }          
  }
  
  if (found) {
    updateHit(currentHit);
    return true;
  } else {
    return false;
  }
}

int CaloSD::getNumberOfHits() { return theHC->entries(); }

CaloG4Hit* CaloSD::createNewHit() {
#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD::CreateNewHit for"
                      << " Unit " << currentID.unitID() 
                      << " " << currentID.depth()
                      << " Edeposit = " << edepositEM << " " << edepositHAD;
  LogDebug("CaloSim") << " primary "    << currentID.trackID()
                      << " time slice " << currentID.timeSliceID()
                      << " For Track  " << theTrack->GetTrackID()
                      << " which is a " <<theTrack->GetDefinition()->GetParticleName()
                      << " of energy "  << theTrack->GetKineticEnergy()/GeV
                      << " " << theTrack->GetMomentum().mag()/GeV
                      << " daughter of part. " << theTrack->GetParentID()
                      << " and created by " ;
  
  if (theTrack->GetCreatorProcess()!=NULL)
    LogDebug("CaloSim") << theTrack->GetCreatorProcess()->GetProcessName() ;
  else 
    LogDebug("CaloSim") << "NO process";
#endif  
  
  CaloG4Hit* aHit;
  if (reusehit.size() > 0) {
    aHit = reusehit[0];
    aHit->setEM(0.);
    aHit->setHadr(0.);
    reusehit.erase(reusehit.begin());
  } else {
    aHit = new CaloG4Hit;
  }
  
  aHit->setID(currentID);
  aHit->setEntry(entrancePoint.x(),entrancePoint.y(),entrancePoint.z());
  aHit->setEntryLocal(entranceLocal.x(),entranceLocal.y(),entranceLocal.z());
  aHit->setPosition(posGlobal.x(),posGlobal.y(),posGlobal.z());
  aHit->setIncidentEnergy(incidentEnergy);
  updateHit(aHit);
  
  storeHit(aHit);
  double etrack = 0;
  if (currentID.trackID() == primIDSaved) { // The track is saved; nothing to be done
  } else if (currentID.trackID() == theTrack->GetTrackID()) {
    etrack= theTrack->GetKineticEnergy();
    //edm::LogInfo("CaloSim") << "CaloSD: set save the track " << currentID.trackID()
    //      << " etrack " << etrack << " eCut " << energyCut << " flag " << forceSave;
    if (etrack >= energyCut || forceSave) {
      TrackInformation* trkInfo = (TrackInformation *)(theTrack->GetUserInformation());
      trkInfo->storeTrack(true);
      trkInfo->putInHistory();
      //      trkInfo->setAncestor();
#ifdef DebugLog
      LogDebug("CaloSim") << "CaloSD: set save the track " << currentID.trackID()
                          << " with Hit";
#endif
    }
  } else {
    TrackWithHistory * trkh = tkMap[currentID.trackID()];
#ifdef DebugLog
    LogDebug("CaloSim") << "CaloSD : TrackwithHistory pointer for " 
                        << currentID.trackID() << " is " << trkh;
#endif
    if (trkh != NULL) {
      etrack = sqrt(trkh->momentum().Mag2());
      if (etrack >= energyCut) {
        trkh->save();
#ifdef DebugLog
        LogDebug("CaloSim") << "CaloSD: set save the track " 
                            << currentID.trackID() << " with Hit";
#endif
      }
    }
  }
  primIDSaved = currentID.trackID();
  if (useMap) totalHits++;
  return aHit;
}  

void CaloSD::updateHit(CaloG4Hit* aHit) {
  if (edepositEM+edepositHAD != 0) {
    aHit->addEnergyDeposit(edepositEM,edepositHAD);
#ifdef DebugLog
    LogDebug("CaloSim") << "CaloSD: Add energy deposit in " << currentID 
                        << " em " << edepositEM/MeV << " hadronic " 
                        << edepositHAD/MeV << " MeV"; 
#endif
  }

  // buffer for next steps:
  previousID = currentID;
}

void CaloSD::resetForNewPrimary(G4ThreeVector point, double energy) { 
  entrancePoint  = point;
  entranceLocal  = setToLocal(entrancePoint, preStepPoint->GetTouchable());
  incidentEnergy = energy;
#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: Incident energy " << incidentEnergy/GeV 
                      << " GeV and" << " entrance point " << entrancePoint 
                      << " (Global) " << entranceLocal << " (Local)";
#endif
}

double CaloSD::getAttenuation(G4Step* aStep, double birk1, double birk2, double birk3) {
  double weight = 1.;
  double charge = aStep->GetPreStepPoint()->GetCharge();

  if (charge != 0. && aStep->GetStepLength() > 0) {
    G4Material* mat = aStep->GetPreStepPoint()->GetMaterial();
    double density = mat->GetDensity();
    double dedx    = aStep->GetTotalEnergyDeposit()/aStep->GetStepLength();
    double rkb     = birk1/density;
    double c       = birk2*rkb*rkb;
    if (std::abs(charge) >= 2.) rkb /= birk3; // based on alpha particle data
    weight = 1./(1.+rkb*dedx+c*dedx*dedx);
#ifdef DebugLog
    LogDebug("CaloSim") << "CaloSD::getAttenuation in " << mat->GetName() 
                        << " Charge " << charge << " dE/dx " << dedx 
                        << " Birk Const " << rkb << ", " << c << " Weight = " 
                        << weight << " dE " << aStep->GetTotalEnergyDeposit();
#endif
  }
  return weight;
}

void CaloSD::update(const BeginOfRun *) {

  
  if(debugCaloSD) cout<<"CaloSD::update BeginOfRun" << csdname << endl; 


  G4ParticleTable * theParticleTable = G4ParticleTable::GetParticleTable();
  G4String particleName;
  emPDG = theParticleTable->FindParticle(particleName="e-")->GetPDGEncoding();
  epPDG = theParticleTable->FindParticle(particleName="e+")->GetPDGEncoding();
  gammaPDG = theParticleTable->FindParticle(particleName="gamma")->GetPDGEncoding();
#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: Particle code for e- = " << emPDG
                      << " for e+ = " << epPDG << " for gamma = " << gammaPDG;
#endif
  initRun();
  runInit = true;
} 

void CaloSD::update(const BeginOfEvent *) {

  if(debugCaloSD) cout<<"CaloSD::update BeginOfEvent " << csdname<< " "<<ncounter1 <<endl;

  if(ncounter1%6==0){
    if(debugCaloSD)cout<<"here_to_cleartree "<< ncounter1 <<endl;
    
    for(int j=-85; j<= 85; j++){
      if(j==0) continue;
      for(int k = 1; k<= 360; k++){
	estepg4EB[j][k].clear();
	tstepg4EB[j][k].clear();
	idstepg4EB[j][k].clear();
	pidstepg4EB[j][k].clear();
	parentidstepg4EB[j][k].clear();
	enterstepg4EB[j][k].clear();
	leavestepg4EB[j][k].clear();

	
	postxstepg4EB[j][k].clear();
	postystepg4EB[j][k].clear();
	postzstepg4EB[j][k].clear();
	
	prexstepg4EB[j][k].clear();
	preystepg4EB[j][k].clear();
	prezstepg4EB[j][k].clear();
	pretstepg4EB[j][k].clear();
	
	preestepg4EB[j][k].clear();

	
      }
    }
    
    for(int j=1; j<=100; j++){
      for(int k=1; k<=100; k++){
	estepg4EEm[j][k].clear();
	tstepg4EEm[j][k].clear();
	estepg4EEp[j][k].clear();
	tstepg4EEp[j][k].clear();

	idstepg4EEm[j][k].clear();
	pidstepg4EEm[j][k].clear();
	parentidstepg4EEm[j][k].clear();
	enterstepg4EEm[j][k].clear();
	leavestepg4EEm[j][k].clear();

	prexstepg4EEm[j][k].clear();
	preystepg4EEm[j][k].clear();
	prezstepg4EEm[j][k].clear();
	pretstepg4EEm[j][k].clear();

	preestepg4EEm[j][k].clear();

	
	
	idstepg4EEp[j][k].clear();
        pidstepg4EEp[j][k].clear();
	parentidstepg4EEp[j][k].clear();
	enterstepg4EEp[j][k].clear();
	leavestepg4EEp[j][k].clear();
	
	postxstepg4EEm[j][k].clear();
	postystepg4EEm[j][k].clear();
	postzstepg4EEm[j][k].clear();
	
	postxstepg4EEp[j][k].clear();
	postystepg4EEp[j][k].clear();
	postzstepg4EEp[j][k].clear();
	
	prexstepg4EEp[j][k].clear();
	preystepg4EEp[j][k].clear();
	prezstepg4EEp[j][k].clear();
	pretstepg4EEp[j][k].clear();

	preestepg4EEp[j][k].clear();

	
	
      }
    }

    eg4EB->clear();
    tg4EB->clear();
    idg4EB->clear();
    pidg4EB->clear();
    parentidg4EB->clear();
    enterg4EB->clear();
    leaveg4EB->clear();

    postxg4EB->clear();
    postyg4EB->clear();
    postzg4EB->clear();
    prexg4EB->clear();
    preyg4EB->clear();
    prezg4EB->clear();
    pretg4EB->clear();
    
    preeg4EB->clear();

    
    eg4EE->clear();
    tg4EE->clear();
    idg4EE->clear();
    pidg4EE->clear();
    parentidg4EE->clear();
    postxg4EE->clear();
    postyg4EE->clear();
    postzg4EE->clear();
    enterg4EE->clear();
    leaveg4EE->clear();
    prexg4EE->clear();
    preyg4EE->clear();
    prezg4EE->clear();
    pretg4EE->clear();
    
    preeg4EE->clear();
    
    
    
  }
  
  ncounter1 ++; 
  if(ncounter1==6){
    ncounter1 = 0;
  }
  

#ifdef DebugLog
  LogDebug("CaloSim")  << "CaloSD: Dispatched BeginOfEvent for " << GetName() 
                       << " !" ;
#endif
  clearHits();
}

void CaloSD::update(const EndOfTrack * trk) {
  int id = (*trk)()->GetTrackID();
  TrackInformation *trkI =(TrackInformation *)((*trk)()->GetUserInformation());
  int lastTrackID = -1;
  if (trkI) lastTrackID = trkI->getIDonCaloSurface();
  if (id == lastTrackID) {
    const TrackContainer * trksForThisEvent = m_trackManager->trackContainer();
    if (trksForThisEvent != NULL) {
      int it = (int)(trksForThisEvent->size()) - 1;
      if (it >= 0) {
        TrackWithHistory * trkH = (*trksForThisEvent)[it];
        if (trkH->trackID() == (unsigned int)(id)) tkMap[id] = trkH;
#ifdef DebugLog
        LogDebug("CaloSim") << "CaloSD: get track " << it << " from "
                            << "Container of size " << trksForThisEvent->size()
                            << " with ID " << trkH->trackID();
      } else {
        LogDebug("CaloSim") << "CaloSD: get track " << it << " from "
                            << "Container of size " << trksForThisEvent->size()
                            << " with no ID";
#endif
      }
    }
  }
}

void CaloSD::update(const ::EndOfEvent * ) {
  int count = 0, wrong = 0;
  bool ok;
  
  slave->ReserveMemory(theHC->entries());

  for (int i=0; i<theHC->entries(); ++i) {
    ok = saveHit((*theHC)[i]);
    ++count;
    if (!ok)  ++wrong;
  }
  
  edm::LogInfo("CaloSim") << "CaloSD: " << GetName() << " store " << count
                          << " hits recorded with " << wrong 
                          << " track IDs not given properly and "
                          << totalHits-count << " hits not passing cuts";
  summarize();

  tkMap.erase (tkMap.begin(), tkMap.end());
}

void CaloSD::clearHits() {  
  if (useMap) hitMap.erase (hitMap.begin(), hitMap.end());
  for (unsigned int i = 0; i<reusehit.size(); ++i) delete reusehit[i];
  std::vector<CaloG4Hit*>().swap(reusehit);
  cleanIndex  = 0;
  previousID.reset();
  primIDSaved = -99;
#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: Clears hit vector for " << GetName() << " " << slave;
#endif
  slave->Initialize();
#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: Initialises slave SD for " << GetName();
#endif
}

void CaloSD::initRun() {}

int CaloSD::getTrackID(G4Track* aTrack) {
  int primaryID = 0;
  forceSave = false;
  TrackInformation* trkInfo=(TrackInformation *)(aTrack->GetUserInformation());
  if (trkInfo) {
    primaryID = trkInfo->getIDonCaloSurface(); 
#ifdef DebugLog
    LogDebug("CaloSim") << "CaloSD: hit update from track Id on Calo Surface " 
                        << trkInfo->getIDonCaloSurface();
#endif   
  } else {
    primaryID = aTrack->GetTrackID();
#ifdef DebugLog
    edm::LogWarning("CaloSim") << "CaloSD: Problem with primaryID **** set by "
                               << "force to TkID **** " << primaryID << " in "
                               << preStepPoint->GetTouchable()->GetVolume(0)->GetName();
#endif
  }
  return primaryID;
}

uint16_t CaloSD::getDepth(G4Step*) { return 0; }

bool CaloSD::filterHit(CaloG4Hit* hit, double time) {
  double emin(eminHit);
  if (hit->getDepth() > 0) emin = eminHitD;
#ifdef DebugLog
  LogDebug("CaloSim") << "Depth " << hit->getDepth() << " Emin = " << emin << " ("
                      << eminHit << ", " << eminHitD << ")";
#endif   
  return ((time <= tmaxHit) && (hit->getEnergyDeposit() > emin));
}

double CaloSD::getResponseWt(G4Track* aTrack) {
  if (meanResponse) {
    TrackInformation * trkInfo = (TrackInformation *)(aTrack->GetUserInformation());
    return meanResponse->getWeight(trkInfo->genParticlePID(), trkInfo->genParticleP());
  } else {
    return 1;
  }
}

void CaloSD::storeHit(CaloG4Hit* hit) {
  if (previousID.trackID()<0) return;
  if (hit == 0) {
    edm::LogWarning("CaloSim") << "CaloSD: hit to be stored is NULL !!";
    return;
  }
  
  theHC->insert(hit);
  if (useMap) hitMap.insert(std::pair<CaloHitID,CaloG4Hit*>(previousID,hit));
}

bool CaloSD::saveHit(CaloG4Hit* aHit) {  
  int tkID;
  bool ok   = true;
  if (m_trackManager) {
    tkID = m_trackManager->giveMotherNeeded(aHit->getTrackID());
    if (tkID == 0) {
      if (m_trackManager->trackExists(aHit->getTrackID())) tkID = (aHit->getTrackID());
      else {
	ok = false;
      }
    }
  } else {
    tkID = aHit->getTrackID();
    ok = false;
  }
  //  edm::LogInfo("CaloSim") << "CalosD: Track ID " << aHit->getTrackID() << " changed to " << tkID << " by SimTrackManager" << " Status " << ok;
#ifdef DebugLog
  LogDebug("CaloSim") << "CalosD: Track ID " << aHit->getTrackID() 
                      << " changed to " << tkID << " by SimTrackManager"
                      << " Status " << ok;
#endif
  double time = aHit->getTimeSlice();
  if (corrTOFBeam) time += correctT;
  slave->processHits(aHit->getUnitID(), aHit->getEM()/GeV, 
                     aHit->getHadr()/GeV, time, tkID, aHit->getDepth());
#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: Store Hit at " << std::hex 
                      << aHit->getUnitID() << std::dec << " " 
                      << aHit->getDepth() << " due to " << tkID 
                      << " in time " << time << " of energy " 
                      << aHit->getEM()/GeV << " GeV (EM) and " 
                      << aHit->getHadr()/GeV << " GeV (Hadr)";
#endif
  return ok;
}

void CaloSD::summarize() {}

void CaloSD::update(const BeginOfTrack * trk) {
  int primary = -1;
  TrackInformation * trkInfo = (TrackInformation *)((*trk)()->GetUserInformation());
  if ( trkInfo->isPrimary() ) primary = (*trk)()->GetTrackID();
  
#ifdef DebugLog
  LogDebug("CaloSim") << "New track: isPrimary " << trkInfo->isPrimary() 
                      << " primary ID = " << primary 
                      << " primary ancestor ID " << primAncestor;
#endif
  
  // update the information if a different primary track ID 
  
  if (primary > 0 && primary != primAncestor) {
    primAncestor = primary;
    
    // clean the hits information
    
    if (theHC->entries()>0) cleanHitCollection();
    
  }
}

void CaloSD::cleanHitCollection() {
  std::vector<CaloG4Hit*>* theCollection = theHC->GetVector();

#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: collection before merging, size = " << theHC->entries();
#endif
  
  selIndex.reserve(theHC->entries()-cleanIndex);
  if ( reusehit.size() == 0 ) reusehit.reserve(theHC->entries()-cleanIndex); 

  // if no map used, merge before hits to have the save situation as a map
  if ( !useMap ) {
    hitvec.swap(*theCollection);
    sort((hitvec.begin()+cleanIndex), hitvec.end(), CaloG4HitLess());
#ifdef DebugLog
    LogDebug("CaloSim") << "CaloSD::cleanHitCollection: sort hits in buffer "
                        << "starting from element = " << cleanIndex;
    for (unsigned int i = 0; i<hitvec.size(); ++i) 
      LogDebug("CaloSim")<<i<<" "<<*hitvec[i];
#endif
    unsigned int i, j;
    CaloG4HitEqual equal;
    for (i=cleanIndex; i<hitvec.size(); ++i) {
      selIndex.push_back(i-cleanIndex);
      int jump = 0;
      for (j = i+1; j <hitvec.size() && equal(hitvec[i], hitvec[j]); ++j) {
        ++jump;
        // merge j to i
        (*hitvec[i]).addEnergyDeposit(*hitvec[j]);
        (*hitvec[j]).setEM(0.);
        (*hitvec[j]).setHadr(0.);
        reusehit.push_back(hitvec[j]);
      }
      i+=jump;
    }
#ifdef DebugLog
    LogDebug("CaloSim") << "CaloSD: cleanHitCollection merge the hits in buffer ";
    for (unsigned int i = 0; i<hitvec.size(); ++i) 
      LogDebug("CaloSim")<<i<<" "<<*hitvec[i];
#endif
    for ( unsigned int i = cleanIndex; i < cleanIndex+selIndex.size(); ++i ) {
      hitvec[i] = hitvec[selIndex[i-cleanIndex]+cleanIndex];
    }
    hitvec.resize(cleanIndex+selIndex.size());
#ifdef DebugLog
    LogDebug("CaloSim") << "CaloSD::cleanHitCollection: remove the merged hits in buffer,"
                        << " new size = " << hitvec.size();
    for (unsigned int i = 0; i<hitvec.size(); ++i) 
      LogDebug("CaloSim")<<i<<" "<<*hitvec[i];
#endif
    hitvec.swap(*theCollection);
    std::vector<CaloG4Hit*>().swap(hitvec);
    selIndex.clear();
    totalHits = theHC->entries();
  }

#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: collection after merging, size = " << theHC->entries();
#endif

  int addhit = 0;

#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: Size of reusehit after merge = " << reusehit.size();
  LogDebug("CaloSim") << "CaloSD: Starting hit selection from index = " << cleanIndex;
#endif
  
  selIndex.reserve(theCollection->size()-cleanIndex);
  for (unsigned int i = cleanIndex; i<theCollection->size(); ++i) {   
    CaloG4Hit* aHit((*theCollection)[i]);
    
    // selection
    
    double time = aHit->getTimeSlice();
    if (corrTOFBeam) time += correctT;
    if (!filterHit(aHit,time)) {
#ifdef DebugLog
      LogDebug("CaloSim") << "CaloSD: dropped CaloG4Hit " << " " << *aHit; 
#endif
      
      // create the list of hits to be reused
      
      reusehit.push_back((*theCollection)[i]);
      ++addhit;
    } else {
      selIndex.push_back(i-cleanIndex);
    }
  }

#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: Size of reusehit after selection = " << reusehit.size()  
                      << " Number of added hit = " << addhit;
#endif
  if (useMap) {
    if ( addhit>0 ) {
      int offset = reusehit.size()-addhit;
      for (int ii = addhit-1; ii>=0; --ii) {
        CaloHitID theID = reusehit[offset+ii]->getID();
        hitMap.erase(theID);
      }
    }
  }
  for (unsigned int j = 0; j<selIndex.size(); ++j) {
    (*theCollection)[cleanIndex+j] = (*theCollection)[cleanIndex+selIndex[j]];
  }

  theCollection->resize(cleanIndex+selIndex.size());
  std::vector<unsigned int>().swap(selIndex);

#ifdef DebugLog
  LogDebug("CaloSim") << "CaloSD: hit collection after selection, size = "
                      << theHC->entries();
  theHC->PrintAllHits();
#endif
    
  cleanIndex = theHC->entries();
}
