// -*- C++ -*-
//
// Package:    RecoAnalyzer
// Class:      RecoAnalyzer
// 
/**\class RecoAnalyzer RecoAnalyzer.cc RECOAnalyzer/RecoAnalyzer/src/RecoAnalyzer.cc

 Description: [one line class summary]

 Implementation:
     [Notes on implementation]
*/
//
// Original Author:  Yong Yang,32 3-A03,+41227679643,
//         Created:  Fri May 13 16:42:13 CEST 2011
// $Id: RecoAnalyzer.cc,v 1.5 2012/02/18 22:04:16 yangyong Exp $
//
//
 
#include "RECOAnalyzer/RecoAnalyzer/interface/RecoAnalyzer.h"



#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectron.h"
#include "DataFormats/EgammaReco/interface/ElectronSeed.h"
#include "DataFormats/GsfTrackReco/interface/GsfTrack.h"
#include "DataFormats/TrackCandidate/interface/TrackCandidateCollection.h"
#include "DataFormats/TrackingRecHit/interface/TrackingRecHitFwd.h"
#include "DataFormats/EgammaCandidates/interface/GsfElectronFwd.h"
#include "DataFormats/BeamSpot/interface/BeamSpot.h"
#include "DataFormats/PatCandidates/interface/Electron.h"
#include "DataFormats/VertexReco/interface/Vertex.h"
#include "DataFormats/VertexReco/interface/VertexFwd.h"

#include "SimDataFormats/CaloHit/interface/PCaloHit.h"
#include "SimDataFormats/CrossingFrame/interface/MixCollection.h"


using namespace std;
using namespace edm;
using namespace reco;



class PtSorter {
public:
  template <class T> bool operator() ( const T& a, const T& b ) {
    return ( a.pt() > b.pt() );
  }
};

class EtSorter {
public:
  template <class T> bool operator() ( const T& a, const T& b ) {
    return ( a.et() > b.et() );
  }
};



//
// constants, enums and typedefs
//

//
// static data member definitions
//

//
// constructors and destructor
//
RecoAnalyzer::RecoAnalyzer(const edm::ParameterSet& iConfig)

{
  
  cout<<"RecoAnalyzer initialzied.. " <<endl;
  
  
  m_muonSrc    = iConfig.getUntrackedParameter<edm::InputTag>("muons",edm::InputTag("muons"));
  
  m_genSrc = iConfig.getUntrackedParameter<edm::InputTag>("genParticles",edm::InputTag("genParticles"));
  m_photonSrc = iConfig.getUntrackedParameter<edm::InputTag>("photon",edm::InputTag("photons"));
  m_electronSrc = iConfig.getUntrackedParameter<edm::InputTag>("electron",edm::InputTag("gsfElectrons"));
  m_PileupSrc = iConfig.getUntrackedParameter<edm::InputTag>("pileup",edm::InputTag("addPileupInfo"));
  
  
  m_hitsProducerTag  = iConfig.getUntrackedParameter<std::string>("hitsProducer","g4SimHits");
  
  
  debug_ = iConfig.getParameter<int> ("debugLevel");
  
  
  outputFile_   = iConfig.getParameter<std::string>("outputFile");
  rootFile_ = new TFile(outputFile_.c_str(),"RECREATE"); // open output file to store root-trees. 

  

  nEventsProcessed  =0;
  nWarning = 0; 
  
}


RecoAnalyzer::~RecoAnalyzer()
{
 
  // do anything here that needs to be done at desctruction time
   // (e.g. close files, deallocate resources etc.)

}


//
// member functions
//

// ------------ method called for each event  ------------
void
RecoAnalyzer::analyze(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
    
  
  runNumber      = iEvent.id().run();
  evtNumber    = iEvent.id().event();
  lumiBlock = iEvent.luminosityBlock();
  bunchX          = iEvent.bunchCrossing();
  orbitNumber       = iEvent.orbitNumber();
  const edm::Timestamp jtime = iEvent.time();
  evtTime = jtime.value() >> 32;
  isRealData = iEvent.isRealData();
  

  if( !isRealData && nEventsProcessed ==0  ){
    
    
    cout<<"additional branch for MC matchinng defined " <<endl; 
  
  ///gen electron, muon,photon
  Analysis->Branch("nGenPht",&nGenPht,"nGenPht/I");
  Analysis->Branch("etaGenPht",etaGenPht,"etaGenPht[nGenPht]/F");
  Analysis->Branch("phiGenPht",phiGenPht,"phiGenPht[nGenPht]/F");
  Analysis->Branch("ptGenPht",ptGenPht,"ptGenPht[nGenPht]/F");
  Analysis->Branch("vxGenPht",vxGenPht,"vxGenPht[nGenPht]/F");
  Analysis->Branch("vyGenPht",vyGenPht,"vyGenPht[nGenPht]/F");
  Analysis->Branch("vzGenPht",vzGenPht,"vzGenPht[nGenPht]/F");
  Analysis->Branch("pidmomGenPht",pidmomGenPht,"pidmomGenPht[nGenPht]/I");
  Analysis->Branch("pidmom2GenPht",pidmom2GenPht,"pidmom2GenPht[nGenPht]/I");
  Analysis->Branch("indmom2GenPht",indmom2GenPht,"indmom2GenPht[nGenPht]/I");
  Analysis->Branch("pidmom3GenPht",pidmom3GenPht,"pidmom3GenPht[nGenPht]/I");
  Analysis->Branch("statusGenPht",statusGenPht,"statusGenPht[nGenPht]/I");
  Analysis->Branch("convGenPht",convGenPht,"convGenPht[nGenPht]/I");
  Analysis->Branch("vtxconvGenPht",vtxconvGenPht,"vtxconvGenPht[nGenPht][3]/F");
  
  Analysis->Branch("nGenEle",&nGenEle,"nGenEle/I");
  Analysis->Branch("etaGenEle",etaGenEle,"etaGenEle[nGenEle]/F");
  Analysis->Branch("phiGenEle",phiGenEle,"phiGenEle[nGenEle]/F");
  Analysis->Branch("ptGenEle",ptGenEle,"ptGenEle[nGenEle]/F");
  Analysis->Branch("pidmomGenEle",pidmomGenEle,"pidmomGenEle[nGenEle]/I");
  Analysis->Branch("pidmom2GenEle",pidmom2GenEle,"pidmom2GenEle[nGenEle]/I");
  Analysis->Branch("pidmom3GenEle",pidmom3GenEle,"pidmom3GenEle[nGenEle]/I");

  Analysis->Branch("barcodemomGenEle",barcodemomGenEle,"barcodemomGenEle[nGenEle]/I");

  Analysis->Branch("statusGenEle",statusGenEle,"statusGenEle[nGenEle]/I");
  Analysis->Branch("vxGenEle",vxGenEle,"vxGenEle[nGenEle]/F");
  Analysis->Branch("vyGenEle",vyGenEle,"vyGenEle[nGenEle]/F");
  Analysis->Branch("vzGenEle",vzGenEle,"vzGenEle[nGenEle]/F");
  Analysis->Branch("chaGenEle",chaGenEle,"chaGenEle[nGenEle]/I");
  

  Analysis->Branch("nGenMu",&nGenMu,"nGenMu/I");
  Analysis->Branch("etaGenMu",etaGenMu,"etaGenMu[nGenMu]/F");
  Analysis->Branch("phiGenMu",phiGenMu,"phiGenMu[nGenMu]/F");
  Analysis->Branch("ptGenMu",ptGenMu,"ptGenMu[nGenMu]/F");
  Analysis->Branch("chaGenMu",chaGenMu,"chaGenMu[nGenMu]/I");
  Analysis->Branch("pidmomGenMu",pidmomGenMu,"pidmomGenMu[nGenMu]/I");
  Analysis->Branch("pidmom2GenMu",pidmom2GenMu,"pidmom2GenMu[nGenMu]/I");
  Analysis->Branch("pidmom3GenMu",pidmom3GenMu,"pidmom3GenMu[nGenMu]/I");
  Analysis->Branch("statusGenMu",statusGenMu,"statusGenMu[nGenMu]/I");
  Analysis->Branch("vxGenMu",vxGenMu,"vxGenMu[nGenMu]/F");
  Analysis->Branch("vyGenMu",vyGenMu,"vyGenMu[nGenMu]/F");
  Analysis->Branch("vzGenMu",vzGenMu,"vzGenMu[nGenMu]/F");
  Analysis->Branch("barcodemomGenMu",barcodemomGenMu,"barcodemomGenMu[nGenMu]/I");
    
  Analysis->Branch("genhiggsm",&genhiggsm,"genhiggsm/F");
  Analysis->Branch("genhiggspt",&genhiggspt,"genhiggspt/F");
  Analysis->Branch("genhiggseta",&genhiggseta,"genhiggseta/F");
  Analysis->Branch("genhiggsphi",&genhiggsphi,"genhiggsphi/F");
  Analysis->Branch("genhiggsvx",&genhiggsvx,"genhiggsvx/F");
  Analysis->Branch("genhiggsvy",&genhiggsvy,"genhiggsvy/F");
  Analysis->Branch("genhiggsvz",&genhiggsvz,"genhiggsvz/F");
  
    ///Analysis->Branch("genhiggsstatus",&genhiggsstatus,"genhiggsstatus/I");
    
  }
  
  
  //clear
  for(int j=-85; j<=85; j++){
    if(j==0) continue; 
    for(int k=1; k<=360; k++){
      eEBsim[j][k].clear();
      tEBsim[j][k].clear();
      bEBsim[j][k].clear();
    }
  }
  for(int j=1;j<=100; j++){
    for(int k=1; k<=100; k++){
      eEEmsim[j][k].clear();
      tEEmsim[j][k].clear();
      bEEmsim[j][k].clear();
      eEEpsim[j][k].clear();
      tEEpsim[j][k].clear();
      bEEpsim[j][k].clear();
    }
  }
  

  
   // Get input
   edm::Handle<CrossingFrame<PCaloHit> > crossingFrame;

   // test access to SimHits
   const std::string barrelHitsName    ( m_hitsProducerTag + "EcalHitsEB" ) ;
   const std::string endcapHitsName    ( m_hitsProducerTag + "EcalHitsEE" ) ;
   const std::string preshowerHitsName ( m_hitsProducerTag + "EcalHitsES" ) ;

   iEvent.getByLabel( "mix",
                     barrelHitsName ,
                     crossingFrame    ) ;

   MixCollection<PCaloHit>* EBHits (
      !crossingFrame.isValid() ? 0 :
      new MixCollection<PCaloHit>( crossingFrame.product() ) ) ;
   
   const bool isEB ( crossingFrame.isValid() &&
                     0 != EBHits             &&
                     EBHits->inRegistry()       ) ;
   
   
   iEvent.getByLabel( "mix",
		      endcapHitsName ,
                     crossingFrame    ) ;
   
   MixCollection<PCaloHit>* EEHits (
				    !crossingFrame.isValid() ? 0 :
				    new MixCollection<PCaloHit>( crossingFrame.product() ) ) ;
   
   const bool isEE ( crossingFrame.isValid() &&
                     0 != EEHits             &&
                     EEHits->inRegistry()       ) ;
   
   
   
   if( isEB ) {
     std::auto_ptr<MixCollection<PCaloHit> >  hits( EBHits ) ;
     
     for( MixCollection<PCaloHit>::MixItr hitItr ( hits->begin() ) ;
	  hitItr != hits->end() ; ++hitItr )       {
       const PCaloHit& hit ( *hitItr ) ;
       const int bunch ( hitItr.bunch() ) ;
       
       if(  std::isnan(hit.time())) continue; 
       const DetId detId ( hit.id() );
       if( ! detId.null() ){
	 EBDetId eb(detId);
	 int x = eb.ieta();
	 int y = eb.iphi();
	 eEBsim[x][y].push_back(hit.energy());
	 tEBsim[x][y].push_back(hit.time());
	 bEBsim[x][y].push_back(bunch);
       }else{
	 cout<<"warning null detIDeb?? "<<endl; 
       }
     }
   }
   

   if( isEE ) {
     std::auto_ptr<MixCollection<PCaloHit> >  hits( EEHits ) ;

     for( MixCollection<PCaloHit>::MixItr hitItr ( hits->begin() ) ;
	  hitItr != hits->end() ; ++hitItr )       {
       const PCaloHit& hit ( *hitItr ) ;
       const int bunch ( hitItr.bunch() ) ;
       
       if(  std::isnan(hit.time())) continue; 
       const DetId detId ( hit.id() );
       if( ! detId.null() ){
	 EEDetId eb(detId);
	 int x = eb.ix();
	 int y = eb.iy();
	 if(eb.zside()==-1){
	   eEEmsim[x][y].push_back(hit.energy());
	   tEEmsim[x][y].push_back(hit.time());
	   bEEmsim[x][y].push_back(bunch);
	 }else{
	   eEEpsim[x][y].push_back(hit.energy());
	   tEEpsim[x][y].push_back(hit.time());
	   bEEpsim[x][y].push_back(bunch);
	 }
	 
       }else{
	 cout<<"warning null detIDee?? "<<endl; 
       }
     }
   }
   


   ///geometry
   
   //bool saveGem = true; 
   bool saveGem = false ;
   
   if( saveGem){
   cout<<"getting geom " <<endl; 
   edm::ESHandle<CaloGeometry> pG;
   iSetup.get<CaloGeometryRecord>().get(pG); 
   const CaloSubdetectorGeometry* geom_eb=pG->getSubdetectorGeometry(DetId::Ecal,EcalBarrel);//EcalBarrel = 1
   const CaloSubdetectorGeometry* geom_ee=pG->getSubdetectorGeometry(DetId::Ecal,EcalEndcap);
   const CaloSubdetectorGeometry* geom_es=pG->getSubdetectorGeometry(DetId::Ecal,EcalPreshower);
   
   cout<<"got geom " <<endl; 


   //preshower
   for(int iz=-1; iz<=1; iz++){
     if(iz==0) continue; 
     int indz = iz<0? 0: 1; 
     for(int plane =1; plane<=2; plane++){
       for(int strip=1; strip<=32; strip++){
	 for(int ix=1; ix<=40; ix++){
	   for(int iy=1; iy<=40; iy++){

	     xESAll[indz][plane-1][strip-1][ix-1][iy-1] = -999;
	     yESAll[indz][plane-1][strip-1][ix-1][iy-1] = -999;
	     zESAll[indz][plane-1][strip-1][ix-1][iy-1] = -999;
	     
	     //ESDetId::validDetId(int istrip, int ixs, int iys, int iplane, int iz) 
	     bool valid = ESDetId::validDetId(strip, ix, iy, plane, iz);
	     if(!valid) continue; 
	     ESDetId det =  ESDetId(strip, ix, iy, plane, iz);
	     const CaloCellGeometry* cell = geom_es->getGeometry(det);
	     GlobalPoint posi = cell->getPosition();
	     //iz, plane, strip, ix, iy
	     //double xESAll[2][2][32][40][40];
	     xESAll[indz][plane-1][strip-1][ix-1][iy-1] = posi.x();
	     yESAll[indz][plane-1][strip-1][ix-1][iy-1] = posi.y();
	     zESAll[indz][plane-1][strip-1][ix-1][iy-1] = posi.z();
	     
	   }
	 }
       }
     }
   }
   
   


      
   for(int x=-85; x<=85; x++){
     if(x==0) continue;
     for(int y=1; y<=360; y++){
       EBDetId det = EBDetId(x,y,EBDetId::ETAPHIMODE);
       const CaloCellGeometry* cell=geom_eb->getGeometry(det);
       GlobalPoint posi = cell->getPosition();
       int ix = x+85; 
       if(x>=1) ix= ix-1; 
       xEBAll[ix][y-1] = posi.x();
       yEBAll[ix][y-1] = posi.y();
       zEBAll[ix][y-1] = posi.z();
       etaEBAll[ix][y-1] = posi.eta();
       phiEBAll[ix][y-1] = posi.phi();
       const CaloCellGeometry::CornersVec& co ( cell->getCorners() ) ;
       int nc = 0; 
       for( CaloCellGeometry::CornersVec::const_iterator point = co.begin();
	    point != co.end() && nc < 8 ; ++point){
	 // if(j==1 && k>=1 && k<=10) cout<<"corner j/k" << j<<"/"<<k<<" pos "<< point->x()<<" "<<point->y()<<" "<<point->z()<<endl; 
	 coxEBAll[ix][y-1][nc] = point->x();
	 coyEBAll[ix][y-1][nc] = point->y();
	 cozEBAll[ix][y-1][nc] = point->z();
	 nc ++; 
       }
     }
   }
   
   
   for( int j=0; j<2; j++){
     int iz = -1;
     if( j==1) iz = 1;
     for( int ix = 0; ix <=100; ix++){
       for( int iy =0; iy <=100; iy++){
	 xEEAll[j][ix][iy]  =0;
	 yEEAll[j][ix][iy] = 0;
	 zEEAll[j][ix][iy] = 0;
	 etaEEAll[j][ix][iy] = 0;
	 phiEEAll[j][ix][iy] = 0;
	 
	 for(int n=0; n<8; n++){
	   coxEEAll[j][ix][iy][n] = 0;
	   coyEEAll[j][ix][iy][n] = 0;
	   cozEEAll[j][ix][iy][n] = 0;
	 }
	 
	 if( ! EEDetId::validDetId(ix,iy,iz) ) continue; 
	 
	 
	 EEDetId det = EEDetId(ix,iy,iz,EEDetId::XYMODE);
	 const CaloCellGeometry* cell = geom_ee->getGeometry(det);
	 GlobalPoint posi = cell->getPosition();
	 xEEAll[j][ix][iy] = posi.x();
	 yEEAll[j][ix][iy] = posi.y();
	 zEEAll[j][ix][iy] = posi.z();
	 etaEEAll[j][ix][iy] = posi.eta();
	 phiEEAll[j][ix][iy] = posi.phi();
	 
	 const CaloCellGeometry::CornersVec& co ( cell->getCorners() ) ;
	 int nc = 0;
	 for( CaloCellGeometry::CornersVec::const_iterator point = co.begin();
	      point != co.end() && nc < 8 ; ++point){
	   
	   ///cout<<"printme: " <<ix<<" "<<iy<<" "<<iz<<" pos "<< posi.x() <<" "<< posi.y() <<" "<< posi.z()  <<" co "<< point->x()<<" "<<point->y()<<" "<<point->z()<<endl; 

	   coxEEAll[j][ix][iy][nc] = point->x();
	   coyEEAll[j][ix][iy][nc] = point->y();
	   cozEEAll[j][ix][iy][nc] = point->z();
	   nc ++;
	 }
	 //       }catch(...){}
	 
       }
     }
   }
   
   }   

   ///fillin
   nsimEB = 0; 
   for(int j=-85; j<=85; j++){
     if(j==0) continue;
     for(int k=1; k<=360; k++){
       int nhits = eEBsim[j][k].size(); 
       if(nhits<1) continue; 
       vector<float> temp;
       vector<float> temp1;
       vector<float> temp2;
       double esum = 0; 
       double tmin = 1E9;
       for(int n=0; n<nhits; n++){
	 temp.push_back(eEBsim[j][k][n]);
	 esum += eEBsim[j][k][n];
	 if(tmin> tEBsim[j][k][n]){
	   tmin = tEBsim[j][k][n];
	 }

	 temp1.push_back(tEBsim[j][k][n]);
	 temp2.push_back(bEBsim[j][k][n]);
       }
       ietasimEB[nsimEB] = j;
       iphisimEB[nsimEB] = k;
       esumsimEB[nsimEB] = esum; 
       tminsimEB[nsimEB] = tmin;

       esimEB->push_back(temp);
       tsimEB->push_back(temp1);
       bsimEB->push_back(temp2);
       
       nsimEB ++; 
     }
   }

   nsimEE = 0; 
   for(int j=1;j<=100; j++){
     for(int k=1; k<=100; k++){
       int nhits = eEEmsim[j][k].size();
       if(nhits<1) continue;
       vector<float> temp;
       vector<float> temp1;
       vector<float> temp2;
       double esum = 0; 
       double tmin = 1E9;
       for(int n=0; n<nhits; n++){
         temp.push_back(eEEmsim[j][k][n]);
	 esum += eEEmsim[j][k][n];
         temp1.push_back(tEEmsim[j][k][n]);
         temp2.push_back(bEEmsim[j][k][n]);

	 if(tmin> tEEmsim[j][k][n]){
           tmin = tEEmsim[j][k][n];
         }


       }
       ixsimEE[nsimEE] = j;
       iysimEE[nsimEE] = k;
       izsimEE[nsimEE] = -1;
       esumsimEE[nsimEE] = esum; 
       tminsimEE[nsimEE] = tmin;

       esimEE->push_back(temp);
       tsimEE->push_back(temp1);
       bsimEE->push_back(temp2);
       nsimEE ++;
     }
   }
   for(int j=1;j<=100; j++){
     for(int k=1; k<=100; k++){
       int nhits = eEEpsim[j][k].size();
       if(nhits<1) continue;
       vector<float> temp;
       vector<float> temp1;
       vector<float> temp2;
       double esum = 0;
       double tmin = 1E9;
       for(int n=0; n<nhits; n++){
         temp.push_back(eEEpsim[j][k][n]);
	 esum += eEEpsim[j][k][n];
         temp1.push_back(tEEpsim[j][k][n]);
         temp2.push_back(bEEpsim[j][k][n]);

	 if(tmin> tEEpsim[j][k][n]){
           tmin = tEEpsim[j][k][n];
         }

       }
       esumsimEE[nsimEE] = esum; 
       tminsimEE[nsimEE] = tmin;
       
       ixsimEE[nsimEE] = j;
       iysimEE[nsimEE] = k;
       izsimEE[nsimEE] = 1;
       esimEE->push_back(temp);
       tsimEE->push_back(temp1);
       bsimEE->push_back(temp2);
       nsimEE ++;
     }
   }
   
   
   
  
  nEventsProcessed ++; 
   
  /////////////// ============== Generator information ================= ////////////////
  partonList.clear();
  if( !isRealData){
        
    
    Handle<GenParticleCollection> genPart;
    try {
      //iEvent.getByLabel( "genParticles", genPart);
      iEvent.getByLabel(m_genSrc,genPart);
      
      nMCpart = 0; 
      for(int i=0; i< int(genPart->size())&& nMCpart<MAXMC; i++){
	const GenParticle & p = (*genPart)[i];
	eMCpart[nMCpart] = p.energy();
	ptotMCpart[nMCpart] = p.p();
	pxMCpart[nMCpart] =p.px();
	pyMCpart[nMCpart] =p.py();
	pzMCpart[nMCpart] =p.pz();
	ptMCpart[nMCpart] =p.pt();
	etMCpart[nMCpart] = p.et();
	etaMCpart[nMCpart] = p.eta();
	phiMCpart[nMCpart] = p.phi();
	mMCpart[nMCpart] = p.mass();
	pidMCpart[nMCpart] = p.pdgId();
	chargeMCpart[nMCpart] = p.charge();
	statusMCpart[nMCpart] = p.status();
	nDauMCpart[nMCpart] = p.numberOfDaughters();
	vtxXMCpart[nMCpart] = p.vx();
	vtxYMCpart[nMCpart] = p.vy();
	vtxZMCpart[nMCpart] = p.vz();

	for(int k= 0; k<3; k++){
	  convPhtMCpart[nMCpart][k] = 0;
	}
	
	pidFDauMCpart[nMCpart] = 0; 
	if(nDauMCpart[nMCpart]>0){ ///for gluon, this is mostly 92
	  pidFDauMCpart[nMCpart] = p.daughter(0)->pdgId(); 
	}
	
	if( p.pdgId() == 25  && (p.status() ==2 || p.status()==62) ){ //two status 3 and 2 and the one with 2 is saved
	  ///pythia8  status==62 decays
	  genhiggsm = p.mass();
	  genhiggspt = p.pt();
	  genhiggseta = p.eta();
	  genhiggsphi = p.phi();
	  genhiggsstatus = p.status();
	  genhiggsvx = p.vx();
	  genhiggsvy = p.vy();
	  genhiggsvz = p.vz();
	  if(debug_>0) cout<<"mcfoundhiggs " << runNumber <<" "<<  evtNumber <<" "<<  i <<" "<< p.status()<<" "<< p.pdgId()<<" "<< p.pt()<<endl; 
	  
	  if(debug_>10) cout<<" higgs dau "<< p.numberOfDaughters() <<" "<<  p.daughter(0)->pdgId() <<" "<<  p.daughter(1)->pdgId()<<endl; 
	  

	}
	
	
	if( nMCpart >= 6){// Don't take into account first 6 particles in generator list, this list is used for quark/gluon matching
	  int pid = abs(p.pdgId()); 
	  int status = p.status();
	  if(  (( pid >=1 && pid <=5) || pid==21) && status ==3  ){
	    partonList.push_back(nMCpart);
	  }
	  else if( pidFDauMCpart[nMCpart] == 91 || pidFDauMCpart[nMCpart] ==92 ){
	    partonList.push_back(nMCpart);
	  }
	}	
	
	nMCpart++;
      }
    }catch(std::exception& ex ){
      nWarning ++; 
      if(nWarning<100) cout<<"genParticles not working.."<<endl;
    }   
    
    
    //find mother index
    for(int i=0; i< int(genPart->size())&& i <MAXMC; i++){
      const Candidate & p = (*genPart)[i];
      //int pid = p.pdgId();
      const Candidate * mom = p.mother();
      if( mom !=NULL){
	pidmomMCpart[i] = mom->pdgId();
	barcodemomMCpart[i] = indexofParticle(mom->px(),mom->pz(),mom->status());
	if(debug_>99) cout<<"indoxeofparticle: " <<  i<<" "<< p.pdgId()<<" "<< mom->pdgId() <<" "<< barcodemomMCpart[i]<<endl;
      }else{
	pidmomMCpart[i] = 0;
	barcodemomMCpart[i] = -1;
      }
    }
    
    if(debug_>9) cout<<" save all gen ele/mu/pht "<<endl; 
    
  
    
    
    //bool runsimTrack = false; //not in AOD
    bool runsimTrack = true; 
    
    if( runsimTrack ){


      // Simulated tracks (i.e. GEANT particles).
      Handle<edm::SimTrackContainer> psimtracks; 
      // Get the associated vertices
      Handle<SimVertexContainer> simvertices;
      try{
	iEvent.getByLabel("g4SimHits", simvertices);
      }catch( std::exception& ex){
	cout<<" SimVertexContainer not working.."<<endl;
      }
      
      try{
	iEvent.getByLabel("g4SimHits", psimtracks);
	
	// Need to check that SimTrackContainer is sorted; otherwise, copy and sort :-(
	std::auto_ptr<SimTrackContainer> simtracksTmp;
	const SimTrackContainer * simtracksSorted = &* psimtracks;
	
	if (!__gnu_cxx::is_sorted(psimtracks->begin(), psimtracks->end(), LessById())) {
	  simtracksTmp.reset(new SimTrackContainer(*psimtracks));
	  std::sort(simtracksTmp->begin(), simtracksTmp->end(), LessById());
	  simtracksSorted = &* simtracksTmp;
	}
	
      
	int nMCpartGen = nMCpart; 
	
	
	map<unsigned int,int> map_simtrackid; 
    
	// loop through all simParticles
	for (SimTrackContainer::const_iterator iM = psimtracks->begin(); 
	     iM != psimtracks->end(); ++iM) {
      
	  // Skip PYTHIA tracks.
	  if (iM->genpartIndex() != -1) continue; 
	
	
	  if( nMCpart >= MAXMC) break; 
	  
	  map_simtrackid.insert(make_pair(iM->trackId(),nMCpart));
      
	  pxMCpart[nMCpart] = iM->momentum().px();
	  pyMCpart[nMCpart] = iM->momentum().py();
	  pzMCpart[nMCpart] = iM->momentum().pz();
	  ptMCpart[nMCpart] = iM->momentum().pt();
	  etaMCpart[nMCpart] = iM->momentum().eta();
	  phiMCpart[nMCpart] = iM->momentum().phi();
	  pidMCpart[nMCpart] = iM->type();
      
	  chargeMCpart[nMCpart] = int(iM->charge());
	  eMCpart[nMCpart] = sqrt(0.000511* 0.000511 + iM->momentum().px()*iM->momentum().px() + iM->momentum().py()*iM->momentum().py() + iM->momentum().pz()*iM->momentum().pz()); 
	
	  statusMCpart[nMCpart] = -99; 
      
	  vtxXMCpart[nMCpart] = -9;
	  vtxYMCpart[nMCpart] = -9;
	  vtxZMCpart[nMCpart] = -9;
	
	  if (! iM->noVertex()) {
	    const SimVertex &vtx = (*simvertices)[iM->vertIndex()];
	    vtxXMCpart[nMCpart] = vtx.position().x(); 
	    vtxYMCpart[nMCpart] = vtx.position().y(); 
	    vtxZMCpart[nMCpart] = vtx.position().z(); 
	  }
      
	  barcodemomMCpart[nMCpart] = -2; 
	  nMCpart ++; 
	
	}
	
	
	vector<int> inde1;  ///+ /
	vector<int> inde2; ///- 
	vector<int> indpht; 
      
	int nSimTrack = 0; 
      
	///loop over again
	// loop through all simParticles
	for (SimTrackContainer::const_iterator iM = psimtracks->begin(); 
	     iM != psimtracks->end(); ++iM) {
      
	  // Skip PYTHIA tracks.
	  if (iM->genpartIndex() != -1) continue; 
	  if( nMCpartGen + nSimTrack >= MAXMC) break; 
	
	
	  barcodemomMCpart[nMCpartGen+nSimTrack] = -2; 
	  
	
	  // find simtrack that has a genParticle match to its parent
	  // Look at the production vertex. If there is no vertex, I can do nothing...
	  if (! iM->noVertex()) {
	    // Pick the vertex (isimtrk.vertIndex() is really an index)
	    const SimVertex &vtx = (*simvertices)[iM->vertIndex()];
	
	    // Now note that vtx.parentIndex() is NOT an index, it's a track id, so I have to search for it 
	    unsigned int idx = vtx.parentIndex();
	    SimTrackContainer::const_iterator it = std::lower_bound(simtracksSorted->begin(), simtracksSorted->end(), idx, LessById());
	    
	    if ((it != simtracksSorted->end()) && (it->trackId() == idx)) { //it is the parent sim track
	      if (it->genpartIndex() != -1) {
	    
		barcodemomMCpart[nMCpartGen+nSimTrack] = it->genpartIndex() -1; 
		if(abs(iM->type())==11 && pidMCpart[it->genpartIndex() -1] == 22) {
		  convPhtMCpart[it->genpartIndex() -1][0] = 1; ///this gen photon converts
		  if( vtxXMCpart[nMCpartGen + nSimTrack] != -999  || vtxYMCpart[nMCpartGen + nSimTrack] != -999  || vtxZMCpart[nMCpartGen + nSimTrack] != -999){
		    if(  convPhtMCpart[it->genpartIndex() -1][1] < 0) convPhtMCpart[it->genpartIndex() -1][1] = nMCpartGen + nSimTrack; ///the first electon of this gen photon converts
		    else if(  convPhtMCpart[it->genpartIndex() -1][2] < 0) convPhtMCpart[it->genpartIndex() -1][2] = nMCpartGen + nSimTrack; ///the second electon of this gen photon converts
		  }
		}
		
		vector<int>::iterator iit = find(indpht.begin(),indpht.end(),it->genpartIndex() -1);
		if( iit == indpht.end()){
		  indpht.push_back(it->genpartIndex() -1); 
		}
		if( iM->type() == 11) inde1.push_back(nMCpartGen+nSimTrack); 
		else inde2.push_back(nMCpartGen+nSimTrack); 
	      }else{
		map<unsigned int,int>::iterator intt = map_simtrackid.find(idx);
		if( intt != map_simtrackid.end()){
		  barcodemomMCpart[nMCpartGen+nSimTrack]  = intt->second;
		}else{
		  barcodemomMCpart[nMCpartGen+nSimTrack] = -3; ////check this if any
		}
	      
	      }
	    } ///found parent in simtrack
	  
	  }
	
	  nSimTrack ++; 
	
	
	} ///end of 2nd loop 
	///now let's try to find for each indpht, it's correponding two electrons index
	vector<int> pht;
	vector<int> e1;
	vector<int> e2;
      

	if(debug_ >= 2) cout<<" indpht_conv: "<< int( indpht.size()) <<endl; 

	for(int j =0; j< int(indpht.size());j++){
	  int k = indpht[j];
      
	  int found1 = -1;
	  int found2 = -1;
	  vector<int>::iterator it1= inde1.begin();
	  for(int j1 =0; j1< int(inde1.size()); j1++){
	    if(k ==   barcodemomMCpart[inde1[j1]]){
	      found1 = inde1[j1];
	      inde1.erase(it1+j1);
	      break;
	    }
	  }
	  vector<int>::iterator it2= inde2.begin();
	  for(int j2 =0; j2< int(inde2.size()); j2++){
	    if(k ==   barcodemomMCpart[inde2[j2]]){
	      found2 = inde2[j2];
	      inde2.erase(it2+j2);
	      break;
	    }
	  }
      
	  convPhtMCpart[k][1] = found1; 
	  convPhtMCpart[k][2] = found2; 
	  
	  if( debug_ >=2) cout<<"convPhtMCpart: "<< k<<" " << found1 <<" "<< found2<<" "<< eMCpart[k] << endl; 
	  
	  if( found1 > 0 && found2 >0) convPhtMCpart[k][0] = 2; /// both electron found in sim
	  else if( found1 >0 && found2 < 0) convPhtMCpart[k][0] = 3; /// only one found in sim 
	  else if( found1 <0 && found2 > 0) convPhtMCpart[k][0] = 3; /// only one found in sim 
	  else if( found1 <0 && found2 < 0) convPhtMCpart[k][0] = -10; /// neither found in sim 
	  
	}
	
      
	
      }catch( std::exception& ex){
	cout<<"SimTrackContainer  not working.."<<endl;
      }
    
    }///end of if runsimTrack
    
      //Save all gen-level electron muon, photon (with pt > 1) 
    nGenPht = 0; 
    nGenMu = 0; 
    nGenEle = 0; 
    for(int i=0; i< int(genPart->size())&& i <MAXMC; i++){
      const Candidate & p = (*genPart)[i];
	  	  
      if( ( abs(p.pdgId()) == 13 ) && nGenMu < MAXGenSaved ){
	etaGenMu[nGenMu] = p.eta();
	phiGenMu[nGenMu] = p.phi();
	ptGenMu[nGenMu] = p.pt();
	pidmomGenMu[nGenMu] = 0; 
	pidmom2GenMu[nGenMu] = 0;
	statusGenMu[nGenMu] = p.status();
	pidGenMu[nGenMu] = p.pdgId();
	vxGenMu[nGenMu] = p.vx();
	vyGenMu[nGenMu] = p.vy();
	vzGenMu[nGenMu] = p.vz();
	chaGenMu[nGenMu] = p.charge();
	
	const Candidate * mom = p.mother();
	if( mom !=NULL){
	  pidmomGenMu[nGenMu] = mom->pdgId();
	  const Candidate * mom2 = mom->mother();
	  if( mom2 !=NULL){
	    pidmom2GenMu[nGenMu] = mom2->pdgId();
	  }else{
	    pidmom2GenMu[nGenMu] = 0;
	  }
	}else{
	  pidmomGenMu[nGenMu] = 0;
	  pidmom2GenMu[nGenMu] = 0;
	  
	}

	if(debug_>1) cout<<"genmuon " << i <<" "<<endl;
	int indmom = getMotherIndex(i);
	if(indmom>=0){
	  pidmom3GenMu[nGenMu] = pidMCpart[indmom];
	  barcodemomGenMu[nGenMu] = indmom;
	}else{
	  pidmom3GenMu[nGenMu] = 0;
	  barcodemomGenMu[nGenMu] = -1;
	}
	if(debug_>1) cout<<"genmuonpidmom3 " << i <<" "<< indmom <<endl; 
	
	nGenMu ++; 
      }
      
      if( ( abs(p.pdgId()) == 11 ) && nGenEle < MAXGenSaved ){
	etaGenEle[nGenEle] = p.eta();
	phiGenEle[nGenEle] = p.phi();
	ptGenEle[nGenEle] = p.pt();
	pidmomGenEle[nGenEle] = 0; 
	pidmom2GenEle[nGenEle] = 0;
	statusGenEle[nGenEle] = p.status();
	pidGenEle[nGenEle] = p.pdgId();
	vxGenEle[nGenEle] = p.vx();
	vyGenEle[nGenEle] = p.vy();
	vzGenEle[nGenEle] = p.vz();
	chaGenEle[nGenEle] = p.charge();
	

	const Candidate * mom = p.mother();
	if( mom !=NULL){
	  pidmomGenEle[nGenEle] = mom->pdgId();
	  const Candidate * mom2 = mom->mother();
	  if( mom2 !=NULL){
	    pidmom2GenEle[nGenEle] = mom2->pdgId();
	  }else{
	    pidmom2GenEle[nGenEle] = 0;
	  }
	}else{
	  pidmomGenEle[nGenEle] = 0;
	}

	if(debug_>1) cout<<"genele " << i <<" "<<endl;
	int indmom = getMotherIndex(i);
	if(indmom>=0){
	  pidmom3GenEle[nGenEle] = pidMCpart[indmom];   ///this is for cross check, trace back all untill find the mother ( not id itself). 
	  barcodemomGenEle[nGenEle] = indmom;
	}else{
	  pidmom3GenEle[nGenEle] = 0;
	  barcodemomGenEle[nGenEle] = -1;
	}
	if(debug_>1) cout<<"genelemom3 " << i <<" "<<endl;
	nGenEle ++; 
      }
      
      
      if(p.pdgId() == 22 ){
	bool tosave = false; 
	if( p.pt()>1) tosave = true; 
	const Candidate * mom = p.mother();
	if( mom !=NULL){
	  if( mom->pdgId() <100) tosave = true; 
	}
	if( tosave  && nGenPht < MAXGenSaved ) {
	  etaGenPht[nGenPht] = p.eta();
	  phiGenPht[nGenPht] = p.phi();
	  ptGenPht[nGenPht] = p.pt();
	  pidmomGenPht[nGenPht] = 0; 
	  pidmom2GenPht[nGenPht] = 0;
	  indmom2GenPht[nGenPht] = -1;
	  statusGenPht[nGenPht] = p.status();
	  vxGenPht[nGenPht] = p.vx();
	  vyGenPht[nGenPht] = p.vy();
	  vzGenPht[nGenPht] = p.vz();
	  
	  int ind = indexofParticle(p.px(),p.pz(),p.status());
	  if(ind>=0){
	    convGenPht[nGenPht] = convPhtMCpart[ind][0];
	    vtxconvGenPht[nGenPht][0] = -999;
	    vtxconvGenPht[nGenPht][1] = -999;
	    vtxconvGenPht[nGenPht][2] = -999;
	    if(convPhtMCpart[ind][1]>0){
	      int e1 =  convPhtMCpart[ind][1];
	      vtxconvGenPht[nGenPht][0] = vtxXMCpart[e1];
	      vtxconvGenPht[nGenPht][1] = vtxYMCpart[e1];
	      vtxconvGenPht[nGenPht][2] = vtxZMCpart[e1];
	    }else if( convPhtMCpart[ind][2]>0){
	      int e1 =  convPhtMCpart[ind][2];
              vtxconvGenPht[nGenPht][0] = vtxXMCpart[e1];
              vtxconvGenPht[nGenPht][1] = vtxYMCpart[e1];
              vtxconvGenPht[nGenPht][2] = vtxZMCpart[e1];
	    }
	    

	    if(debug_==999 ){
	      if( convPhtMCpart[ind][0]==2){
		int e1 = convPhtMCpart[ind][1]; 
		int e2 = convPhtMCpart[ind][2]; 
		cout<<"photon2e "<< evtNumber<<" "<< p.eta()<<" "<<p.phi()<<" "<<p.energy()<<" "<< eMCpart[e1]<<" "<< eMCpart[e2]<<" eta/phi"<< etaMCpart[e1]<<" "<<phiMCpart[e1]<<" "<<etaMCpart[e2]<<" "<<phiMCpart[e2]<<" vtx "<< vtxXMCpart[e1]<<" "<< vtxYMCpart[e1]<<" "<< vtxZMCpart[e1] <<" "<< vtxXMCpart[e2]<<" "<< vtxYMCpart[e2]<<" "<< vtxZMCpart[e2]<<endl; 
		
	      }else{
		cout<<"photon2e "<<evtNumber<<" "<< convPhtMCpart[ind][0] <<endl; 
	      }
	    }


	  }else{
	    convGenPht[nGenPht] = -9; 
	  }
	  
	  
	  const Candidate * mom = p.mother();
	  
	  if(debug_>1) cout<<"genpht " << i <<endl; 
	  
	  if( mom !=NULL){
	    pidmomGenPht[nGenPht] = mom->pdgId();
	    int indmom = barcodemomMCpart[i];
	    
	    if(debug_>1) cout<<" phtonmom " << mom->pdgId() <<" "<< indmom <<" "<< getMotherIndex(indmom) <<endl; 
	    
	    if( indmom >=0){
	      int indmomfirst = getMotherIndex2(indmom); ///trace back to the 1st ele with FSR.
	      indmom2GenPht[nGenPht] = indmomfirst; /// the index, to check if more than one FSR photon
	      int indmom2 = getMotherIndex(indmom); 
	      pidmom2GenPht[nGenPht] = indmom2 >=0 ? pidMCpart[indmom2]: 0; ///track back to electron's mother 
	      if(debug_>1) cout<<"indmofirst " << indmomfirst <<" "<< pidmom2GenPht[nGenPht] <<endl;
	    }else{
	      indmom2GenPht[nGenPht] = -1; 
	      pidmom2GenPht[nGenPht] = 0;
	    }
	    
	  }else{
	    pidmomGenPht[nGenPht] = 0;
	    pidmom2GenPht[nGenPht] = 0;
	    indmom2GenPht[nGenPht] = -1;
	  }

	  if(debug_>1) cout<<"genpht " << i <<" "<<endl;
	  int indmom = getMotherIndex(i);
	  if(indmom>=0){
	    pidmom3GenPht[nGenPht] = pidMCpart[indmom];
	  }else{
	    pidmom3GenPht[nGenPht] = 0; 
	  }
	  if(debug_>1) cout<<"genphtmom3 " << i <<" "<<indmom <<endl; 
	  nGenPht ++; 
	  
	}
	
      }
    }
    
    
    
  }
  /////////// =================END OF Generator information ================ ////////////////////
  
  
  pileupTrueNumInterations = 0;
  pileupBunchX->clear();
  pileupNInteraction->clear();

  //PILE-UP
  if( !isRealData){
    try{
      Handle<std::vector< PileupSummaryInfo > >  PupInfo;
      iEvent.getByLabel(m_PileupSrc, PupInfo);
      std::vector<PileupSummaryInfo>::const_iterator PVI;
      for(PVI = PupInfo->begin(); PVI != PupInfo->end(); ++PVI) {
        if(debug_> 1) std::cout << " Pileup Information: bunchXing, nvtx: " << PVI->getBunchCrossing() << " " << PVI->getPU_NumInteractions() << std::endl;
        pileupBunchX->push_back(PVI->getBunchCrossing());
        pileupNInteraction->push_back(PVI->getPU_NumInteractions());
	pileupTrueNumInterations = PVI->getTrueNumInteractions(); ///Only >= 4_2_8 Fall11MC
      }
    }catch(std::exception& ex ){
      nWarning ++; 
      if(nWarning<100) cout<<"PileupSummaryInfo not working.."<<endl;
    }
  }
  

  
  //edm::Handle<reco::GsfElectronCollection> electrons;
  //edm::Handle<reco::GsfElectronCollection> electrons;
  
  //edm::Handle<edm::View<reco::Candidate> > Electrons ;
  edm::Handle<edm::View<pat::Electron> > Electrons;
  
  try{
    nElectron = 0; 
    iEvent.getByLabel(m_electronSrc,Electrons) ;
    for(edm::View<pat::Electron>::const_iterator ele=Electrons->begin(); ele!=Electrons->end(); ++ele){
      
      electronpt[nElectron] = ele->pt();
      electroneta[nElectron] = ele->eta();
      electronphi[nElectron] = ele->phi();
      electroncharge[nElectron]  = ele->charge();

      electronvertexx[nElectron] = ele->vertex().x();
      electronvertexy[nElectron] = ele->vertex().y();
      electronvertexz[nElectron] = ele->vertex().z();
      
      electronsip[nElectron] = ele->userFloat("sip");
      electronnewID[nElectron] = ele->userInt("newID");
      electronmvaID[nElectron] = ele->userInt("mvaID");
      electronbdtID[nElectron] = ele->userFloat("bdtID");
      
      reco::SuperClusterRef scRef = ele->superCluster();
      electronsceta[nElectron] = scRef->eta();
      electronscphi[nElectron] = scRef->phi();
      reco::GsfTrackRef trackref = ele->gsfTrack();
      if(trackref.isNonnull()){
	const reco::HitPattern& p_inner = trackref->trackerExpectedHitsInner();
	electronExpectedHitsInnernumberOfHits[nElectron] = p_inner.numberOfHits();
	electronExpectedHitsOuternumberOfHits[nElectron] = trackref->trackerExpectedHitsOuter().numberOfHits();
      }else{
	electronExpectedHitsInnernumberOfHits[nElectron] = -99;
	electronExpectedHitsOuternumberOfHits[nElectron] = -99;
      }
      electronpfCombRelIso04EACorr[nElectron] = ele->userFloat("pfCombRelIso04EACorr");
      
      nElectron ++; 
    }
        
  }catch( std::exception& ex){
    nWarning ++; 
    if(nWarning<100) cout<<"electron  not working.."<<endl;
  }
  
  edm::Handle<edm::View<pat::Muon> > Muons;
  try{
    nMuon = 0; 
    iEvent.getByLabel(m_muonSrc,Muons) ;

    for(edm::View<pat::Muon>::const_iterator mu=Muons->begin(); mu!=Muons->end(); ++mu){
      
      muonpt[nMuon] = mu->pt();
      muoneta[nMuon] = mu->eta();
      muonphi[nMuon] = mu->phi();
      muoncharge[nMuon] =  mu->charge();
      muonvx[nMuon] = mu->vx();
      muonvy[nMuon] = mu->vy();
      muonvz[nMuon] = mu->vz();
            
      muonpfMuId[nMuon] = mu->userInt("pfMuId");
      muonpfCombRelIso04EACorr[nMuon] = mu->userFloat("pfCombRelIso04EACorr");
      muonsip[nMuon] = mu->userFloat("sip");
      
      nMuon ++; 
    }
    
  }catch( std::exception& ex){
    nWarning ++; 
    if(nWarning<100) cout<<"muon  not working.."<<endl;
  }
 

  edm::Handle<edm::View<pat::PFParticle> > Photons;
  try{
    nPhoton = 0; 
    iEvent.getByLabel(m_photonSrc,Photons) ;
    
    for(edm::View<pat::PFParticle>::const_iterator pht=Photons->begin(); pht!=Photons->end(); ++pht){
      
      photonpt[nPhoton] = pht->pt();
      photoneta[nPhoton] = pht->eta();
      photonphi[nPhoton] = pht->phi();
      
      photonhasOverlapseleVeto[nPhoton] = pht->hasOverlaps("eleVeto");
      photonhasOverlapsgoodLepIso[nPhoton] = pht->hasOverlaps("goodLepIso");
      
      if(debug_>1) cout<<"photon overlap " << photonhasOverlapseleVeto[nPhoton] <<endl; 
      
      
      nPhoton ++; 
    }
    
  }catch( std::exception& ex){
    nWarning ++; 
    if(nWarning<100)cout<<"photon  not working.."<<endl;
  }
  


  ///Fill the branches
  Analysis->Fill();
  
  
}


// ------------ method called once each job just before starting event loop  ------------
void 
RecoAnalyzer::beginJob()
{


  cout<<" RecoAnalyzer::beginJob( " <<endl; 
  

  
  Analysis = new TTree("Analysis","Reco Analysis");

  bool saveGeo = false;
  //bool saveGeo = true; 

  if(saveGeo){
  
  Analysis->Branch("xEBAll",xEBAll,"xEBAll[170][360]/D");
  Analysis->Branch("yEBAll",yEBAll,"yEBAll[170][360]/D");
  Analysis->Branch("zEBAll",zEBAll,"zEBAll[170][360]/D");
  Analysis->Branch("etaEBAll",etaEBAll,"etaEBAll[170][360]/D");
  Analysis->Branch("phiEBAll",phiEBAll,"phiEBAll[170][360]/D");
  
  Analysis->Branch("coxEBAll",coxEBAll,"coxEBAll[170][360][8]/D");
  Analysis->Branch("coyEBAll",coyEBAll,"coyEBAll[170][360][8]/D");
  Analysis->Branch("cozEBAll",cozEBAll,"cozEBAll[170][360][8]/D");
  

  Analysis->Branch("xEEAll",xEEAll,"xEEAll[2][101][101]/D");
  Analysis->Branch("yEEAll",yEEAll,"yEEAll[2][101][101]/D");
  Analysis->Branch("zEEAll",zEEAll,"zEEAll[2][101][101]/D");
  Analysis->Branch("etaEEAll",etaEEAll,"etaEEAll[2][101][101]/D");
  Analysis->Branch("phiEEAll",phiEEAll,"phiEEAll[2][101][101]/D");

  Analysis->Branch("coxEEAll",coxEEAll,"coxEEAll[2][101][101][8]/D");
  Analysis->Branch("coyEEAll",coyEEAll,"coyEEAll[2][101][101][8]/D");
  Analysis->Branch("cozEEAll",cozEEAll,"cozEEAll[2][101][101][8]/D");

  ///xESAll[2][2][32][40][40];
  Analysis->Branch("xESAll",xESAll,"xESAll[2][2][32][40][40]/D");
  Analysis->Branch("yESAll",yESAll,"yESAll[2][2][32][40][40]/D");
  Analysis->Branch("zESAll",zESAll,"zESAll[2][2][32][40][40]/D");
  
  }


  

  pileupBunchX = new std::vector<short>; pileupBunchX->clear();
  pileupNInteraction = new std::vector<short>; pileupNInteraction->clear();
  //nPileUP
  Analysis->Branch("pileupBunchX","std::vector<short>", &pileupBunchX);
  Analysis->Branch("pileupNInteraction","std::vector<short>", &pileupNInteraction);
  Analysis->Branch("pileupTrueNumInterations", &pileupTrueNumInterations,"pileupTrueNumInterations/F");
  

  //g4simHits
  
  Analysis->Branch("nsimEB",&nsimEB,"nsimEB/I");
  Analysis->Branch("ietasimEB",ietasimEB,"ietasimEB[nsimEB]/I");
  Analysis->Branch("iphisimEB",iphisimEB,"iphisimEB[nsimEB]/I");
  Analysis->Branch("esumsimEB",esumsimEB,"esumsimEB[nsimEB]/F");
  Analysis->Branch("tminsimEB",tminsimEB,"tminsimEB[nsimEB]/F");

  
  esimEB = new std::vector<std::vector<float> >; esimEB->clear();
  tsimEB = new std::vector<std::vector<float> >; tsimEB->clear();
  bsimEB = new std::vector<std::vector<float> >; bsimEB->clear();

  bool saveme = false; 
  //bool saveme = true; 
  if(saveme){
    Analysis->Branch("esimEB", "std::vector<std::vector<float> >", &esimEB);
    Analysis->Branch("tsimEB", "std::vector<std::vector<float> >", &tsimEB);
    Analysis->Branch("bsimEB", "std::vector<std::vector<float> >", &bsimEB);
  }

  Analysis->Branch("nsimEE",&nsimEE,"nsimEE/I");
  Analysis->Branch("ixsimEE",ixsimEE,"ixsimEE[nsimEE]/I");
  Analysis->Branch("iysimEE",iysimEE,"iysimEE[nsimEE]/I");
  Analysis->Branch("izsimEE",izsimEE,"izsimEE[nsimEE]/I");
  Analysis->Branch("esumsimEE",esumsimEE,"esumsimEE[nsimEE]/F");
  Analysis->Branch("tminsimEE",tminsimEE,"tminsimEE[nsimEE]/F");


  esimEE = new std::vector<std::vector<float> >; esimEE->clear();
  tsimEE = new std::vector<std::vector<float> >; tsimEE->clear();
  bsimEE = new std::vector<std::vector<float> >; bsimEE->clear();
  if(saveme){
    Analysis->Branch("esimEE", "std::vector<std::vector<float> >", &esimEE);
    Analysis->Branch("tsimEE", "std::vector<std::vector<float> >", &tsimEE);
    Analysis->Branch("bsimEE", "std::vector<std::vector<float> >", &bsimEE);
  }


  //Event info
  Analysis->Branch("lumiBlock",&lumiBlock,"lumiBlock/I");
  Analysis->Branch("runNumber",&runNumber,"runNumber/I");
  Analysis->Branch("evtNumber",&evtNumber,"evtNumber/I");
  Analysis->Branch("bunchX",&bunchX,"bunchX/I");
  Analysis->Branch("orbitNumber",&orbitNumber,"orbitNumber/I");
  Analysis->Branch("evtTime",&evtTime,"evtTime/I");
  Analysis->Branch("isRealData",&isRealData,"isRealData/I");
    
  /////Photon
  Analysis->Branch("nPhoton",&nPhoton, "nPhoton/I");

  Analysis->Branch("photonpt",photonpt,"photonpt[nPhoton]/F");
  Analysis->Branch("photonenergy",photonenergy,"photonenergy[nPhoton]/F");
  Analysis->Branch("photoneta",photoneta,"photoneta[nPhoton]/F");
  Analysis->Branch("photonphi",photonphi,"photonphi[nPhoton]/F");
  
  Analysis->Branch("photonhasOverlapseleVeto",photonhasOverlapseleVeto,"photonhasOverlapseleVeto[nPhoton]/I");
  Analysis->Branch("photonhasOverlapsgoodLepIso",photonhasOverlapsgoodLepIso,"photonhasOverlapsgoodLepIso[nPhoton]/I");
  
  
  ///electron
  Analysis->Branch("nElectron",&nElectron, "nElectron/I");
  Analysis->Branch("electronpt",electronpt,"electronpt[nElectron]/F");
  Analysis->Branch("electroneta",electroneta,"electroneta[nElectron]/F");
  Analysis->Branch("electronphi",electronphi,"electronphi[nElectron]/F");
  Analysis->Branch("electroncharge",electroncharge,"electroncharge[nElectron]/I");
  
  Analysis->Branch("electronvertexx",electronvertexx,"electronvertexx[nElectron]/F");
  Analysis->Branch("electronvertexy",electronvertexy,"electronvertexy[nElectron]/F");
  Analysis->Branch("electronvertexz",electronvertexz,"electronvertexz[nElectron]/F");
    
  Analysis->Branch("electronsceta",electronsceta,"electronsceta[nElectron]/F");
  Analysis->Branch("electronscphi",electronscphi,"electronscphi[nElectron]/F");

  Analysis->Branch("electronsip",electronsip,"electronsip[nElectron]/F");
  Analysis->Branch("electronnewID",electronnewID,"electronnewID[nElectron]/I");
  Analysis->Branch("electronmvaID",electronmvaID,"electronmvaID[nElectron]/I");  
  Analysis->Branch("electronbdtID",electronbdtID,"electronbdtID[nElectron]/F");  
  Analysis->Branch("electronpfCombRelIso04EACorr",electronpfCombRelIso04EACorr,"electronpfCombRelIso04EACorr[nElectron]/F");
  Analysis->Branch("electronExpectedHitsInnernumberOfHits",electronExpectedHitsInnernumberOfHits,"electronExpectedHitsInnernumberOfHits[nElectron]/I");
  Analysis->Branch("electronExpectedHitsOuternumberOfHits",electronExpectedHitsOuternumberOfHits,"electronExpectedHitsOuternumberOfHits[nElectron]/I");
  
  
  ///Muon
  Analysis->Branch("nMuon",&nMuon,"nMuon/I");
  Analysis->Branch("muonpt",muonpt,"muonpt[nMuon]/F");
  Analysis->Branch("muoncharge",muoncharge,"muoncharge[nMuon]/I");
  Analysis->Branch("muoneta",muoneta,"muoneta[nMuon]/F");
  Analysis->Branch("muonphi",muonphi,"muonphi[nMuon]/F");
  Analysis->Branch("muonpfMuId",muonpfMuId,"muonpfMuId[nMuon]/I");
  Analysis->Branch("muonsip",muonsip,"muonsip[nMuon]/F");
  Analysis->Branch("muonpfCombRelIso04EACorr",muonpfCombRelIso04EACorr,"muonpfCombRelIso04EACorr[nMuon]/F");
  Analysis->Branch("muonvx",muonvx,"muonvx[nMuon]/F");
  Analysis->Branch("muonvy",muonvy,"muonvy[nMuon]/F");
  Analysis->Branch("muonvz",muonvz,"muonvz[nMuon]/F");
  

  cout<<" RecoAnalyzer tree branches defined.. " <<endl; 
  
  

}



////test the particle-to-parton Matching
///retrun index in MCpart array
int RecoAnalyzer::partonMatchingAlgo(float eta, float phi, float dr){
  
  
  int tempParticle = -1;
  int tempPartonHighestPt = -1;
  int tempNearest = -1;
  float maxPt = 0;
  float minDr = 1000;
  ///bool foundPriority = false;
  
  float coneSizeToAssociate = dr; 
  
  if(  nMCpart <1) return -99; 
  if( int(partonList.size()) <1) return -99; 
  
    
  for (int m = 0; m<int(partonList.size()); m++){
    
    int j = partonList[m];
    
    int pid = abs(pidMCpart[j]);
    
    /// if( pid >= 6 && pid != 21 ) continue; 
    float dist = GetDeltaR(eta,etaMCpart[j],phi,phiMCpart[j]);
    
    if( pidFDauMCpart[j] == 91 || pidFDauMCpart[j] == 92){
      
      
      if( dist < coneSizeToAssociate ){
	if( dist < minDr ) {
	  minDr = dist;
	  tempNearest = j;
	}
	
	if( tempParticle == -1 && pid == 4  ) tempParticle = j;
	if(  pid == 5    ) tempParticle = j;
	if( ptMCpart[j] > maxPt ) {
	  maxPt = ptMCpart[j];
	  tempPartonHighestPt = j;
	}
      }
    }
  }
  
  if(tempParticle ==-1) tempParticle = tempPartonHighestPt; 
  
  if( tempParticle <0) return 0; 
  else return tempParticle; 
    
  
}


int RecoAnalyzer::findIndexofTRK(float pt, float eta, float phi){
  
  float drmin = 0.00001; 
  int flag = -1; 
  for(int j=0; j< nTrack; j++){
    if( fabs( pt - trackpt[j]) > 0.00001) continue; 
    float dr = GetDeltaR(eta,tracketa[j],phi,trackphi[j]);
    if( dr <drmin){
      drmin = drmin; 
      flag = j; 
    }
  }
  return flag;
  
}


int RecoAnalyzer::findIndexofSC(float en,float eta,float phi){
  int flag =-1; 
  int nfound =0;
  float err = 0.0001; 
  for( int j=0; j<nSC; j++){
    if( fabs(eSC[j]-en)<err && fabs(etaSC[j]-eta)<err && fabs(phiSC[j]-phi)<err){
      flag = j; 
      nfound++; 
    }
  }
  return flag; 
}


int RecoAnalyzer::getMotherIndex(int j){
  if(j<0) {
    cout<<"no mother. input -1!!!!!!!!"<<endl;
    return -1; 
  }
  
  int indmom = barcodemomMCpart[j];
  int nloop = 0; 
  while( indmom >=0 && pidMCpart[indmom]==pidMCpart[j] ){ //if it is the same pid 
    //int st = statusMCpart[indmom];
    indmom = barcodemomMCpart[indmom];
    nloop ++;
    if( nloop > nMCpart ){
      cout<<" RecoAnalyzer::getMotherIndex deadloop! " << evtNumber <<endl; 
      break; 
    }
  }
  
  return indmom; 
  
  
}


///trace back till the last dauther of the same mother 
///Input .e.g electron. trace back to the first one undergo FSR..
int RecoAnalyzer::getMotherIndex2(int j){
  if(j<0) {
    cout<<"no mother. input -1!!!!!!!!"<<endl;
    return -1; 
  }
  
  int indmom = barcodemomMCpart[j];
  int inddau = j; 
  int nloop = 0; 
  while( indmom >=0 && pidMCpart[indmom] ==pidMCpart[j]){//the same 
    inddau = indmom;
    indmom = barcodemomMCpart[indmom];
    nloop ++; 
    if( nloop > nMCpart){
      cout<<" RecoAnalyzer::getMotherIndex2 deadloop? "<< lumiBlock <<" "<< runNumber  <<" "<< evtNumber <<endl;
      break; 
    }
  }
  return inddau; 
    
}




void RecoAnalyzer::MatchToGenElectron(float eta, float phi, float res[]){
  
  for(int j=0; j<20; j++) res[j] = 1; 
  float drMin = 0.3; 
  int indMin = -1; 

  float drMin1 = 0.3; 
  int indMin1 = -1; 

  
  float drMin2 = 0.3; 
  int indMin2 = -1;

  
  if( nMCpart <1) return; 

  for(int j=1; j<nMCpart; j++){

    float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
    
    if( ptMCpart[j] > 3 && statusMCpart[j] ==1 && abs(pidMCpart[j])==11){ ///gen real electron
      if(dr <drMin){
	drMin = dr; 
	indMin = j; 
      }
    }
    
    if( ptMCpart[j] > 3 && statusMCpart[j] == 1 &&  abs(chargeMCpart[j]) == 1  && abs(pidMCpart[j]) != 11){ //gen charged, not electron  
      if(dr <drMin1){
	drMin1 = dr; 
	indMin1 = j; 
      }
    }
    if(statusMCpart[j] == -99 && barcodemomMCpart[j]>0 && ptMCpart[j] >3 && abs(pidMCpart[j]) == 11){ ///from photon conversion
      if( statusMCpart[barcodemomMCpart[j]] ==1 && pidMCpart[barcodemomMCpart[j]] == 22){
	if( dr < drMin2){
	  drMin2 = dr; 
	  indMin2 = j; 
	}
      }
      
    } 
  }
  

  if( indMin>0){ ///matched to a gen-electron
    res[0] = drMin; 
    res[1] = ptMCpart[indMin];
    res[2] = pidMCpart[indMin];
    int indm = getMotherIndex(indMin); 
    if(indm >=0){
      res[3] = pidMCpart[indm];
    }else{
      res[3] = 0; 
    }
  }
  
  if( indMin1>0){///matched to a gen-charged
    res[4] = drMin1; 
    res[5] = ptMCpart[indMin1];
    res[6] = pidMCpart[indMin1];
    int indm = getMotherIndex(indMin1); 
    if(indm >=0){
      res[7] = pidMCpart[indm];
    }else{
      res[7] = 0; 
    }
  }
  if( indMin2>0){ ///matched to electron from photon conv.
    res[8] = drMin2; 
    res[9] = ptMCpart[indMin2];
    res[10] = pidMCpart[indMin2];
    res[11] = vtxXMCpart[indMin2];
    res[12] = vtxYMCpart[indMin2];
    res[13] = vtxZMCpart[indMin2];
    int indmm = getMotherIndex(barcodemomMCpart[indMin2]);
    if(indmm >=0){
      res[14] = pidMCpart[indmm]; ///photon's mother
    }else{
      res[14] = 0; 
    }
  }
    
}

void RecoAnalyzer::MatchToGenPhoton(float eta, float phi, float res[]){
  
  
  float drMin = 0.2; 
  int indMin = -1; 
  for(int j=0; j<20; j++){
    res[j] = 1; 
  }
  
  if( nMCpart <1) return; 
  
  for(int j=0; j<nMCpart; j++){
    if(statusMCpart[j] !=1) continue;
    if( ptMCpart[j]< 3 ) continue;  ///don't matching to very small pt photon's
    if( abs(pidMCpart[j]) != 22 ) continue; 
    float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);

    if(dr <drMin){
      drMin = dr; 
      indMin = j; 
    }
    
  }
  
  
  if(indMin>=0){
    res[0] = drMin; 
    ///res[1] = ptMCpart[indMin]/sin(2*atan(exp(-etaMCpart[indMin])));
    //res[2] = etaMCpart[indMin];
    //res[3] = phiMCpart[indMin];
    res[1] = ptMCpart[indMin];
    ///mother's index and pid
    ///for qcd events,45->67.
    int indmom = getMotherIndex(indMin);
    //res[5] = indmom;
    if( indmom>=0){
      res[2] = pidMCpart[indmom];
    }else{
      res[2] = 0;
    }
    ///NOTE pythia 8 have photon's mother id with 22 , so change from 1,6, 21 13, 11, to >=1 && <30 
    ///if it's pi0/eta check decay angle of two photon
    ///if it's not from ISR or FSR.
    
    
  }
  
  float drMin1 = 0.2; 
  int indMin1 = -1; 
  
  for(int j=0; j<nMCpart; j++){
    
    if(statusMCpart[j] !=1) continue;
    if( abs(pidMCpart[j]) !=11) continue;  ///matched to electron
    
    if( ptMCpart[j]<10) continue; 
    
    float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
    
    if(dr <drMin1){
      drMin1 = dr; 
      indMin1 = j; 
    }
  }
  
  
  if( indMin1>0){ //matched to a real electron
    res[7] = drMin1; 
    res[8] = ptMCpart[indMin1];
    int indm = getMotherIndex(indMin1); 
    if(indm >=0){
      res[9] = pidMCpart[indm];
    }else{
      res[9] = 0; 
    }
  }
  
  
  if(indMin>=0){
    if( convPhtMCpart[indMin][0] > 0){
      //conv vtx
      res[10] = 1; 
      res[11] = vtxXMCpart[convPhtMCpart[indMin][1]]; 
      res[12] = vtxYMCpart[convPhtMCpart[indMin][1]]; 
      res[13] = vtxZMCpart[convPhtMCpart[indMin][1]]; 
    }
    else {
      res[10] = 0; 
      res[11] = 0; 
      res[12] = 0; 
      res[13] = 0; 
    }
  }else{
    res[10] = 0; 
    res[11] = 0; 
    res[12] = 0; 
    res[13] = 0; 
  }
  
  
  
}



int RecoAnalyzer::indexofParticle(float px, float pz, int status){
  
  float err = 0.00001; 
  
  
  int indMin = -1; 
  float errMin = 0.002; 
  for( int j=0; j<nMCpart; j++){
    if(fabs(px-pxMCpart[j])<err && fabs(pz-pzMCpart[j])<err
       && status == statusMCpart[j]){
      float err1 = fabs(px-pxMCpart[j]) + fabs(pz-pzMCpart[j]);
      if( err1 < errMin){
	errMin = err1; 
	indMin = j; 
      }
    }
  }
  
  return indMin; 
  
  //   if(int(nn.size()) !=1){ ///diquarks or gluon has two copy sometimesss!!! 
  //     cout<<"wrong! indexofParticle: n: "<<int(nn.size())<<endl;
  //     for(int j=0; j<int(nn.size()); j++){
  //       cout<<"  dup j: "<<j<<" ind "<<nn[j]<<endl;
  //     }
  //   }
  //return nn[0]; 
  
}


double RecoAnalyzer::DeltaPhi(double v1, double v2)
{ // Computes the correctly normalized phi difference
  // v1, v2 = phi of object 1 and 2
  double diff = v1 - v2;
  
  while (diff >acos(-1)) diff -= 2*acos(-1);
  while (diff <= -acos(-1)) diff += 2*acos(-1);
  
  return diff; 
}


double RecoAnalyzer::GetDeltaR(double eta1, double eta2, double phi1, double phi2){ 
  // Computes the DeltaR of two objects from their eta and phi values
  
  return sqrt( (eta1-eta2)*(eta1-eta2) + DeltaPhi(phi1, phi2)*DeltaPhi(phi1, phi2) );
  
}


void RecoAnalyzer::MatchToGenMuon(float eta, float phi, float res[]){
  
  float drMin = 0.3; 
  int indMin = -1; 
  
  float drMin1 = 0.3; 
  int indMin1 = -1; 
  
  float drMin2 = 0.3; 
  int indMin2 = -1; 
  
  for(int j=0; j<20; j++){
    res[j] = 1; 
  }

  if( nMCpart <1) return; 
    

  for(int j=6; j<nMCpart; j++){
    
    if( ptMCpart[j] < 3 ) continue; 
        
    ///gen real muon 
    if( abs(pidMCpart[j]) == 13 && statusMCpart[j] ==1 ){
      float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
      if(dr <drMin){
	drMin = dr; 
	indMin = j; 
      }
    }
    
    ///gen real charged non-muon , non-eletron particles
    if(statusMCpart[j] ==1 && abs(chargeMCpart[j]) ==1 && abs(pidMCpart[j]) != 13 && abs(pidMCpart[j]) !=11){
      float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
     if(dr <drMin1){
	drMin1 = dr; 
	indMin1 = j; 
     }
    }
    ////sim real moun from gen
    if(statusMCpart[j]== -99 && abs(pidMCpart[j]) == 13 && barcodemomMCpart[j] >0 ){
      float dr = GetDeltaR(etaMCpart[j],eta,phiMCpart[j],phi);
      if( dr< drMin2){
	drMin2 = dr; 
	indMin2 = j; 
      }
    }
    
  }
  
  
  if(indMin>0){ ///if found a real muon in gen-level
    res[0] = drMin; 
    res[1] = ptMCpart[indMin];
    res[2] = pidMCpart[indMin];

    int indm = getMotherIndex(indMin); 
    if( indm >=0) {
      res[3] = pidMCpart[indm];
    }else{
      res[3] = 0; 
    }
    
  }
  
  if(indMin1>0){ ///if found a charged particle in gen-level
    res[4] = drMin1; 
    res[5] = ptMCpart[indMin1];
    res[6] = pidMCpart[indMin1];
    int indm = getMotherIndex(indMin1);
    if( indm>=0){
      res[7] = pidMCpart[indm];
    }else{
      res[7] = 0; 
    }
    
  }
  
  if( indMin2>0){ ///if found a muon in sim-level
    res[8] = drMin2; 
    res[9] = ptMCpart[indMin2];
    res[10] = pidMCpart[indMin2];
    res[11] = vtxXMCpart[indMin2];
    res[12] = vtxYMCpart[indMin2];
    res[13] = vtxZMCpart[indMin2];
    int indm = getMotherIndex(indMin2);
    if( indm>=0){
      ///res[14] = pidMCpart[barcodemomMCpart[indMin2]]; /// muon's mother, mostly should be pi0/kion decay with staus ==1 
      res[14] =  pidMCpart[indm];
    }else{
      res[14] = 0; 
    }
  }
  
  
}


double RecoAnalyzer::etaTransformation(  float EtaParticle , float Zvertex)  {
  
  //---Definitions
  const float pi = 3.1415927;

  //---Definitions for ECAL
  const float R_ECAL           = 136.5;
  const float Z_Endcap         = 328.0;
  const float etaBarrelEndcap  = 1.479; 
   
  //---ETA correction

  float Theta = 0.0  ; 
  float ZEcal = R_ECAL*sinh(EtaParticle)+Zvertex;

  if(ZEcal != 0.0) Theta = atan(R_ECAL/ZEcal);
  if(Theta<0.0) Theta = Theta+pi ;
  double ETA = - log(tan(0.5*Theta));
         
  if( fabs(ETA) > etaBarrelEndcap )
    {
      float Zend = Z_Endcap ;
      if(EtaParticle<0.0 )  Zend = -Zend ;
      float Zlen = Zend - Zvertex ;
      float RR = Zlen/sinh(EtaParticle); 
      Theta = atan(RR/Zend);
      if(Theta<0.0) Theta = Theta+pi ;
      ETA = - log(tan(0.5*Theta));		      
    } 
  //---Return the result
  return ETA;
  //---end
}



//////////from HiggsAnalysis/HiggsTo2photons/h2gglobe/GeneralFunctions_cc.h
////retrun the index to assced conv..[]
int  RecoAnalyzer::matchPhotonToConversion( int lpho) {

  int result=-99;
  double conv_eta=-999.;
  double conv_phi=-999.;
  
  //float sc_eta  = ((TVector3 *) pho_calopos->At(lpho))->Eta();
  //  float sc_phi  = ((TVector3 *) pho_calopos->At(lpho))->Phi();
  float sc_eta = photonsceta[lpho];
  float sc_phi = photonscphi[lpho];
  

  //TLorentzVector * p4 = (TLorentzVector *) pho_p4->At(lpho);
  // float et = p4->Energy() / cosh(sc_eta);
  //float et = photonenergy[lpho] / cosh(sc_eta);
  
  //  cout << " photon index " << lpho << " eta " <<sc_eta<< " phi " << sc_phi << " et " << et << endl; 
  
  float detaMin=999.;
  float dphiMin=999.;   
  ///float dRMin = 999.;

  float mconv_pt=-999999;
  int iMatch=-1;     
  
  if(debug_ > 10 )  cout << "   LoopAll::matchPhotonToConversion conv_n " << nConv << endl; 
  for(int iconv=0; iconv<nConv; iconv++) {
    
    //vector<TVector3> refittedPairMomentum

    ///refittedPairMomentum.push_back(*((TVector3*) conv_refitted_momentum->At(iconv)));
    
    TVector3 refittedPairMomentum(convrefittedPairMomentumx[iconv],convrefittedPairMomentumy[iconv],convrefittedPairMomentumz[iconv]);
    float conv_pt =  refittedPairMomentum.Pt();
    if (conv_pt < 1 ) continue;    
    if ( convconversionVertexisValid[iconv] == 0 ) continue;
    
    conv_phi  = refittedPairMomentum.Phi();
    float eta  = refittedPairMomentum.Eta();
    conv_eta = etaTransformation(eta, photonconversionzOfPrimaryVertexFromTracks[iconv] ); //transfer eta to eta at ECAL 

    //float dPhi =conv_phi - sc_phi;       
    //double delta_phi = phiNorm (dPhi);
    double delta_phi = DeltaPhi(conv_phi, sc_phi);
    double delta_eta = conv_eta - sc_eta;
    
    //cout << " delta_eta " << delta_eta << " delta_phi " << delta_phi << endl;
    //delta_phi=pow(delta_phi,2);
    //delta_eta=pow(delta_eta,2);
    //float dR = sqrt( delta_phi+delta_eta); 
    
    if ( fabs(delta_eta) < detaMin && fabs(delta_phi) < dphiMin ) {  ///why not using a deltaR matching  ??
      //    if ( dR < dRMin ) {
      detaMin=  fabs(delta_eta);
      dphiMin=  fabs(delta_phi);
      //dRMin=dR;
      iMatch=iconv;
      mconv_pt = conv_pt;
    }
    
  }
  //  cout << " minimized conversion index " << iMatch << " eta " <<conv_eta<< " phi " << conv_phi <<endl; 
  if ( detaMin < 0.1 && dphiMin < 0.1 ) {
    //  if ( dRMin< 0.1 ) {
    if(debug_)    cout << " matched conversion index " << iMatch << " eta " <<conv_eta<< " phi " << conv_phi << " pt " << mconv_pt << endl; 	
    result = iMatch;
  } else {
    result = -1;
  }
  
  return result;
    
}



///from  HggVertexFromConversions::vtxZ(const PhotonInfo & pho)
double RecoAnalyzer::HggVertexFromConversionsvtxZ(int j)
{
  
  bool useAllConvs= true; 
  
  
  double phtConX = 0; 
  double phtConY = 0; 
  double phtConZ = 0; 
  if( useAllConvs && photonmatchToallConv[j] >=0){
    phtConX = convconversionVertexx[photonmatchToallConv[j]];
    phtConY = convconversionVertexy[photonmatchToallConv[j]];
    phtConZ = convconversionVertexz[photonmatchToallConv[j]];
  }else{
    phtConX = photonconversionVertexx[j];
    phtConY = photonconversionVertexy[j];
    phtConZ = photonconversionVertexz[j];
  }
  
  // get the z from conversions
  double deltaX1 =  photoncaloPositionx[j] - phtConX;
  double deltaY1 =  photoncaloPositiony[j] - phtConY;
  double deltaZ1 =  photoncaloPositionz[j] - phtConZ;
  double R1 = sqrt(deltaX1*deltaX1+deltaY1*deltaY1);
  double tantheta = R1/deltaZ1;
  
  double deltaX2 = phtConX-beamSpotX; 
  double deltaY2 = phtConY-beamSpotY; 
  double R2 = sqrt(deltaX2*deltaX2+deltaY2*deltaY2);
  double deltaZ2 = R2/tantheta;
  double higgsZ =  photoncaloPositionz[j]-deltaZ1-deltaZ2;
  return higgsZ;
  
}



///from  HggVertexFromConversions::vtxdZ(const PhotonInfo & pho)
double RecoAnalyzer::HggVertexFromConversionsvtxdZ(int j){
  
  ///// attribute the error depending on the tracker region
  double dz=-99999;
  
  bool useAllConvs= true; 
  
  double phtConX = 0; 
  double phtConY = 0; 
  double phtConZ = 0; 
  if( useAllConvs && photonmatchToallConv[j] >=0){
    phtConX = convconversionVertexx[photonmatchToallConv[j]];
    phtConY = convconversionVertexy[photonmatchToallConv[j]];
    phtConZ = convconversionVertexz[photonmatchToallConv[j]];
  }else{
    phtConX = photonconversionVertexx[j];
    phtConY = photonconversionVertexy[j];
    phtConZ = photonconversionVertexz[j];
  }
  
  float sigmaPix_ =0.06;
  float sigmaTib_ =0.67;
  float sigmaTob_ =2.04;
  float sigmaFwd1_ =0.18;
  float sigmaFwd2_ =0.61;
  float sigmaFwd3_ =0.99;
  
  
  ///if ( pho.iDet() ==1 ) { // barrel
  if ( fabs(photonsceta[j])<1.48 ) { // barrel

    double perp = sqrt( phtConX*phtConX + phtConY*phtConY);
    if ( perp <=15 ) {
      dz=sigmaPix_;
    } else if ( perp > 15 && perp <=60 ) {
      dz=sigmaTib_;
    } else {
      dz=sigmaTob_;
    }

  } else { // endcap
    
    if ( fabs(phtConZ ) <=50 ) {
      dz=sigmaFwd1_;
    } else if (phtConZ > 50 &&  phtConZ <= 100 ) {
      dz=sigmaFwd2_;
    } else {
      dz=sigmaFwd3_;
    }
  }
  
  return dz;
  
}


TLorentzVector RecoAnalyzer::photonp4wrtvertex(int indpht, int indvtx){
  
  TVector3 vPos(vertexx[indvtx],vertexy[indvtx],vertexz[indvtx]);
  TVector3 caloPosition(photoncaloPositionx[indpht],photoncaloPositiony[indpht],photoncaloPositionz[indpht]);
  TVector3 direction = caloPosition - vPos;
  TVector3 p = direction.Unit() * photonenergy[indpht];
  TLorentzVector p4(p.x(),p.y(),p.z(),photonenergy[indpht]);
  return p4;
  
}



////from HiggsAnalysis/HiggsTo2photons/h2gglobe/GeneralFunctions_cc.h
//// LoopAll::vertexSelection(HggVertexAnalyzer & vtxAna, HggVertexFromConversions & vtxAnaFromConv, 
////					  PhotonInfo & pho1, PhotonInfo & pho2, std::vector<std::string> & vtxVarNames)

//Input index of the two photon
void RecoAnalyzer::higgsVertexAnalysis(int indpht1, int indpht2){
  

  ///first build quantities /// inside HggVertexAnalyzer::analyze()
  //analyze(const VertexInfoAdapter & e, const PhotonInfo & p1, const PhotonInfo & p2)
  
  std::vector<unsigned short> vtxTracksBuf;
  std::vector<int> vtxTracksSizeBuf;
  //if( ! e.hasVtxTracks() ) {
  if( false ) {
 //    vtxTracksBuf.resize(nvtx*e.ntracks());
//     vtxTracksSizeBuf.resize(nVertex,0);
//     for(int it=0; it<e.ntracks(); ++it) {
//       int vid = e.tkVtxId(it);
      
//       if( it < 5) cout<<" tkVtxId " << it <<" "<<  vid <<endl; 
      
//       int & ntks = vtxTracksSizeBuf[vid];
//       vtxTracksBuf[ vid*e.ntracks() + ntks ] = it;
//       ++ntks;
//     }
  }
  
  vector<int> preselection_; 
  preselection_.clear();

  bool usehighPurityOnly = false; 
    
  float maxD0Signif = 999; 
  float maxDzSignif = 999; 

  bool rescaleTkPtByError = false; 
  float trackCountThr = 0; 
  
  bool removeTracksInCone = true; 
  float coneSize = 0.05; 
  
  std::vector<float> ptbal_;
  std::vector<float> thrust_;
  std::vector<float> sumpt_;
  std::vector<float> sumpt2_;
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
  
  
  ptbal_.clear(); ptbal_.resize(nVertex,0.);
  thrust_.clear(); thrust_.resize(nVertex,0.);
  sumpt_.clear(); sumpt_.resize(nVertex,0.);
  sumpt2_.clear(); sumpt2_.resize(nVertex,0.);
  sumawy_.clear(); sumawy_.resize(nVertex,0.);
  sumtwd_.clear(); sumtwd_.resize(nVertex,0.);
  sumtrv_.clear(); sumtrv_.resize(nVertex,0.);
  sumweight_.clear(); sumweight_.resize(nVertex,0.);
  ptmax_.clear(); ptmax_.resize(nVertex,0.);
  nchthr_.clear(); nchthr_.resize(nVertex,0.);
  nch_.clear(); nch_.resize(nVertex,0.);
  vtxP_.clear(); vtxP_.resize(nVertex,0.);
  tksPt_.clear(); tksPt_.resize(nVertex, vector<float>(1));
  sphers_.clear(); sphers_.resize(nVertex,TMatrixDSym(3));
  sumpr_.clear(); sumpr_.resize(nVertex,0.);
  spher_.clear(); spher_.resize(nVertex,0.);
  tspher_.clear(); tspher_.resize(nVertex,0.);
  aplan_.clear(); aplan_.resize(nVertex,0.);
  threejetC_.clear(); threejetC_.resize(nVertex,0.);
  fourjetD_.clear(); fourjetD_.resize(nVertex,0.);
	
  diphopt_.clear(); diphopt_.resize(nVertex);
  diPhotonPt_.clear(); diPhotonPt_.resize(nVertex);
  vtxPt_.clear(); vtxPt_.resize(nVertex);
  ptvtx_.clear(); ptvtx_.resize(nVertex);
  diPhotonPz_.clear(); diPhotonPz_.resize(nVertex);
	
  acosA_.clear(); acosA_.resize(nVertex);
  ptasym_.clear(); ptasym_.resize(nVertex);
	
  ptmax3_.clear(); ptmax3_.resize(nVertex);
  thrust_.clear(); thrust_.resize(nVertex);
	
  ptratio_.clear(); ptratio_.resize(nVertex);
  pzasym_.clear(); pzasym_.resize(nVertex);
	
  awytwdasym_.clear(); awytwdasym_.resize(nVertex);
  


  static const float spherPwr_ = 1.5; 
  
  

  std::vector<TLorentzVector> diPhoton_;
  diPhoton_.clear(); diPhoton_.resize(nVertex);
  for(int i=0; i<nVertex; ++i) {
    diPhoton_[i] = photonp4wrtvertex(indpht1,i)
      + photonp4wrtvertex(indpht2,i); 
  }
  
  
  // filling loop over vertexes
  for(int n=0; n<nVertex; n++) {
    
    //foreach vertex loopover all the it's associated  track
    for(int t =0 ; t < vertextrackSize[n]; t++){
      
      float tkWeight = vertex_trkWeight->at(n)[t];
      int indtrk = vertex_trkind->at(n)[t];

      if( indtrk <0 || indtrk >= nTrack){ // 
	continue; 
      }
      
      if( usehighPurityOnly ){
	int qua = trackqualityFlagTracks[indtrk]; 
	
	int highPurityFlag = 3; 
	bool isHP = false; 
	if( ( qua & 1 << highPurityFlag) > 0){
	  isHP = true; 
	}
	if( !isHP) continue; 
      }
      
      if( trackd0[indtrk]/ trackd0Error[indtrk] > maxD0Signif ) continue; 
      if( trackdz[indtrk]/ trackdzError[indtrk] > maxDzSignif ) continue; 
      
      float tkPt = trackpt[indtrk];
      float modpt = tkPt > trackptError[indtrk] ? tkPt - trackptError[indtrk]  : 0.;
      
      TVector3 tkPVec(trackpx[indtrk],trackpy[indtrk],trackpz[indtrk]);
      
      TVector2 tkPtVec(trackpx[indtrk],trackpy[indtrk]);
            
      
      if(rescaleTkPtByError){
	if( modpt ==0) continue; 
	float ptcorr = modpt / tkPt; 
	tkPtVec *= ptcorr;
	tkPt = modpt;
      }
      sumpt2_[n] += tkPtVec.Mod2();
      sumpt_[n] += tkPt;
      if(tkPt > trackCountThr) nchthr_[n] += 1;
      nch_[n] += 1;
      
      
      if ( removeTracksInCone){
	float dr1 = tkPVec.DeltaR( photonp4wrtvertex(indpht1,n).Vect());
	float dr2 = tkPVec.DeltaR( photonp4wrtvertex(indpht2,n).Vect());
	if( dr1 < coneSize || dr2 < coneSize){
	  continue; 
	}
      }
      ptbal_[n] -= tkPtVec * diPhoton_[n].Vect().XYvector().Unit();
      float cosTk = tkPVec.Unit() * diPhoton_[n].Vect().Unit();
      float val = tkPtVec.Mod();
      if ( cosTk < -0.5 )	{
	sumawy_[n] += val;
      } else if ( cosTk > 0.5 ){
	sumtwd_[n] += val;
      } else {
	sumtrv_[n] += val;
      }
      sumweight_[n] += tkWeight;
      vtxP_[n] += tkPVec;
      tksPt_[n].push_back(tkPt);
      
      Float_t p[3] = {0.,0.,0.};
      tkPVec.GetXYZ(p);
      for(int j=3; j--;){
	for(int k=j+1; k--;){
	  (sphers_[n])[j][k] += pow(tkPVec.Mag(),spherPwr_-2.) * p[j]*p[k];
	}
      }
      sumpr_[n] += pow(tkPVec.Mag(),spherPwr_);
      
      
      //if( ( usehighPurityOnly && !e.tkIsHighPurity(tid)  )
      ///|| fabs(e.tkd0(tid,n)/e.tkd0Err(tid,n)) > params_.maxD0Signif 
	    //|| fabs(e.tkdz(tid,n)/e.tkdzErr(tid,n)) > params_.maxDzSignif ) {
	//// std::cerr << "skipping [line" << __LINE__ << "] ";
	//// std::cerr << (params_.highPurityOnly && !e.tkIsHighPurity(tid)  )
	//// 	  << " " 
	//// 	  << (fabs(e.tkd0(tid,n)/e.tkd0Err(tid,n)) > params_.maxD0Signif )
	//// 	  << " " 
	//// 	  << (fabs(e.tkdz(tid,n)/e.tkdzErr(tid,n)) > params_.maxDzSignif) << std::endl;
      
      
    }
    
    sphers_[n] *= 1./sumpr_[n];
    
    TVectorD eigVals(3);
    eigVals = TMatrixDSymEigen(sphers_[n]).GetEigenValues();
    
    spher_[n] = 1.5 * (eigVals[1]+eigVals[2]);
    tspher_[n] = 2. * eigVals[1] / (eigVals[0]+eigVals[1]);
    aplan_[n] = 1.5 * eigVals[2];
    
    threejetC_[n] = 3. * (eigVals[0]*eigVals[1] + eigVals[0]*eigVals[2] + eigVals[1]*eigVals[2]);
    fourjetD_[n] = 27. * eigVals[0]*eigVals[1]*eigVals[2];
    
    
    diPhotonPt_[n] = diPhoton_[n].Vect().XYvector();
    diphopt_[n]    = diPhotonPt_[n].Mod();
    vtxPt_[n]      = vtxP_[n].XYvector();
    ptvtx_[n]      = vtxPt_[n].Mod();
    diPhotonPz_[n]   = diPhoton_[n].Vect().Pz();
    
    sort(tksPt_[n].begin(), tksPt_[n].end(), greater<float>());
    
    acosA_[n] =  	acos(vtxPt_[n].Unit() * diPhotonPt_[n].Unit());
    ptasym_[n] = 	(vtxPt_[n].Mod() - diPhotonPt_[n].Mod())/(vtxPt_[n].Mod() + diPhotonPt_[n].Mod());
    
    ptmax_ [n] = 	tksPt_[n][0];
    ptmax3_[n] = 	accumulate(tksPt_[n].begin(),tksPt_[n].begin() + min(tksPt_[n].size(),(size_t)3), 0.0) ;
    thrust_[n] = 	ptbal_[n]/sumpt_[n];
    
    ptratio_[n] = 	vtxPt_[n].Mod()/diPhotonPt_[n].Mod();
    pzasym_[n] = 	fabs( (vtxP_[n].Pz() - diPhotonPz_[n])/(vtxP_[n].Pz() + diPhotonPz_[n]) );
    
    awytwdasym_[n] = (sumawy_[n]-sumtwd_[n])/(sumawy_[n]+sumtwd_[n]);
    
    
  } //end of loop over each rec vertex
    
  
  bool useAllConvs= true; 
  
  // preselect vertices : all vertices
  std::vector<int> preselAll;
  for(int i=0; i< nVertex; i++) {
    preselAll.push_back(i); 
  }
  
  float zconv = 0; 
  float dzconv = 0;
  std::vector<int> preselConv;
  
  //// VertexAnalysis/src/PhotonInfo.cc:  if (  nTracks_ == 2  &&  convVtxValid_ &&   convVtxChi2Prob_ > 0.0005 )  isAConversion=true;
  //bool  PhotonInfo::isAConversion() {
  // bool isAConversion=false;
  // if (  nTracks_ == 2  &&  convVtxValid_ &&   convVtxChi2Prob_ > 0.0005 )  isAConversion=true;
  // return isAConversion;

  //}

  int iConv1 = useAllConvs ? matchPhotonToConversion(indpht1) : -1;
  int iConv2 = useAllConvs ? matchPhotonToConversion(indpht2) : -1;
  
  //saved iConv1 and iConv2 as well
  photonmatchToallConv[indpht1] = iConv1; 
  photonmatchToallConv[indpht2] = iConv2; 
  
  
  bool ispht1AConversion; 
  if( iConv1 >=0 ){
    ispht1AConversion = convnTracks[iConv1] == 2 && convconversionVertexisValid[iConv1] && convconversionVertexChiSquaredProbability[iConv1] > 0.0005; 
  }else{
    ispht1AConversion = photonconversionnTracks[indpht1] ==2 && photonconversionVertexisValid[indpht1] ==1  && photonconversionChiSquaredProbability[indpht1] > 0.0005; 
  }
  bool ispht2AConversion; 
  
  if( iConv2 >=0 ){
    ispht2AConversion = convnTracks[iConv2] == 2 && convconversionVertexisValid[iConv2] && convconversionVertexChiSquaredProbability[iConv2] > 0.0005; 
  } else{
    ispht2AConversion= photonconversionnTracks[indpht2] ==2 && photonconversionVertexisValid[indpht2] ==1  && photonconversionChiSquaredProbability[indpht2] > 0.0005; 
  }
  if( ispht1AConversion || ispht2AConversion){
    
    if( ispht1AConversion && ! ispht2AConversion){
      //zconv  = vtxAnaFromConv.vtxZ(pho1);
      //dzconv = vtxAnaFromConv.vtxdZ(pho1);
      
      zconv = HggVertexFromConversionsvtxZ(indpht1);
      dzconv = HggVertexFromConversionsvtxdZ(indpht1);
    }
    if( ! ispht1AConversion && ispht2AConversion){
      zconv = HggVertexFromConversionsvtxZ(indpht2);
      dzconv = HggVertexFromConversionsvtxdZ(indpht2);
    }
    if(  ispht1AConversion && ispht2AConversion){
      float z1 = HggVertexFromConversionsvtxZ(indpht1);
      float dz1 = HggVertexFromConversionsvtxdZ(indpht1);
      float z2 = HggVertexFromConversionsvtxZ(indpht2);
      float dz2 = HggVertexFromConversionsvtxdZ(indpht2);
      
      zconv  = sqrt ( 1./(1./dz1/dz1 + 1./dz2/dz2 )*(z1/dz1/dz1 + z2/dz2/dz2) ) ;  // weighted average
      dzconv = sqrt( 1./(1./dz1/dz1 + 1./dz2/dz2)) ;
      
    }
    // preselect vertices : only vertices in a window zconv +/- dzconv
    for(int i=0; i < nVertex ; i++) {
      //TVector3 * vtxpos= (TVector3 *) vtx_std_xyz->At(i);
      //if ( fabs(zconv - vtxpos->Z() ) < dzconv ) 
      if( fabs( zconv - vertexz[i] ) < dzconv){
	preselConv.push_back(i); 
      }
    }
    
  }
  
  ///now preslected vertex 
  vector<int> preselection; 
  if ( preselConv.size()==0 ) 
    preselection = preselAll; 
  else 
    preselection = preselConv; 
  
  
  //std::vector<int> rankprod = vtxAna.rankprod(vtxVarNames);
  //return rankprod;
  
  
  
  
}









float RecoAnalyzer::recHitE( const  DetId id,  const EcalRecHitCollection &recHits )
{
  if ( id == DetId(0) ) {
    return 0;
  } else {
    EcalRecHitCollection::const_iterator it = recHits.find( id );
    if ( it != recHits.end() ) return (*it).energy();
  }
  return 0;
}


float RecoAnalyzer::recHitE( const DetId id, const EcalRecHitCollection & recHits,int di, int dj )
{
  // in the barrel:   di = dEta   dj = dPhi
  // in the endcap:   di = dX     dj = dY
  
  DetId nid;
  if( id.subdetId() == EcalBarrel) nid = EBDetId::offsetBy( id, di, dj );
  else if( id.subdetId() == EcalEndcap) nid = EEDetId::offsetBy( id, di, dj );
  
  return ( nid == DetId(0) ? 0 : recHitE( nid, recHits ) );
}




// ------------ method called to produce the data  ------------
bool RecoAnalyzer::filter(const edm::Event& iEvent, const edm::EventSetup& iSetup)
{
  bool accept = false; 
  
  if (iEvent.isRealData()) {
    // for real data, access the "physics enabled" bit in the L1 GT data
    edm::Handle<L1GlobalTriggerReadoutRecord> h_gtDigis;
    if (not iEvent.getByLabel(m_l1GtRecordInputTag, h_gtDigis)) {
      ////edm::LogWarning(h_gtDigis.whyFailed()->category()) << h_gtDigis.whyFailed()->what();
      // not enough informations to make a decision - reject the event
      return false;
    } else {
      L1GtFdlWord fdlWord = h_gtDigis->gtFdlWord();
      if (fdlWord.physicsDeclared() == 1) 
        accept = true;
    }
  } else {
    // for MC, assume the "physics enabled" bit to be always set
    accept = true;
  }
  
  return accept;
  
  
}







float RecoAnalyzer::recHitApproxEt(  const DetId id,  const EcalRecHitCollection &recHits )
{
  // for the time being works only for the barrel
  if ( id.subdetId() == EcalBarrel ) {
    return recHitE( id, recHits ) / cosh( EBDetId::approxEta( id ) );
  }
  return 0;
}

double RecoAnalyzer::mye2overe9( const DetId id, const EcalRecHitCollection & recHits)
{
  ///////////start calculating e2/e9                                                                                                                                   
  ////http://cmslxr.fnal.gov/lxr/source/RecoLocalCalo/EcalRecAlgos/src/EcalSeverityLevelAlgo.cc#240                                                                    
  // compute e2overe9                                                                                                                                                  
  //                                                                                                                                                                   
  //   | | | |                                                                                                                                                         
  //   +-+-+-+                                                                                                                                                         
  //   | |1|2|                                                                                                                                                         
  //   +-+-+-+                                                                                                                                                         
  //   | | | |                                                                                                                                                         
  //                                                                                                                                                                   
  //   1 - input hit,  2 - highest energy hit in a 3x3 around 1  
  //                                                                                                                                                                     
  //   rechit 1 must have E_t > recHitEtThreshold                                                                                                                        
  //   rechit 2 must have E_t > recHitEtThreshold2                                                                                                                       
  //                                                                                                                                                                     
  //   function returns value of E2/E9 centered around 1 (E2=energy of hits 1+2) if energy of 1>2                                                                        
  //                                                                                                                                                                     
  //   if energy of 2>1 and KillSecondHit is set to true, function returns value of E2/E9 centered around 2                                                              
  //   *provided* that 1 is the highest energy hit in a 3x3 centered around 2, otherwise, function returns 0                                                             

  
  float recHitEtThreshold = 10.0;
  float recHitEtThreshold2 = 1.0;
  bool avoidIeta85=false;
  bool KillSecondHit=true;
  

  if ( id.subdetId() == EcalBarrel ) {

    EBDetId ebId( id );
    // avoid recHits at |eta|=85 where one side of the neighbours is missing                                                                                             
    if ( abs(ebId.ieta())==85 && avoidIeta85) return 0;

    // select recHits with Et above recHitEtThreshold                                                                                                                    
    float e1 = recHitE( id, recHits );
    float ete1=recHitApproxEt( id, recHits );


    // check that rechit E_t is above threshold                                                                                                                          
    if (ete1 < std::min(recHitEtThreshold,recHitEtThreshold2) ) return 0;

    if (ete1 < recHitEtThreshold && !KillSecondHit ) return 0;

    float e2=-1;
    float ete2=0;
    float s9 = 0;

    // coordinates of 2nd hit relative to central hit 
    int e2eta=0;
    int e2phi=0;

    // LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 1                                                                                                                         

    for ( int deta = -1; deta <= +1; ++deta ) {
      for ( int dphi = -1; dphi <= +1; ++dphi ) {

	// compute 3x3 energy                                                                                                                                            
	float etmp=recHitE( id, recHits, deta, dphi );
	s9 += etmp;

	EBDetId idtmp=EBDetId::offsetBy(id,deta,dphi);
	float eapproxet=recHitApproxEt( idtmp, recHits );
	
	// remember 2nd highest energy deposit (above threshold) in 3x3 array                                                                                            
	if (etmp>e2 && eapproxet>recHitEtThreshold2 && !(deta==0 && dphi==0)) {
	  e2=etmp;
	  ete2=eapproxet;
	  e2eta=deta;
	  e2phi=dphi;

	}

      }
    }

    if ( e1 == 0 )  return 0;

    // return 0 if 2nd hit is below threshold                                                                                                                            
    if ( e2 == -1 ) return 0;

    // compute e2/e9 centered around 1st hit                                                                                                                             
    float e2nd=e1+e2;
    float e2e9=0;
    if (s9!=0) e2e9=e2nd/s9;

    // if central hit has higher energy than 2nd hit                                                                                                                     
    //  return e2/e9 if 1st hit is above E_t threshold                                                                                                                   

    if (e1 > e2 && ete1>recHitEtThreshold) return e2e9;

    // if second hit has higher energy than 1st hit                                                                                                                      
    if ( e2 > e1 ) {


      // return 0 if user does not want to flag 2nd hit, or                                                                                                              
      // hits are below E_t thresholds - note here we                                                                                                                    
      // now assume the 2nd hit to be the leading hit.                                                                                                                   

      if (!KillSecondHit || ete2<recHitEtThreshold || ete1<recHitEtThreshold2) {
	return 0;

      }
            
      else {

	// LOOP OVER 3x3 ARRAY CENTERED AROUND HIT 2                                                                                                                     
	float s92nd=0;

	float e2nd_prime=0;
	int e2prime_eta=0;
	int e2prime_phi=0;

	EBDetId secondid=EBDetId::offsetBy(id,e2eta,e2phi);

	for ( int deta = -1; deta <= +1; ++deta ) {
	  for ( int dphi = -1; dphi <= +1; ++dphi ) {
	    // compute 3x3 energy                                                                                                                                        
	    float etmp=recHitE( secondid, recHits, deta, dphi );
	    s92nd += etmp;

	    if (etmp>e2nd_prime && !(deta==0 && dphi==0)) {
	      e2nd_prime=etmp;
	      e2prime_eta=deta;
	      e2prime_phi=dphi;
	    }

	  }
	}

	// if highest energy hit around E2 is not the same as the input hit, return 0;                                                                                   
	if (!(e2prime_eta==-e2eta && e2prime_phi==-e2phi))
	  {
	    return 0;
	  }


	// compute E2/E9 around second hit                                                                                                                               
	float e2e9_2=0;
	if (s92nd!=0) e2e9_2=e2nd/s92nd;

	//   return the value of E2/E9 calculated around 2nd hit                                                                                                         
	return e2e9_2;


      }

    }

  } else if ( id.subdetId() == EcalEndcap ) {
    // only used for EB at the moment                                                                                                                                    
    return 0;
  }
  return 0;
}//double Analyser::mye2overe9( const DetId id, const EcalRecHitCollection & recHits)                                                                                  



// ------------ method called once each job just after ending the event loop  ------------
void 
RecoAnalyzer::endJob() 
{

  rootFile_->cd();
  Analysis->Write();
  
  rootFile_->Write() ;
  rootFile_->Close() ;

}

// ------------ method called when starting to processes a run  ------------
void 
RecoAnalyzer::beginRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a run  ------------
void 
RecoAnalyzer::endRun(edm::Run const&, edm::EventSetup const&)
{
}

// ------------ method called when starting to processes a luminosity block  ------------
void 
RecoAnalyzer::beginLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method called when ending the processing of a luminosity block  ------------
void 
RecoAnalyzer::endLuminosityBlock(edm::LuminosityBlock const&, edm::EventSetup const&)
{
}

// ------------ method fills 'descriptions' with the allowed parameters for the module  ------------
void
RecoAnalyzer::fillDescriptions(edm::ConfigurationDescriptions& descriptions) {
  //The following says we do not know what parameters are allowed so do no validation
  // Please change this to state exactly what you do use, even if it is no parameters
  edm::ParameterSetDescription desc;
  desc.setUnknown();
  descriptions.addDefault(desc);
}

//define this as a plug-in
DEFINE_FWK_MODULE(RecoAnalyzer);
