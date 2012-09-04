#include "RecoLocalCalo/EcalRecProducers/plugins/EcalRecHitWorkerSimple.h"
#include "RecoLocalCalo/EcalRecAlgos/interface/EcalSeverityLevelAlgo.h"
#include "FWCore/Framework/interface/EventSetup.h"
#include "FWCore/Framework/interface/Event.h"
#include "CondFormats/DataRecord/interface/EcalIntercalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeCalibConstantsRcd.h"
#include "CondFormats/DataRecord/interface/EcalADCToGeVConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalTimeOffsetConstantRcd.h"
#include "CondFormats/DataRecord/interface/EcalChannelStatusRcd.h"
#include "CondFormats/EcalObjects/interface/EcalTimeCalibConstants.h"
#include "CalibCalorimetry/EcalLaserCorrection/interface/EcalLaserDbRecord.h"
#include "FWCore/MessageLogger/interface/MessageLogger.h"
#include "CommonTools/Utils/interface/StringToEnumValue.h"

#include <iomanip>
#include <fstream>
#include <iostream>
using namespace std;

EcalRecHitWorkerSimple::EcalRecHitWorkerSimple(const edm::ParameterSet&ps) :
        EcalRecHitWorkerBaseClass(ps)
{
        rechitMaker_ = new EcalRecHitSimpleAlgo();
        v_chstatus_ = ps.getParameter<std::vector<int> >("ChannelStatusToBeExcluded");
	v_DB_reco_flags_ = ps.getParameter<std::vector<int> >("flagsMapDBReco");
        killDeadChannels_ = ps.getParameter<bool>("killDeadChannels");
        laserCorrection_ = ps.getParameter<bool>("laserCorrection");
	EBLaserMIN_ = ps.getParameter<double>("EBLaserMIN");
	EELaserMIN_ = ps.getParameter<double>("EELaserMIN");
	EBLaserMAX_ = ps.getParameter<double>("EBLaserMAX");
	EELaserMAX_ = ps.getParameter<double>("EELaserMAX");


		testNewInterCalibEB_ = ps.getUntrackedParameter<bool>("testNewInterCalibEB",false);
	testNewInterCalibEE_ = ps.getUntrackedParameter<bool>("testNewInterCalibEE",false);
	
	if( testNewInterCalibEE_){
	  pathInterCalibFileEE_ = ps.getUntrackedParameter<std::string>("pathInterCalibFileEE");
	  std::cout<<pathInterCalibFileEE_.c_str()<<std::endl;
          ifstream inputcc(pathInterCalibFileEE_.c_str(),ios::in);
          bool foundfile = true;
          if (inputcc.fail()){
	    std::cout<<"error open file.. " << pathInterCalibFileEE_.c_str()<<std::endl;
            foundfile = false;
          }
	  if( foundfile){
	    int x; 
	    int y; 
	    int iz; 
	    float tmp;
	    double mean = 0;
	    int nread = 0; 
	    while( inputcc.good()){
              inputcc >> iz >>x >>y >> tmp;

	      int izz = iz < 0 ? 0:1;
	      if(tmp<0) tmp=1; 
	      
	      interCalibEndcap[izz][x][y]  = tmp;
	      mean += tmp; 
	      nread ++; 
	      if( nread >= 14648){
		break; 
	      }
	    }
	    mean /= nread; 
	    std::cout<<"mean " << mean <<std::endl;
	    for(int izz = 0; izz < 2; izz++){
	      for(int x = 0; x< 101; x++){
		for(int y = 0; y< 101; y++){
		  interCalibEndcap[izz][x][y] /= mean;
		}
	      }
	    }
	  }
	}
	
	
	if( testNewInterCalibEB_ ){
	  pathInterCalibFile_ = ps.getUntrackedParameter<std::string>("pathInterCalibFile");
	  ////file under lib/slc5...//
	  std::cout<<pathInterCalibFile_.c_str()<<std::endl;
	  ifstream inputcc(pathInterCalibFile_.c_str(),ios::in);
	  bool foundfile = true;
	  if (inputcc.fail()){
	    std::cout<<"error open file.. " << pathInterCalibFile_.c_str()<<std::endl; 
	    foundfile = false; 
	  }
	  
	  if( foundfile){
	    int ieta; 
	    int iphi; 
	    float tmp; 
	    double mean = 0; 
	    int nread = 0; 
	    while( inputcc.good()){
	      inputcc >> ieta >> iphi >> tmp; 
	      if(tmp<0) tmp=1; 
	      convxtalid(iphi,ieta);
	      interCalibBarrel[ieta+85][iphi] = tmp;
	      
	      mean += tmp; 
	      nread ++; 
	      if(nread >= 61200){
		break; 
	      }
	    }
	    mean = mean/nread;
	    std::cout<<"mean " << mean <<std::endl; 
	    for(int j=0; j<170; j++){
	      for(int k=0; k<360; k++){
		interCalibBarrel[j][k] /= mean;
	      }
	    }
	    
	  }
	}

}



void EcalRecHitWorkerSimple::set(const edm::EventSetup& es)
{
        es.get<EcalIntercalibConstantsRcd>().get(ical);
        es.get<EcalTimeCalibConstantsRcd>().get(itime);
        es.get<EcalTimeOffsetConstantRcd>().get(offtime);
        es.get<EcalADCToGeVConstantRcd>().get(agc);
        es.get<EcalChannelStatusRcd>().get(chStatus);
        if ( laserCorrection_ ) es.get<EcalLaserDbRecord>().get(laser);
}


bool
EcalRecHitWorkerSimple::run( const edm::Event & evt,
                const EcalUncalibratedRecHit& uncalibRH,
                EcalRecHitCollection & result )
{
        DetId detid=uncalibRH.id();

        EcalChannelStatusMap::const_iterator chit = chStatus->find(detid);
        EcalChannelStatusCode chStatusCode = 1;
        if ( chit != chStatus->end() ) {
                chStatusCode = *chit;
        } else {
                edm::LogError("EcalRecHitError") << "No channel status found for xtal " 
                        << detid.rawId() 
                        << "! something wrong with EcalChannelStatus in your DB? ";
        }
        if ( v_chstatus_.size() > 0) {
                uint16_t code = chStatusCode.getStatusCode() & 0x001F;
                std::vector<int>::const_iterator res = std::find( v_chstatus_.begin(), v_chstatus_.end(), code );
                if ( res != v_chstatus_.end() ) {
                        return false;
                }
        }

        // find the proper flag for the recHit
        // from a configurable vector
        // (see cfg file for the association)
        uint32_t recoFlag = 0;
        uint16_t statusCode = chStatusCode.getStatusCode() & 0x001F;
        if ( statusCode < v_DB_reco_flags_.size() ) {
                // not very nice...
                recoFlag = v_DB_reco_flags_[ statusCode ];  
        } else {
                edm::LogError("EcalRecHitError") << "Flag " << statusCode 
                        << " in DB exceed the allowed range of " << v_DB_reco_flags_.size();
        }

	float offsetTime = 0; // the global time phase
	const EcalIntercalibConstantMap& icalMap = ical->getMap();  
        if ( detid.subdetId() == EcalEndcap ) {
                rechitMaker_->setADCToGeVConstant( float(agc->getEEValue()) );
		offsetTime = offtime->getEEValue();
        } else {
                rechitMaker_->setADCToGeVConstant( float(agc->getEBValue()) );
		offsetTime = offtime->getEBValue();
        }

        // first intercalibration constants
        EcalIntercalibConstantMap::const_iterator icalit = icalMap.find(detid);
        EcalIntercalibConstant icalconst = 1;
        if( icalit!=icalMap.end() ) {
                icalconst = (*icalit);
        } else {
                edm::LogError("EcalRecHitError") << "No intercalib const found for xtal "
                        << detid.rawId()
                        << "! something wrong with EcalIntercalibConstants in your DB? ";
        }

        // get laser coefficient
        float lasercalib = 1.;
        if ( laserCorrection_ )	lasercalib = laser->getLaserCorrection( detid, evt.time());
	

        // get time calibration coefficient
        const EcalTimeCalibConstantMap & itimeMap = itime->getMap();  
        EcalTimeCalibConstantMap::const_iterator itime = itimeMap.find(detid);
        EcalTimeCalibConstant itimeconst = 0;
        if( itime!=itimeMap.end() ) {
                itimeconst = (*itime);
		  } else {
                edm::LogError("EcalRecHitError") << "No time calib const found for xtal "
                        << detid.rawId()
                        << "! something wrong with EcalTimeCalibConstants in your DB? ";
        }
          
	 
	
	if ( testNewInterCalibEB_ && detid.subdetId() == EcalBarrel ) {
	  EBDetId ebd = EBDetId(detid);
	  int ieta = ebd.ieta();
	  int iphi = ebd.iphi();
	  convxtalid(iphi,ieta);
	  icalconst = interCalibBarrel[ieta+85][iphi];
	}
	if( testNewInterCalibEE_ && detid.subdetId() == EcalEndcap ) {
	  EEDetId ebd = EEDetId(detid);
	  int ieta = ebd.ix();
	  int iphi = ebd.iy();
	  int izz = ebd.zside() <0 ? 0:1;
	  icalconst = interCalibEndcap[izz][ieta][iphi];
	}
	

        // make the rechit and put in the output collection
	if (recoFlag<=EcalRecHit::kLeadingEdgeRecovered || !killDeadChannels_) {
          EcalRecHit myrechit( rechitMaker_->makeRecHit(uncalibRH, icalconst * lasercalib, (itimeconst + offsetTime), /*recoflags_*/ 0) );	
	  if (detid.subdetId() == EcalBarrel && (lasercalib < EBLaserMIN_ || lasercalib > EBLaserMAX_)) myrechit.setFlag(EcalRecHit::kPoorCalib);
	  if (detid.subdetId() == EcalEndcap && (lasercalib < EELaserMIN_ || lasercalib > EELaserMAX_)) myrechit.setFlag(EcalRecHit::kPoorCalib);
	  result.push_back(myrechit);
	}

        return true;
}


void EcalRecHitWorkerSimple::convxtalid(Int_t &nphi,Int_t &neta)
{
  // Barrel only
  // Output nphi 0...359; neta 0...84; nside=+1 (for eta>0), or 0 (for eta<0).
  // neta will be [-85,-1] , or [0,84], the minus sign indicates the z<0 side.
  
  if(neta > 0) neta -= 1;
  if(nphi > 359) nphi=nphi-360;
  
  // final check
  if(nphi >359 || nphi <0 || neta< -85 || neta > 84)
    {
      std::cout <<" unexpected fatal error in HLTEcalResonanceFilter::convxtalid "<<  nphi <<  " " << neta <<  " " <<std::endl;
    }
} //end of convxtalid





EcalRecHitWorkerSimple::~EcalRecHitWorkerSimple(){

  delete rechitMaker_;
}


#include "FWCore/Framework/interface/MakerMacros.h"
#include "RecoLocalCalo/EcalRecProducers/interface/EcalRecHitWorkerFactory.h"
DEFINE_EDM_PLUGIN( EcalRecHitWorkerFactory, EcalRecHitWorkerSimple, "EcalRecHitWorkerSimple" );
