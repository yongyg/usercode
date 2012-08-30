#include "rootheader.h"
#include "testSelection_cluster.h"

#include "common_functions.cc"


float interCalib_preCalib[170][360];
float interCalibEndcap_preCalib[2][101][101];


bool isTestBeamSM(int iSM){
    


  //if( iSM == 1 || iSM == 2  || iSM == 11 || iSM == 15  || iSM == 21 || iSM == 24) return true; 
  
  if( iSM == 1  || iSM == 11 || iSM == 15  || iSM == 24) return true; 
  
  
  else return false; 
    
}

#include "effSigma.C"

#include "gausfit.cc"



#include "foldEndcap.cc"


#include "getCalibConstants.cc"




void checkHowManyDeadCrystalNext(int iz, int j, int k, int flagmap_ietaiphi[2][101][101],int &ncorner, int &nside, int &nflagdeadNextTo){

  
  ///input i, j after convxvalid 
  
  if( ! (iz >=0 && iz<=1) ){
    cout<<"iz 0/1 ? " << iz <<endl; 
    return; 
  }
  
  
  
  int nbad = 0; 
  nside = 0;
  ncorner = 0;




  for(int nieta = -1; nieta <=1; nieta++){
    int jj= j+ nieta; 
    if( jj <=0 || jj >=101) continue; 
    
    for(int niphi = -1; niphi <=1; niphi++){
      
      if(nieta==0 && niphi==0) continue; //itself
      int kk = k + niphi; 
      
      if( validRecHitEndCap[iz][jj][kk] ==0) continue; 
            
      if(flagmap_ietaiphi[iz][jj][kk]>0 ){
	nbad ++;
	int deta = diff_neta(j,jj);
	int dphi = diff_nphi(k,kk);
	
	if( deta ==1 && dphi ==1) ncorner ++;
	else nside ++;
	
      }
    }
    
  }
  
  nflagdeadNextTo = 0; 
   
  if(nbad==0){
    nflagdeadNextTo = 0; 
  }
  else if(nbad==1){
    if(ncorner ==1) {
      nflagdeadNextTo = 1; 
    }else{
      nflagdeadNextTo = 2; 
    }
  }else if(nbad==2){
    if(ncorner ==2) {
      nflagdeadNextTo = 3; 
    }else if( nside ==2) {  //
      nflagdeadNextTo = 4; 
    }else { //NZ
      nflagdeadNextTo = 5; 
    }
  }else if(nbad==3){
    if(ncorner ==3) {
      nflagdeadNextTo = 6; 
    }else if( nside ==3) { ///NZ
      nflagdeadNextTo = 7; 
    }else if( ncorner ==1 && nside ==2) {
      nflagdeadNextTo = 8; 
    }else if( ncorner ==2 && nside ==1) {  //NZ
      nflagdeadNextTo = 9; 
    }
  }else if( nbad==4){
    if(ncorner ==4 ) {
      nflagdeadNextTo = 10; 
    }else if( nside ==4) {
      nflagdeadNextTo = 11; 
    }else if( ncorner ==2 && nside ==2) {  //NZ
      nflagdeadNextTo = 12; 
    }else if( ncorner ==1 && nside ==3) {
      nflagdeadNextTo = 13; 
    }else if( ncorner ==3 && nside ==1) {  //NZ
      nflagdeadNextTo = 14;
    }
  } else if( nbad ==5){
    if( ncorner ==4 ) {
      nflagdeadNextTo = 15;
    }
    else if( nside ==4 ) {
      nflagdeadNextTo = 16;
    }
    else if( nside ==3 && ncorner ==2 ) {
      nflagdeadNextTo = 17;
    }
    else if( nside ==2 && ncorner ==3 ) {
      nflagdeadNextTo = 18;
    }
    
    else{
      cout<<"? " << nbad <<" "<< nside <<" "<< ncorner <<endl; 
    }
  }
  
  else if(nbad == 6){
    if( ncorner ==4 && nside ==2  ) {
      nflagdeadNextTo = 19;
    }
    
    else{
      cout<<"? " << nbad <<" "<< nside <<" "<< ncorner <<endl; 
    }
    
  }

  else{
    cout<<"nbad? "<< nbad <<" "<< iz <<" "<< j<<" "<<k<<" "<< ncorner <<" "<< nside <<endl; 
    exit(1); 
  }
  
}

///input the IC file in calibres/
void make_deadflagEE(char *inputfile, int dataflag ){
  
  for(int j=0; j<2; j++){
    for(int x =0; x<101; x++){
      for(int y=0; y<101; y++){
	validRecHitEndCap[j][x][y] = 0;
      }
    }
  }
  get_xyzEBrechits();
  setEtaRingBoundaryEndcap();
  
  

  readInterCalibEndcap_GR09_V8();
  
  ifstream txtin(inputfile,ios::in);
  
  string sinput = inputfile; 

  int pizEta = 1;
  if( sinput.find("pe2") != string::npos){
    pizEta = 2; 
  }


  cout<<"pizEta " << pizEta <<endl; 
  int nMaxRingIC = 20; 
  if( pizEta==2){
    nMaxRingIC = 30; 
  }
  
  
  int iz; 
  int ix; 
  int iy; 
  int flag; 
  float cc; 
  int flagEE[2][101][101];
  while (txtin.good()){
    txtin>> iz >> ix >> iy >> cc; 
    
    flag = cc <0 ;   //-1 means dead xtal 
    
    int iring = iRingEndCap(2*iz-1,ix,iy); ///input -1/1,  
    
    ///if(iring > nMaxRingIC ) flag = 0; //always good for those rings...
    
    flagEE[iz][ix][iy]= flag; 
  }
  
  
  TString output = TString(Form("crystal_deadflag_ee_dflag%d.txt",dataflag));
  cout<<output<<endl; 
  ofstream txtout(output,ios::out);
  
  for(int iz=0; iz<2; iz++){
    for(int j=0; j<101; j++){
      for(int k=0; k< 101; k++){
        if( validRecHitEndCap[iz][j][k] <1) continue;
	
	int ncorner;
	int nside; 
	int nflagdead; 

	if( flagEE[iz][j][k] ==0){
	  checkHowManyDeadCrystalNext(iz,j,k,flagEE,ncorner,nside,nflagdead);
	  txtout<<2*iz-1 <<" "<< j<<" "<<k<<" "<< flagEE[iz][j][k] <<" "<< ncorner <<" "<<nside <<" "<< nflagdead <<endl;
	}else{
	  txtout<<2*iz-1 <<" "<< j<<" "<<k<<" "<< flagEE[iz][j][k]<<" "<< -1<<" "<< -1 <<" "<<-1 <<endl;
	}
	
      }
    }
  }
  
  
  
}
