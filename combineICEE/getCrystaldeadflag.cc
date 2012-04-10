void getCrystaldeadflagEndcap(){
  
  
  char *file_input = new char[500];

  sprintf(file_input,"crystal_deadflag_run2011amay10.txt");



  cout<<"READING fro deadcystalbarrel file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file.. " << file_input<<endl;
    exit(1);
  }
   
  int iz; 
  int ieta; 
  int iphi; 
  int flag; 
  
  int n = 0; 
  int ndead = 0; 

  int ncorner;
  int nside; 

  int ndeadflag; 
  
  while (inputcc.good()){
    
    inputcc >> iz >> ieta >> iphi >> flag >> ncorner >> nside  >> ndeadflag; 
    if( abs(iz) !=1){
      cout<<"worng ?? " << iz <<endl;
      exit(1);
    }

    int izz = iz < 0 ? 0:1;
    flag_endcap[izz][ieta][iphi] = flag; 
    ndeadcorner_endcap[izz][ieta][iphi] = ncorner; 
    ndeadside_endcap[izz][ieta][iphi] = nside;
    ndeadflag_endcap[izz][ieta][iphi] = ndeadflag; 
    
    
    n ++; 
    if( flag ==1) ndead ++; 
    
    if( n>= 14648) break; 
  }
  
  cout<<"deadflag read: "<< n <<" "<< ndead <<endl; 
  

}


void getCrystaldeadflagEndcap_v1(char *file_input,int ndeadflagg[2][101][101]){
  //  char *file_input = new char[500];
  //sprintf(file_input,"crystal_deadflag_run2011amay10.txt");

  
  cout<<"READING fro deadcystalbarrel file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file.. " << file_input<<endl;
    exit(1);
  }
   
  int iz; 
  int ieta; 
  int iphi; 
  int flag; 
  
  int n = 0; 
  int ndead = 0; 

  int ncorner;
  int nside; 

  int ndeadflag; 
  
  while (inputcc.good()){
    
    inputcc >> iz >> ieta >> iphi >> flag >> ncorner >> nside  >> ndeadflag; 
    if( abs(iz) !=1){
      cout<<"worng ?? " << iz <<endl;
      exit(1);
    }

    int izz = iz < 0 ? 0:1;
    flag_endcap[izz][ieta][iphi] = flag; 
    ndeadcorner_endcap[izz][ieta][iphi] = ncorner; 
    ndeadside_endcap[izz][ieta][iphi] = nside;
    ndeadflagg[izz][ieta][iphi] = ndeadflag; 
    
    
    n ++; 
    if( flag ==1) ndead ++; 
    
    if( n>= 14648) break; 
  }
  
  cout<<"deadflag read: "<< n <<" "<< ndead <<endl; 
  

}



void getCrystaldeadflagBarrel_v1(char *file_input , int ndeadflagg[170][360]){
  
   //  char *file_input = new char[500];

//   //sprintf(file_input,"crystal_deadflag.txt");
//   sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcapiz.tightCut.eg2ps2.runto144114.txt");
//   sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcapiz.run2010b.146428to194442.txt");
//   sprintf(file_input,"crystal_deadflag_alcapiz.run2010b.146428to194442.txtnew"); ///two more crystals set to dead by hand.
  
  
  cout<<"READING fro deadcystalbarrel file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  int ieta; 
  int iphi; 
  int flag; 
  
  int n = 0; 
  int ndead = 0; 

  int ncorner;
  int nside; 

  int ndeadflag; 

  //format
  /// ieta, iphi ( after convxtalid, ieta+85 ) , flag ( 0/1) ncorner, nside, ndeadflag
    
  while (inputcc.good()){

    inputcc >> ieta >> iphi >> flag >> ncorner >> nside  >> ndeadflag; 

    flag_ietaiphi[ieta][iphi] = flag; 
    ndeadcorner_ietaiphi[ieta][iphi] = ncorner; 
    ndeadside_ietaiphi[ieta][iphi] = nside;
    
    ndeadflagg[ieta][iphi] = ndeadflag; 
    
    
    n ++; 
    if( flag ==1) ndead ++; 
    
    if( n>= 61200) break; 
  }
  
  cout<<"deadflag read: "<< n <<" "<< ndead <<endl; 
    
  
  
}

void getCrystaldeadflagBarrel(){
    
  char *file_input = new char[500];

  //sprintf(file_input,"crystal_deadflag.txt");
  sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcapiz.tightCut.eg2ps2.runto144114.txt");
  sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcapiz.run2010b.146428to194442.txt");
  sprintf(file_input,"crystal_deadflag_alcapiz.run2010b.146428to194442.txtnew"); ///two more crystals set to dead by hand.
  
  
  ///2011A
  sprintf(file_input,"crystal_deadflag_alcapiz.run2011a.160405to163869.txt");
  

  
  
  sprintf(file_input,"crystal_deadflag_eb_dflag64.txt");
  

  
  if(doPizEta ==1){
    
//     if(dataflag==25){
//       sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcapiz.looseCut.runto139980.txt");
//     }else if( dataflag==24){
//       sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcapiz.tightCut.runto143193.txt");
//     }else if( dataflag==23){
//       sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcapiz.tightCut.eg2ps2.runto144114.txt");
//     }
//     else if( dataflag ==26 || dataflag == 27 || dataflag == 28){
//       sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcapiz.run2010b.146428to194442.txt");
//     }
//     else{
//       cout<<"warning deadcrystalmap .."<< dataflag<<endl; 
//     }
    

  }else{
    
//     if(dataflag==201){
//       sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcaeta.runto143193.txt");
//     }else if(dataflag==202){
//       sprintf(file_input,"/uscms/home/yongy/work/crystal_deadflag_alcaeta.eg2ps2.runto143193.txt");
//     }
        
//     else{
//       cout<<"warning deadcrystalmap .."<< dataflag<<endl;
//     }
        
  }
  
  cout<<"READING fro deadcystalbarrel file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file barrel.. " << file_input<<endl;
    exit(1);
  }
  
  int ieta; 
  int iphi; 
  int flag; 
  
  int n = 0; 
  int ndead = 0; 

  int ncorner;
  int nside; 

  int ndeadflag; 

  //format
  /// ieta, iphi ( after convxtalid, ieta+85 ) , flag ( 0/1) ncorner, nside, ndeadflag
    
  while (inputcc.good()){

    inputcc >> ieta >> iphi >> flag >> ncorner >> nside  >> ndeadflag; 
    flag_ietaiphi[ieta][iphi] = flag; 
    ndeadcorner_ietaiphi[ieta][iphi] = ncorner; 
    ndeadside_ietaiphi[ieta][iphi] = nside;
    
    ndeadflag_ietaiphi[ieta][iphi] = ndeadflag; 
    
    
    n ++; 
    if( flag ==1) ndead ++; 
    
    if( n>= 61200) break; 
  }
  
  cout<<"deadflag read: "<< n <<" "<< ndead <<endl; 
  
  
}



void getCorrFactorDead(){
  
  for(int j=0; j< 20; j++){
    corrfactorDead[j] = 1; 
  }
  
  
}


///No iteration, just normalize the mass peak of crystals next to dead crystals to the peak of crytsals without dead crystals nearby
// and corretion on the cluster energy 

///maybe apply correction to mpair ; mpair *= fcorr1 * fcorr2 


void getCorrFactorDead_Barrel_CorrectionOnClusterEnergy(){
  
  char *filename = new char[500];
  
  if( applyCorrEta == 1 && applyCorrPhi ==1){
    sprintf(filename,"/uscmst1b_scratch/lpc1/3DayLifetime/yangyong/calib_smear0v1/calibres/deriveCalibConst.testCalibv1.dflag%d.pe%d.cut%d.rmOvlap%d.step%d.method%d.corrEta%d.corrPhi%d.corrSM%d.corrDead%d.precalib%d.vtx%d.encorr%d.evtNot%d.trig%d.root", dataflag,doPizEta, cutset,rmOverlap,11,1,1,0,1,-1, usePreCalib,useVtx,applyEnCorr,-1,trigger);
  }else if( applyCorrEta == 21 && applyCorrPhi ==11){
    sprintf(filename,"/uscmst1b_scratch/lpc1/3DayLifetime/yangyong/calib_smear0v1/calibres/deriveCalibConst.testCalibv1.dflag%d.pe%d.cut%d.rmOvlap%d.step%d.method%d.corrEta%d.corrPhi%d.corrSM%d.corrDead%d.precalib%d.vtx%d.encorr%d.evtNot%d.trig%d.root", dataflag,doPizEta, cutset,rmOverlap,11,1,21,10,1,-1, usePreCalib,useVtx,applyEnCorr,-1,trigger);
  }
  
  //calibres/deriveCalibConst.testCalibv1.dflag19.pe1.cut4.rmOvlap0.step10.method1.corrEta1.corrPhi0.corrSM1.corrDead0.precalib1.vtx0.encorr24.evtNot-1.root
  
  cout<<filename<<endl; 
  
  TFile *ff =new TFile(filename,"read"); 

  TH1F *hh_corr_dead = (TH1F*)ff->Get("hh_corr_deadv1");
  
  for(int b=1; b<= hh_corr_dead->GetNbinsX(); b++){
    corrfactorDead[b-1] = hh_corr_dead->GetBinContent(b); 
    if( b==1 || corrfactorDead[b-1]!=1) cout<<"corr_dead"<<b-1<<" "<< corrfactorDead[b-1]<<endl; 
  }
  cout<<"correction of dead crystal read.."<<endl; 
    
  
}
