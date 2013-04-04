#include "rootheader.h"

#include "testSim.h"

bool debug = 0;

TChain *fChain;
TChain *fChain1;
int entry; 

map<int,map<int, vector<double> > > map_planes_eb; /// ix, iy 6planes 0: 0123, 1: 0145, 2: 1256, 3:2367, 4:0347, 5: 4567

map<int,map<int, map<int, vector<double> > > > map_planes_ee; /// iz , ix, iy 6planes 0: 0123, 1: 0145, 2: 1256, 3:2367, 4:0347, 5: 4567

int coindall[6][5] ={
  {0,1,2,3,0},
  {0,1,5,4,0},
  {1,2,6,5,1},
  {2,3,7,6,2},
  {0,3,7,4,0},
  {4,5,6,7,4}
};



#include "setbranchaddress.cc"

#include "utils.cc"

#include "th3f.cc"

#include "planes.cc"

void testSim(){

  ///sims tep 
  fChain = new TChain("G4SIM");
  fChain->Add("g4simhitsEcal.root");
  fChain1 = new TChain("Analysis");

  setbranchaddress_g4step();

  int N = fChain->GetEntries();

  fChain1->Add("analysis.root");
  setbranchaddress_g4sim();
  
  get_xyzEBrechits();

  TFile *fnew = new TFile("testSime.root","recreate");    
  TTree *mytree = new TTree("Analysis","conv xy in front face");

  int isconv; 
  int indplane; 
  float distplane; 
  float split_new_x;
  float split_new_y;
  float prepos[3];
  float pretime;
  
  float pretimev1[10]; //E>1, E>2, E>3, 4., ...10
  
  int pids[3];  //0: primary id, 1: particle id, 2: parent id
  mytree->Branch("split_new_x",&split_new_x,"split_new_x/F");
  mytree->Branch("split_new_y",&split_new_y,"split_new_y/F");
  mytree->Branch("indplane",&indplane,"indplane/I");
  mytree->Branch("distplane",&distplane,"distplane/F");
  mytree->Branch("prepos",prepos,"prepos[3]/F");
  mytree->Branch("pretime",&pretime,"pretime/F");
  mytree->Branch("pretimev1",pretimev1,"pretimev1[10]/F");
  mytree->Branch("isconv",&isconv,"isconv/I");
  mytree->Branch("pids",pids,"pids[3]/I");
  float ine; 
  mytree->Branch("ine",&ine,"ine/F");
  

  //one in eb
  makeTH1F("pretmin_0",30000,3,6);
  makeTH1F("pretmin_unconv_0",30000,3,6);
  makeTH1F("pretmin_conv_0",30000,3,6);
  
  //one xtal in ee
  makeTH1F("pretmin_1",30000,8,14);
  makeTH1F("pretmin_unconv_1",30000,8,14);
  makeTH1F("pretmin_conv_1",30000,8,14);
  
  double speedoflight = 29.9792458;  //cm /ns
  double p[4];
  bool interset;

  ofstream txtout("testSim.txt",ios::out);
  
  double pzero[3] = {0,0,0};
  
  
  /// calculates here 6planes 0: 0123, 1: 0145, 2: 1256, 3:2367, 4:0347, 5: 4567
  int coind[6][3] ={
    {0,1,2},
    {0,1,4},
    {1,2,5},
    {2,3,6},
    {0,3,4},
    {4,5,6}
  };
  
  for(int k= 0; k<=169; k++){
    for(int j=0; j<=359; j++){
      
      double a[3] = {coxEBAll[k][j][0],coyEBAll[k][j][0],cozEBAll[k][j][0]};
      double b[3] = {coxEBAll[k][j][1],coyEBAll[k][j][1],cozEBAll[k][j][1]};
      double c[3] = {coxEBAll[k][j][2],coyEBAll[k][j][2],cozEBAll[k][j][2]};
      computePlane(a,b,c,p);
      
      double center[3]= {0};
      double cop[8][3];
      for(int n=0; n<8; n++){
	cop[n][0] = coxEBAll[k][j][n];
	cop[n][1] = coyEBAll[k][j][n];
	cop[n][2] = cozEBAll[k][j][n];

	center[0] += 1./8 * coxEBAll[k][j][n];
	center[1] += 1./8 * coyEBAll[k][j][n];
	center[2] += 1./8 * cozEBAll[k][j][n];

	xctEBAll[k][j] += 1./8 * coxEBAll[k][j][n];
	yctEBAll[k][j] += 1./8 * coyEBAll[k][j][n];
	zctEBAll[k][j] += 1./8 * cozEBAll[k][j][n];
      }
      
      /// calculates here 6planes 0: 0123, 1: 0145, 2: 1256, 3:2367, 4:0347, 5: 4567
      for(int n=0; n<6; n++){
	computePlane( cop[coind[n][0]], cop[coind[n][1]],cop[coind[n][2]], p);  
	for(int kk=0; kk<4; kk++){
	  map_planes_eb[k][j].push_back(p[kk]);
	}
      }

    }
  }

  
  
  for(int iz=0; iz<2; iz++){
    for(int k= 1; k<=100; k++){
      for(int j=1; j<=100; j++){
	
	if ( fabs(zEEAll[iz][k][j])<1) continue; //not valid 
	double a[3] = {coxEEAll[iz][k][j][0],coyEEAll[iz][k][j][0],cozEEAll[iz][k][j][0]};
	double b[3] = {coxEEAll[iz][k][j][1],coyEEAll[iz][k][j][1],cozEEAll[iz][k][j][1]};
	double c[3] = {coxEEAll[iz][k][j][2],coyEEAll[iz][k][j][2],cozEEAll[iz][k][j][2]};
	computePlane(a,b,c,p);
	
	double center[3]= {0};
	double cop[8][3];
	for(int n=0; n<8; n++){
	  cop[n][0] = coxEEAll[iz][k][j][n];
	  cop[n][1] = coyEEAll[iz][k][j][n];
	  cop[n][2] = cozEEAll[iz][k][j][n];
	  
	  center[0] += 1./8 * coxEEAll[iz][k][j][n];
	  center[1] += 1./8 * coyEEAll[iz][k][j][n];
	  center[2] += 1./8 * cozEEAll[iz][k][j][n];

	  xctEEAll[iz][k][j] += 1./8 * coxEEAll[iz][k][j][n];
	  yctEEAll[iz][k][j] += 1./8 * coyEEAll[iz][k][j][n];
	  zctEEAll[iz][k][j] += 1./8 * cozEEAll[iz][k][j][n];
	}
	/// calculates here 6planes 0: 0123, 1: 0145, 2: 1256, 3:2367, 4:0347, 5: 4567
	for(int n=0; n<6; n++){
	  computePlane( cop[coind[n][0]], cop[coind[n][1]],cop[coind[n][2]], p);  
	  for(int kk=0; kk<4; kk++){
	    map_planes_ee[iz][k][j].push_back(p[kk]);
	  }
	}
	
      }
    }
  }
  
  float res[10];
  
  cout<<"N " << N <<endl; 

  //  N = 1;
  for(entry =0; entry< N; entry++){    

    if(entry%10==0) cout<<"entry "<< entry <<endl; 

    fChain->GetEntry(entry);
    fChain1->GetEntry(entry);
    
    
    for(int j=0; j<nGenPht; j++){
      float egen = ptGenPht[j]/sin(2*atan(exp(-etaGenPht[j])));
      
      vector<int> vj2 = get5x5CrystalStep(etaGenPht[j],phiGenPht[j],vxGenPht[j],vyGenPht[j],vzGenPht[j]);

      
      if(vj2.size()<1) continue; 
      double pretmin[25];
      double pretminv1[25][10];
      int ind_pretmin[25];
      
      for(int jj=0; jj<25; jj++){
	pretmin[jj] =1E3;
	ind_pretmin[jj] = 0; 
	for(int kk=0; kk<10; kk++){
	  pretminv1[jj][kk] = 1E3; 
	}
      }
      
      if( fabs(etaGenPht[j])<1.5){
	for(int jj=0; jj< int(vj2.size()); jj++){
	  
	  int jmax = vj2[jj];
	  int ieta1 = ietag4EB[jmax];
          int iphi1 = iphig4EB[jmax];
          convxtalid(iphi1,ieta1);
          int ieta = ieta1+85;
          int iphi = getIndphixyzEBAll(iphi1);

	  vector<int> id = idg4EB->at(jmax); 
	  vector<int> pid = pidg4EB->at(jmax); 
	  vector<int> parentid = parentidg4EB->at(jmax); 
	  
	  vector<float> posx = postxg4EB->at(jmax);
          vector<float> posy = postyg4EB->at(jmax);
          vector<float> posz = postzg4EB->at(jmax);

	  
	  vector<float> prex = prexg4EB->at(jmax);
	  vector<float> prey = preyg4EB->at(jmax);
	  vector<float> prez = prezg4EB->at(jmax);
	  vector<float> pret = pretg4EB->at(jmax);
	  vector<float> pree = preeg4EB->at(jmax);
	  vector<int> enter = enterg4EB->at(jmax);
	  
	  for(int k=0; k<int(posx.size()); k++){
	    if(enter[k]==1){
	      if( pretmin[jj] > pret[k]){
		ind_pretmin[jj] = k; 
		pretmin[jj] = pret[k]; 
	      }
	      
	      for(int kk=0; kk<10; kk++){
		if(pree[k]/1000.>1+kk){
		  if( pretminv1[jj][kk] > pret[k]){
		    pretminv1[jj][kk] = pret[k];
		  }
		}
	      }
	    }
	  }
	

  
	  if(ietag4EB[jmax]==8 && iphig4EB[jmax]==9) {
	    
	    int ind = ind_pretmin[jj];
	    double pos[3] = {posx[ind]/10,posy[ind]/10,posz[ind]/10};
	    double pre[3] = {prex[ind]/10,prey[ind]/10,prez[ind]/10};
	    distPreStepPointAllSides(pos, pre, ieta, iphi,res,0,1);
	    
	    distplane = res[0];
	    indplane = int(res[1]+0.1);
	    split_new_x = res[2];
	    split_new_y = res[3];
	    prepos[0] = pre[0];
	    prepos[1] = pre[1];
	    prepos[2] = pre[2];
	    pretime = pretmin[jj];
	    
	    for(int k1=0; k1<10;k1++){
	      pretimev1[k1] = pretminv1[jj][k1];
	    }
	    
	    isconv = convGenPht[j];
	    pids[0] = id[ind];
	    pids[1] = pid[ind];
	    pids[2] = parentid[ind];
            ine = pree[ind];
	    
	    mytree->Fill();
	    
	    fillTH1F("pretmin_0",pretmin[jj]);
	    if(convGenPht[j]==0) {
	      fillTH1F("pretmin_unconv_0",pretmin[jj]);
	    }
	    else {
	      fillTH1F("pretmin_conv_0",pretmin[jj]);
	    }
	    
	  }///for this crystal
	  
	}///all 5x5
	
      }else{//endcap
	

	bool fillxtal1ee = false; 
	
	for(int jj=0; jj< int(vj2.size()); jj++){
	  int jmax = vj2[jj];
	  int ieta = ixg4EE[jmax];
          int iphi = iyg4EE[jmax];
	  int indz = izg4EE[jmax]>0;
	  
	  vector<float> posx = postxg4EE->at(jmax);
          vector<float> posy = postyg4EE->at(jmax);
          vector<float> posz = postzg4EE->at(jmax);
	  
	  vector<float> prex = prexg4EE->at(jmax);
	  vector<float> prey = preyg4EE->at(jmax);
	  vector<float> prez = prezg4EE->at(jmax);
	  vector<int> pid = pidg4EE->at(jmax); 
	  vector<int> id = idg4EE->at(jmax); 
	  vector<int> parentid = parentidg4EE->at(jmax);
	  vector<float> pree = preeg4EE->at(jmax);
	  vector<float> pret = pretg4EE->at(jmax);
	  vector<int> enter = enterg4EE->at(jmax);
	  for(int k=0; k<int(posx.size()); k++){
	    if(enter[k]==1){
              if( pretmin[jj] > pret[k] ){
		ind_pretmin[jj] = k;
                pretmin[jj] = pret[k];
              }
	      for(int kk=0; kk<10; kk++){
                if(pree[k]/1000.>1+kk){
		  if( pretminv1[jj][kk] > pret[k]){
		    pretminv1[jj][kk] = pret[k];
		  }
                }
              }
	    }
	  }//all steps
	  
	  if(ixg4EE[jmax]==20 && iyg4EE[jmax]==50 && izg4EE[jmax] ==1) {
	    fillxtal1ee = true;  //found in 5x5
	    int ind = ind_pretmin[jj];
	    
	    double pos[3] = {posx[ind]/10,posy[ind]/10,posz[ind]/10};
            double pre[3] = {prex[ind]/10,prey[ind]/10,prez[ind]/10};
            distPreStepPointAllSides(pos, pre, ieta,iphi,res,indz,0);
	    
	    distplane = res[0];
	    indplane = int(res[1]+0.1);
	    split_new_x = res[2];
	    split_new_y = res[3];
	    prepos[0] = pre[0];
	    prepos[1] = pre[1];
	    prepos[2] = pre[2];
	    pretime = pretmin[jj];
	    isconv = convGenPht[j];
	    
	    pids[0] = id[ind];
	    pids[1] = pid[ind];
	    pids[2] = parentid[ind];
	    for(int k1=0; k1<10;k1++){
              pretimev1[k1] = pretminv1[jj][k1];
            }
	    
	    ine = pree[ind];

	    mytree->Fill();
	    
            fillTH1F("pretmin_1",pretmin[jj]);
	    if(convGenPht[j]==0){
	      fillTH1F("pretmin_unconv_1",pretmin[jj]);
	    }else {
	      fillTH1F("pretmin_conv_1",pretmin[jj]);
	    }
	  }
	  
	} ///5x5
	
	if( fillxtal1ee == false){
	  //cout<<"not found xtal1ee in 5x5 .. " <<endl; 
	  bool notfilled = false; 

	  for(int k1=0; k1<10; k1++){
	    pretimev1[k1] = 1E3;
	  }
	  
	  for(int j1=0; j1< ng4EE; j1++){
	    
	    if(ixg4EE[j1]==20 && iyg4EE[j1]==50 && izg4EE[j1] ==1) {
	      vector<float> posx = postxg4EE->at(j1);
              vector<float> posy = postyg4EE->at(j1);
              vector<float> posz = postzg4EE->at(j1);

	      vector<float> pret = pretg4EE->at(j1);
	      vector<float> pree = preeg4EE->at(j1);
	      vector<int> enter = enterg4EE->at(j1);
	      vector<float> prex = prexg4EE->at(j1);
	      vector<float> prey = preyg4EE->at(j1);
	      vector<float> prez = prezg4EE->at(j1);
	      vector<int> id = idg4EE->at(j1); 
	      vector<int> pid = pidg4EE->at(j1); 
	      vector<int> parentid = parentidg4EE->at(j1);
	      
	      float tmin = 1E9; 
	      int ind  = 0; 
	      for(int k=0; k<int(pret.size()); k++){
		if(enter[k]==0) continue; 
		if(tmin> pret[k] ){
		  tmin = pret[k];
		  ind = k; 
		}
		for(int k1=0; k1<10;k1++){
		  if( pree[k1]/1000.>1+k1){
		    if( pretimev1[k1] < pret[k]) pretimev1[k1] = pret[k];
		  }
		}
	      }
	      notfilled = true; 
	      fillTH1F("pretmin_1",tmin);
	      if(convGenPht[j]==0){
		fillTH1F("pretmin_unconv_1",tmin);
	      }else {
		fillTH1F("pretmin_conv_1",tmin);
	      }
	      
	      int ieta = ixg4EE[j1];
	      int iphi = iyg4EE[j1];
	      int indz = izg4EE[j1]>0;
	      
	      double pos[3] = {posx[ind]/10,posy[ind]/10,posz[ind]/10};
	      double pre[3] = {prex[ind]/10,prey[ind]/10,prez[ind]/10};
	      distPreStepPointAllSides(pos, pre, ieta,iphi,res,indz,0);
	      
	      distplane = res[0];
	      indplane = int(res[1]+0.1);
	      split_new_x = res[2];
	      split_new_y = res[3];
	      prepos[0] = pre[0];
	      prepos[1] = pre[1];
	      prepos[2] = pre[2];
	      pretime = tmin; 
	      isconv = convGenPht[j];
	      pids[0] = id[ind];
	      pids[1] = pid[ind];
	      pids[2] = parentid[ind];
	      ine = pree[ind];
	      mytree->Fill();
	      
	      
	    }
	  }
	  
	  if( !notfilled){
	    cout<<"not found at all?? "<< entry <<endl; 
	    //return; 
	  }
	}
	
      } //endcap
      
      
    } ///all genPht
    
    
    
  }
  
  

  mytree->Write();
  fnew->Write();
  fnew->Close();
  
}
