


vector<int> run_goodLS; 

map<int,vector<int> > map_runLS; 

vector<int> goodRunList; 

vector<int> crabRunList; 




//get LS ranges
/// t3-susy [data/analysisMX] cat Cert_132440-142664_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt | grep : | grep - | sed s/"-"/" "/g | sed s/":"/" "/g | sed s/"'"/""/g | sed s/","/""/g 


///this is from the offical JSON files of reReco or promptReco
///not exactly the same as crabRunList. 

///update for multiple Json file 

void getLSrangeofEachRuns(vector<string> certfiles){
  
  if(certfiles.size() <=0){
    cout<<" getLSrangeofEachRuns No certfile  " << endl; 
    exit(1);
  }
  
  
  int run; 
  int ls1;
  int ls2;
  
  int prerun = 0; 
  
  vector<int> lsbl; 
  
  for(int n=0; n< int(certfiles.size()); n++){
    
    cout<<"certfile: "<< certfiles[n].c_str() <<endl; 
    ifstream txtin(certfiles[n].c_str(),ios::in);
    
    while( txtin.good()){
    
      txtin >> run >> ls1 >> run >> ls2; 
      
      ///  cout<<run<<" "<< ls1 <<" "<< ls2 <<endl; 
      bool newrun = false; 
      if(run != prerun){
	newrun = true; 
	
      }
      
      vector<int>::iterator it = find( goodRunList.begin(),goodRunList.end(),run);
      if( it == goodRunList.end()){
	goodRunList.push_back(run); 
      }
      
      /// cout<<"run: "<<prerun<<" "<<run <<endl; 
      
    
      if( newrun){
      
	if( int(lsbl.size()) >0){
	  //insert into map
	  map_runLS.insert( make_pair(prerun,lsbl)); 
	}
	
	prerun = run; 
	
	lsbl.clear(); 
	
	
      }
    
      lsbl.push_back(ls1);
      lsbl.push_back(ls2);
      
    
    }
  }
  
  
  map_runLS.insert( make_pair(run,lsbl));
  
  
  cout<<"total runs: "<< goodRunList.size()<<" "<< map_runLS.size()<<endl; 
    
  
  for(int n=0; n< int( goodRunList.size());n++){
    
    map<int,vector<int> > :: const_iterator iterMap = map_runLS.find(goodRunList[n]); 
    
    if( iterMap != map_runLS.end()){
      
      
      vector<int> tmp = iterMap->second; 
      
      for(int j=0; j< int(tmp.size())/2; j++){
	if( n<10 || n >int( goodRunList.size())-10 ) cout<< goodRunList[n]<<" "<< tmp[j*2] <<" "<< tmp[2*j+1]<<endl; 
      }
      
    }else{
      cout<<"wrongcheckmaprunLS.."<<endl; 
      exit(1);
    }
  }
    
  
}



bool checkLumiBlockofRun(){
  
  
  if(dataOrMC==2) return true; 
  
  bool goodLS = false; 
  
  map<int,vector<int> > :: const_iterator iterMap = map_runLS.find(runNumber); 
  
  if( iterMap != map_runLS.end()){
    vector<int> tmp = iterMap->second; 
    
    for(int j=0; j< int(tmp.size())/2; j++){
      if( tmp[j*2] <= lumiBlock && lumiBlock <= tmp[2*j+1]){
	goodLS = true; 
	break; 
      }
    }
    
  }else{
    cout<<"wrong.. no run found.." << runNumber <<endl; 
    //exit(1); 
    return false; 
  }
  
  return goodLS; 
  
  
}
