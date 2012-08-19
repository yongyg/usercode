


vector<int> run_goodLS; 

map<int,vector<int> > map_runLS; 

vector<int> goodRunList; 

vector<int> crabRunList; 


const int MAX_CHARS_PER_LINE = 99999;
const int MAX_TOKENS_PER_LINE = 99999;
const char* const DELIMITER = " ";

//get LS ranges
///cat Cert_132440-142664_7TeV_StreamExpress_Collisions10_CMSSWConfig.txt | grep : | grep - | sed s/"-"/" "/g | sed s/":"/" "/g | sed s/"'"/""/g | sed s/","/""/g 

///update for multiple Json file 

void getLSrangeofEachRuns(vector<string> certfiles){
  
  int run; 
  
  for(int n=0; n< int(certfiles.size()); n++){
    
    cout<<"certfile: "<< certfiles[n].c_str() <<endl; 
    ifstream txtin(certfiles[n].c_str(),ios::in);
    if (txtin.fail()) {
      cout<<" file can not be opened!" <<endl; 
      exit(1);
    }
    

    while( txtin.good()){
      
      char buf[MAX_CHARS_PER_LINE];
      txtin.getline(buf, MAX_CHARS_PER_LINE);

      if(txtin.eof()) break; 
      
      const char* token[MAX_TOKENS_PER_LINE] = {0}; // initialize to 0
      token[0] = strtok(buf, DELIMITER); // first token
      vector<int> lsbl;
      run =  atoi(token[0]);
      if (token[0]){
	for (int k = 1; k < MAX_TOKENS_PER_LINE; k++){
	  token[k] = strtok(0, DELIMITER); // subsequent tokens
	  if (!token[k]) break; // no more tokens
	  int val = atoi(token[k]);
	  lsbl.push_back(val);
	}
      }
      vector<int>::iterator it = find( goodRunList.begin(),goodRunList.end(),run);
      if( it == goodRunList.end()){
	goodRunList.push_back(run);
      }
      if( int(lsbl.size())< 2 || int(lsbl.size()) %2!=0){
	cout<<"wrong format of goodLS file!! " << run <<" "<< lsbl.size()<<endl; 
	exit(1);
      }
      //insert into map
      map_runLS.insert( make_pair(run,lsbl)); 
      
    }
    
  }
  
  
  cout<<"total runs: "<< goodRunList.size()<<" map size "<< map_runLS.size()<<endl; 
  
  for(int n=0; n< int( goodRunList.size());n++){
    
    map<int,vector<int> > :: const_iterator iterMap = map_runLS.find(goodRunList[n]); 
    
    if( iterMap != map_runLS.end()){
      vector<int> tmp = iterMap->second; 
      for(int j=0; j< int(tmp.size())/2; j++){
	if( n<1 || n >int( goodRunList.size())-2 ) cout<< goodRunList[n]<<" "<< tmp[j*2] <<" "<< tmp[2*j+1]<<endl; 
      }
    }else{
      cout<<"wrongcheckmaprunLS.."<<endl; 
      exit(1);
    }
  }
  
  
}



bool checkLumiBlockofRun(){
  
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
    return false; 
  }
  
  return goodLS; 
  
    
}
