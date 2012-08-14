///read initi calibration constants for each crystal
void readInterCalibEndcap_GR09_V8(){

  TString file_input = workingDirectory + "/interCalib_START3X_V25B_endcap_mc.txt";
  cout<<"READING fro inter-calib file "<<file_input<<endl;
  ifstream inputcc(file_input,ios::in);
  if (inputcc.fail()){
    cout<<"error open file endcap.. " << file_input<<endl;
    exit(1);
  }
  
  int eta; 
  int phi; 
  int iz; 
  float c = 1; 
  float c1 = 1; 
  
  int n = 0; 
  
  int isc; 
  
  while (inputcc.good()){
    inputcc >> iz >> eta >> phi >>isc >> c >>c1; 
    int izz = iz < 0? 0:1; 
    validRecHitEndCap[izz][eta][phi] = 1; 
    n ++; 
    
  }
  
  cout<<n<<" constants read. "<<endl;
  
}
