vector<string> l1menu_algonames; 

map<string, int> l1menu_prescale1; ///5E33
map<string, int> l1menu_prescale2; ////7E33

map<string, int> l1menu_count1;
map<string, int> l1menu_count2;


void loadL1SeedsAndPrescale_5e33_7e33(){
  
  ifstream txtin; 
  txtin.open("L1_draft4_v2_edited_v2_simple.txt",ios::in);
  
    
  string lname; 
  int ps1; 
  int ps2; 
  while (txtin.good()){
    txtin>> lname >> ps1 >> ps2;
    if( txtin.eof()) break; 
    
    std::vector<string>::iterator its = find( l1menu_algonames.begin(),
					      l1menu_algonames.end(),
					      lname);
    if( its == l1menu_algonames.end()){
      l1menu_algonames.push_back(lname);
    }else{
      cout<<"warning same L1sees in the menu " << lname.c_str()<<endl; 
    }
    if( ps1 ==0 || ps2 ==0 ){
      cout<<"ps ? " << lname.c_str()<<" "<< ps1 <<" "<< ps2 <<endl; 
      ps1 = 99999;
      ps2 = 99999;
    }
    l1menu_prescale1[lname] = ps1;
    l1menu_prescale2[lname] = ps2;
  }
  cout<<" loadL1SeedsAndPrescale_5e33_7e33  " << l1menu_algonames.size()<<endl; 
  
  
}

void checkMissingL1Prescale(){
  for(unsigned int j=0; j< l1algoName->size() ;j++){
    string l1 = l1algoName->at(j);
    if( l1menu_prescale1[l1] == 0 ){
      l1menu_prescale1[l1] = 9999;
      cout<<" checkMissingL1Prescale1 " << l1.c_str()<<endl; 
    }
    if( l1menu_prescale2[l1] == 0){
      l1menu_prescale2[l1] = 9999;
      cout<<" checkMissingL1Prescale2 " << l1.c_str()<<endl; 
    }
  }
}



void prescale_L1seeds(){
  
  l1bitFired_prescale1->clear();
  l1bitFired_prescale2->clear();
  
  for(unsigned int j=0; j< l1algoName->size() ;j++){
    string l1 = l1algoName->at(j);
    if( isL1PathFired(l1) ){
      l1menu_count1[l1] ++; 

      if( l1menu_prescale1[l1] == 0 || l1menu_prescale2[l1] == 0){
	cout<<" pre-scale 0 ? " << l1.c_str()<<" "<< l1menu_prescale1[l1]<<" "<< l1menu_prescale2[l1]<<endl;
      }
      
      if(  int( l1menu_count1[l1] %  l1menu_prescale1[l1] ) == 0  ){
	l1bitFired_prescale1->push_back(j);
      }
      l1menu_count2[l1] ++; 
      if( int( l1menu_count2[l1] % l1menu_prescale2[l1] ) == 0  ){
	l1bitFired_prescale2->push_back(j);
      }
    }
  }
  
  
}
  
