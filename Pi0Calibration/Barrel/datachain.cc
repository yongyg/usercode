#include "dataChain/datachain_test.cc"


void addDataChain(){
  
  if(dataflag==1){
    datachain_test();
  }
  
  else{
    cout<<"dataflag NA. " <<endl; 
    exit(1);
  }
  
}
