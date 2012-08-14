#include "dataChain/datachain_etaee2012CGRPV39.cc"
#include "dataChain/datachain_pizee2012CGRPV39.cc"

void addDataChain(){
  
   if(dataflag==1){
     //datachain_test();
   }
   else if(dataflag==2){
     datachain_pizee2012CGRPV39();
   }else if(dataflag==3){
     datachain_etaee2012CGRPV39();
   }
  
   else{
    cout<<"dataflag NA. " <<dataflag <<endl; 
    exit(1);
   }
   
}
