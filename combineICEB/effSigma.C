//******************************************//
// author:
//
//  This function evaluates the effective 
//  sigma of a histogram. 
//
//
//  TODO list:
//
//
//
//  Modified by:
//
//******************************************//




//Double_t effSigma(TH1 * hist)

void effSigma(TH1 * hist,double res[])

{

  TAxis *xaxis = hist->GetXaxis();
  Int_t nb = xaxis->GetNbins();
  if(nb < 10) {
    cout << "effsigma: Not a valid histo. nbins = " << nb << endl;
    return ; 
  }
  
  Double_t bwid = xaxis->GetBinWidth(1);
  if(bwid == 0) {
    cout << "effsigma: Not a valid histo. bwid = " << bwid << endl;
    return ; 
  }
  Double_t xmax = xaxis->GetXmax();
  Double_t xmin = xaxis->GetXmin();
  Double_t ave = hist->GetMean();
  Double_t rms = hist->GetRMS();

  Double_t total=0.;
  for(Int_t i=0; i<nb+2; i++) {
    total+=hist->GetBinContent(i);
  }
  //  if(total < 100.) {
  if(total < 50.) {
    cout << "effsigma: Too few entries " << total << endl;
    return ; 
  }
  Int_t ierr=0;
  Int_t ismin=999;
  
  Double_t rlim=0.683*total;
  Int_t nrms=int(rms/(bwid));    // Set scan size to +/- rms
  if(nrms > nb/10) nrms=nb/10; // Could be tuned...

  Double_t widmin=9999999.;
  
  // cout<<"iscan: "<< xmin <<" "<< xmax <<" "<< ave <<" "<< rms <<" "<< rlim <<" "<< nrms <<endl; 
  
  int binlow = 1; 
  int binhigh = nb; 
  
  double mean_min = ave; 
  
  int binlow2 = 1; 
  int binhigh2 = nb; 
  
  
  
  for(Int_t iscan=-nrms;iscan<nrms+1;iscan++) { // Scan window centre
    Int_t ibm=(ave-xmin)/bwid+1+iscan;
    Double_t x=(ibm-0.5)*bwid+xmin;
    Double_t xj=x;
    Double_t xk=x;
    Int_t jbm=ibm;
    Int_t kbm=ibm;
    Double_t bin=hist->GetBinContent(ibm);
    total=bin;

    //  cout<<"ibm: "<< iscan <<" "<< ibm <<" "<< xj<<" "<< xk <<" "<<jbm <<" "<<kbm <<endl; 
    
    double mean = bin * x; 
    
    
    for(Int_t j=1;j<nb;j++){
      if(jbm < nb) {
        jbm++;
        xj+=bwid;
        bin=hist->GetBinContent(jbm);
        total+=bin;
	
        if(total > rlim) break;
      }
      else ierr=1;
      if(kbm > 0) {
        kbm--;
        xk-=bwid;
        bin=hist->GetBinContent(kbm);
        total+=bin;
	
        if(total > rlim) break;
      }
      else ierr=1;
    }
    Double_t dxf=(total-rlim)*bwid/bin;
    Double_t wid=(xj-xk+bwid-dxf)*0.5;

    
    //cout<<"dxf: "<< total <<" "<< rlim <<" "<< bwid <<" "<< bin <<" "<< wid <<" "<< xj<<" "<<xk<<" "<< dxf <<endl; 

    if(wid < widmin) {
      widmin=wid;
      ismin=iscan;
      
      binlow = kbm; 
      binhigh = jbm; 


      binlow2 = ibm - 2*(ibm - kbm);  
      binhigh2 = ibm + 2*(-ibm + jbm);  
      
    }
    
    // cout<<"widmin: "<< wid<<" "<< widmin <<" "<< ismin <<" "<<ibm <<" "<< jbm<<" "<< kbm <<" "<<x <<" "<< xj <<" "<< xk <<" "<<(ibm-0.5)*bwid+xmin<<" "<<
    //" "<< (jbm-0.5)*bwid+xmin <<" "<<(kbm-0.5)*bwid+xmin <<endl; 
    
  }
  if(ismin == nrms || ismin == -nrms) ierr=3;
  if(ierr != 0) cout << "effsigma: Error of type " << ierr << endl;
  
  //return widmin;
  
  double mean = 0; 
  double n =0; 
  for(int b = binlow ; b<= binhigh; b++){
    
    if( b<1 || b> nb) continue; 

    double x= (b-0.5)*bwid+xmin;
    mean += hist->GetBinContent(b) * x;
    n += hist->GetBinContent(b); 
  }
  
  mean /= n; 
  
  double mean2 = 0; 
  double n2 =0; 
  for(int b = binlow2 ; b<= binhigh2; b++){

    if( b<1 || b> nb) continue; 
    
    double x= (b-0.5)*bwid+xmin;
    mean2 += hist->GetBinContent(b) * x;
    n2 += hist->GetBinContent(b); 
  }
  
  mean2 /= n2; 
    

  //cout<<widmin <<" "<<mean <<endl; 
  
  res[0] = mean; 
  res[1] = mean2; 
  
  res[2] = widmin; 
  
  
}
