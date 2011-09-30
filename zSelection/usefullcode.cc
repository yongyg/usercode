
float getMaximumY(TH1 *h1, TH1 *h2){
  return h1->GetMaximum() > h2->GetMaximum() ? h1->GetMaximum() : h2->GetMaximum();
}

int getBinNumberwithLabel(TH1 *h1, string bnam){
  int nbin1 = h1->GetNbinsX();
  
   
  for(int b=1; b<=nbin1; b++){
    if( h1->GetXaxis()->GetBinLabel(b) == bnam){
      return b; 
    }
  }
  cout<<"not found " << bnam.c_str()<<endl; 
  return -1; 
}




void AddHistByQuadrature(TH1 *h1, float xx){
  
  int nbin1 = h1->GetNbinsX();

  
  for(int b=1; b<=nbin1; b++){
    
    float tmp1 =  h1->GetBinContent(b);
    h1->SetBinContent(b,sqrt(tmp1*tmp1 + xx*xx));
    
  }
    
}



void subtractHistByQuadrature(TH1 *h1, float xx){
  
  int nbin1 = h1->GetNbinsX();

  
  for(int b=1; b<=nbin1; b++){
    
    float tmp1 =  h1->GetBinContent(b);
    if(tmp1 > xx){
      h1->SetBinContent(b,sqrt(tmp1*tmp1 - xx*xx));
    }else{
      ///h1->SetBinContent(b,sqrt(tmp1*tmp1 - xx*xx));
    }
  }
    
}




void subtractTwoHistQuadrature(TH1 *h1, TH1 *h2, TH1 *h3){
  
  int nbin1 = h1->GetNbinsX();
  int nbin2 = h2->GetNbinsX();

  if( nbin1 != nbin2) {
    cout<<" subtractTwoHistQuadrature nbin " << nbin1 <<" "<<nbin2 <<endl; 
  }
  
  for(int b=1; b<=nbin1; b++){
    
    float tmp1 =  h1->GetBinContent(b);
    float tmp1err =  h1->GetBinError(b);
    float tmp2 =  h2->GetBinContent(b);
    float tmp2err =  h2->GetBinError(b);
    
    float diff = -1; 
    float diffErr = 0; 
    if( tmp1 >= tmp2){
      diff = sqrt( tmp1*tmp1 - tmp2*tmp2);
      diffErr = sqrt(  pow( tmp1/diff*tmp1err,2) +  pow( tmp2/diff*tmp2err,2) );

      //cout<<"check " << tmp1 <<" "<< tmp2 <<" "<< diff <<endl; 

    }
    h3->SetBinContent(b,diff);
    h3->SetBinError(b,diffErr);
  }
  
  
}




double one_sigma_up(std::vector<double> limits) {
 double v_one_sigma_up = 0;
 sort(limits.begin(),limits.end());
 int samplet_size = limits.size();
 int v_one_sigma_up_index = 0;
 v_one_sigma_up_index = samplet_size * 841 / 1000;
 v_one_sigma_up = limits[v_one_sigma_up_index-1];
 return v_one_sigma_up;
}

double one_sigma_down(std::vector<double> limits) {
 double v_one_sigma_down = 0;
 sort(limits.begin(),limits.end());
 int samplet_size = limits.size();
 int v_one_sigma_down_index = 0;
 v_one_sigma_down_index = samplet_size * 159 / 1000;
 v_one_sigma_down = limits[v_one_sigma_down_index];
 return v_one_sigma_down;
}

double two_sigma_up(std::vector<double> limits) {
 double v_two_sigma_up = 0;
 sort(limits.begin(),limits.end());
 int samplet_size = limits.size();
 int v_two_sigma_up_index = 0;
 v_two_sigma_up_index = samplet_size * 979 / 1000;
 v_two_sigma_up = limits[v_two_sigma_up_index-1];
 return v_two_sigma_up;
}

double two_sigma_down(std::vector<double> limits) {
 double v_two_sigma_down = 0;
 sort(limits.begin(),limits.end());
 int samplet_size = limits.size();
 int v_two_sigma_down_index = 0;
 v_two_sigma_down_index = samplet_size * 21 / 1000;
 v_two_sigma_down = limits[v_two_sigma_down_index];
 return v_two_sigma_down;
}






float getIntegralErrorHistogram(TH1F *h1){
  
  float err = 0; 
  for(int b=1; b<= h1->GetNbinsX(); b++){
    err += pow( h1->GetBinError(b),2); 
    
  }
  return sqrt(err);
}


void getIntegralAndWtErrorHistogram(TH1F *h1, int startBin, float res[]){
  
  float tot = 0; 
  float err = 0; 
  for(int b=startBin; b<= h1->GetNbinsX(); b++){
    err += pow( h1->GetBinError(b),2); 
    tot +=  h1->GetBinContent(b);
  }
  
  res[0] = tot; 
  res[1] = sqrt(err);
  
  
}



void printbinContentHistogram(TH1 *h1){
  
  for(int b=1; b<= h1->GetNbinsX(); b++){
    cout<<"bin: "<< h1->GetBinLowEdge(b)<<"to"<<  h1->GetBinLowEdge(b)+h1->GetBinWidth(b) <<" "<< h1->GetBinContent(b)<<"+/-"<< h1->GetBinError(b)<<endl; 
  }
  
}




void calcRatioBayes(float n1, float n2, float &r, float &rup, float &rlow){
  
  //compute using Bayes                                                                                                                                 
  TH1D h1("dummy1","",1,1,2);
  h1.SetBinContent(1,n1);
  
  TH1D h2("dummy2","",1,1,2);
  h2.SetBinContent(1,n2);
  
  TGraphAsymmErrors g;
  g.BayesDivide(&h1,&h2);
  r = g.GetY()[0];
  rup = g.GetErrorYhigh(0);
  rlow = g.GetErrorYlow(0);
  
}




double ErrorInProduct(double x, double errx, double y, 
                      double erry, double corr) {
   double xFrErr = errx/x;
   double yFrErr = erry/y;
   return sqrt(pow(xFrErr,2) + pow(yFrErr,2) + 2.0*corr*xFrErr*yFrErr)*x*y;
}



void  getErrorInRatio(double x, double errx, double y, 
                      double erry, double res[]) {

  if( x==0 || y==0){
    cout<<"getErrorInRatio x/y" << x <<" "<<y<<endl; 
    exit(1);
  }
  
  res[0] = x/y; 
  res[1] = x/y *sqrt( pow(errx/x,2) + pow(erry/y,2));

}



void get3Max(double x1, double x2, double x3,double res[]){
  double max = x1;
  res[1] = 0;

  if(max<x2) {
    max = x2;
    res[1] = 1;
  }
  if(max<x3){
    max = x3;
    res[1] = 2;
  }

  res[0] = max;
  //  return max;

}



float getMaxOfThree(float x1, float x2, float x3){
  
  float xmax = x1 > x2 ? x1: x2; 
  xmax = xmax > x3 ? xmax: x3; 
  return xmax; 
  
  
}


void generateToyDataFromHist(TH1* h1, TH1 *h2){
  
  int nbinx1 = h1->GetNbinsX();
  int nbinx2 = h2->GetNbinsX();
  
  if(nbinx1 != nbinx2){
    cout<<" generateToyDataFromHist " << nbinx1 <<" "<<  nbinx2 <<endl; 
    exit(1);
  }
  TRandom3 *grand = new TRandom3(0);
  for(int b=1 ; b<= nbinx1; b++){
    float tmp = h1->GetBinContent(b);
    float tmp1 = grand->Poisson(tmp);
    h2->SetBinContent(b,tmp1);
  }
  
  
}


void calcTotalErorrPDF(int npdf, double X0,double X0Err,double XP[100],double XPErr[100],double XM[100],double XMErr[100],double results[]){
  
  
  double dXP = 0; 
  double dXM = 0; 
  double simpleDX = 0; 
  
  double dsimpleDX = 0; 
  double ddXP = 0; 
  double ddXM = 0; 
  
  double res[10];
  
  if(npdf%2 !=0) {
    cout<<" calcTotalErorrPDF npdf: "<<npdf <<endl; 
    exit(1);
  }
    
  
  for( int j=1; j<= npdf/2; j++){
    
    get3Max(XP[j]-X0,XM[j]-X0,0,res);
    double dxP = res[0];
    if( fabs(res[1]-0)<0.01){
      ddXP += 4*pow(XP[j]-X0,2)*(XPErr[j]*XPErr[j] + X0Err*X0Err);
    }else if( fabs(res[1]-1)<0.01){
      ddXP += 4*pow(XM[j]-X0,2)*(XMErr[j]*XMErr[j] + X0Err*X0Err);
    }
    
    get3Max(X0-XP[j],X0-XM[j],0,res);
    double dxM = res[0];
    if( fabs(res[1]-0)<0.01){
      ddXM += 4*pow(XP[j]-X0,2)*(XPErr[j]*XPErr[j] + X0Err*X0Err);
    }else if( fabs(res[1]-1)<0.01){
      ddXM += 4*pow(XM[j]-X0,2)*(XMErr[j]*XMErr[j] + X0Err*X0Err);
    }
    
    dXP += dxP*dxP; 
    dXM += dxM*dxM; 
    
    simpleDX +=  pow(XP[j]-XM[j],2);

    dsimpleDX += pow(XP[j]-XM[j],2)*(XPErr[j]*XPErr[j]+XMErr[j]*XMErr[j]);
    
    
    ///cout<<"checkme: "<<j<<"  "<< dXP <<" "<< dXM <<endl; 
    

  }


  dXP = sqrt(dXP);
  dXM = sqrt(dXM);

  ddXP = sqrt(ddXP);

  if( dXP>0){
    ddXP = ddXP/(2*dXP);
  }
  ddXM = sqrt(ddXM);
  if( dXM >0){
    ddXM = ddXM/(2*dXM);
  }
  
  
  
  
  simpleDX = 0.5* sqrt(simpleDX);
  dsimpleDX =  sqrt(dsimpleDX);
  if(simpleDX >0){
    dsimpleDX = dsimpleDX/(2*simpleDX);
  }
  
  ////asymmetry errors
  results[0] = dXP; 
  results[1] = dXM; 
  ///master formula
  results[2] = simpleDX; 
  
  ////error of above
  results[3] = ddXP; 
  results[4] = ddXM; 
  results[5] = dsimpleDX; 
  
  
}

int mostProb(float b){
  
  int nbest = int(b+1); 
  float prob_max = 0; 
  for(int n = int(b-1) ; n<= int(b+1); n++){
    if( n<0) n =0; 


    float prob = TMath::Poisson(n,b);

    //    cout<<prob<<" "<<n<<" "<<b<<endl;
    
    if( prob > prob_max){
      prob_max = prob; 
      nbest = n; 
    }
  }
  return nbest;
  
}


float getBinContentHistogram(TH1F *hhtmp, float pt){
  
  for(int b=1; b<= hhtmp->GetNbinsX(); b++){
    if( pt >= hhtmp->GetBinLowEdge(b) && pt < hhtmp->GetBinLowEdge(b) +  hhtmp->GetBinWidth(b)){
      return  hhtmp->GetBinContent(b);
    }
  }
  return -999; 
  
  
}






void weightedAeverageTwoHistogram(TH1F *hhtmp1, TH1F *hhtmp2, TH1F *hhres){
  
  if( hhtmp1->GetNbinsX() != hhtmp2->GetNbinsX()){
    cout<<"warning  weightedAeverageTwoHistogram  "<< hhtmp1->GetNbinsX() <<" "<< hhtmp2->GetNbinsX() <<endl; 
    exit(1); 
  }

  for(int b=1; b<= hhtmp1->GetNbinsX(); b++){
    float tmp1 =  hhtmp1->GetBinContent(b);
    float tmp2 =  hhtmp2->GetBinContent(b);
    float tmp1e =  hhtmp1->GetBinError(b);
    float tmp2e =  hhtmp2->GetBinError(b);
    
    float av = 0; 
    float ave = 0; 
    if( tmp1 > 0 && tmp2 >0){
      av = (tmp1/(tmp1e*tmp1e) + tmp2/ (tmp2e*tmp2e)) / ( 1./(tmp1e*tmp1e) + 1./(tmp2e*tmp2e));
      ave =  sqrt(1/ (1./(tmp1e*tmp1e) + 1./(tmp2e*tmp2e)) ); 
    }else if (tmp1 > 0 && tmp2 ==0){
      av = tmp1; 
      ave = tmp1e; 
    }else if( tmp1 ==0 && tmp2 >0){
      av = tmp2; 
      ave = tmp2e; 
    }else{
      cout<<"no data for combine.."<< b<<endl; 
    }
    
    hhres->SetBinContent(b,av);
    hhres->SetBinError(b,ave);
    
  }

}




float getFractionOfHistogram(TH1F *hhtmp, double xmin, double xmax){
  
  double sum = hhtmp->Integral(); 
  if( sum <=0) {
    cout<<" getFractionOfHistogram zero.." << sum << " "<< hhtmp->GetEntries();
    exit(1); 
  }
  
  int nbinsx = hhtmp->GetNbinsX();
  int bin1 =1;  
  int bin2 = nbinsx; 
  
  for(int b=1; b<= nbinsx; b++){
    if( xmin >= hhtmp->GetBinLowEdge(b) && xmin -hhtmp->GetBinLowEdge(b) - hhtmp->GetBinWidth(b) < 1E-6 ){
      bin1 = b; 
    }
    if( xmax > hhtmp->GetBinLowEdge(b) && xmax - hhtmp->GetBinLowEdge(b) - hhtmp->GetBinWidth(b) < 1E-6 ){
      bin2 = b; 
    }
  }
  
  ///cout<<"checkme"<<bin1 <<" "<<bin2 <<" "<< hhtmp->GetBinLowEdge(bin1)  <<" "<< xmax - hhtmp->GetBinLowEdge(bin2) - hhtmp->GetBinWidth(bin2) <<" "<< xmax - hhtmp->GetBinLowEdge(bin2-1) - hhtmp->GetBinWidth(bin2-1)<<endl; 
  
  return hhtmp->Integral(bin1,bin2) / sum; 
  
  
}




void cloneHistogramWithWeight(TH1 *hhtmp1, TH1* hhtmp2, float wt){
  
  int nbins = hhtmp1->GetNbinsX();
  if( nbins != hhtmp2->GetNbinsX()){
    cout<<"error.cloneHistogramWithWeight  "<<nbins <<" "<< hhtmp2->GetNbinsX() << endl; 
    exit(1); 
  }
  
  for(int j=1; j<= nbins+1; j++){ ///bins+1 for overflowbins
    
    float tmp = hhtmp1->GetBinContent(j) * wt; 

    hhtmp2->SetBinContent(j,tmp);

  }
  
  
}




void cloneHistogramWithWeightv1(TH1 *hhtmp1, TH1* hhtmp2, float wt){
  
  int nbins = hhtmp1->GetNbinsX();
  if( nbins >  hhtmp2->GetNbinsX()){
    cout<<"error.cloneHistogramWithWeightv1  "<<nbins <<" "<< hhtmp2->GetNbinsX() << endl; 
    exit(1); 
  }
  
  for(int j=1; j<= nbins; j++){
    
    float tmp = hhtmp1->GetBinContent(j) * wt; 

    hhtmp2->SetBinContent(j,tmp);

  }
  
  
}




///check which bin the value in histogram
int getBinOfHistogram(TH1 *hh, float tmp){
  
  int nbins = hh->GetNbinsX();
  
  if( tmp < hh->GetXaxis()->GetXmin() || tmp >= hh->GetXaxis()->GetXmax()) return -1; 
  
  for(int b = 1; b<= nbins; b++){
    
    if( tmp >= hh->GetBinLowEdge(b) && tmp < hh->GetBinLowEdge(b) +  hh->GetBinWidth(b)) return b; 
    
  }
  
  cout<<"error.. getBinOfHistogram "<<tmp <<endl;
  
  return -1; 
  
  
}





///check which bin the value in histogram
int getBinOfHistogramv1(TH1 *hh, float tmp){
  
  int nbins = hh->GetNbinsX();
  
  if( tmp < hh->GetXaxis()->GetXmin() ) return 1; 
  if( tmp >= hh->GetXaxis()->GetXmax() ) return nbins; 
    
  for(int b = 1; b<= nbins; b++){
    
    if( tmp >= hh->GetBinLowEdge(b) && tmp < hh->GetBinLowEdge(b) +  hh->GetBinWidth(b)) return b; 
    
  }
  
  cout<<"error.. getBinOfHistogram "<<tmp <<endl;
  exit(1);
  
  return -1; 
  
  
}






double BinP(int N, double p, int x1, int x2) {
   double q=p/(1-p); 
   int k=0; 
   double v = 1; 
   double s=0; 
   double tot=0.0;
    while(k<=N) {
       tot=tot+v;
       if(k>=x1 & k<=x2) { s=s+v; }
       if(tot>1e30){s=s/1e30; tot=tot/1e30; v=v/1e30;}
       k=k+1; 
       v=v*q*(N+1-k)/k;
    }
    return s/tot;
}




void ClopperPearsonLimits(double numerator, double denominator, double &ratio,
double &lowerLimit, double &upperLimit, const double CL_low=1.0, 
const double CL_high=1.0) 
{  
  //Confidence intervals are in the units of \sigma.
    

  ratio = numerator/denominator;
   
// first get the lower limit
   if(numerator==0)   lowerLimit = 0.0; 
   else { 
      double v=ratio/2; 
      double vsL=0; 
      double vsH=ratio; 
      double p=CL_low/100;
      while((vsH-vsL)>1e-5) { 
	if(BinP(int(denominator),v,int(numerator),int(denominator))>p) 
         { vsH=v; v=(vsL+v)/2; } 
         else { vsL=v; v=(v+vsH)/2; } 
      }
      lowerLimit = v; 
   }
   
// now get the upper limit
   if(numerator==denominator) upperLimit = 1.0;
   else { 
      double v=(1+ratio)/2; 
      double vsL=ratio; 
      double vsH=1; 
      double p=CL_high/100;
      while((vsH-vsL)>1e-5) { 
	if(BinP(int(denominator),v,0,int(numerator))<p) { vsH=v; v=(vsL+v)/2; } 
         else { vsL=v; v=(v+vsH)/2; } 
      }
      upperLimit = v;
   }
}






