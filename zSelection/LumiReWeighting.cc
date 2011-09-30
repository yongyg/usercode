//// from  PhysicsTools_Utilities_LumiReWeighting_cc


//int nWarning = 0; 

/**
  \class    LumiReWeighting LumiReWeighting.h "PhysicsTools/Utilities/interface/LumiReWeighting.h"
  \brief    Class to provide lumi weighting for analyzers to weight "flat-to-N" MC samples to data

  This class will trivially take two histograms:
  1. The generated "flat-to-N" distributions from a given processing (or any other generated input)
  2. A histogram generated from the "estimatePileup" macro here:

  https://twiki.cern.ch/twiki/bin/view/CMS/LumiCalc#How_to_use_script_estimatePileup

  and produce weights to convert the input distribution (1) to the latter (2).

  \author Salvatore Rappoccio, modified by Mike Hildreth
  
*/

TH1D *weights_; 
double WeightOOTPU_[25][25];

void LumiReWeighting_weightOOT_init() {

  // The following are poisson distributions with different means, where the maximum
  // of the function has been normalized to weight 1.0
  // These are used to reweight the out-of-time pileup to match the in-time distribution.
  // The total event weight is the product of the in-time weight, the out-of-time weight,
  // and a residual correction to fix the distortions caused by the fact that the out-of-time
  // distribution is not flat.

  static double weight_24[25] = {
    0,
    0,
    0,
    0,
    2.46277e-06,
    2.95532e-05,
    0.000104668,
    0.000401431,
    0.00130034,
    0.00342202,
    0.00818132,
    0.0175534,
    0.035784,
    0.0650836,
    0.112232,
    0.178699,
    0.268934,
    0.380868,
    0.507505,
    0.640922,
    0.768551,
    0.877829,
    0.958624,
    0.99939,
    1
  };

  static double weight_23[25] = {
    0,
    1.20628e-06,
    1.20628e-06,
    2.41255e-06,
    1.20628e-05,
    6.39326e-05,
    0.000252112,
    0.000862487,
    0.00244995,
    0.00616527,
    0.0140821,
    0.0293342,
    0.0564501,
    0.100602,
    0.164479,
    0.252659,
    0.36268,
    0.491427,
    0.627979,
    0.75918,
    0.873185,
    0.957934,
    0.999381,
    1,
    0.957738
  };

  static double weight_22[25] = {
    0,
    0,
    0,
    5.88636e-06,
    3.0609e-05,
    0.000143627,
    0.000561558,
    0.00173059,
    0.00460078,
    0.0110616,
    0.0238974,
    0.0475406,
    0.0875077,
    0.148682,
    0.235752,
    0.343591,
    0.473146,
    0.611897,
    0.748345,
    0.865978,
    0.953199,
    0.997848,
    1,
    0.954245,
    0.873688
  };

  static double weight_21[25] = {
    0,
    0,
    1.15381e-06,
    8.07665e-06,
    7.1536e-05,
    0.000280375,
    0.00107189,
    0.00327104,
    0.00809396,
    0.0190978,
    0.0401894,
    0.0761028,
    0.13472,
    0.216315,
    0.324649,
    0.455125,
    0.598241,
    0.739215,
    0.861866,
    0.953911,
    0.998918,
    1,
    0.956683,
    0.872272,
    0.76399
  };
 
 
  static double weight_20[25] = {
    0,
    0,
    1.12532e-06,
    2.58822e-05,
    0.000145166,
    0.000633552,
    0.00215048,
    0.00592816,
    0.0145605,
    0.0328367,
    0.0652649,
    0.11893,
    0.19803,
    0.305525,
    0.436588,
    0.581566,
    0.727048,
    0.8534,
    0.949419,
    0.999785,
    1,
    0.953008,
    0.865689,
    0.753288,
    0.62765
  }; 
  static double weight_19[25] = {
    0,
    0,
    1.20714e-05,
    5.92596e-05,
    0.000364337,
    0.00124994,
    0.00403953,
    0.0108149,
    0.025824,
    0.0544969,
    0.103567,
    0.17936,
    0.283532,
    0.416091,
    0.562078,
    0.714714,
    0.846523,
    0.947875,
    1,
    0.999448,
    0.951404,
    0.859717,
    0.742319,
    0.613601,
    0.48552
  };

  static double weight_18[25] = {
    0,
    3.20101e-06,
    2.88091e-05,
    0.000164319,
    0.000719161,
    0.00250106,
    0.00773685,
    0.0197513,
    0.0443693,
    0.0885998,
    0.159891,
    0.262607,
    0.392327,
    0.543125,
    0.69924,
    0.837474,
    0.943486,
    0.998029,
    1,
    0.945937,
    0.851807,
    0.729309,
    0.596332,
    0.467818,
    0.350434
  };

 
  static double weight_17[25] = {
    1.03634e-06,
    7.25437e-06,
    4.97443e-05,
    0.000340956,
    0.00148715,
    0.00501485,
    0.0143067,
    0.034679,
    0.0742009,
    0.140287,
    0.238288,
    0.369416,
    0.521637,
    0.682368,
    0.828634,
    0.939655,
    1,
    0.996829,
    0.94062,
    0.841575,
    0.716664,
    0.582053,
    0.449595,
    0.331336,
    0.234332
  };

 
  static double weight_16[25] = {
    4.03159e-06,
    2.41895e-05,
    0.000141106,
    0.00081942,
    0.00314565,
    0.00990662,
    0.026293,
    0.0603881,
    0.120973,
    0.214532,
    0.343708,
    0.501141,
    0.665978,
    0.820107,
    0.938149,
    1,
    0.99941,
    0.940768,
    0.837813,
    0.703086,
    0.564023,
    0.42928,
    0.312515,
    0.216251,
    0.14561
  };
 
 
  static double weight_15[25] = {
    9.76084e-07,
    5.07564e-05,
    0.000303562,
    0.00174036,
    0.00617959,
    0.0188579,
    0.047465,
    0.101656,
    0.189492,
    0.315673,
    0.474383,
    0.646828,
    0.809462,
    0.934107,
    0.998874,
    1,
    0.936163,
    0.827473,
    0.689675,
    0.544384,
    0.40907,
    0.290648,
    0.198861,
    0.12951,
    0.0808051
  };
 
 
  static double weight_14[25] = {
    1.13288e-05,
    0.000124617,
    0.000753365,
    0.00345056,
    0.0123909,
    0.0352712,
    0.0825463,
    0.16413,
    0.287213,
    0.44615,
    0.625826,
    0.796365,
    0.930624,
    0.999958,
    1,
    0.934414,
    0.816456,
    0.672939,
    0.523033,
    0.386068,
    0.269824,
    0.180342,
    0.114669,
    0.0698288,
    0.0406496
  };

 
  static double weight_13[25] = {
    2.54296e-05,
    0.000261561,
    0.00167018,
    0.00748083,
    0.0241308,
    0.0636801,
    0.138222,
    0.255814,
    0.414275,
    0.600244,
    0.779958,
    0.92256,
    0.999155,
    1,
    0.927126,
    0.804504,
    0.651803,
    0.497534,
    0.35976,
    0.245834,
    0.160904,
    0.0991589,
    0.0585434,
    0.0332437,
    0.0180159
  };

  static double weight_12[25] = {
    5.85742e-05,
    0.000627706,
    0.00386677,
    0.0154068,
    0.0465892,
    0.111683,
    0.222487,
    0.381677,
    0.5719,
    0.765001,
    0.915916,
    1,
    0.999717,
    0.921443,
    0.791958,
    0.632344,
    0.475195,
    0.334982,
    0.223666,
    0.141781,
    0.0851538,
    0.048433,
    0.0263287,
    0.0133969,
    0.00696683
  };

 
  static double weight_11[25] = {
    0.00015238,
    0.00156064,
    0.00846044,
    0.0310939,
    0.0856225,
    0.187589,
    0.343579,
    0.541892,
    0.74224,
    0.909269,
    0.998711,
    1,
    0.916889,
    0.77485,
    0.608819,
    0.447016,
    0.307375,
    0.198444,
    0.121208,
    0.070222,
    0.0386492,
    0.0201108,
    0.0100922,
    0.00484937,
    0.00222458
  };

  static double weight_10[25] = {
    0.000393044,
    0.00367001,
    0.0179474,
    0.060389,
    0.151477,
    0.302077,
    0.503113,
    0.720373,
    0.899568,
    1,
    0.997739,
    0.909409,
    0.75728,
    0.582031,
    0.415322,
    0.277663,
    0.174147,
    0.102154,
    0.0566719,
    0.0298642,
    0.0147751,
    0.00710995,
    0.00319628,
    0.00140601,
    0.000568796
  };

 
  static double weight_9[25] = {
    0.00093396,
    0.00854448,
    0.0380306,
    0.113181,
    0.256614,
    0.460894,
    0.690242,
    0.888781,
    1,
    0.998756,
    0.899872,
    0.735642,
    0.552532,
    0.382726,
    0.246114,
    0.147497,
    0.0825541,
    0.0441199,
    0.0218157,
    0.0103578,
    0.00462959,
    0.0019142,
    0.000771598,
    0.000295893,
    0.000111529
  };

 
  static double weight_8[25] = {
    0.00240233,
    0.0192688,
    0.0768653,
    0.205008,
    0.410958,
    0.65758,
    0.875657,
    0.999886,
    1,
    0.889476,
    0.711446,
    0.517781,
    0.345774,
    0.212028,
    0.121208,
    0.0644629,
    0.0324928,
    0.0152492,
    0.00673527,
    0.0028547,
    0.00117213,
    0.000440177,
    0.000168471,
    5.80689e-05,
    1.93563e-05
  };

  static double weight_7[25] = {
    0.00617233,
    0.0428714,
    0.150018,
    0.350317,
    0.612535,
    0.856525,
    0.999923,
    1,
    0.87544,
    0.679383,
    0.478345,
    0.303378,
    0.176923,
    0.0950103,
    0.0476253,
    0.0222211,
    0.00972738,
    0.00392962,
    0.0015258,
    0.000559168,
    0.000183928,
    6.77983e-05,
    1.67818e-05,
    7.38398e-06,
    6.71271e-07
  };
 
  static double weight_6[25] = {
    0.0154465,
    0.0923472,
    0.277322,
    0.55552,
    0.833099,
    0.999035,
    1,
    0.855183,
    0.641976,
    0.428277,
    0.256804,
    0.139798,
    0.0700072,
    0.0321586,
    0.0137971,
    0.00544756,
    0.00202316,
    0.000766228,
    0.000259348,
    8.45836e-05,
    1.80362e-05,
    8.70713e-06,
    3.73163e-06,
    6.21938e-07,
    0
  };
 
 
  static double weight_5[25] = {
    0.0382845,
    0.191122,
    0.478782,
    0.797314,
    1,
    0.997148,
    0.831144,
    0.59461,
    0.371293,
    0.205903,
    0.103102,
    0.0471424,
    0.0194997,
    0.00749415,
    0.00273709,
    0.000879189,
    0.000286049,
    0.000102364,
    1.70606e-05,
    3.98081e-06,
    2.27475e-06,
    0,
    0,
    0,
    0
  };
 
 
  static double weight_4[25] = {
    0.0941305,
    0.373824,
    0.750094,
    1,
    0.997698,
    0.800956,
    0.532306,
    0.304597,
    0.152207,
    0.0676275,
    0.0270646,
    0.00975365,
    0.00326077,
    0.00101071,
    0.000301781,
    7.41664e-05,
    1.58563e-05,
    3.58045e-06,
    1.02299e-06,
    0,
    5.11493e-07,
    0,
    0,
    0,
    0
  };
 
 
  static double weight_3[25] = {
    0.222714,
    0.667015,
    1,
    0.999208,
    0.750609,
    0.449854,
    0.224968,
    0.0965185,
    0.0361225,
    0.012084,
    0.00359618,
    0.000977166,
    0.000239269,
    6.29422e-05,
    1.16064e-05,
    1.78559e-06,
    0,
    4.46398e-07,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };
 
  static double weight_2[25] = {
    0.499541,
    0.999607,
    1,
    0.666607,
    0.333301,
    0.13279,
    0.0441871,
    0.0127455,
    0.00318434,
    0.00071752,
    0.000132204,
    2.69578e-05,
    5.16999e-06,
    2.21571e-06,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };
 
  static double weight_1[25] = {
    0.999165,
    1,
    0.499996,
    0.166868,
    0.0414266,
    0.00831053,
    0.00137472,
    0.000198911,
    2.66302e-05,
    2.44563e-06,
    2.71737e-07,
    2.71737e-07,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };
 
  static double weight_0[25] = {
    1,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0,
    0
  };

  //WeightOOTPU_ = {0};

  double* WeightPtr = 0;

  for(int iint = 0; iint<25; ++iint){
    if(iint ==0) WeightPtr = weight_0;
    if(iint ==1) WeightPtr = weight_1;
    if(iint ==2) WeightPtr = weight_2;
    if(iint ==3) WeightPtr = weight_3;
    if(iint ==4) WeightPtr = weight_4;
    if(iint ==5) WeightPtr = weight_5;
    if(iint ==6) WeightPtr = weight_6;
    if(iint ==7) WeightPtr = weight_7;
    if(iint ==8) WeightPtr = weight_8;
    if(iint ==9) WeightPtr = weight_9;
    if(iint ==10) WeightPtr = weight_10;
    if(iint ==11) WeightPtr = weight_11;
    if(iint ==12) WeightPtr = weight_12;
    if(iint ==13) WeightPtr = weight_13;
    if(iint ==14) WeightPtr = weight_14;
    if(iint ==15) WeightPtr = weight_15;
    if(iint ==16) WeightPtr = weight_16;
    if(iint ==17) WeightPtr = weight_17;
    if(iint ==18) WeightPtr = weight_18;
    if(iint ==19) WeightPtr = weight_19;
    if(iint ==20) WeightPtr = weight_20;
    if(iint ==21) WeightPtr = weight_21;
    if(iint ==22) WeightPtr = weight_22;
    if(iint ==23) WeightPtr = weight_23;
    if(iint ==24) WeightPtr = weight_24;

    for(int ibin = 0; ibin<25; ++ibin){
      WeightOOTPU_[iint][ibin] = *(WeightPtr+ibin);
    }
  }

}



void LumiReWeighting_setpileupWeights(TH1D *hdata, TH1D *hmc){
  TH1D *num = (TH1D*)hdata->Clone("datawt");

  cout<<"dd" <<endl;
  weights_ = (TH1D*)hdata->Clone("normwt");
  double total = weights_->Integral();
  cout<<"dd1" <<endl;

  num ->Scale( 1.0/total);
  TH1D *den = (TH1D*)hmc->Clone("mcwt");
  double total1 = den->Integral();
  den ->Scale(1.0/total1);
  ////weights_ ->Divide( den);
  int nbin = weights_->GetNbinsX();
  
  double sum_wt = 0; 
  for(int b=1; b<= nbin; b++){
    if( den->GetBinContent(b) > 0){
      weights_ ->SetBinContent(b, num->GetBinContent(b)/ den->GetBinContent(b));
    }else{
      weights_ ->SetBinContent(b, 0);
    }
    sum_wt += weights_ ->GetBinContent(b) *  den->GetBinContent(b); 
    cout<<"intimepileweight "<< b-1 <<" "<<  weights_ ->GetBinContent(b)<<endl;
  }
  cout<<"sum wt "<< sum_wt <<" "<< weights_->GetMean()<<endl;  
  
  LumiReWeighting_weightOOT_init();
  
}



///Now use this version
vector<double>  LumiReWeighting_getWeights(TH1D *hdata, TH1D *hmc){
  TH1D *num = (TH1D*)hdata->Clone("datawt");
  double total = num->Integral();
  num ->Scale( 1.0/total);
  
  TH1D *den = (TH1D*)hmc->Clone("mcwt");
  double total1 = den->Integral();
  den ->Scale(1.0/total1);
  
  int nbindata = num->GetNbinsX();
  int nbinmc = den->GetNbinsX();
  
  if( nbinmc != nbindata ){
    cout<<" waring! nbinmc != nbindata " << nbinmc <<" "<< nbindata <<endl; 
  }
  
  
  vector<double> wts; 
  int nbin = nbinmc < nbindata ? nbinmc: nbindata; 
  
  for(int b=1; b<= nbin; b++){
    if( den->GetBinContent(b) > 0){
      wts.push_back( num->GetBinContent(b)/ den->GetBinContent(b));
    }else{
      wts.push_back(0);
    }
    cout<<" LumiReWeighting_getWeights " << wts[b-1] <<endl; 
    
  }
  if( nbinmc > nbindata){
    for(int b= nbindata; b < nbinmc; b++){
      wts.push_back(0);
    }
  }
  
  //final check
  if( int(wts.size()) != nbinmc){
    cout<<" problem of seeting pileup wts for MC " <<  wts.size() <<" "<< nbinmc <<endl; 
    exit(1);
  }

  return wts;
  
  
  ///LumiReWeighting_weightOOT_init();
  
}





double LumiReWeighting_weightINT_v1(vector<double> weight){
  
  if( weight.size() <=0){
    cout<<" LumiReWeighting_weightINT_v1 bad weight vector " << endl; 
    exit(1);
  }
  
  int npv = -1; 
  for(unsigned int j =0; j < pileupBunchX->size(); j++){
    int BX = pileupBunchX->at(j);
    if(BX==0){
      npv = pileupNInteraction->at(j);
      break; 
    }
  }
  if(npv < 0) {
    cout << " no in-time beam crossing found\n! " <<endl; 
    return 0;
  }
  if(npv > int(weight.size()) ) {
    cout<<"warning! LumiReWeighting_weightINT_v1 " << npv <<" "<< weight.size() <<endl; 
    return 0; 
  }
  return weight[npv];
  
}




double LumiReWeighting_weightINT_v2(int npv, vector<double> weight){
  
  if( weight.size() <=0){
    cout<<" LumiReWeighting_weightINT_v2 bad weight vector " << endl; 
    exit(1);
  }
  
//   int npv = -1; 
//   for(unsigned int j =0; j < pileupBunchX->size(); j++){
//     int BX = pileupBunchX->at(j);
//     if(BX==0){
//       npv = pileupNInteraction->at(j);
//       break; 
//     }
//   }
  if(npv < 0) {
    cout << " no in-time beam crossing found\n! " <<endl; 
    return 0;
  }
  if(npv > int(weight.size()) ) {
    cout<<"warning! LumiReWeighting_weightINT_v2 " << npv <<" "<< weight.size() <<endl; 
    return 0; 
  }
  return weight[npv];
  
}





double LumiReWeighting_weightINT(){
  
  int npv = -1; 
  for(unsigned int j =0; j < pileupBunchX->size(); j++){
    int BX = pileupBunchX->at(j);
    if(BX==0){
      npv = pileupNInteraction->at(j);
      break; 
    }
  }
  if(npv < 0) cout << " no in-time beam crossing found\n! " <<endl; 

  if( weights_ == NULL){
    cout<<" weights histogram not filled"<<endl;
    return 1; 
  }
  int bin = weights_->GetXaxis()->FindBin( npv );
  double inTimeWeight = weights_->GetBinContent( bin );
  return inTimeWeight; 
  
}


double LumiReWeighting_weightOOT(string dataset) {
    
  if( dataset.find("SIM") == string::npos){
    cout<<" weight not on MC ? " << dataset.c_str() <<endl; 
    exit(1);
  }
  
  
  static double Correct_Weights2011[25] = { // residual correction to match lumi spectrum
    5.30031,
    2.07903,
    1.40729,
    1.27687,
    1.0702,
    0.902094,
    0.902345,
    0.931449,
    0.78202,
    0.824686,
    0.837735,
    0.910261,
    1.01394,
    1.1599,
    1.12778,
    1.58423,
    1.78868,
    1.58296,
    2.3291,
    3.86641,
    0,
    0,
    0,
    0,
    0
  };                        
  
  
  int npv = -1;
  int npv50ns = -1;
  for(unsigned int j =0; j < pileupBunchX->size(); j++){
    int BX = pileupBunchX->at(j);
    if(BX==0){
      npv = pileupNInteraction->at(j);
    }
    if(BX == 1) { 
      npv50ns = pileupNInteraction->at(j);
    }
    
  }
  
  // Note: for the "uncorrelated" out-of-time pileup, reweighting is only done on the 50ns
  // "late" bunch (BX=+1), since that is basically the only one that matters in terms of 
  // energy deposition.  

  if(npv < 0) {
    std::cerr << " no in-time beam crossing found\n! " ;
    std::cerr << " Returning event weight=0\n! ";
    return 0.;
  }
 //  if(npv50ns < 0) {
//     std::cerr << " no out-of-time beam crossing found\n! " ;
//     std::cerr << " Returning event weight=0\n! ";
//     return 0.;
//   }
  
  int bin = weights_->GetXaxis()->FindBin( npv );

  double inTimeWeight = weights_->GetBinContent( bin );
  
  double TotalWeight = 1.0;
  
  
  if( dataset.find("S3") != string::npos){    //if(Reweight_4_2_2p2_) {
    TotalWeight = inTimeWeight * WeightOOTPU_[bin-1][npv50ns] * Correct_Weights2011[bin-1];
  }
  else {
    TotalWeight = inTimeWeight;
  }
  
//   if( TotalWeight < 1E-5 && nWarning < 10 ){
//     cout<<"warnign 0 weight " << TotalWeight <<" "<< inTimeWeight <<" "<< bin <<" "<< npv<<" "<< npv50ns <<" "<< WeightOOTPU_[bin-1][npv50ns] <<" "<<Correct_Weights2011[bin-1] <<endl; 
//     nWarning ++; 
//   }
  
  
  return TotalWeight;
  
}




double LumiReWeighting_weightOOT_givenWeightVector(string dataset, vector<double> puwts) {
  
  if( dataset.find("SIM") == string::npos){
    cout<<" weight not on MC ? " << dataset.c_str() <<endl; 
    exit(1);
  }
  
  
  static double Correct_Weights2011[25] = { // residual correction to match lumi spectrum
    5.30031,
    2.07903,
    1.40729,
    1.27687,
    1.0702,
    0.902094,
    0.902345,
    0.931449,
    0.78202,
    0.824686,
    0.837735,
    0.910261,
    1.01394,
    1.1599,
    1.12778,
    1.58423,
    1.78868,
    1.58296,
    2.3291,
    3.86641,
    0,
    0,
    0,
    0,
    0
  };                        
  
  
  int npv = -1;
  int npv50ns = -1;
  for(unsigned int j =0; j < pileupBunchX->size(); j++){
    int BX = pileupBunchX->at(j);
    if(BX==0){
      npv = pileupNInteraction->at(j);
    }
    if(BX == 1) { 
      npv50ns = pileupNInteraction->at(j);
    }
    
  }
  
  // Note: for the "uncorrelated" out-of-time pileup, reweighting is only done on the 50ns
  // "late" bunch (BX=+1), since that is basically the only one that matters in terms of 
  // energy deposition.  

  if(npv < 0) {
    std::cerr << " no in-time beam crossing found\n! " ;
    std::cerr << " Returning event weight=0\n! ";
    return 0.;
  }
 //  if(npv50ns < 0) {
//     std::cerr << " no out-of-time beam crossing found\n! " ;
//     std::cerr << " Returning event weight=0\n! ";
//     return 0.;
//   }
  
  
  
  double inTimeWeight = puwts[npv];
  
  double TotalWeight = 1.0;
  
  
  
  if( dataset.find("S3") != string::npos){    //if(Reweight_4_2_2p2_) {
    TotalWeight = inTimeWeight * WeightOOTPU_[npv][npv50ns] * Correct_Weights2011[npv];
  }
  else {
    TotalWeight = inTimeWeight;
  }
  
//   if( TotalWeight < 1E-5 && nWarning < 10 ){
//     cout<<"warnign 0 weight " << TotalWeight <<" "<< inTimeWeight <<" "<< bin <<" "<< npv<<" "<< npv50ns <<" "<< WeightOOTPU_[bin-1][npv50ns] <<" "<<Correct_Weights2011[bin-1] <<endl; 
//     nWarning ++; 
//   }
  
  
  return TotalWeight;
  
}

