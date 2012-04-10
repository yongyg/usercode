

TF1* fitgaus(TH1F *hist,float low, float high,float res[]){
  
  //  Float_t rms = hist->GetRMS();
  // Float_t maxbin = hist->GetMaximumBin() * hist->GetBinWidth(1) +hist->GetXaxis()->GetXmin(); 
  // gStyle->SetOptFit(111111);
  
  TF1 *f1 = new TF1("f1","gaus",low,high);
  
  
  
  hist->Fit("f1","RQ");
  Float_t sigma = f1->GetParameter(2); 
  Float_t mean  = f1->GetParameter(1); 
  float sigerr = f1->GetParError(2);
  res[0] = mean; 
  res[1] = f1->GetParError(1); 
  res[2] = sigma; 
  res[3] = sigerr; 
  res[4] = mean; 
  res[5] = f1->GetParError(1);
  res[6] = sigma; 
  
  
  return f1; 

}



TF1* fitgausD(TH1D *hist,float low, float high,float res[]){
  
  //  Float_t rms = hist->GetRMS();
  // Float_t maxbin = hist->GetMaximumBin() * hist->GetBinWidth(1) +hist->GetXaxis()->GetXmin(); 
  // gStyle->SetOptFit(111111);
  
  TF1 *f1 = new TF1("f1","gaus",low,high);

  


  hist->Fit("f1","RQ");
  Float_t sigma = f1->GetParameter(2); 
  Float_t mean  = f1->GetParameter(1); 
  float sigerr = f1->GetParError(2);
  res[0] = mean; 
  res[1] = f1->GetParError(1); 
  res[3] = sigma; 
  res[2] = sigerr; 
  res[4] = mean; 
  res[5] = f1->GetParError(1);
  res[6] = sigma; 
  
  
  return f1; 

}




TF1* fitgauswind(TH1F *hist, Float_t mins, Float_t maxs){
  
  float res[20];
  
  
  Float_t rms = hist->GetRMS();
  Float_t maxbin = hist->GetMaximumBin() * hist->GetBinWidth(1) +hist->GetXaxis()->GetXmin(); 
  ///gStyle->SetOptFit(111111);
  
  TF1 *f1 = new TF1("f1","gaus",maxbin-2*rms, maxbin+2*rms);
  hist->Fit("f1","R0Q");
  Float_t sigma = f1->GetParameter(2); 
  Float_t mean  = f1->GetParameter(1); 
  //  float sigerr = f1->GetParError(2);
  
  for( int ll=0; ll<10; ll++){
    f1 = fitgaus(hist,mean-mins*sigma, mean+maxs*sigma,res);
    if(fabs(sigma-res[6])>0.001){
      sigma = res[6];
      mean = res[4];
    }
    else
      break; 
  }
  
  return f1; 

}


void fitgauswindRefit(TH1F *hist,float res[]){
    
  Float_t rms = hist->GetRMS();
  Float_t maxbin = hist->GetMaximumBin() * hist->GetBinWidth(1) +hist->GetXaxis()->GetXmin(); 
  ///gStyle->SetOptFit(111111);
  
  TF1 *f1 = new TF1("f1","gaus",maxbin-2*rms, maxbin+2*rms);
  hist->Fit("f1","RQ");
  
  Float_t sigma = f1->GetParameter(2); 
  Float_t mean  = f1->GetParameter(1); 
  
  TF1 *f2 = new TF1("f2","gaus",mean-5*sigma, mean+5*sigma); 
  for(int j=0;j<3; j++){
    f2->SetParameter(j,f1->GetParameter(j));
  }
  hist->Fit("f2","RQ");
  
  sigma = f2->GetParameter(2); 
  mean  = f2->GetParameter(1); 
  float sigerr = f2->GetParError(2);
  
  res[0] = mean; 
  res[1] = f2->GetParError(1); 
  res[2] = sigma; 
  res[3] = sigerr; 
  res[4] = mean; 
  res[5] = f2->GetParError(1);
  res[6] = sigma; 
    
}




void fitgauswind2(TH1F *hist, Float_t mins, Float_t maxs,float res[]){
  
  
  Float_t rms = hist->GetRMS();
  Float_t maxbin = hist->GetMaximumBin() * hist->GetBinWidth(1) +hist->GetXaxis()->GetXmin(); 
  ///gStyle->SetOptFit(111111);
  
  TF1 *f1 = new TF1("f1","gaus",maxbin-2*rms, maxbin+2*rms);
  hist->Fit("f1","R0Q");
  Float_t sigma = f1->GetParameter(2); 
  Float_t mean  = f1->GetParameter(1); 
  //  float sigerr = f1->GetParError(2);
  
  for( int ll=0; ll<10; ll++){
    f1 = fitgaus(hist,mean-mins*sigma, mean+maxs*sigma,res);
    if(fabs(sigma-res[6])>0.001){
      sigma = res[6];
      mean = res[4];
    }
    else
      break; 
  }
  
  //  return f1; 

}



void fitgauswindD(TH1D *hist, Float_t mins, Float_t maxs,float res[]){
  
  
  Float_t rms = hist->GetRMS();
  Float_t maxbin = hist->GetMaximumBin() * hist->GetBinWidth(1) +hist->GetXaxis()->GetXmin(); 
  ///gStyle->SetOptFit(111111);
  
  hist->GetXaxis()->SetRangeUser(maxbin-3*rms, maxbin+3*rms);
  
  // cout<<"maxbin: "<< maxbin <<" "<< maxbin-2*rms <<" "<< maxbin+2*rms <<endl; 
  

  TF1 *f1 = new TF1("f1","gaus",maxbin-2*rms, maxbin+2*rms);
  hist->Fit("f1","R0Q");
  Float_t sigma = f1->GetParameter(2); 
  Float_t mean  = f1->GetParameter(1); 
  //  float sigerr = f1->GetParError(2);
  
  //cout<<"maxbin_iter: "<< maxbin <<" "<< maxbin-2*rms <<" "<< maxbin+2*rms <<endl; 

  for( int ll=0; ll<10; ll++){

    // cout<<"ll: "<<ll <<" "<< mean <<" "<< sigma <<endl; 

    f1 = fitgausD(hist,mean-mins*sigma, mean+maxs*sigma,res);
    if(fabs(sigma-res[6])>0.001){
      sigma = res[6];
      mean = res[4];
    }
    else
      break; 
  }
  
  //  return f1; 
  
}



TF1* fitgausI(TH1I *hist,float low, float high,float res[]){
  
  //  Float_t rms = hist->GetRMS();
  // Float_t maxbin = hist->GetMaximumBin() * hist->GetBinWidth(1) +hist->GetXaxis()->GetXmin(); 
  // gStyle->SetOptFit(111111);
  
  TF1 *f1 = new TF1("f1","gaus",low,high);
  hist->Fit("f1","RQ");
  Float_t sigma = f1->GetParameter(2); 
  Float_t mean  = f1->GetParameter(1); 
  float sigerr = f1->GetParError(2);
  res[0] = mean; 
  res[1] = f1->GetParError(1); 
  res[3] = sigma; 
  res[2] = sigerr; 
  res[4] = mean; 
  res[5] = f1->GetParError(1);
  res[6] = sigma; 
  
  
  return f1; 

}

TF1* fitgauswindI(TH1I *hist, Float_t mins, Float_t maxs){
  
  float res[20];
  
  
  Float_t rms = hist->GetRMS();
  Float_t maxbin = hist->GetMaximumBin() * hist->GetBinWidth(1) +hist->GetXaxis()->GetXmin(); 
  ///gStyle->SetOptFit(111111);
  
  TF1 *f1 = new TF1("f1","gaus",maxbin-2*rms, maxbin+2*rms);
  hist->Fit("f1","R0Q");
  Float_t sigma = f1->GetParameter(2); 
  Float_t mean  = f1->GetParameter(1); 
  //  float sigerr = f1->GetParError(2);
  
  for( int ll=0; ll<10; ll++){
    f1 = fitgausI(hist,mean-mins*sigma, mean+maxs*sigma,res);
    if(fabs(sigma-res[6])>0.001){
      sigma = res[6];
      mean = res[4];
    }
    else
      break; 
  }
  
  return f1; 

}

