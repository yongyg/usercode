#include "TH1.h"
#include "TF1.h"
#include "TLegend.h"
#include "TCanvas.h"


#include "TStyle.h"

#include "TLatex.h"

#include "TGaxis.h"


// nice peak  .x pi0_mfitpeak.C ("../pi0anal/pi0ana_punorm.root","minv_spb",0)
// pileup nice .x pi0_mfitpeak.C ("../pi0calib/pi0calibwPU.root","minv_spb",0) 

//void pi0_mfitpeak(char *FileName, char *HistName, Int_t etapi0flag) 
void pi0_mfitpeak(TH1F *mh1, Int_t etapi0flag, float xmin, float xmax, int npol,float res[],int posFlag, const char *dirName, const char *histName, float text_x, float text_y, const char *texName) 
  
{
  TGaxis::SetMaxDigits(3);


  // TVirtualFitter::SetDefaultFitter("Minuit");
  


  // This script attempts to fit any pi0 peak so that the freaking fit function would converge
  // currently background is fitted to a pol4 function and the peak by a gaussian;
  // results are not very nice
  // usage  .x pi0_mfitpeak.C++ ("pi0calib.root","minv_spb")
  // or eg.  .x pi0_mfitpeak.C ("../pi0anal/pi0ana_punorm.root","minv_spb",0)

  gROOT->Reset();
  //  gStyle->SetOptFit();
  //  gStyle->SetOptFit(0);
  //  gStyle->SetOptStat(0);
  //  gStyle->SetOptTitle(0);
Bool_t NOTE=1;
  if(NOTE) gStyle->SetCanvasBorderMode(0);

  gStyle->SetPadTopMargin(0.08);
  gStyle->SetPadBottomMargin(0.12);
  gStyle->SetPadLeftMargin(0.15);
  gStyle->SetPadRightMargin(0.08);

  mh1->GetXaxis()->SetRangeUser(xmin,xmax);
  
  Int_t highx=500;
  TCanvas *c2 = new TCanvas("c2","",200,10,highx,500);


    // cout<<FileName<<" "<<HistName<<endl;

  //     TFile f(FileName);
  //   TH1F *mh1 = (TH1F*) f.Get(HistName);
      mh1->SetMarkerStyle(20);
      mh1->SetMarkerSize(1.);
      mh1->SetStats(0); // 1/0 to set the stat box

   mh1->GetXaxis()->SetTitle("Invariant Mass of Photon Pairs (GeV/c^{2})");


   float binwidth = mh1->GetBinWidth(1);
   
   char *ytitle = new char[100];

   sprintf(ytitle,"Photon Pairs / %4.3f GeV/c^{2}",binwidth);
   

   mh1->GetYaxis()->SetTitle(ytitle);


   mh1->GetXaxis()->SetTitleSize(0.055);
   mh1->GetYaxis()->SetTitleSize(0.055);
   mh1->GetXaxis()->SetLabelSize(0.045);
   mh1->GetYaxis()->SetLabelSize(0.045);
   mh1->GetXaxis()->SetTitleOffset(0.90);
   mh1->GetXaxis()->CenterTitle();
   mh1->GetYaxis()->SetTitleOffset(1.32);
   

   // First work with the histogram and find the peak and fit ranges
   TAxis *xaxis = mh1->GetXaxis();
   Float_t binsiz= xaxis->GetBinCenter(3) - xaxis->GetBinCenter(2);
   Int_t nbins = xaxis->GetNbins(); 
   Float_t nevtperbin0[10000];
   Float_t errorbin0[10000];
   Float_t nevttot;
   Float_t maxbin=0; Int_t nmaxbin=0, nminbord=0, nmaxbord=nbins;
   
   for (Int_t nn=1; nn <= nbins; nn++)
     {
       nevtperbin0[nn] = mh1->GetBinContent(nn); 
       if(nevtperbin0[nn] > maxbin) { maxbin=nevtperbin0[nn]; nmaxbin=nn; }
       errorbin0[nn] = mh1->GetBinError(nn); 
       nevttot+=nevtperbin0[nn];
       if(nevtperbin0[nn] > 0 && nminbord == 0) nminbord=nn; 
       if(nevtperbin0[nn] == 0 && (nn > nminbord +10) && nmaxbord==0 && nminbord > 0) nmaxbord=nn; 
     }
   cout<<"Minbordl "<<nminbord<<" with events: "<<nevtperbin0[nminbord]<<endl;
   cout<<"Maxbordl "<<nmaxbord<<" with events: "<<nevtperbin0[nmaxbord]<<endl;
   nminbord+=0;
   nmaxbord-=0;
   Int_t nmin0=nminbord;
   while(nevtperbin0[nminbord] < nevtperbin0[nmaxbin]*0.025) nminbord++;
   while(nevtperbin0[nmaxbord] < nevtperbin0[nmaxbin]*0.025) nmaxbord--;
   // the above was just to get the info and low/high bins

   // Set the fit range ! This is for total fit !	 
   Float_t fitl=xmin;
     float fith=xmax;
   //     Float_t fitl=0.07, fith=0.2;// this works better for pileup
   //         Float_t fitl=0.08, fith=0.18;// this works even better for pileup
     //  if(etapi0flag == 1)
     // {
     //  fitl=0.35; fith=0.75;
     //}
     

     //   if(fitl < xaxis->GetBinCenter(nmin0)) fitl = xaxis->GetBinCenter(nmin0);
   //if(fith > xaxis->GetBinCenter(nmaxbord)) fith = xaxis->GetBinCenter(nmaxbord);
 
   


   cout<<" fit range "<<fitl<<" -- "<<fith<<endl;

   cout <<"Bin size "<<binsiz<<endl;
   cout<<"Total events "<<nevttot<<endl;
   cout<<"MaxBin "<<nmaxbin<<" with events: "<<nevtperbin0[nmaxbin]<<endl;
   cout<<"Minbord "<<nminbord<<" with events: "<<nevtperbin0[nminbord]<<endl;
   cout<<"Maxbord "<<nmaxbord<<" with events: "<<nevtperbin0[nmaxbord]<<endl;
      mh1->DrawCopy("sep");

      Float_t lowgauss=0.135-4.*0.010;
      Float_t highgauss=0.135+4.*0.010;
      if(etapi0flag == 1)
	{
	  lowgauss=0.55-5.*0.025;
              highgauss=0.55+5.*0.025;
	}
      Int_t nlowgauss=Int_t((lowgauss-xaxis->GetBinCenter(1))/Float_t(binsiz)+0.5);
      Int_t nhighgauss=Int_t((highgauss-xaxis->GetBinCenter(1))/Float_t(binsiz)+0.5);
      cout <<xaxis->GetBinCenter(nlowgauss)<<" "<<xaxis->GetBinCenter(nhighgauss)<<endl;
   // now make the "background" histogram and fit it with p4
      Float_t lowvalgauss=nevtperbin0[nlowgauss];
      Float_t increm=(nevtperbin0[nhighgauss]-nevtperbin0[nlowgauss])/Float_t(nhighgauss-nlowgauss);
      TH1F *hbkg = (TH1F*)mh1->Clone();
      hbkg->SetName("bkg_clone");
      for (Int_t nn=nlowgauss; nn<=nhighgauss; nn++)
	{
	  hbkg->SetBinContent(nn,Float_t(lowvalgauss+(nn-nlowgauss)*increm));
	  hbkg->SetBinError(nn,sqrt(lowvalgauss+(nn-nlowgauss)*increm));
	}
      hbkg->DrawCopy("samesep");
      //      break;
      // Now define the "gaussian" histogram
      TH1F *hgauss = (TH1F*)mh1->Clone();
      hgauss->SetName("gauss_clone");
      hgauss->Sumw2();
      hgauss->Add(mh1,hbkg,1,-1); // if errors are independent Add needs to be used !
       for (Int_t nn=1; nn <= nbins; nn++)
	 {
	   if(hgauss->GetBinContent(nn) < 0.) hgauss->SetBinContent(nn,0.001*nevtperbin0[nmaxbin]);
	   hgauss->SetBinError(nn,sqrt(hgauss->GetBinContent(nn)));
	 }

   // Declare function with wich to fit
       TF1 *g1 = new TF1("g1","gaus",lowgauss,highgauss);
   hgauss->Fit(g1,"R0");
   hgauss->DrawCopy("sep");
   g1->Draw("same");
   //  break;

   char *polff = new char[20];

   sprintf(polff,"pol%d",npol);
   
   TF1 *p4bkg; 
   if(etapi0flag != 1)
     p4bkg   = new TF1("pm2",polff, xaxis->GetBinCenter(nminbord),xaxis->GetBinCenter(nmaxbord));
   else
     p4bkg   = new TF1("pm2",polff, 0.35,0.75);
   
   

   hbkg->Fit(p4bkg,"R0");
   hbkg->DrawCopy("sep");
   p4bkg->SetLineStyle(kDashed);
   p4bkg->Draw("same");
   // break;
   
   
   Double_t par[20],parf[20],errparf[20];
   g1->GetParameters(&par[0]);
   p4bkg->GetParameters(&par[3]);

   char *totff = new char[20];
   
   sprintf(totff,"gaus(0)+pol%d(3)",npol);
   
   
   TF1 *total = new TF1("total",totff,fitl,fith);
   TF1 *p4bkgfin   = new TF1("pm2",polff,fitl,fith);
   total->SetParameters(par);

   if(etapi0flag==0){
     total->SetParLimits(1,0.10,0.15);
     total->SetParLimits(2,0.135*0.06,0.135*0.3);
     
   }else{
     total->SetParLimits(1,0.35,0.65);
   }
   
   
   //  total->FixParameter(1,1.21340e-01); 
   // total->FixParameter(2,2.69780e-02);
   
   
   

   

   mh1->Fit(total,"R0");
   
   cout<<" yield.. "<< total->GetParameter(0) <<"+/- " << total->GetParError(0)<<endl;
   

     total->GetParameters(parf);


     for( Int_t nn=0; nn < 3+npol+1; nn++) errparf[nn]=total->GetParError(nn);
     g1->SetParameters(&parf[0]);
     p4bkgfin->SetParameters(&parf[3]);
     cout <<" Piz Mass = "<<parf[1]*1000.<<" +- "<<errparf[1]*1000.<<
       " Sigma ="<<parf[2]*1000.<<" +- "<<errparf[2]*1000.<<endl;
      cout << " Sigma Rel. = "<< parf[2]/parf[1]<<" +- "
	   << errparf[2]/parf[1] <<endl;

    Float_t int_min=parf[1]-3.*parf[2];
    Float_t int_max=parf[1]+3.*parf[2];
    Float_t sig_peak=g1->Integral(int_min,int_max)/binsiz;
    Float_t bkgd_peak=p4bkgfin->Integral(int_min,int_max)/binsiz;
    cout<<" In +-3. sigma window: Signal="<<sig_peak<<" Bkgd="<<bkgd_peak<<endl;
    Float_t SB=sig_peak/bkgd_peak; Float_t SBerr=SB*(sqrt(1./sig_peak+1./bkgd_peak));
    cout<<" Signal/Bkgd "<<SB<<" +- "<<SBerr<<endl;   

    int_min=parf[1]-2.*parf[2];
    int_max=parf[1]+2.*parf[2];
    sig_peak=g1->Integral(int_min,int_max)/binsiz;
    bkgd_peak=p4bkgfin->Integral(int_min,int_max)/binsiz;
    cout<<" In +-2.sigma window: Signal="<<sig_peak<<" Bkgd="<<bkgd_peak<<endl;
    SB=sig_peak/bkgd_peak; SBerr=SB*(sqrt(1./sig_peak+1./bkgd_peak));
    cout<<" Signal/Bkgd "<<SB<<" +- "<<SBerr<<endl;   
 			

    float S = sig_peak; 
    float B = bkgd_peak; 
    
    int_min=parf[1]-20.*parf[2];
    int_max=parf[1]+20.*parf[2];
    float S_all = g1->Integral(int_min,int_max)/binsiz;
    
    

    float test_sall =  parf[0] * parf[2] * sqrt(acos(-1)) / binsiz; 
    
    cout<<"signal all: "<< S_all << " "<< test_sall <<endl; 
    
    float Serr_all = errparf[0]/ parf[0] * S_all; 
    


    res[0] = S_all; 
    res[1] = Serr_all; 

    res[2] = parf[1]; 
    res[3] = errparf[1];
    
    res[4] = parf[2]; 
    res[5] = errparf[2];

    res[6] = SB; 
    res[7] = SBerr; 
    


    total->SetLineWidth(3);
    
   total->SetLineColor(kBlue);

   p4bkgfin->SetLineWidth(3);
   
   p4bkgfin->SetLineColor(kRed);
   mh1->DrawCopy("sep");

   //   total->SetRange(0.07,0.185);
   //    p4bkgfin->SetRange(0.07,0.185);
   total->Draw("same");
   p4bkgfin->SetLineStyle(kDashed);
   p4bkgfin->Draw("same");   

   TLatex l;
   
   //   l.SetTextSize(0.06);
   l.SetTextSize(0.05);
   
   
   l.SetTextColor(1);
   l.SetNDC();
   
  
   float sigma = parf[2]/ parf[1]*100; 
   float sigmaerr = errparf[2]/ parf[1]*100; 
   char *sigma_name = new char[50]; 
   if(sigmaerr>0.005)
     sprintf(sigma_name,"#sigma = %3.2f #pm %3.2f %% ",sigma,sigmaerr);
   else sprintf(sigma_name,"#sigma = %3.2f %% ",sigma);
   
   if(posFlag==1){
     l.DrawLatex(0.54,0.75,sigma_name);
   }else if( posFlag==2){
     l.DrawLatex(0.54,0.3,sigma_name);

   }
      
      //  sprintf(sigma_name,"S/B = %3.2f #pm %3.2f ",SB,SBerr);
      //
      sprintf(sigma_name,"S = %2.1f #pm %2.1f",S_all,Serr_all);
	      //
      
      // l.DrawLatex(0.5,0.5,sigma_name);
      
      sprintf(sigma_name,"M = %3.1f #pm %3.1f MeV",parf[1]*1000.,errparf[1]*1000);
      if(posFlag==1){
      l.DrawLatex(0.54,0.82,sigma_name);
      }else if( posFlag==2){
	l.DrawLatex(0.54,0.37,sigma_name);

      }
      ///l.DrawLatex(0.169,470.,"d)");


      c2->Modified();
   c2->Update();
   //   c2->SaveAs("nice_pi0.gif");

   
   
   if( text_x >0 && text_y >0){
    TLatex *   tex = new TLatex(text_x, text_y, texName);
    tex->SetNDC();
    tex->SetTextSize(0.06);
    tex->SetLineWidth(2);
    tex->Draw();
  }
  

   char *filename = new char[1000];

   
   sprintf(filename,"%s/%s.gif",dirName,histName);
   c2->Print(filename);
   sprintf(filename,"%s/%s.C",dirName,histName);
   c2->Print(filename);
   
   // .x pi0_mfitpeak.C ("pi0ana_508030.root","minv_spb",0)
   // .x pi0_mfitpeak.C ("pi0ana_npu.root","minv_bkg",0)  
   // .x pi0_mfitpeak.C ("pi0ana_punormv2.root","minv_spb",0)
}


///mean , the maximum bin.
float getMeanEastimateMaximumBin(TH1F *hh){
  
  
  float xmin = hh->GetXaxis()->GetXmin();
  float xmax = hh->GetXaxis()->GetXmax();
  
  float binwidth = hh->GetBinWidth(1);
  
  return xmin + binwidth*(hh->GetMaximumBin()+0.5);
  
  
  
}
