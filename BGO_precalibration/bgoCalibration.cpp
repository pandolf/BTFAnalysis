#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TFile.h"
#include "TTree.h"
#include "TF1.h"
#include "TCanvas.h"



bool useFermiDirac_ = true;


Double_t FermiDirac(Double_t *x, Double_t *par)
 // Fermi-Dirac distribution
{
  return 1./(exp((x[0]-par[0])/par[1])+1.);
}


Double_t FermiDiracPlusGaus(Double_t *x, Double_t *par)
{
  return par[0]/(exp((x[0]-par[1])/par[2])+1.) + par[3]*exp(-0.5*((x[0]-par[4])/par[5])*((x[0]-par[4])/par[5])) ;
}


float getSinglePeakPosition( const std::string& outdir, int hv, const std::string& suffix, int iChannel);


int main() {


  // first get pedestal:
  TFile* file_ped = TFile::Open("../run_precalib_BGO_pedestal_noSource.root");
  TTree* tree_ped = (TTree*)file_ped->Get("eventRawData");

  TH1D* h1_pedestal = new TH1D("ped", "", 200, 100., 200.);
  h1_pedestal->SetXTitle( "ADC Counts" );
  tree_ped->Project("ped", "adcData[0]");


  float pedestal = h1_pedestal->GetMean();

  int hv = 1700;
  std::string suffix = "front";

  std::string suffix_dir(suffix);
  if( useFermiDirac_ ) suffix_dir += "_FD";
  std::string outdir = Form("bgoCalibration_hv%d_%s", hv, suffix_dir.c_str() );
  std::string mkdir_command = "mkdir -p " + outdir;
  system(mkdir_command.c_str());

  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1_pedestal->Draw();
    
  c1->SaveAs( Form("%s/pedestal.eps", outdir.c_str()) );
  c1->SaveAs( Form("%s/pedestal.png", outdir.c_str()) );

  delete c1;


  std::string resultsFileName = Form("%s/results.txt", outdir.c_str());
  ofstream ofs(resultsFileName.c_str());

  ofs << Form("ped\t%.3f", pedestal) << std::endl;

  int nChannels=8;

  for( unsigned i=1; i<nChannels+1; ++i ) 
    ofs << i << "\t" << getSinglePeakPosition(outdir, hv, suffix, i) - pedestal << std::endl;

  ofs.close();

  std::cout << std::endl << "-> Saved results in file: " << resultsFileName << std::endl; 

  return 0;

}



float getSinglePeakPosition( const std::string& outdir, int hv, const std::string& suffix, int iChannel) {

  TFile* file = TFile::Open(Form("../run_precalib_BGO_ch%d_HV%d_%s.root", iChannel, hv, suffix.c_str()));

  TTree* tree = (TTree*)file->Get("eventRawData");

  float adcMax = 1000.;

  std::string name = Form("Channel %d", iChannel);
  TH1D* h1_spectrum_full  = new TH1D(name.c_str(),  "", 1000, 0., adcMax );
  TH1D* h1_spectrum_trunc = new TH1D("spectrum_trunc", "", 70, 300., adcMax );

  h1_spectrum_full->SetXTitle( "ADC Counts" );
  
  tree->Project(name.c_str(),  "adcData[0]");
  tree->Project("spectrum_trunc", "adcData[0]");

  // initialize position to point of maximum of truncated histo:
  float mean = h1_spectrum_trunc->GetBinCenter( h1_spectrum_trunc->GetMaximumBin() );
  float rms  = h1_spectrum_trunc->GetRMS()/2.; // just a guess (maybe 50 is better?)
  float integral  = h1_spectrum_trunc->Integral();

  float nSigma_lo = 1.;
  float nSigma_hi = 2.;

  float fitMin;
  float fitMax;

  TF1* f1_fit;

  if( useFermiDirac_ ) {

    fitMin = 220.;
    fitMax = adcMax;
    if     ( iChannel==1 )  fitMin = 220.;
    else if( iChannel==2 )  fitMin = 220.;
    else if( iChannel==3 )  fitMin = 230.;
    else if( iChannel==4 )  fitMin = 260.;
    else if( iChannel==5 )  fitMin = 220.;
    else if( iChannel==6 )  fitMin = 230.;
    else if( iChannel==7 )  fitMin = 250.;
    else if( iChannel==8 )  fitMin = 240.;

    f1_fit = new TF1("totFunc", FermiDiracPlusGaus, fitMin, fitMax, 6 );

    // FD normalization:
    f1_fit->SetParameter(0,300.);
    //f1_fit->SetParLimits(0,150., 500.);
    
    // 1  p0           3.05564e+02   2.02221e+01   3.31646e-02  -2.73120e-04
    // 2  p1           2.87305e+02   1.05400e+01   1.45260e-02  -5.31457e-04
    // 3  p2           8.16511e+01   1.60388e+00   8.57439e-03  -3.75531e-04
    // 4  p3           2.91747e+02   2.47999e+00   2.88441e-02  -3.04096e-04
    // 5  p4           4.50352e+02   3.33186e-01   4.18731e-03  -2.80640e-03
    // 6  p5           3.75147e+01   2.99857e-01   2.80790e-03  -3.28132e-03

    // FD flex point:
    f1_fit->SetParameter(1,400.);
    //f1_fit->SetParLimits(1,250., 600.);
    
    // FD temp:
    f1_fit->SetParameter(2,80.);
    //f1_fit->SetParLimits(2,0., 150.);
    
    // gauss normalization:
    f1_fit->SetParameter(3,300.);
    //f1_fit->SetParLimits(3,0.2*integral, 2.*integral);
    
    // gauss mean:
    //f1_fit->SetParameter(4,450.);
         if( iChannel==1 ) f1_fit->SetParameter(4, 145.8+192.4);
    else if( iChannel==2 ) f1_fit->SetParameter(4, 145.8+301.3);
    else if( iChannel==3 ) f1_fit->SetParameter(4, 145.8+258.5);
    else if( iChannel==4 ) f1_fit->SetParameter(4, 145.8+290.3);
    else if( iChannel==5 ) f1_fit->SetParameter(4, 145.8+197.8);
    else if( iChannel==6 ) f1_fit->SetParameter(4, 145.8+190.8);
    else if( iChannel==7 ) f1_fit->SetParameter(4, 145.8+225.7);
    else if( iChannel==8 ) f1_fit->SetParameter(4, 145.8+235.2);
    //f1_fit->SetParLimits(4,300., 500.);
    
    // gauss sigma:
    f1_fit->SetParameter(5,40.);
    //f1_fit->SetParLimits(5,10., 100.);
    
    f1_fit->SetLineColor(kRed);

    h1_spectrum_full->Fit( f1_fit, "R" );
    TF1* f1_fit = new TF1("totFunc", FermiDiracPlusGaus, fitMin, adcMax, 6 );

    // FD normalization:
    f1_fit->SetParameter(0,300.);
    //f1_fit->SetParLimits(0,150., 500.);
    
    // 1  p0           3.05564e+02   2.02221e+01   3.31646e-02  -2.73120e-04
    // 2  p1           2.87305e+02   1.05400e+01   1.45260e-02  -5.31457e-04
    // 3  p2           8.16511e+01   1.60388e+00   8.57439e-03  -3.75531e-04
    // 4  p3           2.91747e+02   2.47999e+00   2.88441e-02  -3.04096e-04
    // 5  p4           4.50352e+02   3.33186e-01   4.18731e-03  -2.80640e-03
    // 6  p5           3.75147e+01   2.99857e-01   2.80790e-03  -3.28132e-03

    // FD flex point:
    f1_fit->SetParameter(1,400.);
    //f1_fit->SetParLimits(1,250., 600.);
    
    // FD temp:
    f1_fit->SetParameter(2,80.);
    //f1_fit->SetParLimits(2,0., 150.);
    
    // gauss normalization:
    f1_fit->SetParameter(3,300.);
    //f1_fit->SetParLimits(3,0.2*integral, 2.*integral);
    
    // gauss mean:
    //f1_fit->SetParameter(4,450.);
         if( iChannel==1 ) f1_fit->SetParameter(4, 145.8+192.4);
    else if( iChannel==2 ) f1_fit->SetParameter(4, 145.8+301.3);
    else if( iChannel==3 ) f1_fit->SetParameter(4, 145.8+258.5);
    else if( iChannel==4 ) f1_fit->SetParameter(4, 145.8+290.3);
    else if( iChannel==5 ) f1_fit->SetParameter(4, 145.8+197.8);
    else if( iChannel==6 ) f1_fit->SetParameter(4, 145.8+190.8);
    else if( iChannel==7 ) f1_fit->SetParameter(4, 145.8+225.7);
    else if( iChannel==8 ) f1_fit->SetParameter(4, 145.8+235.2);
    //f1_fit->SetParLimits(4,300., 500.);
    
    // gauss sigma:
    f1_fit->SetParameter(5,40.);

  } else { // iterative gaussian fit

    fitMin = mean-nSigma_lo*rms;
    fitMax = mean+nSigma_hi*rms;

    TF1* f1_gaus = new TF1("gaus0", "gaus", fitMin, fitMax);
    f1_gaus->SetParameter(0, integral);
    f1_gaus->SetParameter(1, mean);
    f1_gaus->SetParameter(2, rms);

    int niter = 4;
    for( unsigned iter=0; iter<niter; ++iter ) {
      std::string option = (iter==niter-1) ? "R" : "RNQ";
      h1_spectrum_full->Fit( f1_gaus, option.c_str() );
      float mean_new = f1_gaus->GetParameter(1);
      float rms_new = f1_gaus->GetParameter(2);
      f1_gaus->SetRange( mean_new-nSigma_lo*rms_new, mean_new+nSigma_hi*rms_new );
    }

  }

  
  f1_fit->SetLineColor(kRed);

  h1_spectrum_full->Fit( f1_fit, "R" );




  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1_spectrum_full->Draw();
    
  c1->SaveAs( Form("%s/fit_ch%d.eps", outdir.c_str(), iChannel) );
  c1->SaveAs( Form("%s/fit_ch%d.png", outdir.c_str(), iChannel) );

  delete c1;


  return f1_fit->GetParameter(4);
  //return f1_gaus->GetParameter(1);

}
  
