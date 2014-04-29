#include <iostream>
#include <string>
#include <cmath>
#include <stdlib.h>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "TMath.h"


int NORDERS;
//double XMAX;


Double_t PMTFunction(Double_t *x, Double_t *par)
{

   float N = par[0];
   float mu = par[1];
   float Q1 = par[2];
   float sigma = par[3];
   float offset = par[4];
   float sigmaoffset = par[5];
   float alpha = par[6];
   float w = par[7];

   float xx = x[0];
   double value = 0.;

   for( unsigned i=1; i<NORDERS; ++i ) {

     //double Qn = offset + (double)(i)*Q1;
     double sigma_n = sqrt( (double)(i)*sigma*sigma + sigmaoffset*sigmaoffset);

     //double poisson = TMath::Poisson( i, mu );
     //double gauss   = TMath::Gaus( xx, Qn, sigma_n );

     //double xxp = xx     - Qn - alpha*sigma_n*sigma_n;
     //double Q0p = offset - Qn - alpha*sigma_n*sigma_n;
     //double bg      = 0.5*alpha * TMath::Exp(-alpha*xxp)* (
     //                     TMath::Erf( abs(Q0p)/(sigma_n*sqrt(2) ) ) +  xxp/abs(xxp) * TMath::Erf( abs(xxp)/(sigma_n*sqrt(2)) ) );
     //value = value + N*( poisson * ( (1.-w)*gauss + w*bg ) );
     value = value + N*(TMath::Poisson( i, mu ) * TMath::Gaus( xx, (double)i*Q1 + offset, sigma_n) );
     //value = value + N*(TMath::Poisson( i, mu ) * TMath::Gaus( xx, (double)i*Q1 + offset, sqrt((double)i)*sigma ));
   }

   return value;

}



//void fitSingleHisto( TFile* file, const std::string& histoName, double xMin );
void fitSingleHisto( TFile* file, const std::string& histoName, double xMin, double xMax );

int main( int argc, char* argv[] ) {


  NORDERS = 5;
  if( argc>1 ) {
    NORDERS = atoi(argv[1]);
    std::cout << "-> NORDER is set to: " << NORDERS << std::endl;
  }

  //XMAX = 250.;
  //if( argc>2 ) {
  //  XMAX = (double)atoi(argv[2]);
  //  std::cout << "-> XMAX is set to: " << XMAX << std::endl;
  //}

  TFile* file = TFile::Open("run_BTF_000001_cosmics_dqmPlots.root");

  fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_0", 138., 185. );
  fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_1", 125., 190. );
  fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_2", 128., 210. );
  fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_3", 137., 200. );

  return 0;

}



void fitSingleHisto( TFile* file, const std::string& histoName, double xMin, double xMax ) {

  TH1D* histo = (TH1D*)file->Get(histoName.c_str());
  float integral = histo->Integral();

  TF1* f1 = new TF1( "func", PMTFunction, xMin, xMax, 8 );
  f1->SetParameter( 0, integral ); //normalization
  f1->SetParameter( 1, 1. ); //poiss mu
  f1->SetParameter( 2, 25. ); //gauss step
  f1->SetParameter( 3, 10. ); //gauss sigma
  f1->SetParameter( 4, 100 ); //offset
  f1->SetParameter( 5, 3. ); //sigmaoffset
  f1->SetParameter( 6, 0.03 ); //alpha
  f1->SetParameter( 7, 0.4 ); //w

  f1->FixParameter( 5, 0. ); //sigmaoffset
  f1->FixParameter( 6, 0. ); //alpha
  f1->FixParameter( 7, 0. ); //w

  f1->SetParLimits( 1, 0.5, 2.5 ); //poiss mu
  f1->SetParLimits( 2, 10., 40. ); //gauss step
  f1->SetParLimits( 3, 3., 12. ); //gauss sigma
  //f1->SetParLimits( 4, 90., 110.); //offset
  //f1->SetParLimits( 5, 0., 8. ); //gauss sigma


  histo->Fit( f1, "R" );

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();

  c1->SetLogy();

  TH2D* h2_axes = new TH2D("axes", "", 10, 100., 350., 10, 9., 7.*histo->GetMaximum() );
  h2_axes->SetXTitle( "ADC Counts" );
  h2_axes->Draw();

  histo->Draw("same");

  std::string canvasName = histoName + ".eps";
  c1->SaveAs( canvasName.c_str() );

  delete c1;
  delete f1;
  delete h2_axes;

}
