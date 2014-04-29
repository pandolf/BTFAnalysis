#include <iostream>
#include <string>

#include "TFile.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"

#include "TMath.h"




Double_t PMTFunction(Double_t *x, Double_t *par)
{

   float mu = par[0];
   float Q1 = par[1];
   float sigma = par[2];

   float xx = x[0];
   double value = 0.;

   for( unsigned i=0; i<4; ++i ) // upt o 4 photoelectrons for now
     value += (TMath::Poisson( i, mu ) * TMath::Gaus( xx, (double)i*Q1, (double)i*sigma ));

   return value;

}



void fitSingleHisto( TFile* file, const std::string& histoName );

int main() {

  TFile* file = TFile::Open("run_BTF_000001_cosmics_dqmPlots.root");

  fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_2" );

  return 0;

}



void fitSingleHisto( TFile* file, const std::string& histoName ) {

  TH1D* histo = (TH1D*)file->Get(histoName.c_str());

  TF1* f1 = new TF1( "func", PMTFunction, 0., 350., 3 );

  histo->Fit( f1, "R" );

  TFile* outfile = TFile::Open("prova.root", "recreate");
  outfile->cd();
  histo->Write();

  outfile->Close();

}
