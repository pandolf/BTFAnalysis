#include <iostream>
#include <string>
#include <stdlib.h>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"

#include "TMath.h"



struct FitResults {

   float ped_mu;
   float ped_mu_err;
   float ped_sigma;
   float ped_sigma_err;

   float mu;
   float mu_err;
   float offset;
   float offset_err;
   float Q1;
   float Q1_err;
   float sigma;
   float sigma_err;

};


int NORDERS;


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



FitResults fitSingleHisto( TFile* file, const std::string& histoName, double pedMin, double pedMax, double xMin, double xMax );













int main( int argc, char* argv[] ) {


  NORDERS = 5;
  if( argc>1 ) {
    NORDERS = atoi(argv[1]);
    std::cout << "-> NORDER is set to: " << NORDERS << std::endl;
  }



  TFile* file = TFile::Open("run_BTF_000001_cosmics_dqmPlots.root");

  FitResults fr_0 = fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_0", 120., 130., 138., 185. );
  FitResults fr_1 = fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_1", 100., 125., 125., 190. );
  FitResults fr_2 = fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_2", 100., 125., 128., 210. );
  FitResults fr_3 = fitSingleHisto( file, "CEF3RAW_cef3RawSpectrum_3", 100., 125., 137., 200. );

  TH1D* h1_Q1 = new TH1D("Q1", "", 4, -0.5, 3.5);
  h1_Q1->SetXTitle( "Channel Number");
  h1_Q1->SetBinContent( 1, fr_0.Q1 );
  h1_Q1->SetBinContent( 2, fr_1.Q1 );
  h1_Q1->SetBinContent( 3, fr_2.Q1 );
  h1_Q1->SetBinContent( 4, fr_3.Q1 );

  h1_Q1->SetBinError( 1, fr_0.Q1_err );
  h1_Q1->SetBinError( 2, fr_1.Q1_err );
  h1_Q1->SetBinError( 3, fr_2.Q1_err );
  h1_Q1->SetBinError( 4, fr_3.Q1_err );

  h1_Q1->SetMarkerStyle(20);
  h1_Q1->SetMarkerSize(1.6);
  //h1_Q1->SetMarkerColor();

  //TH1D* h1_offsetDiff = new TH1D("offsetDiff", "", 4, -0.5, 3.5);
  //h1_offsetDiff->SetXTitle( "Channel Number");
  //h1_offsetDiff->SetBinContent( 1, fr_0.offset - fr_0.ped_mu );
  //h1_offsetDiff->SetBinContent( 2, fr_1.offset - fr_1.ped_mu );
  //h1_offsetDiff->SetBinContent( 3, fr_2.offset - fr_2.ped_mu );
  //h1_offsetDiff->SetBinContent( 4, fr_3.offset - fr_3.ped_mu );

  //h1_offsetDiff->SetBinError( 1, sqrt( fr_0.offset_err*fr_0.offset_err + fr_0.ped_mu_err*fr_0.ped_mu_err) );
  //h1_offsetDiff->SetBinError( 2, sqrt( fr_1.offset_err*fr_1.offset_err + fr_1.ped_mu_err*fr_1.ped_mu_err) );
  //h1_offsetDiff->SetBinError( 3, sqrt( fr_2.offset_err*fr_2.offset_err + fr_2.ped_mu_err*fr_2.ped_mu_err) );
  //h1_offsetDiff->SetBinError( 4, sqrt( fr_3.offset_err*fr_3.offset_err + fr_3.ped_mu_err*fr_3.ped_mu_err) );

  //h1_offsetDiff->SetMarkerStyle(21);
  //h1_offsetDiff->SetMarkerSize(1.6);
  //h1_offsetDiff->SetMarkerColor(kRed);

  TH1D* h1_offset = new TH1D("offset", "", 4, -0.5, 3.5);
  h1_offset->SetXTitle( "Channel Number");
  h1_offset->SetBinContent( 1, fr_0.offset );
  h1_offset->SetBinContent( 2, fr_1.offset );
  h1_offset->SetBinContent( 3, fr_2.offset );
  h1_offset->SetBinContent( 4, fr_3.offset );

  h1_offset->SetBinError( 1, fr_0.offset_err );
  h1_offset->SetBinError( 2, fr_1.offset_err );
  h1_offset->SetBinError( 3, fr_2.offset_err );
  h1_offset->SetBinError( 4, fr_3.offset_err );

  h1_offset->SetMarkerStyle(21);
  h1_offset->SetMarkerSize(1.6);
  h1_offset->SetMarkerColor(kRed);

  TH1D* h1_ped_mu = new TH1D("ped_mu", "", 4, -0.5, 3.5);
  h1_ped_mu->SetXTitle( "Channel Number");
  h1_ped_mu->SetBinContent( 1, fr_0.ped_mu );
  h1_ped_mu->SetBinContent( 2, fr_1.ped_mu );
  h1_ped_mu->SetBinContent( 3, fr_2.ped_mu );
  h1_ped_mu->SetBinContent( 4, fr_3.ped_mu );

  h1_ped_mu->SetBinError( 1, fr_0.ped_mu_err );
  h1_ped_mu->SetBinError( 2, fr_1.ped_mu_err );
  h1_ped_mu->SetBinError( 3, fr_2.ped_mu_err );
  h1_ped_mu->SetBinError( 4, fr_3.ped_mu_err );

  h1_ped_mu->SetMarkerStyle(20);
  h1_ped_mu->SetMarkerSize(1.6);
  h1_ped_mu->SetMarkerColor(38);



  // set the TStyle
  TStyle* style_ = new TStyle("DrawBaseStyle", "");
  style_->SetCanvasColor(0);
  style_->SetPadColor(0);
  style_->SetFrameFillColor(0);
  style_->SetStatColor(0);
  style_->SetOptStat(0);
  style_->SetTitleFillColor(0);
  style_->SetCanvasBorderMode(0);
  style_->SetPadBorderMode(0);
  style_->SetFrameBorderMode(0);
  style_->SetPadBottomMargin(0.12);
  style_->SetPadLeftMargin(0.12);
  style_->cd();

  // For the canvas:
  style_->SetCanvasBorderMode(0);
  style_->SetCanvasColor(kWhite);
  style_->SetCanvasDefH(600); //Height of canvas
  style_->SetCanvasDefW(600); //Width of canvas
  style_->SetCanvasDefX(0); //POsition on screen
  style_->SetCanvasDefY(0);

  // For the Pad:
  style_->SetPadBorderMode(0);
  style_->SetPadColor(kWhite);
  style_->SetPadGridX(false);
  style_->SetPadGridY(false);
  style_->SetGridColor(0);
  style_->SetGridStyle(3);
  style_->SetGridWidth(1);

  // For the frame:
  style_->SetFrameBorderMode(0);
  style_->SetFrameBorderSize(1);
  style_->SetFrameFillColor(0);
  style_->SetFrameFillStyle(0);
  style_->SetFrameLineColor(1);
  style_->SetFrameLineStyle(1);
  style_->SetFrameLineWidth(1);


  // Margins:
  style_->SetPadTopMargin(0.05);
  style_->SetPadBottomMargin(0.15);//0.13);
  style_->SetPadLeftMargin(0.15);//0.16);
  style_->SetPadRightMargin(0.05);//0.02);

  // For the Global title:

  style_->SetOptTitle(0);
  style_->SetTitleFont(42);
  style_->SetTitleColor(1);
  style_->SetTitleTextColor(1);
  style_->SetTitleFillColor(10);
  style_->SetTitleFontSize(0.05);

  // For the axis titles:

  style_->SetTitleColor(1, "XYZ");
  style_->SetTitleFont(42, "XYZ");
  style_->SetTitleSize(0.05, "XYZ");
  style_->SetTitleXOffset(1.15);//0.9);
  style_->SetTitleYOffset(1.4); // => 1.15 if exponents

  // For the axis labels:

  style_->SetLabelColor(1, "XYZ");
  style_->SetLabelFont(42, "XYZ");
  style_->SetLabelOffset(0.007, "XYZ");
  style_->SetLabelSize(0.045, "XYZ");

  // For the axis:

  style_->SetAxisColor(1, "XYZ");
  style_->SetStripDecimals(kTRUE);
  style_->SetTickLength(0.03, "XYZ");
  style_->SetNdivisions(510, "XYZ");
  style_->SetPadTickX(1); // To get tick marks on the opposite side of the frame
  style_->SetPadTickY(1);

  style_->cd();


  TCanvas* c2 = new TCanvas( "c2", "", 600, 600 );

  c2->cd();

  TH2D* h2_axes = new TH2D("axes", "", 10, -0.5, 3.5, 10, 0., 200.);
  h2_axes->SetXTitle( "Channel Number" );
  h2_axes->SetYTitle( "ADC Counts" );
  h2_axes->Draw();

  TLegend* legend = new TLegend( 0.4, 0.7, 0.68, 0.9 );
  legend->SetFillColor(0);
  legend->SetFillStyle(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( h1_Q1, "Q1", "P" );
  legend->AddEntry( h1_offset, "Offset", "P" );
  legend->AddEntry( h1_ped_mu, "Pedestal", "P" );
  legend->Draw("same");

  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(62);
  label_top->AddText("Cosmic Run");
  label_top->Draw("same");

  gStyle->SetErrorX(0);

  h1_Q1->Draw("P same");
  h1_offset->Draw("P same");
  h1_ped_mu->Draw("P same");

  c2->SaveAs("summaryPlot.eps");
  c2->SaveAs("summaryPlot.png");


  return 0;

}



FitResults fitSingleHisto( TFile* file, const std::string& histoName, double pedMin, double pedMax, double xMin, double xMax ) {

  TH1D* histo = (TH1D*)file->Get(histoName.c_str());
  float integral = histo->Integral();

  TF1* f1_ped = new TF1( "ped", "gaus", pedMin, pedMax );
  f1_ped->SetParameter(0, integral);
  f1_ped->SetParameter(1, 110.);
  f1_ped->SetParameter(2, 10.);

  f1_ped->SetLineColor(kRed);

  histo->Fit( f1_ped, "RQN" );

  int nSteps = 4;
  for( unsigned iStep=0; iStep<nSteps; iStep++ ) {

    float ped_mean = f1_ped->GetParameter(1);
    float ped_sigma = f1_ped->GetParameter(2);

    float nSigma = 2.;
    float newMin = ped_mean-nSigma*ped_sigma;
    float newMax = ped_mean+nSigma*ped_sigma;

    f1_ped->SetRange( newMin, newMax );
 
    std::string option = (iStep<(nSteps-1)) ? "RQN" : "RQ+";
    histo->Fit( f1_ped, option.c_str() );

  }




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


  histo->Fit( f1, "R+" );

  TCanvas* c1 = new TCanvas("c1", "c1", 600, 600);
  c1->cd();

  c1->SetLogy();

  TH2D* h2_axes = new TH2D("axes", "", 10, 100., 350., 10, 9., 7.*histo->GetMaximum() );
  h2_axes->SetXTitle( "ADC Counts" );
  h2_axes->Draw();

  histo->Draw("same");

  std::string canvasName = histoName + ".eps";
  c1->SaveAs( canvasName.c_str() );


  FitResults fr;
  fr.ped_mu = f1_ped->GetParameter(1);
  fr.ped_mu_err = f1_ped->GetParError(1);
  fr.ped_sigma = f1_ped->GetParameter(2);
  fr.ped_sigma_err = f1_ped->GetParError(2);

  fr.mu = f1->GetParameter(1);
  fr.mu_err = f1->GetParError(1);
  fr.offset = f1->GetParameter(4);
  fr.offset_err = f1->GetParError(4);
  fr.Q1 = f1->GetParameter(2);
  fr.Q1_err = f1->GetParError(2);
  fr.sigma = f1->GetParameter(3);
  fr.sigma_err = f1->GetParError(3);

  delete c1;
  delete f1;
  delete f1_ped;
  delete h2_axes;

  return fr;

}
