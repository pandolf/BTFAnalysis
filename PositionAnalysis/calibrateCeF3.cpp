#include <iostream>
#include <fstream>

#include "TFile.h"
#include "TTree.h"
#include "TH1D.h"
#include "TF1.h"
#include "TCanvas.h"



void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name );
TF1* fitSingleElectronPeak( const std::string& outputdir, int i, TTree* tree );
TF1* checkTotalResolution( const std::string& outputdir, TTree* tree );


int main() {

  std::string runName = "BTF_92_20140430-020137_beam_uncalib";
  TFile* file = TFile::Open(Form("PosAn_%s.root", runName.c_str()) );
  TTree* tree = (TTree*)file->Get("tree_passedEvents");


  
  std::string outputdir = "CeF3_Calibration_" + runName;
  std::string mkdir_command = "mkdir -p " + outputdir;
  system( mkdir_command.c_str() );

  std::string ofsName = outputdir + "/constants.txt";
  ofstream ofs(ofsName.c_str());

  std::vector<float> calibConstants;
  for( unsigned i=0; i<4; ++i ) {
    TF1* f1 = fitSingleElectronPeak( outputdir, i, tree );
    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    calibConstants.push_back(sigma/mean);
    std::cout << std::endl;
    std::cout << "Channel " << i << std::endl;
    std::cout << "  Mean       : " << mean << std::endl;
    std::cout << "  Sigma      : " << sigma << std::endl;
    std::cout << "  Resolution : " << sigma/mean << std::endl;

    ofs << i << "\t" << mean << std::endl;

  }

  ofs.close();

  checkTotalResolution( outputdir, tree );

  std::cout << "-> Saved constants in: " << ofsName << std::endl;

  return 0;

}



TF1* fitSingleElectronPeak( const std::string& outputdir, int i, TTree* tree ) {

  std::string histoName(Form("h1_%d", i));
  TH1D* h1 = new TH1D(histoName.c_str(), "", 100, 0., 3000.);
  tree->Project( histoName.c_str(), Form("cef3_corr[%d]", i), "(nHodoFibersCorrX==1 && nHodoFibersCorrY==1)");

  TF1* f1 = new TF1( Form("gaus_%d", i), "gaus", 400., 1200.);
  f1->SetParameter(0, 3000.);
  f1->SetParameter(1, 800.);
  f1->SetParameter(2, 150.);

  doSingleFit( h1, f1, outputdir, Form("%d", i) );

  return f1;

}



TF1* checkTotalResolution( const std::string& outputdir, TTree* tree ) {

  std::string histoName("h1_tot");
  TH1D* h1 = new TH1D(histoName.c_str(), "", 500, 0., 12000.);
  tree->Project( histoName.c_str(), "cef3_corr[0]+cef3_corr[1]+cef3_corr[2]+cef3_corr[3]", "(nHodoFibersCorrX==1 && nHodoFibersCorrY==1)");

  TF1* f1 = new TF1("gaus_tot", "gaus", 1600., 4800.);
  f1->SetParameter(0, 3000.);
  f1->SetParameter(1, 3000.);
  f1->SetParameter(2, 600.);

  doSingleFit( h1, f1, outputdir, "tot" );

  std::cout << std::endl;
  std::cout << std::endl;
  std::cout << "Total resolution: " << f1->GetParameter(2)/f1->GetParameter(1) << std::endl;

  return f1;

}


void doSingleFit( TH1D* h1, TF1* f1, const std::string& outputdir, const std::string& name ) {

  h1->Fit( f1, "RQN" );

  int niter = 4.;
  float nSigma = 2.;

  for( unsigned iter=0; iter<niter; iter++ ) {

    float mean  = f1->GetParameter(1);
    float sigma = f1->GetParameter(2);
    float fitMin = mean - nSigma*sigma;
    float fitMax = mean + nSigma*sigma;
    f1->SetRange( fitMin, fitMax );
    if( iter==(niter-1) )
      h1->Fit( f1, "RQN" );
    else
      h1->Fit( f1, "RQ+" );
  }


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  h1->Draw();

  c1->SaveAs( Form("%s/fit_%s.eps", outputdir.c_str(), name.c_str()) );
  c1->SaveAs( Form("%s/fit_%s.png", outputdir.c_str(), name.c_str()) );

  delete c1;

}
