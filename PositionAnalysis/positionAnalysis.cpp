#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"


std::vector<float> getPedestals( const std::string& fileName );



int main( int argc, char* argv[] ) {


  std::string runName = "precalib_BGO_pedestal_noSource";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }

  std::string fileName = "data/run_" + runName + ".root";
  TFile* file = TFile::Open(fileName.c_str());
  if( file==0 ) {
    std::cout << "ERROR! Din't find file " << fileName << std::endl;
    std::cout << "Exiting." << std::endl;
    exit(11);
  }
  TTree* tree = (TTree*)file->Get("eventRawData");


  std::vector<float> pedestals = getPedestals("run_BTF_000002_ped_dqmPlots.root");
  std::cout << std::endl;
  std::cout << "-> Got pedestals: " << std::endl;
  std::cout << " Channel 0: " << pedestals[0] << std::endl;
  std::cout << " Channel 1: " << pedestals[1] << std::endl;
  std::cout << " Channel 2: " << pedestals[2] << std::endl;
  std::cout << " Channel 3: " << pedestals[3] << std::endl;
  std::cout << std::endl;

  UInt_t evtNumber;
  tree->SetBranchAddress( "evtNumber", &evtNumber );
  UInt_t adcData[40];
  tree->SetBranchAddress( "adcData", adcData );
  UInt_t adcBoard[40];
  tree->SetBranchAddress( "adcBoard", adcBoard );
  UInt_t adcChannel[40];
  tree->SetBranchAddress( "adcChannel", adcChannel );


  int nBins = 500;
  float xMax = 40.;

  TH1D* h1_xPos = new TH1D("xPos", "", nBins, -xMax, xMax);
  TH1D* h1_yPos = new TH1D("yPos", "", nBins, -xMax, xMax);

  TH2D* h2_xyPos = new TH2D("xyPos", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH1D* h1_cef3_0   = new TH1D("cef3_0",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_1   = new TH1D("cef3_1",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_2   = new TH1D("cef3_2",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_3   = new TH1D("cef3_3",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_tot = new TH1D("cef3_tot", "", 200, 0., 4.*4097.);
  
  TH1D* h1_cef3_corr_0   = new TH1D("cef3_corr_0",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_corr_1   = new TH1D("cef3_corr_1",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_corr_2   = new TH1D("cef3_corr_2",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_corr_3   = new TH1D("cef3_corr_3",   "", 4097, 0., 4097.);
  TH1D* h1_cef3_corr_tot = new TH1D("cef3_corr_tot", "", 200, 0., 4.*4097.);
  


  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( iEntry % 100 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


    // CeF3 fibres are channels 2-5 of board 1
    float cef3_0(-1.);
    float cef3_1(-1.);
    float cef3_2(-1.);
    float cef3_3(-1.);

    for( unsigned i=0; i<40; ++i ) {
      int board  = adcBoard[i];
      int channel= adcChannel[i];
      if( board!=1 ) continue;
      if( channel==2 ) cef3_0 = adcData[i];
      if( channel==3 ) cef3_1 = adcData[i];
      if( channel==4 ) cef3_2 = adcData[i];
      if( channel==5 ) cef3_3 = adcData[i];
    }

    if( cef3_0<0. ) {
      std::cout << "(Event: " << evtNumber << ") Didn't find data for fibre N.0!" << std::endl;
      continue;
    }

    if( cef3_1<0. ) {
      std::cout << "(Event: " << evtNumber << ") Didn't find data for fibre N.1!" << std::endl;
      continue;
    }

    if( cef3_2<0. ) {
      std::cout << "(Event: " << evtNumber << ") Didn't find data for fibre N.2!" << std::endl;
      continue;
    }

    if( cef3_3<0. ) {
      std::cout << "(Event: " << evtNumber << ") Didn't find data for fibre N.3!" << std::endl;
      continue;
    }

    if( cef3_0>=4095. ) {
      std::cout << "(Event: " << evtNumber << ") Fibre N.0 in overflow!" << std::endl;
      continue;
    }

    if( cef3_1>=4095. ) {
      std::cout << "(Event: " << evtNumber << ") Fibre N.1 in overflow!" << std::endl;
      continue;
    }

    if( cef3_2>=4095. ) {
      std::cout << "(Event: " << evtNumber << ") Fibre N.2 in overflow!" << std::endl;
      continue;
    }

    if( cef3_3>=4095. ) {
      std::cout << "(Event: " << evtNumber << ") Fibre N.3 in overflow!" << std::endl;
      continue;
    }


    // subtract pedestals:
    float cef3_0_corr = cef3_0 - pedestals[0];
    float cef3_1_corr = cef3_1 - pedestals[1];
    float cef3_2_corr = cef3_2 - pedestals[2];
    float cef3_3_corr = cef3_3 - pedestals[3];

    //   0      1
    //          
    //          
    //   2      3


    float xySize = 25.; // in mm
    float chamfer = 2.1; // in mm

    float position = xySize/2. - chamfer/4.;

    float xPosW_0 = cef3_0_corr*(-position);
    float xPosW_1 = cef3_1_corr*(+position);
    float xPosW_2 = cef3_2_corr*(-position);
    float xPosW_3 = cef3_3_corr*(+position);

    float yPosW_0 = cef3_0_corr*(+position);
    float yPosW_1 = cef3_1_corr*(+position);
    float yPosW_2 = cef3_2_corr*(-position);
    float yPosW_3 = cef3_3_corr*(-position);

    float eTot = cef3_0+cef3_1+cef3_2+cef3_3;
    float eTot_corr = cef3_0_corr+cef3_1_corr+cef3_2_corr+cef3_3_corr;

    float xPos = (xPosW_0+xPosW_1+xPosW_2+xPosW_3)/eTot_corr;
    float yPos = (yPosW_0+yPosW_1+yPosW_2+yPosW_3)/eTot_corr;

    h1_cef3_0->Fill( cef3_0 );
    h1_cef3_1->Fill( cef3_1 );
    h1_cef3_2->Fill( cef3_2 );
    h1_cef3_3->Fill( cef3_3 );
std::cout << "eTot: " << eTot << std::endl;
std::cout << "eTot_corr: " << eTot_corr << std::endl;
    h1_cef3_3->Fill( eTot );

    h1_cef3_corr_0->Fill( cef3_0_corr );
    h1_cef3_corr_1->Fill( cef3_1_corr );
    h1_cef3_corr_2->Fill( cef3_2_corr );
    h1_cef3_corr_3->Fill( cef3_3_corr );
    h1_cef3_corr_tot->Fill( eTot_corr );

    h1_xPos->Fill( xPos );
    h1_yPos->Fill( yPos );

    h2_xyPos->Fill( xPos, yPos );

  }


  std::string outfileName = "PosAn_" + runName + ".root";
  TFile* outfile = TFile::Open( outfileName.c_str(), "RECREATE" );
  outfile->cd();

  h1_xPos->Write();
  h1_yPos->Write();
  
  h2_xyPos->Write();

  h1_cef3_0->Write();
  h1_cef3_1->Write();
  h1_cef3_2->Write();
  h1_cef3_3->Write();
  
  h1_cef3_corr_0->Write();
  h1_cef3_corr_1->Write();
  h1_cef3_corr_2->Write();
  h1_cef3_corr_3->Write();
  
  h1_cef3_tot->Write();
  h1_cef3_corr_tot->Write();
  
  
  outfile->Close();

  return 0;

}



std::vector<float> getPedestals( const std::string& fileName ) {

  TFile* file = TFile::Open(fileName.c_str());

  TH1D* h1_ped0 = (TH1D*)file->Get("CEF3RAW_cef3RawSpectrum_0");
  TH1D* h1_ped1 = (TH1D*)file->Get("CEF3RAW_cef3RawSpectrum_1");
  TH1D* h1_ped2 = (TH1D*)file->Get("CEF3RAW_cef3RawSpectrum_2");
  TH1D* h1_ped3 = (TH1D*)file->Get("CEF3RAW_cef3RawSpectrum_3");

  std::vector<float> peds;
  peds.push_back(h1_ped0->GetMean());
  peds.push_back(h1_ped1->GetMean());
  peds.push_back(h1_ped2->GetMean());
  peds.push_back(h1_ped3->GetMean());

  return peds;
  
}
