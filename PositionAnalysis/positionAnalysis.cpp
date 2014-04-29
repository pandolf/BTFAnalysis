#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"





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



  UInt_t evtNumber;
  tree->SetBranchAddress( "evtNumber", &evtNumber );
  UInt_t adcData[40];
  tree->SetBranchAddress( "adcData", adcData );
  UInt_t adcBoard[40];
  tree->SetBranchAddress( "adcBoard", adcBoard );
  UInt_t adcChannel[40];
  tree->SetBranchAddress( "adcChannel", adcChannel );


  TH1D* h1_xPos = new TH1D("xPos", "", 100, -100., 100.);
  TH1D* h1_yPos = new TH1D("yPos", "", 100, -100., 100.);

  TH2D* h2_xyPos = new TH2D("xyPos", "", 100, -100., 100., 100, -100., 100.);

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

    int xPos = (cef3_0+cef3_2) - (cef3_1+cef3_3);
    int yPos = (cef3_0+cef3_1) - (cef3_2+cef3_3);

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
  
  outfile->Close();

  return 0;

}
