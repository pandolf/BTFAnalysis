#include <cstdlib>
#include <fstream>
#include <iostream>
#include <cmath>

#include "TH1D.h"
#include "TH2D.h"
#include "TFile.h"
#include "TTree.h"

#include "fastDQM_CeF3_BTF.h"





std::vector< std::pair<float, float> > getPedestals( const std::string& fileName, const std::string& name );
std::vector< std::pair<float, float> > getPedestalsHodo( const std::string type, const std::string fileName );
std::vector<float> subtractPedestals( std::vector<float> raw, std::vector< std::pair<float, float> > pedestals, float nSigma );
float sumVector( std::vector<float> v );
bool checkVector( std::vector<float> v, float theMax=4095. );
float getMeanposHodo( std::vector<float> hodo, std::vector< std::pair<float, float> > pedestals);



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


  std::vector<std::pair<float, float> > pedestals = getPedestals("run_BTF_000002_ped_dqmPlots.root", "CEF3RAW_cef3RawSpectrum" );
  std::cout << std::endl;
  std::cout << "-> Got pedestals of CeF3: " << std::endl;
  for( unsigned i=0; i<CEF3_CHANNELS; ++i )
    std::cout << " CeF3 Channel " << i << ": " << pedestals[i].first << std::endl;
  std::cout << std::endl;

  std::vector<std::pair<float, float> > pedestals_bgo = getPedestals("run_BTF_000002_ped_dqmPlots.root", "BGORAW_bgoRawSpectrum" );
  std::cout << std::endl;
  std::cout << "-> Got pedestals of BGO: " << std::endl;
  for( unsigned i=0; i<BGO_CHANNELS; ++i )
    std::cout << " BGO Channel " << i << ": " << pedestals_bgo[i].first << std::endl;
  std::cout << std::endl;

  std::vector<std::pair<float, float> > pedestals_hodox = getPedestalsHodo("X", "run_BTF_000002_ped_dqmPlots.root");
  std::vector<std::pair<float, float> > pedestals_hodoy = getPedestalsHodo("Y", "run_BTF_000002_ped_dqmPlots.root");
  std::cout << "-> Got Hodoscope pedestals: " << std::endl;
  std::cout << std::endl;
  for( unsigned i=0; i<HODOX_CHANNELS; ++i )
    std::cout << "Channel " << i << ":  X: " << pedestals_hodox[i].first << " Y: " << pedestals_hodoy[i].first << std::endl;


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

  TH1D* h1_bgo_corr_0   = new TH1D("bgo_corr_0",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_1   = new TH1D("bgo_corr_1",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_2   = new TH1D("bgo_corr_2",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_3   = new TH1D("bgo_corr_3",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_4   = new TH1D("bgo_corr_4",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_5   = new TH1D("bgo_corr_5",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_6   = new TH1D("bgo_corr_6",   "", 4097, 0., 4097.);
  TH1D* h1_bgo_corr_7   = new TH1D("bgo_corr_7",   "", 4097, 0., 4097.);
  
  TH1D* h1_xPos_hodo = new TH1D("xPos_hodo", "", nBins, -xMax, xMax);
  TH1D* h1_yPos_hodo = new TH1D("yPos_hodo", "", nBins, -xMax, xMax);
  TH2D* h2_xyPos_hodo = new TH2D("xyPos_hodo", "", nBins, -xMax, xMax, nBins, -xMax, xMax);

  TH2D* h2_correlation_xPos = new TH2D("correlation_xPos", "", nBins, -xMax, xMax,  nBins, -xMax, xMax);
  TH2D* h2_correlation_yPos = new TH2D("correlation_yPos", "", nBins, -xMax, xMax,  nBins, -xMax, xMax);


  int nentries = tree->GetEntries();

  for( unsigned iEntry=0; iEntry<nentries; ++iEntry ) {

    tree->GetEntry(iEntry);

    if( iEntry % 100 == 0 ) std::cout << "Entry: " << iEntry << " / " << nentries << std::endl;


    // CeF3 fibres
    std::vector<float> cef3;
    for( unsigned i=0; i<CEF3_CHANNELS; ++i ) cef3.push_back(-1.);

    // BGO
    std::vector<float> bgo;
    for( unsigned i=0; i<BGO_CHANNELS; ++i ) bgo.push_back(-1.);

    // hodoX
    std::vector<float> hodox;
    for( unsigned i=0; i<HODOX_CHANNELS; ++i ) hodox.push_back(-1.);

    // hodoY
    std::vector<float> hodoy;
    for( unsigned i=0; i<HODOY_CHANNELS; ++i ) hodoy.push_back(-1.);


    for( unsigned i=0; i<40; ++i ) {
      int board  = adcBoard[i];
      int channel= adcChannel[i];

      if( board==0 ) {
        if     ( channel==(BGO_ADC_START_CHANNEL  ) )  bgo[0] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+1) )  bgo[1] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+2) )  bgo[2] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+3) )  bgo[3] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+4) )  bgo[4] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+5) )  bgo[5] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+6) )  bgo[6] = adcData[i];
        else if( channel==(BGO_ADC_START_CHANNEL+7) )  bgo[7] = adcData[i];
      } else if( board==1 ) {
        if     ( channel==(CEF3_ADC_START_CHANNEL  ) )  cef3[0]  = adcData[i];
        else if( channel==(CEF3_ADC_START_CHANNEL+1) )  cef3[1]  = adcData[i];
        else if( channel==(CEF3_ADC_START_CHANNEL+2) )  cef3[2]  = adcData[i];
        else if( channel==(CEF3_ADC_START_CHANNEL+3) )  cef3[3]  = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL  ) ) hodox[0] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+1) ) hodox[1] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+2) ) hodox[2] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+3) ) hodox[3] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+4) ) hodox[4] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+5) ) hodox[5] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+6) ) hodox[6] = adcData[i];
        else if( channel==(HODOX_ADC_START_CHANNEL+7) ) hodox[7] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL  ) ) hodoy[0] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+1) ) hodoy[1] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+2) ) hodoy[2] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+3) ) hodoy[3] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+4) ) hodoy[4] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+5) ) hodoy[5] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+6) ) hodoy[6] = adcData[i];
        else if( channel==(HODOY_ADC_START_CHANNEL+7) ) hodoy[7] = adcData[i];
      }
    }


    std::vector<float>  bgo_corr = subtractPedestals( bgo , pedestals_bgo, 2. );
    std::vector<float> cef3_corr = subtractPedestals( cef3, pedestals,     2. );

    bool cef3_ok = checkVector(cef3);
    bool cef3_corr_ok = checkVector(cef3_corr);



    // FIRST GET POSITION FROM HODOSCOPE:

    float xPos_hodo = getMeanposHodo(hodox, pedestals_hodox);
    float yPos_hodo = getMeanposHodo(hodoy, pedestals_hodoy);

    h1_xPos_hodo->Fill(xPos_hodo);
    h1_yPos_hodo->Fill(yPos_hodo);

    h2_xyPos_hodo->Fill(xPos_hodo, yPos_hodo);




    if( cef3_ok ) {

      // THEN USE CeF3 DATA:

      float eTot      = sumVector(cef3);
      float eTot_corr = sumVector(cef3_corr);

      h1_cef3_0->Fill( cef3[0] );
      h1_cef3_1->Fill( cef3[1] );
      h1_cef3_2->Fill( cef3[2] );
      h1_cef3_3->Fill( cef3[3] );
      h1_cef3_3->Fill( eTot );

      h1_cef3_corr_0->Fill( cef3_corr[0] );
      h1_cef3_corr_1->Fill( cef3_corr[1] );
      h1_cef3_corr_2->Fill( cef3_corr[2] );
      h1_cef3_corr_3->Fill( cef3_corr[3] );
      h1_cef3_corr_tot->Fill( eTot_corr );
      

      if( cef3_corr_ok ) {

        //   0      1
        //          
        //          
        //   3      2


        float xySize = 25.; // in mm
        float chamfer = 2.1; // in mm

        float position = xySize/2. - chamfer/4.;

        float xPosW_0 = cef3_corr[0]*(-position);
        float xPosW_1 = cef3_corr[1]*(+position);
        float xPosW_2 = cef3_corr[2]*(+position);
        float xPosW_3 = cef3_corr[3]*(-position);

        float yPosW_0 = cef3_corr[0]*(+position);
        float yPosW_1 = cef3_corr[1]*(+position);
        float yPosW_2 = cef3_corr[2]*(-position);
        float yPosW_3 = cef3_corr[3]*(-position);


        float xPos = (xPosW_0+xPosW_1+xPosW_2+xPosW_3)/eTot_corr;
        float yPos = (yPosW_0+yPosW_1+yPosW_2+yPosW_3)/eTot_corr;

        h1_xPos->Fill( xPos );
        h1_yPos->Fill( yPos );

        h2_xyPos->Fill( xPos, yPos );


        // CORRELATIONS BETWEEN CALO AND HODO:
  
        h2_correlation_xPos->Fill( xPos, xPos_hodo );
        h2_correlation_yPos->Fill( yPos, yPos_hodo );

      } // if cef3_ok

    }

  }


  std::cout << "-> Events passing overflow cut: " << h1_xPos->GetEntries() << "/" << nentries << " (" << 100.* h1_xPos->GetEntries()/nentries << "%)" << std::endl;

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
  
  h1_xPos_hodo->Write();
  h1_yPos_hodo->Write();
  h2_xyPos_hodo->Write();

  h2_correlation_xPos->Write();
  h2_correlation_yPos->Write();
  
  outfile->Close();
  std::cout << "-> Histograms saved in: " << outfile->GetName() << std::endl;

  return 0;

}



std::vector< std::pair<float, float> > getPedestals( const std::string& fileName, const std::string& name ) {

  TFile* file = TFile::Open(fileName.c_str());

  std::vector< std::pair<float, float> > peds;
  for( int i=0; i<CEF3_CHANNELS; ++i ) {
    TH1D* h1_ped = (TH1D*)file->Get(Form("%s_%d", name.c_str(), i));
    std::pair<float, float>  thispair;
    thispair.first  = h1_ped->GetMean();
    thispair.second = h1_ped->GetRMS();
    peds.push_back(thispair);
  }

  return peds;
  
}



float getMeanposHodo( std::vector<float> hodo, std::vector<std::pair<float,float> > pedestals ) {

  float mean = 0.;
  float eTot = 0.;
  for( unsigned i=0; i<hodo.size(); ++i ) {
    if( hodo[i] > (pedestals[i].first + 2.* pedestals[i].second) ) {
      mean += (i-3.5);
      eTot += 1.;
    }
  }

  return mean/eTot;

}


std::vector< std::pair<float, float> > getPedestalsHodo( const std::string type, const std::string fileName ) {

  TFile* file = TFile::Open(fileName.c_str());

  std::vector< std::pair<float, float> > peds;
  for( unsigned i=0; i<HODOX_CHANNELS; ++i ) {
    TH1D* h1_ped = (TH1D*)file->Get(Form("HODO%sRAW_hodo%sRawSpectrum_%d", type.c_str(), type.c_str(), i));
    std::pair<float, float>  thispair;
    thispair.first  = h1_ped->GetMean();
    thispair.second = h1_ped->GetRMS();
    peds.push_back(thispair);
  }

  return peds;

}


std::vector<float> subtractPedestals( std::vector<float> raw, std::vector< std::pair<float, float> > pedestals, float nSigma ) {

  std::vector<float> corr;

  for(unsigned i=0; i<raw.size(); ++i ) {

    float iCorr = ( raw[i] > pedestals[i].first + nSigma*pedestals[i].second ) ? (raw[i] - pedestals[i].first) : 0.;
    corr.push_back( iCorr );

  }

  return corr;

}


float sumVector( std::vector<float> v ) {

  float sum=0.;
  for( unsigned i=0; i<v.size(); ++i ) sum += v[i];

  return sum;

}


bool checkVector( std::vector<float> v, float theMax ) {

  bool returnBool = true;

  for( unsigned i=0; i<v.size(); ++i ) {
    if( v[i]<=0. ) return false;
    if( v[i]>=theMax ) return false;
  }

  return returnBool;

}

