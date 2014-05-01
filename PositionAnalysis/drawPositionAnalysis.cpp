#include <iostream>
#include <string>
#include <stdlib.h>
#include <cmath>

#include "TFile.h"
#include "TH1D.h"
#include "TH2D.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "TGraphErrors.h"
#include "TEllipse.h"
#include "TString.h"



std::string runName;


TStyle* setStyle();
TGraphErrors* get_xyCenter( TH2D* h2_xyPos );
void drawSinglePositionPlot( const std::string& outputdir, TFile* file, const std::string& runName, const std::string& suffix );
void drawSinglePlot( const std::string& outputdir, const std::string& saveName, TFile* file, const std::string& name, const std::string& axisName, int nChannels, float xMin=0, float xMax=4095, int rebin=1, bool plotLog=false );
void fitHodoWithBeam( const std::string& outputdir, const std::string& suffix, TH1D* h1, float r, float& pos, float& pos_err );
TPaveText* getLabelTop();
TPaveText* getLabelRun( const std::string& runName, bool top=true );
void getBeamPosition( const std::string& runName, float& beamX, float& beamY );


int main( int argc, char* argv[] ) {

  
  runName = "test_10";
  if( argc>1 ) {
    std::string runName_str(argv[1]);
    runName = runName_str;
  }


  std::string fileName = "PosAn_" + runName + ".root";
  TFile* file = TFile::Open( fileName.c_str() );
  std::cout << "-> Opened file: " << fileName << std::endl;
  
  std::string outputdir = "Plots_" + runName;
  std::string mkdir_command = "mkdir -p " + outputdir;
  system(mkdir_command.c_str());


  TStyle* style = setStyle();
  style->cd();

  drawSinglePositionPlot( outputdir, file, runName, "" );
  drawSinglePositionPlot( outputdir, file, runName, "_singleEle" );

  //drawSinglePlot( outputdir, "cef3_spectrum"      , file, "cef3"     , "ADC Counts", 4, 0., 3500., 10, true );
  //drawSinglePlot( outputdir, "cef3_corr_spectrum" , file, "cef3_corr", "ADC Counts", 4, 0., 3500., 10, true );

  drawSinglePlot( outputdir, "cef3_spectrum_lin"      , file, "cef3"     , "ADC Counts", 4, 0., 3500., 10, false );
  drawSinglePlot( outputdir, "cef3_corr_spectrum_lin" , file, "cef3_corr", "ADC Counts", 4, 0., 3500., 10, false );

  return 0;

}






void drawSinglePositionPlot( const std::string& outputdir, TFile* file, const std::string& runName, const std::string& suffix ) {

  // manually set beam nominal position for some known runs:
  float beamX = -999.;
  float beamY = -999.;
  float beamRX = 4.;
  float beamRY = 2.;
  getBeamPosition( runName, beamX, beamY );


  bool drawBeam = ((beamX>-999.) && (beamY>-999.));
  //bool beamInsideHodo = ((fabs(beamX)<4.) && (fabs(beamY)<4.));


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_xyPos      = (TH2D*)file->Get(Form("xyPos%s", suffix.c_str())); 
  TH2D* h2_xyPos_hodo = (TH2D*)file->Get(Form("xyPos%s_hodo", suffix.c_str())); 
  TH2D* h2_xyPos_bgo  = (TH2D*)file->Get(Form("xyPos%s_bgo", suffix.c_str())); 

  float xySize = 25.;
  float xMax = xySize*3./2.;
  TH2D* h2_axes = new TH2D("axes", "", 10, -xMax, xMax, 10, -xMax, xMax);
  h2_axes->SetXTitle("X Position [mm]");
  h2_axes->SetYTitle("Y Position [mm]");

  h2_axes->Draw();


  TGraphErrors* gr_xyCenter      = get_xyCenter( h2_xyPos );
  TGraphErrors* gr_xyCenter_hodo = get_xyCenter( h2_xyPos_hodo );
  TGraphErrors* gr_xyCenter_bgo  = get_xyCenter( h2_xyPos_bgo  );

  gr_xyCenter->SetMarkerColor(kRed+2);
  gr_xyCenter->SetLineColor(kRed+2);
  gr_xyCenter->SetMarkerStyle(20);
  gr_xyCenter->SetMarkerSize(1.6);

  gr_xyCenter_hodo->SetMarkerColor(kBlack);
  gr_xyCenter_hodo->SetLineColor(kBlack);
  gr_xyCenter_hodo->SetMarkerStyle(20);
  gr_xyCenter_hodo->SetMarkerSize(1.6);

  gr_xyCenter_bgo->SetMarkerColor(kGreen+3);
  gr_xyCenter_bgo->SetLineColor(kGreen+3);
  gr_xyCenter_bgo->SetMarkerStyle(20);
  gr_xyCenter_bgo->SetMarkerSize(1.6);


  // try to fit the hodo points with the expected beam size
  float xPos_hodo_fit, xPos_hodo_fit_err;
  float yPos_hodo_fit, yPos_hodo_fit_err;
  fitHodoWithBeam( outputdir, "X"+suffix, h2_xyPos_hodo->ProjectionX(), beamRX, xPos_hodo_fit, xPos_hodo_fit_err );
  fitHodoWithBeam( outputdir, "Y"+suffix, h2_xyPos_hodo->ProjectionY(), beamRY, yPos_hodo_fit, yPos_hodo_fit_err );
  //TFile* file_prova = TFile::Open("prova.root", "recreate");
  //file_prova->cd();
  //h2_xyPos_hodo->Write();
  //h2_xyPos_hodo->Write();
  //file_prova->Close();
  //exit(11);

  TGraphErrors* gr_xyCenter_hodo_fit = new TGraphErrors(0);
  gr_xyCenter_hodo_fit->SetPoint(0, xPos_hodo_fit, yPos_hodo_fit);
  gr_xyCenter_hodo_fit->SetPointError(0, xPos_hodo_fit_err, yPos_hodo_fit_err);
  gr_xyCenter_hodo_fit->SetMarkerStyle(20);
  gr_xyCenter_hodo_fit->SetMarkerSize(1.6);




  TGraph* gr_beamPos = new TGraph(0);
  gr_beamPos->SetMarkerStyle(24);
  gr_beamPos->SetMarkerSize(4);
  TLegend* legend2 = new TLegend( 0.18, 0.21, 0.52, 0.25 );
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry( gr_beamPos, "Beam Position", "P" );
  if( drawBeam ) 
    legend2->Draw("same");


  TPaveText* label_top = getLabelTop();
  TPaveText* label_run = getLabelRun(runName);
  label_top->Draw("same");
  label_run->Draw("same");


  int lineColor = 17;

  TLine* line_x1 = new TLine( -xMax, -xySize/2., +xMax, -xySize/2. );
  line_x1->SetLineColor(lineColor);
  line_x1->Draw("same");

  TLine* line_x2 = new TLine( -xMax, +xySize/2., +xMax, +xySize/2. );
  line_x2->SetLineColor(lineColor);
  line_x2->Draw("same");

  TLine* line_y1 = new TLine( -xySize/2., -xMax, -xySize/2., +xMax );
  line_y1->SetLineColor(lineColor);
  line_y1->Draw("same");

  TLine* line_y2 = new TLine( +xySize/2., -xMax, +xySize/2., +xMax );
  line_y2->SetLineColor(lineColor);
  line_y2->Draw("same");


  TLegend* legend = new TLegend( 0.75, 0.21, 0.9, 0.39 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( gr_xyCenter, "CeF3", "P" );
  legend->AddEntry( gr_xyCenter_bgo, "BGO", "P" );
  legend->AddEntry( gr_xyCenter_hodo_fit, "Hodo", "P" );
  legend->Draw("same");


  h2_xyPos_hodo->SetMarkerColor(14);
  h2_xyPos->SetMarkerColor(46);
  h2_xyPos_bgo->SetMarkerColor(30);

  float hodoSize = 8.;
  TLine* lineHodo_x1 = new TLine( -hodoSize/2., -hodoSize/2., +hodoSize/2., -hodoSize/2. );
  lineHodo_x1->SetLineColor(kBlack);
  lineHodo_x1->SetLineStyle(2);
  lineHodo_x1->Draw("same");

  TLine* lineHodo_x2 = new TLine( -hodoSize/2., +hodoSize/2., +hodoSize/2., +hodoSize/2. );
  lineHodo_x2->SetLineColor(kBlack);
  lineHodo_x2->SetLineStyle(2);
  lineHodo_x2->Draw("same");

  TLine* lineHodo_y1 = new TLine( -hodoSize/2., -hodoSize/2., -hodoSize/2., +hodoSize/2. );
  lineHodo_y1->SetLineColor(kBlack);
  lineHodo_y1->SetLineStyle(2);
  lineHodo_y1->Draw("same");

  TLine* lineHodo_y2 = new TLine( +hodoSize/2., -hodoSize/2., +hodoSize/2., +hodoSize/2. );
  lineHodo_y2->SetLineColor(kBlack);
  lineHodo_y2->SetLineStyle(2);
  lineHodo_y2->Draw("same");


  h2_xyPos_bgo->Draw("same");
  h2_xyPos->Draw("same");
  h2_xyPos_hodo->Draw("same");


  TEllipse* beamPos = new TEllipse( beamX, beamY, beamRX, beamRY );
  beamPos->SetLineColor(kBlack);
  beamPos->SetFillStyle(0);
  beamPos->Draw("same"); 

  //gr_xyCenter_hodo->Draw("p same");  // now using hodo_fit
  gr_xyCenter_hodo_fit->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  gr_xyCenter->Draw("p same");


  c1->SaveAs(Form("%s/xyPos%s.eps", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s.png", outputdir.c_str(), suffix.c_str()) );

  c1->Clear();

  xMax = xySize/2.;
  TH2D* h2_axes_zoom = new TH2D("axes_zoom", "", 10, -xMax, xMax, 10, -xMax, xMax);
  h2_axes_zoom->SetXTitle("X Position [mm]");
  h2_axes_zoom->SetYTitle("Y Position [mm]");
  h2_axes_zoom->Draw();

  drawBeam = ((fabs(beamX)<xMax) && (fabs(beamY)<xMax));
  if( drawBeam )
    legend2->Draw("same");
  label_top->Draw("same");
  label_run->Draw("same");
  legend->Draw("same");

  h2_xyPos_bgo->Draw("same");
  h2_xyPos_hodo->Draw("same");
  h2_xyPos->Draw("same");


  lineHodo_x1->Draw("same");
  lineHodo_x2->Draw("same");
  lineHodo_y1->Draw("same");
  lineHodo_y2->Draw("same");

  beamPos->Draw("same"); 

  //gr_xyCenter_hodo->Draw("p same"); // now using hodo_fit
  gr_xyCenter_hodo_fit->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  gr_xyCenter->Draw("p same");


  c1->SaveAs(Form("%s/xyPos%s_zoom.eps", outputdir.c_str(), suffix.c_str()) );
  c1->SaveAs(Form("%s/xyPos%s_zoom.png", outputdir.c_str(), suffix.c_str()) );

  delete c1;
  delete h2_axes;
  delete h2_axes_zoom;
  delete legend;

}




void drawSinglePlot( const std::string& outputdir, const std::string& saveName, TFile* file, const std::string& name, const std::string& axisName, int nChannels, float xMin, float xMax, int rebin, bool plotLog ) {


  std::vector<int> colors;
  colors.push_back( 46 );
  colors.push_back( 38 );
  colors.push_back( 30 );
  colors.push_back( 42 );
  colors.push_back( 29 );
  colors.push_back( kBlack );
  colors.push_back( kGreen );
  colors.push_back( kBlue  );

  float yMax_leg = 0.8;
  float yMin_leg = yMax_leg-0.05*nChannels;
  TLegend* legend = new TLegend( 0.6, yMin_leg, 0.9, yMax_leg );
  legend->SetTextSize(0.035);
  legend->SetFillColor(0);

  float yMax;
  std::vector<TH1D*> histos;

  for(unsigned i=0; i<nChannels; ++i ) {

    TH1D* h1 = (TH1D*)file->Get(Form("%s_%d", name.c_str(), i));
    h1->SetLineColor( colors[i] );
    h1->SetLineWidth( 2 );
    h1->Rebin(rebin);
  
    histos.push_back(h1);

    legend->AddEntry( h1, Form("Channel %d", i), "L" );


    float thisMax = h1->GetMaximum();
    if( thisMax>yMax ) yMax = thisMax;

  }

  if( !plotLog ) yMax = 3000.;

  TCanvas* c1 = new TCanvas( "c1_new", "", 600, 600 );
  c1->cd();
  if( plotLog ) c1->SetLogy();


  float yScaleFactor = 1.3;
  //float yScaleFactor = (plotLog) ? 5.  : 1.3;
  float yMin         = (plotLog) ? 10. :  0.;

  TH2D* h2_axes = new TH2D("axes_new", "", 10, xMin, xMax, 10, yMin, yMax*yScaleFactor );
  h2_axes->SetXTitle( axisName.c_str() );
  h2_axes->SetYTitle( "Entries" );

  h2_axes->Draw();

  legend->Draw("Same");

  TPaveText* label_top = getLabelTop();
  TPaveText* label_run = getLabelRun(runName, !plotLog);
  label_top->Draw("same");
  label_run->Draw("same");

  for(unsigned i=0; i<nChannels; ++i ) 
    histos[i]->Draw("same");

  gPad->RedrawAxis();

  c1->SaveAs( Form("%s/%s.eps", outputdir.c_str(), saveName.c_str()) );
  c1->SaveAs( Form("%s/%s.png", outputdir.c_str(), saveName.c_str()) );


  c1->Clear();

  // and now draw also individual channels:
  for( unsigned i=0; i<nChannels; ++i ) {
  
    h2_axes->Draw();
    label_top->Draw("Same");
    label_run->Draw("Same");

    histos[i]->SetLineColor(46);
    histos[i]->Draw("same");

    TPaveText* labelChannel = new TPaveText( 0.65, 0.8, 0.9, 0.85, "brNDC" );
    labelChannel->SetTextSize( 0.035 );
    labelChannel->SetFillColor( 0 );
    labelChannel->AddText( Form("Channel %d", i) );
    labelChannel->Draw("same");
  
    c1->SaveAs( Form("%s/%s_%d.eps", outputdir.c_str(), saveName.c_str(), i) );
    c1->SaveAs( Form("%s/%s_%d.png", outputdir.c_str(), saveName.c_str(), i) );

  }


  delete c1;
  delete h2_axes;
  delete legend;

}









TStyle* setStyle() {

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

  return style_;

}


TGraphErrors* get_xyCenter( TH2D* h2_xyPos ) {

  float x     = h2_xyPos->ProjectionX()->GetMean();
  float x_err = h2_xyPos->ProjectionX()->GetRMS();

  float y     = h2_xyPos->ProjectionY()->GetMean();
  float y_err = h2_xyPos->ProjectionY()->GetRMS();

  TGraphErrors* gr_point = new TGraphErrors(0);
  gr_point->SetPoint( 0, x, y );
  gr_point->SetPointError( 0, x_err, y_err );

  return gr_point;

}





TPaveText* getLabelTop() {

  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(62);
  label_top->AddText("500 MeV Electron Beam");

  return label_top;

}


TPaveText* getLabelRun( const std::string& runName, bool top ) {

  float yMin = (top) ? 0.86 : 0.18;
  float yMax = (top) ? 0.91 : 0.23;
  TPaveText* label_run = new TPaveText(0.34,yMin,0.9,yMax, "brNDC");
  label_run->SetFillColor(kWhite);
  label_run->SetTextSize(0.033);
  label_run->SetTextAlign(11); // align right
  label_run->SetTextFont(42);
  label_run->AddText(Form("Run %s", runName.c_str()));

  return label_run;

}



void fitHodoWithBeam( const std::string& outputdir, const std::string& suffix, TH1D* h1, float r, float& pos, float& pos_err ) {

  h1->Rebin(5);
  float nentries =  h1->GetEntries();

  TF1* f1_gaus = new TF1( "gaus_hodo", "gaus" );
  f1_gaus->SetRange(-4., 4.);
  //f1_gaus->SetParameter(0, nentries);
  f1_gaus->FixParameter(0, nentries);
  f1_gaus->SetParameter(1, 0.);

  //f1_gaus->SetParLimits(0, 0.01*nentries, nentries);
  f1_gaus->SetParLimits(1, -40., 40.);
  f1_gaus->FixParameter(2, r );

  h1->Fit(f1_gaus, "RNL" );
  TCanvas* c1 = new TCanvas("c1_temp", "", 600, 600);
  c1->cd();
  h1->SetLineColor(kRed);
  h1->SetLineWidth(2);
  h1->SetXTitle("Position [mm]");
  h1->SetYTitle("Hits");
  h1->Draw("");
  f1_gaus->SetRange(-40., 40.);
  f1_gaus->Draw("same");
  c1->SaveAs(Form("%s/tmpFit%s.eps", outputdir.c_str(), suffix.c_str()));
  c1->SaveAs(Form("%s/tmpFit%s.png", outputdir.c_str(), suffix.c_str()));
  delete c1;

  pos = f1_gaus->GetParameter(1);
  pos_err = f1_gaus->GetParameter(2);

  delete f1_gaus;

}





void getBeamPosition( const std::string& runName, float& beamX, float& beamY ) {

  TString runName_tstr(runName.c_str());


  if( runName == "BTF_94_20140430-073300_beam" ) {
    beamX = -3.;
    beamY = +3.;
  } else if( runName == "BTF_96_20140430-083733_beam" ) {
    beamX = -6.;
    beamY = +6.;
  } else if( runName == "BTF_98_20140430-092026_beam" ) {
    beamX = -9.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_144_2014") ) {
    beamX = +12.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_145_2014") ) {
    beamX = +10.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_146_2014") ) {
    beamX = +8.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_147_2014") ) {
    beamX = +6.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_148_2014") ) {
    beamX = +4.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_149_2014") ) {
    beamX = +2.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_150_2014") ) {
    beamX = +0.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_151_2014") ) {
    beamX = -2.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_152_2014") ) {
    beamX = -4.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_153_2014") ) {
    beamX = -6.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_154_2014") ) {
    beamX = -8.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_155_2014") ) {
    beamX = -10.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_156_2014") ) {
    beamX = -12.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_157_2014") ) {
    beamX = +0.;
    beamY = +8.;
  } else if( runName_tstr.BeginsWith("BTF_158_2014") ) {
    beamX = +0.;
    beamY = +6.;
  } else if( runName_tstr.BeginsWith("BTF_159_2014") ) {
    beamX = +0.;
    beamY = +4.;
  } else if( runName_tstr.BeginsWith("BTF_160_2014") ) {
    beamX = +0.;
    beamY = +2.;
  } else if( runName_tstr.BeginsWith("BTF_161_2014") ) {
    beamX = +0.;
    beamY = +0.;
  } else if( runName_tstr.BeginsWith("BTF_162_2014") ) {
    beamX = +0.;
    beamY = -2.;
  } else if( runName_tstr.BeginsWith("BTF_163_2014") ) {
    beamX = +0.;
    beamY = -4.;
  } else if( runName_tstr.BeginsWith("BTF_164_2014") ) {
    beamX = +0.;
    beamY = -6.;
  } else if( runName_tstr.BeginsWith("BTF_165_2014") ) {
    beamX = +0.;
    beamY = -8.;
  } else {
    beamX = -999.;
    beamY = -999.;
  }



}

