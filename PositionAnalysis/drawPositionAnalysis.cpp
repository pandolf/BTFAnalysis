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



TStyle* setStyle();
TGraphErrors* get_xyCenter( TH2D* h2_xyPos );


int main( int argc, char* argv[] ) {

  
  std::string runName = "test_10";
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


  // manually set beam nominal position for some known runs:
  float beamX = -999.;
  float beamY = -999.;
  float beamRX = 4.;
  float beamRY = 2.;
  if( runName == "BTF_94_20140430-073300_beam" ) {
    beamX = -3.;
    beamY = +3.;
  } else if( runName == "BTF_96_20140430-083733_beam" ) {
    beamX = -6.;
    beamY = +6.;
  } else if( runName == "BTF_98_20140430-092026_beam" ) {
    beamX = -9.;
    beamY = +9.;
  }

  bool drawBeam = ((beamX>-999.) && (beamY>-999.));
  //bool beamInsideHodo = ((fabs(beamX)<4.) && (fabs(beamY)<4.));


  TCanvas* c1 = new TCanvas("c1", "", 600, 600);
  c1->cd();

  TH2D* h2_xyPos = (TH2D*)file->Get("xyPos"); 
  TH2D* h2_xyPos_hodo = (TH2D*)file->Get("xyPos_hodo"); 
  TH2D* h2_xyPos_bgo = (TH2D*)file->Get("xyPos_bgo"); 

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




  TLegend* legend = new TLegend( 0.75, 0.21, 0.9, 0.39 );
  legend->SetFillColor(0);
  legend->SetTextSize(0.038);
  legend->AddEntry( gr_xyCenter, "CeF3", "P" );
  legend->AddEntry( gr_xyCenter_bgo, "BGO", "P" );
  legend->AddEntry( gr_xyCenter_hodo, "Hodo", "P" );
  legend->Draw("same");


  TGraph* gr_beamPos = new TGraph(0);
  gr_beamPos->SetMarkerStyle(24);
  gr_beamPos->SetMarkerSize(4);
  TLegend* legend2 = new TLegend( 0.18, 0.21, 0.52, 0.25 );
  legend2->SetFillColor(0);
  legend2->SetTextSize(0.035);
  legend2->AddEntry( gr_beamPos, "Beam Position", "P" );
  if( drawBeam ) 
    legend2->Draw("same");


  TPaveText* label_top = new TPaveText(0.4,0.953,0.975,0.975, "brNDC");
  label_top->SetFillColor(kWhite);
  label_top->SetTextSize(0.038);
  label_top->SetTextAlign(31); // align right
  label_top->SetTextFont(62);
  label_top->AddText("500 MeV Electron Beam");
  label_top->Draw("same");

  TPaveText* label_run = new TPaveText(0.34,0.86,0.9,0.91, "brNDC");
  label_run->SetFillColor(kWhite);
  label_run->SetTextSize(0.033);
  label_run->SetTextAlign(11); // align right
  label_run->SetTextFont(42);
  label_run->AddText(Form("Run %s", runName.c_str()));
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

  gr_xyCenter_hodo->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  gr_xyCenter->Draw("p same");



  c1->SaveAs(Form("%s/xyPos.eps", outputdir.c_str()) );
  c1->SaveAs(Form("%s/xyPos.png", outputdir.c_str()) );

  c1->Clear();

  xMax = xySize/2.;
  TH2D* h2_axes_zoom = new TH2D("axes_zoom", "", 10, -xMax, xMax, 10, -xMax, xMax);
  h2_axes_zoom->SetXTitle("X Position [mm]");
  h2_axes_zoom->SetYTitle("Y Position [mm]");
  h2_axes_zoom->Draw();

  legend->Draw("same");
  drawBeam = ((fabs(beamX)<xMax) && (fabs(beamY)<xMax));
  if( drawBeam )
    legend2->Draw("same");
  label_top->Draw("same");
  label_run->Draw("same");


  h2_xyPos_bgo->Draw("same");
  h2_xyPos_hodo->Draw("same");
  h2_xyPos->Draw("same");

  lineHodo_x1->Draw("same");
  lineHodo_x2->Draw("same");
  lineHodo_y1->Draw("same");
  lineHodo_y2->Draw("same");

  beamPos->Draw("same"); 

  gr_xyCenter_hodo->Draw("p same");
  gr_xyCenter_bgo->Draw("p same");
  gr_xyCenter->Draw("p same");


  c1->SaveAs(Form("%s/xyPos_zoom.eps", outputdir.c_str()) );
  c1->SaveAs(Form("%s/xyPos_zoom.png", outputdir.c_str()) );

  return 0;

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
