//
//  Style, based on a style file from BaBar
//

#include <iostream>
#include "Style.h"
#include "TROOT.h"

void SetQQbarStyle ()
{
  static TStyle* qqbarStyle = 0;
  //  std::cout << "\nApplying QQBAR style settings...\n" << std::endl ;
  if ( qqbarStyle==0 ) qqbarStyle = QQbarStyle();
  gROOT->SetStyle("QQBAR");
  gROOT->ForceStyle();
}

TStyle* QQbarStyle() 
{
  TStyle *qqbarStyle = new TStyle("QQBAR","QQbar style");

  // use plain black on white colors
  Int_t icol=0; // WHITE
  qqbarStyle->SetFrameBorderMode(icol);
  qqbarStyle->SetFrameFillColor(icol);
  qqbarStyle->SetCanvasBorderMode(icol);
  qqbarStyle->SetCanvasColor(icol);
  qqbarStyle->SetPadBorderMode(icol);
  qqbarStyle->SetPadColor(icol);
  qqbarStyle->SetStatColor(icol);
  //qqbarStyle->SetFillColor(icol); // don't use: white fill color for *all* objects

  // set the paper & margin sizes
  qqbarStyle->SetPaperSize(20,26);

  // set margin sizes
  qqbarStyle->SetPadTopMargin(0.05);
  qqbarStyle->SetPadRightMargin(0.05);
  qqbarStyle->SetPadBottomMargin(0.16);
  qqbarStyle->SetPadLeftMargin(0.16);

  // set title offsets (for axis label)
  //  qqbarStyle->SetTitleXOffset(1.4);
  //qqbarStyle->SetTitleYOffset(1.4);

  // use large fonts
  //Int_t font=72; // Helvetica italics
  Int_t font=42; // Helvetica
  Double_t tsize=0.04;
  qqbarStyle->SetTextFont(font);

  qqbarStyle->SetTitleX(0.2);
  qqbarStyle->SetTitleY(0.92);
  
  qqbarStyle->SetLabelFont(42,"xyz");
  qqbarStyle->SetLabelSize(0.04,"xyz");
  qqbarStyle->SetLabelOffset(0.015,"xyz");
  qqbarStyle->SetTitleFont(42,"xyz");
  qqbarStyle->SetTitleSize(0.05,"xyz");
  qqbarStyle->SetTitleOffset(1.35,"yz");
  qqbarStyle->SetTitleOffset(1.1,"x");
  
  /*qqbarStyle->SetTextSize(tsize);
  qqbarStyle->SetLabelFont(font,"x");
  qqbarStyle->SetTitleFont(font,"x");
  qqbarStyle->SetLabelFont(font,"y");
  qqbarStyle->SetTitleFont(font,"y");
  qqbarStyle->SetLabelFont(font,"z");
  qqbarStyle->SetTitleFont(font,"z");
  
  qqbarStyle->SetLabelSize(tsize,"x");
  qqbarStyle->SetTitleSize(tsize,"x");
  qqbarStyle->SetLabelSize(tsize,"y");
  qqbarStyle->SetTitleSize(tsize,"y");
  qqbarStyle->SetLabelSize(tsize,"z");
  qqbarStyle->SetTitleSize(tsize,"z");
  */
  
  // use bold lines and markers
  qqbarStyle->SetMarkerStyle(20);
  qqbarStyle->SetMarkerSize(0.8);
  qqbarStyle->SetHistLineWidth(2.);
  qqbarStyle->SetLineStyleString(2,"[12 12]"); // postscript dashes

  // get rid of X error bars 
  //qqbarStyle->SetErrorX(0.001);
  // get rid of error bar caps
  qqbarStyle->SetEndErrorSize(0.);

  // do not display any of the standard histogram decorations
  qqbarStyle->SetOptTitle(0);
  //qqbarStyle->SetOptStat(1111);
  qqbarStyle->SetOptStat(0);
  //qqbarStyle->SetOptFit(1111);
  qqbarStyle->SetOptFit(0);

  // put tick marks on top and RHS of plots
  qqbarStyle->SetPadTickX(1);
  qqbarStyle->SetPadTickY(1);

  return qqbarStyle;

}

