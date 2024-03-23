#include "Labels.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPave.h"
#include "TMarker.h"

/*void QQBARLabel(Double_t x,Double_t y,TString text,Color_t color)
{
  TLatex l;
  //l.SetTextAlign(12);
  l.SetTextSize(0.065);
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextColor(color);

  double delx = 1200*gPad->GetWh()/(1000*gPad->GetWw());
  
  l.DrawLatex(x,y,"ILD");
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextSize(0.05);
    p.SetTextFont(42);
    p.SetTextColor(color);
    p.DrawLatex(x+0.1,y,text);
  }
  }*/

void QQBARLabel(Double_t x,Double_t y,TString text,Color_t color)
{
  TLatex l;
  //l.SetTextAlign(12);
  l.SetTextSize(0.065);
  l.SetNDC();
  l.SetTextFont(42);
  l.SetTextColor(color);

  double delx = 1200*gPad->GetWh()/(1000*gPad->GetWw());
  
  l.DrawLatex(x,y,"ILD");
  if (text) {
    TLatex p;
    p.SetNDC();
    p.SetTextSize(0.03);
    p.SetTextFont(50);
    p.SetTextColor(kRed-3);
    p.DrawLatex(x-0.132,y-0.03,text);
  }
  }

void QQBARLabel2(Double_t x,Double_t y,TString text,Color_t color, Double_t textsize, Double_t angle)
{
  
  TLatex p;
  p.SetNDC();
  p.SetTextSize(textsize);
  p.SetTextFont(52);
  p.SetTextColor(color);
  p.SetTextAngle(angle);
  p.DrawLatex(x,y,text);
  
}

void Labels(TString energy, TString lumi, TString pol, TString cat){
  // 250: 20 y 36
  // 500 old: 47, 47
  QQBARLabel(0.86,0.952,"");
  if(pol=="eL_pR"){
    QQBARLabel2(0.1,0.02, "e_{L}^{-}e_{R}^{+} #rightarrow q#bar{q}, q=udscb, "+energy+" GeV, "+lumi+" fb^{-1}, cat. "+cat,kGray+2);
  }
  if(pol=="eR_pL"){
    QQBARLabel2(0.1,0.02, "e_{R}^{-}e_{L}^{+} #rightarrow q#bar{q}, q=udscb, "+energy+" GeV, "+lumi+" fb^{-1}, cat. "+cat,kGray+2);
  }
}



