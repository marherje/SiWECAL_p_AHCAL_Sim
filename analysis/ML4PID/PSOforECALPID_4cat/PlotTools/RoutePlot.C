#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <fstream>
#include <map>
#include "TChain.h"
#include "TGraph.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine.h"
#include "TAxis3D.h"
#include "TArrow.h"
#include "TMath.h"
#include "TPaveText.h"

#include "TFile.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include <string>
#include <stdlib.h>


void RoutePlot(){
gStyle->SetTitleX(0.21);
gStyle->SetTitleY(0.98);

  const int Iteration=5;
  int nParticles=15;
  
  TGraph* emptygraph = new TGraph(1);
  emptygraph->SetPoint(0,1000,0.005);
  
  std::vector< std::vector<TGraph*> > graphs;
  std::vector<TPolyLine*> lines;
  std::vector<TGraph*> FirstPoints;
  std::vector<TGraph*> LastPoints;
  std::vector< std::vector<TArrow*> > arrows;
  

  //get max an min ROC
  Int_t bestPartIteration[Iteration];
  Double_t bestROCIteration[Iteration];
  for(int o=0;o<Iteration;o++){bestROCIteration[o]=0.0;}
  Double_t maxROC=0.0;
  Double_t minROC=1.0;
  Double_t bestROCGlobal=0.0;
  for(int i=0;i<nParticles;i++){
      if(!(i==4 or i==5 or i==6 or i==13))continue;

//     std::cout<<i<<std::endl;
    std::string buffer="";
    std::stringstream buffer2;
    buffer2<<i;
//     std::cout<<buffer2.str()<<" "<<buffer2.str()<<std::endl;
//     buffer2>>buffer;
    std::string filename="PartForPlot/Particle";
    filename+=buffer2.str();
    filename+="/ParticleRoute.txt";
//     std::cout<<filename.c_str()<<std::endl;
    std::ifstream inputfile(filename.c_str());
//     std::ifstream inputfile("Particles/Particle1/ParticleRoute.txt");

    int k=0;
    double nTrees=0.0;
    double shrinkage=0.0;
    double bagging=0.0;
    double cuts=0.0;
    double ROC=0.0;
    double KS=0.0;
    double depth=0.0;
    bool readline=true;
    while(readline){
      inputfile>>ROC;
      inputfile>>KS;
      inputfile>>nTrees;
      inputfile>>shrinkage;
      inputfile>>bagging;
      inputfile>>cuts;
//       inputfile>>depth;
      std::string dump="";
      int j=0;
      do{
        inputfile>>dump;
//         std::cout<<dump<<std::endl;
        j++;
//         if(dump=="--Next--")std::cout<<"here"<<std::endl;
//         if(dump=="--Next--\n")std::cout<<"too"<<std::endl;
//         if(dump=="")std::cout<<"never"<<std::endl;

      }while(dump!="--Next--" and dump!="--Next--\n" and dump!="");
//       std::cout<<"Data "<<nTrees<<" "<<shrinkage<<" "<<ROC<<std::endl;
      if(dump!=""){
        if(ROC>maxROC)maxROC=ROC;
        if(ROC<minROC and ROC!=0)minROC=ROC;
      
      }
      k++;
//      inputfile>>dump;
      if(inputfile.eof() or dump=="" or k==Iteration)readline=false;
    }

    inputfile.close();
    std::cout<<"NPoints: "<<k<<std::endl;
  }
  std::cout<<minROC<<"  "<<maxROC<<std::endl;
  
  //get marker scaling;
  Double_t maxSize=2.4;
  Double_t minSize=0.3;
  minROC=0.820;
  Double_t m = (maxSize-minSize)/(maxROC-minROC);
  Double_t abschnitt = maxSize-m*maxROC;
  std::cout<<m<<" "<<abschnitt<<std::endl;
  std::cout<<minSize<<" "<<maxSize<<std::endl;

      
  //read points
  for(int i=0;i<nParticles;i++){
    if(!(i==4 or i==5 or i==6 or i==13))continue;

    std::cout<<i<<std::endl;
    std::string buffer="";
    std::stringstream buffer2;
    buffer2<<i;
//     std::cout<<buffer2.str()<<" "<<buffer2.str()<<std::endl;
//     buffer2>>buffer;
    std::string filename="PartForPlot/Particle";
    filename+=buffer2.str();
    filename+="/ParticleRoute.txt";
    std::cout<<filename.c_str()<<std::endl;
    std::ifstream inputfile(filename.c_str());
//     std::ifstream inputfile("Particles/Particle1/ParticleRoute.txt");

    std::vector< TGraph*> buffgraph;
    graphs.push_back( buffgraph);
    lines.push_back(new TPolyLine);
    FirstPoints.push_back(new TGraph);
    LastPoints.push_back(new TGraph);
    std::vector< TArrow*> buffarrow;
    arrows.push_back(buffarrow);


    int k=0;
    double nTrees=0.0;
    double shrinkage=0.0;
    double bagging=0.0;
    double cuts=0.0;
    double ROC=0.0;
    double KS=0.0;
    double depth=0.0;
    
    bool readline=true;
    while(readline){
      inputfile>>ROC;
      inputfile>>KS;
      inputfile>>nTrees;
      inputfile>>shrinkage;
      inputfile>>bagging;
      inputfile>>cuts;
//       inputfile>>depth;
      std::string dump="";
      int j=0;
      do{
        inputfile>>dump;
//         std::cout<<dump<<std::endl;
        j++;
//         if(dump=="--Next--")std::cout<<"here"<<std::endl;
//         if(dump=="--Next--\n")std::cout<<"too"<<std::endl;
//         if(dump=="")std::cout<<"never"<<std::endl;

      }while(dump!="--Next--" and dump!="--Next--\n" and dump!="");
//       std::cout<<"Data "<<nTrees<<" "<<shrinkage<<" "<<ROC<<std::endl;
      if(dump!=""){
        graphs.back().push_back(new TGraph);
        graphs.back().back()->SetPoint(0,nTrees,shrinkage);
        if(ROC==0)graphs.back().back()->SetMarkerStyle(3);
        else{
          graphs.back().back()->SetMarkerStyle(8);
//           std::cout<<TMath::Exp(ROC)<<std::endl;
          Double_t size=m*ROC+abschnitt;
          graphs.back().back()->SetMarkerSize(size); 
          std::cout<<ROC<<" "<<size<<std::endl;

        }
        if(k==0)graphs.back().back()->SetMarkerStyle(34);
        lines.back()->SetPoint(k,nTrees,shrinkage);
        
        
           
          
          //if(k==0)FirstPoints.back()->SetPoint(0,nTrees,shrinkage);
      
      }
      if(ROC>=bestROCIteration[k]){
      bestROCIteration[k]=ROC;
      bestPartIteration[k]=i;
//       if(ROC>=bestROCGlobal)bestROCGlobal=ROC;
//       if(bestROCGlobal>=bestROCIteration[k])bestROCIteration[k]=bestROCGlobal;
      }
      
      k++;
//      inputfile>>dump;
      if(inputfile.eof() or dump=="" or k==Iteration)readline=false;
    }
    
//     LastPoints.back()->SetPoint(0,nTrees,bagging);

    inputfile.close();
    std::cout<<"NPoints: "<<k<<std::endl;
    
    //get velocities
    for(size_t l=1;l<graphs.back().size()-1;l++){
//       std::cout<<graphs.size()<<std::endl;
      Double_t PrevTrees=0.001;
      Double_t PrevShrinkage=0.001;
      graphs.back().at(l-1)->GetPoint(0,PrevTrees,PrevShrinkage);
//       std::cout<<PrevTrees<<" "<<PrevShrinkage<<std::endl;
      Double_t CurrTrees;
      Double_t CurrShrinkage;
      graphs.back().at(l)->GetPoint(0,CurrTrees,CurrShrinkage);
      Double_t NextTrees;
      Double_t NextShrinkage;
      graphs.back().at(l+1)->GetPoint(0,NextTrees,NextShrinkage);
      Double_t velTree=CurrTrees-PrevTrees;
      Double_t velShrinkage=CurrShrinkage-PrevShrinkage;
      Double_t nvelTree=NextTrees-CurrTrees;
      Double_t nvelShrinkage=NextShrinkage-CurrShrinkage;
      Double_t forceTrees=nvelTree-0.5*velTree;
      Double_t forceShrinkage=nvelShrinkage-0.5*velShrinkage;
      arrows.back().push_back(new TArrow(CurrTrees,CurrShrinkage,CurrTrees+forceTrees,CurrShrinkage+forceShrinkage,0.02,">"));
      
    }
    
  }//end input loop
    
    //Get velocities
    
    
    TCanvas* c = new TCanvas("c","c",800,600);
    TLegend* leg= new TLegend(0.5,0.65,0.7,0.85);
    leg->SetFillColor(0);
    leg->SetTextFont(42);
    std::vector<TGraph*> leggraphs;

    for(size_t i=0;i<graphs.size();i++){
      leggraphs.push_back(new TGraph);
      leggraphs.back()->SetMarkerColor(1+i);
      leggraphs.back()->SetMarkerStyle(8);
      leggraphs.back()->SetMarkerSize(1);
      if(i==2)leggraphs.back()->SetMarkerColor(6);
      leg->AddEntry(leggraphs.back(), Form("particle %i",i+1),"p");
    }
    TLegend* legMarker= new TLegend(0.7,0.55,0.9,0.85);
    legMarker->SetFillColor(0);
    legMarker->SetFillStyle(0);
    legMarker->SetTextFont(42);
    TLegend* legCMS= new TLegend(0.07,0.75,0.39,0.85);
    legCMS->SetFillColor(0);
    legCMS->SetFillStyle(0);
    legCMS->SetBorderSize(0);
    legCMS->SetTextFont(42);
    legCMS->AddEntry((TObject*)0,"CMS private work","");
    TPaveText* legExplain= new TPaveText(0.23,0.84,0.73,0.94,"NDC");
    legExplain->SetFillColor(3);
    legExplain->SetFillStyle(0);
    legExplain->SetBorderSize(0);
    legExplain->AddText("BDTs trained and tested with t#bar{t}H, H#rightarrowb#bar{b} and t#bar{t} MC-Samples");
    legExplain->SetTextFont(42);
    
    TGraph* initMarker = new TGraph();
    initMarker->SetMarkerStyle(34);
    initMarker->SetMarkerSize(1.5);
    legMarker->AddEntry(initMarker,"init. pos.","p");
    
    TArrow* legArrow = new TArrow(100,0.01,200,0.02,0.02,">");
    TArrow* legArrow2 = new TArrow(1715,0.0334,1793,0.0334,0.02,">");
    
    legArrow->SetLineColor(0);
    legArrow->SetFillColor(0);
    legArrow2->SetLineWidth(2);
    legArrow2->SetAngle(40);
    legMarker->AddEntry(legArrow,"social pressure","");
//         legMarker->AddEntry((TObject*)0,"","");

    TGraph* sizegraph1 = new TGraph();
    sizegraph1->SetMarkerStyle(8);
    Double_t size2=m*0.82+abschnitt;
    sizegraph1->SetMarkerSize(size2);
    legMarker->AddEntry(sizegraph1, "A_{ROC} = 0.82","p");
    TGraph* sizegraph2 = new TGraph();
    sizegraph2->SetMarkerStyle(8);
    size2=m*0.83+abschnitt;
    sizegraph2->SetMarkerSize(size2);
    legMarker->AddEntry(sizegraph2, "A_{ROC} = 0.83","p");
    TGraph* sizegraph3 = new TGraph();
    sizegraph3->SetMarkerStyle(8);
    size2=m*0.84+abschnitt;
    sizegraph3->SetMarkerSize(size2);
    legMarker->AddEntry(sizegraph3, "A_{ROC} = 0.84","p");
    
    
    
    for(size_t i=0;i<graphs.size();i++){

      if(i==0){
        c->cd();
        
        
        emptygraph->Draw("AP");
        c->Update();
        emptygraph->GetXaxis()->SetLimits(100,2200);
        emptygraph->GetXaxis()->SetRangeUser(100,2200);

        emptygraph->GetYaxis()->SetLimits(0.0,0.038);
        emptygraph->GetYaxis()->SetRangeUser(0.0,0.038);
        emptygraph->SetTitle("particle paths in the n_{Trees} - shrinkage plane");
        emptygraph->GetXaxis()->SetTitle("n_{Trees}");
        emptygraph->GetYaxis()->SetTitle("shrinkage");
        emptygraph->GetYaxis()->SetTitleOffset(1.5);
        emptygraph->GetXaxis()->SetTitleOffset(1.2);
       std::cout<<emptygraph->GetXaxis()->GetTitleFont()<<std::endl;

	c->SetTopMargin(0.15);
        c->Update();
        for(size_t j=0;j<graphs.at(i).size();j++){
//           graphs.at(i).at(j)->SetMarkerStyle(3);
          graphs.at(i).at(j)->SetMarkerColor(1+i);
          graphs.at(i).at(j)->SetLineColor(1+i);
//         graphs.at(i)->GetXaxis()->SetLimits(200.0,2500.0);
//         graphs.at(i)->GetYaxis()->SetLimits(0.1,0.9);
//         graphs.at(i)->GetZaxis()->SetRangeUser(0.0,1.0);
          graphs.at(i).at(j)->Draw("P SAME");
        c->Update();
        }
        for(size_t j=0;j<arrows.at(i).size();j++){
          arrows.at(i).at(j)->SetLineColor(1+i);
          arrows.at(i).at(j)->SetFillColor(1+i);
          arrows.at(i).at(j)->SetAngle(40);
          arrows.at(i).at(j)->SetLineWidth(1.9);
          arrows.at(i).at(j)->Draw("");
          c->Update();
        }
        lines.at(i)->SetLineColor(1+i);
        lines.at(i)->SetLineStyle(7);
        lines.at(i)->Draw();
        c->Update();
        
//         FirstPoints.at(i)->SetMarkerStyle(20);
//         FirstPoints.at(i)->SetMarkerColor(1+i);
//         FirstPoints.at(i)->Draw("P SAME");
//         LastPoints.at(i)->SetMarkerStyle(21);
//         LastPoints.at(i)->SetMarkerColor(1+i);
//         LastPoints.at(i)->Draw("P SAME");
//         c->Update();
//         LastPoints.at(i)->GetXaxis()->SetRangeUser(200.0,2500.0);
//         LastPoints.at(i)->GetYaxis()->SetRangeUser(0.0001,0.05);
//         LastPoints.at(i)->GetZaxis()->SetRangeUser(0.0,1.0);
//         LastPoints.at(i)->Draw("P SAME");

//         graphs.at(i)->Draw("P SAME");
        c->Update();
//         TObject* view = c->GetView3D();
//         TAxis3D *axis = TAxis3D::GetPadAxis();
//         std::cout<<view<<std::endl;
//         TAxis3D::ToggleRulers();     // To pop axice down
//         axis->SetLabelColor(kBlue); // Paint the axice labels with blue color
//         axis->SetAxisColor(kRed);   // Paint the axice itself with blue color
//         TAxis3D::ToggleRulers(); 
//         axis->Paint();// To pop axice up
//         c->Update();
        
      }
      else{
        c->cd();
        for(size_t j=0;j<graphs.at(i).size();j++){
//           graphs.at(i).at(j)->SetMarkerStyle(3);
          graphs.at(i).at(j)->SetMarkerColor(1+i);
          graphs.at(i).at(j)->SetLineColor(1+i);
          if(i==2){
          graphs.at(i).at(j)->SetMarkerColor(6);
          graphs.at(i).at(j)->SetLineColor(6);
          }
//         graphs.at(i)->GetXaxis()->SetRangeUser(200.0,2500.0);
//         graphs.at(i)->GetYaxis()->SetRangeUser(0.0001,0.05);
//         graphs.at(i)->GetZaxis()->SetRangeUser(0.0,1.0);
          graphs.at(i).at(j)->Draw("SAME P");
        c->Update();
        }
        for(size_t j=0;j<arrows.at(i).size();j++){
        arrows.at(i).at(j)->SetLineColor(1+i);
        arrows.at(i).at(j)->SetFillColor(1+i);
        if(i==2){
        arrows.at(i).at(j)->SetLineColor(6);
        arrows.at(i).at(j)->SetFillColor(6);
        }
        arrows.at(i).at(j)->SetAngle(40);
        arrows.at(i).at(j)->SetLineWidth(1.9);


        arrows.at(i).at(j)->Draw("");
        c->Update();
        }
        lines.at(i)->SetLineColor(1+i);
        lines.at(i)->SetLineStyle(7);
        if(i==2)lines.at(i)->SetLineColor(6);

        lines.at(i)->Draw();
        c->Update();
        
// d        FirstPoints.at(i)->SetMarkerStyle(20);
//     d    FirstPoints.at(i)->SetMarkerColor(1+i);
//         FirstPoints.at(i)->Draw("P Same");
//         LastPoints.at(i)->SetMarkerStyle(21);
//         LastPoints.at(i)->SetMarkerColor(1+i);
//         LastPoints.at(i)->Draw("P Same");
        c->Update();
      }
      
    }
    legArrow2->Draw("");
    leg->Draw();
    legMarker->Draw();
    legCMS->Draw();
    legExplain->Draw();
    c->Update();
    TString outfile="ParticleRoute_";
    outfile+=Iteration;
    outfile+="eng.eps";
    c->SaveAs(outfile);
    
    
    for(int p=0;p<Iteration;p++){
    if(p>0 and bestROCIteration[p-1]>bestROCIteration[p])bestROCIteration[p]=bestROCIteration[p-1];

    std::cout<<p<<" "<<bestROCIteration[p]<<" "<<bestPartIteration[p]<<std::endl;
    
    }
}
