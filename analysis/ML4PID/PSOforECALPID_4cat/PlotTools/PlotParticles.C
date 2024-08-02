#include <vector>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <fstream>
#include <fstream>
#include <map>
#include "TChain.h"
#include "TGraph2D.h"
#include "TH2F.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TStopwatch.h"
#include "TLegend.h"
#include "TLatex.h"
#include "TLine.h"
#include "TPolyLine3D.h"
#include "TAxis3D.h"

#include "TFile.h"
#include "TVector3.h"
#include "TVector2.h"
#include "TLorentzVector.h"
#include <string>
#include <stdlib.h>


void PlotParticles(){
  int nParticles=15;
  
  std::vector<TGraph2D*> graphs;
  std::vector<TPolyLine3D*> lines;
  std::vector<TGraph2D*> FirstPoints;
  std::vector<TGraph2D*> LastPoints;

  for(int i=0;i<nParticles;i++){
    std::cout<<i<<std::endl;
    std::string buffer="";
    std::stringstream buffer2;
    buffer2<<i;
    std::cout<<buffer2.str()<<" "<<buffer2.str()<<std::endl;
//     buffer2>>buffer;
    std::string filename="Particles/Particle";
    filename+=buffer2.str();
    filename+="/ParticleRoute.txt";
    std::cout<<filename.c_str()<<std::endl;
    std::ifstream inputfile(filename.c_str());
//     std::ifstream inputfile("Particles/Particle1/ParticleRoute.txt");

    graphs.push_back(new TGraph2D);
    lines.push_back(new TPolyLine3D);
    FirstPoints.push_back(new TGraph2D);
    LastPoints.push_back(new TGraph2D);


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
        std::cout<<dump<<std::endl;
        j++;
        if(dump=="--Next--")std::cout<<"here"<<std::endl;
        if(dump=="--Next--\n")std::cout<<"too"<<std::endl;
        if(dump=="")std::cout<<"never"<<std::endl;

      }while(dump!="--Next--" and dump!="--Next--\n" and dump!="");
      std::cout<<"Data "<<nTrees<<" "<<shrinkage<<" "<<ROC<<std::endl;
      if(dump!=""){graphs.back()->SetPoint(k,nTrees,shrinkage,ROC);
        lines.back()->SetPoint(k,nTrees,shrinkage,ROC);
        if(k==0)FirstPoints.back()->SetPoint(0,nTrees,shrinkage,ROC);
      
      }
      k++;
//      inputfile>>dump;
      if(inputfile.eof() or dump=="")readline=false;
    }
    
    LastPoints.back()->SetPoint(0,nTrees,shrinkage,ROC);

    inputfile.close();
    std::cout<<"NPoints: "<<k<<std::endl;
    }
    
    TCanvas* c = new TCanvas("c","c",800,600);

    for(size_t i=0;i<graphs.size();i++){

      if(i==0){
        c->cd();
        
        graphs.at(i)->SetMarkerStyle(3);
        graphs.at(i)->SetMarkerColor(1+i);
        graphs.at(i)->SetLineColor(1+i);
        graphs.at(i)->GetXaxis()->SetLimits(200.0,2500.0);
        graphs.at(i)->GetYaxis()->SetLimits(0.0001,0.05);
        graphs.at(i)->GetZaxis()->SetRangeUser(0.0,1.0);
        graphs.at(i)->Draw("P");
        c->Update();
        
        lines.at(i)->SetLineColor(1+i);
        lines.at(i)->Draw();
        c->Update();
        FirstPoints.at(i)->SetMarkerStyle(20);
        FirstPoints.at(i)->SetMarkerColor(1+i);
        FirstPoints.at(i)->Draw("P SAME");
        LastPoints.at(i)->SetMarkerStyle(21);
        LastPoints.at(i)->SetMarkerColor(1+i);
        LastPoints.at(i)->Draw("P SAME");
        c->Update();
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

        graphs.at(i)->SetMarkerStyle(3);
        graphs.at(i)->SetMarkerColor(1+i);
        graphs.at(i)->SetLineColor(1+i);
//         graphs.at(i)->GetXaxis()->SetRangeUser(200.0,2500.0);
//         graphs.at(i)->GetYaxis()->SetRangeUser(0.0001,0.05);
//         graphs.at(i)->GetZaxis()->SetRangeUser(0.0,1.0);
        graphs.at(i)->Draw("SAME P");
        c->Update();
        lines.at(i)->SetLineColor(1+i);
        lines.at(i)->Draw();
        c->Update();
        FirstPoints.at(i)->SetMarkerStyle(20);
        FirstPoints.at(i)->SetMarkerColor(1+i);
        FirstPoints.at(i)->Draw("P Same");
        LastPoints.at(i)->SetMarkerStyle(21);
        LastPoints.at(i)->SetMarkerColor(1+i);
        LastPoints.at(i)->Draw("P Same");
        c->Update();
      }
      
    }
  
}