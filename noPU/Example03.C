//
// Plots one waveform. To run it:
// > root -l Example03.C
//


#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TRandom.h>
#include <iostream>


void Example03()
{

  TCanvas *c = new TCanvas("c", "c", 300,300);
  //TFile *file = new TFile("data/waveform_signal_10GeV_eta_0.0_pu_140.root");
  //TFile *file = new TFile("output_example2_noPU.root");
  //TFile *file = new TFile("outputPSWithPU.root");
  TFile *file = new TFile("output_10.000000.root");
  //TFile *file = new TFile("output_example2_PU50.root");

  int    nWF;
  double waveform[500];
  TTree *tree = (TTree*)file->Get("Waveforms");
  tree->SetBranchAddress("nWF",      &nWF);
  tree->SetBranchAddress("waveform", waveform);

  //tree->GetEntry(10);
  tree->GetEntry(0);
  
  TGraph *gr = new TGraph();
  for(int i=0; i<nWF; i++){
    gr->SetPoint(i, i, waveform[i]);
  }

  //TCanvas *c = new TCanvas("c",500,500);
  gr->SetLineWidth(3);
  gr->Draw("AL");

  c->SaveAs("plot_ps.C");
  c->SaveAs("plot_ps.png");

  /*
  ///Draw another one 
  TFile *file1 = new TFile("outputPS.root");

  TTree *tree1 = (TTree*)file1->Get("Waveforms");
  tree1->SetBranchAddress("nWF",      &nWF);
  tree1->SetBranchAddress("waveform", waveform);

  tree1->GetEntry(10);
  
  TGraph *gr1 = new TGraph();
  for(int i=0; i<nWF; i++){
    gr1->SetPoint(i, i, waveform[i]);
  }
  
  gr1->SetLineWidth(3);
  gr1->SetLineColor(2);
  gr1->Draw("Lsame");

  c->Modified();
  c->Update();
  c->SaveAs("ps_withPU_noPU.png");
  */

}
