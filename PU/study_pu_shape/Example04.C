//
// Takes waveforms from the "file".
// Creates NSAMPLES starting with IDSTART time with the step of NFREQ ns
// Applies noise (correlated) for each sample
// Stores samples and true in-time amplitude
// To run:
// > root -l -q Example04.C+
//


#include "Pulse.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TProfile.h>
#include <TF1.h>
#include <TGraph.h>
#include <TRandom.h>
#include <TMath.h>
#include <iostream>



void Example04()
{
  Pulse pSh;

  //TString filenameOutput = "output.root"; 
  //TString filenameOutput = "/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/outputPSWithPU_withNoise.root"; 
  //TString filenameOutput = Form("/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/noise_0p15/0p2GeV/outputPSWithPU_withNoise.root"); 
  TString filenameOutput = Form("/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/noise_0p15/0p15GeV/outputPSWithPU_withNoise.root"); 

  cout<<"output file name is "<<filenameOutput.Data()<<endl;
  // Noise level (GeV)
  //double sigmaNoise = 0.044;
  //double sigmaNoise = 0.060;
  double sigmaNoise = 0.15;
  

  // input Waveforms

  //TFile *file = new TFile("data/waveform_signal_10GeV_pu_0.root");
  //TFile *file = new TFile("output_example2_noPU.root");
  //TString fileName = Form("/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/noise_0p06/1GeV/outputPSWithPU.root");
  TString fileName = Form("/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/outputPSWithPU_0.150000.root");
  TFile *file = new TFile(fileName);

  cout<<"input file name is "<<fileName.Data()<<endl;

  int    BX0;
  int    nWF;
  double waveform[WFLENGTH], PUOnlywaveform[WFLENGTH];
  double energyPU[NBXTOTAL];
  double signalTruth;
  TTree *tree = (TTree*)file->Get("Waveforms");
  tree->SetBranchAddress("nWF",          &nWF);
  tree->SetBranchAddress("waveform",     waveform);
  tree->SetBranchAddress("PUOnlywaveform",     PUOnlywaveform);
  tree->SetBranchAddress("BX0",          &BX0);
  tree->SetBranchAddress("signalTruth",  &signalTruth);
  tree->SetBranchAddress("energyPU" ,    energyPU);



  // output samples
  
  int nSmpl = NSAMPLES;
  int nFreq = NFREQ;
  double samples[NSAMPLES], samplesNoise[NSAMPLES], PUOnlywaveform_[NSAMPLES];
  double ysamples[1], ysamplesNoise[1];

  double amplitudeTruth;
  TFile *fileOut = new TFile(filenameOutput.Data(),"recreate");
  TTree *treeOut = new TTree("Samples", "");

  treeOut->Branch("nSmpl",             &nSmpl,               "nSmpl/I");
  treeOut->Branch("nFreq",             &nFreq,               "nFreq/I");
  treeOut->Branch("amplitudeTruth",    &amplitudeTruth,      "amplitudeTruth/D");
  treeOut->Branch("samples",           samples,              "samples[nSmpl]/D");
  treeOut->Branch("samplesNoise",     samplesNoise, "samplesNoise[nSmpl]/D");

  treeOut->Branch("ysamples",     ysamples, "ysamples[1]/D");
  treeOut->Branch("ysamplesNoise",     ysamplesNoise, "ysamplesNoise[1]/D");
  treeOut->Branch("PUOnlywaveform",     PUOnlywaveform_, "PUOnlywaveform[nSmpl]/D");

  int nentries = tree->GetEntries();
  for(int ievt=0; ievt<nentries; ievt++){

    double samplesUncorrelated[NSAMPLES];
    
    for(int i=0; i<NSAMPLES; ++i){
      samplesUncorrelated[i] = rnd.Gaus(0,1);
    }
    
    // Noise correlations
    for(int i=0; i<NSAMPLES; ++i){
      samples[i]=0;
      for(int j=0; j<NSAMPLES; ++j){
	samples[i] += pSh.cholesky(i,j) * samplesUncorrelated[j];
      }
    }

    double maxNoiseSample = -99;

    for(int i=0; i<NSAMPLES; ++i){
      samples[i]   *= sigmaNoise;
      int index = IDSTART + i * NFREQ;
      samplesNoise[i] = samples[i] + PUOnlywaveform[index];
      PUOnlywaveform_[i] = PUOnlywaveform[index];
      //cout<<"isample : Samples : PU : samplesNoise "<<i<<" "<<samples[i]<<" "<<PUOnlywaveform[index]<<" "<<samplesNoise[i]<<endl;

      if( fabs(samplesNoise[i]) > maxNoiseSample)
	maxNoiseSample = fabs(samplesNoise[i]);
    }
    ysamplesNoise[0] = 0;

    ///Normalize wrt max sample
    //double maxNoiseSample = TMath::MaxElement(NSAMPLES, samplesNoise);
    for(int i=0; i<NSAMPLES; ++i){
      //cout<<" Before normalizing, samplesNoise of "<<i<<"th sample is "<<samplesNoise[i]<<endl;
      samplesNoise[i] = samplesNoise[i]/maxNoiseSample;
      //cout<<"samplesNoise of "<<i<<"th sample is "<<samplesNoise[i]<<endl;
      if(samplesNoise[i] > 1) cout<<"WARNING!!! sample noise is > 1"<<endl;
    }

    // add signal and pileup
    //SJ comment: It seems more like adding waveform with PU WITH Noise
    tree->GetEntry(ievt);
    double maxSample = -99;

    for(int i=0; i<NSAMPLES; ++i){
      int index = IDSTART + i * NFREQ;
      samples[i]   += waveform[index];
      if( fabs(samples[i]) > maxSample)
	maxSample = fabs(samples[i]);
    }    
    ysamples[0] = 1;

    /// Normalize wrt max sample
    //double maxSample = TMath::MaxElement(NSAMPLES, samples);
    for(int i=0; i<NSAMPLES; ++i){
      samples[i] = samples[i]/maxSample;
      if(samples[i] > 1) cout<<"WARNING!!! sample EM is > 1"<<endl;
    }


    // true amplitude = in-time pileup + signal
    amplitudeTruth = signalTruth + energyPU[BX0];

    treeOut->Fill();
  }
  
  treeOut->Write();
  fileOut->Close();
  file->Close();

}
