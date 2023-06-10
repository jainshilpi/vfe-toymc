//
// Takes existing "empty" file. Uses pileup PDFs from "empty"
// file. Creates TTree "Events" and fills each of BX in a bunch train
// with the following info:
//  - number of minbias interactions
//  - total energy from PU
//  - MC truth for signal amplitude
//  - waveform of PU + Signal
//
// Also, picks randomly the in-time BX
//
// suggestions to run:
// - change inputs (see below)
// - run command:
//   > root -l -q Example02.C+
// 
//

#include "Pulse.h"
#include <TFile.h>
#include <TTree.h>
#include <TH1.h>
#include <TF1.h>
#include <TGraph.h>
#include <TRandom.h>
#include <iostream>
#include <sys/stat.h> /// for creating any directory using mkdir

void createDir(string dirPath){
  
  int check = mkdir(dirPath.c_str(), S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
  // check if directory is created or not
  if (!check)
    printf("Directory created\n");
  else {
    printf("Unable to create directory\n");
    exit(1);
  }

}


void Example02()
{

  // make sure these inputs are what you really want
  //const float   signalAmplitude = 10.0;
  //const float   signalAmplitude = 1.0;
  //const float   signalAmplitude = 0.2;
  const float   signalAmplitude = 0.15;
  
  //const TString fileInput       = "data/EmptyFileCRRC43.root";
  const TString fileInput       = "pulseShapeEB.root";
  //const TString fileOutput      = "output.root";
  //string subDir = Form("%d",signalAmplitude);
  const TString fileOutput      = "/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/outputPSWithoutPU";
  //const int     nPU             = 0;
  const int     nPU             = 50;
  //const int     nEventsTotal    = 100;
  const int     nEventsTotal    = 100000;
  //const float   eta             = 0.0;
  const float   eta             = 1.4;
  
  TFile *file = new TFile(fileInput.Data());

  // Get PDF for pileup

  int indx = 10 * fabs(eta) / 0.1;
  if( indx < 0 )  indx = 0;
  if( indx > 13 ) indx = 13;
  char hname[120];
  sprintf(hname,"PileupPDFs/pupdf_%d",indx);
  TH1D *pupdf = (TH1D*)file->Get(hname);
  pupdf->SetDirectory(0);
  
  // Get the Pulse Shape

  Pulse pSh;
   
  // Output file will be created

  TFile *fileOut = new TFile(Form("%s_%f.root",fileOutput.Data(),signalAmplitude),"recreate");
  TTree *tree = new TTree("Waveforms", "");

  // Add branches

  int    BX0;
  int    nBX    = NBXTOTAL;
  int    nMinBias[NBXTOTAL];
  double energyPU[NBXTOTAL];
  int    nWF    = WFLENGTH;
  double waveform[WFLENGTH];
  double signalTruth;
  
  tree->Branch("BX0",         &BX0,         "BX0/I");
  tree->Branch("nBX",         &nBX,         "nBX/I");
  tree->Branch("nWF",         &nWF,         "nWF/I");
  tree->Branch("nMinBias",    nMinBias,     "nMinBias[nBX]/I");
  tree->Branch("energyPU",    energyPU,     "energyPU[nBX]/D");
  tree->Branch("waveform",    waveform,     "waveform[nWF]/D");
  tree->Branch("signalTruth", &signalTruth, "signalTruth/D");
  
  for(int ievt = 0; ievt < nEventsTotal; ievt++){
    for(int ibx = 0; ibx < nBX; ibx++){

      // number of min-bias interactions in each bunch crossing
      nMinBias[ibx] = rnd.Poisson(nPU);

      // total energy per BX
      energyPU[ibx] = 0.;
      //std::cout<<"ibx and nMinBias[ibx] "<<ibx<<" "<<nMinBias[ibx]<<std::endl;
      for(int imb = 0; imb < nMinBias[ibx]; imb++){
	energyPU[ibx] += pow(10., pupdf->GetRandom()); 
	//std::cout<<"ibx is "<< imb <<" Minbias is "<<nMinBias[ibx]<<"energy added to PU "<<pow(10., pupdf->GetRandom())<<" and energy in PU "<<energyPU[ibx]<<std::endl; 
      }

      // pick in-time BX
      
      BX0 = int(nBX * rnd.Rndm());
      //BX0 = 10;
      //std::cout<<"nBX : iev : ibx : BX0 "<<nBX<<" "<<ievt<<" "<<BX0<<std::endl;
    }

    // waveform, initialize
    
    for(int iwf = 0; iwf < nWF; iwf++){
      waveform[iwf] = 0.;
    }

    // add pileup to the waveform
    // time window is nWF ns wide and is centered at BX0

    ///remove PU
    /*
    int ibxMax = min( BX0+11, nBX );
    for(int ibx = 0; ibx < ibxMax; ibx++){
      for(int iwf = 0; iwf < nWF; iwf++){
	double t = (BX0 - ibx) * 25. + iwf - (nWF / 2);
	waveform[iwf] += energyPU[ibx] * pSh.fShape(t);
	//std::cout<<"BX0 : "<<BX0<<" ibx : "<<ibx<<" ,iwf : "<<iwf<<" ,t : "<<t<<std::endl;
      }
    }
    */

    // save MC truth for signal
    
    signalTruth = signalAmplitude;

    // add signal to the waveform

    for(int iwf = 0; iwf < nWF; iwf++){
      waveform[iwf] += signalTruth * pSh.fShape(iwf - (nWF / 2));
    }

    tree->Fill();
  }

  tree->Write();
  fileOut->Close();
  file->Close();
}

