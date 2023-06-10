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


void Example02()
{

  //const float   signalAmplitude = 10.0;
  //const float   signalAmplitude = 1.0;
  //const float   signalAmplitude = 0.2;
  const float   signalAmplitude = 0.15;

  // make sure these inputs are what you really want
  
  //const TString fileInput       = "data/EmptyFileCRRC43.root";
  const TString fileInput       = "pulseShapeEB.root";
  //const TString fileOutput      = "output.root";
  const TString fileOutput      = "/eos/user/s/shilpi/SWAN_projects/ECAL_noise_EM_discrimination/data/outputPSWithPU_PUstudy";
  //const int     nPU             = 0;
  const int     nPU             = 50;
  const int     nEventsTotal    = 100;
  //const int     nEventsTotal    = 100000;
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
  TFile *fileOut = new TFile(Form("%s.root",fileOutput.Data()),"recreate");
  TTree *tree = new TTree("Waveforms", "");

  // Add branches

  int    BX0;
  int    nBX    = NBXTOTAL;
  int    nMinBias[NBXTOTAL];
  double energyPU[NBXTOTAL];
  int    nWF    = WFLENGTH;
  double waveform[WFLENGTH];
  double PUOnlywaveform[WFLENGTH];
  double signalTruth;

  const int NCONT = 10;
  Double_t PUOnlywaveformNCont[NCONT][WFLENGTH];
  
  Double_t PUOnlywaveform0Cont[WFLENGTH], PUOnlywaveform1Cont[WFLENGTH], PUOnlywaveform2Cont[WFLENGTH], PUOnlywaveform3Cont[WFLENGTH], PUOnlywaveform4Cont[WFLENGTH], PUOnlywaveform5Cont[WFLENGTH], PUOnlywaveform6Cont[WFLENGTH], PUOnlywaveform7Cont[WFLENGTH], PUOnlywaveform8Cont[WFLENGTH], PUOnlywaveform9Cont[WFLENGTH]  ;


  int ncontributions = NCONT;

  tree->Branch("BX0",         &BX0,         "BX0/I");
  tree->Branch("nBX",         &nBX,         "nBX/I");
  tree->Branch("nWF",         &nWF,         "nWF/I");
  tree->Branch("nMinBias",    nMinBias,     "nMinBias[nBX]/I");
  tree->Branch("energyPU",    energyPU,     "energyPU[nBX]/D");
  tree->Branch("waveform",    waveform,     "waveform[nWF]/D");
  tree->Branch("signalTruth", &signalTruth, "signalTruth/D");
  tree->Branch("PUOnlywaveform",    PUOnlywaveform,     "PUOnlywaveform[nWF]/D");

  tree->Branch("ncontributions",         &ncontributions,         "ncontributions/I");
  tree->Branch("PUOnlywaveformNCont",    PUOnlywaveformNCont,     "PUOnlywaveformNCont[10][nWF]/D");
  
  tree->Branch("PUOnlywaveform0Cont",    PUOnlywaveform0Cont,     "PUOnlywaveform0Cont[nWF]/D");
  tree->Branch("PUOnlywaveform1Cont",    PUOnlywaveform1Cont,     "PUOnlywaveform1Cont[nWF]/D");
  tree->Branch("PUOnlywaveform2Cont",    PUOnlywaveform2Cont,     "PUOnlywaveform2Cont[nWF]/D");
  tree->Branch("PUOnlywaveform3Cont",    PUOnlywaveform3Cont,     "PUOnlywaveform3Cont[nWF]/D");
  tree->Branch("PUOnlywaveform4Cont",    PUOnlywaveform4Cont,     "PUOnlywaveform4Cont[nWF]/D");
  tree->Branch("PUOnlywaveform5Cont",    PUOnlywaveform5Cont,     "PUOnlywaveform5Cont[nWF]/D");
  tree->Branch("PUOnlywaveform6Cont",    PUOnlywaveform6Cont,     "PUOnlywaveform6Cont[nWF]/D");
  tree->Branch("PUOnlywaveform7Cont",    PUOnlywaveform7Cont,     "PUOnlywaveform7Cont[nWF]/D");
  tree->Branch("PUOnlywaveform8Cont",    PUOnlywaveform8Cont,     "PUOnlywaveform8Cont[nWF]/D");
  tree->Branch("PUOnlywaveform9Cont",    PUOnlywaveform9Cont,     "PUOnlywaveform9Cont[nWF]/D");

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
      PUOnlywaveform0Cont[iwf] = 0.;
      PUOnlywaveform1Cont[iwf] = 0.;
      PUOnlywaveform2Cont[iwf] = 0.;
      PUOnlywaveform3Cont[iwf] = 0.;
      PUOnlywaveform4Cont[iwf] = 0.;
      PUOnlywaveform5Cont[iwf] = 0.;
      PUOnlywaveform6Cont[iwf] = 0.;
      PUOnlywaveform7Cont[iwf] = 0.;
      PUOnlywaveform8Cont[iwf] = 0.;
      PUOnlywaveform9Cont[iwf] = 0.;

    }

    // add pileup to the waveform
    // time window is nWF ns wide and is centered at BX0

    ///remove PU
    int ibxMax = min( BX0+11, nBX );
    for(int ibx = 0; ibx < ibxMax; ibx++){
      for(int iwf = 0; iwf < nWF; iwf++){
	double t = (BX0 - ibx) * 25. + iwf - (nWF / 2);
	waveform[iwf] += energyPU[ibx] * pSh.fShape(t);
	//waveform2D[ibx][iwf] = energyPU[ibx] * pSh.fShape(t); 
	//std::cout<<"BX0 : "<<BX0<<" ibx : "<<ibx<<" ,iwf : "<<iwf<<" ,t : "<<t<<std::endl;
	/// contribution from its BX only +/- N BX to study the shape of a PU waveform
	for(int icont=0; icont<=NCONT; icont++){
	  if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	    //cout<<"ibx : icont : BX0 "<<ibx<<" "<<icont<<" "<<BX0<<endl;
	    PUOnlywaveformNCont[icont][iwf] = (Double32_t)(energyPU[ibx] * pSh.fShape(t));
	    //PUOnlywaveformNCont[icont][iwf] = 0.001;
	    //PUOnlywaveformNCont[icont][iwf] = 1;
	  }
	}//
	

	///above does not work somehow so just store 10 variables for this
	int icont = 0;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform0Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 1;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform1Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 2;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform2Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 3;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform3Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 4;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform4Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 5;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform5Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 6;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform6Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 7;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform7Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 8;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform8Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

	icont = 9;
	if( (ibx>=BX0-icont && BX0-icont>0 && ibx<=BX0+icont) ){
	  PUOnlywaveform9Cont[iwf] += (Double_t)(energyPU[ibx] * pSh.fShape(t));
	}

      }//for(int iwf = 0; iwf < nWF; iwf++)
    }//for(int ibx = 0; ibx < ibxMax; ibx++)
    
    /*
    for(int icont=0; icont<=NCONT; icont++){
      for(int iwf = 0; iwf < nWF; iwf++){
	cout<<"PUOnlywaveformNCont "<<PUOnlywaveformNCont[icont][iwf]<<endl;
      }
    }
    */	

    for(int iwf = 0; iwf < nWF; iwf++){
      PUOnlywaveform[iwf] = waveform[iwf];
    }
	

    // save MC truth for signal
    
    signalTruth = signalAmplitude;

    // add signal to the waveform

    for(int iwf = 0; iwf < nWF; iwf++){

      waveform[iwf] += signalTruth * pSh.fShape(iwf - (nWF / 2));
    }

    tree->Fill();
  }

  fileOut->cd();
  tree->Write();
  fileOut->Close();
  file->Close();
}

