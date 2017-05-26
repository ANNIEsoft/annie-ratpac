// standard library includes
#include <iostream>
#include <vector>

// ROOT includes
#include "TFile.h"
#include "TTree.h"
#include "TMath.h"
#include "Rtypes.h"
#include "TRandom3.h"

// RAT-PAC includes
#include "RAT/DS/Root.hh"
#include "RAT/DS/MC.hh"
#include "RAT/DS/PMT.hh"
#include "RAT/DSReader.hh"

int rat2reco(char *filename) {

  // create a random number generator
  gRandom = new TRandom3();

  // Load the files
  TFile inputfile(filename,"READ");
  RAT::DSReader *dsReader;
  RAT::DS::Root   *ds;
  dsReader = new RAT::DSReader(inputfile.GetName());
  ULong64_t NbEntries = dsReader->GetTotal();

  // some values about your waveform
  Float_t sigma_noise = 1.*0.001; // in V
  Float_t mean_noise = 0.0*0.001; // in V
  Float_t pulse_shape[6] = {5*0.001/35.,10*0.001/35.,8*0.001/35.,6*0.001/35.,4*0.001/35.,2*0.001/35.};
  Int_t pulse_location[3] = {10000,10010,10013}; // where you want your pulses to be located
  const Int_t nb_samples = 40000; // nb points in your waveform

  Float_t waveform[nb_samples], waveform_samples[nb_samples];

  // Create output file
  TFile outputfile(Form("rat2reco_%s", inputfile.GetName() ),"RECREATE");

  // Tree and branhces
  TTree *PMTData = new TTree("PMTData","PMT Data tree");
  long LastSync ,StartCount, TriggerCount;
  Int_t SequenceID, StartTimeSec, StartTimeNSec, TriggerNumber, CardID, Channel, BufferSize, Trigger, PMTID, PMTx, PMTy, PMTz, PMTf;
  int Rate;
  Float_t Data[nb_samples];
  PMTData->Branch("LastSync",&LastSync,"LastSync/l");// no idea, set to 0
  PMTData->Branch("SequenceID",&SequenceID,"SequenceID/I");// no idea, set to 0
  PMTData->Branch("StartTimeSec",&StartTimeSec,"StartTimeSec/I");// no idea, set to 0
  PMTData->Branch("StartTimeNSec",&StartTimeNSec,"StartTimeNSec/I");// no idea, set to 0
  PMTData->Branch("StartCount",&StartCount,"StartCount/l");// no idea, set to 0
  PMTData->Branch("TriggerNumber",&TriggerNumber,"TriggerNumber/I");// no idea, set to 1
  PMTData->Branch("TriggerCount",&TriggerCount,"TriggerCount/l");// no idea, set to 0
  PMTData->Branch("CardID",&CardID,"CardID/I");
  PMTData->Branch("Channel",&Channel,"Channel/I");
  PMTData->Branch("Rate",&Rate,"Rate/i"); // no idea, set to 0
  PMTData->Branch("BufferSize",&BufferSize,"BufferSize/I");
  PMTData->Branch("Trigger",&Trigger,"Trigger/I"); //event number
  PMTData->Branch("Data",Data,"Data[BufferSize]/F");
  PMTData->Branch("PMTID",&PMTID,"PMTID/I");
  PMTData->Branch("PMTf",&PMTf,"PMTf/I");
  PMTData->Branch("PMTx",&PMTx,"PMTx/I");
  PMTData->Branch("PMTy",&PMTy,"PMTy/I");
  PMTData->Branch("PMTz",&PMTz,"PMTz/I");


  // Variables
  Bool_t broken_tube = false;
  std::vector<Int_t> hitPMT;
  Bool_t ncv_hit;

  for(Int_t i = 0; i < nb_samples; ++i) {
    waveform_samples[i] = i;
  }

  // PMT grid coordinates (not in ratpac x,y,z)
  int pmt_x_array_run1[] = {0,0,0,1,2,3,4,5,6,7,7,7,1,2,2,1,3,4,4,3,5,6,6,5,1,2,2,1,3,4,4,3,5,6,6,5,1,2,2,1,3,4,4,3,5,6,6,5,0,0,0,1,2,3,4,5,6,7,7,7};
  int pmt_z_array_run1[] = {3,2,1,0,0,0,0,0,0,1,2,3,1,1,2,2,1,1,2,2,1,1,2,2,3,3,4,4,3,3,4,4,3,3,4,4,5,5,6,6,5,5,6,6,5,5,6,6,4,5,6,7,7,7,7,7,7,6,5,4};
  int pmt_card_array_run1[] = {3,3,3,3,4,4,4,4,5,5,5,5,6,6,6,6,8,8,8,8,9,9,9,9,10,10,10,10,11,11,11,11,13,13,13,13,14,14,14,14,15,15,15,15,16,16,16,16,18,18,18,18,19,19,19,19,20,20,20,20};
  int pmt_channel_array_run1[] = {0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3,0,1,2,3};

  ULong64_t num_entries = dsReader->GetTotal();
  for(ULong64_t entry = 0; entry < num_entries; ++entry) {
    ds = dsReader->GetEvent(entry);
    if (entry % 10 == 0) std::cout << "Event " << entry << '\n';

    // Some initilizations
    broken_tube = false; hitPMT.clear(); ncv_hit = false;

//     cout << "New evt -- " << endl;

    for( Int_t jPMT = 0; jPMT < ds->GetMC()->GetMCPMTCount(); ++jPMT ){
     	if (ds->GetMC()->GetMCPMT(jPMT)->GetID()+1 == 61 || ds->GetMC()->GetMCPMT(jPMT)->GetID()+1 == 62){
	 ncv_hit = true;
	 break;
	}
    }

    if(ncv_hit) {
    for( Int_t jPMT = 0; jPMT < ds->GetMC()->GetMCPMTCount(); ++jPMT ){

      if (ds->GetMC()->GetMCPMT(jPMT)->GetID()+1 == 6 ||
	ds->GetMC()->GetMCPMT(jPMT)->GetID()+1 == 49 ||
	ds->GetMC()->GetMCPMT(jPMT)->GetID()+1 == 19 ||
	ds->GetMC()->GetMCPMT(jPMT)->GetID()+1 == 37) {
	continue;
	}

// 	 	cout << "Hit: " << ds->GetMC()->GetMCPMT(jPMT)->GetID() << endl;

	BufferSize = nb_samples;
      LastSync = 0;
      StartCount = 0;
      StartTimeSec = 0;
      StartTimeNSec = 0;
      TriggerNumber = 0;
      TriggerCount = 0;
      Rate = 0;
      Trigger = entry;
      SequenceID = entry;

      if (ds->GetMC()->GetMCPMT(jPMT)->GetID()+1 == 61){
// 	 		cout << " NCV" << endl;
	PMTID = 6;
	PMTx = pmt_x_array_run1[5];
	PMTz = pmt_z_array_run1[5];
	CardID = pmt_card_array_run1[5];
	Channel = pmt_channel_array_run1[5];
      } else if (ds->GetMC()->GetMCPMT(jPMT)->GetID()+1 == 62){
// 	 		cout << " NCV" << endl;
	PMTID = 49;
	PMTx = pmt_x_array_run1[48];
	PMTz = pmt_z_array_run1[48];
	CardID = pmt_card_array_run1[48];
	Channel = pmt_channel_array_run1[48];
      } else {
	PMTID = ds->GetMC()->GetMCPMT(jPMT)->GetID()+1;
	PMTx = pmt_x_array_run1[ds->GetMC()->GetMCPMT(jPMT)->GetID()];
	PMTz = pmt_z_array_run1[ds->GetMC()->GetMCPMT(jPMT)->GetID()];
	CardID = pmt_card_array_run1[ds->GetMC()->GetMCPMT(jPMT)->GetID()];
	Channel = pmt_channel_array_run1[ds->GetMC()->GetMCPMT(jPMT)->GetID()];
      }

      hitPMT.push_back(PMTID-1);

      for(Int_t i = 0; i < nb_samples; ++i) {
	Data[i] = gRandom->Gaus(mean_noise,sigma_noise);
      }


      for( Int_t iPhot = 0; iPhot < ds->GetMC()->GetMCPMT(jPMT)->GetMCPhotonCount(); ++iPhot ){ // loop on photons that generated a PE

	if (TMath::FloorNint(ds->GetMC()->GetMCPMT(jPMT)->GetMCPhoton(ds->GetMC()->GetMCPMT(jPMT)->GetMCPhotonCount()-1)->GetHitTime()*0.5 - ds->GetMC()->GetMCPMT(jPMT)->GetMCPhoton(0)->GetHitTime()*0.5) < 25000) {
	  for(Int_t i = 0; i < 6; ++i) {
	    Data[10000 + TMath::FloorNint(ds->GetMC()->GetMCPMT(jPMT)->GetMCPhoton(ds->GetMC()->GetMCPMT(jPMT)->GetMCPhotonCount()-1)->GetHitTime()*0.5 - ds->GetMC()->GetMCPMT(jPMT)->GetMCPhoton(0)->GetHitTime()*0.5)] += pulse_shape[i]*ds->GetMC()->GetMCPMT(jPMT)->GetMCPhoton(iPhot)->GetCharge();
	  }
	}


      }
      PMTData->Fill();
    }

    // PMT loop
    for( Int_t iPMT = 0; iPMT < 64; ++iPMT ){
      if(!ncv_hit) { break;}
      if (find(hitPMT.begin(), hitPMT.end(), iPMT) == hitPMT.end()) {
// 	 	cout << "No hit: " << iPMT << endl;
	for(Int_t i = 0; i < nb_samples; ++i) {
	  Data[i] = gRandom->Gaus(mean_noise,sigma_noise);
	}
	if (iPMT >= 60) {
	  PMTID = 100;
	  PMTx = -10;
	  PMTz = -10;
	  PMTy = -10;
	  if (iPMT == 60) {
	    CardID = 21;
	    Channel = 0;
	  }
	  if (iPMT == 61) {
	    CardID = 21;
	    Channel = 1;
	  }
	  if (iPMT == 62) {
	    CardID = 21;
	    Channel = 2;
	  }
	  if (iPMT == 63) {
	    CardID = 21;
	    Channel = 3;
	  }
	} else {
	  PMTID = iPMT+1;
	  PMTx = pmt_x_array_run1[iPMT];
	  PMTz = pmt_z_array_run1[iPMT];
	  CardID = pmt_card_array_run1[iPMT];
	  Channel = pmt_channel_array_run1[iPMT];
	  PMTy = 0;
	}
	PMTf = 0;

	PMTData->Fill();
      }

    }
    }

  }

  outputfile.Write();
  outputfile.Close();
  inputfile.Close();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Please specify an input file name\n";
    return 1;
  }

  rat2reco(argv[1]);

  return 0;
}
