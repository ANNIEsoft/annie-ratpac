// standard library includes
#include <iostream>
#include <memory>
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

class PMTInfo {
  public:
    PMTInfo(int x, int y, int z, int card, int channel)
      : x_(x), y_(y), z_(z), card_(card), channel_(channel) {}
    int x() const { return x_; }
    int y() const { return y_; }
    int z() const { return z_; }
    int card() const { return card_; }
    int channel() const { return channel_; }

  protected:
    int x_;
    int y_;
    int z_;
    int card_;
    int channel_;
};

const std::map<int, PMTInfo> pmt_id_to_info = {
  // { ID, { x, y, z, card, channel } }
  { 1, { 0, 0, 3, 3, 0 } },
  { 2, { 0, 0, 2, 3, 1 } },
  { 3, { 0, 0, 1, 3, 2 } },
  { 4, { 1, 0, 0, 3, 3 } },
  { 5, { 2, 0, 0, 4, 0 } },
  { 6, { 3, 0, 0, 4, 1 } },
  { 7, { 4, 0, 0, 4, 2 } },
  { 8, { 5, 0, 0, 4, 3 } },
  { 9, { 6, 0, 0, 5, 0 } },
  { 10, { 7, 0, 1, 5, 1 } },
  { 11, { 7, 0, 2, 5, 2 } },
  { 12, { 7, 0, 3, 5, 3 } },
  { 13, { 1, 0, 1, 6, 0 } },
  { 14, { 2, 0, 1, 6, 1 } },
  { 15, { 2, 0, 2, 6, 2 } },
  { 16, { 1, 0, 2, 6, 3 } },
  { 17, { 3, 0, 1, 8, 0 } },
  { 18, { 4, 0, 1, 8, 1 } },
  { 19, { 4, 0, 2, 8, 2 } },
  { 20, { 3, 0, 2, 8, 3 } },
  { 21, { 5, 0, 1, 9, 0 } },
  { 22, { 6, 0, 1, 9, 1 } },
  { 23, { 6, 0, 2, 9, 2 } },
  { 24, { 5, 0, 2, 9, 3 } },
  { 25, { 1, 0, 3, 10, 0 } },
  { 26, { 2, 0, 3, 10, 1 } },
  { 27, { 2, 0, 4, 10, 2 } },
  { 28, { 1, 0, 4, 10, 3 } },
  { 29, { 3, 0, 3, 11, 0 } },
  { 30, { 4, 0, 3, 11, 1 } },
  { 31, { 4, 0, 4, 11, 2 } },
  { 32, { 3, 0, 4, 11, 3 } },
  { 33, { 5, 0, 3, 13, 0 } },
  { 34, { 6, 0, 3, 13, 1 } },
  { 35, { 6, 0, 4, 13, 2 } },
  { 36, { 5, 0, 4, 13, 3 } },
  { 37, { 1, 0, 5, 14, 0 } },
  { 38, { 2, 0, 5, 14, 1 } },
  { 39, { 2, 0, 6, 14, 2 } },
  { 40, { 1, 0, 6, 14, 3 } },
  { 41, { 3, 0, 5, 15, 0 } },
  { 42, { 4, 0, 5, 15, 1 } },
  { 43, { 4, 0, 6, 15, 2 } },
  { 44, { 3, 0, 6, 15, 3 } },
  { 45, { 5, 0, 5, 16, 0 } },
  { 46, { 6, 0, 5, 16, 1 } },
  { 47, { 6, 0, 6, 16, 2 } },
  { 48, { 5, 0, 6, 16, 3 } },
  { 49, { 0, 0, 4, 18, 0 } },
  { 50, { 0, 0, 5, 18, 1 } },
  { 51, { 0, 0, 6, 18, 2 } },
  { 52, { 1, 0, 7, 18, 3 } },
  { 53, { 2, 0, 7, 19, 0 } },
  { 54, { 3, 0, 7, 19, 1 } },
  { 55, { 4, 0, 7, 19, 2 } },
  { 56, { 5, 0, 7, 19, 3 } },
  { 57, { 6, 0, 7, 20, 0 } },
  { 58, { 7, 0, 6, 20, 1 } },
  { 59, { 7, 0, 5, 20, 2 } },
  { 60, { 7, 0, 4, 20, 3 } },
};

constexpr int NCV_PMT_1_ID = 6;
constexpr int NCV_PMT_2_ID = 49;

int rat2reco(const std::string& file_name) {

  // create a random number generator
  gRandom = new TRandom3();

  // Load the files
  TFile input_file(file_name.c_str(), "read");

  std::unique_ptr<RAT::DSReader> dsReader
    = std::make_unique<RAT::DSReader>(input_file.GetName());

  RAT::DS::Root* ds;
  ULong64_t NbEntries = dsReader->GetTotal();

  // some values about your waveform
  Float_t sigma_noise = 1.*0.001; // in V
  Float_t mean_noise = 0.0*0.001; // in V
  Float_t pulse_shape[6] = {5*0.001/35.,10*0.001/35.,8*0.001/35.,6*0.001/35.,4*0.001/35.,2*0.001/35.};
  Int_t pulse_location[3] = {10000,10010,10013}; // where you want your pulses to be located
  const Int_t nb_samples = 40000; // nb points in your waveform
  Int_t trigger_offset = 100; // amount of samples you want to keep free at the beginning of your window

  // Create output file
  TFile output_file(Form("rat2reco_%s", input_file.GetName() ),"RECREATE");

  // Tree and branhces
  TTree* PMTData = new TTree("PMTData","PMT Data tree");
  long LastSync ,StartCount, TriggerCount;
  Int_t SequenceID, StartTimeSec, StartTimeNSec, TriggerNumber, CardID;
  Int_t Channel, BufferSize, Trigger, PMTID, PMTx, PMTy, PMTz, PMTf;
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

  // Branch variables that are the same for every entry
  BufferSize = nb_samples;
  LastSync = 0;
  StartCount = 0;
  StartTimeSec = 0;
  StartTimeNSec = 0;
  TriggerNumber = 0;
  TriggerCount = 0;
  Rate = 0;
  PMTf = 0;

  // Variables
  std::set<int> hit_PMT_indices;
  bool ncv_hit;

  ULong64_t num_entries = dsReader->GetTotal();
  for(ULong64_t entry = 0; entry < num_entries; ++entry) {
    RAT::DS::Root* ds = dsReader->GetEvent(entry);
    auto* mc = ds->GetMC();

    // Use the entry number as the trigger and sequence ID for rat2reco fake
    // data
    Trigger = entry;
    SequenceID = entry;

    if (entry % 10 == 0) std::cout << "Event " << entry << '\n';

    // Some initializations
    hit_PMT_indices.clear();
    ncv_hit = false;

    // cout << "New evt -- " << endl;
    int num_mc_pmts = mc->GetMCPMTCount();
    for ( int jPMT = 0; jPMT < num_mc_pmts; ++jPMT ) {
      int mc_pmt_id = mc->GetMCPMT(jPMT)->GetID();
      hit_PMT_indices.insert(jPMT);
      if (!ncv_hit && (mc_pmt_id == 60 || mc_pmt_id == 61)) ncv_hit = true;
    }

    if (ncv_hit) {
      for ( int iPMT = 0; iPMT < 64; ++iPMT ) {

        RAT::DS::MCPMT* mc_pmt = nullptr;
        if (hit_PMT_indices.count(iPMT)) {
          mc_pmt = mc->GetMCPMT(iPMT);
          int mc_pmt_id = mc_pmt->GetID();

          // Skip the two channels that have been replaced
          // by the NCV PMTs (RAT-PAC IDs 5 and 48). Also
          // skip the two channels (RAT-PAC IDs 18 and 36)
          // that have been replaced by the neutron source
          // trigger PMT and the cosmic trigger input,
          // respectively.
          if (mc_pmt_id == 5 || mc_pmt_id == 48 ||
            mc_pmt_id == 18 || mc_pmt_id == 36)
          {
            continue;
          }

          // NCV PMTs
          if (mc_pmt_id == 60) PMTID = NCV_PMT_1_ID;
          else if (mc_pmt_id == 61) PMTID = NCV_PMT_2_ID;
          // MC IDs are zero-based, data IDs are one-based
          else PMTID = mc_pmt_id + 1;
        }
        else PMTID = iPMT + 1;

        // Get special channel PMT info
        // TODO: clean this up (change map to accomodate these channels?)
        if (iPMT >= 60) {
          PMTID = 100;
          PMTx = -10;
          PMTy = -10;
          PMTz = -10;
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
        }
        else {
          const PMTInfo& info = pmt_id_to_info.at(PMTID);
          PMTx = info.x();
          PMTy = info.y();
          PMTz = info.z();
          CardID = info.card();
          Channel = info.channel();
        }

        // Initialize waveform for this PMT using random noise
        for (int i = 0; i < nb_samples; ++i) {
          Data[i] = gRandom->Gaus(mean_noise, sigma_noise);
        }

        // If there were MC photon hits on this PMT, add those to the waveform
        if (mc_pmt) {
          int num_photons = mc_pmt->GetMCPhotonCount();
          for ( int iPhot = 0; iPhot < num_photons; ++iPhot ) {
            // loop on photons that generated a PE
            const auto* mc_photon = mc_pmt->GetMCPhoton(iPhot);
            double hit_time = mc_photon->GetHitTime(); // ns

            // Factor of 1/2 converts from ns to samples (2 ns resolution)
            int hit_sample = TMath::FloorNint(0.5*hit_time);
            double hit_charge = mc_photon->GetCharge();

            if (hit_sample < nb_samples - trigger_offset) {
              for (int i = 0; i < 6; ++i) {
                Data[trigger_offset + hit_sample] += pulse_shape[i] * hit_charge;
              }
            }
          }
        }

        // Save the data for the current PMT to the tree
        PMTData->Fill();
      }
    }
  }

  output_file.Write();
  output_file.Close();
  input_file.Close();
}

int main(int argc, char* argv[]) {
  if (argc < 2) {
    std::cerr << "Please specify an input file name\n";
    return 1;
  }

  rat2reco(argv[1]);

  return 0;
}
