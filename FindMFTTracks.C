using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
std::vector<MFTTrack> mMFTTracks;
std::vector<MFTCluster> mMFTClusters;
std::vector<int> mtrackExtClsIDs;

//_________________________________________________________________________________________________


TH2D *histPtVsEtaMFTTracks=0, *histPhiVsEtaMFTTracks=0, *histPhiVsPtMFTTracks=0, *histZvtxVsEtaMFTTracks=0, *histRVsZMFTTracks=0;


void BookHistos();


void loadMFTTracks()
{
  // Load all mft tracks

  std::string trkFile = "mfttracks.root";
  TFile* trkFileIn = new TFile(trkFile.c_str());
  TTree* mftTrackTree = (TTree*)trkFileIn->Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackMFTVec, *trackMFTVecP = &trackMFTVec;
  mftTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);

  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  mftTrackTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);

  mftTrackTree->GetEntry(0);
  mMFTTracks.swap(trackMFTVec);
  mtrackExtClsIDs.swap(trackExtClsVec);
}

void FindMFTTracks(const Char_t *ofname = "reconstructed_outputfile.root")
{
  loadMFTTracks();   // Loads mfttracks into std::vector<o2::mft::TrackMFT> mMFTTracks;
  BookHistos();

  double phi=0;

  TFile of_histo_reconstructed(ofname, "RECREATE");

  // Print track parameters and cluster positions for all mft tracks
  for (auto& mftTrack : mMFTTracks)
  {//loop over the MFT tracks
    auto ncls = mftTrack.getNumberOfPoints();
    auto offset = mftTrack.getExternalClusterIndexOffset();
    //std::cout << "\nMFT Track has " << ncls << " clusters:\n";

    //std::cout << " Parameters: Eta = " << mftTrack.getEta() << " Phi =" << mftTrack.getPhi() << std::endl;


    histPtVsEtaMFTTracks->Fill(TMath::Abs(mftTrack.getEta()),mftTrack.getPt());
    if (mftTrack.getPhi()<0)
    {
      phi=mftTrack.getPhi()+2*TMath::pi();
    }
    else
    {
      phi=mftTrack.getPhi();
    }

    histPhiVsEtaMFTTracks->Fill(TMath::Abs(mftTrack.getEta()),phi);
    histPhiVsPtMFTTracks->Fill(mftTrack.getPt(),phi);
  }

  TCanvas *cnvPtVsEta = new TCanvas("cnvPtVsEta","cnvPtVsEta",600,600);
  cnvPtVsEta->SetGridx();
  cnvPtVsEta->SetGridy();
  histPtVsEtaMFTTracks -> Draw("colz");

  TCanvas *cnvPhiVsEta = new TCanvas("cnvPhiVsEta","cnvPhiVsEta",600,600);
  cnvPhiVsEta->SetGridx();
  cnvPhiVsEta->SetGridy();
  histPhiVsEtaMFTTracks -> Draw("colz");

  TCanvas *cnvPhiVsPt = new TCanvas("cnvPhiVsPt","cnvPhiVsPt",600,600);
  cnvPhiVsPt->SetGridx();
  cnvPhiVsPt->SetGridy();
  histPhiVsPtMFTTracks -> Draw("colz");

  of_histo_reconstructed.cd();

  histPtVsEtaMFTTracks->Write();
  histPhiVsEtaMFTTracks->Write();
  histPhiVsPtMFTTracks->Write();
  //histZvtxVsEtaMFTTracks->Write();
  //histRVsZMFTTracks->Write();

  of_histo_reconstructed.Close();
}

void BookHistos()
{
  float zVtxBinLimits[] = {-15., -10., -5., -3., -1., 1., 3., 5., 10., 15.};
  float xEtaLim[100]={};//100 is the number of bins for eta

  //initialisation of the array of limits for eta, xEtaLim
  float xEta=1.5; //eta min here

  for (int i = 0 ; i < sizeof(xEtaLim)/sizeof(xEtaLim[0]) ; i++)
  {
    xEtaLim[i]=xEta;
    xEta+=0.03;
  }



  histPtVsEtaMFTTracks = new TH2D("histPtVsEtaMFTTracks", "Pt Vs Eta of trackable MFT tracks", 50, 1.0, 4.5, 100, 0., 10.);
  //histPtVsEtaMFTTracks->Sumw2();
  histPtVsEtaMFTTracks->SetXTitle("#eta of trackable MFT tracks");
  histPtVsEtaMFTTracks->SetYTitle("#p_{T} of trackable MFT tracks");

  histPhiVsEtaMFTTracks = new TH2D("histPhiVsEtaMFTTracks", "", 50, 1.0, 4.5, 100, 0., 2*TMath::Pi());
  //histPhiVsEtaMFTTracks->Sumw2();
  histPhiVsEtaMFTTracks->SetXTitle("#eta of trackable MFT tracks");
  histPhiVsEtaMFTTracks->SetYTitle("#phi of trackable MFT tracks");

  histPhiVsPtMFTTracks = new TH2D("histPhiVsPtMFTTracks", "", 100, 0., 10., 100, 0., 2*TMath::Pi());
  //histPhiVsPtMFTTracks->Sumw2();
  histPhiVsPtMFTTracks->SetXTitle("#p_{T} of trackable MFT tracks");
  histPhiVsPtMFTTracks->SetYTitle("#phi of trackable MFT tracks");

  //histZvtxVsEtaMFTTracks = new TH2D("histZvtxVsEtaMFTTracks", "", sizeof(xEtaLim)/sizeof(xEtaLim[0]), xEtaLim, (sizeof(zVtxBinLimits) / sizeof(zVtxBinLimits[0]) - 1), zVtxBinLimits);
  histZvtxVsEtaMFTTracks = new TH2D("histZvtxVsEtaMFTTracks", "", 50, 1.0, 4.5, 15, -15, 15);
  histZvtxVsEtaMFTTracks->Sumw2();
  histZvtxVsEtaMFTTracks->SetXTitle("#eta of trackable MFT tracks");
  histZvtxVsEtaMFTTracks->SetYTitle("#z_{Vtx} of trackable MFT tracks");

  histRVsZMFTTracks = new TH2D("histRVsZMFTTracks", "", 400, -80., 20., 100, 0., 20.);
  histRVsZMFTTracks->Sumw2();
  histRVsZMFTTracks->SetXTitle("z origin of trackable MFT tracks");
  histRVsZMFTTracks->SetYTitle("R radius of origin of trackable MFT tracks");

}
