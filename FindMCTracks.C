using o2::MCTrack;


TH2D *histPtVsEtaMCTracks=0, *histPhiVsEtaMCTracks=0, *histPhiVsPtMCTracks=0, *histZvtxVsEtaMCTracks=0, *histRVsZMCTracks=0;

void BookHistos();

void FindMCTracks(const Char_t *ofname = "generated_outputfile.root", const Char_t *kineFileName = "o2sim_Kine.root")
{
  BookHistos();

  // MC tracks
  TFile kineFile(kineFileName);//o2sim_Kine.root contains the kinematic of generated MCTracks
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  TFile of_histo_generated(ofname, "RECREATE");

  Int_t nEvents = kineTree->GetEntries();
  double zVtx=0, pt=0, eta=0, phi=0, z=0, R=0; //Of the MCTrack

  for (Int_t event = 0; event < nEvents ; event++)
  {
    kineTree->GetEntry(event);
    for (Int_t trkID = 0 ; trkID < eventHeader->getMCEventStats().getNKeptTracks(); trkID++)
    {
      MCTrack* thisTrack =  &(mcTrkVec)[trkID];
      zVtx = eventHeader->GetZ();//YES ?? check with Antonio

      pt = thisTrack->GetPt();
      eta = TMath::Abs(thisTrack->GetEta());
      phi = thisTrack->GetPhi();
      z = thisTrack->GetStartVertexCoordinatesZ();
      R = sqrt(pow(thisTrack->GetStartVertexCoordinatesX(),2)+pow(thisTrack->GetStartVertexCoordinatesY(),2));

      histPtVsEtaMCTracks->Fill(eta,pt);
      histPhiVsEtaMCTracks->Fill(eta,phi);
      histPhiVsPtMCTracks->Fill(pt,phi);
      histZvtxVsEtaMCTracks->Fill(eta,zVtx);
      histRVsZMCTracks->Fill(z,R);
    }
  }

  TCanvas *cnvPtVsEta = new TCanvas("cnvPtVsEta","cnvPtVsEta",600,600);
  cnvPtVsEta->SetGridx();
  cnvPtVsEta->SetGridy();
  histPtVsEtaMCTracks -> Draw("colz");

  TCanvas *cnvPhiVsEta = new TCanvas("cnvPhiVsEta","cnvPhiVsEta",600,600);
  cnvPhiVsEta->SetGridx();
  cnvPhiVsEta->SetGridy();
  histPhiVsEtaMCTracks -> Draw("colz");

  TCanvas *cnvPhiVsPt = new TCanvas("cnvPhiVsPt","cnvPhiVsPt",600,600);
  cnvPhiVsPt->SetGridx();
  cnvPhiVsPt->SetGridy();
  histPhiVsPtMCTracks -> Draw("colz");



  of_histo_generated.cd();

  histPtVsEtaMCTracks->Write();
  histPhiVsEtaMCTracks->Write();
  histPhiVsPtMCTracks->Write();
  histZvtxVsEtaMCTracks->Write();
  histRVsZMCTracks->Write();

  of_histo_generated.Close();
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



  histPtVsEtaMCTracks = new TH2D("histPtVsEtaMCTracks", "Pt Vs Eta of all MC tracks", 100, 1.0, 4.5, 200, 0., 10.);
  //histPtVsEtaMCTracks->Sumw2();
  histPtVsEtaMCTracks->SetXTitle("#eta of all MC tracks");
  histPtVsEtaMCTracks->SetYTitle("p_{T} of all MC tracks");

  histPhiVsEtaMCTracks = new TH2D("histPhiVsEtaMCTracks", "", 100, 1.0, 4.5, 200, 0., 2*TMath::Pi());
  //histPhiVsEtaMCTracks->Sumw2();
  histPhiVsEtaMCTracks->SetXTitle("#eta of all MC tracks");
  histPhiVsEtaMCTracks->SetYTitle("#phi of all MC tracks");

  histPhiVsPtMCTracks = new TH2D("histPhiVsPtMCTracks", "", 200, 0., 10., 200, 0., 2*TMath::Pi());
  //histPhiVsPtMCTracks->Sumw2();
  histPhiVsPtMCTracks->SetXTitle("p_{T} of all MC tracks");
  histPhiVsPtMCTracks->SetYTitle("#phi of all MC tracks");

  //histZvtxVsEtaMCTracks = new TH2D("histZvtxVsEtaMCTracks", "", sizeof(xEtaLim)/sizeof(xEtaLim[0]), xEtaLim, (sizeof(zVtxBinLimits) / sizeof(zVtxBinLimits[0]) - 1), zVtxBinLimits);
  histZvtxVsEtaMCTracks = new TH2D("histZvtxVsEtaMCTracks", "", 100, 1.0, 4.5, 15, -15, 15);
  histZvtxVsEtaMCTracks->Sumw2();
  histZvtxVsEtaMCTracks->SetXTitle("#eta of all MC tracks");
  histZvtxVsEtaMCTracks->SetYTitle("z_{Vtx} of all MC tracks");

  histRVsZMCTracks = new TH2D("histRVsZMCTracks", "", 400, -80., 20., 100, 0., 20.);
  histRVsZMCTracks->Sumw2();
  histRVsZMCTracks->SetXTitle("z origin of all MC tracks");
  histRVsZMCTracks->SetYTitle("R radius of origin of all MC tracks");

}
