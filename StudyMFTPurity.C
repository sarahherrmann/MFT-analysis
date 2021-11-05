using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
using o2::MCTrack;
using o2::itsmft::CompClusterExt;

std::vector<MFTTrack> mMFTTracks;
std::vector<MFTCluster> mMFTClusters;
std::vector<int> mtrackExtClsIDs;

TH2D *histPhiRecVsPhiGen=0, *histEtaRecVsEtaGen=0, *histPhiVsEtaPurityTrue=0, *histPhiVsEtaPurityRec=0;

enum typeOfTracks {kGen, kRecoTrue, kReco, kTypeOfTracks};
std::vector<string> nameOfTracks = {"Gen","RecTrue","Rec"};//of size kTypeOfTracks

void BookHistos();

//StudyMFTPurity("output.root","./mcarchive/tf10/sgn_10_Kine.root","./mcarchive/tf10/mftclusters.root", "./mcarchive/tf10/mfttracks.root")
std::pair<long long int, int> findEntryWithLargestValue(std::map<long long int, int> sampleMap);
//______________________________________________________________________________

void StudyMFTPurity(const Char_t *ofname = "outputfile_studyTracks.root", const Char_t *kineFileName = "o2sim_Kine.root", const Char_t *clusterFileName = "mftclusters.root", const Char_t *trkFileName = "mfttracks.root")
{

  BookHistos();
  TFile of(ofname, "RECREATE");//output file

  double eta[kTypeOfTracks]={0}, phi[kTypeOfTracks]={0}, pt[kTypeOfTracks]={0};
  double zVtx[kTypeOfTracks]={0}, z[kTypeOfTracks-1]={0}, R[kTypeOfTracks-1]={0}; //Of the MCTrack

  // Load the file containing MFT tracks
  // Load the tree branches "MFTTrack","MFTTrackMCTruth" and "MFTTrackClusIdx"
  TFile trkFile(trkFileName);
  TTree *trkTree = (TTree*)trkFile.Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackVec, *trackVecP = &trackVec;
  trkTree->SetBranchAddress("MFTTrack", &trackVecP);
  std::vector<o2::MCCompLabel>* trkLabels = nullptr;
  if (trkTree->GetBranch("MFTTrackMCTruth"))
  {
    trkTree->SetBranchAddress("MFTTrackMCTruth", &trkLabels);
  }
  else
  {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  trkTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);

  trkTree->GetEntry(0);//the trkTree has only one Entry


  //______________________________________________________________________________
  // Cluster pattern dictionary
  std::string dictfile = "MFTdictionary.bin";
  o2::itsmft::TopologyDictionary dict;
  std::ifstream file(dictfile.c_str());
  if (file.good()) {
    printf("Running with dictionary: %s \n", dictfile.c_str());
    dict.readBinaryFile(dictfile);
  } else {
    printf("Can not run without dictionary !\n");
    return;
  }

  //Load the file containing clusters
  //Initialize cluster tree
  TFile clusterFile(clusterFileName);
  TTree* clsTree = (TTree*)clusterFile.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);

  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clsLabels = nullptr; // This goes to global variables
  if (clsTree->GetBranch("MFTClusterMCTruth"))
  {
    clsTree->SetBranchAddress("MFTClusterMCTruth", &clsLabels);
  }
  else
  {
    printf("No Monte-Carlo information in this file\n");
    return;
  }
  clsTree -> GetEntry(0);

  // MC tracks file: initializing generated kinematics tree
  TFile kineFile(kineFileName);//kineFileName contains the kinematic of generated MCTracks
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  int srcID, trkID, evnID;
  bool fake;
  Int_t iTrack = 0;
  for (auto &track : trackVec)
  {
    auto ncls = track.getNumberOfPoints();
    auto offset = track.getExternalClusterIndexOffset();
    std::map<long long int, int> mcLabels;
    for (int icls = 0; icls < ncls; ++icls)
    {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto cluster = clsVec[clsEntry];

      for (int ilabel = 0; ilabel < (clsLabels->getLabels(clsEntry)).size(); ++ilabel)
      {
        auto& clsLabel = (clsLabels->getLabels(clsEntry))[ilabel];
        clsLabel.get(trkID, evnID, srcID, fake);
        if (!clsLabel.isNoise())
        {
  	       mcLabels[clsLabel.getRawValue()] += 1;
        }
      }


    }//end of cluster loop 1

    std::pair<long long int, int> entryWithMaxValue = findEntryWithLargestValue(mcLabels);



    Int_t thisEvnID = -1, thisSrcID = -1, thisTrkID = -1, thisEventIDLabel = -1, nPoints = track.getNumberOfPoints();
    bool thisFake=false;

    //cout << "Entry with highest value: "<< entryWithMaxValue.first << " = "<< entryWithMaxValue.second << endl;
    //printf("Number of points for this track = %d\n", nPoints);
    for (int icls = 0; icls < ncls; ++icls)
    {
      auto clsEntry = trackExtClsVec[offset + icls];
      auto cluster = clsVec[clsEntry];
      for (int ilabel = 0; ilabel < (clsLabels->getLabels(clsEntry)).size(); ++ilabel)
      {
        auto& clsLabel = (clsLabels->getLabels(clsEntry))[ilabel];
        clsLabel.get(trkID, evnID, srcID, fake);



        if ((!clsLabel.isNoise()) && (clsLabel.getRawValue()==entryWithMaxValue.first))
        {
  	       if (((Float_t)(mcLabels[clsLabel.getRawValue()]) / (Float_t)(nPoints)) >= 0.8)
           { // Must have at least 80% of its clusters from the same MC Track
  	          thisTrkID = trkID;
  	          thisSrcID = srcID;
  	          thisEvnID = evnID;
  	          thisEventIDLabel = icls;
              thisFake = fake;

  	       }
        }
      }

    }//end of cluster loop 2

    //after that, if thisTrkID, thisSrcID, thisEvnID are still equal to -1, it means the track isn't considered a true track !

    if (thisTrkID==-1)
    {
      printf("NO This track number %d in the reco track tree is not a true track\n", iTrack);
    }
    else
    {
      printf("YES this track number %d in the reco track tree of trkID= %d is true\n", iTrack, thisTrkID);
      kineTree->GetEntry(thisEvnID);

      MCTrack* thisTrack =  &(mcTrkVec)[thisTrkID];
      //zVtx[kGen] = eventHeader->GetZ();

      //pt[kGen] = thisTrack->GetPt();
      eta[kGen] = -1*thisTrack->GetEta();
      phi[kGen] = thisTrack->GetPhi();

      eta[kRecoTrue]=-1*track.getEta();

      if (track.getPhi()>=TMath::Pi()/2)
      {
        phi[kRecoTrue]=-track.getPhi()+TMath::Pi()/2+2*TMath::Pi();
      }
      else
      {
        phi[kRecoTrue]=-track.getPhi()+TMath::Pi()/2;
      }

      printf("phigen=%f, phirec=%f\n", phi[kGen], phi[kRecoTrue]);
      histPhiRecVsPhiGen->Fill(phi[kGen], phi[kRecoTrue]);
      histEtaRecVsEtaGen->Fill(eta[kGen], eta[kRecoTrue]);
      histPhiVsEtaPurityTrue->Fill(eta[kRecoTrue], phi[kRecoTrue]);

    }

    eta[kReco]=-1*track.getEta();

    if (track.getPhi()>=TMath::Pi()/2)
    {
      phi[kReco]=-track.getPhi()+TMath::Pi()/2+2*TMath::Pi();
    }
    else
    {
      phi[kReco]=-track.getPhi()+TMath::Pi()/2;
    }

    histPhiVsEtaPurityRec->Fill(eta[kReco], phi[kReco]);
    iTrack++;
  }

  auto cnv_phi = new TCanvas("cnv_phi","",800,800);
  histPhiRecVsPhiGen->Draw("colz");

  auto cnv_eta = new TCanvas("cnv_eta","",800,800);
  histEtaRecVsEtaGen->Draw("colz");

  auto cnv_RecoTrue = new TCanvas("cnv_RecoTrue","",800,800);
  histPhiVsEtaPurityTrue->Draw("colz");

  auto cnv_AllReco = new TCanvas("cnv_AllReco","",800,800);
  histPhiVsEtaPurityRec->Draw("colz");

  of.cd();
  histPhiRecVsPhiGen     ->Write();
  histEtaRecVsEtaGen     ->Write();
  histPhiVsEtaPurityTrue ->Write();
  histPhiVsEtaPurityRec  ->Write();
  of.Close();



}//end of StudyMFTPurity()

// Function to find the Entry
// with largest Value in a Map
std::pair<long long int, int> findEntryWithLargestValue(std::map<long long int, int> sampleMap)
{

	// Reference variable to help find
	// the entry with the highest value
	std::pair<long long int, int> entryWithMaxValue = std::make_pair(0, 0);

	// Iterate in the map to find the required entry
	std::map<long long int, int>::iterator currentEntry;
	for (currentEntry = sampleMap.begin(); currentEntry != sampleMap.end(); ++currentEntry)
  {

		// If this entry's value is more
		// than the max value
		// Set this entry as the max
		if (currentEntry->second > entryWithMaxValue.second)
    {

			entryWithMaxValue= std::make_pair(currentEntry->first,currentEntry->second);
		}
	}

	return entryWithMaxValue;
}

void BookHistos()
{

  //histPhiRecVsPhiGen = new TH2D("histPhiRecVsPhiGen", "Phi Rec Vs Phi Gen of true reco tracks ", 12, -TMath::Pi(), TMath::Pi(), 12, -TMath::Pi(), TMath::Pi());
  histPhiRecVsPhiGen = new TH2D("histPhiRecVsPhiGen", "Phi Rec Vs Phi Gen of true reco tracks ", 24, 0., 2*TMath::Pi(), 24, 0., 2*TMath::Pi());

  histPhiRecVsPhiGen->SetXTitle((string("#phi of ")+nameOfTracks[kGen]).c_str());
  histPhiRecVsPhiGen->SetYTitle((string("#phi of ")+nameOfTracks[kRecoTrue]).c_str());
  histPhiRecVsPhiGen->Sumw2();

  histEtaRecVsEtaGen = new TH2D("histEtaRecVsEtaGen", "Eta Rec Vs Eta Gen of true reco tracks ", 35, 1.0, 4.5, 35, 1.0, 4.5);
  histEtaRecVsEtaGen->SetXTitle((string("#eta of ")+nameOfTracks[kGen]).c_str());
  histEtaRecVsEtaGen->SetYTitle((string("#eta of ")+nameOfTracks[kRecoTrue]).c_str());
  histEtaRecVsEtaGen->Sumw2();

  histPhiVsEtaPurityTrue = new TH2D("histPhiVsEtaPurityTrue", "Phi Vs Eta of true reco tracks ", 35, 1.0, 4.5, 24, 0., 2*TMath::Pi());
  histPhiVsEtaPurityTrue->SetXTitle((string("#eta of ")+nameOfTracks[kRecoTrue]).c_str());
  histPhiVsEtaPurityTrue->SetYTitle((string("#phi of ")+nameOfTracks[kRecoTrue]).c_str());
  histPhiVsEtaPurityTrue->Sumw2();

  histPhiVsEtaPurityRec = new TH2D("histPhiVsEtaPurityRec", "Phi Vs Eta of all reco tracks ", 35, 1.0, 4.5, 24, 0., 2*TMath::Pi());
  histPhiVsEtaPurityRec->SetXTitle((string("#eta of ")+nameOfTracks[kReco]).c_str());
  histPhiVsEtaPurityRec->SetYTitle((string("#phi of ")+nameOfTracks[kReco]).c_str());
  histPhiVsEtaPurityRec->Sumw2();



}
