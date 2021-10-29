using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
using o2::MCTrack;
using o2::itsmft::CompClusterExt; //useful
o2::itsmft::ChipMappingMFT mftChipMapper;
std::vector<MFTTrack> mMFTTracks;
std::vector<MFTCluster> mMFTClusters;

TH1D *histNoccurencesPerMCLabel=0;
TH2D *histPtVsEtaTrackables=0, *histPhiVsEtaTrackables=0, *histPhiVsPtTrackables=0, *histZvtxVsEtaTrackables=0, *histRVsZTrackables=0;


void BookHistos();

void FindTrackableTracks(const Char_t *ofname = "trackable_outputfile.root", const Char_t *kineFileName = "o2sim_Kine.root", const Char_t *clusterFileName = "mftclusters.root")
{
  BookHistos();


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

  // open output file
  TFile of_histo_trackables(ofname, "RECREATE");


  TFile clusterFile(clusterFileName);
  TTree* clsTree = (TTree*)clusterFile.Get("o2sim");
  std::vector<CompClusterExt> clsVec, *clsVecP = &clsVec;
  clsTree->SetBranchAddress("MFTClusterComp", &clsVecP);

  o2::dataformats::MCTruthContainer<o2::MCCompLabel>* clsLabels = nullptr; // This goes to global variables
    if (clsTree->GetBranch("MFTClusterMCTruth"))
    { // This goes to LoadClusters
      clsTree->SetBranchAddress("MFTClusterMCTruth", &clsLabels);
    }
    else
    {
      printf("No Monte-Carlo information in this file\n");
      return;
    }

    // MC tracks
    TFile kineFile(kineFileName);//o2sim_Kine.root contains the kinematic of generated MCTracks
    TTree *kineTree = (TTree*)kineFile.Get("o2sim");
    std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
    kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);
    o2::dataformats::MCEventHeader* eventHeader = nullptr;
    kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);


  int nEntries = clsTree->GetEntries();
  //printf("Number of entries in clusters tree %d \n", nEntries);

  clsTree -> GetEntry(0);//clsTree instead
  int nClusters = (int)clsVec.size(); // Number of mft hits in this event --NEEDS TO BE THE NB OF CLUSTERS
    //std::cout << "Event " << event << " has " << eventHeader->getMCEventStats().getNKeptTracks() << " tracks and " << nMFTHits << " hits\n";
    //printf("Number of clusters detected = %d\n", nClusters);

    int nMFTDisksHasClusters = 0;
    int nMFTTrackable=0;
    int trkID=0, evnID=0, srcID=0;
    int index=-1;

    double zVtx=0, pt=0, eta=0, phi=0, z=0, R=0; //Of the MCTrack

    bool fake= false;

    o2::MCCompLabel mcLabel;

    std::vector<std::array<bool,5>> mcLabelHasClustersInMFTDisks(nClusters,{0,0,0,0,0});//taille ?
    std::vector<long long int> mcLabelTableRaw; //contains the raw values of the mcLabel which have at least one cluster
    std::vector<o2::MCCompLabel> mcLabelTable; //contains the mcLabel which have at least one cluster
    std::vector<int> nbLabelOcurrences;
    std::vector<long long int>::iterator it; //iterator to find the index where trkID is in the mcLabelTableRaw



      for (int icls = 0; icls < nClusters; ++icls)
      {
        auto clsEntry = icls;//? oui ?
        auto cluster = clsVec[clsEntry];
        auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];//WE SHOULD LOOP OVER LABELS
        clsLabel.get(trkID, evnID, srcID, fake);
        //printf("icls =%d , label = %llu, trackID = %d\n", icls, clsLabel.getRawValue(), trkID);
        auto clsLayer =  mftChipMapper.chip2Layer(cluster.getChipID());
        int clsMFTdiskID = clsLayer/2; //entier pour root
        it=std::find(mcLabelTableRaw.begin(), mcLabelTableRaw.end(), clsLabel.getRawValue());
          if (it != mcLabelTableRaw.end())//mcLabel is already in the array mcLabelTableRaw
          {
            index=std::distance(mcLabelTableRaw.begin(), it);//at THIS index
            nbLabelOcurrences[index]+=1;
            //printf("#####icls =%d, mcLabel=%llu , nOcurrences =%d\n", icls, clsLabel.getRawValue(), nbLabelOcurrences[index]);
          }
          else //mcLabel is not yet in mcLabelTableRaw
          {
            index=mcLabelTableRaw.size();
            mcLabelTableRaw.push_back(clsLabel.getRawValue()); //we add it to the tables
            mcLabelTable.push_back(clsLabel);
            nbLabelOcurrences.push_back(1);
          }
        //printf("index = %d,diskID =%d\n", index, clsMFTdiskID);

        mcLabelHasClustersInMFTDisks[index][clsMFTdiskID]=true;

      }
      //printf("Number of mcLabels stored in the mcLabelTableRaw : %lu\n", mcLabelTableRaw.size());

      for (auto ilabel = 0 ; ilabel < mcLabelTableRaw.size() ; ilabel++)//ilabel is the index corresponding to each mcLabel
      {
        //printf("ilabel= %d\n", ilabel);
        nMFTDisksHasClusters = 0;//has to be put in a separate IsTrackTrackable() method ?
        mcLabel=mcLabelTable[ilabel];//this is a real MC label, not only the raw value
        mcLabel.get(trkID, evnID, srcID, fake);
        kineTree->GetEntry(evnID);

        if (!mcLabel.isNoise())//the Labels corresponding to noise don't give tracks in the MCTrack branch
        {
          MCTrack* thisTrack =  &(mcTrkVec)[trkID];
          zVtx = eventHeader->GetZ();//YES ?? check with Antonio

          pt = thisTrack->GetPt();
          eta = TMath::Abs(thisTrack->GetEta());
          phi = thisTrack->GetPhi();
          z = thisTrack->GetStartVertexCoordinatesZ();
          R = sqrt(pow(thisTrack->GetStartVertexCoordinatesX(),2)+pow(thisTrack->GetStartVertexCoordinatesY(),2));

          //printf("zVtx =%f, pt=%f, eta=%f, phi=%f, R=%f\n", zVtx, pt, eta, phi, R);



          for(auto disk: {0,1,2,3,4}) nMFTDisksHasClusters+= int(mcLabelHasClustersInMFTDisks[ilabel][disk]);
          if(nMFTDisksHasClusters>=4)// if IsTrackTrackable(mcLabel)
          {   //Track is trackable if has left at least 1 cluster on at least 4 different disks
            nMFTTrackable++;
            //std::cout << "So far we have " << nMFTTrackable << " mftlabels identifying trackable MFT tracks " << std::endl;

            histPtVsEtaTrackables->Fill(eta,pt);
            histPhiVsEtaTrackables->Fill(eta,phi);
            histPhiVsPtTrackables->Fill(pt,phi);
            histZvtxVsEtaTrackables->Fill(eta,zVtx);
            histRVsZTrackables->Fill(z,R);
          }


        }

        histNoccurencesPerMCLabel->Fill(nbLabelOcurrences[ilabel]);
      }





      of_histo_trackables.cd();

      histNoccurencesPerMCLabel->Write();
      histPtVsEtaTrackables->Write();
      histPhiVsEtaTrackables->Write();
      histPhiVsPtTrackables->Write();
      histZvtxVsEtaTrackables->Write();
      histRVsZTrackables->Write();

      of_histo_trackables.Close();



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


  histNoccurencesPerMCLabel = new TH1D("histNoccurencesPerMCLabel", "",100 , 0., 100.);
  histNoccurencesPerMCLabel->Sumw2();
  histNoccurencesPerMCLabel->SetXTitle("Number of occurences");
  histNoccurencesPerMCLabel->SetYTitle("Number of labels having this occurence");

  histPtVsEtaTrackables = new TH2D("histPtVsEtaTrackables", "Pt Vs Eta of trackable MFT tracks", 100, 1.0, 4.5, 200, 0., 10.);
  //histPtVsEtaTrackables->Sumw2();
  histPtVsEtaTrackables->SetXTitle("#eta of trackable MFT tracks");
  histPtVsEtaTrackables->SetYTitle("#p_{T} of trackable MFT tracks");

  histPhiVsEtaTrackables = new TH2D("histPhiVsEtaTrackables", "", 100, 1.0, 4.5, 200, 0., 2*TMath::Pi());
  //histPhiVsEtaTrackables->Sumw2();
  histPhiVsEtaTrackables->SetXTitle("#eta of trackable MFT tracks");
  histPhiVsEtaTrackables->SetYTitle("#phi of trackable MFT tracks");

  histPhiVsPtTrackables = new TH2D("histPhiVsPtTrackables", "", 200, 0., 10., 200, 0., 2*TMath::Pi());
  //histPhiVsPtTrackables->Sumw2();
  histPhiVsPtTrackables->SetXTitle("#p_{T} of trackable MFT tracks");
  histPhiVsPtTrackables->SetYTitle("#phi of trackable MFT tracks");

  //histZvtxVsEtaTrackables = new TH2D("histZvtxVsEtaTrackables", "", sizeof(xEtaLim)/sizeof(xEtaLim[0]), xEtaLim, (sizeof(zVtxBinLimits) / sizeof(zVtxBinLimits[0]) - 1), zVtxBinLimits);
  histZvtxVsEtaTrackables = new TH2D("histZvtxVsEtaTrackables", "", 100, 1.0, 4.5, 15, -15, 15);
  histZvtxVsEtaTrackables->Sumw2();
  histZvtxVsEtaTrackables->SetXTitle("#eta of trackable MFT tracks");
  histZvtxVsEtaTrackables->SetYTitle("#z_{Vtx} of trackable MFT tracks");

  histRVsZTrackables = new TH2D("histRVsZTrackables", "", 400, -80., 20., 100, 0., 20.);
  histRVsZTrackables->Sumw2();
  histRVsZTrackables->SetXTitle("z origin of trackable MFT tracks");
  histRVsZTrackables->SetYTitle("R radius of origin of trackable MFT tracks");

}
