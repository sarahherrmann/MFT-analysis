using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
using o2::itsmft::CompClusterExt; //useful
o2::itsmft::ChipMappingMFT mftChipMapper;
std::vector<MFTTrack> mMFTTracks;
std::vector<MFTCluster> mMFTClusters;

TH1D *histNoccurencesPerMCLabel=0;

void BookHistos();

void FindTrackableTracks(const Char_t *kineFileName = "o2sim_Kine.root", const Char_t *clusterFileName = "mftclusters.root")
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



  int nMFTTrackable;
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
  printf("Number of entries in clusters tree %d \n", nEntries);

  clsTree -> GetEntry(0);//clsTree instead
  int nClusters = (int)clsVec.size(); // Number of mft hits in this event --NEEDS TO BE THE NB OF CLUSTERS
    //std::cout << "Event " << event << " has " << eventHeader->getMCEventStats().getNKeptTracks() << " tracks and " << nMFTHits << " hits\n";
    printf("Number of clusters detected = %d\n", nClusters);

    int nMFTDisksHasClusters = 0;
    int trkID=0, evnID=0, srcID=0;
    bool fake= false;
    long long int mcLabel = 0;
    std::vector<std::array<bool,5>> mcLabelHasClustersInMFTDisks(nClusters,{0,0,0,0,0});//taille ?
    std::vector<long long int> mcLabelTable; //contains mcLabel and the number of occurences of this mcLabel
    std::vector<int> nbLabelOcurrences;
    std::vector<long long int>::iterator it; //iterator to find the index where trkID is in the mcLabelTable
    int index=-1;


      for (int icls = 0; icls < nClusters; ++icls)
      {
        auto clsEntry = icls;//? oui ?
        auto cluster = clsVec[clsEntry];
        auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];//1er label seulement ?? faire une boucle sur les labels ?
        clsLabel.get(trkID, evnID, srcID, fake);
        //printf("icls =%d , label = %llu, trackID = %d\n", icls, clsLabel.getRawValue(), trkID);
        auto clsLayer =  mftChipMapper.chip2Layer(cluster.getChipID());
        int clsMFTdiskID = clsLayer/2; //entier pour root
        it=std::find(mcLabelTable.begin(), mcLabelTable.end(), clsLabel.getRawValue());
          if (it != mcLabelTable.end())//trkID is already in the array mcLabelTable
          {
            index=std::distance(mcLabelTable.begin(), it);//at THIS index
            nbLabelOcurrences[index]+=1;
            printf("#####icls =%d, mcLabel=%llu , nOcurrences =%d\n", icls, clsLabel.getRawValue(), nbLabelOcurrences[index]);
          }
          else //trkID is not yet in mcLabelTable
          {
            index=mcLabelTable.size();
            mcLabelTable.push_back(clsLabel.getRawValue()); //we add it to the table
            nbLabelOcurrences.push_back(1);
          }
        printf("index = %d,diskID =%d\n", index, clsMFTdiskID);

        mcLabelHasClustersInMFTDisks[index][clsMFTdiskID]=true;

      }
      printf("Number of mcLabels stored in the mcLabelTable : %lu\n", mcLabelTable.size());

      for (auto ilabel = 0 ; ilabel < mcLabelTable.size() ; ilabel++)//ilabel is the index corresponding to each mcLabel
      {

        nMFTDisksHasClusters = 0;//has to be put in a separate IsTrackTrackable() method ?
        //mcLabel=mcLabelTable[ilabel];
        for(auto disk: {0,1,2,3,4}) nMFTDisksHasClusters+= int(mcLabelHasClustersInMFTDisks[mcLabel][disk]);
        /*if(nMFTDisksHasClusters>=4)// if IsTrackTrackable(mcLabel)
        {   //Track is trackable if has left at least 1 cluster on at least 4 different disks
          nMFTTrackable++;
          std::cout << "So far we have " << nMFTTrackable << " mftlabels identifying trackable MFT tracks " << std::endl;
          MFTTrackablesEtaZ->Fill(z,eta);
          TrackablepT->Fill(pt);
          Trackablep->Fill(p);
          TrackableEta->Fill(eta);
        }*/


        histNoccurencesPerMCLabel->Fill(nbLabelOcurrences[ilabel]);
      }


      histNoccurencesPerMCLabel->SetXTitle("Number of occurences");
      histNoccurencesPerMCLabel->SetYTitle("Number of labels having this occurence");

      histNoccurencesPerMCLabel->Draw();
}

void BookHistos()
{
  histNoccurencesPerMCLabel = new TH1D("histNOccurencesPerMCLabel", "",100 , 0., 100.);
  //histNoccurencesPerMCLabel->Sumw2();
}
