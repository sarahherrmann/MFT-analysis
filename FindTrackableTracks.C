using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
o2::itsmft::ChipMappingMFT mftChipMapper;
std::vector<MFTTrack> mMFTTracks;
std::vector<MFTCluster> mMFTClusters;
std::vector<int> mtrackExtClsIDs;

void FindTrackableTracks()
{

  using o2::itsmft::CompClusterExt;

  constexpr float DefClusErrorRow = 26.88e-4 * 0.5;
  constexpr float DefClusErrorCol = 29.24e-4 * 0.5;
  constexpr float DefClusError2Row = DefClusErrorRow * DefClusErrorRow;
  constexpr float DefClusError2Col = DefClusErrorCol * DefClusErrorCol;

  // Geometry and matrix transformations
  std::string inputGeom = "o2sim_geometry.root";
  o2::base::GeometryManager::loadGeometry(inputGeom);
  auto gman = o2::mft::GeometryTGeo::Instance();
  gman->fillMatrixCache(
    o2::math_utils::bit2Mask(o2::math_utils::TransformType::L2G));

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
  TFile fileC("~/Documents/DOCTORAT/simpp_10Ev_2/mftclusters.root");
  TTree* clsTree = (TTree*)fileC.Get("o2sim");
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
    //clsLabel.get(trkID, evnID, srcID, fake);

  int nEntries = clsTree->GetEntries();
  printf("Number of entries in clusters tree %d \n", nEntries);

  clsTree -> GetEntry(0);//clsTree instead
  int nClusters = (int)clsVec.size(); // Number of mft hits in this event --NEEDS TO BE THE NB OF CLUSTERS
    //std::cout << "Event " << event << " has " << eventHeader->getMCEventStats().getNKeptTracks() << " tracks and " << nMFTHits << " hits\n";
    printf("Number of clusters detected = %d\n", nClusters);
    int trkID=0, evnID=0, srcID=0;
    bool fake= false;
    std::vector<std::array<bool,5>> mcLabelHasClustersInMFTDisks(nClusters,{0,0,0,0,0});//taille ?
    std::vector<int> trkIDTable;
    std::vector<int>::iterator it; //iterator to find the index where trkID is in the trkIDTable
    int index=-1;


      for (int icls = 0; icls < nClusters; ++icls)
      {
        auto clsEntry = icls;//? oui ?
        auto cluster = clsVec[clsEntry];
        auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];//1er label seulement ?? faire une boucle sur les labels ?
        clsLabel.get(trkID, evnID, srcID, fake);
        printf("icls =%d , label = %llu, trackID = %d\n", icls, clsLabel.getRawValue(), trkID);
        auto clsLayer =  mftChipMapper.chip2Layer(cluster.getChipID());
        int clsMFTdiskID = clsLayer/2; //entier pour root
        it=std::find(trkIDTable.begin(), trkIDTable.end(), trkID);
          if (it != trkIDTable.end())//trkID is already in the array trkIDTable
          {
            index=std::distance(trkIDTable.begin(), it);//at THIS index
          }
          else //trkID is not yet in trkIDTable
          {
            index=trkIDTable.size();
            trkIDTable.push_back(trkID); //we add it to the table
          }
        printf("index = %d, trackID = %d, diskID =%d\n", index, trkID, clsMFTdiskID);
        mcLabelHasClustersInMFTDisks[index][clsMFTdiskID]=true;

      }
      printf("Number of trkID stored in the trkIDTable : %lu\n", trkIDTable.size());
      /*for (auto& clabel : clsLabels) // sur tous les MCLabel A REVOIR
      {
        int nMFTDisksHasHits = 0;//has clusters instead of hits now
          for(auto disk: {0,1,2,3,4}) nMFTDisksHasHits+= int(mcLabelHasClustersInMFTDisks[clabel][disk]);
          if(nMFTDisksHasHits>=4)
          {   //Track is trackable if has left hits on at least 4 disks
            nMFTTrackable++;
            std::cout << "So far we have " << nMFTTrackable << " mftlabels identifying trackable MFT tracks " << std::endl;
            MFTTrackablesEtaZ->Fill(z,eta);
            TrackablepT->Fill(pt);
            Trackablep->Fill(p);
            TrackableEta->Fill(eta);
          }
      }*/
      auto& mftTrack : mMFTTracks
}
