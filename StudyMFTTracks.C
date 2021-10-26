using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
using o2::MCTrack;
using o2::itsmft::CompClusterExt; //useful
o2::itsmft::ChipMappingMFT mftChipMapper;
std::vector<MFTTrack> mMFTTracks;
std::vector<MFTCluster> mMFTClusters;
std::vector<int> mtrackExtClsIDs;

enum typeOfTracks {kGen, kTrackable, kReco, kNumberOfTracks};

std::vector<string> nameOfTracks = {"MCTracks","Trackables","MFTTracks"};


std::vector<TH2D*> histPhiVsEta(kNumberOfTracks);
std::vector<TH2D*> histPtVsEta(kNumberOfTracks); //even though pt not much sense for efficiency
std::vector<TH2D*> histPhiVsPt(kNumberOfTracks);
std::vector<TH2D*> histZvtxVsEta(kNumberOfTracks-1);
std::vector<TH2D*> histRVsZ(kNumberOfTracks-1);

//-----------------------------------------------------------------------------------------------

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

void StudyMFTTracks(const Char_t *ofname = "outputfile_studyTracks.root", const Char_t *kineFileName = "o2sim_Kine.root", const Char_t *clusterFileName = "mftclusters.root")
{
  BookHistos();
  loadMFTTracks();
  TFile of(ofname, "RECREATE");//output file

  // MC tracks file
  TFile kineFile(kineFileName);//o2sim_Kine.root contains the kinematic of generated MCTracks
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);




  double eta[kNumberOfTracks]={0}, phi[kNumberOfTracks]={0}, pt[kNumberOfTracks]={0};
  double zVtx[2]={0}, z[2]={0}, R[2]={0}; //Of the MCTrack


  //-------------------------GENERATED TRACKS-----------------------------------

  Int_t nEvents = kineTree->GetEntries();

  for (Int_t event = 0; event < nEvents ; event++)
  {
    kineTree->GetEntry(event);
    for (Int_t trkID = 0 ; trkID < eventHeader->getMCEventStats().getNKeptTracks(); trkID++)
    {
      MCTrack* thisTrackGen =  &(mcTrkVec)[trkID];
      zVtx[kGen] = eventHeader->GetZ();//YES ?? check with Antonio

      pt[kGen] = thisTrackGen->GetPt();
      eta[kGen] = TMath::Abs(thisTrackGen->GetEta());
      phi[kGen] = thisTrackGen->GetPhi();
      z[kGen] = thisTrackGen->GetStartVertexCoordinatesZ();
      R[kGen] = sqrt(pow(thisTrackGen->GetStartVertexCoordinatesX(),2)+pow(thisTrackGen->GetStartVertexCoordinatesY(),2));

      histPtVsEta[kGen]->Fill(eta[kGen],pt[kGen]);
      histPhiVsEta[kGen]->Fill(eta[kGen],phi[kGen]);
      histPhiVsPt[kGen]->Fill(pt[kGen],phi[kGen]);
      histZvtxVsEta[kGen]->Fill(eta[kGen],zVtx[kGen]);
      histRVsZ[kGen]->Fill(z[kGen],R[kGen]);
      //end of loop on trkIDs
    }
    //end of loop on events
  }


  //-----------------------------TRACKABLES-------------------------------------
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

  // Cluster file
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

    //------------------We look for trackables
    int nEntries = clsTree->GetEntries();
    //printf("Number of entries in clusters tree %d \n", nEntries);

    clsTree -> GetEntry(0);//clsTree instead
    int nClusters = (int)clsVec.size(); // Number of mft hits in this event --NEEDS TO BE THE NB OF CLUSTERS
      //std::cout << "Event " << event << " has " << eventHeader->getMCEventStats().getNKeptTracks() << " tracks and " << nMFTHits << " hits\n";
      //printf("Number of clusters detected = %d\n", nClusters);

      int nMFTDisksHasClusters = 0;
      int nMFTTrackable=0;
      int trackID=0, evnID=0, srcID=0;
      int index=-1;
      bool fake= false;

      o2::MCCompLabel mcLabel;

      std::vector<std::array<bool,5>> mcLabelHasClustersInMFTDisks(nClusters,{0,0,0,0,0});//taille ?
      std::vector<long long int> mcLabelTableRaw; //contains the raw values of the mcLabel which have at least one cluster
      std::vector<o2::MCCompLabel> mcLabelTable; //contains the mcLabel which have at least one cluster
      std::vector<int> nbLabelOcurrences;
      std::vector<long long int>::iterator it; //iterator to find the index where trackID is in the mcLabelTableRaw



        for (int icls = 0; icls < nClusters; ++icls)
        {
          auto clsEntry = icls;//? oui ?
          auto cluster = clsVec[clsEntry];
          auto& clsLabel = (clsLabels->getLabels(clsEntry))[0];//WE SHOULD LOOP OVER LABELS
          clsLabel.get(trackID, evnID, srcID, fake);
          //printf("icls =%d , label = %llu, trackID = %d\n", icls, clsLabel.getRawValue(), trackID);
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
          mcLabel.get(trackID, evnID, srcID, fake);
          kineTree->GetEntry(evnID);

          if (!mcLabel.isNoise())//the Labels corresponding to noise don't give tracks in the MCTrack branch
          {
            MCTrack* thisTrack =  &(mcTrkVec)[trackID];
            zVtx[kTrackable] = eventHeader->GetZ();//YES ?? check with Antonio

            pt[kTrackable] = thisTrack->GetPt();
            eta[kTrackable] = TMath::Abs(thisTrack->GetEta());
            phi[kTrackable] = thisTrack->GetPhi();
            z[kTrackable] = thisTrack->GetStartVertexCoordinatesZ();
            R[kTrackable] = sqrt(pow(thisTrack->GetStartVertexCoordinatesX(),2)+pow(thisTrack->GetStartVertexCoordinatesY(),2));

            //printf("zVtx =%f, pt=%f, eta=%f, phi=%f, R=%f\n", zVtx, pt, eta, phi, R);



            for(auto disk: {0,1,2,3,4}) nMFTDisksHasClusters+= int(mcLabelHasClustersInMFTDisks[ilabel][disk]);
            if(nMFTDisksHasClusters>=4)// if IsTrackTrackable(mcLabel)
            {   //Track is trackable if has left at least 1 cluster on at least 4 different disks
              nMFTTrackable++;
              //std::cout << "So far we have " << nMFTTrackable << " mftlabels identifying trackable MFT tracks " << std::endl;

              histPtVsEta[kTrackable]->Fill(eta[kTrackable],pt[kTrackable]);
              histPhiVsEta[kTrackable]->Fill(eta[kTrackable],phi[kTrackable]);
              histPhiVsPt[kTrackable]->Fill(pt[kTrackable],phi[kTrackable]);
              histZvtxVsEta[kTrackable]->Fill(eta[kTrackable],zVtx[kTrackable]);
              histRVsZ[kTrackable]->Fill(z[kTrackable],R[kTrackable]);
            }


          }
          //____________________________________________________________________SEE WHAT TO DO HERE
          //histNoccurencesPerMCLabel->Fill(nbLabelOcurrences[ilabel]);
        }


        //----------------------------RECONSTRUCTED TRACKS
        for (auto& mftTrack : mMFTTracks)
        {//loop over the MFT tracks
          auto ncls = mftTrack.getNumberOfPoints();
          auto offset = mftTrack.getExternalClusterIndexOffset();
          //std::cout << "\nMFT Track has " << ncls << " clusters:\n";

          //std::cout << " Parameters: Eta = " << mftTrack.getEta() << " Phi =" << mftTrack.getPhi() << std::endl;


          histPtVsEta[kReco]->Fill(TMath::Abs(mftTrack.getEta()),mftTrack.getPt());
          if (mftTrack.getPhi()<0)
          {
            phi[kReco]=mftTrack.getPhi()+2*TMath::Pi();
          }
          else
          {
            phi[kReco]=mftTrack.getPhi();
          }

          histPhiVsEta[kReco]->Fill(TMath::Abs(mftTrack.getEta()),phi[kReco]);
          histPhiVsPt[kReco]->Fill(mftTrack.getPt(),phi[kReco]);
        }
        //----------------------END OF RECO TRACKS

//Write everything in one output file
of.cd();
for (int i = 0; i < kNumberOfTracks ; i++)
  {
    histPhiVsEta[i]->Write();
    histPtVsEta[i]->Write();
    histPhiVsPt[i]->Write();


    if (i < kNumberOfTracks-1)//information only available for generated and trackable tracks
    {
      histZvtxVsEta[i]->Write();
      histRVsZ[i]->Write();
    }
  }
of.Close();



//end of StudyMFTTracks
}

void BookHistos()
{
  for (int i = 0; i < kNumberOfTracks ; i++)
  {
    //histPhiVsEta
    histPhiVsEta[i] = new TH2D((string("histPhiVsEta")+nameOfTracks[i]).c_str(), (string("Phi Vs Eta of ")+nameOfTracks[i]).c_str(), 50, 1.0, 4.5, 100, 0., 2*TMath::Pi());
    histPhiVsEta[i]->SetXTitle((string("#eta of ")+nameOfTracks[i]).c_str());
    histPhiVsEta[i]->SetYTitle((string("#phi of ")+nameOfTracks[i]).c_str());


    //histPtVsEta
    histPtVsEta[i] = new TH2D((string("histPtVsEta")+nameOfTracks[i]).c_str(), (string("Pt Vs Eta of ")+nameOfTracks[i]).c_str(), 50, 1.0, 4.5, 100, 0., 10.);
    histPtVsEta[i]->SetXTitle((string("#eta of ")+nameOfTracks[i]).c_str());
    histPtVsEta[i]->SetYTitle((string("p_{T} of ")+nameOfTracks[i]).c_str());

    //histPhiVsPt
    histPhiVsPt[i] = new TH2D((string("histPhiVsPt")+nameOfTracks[i]).c_str(), (string("Phi Vs Pt of ")+nameOfTracks[i]).c_str(), 100, 0., 10., 100, 0., 2*TMath::Pi());
    histPhiVsPt[i]->SetXTitle((string("p_{T} of ")+nameOfTracks[i]).c_str());
    histPhiVsPt[i]->SetYTitle((string("#phi of ")+nameOfTracks[i]).c_str());

    if (i < kNumberOfTracks-1)//information only available for generated and trackable tracks
    {
      //histZvtxVsEta
      histZvtxVsEta[i] = new TH2D((string("histZvtxVsEta")+nameOfTracks[i]).c_str(), (string("Z_{vtx} Vs Eta of ")+nameOfTracks[i]).c_str(), 50, 1.0, 4.5, 15, -15, 15);
      histZvtxVsEta[i]->SetXTitle((string("#eta of ")+nameOfTracks[i]).c_str());
      histZvtxVsEta[i]->SetYTitle((string("z_{vtx} of ")+nameOfTracks[i]).c_str());

      //histRVsZ
      histRVsZ[i] = new TH2D((string("histRVsZ")+nameOfTracks[i]).c_str(), (string("R Vs Z of ")+nameOfTracks[i]).c_str(), 400, -80., 20., 130, 0., 30.);
      histRVsZ[i]->SetXTitle((string("z origin of ")+nameOfTracks[i]).c_str());
      histRVsZ[i]->SetYTitle((string("R radius of origin of ")+nameOfTracks[i]).c_str());

    }
  }
}
