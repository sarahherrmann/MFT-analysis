using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
using o2::MCTrack;
using o2::itsmft::CompClusterExt;
o2::itsmft::ChipMappingMFT mftChipMapper;
std::vector<MFTTrack> mMFTTracks;
std::vector<MFTCluster> mMFTClusters;
std::vector<int> mtrackExtClsIDs;

enum typeOfTracks {kGen, kTrackable, kReco, kTypeOfTracks};

std::vector<string> nameOfTracks = {"Gen","Trackable","Rec"};//of size kTypeOfTracks


std::vector<TH2D*> histPhiVsEta(kTypeOfTracks);
std::vector<TH2D*> histPtVsEta(kTypeOfTracks); //even though pt not much sense for efficiency
std::vector<TH2D*> histPhiVsPt(kTypeOfTracks);
std::vector<TH2D*> histZvtxVsEta(kTypeOfTracks);
std::vector<TH2D*> histRVsZ(kTypeOfTracks-1);//-1 because only for gen and trackable tracks

bool withMC = false;

std::vector<std::array<bool,5>> mcLabelHasClustersInMFTDisks;
int nCrossedDisksPerLabel = 0;

//-----------------------------------------------------------------------------------------------

void BookHistos();
void loadMFTTracks(const Char_t *recoFileName = "mfttracks.root");
bool IsTrackTrackable(std::array<bool,5> hasClusterInMFTDisks);

//"./mcarchive/tf10/sgn_10_Kine.root"
void StudyMFTTracks(const Char_t *ofname = "outputfile_studyTracks.root", const Char_t *kineFileName = "o2sim_Kine.root", const Char_t *clusterFileName = "mftclusters.root", const Char_t *recoFileName = "mfttracks.root")
{
  BookHistos();
  loadMFTTracks(recoFileName);


  // MC tracks file: initializing generated kinematics tree
  TFile kineFile(kineFileName);//kineFileName contains the kinematic of generated MCTracks
  //if (kineFile.IsOpen())
  TTree *kineTree = (TTree*)kineFile.Get("o2sim");
  std::vector<o2::MCTrack> mcTrkVec, *mcTrkVecP = &mcTrkVec;
  kineTree->SetBranchAddress("MCTrack",&mcTrkVecP);
  o2::dataformats::MCEventHeader* eventHeader = nullptr;
  kineTree->SetBranchAddress("MCEventHeader.", &eventHeader);

  withMC=true;





  double eta[kTypeOfTracks]={0}, phi[kTypeOfTracks]={0}, pt[kTypeOfTracks]={0};
  double zVtx[kTypeOfTracks]={0}, z[kTypeOfTracks-1]={0}, R[kTypeOfTracks-1]={0}; //Of the MCTrack


  //-------------------------GENERATED TRACKS-----------------------------------

  Int_t nEvents = kineTree->GetEntries();

  for (Int_t event = 0; event < nEvents ; event++)
  {
    kineTree->GetEntry(event);
    for (Int_t trkID = 0 ; trkID < eventHeader->getMCEventStats().getNKeptTracks(); trkID++)
    {
      MCTrack* thisTrackGen =  &(mcTrkVec)[trkID];
      zVtx[kGen] = eventHeader->GetZ();//position of the primary vertex on the z axis

      pt[kGen]   = thisTrackGen->GetPt();
      eta[kGen]  = thisTrackGen->GetEta()*-1;
      phi[kGen]  = thisTrackGen->GetPhi();
      z[kGen]    = thisTrackGen->GetStartVertexCoordinatesZ();
      R[kGen]    = sqrt(pow(thisTrackGen->GetStartVertexCoordinatesX(),2)+pow(thisTrackGen->GetStartVertexCoordinatesY(),2));

      if ((R[kGen]<0.2) && (TMath::Abs(z[kGen])<20))//coupure sur l'origine de la trace
      {
        histPtVsEta[kGen]  ->Fill(eta[kGen],pt[kGen]);
        histPhiVsEta[kGen] ->Fill(eta[kGen],phi[kGen]);
        histPhiVsPt[kGen]  ->Fill(pt[kGen],phi[kGen]);
        histZvtxVsEta[kGen]->Fill(eta[kGen],zVtx[kGen]);
        histRVsZ[kGen]     ->Fill(z[kGen],R[kGen]);
      }


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

  // Cluster file: initializing cluster tree
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

  clsTree -> GetEntry(0);//clsTree instead
  int nClusters = (int)clsVec.size(); // Number of cluster

  int nCrossedDisksPerLabel = 0;//number of MFTdisks containing a cluster given from a mclabel
  int trackID=0, evnID=0, srcID=0;
  int index=-1;
  bool fake= false;

  o2::MCCompLabel mcLabel;


  std::array<bool,5> boolClusterInMFTDisks={0,0,0,0,0};
  std::vector<long long int> mcLabelTableRaw; //contains the raw values of the mcLabel which have at least one cluster
  std::vector<o2::MCCompLabel> mcLabelTable; //contains the mcLabel which have at least one cluster
  std::vector<int> nClusterPerLabel;
  std::vector<long long int>::iterator it; //iterator for the mcLabelTableRaw array



        for (int icls = 0; icls < nClusters; ++icls)
        {
          auto cluster = clsVec[icls];
          auto labelSize=(clsLabels->getLabels(icls)).size();
          for (auto il = 0 ; il < labelSize ; il++)
          {
            auto& clsLabel = (clsLabels->getLabels(icls))[il];
            clsLabel.get(trackID, evnID, srcID, fake);


            auto clsLayer =  mftChipMapper.chip2Layer(cluster.getChipID());
            int clsMFTdiskID = clsLayer/2; //entier pour root
            it=std::find(mcLabelTableRaw.begin(), mcLabelTableRaw.end(), clsLabel.getRawValue());
            if (it != mcLabelTableRaw.end())//mcLabel is already in the array mcLabelTableRaw
              {
                index=std::distance(mcLabelTableRaw.begin(), it);//at THIS index
                nClusterPerLabel[index]+=1;
                //printf("#####icls =%d, mcLabel=%llu , nOcurrences =%d\n", icls, clsLabel.getRawValue(), nClusterPerLabel[index]);
              }
            else //mcLabel is not yet in mcLabelTableRaw
              {
                index=mcLabelTableRaw.size();
                mcLabelTableRaw.push_back(clsLabel.getRawValue()); //we add it to the tables
                mcLabelTable.push_back(clsLabel);
                nClusterPerLabel.push_back(1);
                mcLabelHasClustersInMFTDisks.push_back(boolClusterInMFTDisks);

              }
            //printf("index = %d,diskID =%d\n", index, clsMFTdiskID);

            mcLabelHasClustersInMFTDisks[index][clsMFTdiskID]=true;
          //end of the loop on mcLabels for the cluster icls
          }

        //end of the loop on clusters
        }
        //printf("Number of mcLabels stored in the mcLabelTableRaw : %lu\n", mcLabelTableRaw.size());

        for (auto ilabel = 0 ; ilabel < mcLabelTableRaw.size() ; ilabel++)//ilabel is the index corresponding to each mcLabel
        {
          //printf("ilabel= %d\n", ilabel);
          nCrossedDisksPerLabel = 0;//has to be put in a separate IsTrackTrackable() method ?
          mcLabel=mcLabelTable[ilabel];//this is a real MC label, not only the raw value
          mcLabel.get(trackID, evnID, srcID, fake);
          kineTree->GetEntry(evnID);

          if (!mcLabel.isNoise())//the Labels corresponding to noise don't give tracks in the MCTrack branch
          {
            MCTrack* thisTrack =  &(mcTrkVec)[trackID];
            zVtx[kTrackable] = eventHeader->GetZ();

            pt[kTrackable] = thisTrack->GetPt();
            eta[kTrackable] = -1*thisTrack->GetEta();
            phi[kTrackable] = thisTrack->GetPhi();
            z[kTrackable] = thisTrack->GetStartVertexCoordinatesZ();
            R[kTrackable] = sqrt(pow(thisTrack->GetStartVertexCoordinatesX(),2)+pow(thisTrack->GetStartVertexCoordinatesY(),2));

            //printf("zVtx =%f, pt=%f, eta=%f, phi=%f, R=%f\n", zVtx, pt, eta, phi, R);



            if (IsTrackTrackable(mcLabelHasClustersInMFTDisks[ilabel]))
            {   //Track is trackable if has left at least 1 cluster on at least 4 different disks

              histPtVsEta[kTrackable]  ->Fill(eta[kTrackable],pt[kTrackable]);
              histPhiVsEta[kTrackable] ->Fill(eta[kTrackable],phi[kTrackable]);
              histPhiVsPt[kTrackable]  ->Fill(pt[kTrackable],phi[kTrackable]);
              histZvtxVsEta[kTrackable]->Fill(eta[kTrackable],zVtx[kTrackable]);
              histRVsZ[kTrackable]     ->Fill(z[kTrackable],R[kTrackable]);
            }


          }
          //____________________________________________________________________SEE WHAT TO DO HERE
          //histNoccurencesPerMCLabel->Fill(nClusterPerLabel[ilabel]);
        }


        //----------------------------RECONSTRUCTED TRACKS
        for (auto& mftTrack : mMFTTracks)
        {//loop over the MFT tracks
          auto ncls = mftTrack.getNumberOfPoints();
          eta[kReco] = -1*mftTrack.getEta();
          pt[kReco] = mftTrack.getPt();


          if (mftTrack.getPhi()>=TMath::Pi()/2)
          {
            phi[kReco]=-mftTrack.getPhi()+TMath::Pi()/2+2*TMath::Pi();
          }
          else
          {
            phi[kReco]=-mftTrack.getPhi()+TMath::Pi()/2;
          }

          histPtVsEta[kReco]->Fill(eta[kReco],pt[kReco]);
          histPhiVsEta[kReco]->Fill(eta[kReco],phi[kReco]);
          histPhiVsPt[kReco]->Fill(pt[kReco],phi[kReco]);
        }
        //----------------------END OF RECO TRACKS

TFile of(ofname, "RECREATE");//output file
//Write everything in one output file
of.cd();
for (int i = 0; i < kTypeOfTracks ; i++)
  {
    histPhiVsEta[i]->Write();
    histPtVsEta[i]->Write();
    histPhiVsPt[i]->Write();


    if (i < kTypeOfTracks-1)//information only available for generated and trackable tracks
    {
      histZvtxVsEta[i]->Write();
      histRVsZ[i]->Write();
    }
  }
of.Close();



//end of StudyMFTTracks
}

void loadMFTTracks(const Char_t *recoFileName = "mfttracks.root")
{
  // Load all mft tracks

  TFile* trkFileIn = new TFile(recoFileName);
  TTree* mftTrackTree = (TTree*)trkFileIn->Get("o2sim");
  std::vector<o2::mft::TrackMFT> trackMFTVec, *trackMFTVecP = &trackMFTVec;
  mftTrackTree->SetBranchAddress("MFTTrack", &trackMFTVecP);

  std::vector<int> trackExtClsVec, *trackExtClsVecP = &trackExtClsVec;
  mftTrackTree->SetBranchAddress("MFTTrackClusIdx", &trackExtClsVecP);

  mftTrackTree->GetEntry(0);
  mMFTTracks.swap(trackMFTVec);
  mtrackExtClsIDs.swap(trackExtClsVec);
}

bool IsTrackTrackable(std::array<bool,5> hasClusterInMFTDisks)
{
  nCrossedDisksPerLabel = 0;
  for(auto disk: {0,1,2,3,4}) nCrossedDisksPerLabel+= int(hasClusterInMFTDisks[disk]);
  if(nCrossedDisksPerLabel>=4)
  {
    return true;
  }
  else
  {
    return false;
  }
}



void BookHistos()
{
  for (int i = 0; i < kTypeOfTracks ; i++)
  {
    //histPhiVsEta
    histPhiVsEta[i] = new TH2D((string("histPhiVsEta")+nameOfTracks[i]).c_str(), (string("Phi Vs Eta of ")+nameOfTracks[i]).c_str(), 35, 1.0, 4.5, 24, 0., 2*TMath::Pi());
    histPhiVsEta[i]->SetXTitle((string("#eta of ")+nameOfTracks[i]).c_str());
    histPhiVsEta[i]->SetYTitle((string("#phi of ")+nameOfTracks[i]).c_str());
    histPhiVsEta[i]->Sumw2();

    //histPtVsEta
    histPtVsEta[i] = new TH2D((string("histPtVsEta")+nameOfTracks[i]).c_str(), (string("Pt Vs Eta of ")+nameOfTracks[i]).c_str(), 35, 1.0, 4.5, 40, 0., 10.);
    histPtVsEta[i]->SetXTitle((string("#eta of ")+nameOfTracks[i]).c_str());
    histPtVsEta[i]->SetYTitle((string("p_{T} (GeV/c) of ")+nameOfTracks[i]).c_str());
    histPtVsEta[i]->Sumw2();

    //histPhiVsPt
    histPhiVsPt[i] = new TH2D((string("histPhiVsPt")+nameOfTracks[i]).c_str(), (string("Phi Vs Pt of ")+nameOfTracks[i]).c_str(), 40, 0., 10., 24, 0., 2*TMath::Pi());
    histPhiVsPt[i]->SetXTitle((string("p_{T} (GeV/c) of ")+nameOfTracks[i]).c_str());
    histPhiVsPt[i]->SetYTitle((string("#phi of ")+nameOfTracks[i]).c_str());
    histPhiVsPt[i]->Sumw2();

    //histZvtxVsEta
    histZvtxVsEta[i] = new TH2D((string("histZvtxVsEta")+nameOfTracks[i]).c_str(), (string("Z_{vtx} Vs Eta of ")+nameOfTracks[i]).c_str(), 35, 1.0, 4.5, 15, -15, 15);
    histZvtxVsEta[i]->SetXTitle((string("#eta of ")+nameOfTracks[i]).c_str());
    histZvtxVsEta[i]->SetYTitle((string("z_{vtx} (cm) of ")+nameOfTracks[i]).c_str());
    histZvtxVsEta[i]->Sumw2();

    if (i < kTypeOfTracks-1)//information only available for generated and trackable tracks
    {
      //histRVsZ
      histRVsZ[i] = new TH2D((string("histRVsZ")+nameOfTracks[i]).c_str(), (string("R Vs Z of ")+nameOfTracks[i]).c_str(), 400, -80., 20., 400, 0., 80.);
      histRVsZ[i]->SetXTitle((string("z (cm) origin of ")+nameOfTracks[i]).c_str());
      histRVsZ[i]->SetYTitle((string("R (cm) radius of origin of ")+nameOfTracks[i]).c_str());
      histRVsZ[i]->Sumw2();

    }
  }
}
