using MFTCluster = o2::BaseCluster<float>;
using MFTTrack = o2::mft::TrackMFT;
using o2::MCTrack;
using o2::itsmft::CompClusterExt;
o2::itsmft::ChipMappingMFT mftChipMapper;
std::vector<MFTTrack> mMFTTracks;
std::vector<MFTCluster> mMFTClusters;
std::vector<int> mtrackExtClsIDs;

enum typeOfTracks {kGen, kTrackable, kReco, kRecoTrue, kTypeOfTracks};

std::vector<string> nameOfTracks = {"Gen","Trackable","Rec","Rec and True"};//of size kTypeOfTracks

TH2D *histPhiRecVsPhiGen=0, *histEtaRecVsEtaGen=0;

std::vector<TH2D*> histPhiVsEta(kTypeOfTracks);
std::vector<TH2D*> histPtVsEta(kTypeOfTracks); //even though pt not much sense for efficiency
std::vector<TH2D*> histPhiVsPt(kTypeOfTracks);
std::vector<TH2D*> histZvtxVsEta(kTypeOfTracks);
std::vector<TH2D*> histRVsZ(kTypeOfTracks-2);//-2 because only for gen and trackable tracks



std::vector<std::array<bool,5>> mcLabelHasClustersInMFTDisks;
int nCrossedDisksPerLabel = 0;



//-----------------------------------------------------------------------------------------------

TH1D *histCompteurDeMCLabel=0, *histPtOfPrimary=0, *histPOfPrimary=0, *histEtaOfVeryLowPPrimary=0;
TH1D *histPxOfVeryLowPPrimary=0, *histPyOfVeryLowPPrimary=0, *histPzOfVeryLowPPrimary=0;
TH1I *histPdgOfVeryLowPPrimary=0;
//---------------------------------------------------------------------------------------

void BookHistos();
void loadMFTTracks(const Char_t *recoFileName = "mfttracks.root");
bool IsTrackTrackable(std::array<bool,5> hasClusterInMFTDisks);
std::pair<uint64_t, int> findEntryWithLargestValue(std::map<uint64_t, int> sampleMap);



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





  double eta[kTypeOfTracks]={0}, phi[kTypeOfTracks]={0}, pt[kTypeOfTracks]={0};
  double zVtx[kTypeOfTracks]={0}, z[kTypeOfTracks-1]={0}, R[kTypeOfTracks-1]={0}; //Of the MCTrack


  //-------------------------GENERATED TRACKS-----------------------------------

  Int_t nEvents = kineTree->GetEntries();
  double charge;
  long pdgCode;

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
      pdgCode = thisTrackGen->GetPdgCode();
      if (TDatabasePDG::Instance()->GetParticle(pdgCode))
      {
        charge = TMath::Abs(TDatabasePDG::Instance()->GetParticle(pdgCode)->Charge());
      }
      else
      {
        //printf("strange pdgCode = %ld\n", pdgCode);
        charge = 0;
      }

      if ((thisTrackGen->isPrimary()) && (charge>0.1))//coupure sur l'origine de la trace et particule chargée
      {
        histPtVsEta[kGen]  ->Fill(eta[kGen],pt[kGen]);
        histPhiVsEta[kGen] ->Fill(eta[kGen],phi[kGen]);
        histPhiVsPt[kGen]  ->Fill(pt[kGen],phi[kGen]);
        histZvtxVsEta[kGen]->Fill(eta[kGen],zVtx[kGen]);
        histRVsZ[kGen]     ->Fill(z[kGen],R[kGen]);
      }
      if(thisTrackGen->isPrimary())
      {
        histPtOfPrimary ->Fill(pt[kGen]);
        histPOfPrimary ->Fill(thisTrackGen->GetP());
        if (TMath::Abs(pdgCode)==13){printf("Un muon ! %ld\n", pdgCode);}
        //if((pt[kGen]<0.005)&&(thisTrackGen->GetEta()<1000000000000000019884624838655.0))
        if((pt[kGen]<0.005)&&(charge>0.1))
        {
          histEtaOfVeryLowPPrimary->Fill(thisTrackGen->GetEta());
          histPxOfVeryLowPPrimary->Fill(thisTrackGen->GetStartVertexMomentumX());
          histPyOfVeryLowPPrimary->Fill(thisTrackGen->GetStartVertexMomentumY());
          histPzOfVeryLowPPrimary->Fill(thisTrackGen->GetStartVertexMomentumZ());
          histPdgOfVeryLowPPrimary->Fill(thisTrackGen->GetPdgCode());
        }
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
  std::vector<uint64_t> mcLabelTableRaw; //contains the raw values of the mcLabel which have at least one cluster
  std::vector<o2::MCCompLabel> mcLabelTable; //contains the mcLabel which have at least one cluster
  std::vector<int> nClusterPerLabel;
  std::vector<uint64_t>::iterator it; //iterator for the mcLabelTableRaw array



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
        int srcIDTrue, trkIDTrue, evnIDTrue;
        fake = false;
        Int_t iTrack = 0;
        for (auto& mftTrack : mMFTTracks)
        {//loop over the MFT tracks

          auto ncls = mftTrack.getNumberOfPoints();
          auto offset = mftTrack.getExternalClusterIndexOffset();
          std::map<uint64_t, int> mcLabels;
          for (int icls = 0; icls < ncls; ++icls)//cluster loop 1
          {
            auto clsEntry = mtrackExtClsIDs[offset + icls];
            auto cluster = clsVec[clsEntry];

            for (int ilabel = 0; ilabel < (clsLabels->getLabels(clsEntry)).size(); ++ilabel)
            {
              auto& clsLabel = (clsLabels->getLabels(clsEntry))[ilabel];
              clsLabel.get(trkIDTrue, evnIDTrue, srcIDTrue, fake);
              if (!clsLabel.isNoise())
              {
        	       mcLabels[clsLabel.getRawValue()] += 1;
              }
            }


          }//end of cluster loop 1

          histCompteurDeMCLabel->Fill(mcLabels.size());

          std::pair<uint64_t, int> entryWithMaxValue = findEntryWithLargestValue(mcLabels);



          Int_t thisEvnID = -1, thisSrcID = -1, thisTrkID = -1, thisEventIDLabel = -1, nPoints = mftTrack.getNumberOfPoints();
          bool thisFake=false;

          //cout << "Entry with highest value: "<< entryWithMaxValue.first << " = "<< entryWithMaxValue.second << endl;
          //printf("Number of points for this track = %d\n", nPoints);
          for (int icls = 0; icls < ncls; ++icls)
          {
            auto clsEntry = mtrackExtClsIDs[offset + icls];
            auto cluster = clsVec[clsEntry];
            for (int ilabel = 0; ilabel < (clsLabels->getLabels(clsEntry)).size(); ++ilabel)
            {
              auto& clsLabel = (clsLabels->getLabels(clsEntry))[ilabel];
              clsLabel.get(trkIDTrue, evnIDTrue, srcIDTrue, fake);



              if ((!clsLabel.isNoise()) && (clsLabel.getRawValue()==entryWithMaxValue.first))
              {
        	       if (((Float_t)(mcLabels[clsLabel.getRawValue()]) / (Float_t)(nPoints)) >= 0.8)
                 { // Must have at least 80% of its clusters from the same MC Track
        	          thisTrkID = trkIDTrue;
        	          thisSrcID = srcIDTrue;
        	          thisEvnID = evnIDTrue;
        	          thisEventIDLabel = icls;
                    thisFake = fake;

        	       }
              }
            }

          }//end of cluster loop 2

          //after that, if thisTrkID, thisSrcID, thisEvnID are still equal to -1, it means the track isn't considered a true track !

          if (thisTrkID==-1)
          {
            //printf("NO This track number %d in the reco track tree is not a true track\n", iTrack);
          }
          else
          {
            //printf("YES this track number %d in the reco track tree of trkIDTrue= %d is true\n", iTrack, thisTrkID);
            kineTree->GetEntry(thisEvnID);

            MCTrack* thisTrack =  &(mcTrkVec)[thisTrkID];
            //zVtx[kGen] = eventHeader->GetZ();

            //pt[kGen] = thisTrack->GetPt();
            eta[kGen] = -1*thisTrack->GetEta();
            phi[kGen] = thisTrack->GetPhi();

            eta[kRecoTrue]=-1*mftTrack.getEta();

            if (mftTrack.getPhi()>=TMath::Pi()/2)
            {
              phi[kRecoTrue]=-mftTrack.getPhi()+TMath::Pi()/2+2*TMath::Pi();
            }
            else
            {
              phi[kRecoTrue]=-mftTrack.getPhi()+TMath::Pi()/2;
            }

            //printf("phigen=%f, phirec=%f\n", phi[kGen], phi[kRecoTrue]);
            histPhiRecVsPhiGen->Fill(phi[kGen], phi[kRecoTrue]);
            histEtaRecVsEtaGen->Fill(eta[kGen], eta[kRecoTrue]);
            histPhiVsEta[kRecoTrue]->Fill(eta[kRecoTrue], phi[kRecoTrue]);

          }
          z[kReco]=mftTrack.getZ();
          eta[kReco]=-1*mftTrack.getEta();

          if (mftTrack.getPhi()>=TMath::Pi()/2)
          {
            phi[kReco]=-mftTrack.getPhi()+TMath::Pi()/2+2*TMath::Pi();
          }
          else
          {
            phi[kReco]=-mftTrack.getPhi()+TMath::Pi()/2;
          }

          histPhiVsEta[kReco]->Fill(eta[kReco], phi[kReco]);
          iTrack++;


          //_______end of purity code
        }
        //----------------------END OF RECO TRACKS

  TFile of(ofname, "RECREATE");//output file
  //Write everything in one output file
  of.cd();

  histCompteurDeMCLabel->Write();
  histPtOfPrimary->Write();
  histPOfPrimary->Write();

  histEtaOfVeryLowPPrimary->Write();
  histPxOfVeryLowPPrimary->Write();
  histPyOfVeryLowPPrimary->Write();
  histPzOfVeryLowPPrimary->Write();
  histPdgOfVeryLowPPrimary->Write();

  histPhiRecVsPhiGen->Write();
  histEtaRecVsEtaGen->Write();

  for (int i = 0; i < kTypeOfTracks ; i++)
    {
      histPhiVsEta[i]->Write();
      histPtVsEta[i] ->Write();
      histPhiVsPt[i] ->Write();


      if (i < kTypeOfTracks-2)//information only available for generated and trackable tracks
      {
        histZvtxVsEta[i]->Write();
        histRVsZ[i]     ->Write();
      }
    }
  of.Close();



//end of StudyMFTTracks
}


//______________________________________________________________________________
//______________________________________________________________________________


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

//______________________________________________________________________________

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
}//end of IsTrackTrackable



//______________________________________________________________________________


// Function to find the Entry
// with largest Value in a Map
std::pair<uint64_t, int> findEntryWithLargestValue(std::map<uint64_t, int> sampleMap)
{

	// Reference variable to help find
	// the entry with the highest value
	std::pair<uint64_t, int> entryWithMaxValue = std::make_pair(0, 0);

	// Iterate in the map to find the required entry
	std::map<uint64_t, int>::iterator currentEntry;
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


//______________________________________________________________________________

void BookHistos()
{

  histPhiRecVsPhiGen = new TH2D("histPhiRecVsPhiGen", "Phi Rec Vs Phi Gen of true reco tracks ", 24, 0., 2*TMath::Pi(), 24, 0., 2*TMath::Pi());
  histPhiRecVsPhiGen->SetXTitle((string("#phi of ")+nameOfTracks[kGen]).c_str());
  histPhiRecVsPhiGen->SetYTitle((string("#phi of ")+nameOfTracks[kRecoTrue]).c_str());
  histPhiRecVsPhiGen->Sumw2();

  histEtaRecVsEtaGen = new TH2D("histEtaRecVsEtaGen", "Eta Rec Vs Eta Gen of true reco tracks ", 35, 1.0, 4.5, 35, 1.0, 4.5);
  histEtaRecVsEtaGen->SetXTitle((string("#eta of ")+nameOfTracks[kGen]).c_str());
  histEtaRecVsEtaGen->SetYTitle((string("#eta of ")+nameOfTracks[kRecoTrue]).c_str());
  histEtaRecVsEtaGen->Sumw2();

  histCompteurDeMCLabel = new TH1D("histCompteurDeMCLabel","Compteur de MCLabel", 100, 0, 10);
  histCompteurDeMCLabel->SetXTitle("nb de mcLabel par trace");

  histPtOfPrimary = new TH1D("histPtOfPrimary","Transverse momentum of primary tracks", 1000, 0, 1.0);
  histPtOfPrimary->SetXTitle("p_{T} (GeV/c)");

  histPOfPrimary = new TH1D("histPOfPrimary","Momentum of primary tracks", 1000, 0, 1.0);
  histPOfPrimary->SetXTitle("p (GeV/c)");

  histEtaOfVeryLowPPrimary = new TH1D("histEtaOfVeryLowPPrimary","#eta of very low momentum primary tracks", 1000,-20, 20);
  histEtaOfVeryLowPPrimary->SetXTitle("#eta");

  histPxOfVeryLowPPrimary = new TH1D("histPxOfVeryLowPPrimary","#p_x of very low momentum primary tracks", 1000,-10, 10);
  histPxOfVeryLowPPrimary->SetXTitle("#p_{x}");

  histPyOfVeryLowPPrimary = new TH1D("histPyOfVeryLowPPrimary","#p_x of very low momentum primary tracks", 1000,-10, 10);
  histPyOfVeryLowPPrimary->SetXTitle("#p_{y}");

  histPzOfVeryLowPPrimary = new TH1D("histPzOfVeryLowPPrimary","#p_x of very low momentum primary tracks", 2000,-900, 900);
  histPzOfVeryLowPPrimary->SetXTitle("#p_{z}");

  histPdgOfVeryLowPPrimary = new TH1I("histPdgOfVeryLowPPrimary","#p_x of very low momentum primary tracks", 4000,-2000, 5000);
  histPdgOfVeryLowPPrimary->SetXTitle("pdgCode");

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

    if (i < kTypeOfTracks-2)//information only available for generated and trackable tracks
    {
      //histRVsZ
      histRVsZ[i] = new TH2D((string("histRVsZ")+nameOfTracks[i]).c_str(), (string("R Vs Z of ")+nameOfTracks[i]).c_str(), 400, -80., 20., 400, 0., 80.);
      histRVsZ[i]->SetXTitle((string("z (cm) origin of ")+nameOfTracks[i]).c_str());
      histRVsZ[i]->SetYTitle((string("R (cm) radius of origin of ")+nameOfTracks[i]).c_str());
      histRVsZ[i]->Sumw2();

    }
  }
}
