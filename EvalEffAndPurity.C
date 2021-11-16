//Macro EvalEffAndPurity.C evaluates the purity, the acceptance and the efficiency of the MFT
//example of how to use:
//EvalEffAndPurity("./08-11-2021-18:00/histoResults_08-11-2021-18:00/outputfile_studyTracks_merged.root")

TH2D *hPhiVsEtaTrackableInTmp=0, *hPhiVsEtaGenInTmp=0, *hPhiVsEtaEfficiency=0, *hPhiVsEtaAcceptance=0;

TH2D *hPhiVsEtaRecoTrueInTmp=0, *hPhiVsEtaRecoInTmp=0, *hPhiVsEtaPurity=0;

enum typeOfTracks {kGen, kTrackable, kReco, kRecoTrue, kTypeOfTracks};

std::vector<string> nameOfTracks = {"Gen","Trackable","Rec","Rec and True"};//of size kTypeOfTracks

void EvalEffAndPurity(const Char_t *ifname = "outputfile_studyTracks_merged.root", const Char_t *ofnameStart = "EffAccPurity")
{
  TFile *fileInTmp = new TFile(ifname);

  const Char_t *ofname = (string(ofnameStart)+".root").c_str();
  TFile saveFile(ofname,"recreate");

  //-------------------------------PURITY---------------------------------------
  hPhiVsEtaRecoTrueInTmp = (TH2D*) fileInTmp->Get((string("histPhiVsEta")+nameOfTracks[kRecoTrue]).c_str());
  hPhiVsEtaRecoInTmp = (TH2D*) fileInTmp->Get((string("histPhiVsEta")+nameOfTracks[kReco]).c_str());


  auto cnv_Purity = new TCanvas("cnv_Purity","",800,800);
  hPhiVsEtaPurity = (TH2D*)hPhiVsEtaRecoTrueInTmp->Clone("hPhiVsEtaPurity");

  hPhiVsEtaPurity->SetTitle("MFT tracking purity as a function of #eta and #phi;#eta;#phi");
  hPhiVsEtaPurity->Divide(hPhiVsEtaRecoInTmp);
  cnv_Purity->SetLogz();
  hPhiVsEtaPurity->Draw("COLZ");
  cnv_Purity->Update();
  cnv_Purity->Draw();

  //----------------------------ACCEPTANCE--------------------------------------

  hPhiVsEtaTrackableInTmp = (TH2D*) fileInTmp->Get((string("histPhiVsEta")+nameOfTracks[kTrackable]).c_str());
  hPhiVsEtaRecoInTmp = (TH2D*) fileInTmp->Get((string("histPhiVsEta")+nameOfTracks[kReco]).c_str());
  hPhiVsEtaGenInTmp = (TH2D*) fileInTmp->Get((string("histPhiVsEta")+nameOfTracks[kGen]).c_str());


  hPhiVsEtaEfficiency = (TH2D*)hPhiVsEtaRecoInTmp->Clone("hPhiVsEtaEfficiency");
  hPhiVsEtaEfficiency->SetTitle("MFT tracking efficiency as a function of #eta and #phi;#eta;#phi");
  hPhiVsEtaEfficiency->Divide(hPhiVsEtaTrackableInTmp);


  auto cnv_Efficiency = new TCanvas("cnv_Efficiency","",800,800);
  hPhiVsEtaEfficiency->Draw("COLZ");
  cnv_Efficiency->Update();
  cnv_Efficiency->Draw();

  //---------------------------ACCEPTANCE---------------------------------------

  hPhiVsEtaAcceptance = (TH2D*)hPhiVsEtaTrackableInTmp->Clone("hPhiVsEtaAcceptance");
  hPhiVsEtaAcceptance->SetTitle("MFT Acceptance as a function of #eta and #phi;#eta;#phi");
  hPhiVsEtaAcceptance->Divide(hPhiVsEtaGenInTmp);



  auto cnv_Acceptance = new TCanvas("cnv_Acceptance","",800,800);
  hPhiVsEtaAcceptance->Draw("colz");
  cnv_Acceptance->SetLogz();
  cnv_Acceptance->Update(); // needed to get the painted histogram

  cnv_Acceptance->Draw();


  gStyle->SetOptStat(0);

  auto cnv_ToPDF = new TCanvas("cnv_ToPDF","",800,800);
  hPhiVsEtaPurity->Draw("COLZ");
  cnv_ToPDF->Print((string(ofnameStart)+".pdf"+"(").c_str()); //write canvas and keep the pdf file open
  hPhiVsEtaEfficiency->Draw("COLZ");
  cnv_ToPDF->Print((string(ofnameStart)+".pdf").c_str());
  hPhiVsEtaAcceptance->Draw("colz");
  cnv_ToPDF->Print((string(ofnameStart)+".pdf").c_str());

  hPhiVsEtaRecoTrueInTmp->Draw("COLZ");
  cnv_ToPDF->Print((string(ofnameStart)+".pdf").c_str()); //canvas is added to the pdf file
  hPhiVsEtaRecoInTmp->Draw("COLZ");
  cnv_ToPDF->Print((string(ofnameStart)+".pdf").c_str());
  hPhiVsEtaTrackableInTmp->Draw("COLZ");
  cnv_ToPDF->Print((string(ofnameStart)+".pdf").c_str());
  hPhiVsEtaGenInTmp->Draw("COLZ");
  cnv_ToPDF->Print((string(ofnameStart)+".pdf"+")").c_str()); //canvas is added to file and pdf file is closed


  saveFile.cd();
  hPhiVsEtaPurity->Write();
  hPhiVsEtaEfficiency->Write();
  hPhiVsEtaAcceptance->Write();
  saveFile.Close();
}
