
void ProcessFile(std::string filename){

  TFile* f = TFile::Open(filename.c_str());
  
  TTree* tree = (TTree*)f->Get("Hits");

  Int_t Evn,ParticleID,ParticleTag,ParentID,VolumeNumber;
  Double_t x_hits,y_hits,z_hits,VolumeTraslX,VolumeTraslY,VolumeTraslZ,EnergyDeposit;

  Char_t Nucleus[10];

  Char_t ParticleName[10];

  tree->SetBranchAddress("EventNumber",&Evn);
  tree->SetBranchAddress("ParticleID",&ParticleID);
  tree->SetBranchAddress("ParticleTag",&ParticleTag);
  tree->SetBranchAddress("ParticleName",&ParticleName);
  tree->SetBranchAddress("ParentID",&ParentID);
  tree->SetBranchAddress("VolumeNumber",&VolumeNumber);
  tree->SetBranchAddress("x_hits",&x_hits);
  tree->SetBranchAddress("y_hits",&y_hits);
  tree->SetBranchAddress("z_hits",&z_hits);
  tree->SetBranchAddress("EnergyDeposit",&EnergyDeposit);
  tree->SetBranchAddress("VolumeNumber",&VolumeNumber);
  tree->SetBranchAddress("VolumeTraslX",&VolumeTraslX);
  tree->SetBranchAddress("VolumeTraslY",&VolumeTraslY);
  tree->SetBranchAddress("VolumeTraslZ",&VolumeTraslZ);
  tree->SetBranchAddress("Nucleus",&Nucleus);

  TFile* fout= new TFile("test.root","recreate");

  fout->cd();
  
  Int_t Out_evNumber=0;
  std::string Out_PartName;
  double_t Out_TotalEdep=0;
  Int_t Out_VolNum;
  Int_t Out_PartID;
  Int_t Out_PartentID;
  
  TTree* outTree = new TTree("elabHits","elabHits");

  outTree->Branch("evNumber",&Out_evNumber);
  outTree->Branch("PartName",&Out_PartName);
  outTree->Branch("Edep",&Out_TotalEdep);
  outTree->Branch("PartID",&Out_PartID);
  outTree->Branch("ParentID",&Out_ParentID);
  outTree->Branch("VolNum",&Out_VolNum);
    
  tree->GetEntry(0);
  
  Out_evNumber=Evn;
  Out_PartName=ParticleName;
  Out_VolNum=VolumeNumber;
  
  for(int i=0;i<4000;i++){
    tree->GetEntry(i);

    //if neutrino CONTINUE
    if(strcmp(ParticleName,"anti_nu_e")==0 || strcmp(ParticleName,"anti_nu_e")==0){
      continue;
    }else if(Evn == Out_evNumber && Out_PartName==ParticleName && Out_VolNum==VolumeNumber){
      Out_TotalEdep+=EnergyDeposit;
    } else {
      outTree->Fill();
      Out_evNumber=Evn;
      Out_PartName=ParticleName;
      Out_VolNum=VolumeNumber;
      Out_TotalEdep=0;
    }
    
  }//chiudo for

  outTree->Write();
  fout->Save();
  fout->Close();

}
