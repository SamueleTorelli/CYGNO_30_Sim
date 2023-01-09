#include "Detector.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"
#include "G4RunManager.hh"
#include "G4AnalysisManager.hh"


SensitiveDetector::SensitiveDetector(G4String name) :
  G4VSensitiveDetector(name)
{}

SensitiveDetector::~SensitiveDetector()
{}

G4bool SensitiveDetector::ProcessHits(G4Step * aStep, G4TouchableHistory* Rohist)
{

  G4Track* track = aStep->GetTrack();

  G4StepPoint* preStepPoint = aStep->GetPreStepPoint();
  G4StepPoint* postStepPoint = aStep->GetPostStepPoint();

  G4ThreeVector posParticle = preStepPoint->GetPosition();

  G4String particleName = track->GetParticleDefinition()->GetParticleName();
  G4int particleID = track->GetTrackID();
  G4double EdepStep = aStep->GetTotalEnergyDeposit();
  G4int VolumeCopyNumber = track->GetVolume()->GetCopyNo();
  G4int particleParentID = track->GetParentID();
  
  G4int particleTag=-1;

  if(particleName == "e-"){
    particleTag=1;
  } else if(particleName == "e+"){
    particleTag=2;
  } else if(particleName == "gamma"){
    particleTag=0;
  } else{
    particleTag=-1;
  }
  
  G4cout << "position of: " << particleName <<" " << track->GetTrackID() << "  is:  "<< posParticle << "Energy deposited:  " << EdepStep << "  in volume:  " << VolumeCopyNumber << " ParentID: "  << track->GetParentID()  << G4endl;

  G4int evt = G4RunManager::GetRunManager()->GetCurrentEvent()->GetEventID();
  
  G4AnalysisManager* AnalysisManager = G4AnalysisManager::Instance(); 

  AnalysisManager->FillNtupleIColumn(0,evt);
  AnalysisManager->FillNtupleSColumn(1,particleName);
  AnalysisManager->FillNtupleIColumn(2,particleID);
  AnalysisManager->FillNtupleIColumn(3,particleTag);
  AnalysisManager->FillNtupleIColumn(4,particleParentID);
  AnalysisManager->FillNtupleDColumn(5,posParticle[0]);
  AnalysisManager->FillNtupleDColumn(6,posParticle[1]);
  AnalysisManager->FillNtupleDColumn(7,posParticle[2]);
  AnalysisManager->FillNtupleDColumn(8,EdepStep);
  AnalysisManager->FillNtupleIColumn(9,VolumeCopyNumber);
  
  AnalysisManager->AddNtupleRow(0);


}
