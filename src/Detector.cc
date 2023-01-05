#include "Detector.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "G4Track.hh"
#include "G4StepPoint.hh"
#include "G4ThreeVector.hh"

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

  G4ThreeVector posElectron = preStepPoint->GetPosition();

  G4cout << "position of: " << track->GetParticleDefinition()->GetParticleName() <<" " << track->GetTrackID() << "  is:  "<< track->GetPosition() << "Energy deposited:  " << aStep->GetTotalEnergyDeposit() << "  in volume:  " << track->GetVolume()->GetName()  << G4endl;

  
}
