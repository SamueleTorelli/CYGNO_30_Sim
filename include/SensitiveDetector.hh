#ifndef DETECTOR_HH
#define DETECTOR_HH

#include "G4VSensitiveDetector.hh"
#include "G4Step.hh"
#include "G4TouchableHistory.hh"
#include "RunAction.hh"

class SensitiveDetector : public G4VSensitiveDetector
{

public:
  SensitiveDetector(G4String);
  ~SensitiveDetector();

private:
  virtual G4bool ProcessHits(G4Step *, G4TouchableHistory*);

};


#endif