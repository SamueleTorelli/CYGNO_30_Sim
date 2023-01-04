//
// ********************************************************************
// * License and Disclaimer                                           *
// *                                                                  *
// * The  Geant4 software  is  copyright of the Copyright Holders  of *
// * the Geant4 Collaboration.  It is provided  under  the terms  and *
// * conditions of the Geant4 Software License,  included in the file *
// * LICENSE and available at  http://cern.ch/geant4/license .  These *
// * include a list of copyright holders.                             *
// *                                                                  *
// * Neither the authors of this software system, nor their employing *
// * institutes,nor the agencies providing financial support for this *
// * work  make  any representation or  warranty, express or implied, *
// * regarding  this  software system or assume any liability for its *
// * use.  Please see the license in the file  LICENSE  and URL above *
// * for the full disclaimer and the limitation of liability.         *
// *                                                                  *
// * This  code  implementation is the result of  the  scientific and *
// * technical work of the GEANT4 collaboration.                      *
// * By using,  copying,  modifying or  distributing the software (or *
// * any work based  on the software)  you  agree  to acknowledge its *
// * use  in  resulting  scientific  publications,  and indicate your *
// * acceptance of all terms of the Geant4 Software license.          *
// ********************************************************************
//
/// \file PrimaryGeneratorAction.cc
/// \brief Implementation of the PrimaryGeneratorAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
#include "PrimaryGeneratorAction.hh"
#include "G4LogicalVolume.hh"
#include "G4Event.hh"
#include "G4ParticleTable.hh"
#include "G4IonTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4Geantino.hh"
#include "G4ParticleTable.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4VSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::PrimaryGeneratorAction(DetectorConstruction* Detector)
  : G4VUserPrimaryGeneratorAction(),
    fParticleGun(0),
    fDetector(Detector)
{
  
  G4int n_particle = 1;
  fParticleGun  = new G4ParticleGun(n_particle);
  
  G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
  G4ParticleDefinition* particle = particleTable->FindParticle("geantino");  
    
  double randNum = G4UniformRand();

  G4cout << "myTest " << fDetector->GetCathodesVolumes()->GetCopyNo() << G4endl;
  G4cout << "myTest2 " << fDetector->GetCathodesVolumes()->GetMultiplicity() << G4endl;
  G4cout << "myTest3 " <<  randNum << G4endl;
  
  fParticleGun->SetParticleEnergy(0*eV);
  //fParticleGun->SetParticlePosition(GetPointOnDetectorElement("Cathodes"));
  fParticleGun->SetParticleMomentumDirection(G4ThreeVector(1.,0.,0.));          
  fParticleGun->SetParticleDefinition(particle);

  G4cout << "here" << G4endl;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void PrimaryGeneratorAction::GeneratePrimaries(G4Event* anEvent)
{

  if (fParticleGun->GetParticleDefinition() == G4Geantino::Geantino()) {  
    G4int Z = 29, A = 64;
    G4double ionCharge   = 0.*eplus;
    G4double excitEnergy = 0.*keV;
    
    G4ParticleDefinition* ion
      = G4IonTable::GetIonTable()->GetIon(Z,A,excitEnergy);
    fParticleGun->SetParticleDefinition(ion);
    fParticleGun->SetParticleCharge(ionCharge);
  }    
  
  //create vertex
  //   
  fParticleGun->GeneratePrimaryVertex(anEvent);
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


G4ThreeVector PrimaryGeneratorAction::GetPointOnDetectorElement(G4String El){

  std::vector<G4String> Elements;
    
  if(El == "Cathodes"){
    Elements = fDetector->GetCathodesList();
  } else if (El == "GEMs"){
    Elements = fDetector->GetGEMsLists();
  } else if (El == "Rings"){
    Elements = fDetector->GetRingsList();
  } else if (El == "Vessel"){
    Elements.push_back("0");
  } else if (El == "Lens"){
    Elements = fDetector->GetLensList();
  } else if (El == "Sensors"){
    Elements = fDetector->GetSensorsList();
  }

  /*  for(int i = 0 ; i< Elements.size();i++){
    G4cout <<"element " << i << " "  <<Elements[i] << G4endl;
    
    } */ 

  if(Elements[0]!="0"){

    G4int min = 0;
    G4int max = Elements.size();
    
    G4int nEl = min + (int)(G4UniformRand() * (max - min));

    G4cout<< "ElementNumber_____ " << nEl << G4endl;
    G4cout<< "Element " << Elements[nEl] << G4endl;
    
    G4VPhysicalVolume* vol = fDetector->GetVolumeStored()->GetVolume(Elements[nEl]);
    
    G4ThreeVector PointOnSurface = vol->GetLogicalVolume()->GetSolid()->GetPointOnSurface();

    G4ThreeVector TranslationVolume = vol->GetObjectTranslation();

    G4ThreeVector Point = PointOnSurface + TranslationVolume;
    
    G4cout <<"PointOnSurface " << PointOnSurface << G4endl;
    G4cout <<"TranslationVolume " << TranslationVolume << G4endl;
    
    return Point;
    
  } else {
    
    G4VPhysicalVolume* vol = fDetector->GetVessel();
    
    G4ThreeVector PointOnSurface = vol->GetLogicalVolume()->GetSolid()->GetPointOnSurface();
    G4ThreeVector TranslationVolume = vol->GetObjectTranslation();
    
    G4ThreeVector Point = PointOnSurface + TranslationVolume;

    return Point;
    
  }

  
}
