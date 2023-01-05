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
/// \file DetectorConstruction.hh
/// \brief Definition of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo...... 

#ifndef DetectorConstruction_h
#define DetectorConstruction_h 1

#include "G4VUserDetectorConstruction.hh"
#include "G4VPhysicalVolume.hh"
#include "G4PhysicalVolumeStore.hh"
#include "globals.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
class G4VPhysicalVolume;

class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:
  
    DetectorConstruction();
   ~DetectorConstruction();

    virtual     
    G4VPhysicalVolume* Construct();
                        
  G4double GetWorldSizeX() {return fWorldSize_x;};
  G4double GetWorldSizeY() {return fWorldSize_y;};
  G4double GetWorldSizeZ() {return fWorldSize_z;}; 

  G4VPhysicalVolume* GetCathodesVolumes() {return fPhysicalCathodes;}
  G4VPhysicalVolume* GetGEMVolumesPlus() {return fPhysicGEMsPlus;}
  G4VPhysicalVolume* GetGEMVolumesMinus() {return fPhysicGEMsMinus;}
  G4VPhysicalVolume* GetRingsVolumesPlus() {return fPhysicRingsPlus;}
  G4VPhysicalVolume* GetRingsVolumesMinus() {return fPhysicRingsMinus;}
  G4VPhysicalVolume* GetVessel() {return fPhysicVessel;} 

  G4PhysicalVolumeStore* GetVolumeStored() {return fPhysVolStore;}

  std::vector<G4String> GetCathodesList() {return fListCathodes;}
  std::vector<G4String> GetGEMsLists() {return fListGEMs;}
  std::vector<G4String> GetRingsList() {return fListRings;}
  std::vector<G4String> GetLensList() {return fListLens;}
  std::vector<G4String> GetSensorsList() {return fListSensors;}
  
  private:
  
  G4double fWorldSize_x;
  G4double fWorldSize_y;
  G4double fWorldSize_z;

  std::vector<G4String> fListCathodes;
  std::vector<G4String> fListGEMs;
  std::vector<G4String> fListRings;
  std::vector<G4String> fListLens;
  std::vector<G4String> fListSensors;
  std::vector<G4String> fListDetector;
  
  G4VPhysicalVolume* fPhysicalCathodes;
  G4VPhysicalVolume* fPhysicGEMsPlus;
  G4VPhysicalVolume* fPhysicGEMsMinus;
  G4VPhysicalVolume* fPhysicRingsPlus;
  G4VPhysicalVolume* fPhysicRingsMinus;
  G4VPhysicalVolume* fPhysicVessel;
  G4VPhysicalVolume* fPhysicLensPlus;
  G4VPhysicalVolume* fPhysicLensMinus;
  G4VPhysicalVolume* fPhysicSensorsPlus;
  G4VPhysicalVolume* fPhysicSensorsMinus;


  G4PhysicalVolumeStore* fPhysVolStore;

  G4LogicalVolume* fLogicalGasVolume;

  virtual void ConstructSDandField();
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

