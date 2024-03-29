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
#include "SensitiveDetector.hh"
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
  G4VPhysicalVolume* GetRingsSupportPlus() {return fPhysicRingsSupportPlus;}
  G4VPhysicalVolume* GetRingsSupportMinus() {return fPhysicRingsSupportMinus;}
  G4VPhysicalVolume* GetResistorsPlus() {return fPhysicResistorsPlus;}
  G4VPhysicalVolume* GetResistorsMinus() {return fPhysicResistorsMinus;}
  G4VPhysicalVolume* GetRingStripsPlus() {return fPhysicRingStripsPlus;}
  G4VPhysicalVolume* GetRingStripsMinus() {return fPhysicRingStripsMinus;}
  G4VPhysicalVolume* GetVessel() {return fPhysicVessel;} 

  G4double GetCathodeWidth() {return fCathodeWidth;}
  G4double GetGEMOuterWidth() {return fGEMOuterWidth;}
  G4double GetGEMInnerWidth() {return fGEMCoreWidth;}
  G4double GetRingsSupportWidth() {return fRingSupportWidth;}
  G4double GetRingStripsWidth() {return fRingStripWidth;}
  G4double GetResistorWidth() {return fResistorWidth;}
  G4double GetLensWidth() {return fLensWidth;}
  G4double GetSensorWidth() {return fSensorWidth;}
  G4double GetVesselWidth() {return fVesselWidth;}
  
  G4PhysicalVolumeStore* GetVolumeStored() {return fPhysVolStore;}

  std::vector<G4String> GetCathodesList() {return fListCathodes;}
  std::vector<G4String> GetGEMsOuterLists() {return fListGEMsOuter;}
  std::vector<G4String> GetGEMsInnerLists() {return fListGEMsCore;}
  std::vector<G4String> GetRingsSupportList() {return fListSupportRings;}
  std::vector<G4String> GetRingStripstList() {return fListRingStrips;}
  std::vector<G4String> GetResistorList() {return fListResistors;}
  std::vector<G4String> GetLensList() {return fListLens;}
  std::vector<G4String> GetSensorsList() {return fListSensors;}

  SensitiveDetector* GetSensitiveDetector(){return fSensitiveDetector;}
  
  private:
  
  G4double fWorldSize_x;
  G4double fWorldSize_y;
  G4double fWorldSize_z;

  G4double fCathodeWidth;
  G4double fGEMOuterWidth;
  G4double fGEMCoreWidth;
  G4double fRingSupportWidth;
  G4double fRingStripWidth;
  G4double fResistorWidth;
  G4double fLensWidth;
  G4double fSensorWidth;
  G4double fVesselWidth;
  
  std::vector<G4String> fListCathodes;
  std::vector<G4String> fListGEMsOuter;
  std::vector<G4String> fListGEMsCore;
  std::vector<G4String> fListSupportRings;
  std::vector<G4String> fListRingStrips;
  std::vector<G4String> fListResistors;
  std::vector<G4String> fListLens;
  std::vector<G4String> fListSensors;
  std::vector<G4String> fListDetector;
  
  G4VPhysicalVolume* fPhysicalCathodes;
  G4VPhysicalVolume* fPhysicGEMsPlus;
  G4VPhysicalVolume* fPhysicGEMsMinus;
  G4VPhysicalVolume* fPhysicGEMsCorePlus;
  G4VPhysicalVolume* fPhysicGEMsCoreMinus;
  G4VPhysicalVolume* fPhysicRingsSupportPlus;
  G4VPhysicalVolume* fPhysicRingsSupportMinus;
  G4VPhysicalVolume* fPhysicRingStripsPlus;
  G4VPhysicalVolume* fPhysicRingStripsMinus;
  G4VPhysicalVolume* fPhysicResistorsMinus;
  G4VPhysicalVolume* fPhysicResistorsPlus;
  G4VPhysicalVolume* fPhysicVessel;
  G4VPhysicalVolume* fPhysicLensPlus;
  G4VPhysicalVolume* fPhysicLensMinus;
  G4VPhysicalVolume* fPhysicSensorsPlus;
  G4VPhysicalVolume* fPhysicSensorsMinus;
  
  
  G4PhysicalVolumeStore* fPhysVolStore;

  G4LogicalVolume* fLogicalGasVolume;

  SensitiveDetector* fSensitiveDetector;
  
  virtual void ConstructSDandField();
  
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

