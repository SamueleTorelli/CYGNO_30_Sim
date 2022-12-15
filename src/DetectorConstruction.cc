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
/// \file DetectorConstruction.cc
/// \brief Implementation of the DetectorConstruction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4PVReplica.hh"
#include "G4VPVParameterisation.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"
#include "G4SubtractionSolid.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():G4VUserDetectorConstruction()
{
  fWorldSize_x = 40*m;
  fWorldSize_y = 1*m;
  fWorldSize_z = 3*m;

  fCathodeSize_x = 50*cm;
  fCathodeSize_y = 80*cm;
  fCathodeSize_z = 0.0009*mm;

  fGEMSize_x = 50*cm;
  fGEMSize_y = 80*cm;
  fGEMSize_z = 0.06*mm;
  
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::~DetectorConstruction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* DetectorConstruction::Construct()
{
  //
  // define a material
  //

  G4NistManager* nist = G4NistManager::Instance();
  
  G4Material* Air =
  nist->FindOrBuildMaterial("G4_AIR"); 

  G4Material* Alluminium =
  nist->FindOrBuildMaterial("G4_Al");

  G4Material* Copper =
  nist->FindOrBuildMaterial("G4_Cu");

  //
  //defining PMMA
  //

  std::vector<G4int> natoms;
  std::vector<G4String> elements;

  elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);

  G4double density = 1.190*g/cm3;

  G4Material* PMMA = nist->ConstructNewMaterial("PMMA", elements, natoms, density);
  
  //     
  // World
  //
  
  G4Box*  
  solidWorld = new G4Box("World",                          //its name
                   fWorldSize_x/2,fWorldSize_y/2,fWorldSize_z/2);//its size
                   
  G4LogicalVolume*                         
  logicWorld = new G4LogicalVolume(solidWorld,             //its solid
                                   Air,                    //its material
                                   "World");               //its name
  G4VPhysicalVolume*                                   
  physiWorld = new G4PVPlacement(0,                      //no rotation
                                 G4ThreeVector(),        //at (0,0,0)
                                 logicWorld,             //its logical volume
                                 "World",                //its name
                                 0,                      //its mother  volume
                                 false,                  //no boolean operation
                                 0);                     //copy number


  //
  //Cathode
  //

  G4Colour AlColor(0.69, 0.77, 1.00);
  G4VisAttributes* cathodeVisAttributes = new G4VisAttributes(AlColor);
  cathodeVisAttributes->SetForceSolid(true);

  
  G4Box*
    solidCathode = new G4Box("Cathode",
			fCathodeSize_x/2,fCathodeSize_y/2,fCathodeSize_z);

  G4LogicalVolume*
    logicCathode = new G4LogicalVolume(solidCathode,
				       Alluminium,
				       "Cathode");
  
  logicCathode->SetVisAttributes(cathodeVisAttributes);
    
  G4double detectorSpace = 1.*cm;
     
  for(G4int i=-37;i<38;i++){
  
    G4VPhysicalVolume*
      physicCathode = new G4PVPlacement(0,
					G4ThreeVector(i*(fCathodeSize_x+detectorSpace),0,0),
					logicCathode,
					"Cathode",
					logicWorld,
					false,
					i+37);
    
    std::cout << "Cathode index: " << i+37 << std::endl;
  }


  
  
  //
  //GEMs
  //
  

  G4Colour CuColor(0.45,0.25,0.0,0.9);
  G4VisAttributes* GEMVisAttributes = new G4VisAttributes(CuColor);
  GEMVisAttributes->SetForceSolid(true);
  
  G4Box*
    solidGEM = new G4Box("GEM",
			fGEMSize_x/2,fGEMSize_y/2,fGEMSize_z);

  G4LogicalVolume*
    logicGEM = new G4LogicalVolume(solidGEM,
				   Copper,
				   "Cathode");

  logicGEM->SetVisAttributes(GEMVisAttributes);
  
  G4double GEMGap = 0.2*cm;
  G4double GEMDistanceFromCathode = 50*cm;

  for(G4int i =-37;i<38;i++){

    //GEMs on positive side of z axis

    for(G4int j=0;j<3;j++){
      G4VPhysicalVolume* 
	physicGEMs = new G4PVPlacement(0,
				       G4ThreeVector(i*(fCathodeSize_x+detectorSpace),0,GEMDistanceFromCathode+j*GEMGap),
				       logicGEM,
				       "GEMs",
				       logicWorld,
				       false,
				       (37+i)*10+abs(j+1));
      std::cout << "North: " << (37+i)*10+abs(j+1) << std::endl;
    }//chiudo for su j
    
    //GEMs on negative side of z axis
    
    for(G4int j=0;j>-3;j--){
      G4VPhysicalVolume* 
	physicGEMs = new G4PVPlacement(0,
				       G4ThreeVector(i*(fCathodeSize_x+detectorSpace),0,-1*GEMDistanceFromCathode+j*GEMGap),
				       logicGEM,
				       "GEMs",
				       logicWorld,
				       false,
				       -1*(37+i)*10-abs(j-1));
      std::cout << "South: " << -1*(37+i)*10-abs(j-1) << std::endl; 
    }//chiudo for su j
    
  }//chiudo for su i


  //
  //Field ring
  //

  G4Colour CaptonColor(0.2,1.0,0.0,0.2);
  G4VisAttributes* RingsVisAttributes = new G4VisAttributes(CaptonColor);
  RingsVisAttributes->SetForceSolid(true);

  
  G4double Ring_z = 1*cm; //depth of the ring in Z
  G4double Ring_width = 0.2*cm; // width of the ring in x-y
  G4int NRings = 32;

  G4double ringspacing = (GEMDistanceFromCathode-NRings*Ring_z)/(NRings+1);
  
  G4Box*
    outerShapeRing = new G4Box("RingOuterShape",
			       fGEMSize_x/2,fGEMSize_y/2,Ring_z/2);

  G4Box*
    innerShapeRing = new G4Box("RingOuterShape",
			 fGEMSize_x/2-Ring_width/2,fGEMSize_y/2-Ring_width/2,Ring_z/2+0.3*cm);

  G4VSolid* solidRing =
    new G4SubtractionSolid("Ring",outerShapeRing,innerShapeRing);

  G4LogicalVolume*
    logicRing = new G4LogicalVolume(solidRing,
				   Copper,
				   "Ring");

  logicRing->SetVisAttributes(RingsVisAttributes);

  for(G4int i=-37;i<38;i++){

    for(G4int j=0;j<NRings;j++){

      //Rings on positive side of z axis
      
      G4VPhysicalVolume* 
	physicRings = new G4PVPlacement(0,
					G4ThreeVector(i*(fCathodeSize_x+detectorSpace),0,(j+1)*ringspacing+0.5*Ring_z+j*Ring_z),
					logicRing,
					"Rings",
					logicWorld,
					false,
					(37+i)*100+j);
      
    }

    //Rings on positive side of z axis
    
    for(G4int j=0;j<NRings;j++){
    
      G4VPhysicalVolume* 
	physicRings = new G4PVPlacement(0,
					G4ThreeVector(i*(fCathodeSize_x+detectorSpace),0,-1*((j+1)*ringspacing+0.5*Ring_z+j*Ring_z)),
					logicRing,
					"Rings",
					logicWorld,
					false,
					-1*((37+i)*100+j));
      
    }
  }


  //
  //PMMA Vessel
  //

  G4double Vesselwidth = 1*cm;

  G4double VesselSize_x_outer = (fCathodeSize_x/2 + 37*(+fCathodeSize_x + detectorSpace) + 2*cm);
  G4double VesselSize_y_outer = 84*cm/2;
  G4double VesselSize_z_outer = 104*cm/2;

  G4double VesselSize_x_inner = VesselSize_x_outer-Vesselwidth;
  G4double VesselSize_y_inner = VesselSize_y_outer-Vesselwidth;
  G4double VesselSize_z_inner = VesselSize_z_outer-Vesselwidth;

  G4Colour VesselColor(1,1,1,0.1);
  G4VisAttributes* VesselVisAttributes = new G4VisAttributes(VesselColor);
  VesselVisAttributes->SetForceSolid(true);


  G4Box*
    outerShapeVessel = new G4Box("VesselOuterShape",
				 VesselSize_x_outer,VesselSize_y_outer,VesselSize_z_outer);

  G4Box*
    innerShapeVessel = new G4Box("VesselOuterShape",
				 VesselSize_x_inner,VesselSize_y_inner,VesselSize_z_inner);
  
  G4VSolid* solidVessel =
    new G4SubtractionSolid("Vessel",outerShapeVessel,innerShapeVessel);

  G4LogicalVolume*
    logicVessel = new G4LogicalVolume(solidVessel,
				   PMMA,
				   "Vessel");

  logicVessel->SetVisAttributes(VesselVisAttributes);
  
  G4VPhysicalVolume* 
    physicVessel = new G4PVPlacement(0,
				     G4ThreeVector(0,0,0),
				     logicVessel,
				     "Vessel",
				     logicWorld,
				     false,
				     0);
  
  
  //
  //always return the physical World
  //  
  
  return physiWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
