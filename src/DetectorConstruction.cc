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
#include "G4PhysicalVolumeStore.hh"
#include "G4Tubs.hh"
#include "SensitiveDetector.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():G4VUserDetectorConstruction()
{
  fWorldSize_x = 14*m;
  fWorldSize_y = 3*m;
  fWorldSize_z = 3*m;  

  fListCathodes.clear();
  fListGEMs.clear();
  fListRings.clear();
  fListLens.clear();
  fListSensors.clear();
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
  //defining Glass
  //

  G4double GlassDensity = 2.65*g/cm3;
  
  G4Material* Glass = new G4Material("Silicon Dioxide", GlassDensity, 2);

  Glass->AddElement(nist->FindOrBuildElement("Si"), 1);
  Glass->AddElement(nist->FindOrBuildElement("O"), 2);
  
  //
  //defining PMMA
  //

  std::vector<G4int> natoms;
  std::vector<G4String> elements;

  elements.push_back("C");     natoms.push_back(5);
  elements.push_back("H");     natoms.push_back(8);
  elements.push_back("O");     natoms.push_back(2);

  G4double PMMADensity = 1.190*g/cm3;

  G4Material* PMMA = nist->ConstructNewMaterial("PMMA", elements, natoms, PMMADensity);

  //
  //defining sensor material
  //

  G4Material* Silicon = nist->FindOrBuildMaterial("G4_Si");    

  //
  //defining detector gas mixture
  //

  G4double aHe = 4.002602*g/mole;
  G4Element* elHe = new G4Element("Helium","He", 2, aHe);
  
  G4double aC = 12.0107*g/mole;
  G4Element* elC = new G4Element("Carbon", "C", 6, aC);
  
  G4double aF=18.998*g/mole;
  G4Element* elF = new G4Element("Flourine"  ,"F" , 9., aF);

  G4double He_frac = 0.6;
  G4double CF4_frac = 0.4;
  
  G4double densityHe = 162.488*He_frac*g/m3;
  G4double pressureHe = 1*He_frac*atmosphere;
  G4double temperatureHe = 300*kelvin;
  G4Material* He_gas = new G4Material("He_gas", densityHe, 1, kStateGas, temperatureHe, pressureHe);
  He_gas->AddElement(elHe, 1);

  //CF4_gas
  G4double densityCF4 = 3574.736*CF4_frac*g/m3;
  G4double pressureCF4 = 1*CF4_frac*atmosphere;
  G4double temperatureCF4 = 300*kelvin;
  G4Material* CF4_gas = new G4Material("CF4_gas", densityCF4, 2, kStateGas, temperatureCF4, pressureCF4);
  CF4_gas->AddElement(elC, 1);
  CF4_gas->AddElement(elF, 4);

  //CYGNO_gas
  G4double densityMix = He_gas->GetDensity()+CF4_gas->GetDensity();
  G4double pressureMix = He_gas->GetPressure()+CF4_gas->GetPressure();
  G4double temperatureMix = 300*kelvin;
  G4Material* CYGNO_gas = new G4Material("CYGNO_gas", densityMix, 2, kStateGas, temperatureMix, pressureMix);
  CYGNO_gas->AddMaterial(He_gas, He_gas->GetDensity()/densityMix*100*perCent);
  CYGNO_gas->AddMaterial(CF4_gas, CF4_gas->GetDensity()/densityMix*100*perCent);

  
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

  G4double CathodeSize_x = 50*cm;
  G4double CathodeSize_y = 80*cm;
  G4double CathodeSize_z = 0.0009*mm;

  fCathodeWidth=CathodeSize_z;
  
  G4Colour AlColor(0.69, 0.77, 1.00);
  G4VisAttributes* cathodeVisAttributes = new G4VisAttributes(AlColor);
  cathodeVisAttributes->SetForceSolid(true);

  
  G4Box*
    solidCathode = new G4Box("Cathode",
			CathodeSize_x/2,CathodeSize_y/2,CathodeSize_z);

  G4LogicalVolume*
    logicCathode = new G4LogicalVolume(solidCathode,
				       Alluminium,
				       "Cathode");
  
  logicCathode->SetVisAttributes(cathodeVisAttributes);
    
  G4double detectorSpace = 0.3*cm;
     
  for(G4int i=-12;i<13;i++){
    for(G4int j=-1;j<2;j++){
      
      fPhysicalCathodes = new G4PVPlacement(0,
					    G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace),0),
					    logicCathode,
					    "Cathode_"+std::to_string( (j+2)*100+(i+12) ),
					    logicWorld,
					    false,
					    i+37);
      
      fListCathodes.push_back("Cathode_"+std::to_string( (j+2)*100+(i+12)) );
      
      //std::cout << "Cathode index: " << "Cathode_"+std::to_string( (j+2)*100+(i+12) ) << "\n";
      
    }
  }
  
  //
  //GEMs
  //
  
  G4double GEMSize_x = 50*cm;
  G4double GEMSize_y = 80*cm;
  G4double GEMSize_z = 0.06*mm;

  fGEMWidth = GEMSize_z;

  G4Colour CuColor(0.45,0.25,0.0,0.9);
  G4VisAttributes* GEMVisAttributes = new G4VisAttributes(CuColor);
  GEMVisAttributes->SetForceSolid(true);
  
  G4Box*
    solidGEM = new G4Box("GEM",
			GEMSize_x/2,GEMSize_y/2,GEMSize_z);

  G4LogicalVolume*
    logicGEM = new G4LogicalVolume(solidGEM,
				   Copper,
				   "GEM");

  logicGEM->SetVisAttributes(GEMVisAttributes);
  
  G4double GEMGap = 0.2*cm;
  G4double GEMDistanceFromCathode = 50*cm;

  for(G4int i =-12;i<13;i++){

    //GEMs on positive side of z axis

    for(G4int j=0;j<3;j++){
      
      for(G4int k=-1;k<2;k++){
	
	fPhysicGEMsPlus = new G4PVPlacement(0,
					    G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace),GEMDistanceFromCathode+j*GEMGap),
					    logicGEM,
					    "GEM_"+std::to_string((i+13)*100+j*10+k+1),
					    logicWorld,
					    false,
					    (i+13)*100+j*10+k+1
					    );
	
	fListGEMs.push_back("GEM_"+std::to_string((i+13)*100+j*10+k+1));
	
	//std::cout << "North GEM: " << (i+13)*100+j*10+k+1 << "\n";

      }//chiudo for su k
      
    }//chiudo for su j


    
    //GEMs on negative side of z axis
    
    for(G4int j=0;j>-3;j--){
      
      for(G4int k=-1;k<2;k++){
	
	fPhysicGEMsMinus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace),-1*GEMDistanceFromCathode+j*GEMGap),
					     logicGEM,
					     "GEM_"+std::to_string((i+13)*100+j*10+k+4 ),
					     logicWorld,
					     false,
					     (i+13)*100+j*10+k+4
					     );
	
	fListGEMs.push_back("GEM_"+std::to_string((i+13)*100+j*10+k+4 ));
	
	//std::cout << "South GEM: " << (i+13)*100+j*10+k+4 << "\n"; 
      }//chiudo for su j

    }//chiudo for su k

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

  fRingWidth=Ring_width;
  
  G4Box*
    outerShapeRing = new G4Box("RingOuterShape",
			       GEMSize_x/2,GEMSize_y/2,Ring_z/2);

  G4Box*
    innerShapeRing = new G4Box("RingOuterShape",
			 GEMSize_x/2-Ring_width/2,GEMSize_y/2-Ring_width/2,Ring_z/2+0.3*cm);

  G4VSolid* solidRing =
    new G4SubtractionSolid("Ring",outerShapeRing,innerShapeRing);
  
  G4LogicalVolume*
    logicRing = new G4LogicalVolume(solidRing,
				   Copper,
				   "Ring");

  logicRing->SetVisAttributes(RingsVisAttributes);

  for(G4int i=-12;i<13;i++){

    for(G4int j=0;j<NRings;j++){

      //Rings on positive side of z axis
      for(G4int k=-1;k<2;k++){

	fPhysicRingsPlus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace),(j+1)*ringspacing+0.5*Ring_z+j*Ring_z),
					     logicRing,
					     "Ring_"+std::to_string((i+13)*1000+(j+10)*10+k+1),
					     logicWorld,
					     false,
					     (i+13)*1000+(j+10)*10+k+1
					     );
	
	fListRings.push_back("Ring_"+std::to_string( (i+13)*1000+(j+10)*10+k+1) );

	//std::cout << "Nord Rings: " << (i+13)*1000+(j+10)*10+k+1 << "\n";
	
      }//chiudo for k 
    }//chiudo for j
    
    //Rings on positive side of z axis
    
    for(G4int j=0;j<NRings;j++){
      
      for(G4int k=-1;k<2;k++){ 

	fPhysicRingsMinus = new G4PVPlacement(0,
					      G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace),-1*((j+1)*ringspacing+0.5*Ring_z+j*Ring_z)),
					      logicRing,
					      "Ring_"+std::to_string((i+13)*1000+(j+10)*10+k+4),
					      logicWorld,
					      false,
					      (i+13)*1000+(j+10)*10+k+4
					      );
      
	fListRings.push_back("Ring_"+std::to_string( (i+13)*1000+(j+10)*10+k+4) );
	
	//std::cout << "South Rings: " << (i+13)*1000+(j+10)*10+k+4 << "\n"; 

      }//chiudo for k
    }//chiudo for j
  }
  
  
  //
  //PMMA Vessel
  //

  G4double Vesselwidth = 1*cm;

  G4double VesselSize_x_outer = (CathodeSize_x/2 + 12*(+CathodeSize_x + detectorSpace) + 2*Vesselwidth);
  G4double VesselSize_y_outer = (CathodeSize_y/2 + (+CathodeSize_y + detectorSpace) + 2*Vesselwidth);
  G4double VesselSize_z_outer = (GEMDistanceFromCathode+2*GEMGap+3*GEMSize_z+ 2*Vesselwidth);

  G4double VesselSize_x_inner = VesselSize_x_outer-Vesselwidth;
  G4double VesselSize_y_inner = VesselSize_y_outer-Vesselwidth;
  G4double VesselSize_z_inner = VesselSize_z_outer-Vesselwidth;

  G4Colour VesselColor(1,1,1,0.1);
  G4VisAttributes* VesselVisAttributes = new G4VisAttributes(VesselColor);
  VesselVisAttributes->SetForceSolid(true);

  fVesselWidth=Vesselwidth;

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
  
  
  fPhysicVessel = new G4PVPlacement(0,
				    G4ThreeVector(0,0,0),
				    logicVessel,
				    "Vessel",
				    logicWorld,
				    true,
				    0
				    );

  //
  //Camera lenses
  //

  G4double LensDiameter = 1*cm;
  G4double LensThickness = 2*mm;
  G4double LensDistanceFromGEMs = 57.6*cm;
  G4double VerticalLensSpacing = 38*cm;

  fLensWidth=LensThickness;
  
  G4Tubs* solidLens = new G4Tubs("Lens", 0*mm, LensDiameter, LensThickness/2, 0, 360*deg);
  
  G4LogicalVolume*
    logicLens = new G4LogicalVolume(solidLens,
				    Glass,
				    "Lens"
				    );

  
  for(G4int i=-12;i<13;i++){
    
    for(G4int j=-1;j<2;j++){

      for(G4float k=-0.5;k<1;k++){

	fPhysicLensPlus = new G4PVPlacement(0,
					    G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace)+k*VerticalLensSpacing/2, GEMDistanceFromCathode+2*GEMGap+3*GEMSize_z + LensDistanceFromGEMs ),
					    logicLens,
					    "Lens_"+std::to_string((i+13)*1000+(j+1)*10+(k+0.5)),
					    logicWorld,
					    false,
					    (i+13)*1000+(j+1)*10+(k+0.5)
					    );

	fListLens.push_back("Lens_"+std::to_string((i+13)*1000+(j+1)*10+(k+0.5)));
	
      }//chiudo for k
      
    }//chiudo for j
    
    for(G4int j=-1;j<2;j++){
      
      for(G4float k=-0.5;k<1;k++){
	
	fPhysicLensMinus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace)+k*VerticalLensSpacing/2, -1*(GEMDistanceFromCathode+2*GEMGap+3*GEMSize_z + LensDistanceFromGEMs) ),
					    logicLens,
					    "Lens_"+std::to_string((i+13)*1000+(j+1+4)*10+(k+0.5)),
					    logicWorld,
					    false,
					    (i+13)*1000+(j+1)*10+(k+0.5)
					    );

	fListLens.push_back("Lens_"+std::to_string((i+13)*1000+(j+1+4)*10+(k+0.5)));
	
      }//chiudo for k
      
    }//chiudo for j

    
  }//chiudo for i

  
  //
  //Camera sensors
  //

  G4double SensorSize_x = 10.6*mm;
  G4double SensorSize_y = 18.8*mm;
  G4double SensorSize_z = 1*mm;
  G4double SensorDistanceFromLens = 6 * cm; 

  fSensorWidth=SensorSize_z;
  
  G4Box*
    solidSensor = new G4Box("Sensor",
			       SensorSize_x/2,SensorSize_y/2,SensorSize_z/2);

  G4LogicalVolume*
    logicSensor = new G4LogicalVolume(solidSensor,
				    Silicon,
				    "Sensor"
				    );


  for(G4int i=-12;i<13;i++){
    
    for(G4int j=-1;j<2;j++){
      
      for(G4float k=-0.5;k<1;k++){
	
	fPhysicSensorsPlus = new G4PVPlacement(0,
					    G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace)+k*VerticalLensSpacing/2, GEMDistanceFromCathode+2*GEMGap+3*GEMSize_z + LensDistanceFromGEMs + SensorDistanceFromLens ),
					    logicSensor,
					    "Sensor_"+std::to_string((i+13)*1000+(j+1)*10+(k+0.5)),
					    logicWorld,
					    false,
					    (i+13)*1000+(j+1)*10+(k+0.5)
					    );

	fListSensors.push_back("Sensor_"+std::to_string((i+13)*1000+(j+1)*10+(k+0.5)));
	
      }//chiudo for k
      
    }//chiudo for j
    
    for(G4int j=-1;j<2;j++){
      
      for(G4float k=-0.5;k<1;k++){
	
	fPhysicSensorsMinus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace)+k*VerticalLensSpacing/2, -1*(GEMDistanceFromCathode+2*GEMGap+3*GEMSize_z + LensDistanceFromGEMs + SensorDistanceFromLens) ),
					    logicSensor,
					    "Sensor_"+std::to_string((i+13)*1000+(j+1+4)*10+(k+0.5)),
					    logicWorld,
					    false,
					    (i+13)*1000+(j+1)*10+(k+0.5)
					    );

	fListSensors.push_back("Sensor_"+std::to_string((i+13)*1000+(j+1+4)*10+(k+0.5)));
	
      }//chiudo for k
      
    }//chiudo for j

    
  }//chiudo for i


  //
  //CF4 sensitive volume
  //


  G4Colour GasColor(0.0,0.0,1.0,0.1);
  G4VisAttributes* GasVisAttributes = new G4VisAttributes(GasColor);
  GasVisAttributes->SetForceSolid(true);
  
  G4Box* solidGasVolume = new G4Box("GasVolume",CathodeSize_x/2,CathodeSize_y/2,GEMDistanceFromCathode/2);

  fLogicalGasVolume = new G4LogicalVolume(solidGasVolume,
					  CYGNO_gas,
					  "GasVolume"
					  );

  fLogicalGasVolume->SetVisAttributes(GasVisAttributes);

  G4int counter = 0;
  
  for(G4int i=-12;i<13;i++){  
    for(G4int j=-1;j<2;j++){
      
      G4VPhysicalVolume* GasVolumePlus = new G4PVPlacement(0,
							   G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace),GEMDistanceFromCathode/2),
							   fLogicalGasVolume,
							   "GasVolume_"+std::to_string(counter),
							   logicWorld,
							   false,
							   (i+13)*100+(j+1)
							   );

      fListDetector.push_back("GasVolume_"+std::to_string(counter));
      counter++;
      
    }//chiudo for j
  }//chiudo for i
  
  for(G4int i=-12;i<13;i++){ 
    for(G4int j=-1;j<2;j++){
      
      G4VPhysicalVolume* GasVolumePlus = new G4PVPlacement(0,
							   G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace),-GEMDistanceFromCathode/2),
							   fLogicalGasVolume,
							   "GasVolume_"+std::to_string(counter),
							   logicWorld,
							   false,
							   (i+13)*100+(j+1)
							   );
      
      fListDetector.push_back("GasVolume_"+std::to_string(counter));
      counter++;
     
    }//chiudo for j
  }//chiudo secondo for i
  
  //
  //Creating the physicalvolumestore
  //
  
  fPhysVolStore = G4PhysicalVolumeStore::GetInstance();
  
  
  //
  //always return the physical World
  //  
  
  return physiWorld;
}


void DetectorConstruction::ConstructSDandField()
{

  SensitiveDetector* SensDet = new SensitiveDetector("SensitiveDetector");  
  fLogicalGasVolume->SetSensitiveDetector(SensDet);
  
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
