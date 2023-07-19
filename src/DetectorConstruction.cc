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
#include <G4UnitsTable.hh>

#include "SensitiveDetector.hh"


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

DetectorConstruction::DetectorConstruction():G4VUserDetectorConstruction()
{
  fWorldSize_x = 14*m;
  fWorldSize_y = 3*m;
  fWorldSize_z = 3*m;  

  fListCathodes.clear();
  fListGEMsOuter.clear();
  fListGEMsCore.clear();
  fListSupportRings.clear();
  fListRingStrips.clear();
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
  //definition of Al2O3
  //

  G4Element* elAl = nist->FindOrBuildElement("Al");
  G4Element* elO = nist->FindOrBuildElement("O");
  
  G4Material* al2o3 = new G4Material("Al2O3", 3.97 * g/cm3, 2);
  al2o3->AddElement(elAl, 2);
  al2o3->AddElement(elO, 3);
    
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

  std::map<G4String,G4double> MassMap={
    {"Cathodes",0},
    {"SupportRings",0},
    {"RingStrips",0},
    {"Resistor",0},
    {"GEMsOuter",0},
    {"GEMsCore",0},
    {"Vessel",0},
    {"Lens",0},
    {"Sensors",0}
  };
  
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
  G4double CathodeSize_z = 0.05*cm;

  fCathodeWidth=CathodeSize_z;
  
  G4Colour CathodeColor(0.6, 0.4, 0.2,0.0);
  G4VisAttributes* cathodeVisAttributes = new G4VisAttributes(CathodeColor);
  cathodeVisAttributes->SetForceSolid(true);

  
  G4Box*
    solidCathode = new G4Box("Cathode",
			CathodeSize_x/2,CathodeSize_y/2,CathodeSize_z);

  G4LogicalVolume*
    logicCathode = new G4LogicalVolume(solidCathode,
				       Copper,
				       "Cathode");
  
  logicCathode->SetVisAttributes(cathodeVisAttributes);
    
  G4double detectorSpace = 0.4*cm; // original 0.3*cm

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
      MassMap["Cathodes"]+=logicCathode->GetMass();
      //std::cout << "Cathode index: " << "Cathode_"+std::to_string( (j+2)*100+(i+12) ) << "\n";
      
    }
  }
  
  //
  //GEMs outer
  //
  
  G4double GEMOuterSize_x = 50*cm;
  G4double GEMOuterSize_y = 80*cm;
  G4double GEMOuterSize_z = 0.05*mm;

  G4double GEMCoreSize_xSub = 50.01*cm;
  G4double GEMCoreSize_ySub = 80.01*cm;
  G4double GEMCoreSize_zSub = 0.005*mm;

  G4double GEMGap = 0.2*cm;
  G4double GEMDistanceFromCathode = 50*cm;

  
  fGEMOuterWidth = (GEMOuterSize_z-GEMCoreSize_zSub)/2;

  G4Colour CuColor(0.45,0.25,0.0,0.0);
  G4VisAttributes* GEMVisAttributes = new G4VisAttributes(CuColor);
  GEMVisAttributes->SetForceSolid(true);
  
  G4Box*
    solidOuterGEM = new G4Box("GEMOuter",
			GEMOuterSize_x/2,GEMOuterSize_y/2,GEMOuterSize_z/2);
  G4Box*
    solidCoreGEMSub = new G4Box("GEMCore",
			GEMCoreSize_xSub/2,GEMCoreSize_ySub/2,GEMCoreSize_zSub/2);

  G4VSolid* solidGEM =
    new G4SubtractionSolid("GEMSubSolid",solidOuterGEM,solidCoreGEMSub,0,G4ThreeVector( 0,0,-(GEMCoreSize_zSub/2)+0.095*(GEMCoreSize_zSub/2) )); //Translation inserted by hand due to some bug in G4SubtractionSolid
  
  G4LogicalVolume*
    logicGEM = new G4LogicalVolume(solidGEM,
				   Copper,
				   "logicGEM");

  logicGEM->SetVisAttributes(GEMVisAttributes);
  

  for(G4int i =-12;i<13;i++){

    //GEMs on positive side of z axis

    for(G4int j=0;j<3;j++){
      
      for(G4int k=-1;k<2;k++){
	
	fPhysicGEMsPlus = new G4PVPlacement(0,
					    G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace),CathodeSize_z/2+GEMDistanceFromCathode+j*GEMGap+GEMOuterSize_z/2),
					    logicGEM,
					    "GEM_"+std::to_string((i+13)*100+j*10+k+1),
					    logicWorld,
					    false,
					    (i+13)*100+j*10+k+1
					    );
	
	fListGEMsOuter.push_back("GEM_"+std::to_string((i+13)*100+j*10+k+1));
	MassMap["GEMsOuter"]+=logicGEM->GetMass(); 
	//std::cout << "North GEM: " << (i+13)*100+j*10+k+1 << "\n";

      }//chiudo for su k
      
    }//chiudo for su j


    
    //GEMs on negative side of z axis
    
    for(G4int j=0;j>-3;j--){
      
      for(G4int k=-1;k<2;k++){
	
	fPhysicGEMsMinus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace),-1*(CathodeSize_z/2+GEMDistanceFromCathode)+j*GEMGap-GEMOuterSize_z/2),
					     logicGEM,
					     "GEM_"+std::to_string((i+13)*100+j*10+k+4 ),
					     logicWorld,
					     false,
					     (i+13)*100+j*10+k+4
					     );
	
	fListGEMsOuter.push_back("GEM_"+std::to_string((i+13)*100+j*10+k+4 ));
	MassMap["GEMsOuter"]+=logicGEM->GetMass();  
	//std::cout << "South GEM: " << (i+13)*100+j*10+k+4 << "\n"; 
      }//chiudo for su j

    }//chiudo for su k

  }//chiudo for su i




  
  //
  //GEMs core
  //

  
  
  G4Colour PMMAColor(1,1,1,0.0);
  G4VisAttributes* PMMAVisAttributes = new G4VisAttributes(PMMAColor);
  PMMAVisAttributes->SetForceSolid(true);
  
  G4double GEMCoreSize_x = 50.0*cm;
  G4double GEMCoreSize_y = 80.0*cm;
  G4double GEMCoreSize_z = 0.029*mm;

  fGEMCoreWidth=GEMCoreSize_z/2;
  
  G4Box*
    solidCoreGEM = new G4Box("GEMInnner",
				GEMCoreSize_x/2,GEMCoreSize_y/2,GEMCoreSize_z/2);

  G4LogicalVolume*
    logicCoreGEM = new G4LogicalVolume(solidCoreGEM,
				   PMMA,
				   "GEMInnerLogical");
  
  logicCoreGEM->SetVisAttributes(PMMAVisAttributes);

  for(G4int i =-12;i<13;i++){

    //GEMs on positive side of z axis

    for(G4int j=0;j<3;j++){
      
      for(G4int k=-1;k<2;k++){
	
	fPhysicGEMsCorePlus = new G4PVPlacement(0,
					    G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace),CathodeSize_z/2+GEMDistanceFromCathode+j*GEMGap+GEMOuterSize_z/2),
					    logicCoreGEM,
					    "GEMCore_"+std::to_string((i+13)*100+j*10+k+1),
					    logicWorld,
					    false,
					    (i+13)*100+j*10+k+1
					    );
	
	fListGEMsCore.push_back("GEMCore_"+std::to_string((i+13)*100+j*10+k+1));
	MassMap["GEMsCore"]+=logicCoreGEM->GetMass(); 
	//std::cout << "North GEM: " << (i+13)*100+j*10+k+1 << "\n";
	
      }//chiudo for su k
      
    }//chiudo for su j


    
    //GEMs on negative side of z axis
    
    for(G4int j=0;j>-3;j--){
      
      for(G4int k=-1;k<2;k++){
	
	fPhysicGEMsCoreMinus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace),-1*(CathodeSize_z/2+GEMDistanceFromCathode)+j*GEMGap-GEMOuterSize_z/2),
					     logicCoreGEM,
					     "GEMCore_"+std::to_string((i+13)*100+j*10+k+4 ),
					     logicWorld,
					     false,
					     (i+13)*100+j*10+k+4
					     );
	
	fListGEMsCore.push_back("GEMCore_"+std::to_string((i+13)*100+j*10+k+4 ));
	MassMap["GEMsCore"]+=logicGEM->GetMass();  
	//std::cout << "South GEM: " << (i+13)*100+j*10+k+4 << "\n"; 
	
	}//chiudo for su j

    }//chiudo for su k

  }//chiudo for su i



  
  //
  //Old field ring
  //

  /*
  
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
			       GEMOuterSize_x/2,GEMOuterSize_y/2,Ring_z/2);

  G4Box*
    innerShapeRing = new G4Box("RingOuterShape",
			 GEMOuterSize_x/2-Ring_width/2,GEMOuterSize_y/2-Ring_width/2,Ring_z/2+0.3*cm);
			 
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
	MassMap["Rings"]+=logicRing->GetMass();  
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
	MassMap["Rings"]+=logicRing->GetMass(); 
	//std::cout << "South Rings: " << (i+13)*1000+(j+10)*10+k+4 << "\n"; 

      }//chiudo for k
    }//chiudo for j
  }

  */


  //
  //RingSupport
  //

  G4Colour SupportRingColor(1.0,1.0,0.0,0.4);
  G4VisAttributes* SupportRingsVisAttributes = new G4VisAttributes(SupportRingColor);
  SupportRingsVisAttributes->SetForceSolid(true);
  
  G4double RingSupportWidth = 0.075*mm;
  G4double RingStripWidth = 0.035*mm;

  fRingSupportWidth=RingSupportWidth;
  
  G4Box*
    outerRingSupport = new G4Box("outerRingSupport",
				 (GEMOuterSize_x/2+RingSupportWidth+RingStripWidth),(GEMOuterSize_y/2+RingSupportWidth+RingStripWidth),GEMDistanceFromCathode/2);

  G4Box*
    innerRingSupport = new G4Box("innerRingSupport",
				 (GEMOuterSize_x/2+RingStripWidth),(GEMOuterSize_y/2+RingStripWidth),GEMDistanceFromCathode/2+0.5*cm);

  

  G4VSolid* solidRingSupport =
    new G4SubtractionSolid("solidRingSupport",outerRingSupport,innerRingSupport );
  
  G4LogicalVolume*
    logicRingSupport = new G4LogicalVolume(solidRingSupport,
					   PMMA,
					   "RingSupport");
  
  logicRingSupport->SetVisAttributes(SupportRingsVisAttributes);
  
  for(G4int i=-12;i<13;i++){  
    for(G4int j=-1;j<2;j++){
      
      fPhysicRingsSupportPlus = new G4PVPlacement(0,
						  G4ThreeVector(i*(CathodeSize_x+detectorSpace)-RingStripWidth,j*(CathodeSize_y+detectorSpace)-RingStripWidth/2,CathodeSize_z/2+GEMDistanceFromCathode/2),
						  logicRingSupport,
						  "RingSupport_"+std::to_string((j+2)*100+(i+12)),
						  logicWorld,
						  false,
						  (j+2)*100+(i+12)
						  );
      
      fListSupportRings.push_back("RingSupport_"+std::to_string((j+2)*100+(i+12)));
      MassMap["SupportRings"]+=logicRingSupport->GetMass();
            
    }//chiudo for j
  }//chiudo for i

  for(G4int i=-12;i<13;i++){  
    for(G4int j=-1;j<2;j++){
      
      fPhysicRingsSupportMinus = new G4PVPlacement(0,
						   G4ThreeVector(i*(CathodeSize_x+detectorSpace)-RingStripWidth,j*(CathodeSize_y+detectorSpace)-RingStripWidth/2,-1*(CathodeSize_z/2+GEMDistanceFromCathode/2)),
						  logicRingSupport,
						  "RingSupport_"+std::to_string((j+6)*100+(i+12)),
						  logicWorld,
						  false,
						  (j+6)*100+(i+12)
						  );
      
      fListSupportRings.push_back("RingSupport_"+std::to_string((j+6)*100+(i+12)));
      MassMap["SupportRings"]+=logicRingSupport->GetMass();
            
    }//chiudo for j
  }//chiudo for i



  


  

  //
  //RingStrip
  //

  G4Colour StripRingColor(0.5,0.5,0.5,0.4);
  G4VisAttributes* StripRingsVisAttributes = new G4VisAttributes(StripRingColor);
  StripRingsVisAttributes->SetForceSolid(true);

  G4double Ring_z = 5.5*cm; //depth of the ring in Z
  G4int NRings = 4;
  
  G4double RingSpacing = (GEMDistanceFromCathode-NRings*Ring_z)/(NRings+1);

  fRingStripWidth=RingStripWidth;
  
  G4Box*
    outerRingStrip = new G4Box("outerRingStrip",
				 (GEMOuterSize_x/2+RingStripWidth),(GEMOuterSize_y/2+RingStripWidth),Ring_z/2);

  G4Box*
    innerRingStrip = new G4Box("innerRingStrip",
				 (GEMOuterSize_x)/2,(GEMOuterSize_y)/2,Ring_z/2+0.5*cm);


  G4VSolid* solidRingStrip =
    new G4SubtractionSolid("solidRingStrip",outerRingStrip,innerRingStrip,0,G4ThreeVector(-RingStripWidth/2,-RingStripWidth/2,0) );
  
  G4LogicalVolume*
    logicRingStrip = new G4LogicalVolume(solidRingStrip,
					   Copper,
					   "RingStrip");
  
  logicRingStrip->SetVisAttributes(StripRingsVisAttributes);
  
  for(G4int i=-12;i<13;i++){
    
    for(G4int j=0;j<NRings;j++){
      
      //Rings on positive side of z axis
      for(G4int k=-1;k<2;k++){
	
	fPhysicRingStripsPlus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace)-0.26*RingStripWidth,k*(CathodeSize_y+detectorSpace),(j+1)*RingSpacing+0.5*Ring_z+j*Ring_z),
					     logicRingStrip,
					     "RingStrip_"+std::to_string((i+13)*1000+(j+10)*10+k+1),
					     logicWorld,
					     false,
					     (i+13)*1000+(j+10)*10+k+1
					     );
	
	fListRingStrips.push_back("RingStrip_"+std::to_string( (i+13)*1000+(j+10)*10+k+1) );
	MassMap["RingStrips"]+=logicRingStrip->GetMass();  
	
	
      }//chiudo for k 
    }//chiudo for j

    for(G4int j=0;j<NRings;j++){
      
      for(G4int k=-1;k<2;k++){ 

	fPhysicRingStripsMinus = new G4PVPlacement(0,
					      G4ThreeVector(i*(CathodeSize_x+detectorSpace)-0.26*RingStripWidth,k*(CathodeSize_y+detectorSpace),-1*((j+1)*RingSpacing+0.5*Ring_z+j*Ring_z)),
					      logicRingStrip,
					      "RingStrip_"+std::to_string((i+13)*1000+(j+10)*10+k+4),
					      logicWorld,
					      false,
					      (i+13)*1000+(j+10)*10+k+4
					      );
      
	fListRingStrips.push_back("RingStrip_"+std::to_string( (i+13)*1000+(j+10)*10+k+4) );
	MassMap["RingStrips"]+=logicRingStrip->GetMass(); 
	//std::cout << "South Rings: " << (i+13)*1000+(j+10)*10+k+4 << "\n"; 

      }//chiudo for k
    }//chiudo for j

    
  }//chiudo for i

 

  //
  //SMD Resistors
  //
  
  G4double ResistorSize_x = 1.6*mm;
  G4double ResistorSize_y = 0.55*mm;
  G4double ResistorSize_z = 3.2*mm;

  fResistorWidth=ResistorSize_y;
  
  G4Colour ResistorColor(0.0,0.0,0.0,0.3);
  G4VisAttributes* ResistorVisAttributes = new G4VisAttributes(ResistorColor);
  ResistorVisAttributes->SetForceSolid(true);
  
  G4Box*
    solidResistor = new G4Box("resistorShape",
			      ResistorSize_x/2,ResistorSize_y/2,ResistorSize_z/2);

  G4LogicalVolume*
    logicResistor = new G4LogicalVolume(solidResistor,
					 al2o3,
					 "Resistor");
  

  for(G4int i=-12;i<13;i++){
    
    for(G4int j=0;j<NRings+1;j++){
      
      //Rings on positive side of z axis
      for(G4int k=-1;k<2;k++){
	
	fPhysicResistorsPlus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace)+CathodeSize_y/2+RingSupportWidth+RingStripWidth+ResistorSize_y/2,(j)*RingSpacing+0.5*Ring_z+j*Ring_z),
					     logicResistor,
					     "Resistor_"+std::to_string((i+13)*1000+(j+10)*10+k+1),
					     logicWorld,
					     false,
					     (i+13)*1000+(j+10)*10+k+1
					     );
	
	fListResistors.push_back("Resistor_"+std::to_string( (i+13)*1000+(j+10)*10+k+1) );
	MassMap["Resistors"]+=logicResistor->GetMass();  
	
	
      }//chiudo for k 
    }//chiudo for j

    for(G4int j=0;j<NRings+1;j++){
      
      for(G4int k=-1;k<2;k++){ 

	fPhysicResistorsMinus = new G4PVPlacement(0,
					      G4ThreeVector(i*(CathodeSize_x+detectorSpace),k*(CathodeSize_y+detectorSpace)+CathodeSize_y/2+RingSupportWidth+RingStripWidth+ResistorSize_y/2,-1*((j)*RingSpacing+0.5*Ring_z+j*Ring_z)),
					      logicResistor,
					      "Resistor_"+std::to_string((i+13)*1000+(j+10)*10+k+4),
					      logicWorld,
					      false,
					      (i+13)*1000+(j+10)*10+k+4
					      );
      
	fListResistors.push_back("Resistor_"+std::to_string( (i+13)*1000+(j+10)*10+k+4) );
	MassMap["Resistors"]+=logicRingStrip->GetMass(); 
	//std::cout << "South Rings: " << (i+13)*1000+(j+10)*10+k+4 << "\n"; 

      }//chiudo for k
    }//chiudo for j

    
  }//chiudo for i








  
  
  
  //
  //PMMA Vessel
  //

  G4double Vesselwidth = 0.5*cm;

  G4double VesselSize_x_outer = (CathodeSize_x/2 + 12*(+CathodeSize_x + detectorSpace) + 2*Vesselwidth);
  G4double VesselSize_y_outer = (CathodeSize_y/2 + (+CathodeSize_y + detectorSpace) + 2*Vesselwidth);
  G4double VesselSize_z_outer = (CathodeSize_z/2+GEMDistanceFromCathode+2*GEMGap+3*GEMOuterSize_z+ 2*Vesselwidth);

  //G4cout << G4BestUnit(2*VesselSize_x_outer,"Length") << "\t" << G4BestUnit(2*VesselSize_y_outer,"Length") << "\t" <<G4BestUnit(2*VesselSize_z_outer,"Length") << G4endl;
  
  G4double VesselSize_x_inner = VesselSize_x_outer-Vesselwidth;
  G4double VesselSize_y_inner = VesselSize_y_outer-Vesselwidth;
  G4double VesselSize_z_inner = VesselSize_z_outer-Vesselwidth;

  //G4cout << G4BestUnit(2*VesselSize_x_inner,"Length") << "\t" << G4BestUnit(2*VesselSize_y_inner,"Length") << "\t" <<G4BestUnit(2*VesselSize_z_inner,"Length") << G4endl;
  
  
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
				      Copper,
				      "Vessel");
  
  logicVessel->SetVisAttributes(PMMAVisAttributes);
  
  
  fPhysicVessel = new G4PVPlacement(0,
				    G4ThreeVector(0,0,0),
				    logicVessel,
				    "Vessel",
				    logicWorld,
				    true,
				    0
				    );

  MassMap["Vessel"]+=logicVessel->GetMass();
  
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
					    G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace)+k*VerticalLensSpacing/2, GEMDistanceFromCathode+2*GEMGap+3*GEMOuterSize_z + LensDistanceFromGEMs ),
					    logicLens,
					    "Lens_"+std::to_string((i+13)*1000+(j+1)*10+(k+0.5)),
					    logicWorld,
					    false,
					    (i+13)*1000+(j+1)*10+(k+0.5)
					    );

	fListLens.push_back("Lens_"+std::to_string((i+13)*1000+(j+1)*10+(k+0.5)));
	MassMap["Lens"]+=logicLens->GetMass(); 
      }//chiudo for k
      
    }//chiudo for j
    
    for(G4int j=-1;j<2;j++){
      
      for(G4float k=-0.5;k<1;k++){
	
	fPhysicLensMinus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace)+k*VerticalLensSpacing/2, -1*(GEMDistanceFromCathode+2*GEMGap+3*GEMOuterSize_z + LensDistanceFromGEMs) ),
					    logicLens,
					    "Lens_"+std::to_string((i+13)*1000+(j+1+4)*10+(k+0.5)),
					    logicWorld,
					    false,
					    (i+13)*1000+(j+1)*10+(k+0.5)
					    );

	fListLens.push_back("Lens_"+std::to_string((i+13)*1000+(j+1+4)*10+(k+0.5)));
	MassMap["Lens"]+=logicLens->GetMass(); 
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
					    G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace)+k*VerticalLensSpacing/2, GEMDistanceFromCathode+2*GEMGap+3*GEMOuterSize_z + LensDistanceFromGEMs + SensorDistanceFromLens ),
					    logicSensor,
					    "Sensor_"+std::to_string((i+13)*1000+(j+1)*10+(k+0.5)),
					    logicWorld,
					    false,
					    (i+13)*1000+(j+1)*10+(k+0.5)
					    );

	fListSensors.push_back("Sensor_"+std::to_string((i+13)*1000+(j+1)*10+(k+0.5)));
	MassMap["Sensors"]+=logicSensor->GetMass(); 
      }//chiudo for k
      
    }//chiudo for j
    
    for(G4int j=-1;j<2;j++){
      
      for(G4float k=-0.5;k<1;k++){
	
	fPhysicSensorsMinus = new G4PVPlacement(0,
					     G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace)+k*VerticalLensSpacing/2, -1*(GEMDistanceFromCathode+2*GEMGap+3*GEMOuterSize_z + LensDistanceFromGEMs + SensorDistanceFromLens) ),
					    logicSensor,
					    "Sensor_"+std::to_string((i+13)*1000+(j+1+4)*10+(k+0.5)),
					    logicWorld,
					    false,
					    (i+13)*1000+(j+1)*10+(k+0.5)
					    );

	fListSensors.push_back("Sensor_"+std::to_string((i+13)*1000+(j+1+4)*10+(k+0.5)));
	MassMap["Sensors"]+=logicSensor->GetMass();
      }//chiudo for k
      
    }//chiudo for j

    
  }//chiudo for i


  //
  //CF4 sensitive volume
  //


  G4Colour GasColor(0.0,0.0,1.0,0.2);
  G4VisAttributes* GasVisAttributes = new G4VisAttributes(GasColor);
  GasVisAttributes->SetForceSolid(true);
  
  G4Box* solidGasVolume = new G4Box("GasVolume",CathodeSize_x/2,CathodeSize_y/2,(GEMDistanceFromCathode)/2);

  fLogicalGasVolume = new G4LogicalVolume(solidGasVolume,
					  CYGNO_gas,
					  "GasVolume"
					  );

  fLogicalGasVolume->SetVisAttributes(GasVisAttributes);

  G4int counter = 0;
  
  for(G4int i=-12;i<13;i++){  
    for(G4int j=-1;j<2;j++){
      
      G4VPhysicalVolume* GasVolumePlus = new G4PVPlacement(0,
							   G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace),CathodeSize_z/2+GEMDistanceFromCathode/2),
							   fLogicalGasVolume,
							   "GasVolume_"+std::to_string(counter),
							   logicWorld,
							   false,
							   counter
							   );

      fListDetector.push_back("GasVolume_"+std::to_string(counter));

      G4cout  << "detN " << counter << "\t coord " << i*(CathodeSize_x+detectorSpace)  << " " <<j*(CathodeSize_y+detectorSpace) << " " << CathodeSize_z/2+GEMDistanceFromCathode/2<< "\n";

      counter++;

      
    }//chiudo for j
  }//chiudo for i
  
  for(G4int i=-12;i<13;i++){ 
    for(G4int j=-1;j<2;j++){
      
      G4VPhysicalVolume* GasVolumeMinus = new G4PVPlacement(0,
							    G4ThreeVector(i*(CathodeSize_x+detectorSpace),j*(CathodeSize_y+detectorSpace),-(CathodeSize_z/2+GEMDistanceFromCathode/2)),
							   fLogicalGasVolume,
							   "GasVolume_"+std::to_string(counter),
							   logicWorld,
							   false,
							   counter
							   );
      
      fListDetector.push_back("GasVolume_"+std::to_string(counter));

      G4cout  << "detN " << counter << "\t coord " << i*(CathodeSize_x+detectorSpace)  << " " <<j*(CathodeSize_y+detectorSpace) << " " << -(CathodeSize_z/2+GEMDistanceFromCathode/2)<< "\n";
      
      counter++;

    }//chiudo for j
  }//chiudo secondo for i
  
  //
  //Creating the physicalvolumestore
  //
  
  fPhysVolStore = G4PhysicalVolumeStore::GetInstance();
  
  for (auto& el: MassMap) {
    std::cout << "Mass of: " << el.first << " = " << G4BestUnit(el.second,"Mass") << "\n";
  }
  //
  //always return the physical World
  //  
  
  return physiWorld;
}


void DetectorConstruction::ConstructSDandField()
{

  fSensitiveDetector = new SensitiveDetector("SensitiveDetector");  
  fLogicalGasVolume->SetSensitiveDetector(fSensitiveDetector);
  
  
}


//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
