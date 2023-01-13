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
/// \file EventAction.cc
/// \brief Implementation of the EventAction class
//
// 
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "EventAction.hh"
#include "HistoManager.hh"
#include "Run.hh"
#include "PrimaryGeneratorAction.hh" 

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4GenericMessenger.hh"

#include <iomanip>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::EventAction(PrimaryGeneratorAction* primary)
:G4UserEventAction(),
 fDecayChain(),
 fEvisTot(0.),
 fPrimary(primary)
{
  // Set default print level 
  G4RunManager::GetRunManager()->SetPrintProgress(10000);

  fEventMessenger = new G4GenericMessenger(this, "/detector/","Element of the detector");
  fEventMessenger->DeclareProperty("RadElement", fElement, "Select radioactive Element");
  //fElement = "Vessel";
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

EventAction::~EventAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::BeginOfEventAction(const G4Event*)
{
 fDecayChain = " ";
 fEvisTot = 0.;
 
 fPrimary->GetParticleGun()->SetParticlePosition(fPrimary->GetPointOnDetectorElement( fElement ));

 //fPrimary->GetParticleGun()->GeneratePrimaryVertex(); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void EventAction::EndOfEventAction(const G4Event* evt)
{
 G4int evtNb = evt->GetEventID(); 
 G4int printProgress = G4RunManager::GetRunManager()->GetPrintProgress();
 //printing survey
 //

 if(evtNb%1000==0){
   G4cout << evtNb << G4endl;
 }
 
 if (evtNb%printProgress == 0) 
   G4cout << "    End of event. Decay chain:" << fDecayChain 
          << G4endl << G4endl;
 
 //total visible energy
 /*G4AnalysisManager::Instance()->FillH1(9, fEvisTot);
 Run* run 
  = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
 run->EvisEvent(fEvisTot);
 */
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......


