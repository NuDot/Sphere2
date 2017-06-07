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
//
// $Id$
//
// 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Sphere1TrackingActionMessenger.hh"

#include "Sphere1TrackingAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1TrackingActionMessenger::Sphere1TrackingActionMessenger(Sphere1TrackingAction* tTrack)
:tTrackingAction(tTrack)
{
  Sphere1Dir = new G4UIdirectory("/Sphere1/"); 
  Sphere1Dir->SetGuidance("UI commands of this example"); 
  
  trackDir = new G4UIdirectory("/Sphere1/track/"); 
  trackDir->SetGuidance("TrackingAction control"); 

  killCmd = new G4UIcmdWithAString("/Sphere1/track/killIon", this); 
  killCmd->SetGuidance("set ion to be killed in Bi-214 decay chain"); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1TrackingActionMessenger::~Sphere1TrackingActionMessenger()
{
  delete killCmd;
  delete trackDir; 
  delete Sphere1Dir; 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1TrackingActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{       
  if( command == killCmd ){ 
    {tTrackingAction->SetKillIon(newValue);}
  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
