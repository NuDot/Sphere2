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
// $Id$
//
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#include "Sphere1StackingAction.hh"

#include "G4ParticleDefinition.hh"
#include "G4ParticleTypes.hh"
#include "G4Track.hh"
#include "G4ios.hh"
//#include "G4VProcess.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1StackingAction::Sphere1StackingAction()
: gammaCounter(0)
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1StackingAction::~Sphere1StackingAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4ClassificationOfNewTrack
Sphere1StackingAction::ClassifyNewTrack(const G4Track * aTrack)
{
  if(aTrack->GetDefinition() == G4OpticalPhoton::OpticalPhotonDefinition())
  { // particle is optical photon
    
    if(aTrack->GetParentID()>0)
    { // particle is secondary
      gammaCounter++;
      //get wavelength of optical photon
      G4double photon_energy=aTrack->GetDynamicParticle()->GetKineticEnergy();
      G4double photon_wavelength=1.2398*0.001/photon_energy;
      photonWavelengths.push_back(photon_wavelength);

     }
  }
  return fUrgent;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1StackingAction::NewStage()
{
  

  G4cout << "Number of optical photons produced in this event: " << gammaCounter << G4endl;

  /*
  //write wavelengths to file
  //?
  //photonWavelengths_file.open("output/wavelengths.txt",std::ios::app);

  //print a list of wavelengths of the produced optical photons.
  //G4cout << "List of the photon wavelengths: " << G4endl;
  
  for (unsigned int iii=0; iii<photonWavelengths.size(); iii++) 
    {
      //G4cout << photonWavelengths[iii] << G4endl;
      photonWavelengths_file << photonWavelengths[iii] << "\n";
    }
  //G4cout << "--------------------------------" << G4endl;
  

  //write to file
  //photonWavelengths_file.close();
  */

  return;
       
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1StackingAction::PrepareNewEvent()
{ gammaCounter = 0;
  photonWavelengths.clear();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
