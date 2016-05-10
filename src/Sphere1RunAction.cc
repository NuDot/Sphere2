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
/// \file Sphere1RunAction.cc
/// \brief Implementation of the Sphere1RunAction class

#include "Sphere1RunAction.hh"
#include "Sphere1PrimaryGeneratorAction.hh"
#include "Sphere1EventAction.hh"
#include "Sphere1SteppingAction.hh"
  // use of other actions 
  // - primary generator: to get info for printing about the primary
  // - event action: to get and reset accumulated energy sums
  // - stepping action: to get info about accounting volume 

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UnitsTable.hh"
#include "G4SystemOfUnits.hh"

#include "Randomize.hh"
#include <time.h>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1RunAction::Sphere1RunAction()
: G4UserRunAction()
{
  // add new units for dose
  //

  /* 
  const G4double milligray = 1.e-3*gray;
  const G4double microgray = 1.e-6*gray;
  const G4double nanogray  = 1.e-9*gray;  
  const G4double picogray  = 1.e-12*gray;
   
  new G4UnitDefinition("milligray", "milliGy" , "Dose", milligray);
  new G4UnitDefinition("microgray", "microGy" , "Dose", microgray);
  new G4UnitDefinition("nanogray" , "nanoGy"  , "Dose", nanogray);
  new G4UnitDefinition("picogray" , "picoGy"  , "Dose", picogray);        
  */

}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1RunAction::~Sphere1RunAction()
{}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1RunAction::BeginOfRunAction(const G4Run* aRun)
{ 

//AE: random seed added for generating large samples
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine); //selection of a random engine
  G4long seed=time(0); //returns time in seconds as an integer
  CLHEP::HepRandom::setTheSeed(seed);//changes the seed of the random engine
  CLHEP::HepRandom::showEngineStatus();//shows the actual seed

  G4cout << "### Run " << aRun->GetRunID() << " start." << G4endl;

  //inform the runManager to save random number seed
  //?
  G4RunManager::GetRunManager()->SetRandomNumberStore(false);
    
  //initialize event cumulative quantities
  Sphere1EventAction::Instance()->Reset();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1RunAction::EndOfRunAction(const G4Run* aRun)
{
  G4int nofEvents = aRun->GetNumberOfEvent();
  if (nofEvents == 0) return;
  
  // Compute dose
  //
  G4double energySum  = Sphere1EventAction::Instance()->GetEnergySum();
  G4double energy2Sum = Sphere1EventAction::Instance()->GetEnergy2Sum();
  G4double rms = energy2Sum - energySum*energySum/nofEvents;
  if (rms > 0.) rms = std::sqrt(rms); else rms = 0.;

  /*
  G4double mass = Sphere1SteppingAction::Instance()->GetVolume()->GetMass();
  G4double dose = energySum/mass;
  G4double rmsDose = rms/mass;
  */

  // Run conditions
  //
/*
  const G4ParticleGun* particleGun 
    = Sphere1PrimaryGeneratorAction::Instance()->GetParticleGun();
  G4String particleName 
    = particleGun->GetParticleDefinition()->GetParticleName();                       
  G4double particleEnergy = particleGun->GetParticleEnergy();
        
  // Print
  //  
  G4cout
     << "\n--------------------End of Run------------------------------\n"
     << " The run consists of " << nofEvents << " "<< particleName << " of "
     <<   G4BestUnit(particleEnergy,"Energy")      
     << "\n Dose in scoring volume " 
     << Sphere1SteppingAction::Instance()->GetVolume()->GetName() << " : " 
     << G4BestUnit(energySum,"Energy") 
     << " +- "                   << rms
     << "\n------------------------------------------------------------\n"
     << G4endl;
*/
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
