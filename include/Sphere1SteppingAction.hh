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
/// \file Sphere1SteppingAction.hh
/// \brief Definition of the Sphere1SteppingAction class

#ifndef Sphere1SteppingAction_h
#define Sphere1SteppingAction_h 1

#include "G4UserSteppingAction.hh"
#include "globals.hh"
#include <iostream>
#include <fstream>

#include "EventStructure.hh"

class G4LogicalVolume;

/// Stepping action class
/// 
/// It holds data member fEnergy for accumulating the energy deposit
/// in a selected volume step by step.
/// The selected volume is set from  the detector construction via the  
/// SetVolume() function. The accumulated energy deposit is reset for each 
/// new event via the Reset() function from the event action.

class Sphere1SteppingAction : public G4UserSteppingAction
{
  public:
    Sphere1SteppingAction(event*);
    virtual ~Sphere1SteppingAction();

    // static access method
    static Sphere1SteppingAction* Instance();

    // method from the base class
    virtual void UserSteppingAction(const G4Step*);

    // reset accumulated hits
    void Reset();


    // apply the transit time spread when detecting photons with a PMT
    G4double ApplyTransitTimeSpread(G4double);      

    //apply a time smearing due to the vertex reco uncertainty
    G4double ApplyVertexRecoSpread(G4double);

    // apply the PMT coverage
    G4bool ApplyCoverage();
 
    //correct the hit time
    G4double CorrectHitTime(G4double, G4double, G4double, G4double); 

    // set methods
    void SetVolume(G4LogicalVolume* volume) { fVolume = volume; }
  
    // get methods
    G4LogicalVolume* GetVolume() const { return fVolume; }
    G4double GetEnergy() const { return fEnergy; }
    G4int GetPE() const { return fPE;}
    G4int GetForwardHits() const { return fForwardHit; }
    G4int GetBackwardHits() const { return fBackwardHit; }
  
  private:
    static Sphere1SteppingAction* fgInstance;  

    event* sEvt;
  
    G4LogicalVolume* fVolume;
    G4double  fEnergy;
    G4int fPE;
    G4int fForwardHit;
    G4int fBackwardHit;
    std::ofstream hit_positions_file;
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif
