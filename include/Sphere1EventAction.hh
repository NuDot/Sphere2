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
/// \file Sphere1EventAction.hh
/// \brief Definition of the Sphere1EventAction class

#ifndef Sphere1EventAction_h
#define Sphere1EventAction_h 1

#include "G4UserEventAction.hh"
#include "globals.hh"

#include <iostream>
#include <fstream>

#include "EventStructure.hh"
#include "TTree.h"

class Sphere1SteppingAction;

/// Event action class
///
/// It holds data member fEnergySum and fEnergy2Sum for accumulating 
/// the event energy deposit its square event by event.
/// These data are then used in the run action to compute the dose.
/// The accumulated energy and energy square sums are reset for each 
/// new run via the Reset() function from the run action.

class Sphere1EventAction : public G4UserEventAction
{
  public:
    Sphere1EventAction(TTree*, event*);
    virtual ~Sphere1EventAction();
    
    // static access method
    static Sphere1EventAction* Instance();

    virtual void BeginOfEventAction(const G4Event* event);
    virtual void EndOfEventAction(const G4Event* event);

    void Reset();

    // get methods
    G4double GetEnergySum() const { return fEnergySum; }
    G4double GetEnergy2Sum() const { return fEnergy2Sum; }

    G4double GetTrueVertex_X() { return true_vertex[0]; }
    G4double GetTrueVertex_Y() { return true_vertex[1]; }
    G4double GetTrueVertex_Z() { return true_vertex[2]; }

    G4double GetFakeRecoVertex_X() { return fake_reco_vertex[0]; }
    G4double GetFakeRecoVertex_Y() { return fake_reco_vertex[1]; }
    G4double GetFakeRecoVertex_Z() { return fake_reco_vertex[2]; }
 
    void SetReconstructedVertex(G4double,G4double,G4double);
    G4int GetCopyOfEventID() {return copy_of_eventID; }
         
  private:
    event* eEvt;
    TTree* eTree;
    static Sphere1EventAction* fgInstance;  

    G4int     fPrintModulo;
    G4double  fEnergySum;
    G4double  fEnergy2Sum;

    G4double true_vertex[3];
    G4double fake_reco_vertex[3];

    std::ofstream events_file;

    void CollectEventInformation(const G4Event* event);
    G4int copy_of_eventID;

};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

    
