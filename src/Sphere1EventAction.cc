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
/// \file Sphere1EventAction.cc
/// \brief Implementation of the Sphere1EventAction class

#include "Sphere1EventAction.hh"

#include "Sphere1RunAction.hh"
#include "Sphere1SteppingAction.hh"
  // use of stepping action to get and reset accumulated energy  

#include "G4RunManager.hh"
#include "G4Event.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1EventAction* Sphere1EventAction::fgInstance = 0;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1EventAction* Sphere1EventAction::Instance()
{
// Static acces function via G4RunManager 

  return fgInstance;
}      

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1EventAction::Sphere1EventAction(TTree* fTree, event* fEvt)
: G4UserEventAction(),
  fPrintModulo(1),
  fEnergySum(0.),
  fEnergy2Sum(0.)
{ 
  fgInstance = this;
  
  eTree=fTree;
  eEvt=fEvt;
  //open file for the event wise information:
//AE  events_file.open("output/events.txt",std::ios::app);
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1EventAction::~Sphere1EventAction()
{ 
  fgInstance = 0;
//AE  events_file.close();
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1EventAction::BeginOfEventAction(const G4Event* event)
{
  eEvt->t0=0.0;

  eEvt->N_epg_i=0;
  eEvt->N_epg_f=0;  
  eEvt->isB10_718=0;
  eEvt->isB10_1740=0;
  eEvt->C10_0_ti=0;
  eEvt->nue_ti=0;
  eEvt->B10_718_ti=0;
  eEvt->B10_1740_ti=0;
  eEvt->N_phot=0;
  eEvt->edep=0;
  eEvt->edep_cor=0;

  G4int eventNb = event->GetEventID();
  if (eventNb%fPrintModulo == 0) { 
    G4cout << "\n---> Begin of event: " << eventNb << G4endl;
  }

  copy_of_eventID = event->GetEventID();

  //draw a random vertex reconstruction point around the true vertex (this is then used to "correct" the hit times). This accounts for the fact that in reality the vertex reconstruction is not perfect. It mimics event by event the additional time smearing in the hits due to vertex reconstruction resolution
  
  true_vertex[0]=event->GetPrimaryVertex()->GetX0();
  true_vertex[1]=event->GetPrimaryVertex()->GetY0();
  true_vertex[2]=event->GetPrimaryVertex()->GetZ0();

  SetReconstructedVertex(true_vertex[0],true_vertex[1],true_vertex[2]);

  G4int Nprimaries = event->GetPrimaryVertex()->GetNumberOfParticle();
  for(G4int i=0;i!=Nprimaries;++i)
  {
    G4PrimaryParticle* primary = event->GetPrimaryVertex()->GetPrimary(i); 
    //get the particle direction (DEBUG of HEPevt interface)
    G4cout << "----------------------The primary particle # "<<i+1<<" information: " << G4endl; 
    G4cout << "The PDG code is " << primary->GetPDGcode() << G4endl;
    G4cout << "The mass is " << primary->GetMass() << G4endl;
    G4cout << "The charge is " << primary->GetCharge() << G4endl;
    //G4cout << "The momentum is " << event->GetPrimaryVertex()->GetPrimary()->GetMomentum() << G4endl;
    G4cout << "The momentum is (" << primary->GetPx() << "," << primary->GetPy() << "," <<primary->GetPz() << ")." << G4endl;  
    G4cout << "The polarization is " << primary->GetPolarization() << G4endl;
  }
   // Reset accounted energy in stepping action
  Sphere1SteppingAction::Instance()->Reset();
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1EventAction::EndOfEventAction(const G4Event* event)
{
  eTree->Fill();
  eEvt->evt_num++;
  // accumulate statistics
  G4double energy = Sphere1SteppingAction::Instance()->GetEnergy();
  fEnergySum  += energy;
  fEnergy2Sum += energy*energy;
 

  CollectEventInformation(event);
}  

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1EventAction::Reset()
{
  //reset cumulative quantities
  //
  fEnergySum = 0.;
  fEnergy2Sum = 0.;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
//This function gives a fake reconstructed vertex (just draws randomly a reconstructed vertex around the true vertex)
void Sphere1EventAction::SetReconstructedVertex(G4double true_x,G4double true_y,G4double true_z)
{
  //?
  fake_reco_vertex[0]=true_x + G4RandGauss::shoot(0.,120.);  // hardcoded 120mm vertex resolution in each direction
  fake_reco_vertex[1]=true_y + G4RandGauss::shoot(0.,120.);
  fake_reco_vertex[2]=true_z + G4RandGauss::shoot(0.,120.);  
  G4cout << "The true vertex is (x,y,z)= (" << true_x << "," << true_y << "," << true_z << ")." << G4endl;
  G4cout << "The fake reco vertex is (x,y,z)= (" << fake_reco_vertex[0] <<"," << fake_reco_vertex[1] << "," <<fake_reco_vertex[2] << ")." << G4endl;
}
//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1EventAction::CollectEventInformation(const G4Event* event)
{

  //get number of primaries within the given event
  G4int Nprimaries = event->GetPrimaryVertex()->GetNumberOfParticle();
  G4double pos_x = event->GetPrimaryVertex()->GetX0();
  G4double pos_y = event->GetPrimaryVertex()->GetY0();
  G4double pos_z = event->GetPrimaryVertex()->GetZ0();

  //for now only support up to two primary particles (1: elastic nuebar e- scattering, 2: double beta decay)
  //particle 1:  
  G4double px1=event->GetPrimaryVertex()->GetPrimary(0)->GetPx();
  G4double py1=event->GetPrimaryVertex()->GetPrimary(0)->GetPy();
  G4double pz1=event->GetPrimaryVertex()->GetPrimary(0)->GetPz();
  G4double momentum1=sqrt(px1*px1+py1*py1+pz1*pz1);
  G4double dir_x1=px1/momentum1;
  G4double dir_y1=py1/momentum1;
  G4double dir_z1=pz1/momentum1;
  G4double mass1=event->GetPrimaryVertex()->GetPrimary(0)->GetMass();
  G4double Etot1=sqrt(momentum1*momentum1+mass1*mass1);
  G4double Ekin1=Etot1-mass1;
  G4int pdg1=event->GetPrimaryVertex()->GetPrimary(0)->GetPDGcode();
  //particle 2:
  G4int pdg2;
  G4double Ekin2;
  G4double dir_x2;
  G4double dir_y2;
  G4double dir_z2;
  if(event->GetPrimaryVertex()->GetPrimary(1))
    {		
      G4double px2=event->GetPrimaryVertex()->GetPrimary(1)->GetPx();
      G4double py2=event->GetPrimaryVertex()->GetPrimary(1)->GetPy();
      G4double pz2=event->GetPrimaryVertex()->GetPrimary(1)->GetPz();
      G4double momentum2=sqrt(px2*px2+py2*py2+pz2*pz2);
      dir_x2=px2/momentum2;
      dir_y2=py2/momentum2;
      dir_z2=pz2/momentum2;
      G4double mass2=event->GetPrimaryVertex()->GetPrimary(1)->GetMass();
      G4double Etot2=sqrt(momentum2*momentum2+mass2*mass2);
      Ekin2=Etot2-mass2;
      pdg2=event->GetPrimaryVertex()->GetPrimary(1)->GetPDGcode();
    }
  else
    {
      pdg2=0;
      Ekin2=0;
      dir_x2=0;
      dir_y2=0;
      dir_z2=0;
    }

  //get PEs
  G4int Npe = Sphere1SteppingAction::Instance()->GetPE();

  // get the info from the PEs which pass a hard-coded time cut: PEs in forward/backward direction
  G4int forward_hits = Sphere1SteppingAction::Instance()->GetForwardHits();
  G4int backward_hits = Sphere1SteppingAction::Instance()->GetBackwardHits();

  //G4cout << "Event number " << event->GetEventID() << ": The number of hits detected in forward direction (after the time cut) is " << forward_hits << " and the number of backward events is " << backward_hits << G4endl; 

  //fill event file
//AE  events_file << event->GetEventID() << " " << Nprimaries << " " << pos_x << " " << pos_y << " " << pos_z << " " << fake_reco_vertex[0] << " " << fake_reco_vertex[1] << " " << fake_reco_vertex[2] << " " << pdg1 << " " << Ekin1 << " " << dir_x1 << " " << dir_y1 << " " << dir_z1 << " " << pdg2 << " " << Ekin2 << " " << dir_x2 << " " << dir_y2 << " " << dir_z2 << " " << Npe << " " << forward_hits << " " << backward_hits << "\n";
  
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

