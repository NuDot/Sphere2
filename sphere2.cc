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
/// \file sphere1.cc
/// \brief Main program of the Sphere1 application

//?
//#include "QGSP_BERT_HP.hh" //***
//#include "G4OpticalPhysics.hh" //***
//#include "G4OpticalProcessIndex.hh" //***
//#include "G4VModularPhysicsList.hh" //***
#include "Sphere1DetectorConstruction.hh"
#include "Sphere1PrimaryGeneratorAction.hh"
#include "Sphere1RunAction.hh"
#include "Sphere1EventAction.hh"
#include "Sphere1TrackingAction.hh"
#include "Sphere1SteppingAction.hh"
#include "Sphere1StackingAction.hh"
#include "Sphere1PhysicsList.hh" 

#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "QBBC.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

#ifdef G4UI_USE
#include "G4UIExecutive.hh"
#endif

#include "Randomize.hh"

#include "EventStructure.hh"

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include <iostream>
#include <string>

event Ev; 

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

int main(int argc,char** argv)
{
  std::string mystr;
  std::cout << "Enter the filename:";
  getline (std::cin, mystr);
  const char * filename = mystr.c_str();
  TROOT root("", "");
  TFile* f;
  //f = TFile::Open("sph_out_bkgC10_rndVtx_3p0mSphere_1.root", "recreate");
//  f = TFile::Open("sph_out_promptC10_2p529MeV_center_5.root", "recreate");
//  f = TFile::Open("sph_out_1gamma_0p718MeV_center_1k_test3.root", "recreate");
  f = TFile::Open(filename, "recreate");
//    f = TFile::Open("sph_out_topology180_center_NoMultScat_100.root", "recreate");

  TTree* epgTree = new TTree("epgTree", "epgTree");
  epgTree->Branch("evt_num",&Ev.evt_num,"evt_num/I");
  epgTree->Branch("C10dT",&Ev.C10dT,"C10dT/F");
  epgTree->Branch("B10dT",&Ev.B10dT,"B10dT/F");
  epgTree->Branch("t0",&Ev.t0,"t0/D");

  epgTree->Branch("isB10_718",&Ev.isB10_718,"isB10_718/I");
  epgTree->Branch("isB10_1740",&Ev.isB10_1740,"isB10_1740/I");
  epgTree->Branch("C10_0_ti",&Ev.C10_0_ti,"C10_0_ti/D");
  epgTree->Branch("nue_ti",&Ev.nue_ti,"nue_ti/D");
  epgTree->Branch("B10_0_ti",&Ev.B10_0_ti,"B10_0_ti/D");
  epgTree->Branch("B10_718_ti",&Ev.B10_718_ti,"B10_718_ti/D");
  epgTree->Branch("B10_1740_ti",&Ev.B10_1740_ti,"B10_1740_ti/D");

  epgTree->Branch("N_epg_i",&Ev.N_epg_i,"N_epg_i/I");
  epgTree->Branch("epg_pIDi",&Ev.epg_pIDi,"epg_pIDi[N_epg_i]/I");
  epgTree->Branch("epg_tIDi",&Ev.epg_tIDi,"epg_tIDi[N_epg_i]/I");
  epgTree->Branch("epg_qi",&Ev.epg_qi,"epg_qi[N_epg_i]/F");
  epgTree->Branch("epg_ti",&Ev.epg_ti,"epg_ti[N_epg_i]/D");
  epgTree->Branch("epg_xi",&Ev.epg_xi,"epg_xi[N_epg_i]/F");
  epgTree->Branch("epg_yi",&Ev.epg_yi,"epg_yi[N_epg_i]/F");
  epgTree->Branch("epg_zi",&Ev.epg_zi,"epg_zi[N_epg_i]/F");
  epgTree->Branch("epg_pxi",&Ev.epg_pxi,"epg_pxi[N_epg_i]/F");
  epgTree->Branch("epg_pyi",&Ev.epg_pyi,"epg_pyi[N_epg_i]/F");
  epgTree->Branch("epg_pzi",&Ev.epg_pzi,"epg_pzi[N_epg_i]/F");
  epgTree->Branch("epg_ei",&Ev.epg_ei,"epg_ei[N_epg_i]/F");
  epgTree->Branch("epg_ti",&Ev.epg_ti,"epg_ti[N_epg_i]/D");
  epgTree->Branch("epg_isOPi",&Ev.epg_isOPi,"epg_isOPi[N_epg_i]/I");
  epgTree->Branch("epg_isChei",&Ev.epg_isChei,"epg_isChei[N_epg_i]/I");

  epgTree->Branch("N_epg_f",&Ev.N_epg_f,"N_epg_f/I");
  epgTree->Branch("epg_pIDf",&Ev.epg_pIDf,"epg_pIDf[N_epg_f]/I");
  epgTree->Branch("epg_tIDf",&Ev.epg_tIDf,"epg_tIDf[N_epg_f]/I");
  epgTree->Branch("epg_qf",&Ev.epg_qf,"epg_qf[N_epg_f]/F");
  epgTree->Branch("epg_tf",&Ev.epg_tf,"epg_tf[N_epg_f]/D");
  epgTree->Branch("epg_ltf",&Ev.epg_ltf,"epg_ltf[N_epg_f]/D");
  epgTree->Branch("epg_xf",&Ev.epg_xf,"epg_xf[N_epg_f]/F");
  epgTree->Branch("epg_yf",&Ev.epg_yf,"epg_yf[N_epg_f]/F");
  epgTree->Branch("epg_zf",&Ev.epg_zf,"epg_zf[N_epg_f]/F");
  epgTree->Branch("epg_pxf",&Ev.epg_pxf,"epg_pxf[N_epg_f]/F");
  epgTree->Branch("epg_pyf",&Ev.epg_pyf,"epg_pyf[N_epg_f]/F");
  epgTree->Branch("epg_pzf",&Ev.epg_pzf,"epg_pzf[N_epg_f]/F");
  epgTree->Branch("epg_ef",&Ev.epg_ef,"epg_ef[N_epg_f]/F");
  epgTree->Branch("epg_lf",&Ev.epg_lf,"epg_lf[N_epg_f]/F");
  epgTree->Branch("epg_isOPf",&Ev.epg_isOPf,"epg_isOPf[N_epg_f]/I");
  epgTree->Branch("epg_isChef",&Ev.epg_isChef,"epg_isChef[N_epg_f]/I");

  epgTree->Branch("edep",&Ev.edep,"edep/F");
  epgTree->Branch("edep_cor",&Ev.edep_cor,"edep_cor/F");

  epgTree->Branch("N_phot",&Ev.N_phot,"N_phot/I");
  epgTree->Branch("t_start",&Ev.t_start,"t_start[N_phot]/F");
  epgTree->Branch("x_start",&Ev.x_start,"x_start[N_phot]/F");
  epgTree->Branch("y_start",&Ev.y_start,"y_start[N_phot]/F");
  epgTree->Branch("z_start",&Ev.z_start,"z_start[N_phot]/F");
  epgTree->Branch("x_hit",&Ev.x_hit,"x_hit[N_phot]/F");
  epgTree->Branch("y_hit",&Ev.y_hit,"y_hit[N_phot]/F");
  epgTree->Branch("z_hit",&Ev.z_hit,"z_hit[N_phot]/F");
  epgTree->Branch("x_vtx",&Ev.x_vtx,"x_vtx[N_phot]/F");
  epgTree->Branch("y_vtx",&Ev.y_vtx,"y_vtx[N_phot]/F");
  epgTree->Branch("z_vtx",&Ev.z_vtx,"z_vtx[N_phot]/F");
  epgTree->Branch("cos_theta",&Ev.cos_theta,"cos_theta[N_phot]/F"); 
  epgTree->Branch("photon_wavelength",&Ev.photon_wavelength,"photon_wavelength[N_phot]/F");
  epgTree->Branch("true_time",&Ev.true_time,"true_time[N_phot]/F");
  epgTree->Branch("PE_creation",&Ev.PE_creation,"PE_creation[N_phot]/I");
  epgTree->Branch("PE_time",&Ev.PE_time,"PE_time[N_phot]/F");
  epgTree->Branch("detector_coverage_included",&Ev.detector_coverage_included,"detector_coverage_included[N_phot]/I");
  epgTree->Branch("true_time_corrected",&Ev.true_time_corrected,"true_time_corrected[N_phot]/F");
  epgTree->Branch("PE_time_corrected",&Ev.PE_time_corrected,"PE_time_corrected[N_phot]/F");
//  epgTree->Branch("cos_theta_reco",&Ev.cos_theta_reco,"cos_theta_reco[N_phot]/F");
//  epgTree->Branch("theta_reco",&Ev.,"[N_phot]/F");
//  epgTree->Branch("",&Ev.,"[N_phot]/F");
//  epgTree->Branch("",&Ev.,"[N_phot]/F");
  epgTree->Branch("process",&Ev.process,"process[N_phot]/I");
  
  epgTree->Branch("trueVtxX",&Ev.trueVtxX,"trueVtxX/D");
  epgTree->Branch("trueVtxY",&Ev.trueVtxY,"trueVtxY/D");
  epgTree->Branch("trueVtxZ",&Ev.trueVtxZ,"trueVtxZ/D");
  epgTree->Branch("trurDirX",&Ev.trueDirX,"trueDirX/D");
  epgTree->Branch("trurDirY",&Ev.trueDirY,"trueDirY/D");
  epgTree->Branch("trurDirZ",&Ev.trueDirZ,"trueDirZ/D");

  epgTree->Branch("gs1_r",&Ev.gs1_r,"gs1_r/F");
  epgTree->Branch("gs1_e",&Ev.gs1_e,"gs1_e/F");
  epgTree->Branch("gs2_r",&Ev.gs2_r,"gs2_r/F");
  epgTree->Branch("gs2_e",&Ev.gs2_e,"gs2_e/F");
  epgTree->Branch("gs3_r",&Ev.gs3_r,"gs3_r/F");
  epgTree->Branch("gs3_e",&Ev.gs3_e,"gs3_e/F");
  epgTree->Branch("gs4_r",&Ev.gs4_r,"gs4_r/F");
  epgTree->Branch("gs4_e",&Ev.gs4_e,"gs4_e/F");

  epgTree->Branch("primary_track_length",&Ev.primary_track_length,"primary_track_length/F");
  epgTree->Branch("primary_track_time",&Ev.primary_track_time,"primary_track_time/F");

  // Choose the Random engine
  //
  CLHEP::HepRandom::setTheEngine(new CLHEP::RanecuEngine);
  long seed=CLHEP::HepRandom::getTheSeed();
  G4cout << "The seed of the RNG is: " << seed << " " << G4endl;
  //? seed
  CLHEP::HepRandom::setTheSeed(58); //18April: 58
  G4cout << "The seed of the RNG is: " << seed << " " << G4endl;
  
  // Construct the default run manager
  //
  G4RunManager * runManager = new G4RunManager;
  G4cout << "runManager allocated!" << G4endl;

  // Set mandatory initialization classes
  //
  // Detector construction
  runManager->SetUserInitialization(new Sphere1DetectorConstruction());


  
  // Physics list: switch here by (un)commenting
  // 1) QBBC Physics list
  // ?
  //G4VModularPhysicsList* physicsList = new QBBC;
  //G4VModularPhysicsList* physicsList = new QGSP_BERT_HP; //or FTFP_BERT, //*** 
  //G4OpticalPhysics* opticalPhysics = new G4OpticalPhysics();  //***
  //RegisterPhysics( opticalPhysics ); //***
  // or 2) Physics List copied from ExN06
  G4VUserPhysicsList* physicsList = new Sphere1PhysicsList; 
  G4cout << "physicsList allocated!" << G4endl;


  physicsList->SetVerboseLevel(1);
  runManager->SetUserInitialization(physicsList);
  G4cout << "physicsList set as user initialization!" << G4endl;
    
  // Primary generator action
  runManager->SetUserAction(new Sphere1PrimaryGeneratorAction(&Ev));
  G4cout << "PrimaryGeneratorAction set!" << G4endl;

  // Set user action classes
  //
  // Stepping action
  runManager->SetUserAction(new Sphere1SteppingAction(&Ev));     
  
  runManager->SetUserAction(new Sphere1TrackingAction(&Ev));

  // Event action
  runManager->SetUserAction(new Sphere1EventAction(epgTree, &Ev));

  // Run action
  runManager->SetUserAction(new Sphere1RunAction());
    
  // Stacking action
  runManager->SetUserAction(new Sphere1StackingAction());
 
  G4cout << "All user actions set!" << G4endl;

  // Initialize G4 kernel
  //
  runManager->Initialize();
  G4cout << "runManager Initialized!" << G4endl;
  
#ifdef G4VIS_USE
  // Initialize visualization
  G4VisManager* visManager = new G4VisExecutive;
  // G4VisExecutive can take a verbosity argument - see /vis/verbose guidance.
  // G4VisManager* visManager = new G4VisExecutive("Quiet");
  visManager->Initialize();
#endif

  // Get the pointer to the User Interface manager
  G4UImanager* UImanager = G4UImanager::GetUIpointer();

  if (argc!=1) {
    // batch mode
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    UImanager->ApplyCommand(command+fileName);
  }
  else {
    // interactive mode : define UI session
#ifdef G4UI_USE
    G4UIExecutive* ui = new G4UIExecutive(argc, argv,"csh");
#ifdef G4VIS_USE
    UImanager->ApplyCommand("/control/execute init_vis.mac"); 
#else
    UImanager->ApplyCommand("/control/execute init.mac"); 
#endif
    ui->SessionStart();
    delete ui;
#endif
  }

  // Job termination
  // Free the store: user actions, physics_list and detector_description are
  // owned and deleted by the run manager, so they should not be deleted 
  // in the main() program !

  f->cd();
  epgTree->Write();
  f->Close();
  
#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;

  return 0;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo.....
