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
/// \file Sphere1DetectorConstruction.cc
/// \brief Implementation of the Sphere1DetectorConstruction class

#include "Sphere1DetectorConstruction.hh"
#include "Sphere1SteppingAction.hh"
   // use of stepping action to set the accounting volume

#include "G4RunManager.hh"
#include "G4NistManager.hh"
#include "G4Box.hh"
#include "G4Cons.hh"
#include "G4Orb.hh"
#include "G4Sphere.hh"
#include "G4Trd.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4SystemOfUnits.hh"
#include "G4LogicalBorderSurface.hh"
#include "G4OpticalSurface.hh"

#include <fstream>
#include <sstream>

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1DetectorConstruction::Sphere1DetectorConstruction()
: G4VUserDetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1DetectorConstruction::~Sphere1DetectorConstruction()
{ }

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

G4VPhysicalVolume* Sphere1DetectorConstruction::Construct()
{ 

  G4double a, z, density, fractionmass;
  G4int nelements;
 
  // Get nist material manager
  G4NistManager* nist = G4NistManager::Instance();
  
  // Sphere parameters
  //
  G4double sphere_inner_radius = 6500.000;//650*cm;
  G4double sphere_outer_radius = 6500.001;//660*cm;
  G4Material* sphere_mat = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
   
  // Option to switch on/off checking of volumes overlaps
  //
  G4bool checkOverlaps = true;

  //     
  // World
  //
  G4double world_sizeX = 2000*cm;
  G4double world_sizeY = 2000*cm;
  G4double world_sizeZ = 2000*cm;

  //define materials
  //first specify the needed elements
  G4Element* C = new G4Element("Carbon", "C", z=6 , a=12.0107*g/mole);
  G4Element* H = new G4Element("Hydrogen", "H", z=1 , a=1.00794*g/mole);
  G4Element* O = new G4Element("Oxygen"  , "O", z=8 , a=15.9994*g/mole);
  G4Element* N = new G4Element("Nitrogen"  , "N", z=7 , a=14.0067*g/mole);
  //now define some molecules or just use water: 
//  G4Material* world_mat = nist->FindOrBuildMaterial("G4_Galactic");//("G4_WATER");
  
  //Toluene  
  //G4Material* world_mat = new G4Material("Toluene", density=0.865*g/cm3, nelements=2);
  //world_mat->AddElement(C,natoms=7);
  //world_mat->AddElement(H,natoms=8);

  //KamLAND scintillator (80% Dodecane, 20% PC, 1.52g/l PPO)
  G4Material* world_mat = new G4Material("KamLAND_scintillator",density=0.7752*g/cm3,nelements=4); //density: 0.8*0.75 + 0.2*0.876
  
  world_mat->AddElement(C,fractionmass=85.8090*perCent);
  world_mat->AddElement(H,fractionmass=14.1644*perCent);
  world_mat->AddElement(N,fractionmass=0.0124*perCent);
  world_mat->AddElement(O,fractionmass=0.0142*perCent);


  //------------ Generate & Add Material Properties Tables ------------
  // read the indices of refraction from file

  //?
  //std::string refraction_file="data/refraction_DC_Target.txt";
  std::string refraction_file="data/KamLAND/IndexOfRefraction_KamLAND_Calculated_smooth.txt";
  //std::string refraction_file="data/KamLAND/refraction_KamLAND_2.txt";
  ReadDataFile(refraction_file);

  G4MaterialPropertiesTable* myMPT1 = new G4MaterialPropertiesTable();
  myMPT1->AddProperty("RINDEX",Ephoton_vector.data(),value_vector.data(),Ephoton_vector.size());

  
  /*
  //debug: 
  const G4int NUMENTRIES = 17;
  G4double ppckov[NUMENTRIES] = {6.27105e-06*MeV,5.99528e-06*MeV,5.88034e-06*MeV,5.716e-06*MeV,5.43736e-06*MeV,5.14118e-06*MeV,4.78062e-06*MeV,4.34087e-06*MeV,4.15e-06*MeV,3.91184e-06*MeV,3.65546e-06*MeV,3.40206e-06*MeV,3.12481e-06*MeV,2.81639e-06*MeV,2.6125e-06*MeV,2.09863e-06*MeV,1.54975e-06*MeV};
  //G4double rindex[NUMENTRIES] = {1.773,1.751,1.72,1.684,1.637,1.586,1.547,1.514,1.501,1.487,1.473,1.468,1.462,1.455,1.448,1.438,1.43};
  G4double rindex[NUMENTRIES] = {1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462,1.462};
  myMPT1->AddProperty("RINDEX",ppckov,rindex,NUMENTRIES);
  */


  ClearDataVectors(); 


  //absorption length data
  //?
  //std::string absorption_file="data/attlength_DC_Target.txt";
  std::string absorption_file="data/KamLAND/attlength_KamLAND.txt";
  //std::string absorption_file="data/KamLAND/attlength_infinite.txt";
  ReadDataFile(absorption_file);
  
  //input data is in meters
  for(unsigned int iii=0;iii<wavelengths_vector.size()-1;iii++)
    {
      value_vector.at(iii)=value_vector.at(iii)*m;
    }

  myMPT1->AddProperty("ABSLENGTH",Ephoton_vector.data(),value_vector.data(),Ephoton_vector.size());
 
  ClearDataVectors();

  //scintillation_waveforms
  //?
  //std::string scintspec_file="data/PPO+Bis-MSB_em.txt";
  std::string scintspec_file="data/KamLAND/scintspec2_KamLAND.txt";
  ReadDataFile(scintspec_file);

  // Here, we have to account for the fact, that a sampling over lambda is not the same as the sampling over energy which is done by Geant4
  // We have our relative scintillation spectrum in wavelength. If we don't modify it, higher energies will get too much weight since our 
  // values are given equidistant in lambda which is not equidistant in energy. We can correct for that as follows: 
  
  //normalize the spectrum to the value at the lowest wavelength
  G4double edist_zero = Ephoton_vector.at(1) - Ephoton_vector.at(0);
  for(unsigned int iii=0;iii<wavelengths_vector.size()-1;iii++)
    {
      G4double energydistance = Ephoton_vector.at(iii+1) - Ephoton_vector.at(iii);
      G4double norm_factor = edist_zero/energydistance;
      value_vector.at(iii)*=norm_factor; 
    } 

  //same scintillation spectrum for the fast and the slow component (the measurement of the spectrum
  // contains both fast and slow components). 
  //?
  myMPT1->AddProperty("FASTCOMPONENT",Ephoton_vector.data(), value_vector.data(), Ephoton_vector.size());
  //?
  myMPT1->AddProperty("SLOWCOMPONENT",Ephoton_vector.data(), value_vector.data(), Ephoton_vector.size());

  ClearDataVectors();

  //Light yield and resolution scale 
  //?
  //myMPT1->AddConstProperty("SCINTILLATIONYIELD",6000./MeV); //DC like
  myMPT1->AddConstProperty("SCINTILLATIONYIELD",9030.5/MeV); //KamLAND
  //myMPT1->AddConstProperty("SCINTILLATIONYIELD",0./MeV); //scintillation switched off
  //?
  myMPT1->AddConstProperty("RESOLUTIONSCALE",1.0);

  //time constants and weights
  //?
  //myMPT1->AddConstProperty("FASTTIMECONSTANT", 2.63*ns); // DC value
  //myMPT1->AddConstProperty("SLOWTIMECONSTANT",9.69*ns); // DC value
  myMPT1->AddConstProperty("FASTTIMECONSTANT", 6.9*ns); // KamLAND: 6.9
  myMPT1->AddConstProperty("SLOWTIMECONSTANT",8.8*ns); // KamLAND: 8.8
  //?
  myMPT1->AddConstProperty("FASTSCINTILLATIONRISETIME",1.0*ns); // 1.0ns typical scale
  myMPT1->AddConstProperty("SLOWSCINTILLATIONRISETIME",1.0*ns); // 1.0ns typical scale
  
//  myMPT1->AddConstProperty("FASTSCINTILLATIONRISETIME",5.0*ns); // 1.0ns typical scale
//  myMPT1->AddConstProperty("SLOWSCINTILLATIONRISETIME",5.0*ns); // 1.0ns typical scale

  //?
  //myMPT1->AddConstProperty("YIELDRATIO",0.737);  //DC value renormalized for only 2 components
  myMPT1->AddConstProperty("YIELDRATIO",0.87);  //KamLAND: 0.87

  world_mat->SetMaterialPropertiesTable(myMPT1);

  //Birks constant kB
  //world_mat->GetIonisation()->SetBirksConstant(0.202*mm/MeV); //not as important now for us, use DC value although maybe the step sizes etc. are different than for DCGLG4sim
  //?
  world_mat->GetIonisation()->SetBirksConstant(0.106*mm/MeV); // KamLAND
 
  //Volumes
  //----------------------------------------------------------------------
//  G4Box* solidWorld =    
//    new G4Box("World",                       //its name
//       0.5*world_sizeX, 0.5*world_sizeY, 0.5*world_sizeZ);     //its size

    G4Sphere* solidWorld = new G4Sphere("SteelSphere",                    //its name
		 0,
		 sphere_outer_radius,
		 0.,2*pi,   //phi range
		 0.,pi);    //theta range
 
      
  G4LogicalVolume* logicWorld =                         
    new G4LogicalVolume(solidWorld,          //its solid
                        world_mat,           //its material
                        "World");            //its name
                                   
  G4VPhysicalVolume* physWorld = 
    new G4PVPlacement(0,                     //no rotation
                      G4ThreeVector(),       //at (0,0,0)
                      logicWorld,            //its logical volume
                      "World",               //its name
                      0,                     //its mother  volume
                      false,                 //no boolean operation
                      0,                     //copy number
                      checkOverlaps);        //overlaps checking
                     
  //     
  // Sphere
  //  
  G4Sphere* solidSphere =    
    new G4Sphere("SteelSphere",                    //its name
		 sphere_inner_radius,
		 sphere_outer_radius,
		 0.,2*pi,   //phi range
		 0.,pi);    //theta range
      
  G4LogicalVolume* logicSphere =                         
    new G4LogicalVolume(solidSphere,            //its solid
                        sphere_mat,             //its material
                        "SteelSphere");         //its name
               
  G4VPhysicalVolume* physicalSphere = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicSphere,                //its logical volume
                    "Sphere",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
 

  //Surfaces: Steel Sphere inner surface
  G4OpticalSurface* SphereInnerSurface = new G4OpticalSurface("SphereSurface"); 
  SphereInnerSurface->SetType(dielectric_metal); 
  SphereInnerSurface->SetFinish(polished);
  SphereInnerSurface->SetModel(unified);

  //G4LogicalBorderSurface* logialSurfaceSphere = 
  new G4LogicalBorderSurface("WaterSurface",physWorld,physicalSphere,SphereInnerSurface);  

  G4MaterialPropertiesTable* myST1 = new G4MaterialPropertiesTable();

  //reflectivity of steel surface: switched off for now
  //?
  const G4int num = 2;
  G4double Ephoton[num] = {1.5*eV, 6.2*eV};
  G4double Reflectivity[num] = {0.,0.};  // Photons are either reflected or absorbed
   
  myST1->AddProperty("REFLECTIVITY", Ephoton, Reflectivity, num);


  //QE for the surface from a file. 
  //?
  std::string QE_file="data/QE_PMT_DC.txt";
  ReadDataFile(QE_file);

  //efficiency: get this from file (QE of metal sphere)
  //old
  //myST1->AddProperty("EFFICIENCY", Ephoton, Efficiency, num);               
  //new
  myST1->AddProperty("EFFICIENCY", Ephoton_vector.data(), value_vector.data(), Ephoton_vector.size()); 

  ClearDataVectors();
  

  SphereInnerSurface->SetMaterialPropertiesTable(myST1);  
/*
  //     
  // InnerSphere
  //  
  G4Sphere* solidInnerSphere =    
    new G4Sphere("InnerSphere",                    //its name
		 0,
		 sphere_inner_radius,
		 0.,2*pi,   //phi range
		 0.,pi);    //theta range
      
  G4LogicalVolume* logicInneSphere =                         
    new G4LogicalVolume(solidInnerSphere,            //its solid
                        inner_sphere_mat,             //its material
                        "InnerSphere");         //its name
               
  G4VPhysicalVolume* physicalInnerSphere = new G4PVPlacement(0,                       //no rotation
                    G4ThreeVector(),         //at (0,0,0)
                    logicInnerSphere,                //its logical volume
                    "InnerSphere",              //its name
                    logicWorld,              //its mother  volume
                    false,                   //no boolean operation
                    0,                       //copy number
                    checkOverlaps);          //overlaps checking
*/ 


  // Set scoring volume to stepping action 
  // (where we will account energy deposit)
  //
  Sphere1SteppingAction* steppingAction = Sphere1SteppingAction::Instance(); 
  ////steppingAction->SetVolume(logicShape1);
  steppingAction->SetVolume(logicSphere);


  //
  //always return the physical World
  //
  return physWorld;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1DetectorConstruction::ReadDataFile(std::string input_file)
{
  std::ifstream ifs;

  ifs.open(input_file.c_str());
   if (!ifs.good()) {
     G4cout << "Could not open " << input_file <<G4endl;
    }

  std::string line;
  while (ifs.good())
    {
      getline( ifs, line );
      if ( ifs.fail() )
        break;
      // put the line in an istringstream for convenient parsing
      std::istringstream lineStream(line);
      G4double value,value2;
      lineStream >> value >> value2;
      wavelengths_vector.push_back(value);
      value_vector.push_back(value2);
      Ephoton_vector.push_back(1239.8*eV/value); 
    }

  G4cout << "For the file " << input_file << " the number of entries which were read is: " <<  wavelengths_vector.size() << G4endl;


  for(unsigned int iii=0;iii<wavelengths_vector.size();iii++)
    {
      G4cout << "iii= " << iii << " " << wavelengths_vector.at(iii) << " "<< Ephoton_vector.at(iii) << " " <<value_vector.at(iii) << G4endl;
    }

 ifs.close();

 return;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1DetectorConstruction::ClearDataVectors()
{
  wavelengths_vector.clear();
  value_vector.clear();
  Ephoton_vector.clear();

  return;
}

