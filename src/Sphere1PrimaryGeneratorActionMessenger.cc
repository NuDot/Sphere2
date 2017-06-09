#include "Sphere1PrimaryGeneratorActionMessenger.hh"

#include "Sphere1PrimaryGeneratorAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1PrimaryGeneratorActionMessenger::Sphere1PrimaryGeneratorActionMessenger(Sphere1PrimaryGeneratorAction* gGen)
:gPrimaryGeneratorAction(gGen)
{
  Sphere1Dir = new G4UIdirectory("/Sphere1/");
  Sphere1Dir->SetGuidance("UI commands of this example");

  genDir = new G4UIdirectory("/Sphere1/gen/");
  genDir->SetGuidance("PrimaryGeneratorAction control");

  vtxCmd = new G4UIcmdWithAString("/Sphere1/gen/setTrueVtx", this);
  vtxCmd->SetGuidance("set whether decay events are centered or isotropic");

  //nuCmd = new G4UIcmdWithAString("/Sphere1/gen/setNeutrinos", this); 
  //nuCmd->SetGuidance("set whether double beta decays are neutrinoless or not"); 
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

Sphere1PrimaryGeneratorActionMessenger::~Sphere1PrimaryGeneratorActionMessenger()
{
  delete vtxCmd;
  //delete nuCmd; 
  delete genDir;
  delete Sphere1Dir;
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

void Sphere1PrimaryGeneratorActionMessenger::SetNewValue(G4UIcommand* command,
                                          G4String newValue)
{
  if( command == vtxCmd ){
    {gPrimaryGeneratorAction->SetTrueVtx(newValue);}
  }
  //if( command == nuCmd ){ 
    //{gPrimaryGeneratorAction->SetNeutrinos(newValue);}
//  }
}

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......
