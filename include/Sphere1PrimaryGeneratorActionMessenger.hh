#ifndef Sphere1PrimaryGeneratorActionMessenger_h
#define Sphere1PrimaryGeneratorActionMessenger_h 1

#include "globals.hh"
#include "G4UImessenger.hh"

class Sphere1PrimaryGeneratorAction;
class G4UIdirectory;
class G4UIcmdWithAString;

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

class Sphere1PrimaryGeneratorActionMessenger: public G4UImessenger
{
  public:
    Sphere1PrimaryGeneratorActionMessenger(Sphere1PrimaryGeneratorAction* );
   ~Sphere1PrimaryGeneratorActionMessenger();

  public:
    void SetNewValue(G4UIcommand * , G4String);

  private:
    Sphere1PrimaryGeneratorAction*     gPrimaryGeneratorAction;

    G4UIdirectory*        Sphere1Dir;
    G4UIdirectory*        genDir;
    G4UIcmdWithAString*   vtxCmd;
    G4UIcmdWithAString*   nuCmd; 
};

//....oooOO0OOooo........oooOO0OOooo........oooOO0OOooo........oooOO0OOooo......

#endif

