// RunMessenger
// 2012.11.06 Yasuhiro Togano

#include "RunMessenger.hh"
#include "RunAction.hh"
#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4ios.hh"
#include "globals.hh"

// Constructor ==============================================================
RunMessenger::RunMessenger(RunAction* RA)
  : fRunAction(RA)
{
  fRootDir = new G4UIdirectory("/Gamma/root/");
  fRootDir->SetGuidance("root file control");

  fSetRootFileNameCmd = new G4UIcmdWithAString("/Gamma/root/setFileName", 
					       this);
  fSetRootFileNameCmd->SetGuidance("set name of root file");

}

// Destructor ===============================================================
RunMessenger::~RunMessenger(){
  delete fSetRootFileNameCmd;
}

// SetNewValue ==============================================================
void RunMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if(command == fSetRootFileNameCmd){
    fRunAction->SetRootFileName(newValue);
  }

}


