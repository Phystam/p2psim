// PhysicsListMessenger based on TestEm5
// 2012.10.29 Yasuhiro Togano

#include "PhysicsListMessenger.hh"
#include "PhysicsList.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"

// Constructor ===============================================================
PhysicsListMessenger::PhysicsListMessenger(PhysicsList *physlist)
  : fPhysicsList(physlist)
{
  fPhysDir = new G4UIdirectory("/Gamma/phys/");
  fPhysDir->SetGuidance("physics list commands");
  
  fListCmd = new G4UIcmdWithAString("/Gamma/phys/addPhysics",this);
  fListCmd->SetGuidance("Add modula physics list.");
  fListCmd->SetParameterName("PList",false);
  fListCmd->AvailableForStates(G4State_PreInit);

  fGammaCutCmd = new G4UIcmdWithADoubleAndUnit("/Gamma/phys/SetGCut",this);//単位付きの変数 Double型
  fGammaCutCmd->SetGuidance("Set gamma cut.");
  fGammaCutCmd->SetParameterName("Gcut",false);
  fGammaCutCmd->SetUnitCategory("Length");
  fGammaCutCmd->SetRange("Gcut>0.0");//下限値
  fGammaCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fElectCutCmd = new G4UIcmdWithADoubleAndUnit("/Gamma/phys/SetECut",this);
  fElectCutCmd->SetGuidance("Set electron cut.");
  fElectCutCmd->SetParameterName("Ecut",false);
  fElectCutCmd->SetUnitCategory("Length");
  fElectCutCmd->SetRange("Ecut>0.0");
  fElectCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fProtoCutCmd = new G4UIcmdWithADoubleAndUnit("/Gamma/phys/SetPCut",this);
  fProtoCutCmd->SetGuidance("Set positron cut.");
  fProtoCutCmd->SetParameterName("Pcut",false);
  fProtoCutCmd->SetUnitCategory("Length");
  fProtoCutCmd->SetRange("Pcut>0.0");
  fProtoCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fAllCutCmd = new G4UIcmdWithADoubleAndUnit("/Gamma/phys/SetCuts",this);
  fAllCutCmd->SetGuidance("Set All cuts.");
  fAllCutCmd->SetParameterName("cut",false);
  fAllCutCmd->SetUnitCategory("Length");
  fAllCutCmd->SetRange("cut>0.0");
  fAllCutCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

}

// Destructor ===============================================================
PhysicsListMessenger::~PhysicsListMessenger()
{
  delete fListCmd;
  delete fGammaCutCmd;
  delete fElectCutCmd;
  delete fProtoCutCmd;
  delete fAllCutCmd;
  delete fPhysDir;
}

// SetNewValue ==============================================================
void PhysicsListMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  if(command == fListCmd){
    fPhysicsList->AddPhysicsList(newValue);
  }

  if(command == fGammaCutCmd){ 
    fPhysicsList->SetCutForGamma(fGammaCutCmd->GetNewDoubleValue(newValue));
  }

  if(command == fElectCutCmd){ 
    fPhysicsList->SetCutForElectron(fElectCutCmd->GetNewDoubleValue(newValue));
  }

  if(command == fProtoCutCmd){
    fPhysicsList->SetCutForPositron(fProtoCutCmd->GetNewDoubleValue(newValue));
  }

  if(command == fAllCutCmd){
    G4double cut = fAllCutCmd->GetNewDoubleValue(newValue);
    fPhysicsList->SetCutForGamma(cut);
    fPhysicsList->SetCutForElectron(cut);
    fPhysicsList->SetCutForPositron(cut);
  }

}
