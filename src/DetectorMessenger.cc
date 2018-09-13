// DetectorMessenger.cc
// 2012.11.02 Yasuhiro Togano

#include "DetectorMessenger.hh"
#include "DetectorConstruction.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithoutParameter.hh"
#include "G4UIcmdWith3VectorAndUnit.hh"

// Constructor =============================================================
DetectorMessenger::DetectorMessenger(DetectorConstruction *det)
  : G4UImessenger(), fDetectorConstruction(det)
{
  fDirectory = new G4UIdirectory("/Gamma/");
  fDirectory->SetGuidance("UI commands for Gamma");

  fDetDirectory = new G4UIdirectory("/Gamma/det/");
  fDetDirectory->SetGuidance("Detector Construction control");
  
  fCrystalMaterialCmd = new G4UIcmdWithAString("/Gamma/det/setCryMaterial",
					       this);
  fCrystalMaterialCmd->SetGuidance("Set Material of crystal. BGO or CsI");
  fCrystalMaterialCmd->SetParameterName("choice",false);
  fCrystalMaterialCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fCrystalInRadiusCmd = new G4UIcmdWithADoubleAndUnit("/Gamma/det/setInRad", 
						      this);
  fCrystalInRadiusCmd->SetGuidance("Set inner radius sphere array");
  fCrystalInRadiusCmd->SetParameterName("radius",false);
  fCrystalInRadiusCmd->SetRange("radius>=0.");
  fCrystalInRadiusCmd->SetUnitCategory("Length");
  fCrystalInRadiusCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fArrayPosCmd = new G4UIcmdWith3VectorAndUnit("/Gamma/det/setArrayPos",this);
  fArrayPosCmd->SetGuidance("set global offset of center position of array");
  fArrayPosCmd->SetParameterName("posX","posY","posZ",false,false);
  fArrayPosCmd->SetUnitCategory("Length");
  fArrayPosCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fLaBrPosCmd = new G4UIcmdWith3VectorAndUnit("/Gamma/det/setLaBrPos",this);
  fLaBrPosCmd->SetGuidance("set global offset of center position of LaBr");
  fLaBrPosCmd->SetParameterName("posX","posY","posZ",false,false);
  fLaBrPosCmd->SetUnitCategory("Length");
  fLaBrPosCmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fUpdateCmd = new G4UIcmdWithoutParameter("/Gamma/det/update",this);
  fUpdateCmd->SetGuidance("Update crystal geometory.");
  fUpdateCmd->SetGuidance("This command MUST be applied before \"beamOn\" ");
  fUpdateCmd->SetGuidance("if you changed geometrical value(s).");
  fUpdateCmd->AvailableForStates(G4State_Idle);

}

// Destructor ===============================================================
DetectorMessenger::~DetectorMessenger()
{
  delete fDirectory;
  delete fDetDirectory;
  delete fCrystalMaterialCmd;
  delete fCrystalInRadiusCmd;
  delete fArrayPosCmd;
  delete fLaBrPosCmd;
  delete fUpdateCmd;
}

void DetectorMessenger::SetNewValue(G4UIcommand *command, G4String newValue)
{
  if(command == fCrystalMaterialCmd){
    fDetectorConstruction->SetDetMaterialName(newValue);
  }
  if(command == fCrystalInRadiusCmd){
    fDetectorConstruction->
      SetInRadius(fCrystalInRadiusCmd->GetNewDoubleValue(newValue));
  }
  if(command == fArrayPosCmd){
    fDetectorConstruction->
      SetArrayPos(fArrayPosCmd->GetNew3VectorValue(newValue));
  }
  if(command == fLaBrPosCmd){
    fDetectorConstruction->
      SetLaBrPos(fLaBrPosCmd->GetNew3VectorValue(newValue));
  }
  if(command == fUpdateCmd){
    fDetectorConstruction->UpdateGeometry();
  }

}

