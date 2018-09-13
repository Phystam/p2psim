//PGAMessenger.cc (PrimaryGeneratorActionMessenger.cc)

#include "PGAMessenger.hh"

// Constructor-===============================================================
PGAMessenger::PGAMessenger(PrimaryGeneratorAction *pga)
  : G4UImessenger(), fPGA(pga)
{
  fDirectory = new G4UIdirectory("/Gamma/");
  fDirectory->SetGuidance("UI commands for Gamma");

  fPGADirectory = new G4UIdirectory("/Gamma/primary/");
  fPGADirectory->SetGuidance("Initial conditions of the beam and gamma");

  fNumGammaCmd = new G4UIcmdWithAnInteger("/Gamma/primary/setNumGammas",this);
  fNumGammaCmd->SetGuidance("Set the number of Gamma generated in a event");
  fNumGammaCmd->SetParameterName("Ngamma",true);
  fNumGammaCmd->SetDefaultValue(1);
  fNumGammaCmd->SetRange("0<Ngamma<=5");
  fNumGammaCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
  
  fEgamma1Cmd = new G4UIcmdWithADoubleAndUnit("/Gamma/primary/setEgamma1",this);
  fEgamma1Cmd->SetGuidance("Set Egamma1");
  fEgamma1Cmd->SetParameterName("Egamma1",true);
  fEgamma1Cmd->SetDefaultValue(1.0);
  fEgamma1Cmd->SetRange("Egamma1>0");
  fEgamma1Cmd->SetUnitCategory("Energy");
  fEgamma1Cmd->SetUnitCandidates("MeV");
  fEgamma1Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEgamma2Cmd = new G4UIcmdWithADoubleAndUnit("/Gamma/primary/setEgamma2",this);
  fEgamma2Cmd->SetGuidance("Set Egamma2");
  fEgamma2Cmd->SetParameterName("Egamma2",false);
  fEgamma2Cmd->SetDefaultValue(1.0);
  fEgamma2Cmd->SetRange("Egamma2>0");
  fEgamma2Cmd->SetUnitCategory("Energy");
  fEgamma2Cmd->SetUnitCandidates("MeV");
  fEgamma2Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEgamma3Cmd = new G4UIcmdWithADoubleAndUnit("/Gamma/primary/setEgamma3",this);
  fEgamma3Cmd->SetGuidance("Set Egamma3");
  fEgamma3Cmd->SetParameterName("Egamma3",false);
  fEgamma3Cmd->SetDefaultValue(1.0);
  fEgamma3Cmd->SetRange("Egamma3>0");
  fEgamma3Cmd->SetUnitCategory("Energy");
  fEgamma3Cmd->SetUnitCandidates("MeV");
  fEgamma3Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEgamma4Cmd = new G4UIcmdWithADoubleAndUnit("/Gamma/primary/setEgamma4",this);
  fEgamma4Cmd->SetGuidance("Set Egamma4");
  fEgamma4Cmd->SetParameterName("Egamma4",false);
  fEgamma4Cmd->SetDefaultValue(1.0);
  fEgamma4Cmd->SetRange("Egamma4>0");
  fEgamma4Cmd->SetUnitCategory("Energy");
  fEgamma4Cmd->SetUnitCandidates("MeV");
  fEgamma4Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fEgamma5Cmd = new G4UIcmdWithADoubleAndUnit("/Gamma/primary/setEgamma5",this);
  fEgamma5Cmd->SetGuidance("Set Egamma5");
  fEgamma5Cmd->SetParameterName("Egamma5",false);
  fEgamma5Cmd->SetDefaultValue(1.0);
  fEgamma5Cmd->SetRange("Egamma5>0");
  fEgamma5Cmd->SetUnitCategory("Energy");
  fEgamma5Cmd->SetUnitCandidates("MeV");
  fEgamma5Cmd->AvailableForStates(G4State_PreInit, G4State_Idle);

  fBeamBetaCmd = new G4UIcmdWithADouble("/Gamma/primary/setBeamBeta",this);
  fBeamBetaCmd->SetGuidance("Set beta of the beam");
  fBeamBetaCmd->SetParameterName("beta",true);
  fBeamBetaCmd->SetDefaultValue(0.6);
  fBeamBetaCmd->SetRange("0<=beta<=1");
  fBeamBetaCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

// Destructor-===============================================================
PGAMessenger::~PGAMessenger(){
  delete fDirectory;
  delete fPGADirectory;
  delete fNumGammaCmd;
  delete fEgamma1Cmd;
  delete fEgamma2Cmd;
  delete fEgamma3Cmd;
  delete fEgamma4Cmd;
  delete fEgamma5Cmd;
  delete fBeamBetaCmd;
}

//===========================================================================
void PGAMessenger::SetNewValue(G4UIcommand *cmd, G4String newValue){
  if(cmd == fNumGammaCmd){
    fPGA->SetNumGamma(fNumGammaCmd->GetNewIntValue(newValue));
  }
  if(cmd == fEgamma1Cmd){
    fPGA->SetEgamma1(fEgamma1Cmd->GetNewDoubleValue(newValue));
  }
  if(cmd == fEgamma2Cmd){
    fPGA->SetEgamma2(fEgamma2Cmd->GetNewDoubleValue(newValue));
  }
  if(cmd == fEgamma3Cmd){
    fPGA->SetEgamma3(fEgamma3Cmd->GetNewDoubleValue(newValue));
  }
  if(cmd == fEgamma4Cmd){
    fPGA->SetEgamma4(fEgamma4Cmd->GetNewDoubleValue(newValue));
  }
  if(cmd == fEgamma5Cmd){
    fPGA->SetEgamma5(fEgamma5Cmd->GetNewDoubleValue(newValue));
  }
  if(cmd == fBeamBetaCmd){
    fPGA->SetBeamBeta(fBeamBetaCmd->GetNewDoubleValue(newValue));
  }
}  

    
