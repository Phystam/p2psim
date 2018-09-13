// PrimaryGeneratorActionMessenger.hh
#ifndef PGAMESSENGER_HH
#define PGAMESSENGER_HH 1

#include "globals.hh"
#include "G4UImessenger.hh"

#include "G4UIdirectory.hh"
#include "G4UIcmdWithAString.hh"
#include "G4UIcmdWithADoubleAndUnit.hh"
#include "G4UIcmdWithADouble.hh"
#include "G4UIcmdWithAnInteger.hh"

#include "PrimaryGeneratorAction.hh"

class PGAMessenger : public G4UImessenger {

private:
  PrimaryGeneratorAction *fPGA;
  G4UIdirectory *fDirectory;
  G4UIdirectory *fPGADirectory;
  G4UIcmdWithAnInteger *fNumGammaCmd;
  G4UIcmdWithADoubleAndUnit *fEgamma1Cmd;
  G4UIcmdWithADoubleAndUnit *fEgamma2Cmd;
  G4UIcmdWithADoubleAndUnit *fEgamma3Cmd;
  G4UIcmdWithADoubleAndUnit *fEgamma4Cmd;
  G4UIcmdWithADoubleAndUnit *fEgamma5Cmd;
  G4UIcmdWithADouble *fBeamBetaCmd;

public:
  PGAMessenger(PrimaryGeneratorAction* pga);
  virtual ~PGAMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);
};

#endif
