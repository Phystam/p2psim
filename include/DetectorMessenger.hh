// DetectorMessenger.hh
// Control the detector properties from interactive mode.
// 2012.11.02 Yasuhiro Togano
#ifndef DETECTORMESSENGER_H
#define DETECTORMESSENGER_H 1

#include "globals.hh"
#include "G4UImessenger.hh"

class DetectorConstruction;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;
class G4UIcmdWithoutParameter;
class G4UIcmdWith3VectorAndUnit;

class DetectorMessenger : public G4UImessenger
{
private:
  DetectorConstruction* fDetectorConstruction;

  G4UIdirectory* fDirectory;
  G4UIdirectory* fDetDirectory;
  G4UIcmdWithAString* fCrystalMaterialCmd;
  G4UIcmdWithADoubleAndUnit* fCrystalInRadiusCmd;
  G4UIcmdWith3VectorAndUnit* fArrayPosCmd;
  G4UIcmdWithoutParameter* fUpdateCmd;
  G4UIcmdWith3VectorAndUnit* fLaBrPosCmd;

public:
  DetectorMessenger(DetectorConstruction* det);
  virtual ~DetectorMessenger();

  virtual void SetNewValue(G4UIcommand*, G4String);
};

#endif


