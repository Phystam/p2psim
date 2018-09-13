// PhysicsListMessenger.hh
// based on the TestEm5
// 2012.10.29 Yasuhiro Togano

#ifndef PHYSICSLISTMESSENGER_H
#define PHYSICSLISTMESSENGER_H 1

#include "G4UImessenger.hh"
#include "globals.hh"

class PhysicsList;
class G4UIdirectory;
class G4UIcmdWithAString;
class G4UIcmdWithADoubleAndUnit;

class PhysicsListMessenger : public G4UImessenger
{
private:
  PhysicsList *fPhysicsList;
  
  G4UIdirectory *fPhysDir;
  G4UIcmdWithAString *fListCmd;
  G4UIcmdWithADoubleAndUnit *fGammaCutCmd;
  G4UIcmdWithADoubleAndUnit *fElectCutCmd;
  G4UIcmdWithADoubleAndUnit *fProtoCutCmd;
  G4UIcmdWithADoubleAndUnit *fAllCutCmd;

public:
  PhysicsListMessenger(PhysicsList*);
  ~PhysicsListMessenger();

  void SetNewValue(G4UIcommand*, G4String);

};

#endif
