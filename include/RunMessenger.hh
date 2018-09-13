// RunMessenger to change the root file name
// 2012.11.06 Yasuhiro Togano

#ifndef RUNMESSENGER_H
#define RUNMESSENGER_H 1

#include "globals.hh"
#include "G4UImessenger.hh"
#include "G4ios.hh"

class RunAction;
class G4UIdirectory;
class G4UIcmdWithAString;

class RunMessenger : public G4UImessenger
{
public:
  RunMessenger(RunAction*);
  ~RunMessenger();

  void SetNewValue(G4UIcommand*, G4String);

private:
  RunAction* fRunAction;
  G4UIdirectory* fRootDir;
  G4UIcmdWithAString * fSetRootFileNameCmd;
  
};
#endif


