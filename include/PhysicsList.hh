// PhysicsList based on the TestEm5
// 2012/10/29 Yasuhiro Togano

#ifndef PHYSICSLIST_H
#define PHYSICSLIST_H 1 

#include "G4VModularPhysicsList.hh"
#include "globals.hh"

class G4VPhysicsConstructor;
class PhysicsListMessenger;


class PhysicsList : public G4VModularPhysicsList
{
public:
  PhysicsList();
  virtual ~PhysicsList();

  void ConstructParticle();
  void AddPhysicsList(const G4String& name);
  
  void ConstructProcess();
  void AddDecay();
  void AddStepMax();
  
  void SetCuts();
  void SetCutForGamma(G4double);
  void SetCutForElectron(G4double);
  void SetCutForPositron(G4double);

private:
  PhysicsListMessenger* fPhysMessenger;
  
  G4String fEMName;
  G4VPhysicsConstructor *fEMPhysicsList;
  
  G4double fCutForGamma;
  G4double fCutForElectron;
  G4double fCutForPositron;

};
#endif
