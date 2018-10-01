//EventAction
// 2012.10.16 Yasuhiro Togano

#ifndef EVENTACTION_H
#define EVENTACTION_H 1

#include "G4UserEventAction.hh"
#include "CalorHit.hh"
#include "LaBrHit.hh"

#include "globals.hh"

class EventAction : public G4UserEventAction
{
private:
  G4int fPrintModulo;
  G4int CalorCollectionID;
  G4int SiCollectionID;
  G4int LaBrCollectionID;
  G4double resopara[100][2];

public:
  EventAction();
  virtual ~EventAction();
  
  // methods from base class
  virtual void BeginOfEventAction(const G4Event* anEvent);
  virtual void EndOfEventAction(const G4Event* anEvent);

  // set
  void SetPrintModulo(G4int value){fPrintModulo = value;}

  // Get detector resolution
  G4double GetDetResolution(G4double, G4int);
  G4double GetLaBrResolution(G4double);
  
};

#endif
