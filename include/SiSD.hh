// Definition of the CalorimeterSD class
// 2012.10.16 Yasuhiro Togano
// almost the copy of B4cCalorimeterSD.hh from Kobayashi-kun
// Definition of the SiSD class



#ifndef SISD_H
#define SISD_H 1

#include "G4VSensitiveDetector.hh"
#include "SiHit.hh" 

#include <vector>

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"

// Sensitive Detector class 
// Initialize --> Create a hit for each calorimetor. 
// ProcessHits --> physics quantities like time, energy loss are accounted.
//                 This function is called by G4 kernel at each step.

class SiSD : public G4VSensitiveDetector
{
public:
  SiSD(const G4String& name,
		const G4String& hitsCollectionName,
		G4int nbofCells);
  virtual ~SiSD();

  //methods from base class
  virtual void Initialize(G4HCofThisEvent *hitCollection);
  virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory *history);
  virtual void EndOfEvent(G4HCofThisEvent *hitCollection);

private:
  SiHitsCollection *fHitsCollection;
  G4int fNbofCells;
  G4int* HitID;
};

#endif
