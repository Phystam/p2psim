// Definition of the CalorimeterSD class
// 2012.10.16 Yasuhiro Togano
// almost the copy of B4cCalorimeterSD.hh from Kobayashi-kun

#ifndef CALORIMETERSD_H
#define CALORIMETERSD_H 1

#include "G4VSensitiveDetector.hh"
#include "CalorHit.hh" 

#include <vector>

#include "G4Step.hh"
#include "G4HCofThisEvent.hh"

// Sensiive Detector class 
// Initialize --> Create a hit for each calorimetor. 
// ProcessHits --> physics quantities like time, energy loss are accounted.
//                 This function is called by G4 kernel at each step.

class CalorimeterSD : public G4VSensitiveDetector
{
public:
  CalorimeterSD(const G4String& name,
		const G4String& hitsCollectionName,
		G4int nbofCells);
  virtual ~CalorimeterSD();

  //methods from base class
  virtual void Initialize(G4HCofThisEvent *hitCollection);
  virtual G4bool ProcessHits(G4Step *step, G4TouchableHistory *history);
  virtual void EndOfEvent(G4HCofThisEvent *hitCollection);

private:
  CalorHitsCollection *fHitsCollection;
  G4int fNbofCells;
  G4int* HitID;
};

#endif
