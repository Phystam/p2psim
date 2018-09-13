// CalorHit.hh
// 2012.10.16 Yasuhiro Togano
// originated from the B4cCalorHit.hh by Kobayashi-kun

#ifndef CALORHIT_H
#define CALORHIT_H 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

// Calorimetor Hit class
// It defines data members to store the energy deposit and track lengths
// of charged particles in a selected detector volume;
// fEdep, fTrackLength

class CalorHit : public G4VHit
{
private:
  G4double fEdep;
  G4double fTrackLength;
  G4double fGlobalTime;
  G4ThreeVector fHitPosition;
public:
  CalorHit();
  CalorHit(const CalorHit&);
  virtual ~CalorHit();
  
  // operators
  const CalorHit& operator=(const CalorHit&);
  G4int operator==(const CalorHit&) const;

  inline void* operator new(size_t);
  inline void operator delete(void*);

  // methods from base class
  virtual void Draw(){};
  virtual void Print();

  // methods to handle data
  void Add(G4double de, G4double dl);

  // Gets! 
  inline G4double GetEdep() const {return fEdep;}
  inline G4double GetTrackLength() const {return fTrackLength;}
  inline G4ThreeVector GetHitPosition() const {return fHitPosition;}
  inline G4double GetGlobalTime() const {return fGlobalTime;}

  // Set
  inline void SetHitPosition(G4ThreeVector pos){fHitPosition = pos;}
  inline void SetGlobalTime(G4double time){fGlobalTime = time;}

};
//=========================================================================
typedef G4THitsCollection<CalorHit> CalorHitsCollection;
extern G4Allocator<CalorHit> CalorHitAllocator;
//=========================================================================
inline void *CalorHit::operator new(size_t){
  void *hit;
  hit = (void *) CalorHitAllocator.MallocSingle();
  return hit;
}

inline void CalorHit::operator delete(void *hit){
  CalorHitAllocator.FreeSingle((CalorHit*) hit);
}

inline void CalorHit::Add(G4double de, G4double dl){
  fEdep += de;
  fTrackLength += dl;
}

#endif
