// CalorHit.hh
// 2012.10.16 Yasuhiro Togano
// originated from the B4cCalorHit.hh by Kobayashi-kun

#ifndef LABRHIT_H
#define LABRHIT_H 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

// Calorimetor Hit class
// It defines data members to store the energy deposit and track lengths
// of charged particles in a selected detector volume;
// fEdep, fTrackLength

class LaBrHit : public G4VHit
{
private:
  G4double fEdep;
  G4double fTrackLength;
  G4double fGlobalTime;
  G4ThreeVector fHitPosition;
public:
  LaBrHit();
  LaBrHit(const LaBrHit&);
  virtual ~LaBrHit();
  
  // operators
  const LaBrHit& operator=(const LaBrHit&);
  G4int operator==(const LaBrHit&) const;

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
typedef G4THitsCollection<LaBrHit> LaBrHitsCollection;
extern G4Allocator<LaBrHit> LaBrHitAllocator;
//=========================================================================
inline void *LaBrHit::operator new(size_t){
  void *hit;
  hit = (void *) LaBrHitAllocator.MallocSingle();
  return hit;
}

inline void LaBrHit::operator delete(void *hit){
  LaBrHitAllocator.FreeSingle((LaBrHit*) hit);
}

inline void LaBrHit::Add(G4double de, G4double dl){
  fEdep += de;
  fTrackLength += dl;
  //G4cout << "Add is okay" << G4endl;
}

#endif
