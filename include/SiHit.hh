// CalorHit.hh
// 2012.10.16 Yasuhiro Togano
// originated from the B4cSiHit.hh by Kobayashi-kun

#ifndef SIHIT_H
#define SIHIT_H 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

// Siimetor Hit class
// It defines data members to store the energy deposit and track lengths
// of charged particles in a selected detector volume;
// fEdep, fTrackLength

class SiHit : public G4VHit
{
private:
  G4double fEdep;
  G4double fTrackLength;
  G4double fGlobalTime;
  G4ThreeVector fHitPosition;
public:
  SiHit();
  SiHit(const SiHit&);
  virtual ~SiHit();
  
  // operators
  const SiHit& operator=(const SiHit&);
  G4int operator==(const SiHit&) const;

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
typedef G4THitsCollection<SiHit> SiHitsCollection;
extern G4Allocator<SiHit> SiHitAllocator;
//=========================================================================
inline void *SiHit::operator new(size_t){
  void *hit;
  hit = (void *) SiHitAllocator.MallocSingle();
  return hit;
}

inline void SiHit::operator delete(void *hit){
  SiHitAllocator.FreeSingle((SiHit*) hit);
}

inline void SiHit::Add(G4double de, G4double dl){
  fEdep += de;
  fTrackLength += dl;
}

#endif
