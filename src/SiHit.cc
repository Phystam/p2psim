// SiHit.cc
// 2012.10.16 Yasuhiro Togano
// Based on the B4cSiHit.cc by Kobayashi-kun

#include "SiHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4Allocator<SiHit> SiHitAllocator;


// constructor
SiHit::SiHit()
  : G4VHit(), fEdep(0.), fTrackLength(0.)
{}

// destructor
SiHit::~SiHit(){}

SiHit::SiHit(const SiHit& right)
  : G4VHit()
{
  fEdep = right.fEdep;
  fTrackLength = right.fTrackLength;
}

const SiHit& SiHit::operator=(const SiHit& right)
{
  fEdep = right.fEdep;
  fTrackLength = right.fTrackLength;

  return *this;
}


G4int SiHit::operator==(const SiHit& right) const
{ 
  return ( this ==&right ) ? 1: 0;
}

void SiHit::Print(){
  G4cout << "Edep: "
	 << std::setw(7) << G4BestUnit(fEdep,"Energy") << G4endl
	 << "TrackLength: "
	 << std::setw(7) << G4BestUnit(fTrackLength,"Length") << G4endl
	 << "HitPosition: "
	 << G4BestUnit(fHitPosition,"Length") << G4endl;
}



