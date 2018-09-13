// LaBrHit.cc
// 2012.10.16 Yasuhiro Togano
// Based on the B4cCalorHit.cc by Kobayashi-kun

#include "LaBrHit.hh"
#include "G4UnitsTable.hh"
#include "G4VVisManager.hh"
#include "G4Circle.hh"
#include "G4Colour.hh"
#include "G4VisAttributes.hh"

#include <iomanip>

G4Allocator<LaBrHit> LaBrHitAllocator;


// constructor
LaBrHit::LaBrHit()
  : G4VHit(), fEdep(0.), fTrackLength(0.)
{}

// destructor
LaBrHit::~LaBrHit(){}

LaBrHit::LaBrHit(const LaBrHit& right)
  : G4VHit()
{
  fEdep = right.fEdep;
  fTrackLength = right.fTrackLength;
}

const LaBrHit& LaBrHit::operator=(const LaBrHit& right)
{
  fEdep = right.fEdep;
  fTrackLength = right.fTrackLength;

  return *this;
}


G4int LaBrHit::operator==(const LaBrHit& right) const
{ 
  return ( this ==&right ) ? 1: 0;
}

void LaBrHit::Print(){
  G4cout << "Edep: "
	 << std::setw(7) << G4BestUnit(fEdep,"Energy") << G4endl
	 << "TrackLength: "
	 << std::setw(7) << G4BestUnit(fTrackLength,"Length") << G4endl
	 << "HitPosition: "
	 << G4BestUnit(fHitPosition,"Length") << G4endl;
}



