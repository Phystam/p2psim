#ifndef LABR3ARRAY_H
#define LABR3ARRAY_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4Box.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4Sphere.hh"
#include "G4Trap.hh"
#include "G4Polyhedra.hh"
#include "G4IntersectionSolid.hh"
#include "G4SubtractionSolid.hh"
#include "G4UnionSolid.hh"
#include "G4UserLimits.hh"

#include "G4GeometryManager.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "G4NistManager.hh"
#include "G4UnitsTable.hh"

#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "Randomize.hh"
#include "G4ios.hh"
#include <iostream>
#include "TMath.h"

#include "TDetector.hh"
#include "LaBrSD.hh"

class DetectorConstruction;

class TLaBr3Array : public TDetector {
public:
  TLaBr3Array(){};
  TLaBr3Array(const DetectorConstruction *aworld);
  ~TLaBr3Array(){};
  
  virtual void Spec() const;
  virtual void Build();
  
  // Get functions;
  G4int GetNumberOfCrystals() const {return fNumCrystal;}
  G4ThreeVector GetArrayPosition() const { return fArrayPosition;}
  
  // Set functions
  void SetArrayPosition(G4ThreeVector vec){fArrayPosition = vec;}
  void SetHousing(G4bool doHousing, G4double front, G4double side,
		  G4bool transparent);
  void SetInsulation(G4bool doInsulation, G4double front, G4double side, 
		     G4bool transparent);
  
private:
  G4int fNumCrystal;
  G4ThreeVector fArrayPosition;
  LaBrSD *fArraySD;
  G4Region *fArrayRegion;
  FILE *fGeomFile;
  G4double fCrystalDiameter;
  G4double fCrystalLength;
  G4double fHousingThicknessFront;
  G4double fHousingThicknessSide;
  G4double fInsulationThicknessFront;
  G4double fInsulationThicknessSide;

  G4bool fHousing;
  G4bool fTransparentHousing;
  G4bool fInsulation;
  G4bool fTransparentInsulation;

  G4LogicalVolume *logicLaBr3House;
  G4LogicalVolume *logicLaBr3Insulation;
  
};

#endif



  
