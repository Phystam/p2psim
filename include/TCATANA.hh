// Class to construct the sphere gamma detector 
// 2012.11.8 Yasuhiro Togano

#ifndef TCATANA_H
#define TCATANA_H 1

#include "TDetector.hh"
#include "globals.hh"
#include "CalorimeterSD.hh"
#include "G4ThreeVector.hh"

class G4Trd;
//class G4Trap;
//class G4SubtractionSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;
class G4VisAttributes;
class DetectorConstruction;

class TCATANA : public TDetector {
public:
  TCATANA(){};
  TCATANA(const DetectorConstruction *aworld);
  ~TCATANA(){};

  virtual void Spec() const;
  virtual void Build();

  // Get
  G4int GetNbOfCrystal() const {return fNumCrystal;}
  G4double GetInnerRadius() const {return fInRadius;}
  G4String GetMaterialName() const {return fMaterialName;}
  G4ThreeVector GetPositionVector() const {return fPosArray;}


  //Sets
  void SetInnerRadius(G4double radius){fInRadius = radius;}
  void SetMaterialName(G4String name){fMaterialName = name;}
  void SetPositionVector(G4ThreeVector vec){fPosArray = vec;}

private:
  CalorimeterSD *fArraySD;
  G4Region *fArrayRegion;

  G4double  fAngCoveragePerCrystal;
  G4ThreeVector fPosArray; // center of SphereArray
  G4double fInRadius; // innner radius
  G4double fCrystalThickness;
  G4int fNumCrystal;
  G4String fMaterialName; 
  G4double fBeamAcceptance;
};
#endif
