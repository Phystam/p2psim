// DetectorConstruction.hh
// DetectorConstruction class is inherent the class G4VUserDetectorConstruction
// 

#ifndef DETECTOR_CONSTRUCTION_H
#define DETECTOR_CONSTRUCTION_H 1

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "globals.hh"
#include "G4PhysicalVolumeStore.hh"
#include "G4LogicalVolumeStore.hh"
#include "G4SolidStore.hh"

#include "DetectorMessenger.hh"

// detector or material to be used in the DetectorConstruction
#include "TCATANA.hh"
#include "TLaBr3Array.hh"

class G4Box;
class G4LogicalVolume;
class G4VPhysicalVolume;

// class definition
class DetectorConstruction : public G4VUserDetectorConstruction
{
private:
  // for experimental hall
  G4Box *solidExpHall;
  G4LogicalVolume *logicExpHall;
  G4VPhysicalVolume *physiExpHall;
  G4ThreeVector size_ExpHall;

  TCATANA *fCATANA;
  TLaBr3Array *fLaBr3Array;
  
  DetectorMessenger *fDetMessenger;
  
  // buffer for Crystal
  G4double fInRadius;//, fCrystalThickness, fBeamAcceptance;
  G4String fDetMaterialName;
  G4ThreeVector fArrayPos;
  G4ThreeVector fLaBrPos;
public:
  DetectorConstruction();
  ~DetectorConstruction();

  virtual G4VPhysicalVolume *Construct();
  void UpdateGeometry();

  // Gets
  G4VPhysicalVolume *GetExpHall() const {return physiExpHall;};

  // Set functions
  void SetInRadius(G4double value){fInRadius = value;}
  void SetDetMaterialName(G4String name){fDetMaterialName = name;}
  void SetArrayPos(G4ThreeVector vec){fArrayPos = vec;}

  void SetLaBrPos(G4ThreeVector vec){fLaBrPos = vec;}
};

#endif

