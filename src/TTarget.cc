#include "TTarget.hh"
#include "DetectorConstruction.hh"
#include "CalorimeterSD.hh"

#include "G4SDManager.hh"
#include "G4PVPlacement.hh"
#include "G4Region.hh"
#include "globals.hh"
#include "G4UnitsTable.hh"

#include "G4Box.hh"
#include "G4Trd.hh"
#include "G4Trap.hh"
//#include "G4UnionSolid.hh"                                                                                                                         
#include "G4SubtractionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4VisAttributes.hh"
#include "G4ThreeVector.hh"

#include <fstream>
#include "G4SystemOfUnits.hh"


//#define DEBUG                                                                                                                                      

#ifdef DEBUG
#endif
