// TDetector.cc
// Base Class for detector description

#include "TDetector.hh"
#include "DetectorConstruction.hh"

TDetector::TDetector() : world(0)
{
  materialMgr = MaterialManager::GetPointer();
}

TDetector::TDetector( const DetectorConstruction *aworld)
  : world(aworld)
{
  materialMgr = MaterialManager::GetPointer();
}

