// TDetector.hh
// Base class for user detector description
// virtual method is used for two function below.
// Build()
// Spec()

#ifndef TDETECTOR_H
#define TDETECTOR_H 1

#include "MaterialManager.hh"

class DetectorConstruction;

class TDetector{
protected:
  MaterialManager *materialMgr;
  const DetectorConstruction *world;

public:
  TDetector();
  TDetector( const DetectorConstruction *aworld);
  virtual ~TDetector() {}

  virtual void Build() = 0;
  virtual void Spec() const = 0;

};
#endif

