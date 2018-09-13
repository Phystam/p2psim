// MaterialManager.hh
// 2012.10.15 Yasuhiro Togano

#ifndef MATERIAL_MANAGER_H
#define MATERIAL_MANAGER_H 1

#include "G4Material.hh"
#include "G4Element.hh"

//------------------------
// Class definition
//------------------------

class MaterialManager{
private:
  static MaterialManager* ptrMaterialManager;
  // vectors of elements
  const G4ElementTable *elementTable;
  const G4MaterialTable *materialTable;
  // definition of this func. is in icc file
  void DefineMaterials();

protected:
  MaterialManager();

public:
  ~MaterialManager();

  static MaterialManager *GetPointer();

  G4Element *GetElement( const G4String& name) const;
  G4Material *GetMaterial( const G4String& name) const;

};

// inline function

inline MaterialManager *MaterialManager::GetPointer()
{
  if(!ptrMaterialManager) ptrMaterialManager = new MaterialManager;

  return ptrMaterialManager;
}

#endif

