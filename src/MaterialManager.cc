// MaterialManager.cc 
// 2012.10.15 Yasuhiro Togano

#include "MaterialManager.hh"
#include "G4SystemOfUnits.hh"

// initialize the material pointer
MaterialManager *MaterialManager::ptrMaterialManager = 0;

// Definition of Constructor
MaterialManager::MaterialManager(){
  
  elementTable = G4Element::GetElementTable();
  materialTable = G4Material::GetMaterialTable();
  DefineMaterials(); // definition of this func. is in icc file!

  // list of elements and Materials

  G4cout << *elementTable <<G4endl;
  G4cout << *materialTable <<G4endl;

}

//destructor
MaterialManager::~MaterialManager(){}

// GetElement

G4Element* MaterialManager::GetElement(const G4String& name) const
{
  
  G4Element *element = 0;

  G4int i;

  for(i=0; elementTable->size();i++){
    G4Element *aElement = (*elementTable)[i];
    G4String aName = aElement->GetName();

    if(aName==name){
      element = aElement;
      break;
    }
  }

  if(!element){
    G4String errorMsg = "Element named <" + name + "> is not defined.";
    G4cout << errorMsg << G4endl;
    exit(1);
  }

  return element;

}

//GetMaterial
G4Material* MaterialManager::GetMaterial(const G4String& name) const
{
  
  G4Material *material = 0;

  G4int i;

  for(i=0; materialTable->size();i++){
    G4Material *aMaterial = (*materialTable)[i];
    G4String aName = aMaterial->GetName();

    if(aName==name){
      material = aMaterial;
      break;
    }
  }



  if(!material){
    G4String errorMsg = "Material named <" + name + "> is not defined.";
    G4cout << errorMsg << G4endl;
    exit(1);
  }

  return material;

}

// definition of "DefineMaterials" function is in this file
#include "Materials.icc"

