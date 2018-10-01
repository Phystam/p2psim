// SiSD class 
// 2012.10.16 Yasuhiro Togano
// Almost the copy of B4cSiSD.cc from Kobayashi-kun

#include "SiSD.hh"

#include "G4HCofThisEvent.hh"
#include "G4ThreeVector.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"

//#define DEBUG

// Constructor
SiSD::SiSD(
			     const G4String& name,
			     const G4String& hitsCollectionName,
			     G4int nbofCells)
  : G4VSensitiveDetector(name),
    fHitsCollection(0),
    fNbofCells(nbofCells)
{
  collectionName.insert(hitsCollectionName);
  HitID = new G4int[500];
}

// Destructor
SiSD::~SiSD(){
  delete [] HitID;
}

void SiSD::Initialize(G4HCofThisEvent *hcte)
{
  // Create hits collection
  fHitsCollection = 
    new SiHitsCollection(SensitiveDetectorName, collectionName[0]);

  // Add this collection to hcte
  G4int hcID = 
    G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  hcte->AddHitsCollection(hcID, fHitsCollection);

  // Create hits
  //fNbofCells for cells 
  for (G4int i=0;i<fNbofCells; i++){
    fHitsCollection->insert(new SiHit());
    HitID[i] = -1;
  }
}

// ProcessHits =======================
G4bool SiSD::ProcessHits(G4Step *step, G4TouchableHistory*){
  // energy deposit
  G4double edep = step -> GetTotalEnergyDeposit();

  // step length
  G4double stpl = 0.;
  if(step->GetTrack()->GetDefinition()->GetPDGCharge()!=0.){
    stpl = step->GetStepLength();
  }
  
  if( edep==0. && stpl ==0.) return false;
  
  //Get detector id
  G4int detid = step->GetPreStepPoint()->GetTouchableHandle()->GetCopyNumber();

  SiHit *hit = (*fHitsCollection)[detid];
  if(!hit){
    G4cerr << "Cannot access hit " << G4endl;
    exit(1);
  }

  if(HitID[detid] == -1){
    hit->SetHitPosition(step->GetTrack()->GetPosition());
    hit->SetGlobalTime(step->GetTrack()->GetGlobalTime());
  }

  HitID[detid]++;

#ifdef DEBUG
//  G4cout << "HitID: " << HitID << ", Edep: " << edep 
//	 << ", HitPos: " << step->GetTrack()->GetPosition() << G4endl;
  G4cout << "HitID: " << HitID[detid] << ", detid: " << detid 
	 << ", Edep: " << edep
	 <<", time: " << step->GetTrack()->GetGlobalTime() << G4endl;
#endif
  // add values
  hit->Add(edep, stpl);
  //hit->Print();
  return true;
}

void SiSD::EndOfEvent(G4HCofThisEvent*)
{
  if(verboseLevel>1){
    G4int nofHits = fHitsCollection->entries();
    G4cout << "\n---------->Hits Collection: in this event they are "
	   << nofHits
	   << " hits in the Si :" << G4endl;
    for(G4int i=0;i<nofHits;i++){
      (*fHitsCollection)[i]->Print();
    }
  }
}


