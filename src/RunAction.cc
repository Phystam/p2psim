// RunAction
//
#include "RunAction.hh"
#include "RunMessenger.hh"

#include "G4Run.hh"
#include "G4RunManager.hh"
#include "G4UImanager.hh"
#include "G4VVisManager.hh"
#include "Randomize.hh"
#include "G4UnitsTable.hh"

#include <fstream>
#include <time.h>

// Constructor ============================================================
RunAction::RunAction()
  : G4UserRunAction(), fRootFileName("gamma.root")
{
  fRunMessenger = new RunMessenger(this);
  std::ifstream ifs;
  ifs.open("./geometry/CATANA/CATANAGeometryOUT.txt");
  if(!ifs){ 
    G4cerr << "Error! Cannot open geometry.dat at RunAction" << G4endl;
    G4cout << "Could not find the geometry file" << G4endl;
    return;
  }
  else{
    G4int index = 0;
    while(!ifs.eof()){
      G4double detID;
      ifs >> detID >> DetX[index] >> DetY[index] >> DetZ[index];
      index++;
    }
    ifs.close();
    G4cout << "CATANA Geometry file read in at RunAction" << G4endl;
  }

  std::ifstream ifs2;
  ifs2.open("./geometry/LaBr3/LaBr3GeometryOUT.txt");
  if(!ifs2){ 
    G4cerr << "Error! Cannot open geometryLaBr.dat at RunAction" << G4endl;
    G4cout << "Could not find the LaBrgeometry file" << G4endl;
    return;
  }
  else{
    G4int index = 0;
    while(!ifs2.eof()){
      G4double detID;
      ifs2 >> detID >> DetXLaBr[index] >> DetYLaBr[index] >> DetZLaBr[index];
      index++;
    }
    ifs2.close();
    G4cout << "LaBr Geometry file read in at RunAction" << G4endl;
  }

}

// Destructor =============================================================
RunAction::~RunAction()
{}

// BeginOfRunAction =======================================================
void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  G4int seed;
  seed = (G4int)time(NULL);
  CLHEP::HepRandom::setTheSeed(seed);

  if(G4VVisManager::GetConcreteInstance()){
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/scene/notifyHandlers");
  }

#ifdef __OUTPUT_ROOTFILE__
  fOutputRootFile = new TFile(fRootFileName,"RECREATE");
  t1 = new TTree("t1","CATANA");

  DefineBranch();
  tRunNumber = aRun -> GetRunID();
#endif

  G4cout << "#### Run " << aRun->GetRunID() << " start! ##### " << G4endl;
}
// EndOfRunAction =========================================================
void RunAction::EndOfRunAction(const G4Run *aRun)
{
  G4cout << "##### Run " << aRun->GetRunID() << " end! #####" << G4endl;

  if(G4VVisManager::GetConcreteInstance()){
    G4UImanager::GetUIpointer()->ApplyCommand("/vis/viewerupdate");
  }

#ifdef __OUTPUT_ROOTFILE__
  t1->Print();
  t1->Write();
  fOutputRootFile->Close();
  G4cout << G4endl << "Output file Name: " << fRootFileName << G4endl << G4endl;
#endif
}


#ifdef __OUTPUT_ROOTFILE__
void RunAction::DefineBranch()
{
  // information from Run and Event
  t1->Branch("RunNumber",&tRunNumber,"RunNumber/I");
  t1->Branch("EventID",&tEventID,"EventID/I");
  // information from the generated events
  t1->Branch("EGammaCM",tEGammaCM,"EGammaCM[40]/D");
  t1->Branch("EGamma",tEGamma,"EGamma[40]/D");
  t1->Branch("PosGamma",tPosGamma,"PosGamma[40][3]/D");
  t1->Branch("DirGamma",tDirGamma,"DirGamma[40][3]/D");
  t1->Branch("NumGamma",&tNumGamma,"NumGamma/I");

  // information from the detector
  t1->Branch("NumHit",&tNumHitDet, "NumHit/I");
  t1->Branch("Edep",&tEdep);
  t1->Branch("Edep_sm",&tEdep_sm);
  t1->Branch("DetID",&tDetID);
  t1->Branch("TrackLength",&tTrackLength);
  t1->Branch("DetX",&tDetX);
  t1->Branch("DetY",&tDetY);
  t1->Branch("DetZ",&tDetZ);
  t1->Branch("Time",&tTime);
  t1->Branch("TotalEdep",&tTotalEdep,"TotalEdep/D");
  t1->Branch("TotalTrackLength",&tTotalTrackLength,"TotalTrackLength/D");
  //
  t1->Branch("HitX",&tHitX);
  t1->Branch("HitY",&tHitY);
  t1->Branch("HitZ",&tHitZ);

  t1->Branch("SourceX",&tSourceX);
  t1->Branch("SourceY",&tSourceY);
  t1->Branch("SourceZ",&tSourceZ);


  t1->Branch("TrackerPos","TVector3",&tTrackerPos);
  /*
  t1->Branch("NumHitLaBr",&tNumHitLaBr, "NumHitLaBr/I");
  t1->Branch("EdepLaBr",&tEdepLaBr);
  t1->Branch("Edep_smLaBr",&tEdep_smLaBr);
  t1->Branch("DetIDLaBr",&tDetIDLaBr);
  t1->Branch("TrackLengthLaBr",&tTrackLengthLaBr);
  t1->Branch("DetXLaBr",&tDetXLaBr);
  t1->Branch("DetYLaBr",&tDetYLaBr);
  t1->Branch("DetZLaBr",&tDetZLaBr);
  t1->Branch("TimeLaBr",&tTimeLaBr);
  t1->Branch("TotalEdepLaBr",&tTotalEdepLaBr,"TotalEdepLaBr/D");
  t1->Branch("TotalTrackLengthLaBr",&tTotalTrackLengthLaBr,
	     "TotalTrackLengthLaBr/D");
  */

  t1->Branch("BeamBeta",&tBeamBeta,"BeamBeta/D");
 
}

#endif
