// EventAction.cc
// 2012.10.16 Yasuhiro Togano

#include "EventAction.hh"
#include "CalorimeterSD.hh"
#include "CalorHit.hh"
#include "LaBrSD.hh"
#include "LaBrHit.hh"
#include "SiSD.hh"
#include "SiHit.hh"
#include "RunAction.hh"

#include "G4RunManager.hh"
#include "G4Event.hh"
#include "G4SDManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4UnitsTable.hh"

#include "Randomize.hh"
#include <iomanip>
#include <vector>

#include "G4SystemOfUnits.hh"

#define Det_E_Threshold 0.1 // MeV

//#define DEBUG

// Constructor --=================================================
EventAction::EventAction()
  : G4UserEventAction(), fPrintModulo(10000), CalorCollectionID(-1), SiCollectionID(-1),
    LaBrCollectionID(-1)
{
  std::ifstream ifs;
  ifs.open("prm/GetReso_para_20171126.prm");

  for(G4int i=0;i<100;i++){
    for(G4int j=0;j<2;j++){
      ifs >> resopara[i][j];
      }
  }
}

// Destructor --==================================================
EventAction::~EventAction(){}
// ===============================================================
void EventAction::BeginOfEventAction(const G4Event *anEvent)
{

  G4int eventID = anEvent->GetEventID();
  if(eventID%fPrintModulo == 0) {
    G4cout << "\n---> Begin of event: " << eventID << G4endl;
  }
}

// ===============================================================
void EventAction::EndOfEventAction(const G4Event *anEvent)
{
  //G4cout << " ****** Start EndOfEventAction ******" << G4endl;
  G4int ndetCalor = 0;
  G4int ndetSi = 0;
  G4double totalEdep = 0.;
  G4double totalTrackLength = 0.;

  G4SDManager *SDman = G4SDManager::GetSDMpointer();
  G4RunManager *runManager = G4RunManager::GetRunManager();
  RunAction *runAction = (RunAction *) runManager->GetUserRunAction();

  TTree *t1 = runAction->GetTTreePtr();

  G4int eventID = anEvent->GetEventID();

  if(CalorCollectionID == -1){
    CalorCollectionID = SDman -> GetCollectionID("arrayCollection");
  }

  if(SiCollectionID == -1){
    SiCollectionID = SDman -> GetCollectionID("SiCollection");
    //    std::cout<<"SiCollectionID"<<SiCollectionID<<std::endl;
  }

  /*
  if(LaBrCollectionID == -1){
    LaBrCollectionID = SDman->GetCollectionID("LaBrCollection");
  }
  */

  G4HCofThisEvent *hcte = anEvent->GetHCofThisEvent();
  CalorHitsCollection *CalorHC = NULL;
  SiHitsCollection *SiHC = NULL;
  LaBrHitsCollection *LaBrHC = NULL;

  if(hcte){
    CalorHC = (CalorHitsCollection*)(hcte->GetHC(CalorCollectionID));
    SiHC = (SiHitsCollection*)(hcte->GetHC(SiCollectionID));
    LaBrHC = (LaBrHitsCollection*)(hcte->GetHC(LaBrCollectionID));
  }

  std::vector<G4int> detID;
  std::vector<G4double> edep;
  std::vector<G4double> trackLength;
  std::vector<G4ThreeVector> hitpos;
  std::vector<G4double> time;
  G4int numHit=0;

  if(CalorHC){
    ndetCalor = CalorHC -> entries();// ndetCalor==100 : number of CsI detector
  }

  if(SiHC){
    ndetSi = SiHC -> entries();// ndetSi==6 : number of Si detector
    //    std::cout<<ndetSi<<std::endl;
  }

  for(G4int i=0;i<ndetCalor;i++){
    G4double dedep = (*CalorHC)[i]->GetEdep();

    if(dedep>0){
      numHit++;
      detID.push_back(i+1);
      edep.push_back(dedep);
      trackLength.push_back((*CalorHC)[i]->GetTrackLength());
      hitpos.push_back((*CalorHC)[i]->GetHitPosition());
      time.push_back((*CalorHC)[i]->GetGlobalTime());
      totalEdep += dedep;
      totalTrackLength += trackLength[i];
    }
  }

#ifdef DEBUG
  G4cout << "Calor" << G4endl;
  G4cout << "Energy deposit: " << totalEdep << " MeV " << G4endl;
  G4cout << "Track Length: " << totalTrackLength << G4endl;
  G4cout << "Entries: " << ndetCalor << ", numHit: " << numHit << G4endl;
  G4cout << "numHit: " << numHit << G4endl;
  G4cout << "======================================" << G4endl;
#endif

  G4int nLaBr = 0;
  /*
  if(LaBrHC){
    nLaBr = LaBrHC -> entries();
  }
  */

  G4double totalEdep_LaBr = 0.;
  G4double totalTrackLength_LaBr = 0.;

  G4int numHitLaBr=0;
  std::vector<G4int> detID_LaBr;
  std::vector<G4double> edep_LaBr;
  std::vector<G4double> trackLength_LaBr;
  std::vector<G4ThreeVector> hitpos_LaBr;
  std::vector<G4double> time_LaBr;

  for(G4int i=0;i<nLaBr;i++){
    G4double dedep = (*LaBrHC)[i]->GetEdep();
    if(dedep>0){
      numHitLaBr++;
      detID_LaBr.push_back(i+1);
      edep_LaBr.push_back(dedep);
      trackLength_LaBr.push_back((*LaBrHC)[i]->GetTrackLength());
      hitpos_LaBr.push_back((*LaBrHC)[i]->GetHitPosition());
      time_LaBr.push_back((*LaBrHC)[i]->GetGlobalTime());
      totalEdep_LaBr += dedep;
      totalTrackLength_LaBr += trackLength_LaBr[i];
    }
  }

#ifdef DEBUG
  G4cout << "LaBr" << G4endl;
  G4cout << "Energy deposit: " << totalEdep_LaBr << " MeV " << G4endl;
  G4cout << "Track Length: " << totalTrackLength_LaBr << G4endl;
  G4cout << "Entries: " << nLaBr << ", numHit: " << numHit << G4endl;
  G4cout << "numHit: " << numHitLaBr << G4endl;
  G4cout << "======================================" << G4endl;
#endif

  //Si
  std::vector<G4int> detIDSi;
  std::vector<G4double> edepSi;
  std::vector<G4double> trackLengthSi;
  std::vector<G4double> hitposXSi;
  std::vector<G4double> hitposYSi;
  std::vector<G4double> hitposZSi;
  std::vector<TVector3> hitposSi;
  std::vector<G4double> timeSi;
  G4int numHitSi=0;

  for(G4int i=0;i<ndetSi;i++){
    //    G4double dep = (*SiHC)[i]->GetEdep();
    G4double dedepSi = (*SiHC)[i]->GetEdep();
    if(dedepSi>0){
      numHitSi++;
      detIDSi.push_back(i+1);
      edepSi.push_back(dedepSi);
      trackLengthSi.push_back((*SiHC)[i]->GetTrackLength());
      G4ThreeVector hitposg4=(*SiHC)[i]->GetHitPosition();
      //      TVector3 hitposroot(hitposg4.x()/mm,hitposg4.y()/mm,hitposg4.z()/mm);
      TVector3 hitposroot(hitposg4.x()/mm,hitposg4.y()/mm,hitposg4.z()/mm);
      hitposXSi.push_back(hitposroot.X());
      hitposYSi.push_back(hitposroot.Y());
      hitposZSi.push_back(hitposroot.Z());
      time.push_back((*SiHC)[i]->GetGlobalTime());
    }
  }

#ifdef DEBUG
  G4cout << "Si" << G4endl;
  G4cout << "Energy deposit: " << totalEdepSi/MeV << " MeV " << G4endl;
  G4cout << "Track Length: " << totalTrackLength/mm << G4endl;
  G4cout << "Entries: " << ndetSi << ", numHit: " << numHit << G4endl;
  G4cout << "numHit: " << numHitSi << G4endl;
  G4cout << "======================================" << G4endl;
#endif


  //fill Tree
#ifdef __OUTPUT_ROOTFILE__
  G4double val=-1000.;
  //G4ThreeVector tpos, tdir;
  G4int intval=-1000.;

  //tpos.set(0.,0.,0.);

  // EventID
  // BeamPosition and BeamDirection are filled at PrimaryGeneratorAction
  intval = eventID;
  runAction->SettEventID(&intval);
  intval = 0;
  // Detectors

  intval = numHit;
  runAction->SettNumHitDet(&intval);

  val = totalEdep/MeV;
  runAction->SettTotalEdep(&val);

  val = totalTrackLength/cm;
  runAction->SettTotalTrackLength(&val);

  runAction->ClearVectors();
  //CATANA
  for(G4int i=0;i<numHit;i++){
    if(edep[i]<Det_E_Threshold) continue;
    //    G4double sigma = GetDetResolution(edep[i],detID[i]);
    G4double sigma = 0.0456*sqrt(edep[i]);
    
    G4double edep2 = edep[i]>0 ? CLHEP::RandGauss::shoot(edep[i],sigma) : 0;
    runAction->SettEdep(edep[i]);
    runAction->SettEdep_sm(edep2);
    runAction->SettTrackLength(trackLength[i]);
    runAction->SettDetID(detID[i]);
    runAction->SettDetPosition(detID[i]);
    runAction->SettTime(time[i]);

    //
    runAction->SettHitPos(hitpos[i]);
  }

  // LaBr---------------
  intval = numHitLaBr;
  runAction->SettNumHitLaBr(&intval);

  val = totalEdep_LaBr/MeV;
  runAction->SettTotalEdepLaBr(&val);

  val = totalTrackLength_LaBr/cm;
  runAction->SettTotalTrackLengthLaBr(&val);

  for(G4int i=0;i<numHitLaBr;i++){
    if(edep_LaBr[i]<Det_E_Threshold) continue;
    G4double sigma = GetLaBrResolution(edep_LaBr[i]);
    G4double edep2 = CLHEP::RandGauss::shoot(edep_LaBr[i],sigma);
    runAction->SettEdepLaBr(edep_LaBr[i]);
    runAction->SettEdep_smLaBr(edep2);
    runAction->SettTrackLengthLaBr(trackLength_LaBr[i]);
    runAction->SettDetIDLaBr(detID_LaBr[i]);
    runAction->SettDetPositionLaBr(detID_LaBr[i]);
    runAction->SettTimeLaBr(time_LaBr[i]);
  }

  //Si--------------
  intval = numHitSi;
  runAction->SettNumHitSi(&intval);

  // val = totalEdepSi/MeV;
  // runAction->SettTotalEdepSi(&val);

  // val = totalTrackLengthSi/cm;
  // runAction->SettTotalTrackLengthSi(&val);

   for(G4int i=0;i<numHitSi;i++){
    if(edepSi[i]<Det_E_Threshold) continue;
    //    G4double sigma = GetSiResolution(edep_Si[i]);
    runAction->SettDetIDSi(detIDSi[i]);
    runAction->SettDetPosXSi(hitposXSi);
    runAction->SettDetPosYSi(hitposYSi);
    runAction->SettDetPosZSi(hitposZSi);
    
    //    runAction->SettDetPositionSi(detIDSi[i]);
    // runAction->SettTrackerPos(hitposSi[i]);
    // runAction->SettTrackerPos_sm(hitposSi[i]);
  }
 

  t1 -> Fill();
  
#endif

  if(eventID%fPrintModulo ==0){
    G4cout << "###### End of Event: " << eventID << G4endl;
  }
}

// ===============================================================
G4double EventAction::GetDetResolution(G4double edep, G4int ID){
  G4double sigma;

  //  G4double a = 1.5299;// original value by Togano-san
  //  G4double b = 0.5705;
  //  G4double a = 1.48285;// 2017/11/24 
  //  G4double b = 0.5818;
  //  sigma = a/2.35 * pow(edep*1000,b);


  G4double a = resopara[ID-1][0]; 
  G4double b = resopara[ID-1][1];
  sigma = a*pow(edep*1000,b);

  
  return sigma/1000;//keV->MeV
}

// ===============================================================
G4double EventAction::GetLaBrResolution(G4double edep){
  G4double sigma;

  // Parameters are taken from A. Giaz et al, NIMA 729, 910 (2013).
  G4double a = 400.;
  G4double b = 0.625;
  G4double c = 28E-6;
				    
  sigma = sqrt(a + b*edep + c*edep*edep); // given in keV

  return sigma/1000;//keV->MeV
}
