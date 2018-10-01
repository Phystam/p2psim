// RunAction
// 2012.10.16 Yasuhiro Togano

#ifndef RUNACTION_H 
#define RUNACTION_H 1

#include "G4UserRunAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"

#include "G4Run.hh"

#ifdef __OUTPUT_ROOTFILE__
#include "RunMessenger.hh"
#include "TFile.h"
#include "TTree.h"
#include "TVector3.h"
#include "TLorentzVector.h"

#endif
class TVector3;

class RunAction : public G4UserRunAction
{
private:
  G4int RumNumber;
  RunMessenger *fRunMessenger;
public:
  RunAction();
  virtual ~RunAction();

  // virtual methods from base class
  virtual void BeginOfRunAction(const G4Run*);
  virtual void EndOfRunAction(const G4Run*);

  // Get
  G4int GetRunNumber() const {return tRunNumber;}
  
#ifdef __OUTPUT_ROOTFILE__
private:
  TFile *fOutputRootFile;
  TTree *t1;
  G4String fRootFileName;

  // information from the event and run
  G4int tRunNumber, tEventID;

  // Information from the generated events
  G4double tPosGamma[40][3], tDirGamma[40][3];
  G4double tEGamma[40], tEGammaCM[40];
  G4int tNumGamma;

  // Information from Detectors
  G4int tNumHitDet;
  G4double tTotalEdep, tTotalTrackLength;
  G4double tHitPos[3];
  G4double tBeamBeta;
  G4double DetX[402], DetY[402], DetZ[402];
  //G4double Angle[14];
  G4double DetXLaBr[8],DetYLaBr[8],DetZLaBr[8];
  G4double DetXSi[3],DetYSi[3],DetZSi[3];

  
  std::vector<G4int> tDetID;
  std::vector<G4double> tEdep;
  std::vector<G4double> tTrackLength;
  std::vector<G4double> tTime;
  std::vector<G4double> tDetX, tDetY, tDetZ;
  std::vector<G4double> tEdep_sm;

  std::vector<G4double> tHitX, tHitY, tHitZ;
  //  std::vector<G4double> tSourceX, tSourceY, tSourceZ;
  G4double tSourceX, tSourceY, tSourceZ;

  G4int tNumHitLaBr;
  G4double tTotalEdepLaBr, tTotalTrackLengthLaBr;
  std::vector<G4int> tDetIDLaBr;
  std::vector<G4double> tEdepLaBr;
  std::vector<G4double> tTrackLengthLaBr;
  std::vector<G4double> tTimeLaBr;
  std::vector<G4double> tDetXLaBr, tDetYLaBr, tDetZLaBr;
  std::vector<G4double> tEdep_smLaBr;


  G4int tNumHitSi;
  G4double tTotalEdepSi, tTotalTrackLengthSi;
  std::vector<G4int> tDetIDSi;
  std::vector<G4double> tEdepSi;
  std::vector<G4double> tTrackLengthSi;
  std::vector<G4double> tTimeSi;
  std::vector<G4double> tDetXSi, tDetYSi, tDetZSi;
  std::vector<G4double> tEdep_smSi;
  TVector3 tTrackerPos;
  TVector3 tTrackerPos_sm;

public:
  void SetRootFileName(G4String name){fRootFileName = name;}
    
  // Define Branches for tree
  void DefineBranch();

  inline TFile* GetOutputRootFilePtr() const {return fOutputRootFile;}
  inline TTree *GetTTreePtr() const{return t1;}

  void SettRunNumber(G4int *val){tRunNumber = *val;}
  void SettEventID(G4int *val){tEventID = *val;}

  void SettPosGamma(G4ThreeVector val, G4int i){
    tPosGamma[i][0] = val.x();
    tPosGamma[i][1] = val.y();
    tPosGamma[i][2] = val.z();
    //        G4cout << "X: " <<tPosGamma[i][0]<< "  Y: " <<tPosGamma[i][1]
    //    	   << "  Z: " <<tPosGamma[i][2] << G4endl;
    tSourceX=val.x();
    tSourceY=val.y();
    tSourceZ=val.z();

    //    G4cout << tSourceX << G4endl;

  }
  void SettDirGamma(G4ThreeVector val, G4int i){
    tDirGamma[i][0] = val.x();
    tDirGamma[i][1] = val.y();
    tDirGamma[i][2] = val.z();
  }
  void SettEGamma(G4double val, G4int i){
    tEGamma[i] = val;
  }
  void SettEGammaCM(G4double val, G4int i){
    tEGammaCM[i] = val;
  }
  void SettNumGamma(G4int *val){tNumGamma = *val;}

  void SettTotalEdep(G4double *val){tTotalEdep = *val;}
  void SettTotalTrackLength(G4double *val){tTotalTrackLength = *val;}
  void SettHitPos(G4ThreeVector val){
    tHitPos[0] = val.x();
    tHitPos[1] = val.y();
    tHitPos[2] = val.z();
    tHitX.push_back(val.x());
    tHitY.push_back(val.y());
    tHitZ.push_back(val.z());
  }
  void SettBeamBeta(G4double *val){tBeamBeta = *val;}
  void SettNumHitDet(G4int *val){tNumHitDet = *val;}
  void SettEdep(G4double val){tEdep.push_back(val);}
  void SettEdep_sm(G4double val){tEdep_sm.push_back(val);}
  void SettTrackLength(G4double val){tTrackLength.push_back(val);}
  void SettDetID(G4int val){tDetID.push_back(val);}
  void SettDetPosition(G4int DetectorID){
    tDetX.push_back(DetX[DetectorID-1]);
    tDetY.push_back(DetY[DetectorID-1]);
    tDetZ.push_back(DetZ[DetectorID-1]);
  }
  void SettTime(G4double val){tTime.push_back(val);}

  //
  void SettHitX(G4int val){tHitX.push_back(val);}
  //

  //LaBr
  void SettEdepLaBr(G4double val){tEdepLaBr.push_back(val);}
  void SettEdep_smLaBr(G4double val){tEdep_smLaBr.push_back(val);}
  void SettTrackLengthLaBr(G4double val){tTrackLengthLaBr.push_back(val);}
  void SettDetIDLaBr(G4int val){tDetIDLaBr.push_back(val);}
  void SettDetPositionLaBr(G4int DetectorID){
    tDetXLaBr.push_back(DetXLaBr[DetectorID-1]);
    tDetYLaBr.push_back(DetYLaBr[DetectorID-1]);
    tDetZLaBr.push_back(DetZLaBr[DetectorID-1]);
  }
  void SettTimeLaBr(G4double val){tTimeLaBr.push_back(val);}
  void SettNumHitLaBr(G4int *val){tNumHitLaBr = *val;}
  void SettTotalEdepLaBr(G4double *val){tTotalEdepLaBr = *val;}
  void SettTotalTrackLengthLaBr(G4double *val){tTotalTrackLengthLaBr = *val;}

  //Si
  void SettEdepSi(G4double val){tEdepSi.push_back(val);}
  void SettEdep_smSi(G4double val){tEdep_smSi.push_back(val);}
  void SettTrackLengthSi(G4double val){tTrackLengthSi.push_back(val);}
  void SettDetIDSi(G4int val){tDetIDSi.push_back(val);}
  void SettDetPositionSi(G4int DetectorID){
    tDetXSi.push_back(DetXSi[DetectorID-1]);
    tDetYSi.push_back(DetYSi[DetectorID-1]);
    tDetZSi.push_back(DetZSi[DetectorID-1]);
  }
  void SettDetPosXSi(std::vector<G4double> vx){
    tDetXSi=vx;
  }
  void SettDetPosYSi(std::vector<G4double> vy){
    tDetYSi=vy;
  }
  void SettDetPosZSi(std::vector<G4double> vz){
    tDetZSi=vz;
  }
  void SettTimeSi(G4double val){tTimeSi.push_back(val);}
  void SettNumHitSi(G4int *val){tNumHitSi = *val;}
  void SettTotalEdepSi(G4double *val){tTotalEdepSi = *val;}
  void SettTotalTrackLengthSi(G4double *val){tTotalTrackLengthSi = *val;}
  void SettTrackerPos(TVector3 vec){//tTrackerPos.push_back(&vec);}
    tTrackerPos=vec;}
  void SettTrackerPos_sm(TVector3 vec){//tTrackerPos_sm.push_back(&vec);}
    tTrackerPos_sm=vec;}
  void ClearVectors(){
    // tTrackerPos.clear();
    // tTrackerPos_sm.clear();
    tDetID.clear();
    tTrackLength.clear();
    tEdep.clear();
    tEdep_sm.clear();
    tTime.clear();
    tDetX.clear();
    tDetY.clear();
    tDetZ.clear();
    tHitX.clear();
    tHitY.clear();
    tHitZ.clear();
    //    tSourceX.clear();
    //    tSourceY.clear();
    //    tSourceZ.clear();

    tDetIDLaBr.clear();
    tTrackLengthLaBr.clear();
    tEdepLaBr.clear();
    tEdep_smLaBr.clear();
    tTimeLaBr.clear();
    tDetXLaBr.clear();
    tDetYLaBr.clear();
    tDetZLaBr.clear();

    //Si
    tDetIDSi.clear();
    tTrackLengthSi.clear();
    tEdepSi.clear();
    tEdep_smSi.clear();
    tTimeSi.clear();
    tDetXSi.clear();
    tDetYSi.clear();
    tDetZSi.clear();
}

#endif
};

#endif
