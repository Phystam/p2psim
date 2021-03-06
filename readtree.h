//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Thu Sep 20 21:14:56 2018 by ROOT version 5.34/30
// from TTree t1/CATANA
// found on file: root/simucat.root
//////////////////////////////////////////////////////////

#ifndef readtree_h
#define readtree_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.
#include <vector>
#include <vector>
#include <TVector3.h>
#include <TObject.h>

// Fixed size dimensions of array or collections stored in the TTree if any.

class readtree {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Int_t           RunNumber;
   Int_t           EventID;
   Double_t        EGammaCM[40];
   Double_t        EGamma[40];
   Double_t        PosGamma[40][3];
   Double_t        DirGamma[40][3];
   Int_t           NumGamma;
   Int_t           NumHit;
   vector<double>  *Edep;
   vector<double>  *Edep_sm;
   vector<int>     *DetID;
   vector<double>  *TrackLength;
   vector<double>  *DetX;
   vector<double>  *DetY;
   vector<double>  *DetZ;
   vector<double>  *Time;
   Double_t        TotalEdep;
   Double_t        TotalTrackLength;
   vector<double>  *HitX;
   vector<double>  *HitY;
   vector<double>  *HitZ;
   Double_t        SourceX;
   Double_t        SourceY;
   Double_t        SourceZ;
 //TVector3        *TrackerPos;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Double_t        fX;
   Double_t        fY;
   Double_t        fZ;
 //TVector3        *TrackerPos_sm;
   UInt_t          fUniqueID;
   UInt_t          fBits;
   Double_t        fX;
   Double_t        fY;
   Double_t        fZ;
   Double_t        BeamBeta;

   // List of branches
   TBranch        *b_RunNumber;   //!
   TBranch        *b_EventID;   //!
   TBranch        *b_EGammaCM;   //!
   TBranch        *b_EGamma;   //!
   TBranch        *b_PosGamma;   //!
   TBranch        *b_DirGamma;   //!
   TBranch        *b_NumGamma;   //!
   TBranch        *b_NumHit;   //!
   TBranch        *b_Edep;   //!
   TBranch        *b_Edep_sm;   //!
   TBranch        *b_DetID;   //!
   TBranch        *b_TrackLength;   //!
   TBranch        *b_DetX;   //!
   TBranch        *b_DetY;   //!
   TBranch        *b_DetZ;   //!
   TBranch        *b_Time;   //!
   TBranch        *b_TotalEdep;   //!
   TBranch        *b_TotalTrackLength;   //!
   TBranch        *b_HitX;   //!
   TBranch        *b_HitY;   //!
   TBranch        *b_HitZ;   //!
   TBranch        *b_SourceX;   //!
   TBranch        *b_SourceY;   //!
   TBranch        *b_SourceZ;   //!
   TBranch        *b_TrackerPos_fUniqueID;   //!
   TBranch        *b_TrackerPos_fBits;   //!
   TBranch        *b_TrackerPos_fX;   //!
   TBranch        *b_TrackerPos_fY;   //!
   TBranch        *b_TrackerPos_fZ;   //!
   TBranch        *b_TrackerPos_sm_fUniqueID;   //!
   TBranch        *b_TrackerPos_sm_fBits;   //!
   TBranch        *b_TrackerPos_sm_fX;   //!
   TBranch        *b_TrackerPos_sm_fY;   //!
   TBranch        *b_TrackerPos_sm_fZ;   //!
   TBranch        *b_BeamBeta;   //!

   readtree(TTree *tree=0);
   virtual ~readtree();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop();
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef readtree_cxx
readtree::readtree(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
   if (tree == 0) {
      TFile *f = (TFile*)gROOT->GetListOfFiles()->FindObject("root/simucat.root");
      if (!f || !f->IsOpen()) {
         f = new TFile("root/simucat.root");
      }
      f->GetObject("t1",tree);

   }
   Init(tree);
}

readtree::~readtree()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t readtree::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t readtree::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void readtree::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set object pointer
   Edep = 0;
   Edep_sm = 0;
   DetID = 0;
   TrackLength = 0;
   DetX = 0;
   DetY = 0;
   DetZ = 0;
   Time = 0;
   HitX = 0;
   HitY = 0;
   HitZ = 0;
   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("RunNumber", &RunNumber, &b_RunNumber);
   fChain->SetBranchAddress("EventID", &EventID, &b_EventID);
   fChain->SetBranchAddress("EGammaCM", EGammaCM, &b_EGammaCM);
   fChain->SetBranchAddress("EGamma", EGamma, &b_EGamma);
   fChain->SetBranchAddress("PosGamma", PosGamma, &b_PosGamma);
   fChain->SetBranchAddress("DirGamma", DirGamma, &b_DirGamma);
   fChain->SetBranchAddress("NumGamma", &NumGamma, &b_NumGamma);
   fChain->SetBranchAddress("NumHit", &NumHit, &b_NumHit);
   fChain->SetBranchAddress("Edep", &Edep, &b_Edep);
   fChain->SetBranchAddress("Edep_sm", &Edep_sm, &b_Edep_sm);
   fChain->SetBranchAddress("DetID", &DetID, &b_DetID);
   fChain->SetBranchAddress("TrackLength", &TrackLength, &b_TrackLength);
   fChain->SetBranchAddress("DetX", &DetX, &b_DetX);
   fChain->SetBranchAddress("DetY", &DetY, &b_DetY);
   fChain->SetBranchAddress("DetZ", &DetZ, &b_DetZ);
   fChain->SetBranchAddress("Time", &Time, &b_Time);
   fChain->SetBranchAddress("TotalEdep", &TotalEdep, &b_TotalEdep);
   fChain->SetBranchAddress("TotalTrackLength", &TotalTrackLength, &b_TotalTrackLength);
   fChain->SetBranchAddress("HitX", &HitX, &b_HitX);
   fChain->SetBranchAddress("HitY", &HitY, &b_HitY);
   fChain->SetBranchAddress("HitZ", &HitZ, &b_HitZ);
   fChain->SetBranchAddress("SourceX", &SourceX, &b_SourceX);
   fChain->SetBranchAddress("SourceY", &SourceY, &b_SourceY);
   fChain->SetBranchAddress("SourceZ", &SourceZ, &b_SourceZ);
   fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_TrackerPos_fUniqueID);
   fChain->SetBranchAddress("fBits", &fBits, &b_TrackerPos_fBits);
   fChain->SetBranchAddress("fX", &fX, &b_TrackerPos_fX);
   fChain->SetBranchAddress("fY", &fY, &b_TrackerPos_fY);
   fChain->SetBranchAddress("fZ", &fZ, &b_TrackerPos_fZ);
//    fChain->SetBranchAddress("fUniqueID", &fUniqueID, &b_TrackerPos_sm_fUniqueID);
//    fChain->SetBranchAddress("fBits", &fBits, &b_TrackerPos_sm_fBits);
//    fChain->SetBranchAddress("fX", &fX, &b_TrackerPos_sm_fX);
//    fChain->SetBranchAddress("fY", &fY, &b_TrackerPos_sm_fY);
//    fChain->SetBranchAddress("fZ", &fZ, &b_TrackerPos_sm_fZ);
   fChain->SetBranchAddress("BeamBeta", &BeamBeta, &b_BeamBeta);
   Notify();
}

Bool_t readtree::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void readtree::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t readtree::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef readtree_cxx
