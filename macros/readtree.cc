#ifndef __CINT__

#include <iostream>
#include <TMath.h>
#include <TFile.h>
#include <TTree.h>
#include <TF1.h>
#include <TH1.h>
#include <TH2.h>
#include <vector>
#include <TVector3.h>
#endif
using namespace std;
//class vector<double>;
vector<Double_t> Edep;
vector<Double_t> Edep_sm;
vector<Double_t> DetX;
vector<Double_t> DetY;
vector<Double_t> DetZ;
vector<Double_t>* DetX_p=&DetX;
vector<Double_t>* DetY_p=&DetY;
vector<Double_t>* DetZ_p=&DetZ;
vector<Double_t>* Edep_p=&Edep;
vector<Double_t>* Edep_sm_p=&Edep_sm;
Double_t SourceX,SourceY,SourceZ;
Double_t BeamBeta;
void readtree(){
  TFile* file = new TFile("root/simucat.root","readonly");
  TTree* tree = (TTree*)file->Get("t1");
  tree->SetBranchAddress("Edep",&Edep_p);
  tree->SetBranchAddress("Edep_sm",&Edep_sm_p);
  tree->SetBranchAddress("SourceX",&SourceX);
  tree->SetBranchAddress("SourceY",&SourceY);
  tree->SetBranchAddress("SourceZ",&SourceZ);

  tree->SetBranchAddress("DetX",&DetX_p);
  tree->SetBranchAddress("DetY",&DetY_p);
  tree->SetBranchAddress("DetZ",&DetZ_p);
  tree->SetBranchAddress("BeamBeta",&BeamBeta);
  Int_t entries = tree->GetEntries();
  TH1* hEdep=new TH1D("hEdep","Edep",400,0,2000);
  TH1* hEdep_sm=new TH1D("hEdep_sm","Edep Res",400,0,2000);
  TH1* hEdepcor=new TH1D("hEdepcor","Edep Cor",400,0,2000);
  TH1* hEdep_smcor=new TH1D("hEdep_smcor","Edep RESCOR",400,0,2000);
  hEdep_smcor->Sumw2();
  for (int i=0;i<entries;i++){
    tree->GetEntry(i);
    //    BeamBeta-=0.01;
    for(int j=0;j<(int)Edep.size();j++){
      
      TVector3 detpos(DetX[j],DetY[j],DetZ[j]);
      //      TVector3 tgtpos(SourceX,SourceY,SourceZ-10.);
      TVector3 tgtpos(SourceX,SourceY,SourceZ);
      TVector3 beam(0,0,1);
      TVector3 gamma_dir=detpos-tgtpos;
      Double_t theta=beam.Angle(gamma_dir);
      Double_t cost=cos(theta);
      Double_t DoppCorEnergy=Edep[j]*(1.-BeamBeta*cost)/TMath::Sqrt((1.-pow(BeamBeta,2.)));
      Double_t DoppCorEnergy_sm=Edep_sm[j]*(1.-BeamBeta*cost)/TMath::Sqrt((1.-pow(BeamBeta,2.)));
      hEdep->Fill(Edep[j]*1000.);
      hEdep_sm->Fill(Edep_sm[j]*1000.);
      hEdepcor->Fill(DoppCorEnergy*1000.);
      hEdep_smcor->Fill(DoppCorEnergy_sm*1000.,1./entries);
    }
  }
}
    
