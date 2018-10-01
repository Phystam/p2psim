#define readtree3_cxx
#include "readtree3.h"
#include <TH2.h>
#include <TH1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <vector>
#include <TVector3.h>
#include <TLorentzVector.h>
#include <iostream>
#include <TRandom3.h>
using namespace std;
void readtree3::Loop()
{
//   In a ROOT session, you can do:
//      Root > .L readtree3.C
//      Root > readtree3 t
//      Root > t.GetEntry(12); // Fill t data members with entry number 12
//      Root > t.Show();       // Show values of entry 12
//      Root > t.Show(16);     // Read and show values of entry 16
//      Root > t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch

  Double_t amu=931.494;
  Double_t KELi=200*11.0421;
  Double_t A=11;
  Double_t Exe=10;
  Double_t MA=11.0421*amu;
  Double_t Mp=938.279;
  Double_t MB  = 10.0517*amu;

  //resolution
  Double_t LiE_sigma=1*11;//MeV
  Double_t Si_res=0.4/sqrt(12);//mm
  Double_t reactionpos_res=0.2;//mm

  TH1* hmass=new TH1D("hmass","Mass He10",100,0,20);
  TH1* hmass_g=new TH1D("hmass_g","Mass He10 gated 50MeV<",100,0,20);
  TH1* hmass_ref=new TH1D("hmass_ref","Mass He10 ref",100,0,20);
  TH1* hprod=new TH2D("hprod","Mass He10",100,0,1,100,0,1);
  TH1* hKE=new TH2D("hKE","Mass He10",1000,0,400,1000,0,400);


   if (fChain == 0) return;

   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;

      //反応位置
      TVector3 reactionpos(PosGamma[0][0],PosGamma[0][1],PosGamma[0][2]);
      reactionpos.SetX(reactionpos.X()+gRandom->Gaus(0,reactionpos_res));
      reactionpos.SetY(reactionpos.Y()+gRandom->Gaus(0,reactionpos_res));

      vector<TVector3> Sivects;
      if((int)DetXSi->size()!=2){ //2こ以外は除外
	continue;
      }
      for(int i=0;i<(int)DetXSi->size();i++){
	//	cout<<(*DetXSi)[i]<<endl;
	TVector3 vect;
	vect.SetXYZ((*DetXSi)[i], (*DetYSi)[i], (*DetZSi)[i]);
	cout<<vect.X()<<endl;
	vect.SetX(vect.X()+gRandom->Gaus(0,Si_res));
	vect.SetY(vect.Y()+gRandom->Gaus(0,Si_res));
	vect.SetZ(vect.Z()+gRandom->Gaus(0,Si_res));
	Sivects.push_back(vect-reactionpos);
      }

      vector<TVector3> Catvects;
      if(DetX->size()!=2){ //2こ以外は除外
	continue;
      }
      bool OK_cat=true;
      for(int i=0;i<(int)DetX->size();i++){
	if((*Edep)[i]<1){
	  OK_cat=false;
	}
	TVector3 vect((*DetX)[i],(*DetY)[i],(*DetZ)[i]);

	Catvects.push_back(vect-reactionpos);
      }
      if(!OK_cat){
	continue;
      }
      //これで2こずつが担保

      //      cout<<"---"<<endl;
      Int_t index_si[]={-1,-1};
      Int_t index_cat[]={-1,-1};
      Int_t counter=0;
      for(int i=0;i<2;i++){
	for(int j=0;j<2;j++){
	  TVector3 sivectu=Sivects[i].Unit();
	  TVector3 catvectu=Catvects[j].Unit();
	  Double_t prod = sivectu.Dot(catvectu);
	  //	  cout<<prod<<endl;
	  if(prod>0.9){
	    index_si[counter] = i;
	    index_cat[counter] = j;
	    counter++;
	  }
	}
      }


      // TLorentzVectorの構成

      TLorentzVector proton[2];
      Double_t prod[2];
      Double_t KEp[2];
      for(int i=0;i<2;i++){
	Double_t KE= (*Edep_sm)[index_cat[i]];
	KEp[i]=KE;
	Double_t P=sqrt(KE*(KE+2.*Mp));
	TVector3 Pvect=Sivects[index_si[i]].Unit()*P;
	proton[index_si[i]].SetVectM(Pvect,Mp);
	TVector3 sivectu=Sivects[index_si[i]].Unit();
	TVector3 catvectu=Catvects[index_cat[i]].Unit();

	prod[i]=sivectu.Dot(catvectu);
      }

      hprod->Fill(prod[0],prod[1]);
      hKE->Fill(KEp[0],KEp[1]);
      //incoming beam
      TLorentzVector Li;
      {
	Double_t P=sqrt(KELi*(KELi+2.*MA));
	TVector3 Pvect(0,0,P);
	Li.SetVectM(Pvect,MA);
      }
      TLorentzVector p0(0,0,0,Mp);

      TLorentzVector Missing=Li+p0-proton[0]-proton[1];
      Double_t Erel=Missing.M()-MB;
      //      cout<<proton[0].E()-Mp<<" "<<proton[1].E()-Mp<<endl;

      hmass->Fill(Erel);

      if(proton[0].E()-Mp>50 && proton[1].E()-Mp>50){
	hmass_g->Fill(Erel);
      }

      //ref

      TLorentzVector proton_ref[2];
      {
	Double_t PGamma = sqrt(EGamma[0]*(EGamma[0]+2.*Mp));
	TVector3 vect(DirGamma[0][0],DirGamma[0][1],DirGamma[0][2]);
	vect = vect*PGamma;
	proton_ref[0].SetVectM(vect,Mp);
      }
      {
	Double_t PGamma = sqrt(EGamma[1]*(EGamma[1]+2.*Mp));
	TVector3 vect(DirGamma[1][0],DirGamma[1][1],DirGamma[1][2]);
	vect = vect*PGamma;
	proton_ref[1].SetVectM(vect,Mp);
      }
      TLorentzVector Missing_true=Li+p0-proton_ref[0]-proton_ref[1];
      Double_t Erel2=Missing_true.M()-MB;
      hmass_ref->Fill(Erel2);
      

   }
}
