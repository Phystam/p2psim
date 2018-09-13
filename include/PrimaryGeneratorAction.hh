// Primary Generator Action 
// 2012.10.16 Yasuhiro Togano

#ifndef PRIMARYGENERATORACTION_H
#define PRIMARYGENERATORACTION_H 1

#include "G4VUserPrimaryGeneratorAction.hh"
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4LorentzVector.hh"
#include "TChain.h"
#include "TMath.h"
#include "TLorentzVector.h"
class G4ParticleGun;
class G4Event;
class PGAMessenger;
class TChain;
class TLorentzVector;
// Primary Generator action
// Define the primary event to be generated.
// If "/run/BeamOn" command was used, the primary event described 
// in this class will be generated.

class PrimaryGeneratorAction : public G4VUserPrimaryGeneratorAction
{
private:
  G4ParticleGun *fParticleGun;
  G4ParticleGun *fParticleGun_p[2];
  
  G4int fNumGamma;
  G4double fEgamma[5];
  G4double fBeamBeta;
  
  PGAMessenger *fpgaMessenger;
  TChain* ftree;



public:
  PrimaryGeneratorAction();
  virtual ~PrimaryGeneratorAction();
  virtual void GeneratePrimaries(G4Event *anEvent);

  G4LorentzVector DopplerShift(G4double energyCM, 
			       G4ThreeVector gRayDirectionCMVector,
			       G4ThreeVector betaVector);

  void SetNumGamma(G4int n){fNumGamma = n;}
  void SetEgamma1(G4double e){fEgamma[0] = e;}
  void SetEgamma2(G4double e){fEgamma[1] = e;}
  void SetEgamma3(G4double e){fEgamma[2] = e;}
  void SetEgamma4(G4double e){fEgamma[3] = e;}
  void SetEgamma5(G4double e){fEgamma[4] = e;}
  void SetBeamBeta(G4double b){fBeamBeta = b;}

  bool GenerateOneEvent(TLorentzVector& vect1,TLorentzVector& vect2);

  //variables for GenerateOneEvent

   G4double amu;
  // double ENERGY 	 = 397.6827; // Beam energy (MeV/u)
   G4double ENERGY;
  //  int		 A 				 = 12;			 // Mass number of the incident nucleus A
   int A;
   G4double Exe;
   G4double MOM_SIGMA;
   bool   ISOTROPIC;
  
  //  double MA  = 11174.950;			 // Nuclear mass of incident A nucleus(MeV/cｲ)
   G4double MA;
  //  double MB  = 10252.628 + Exe;// Nuclear mass of the residual fragment B (MeV/cｲ)
   G4double MB;
   G4double Ma;
   G4double Mi;
  
  //Constants
   G4double PI;
Double_t theta_1, theta_2, theta_B;		
Double_t phi_1,   phi_2,   phi_B;		
Double_t P1x,     P1y,     P1z;
Double_t P2x,     P2y,     P2z;
Double_t PBx,     PBy,     PBz_lab,  PBz_rf;
Double_t E1,      E2, 	   EB;
Double_t th1_cm,  th2_cm, phi1_cm, phi2_cm;
Double_t P1_cm,   P2_cm;
Double_t Moff, Mandelstam_T, Opang, Dif_phi;
Double_t Mandelstam_S;
TLorentzVector P1LL,P2LL,PALL,PBLL,MissingMass;
Double_t ReconstructedEx;

  double Tkin;
  double PA ;
  double EA ;
  double bA ;
  double gA ;
  double S_first ;

struct cm_values
{
	//Internal cluster
	double	e_clust;
	double	p_clust;
	double	theta_clust;
	//Scattered particle
	double	e_scat;
	double	p_scat;
	double	theta_scat;
	//indicates satisfactory kinematics (i.e. energy & momentum conservation)
	bool 	good;
	double T;
};

  TVector3 DREHUNG(TVector3 v1,TVector3 v2);
  double CINEMA(double x,double y,double z);
  PrimaryGeneratorAction::cm_values CENMASS(double s,double m2off,double m1,double m2,bool isotropic);
  double momentum_CM(double TLAB,double M1,double M2);
  std::pair<double, double> Lorentz(double g,double b,double e,double p);
  double get_T(double sm, double max);



};




#endif
