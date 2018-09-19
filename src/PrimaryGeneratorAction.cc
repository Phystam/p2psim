// class PrimaryGeneratorAction
// 2012.10.16 Yasuhiro Togano
#include <iostream>
#include "PrimaryGeneratorAction.hh"
#include "PGAMessenger.hh"
#include "DetectorConstruction.hh"
#include "RunAction.hh"

#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ParticleGun.hh"
#include "G4ParticleTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4RandomDirection.hh"
#include "Randomize.hh"
#include "G4LorentzVector.hh"
#include "G4ThreeVector.hh"

#include "TGenPhaseSpace.h"
#include "TVector3.h"
#include "TRandom3.h"
#include "G4SystemOfUnits.hh"
#include "TChain.h"


// C++ libs
#include <cstdio>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <stdlib.h>
//#include <sys/time.h>
#include <cstdarg>
#include <vector>
#include <list> //added
//#include <getopt.h>

// ROOT libs
#include "Riostream.h" // added
#include <TApplication.h>
#include <TCanvas.h>
#include <TChain.h>
#include <TCut.h>
#include <TCutG.h>
#include <TDirectory.h>
#include <TF1.h>
#include <TFile.h>
#include <TGraph.h>
#include <TGaxis.h>
#include <TGraphErrors.h>
#include <TH1.h>
#include <TH2.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TSystem.h>
#include <THStack.h>
#include <TEllipse.h>
#include <TPRegexp.h>
#include <TString.h>
#include <TProfile.h>
#include <TRandom3.h>
#include <TText.h>
#include <TLegend.h>
#include <TLatex.h>
#include "simconst.hh"
//#define DEBUG

// Constructor -------------------------------------------------------------
PrimaryGeneratorAction::PrimaryGeneratorAction()
  : G4VUserPrimaryGeneratorAction(), fParticleGun(0),
    fNumGamma(1), fBeamBeta(0.6)
{

  // fEgamma[0] = 0.662 * MeV;
  // for(G4int i=1;i<5;i++){
  //   fEgamma[i] = 0.;
  // }
  G4int nofParticles = 1;
  
  fParticleGun = new G4ParticleGun(nofParticles);
  // fParticleGun_p[0] = new G4ParticleGun();
  // fParticleGun_p[1] = new G4ParticleGun();
  

  fpgaMessenger = new PGAMessenger(this);
  ftree=new TChain("tree","chain of tree");
  ftree->Add("event_generator/quasi.root");


   amu = 931.494*MeV;// MeV
  //const double ENERGY 	 = 397.6827; // Beam energy (MeV/u)
   ENERGY 	 = 200*MeV; // Beam energy (MeV/u)
  // const int		 A 				 = 12;			 // Mass number of the incident nucleus A
   A = 11;			 // Mass number of the incident nucleus A
   Exe = 10.00000*MeV;  // Residual excitation energy (MeV) (change it for deeply bound states)
   MOM_SIGMA = 105.00*MeV;   // Internal momentum spread (Gauss)
     ISOTROPIC = true;   	 // proton-proton CM scattering
  
  // const double MA  = 11174.950;			 // Nuclear mass of incident A nucleus(MeV/cｲ)
   MA  = 11.0421*amu;			 // Nuclear mass of incident A nucleus(MeV/cｲ)
  // const double MB  = 10252.628 + Exe;// Nuclear mass of the residual fragment B (MeV/cｲ)
   MB  = 10.0517*amu + Exe;// Nuclear mass of the residual fragment B (MeV/cｲ)
   Ma  = 938.279;  			 // Mass of the knocked-out nucleon (MeV/cｲ)
   Mi  = 938.279;  			 // Mass of the scattered (target) nucleon (MeV/cｲ)
  
  //Constants
   PI  = TMath::Pi();

//
   Tkin = ENERGY			* A;					 // Total kinetic energy (MeV) of the projectile
   PA = sqrt(Tkin*(Tkin + 2*MA)); 	 // Total 3-momentum of the beam
   EA = sqrt(MA*MA + PA*PA);      	 // Total energy of the beam
   bA = -PA/EA;		      						 // Beta of the beam
   gA = 1/sqrt(1-bA*bA);	      		 // Gamma of the beam
   S_first = (EA+Mi)*(EA+Mi) - PA*PA;// Invariant mass (A+i) (Mandelstam S)

}

// Destructor --------------------------------------------------------------
PrimaryGeneratorAction::~PrimaryGeneratorAction()
{
  delete fParticleGun;
  // delete fParticleGun_p[0];
  // delete fParticleGun_p[1];
  delete fpgaMessenger;
}

void PrimaryGeneratorAction::GeneratePrimaries(G4Event *anEvent)
{
  //===========================================
  //コード移植
  // TLorentzVector Li11;
  // G4double pmass = 938.27*MeV;
  // G4double amu = 931.494*MeV;
  // G4double Li11mass = 11*amu + 40.72*MeV;
  // G4double He10mass = 10*amu + 50*MeV;
  // {
  //   G4double T=gRandom->Gaus(0.20,0.05)*11*GeV;//GeV
  //   G4double E=T+Li11mass;
  //   G4double p=sqrt(pow(E,2)-pow(Li11mass,2));
  //   TVector3 vect(0,0,p);//GeV
  //   Li11.SetVectM(vect,Li11mass);
  // }
  // //斜めに置いてある標的に入射する
  //反応位置は z=tan(angle)*x 平面上にある
  //  G4double theta = 45*deg;
  G4double theta = 0*deg;
  G4double x = gRandom->Gaus(0,10)*mm;//mean 0mm, Sigma 10mm
  G4double y = gRandom->Gaus(0,10)*mm;
  G4double z = tan(theta)*x*mm+TARGETPOS;
  //zについては厚さ分乱数を振っておく
  G4double thickness = 0.1*mm;
  G4double length = (gRandom->Uniform(-thickness/2./cos(theta),thickness/2./cos(theta)))*mm;
  z+=length;// -厚さ半分から+厚さ半分
  length +=thickness/2./cos(theta);//0から+厚さ
  TVector3 hitpos(x,y,z);

  // //10Heの質量を計算
  // //  G4double mass10He = f10Hemass + Ex/1000.;//GeV
  // TLorentzVector target(0,0,0,pmass);
  // G4double mass_array[]={pmass,pmass,He10mass};
  // TLorentzVector totalcoming=Li11+target;
  // TGenPhaseSpace* decay = new TGenPhaseSpace;
  // decay->SetDecay(totalcoming,3,mass_array);
  // decay->Generate();
TLorentzVector proton[2];
  // proton[0]=*decay->GetDecay(0);
  // proton[1]=*decay->GetDecay(1);

  //Treeから読み込み
  

  //移植
 TLorentzVector proton1,proton2;

 while(!GenerateOneEvent(proton1,proton2)){
 }
 proton[0]=proton1;
 proton[1]=proton2;

 // std::cout << proton[1].Px()<<" "<<proton[1].Py() << " " << proton[1].Pz() <<" "<< proton[1].E()-proton[1].M() << std::endl;

  //===========================================

  G4ParticleTable *particleTable = G4ParticleTable::GetParticleTable();
  //particle
  G4ParticleDefinition* particle_def=particleTable->FindParticle("proton");
  std::cout<<particle_def<<std::endl;
  // fParticleGun_p[0]->SetParticleDefinition(particle_def);
  // fParticleGun_p[1]->SetParticleDefinition(particle_def);
   fParticleGun->SetParticleDefinition(particle_def);
 
  //_______________
  G4ThreeVector beamDir[2], beamPos;
  // direction of beam in xyz space
  beamDir[0].set(proton[0].Px()*MeV,proton[0].Py()*MeV,proton[0].Pz()*MeV);
  beamDir[1].set(proton[1].Px()*MeV,proton[1].Py()*MeV,proton[1].Pz()*MeV);

  beamPos = G4ThreeVector(x,y,z);//pos0


  G4RunManager *runManager = G4RunManager::GetRunManager();
  RunAction *runAction = (RunAction*)runManager->GetUserRunAction();

  //    G4double betaBeam = fBeamBeta;
  G4double gammaE[5];
  Int_t numGamma = fNumGamma;

  // Gamma Event Generation ===================================================

  for(G4int i=0;i<2;i++){
    fParticleGun->SetParticlePosition(beamPos);
    fParticleGun->SetParticleMomentumDirection(beamDir[i].unit());
    fParticleGun->SetParticleEnergy((proton[i].E()-proton[i].M())*MeV);
    //    fParticleGun->SetParticleEnergy(100.*MeV);
    fParticleGun -> GeneratePrimaryVertex(anEvent);
    std::cout << proton[i].Px()<<" "<<proton[i].Py() << " " << proton[i].Pz() <<" "<< proton[i].E()-proton[i].M() << std::endl;     
    //    runAction->SettEGammaCM(gammaE[i],i);
    runAction->SettEGamma((proton[i].E()-proton[i].M())*MeV, i);
    runAction->SettPosGamma(beamPos, i);
    runAction->SettDirGamma(beamDir[i].unit(),i);
  }

  // for(G4int i=0;i<numGamma;i++){
  //   G4ThreeVector gammaCMDir = G4RandomDirection().unit();
  //   G4LorentzVector gammaLAB = DopplerShift(gammaE[i], gammaCMDir,
  // 					    betaBeam*beamDir);
  //   G4double energygammaLAB = gammaLAB.e();
  //   G4ThreeVector dir_gammaLAB = gammaLAB.vect();

  //   fParticleGun -> SetParticlePosition(beamPos);
  //   fParticleGun -> SetParticleMomentumDirection(dir_gammaLAB.unit());
  //   fParticleGun -> SetParticleEnergy(energygammaLAB);
  //   fParticleGun -> GeneratePrimaryVertex(anEvent);
    
  //   runAction->SettEGammaCM(gammaE[i],i);
  //   runAction->SettEGamma(energygammaLAB, i);
  //   runAction->SettPosGamma(beamPos, i);
  //   runAction->SettDirGamma(dir_gammaLAB.unit(),i);

  // }

  // runAction->SettBeamBeta(&betaBeam);
  // runAction->SettNumGamma(&numGamma);
}

// DopplerShift ==============================================================
G4LorentzVector PrimaryGeneratorAction::DopplerShift(G4double energyCM, 
				     G4ThreeVector gRayDirectionCMVector,
				     G4ThreeVector BeambetaVector){

  G4double gRayCMMomentum = energyCM;
  G4ThreeVector gRayCMMomentumVector;
  gRayCMMomentumVector = gRayCMMomentum * gRayDirectionCMVector.unit();

  G4LorentzVector gammaCM(gRayCMMomentumVector, energyCM);
  G4LorentzVector gammaLAB = gammaCM.boost(BeambetaVector);

  return gammaLAB;
}

bool PrimaryGeneratorAction::GenerateOneEvent(TLorentzVector& outvect1,TLorentzVector& outvect2){
  TVector3 Pa;
  Pa.SetX(gRandom->Gaus(0.,MOM_SIGMA));
  Pa.SetY(gRandom->Gaus(0.,MOM_SIGMA));
  Pa.SetZ(gRandom->Gaus(0.,MOM_SIGMA));
  
  //------------ Internal momentum of the residual-----------------
  TVector3 PB;
  PB.SetXYZ( (-Pa.X()), (-Pa.Y()), (-Pa.Z()));
  //Tree variables for the fragmnent
  PBx=PB.X();	PBy=PB.Y(); PBz_rf = PB.Z();
  
  //Off-shell mass of the bound nucleon from the energy conservation
  //in the virtual dissociation A->B+a
  double rrtt = MA*MA + MB*MB - 2*MA * sqrt(MB*MB + Pa.Mag2());
  if(rrtt<=0){
    cout<<"\nERROR off-shell mass!!\n";//non-zero and real off-shell mass
    return false;	
  }
  double Ma_off = sqrt(rrtt);
  
  //Total energies of "a" and "B" in the restframe of "A"
  double Ea  = sqrt(Ma_off*Ma_off + Pa.Mag2()); 
  double EB  = sqrt(MB*MB + PB.Mag2()); 
  
  //------- Lorentz transformations into laboratory system ---------
  std::pair<double, double> lora = Lorentz(gA,bA,Ea,Pa.Z());
  double EaL = lora.first; //cluster energy in lab
  Pa.SetZ(lora.second); 	 //cluster Pz in lab
  
  std::pair<double, double> lorB = Lorentz(gA,bA,EB,PB.Z());
  double EBL = lorB.first; //energy of the residual B in lab
  PB.SetZ(lorB.second);    //Pz of the residual B in lab
  
  //---------- Generating CM scattering process ----------
  double S = Ma_off*Ma_off + Mi*Mi + 2*Mi*EaL; //Mandelstam invariant
  Mandelstam_S = S;//filling tree variable
  //Now generate CM scattering kinematics
  PrimaryGeneratorAction::cm_values CM = CENMASS(S,Ma_off,Mi,Ma,ISOTROPIC);
  if(!CM.good) return false;//non-physical output
  
  TVector3 P1cm(0.,0.,1.), P2cm(0.,0.,1.);
  double phi_rand = gRandom->Uniform(-PI,PI);
  
  P2cm.SetMag(CM.p_clust);
  P2cm.SetTheta(CM.theta_clust);
  P2cm.SetPhi(phi_rand);
  
  P1cm.SetX(-P2cm.X());
  P1cm.SetY(-P2cm.Y());
  P1cm.SetZ(-P2cm.Z());
  
  //------- Calculate realtive to the direction of the quasi-particle (cluster) --------
  double beta_cm = -Pa.Mag()/(EaL+Mi);
  double gamma_cm = 1/sqrt(1-beta_cm*beta_cm);
  
  std::pair<double, double> lora1 = Lorentz(gamma_cm,beta_cm,CM.e_scat,P1cm.Z());	
  std::pair<double, double> lora2 = Lorentz(gamma_cm,beta_cm,CM.e_clust,P2cm.Z());
  
  P1cm.SetZ(lora1.second);
  P2cm.SetZ(lora2.second);
  
  //-------- Rotating back to the beam direction -----------
  TVector3 P1L = DREHUNG(P1cm,Pa);
  TVector3 P2L = DREHUNG(P2cm,Pa);
  
  outvect1.SetVectM(P1L,Mi);
  outvect2.SetVectM(P2L,Ma);
  return true;

}


TVector3 PrimaryGeneratorAction::DREHUNG(TVector3 v1,TVector3 v2) 
{
	double CT = v2.Z()/v2.Mag(); // cos(theta) of v2 wrt. Z-axis
	double ST = sqrt(1-CT*CT);   // sin(theta)
	double CF = v2.X()/v2.Mag()/ST;
	double SF = v2.Y()/v2.Mag()/ST;

	TVector3 v3;
	double _v3x =  v1.X()*CT*CF - v1.Y()*SF + v1.Z()*ST*CF;
	double _v3y =  v1.X()*CT*SF + v1.Y()*CF + v1.Z()*ST*SF;
	double _v3z = -v1.X()*ST   +  v1.Z()*CT;
	v3.SetXYZ(_v3x,_v3y,_v3z);
	return v3;
}

//Kinematical function
double PrimaryGeneratorAction::CINEMA(double x,double y,double z)
{	
	double lambda = x*x + y*y + z*z - 2*x*y - 2*x*z - 2*y*z;
	return lambda;
}

// Calculate elastic scattering kinematics in CM-system (1-target proton, 2-cluster)
PrimaryGeneratorAction::cm_values PrimaryGeneratorAction::CENMASS(double s,double m2off,double m1,double m2,bool isotropic)
{
	PrimaryGeneratorAction::cm_values output;
	output.good=false;
	double X = s;
	double Y = m2off*m2off;
	double Z = m1*m1;
	double sqrs = sqrt(s);

	// Kinematics before the scattering process
	// (with one off-shell mass)
	double p2_off = sqrt(CINEMA(X,Y,Z))/2/sqrs;
	double p1_off = p2_off;
	// CM energies
	double e1_off = (s+Z-Y)/2/sqrs;
	double e2_off = (s+Y-Z)/2/sqrs;

	// Now take the real masses (after scattering)
	Y = m2*m2;  Z = m1*m1;
	//And check whether the kinematical function is ok
	//for this specific kinematical case
	double ERROR_CI = CINEMA(X,Y,Z);
	if(ERROR_CI <= 0.){
		//cout << "\nERROR!!! Kinematical function is negative!";
		return output;
	}

	// Kinematics after the scattering process
	// (with all real masses)
	double p2 = sqrt(CINEMA(X,Y,Z))/2/sqrs;
	double p1 = p2;
	double e1 = (s+Z-Y)/2/sqrs;
	double e2 = (s+Y-Z)/2/sqrs;

	// Let's consider momentum transfer <t> from the
	// target particle 1 to the cluster 2
	double tmax = 2*(m1*m1 - e1_off*e1 - p1_off*p1);//COSINE=(-1)
	double tmin = 2*(m1*m1 - e1_off*e1 + p1_off*p1);//COSINE=(1)
	//cout << "\n\n Tmax = " << tmax;
	//cout << "\n Tmin = " << tmin;
	//cout << "\n Mandels = " << X;

	double t;
	// Generate random momentum transfer for this kinematical case
	if(!isotropic){ t = PrimaryGeneratorAction::get_T(s,tmax); }			//Using parameterized cross sections
	else{ t = gRandom->Uniform(0.,3000000.) * (-1); } //Isotropic scattering

	//double COSINE = (t - m2off*m2off - m2*m2 + 2*e2_off*e2)/(2*p2_off*p2);
	double COSINE = (t - 2*m1*m1 + 2*e1_off*e1)/(2*p1_off*p1);
	if(fabs(COSINE)>=1){//momentum transfer out of range
		//cout << "\nERROR! Scattering cosine is larger than 1";
		return output;
	}

	//CM scattering angles
	double theta1 = acos(COSINE);
	double theta2 = PI - theta1;

	output.e_clust = e2;
	output.p_clust = p2;
	output.theta_clust = theta2;

	output.e_scat = e1;
	output.p_scat = p1;
	output.theta_scat = theta1;

	output.T = t;
	output.good = true;

	return output;
}

// Calculate 3-momentum in CM system of two particles M1 and M2
// when M1 has kinetic energy TLAB and M2 is at rest
double PrimaryGeneratorAction::momentum_CM(double TLAB,double M1,double M2)
{
	//Particle M2 is assumed to be in rest
	double PLAB = sqrt(TLAB*(TLAB + 2*M1));//  Total 3-momentum of an incident particle in Lab
	double ELAB = sqrt(PLAB*PLAB + M1*M1); //  Total energy of an incident particle in lab
	double SLAB = M1*M1 + M2*M2 + 2*M2*ELAB;// Mandelstam invariant S in lab
	double PCM = PLAB * M2 / sqrt(SLAB); // Momentum of both particles in CM frame
	return PCM;
}

std::pair<double, double> PrimaryGeneratorAction::Lorentz(double g,double b,double e,double p)
{
	double eL = g*e - g*b*p;
	double pL = g*p - g*b*e;
	return std::make_pair(eL, pL);
}

// Returns a random value of mandelstam T (in (MeV/c)ｲ units)
// distributed according to the parameterized proton-proton
// invarant cross section. Pass as a parameter "sm" the 
// Mandelstam variable S (in MeVｲ)
// and the maximum possible momentum transfer  
double PrimaryGeneratorAction::get_T(double sm, double max)
{
	//TRandom1 rand;
	//rand.SetSeed(0);

	double Tmax = max*0.000001;//convert to GeVｲ units
	//double Tmin = min*0.000001;

	//double Tmax = -2*pCM*pCM*(1 - cos(PI))*0.000001; //in (GeV/c)ｲ
	Double_t rr = gRandom->Uniform(-1.,1.);// to randomize wrt 90 degrees
	double mandels = sm * 0.000001; // in GeVｲ
	//cout << "\nMandelstam S = " << mandels << "\t Tmax/2 = " << Tmax/2 << "\t Random: " << rr;

	//Probability function from the parameterization
	TF1 * foo = new TF1("foo","[0]*exp(x*[1])*(1+0.02*exp((-6)*x))",Tmax/2,0);
	double c = 0.; 
	if(mandels<=4.79) c = -3283.75 + 3064.11*mandels - 1068.44 *mandels*mandels + 164.844*pow(mandels,3) - 9.48152*pow(mandels,4);
	else if(mandels>4.79) c = -776.822 + 586.016*mandels - 175.347 *mandels*mandels + 26.1823*pow(mandels,3) - 1.94889*pow(mandels,4) + 0.0578352*pow(mandels,5);

	foo->FixParameter(0, 25.); //normalization constant (could be anything)
	foo->FixParameter(1, c);	

	double Trand = foo->GetRandom(Tmax/2,0.); //from 90 to 0 degrees
	if(rr>0) Trand = Tmax - Trand; // symmetrization relative to 90 degrees
	//cout << "\n Tmax/2 = " << Tmax/2 *  1000000;
	//cout << "\n Random T = " << Trand*1000000; 
	delete foo;
	return (Trand*1000000); // returning value in MeVｲ
}
