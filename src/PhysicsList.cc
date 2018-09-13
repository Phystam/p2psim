// PhysicsList.cc
// Based on the TestEM5
// 2012.10.29 Yasuhiro Togano

#include "PhysicsList.hh"
#include "PhysicsListMessenger.hh"
#include "PhysListEmStandard.hh"

#include "G4EmStandardPhysics.hh"
#include "G4EmLivermorePhysics.hh"

#include "G4Decay.hh"
#include "StepMax.hh"

#include "G4UnitsTable.hh"
#include "G4ParticleDefinition.hh"
#include "G4ProcessManager.hh"

// Boson
#include "G4ChargedGeantino.hh"
#include "G4Geantino.hh"
#include "G4Gamma.hh"
#include "G4OpticalPhoton.hh"

// leptons
#include "G4MuonPlus.hh"
#include "G4MuonMinus.hh"
#include "G4NeutrinoMu.hh"
#include "G4AntiNeutrinoMu.hh"

#include "G4Electron.hh"
#include "G4Positron.hh"
#include "G4NeutrinoE.hh"
#include "G4AntiNeutrinoE.hh"

// Hadrons
#include "G4MesonConstructor.hh"
#include "G4BaryonConstructor.hh"
#include "G4IonConstructor.hh"

#include "G4SystemOfUnits.hh"

// Constructor =============================================================
PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  fPhysMessenger = new PhysicsListMessenger(this);

  // EM physics
  fEMName = G4String("local");
  fEMPhysicsList = new PhysListEmStandard(fEMName);

  defaultCutValue = 1.*mm;
  
  fCutForGamma = defaultCutValue;
  fCutForElectron = defaultCutValue;
  fCutForPositron = defaultCutValue;
  
  SetVerboseLevel(1);
}

// Destructor =================================================================
PhysicsList::~PhysicsList()
{
  delete fEMPhysicsList;
  delete fPhysMessenger;
}

// ConstructParticle =========================================================
void PhysicsList::ConstructParticle()
{
  // Geantinos
  G4Geantino::GeantinoDefinition();
  G4ChargedGeantino::ChargedGeantinoDefinition();
  
  // gamma
  G4Gamma::GammaDefinition();
  
  // leptons
  G4Electron::ElectronDefinition();
  G4Positron::PositronDefinition();
  G4MuonPlus::MuonPlusDefinition();
  G4MuonMinus::MuonMinusDefinition();

  G4NeutrinoE::NeutrinoEDefinition();
  G4AntiNeutrinoE::AntiNeutrinoEDefinition();
  G4NeutrinoMu::NeutrinoMuDefinition();
  G4AntiNeutrinoMu::AntiNeutrinoMuDefinition();  

  // mesons
  G4MesonConstructor mConstructor;
  mConstructor.ConstructParticle();

  // barions
  G4BaryonConstructor bConstructor;
  bConstructor.ConstructParticle();

  // ions
  G4IonConstructor iConstructor;
  iConstructor.ConstructParticle();
}

// ConstructProcess ==========================================================
void PhysicsList::ConstructProcess()
{
  AddTransportation();
  fEMPhysicsList->ConstructProcess();
  AddDecay();
  AddStepMax();
}

// AddDecay ==================================================================
void PhysicsList::AddDecay()
{
  // Add decay process
  G4Decay *decayProcess = new G4Decay();
  
  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();

    if(decayProcess->IsApplicable(*particle) && !particle->IsShortLived()){
      pmanager->AddProcess(decayProcess);
      // set ordering
      pmanager->SetProcessOrdering(decayProcess, idxPostStep);
      pmanager->SetProcessOrdering(decayProcess, idxAtRest);
    }
  }
}

// AddStepMax ================================================================
//Step limitation seen as a process
void PhysicsList::AddStepMax()
{
  StepMax *stepMaxProcess = new StepMax();

  theParticleIterator->reset();
  while( (*theParticleIterator)() ){
    G4ParticleDefinition *particle = theParticleIterator->value();
    G4ProcessManager *pmanager = particle->GetProcessManager();
    if(stepMaxProcess->IsApplicable(*particle) && !particle->IsShortLived()){
      pmanager->AddDiscreteProcess(stepMaxProcess);
    }
  }
}

// AddPhysicsList ============================================================
void PhysicsList::AddPhysicsList(const G4String& name)
{
  if(verboseLevel>-1){
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">" << G4endl;
  }

  if(name == fEMName) return;
  
  if(name == "local"){
    fEMName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new PhysListEmStandard(name);
  }
  else if(name == "emlivermore"){
    fEMName = name;
    delete fEMPhysicsList;
    fEMPhysicsList = new G4EmLivermorePhysics();
  }
  else {
    G4cout << "PhysicsList::AddPhysicsList: <" << name << ">"
	   << "is not defined"
	   << G4endl;
  }
}

// SetCuts ===================================================================
void PhysicsList::SetCuts()
{
  if(verboseLevel>0){
    G4cout << "PhysicsList::SetCuts:";
    G4cout << "CutLength : " << G4BestUnit(defaultCutValue,"Length") << G4endl;
  }
  
  // set cuts for EM process
  // Always this order!: gamma->e+->e-
  // this is because internal process inside e+/e- need cut value for gamma.
  SetCutValue(fCutForGamma,"gamma");
  SetCutValue(fCutForElectron,"e-");
  SetCutValue(fCutForPositron,"e+");

  if(verboseLevel>0) DumpCutValuesTable();
}


// SetCutForGamma ============================================================
void PhysicsList::SetCutForGamma(G4double cut)
{
  fCutForGamma = cut;
  SetParticleCuts(fCutForGamma, G4Gamma::Gamma());
}

// SetCutForElectron =========================================================
void PhysicsList::SetCutForElectron(G4double cut)
{
  fCutForElectron = cut;
  SetParticleCuts(fCutForElectron, G4Electron::Electron());
}

// SetCutForPositron =========================================================
void PhysicsList::SetCutForPositron(G4double cut)
{
  fCutForPositron = cut;
  SetParticleCuts(fCutForPositron, G4Positron::Positron());
}
