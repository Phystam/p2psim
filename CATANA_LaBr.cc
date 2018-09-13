// Gammacc
// 2012.10.18 Yasuhiro Togano
// main function for TestCsI simulation program

#include "G4RunManager.hh"
#include "G4UIterminal.hh"
#include "G4UImanager.hh"
#include "G4UItcsh.hh"

#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "RunAction.hh"
#include "EventAction.hh"

#ifdef G4VIS_USE
#include "G4VisExecutive.hh"
#endif

int main(int argc, char** argv)
{
  // Construct the default run manager
  G4RunManager *runManager = new G4RunManager;

  // Set mandatory initialization classes
  DetectorConstruction *detector = new DetectorConstruction;
  runManager -> SetUserInitialization(detector);
  runManager -> SetUserInitialization(new PhysicsList);
  G4UIsession *session = 0;

  // set user action classes
  runManager->SetUserAction(new PrimaryGeneratorAction);
  runManager->SetUserAction(new RunAction);
  runManager->SetUserAction(new EventAction);

  if(argc == 1){
#ifdef G4UI_USE_XM
    session = new G4UIXm(argc,argv);
#else
    session  = new G4UIterminal(new G4UItcsh);
#endif
  }

#ifdef G4VIS_USE
  //Initialze visualization
  G4VisManager* visManager = new G4VisExecutive;
  visManager->Initialize();
#endif

  // get the pointer to the UserInterface manager
  G4UImanager *UI = G4UImanager::GetUIpointer();

  if(session){
    // set cuts etc... relating the Physics Processes
    UI -> ApplyCommand("/control/execute init.mac");

    UI -> ApplyCommand("/control/execute vis.mac");
#ifdef G4USE_XM
    UI -> ApplyCommand("/control/execute gui.mac");
#endif
    session -> SessionStart();
    delete session;
  }
  else{
    G4cout << "ELSE" << G4endl;
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    
    UI -> ApplyCommand(command+fileName);
  }

#ifdef G4VIS_USE
  delete visManager;
#endif
  delete runManager;
  return 0;
}

