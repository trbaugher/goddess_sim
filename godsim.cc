/* ************************************************
 * GEANT4 VCGLIB/CAD INTERFACE - basic example
 *
 * File:      cadmesh_example.cc
 *
 * Author:    Christopher M Poole,
 * Email:     mail@christopherpoole.net
 *
 * Date:      20th March, 2011
 **************************************************/

// USER //
#include "DetectorConstruction.hh"
#include "PhysicsList.hh"
#include "PrimaryGeneratorAction.hh"
#include "EventAction.hh"
#include "RunAction.hh"

// GEANT4 //
#include "G4RunManager.hh"
#include "G4UImanager.hh"
//#include "G4UIterminal.hh"
//#include "G4UItcsh.hh"
#include "G4VisExecutive.hh"
//#include "FTFP_BERT.hh"
#include "G4UIExecutive.hh"

std::ofstream outFile;

int main(int argc,char** argv)
{
  G4RunManager* run_manager = new G4RunManager;

  G4VUserDetectorConstruction* detector_construction = new DetectorConstruction;
  run_manager->SetUserInitialization(detector_construction);

  G4VUserPhysicsList* physics_list = new PhysicsList;
  run_manager->SetUserInitialization(physics_list);
  // use FTFP_BERT for now
  //run_manager->SetUserInitialization(new FTFP_BERT);

  G4VUserPrimaryGeneratorAction* primary_generator = new PrimaryGeneratorAction;
  run_manager->SetUserAction(primary_generator);


  EventAction* eventAction = new EventAction();
  run_manager->SetUserAction(eventAction); 
  RunAction* runAction = new RunAction();
  run_manager->SetUserAction(runAction); 

  run_manager->Initialize();
  
  G4VisManager* vis_manager = new G4VisExecutive;
  vis_manager->Initialize();

  G4UImanager * ui_manager = G4UImanager::GetUIpointer();
  if (argc > 1) {
    G4String command = "/control/execute ";
    G4String fileName = argv[1];
    ui_manager->ApplyCommand(command+fileName);
  } else{
    G4UIExecutive* ui = new G4UIExecutive(argc, argv);
    ui_manager->ApplyCommand("/control/execute macros/vis.mac"); 
    ui->SessionStart();
    delete ui; 
  }
  
  return 0;
}


