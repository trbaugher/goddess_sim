#include "RunAction.hh"
#include "G4Run.hh"
#include "G4String.hh"
#include "G4UIManager.hh"
#include "G4VVisManager.hh"
#include "G4ios.hh"
#include "RunActionMessenger.hh"

extern std::ofstream outFile;

RunAction::RunAction()
{
  name = "";
  runMessenger = new RunActionMessenger(this);
}

RunAction::~RunAction()
{
  delete runMessenger;
}

void RunAction::BeginOfRunAction(const G4Run* aRun)
{
  // file default name is run  number
  if (name == "") {
    name = "default.out";
  }

  G4cout << "Saving in " << name << G4endl;
  outFile.open(name);

}

void RunAction::EndOfRunAction(const G4Run* aRun)
{
  outFile.close();
  G4cout << "Closed output file " << name << G4endl;
}

//void RunAction::SetFileName(G4String fname)
//{
//  name = fname;
//}
