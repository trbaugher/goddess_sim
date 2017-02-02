#include "RunActionMessenger.hh"
#include "RunAction.hh"
#include "G4UICmdWithAString.hh"

RunActionMessenger::RunActionMessenger(RunAction* runAction)
  :theRunAction(runAction)
{
  filenameCmd = new G4UIcmdWithAString("/run/FileName", this);
  filenameCmd->SetGuidance("Set output file name for writing");
  filenameCmd->SetParameterName("fname", true);
  filenameCmd->SetDefaultValue("default.out");
  filenameCmd->AvailableForStates(G4State_PreInit, G4State_Idle);
}

RunActionMessenger::~RunActionMessenger()
{
  delete filenameCmd;
}

void RunActionMessenger::SetNewValue(G4UIcommand* command, G4String newValue)
{
  if (command == filenameCmd) {
    theRunAction->SetFileName(newValue);
  }
}

