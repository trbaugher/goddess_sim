#ifndef RunActionMessenger_h
#define RunActionMessenger_h 1

#include "G4UImessenger.hh"
#include "globals.hh"

class RunAction;
class G4UIcmdWithAString;

class RunActionMessenger : public G4UImessenger
{
  public:

  RunActionMessenger(RunAction*);
  ~RunActionMessenger();
  void SetNewValue(G4UIcommand*, G4String);

  private:

  RunAction* theRunAction;
  G4UIcmdWithAString* filenameCmd;
};
#endif
