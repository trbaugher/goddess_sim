#ifndef RunAction_h
#define RunAction_h 1

#include "G4UserRunAction.hh"
#include "globals.hh"

class G4Run;
class RunActionMessenger;
class RunAction : public G4UserRunAction
{
  public:

  RunAction();
  ~RunAction();

  void BeginOfRunAction(const G4Run* aRun);
  void EndOfRunAction(const G4Run* aRun);
  inline void SetFileName(G4String fname) {name = fname;};
  //void SetFileName(G4String fname);
  inline G4String GetFileName() {return name;};

  private:

  G4String  name;
  RunActionMessenger* runMessenger;
};
#endif
