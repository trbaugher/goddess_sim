#ifndef GoddessSD_h
#define GoddessSD_h 1

#include "G4VSensitiveDetector.hh"
#include "globals.hh"

class DetectorConstruction;
class G4HCofThisEvent;
class G4Step;
#include "GoddessHit.hh"

class GoddessSD : public G4VSensitiveDetector
{
  public:

  GoddessSD(G4String);
  ~GoddessSD();
  void Initialize(G4HCofThisEvent*);
  G4bool ProcessHits(G4Step* aStep, G4TouchableHistory* ROHist);
  void EndOfEvent(G4HCofThisEvent*);
  //void clear();
  //void DrawAll();
  //void PrintAll();

  private:
  GoddessHitsCollection* GoddessCollection;
  DetectorConstruction* Detector;
  G4int (*GS_hit);
  G4int (*BGO_hit);
  G4int (*sil_hit_E);
  G4int (*sil_hit_dE);
  G4int n_det_per_ring;
  G4int n_sil_E;
  G4int n_sil_dE;
  G4int n_GS;

};

#endif
