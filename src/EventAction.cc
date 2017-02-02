#include "EventAction.hh"
#include "GoddessHit.hh"
#include "G4SystemOfUnits.hh"
#include "G4Event.hh"
#include "G4EventManager.hh"
#include "G4HCofThisEvent.hh"
#include "G4VHitsCollection.hh"
#include "G4SDManager.hh"
#include "G4UIManager.hh"
#include "G4ios.hh"
#include "G4UnitsTable.hh"
#include "Randomize.hh"

extern std::ofstream outFile;

EventAction::EventAction()
  :goddessCollID(-1)
{
}

EventAction::~EventAction()
{
}

void EventAction::BeginOfEventAction(const G4Event* evt)
{
  G4int evtNb= evt->GetEventID();
  if (evtNb%1000 == 0) {
    G4cout << "Event number " << evtNb << G4endl;
  }

  G4SDManager* SDman = G4SDManager::GetSDMpointer();

  if (goddessCollID == -1) {
    goddessCollID = SDman->GetCollectionID("GoddessCollection");
  }
}

void EventAction::EndOfEventAction(const G4Event* evt)
{
  G4HCofThisEvent* HCE = evt->GetHCofThisEvent();
  GoddessHitsCollection* SHC = 0;
  G4int eventID = evt->GetEventID();

  if (HCE) {
    SHC = (GoddessHitsCollection*)(HCE->GetHC(goddessCollID));
    if (SHC) {
      G4int nHit = SHC->entries();
      if (nHit > 0) {
        //G4cout << "Number of goddess hits in this event " << nHit << G4endl;
        //G4cout << eventID << " " << nHit << G4endl;
        outFile << eventID << " " << nHit << G4endl;
        for (G4int i = 0; i < nHit; i++) {
          G4int detID = (*SHC)[i]->GetID();
          G4double edep = (*SHC)[i]->GetEdep();
          G4double t = (*SHC)[i]->GetTime();
          G4ThreeVector pos = (*SHC)[i]->GetPosition();
          //G4cout << detID << " " << edep << " " << pos.x() << " " << pos.y() << " " << pos.z() << G4endl;
          outFile << detID << " " << edep << " " << t << " " << pos.x() << " " << pos.y() << " " << pos.z() << G4endl;
        }
      }
    }
  }
}
        
