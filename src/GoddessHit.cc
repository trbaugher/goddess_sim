#include "GoddessHit.hh"

G4Allocator<GoddessHit> GoddessHitAllocator;

GoddessHit::GoddessHit()
{
  edep = 0.0;
  pos = G4ThreeVector(0.0, 0.0, 0.0);
  detector_id = -1;
}

GoddessHit::~GoddessHit()
{;}

GoddessHit::GoddessHit(const GoddessHit& right)
  :G4VHit()
{
  edep = right.edep;
  pos = right.pos;
  detector_id = right.detector_id;
}

const GoddessHit& GoddessHit::operator=(const GoddessHit& right)
{
  edep = right.edep;
  pos = right.pos;
  detector_id = right.detector_id;
  return *this;
}

int GoddessHit::operator==(const GoddessHit& right) const
{
  return ((edep == right.edep) && (pos == right.pos) && (detector_id==right.detector_id));
}

void GoddessHit::Draw()
{;}

void GoddessHit::Print()
{ 
  G4cout << "Printing a GoddessHit\n";
  G4cout << "Detector: " << detector_id << G4endl;
  G4cout << "Position: " << pos << G4endl;
  G4cout << "Energy:   " << edep << G4endl;
  G4cout << "Time:     " << time << G4endl;
}
