#ifndef GoddessHit_h
#define GoddessHit_h 1

#include "G4VHit.hh"
#include "G4THitsCollection.hh"
#include "G4Allocator.hh"
#include "G4ThreeVector.hh"

class GoddessHit : public G4VHit
{
  public:
  GoddessHit();
  ~GoddessHit();
  GoddessHit (const GoddessHit&);
  const GoddessHit& operator=(const GoddessHit&);
  G4int operator==(const GoddessHit&) const;
  inline void* operator new(size_t);
  inline void operator delete(void*);

  void Draw();
  void Print();

  inline void AddEdep(G4double e) {edep += e;};
  inline void SetPosition(G4ThreeVector p) {pos = p;};
  inline void SetID(G4int id) {detector_id = id;};
  inline void SetTime(G4double t) {time = t;};
  inline G4double GetEdep() {return edep;};
  inline G4ThreeVector GetPosition() {return pos;};
  inline G4int GetID() {return detector_id;};
  inline G4double GetTime() {return time;};

  private:

  G4double edep;
  G4double time;
  G4ThreeVector pos;
  G4int detector_id;
};

typedef G4THitsCollection<GoddessHit> GoddessHitsCollection;

extern G4Allocator<GoddessHit> GoddessHitAllocator;

inline void* GoddessHit::operator new(size_t)
{
  void* aHit;
  aHit = (void*) GoddessHitAllocator.MallocSingle();
  return aHit;
}

inline void GoddessHit::operator delete(void* aHit)
{
  GoddessHitAllocator.FreeSingle((GoddessHit*) aHit);
}

#endif
