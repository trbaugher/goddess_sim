#include "G4RunManager.hh"
#include "GoddessSD.hh"
#include "GoddessHit.hh"
#include "DetectorConstruction.hh"
#include "G4SystemOfUnits.hh"
#include "G4VPhysicalVolume.hh"
#include "G4Step.hh"
#include "G4VTouchable.hh"
#include "G4SDManager.hh"
#include "G4ios.hh"
#include "G4HCofThisEvent.hh"

GoddessSD::GoddessSD(G4String name):G4VSensitiveDetector(name)
{
  //G4cout << "SD Constructor\n";
  G4RunManager* runManager = G4RunManager::GetRunManager();
  Detector = (DetectorConstruction*)(runManager->GetUserDetectorConstruction());
  n_det_per_ring = Detector->GetNperRing();
  n_sil_E = 2*n_det_per_ring;
  n_sil_dE = 2*n_det_per_ring;
  n_GS = 111;

  //G4cout << "n_sil_dE " << n_sil_dE << G4endl;
  //G4cout << "n_sil_E " << n_sil_E << G4endl;
  sil_hit_dE = new G4int[n_sil_dE];
  sil_hit_E = new G4int[n_sil_E];
  GS_hit = new G4int[n_GS];
  BGO_hit = new G4int[n_GS];
  collectionName.insert("GoddessCollection");
}

GoddessSD::~GoddessSD()
{
  delete [] sil_hit_dE;
  delete [] sil_hit_E;
  delete [] GS_hit;
  delete [] BGO_hit;
}

void GoddessSD::Initialize(G4HCofThisEvent*)
{
  //G4cout << "Initializing SD\n";
  GoddessCollection = new GoddessHitsCollection(SensitiveDetectorName, collectionName[0]);
  for (G4int i = 0; i < n_sil_E; i++) {
    sil_hit_E[i] = -1;
  }
  for (G4int i = 0; i < n_sil_dE; i++) {
    sil_hit_dE[i] = -1;
  }
  for (G4int i = 0; i < n_GS; i++) {
    GS_hit[i] = -1;
    BGO_hit[i] = -1;
  }
}

G4bool GoddessSD::ProcessHits(G4Step* aStep, G4TouchableHistory* ROhist)
{
  //G4cout << "ProcessHit\n";
  G4double edep = aStep->GetTotalEnergyDeposit();
  //G4cout << edep/keV << G4endl;
  if (edep == 0.0) {
    return false;
  }

  // touchable history to get the physical detector geometry
  G4TouchableHistory* theTouchable = (G4TouchableHistory*)(aStep->GetPreStepPoint()->GetTouchable());
  G4VPhysicalVolume* theDetector = theTouchable->GetVolume();
  //G4cout << "name " << theDetector->GetName() << G4endl;
  G4int detector_id = theDetector->GetCopyNo();
  //G4cout << "detector_id " << detector_id << G4endl<<G4endl;
  //G4cout << "edep " << edep/keV << " keV\n";

  /* need this???
  // ReadOut Geometry will be used for the strip ... start at detector level here, 
  // making this section redundant
  G4VPhysicalVolume* ROdetector = 0;
  ROdetector = ROhist->GetVolume();
  G4String ROname = ROdetector->GetName();
  //G4cout << "ROname " << ROname << G4endl;
  G4int RO_detector_id = ROdetector->GetCopyNo();
  //G4cout << "RO_detector_id " << RO_detector_id << G4endl;
 
  // strip-> detector
  G4VPhysicalVolume *ROdetector_1 = ROhist->GetVolume(1);
  //G4cout << "ROdetector_1 " << ROdetector_1 << G4endl;
  G4int RO_detector_id_1 = ROdetector_1->GetCopyNo();
  G4cout << "RO_detector_id_1 " << RO_detector_id_1 << G4endl;
  G4String ROname_1 = ROdetector_1->GetName();
  G4cout << "ROname_1 " << ROname_1 << G4endl;
  */



  // detector number 0-110 are GammaSphere detectors
  if ((detector_id >= 1) && (detector_id <= 110)) {
    //G4cout << "GS HIT\n";
    if (GS_hit[detector_id] == -1) { // new hit
      //G4cout << "new\n";
      GoddessHit *goddessHit = new GoddessHit;
      goddessHit->SetID(detector_id);
      goddessHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());
      goddessHit->AddEdep(edep);
      goddessHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
      GS_hit[detector_id] = GoddessCollection->insert(goddessHit)-1;
    } else { // old hit
      //G4cout << "old\n";
      (*GoddessCollection)[GS_hit[detector_id]]->AddEdep(edep);
    }
  } else if ((detector_id >= 401) && (detector_id <= 510)) {
    if (BGO_hit[detector_id-400] == -1) {
      GoddessHit *goddessHit = new GoddessHit;
      goddessHit->SetID(detector_id);
      goddessHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());
      goddessHit->AddEdep(edep);
      goddessHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
      BGO_hit[detector_id-400] = GoddessCollection->insert(goddessHit)-1;
    } else { // old hit
      //G4cout << "old\n";
      (*GoddessCollection)[BGO_hit[detector_id-400]]->AddEdep(edep);
    }
  } else if ((detector_id >= 200) && (detector_id <= 223)) {
    //G4cout << "E detector ...\n";
    if (sil_hit_E[detector_id-200] == -1) { // new hit
      GoddessHit *goddessHit = new GoddessHit;
      goddessHit->SetID(detector_id);
      goddessHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());
      goddessHit->AddEdep(edep);
      goddessHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
      sil_hit_E[detector_id-200] = GoddessCollection->insert(goddessHit)-1;
    } else { // old hit
      (*GoddessCollection)[sil_hit_E[detector_id-200]]->AddEdep(edep);
    }
  } else if ((detector_id >= 300) && (detector_id <= 323)) {
    //G4cout << "dE detector ...\n";
    if (sil_hit_dE[detector_id-300] == -1) { // new hit
      GoddessHit *goddessHit = new GoddessHit;
      goddessHit->SetID(detector_id);
      goddessHit->SetPosition(aStep->GetPreStepPoint()->GetPosition());
      goddessHit->AddEdep(edep);
      goddessHit->SetTime(aStep->GetPreStepPoint()->GetGlobalTime());
      sil_hit_dE[detector_id-300] = GoddessCollection->insert(goddessHit)-1;
    } else { // old hit
      //G4cout << sil_hit_dE[detector_id-100] << " -> old hit\n";
      (*GoddessCollection)[sil_hit_dE[detector_id-300]]->AddEdep(edep);
    }
  } else {
    G4cerr << "Problem! Unrecognized detector number " << detector_id << G4endl;
  }

  
  return true; 
}

void GoddessSD::EndOfEvent(G4HCofThisEvent* HCE)
{
  //G4cout << "EndOfEvent\n";
  static G4int HCID = -1;
  if (HCID<0) {
    HCID = G4SDManager::GetSDMpointer()->GetCollectionID(collectionName[0]);
  }
  HCE->AddHitsCollection(HCID, GoddessCollection);

  for (G4int i = 0; i < n_sil_E; i++) {
    sil_hit_E[i] = -1;
  }
  for (G4int i = 0; i < n_sil_dE; i++) {
    sil_hit_dE[i] = -1;
  }
  for (G4int i = 0; i < n_GS; i++) {
    GS_hit[i] = -1;
    BGO_hit[i] = -1;
  }
}
