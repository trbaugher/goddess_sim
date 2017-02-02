

#ifndef DetectorConstruction_H
#define DetectorConstruction_H 1

// GEANT4 //
class G4VSolid;
class G4LogicalVolume;
class G4VPhysicalVolume;

#include "G4Material.hh"
#include "G4Element.hh"
#include "G4ThreeVector.hh"
#include "G4VUserDetectorConstruction.hh"
#include "GoddessSD.hh"

class GoddessSD;


class DetectorConstruction : public G4VUserDetectorConstruction
{
  public:

  DetectorConstruction();
  ~DetectorConstruction();

  inline G4int GetNperRing() {return n_sil_per_ring;};

  G4VPhysicalVolume* Construct();

  private:
  void DefineMaterials();

  G4VSolid * world_solid;
  G4LogicalVolume* world_logical;
  G4VPhysicalVolume* world_physical;

  G4VSolid * sil_E_solid;
  G4LogicalVolume* sil_E_logical;
  G4VPhysicalVolume* sil_E_physical;

  G4VSolid * sil_dE_solid;
  G4LogicalVolume* sil_dE_logical;
  G4VPhysicalVolume* sil_dE_physical;
  
  G4VSolid * detector_mounts_solid;
  G4LogicalVolume * detector_mounts_logical;
  G4VPhysicalVolume * detector_mounts_physical;

  G4VSolid * ORRUBA_frame_solid;
  G4LogicalVolume * ORRUBA_frame_logical;
  G4VPhysicalVolume * ORRUBA_frame_physical;

  G4VSolid * scattering_chamber_solid;
  G4LogicalVolume * scattering_chamber_logical;
  G4VPhysicalVolume * scattering_chamber_physical;

  G4VSolid * ge_crystal_solid;
  G4LogicalVolume * ge_crystal_logical;
  G4VPhysicalVolume * ge_crystal_physical;

  G4VSolid* BGOTypeB_alu_solid;
  G4LogicalVolume* BGOTypeB_alu_logical;
  G4VPhysicalVolume* BGOTypeB_alu_physical;

  G4VSolid* BGOTypeB_solid;
  G4LogicalVolume* BGOTypeB_logical;
  G4VPhysicalVolume* BGOTypeB_physical;

  G4VSolid* BGOTypeC_solid;
  G4LogicalVolume* BGOTypeC_logical;
  G4VPhysicalVolume* BGOTypeC_physical;

  G4VSolid* BGOTypeD_solid;
  G4LogicalVolume* BGOTypeD_logical;
  G4VPhysicalVolume* BGOTypeD_physical;

  G4VSolid* HM_B_solid;
  G4LogicalVolume* HM_B_logical;
  G4VPhysicalVolume* HM_B_physical;

  G4VSolid* HM_C_solid;
  G4LogicalVolume* HM_C_logical;
  G4VPhysicalVolume* HM_C_physical;

  G4VSolid* HM_D_solid;
  G4LogicalVolume* HM_D_logical;
  G4VPhysicalVolume* HM_D_physical;

  //G4VSolid* endcaps_solid;
  //G4LogicalVolume* endcaps_logical;
  //G4VPhysicalVolume* endcaps_physical;

  G4VSolid* absorber1_solid;
  G4LogicalVolume* absorber1_logical;
  G4VPhysicalVolume* absorber1_physical;

  G4VSolid* absorber2_solid;
  G4LogicalVolume* absorber2_logical;
  G4VPhysicalVolume* absorber2_physical;

  G4VSolid* flowerpot_solid;
  G4LogicalVolume* flowerpot_logical;
  G4VPhysicalVolume* flowerpot_physical;

  G4VSolid* slewing_ring_solid;
  G4LogicalVolume* slewing_ring_logical;
  G4VPhysicalVolume* slewing_ring_physical;

  G4VSolid* gammasphere_shell_solid;
  G4LogicalVolume* gammasphere_shell_logical;
  G4VPhysicalVolume* gammasphere_shell_physical;

  G4VSolid* flange_solid;
  G4LogicalVolume* flange_logical;
  G4VPhysicalVolume* flange_physical;

	G4VSolid* carbon_disk_solid;
	G4LogicalVolume* carbon_disk_logical;
	G4VPhysicalVolume* carbon_disk_physical;

	G4VSolid* backplug_solid;
	G4LogicalVolume* backplug_logical;
	G4VPhysicalVolume* backplug_physical;



  // detector parameters
  G4int n_sil_per_ring;
  G4double sil_z_offset;
  G4double sil_length;
  G4double sil_width;
  G4double sil_E_thick;
  G4double sil_dE_thick;
  G4double sil_E_radius;
  G4double sil_dE_radius;

  G4double world_size;

  G4double absorber1_thickness;
  G4double absorber2_thickness;

  G4double ge_inner_diameter;
  G4double ge_outer_diameter;
  G4double ge_length;
  G4double ge_taper_diameter;
  G4double ge_taper_length;

  G4double heavimet_thickness;

  G4double BGO_length;

  G4double GS_ge_radius;
  G4double GS_BGO_radius;

	G4double carbon_disk_diameter;
	G4double carbon_disk_thickness;

	G4double backplug_diameter;
	G4double backplug_thickness;

  G4Material *HM_mat;
  G4Material *sil_mat;
  G4Material *BGO_mat;
  G4Material *world_mat;
  G4Material *detector_mounts_mat;
  G4Material* ge_mat;
  G4Material* chamber_mat;
  G4Material* ORRUBA_frame_mat;
  G4Material* endcaps_mat;
  G4Material* flowerpot_mat;
  G4Material* gammasphere_shell_mat;
  G4Material* slewing_ring_mat;
  G4Material* flange_mat;
  G4Material* absorber2_mat;
  G4Material* absorber1_mat;
	G4Material* carbon_disk_mat;
	G4Material* alu_mat;

  GoddessSD* goddessSD;

  G4String GSfilename;
};

#endif

