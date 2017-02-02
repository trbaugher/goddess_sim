/* ************************************************
 * GEANT4 VCGLIB/CAD INTERFACE - basic example
 *
 * File:      DetectorConstruction.cc
 *
 * Author:    Christopher M Poole,
 * Email:     mail@christopherpoole.net
 *
 * Date:      20th March, 2011
 **************************************************/

// USER //
#include "DetectorConstruction.hh"

// CADMESH //
#include "CADMesh.hh"

// GEANT4 //
#include "globals.hh"
#include "G4ThreeVector.hh"
#include "G4Transform3D.hh"

#include "G4Box.hh"
#include "G4Sphere.hh"
#include "G4Tubs.hh"
#include "G4Cons.hh"
#include "G4UnionSolid.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4AssemblyVolume.hh"
#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4SDManager.hh"
#include "GoddessSD.hh"
#include "G4Material.hh"
#include "G4Isotope.hh"
#include "G4Element.hh"
#include "G4VisAttributes.hh"

#include "G4PhysicalConstants.hh"
#include "G4NistManager.hh"

DetectorConstruction::DetectorConstruction()
  // I would really like to know what this construction is ... but it is needed to avois segfault 11
  :world_solid(0), world_logical(0), world_physical(0),
   sil_E_solid(0), sil_E_logical(0), sil_E_physical(0),
   sil_dE_solid(0), sil_dE_logical(0), sil_dE_physical(0),
   ge_crystal_solid(0), ge_crystal_logical(0), ge_crystal_physical(0),
   BGOTypeB_alu_solid(0), BGOTypeB_alu_logical(0), BGOTypeB_alu_physical(0),
   BGOTypeB_solid(0), BGOTypeB_logical(0), BGOTypeB_physical(0),
   BGOTypeC_solid(0), BGOTypeC_logical(0), BGOTypeC_physical(0),
   BGOTypeD_solid(0), BGOTypeD_logical(0), BGOTypeD_physical(0),
	 carbon_disk_solid(0), carbon_disk_logical(0), carbon_disk_physical(0),
	 backplug_solid(0), backplug_logical(0), backplug_physical(0),
	 //coldfinger_solid(0), coldfinger_logical(0), coldfinger_physical(0),
   goddessSD(0)
{
  
  DefineMaterials();
	

  n_sil_per_ring = 12;

  sil_z_offset = 40.843*mm;
  sil_width = 40.*mm;
  sil_length = 75.0*mm;
  sil_E_thick = 1.0*mm;
  sil_dE_thick = 0.065*mm;
  sil_E_radius = 101.10*mm;
  sil_dE_radius = 96.17*mm;

  ge_outer_diameter = 72.0*mm;
  ge_inner_diameter = 0.0*mm;
  ge_length = 76.6*mm;
  // IY Lee paper: crystal taper 2 cm with half angle 7.5 deg = 2.63 mm
  ge_taper_length = 20.00*mm;
  ge_taper_diameter = (ge_outer_diameter-ge_taper_length*tan(7.5*3.14159/180.0))*mm;

	carbon_disk_diameter = ge_outer_diameter;
	carbon_disk_thickness = 5.0*mm;
	
	backplug_diameter = ge_outer_diameter;
	backplug_thickness = 40.0*mm;
	
	//coldfinger_length = 30.0*cm;
	//coldfinger_diameter
	//coldfinger_offset

  heavimet_thickness = 30.5562*mm; // 30.5563 mm is the Z offset from 0,0,0 in the drawings

  absorber1_thickness = 0.127*mm;
  absorber2_thickness = 0.127*mm;

  // focal point to 
  GS_ge_radius = 250.0*mm;
  GS_BGO_radius = 215.0*mm;

  world_size = 2.0*m;

  GSfilename = "GSAngles.txt";
	
  
}

DetectorConstruction::~DetectorConstruction()
{
}
void DetectorConstruction::DefineMaterials()
{
  // Make sure to call this function!!!  

  // define materials
  G4NistManager * nist_manager = G4NistManager::Instance();
  G4Material *silicon = nist_manager->FindOrBuildMaterial("G4_Si");
  G4Material *germanium = nist_manager->FindOrBuildMaterial("G4_Ge");
  G4Material* air = nist_manager->FindOrBuildMaterial("G4_AIR");
  G4Material* galactic = nist_manager->FindOrBuildMaterial("G4_Galactic");
  G4Material* BGO = nist_manager->FindOrBuildMaterial("G4_BGO");
  G4Material* aluminium = nist_manager->FindOrBuildMaterial("G4_Al");
  G4Material* copper = nist_manager->FindOrBuildMaterial("G4_Cu");
  G4Material* tantalum = nist_manager->FindOrBuildMaterial("G4_Ta");
	G4Material* carbon = nist_manager->FindOrBuildMaterial("G4_C");

  G4Material* vacuum = new G4Material("vacuum", 1.0e-5*g/cm3, 1, kStateGas,STP_Temperature,5.0e-9*bar);
  vacuum->AddMaterial(air,1.);

  G4Element* elTa = nist_manager->FindOrBuildElement(73);
  G4Element* elCu = nist_manager->FindOrBuildElement(29);
  G4Element* elNi = nist_manager->FindOrBuildElement(28);
  G4Element* elSi = nist_manager->FindOrBuildElement(14);
  G4Element* elO = nist_manager->FindOrBuildElement(8);
  G4Element* elC = nist_manager->FindOrBuildElement(6);
  G4Element* elH = nist_manager->FindOrBuildElement(1);
  G4Element* elCr = nist_manager->FindOrBuildElement(24);
  G4Element* elMn = nist_manager->FindOrBuildElement(25);
  G4Element* elFe = nist_manager->FindOrBuildElement(26);

  G4Material* G10 = new G4Material("G10", 1.7*g/cm3, 4);
  G10->AddElement(elSi, 1);
  G10->AddElement(elO, 2);
  G10->AddElement(elC, 3);
  G10->AddElement(elH, 4);

  G4Material* heavimet = new G4Material("heavimet", 19.0*g/cm3, 3);
  heavimet->AddElement(elTa, 0.801);
  heavimet->AddElement(elCu, 0.129);
  heavimet->AddElement(elNi, 0.070);

  G4Material* stainless = new G4Material("stainless", 8.06*g/cm3, 6);
  stainless->AddElement(elC, 0.001);
  stainless->AddElement(elSi, 0.007);
  stainless->AddElement(elCr, 0.18);
  stainless->AddElement(elMn, 0.01);
  stainless->AddElement(elFe, 0.712);
  stainless->AddElement(elNi, 0.09);

	// materials 
  flowerpot_mat = aluminium;
  gammasphere_shell_mat = aluminium;
  world_mat = air;
  ge_mat = germanium;
  BGO_mat = BGO;
  chamber_mat = aluminium;
	alu_mat = aluminium;

  HM_mat = heavimet;
  absorber1_mat = copper;
  absorber2_mat = tantalum;

  flange_mat = stainless;
	carbon_disk_mat = carbon;
  detector_mounts_mat = G10;
  ORRUBA_frame_mat = aluminium;
  sil_mat = silicon;
  endcaps_mat = aluminium;
  slewing_ring_mat = aluminium;

  G4cout << *(G4Material::GetMaterialTable()) << G4endl;

  
}


G4VPhysicalVolume* DetectorConstruction::Construct()
{
  
  // +y-axis is up, +z-axis is beam direction, +x axis is beam-left
  G4int i;
  G4int det_num = 0;// continuous detector numbering index
  G4ThreeVector position(0.0, sil_E_radius, -sil_z_offset);
  G4RotationMatrix rotation = G4RotationMatrix();
  G4Transform3D transform = G4Transform3D(rotation, position); 
  G4double phi = 2*3.14159/n_sil_per_ring;

  // define world
  world_solid = new G4Box("world", world_size, world_size, world_size);
  world_logical = new G4LogicalVolume(world_solid, world_mat,"world");
  world_physical = new G4PVPlacement(0, G4ThreeVector(), world_logical, 
                                     "world", 0, false, 0, false);

  // create silicon 1000 um E detectors
  sil_E_solid = new G4Box("sil_E", sil_width/2, sil_E_thick/2, sil_length/2);
  sil_E_logical = new G4LogicalVolume(sil_E_solid, sil_mat, "sil_E");

  // create dE detectors
  sil_dE_solid = new G4Box("sil_dE", sil_width/2, sil_dE_thick/2, sil_length/2);
  sil_dE_logical = new G4LogicalVolume(sil_dE_solid, sil_mat, "sil_dE");

  // backward E detectors (det_num=0-11 clockwise from top (y-axis))
  det_num = 200;
  for (i = 0; i < n_sil_per_ring; i++) {
    //G4cout << "Creating detector # " << det_num << G4endl;
    sil_E_physical = new G4PVPlacement(transform, sil_E_logical, "sil_E", world_logical, false, det_num, false);
    rotation.rotateZ(phi);
    position.rotateZ(phi);
    transform = G4Transform3D(rotation, position);
    det_num++;
  }

  // forward E detectors (det_num=12-23 clockwise from top (y-axis))
  position.setZ(sil_z_offset);
  transform = G4Transform3D(rotation, position);
  for (i = 0; i < n_sil_per_ring; i++) {
    //G4cout << "Creating detector # " << det_num << G4endl;
    sil_E_physical = new G4PVPlacement(transform, sil_E_logical, "sil_E", world_logical, false, det_num, false);
    rotation.rotateZ(phi);
    position.rotateZ(phi);
    transform = G4Transform3D(rotation, position);
    det_num++;
  }

  det_num = 300;
  // backward dE detectors (det_num 100-111 clockwise from top (y-axis)
  position.setY(sil_dE_radius);
  position.setZ(-sil_z_offset);
  transform = G4Transform3D(rotation, position);
  for (i = 0; i < n_sil_per_ring; i++) {
    //G4cout << "Creating detector # " << det_num << G4endl;
    sil_dE_physical = new G4PVPlacement(transform, sil_dE_logical, "sil_dE", world_logical, false, det_num, false);
    rotation.rotateZ(phi);
    position.rotateZ(phi);
    transform = G4Transform3D(rotation, position);
    det_num++;
  }

  // forward dE detectors (det_num 112-123 clockwise from top (y-axis)
  position.setZ(sil_z_offset);
  transform = G4Transform3D(rotation, position);
  for (i = 0; i < n_sil_per_ring; i++) {
    //G4cout << "Creating detector # " << det_num << G4endl;
    sil_dE_physical = new G4PVPlacement(transform, sil_dE_logical, "sil_dE", world_logical, false, det_num, false);
    rotation.rotateZ(phi);
    position.rotateZ(phi);
    transform = G4Transform3D(rotation, position);
    det_num++;
  }

  G4double stlMaxX;
  G4double stlMaxY;
  G4double stlMaxZ;
  G4double stlMinX;
  G4double stlMinY;
  G4double stlMinZ;
  G4ThreeVector offset;
  // GammaSphere
  // BGO type "B"
  // need to be careful with ASYMMETRIC solids ... measured from CAD
  stlMaxX= 89.77 ;stlMaxY= 76.24 ;stlMaxZ= 245.9482 ;
  BGO_length = stlMaxZ-heavimet_thickness;
  offset = G4ThreeVector(-stlMaxX, -stlMaxY, -heavimet_thickness);
  CADMesh *CADTypeB_alu = new CADMesh("models/BGO_det_B_alu.stl", "STL", mm, offset, false);
  BGOTypeB_alu_solid = CADTypeB_alu->TessellatedMesh();
  BGOTypeB_alu_logical = new G4LogicalVolume(BGOTypeB_alu_solid, alu_mat, "BGOtypeB_alu", 0, 0, 0);

  stlMaxX= 88.58 ;stlMaxY= 75.22 ;stlMaxZ= 244.51 ;
  BGO_length = stlMaxZ-heavimet_thickness;
  offset = G4ThreeVector(-stlMaxX, -stlMaxY, -heavimet_thickness-0.02);
  CADMesh *CADTypeB = new CADMesh("models/BGO_det_B_scaled.stl", "STL", mm, offset, false);
  BGOTypeB_solid = CADTypeB->TessellatedMesh();
  BGOTypeB_logical = new G4LogicalVolume(BGOTypeB_solid, BGO_mat, "BGOtypeB", 0, 0, 0);

  // BGO type "C"
  // need to be careful with ASYMMETRIC solids ... measured from CAD
  stlMaxX= 89.80 ;stlMaxY= 76.55 ;stlMaxZ= 245.9482 ;
  offset = G4ThreeVector(-stlMaxX, -stlMaxY, -heavimet_thickness);
  CADMesh *CADTypeC = new CADMesh("models/BGO_det_C.stl", "STL", mm, offset, false);
  BGOTypeC_solid = CADTypeC->TessellatedMesh();
  BGOTypeC_logical = new G4LogicalVolume(BGOTypeC_solid, BGO_mat, "BGOtypeC", 0, 0, 0);
  
  // BGO type "D"
  // import solid bgo shape from CAD
  // need to be careful with ASYMMETRIC solids ... measured from CAD
  stlMaxX= 89.82 ;stlMaxY= 79.35 ;stlMaxZ= 245.9482 ;
  offset = G4ThreeVector(-stlMaxX, -stlMaxY, -heavimet_thickness);
  CADMesh *CADTypeD = new CADMesh("models/BGO_det_D.stl", "STL", mm, offset, false);
  BGOTypeD_solid = CADTypeD->TessellatedMesh();
  BGOTypeD_logical = new G4LogicalVolume(BGOTypeD_solid, BGO_mat, "BGOtypeD", 0, 0, 0);

  // heavimet collimators - 3 types. Assymetric so need to measure offsets in autocad
  // B
  stlMaxX = 44.82; stlMaxY = 38.06; stlMaxZ = 0.0;
  offset = G4ThreeVector(-stlMaxX, -stlMaxY, -stlMaxZ);
  CADMesh *CADHM_B = new CADMesh("models/HM_det_B.stl", "STL", mm, offset, false);
  HM_B_solid = CADHM_B->TessellatedMesh();
  HM_B_logical = new G4LogicalVolume(HM_B_solid, HM_mat, "HM_B", 0, 0, 0);

  // C
  stlMaxX = 44.82; stlMaxY = 38.06; stlMaxZ = 0.0;
  offset = G4ThreeVector(-stlMaxX, -stlMaxY, -stlMaxZ);
  CADMesh *CADHM_C = new CADMesh("models/HM_det_C.stl", "STL", mm, offset, false);
  HM_C_solid = CADHM_C->TessellatedMesh();
  HM_C_logical = new G4LogicalVolume(HM_C_solid, HM_mat, "HM_C", 0, 0, 0);

  // D
  stlMaxX = 44.82; stlMaxY = 38.06; stlMaxZ = 0.0;
  offset = G4ThreeVector(-stlMaxX, -stlMaxY, -stlMaxZ);
  CADMesh *CADHM_D = new CADMesh("models/HM_det_D.stl", "STL", mm, offset, false);
  HM_D_solid = CADHM_D->TessellatedMesh();
  HM_D_logical = new G4LogicalVolume(HM_D_solid, HM_mat, "HM_D", 0, 0, 0);

  G4Tubs *ge_cyl = new G4Tubs("ge_cyl", ge_inner_diameter/2.0, ge_outer_diameter/2.0, (ge_length-ge_taper_length)/2.0, 0.0*deg, 360.0*deg);
  G4Cons *ge_taper = new G4Cons("ge_taper", 0., ge_taper_diameter/2.0, 0.0, ge_outer_diameter/2.0, ge_taper_length/2.0, 0.0, 360.0*degree);
  position = G4ThreeVector(0.0, 0.0,-ge_length/2);
  rotation = G4RotationMatrix();
  transform = G4Transform3D(rotation, position);
  ge_crystal_solid = new G4UnionSolid("GSCrystal", ge_cyl, ge_taper, transform);
  ge_crystal_logical = new G4LogicalVolume(ge_crystal_solid, ge_mat, "GSCrystal");

	carbon_disk_solid = new G4Tubs("carbon_disk", 0, carbon_disk_diameter/2.0, carbon_disk_thickness/2.0, 0.0*deg, 360.0*deg);
	carbon_disk_logical = new G4LogicalVolume(carbon_disk_solid, carbon_disk_mat, "carbon_disk");


  // absorbers
  absorber1_solid = new G4Tubs("absorber1", 0.0, ge_taper_diameter/2.0, absorber1_thickness/2.0, 0.0*degree, 360.0*degree);
  absorber2_solid = new G4Tubs("absorber2", 0.0, ge_taper_diameter/2.0, absorber2_thickness/2.0, 0.0*degree, 360.0*degree);
  absorber1_logical = new G4LogicalVolume(absorber1_solid, absorber1_mat, "absorber1", 0, 0, 0);
  absorber2_logical = new G4LogicalVolume(absorber2_solid, absorber2_mat, "absorber2", 0, 0, 0);

	// BGO backplugs
	stlMinX= 0.10563216 ;stlMinY= 0.026418482 ;stlMinZ= 1e-06 ;
	stlMaxX= 72.000001 ;stlMaxY= 71.973584 ;stlMaxZ= 40.000001 ;
	offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), -((stlMaxZ-stlMinZ)/2.0+stlMinZ));
 	CADMesh *CAD_backplug = new CADMesh("models/backplug.stl","STL", mm, offset, false);
	backplug_solid = CAD_backplug->TessellatedMesh();
	backplug_logical = new G4LogicalVolume(backplug_solid, BGO_mat, "backplug", 0, 0, 0);

  // parse stupid Gammasphere angle file for detector angles
  G4double GStheta;
  G4double GSphi[5];
  G4double rotPhi[5];
  G4int GSnum[5];
  G4String GStype[5];
  std::ifstream GSfile;
  GSfile.open(GSfilename);
  G4String junk;
  // loop over inputfile
  //while (GSfile >> GStheta >> GSnum[0] >> GSnum[1] >> GSnum[2]>> GSnum[3] >> GSnum[4]) {
  while (GSfile >> GStheta) {
    GSfile >> GSnum[0] >> GSnum[1] >> GSnum[2]>> GSnum[3] >> GSnum[4];
    GSfile >> GSphi[0] >> GSphi[1] >> GSphi[2]>> GSphi[3] >> GSphi[4];
    GSfile >> GStype[0] >> GStype[1] >> GStype[2]>> GStype[3] >> GStype[4];
    GSfile >> rotPhi[0] >> rotPhi[1] >> rotPhi[2]>> rotPhi[3] >> rotPhi[4];

    getline(GSfile, junk);
    for (i = 0; i < 5; i++) {
      position.setMag(GS_ge_radius+(ge_length+ge_taper_length)/2.0);
      position.setTheta(GStheta*degree);
      position.setPhi(GSphi[i]*degree);

      rotation=G4RotationMatrix::IDENTITY; 
      rotation.rotateY(GStheta*degree);
      rotation.rotateZ(GSphi[i]*degree);
      transform = G4Transform3D(rotation, position);

      //if((position.getPhi()/degree <  90.0*degree)) {
      ge_crystal_physical = new G4PVPlacement(transform, ge_crystal_logical, "GSCrystal", world_logical, false, GSnum[i], false);

      position.setMag(GS_ge_radius-absorber2_thickness-absorber1_thickness/2.0);
      transform = G4Transform3D(rotation, position);
      absorber1_physical = new G4PVPlacement(transform, absorber1_logical, "absorber1", world_logical, false, 0, false);
      position.setMag(GS_ge_radius-absorber2_thickness/2.0);
      transform = G4Transform3D(rotation, position);
      absorber2_physical = new G4PVPlacement(transform, absorber2_logical, "absorber2", world_logical, false, 0, false);
      //}

			// bgo backplug
			position.setMag(GS_ge_radius+ge_length+carbon_disk_thickness+backplug_thickness/2.0);
			transform = G4Transform3D(rotation, position);
			backplug_physical = new G4PVPlacement(transform, backplug_logical, "backplug", world_logical, false, GSnum[i]+600, false);

			// 
			position.setMag(GS_ge_radius + ge_length + carbon_disk_thickness/2);
			transform = G4Transform3D(rotation, position);
			carbon_disk_physical = new G4PVPlacement(transform, carbon_disk_logical, "carbon_disk", world_logical, false, 0, false);

      //position.setMag(GS_BGO_radius+BGO_length/2.0);
      position.setMag(GS_BGO_radius);
      position.setTheta(GStheta*degree);
      position.setPhi(GSphi[i]*degree);

      rotation=G4RotationMatrix::IDENTITY; 
      rotation.rotateY(GStheta*degree);
      rotation.rotateZ(GSphi[i]*degree);
      rotation.setPhi(rotPhi[i]*degree);
      transform = G4Transform3D(rotation, position);

      //if((position.getPhi()/degree <  90.0*degree)) {
      if (GStype[i][0] == 'B') {
        BGOTypeB_alu_physical = new G4PVPlacement(transform, BGOTypeB_alu_logical, "BGOTypeB_alu", world_logical, false, 0, false);
        BGOTypeB_physical = new G4PVPlacement(transform, BGOTypeB_logical, "BGOTypeB", world_logical, false, GSnum[i]+400, false);
      } else if (GStype[i][0] == 'C') {
        BGOTypeC_physical = new G4PVPlacement(transform, BGOTypeC_logical, "BGOTypeC", world_logical, false, GSnum[i]+400, false);
      } else if (GStype[i][0] == 'D') {
        BGOTypeD_physical = new G4PVPlacement(transform, BGOTypeD_logical, "BGOTypeD", world_logical, false, GSnum[i]+400, false);
      } else {
        G4cerr << "unknown BGO type " << GStype[i] << G4endl;
      }

      position.setMag(GS_BGO_radius-heavimet_thickness);
      transform = G4Transform3D(rotation, position);
      if (GStype[i][0] == 'B') {
        HM_B_physical = new G4PVPlacement(transform, HM_B_logical, "HM_B", world_logical, false, 0, false);
      } else if (GStype[i][0] == 'C') {
        HM_C_physical = new G4PVPlacement(transform, HM_C_logical, "HM_C", world_logical, false, 0, false);
      } else if (GStype[i][0] == 'D') {
        HM_D_physical = new G4PVPlacement(transform, HM_D_logical, "HM_D", world_logical, false, 0, false);
      } else {
        G4cerr << "unknown BGO type " << GStype[i] << G4endl;
      }
      //}
    }
  }

  // Sensitive detectors
  G4SDManager * SDman = G4SDManager::GetSDMpointer();

  if (!goddessSD) {
    goddessSD = new GoddessSD("Goddess_SD");
    SDman->AddNewDetector(goddessSD);
  }

  // needed?
  //G4String ROGeometryName = "SiliconROGeometry";
  //G4VReadOutGeometry* siliconRO = new SiliconROGeometry(ROGemoetryName);
  //siliconRO->BuildROGeometry();
  //siliconSD->SetROGeometry(siliconRO);

  ge_crystal_logical->SetSensitiveDetector(goddessSD);
  BGOTypeB_logical->SetSensitiveDetector(goddessSD);
  BGOTypeC_logical->SetSensitiveDetector(goddessSD);
  BGOTypeD_logical->SetSensitiveDetector(goddessSD);


  // load non-sensitive components with CADMESH

  // Load CAD file as tessellated solid. Program will segfault if file does not exist ...

  // detector mounts
  stlMinX= 1e-06 ;stlMinY= 1e-06 ;stlMinZ= 1e-06 ;
  stlMaxX= 215.53466 ;stlMaxY= 215.53466 ;stlMaxZ= 213.36 ;
  offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), -((stlMaxZ-stlMinZ)/2.0+stlMinZ));
  CADMesh *detector_mounts = new CADMesh("models/detector_mounts.stl", "STL", mm, offset, false);
  detector_mounts_solid = detector_mounts->TessellatedMesh();
  detector_mounts_logical = new G4LogicalVolume(detector_mounts_solid, detector_mounts_mat, "detector_mounts", 0, 0, 0);
  detector_mounts_physical = new G4PVPlacement(0, G4ThreeVector(), detector_mounts_logical, "detector_mounts", world_logical, false, 0, false);

  // ORRUBA frame
  stlMinX= 1e-06 ;stlMinY= 1e-06 ;stlMinZ= 1e-06 ;
  stlMaxX= 245.80749 ;stlMaxY= 245.8072 ;stlMaxZ= 228.6 ;
  offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), -((stlMaxZ-stlMinZ)/2.0+stlMinZ));
  CADMesh* ORRUBA_frame = new CADMesh("models/ORRUBA_frame.stl", "STL", mm, offset, false);
  ORRUBA_frame_solid = ORRUBA_frame->TessellatedMesh();
  ORRUBA_frame_logical = new G4LogicalVolume(ORRUBA_frame_solid, ORRUBA_frame_mat, "ORRUBA_frame", 0, 0, 0);
  ORRUBA_frame_physical = new G4PVPlacement(0, G4ThreeVector(), ORRUBA_frame_logical, "ORRUBA_frame", world_logical, false, 0, false);

  // endcap detectors
  //stlMaxX= 174.0 ;stlMaxY= 174.0 ;stlMaxZ= 177.59437 ;
  //offset = G4ThreeVector(-stlMaxX/2.0, -stlMaxY/2.0, -stlMaxZ/2.0);
  //CADMesh* endcaps = new CADMesh("models/endcaps.stl", "STL", mm, offset, false);
  //endcaps_solid = endcaps->TessellatedMesh();
  //endcaps_logical = new G4LogicalVolume(endcaps_solid, endcaps_mat, "endcaps", 0, 0, 0);
  //endcaps_physical = new G4PVPlacement(0, G4ThreeVector(), endcaps_logical, "endcaps", world_logical, false, 0);

  // scattering chamber
  stlMinX= 0.60727814 ;stlMinY= 1.2124811 ;stlMinZ= 9.9999999e-07 ;
  stlMaxX= 354.99272 ;stlMaxY= 354.38752 ;stlMaxZ= 346.89204 ;
  offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), -((stlMaxZ-stlMinZ)/2.0+stlMinZ));
  CADMesh* scattering_chamber = new CADMesh("models/scattering_chamber.stl", "STL", mm, offset, false);
  scattering_chamber_solid = scattering_chamber->TessellatedMesh();
  scattering_chamber_logical = new G4LogicalVolume(scattering_chamber_solid, chamber_mat, "scattering_chamber", 0, 0, 0);
  scattering_chamber_physical = new G4PVPlacement(0, G4ThreeVector(), scattering_chamber_logical, "scattering_chamber", world_logical, false, 0, false);

  // flowerpot
  stlMinX= 1e-06 ;stlMinY= 1e-06 ;stlMinZ= 228.6 ;
  stlMaxX= 431.8 ;stlMaxY= 431.8 ;stlMaxZ= 568.325 ;
  offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), 0.0);
  CADMesh* flowerpot = new CADMesh("models/flowerpot.stl", "STL", mm, offset, false);
  flowerpot_solid = flowerpot->TessellatedMesh();
  flowerpot_logical = new G4LogicalVolume(flowerpot_solid, chamber_mat, "flowerpot", 0, 0, 0);
  flowerpot_physical = new G4PVPlacement(0, G4ThreeVector(), flowerpot_logical, "flowerpot", world_logical, false, 0, false);

  // slewing ring
  stlMinX= 0.064685807 ;stlMinY= 0.064685807 ;stlMinZ= 103.3071 ;
  stlMaxX= 184.93532 ;stlMaxY= 184.93532 ;stlMaxZ= 137.3071 ;
  offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), 0);
  CADMesh* slewing_ring = new CADMesh("models/slewing_ring.stl", "STL", mm, offset, false);
  slewing_ring_solid = slewing_ring->TessellatedMesh();
  slewing_ring_logical = new G4LogicalVolume(slewing_ring_solid, slewing_ring_mat, "slewing_ring", 0,0,0);
  slewing_ring_physical = new G4PVPlacement(0, G4ThreeVector(), slewing_ring_logical, "slewing_ring", world_logical, false, 0, false);

  // flange
  stlMinX= 1e-06 ;stlMinY= 1e-06 ;stlMinZ= 137.3071 ;
  stlMaxX= 121.92 ;stlMaxY= 121.92 ;stlMaxZ= 218.1606 ;
  offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), 0);
  CADMesh* flange = new CADMesh("models/flange.stl", "STL", mm, offset, false);
  flange_solid = flange->TessellatedMesh();
  flange_logical = new G4LogicalVolume(flange_solid, flange_mat, "flange", 0,0,0);
  flange_physical = new G4PVPlacement(0, G4ThreeVector(), flange_logical, "flange", world_logical, false, 0, false);

  // gammasphere shell
  stlMinX= 9.9999999e-07 ;stlMinY= 0.71656709 ;stlMinZ= 1e-06 ;
  stlMaxX= 1538.4177 ;stlMaxY= 1548.6834 ;stlMaxZ= 1549.4 ;
  offset = G4ThreeVector(-((stlMaxX-stlMinX)/2.0+stlMinX), -((stlMaxY-stlMinY)/2.0+stlMinY), -((stlMaxZ-stlMinZ)/2.0+stlMinZ));
  CADMesh* gammasphere_shell = new CADMesh("models/gammasphere_shell.stl", "STL", mm, offset, false);
  gammasphere_shell_solid = gammasphere_shell->TessellatedMesh();
  gammasphere_shell_logical = new G4LogicalVolume(gammasphere_shell_solid, gammasphere_shell_mat, "gammasphere_shell", 0, 0, 0);
  gammasphere_shell_physical = new G4PVPlacement(0, G4ThreeVector(), gammasphere_shell_logical, "gammasphere_shell", world_logical, false, 0, false);

  G4VisAttributes *yellow = new G4VisAttributes(G4Color(1.0, 1.0, 0.0));
  yellow->SetForceSolid(true);
  G4VisAttributes *gray = new G4VisAttributes(G4Color(0.5, 0.5, 0.5));
  gray->SetForceSolid(true);
  G4VisAttributes* white = new G4VisAttributes(G4Color(1.0, 1.0, 1.0));
  white->SetForceSolid(true);
  G4VisAttributes *red = new G4VisAttributes(G4Color(1.0, 0.0, 0.0));
  red->SetForceSolid(true);
  G4VisAttributes* green = new G4VisAttributes(G4Color(0.0, 1.0, 0.0));
  green->SetForceSolid(true);
  G4VisAttributes * blue = new G4VisAttributes(G4Color(0.0, 0.0, 1.0));
  blue->SetForceSolid(true);
  G4VisAttributes * orange = new G4VisAttributes(G4Color(1.0, 0.647, 0.0));
  orange->SetForceSolid(true);
  G4VisAttributes *invisible = new G4VisAttributes();
  invisible->SetVisibility(false);

  //world_logical->             SetVisAttributes(invisible);
  scattering_chamber_logical->SetVisAttributes(invisible);
  detector_mounts_logical->   SetVisAttributes(invisible);
  ge_crystal_logical->        SetVisAttributes(invisible);
  HM_B_logical->              SetVisAttributes(invisible);
  BGOTypeB_alu_logical->          SetVisAttributes(green);
  BGOTypeB_logical->          SetVisAttributes(red);
  BGOTypeC_logical->          SetVisAttributes(orange);
  HM_C_logical->              SetVisAttributes(invisible);
  BGOTypeD_logical->          SetVisAttributes(blue);
  HM_D_logical->              SetVisAttributes(invisible);
  sil_dE_logical->            SetVisAttributes(invisible);
  sil_E_logical->             SetVisAttributes(invisible);
  ORRUBA_frame_logical->      SetVisAttributes(invisible);
  flowerpot_logical->         SetVisAttributes(invisible);
  absorber1_logical->         SetVisAttributes(invisible);
  absorber2_logical->         SetVisAttributes(invisible);
  flange_logical->            SetVisAttributes(invisible);
  slewing_ring_logical->      SetVisAttributes(invisible);
  gammasphere_shell_logical->SetVisAttributes(invisible);
	backplug_logical->          SetVisAttributes(invisible);
	carbon_disk_logical->       SetVisAttributes(invisible);
  
  //endcaps_logical->SetVisAttributes(white);

  return world_physical;
}

