/* ************************************************
 * GEANT4 VCGLIB/CAD INTERFACE - basic example
 *
 * File:      PhysicsList.cc
 *
 * Author:    Christopher M Poole,
 * Email:     mail@christopherpoole.net
 *
 * Date:      20th March, 2011
 **************************************************/

// USER//
#include "PhysicsList.hh"

// GEANT4 //
#include "G4EmStandardPhysics.hh"
#include "G4RadioactiveDecayPhysics.hh"
#include "G4DecayPhysics.hh"


PhysicsList::PhysicsList() : G4VModularPhysicsList()
{
  RegisterPhysics(new G4EmStandardPhysics());
  RegisterPhysics(new G4DecayPhysics());
  RegisterPhysics(new G4RadioactiveDecayPhysics());
}

PhysicsList::~PhysicsList()
{;}

void PhysicsList::ConstructParticle()
{
  G4VModularPhysicsList::ConstructParticle();
}

void PhysicsList::ConstructProcess()
{
  G4VModularPhysicsList::ConstructProcess();
}


void PhysicsList::SetCuts()
{
  SetCutValue(defaultCutValue, "gamma");
  SetCutValue(defaultCutValue, "e-");
  SetCutValue(defaultCutValue, "e+");
}

