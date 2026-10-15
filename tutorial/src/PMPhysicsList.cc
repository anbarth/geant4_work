#include "PMPhysicsList.hh"

// Constructor
PMPhysicsList::PMPhysicsList()
{
    // EM physics
    RegisterPhysics(new G4EmStandardPhysics());
}

// Destructor (nothing to destroy)
PMPhysicsList::~PMPhysicsList()
{
    
}