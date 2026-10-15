#include "PrimaryGenerator.hh"
#include "DetectorConstruction.hh"

#include "G4ParticleGun.hh"
#include "G4IonTable.hh"
#include "G4Event.hh"
#include "G4RunManager.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include "Randomize.hh"
#include "G4Alpha.hh"

#include <cmath>

PrimaryGenerator::PrimaryGenerator()
{
    fGun = new G4ParticleGun(1);
    fDetector = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
}

PrimaryGenerator::~PrimaryGenerator()
{
    delete fGun;
}

void PrimaryGenerator::GeneratePrimaries(G4Event* event)
{
    // The ion table is only usable once physics is initialised, so set up the
    // ion on the first call rather than in the constructor.
    if (!fInitialised) {
        G4ParticleDefinition* sr90 = G4IonTable::GetIonTable()->GetIon(38, 90, 0.);
        fGun->SetParticleDefinition(sr90);

        //G4ParticleTable* particleTable = G4ParticleTable::GetParticleTable();
        //G4ParticleDefinition* particle = particleTable->FindParticle("e-");
        //fGun->SetParticleDefinition(particle);

        fGun->SetParticleEnergy(0.); // ions are at rest
        fGun->SetParticleMomentumDirection(G4ThreeVector(0, 0, 0));
        fInitialised = true;
    }

    // Uniform position in the active disc
    G4double R = fDetector->GetSourceRadius();
    G4double t = fDetector->GetSourceThickness();
    G4double r = R * std::sqrt(G4UniformRand());
    G4double phi = 2. * CLHEP::pi * G4UniformRand();
    G4double z = (G4UniformRand() - 0.5) * t;

    G4ThreeVector pos = fDetector->GetSourceCentre()
                      + G4ThreeVector(r * std::cos(phi), r * std::sin(phi), z);
    fGun->SetParticlePosition(pos);
    fGun->GeneratePrimaryVertex(event);
}
