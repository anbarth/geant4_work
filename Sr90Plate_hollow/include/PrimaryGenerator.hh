#ifndef PrimaryGenerator_hh
#define PrimaryGenerator_hh

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;
class DetectorConstruction;

// Each event = one Sr-90 nucleus at rest, placed uniformly in the active disc.
// Radioactive decay then produces the Sr90 beta, the Y90 daughter, and its beta.
class PrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
    PrimaryGenerator();
    ~PrimaryGenerator() override;

    void GeneratePrimaries(G4Event*) override;

private:
    G4ParticleGun* fGun = nullptr;
    const DetectorConstruction* fDetector = nullptr;
    bool fInitialised = false;
};

#endif
