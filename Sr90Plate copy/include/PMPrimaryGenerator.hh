#ifndef PMPrimaryGenerator_hh
#define PMPrimaryGenerator_hh

#include "G4VUserPrimaryGeneratorAction.hh"

class G4ParticleGun;
class G4Event;
class PMDetectorConstruction;

// Each event = one Sr-90 nucleus at rest, placed uniformly in the active disc.
// Radioactive decay then produces the Sr90 beta, the Y90 daughter, and its beta.
class PMPrimaryGenerator : public G4VUserPrimaryGeneratorAction
{
public:
    PMPrimaryGenerator();
    ~PMPrimaryGenerator() override;

    void GeneratePrimaries(G4Event*) override;

private:
    G4ParticleGun* fGun = nullptr;
    const PMDetectorConstruction* fDetector = nullptr;
    bool fInitialised = false;
};

#endif
