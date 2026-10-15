#ifndef PMRunAction_hh
#define PMRunAction_hh

#include "G4UserRunAction.hh"

class PMDetectorConstruction;

class PMRunAction : public G4UserRunAction
{
public:
    explicit PMRunAction(const PMDetectorConstruction* det);
    ~PMRunAction() override = default;

    G4Run* GenerateRun() override;
    void EndOfRunAction(const G4Run* run) override;

private:
    const PMDetectorConstruction* fDetector;

    // Used only to convert dose per decay into a dose rate in the printout
    static constexpr double kActivity_Bq = 1.0e6;
};

#endif
