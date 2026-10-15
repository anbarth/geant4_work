#ifndef PMSteppingAction_hh
#define PMSteppingAction_hh

#include "G4UserSteppingAction.hh"

class PMEventAction;
class PMDetectorConstruction;

// Adds up the energy deposited in the liquid of each well
class PMSteppingAction : public G4UserSteppingAction
{
public:
    PMSteppingAction(const PMDetectorConstruction* det, PMEventAction* evt);
    ~PMSteppingAction() override = default;

    void UserSteppingAction(const G4Step*) override;

private:
    const PMDetectorConstruction* fDetector;
    PMEventAction* fEventAction;
};

#endif
