#ifndef SteppingAction_hh
#define SteppingAction_hh

#include "G4UserSteppingAction.hh"

class EventAction;
class DetectorConstruction;

// Adds up the energy deposited in the liquid (and cube) of each well, and histograms the electrons entering the cube
class SteppingAction : public G4UserSteppingAction
{
public:
    SteppingAction(const DetectorConstruction* det, EventAction* evt);
    ~SteppingAction() override = default;

    void UserSteppingAction(const G4Step*) override;

private:
    const DetectorConstruction* fDetector;
    EventAction* fEventAction;
};

#endif
