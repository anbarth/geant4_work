#ifndef SteppingAction_hh
#define SteppingAction_hh

#include "G4UserSteppingAction.hh"

class EventAction;
class DetectorConstruction;

// Records, for the electrons in the cell layer, the track length and the energy-weighted track length per voxel
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
