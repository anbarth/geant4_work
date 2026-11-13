#ifndef EventAction_hh
#define EventAction_hh

#include "G4UserEventAction.hh"
#include "globals.hh"
#include <vector>

class DetectorConstruction;

class EventAction : public G4UserEventAction
{
public:
    explicit EventAction(const DetectorConstruction* det);
    ~EventAction() override = default;

    void BeginOfEventAction(const G4Event*) override;
    void EndOfEventAction(const G4Event*) override;

    void AddEdep(G4int well, G4double edep) { fEdep[well] += edep; }
    // voxel = global voxel index over all voxelised wells. Only the voxels hit in this event are
    // remembered, so the cost per event does not depend on the number of voxels.
    void AddVoxelEdep(G4int voxel, G4double edep)
    {
        if (fVoxelEdep[voxel] == 0.) fTouched.push_back(voxel);   // edep > 0, so 0 means untouched
        fVoxelEdep[voxel] += edep;
    }

private:
    const DetectorConstruction* fDetector;
    std::vector<G4double> fEdep;        // per well, this event
    std::vector<G4double> fVoxelEdep;   // per voxel, this event (empty if voxel scoring is off)
    std::vector<G4int> fTouched;        // voxels with non-zero deposit in this event
};

#endif
