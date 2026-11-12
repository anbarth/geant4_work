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

    // voxel = global voxel index over all voxelised wells. length = step length, ekin = mean kinetic
    // energy over the step. Only the voxels hit in this event are remembered, so the cost per event
    // does not depend on the number of voxels.
    void AddTrack(G4int voxel, G4double length, G4double ekin)
    {
        if (fLength[voxel] == 0.) fTouched.push_back(voxel);   // length > 0, so 0 means untouched
        fLength[voxel] += length;
        fEnergyLength[voxel] += ekin * length;
    }

private:
    const DetectorConstruction* fDetector;
    std::vector<G4double> fLength;         // per voxel, this event: sum of l
    std::vector<G4double> fEnergyLength;   // per voxel, this event: sum of E*l
    std::vector<G4int> fTouched;           // voxels crossed by an electron in this event
};

#endif
