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

    // bin 0..nBins-1 = energy bins, nBins = overflow. Every call also counts in the total.
    void AddSpectrumCount(G4int bin)
    {
        const G4int total = static_cast<G4int>(fSpecCounts.size()) - 1;
        if (fSpecCounts[bin] == 0) fSpecTouched.push_back(bin);
        ++fSpecCounts[bin];
        if (fSpecCounts[total] == 0) fSpecTouched.push_back(total);
        ++fSpecCounts[total];
    }

private:
    const DetectorConstruction* fDetector;
    std::vector<G4double> fEdep;        // per well, this event
    std::vector<G4double> fVoxelEdep;   // per voxel, this event (empty if voxel scoring is off)
    std::vector<G4int> fTouched;        // voxels with non-zero deposit in this event
    std::vector<G4double> fSpecCounts;  // spectrum bins + overflow + total, this event
    std::vector<G4int> fSpecTouched;    // spectrum entries that are non-zero in this event
};

#endif
