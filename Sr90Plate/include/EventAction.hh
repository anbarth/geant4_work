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
    void AddVoxelEdep(G4int voxel, G4double edep) { fVoxelEdep[voxel] += edep; }

private:
    const DetectorConstruction* fDetector;
    std::vector<G4double> fEdep;        // per well, this event
    std::vector<G4double> fVoxelEdep;   // per voxel of the chosen well, this event (empty if off)
};

#endif
