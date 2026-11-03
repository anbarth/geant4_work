#ifndef Run_hh
#define Run_hh

#include "G4Run.hh"
#include "globals.hh"
#include <vector>

// Run object holding, for the whole active layer and for each voxel, the sum of the event energy
// deposits and the sum of their squares (for the statistical error). Merged across threads.
class Run : public G4Run
{
public:
    explicit Run(G4int nVoxels);
    ~Run() override = default;

    // touched = voxels with a non-zero deposit in this event (voxelEdep is indexed by voxel)
    void AddEvent(G4double layerEdep, const std::vector<G4int>& touched,
                  const std::vector<G4double>& voxelEdep);
    void Merge(const G4Run* run) override;

    G4double GetLayerSum() const { return fLayerSum; }
    G4double GetLayerSum2() const { return fLayerSum2; }
    const std::vector<G4double>& GetVoxelSum() const { return fVoxSum; }
    const std::vector<G4double>& GetVoxelSum2() const { return fVoxSum2; }

private:
    G4double fLayerSum = 0.;
    G4double fLayerSum2 = 0.;
    std::vector<G4double> fVoxSum;
    std::vector<G4double> fVoxSum2;
};

#endif
