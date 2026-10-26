#ifndef Run_hh
#define Run_hh

#include "G4Run.hh"
#include "globals.hh"
#include <vector>

// Run object holding, for each well, the sum of event energy deposits and the
// sum of their squares (for the statistical error). Merged across threads.
// The same is kept for each voxel of the chosen well when voxel scoring is on.
class Run : public G4Run
{
public:
    explicit Run(G4int nWells, G4int nVoxels = 0);
    ~Run() override = default;

    void AddEvent(const std::vector<G4double>& edepPerWell);
    void AddVoxelEvent(const std::vector<G4double>& edepPerVoxel);
    void Merge(const G4Run* run) override;

    const std::vector<G4double>& GetSum() const { return fSum; }
    const std::vector<G4double>& GetSum2() const { return fSum2; }
    const std::vector<G4double>& GetVoxelSum() const { return fVoxSum; }
    const std::vector<G4double>& GetVoxelSum2() const { return fVoxSum2; }

private:
    std::vector<G4double> fSum;
    std::vector<G4double> fSum2;
    std::vector<G4double> fVoxSum;
    std::vector<G4double> fVoxSum2;
};

#endif
