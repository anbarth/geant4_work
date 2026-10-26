#include "Run.hh"

Run::Run(G4int nWells, G4int nVoxels)
    : G4Run(), fSum(nWells, 0.), fSum2(nWells, 0.), fVoxSum(nVoxels, 0.), fVoxSum2(nVoxels, 0.) {}

void Run::AddVoxelEvent(const std::vector<G4double>& edepPerVoxel)
{
    for (std::size_t i = 0; i < fVoxSum.size(); ++i) {
        fVoxSum[i] += edepPerVoxel[i];
        fVoxSum2[i] += edepPerVoxel[i] * edepPerVoxel[i];
    }
}

void Run::AddEvent(const std::vector<G4double>& edepPerWell)
{
    for (std::size_t i = 0; i < fSum.size(); ++i) {
        fSum[i] += edepPerWell[i];
        fSum2[i] += edepPerWell[i] * edepPerWell[i];
    }
}

void Run::Merge(const G4Run* run)
{
    const auto* other = static_cast<const Run*>(run);
    for (std::size_t i = 0; i < fSum.size(); ++i) {
        fSum[i] += other->fSum[i];
        fSum2[i] += other->fSum2[i];
    }
    for (std::size_t i = 0; i < fVoxSum.size(); ++i) {
        fVoxSum[i] += other->fVoxSum[i];
        fVoxSum2[i] += other->fVoxSum2[i];
    }
    G4Run::Merge(run);
}
