#include "Run.hh"

Run::Run(G4int nVoxels) : G4Run(), fVoxSum(nVoxels, 0.), fVoxSum2(nVoxels, 0.) {}

void Run::AddEvent(G4double layerEdep, const std::vector<G4int>& touched,
                   const std::vector<G4double>& voxelEdep)
{
    fLayerSum += layerEdep;
    fLayerSum2 += layerEdep * layerEdep;
    for (G4int i : touched) {
        fVoxSum[i] += voxelEdep[i];
        fVoxSum2[i] += voxelEdep[i] * voxelEdep[i];
    }
}

void Run::Merge(const G4Run* run)
{
    const auto* other = static_cast<const Run*>(run);
    fLayerSum += other->fLayerSum;
    fLayerSum2 += other->fLayerSum2;
    for (std::size_t i = 0; i < fVoxSum.size(); ++i) {
        fVoxSum[i] += other->fVoxSum[i];
        fVoxSum2[i] += other->fVoxSum2[i];
    }
    G4Run::Merge(run);
}
