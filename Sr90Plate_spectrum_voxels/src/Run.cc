#include "Run.hh"

Run::Run(G4int nVoxels)
    : G4Run(), fL(nVoxels, 0.), fS(nVoxels, 0.), fLL(nVoxels, 0.), fSS(nVoxels, 0.), fLS(nVoxels, 0.) {}

void Run::AddVoxelEvent(const std::vector<G4int>& touched, const std::vector<G4double>& length,
                        const std::vector<G4double>& energyLength)
{
    for (G4int i : touched) {
        fL[i] += length[i];
        fS[i] += energyLength[i];
        fLL[i] += length[i] * length[i];
        fSS[i] += energyLength[i] * energyLength[i];
        fLS[i] += length[i] * energyLength[i];
    }
}

void Run::Merge(const G4Run* run)
{
    const auto* other = static_cast<const Run*>(run);
    for (std::size_t i = 0; i < fL.size(); ++i) {
        fL[i] += other->fL[i];
        fS[i] += other->fS[i];
        fLL[i] += other->fLL[i];
        fSS[i] += other->fSS[i];
        fLS[i] += other->fLS[i];
    }
    G4Run::Merge(run);
}
