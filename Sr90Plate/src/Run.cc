#include "Run.hh"

Run::Run(G4int nWells) : G4Run(), fSum(nWells, 0.), fSum2(nWells, 0.) {}

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
    G4Run::Merge(run);
}
