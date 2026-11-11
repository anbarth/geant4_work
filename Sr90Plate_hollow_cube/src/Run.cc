#include "Run.hh"

Run::Run(G4int nWells, G4int nSpectrumBins)
    : G4Run(), fSum(nWells, 0.), fSum2(nWells, 0.),
      fSpecSum(nSpectrumBins > 0 ? nSpectrumBins + 2 : 0, 0.),
      fSpecSum2(nSpectrumBins > 0 ? nSpectrumBins + 2 : 0, 0.) {}

void Run::AddSpectrumEvent(const std::vector<G4int>& touched, const std::vector<G4double>& counts)
{
    for (G4int i : touched) {
        fSpecSum[i] += counts[i];
        fSpecSum2[i] += counts[i] * counts[i];
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
    for (std::size_t i = 0; i < fSpecSum.size(); ++i) {
        fSpecSum[i] += other->fSpecSum[i];
        fSpecSum2[i] += other->fSpecSum2[i];
    }
    G4Run::Merge(run);
}
