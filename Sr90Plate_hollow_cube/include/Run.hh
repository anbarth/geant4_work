#ifndef Run_hh
#define Run_hh

#include "G4Run.hh"
#include "globals.hh"
#include <vector>

// Run object holding, for each well, the sum of event energy deposits and the
// sum of their squares (for the statistical error). Merged across threads.
// The same is kept for the electron spectrum entering the cube.
class Run : public G4Run
{
public:
    explicit Run(G4int nWells, G4int nSpectrumBins = 0);
    ~Run() override = default;

    void AddEvent(const std::vector<G4double>& edepPerWell);
    // touched = spectrum entries that are non-zero in this event (counts is indexed like the spectrum:
    // bins, then the overflow bin, then the total number of electrons)
    void AddSpectrumEvent(const std::vector<G4int>& touched, const std::vector<G4double>& counts);
    void Merge(const G4Run* run) override;

    const std::vector<G4double>& GetSum() const { return fSum; }
    const std::vector<G4double>& GetSum2() const { return fSum2; }
    const std::vector<G4double>& GetSpectrumSum() const { return fSpecSum; }
    const std::vector<G4double>& GetSpectrumSum2() const { return fSpecSum2; }

private:
    std::vector<G4double> fSum;
    std::vector<G4double> fSum2;
    std::vector<G4double> fSpecSum;    // bins + overflow + total
    std::vector<G4double> fSpecSum2;
};

#endif
