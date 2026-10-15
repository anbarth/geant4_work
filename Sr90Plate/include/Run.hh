#ifndef Run_hh
#define Run_hh

#include "G4Run.hh"
#include "globals.hh"
#include <vector>

// Run object holding, for each well, the sum of event energy deposits and the
// sum of their squares (for the statistical error). Merged across threads.
class Run : public G4Run
{
public:
    explicit Run(G4int nWells);
    ~Run() override = default;

    void AddEvent(const std::vector<G4double>& edepPerWell);
    void Merge(const G4Run* run) override;

    const std::vector<G4double>& GetSum() const { return fSum; }
    const std::vector<G4double>& GetSum2() const { return fSum2; }

private:
    std::vector<G4double> fSum;
    std::vector<G4double> fSum2;
};

#endif
