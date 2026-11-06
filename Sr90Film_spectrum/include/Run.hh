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
    explicit Run(G4int nVoxels, G4int nSpectrumBins = 0);
    ~Run() override = default;

    // touched = voxels with a non-zero deposit in this event (voxelEdep is indexed by voxel)
    void AddEvent(G4double layerEdep, const std::vector<G4int>& touched,
                  const std::vector<G4double>& voxelEdep);
    // touched = spectrum entries that are non-zero in this event (counts is indexed like the spectrum:
    // bins, then the overflow bin, then the total number of electrons)
    void AddSpectrumEvent(const std::vector<G4int>& touched, const std::vector<G4double>& counts);
    void Merge(const G4Run* run) override;

    G4double GetLayerSum() const { return fLayerSum; }
    G4double GetLayerSum2() const { return fLayerSum2; }
    const std::vector<G4double>& GetVoxelSum() const { return fVoxSum; }
    const std::vector<G4double>& GetVoxelSum2() const { return fVoxSum2; }
    const std::vector<G4double>& GetSpectrumSum() const { return fSpecSum; }
    const std::vector<G4double>& GetSpectrumSum2() const { return fSpecSum2; }

private:
    G4double fLayerSum = 0.;
    G4double fLayerSum2 = 0.;
    std::vector<G4double> fVoxSum;
    std::vector<G4double> fVoxSum2;
    std::vector<G4double> fSpecSum;    // bins + overflow + total
    std::vector<G4double> fSpecSum2;
};

#endif
