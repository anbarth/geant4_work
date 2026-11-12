#ifndef Run_hh
#define Run_hh

#include "G4Run.hh"
#include "globals.hh"
#include <vector>

// Run object holding, for each voxel, the per-event sums needed for the track-length-weighted mean
// electron energy and its statistical error:
//   L = sum of track lengths in the event,  S = sum of (energy x track length) in the event
// and over all events: sum L, sum S, sum L^2, sum S^2, sum L*S. Merged across threads.
class Run : public G4Run
{
public:
    explicit Run(G4int nVoxels);
    ~Run() override = default;

    // touched = voxels crossed by an electron in this event (length / energyLength are indexed by voxel)
    void AddVoxelEvent(const std::vector<G4int>& touched, const std::vector<G4double>& length,
                       const std::vector<G4double>& energyLength);
    void Merge(const G4Run* run) override;

    const std::vector<G4double>& GetL() const { return fL; }
    const std::vector<G4double>& GetS() const { return fS; }
    const std::vector<G4double>& GetLL() const { return fLL; }
    const std::vector<G4double>& GetSS() const { return fSS; }
    const std::vector<G4double>& GetLS() const { return fLS; }

private:
    std::vector<G4double> fL, fS, fLL, fSS, fLS;
};

#endif
