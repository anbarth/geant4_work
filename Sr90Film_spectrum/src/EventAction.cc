#include "EventAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"

EventAction::EventAction(const DetectorConstruction* det) : fDetector(det) {}

void EventAction::BeginOfEventAction(const G4Event*)
{
    // Sized here rather than in the constructor: the action may be created before the geometry is set up
    const std::size_t nVox = static_cast<std::size_t>(fDetector->GetNumVoxels());
    if (fVoxelEdep.size() != nVox) {
        fVoxelEdep.assign(nVox, 0.);
        fTouched.clear();
    }
    const std::size_t nSpec = static_cast<std::size_t>(fDetector->GetNumSpectrumBins()) + 2;
    if (fSpecCounts.size() != nSpec) {
        fSpecCounts.assign(nSpec, 0.);
        fSpecTouched.clear();
    }
    for (G4int b : fSpecTouched) fSpecCounts[b] = 0.;
    fSpecTouched.clear();
    // Reset only what the previous event touched
    for (G4int v : fTouched) fVoxelEdep[v] = 0.;
    fTouched.clear();
    fLayerEdep = 0.;
}

void EventAction::EndOfEventAction(const G4Event*)
{
    auto* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->AddEvent(fLayerEdep, fTouched, fVoxelEdep);
    if (!fSpecTouched.empty()) run->AddSpectrumEvent(fSpecTouched, fSpecCounts);
}
