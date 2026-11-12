#include "EventAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"

EventAction::EventAction(const DetectorConstruction* det) : fDetector(det) {}

void EventAction::BeginOfEventAction(const G4Event*)
{
    // Sized here rather than in the constructor: the action may be created before the geometry is set up
    const std::size_t nVox = static_cast<std::size_t>(fDetector->GetNumVoxels());
    if (fLength.size() != nVox) {
        fLength.assign(nVox, 0.);
        fEnergyLength.assign(nVox, 0.);
        fTouched.clear();
    }
    // Reset only what the previous event touched
    for (G4int v : fTouched) {
        fLength[v] = 0.;
        fEnergyLength[v] = 0.;
    }
    fTouched.clear();
}

void EventAction::EndOfEventAction(const G4Event*)
{
    if (fTouched.empty()) return;
    auto* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->AddVoxelEvent(fTouched, fLength, fEnergyLength);
}
