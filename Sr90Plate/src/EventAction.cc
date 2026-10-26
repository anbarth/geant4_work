#include "EventAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"

#include <algorithm>

EventAction::EventAction(const DetectorConstruction* det)
    : fDetector(det), fEdep(DetectorConstruction::kWells, 0.) {}

void EventAction::BeginOfEventAction(const G4Event*)
{
    std::fill(fEdep.begin(), fEdep.end(), 0.);
    // Sized here rather than in the constructor: the action may be created before the
    // geometry is set up. Empty when voxel scoring is off.
    fVoxelEdep.assign(fDetector->GetNumVoxels(), 0.);
}

void EventAction::EndOfEventAction(const G4Event*)
{
    auto* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->AddEvent(fEdep);
    if (!fVoxelEdep.empty()) run->AddVoxelEvent(fVoxelEdep);
}
