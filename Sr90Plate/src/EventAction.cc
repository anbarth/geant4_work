#include "EventAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4RunManager.hh"

#include <algorithm>

EventAction::EventAction() : fEdep(DetectorConstruction::kWells, 0.) {}

void EventAction::BeginOfEventAction(const G4Event*)
{
    std::fill(fEdep.begin(), fEdep.end(), 0.);
}

void EventAction::EndOfEventAction(const G4Event*)
{
    auto* run = static_cast<Run*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->AddEvent(fEdep);
}
