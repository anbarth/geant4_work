#include "PMEventAction.hh"
#include "PMRun.hh"
#include "PMDetectorConstruction.hh"

#include "G4RunManager.hh"

#include <algorithm>

PMEventAction::PMEventAction() : fEdep(PMDetectorConstruction::kWells, 0.) {}

void PMEventAction::BeginOfEventAction(const G4Event*)
{
    std::fill(fEdep.begin(), fEdep.end(), 0.);
}

void PMEventAction::EndOfEventAction(const G4Event*)
{
    auto* run = static_cast<PMRun*>(G4RunManager::GetRunManager()->GetNonConstCurrentRun());
    run->AddEvent(fEdep);
}
