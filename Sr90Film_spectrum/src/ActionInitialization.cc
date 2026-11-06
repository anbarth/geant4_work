#include "ActionInitialization.hh"
#include "DetectorConstruction.hh"
#include "PrimaryGenerator.hh"
#include "RunAction.hh"
#include "EventAction.hh"
#include "SteppingAction.hh"

#include "G4RunManager.hh"

void ActionInitialization::BuildForMaster() const
{
    auto* det = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    SetUserAction(new RunAction(det));
}

void ActionInitialization::Build() const
{
    auto* det = static_cast<const DetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    SetUserAction(new PrimaryGenerator());
    SetUserAction(new RunAction(det));

    auto* eventAction = new EventAction(det);
    SetUserAction(eventAction);
    SetUserAction(new SteppingAction(det, eventAction));
}
