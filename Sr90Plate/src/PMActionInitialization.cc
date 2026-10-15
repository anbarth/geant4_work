#include "PMActionInitialization.hh"
#include "PMDetectorConstruction.hh"
#include "PMPrimaryGenerator.hh"
#include "PMRunAction.hh"
#include "PMEventAction.hh"
#include "PMSteppingAction.hh"

#include "G4RunManager.hh"

void PMActionInitialization::BuildForMaster() const
{
    auto* det = static_cast<const PMDetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());
    SetUserAction(new PMRunAction(det));
}

void PMActionInitialization::Build() const
{
    auto* det = static_cast<const PMDetectorConstruction*>(
        G4RunManager::GetRunManager()->GetUserDetectorConstruction());

    SetUserAction(new PMPrimaryGenerator());
    SetUserAction(new PMRunAction(det));

    auto* eventAction = new PMEventAction();
    SetUserAction(eventAction);
    SetUserAction(new PMSteppingAction(det, eventAction));
}
