#include "PMSteppingAction.hh"
#include "PMEventAction.hh"
#include "PMDetectorConstruction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

PMSteppingAction::PMSteppingAction(const PMDetectorConstruction* det, PMEventAction* evt)
    : fDetector(det), fEventAction(evt) {}

void PMSteppingAction::UserSteppingAction(const G4Step* step)
{
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return;

    const G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
    if (touch->GetVolume()->GetLogicalVolume() != fDetector->GetLiquidLV()) return;

    // Depth 0 = liquid (copy no. 0); depth 1 = the well it sits in (copy no. = well index)
    const G4int well = touch->GetCopyNumber(1);
    fEventAction->AddEdep(well, edep);
}
