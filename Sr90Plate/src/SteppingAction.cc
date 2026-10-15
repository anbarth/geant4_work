#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"

SteppingAction::SteppingAction(const DetectorConstruction* det, EventAction* evt)
    : fDetector(det), fEventAction(evt) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return;

    const G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
    if (touch->GetVolume()->GetLogicalVolume() != fDetector->GetLiquidLV()) return;

    // Depth 0 = liquid (copy no. 0); depth 1 = the well it sits in (copy no. = well index)
    const G4int well = touch->GetCopyNumber(1);
    fEventAction->AddEdep(well, edep);
}
