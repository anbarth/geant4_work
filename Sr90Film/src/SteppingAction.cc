#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

SteppingAction::SteppingAction(const DetectorConstruction* det, EventAction* evt)
    : fDetector(det), fEventAction(evt) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return;

    const G4StepPoint* pre = step->GetPreStepPoint();
    if (!fDetector->IsActiveLayer(pre->GetPhysicalVolume()->GetLogicalVolume())) return;

    fEventAction->AddLayerEdep(edep);

    // Bin the deposit by (x, y). The mid-point of the step is used; steps in a 25 um layer are far
    // shorter than a voxel. The film has no rotation, so global x, y are used directly.
    const G4ThreeVector mid = 0.5 * (pre->GetPosition() + step->GetPostStepPoint()->GetPosition());
    const G4int voxel = fDetector->GetVoxelIndex(mid.x(), mid.y());
    if (voxel >= 0) fEventAction->AddVoxelEdep(voxel, edep);
}
