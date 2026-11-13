#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"

SteppingAction::SteppingAction(const DetectorConstruction* det, EventAction* evt)
    : fDetector(det), fEventAction(evt) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return;

    const G4TouchableHandle touch = step->GetPreStepPoint()->GetTouchableHandle();
    const G4LogicalVolume* lv = touch->GetVolume()->GetLogicalVolume();
    const bool scored = fDetector->IsScoringVolume(lv);
    const bool inCellLayer = (lv == fDetector->GetCellLayerLV());
    if (!scored && !inCellLayer) return;

    // Depth 0 = cell layer or liquid (copy no. 0); depth 1 = the well it sits in (copy no. = well index)
    const G4int well = touch->GetCopyNumber(1);
    if (scored) fEventAction->AddEdep(well, edep);

    // Voxel scoring: bin the deposit in the cell layer of the chosen well by (x, y) position.
    // The mid-point of the step is used; steps in the thin layer are much shorter than a voxel.
    if (fDetector->VoxelScoring() && inCellLayer && well == fDetector->GetVoxelWell()) {
        const G4ThreeVector mid = 0.5 * (step->GetPreStepPoint()->GetPosition()
                                         + step->GetPostStepPoint()->GetPosition());
        const G4ThreeVector local = mid - fDetector->GetVoxelWellCentre();   // no rotations in the geometry
        const G4int voxel = fDetector->GetVoxelIndex(local.x(), local.y());
        if (voxel >= 0) fEventAction->AddVoxelEdep(voxel, edep);
    }
}
