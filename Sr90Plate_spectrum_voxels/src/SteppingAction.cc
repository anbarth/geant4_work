#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4TouchableHandle.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4Electron.hh"

SteppingAction::SteppingAction(const DetectorConstruction* det, EventAction* evt)
    : fDetector(det), fEventAction(evt) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    // Electrons only (primary beta electrons and secondary delta electrons; positrons are not counted)
    if (step->GetTrack()->GetDefinition() != G4Electron::Definition()) return;

    const G4double length = step->GetStepLength();
    if (length <= 0.) return;

    const G4StepPoint* pre = step->GetPreStepPoint();
    const G4TouchableHandle touch = pre->GetTouchableHandle();
    const G4LogicalVolume* lv = touch->GetVolume()->GetLogicalVolume();
    if (!fDetector->IsCellLayer(lv)) return;

    // Depth 0 = cell layer (copy no. 0); depth 1 = the well it sits in (copy no. = well index)
    const G4int well = touch->GetCopyNumber(1);
    const G4int slot = fDetector->GetVoxelSlot(well);
    if (slot < 0) return;

    // The step is binned by its mid-point (steps in the thin layer are much shorter than a voxel)
    const G4ThreeVector mid = 0.5 * (pre->GetPosition() + step->GetPostStepPoint()->GetPosition());
    const G4ThreeVector local = mid - fDetector->GetWellCentre(well);   // no rotations in the geometry
    const G4int voxel = fDetector->GetVoxelIndex(local.x(), local.y());
    if (voxel < 0) return;

    // Kinetic energy averaged over the step
    const G4double ekin = 0.5 * (pre->GetKineticEnergy() + step->GetPostStepPoint()->GetKineticEnergy());
    fEventAction->AddTrack(slot * fDetector->GetVoxelsPerWell() + voxel, length, ekin);
}
