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

#include <cmath>

SteppingAction::SteppingAction(const DetectorConstruction* det, EventAction* evt)
    : fDetector(det), fEventAction(evt) {}

void SteppingAction::UserSteppingAction(const G4Step* step)
{
    const G4StepPoint* pre = step->GetPreStepPoint();

    // ---- Electron spectrum: electrons crossing the spectrum plane, moving away from the source ----
    {
        const G4StepPoint* post = step->GetPostStepPoint();
        const G4double z0 = fDetector->GetSpectrumPlaneZ();
        const G4int dir = fDetector->GetSpectrumDirection();
        // The step starts strictly before the plane and ends on it or beyond it. A particle that
        // stops exactly on the plane is counted by the step that reaches it, never twice.
        const G4double zPre = (pre->GetPosition().z() - z0) * dir;
        const G4double zPost = (post->GetPosition().z() - z0) * dir;
        if (zPre < 0. && zPost >= 0. && step->GetTrack()->GetDefinition() == G4Electron::Definition()) {
            const G4double f = zPre / (zPre - zPost);   // fraction of the step before the plane
            const G4ThreeVector pos = pre->GetPosition() + f * (post->GetPosition() - pre->GetPosition());
            const G4ThreeVector axis = fDetector->GetSourceCentre();
            if (std::hypot(pos.x() - axis.x(), pos.y() - axis.y()) <= fDetector->GetSpectrumRadius()) {
                // Kinetic energy at the plane, interpolated along the step
                const G4double ekin = pre->GetKineticEnergy() + f * (post->GetKineticEnergy() - pre->GetKineticEnergy());
                const G4int nBins = fDetector->GetNumSpectrumBins();
                G4int bin = static_cast<G4int>(std::floor(ekin / fDetector->GetSpectrumBinWidth()));
                if (bin >= nBins) bin = nBins;   // overflow bin
                fEventAction->AddSpectrumCount(bin);
            }
        }
    }

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

    // Voxel scoring: bin the deposit in the cell layer of the chosen wells by (x, y) position.
    // The mid-point of the step is used; steps in the thin layer are much shorter than a voxel.
    if (fDetector->VoxelScoring() && inCellLayer) {
        const G4int slot = fDetector->GetVoxelSlot(well);
        if (slot >= 0) {
            const G4ThreeVector mid = 0.5 * (step->GetPreStepPoint()->GetPosition()
                                             + step->GetPostStepPoint()->GetPosition());
            const G4ThreeVector local = mid - fDetector->GetWellCentre(well);   // no rotations in the geometry
            const G4int voxel = fDetector->GetVoxelIndex(local.x(), local.y());
            if (voxel >= 0)
                fEventAction->AddVoxelEdep(slot * fDetector->GetVoxelsPerWell() + voxel, edep);
        }
    }
}
