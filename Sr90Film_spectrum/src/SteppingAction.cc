#include "SteppingAction.hh"
#include "EventAction.hh"
#include "DetectorConstruction.hh"

#include "G4Step.hh"
#include "G4StepPoint.hh"
#include "G4LogicalVolume.hh"
#include "G4VPhysicalVolume.hh"
#include "G4ThreeVector.hh"
#include "G4Track.hh"
#include "G4Electron.hh"
#include "G4SystemOfUnits.hh"

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

    if (!fDetector->IsActiveLayer(pre->GetPhysicalVolume()->GetLogicalVolume())) return;

    fEventAction->AddLayerEdep(edep);

    // Bin the deposit by (x, y). The mid-point of the step is used; steps in a 25 um layer are far
    // shorter than a voxel. The film has no rotation, so global x, y are used directly.
    const G4ThreeVector mid = 0.5 * (pre->GetPosition() + step->GetPostStepPoint()->GetPosition());
    const G4int voxel = fDetector->GetVoxelIndex(mid.x(), mid.y());
    if (voxel >= 0) fEventAction->AddVoxelEdep(voxel, edep);
}
