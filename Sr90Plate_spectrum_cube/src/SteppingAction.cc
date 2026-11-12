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

    // ---- Electron spectrum: electrons entering the cube from outside ----
    {
        const G4StepPoint* post = step->GetPostStepPoint();
        // At a volume boundary the post-step point already belongs to the volume being entered
        if (post->GetStepStatus() == fGeomBoundary && post->GetPhysicalVolume() &&
            fDetector->IsCube(post->GetPhysicalVolume()->GetLogicalVolume()) &&
            !fDetector->IsCube(pre->GetPhysicalVolume()->GetLogicalVolume()) &&
            step->GetTrack()->GetDefinition() == G4Electron::Definition() &&
            (!fDetector->CubeBottomFaceOnly() ||
             std::abs(post->GetPosition().z() - fDetector->GetCubeBottomZ()) < 1.0e-5 * mm)) {
            const G4int nBins = fDetector->GetNumSpectrumBins();
            G4int bin = static_cast<G4int>(std::floor(post->GetKineticEnergy() / fDetector->GetSpectrumBinWidth()));
            if (bin >= nBins) bin = nBins;   // overflow bin
            fEventAction->AddSpectrumCount(bin);
        }
    }

    const G4double edep = step->GetTotalEnergyDeposit();
    if (edep <= 0.) return;

    const G4TouchableHandle touch = pre->GetTouchableHandle();
    const G4LogicalVolume* lv = touch->GetVolume()->GetLogicalVolume();

    // The cube is made of the same material as the liquid, so its deposits are part of the liquid of its well
    if (fDetector->IsCube(lv)) {
        fEventAction->AddEdep(fDetector->GetCubeWellIndex(), edep);
        return;
    }
    if (!fDetector->IsScoringVolume(lv)) return;

    // Depth 0 = the liquid; depth 1 = the well it sits in (copy no. = well index)
    fEventAction->AddEdep(touch->GetCopyNumber(1), edep);
}
