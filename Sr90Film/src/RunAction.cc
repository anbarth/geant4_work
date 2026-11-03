#include "RunAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Run.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

RunAction::RunAction(const DetectorConstruction* det) : fDetector(det) {}

G4Run* RunAction::GenerateRun()
{
    return new Run(fDetector->GetNumVoxels());
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    // Only the master (or the sequential run manager) reports merged results
    if (!IsMaster()) return;

    const G4int nEvents = run->GetNumberOfEvent();
    if (nEvents == 0) return;

    const auto* myRun = static_cast<const Run*>(run);
    const auto& vsum = myRun->GetVoxelSum();
    const auto& vsum2 = myRun->GetVoxelSum2();
    const G4int nVox = fDetector->GetNumVoxels();
    const G4int nx = fDetector->GetNumVoxelsX();
    const G4double voxMass = fDetector->GetVoxelMass();
    const G4double layerMass = fDetector->GetActiveLayerMass();

    // Relative error on a sum of per-event deposits: sqrt(sum2 - sum^2/N) / sum, N = all events
    auto relError = [nEvents](G4double sum, G4double sum2) {
        if (sum <= 0.) return 0.;
        return std::sqrt(std::max(0., sum2 - sum * sum / nEvents)) / sum;
    };

    std::ofstream csv("voxel_dose.csv");
    csv << "ix,iy,x_mm,y_mm,area_mm2,edep_per_decay_MeV,dose_per_decay_Gy,rel_error,"
           "dose_rate_Gy_per_h_at_1MBq\n"
        << std::setprecision(6) << std::scientific;

    G4double maxRate = 0., maxX = 0., maxY = 0.;
    for (G4int i = 0; i < nVox; ++i) {
        const G4double meanEdep = vsum[i] / nEvents;                         // per decay
        const G4double doseGy = (meanEdep / voxMass) / gray;                 // per decay
        const G4double rate = doseGy * kActivity_Bq * 3600.;                 // Gy/h at 1 MBq
        G4double x, y;
        fDetector->GetVoxelCentre(i, x, y);
        csv << i % nx << "," << i / nx << "," << x / mm << "," << y / mm << ","
            << fDetector->GetVoxelArea() / mm2 << ","
            << meanEdep / MeV << "," << doseGy << "," << relError(vsum[i], vsum2[i]) << "," << rate << "\n";
        if (rate > maxRate) { maxRate = rate; maxX = x; maxY = y; }
    }

    const G4double layerMean = myRun->GetLayerSum() / nEvents;               // per decay
    const G4double layerRate = (layerMean / layerMass) / gray * kActivity_Bq * 3600.;

    G4cout << "\n--------------------- radiochromic film dose summary ---------------------\n"
           << " Sr-90 decays simulated : " << nEvents << "\n"
           << " Active layer           : " << G4BestUnit(fDetector->GetActiveThickness(), "Length")
           << " thick, mass " << G4BestUnit(layerMass, "Mass") << "\n"
           << " Voxels                 : " << fDetector->GetNumVoxelsX() << " x " << fDetector->GetNumVoxelsY()
           << " of " << fDetector->GetVoxelSize() / mm << " x " << fDetector->GetVoxelSize() / mm
           << " mm (mass " << G4BestUnit(voxMass, "Mass") << " each)\n"
           << std::setprecision(4) << std::scientific
           << " Mean dose rate, whole layer : " << layerRate * 1000. << " mGy/h at " << kActivity_Bq
           << " Bq (rel. error " << std::setprecision(2) << relError(myRun->GetLayerSum(), myRun->GetLayerSum2())
           << ")\n"
           << std::setprecision(4)
           << " Highest voxel dose rate        : " << maxRate * 1000. << " mGy/h at voxel centre ("
           << maxX / mm << ", " << maxY / mm << ") mm from the film centre\n"
           << std::defaultfloat
           << " Voxel results written to voxel_dose.csv\n"
           << "-----------------------------------------------------------------------\n"
           << G4endl;
}
