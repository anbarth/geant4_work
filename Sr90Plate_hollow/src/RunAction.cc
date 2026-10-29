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
#include <string>

RunAction::RunAction(const DetectorConstruction* det) : fDetector(det) {}

G4Run* RunAction::GenerateRun()
{
    return new Run(DetectorConstruction::kWells, fDetector->GetNumVoxels());
}

void RunAction::WriteVoxelResults(const Run* run, G4int nEvents) const
{
    const auto& vsum = run->GetVoxelSum();
    const auto& vsum2 = run->GetVoxelSum2();
    const G4int nVox = fDetector->GetNumVoxels();               // all voxelised wells
    const G4int perWell = fDetector->GetVoxelsPerWell();
    if (nVox == 0 || static_cast<G4int>(vsum.size()) != nVox) return;

    const char rowName[] = "ABCDEF";
    auto wellName = [&](G4int well) {
        return std::string(1, rowName[well / DetectorConstruction::kCols])
               + std::to_string(well % DetectorConstruction::kCols + 1);
    };

    std::ofstream csv("voxel_dose.csv");
    csv << "well,ix,iy,x_mm,y_mm,area_mm2,edep_per_decay_MeV,dose_per_decay_Gy,rel_error,"
           "dose_rate_Gy_per_h_at_1MBq\n"
        << std::setprecision(6) << std::scientific;

    G4int nWritten = 0;
    G4double maxRate = 0.;   // Gy/h at 1 MBq
    G4int maxWell = -1;
    for (G4int i = 0; i < nVox; ++i) {
        const G4int slot = i / perWell;
        const G4int v = i % perWell;                                 // voxel index within its well
        const G4double mass = fDetector->GetVoxelMass(v);
        if (mass <= 0.) continue;   // voxel completely outside the well

        const G4double meanEdep = vsum[i] / nEvents;                          // per decay
        const G4double doseGy = (meanEdep / mass) / gray;                     // per decay
        const G4double rate = doseGy * kActivity_Bq * 3600.;                  // Gy/h

        G4double rel = 0.;
        if (vsum[i] > 0.) {
            const G4double var = std::max(0., vsum2[i] - vsum[i] * vsum[i] / nEvents);
            rel = std::sqrt(var) / vsum[i];
        }

        const G4int well = fDetector->GetVoxelWellIndex(slot);
        G4double x, y;
        fDetector->GetVoxelCentre(v, x, y);
        csv << wellName(well) << ","
            << v % fDetector->GetVoxelsPerSide() << "," << v / fDetector->GetVoxelsPerSide() << ","
            << x / mm << "," << y / mm << "," << fDetector->GetVoxelArea(v) / mm2 << ","
            << meanEdep / MeV << "," << doseGy << "," << rel << "," << rate << "\n";
        if (rate > maxRate) { maxRate = rate; maxWell = well; }
        ++nWritten;
    }

    G4cout << " Voxel scoring: " << fDetector->GetNumVoxelWells() << " well(s), " << nWritten
           << " voxels of " << fDetector->GetVoxelSize() / mm << " x " << fDetector->GetVoxelSize() / mm
           << " x " << fDetector->GetCellLayerThickness() / mm << " mm\n";
    if (maxWell >= 0)
        G4cout << " Highest voxel dose rate: " << std::setprecision(3) << std::scientific
               << maxRate * 1000. << " mGy/h at 1 MBq (well " << wellName(maxWell) << ")\n";
    G4cout << " Written to voxel_dose.csv\n" << std::defaultfloat;
}

void RunAction::EndOfRunAction(const G4Run* run)
{
    // Only the master (or the sequential run manager) reports merged results
    if (!IsMaster()) return;

    const G4int nEvents = run->GetNumberOfEvent();
    if (nEvents == 0) return;

    const auto* myRun = static_cast<const Run*>(run);
    const auto& sum = myRun->GetSum();
    const auto& sum2 = myRun->GetSum2();

    const G4double mass = fDetector->GetScoredMass();   // mass of the scored region in ONE well
    const char rowName[] = "ABCDEF";

    std::ofstream csv("well_dose.csv");
    csv << "well,row,col,edep_per_decay_MeV,dose_per_decay_Gy,rel_error,"
           "dose_rate_Gy_per_h_at_1MBq\n";

    G4cout << "\n--------------------- 48-well plate dose summary ---------------------\n"
           << " Sr-90 decays simulated : " << nEvents << "\n"
           << " Scored region         : ";
    if (fDetector->ScoreCellLayerOnly())
        G4cout << "cell layer only (" << G4BestUnit(fDetector->GetCellLayerThickness(), "Length")
               << " at the well bottom)";
    else
        G4cout << "whole liquid column";
    G4cout << "\n"
           << " Scored mass per well   : " << G4BestUnit(mass, "Mass") << "\n"
           << " Dose rate column assumes " << kActivity_Bq << " Bq of Sr-90 "
           << "(Y-90 in equilibrium)\n\n"
           << " Dose rate [mGy/h] at that activity (rows A-F, columns 1-8):\n\n      ";
    for (G4int c = 0; c < DetectorConstruction::kCols; ++c)
        G4cout << std::setw(10) << c + 1;
    G4cout << "\n";

    for (G4int r = 0; r < DetectorConstruction::kRows; ++r) {
        G4cout << "   " << rowName[r] << "  ";
        for (G4int c = 0; c < DetectorConstruction::kCols; ++c) {
            const G4int i = r * DetectorConstruction::kCols + c;

            const G4double meanEdep = sum[i] / nEvents;                       // per decay
            const G4double dosePerDecay = meanEdep / mass;                    // G4 units
            const G4double dosePerDecayGy = dosePerDecay / gray;
            const G4double rate = dosePerDecayGy * kActivity_Bq * 3600.;      // Gy/h
            // (note that k is decays/s, so the 3600 converts to hours)

            // relative error on the sum: sqrt(sum2 - sum^2/N) / sum
            G4double rel = 0.;
            if (sum[i] > 0.) {
                // var = standard error of the TOTAL energy deposited
                const G4double var = std::max(0., sum2[i] - sum[i] * sum[i] / nEvents);
                // rel = fractional error on sum(energy), which should be the same as fractional error on other quantities
                rel = std::sqrt(var) / sum[i];
            }

            G4cout << std::setw(10) << std::setprecision(3) << std::scientific
                   << rate * 1000.;
            csv << rowName[r] << c + 1 << "," << r << "," << c << ","
                << std::setprecision(6) << std::scientific
                << meanEdep / MeV << "," << dosePerDecayGy << "," << rel << ","
                << rate << "\n";
        }
        G4cout << "\n";
    }
    G4cout << std::defaultfloat
           << "\n Per-well results written to well_dose.csv\n";
    if (fDetector->VoxelScoring()) WriteVoxelResults(myRun, nEvents);
    G4cout << "-----------------------------------------------------------------------\n"
           << G4endl;
}
