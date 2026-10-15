#include "PMRunAction.hh"
#include "PMRun.hh"
#include "PMDetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4UnitsTable.hh"
#include "G4Run.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>

PMRunAction::PMRunAction(const PMDetectorConstruction* det) : fDetector(det) {}

G4Run* PMRunAction::GenerateRun()
{
    return new PMRun(PMDetectorConstruction::kWells);
}

void PMRunAction::EndOfRunAction(const G4Run* run)
{
    // Only the master (or the sequential run manager) reports merged results
    if (!IsMaster()) return;

    const G4int nEvents = run->GetNumberOfEvent();
    if (nEvents == 0) return;

    const auto* pmRun = static_cast<const PMRun*>(run);
    const auto& sum = pmRun->GetSum();
    const auto& sum2 = pmRun->GetSum2();

    const G4double mass = fDetector->GetLiquidMass();   // mass of liquid in ONE well
    const char rowName[] = "ABCDEF";

    std::ofstream csv("well_dose.csv");
    csv << "well,row,col,edep_per_decay_MeV,dose_per_decay_Gy,rel_error,"
           "dose_rate_Gy_per_h_at_1MBq\n";

    G4cout << "\n--------------------- 48-well plate dose summary ---------------------\n"
           << " Sr-90 decays simulated : " << nEvents << "\n"
           << " Liquid mass per well   : " << G4BestUnit(mass, "Mass") << "\n"
           << " Dose rate column assumes " << kActivity_Bq << " Bq of Sr-90 "
           << "(Y-90 in equilibrium)\n\n"
           << " Dose rate [mGy/h] at that activity (rows A-F, columns 1-8):\n\n      ";
    for (G4int c = 0; c < PMDetectorConstruction::kCols; ++c)
        G4cout << std::setw(10) << c + 1;
    G4cout << "\n";

    for (G4int r = 0; r < PMDetectorConstruction::kRows; ++r) {
        G4cout << "   " << rowName[r] << "  ";
        for (G4int c = 0; c < PMDetectorConstruction::kCols; ++c) {
            const G4int i = r * PMDetectorConstruction::kCols + c;

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
           << "\n Per-well results written to well_dose.csv\n"
           << "-----------------------------------------------------------------------\n"
           << G4endl;
}
