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
    return new Run(DetectorConstruction::kWells, fDetector->GetNumSpectrumBins());
}

void RunAction::WriteSpectrum(const Run* run, G4int nEvents) const
{
    const auto& sum = run->GetSpectrumSum();
    const auto& sum2 = run->GetSpectrumSum2();
    const G4int nBins = fDetector->GetNumSpectrumBins();
    if (nBins <= 0 || static_cast<G4int>(sum.size()) != nBins + 2) return;
    const G4double w = fDetector->GetSpectrumBinWidth();

    // Relative error of a per-decay count: sqrt(sum2 - sum^2/N) / sum, N = all events (zeros included)
    auto relError = [nEvents](G4double s, G4double s2) {
        if (s <= 0.) return 0.;
        return std::sqrt(std::max(0., s2 - s * s / nEvents)) / s;
    };

    std::ofstream csv("spectrum.csv");
    csv << "e_low_MeV,e_high_MeV,e_mid_MeV,counts_per_decay,rel_error,counts_per_decay_per_MeV\n"
        << std::setprecision(6) << std::scientific;
    G4double meanE = 0., inRange = 0.;
    for (G4int i = 0; i < nBins; ++i) {
        const G4double perDecay = sum[i] / nEvents;
        const G4double mid = (i + 0.5) * w;
        csv << i * w / MeV << "," << (i + 1) * w / MeV << "," << mid / MeV << "," << perDecay << ","
            << relError(sum[i], sum2[i]) << "," << perDecay / (w / MeV) << "\n";
        meanE += perDecay * mid;
        inRange += perDecay;
    }
    const G4int total = nBins + 1, over = nBins;
    G4cout << " Electrons entering the " << fDetector->GetCubeSize() / mm << " mm cube in well "
           << fDetector->GetCubeWellName() << (fDetector->CubeBottomFaceOnly() ? " (bottom face only): " : " (all faces): ")
           << std::setprecision(4) << std::scientific
           << sum[total] / nEvents << " per decay (rel. error " << std::setprecision(2)
           << relError(sum[total], sum2[total]) << ")\n"
           << std::defaultfloat << std::setprecision(4);
    if (inRange > 0.)
        G4cout << " Mean kinetic energy (from the histogram): " << meanE / inRange / keV << " keV\n";
    if (sum[over] > 0.)
        G4cout << " Warning: " << sum[over] / nEvents << " electrons per decay are above "
               << nBins * w / MeV << " MeV (overflow bin, not in spectrum.csv)\n";
    G4cout << " Spectrum written to spectrum.csv\n";
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
           << " Scored region         : liquid column of each well\n"
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
    WriteSpectrum(myRun, nEvents);
    G4cout << "-----------------------------------------------------------------------\n"
           << G4endl;
}
