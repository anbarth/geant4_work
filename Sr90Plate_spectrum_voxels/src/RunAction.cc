#include "RunAction.hh"
#include "Run.hh"
#include "DetectorConstruction.hh"

#include "G4SystemOfUnits.hh"
#include "G4Run.hh"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <string>

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
    const auto& sL = myRun->GetL();
    const auto& sS = myRun->GetS();
    const auto& sLL = myRun->GetLL();
    const auto& sSS = myRun->GetSS();
    const auto& sLS = myRun->GetLS();

    const G4int nVox = fDetector->GetNumVoxels();
    const G4int perWell = fDetector->GetVoxelsPerWell();
    if (nVox == 0 || static_cast<G4int>(sL.size()) != nVox) return;

    const char rowName[] = "ABCDEF";
    auto wellName = [&](G4int well) {
        return std::string(1, rowName[well / DetectorConstruction::kCols])
               + std::to_string(well % DetectorConstruction::kCols + 1);
    };

    // Mean energy of the electrons in a voxel = track-length-weighted mean (i.e. the mean of the fluence
    // spectrum):  R = sum(E*l) / sum(l)  over all electron steps in the voxel, all events.
    // Event totals S_e = sum(E*l), L_e = sum(l) are the independent samples. With N events (zeros included)
    //   m_S = sum S_e / N,  m_L = sum L_e / N,  R = m_S / m_L
    //   var(R) = ( s_SS + R^2 s_LL - 2 R s_LS ) / ( N m_L^2 )      (ratio estimator, s = sample covariances)
    // Statistical error only.
    const G4double N = nEvents;

    std::ofstream csv("voxel_energy.csv");
    csv << "well,ix,iy,x_mm,y_mm,area_mm2,mean_energy_keV,error_keV,track_length_per_decay_mm,"
           "fluence_per_decay_per_mm2\n"
        << std::setprecision(6) << std::scientific;

    G4int nWritten = 0, nEmpty = 0;
    G4double eMin = 1e30, eMax = 0.;
    for (G4int i = 0; i < nVox; ++i) {
        const G4int slot = i / perWell;
        const G4int v = i % perWell;
        const G4double area = fDetector->GetVoxelArea(v);
        if (area <= 0.) continue;   // voxel completely outside the well

        const G4double mL = sL[i] / N;
        const G4double mS = sS[i] / N;
        G4double meanE = 0., err = 0.;
        if (sL[i] > 0.) {
            meanE = mS / mL;
            const G4double sSSv = sSS[i] / N - mS * mS;
            const G4double sLLv = sLL[i] / N - mL * mL;
            const G4double sLSv = sLS[i] / N - mS * mL;
            const G4double var = (sSSv + meanE * meanE * sLLv - 2. * meanE * sLSv) / (N * mL * mL);
            err = std::sqrt(std::max(0., var));
            eMin = std::min(eMin, meanE);
            eMax = std::max(eMax, meanE);
        } else {
            ++nEmpty;
        }

        const G4double trackPerDecay = mL;   // total electron track length per decay in the voxel
        const G4double fluence = trackPerDecay / (area * fDetector->GetCellLayerThickness());

        const G4int well = fDetector->GetVoxelWellIndex(slot);
        G4double x, y;
        fDetector->GetVoxelCentre(v, x, y);
        csv << wellName(well) << ","
            << v % fDetector->GetVoxelsPerSide() << "," << v / fDetector->GetVoxelsPerSide() << ","
            << x / mm << "," << y / mm << "," << area / mm2 << ","
            << meanE / keV << "," << err / keV << "," << trackPerDecay / mm << ","
            << fluence * mm2 << "\n";
        ++nWritten;
    }

    G4cout << "\n----------------- Mean electron energy per voxel -----------------\n"
           << " Sr-90 decays simulated : " << nEvents << "\n"
           << " Voxels written         : " << nWritten << " (" << fDetector->GetNumVoxelWells()
           << " well(s), " << fDetector->GetVoxelSize() / mm << " x " << fDetector->GetVoxelSize() / mm
           << " x " << fDetector->GetCellLayerThickness() / mm << " mm)\n"
           << " Voxels with no electron: " << nEmpty << " (mean energy written as 0)\n";
    if (nWritten > nEmpty)
        G4cout << std::setprecision(4) << " Mean energy range      : " << eMin / keV << " - " << eMax / keV
               << " keV\n";
    G4cout << " Written to voxel_energy.csv\n"
           << "------------------------------------------------------------------\n" << G4endl;
}
