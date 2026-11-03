#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"
#include "G4Exception.hh"

#include <algorithm>
#include <string>
#include <vector>

namespace {
struct MassFraction {
    const char* symbol;
    G4double percent;
};

// Builds a material from mass percentages. If the percentages do not add up to 100 they are
// normalised (with a warning), so the density you give is the density you get.
G4Material* MakeMaterial(const std::string& name, G4double density, const std::vector<MassFraction>& comp)
{
    auto* nist = G4NistManager::Instance();
    G4double total = 0.;
    for (const auto& c : comp) total += c.percent;
    if (std::abs(total - 100.) > 0.01) {
        G4ExceptionDescription ed;
        ed << "Mass percentages of " << name << " add up to " << total
           << " %, not 100 %. They are normalised to 100 % (each fraction divided by " << total / 100. << ").";
        G4Exception("DetectorConstruction::Construct", "Mat001", JustWarning, ed);
    }
    auto* mat = new G4Material(name, density, static_cast<G4int>(comp.size()));
    for (const auto& c : comp)
        mat->AddElement(nist->FindOrBuildElement(c.symbol), c.percent / total);
    return mat;
}
}  // namespace

G4int DetectorConstruction::GetVoxelIndex(G4double x, G4double y) const
{
    const G4int nx = GetNumVoxelsX(), ny = GetNumVoxelsY();
    const G4double lx = x - fFilmCentreX + fFilmX / 2;   // 0 .. fFilmX across the film
    const G4double ly = y - fFilmCentreY + fFilmY / 2;
    G4int ix = static_cast<G4int>(std::floor(lx / fVoxelSize));
    G4int iy = static_cast<G4int>(std::floor(ly / fVoxelSize));
    // A point exactly on the far edge belongs to the last voxel
    if (ix == nx && lx <= fFilmX + 1e-9) ix = nx - 1;
    if (iy == ny && ly <= fFilmY + 1e-9) iy = ny - 1;
    if (ix < 0 || ix >= nx || iy < 0 || iy >= ny) return -1;
    return iy * nx + ix;
}

void DetectorConstruction::GetVoxelCentre(G4int index, G4double& xRel, G4double& yRel) const
{
    const G4int nx = GetNumVoxelsX();
    xRel = -fFilmX / 2 + ((index % nx) + 0.5) * fVoxelSize;
    yRel = -fFilmY / 2 + ((index / nx) + 0.5) * fVoxelSize;
}

G4double DetectorConstruction::GetVoxelMass() const
{
    if (!fActiveLV) return 0.;
    return GetVoxelArea() * fActiveThickness * fActiveLV->GetMaterial()->GetDensity();
}

G4double DetectorConstruction::GetActiveLayerMass() const
{
    return fActiveLV ? fActiveLV->GetMass() : 0.;
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    // ---------- Materials ----------
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* silver = nist->FindOrBuildMaterial("G4_Ag");
    G4Material* strontium = nist->FindOrBuildMaterial("G4_Sr");

    // Radiochromic film. Mass percentages; the active layer's add up to 101.96 %, so
    // MakeMaterial normalises them (and warns).
    // Maryam's numbers
    //G4Material* activeMat = MakeMaterial("RCF_Active", 0.95 * g / cm3,
    //    {{"H", 10.8}, {"Li", 0.65}, {"C", 52.98}, {"N", 0.86}, {"O", 29.0},
    //     {"Na", 0.36}, {"Al", 6.27}, {"S", 0.49}, {"Cl", 0.55}});
    //G4Material* baseMat = MakeMaterial("RCF_Base", 1.35 * g / cm3,
    //    {{"H", 4.2}, {"C", 62.5}, {"O", 33.3}});
    // Palmer et al (with 0.1% tweaks to make 100%)
    G4Material* activeMat = MakeMaterial("RCF_Active", 1.35 * g / cm3,
        {{"H", 57.1}, {"Li", 0.6}, {"C", 28.5}, {"N", 0.4}, {"O", 11.7},
         {"Na", 0.1}, {"Al", 1.4}, {"S", 0.1}, {"Cl", 0.1}});
    G4Material* baseMat = MakeMaterial("RCF_Base", 1.35 * g / cm3,
            {{"H", 36.4}, {"C", 45.5}, {"O", 18.2}});

    // ---------- Checks ----------
    const G4int nx = GetNumVoxelsX(), ny = GetNumVoxelsY();
    if (fVoxelSize <= 0. || std::abs(nx * fVoxelSize - fFilmX) > 1e-6 * mm ||
        std::abs(ny * fVoxelSize - fFilmY) > 1e-6 * mm) {
        G4Exception("DetectorConstruction::Construct", "Vox001", FatalException,
                    "fFilmX and fFilmY must be whole multiples of fVoxelSize.");
    }
    if (fGap < 0.) {
        G4Exception("DetectorConstruction::Construct", "Src001", FatalException,
                    "fGap must be >= 0 (the source would overlap the film).");
    }
    if (std::hypot(fSourceX - fFilmCentreX, fSourceY - fFilmCentreY) + fSourceRadius >
        0.5 * std::min(fFilmX, fFilmY)) {
        G4Exception("DetectorConstruction::Construct", "Src002", JustWarning,
                    "The active disc of the source extends beyond the film.");
    }

    // ---------- World ----------
    G4double worldX = 100. * mm, worldY = 100. * mm, worldZ = 60. * mm;  // half-lengths
    auto* solidWorld = new G4Box("World", worldX, worldY, worldZ);
    auto* logicWorld = new G4LogicalVolume(solidWorld, air, "World");
    auto* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld,
                                        "World", nullptr, false, 0, true);
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ---------- Film: top base layer, active layer, bottom base layer; top of the film at z = 0 ----------
    const G4double zTopBase = -fTopBaseThickness / 2;
    const G4double zActive = -fTopBaseThickness - fActiveThickness / 2;
    const G4double zBottomBase = -fTopBaseThickness - fActiveThickness - fBaseThickness / 2;
    const G4VisAttributes baseVis(G4Colour(0.5, 0.7, 0.95, 1.0));

    auto* solidTopBase = new G4Box("BaseTop", fFilmX / 2, fFilmY / 2, fTopBaseThickness / 2);
    auto* logicTopBase = new G4LogicalVolume(solidTopBase, baseMat, "BaseTop");
    logicTopBase->SetVisAttributes(baseVis);
    new G4PVPlacement(nullptr, G4ThreeVector(fFilmCentreX, fFilmCentreY, zTopBase),
                      logicTopBase, "BaseTop", logicWorld, false, 0, true);

    auto* solidActive = new G4Box("Active", fFilmX / 2, fFilmY / 2, fActiveThickness / 2);
    fActiveLV = new G4LogicalVolume(solidActive, activeMat, "Active");
    fActiveLV->SetVisAttributes(G4VisAttributes(G4Colour(0.95, 0.6, 0.15, 1.0)));
    new G4PVPlacement(nullptr, G4ThreeVector(fFilmCentreX, fFilmCentreY, zActive),
                      fActiveLV, "Active", logicWorld, false, 0, true);

    auto* solidBottomBase = new G4Box("BaseBottom", fFilmX / 2, fFilmY / 2, fBaseThickness / 2);
    auto* logicBottomBase = new G4LogicalVolume(solidBottomBase, baseMat, "BaseBottom");
    logicBottomBase->SetVisAttributes(baseVis);
    new G4PVPlacement(nullptr, G4ThreeVector(fFilmCentreX, fFilmCentreY, zBottomBase),
                      logicBottomBase, "BaseBottom", logicWorld, false, 0, true);

    // ---------- Sr-90 applicator above the film, window facing down ----------
    auto* solidHousing = new G4Tubs("Housing", 0., fHousingRadius, fHousingHeight / 2, 0., 360. * deg);
    auto* logicHousing = new G4LogicalVolume(solidHousing, silver, "Housing");
    logicHousing->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.6, 0.5)));

    // Bottom face of the housing is fGap above the top of the film (z = 0)
    G4double zHousing = fGap + fHousingHeight / 2;
    new G4PVPlacement(nullptr, G4ThreeVector(fSourceX, fSourceY, zHousing), logicHousing,
                      "Housing", logicWorld, false, 0, true);

    auto* solidSource = new G4Tubs("Source", 0., fSourceRadius, fSourceThickness / 2, 0., 360. * deg);
    auto* logicSource = new G4LogicalVolume(solidSource, strontium, "Source");
    logicSource->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.8, 0.0)));

    // Active disc embedded in the housing, `fWindow` above the bottom face of the housing
    G4double zSourceLocal = -fHousingHeight / 2 + fWindow + fSourceThickness / 2;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zSourceLocal), logicSource,
                      "Source", logicHousing, false, 0, true);

    fSourceCentre = G4ThreeVector(fSourceX, fSourceY, zHousing + zSourceLocal);

    return physWorld;
}
