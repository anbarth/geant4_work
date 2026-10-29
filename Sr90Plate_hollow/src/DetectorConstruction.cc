#include "DetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4SubtractionSolid.hh"
#include "G4Exception.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <string>

namespace {
// Area of the overlap between the rectangle [x0,x1] x [y0,y1] and a disc of radius R centred on
// the origin: integrate the length of the chord inside the rectangle over x (midpoint rule).
G4double RectDiscOverlap(G4double x0, G4double x1, G4double y0, G4double y1, G4double R)
{
    const G4double a = std::max(x0, -R);
    const G4double b = std::min(x1, R);
    if (b <= a) return 0.;
    const G4int n = 2000;
    const G4double dx = (b - a) / n;
    G4double area = 0.;
    for (G4int k = 0; k < n; ++k) {
        const G4double x = a + (k + 0.5) * dx;
        const G4double h = std::sqrt(std::max(0., R * R - x * x));   // half chord at this x
        area += std::max(0., std::min(y1, h) - std::max(y0, -h)) * dx;
    }
    return area;
}
}  // namespace

DetectorConstruction::DetectorConstruction()
{
    if (fVoxelScoring) SetupVoxels();
}

void DetectorConstruction::SetupVoxels()
{
    if (fVoxelSize <= 0.) {
        G4Exception("DetectorConstruction::SetupVoxels", "Vox002", FatalException,
                    "fVoxelSize must be positive.");
    }

    // Wells to voxelise: names like "C5" (row letter A-F, column 1-8), or "ALL"
    fVoxelSlot.assign(kWells, -1);
    fVoxelWellIndex.clear();
    auto addWell = [&](G4int well) {
        if (fVoxelSlot[well] >= 0) return;   // listed twice
        fVoxelSlot[well] = static_cast<G4int>(fVoxelWellIndex.size());
        fVoxelWellIndex.push_back(well);
    };
    for (const std::string& raw : fVoxelWells) {
        std::string name;
        for (char ch : raw) {
            if (!std::isspace(static_cast<unsigned char>(ch)))
                name += static_cast<char>(std::toupper(static_cast<unsigned char>(ch)));
        }
        if (name == "ALL") {
            for (G4int w = 0; w < kWells; ++w) addWell(w);
            continue;
        }
        bool ok = name.size() >= 2 && name.size() <= 3 && name[0] >= 'A' && name[0] < 'A' + kRows;
        G4int col = 0;
        for (std::size_t k = 1; ok && k < name.size(); ++k)
            ok = std::isdigit(static_cast<unsigned char>(name[k])) != 0;
        if (ok) {
            col = std::stoi(name.substr(1));
            ok = col >= 1 && col <= kCols;
        }
        if (!ok) {
            const std::string msg = "Bad well name '" + raw + "' in fVoxelWells: use e.g. \"C5\" "
                                    "(row A-F, column 1-8) or \"ALL\".";
            G4Exception("DetectorConstruction::SetupVoxels", "Vox001", FatalException, msg.c_str());
        }
        addWell((name[0] - 'A') * kCols + (col - 1));
    }
    if (fVoxelWellIndex.empty()) {
        G4Exception("DetectorConstruction::SetupVoxels", "Vox003", FatalException,
                    "fVoxelScoring is on but fVoxelWells is empty.");
    }

    // Square grid centred on the well axis, just big enough to cover the well
    fNVox = static_cast<G4int>(std::ceil(fWellDiameter / fVoxelSize - 1e-9));
    const G4double half = fNVox * fVoxelSize / 2;
    const G4double R = fWellDiameter / 2;

    fVoxelArea.assign(fNVox * fNVox, 0.);
    for (G4int iy = 0; iy < fNVox; ++iy) {
        for (G4int ix = 0; ix < fNVox; ++ix) {
            const G4double x0 = -half + ix * fVoxelSize;
            const G4double y0 = -half + iy * fVoxelSize;
            fVoxelArea[iy * fNVox + ix] = RectDiscOverlap(x0, x0 + fVoxelSize, y0, y0 + fVoxelSize, R);
        }
    }
}

G4ThreeVector DetectorConstruction::GetWellCentre(G4int well) const
{
    const G4int row = well / kCols;
    const G4int col = well % kCols;
    return G4ThreeVector((col - (kCols - 1) / 2.0) * fWellPitch,
                         ((kRows - 1) / 2.0 - row) * fWellPitch, 0.);
}

G4int DetectorConstruction::GetVoxelIndex(G4double xLocal, G4double yLocal) const
{
    const G4double half = fNVox * fVoxelSize / 2;
    const G4int ix = static_cast<G4int>(std::floor((xLocal + half) / fVoxelSize));
    const G4int iy = static_cast<G4int>(std::floor((yLocal + half) / fVoxelSize));
    if (ix < 0 || ix >= fNVox || iy < 0 || iy >= fNVox) return -1;
    return iy * fNVox + ix;
}

void DetectorConstruction::GetVoxelCentre(G4int index, G4double& xLocal, G4double& yLocal) const
{
    const G4double half = fNVox * fVoxelSize / 2;
    xLocal = -half + ((index % fNVox) + 0.5) * fVoxelSize;
    yLocal = -half + ((index / fNVox) + 0.5) * fVoxelSize;
}

G4double DetectorConstruction::GetVoxelMass(G4int index) const
{
    if (!fCellLayerLV) return 0.;
    return fVoxelArea[index] * fCellLayerThickness * fCellLayerLV->GetMaterial()->GetDensity();
}

G4double DetectorConstruction::GetScoredMass() const
{
    if (!fCellLayerLV || !fLiquidLV) return 0.;
    const G4double layerMass = fCellLayerLV->GetMass();
    return fScoreCellLayerOnly ? layerMass : layerMass + fLiquidLV->GetMass();
}

G4VPhysicalVolume* DetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    // ---------- Materials ----------
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* polystyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
    G4Material* medium = nist->FindOrBuildMaterial("G4_WATER");   // culture medium ~ water
    G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");
    G4Material* silver = nist->FindOrBuildMaterial("G4_Ag");

    // Active layer: strontium titanate (SrTiO3), a common Sr-90 carrier
    G4Element* elSr = nist->FindOrBuildElement("Sr");
    G4Element* elTi = nist->FindOrBuildElement("Ti");
    G4Element* elO = nist->FindOrBuildElement("O");
    auto* srtio3 = new G4Material("SrTiO3", 5.12 * g / cm3, 3);
    srtio3->AddElement(elSr, 1);
    srtio3->AddElement(elTi, 1);
    srtio3->AddElement(elO, 3);

    // Alternative active layer: pure Sr
    G4Material* strontium = nist->FindOrBuildMaterial("G4_Sr");

    // ---------- World ----------
    G4double worldX = 100. * mm, worldY = 100. * mm, worldZ = 60. * mm;  // half-lengths
    auto* solidWorld = new G4Box("World", worldX, worldY, worldZ);
    auto* logicWorld = new G4LogicalVolume(solidWorld, air, "World");
    auto* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld,
                                        "World", nullptr, false, 0, true);
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ---------- Plate: a rectangular bottom sheet plus 48 separate hollow tubes standing on it ----------
    // The space between the tubes is just air (the world material).
    const G4double wellDepth = fPlateZ - fFloor;          // height of each tube
    const G4double wellR = fWellDiameter / 2;
    const G4double tubeOuterR = wellR + fWallThickness;
    if (fWallThickness <= 0. || 2 * tubeOuterR > fWellPitch) {
        G4Exception("DetectorConstruction::Construct", "Plate001", FatalException,
                    "fWallThickness must be > 0 and fWellDiameter + 2*fWallThickness must not exceed fWellPitch.");
    }

    auto* solidBase = new G4Box("PlateBase", fPlateX / 2, fPlateY / 2, fFloor / 2);
    auto* logicBase = new G4LogicalVolume(solidBase, polystyrene, "PlateBase");
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -fPlateZ / 2 + fFloor / 2), logicBase,
                      "PlateBase", logicWorld, false, 0, true);
    logicBase->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.8, 0.5)));

    auto* solidTube = new G4Tubs("WellTube", wellR, tubeOuterR, wellDepth / 2, 0., 360. * deg);
    auto* logicTube = new G4LogicalVolume(solidTube, polystyrene, "WellTube");
    logicTube->SetVisAttributes(G4VisAttributes(G4Colour(0.6, 0.6, 0.8, 0.5)));

    // ---------- One well: air cylinder filling the inside of a tube + liquid at the bottom ----------
    auto* solidWell = new G4Tubs("Well", 0., wellR, wellDepth / 2, 0., 360. * deg);
    auto* logicWell = new G4LogicalVolume(solidWell, air, "Well");
    // Drawn as a wireframe so the 48 well outlines are visible
    G4VisAttributes wellVis(G4Colour(0.2, 0.2, 0.2, 1.0));
    wellVis.SetForceWireframe(true);
    logicWell->SetVisAttributes(wellVis);

    if (fCellLayerThickness <= 0. || fCellLayerThickness >= fFillHeight) {
        G4Exception("DetectorConstruction::Construct", "Cell001", FatalException,
                    "fCellLayerThickness must be > 0 and smaller than fFillHeight.");
    }

    // The liquid is split in two stacked cylinders of the same material: a thin cell layer
    // on the floor of the well, and the rest of the liquid above it. Splitting does not
    // change the physics; it only lets us score the layer on its own.
    const G4double liquidHeight = fFillHeight - fCellLayerThickness;

    auto* solidCell = new G4Tubs("CellLayer", 0., fWellDiameter / 2, fCellLayerThickness / 2, 0., 360. * deg);
    fCellLayerLV = new G4LogicalVolume(solidCell, medium, "CellLayer");
    fCellLayerLV->SetVisAttributes(G4VisAttributes(G4Colour(0.2, 0.8, 0.3, 1.0)));

    auto* solidLiquid = new G4Tubs("Liquid", 0., fWellDiameter / 2, liquidHeight / 2, 0., 360. * deg);
    fLiquidLV = new G4LogicalVolume(solidLiquid, medium, "Liquid");
    fLiquidLV->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.3, 0.3, 1.0)));

    // Cell layer sits on the floor of the well cavity (local z of the well cavity),
    // with the rest of the liquid stacked on top of it
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -wellDepth / 2 + fCellLayerThickness / 2),
                      fCellLayerLV, "CellLayer", logicWell, false, 0, true);
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -wellDepth / 2 + fCellLayerThickness + liquidHeight / 2),
                      fLiquidLV, "Liquid", logicWell, false, 0, true);

    // ---------- Place the 48 wells ----------
    G4double zWell = -fPlateZ / 2 + fFloor + wellDepth / 2;
    for (G4int row = 0; row < kRows; ++row) {
        for (G4int col = 0; col < kCols; ++col) {
            G4double x = (col - (kCols - 1) / 2.0) * fWellPitch;
            G4double y = ((kRows - 1) / 2.0 - row) * fWellPitch;
            G4int index = row * kCols + col;
            new G4PVPlacement(nullptr, G4ThreeVector(x, y, zWell), logicWell, "Well",
                              logicWorld, false, index, true);
            new G4PVPlacement(nullptr, G4ThreeVector(x, y, zWell), logicTube, "WellTube",
                              logicWorld, false, index, true);
        }
    }

    // ---------- Lid: flat top above the plate, with a skirt around the outside ----------
    // One solid: a box (top + skirt) with a rectangular cavity cut out from below.
    // The air between lid and plate is just the world material.
    if (fHasLid) {
        if (fGap < fLidAirGap + fLidThickness) {
            G4Exception("DetectorConstruction::Construct", "Lid001", FatalException,
                        "fGap is smaller than the height of the lid top above the plate: "
                        "the source housing would overlap the lid.");
        }

        const G4double eps = 1.0 * mm;   // cavity sticks out below the skirt: no coincident faces
        const G4double innerX = fPlateX + 2 * fLidClearance;
        const G4double innerY = fPlateY + 2 * fLidClearance;
        const G4double outerX = innerX + 2 * fLidWall;
        const G4double outerY = innerY + 2 * fLidWall;
        const G4double lidHeight = fLidThickness + fLidSkirtHeight;

        auto* solidLidOuter = new G4Box("LidOuter", outerX / 2, outerY / 2, lidHeight / 2);
        auto* solidLidCavity = new G4Box("LidCavity", innerX / 2, innerY / 2,
                                         (fLidSkirtHeight + eps) / 2);
        // Cavity top is flush with the underside of the flat top (lid-local z = (skirt - thickness)/2)
        auto* solidLid = new G4SubtractionSolid("Lid", solidLidOuter, solidLidCavity, nullptr,
                                                G4ThreeVector(0, 0, -(fLidThickness + eps) / 2));
        auto* logicLid = new G4LogicalVolume(solidLid, polystyrene, "Lid");
        G4VisAttributes lidVis(G4Colour(0.7, 0.85, 1.0, 0.25));
        logicLid->SetVisAttributes(lidVis);

        // Underside of the flat top sits fLidAirGap above the top of the plate
        G4double zLid = fPlateZ / 2 + fLidAirGap - (fLidSkirtHeight - fLidThickness) / 2;
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zLid), logicLid, "Lid", logicWorld, false, 0, true);
    }

    // ---------- Sr-90 applicator above the plate, window facing down ----------
    auto* solidHousing = new G4Tubs("Housing", 0., fHousingRadius, fHousingHeight / 2, 0., 360. * deg);
    //auto* logicHousing = new G4LogicalVolume(solidHousing, steel, "Housing");
    auto* logicHousing = new G4LogicalVolume(solidHousing, silver, "Housing");
    G4VisAttributes housingVis(G4Colour(0.6, 0.6, 0.6, 0.5));
    //housingVis.SetForceWireframe(true);
    logicHousing->SetVisAttributes(housingVis);

    // G4double zHousing = fPlateZ / 2 + fGap + fHousingHeight / 2;
    G4double zHousing = -fPlateZ / 2 - fHousingHeight / 2;
    new G4PVPlacement(nullptr, G4ThreeVector(fSourceX, fSourceY, zHousing), logicHousing,
                      "Housing", logicWorld, false, 0, true);

    auto* solidSource = new G4Tubs("Source", 0., fSourceRadius, fSourceThickness / 2, 0., 360. * deg);
    //auto* logicSource = new G4LogicalVolume(solidSource, srtio3, "Source");
    auto* logicSource = new G4LogicalVolume(solidSource, strontium, "Source");
    logicSource->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.8, 0.0)));

    // Active disc embedded in the steel, `fWindow` above the bottom face of the housing
    G4double zSourceLocal = -fHousingHeight / 2 + fWindow + fSourceThickness / 2;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zSourceLocal), logicSource,
                      "Source", logicHousing, false, 0, true);

    fSourceCentre = G4ThreeVector(fSourceX, fSourceY, zHousing + zSourceLocal);

    return physWorld;
}
