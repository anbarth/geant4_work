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
#include <vector>

G4double DetectorConstruction::GetScoredMass() const
{
    // Mass of the liquid of ONE well (the cube has the same material, so the well that holds it is the same)
    return fLiquidLV ? fLiquidLV->GetMass() : 0.;
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

    auto* solidLiquid = new G4Tubs("Liquid", 0., wellR, fFillHeight / 2, 0., 360. * deg);
    fLiquidLV = new G4LogicalVolume(solidLiquid, medium, "Liquid");
    fLiquidLV->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.3, 0.3, 1.0)));

    // Liquid sits on the floor of the well cavity (local z of the well cavity)
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -wellDepth / 2 + fFillHeight / 2),
                      fLiquidLV, "Liquid", logicWell, false, 0, true);

    G4LogicalVolume* logicWellC5 = nullptr;   // the well that holds the cube

    // ---------- The well that holds the cube ----------
    // The cube must exist in one well only, but the Well and Liquid volumes are shared by all 48 wells,
    // so this well gets its own copies: a Well and a Liquid (same solid and material) with the cube inside it.
    {
        const std::string name = fCubeWell;
        if (name.size() < 2 || !std::isalpha(static_cast<unsigned char>(name[0])) ||
            !std::isdigit(static_cast<unsigned char>(name[1]))) {
            G4Exception("DetectorConstruction::Construct", "Cube001", FatalException,
                        "fCubeWell must look like C5 (rows A-F, columns 1-8).");
        }
        const G4int row = std::toupper(static_cast<unsigned char>(name[0])) - 'A';
        const G4int col = std::stoi(name.substr(1)) - 1;
        if (row < 0 || row >= kRows || col < 0 || col >= kCols) {
            G4Exception("DetectorConstruction::Construct", "Cube001", FatalException,
                        "fCubeWell must look like C5 (rows A-F, columns 1-8).");
        }
        fCubeWellIndex = row * kCols + col;

        const G4double half = fCubeSize / 2;
        if (fCubeBottomHeight < 0. || fCubeBottomHeight + fCubeSize > fFillHeight ||
            std::hypot(std::abs(fCubeX) + half, std::abs(fCubeY) + half) > wellR) {
            G4Exception("DetectorConstruction::Construct", "Cube002", FatalException,
                        "The cube must lie inside the liquid and inside the well wall: "
                        "check fCubeSize, fCubeX/Y and fCubeBottomHeight.");
        }

        // Cube centre in the well frame (z = 0 at the middle of the well height)
        const G4double zCubeWell = -wellDepth / 2 + fCubeBottomHeight + half;
        const G4double zLiquidWell = -wellDepth / 2 + fFillHeight / 2;
        fCubeBottomZ = -fPlateZ / 2 + fFloor + fCubeBottomHeight;   // world frame

        auto* solidCube = new G4Box("Cube", half, half, half);
        fCubeLV = new G4LogicalVolume(solidCube, medium, "Cube");
        fCubeLV->SetVisAttributes(G4VisAttributes(G4Colour(0.1, 0.3, 1.0, 1.0)));

        fLiquidCubeWellLV = new G4LogicalVolume(solidLiquid, medium, "LiquidCubeWell");
        fLiquidCubeWellLV->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.3, 0.3, 0.3)));
        new G4PVPlacement(nullptr, G4ThreeVector(fCubeX, fCubeY, zCubeWell - zLiquidWell), fCubeLV, "Cube",
                          fLiquidCubeWellLV, false, 0, true);

        logicWellC5 = new G4LogicalVolume(solidWell, air, "WellCube");
        logicWellC5->SetVisAttributes(wellVis);
        new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zLiquidWell), fLiquidCubeWellLV, "Liquid",
                          logicWellC5, false, 0, true);
    }

    // ---------- Place the 48 wells ----------
    G4double zWell = -fPlateZ / 2 + fFloor + wellDepth / 2;
    for (G4int row = 0; row < kRows; ++row) {
        for (G4int col = 0; col < kCols; ++col) {
            G4double x = (col - (kCols - 1) / 2.0) * fWellPitch;
            G4double y = ((kRows - 1) / 2.0 - row) * fWellPitch;
            G4int index = row * kCols + col;
            new G4PVPlacement(nullptr, G4ThreeVector(x, y, zWell),
                              index == fCubeWellIndex ? logicWellC5 : logicWell, "Well",
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

    // ---------- Sr-90 source: strontium disc between two silver discs, between two steel discs ----------
    // Five stacked cylinders with touching faces. From the bottom: steel, silver, strontium, silver, steel.
    // Only the strontium disc contains decaying nuclei (see PrimaryGenerator); the silver and steel
    // are inactive and only absorb / scatter the betas.
    struct Layer {
        const char* name;
        G4Material* material;
        G4double thickness;
        G4double radius;
        G4Colour colour;
    };
    const std::vector<Layer> stack = {
        {"SteelBottom",  steel,    fSteelThickness,  fHousingRadius, G4Colour(0.45, 0.45, 0.5, 0.6)},
        {"SilverBottom", silver,   fSilverThickness, fHousingRadius, G4Colour(0.85, 0.85, 0.9, 0.6)},
        {"Source",       strontium, fSourceThickness, fSourceRadius,  G4Colour(1.0, 0.8, 0.0, 1.0)},
        {"SilverTop",    silver,   fSilverThickness, fHousingRadius, G4Colour(0.85, 0.85, 0.9, 0.6)},
        {"SteelTop",     steel,    fSteelThickness,  fHousingRadius, G4Colour(0.45, 0.45, 0.5, 0.6)}};
    const std::size_t kSourceIndex = 2;

    G4double stackHeight = 0.;
    for (const auto& l : stack) stackHeight += l.thickness;

    G4double zLayer = -fPlateZ / 2 - stackHeight;   // top face of the stack touches the bottom of the plate
    for (std::size_t i = 0; i < stack.size(); ++i) {
        const Layer& l = stack[i];
        auto* solid = new G4Tubs(l.name, 0., l.radius, l.thickness / 2, 0., 360. * deg);
        auto* logic = new G4LogicalVolume(solid, l.material, l.name);
        logic->SetVisAttributes(G4VisAttributes(l.colour));
        const G4double zCentre = zLayer + l.thickness / 2;
        new G4PVPlacement(nullptr, G4ThreeVector(fSourceX, fSourceY, zCentre), logic, l.name,
                          logicWorld, false, 0, true);
        if (i == kSourceIndex) fSourceCentre = G4ThreeVector(fSourceX, fSourceY, zCentre);
        zLayer += l.thickness;
    }

    return physWorld;
}
