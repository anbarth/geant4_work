#include "PMDetectorConstruction.hh"

#include "G4NistManager.hh"
#include "G4Material.hh"
#include "G4Element.hh"
#include "G4Box.hh"
#include "G4Tubs.hh"
#include "G4LogicalVolume.hh"
#include "G4PVPlacement.hh"
#include "G4VisAttributes.hh"
#include "G4Colour.hh"

PMDetectorConstruction::PMDetectorConstruction() {}

G4double PMDetectorConstruction::GetLiquidMass() const
{
    return fLiquidLV ? fLiquidLV->GetMass() : 0.;
}

G4VPhysicalVolume* PMDetectorConstruction::Construct()
{
    auto* nist = G4NistManager::Instance();

    // ---------- Materials ----------
    G4Material* air = nist->FindOrBuildMaterial("G4_AIR");
    G4Material* polystyrene = nist->FindOrBuildMaterial("G4_POLYSTYRENE");
    G4Material* medium = nist->FindOrBuildMaterial("G4_WATER");   // culture medium ~ water
    G4Material* steel = nist->FindOrBuildMaterial("G4_STAINLESS-STEEL");

    // Active layer: strontium titanate (SrTiO3), a common Sr-90 carrier
    G4Element* elSr = nist->FindOrBuildElement("Sr");
    G4Element* elTi = nist->FindOrBuildElement("Ti");
    G4Element* elO = nist->FindOrBuildElement("O");
    auto* srtio3 = new G4Material("SrTiO3", 5.12 * g / cm3, 3);
    srtio3->AddElement(elSr, 1);
    srtio3->AddElement(elTi, 1);
    srtio3->AddElement(elO, 3);

    // ---------- World ----------
    G4double worldX = 100. * mm, worldY = 100. * mm, worldZ = 60. * mm;  // half-lengths
    auto* solidWorld = new G4Box("World", worldX, worldY, worldZ);
    auto* logicWorld = new G4LogicalVolume(solidWorld, air, "World");
    auto* physWorld = new G4PVPlacement(nullptr, G4ThreeVector(), logicWorld,
                                        "World", nullptr, false, 0, true);
    logicWorld->SetVisAttributes(G4VisAttributes::GetInvisible());

    // ---------- Plate body (polystyrene block, wells are placed inside) ----------
    auto* solidPlate = new G4Box("Plate", fPlateX / 2, fPlateY / 2, fPlateZ / 2);
    auto* logicPlate = new G4LogicalVolume(solidPlate, polystyrene, "Plate");
    new G4PVPlacement(nullptr, G4ThreeVector(), logicPlate, "Plate", logicWorld, false, 0, true);
    G4VisAttributes plateVis(G4Colour(0.6, 0.6, 0.8, 0.5));
    //plateVis.SetForceWireframe(true);
    logicPlate->SetVisAttributes(plateVis);

    // ---------- One well: air cylinder (open cavity) + liquid at the bottom ----------
    G4double wellDepth = fPlateZ - fFloor;
    auto* solidWell = new G4Tubs("Well", 0., fWellDiameter / 2, wellDepth / 2, 0., 360. * deg);
    auto* logicWell = new G4LogicalVolume(solidWell, air, "Well");
    // Drawn as a wireframe so the 48 well outlines are visible
    G4VisAttributes wellVis(G4Colour(0.2, 0.2, 0.2, 1.0));
    wellVis.SetForceWireframe(true);
    logicWell->SetVisAttributes(wellVis);

    auto* solidLiquid = new G4Tubs("Liquid", 0., fWellDiameter / 2, fFillHeight / 2, 0., 360. * deg);
    fLiquidLV = new G4LogicalVolume(solidLiquid, medium, "Liquid");
    fLiquidLV->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.3, 0.3, 1.0)));

    // Liquid sits on the floor of the well cavity (local z of the well cavity)
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, -wellDepth / 2 + fFillHeight / 2),
                      fLiquidLV, "Liquid", logicWell, false, 0, true);

    // ---------- Place the 48 wells ----------
    G4double zWell = -fPlateZ / 2 + fFloor + wellDepth / 2;
    for (G4int row = 0; row < kRows; ++row) {
        for (G4int col = 0; col < kCols; ++col) {
            G4double x = (col - (kCols - 1) / 2.0) * fWellPitch;
            G4double y = ((kRows - 1) / 2.0 - row) * fWellPitch;
            G4int index = row * kCols + col;
            new G4PVPlacement(nullptr, G4ThreeVector(x, y, zWell), logicWell, "Well",
                              logicPlate, false, index, true);
        }
    }

    // ---------- Sr-90 applicator above the plate, window facing down ----------
    auto* solidHousing = new G4Tubs("Housing", 0., fHousingRadius, fHousingHeight / 2, 0., 360. * deg);
    //auto* logicHousing = new G4LogicalVolume(solidHousing, steel, "Housing");
    auto* logicHousing = new G4LogicalVolume(solidHousing, air, "Housing");
    G4VisAttributes housingVis(G4Colour(0.6, 0.6, 0.6, 0.5));
    //housingVis.SetForceWireframe(true);
    logicHousing->SetVisAttributes(housingVis);

    G4double zHousing = fPlateZ / 2 + fGap + fHousingHeight / 2;
    new G4PVPlacement(nullptr, G4ThreeVector(fSourceX, fSourceY, zHousing), logicHousing,
                      "Housing", logicWorld, false, 0, true);

    auto* solidSource = new G4Tubs("Source", 0., fSourceRadius, fSourceThickness / 2, 0., 360. * deg);
    auto* logicSource = new G4LogicalVolume(solidSource, srtio3, "Source");
    logicSource->SetVisAttributes(G4VisAttributes(G4Colour(1.0, 0.8, 0.0)));

    // Active disc embedded in the steel, `fWindow` above the bottom face of the housing
    G4double zSourceLocal = -fHousingHeight / 2 + fWindow + fSourceThickness / 2;
    new G4PVPlacement(nullptr, G4ThreeVector(0, 0, zSourceLocal), logicSource,
                      "Source", logicHousing, false, 0, true);

    fSourceCentre = G4ThreeVector(fSourceX, fSourceY, zHousing + zSourceLocal);

    return physWorld;
}
