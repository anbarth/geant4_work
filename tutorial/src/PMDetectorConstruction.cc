#include "PMDetectorConstruction.hh"

PMDetectorConstruction::PMDetectorConstruction()
{

}

PMDetectorConstruction::~PMDetectorConstruction()
{
    
}

// the asterisk means that this function returns only a pointer to a g4vphysicalvolume
// apparently the asterisk could be instead attached to the return type, but i wont mess with anything for now
G4VPhysicalVolume *PMDetectorConstruction::Construct()
{
    // volumes must be separated or else many things break
    // G4 kindly includes a check for this
    G4bool checkOverlaps = true;

    // the asterisks declare nist as a pointer to the instantiated object
    G4NistManager *nist = G4NistManager::Instance();
    // nist->foo() calls the member fxn foo() for the object nist points to
    G4Material *worldMat = nist->FindOrBuildMaterial("G4_AIR");
    G4Material *leadMat = nist->FindOrBuildMaterial("G4_Pb");
    G4Material *detMat = nist->FindOrBuildMaterial("G4_SODIUM_IODIDE");

    // dimensions of the world, with m being meter
    G4double xWorld = 1. * m;
    G4double yWorld = 1. * m;
    G4double zWorld = 1. * m;

    // set up the world volume
    // solid volume is just a volume of space
    G4Box *solidWorld = new G4Box("solidWorld", 0.5*xWorld, 0.5*yWorld, 0.5*zWorld);
    // logical volume is space+properties (in this case, air)
    G4LogicalVolume *logicWorld = new G4LogicalVolume(solidWorld, worldMat, "logicalWorld");
    // arguments: rotation (0=none), position, ... other stuff
    G4VPhysicalVolume *physWorld = new G4PVPlacement(0, G4ThreeVector(0.,0.,0.), logicWorld, "physWorld", 0, false, 0, checkOverlaps);

    // make a lead absorber
    G4double leadThickeness = 2. * mm;
    G4double leadSize = 10. * cm;
    G4Box *solidLead = new G4Box("solidLead", 0.5*leadSize, 0.5*leadSize, 0.5*leadThickeness);
    G4LogicalVolume *logicLead = new G4LogicalVolume(solidLead, leadMat, "logicLead");
    G4VPhysicalVolume *physLead = new G4PVPlacement(0, G4ThreeVector(0., 0., 5. * cm), logicLead, "physLead", logicWorld, false, checkOverlaps);

    // set the appearance of the lead absorber
    G4VisAttributes *leadVisAtt = new G4VisAttributes(G4Color(1.0, 0.0, 0.0, 0.5));
    leadVisAtt->SetForceSolid(true);
    logicLead->SetVisAttributes(leadVisAtt);

    G4double detectorSize = 10.0 * cm;
    G4Box *solidDetector = new G4Box("solidDetector", 0.5*detectorSize, 0.5*detectorSize, 0.5*detectorSize);
    logicDetector = new G4LogicalVolume(solidDetector, detMat, "logicDetector");
    G4VPhysicalVolume *physDetector = new G4PVPlacement(0, G4ThreeVector(0., 0., 10.5 * cm), logicDetector, "physDetector", logicWorld, false, checkOverlaps);

    // appearance for detector
    G4VisAttributes *detVisAtt = new G4VisAttributes(G4Color(1.0, 1.0, 0.0, 0.5));
    detVisAtt->SetForceSolid(true);
    logicDetector->SetVisAttributes(detVisAtt);

    return physWorld;

}

void PMDetectorConstruction::ConstructSDandField()
{
    PMSensitiveDetector *sensDet = new PMSensitiveDetector("sensitive detector");
    logicDetector->SetSensitiveDetector(sensDet);

    // the detector will only be "set off" at the right moments if you
    // "register" it with the sensitive detector manager
    G4SDManager::GetSDMpointer()->AddNewDetector(sensDet);
}