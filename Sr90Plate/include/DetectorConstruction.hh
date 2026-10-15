#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"

class G4LogicalVolume;
class G4VPhysicalVolume;

// Geometry: a 48-well plate (6 rows A-F x 8 columns 1-8) centred at the origin,
// with a Sr-90/Y-90 applicator held above it, window facing down (-z).
//
// Well copy numbers (of the "air well" volume) run 0..47 as  index = row*8 + col
//   row 0 = A (top, +y) ... row 5 = F;   col 0 = 1 (left, -x) ... col 7 = 8
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction();
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;

    static constexpr G4int kRows = 6;
    static constexpr G4int kCols = 8;
    static constexpr G4int kWells = kRows * kCols;

    G4LogicalVolume* GetLiquidLV() const { return fLiquidLV; }
    G4double GetLiquidMass() const;

    // Active disc of the source (for primary generation)
    G4double GetSourceRadius() const { return fSourceRadius; }
    G4double GetSourceThickness() const { return fSourceThickness; }
    G4ThreeVector GetSourceCentre() const { return fSourceCentre; }

private:
    // ---- Edit these to change the setup ----
    // Plate (standard SBS footprint)
    G4double fPlateX = 127.76 * mm;
    G4double fPlateY = 85.48 * mm;
    G4double fPlateZ = 17.0 * mm;       // plate height (no lid)
    G4double fFloor = 1.0 * mm;         // well floor thickness
    G4double fWellDiameter = 11.0 * mm;
    G4double fWellPitch = 13.0 * mm;
    G4double fFillHeight = 2.0 * mm;    // liquid depth 

    // Source applicator
    G4double fSourceRadius = 5.0 * mm;       // active disc radius
    G4double fSourceThickness = 0.1 * mm;    // active disc thickness
    G4double fHousingRadius = 12.5 * mm;
    G4double fHousingHeight = 8.0 * mm;
    G4double fWindow = 0.1 * mm;             // steel between active disc and outside
    G4double fGap = 30.0 * mm;               // source window to top of plate
    G4double fSourceX = 0.0 * mm;            // lateral offset of source
    G4double fSourceY = 0.0 * mm;
    // ----------------------------------------

    G4LogicalVolume* fLiquidLV = nullptr;
    G4ThreeVector fSourceCentre;
};

#endif
