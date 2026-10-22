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
    G4LogicalVolume* GetCellLayerLV() const { return fCellLayerLV; }

    // True if energy deposited in this logical volume counts towards the dose
    // (the cell layer only, or cell layer + liquid above it, see fScoreCellLayerOnly)
    G4bool IsScoringVolume(const G4LogicalVolume* lv) const
    {
        return lv == fCellLayerLV || (!fScoreCellLayerOnly && lv == fLiquidLV);
    }
    G4bool ScoreCellLayerOnly() const { return fScoreCellLayerOnly; }
    G4double GetCellLayerThickness() const { return fCellLayerThickness; }

    // Mass of the scored region in ONE well (what the dose is divided by)
    G4double GetScoredMass() const;

    // Active disc of the source (for primary generation)
    G4double GetSourceRadius() const { return fSourceRadius; }
    G4double GetSourceThickness() const { return fSourceThickness; }
    G4ThreeVector GetSourceCentre() const { return fSourceCentre; }

private:
    // ---- Edit these to change the setup ----
    // Plate (standard SBS footprint)
    G4double fPlateX = 127.76 * mm;
    G4double fPlateY = 85.48 * mm;
    G4double fPlateZ = 17.0 * mm;       // plate height (the lid is separate, see below)
    G4double fFloor = 1.0 * mm;         // well floor thickness
    G4double fWellDiameter = 11.0 * mm;
    G4double fWellPitch = 13.0 * mm;
    G4double fFillHeight = 2.0 * mm;    // liquid depth (includes the cell layer)

    // Scoring region. The bottom `fCellLayerThickness` of the liquid is a separate volume
    // standing in for the adherent cell monolayer (water-equivalent).
    //   true  : dose = energy deposited in that layer only, divided by its mass
    //   false : dose = energy deposited in the whole liquid column, divided by its mass
    G4bool fScoreCellLayerOnly = true;
    G4double fCellLayerThickness = 50.0 * um;

    // Lid (polystyrene): flat top held above the plate, plus a skirt around the outside.
    // Typical values - set them to match the datasheet of your plate.
    G4bool fHasLid = true;              // false = no lid at all
    G4double fLidThickness = 1.0 * mm;  // thickness of the flat top
    G4double fLidAirGap = 1.0 * mm;     // air between plate top and underside of the lid
    G4double fLidSkirtHeight = 4.0 * mm; // skirt length below the underside of the lid top (0 = no skirt)
    G4double fLidWall = 1.0 * mm;       // skirt wall thickness
    G4double fLidClearance = 0.5 * mm;  // air gap between plate side and skirt

    // Source applicator
    G4double fSourceRadius = 8.2 * mm;       // active disc radius
    G4double fSourceThickness = 1.2 * um;    // active disc thickness
    G4double fHousingRadius = 8.2 * mm;
    G4double fHousingHeight = 101.2 * um;
    G4double fWindow = 50.0 * um;             // thickness of housing between active disc and outside
    G4double fGap = 30.0 * mm;               // source window to top of plate (NOT to the lid)
    G4double fSourceX = fWellPitch/2;            // lateral offset of source
    G4double fSourceY = fWellPitch/2;
    // ----------------------------------------

    G4LogicalVolume* fLiquidLV = nullptr;      // liquid above the cell layer
    G4LogicalVolume* fCellLayerLV = nullptr;   // bottom layer of the liquid
    G4ThreeVector fSourceCentre;
};

#endif
