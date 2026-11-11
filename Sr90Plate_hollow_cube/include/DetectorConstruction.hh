#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>
#include <string>
#include <vector>

class G4LogicalVolume;
class G4VPhysicalVolume;

// Geometry: a 48-well plate (6 rows A-F x 8 columns 1-8) centred at the origin,
// built as a rectangular bottom sheet with 48 separate hollow tubes (walls fWallThickness) standing on it,
// with a Sr-90/Y-90 source stack below it. A 1 mm cube in one well (fCubeWell) is used to measure the
// spectrum of the electrons entering it.
//
// Well copy numbers (of the "air well" volume) run 0..47 as  index = row*8 + col
//   row 0 = A (top, +y) ... row 5 = F;   col 0 = 1 (left, -x) ... col 7 = 8
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;

    static constexpr G4int kRows = 6;
    static constexpr G4int kCols = 8;
    static constexpr G4int kWells = kRows * kCols;

    G4LogicalVolume* GetLiquidLV() const { return fLiquidLV; }

    // True if energy deposited in this logical volume counts towards the dose (the liquid of each well;
    // the cube is part of the liquid of its well, see the stepping action)
    G4bool IsScoringVolume(const G4LogicalVolume* lv) const { return lv == fLiquidLV || lv == fLiquidCubeWellLV; }

    // Mass of the scored region in ONE well (what the dose is divided by)
    G4double GetScoredMass() const;

    // Active disc of the source (for primary generation)
    G4double GetSourceRadius() const { return fSourceRadius; }
    G4double GetSourceThickness() const { return fSourceThickness; }
    G4ThreeVector GetSourceCentre() const { return fSourceCentre; }

    // ---- Electron spectrum entering a small cube in one well ----
    // A cube of the same material as the liquid (so the physics is unchanged) sits inside the liquid of
    // one well. Electrons that enter it from outside are histogrammed by kinetic energy.
    G4bool IsCube(const G4LogicalVolume* lv) const { return lv != nullptr && lv == fCubeLV; }
    G4int GetCubeWellIndex() const { return fCubeWellIndex; }          // row*kCols + col
    G4double GetCubeBottomZ() const { return fCubeBottomZ; }           // z of the bottom face, world frame
    G4bool CubeBottomFaceOnly() const { return fCubeBottomFaceOnly; }
    G4double GetCubeSize() const { return fCubeSize; }
    std::string GetCubeWellName() const { return fCubeWell; }
    G4int GetNumSpectrumBins() const { return static_cast<G4int>(std::lround(fSpectrumEmax / fSpectrumBinWidth)); }
    G4double GetSpectrumBinWidth() const { return fSpectrumBinWidth; }

private:
    // ---- Edit these to change the setup ----
    // Plate (standard SBS footprint)
    G4double fPlateX = 127.76 * mm;
    G4double fPlateY = 85.48 * mm;
    G4double fPlateZ = 17.0 * mm;       // total plate height (the lid is separate, see below)
    G4double fFloor = 1.0 * mm;         // thickness of the rectangular bottom sheet (= well floor)
    G4double fWellDiameter = 11.0 * mm; // INNER diameter of each well tube
    G4double fWallThickness = 0.95 * mm; // wall of each well tube. Outer diameter = fWellDiameter +
                                         // 2*fWallThickness must not exceed fWellPitch (13 mm); 0.95 mm
                                         // leaves a 0.1 mm air gap so neighbouring tubes do not touch
    G4double fWellPitch = 13.0 * mm;
    G4double fFillHeight = 5.0 * mm;    // liquid depth

    // Lid (polystyrene): flat top held above the plate, plus a skirt around the outside.
    // Typical values - set them to match the datasheet of your plate.
    G4bool fHasLid = true;              // false = no lid at all
    G4double fLidThickness = 1.0 * mm;  // thickness of the flat top
    G4double fLidAirGap = 1.0 * mm;     // air between plate top and underside of the lid
    G4double fLidSkirtHeight = 4.0 * mm; // skirt length below the underside of the lid top (0 = no skirt)
    G4double fLidWall = 1.0 * mm;       // skirt wall thickness
    G4double fLidClearance = 0.5 * mm;  // air gap between plate side and skirt

    // Source: strontium disc, then a silver disc on each side, then a stainless-steel disc on each side
    G4double fSourceRadius = 4.1 * mm;       // strontium disc radius
    G4double fSourceThickness = 1.2 * um;    // strontium disc thickness
    G4double fSilverThickness = 50.0 * um;   // each of the two silver discs
    G4double fSteelThickness = 50.0 * um;    // each of the two outer stainless-steel discs  (placeholder value)
    G4double fHousingRadius = 4.1 * mm;      // radius of the silver and steel discs
    G4double fGap = 30.0 * mm;               // source window to top of plate (NOT to the lid)
    G4double fSourceX = fWellPitch * 0.0;
    G4double fSourceY = fWellPitch * 0.0;

    // Cube that the spectrum is measured on (a box of the same material as the liquid, inside the liquid)
    std::string fCubeWell = "C5";             // well that contains the cube
    G4double fCubeSize = 1.0 * mm;            // edge length
    G4double fCubeX = 0.0 * mm;               // centre of the cube relative to the well axis
    G4double fCubeY = 0.0 * mm;
    G4double fCubeBottomHeight = 0.0 * mm;    // height of the bottom face above the floor of the well
                                              // (0 = resting on the floor; the cube must stay inside the liquid)
    G4bool fCubeBottomFaceOnly = false;        // true: count only electrons entering through the bottom face, i.e.
                                              // coming from below (the source side); false: through any face

    // Spectrum binning
    G4double fSpectrumEmax = 2.5 * MeV;      // Y-90 end point is 2.28 MeV; higher energies go in an overflow bin
    G4double fSpectrumBinWidth = 5.0 * keV;
    // ----------------------------------------

    G4LogicalVolume* fLiquidLV = nullptr;      // liquid of the wells
    G4ThreeVector fSourceCentre;
    G4LogicalVolume* fCubeLV = nullptr;
    G4LogicalVolume* fLiquidCubeWellLV = nullptr;   // the liquid of the well with the cube (a copy of the liquid)
    G4int fCubeWellIndex = -1;
    G4double fCubeBottomZ = 0.;
};

#endif
