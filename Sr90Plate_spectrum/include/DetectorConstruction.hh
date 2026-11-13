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

    // ---- Voxel scoring: cell layer of the chosen wells split into fVoxelSize x fVoxelSize voxels ----
    // The geometry is not changed: energy deposited in the cell layer of those wells is binned
    // by (x, y) position in the stepping action. Voxels cut by the well wall have a smaller
    // mass (only the area inside the well counts).
    // Each chosen well is a "slot" 0..N-1 in the order of fVoxelWells. A voxel of the well in
    // slot s with per-well index v (0 .. GetVoxelsPerWell()-1) has the global index
    // s * GetVoxelsPerWell() + v.
    G4bool VoxelScoring() const { return fVoxelScoring; }
    G4int GetNumVoxelWells() const { return static_cast<G4int>(fVoxelWellIndex.size()); }
    G4int GetVoxelWellIndex(G4int slot) const { return fVoxelWellIndex[slot]; }   // well index = row*kCols + col
    G4int GetVoxelSlot(G4int well) const { return fVoxelSlot.empty() ? -1 : fVoxelSlot[well]; }   // -1: not voxelised
    G4int GetVoxelsPerSide() const { return fNVox; }
    G4int GetVoxelsPerWell() const { return fNVox * fNVox; }
    G4int GetNumVoxels() const { return GetNumVoxelWells() * GetVoxelsPerWell(); }   // all chosen wells
    G4double GetVoxelSize() const { return fVoxelSize; }
    G4ThreeVector GetWellCentre(G4int well) const;   // x, y of a well in the world frame
    // Per-well voxel index containing (x, y) given relative to the well axis; -1 if outside the grid
    G4int GetVoxelIndex(G4double xLocal, G4double yLocal) const;
    // Voxel centre relative to the well axis (per-well index)
    void GetVoxelCentre(G4int index, G4double& xLocal, G4double& yLocal) const;
    G4double GetVoxelArea(G4int index) const { return fVoxelArea[index]; }   // part inside the well
    G4double GetVoxelMass(G4int index) const;

    // Active disc of the source (for primary generation)
    G4double GetSourceRadius() const { return fSourceRadius; }
    G4double GetSourceThickness() const { return fSourceThickness; }
    G4ThreeVector GetSourceCentre() const { return fSourceCentre; }

    // ---- Electron spectrum in a plane near the steel filter ----
    // Electrons crossing the plane (a disc of radius fSpectrumRadius around the source axis, perpendicular
    // to z) in the direction away from the source are histogrammed by kinetic energy. The plane is
    // virtual: it is not a volume, so it can be put anywhere, also inside other materials.
    G4double GetSpectrumPlaneZ() const { return fSpectrumPlaneZ; }
    G4int GetSpectrumDirection() const { return fSpectrumDirection; }   // +1: electrons moving to +z are counted
    G4double GetSpectrumRadius() const { return fSpectrumRadius; }
    G4int GetNumSpectrumBins() const { return static_cast<G4int>(std::lround(fSpectrumEmax / fSpectrumBinWidth)); }
    G4double GetSpectrumBinWidth() const { return fSpectrumBinWidth; }

private:
    // ---- Edit these to change the setup ----
    // Plate (standard SBS footprint)
    G4double fPlateX = 127.63 * mm;
    G4double fPlateY = 85.47 * mm;
    G4double fPlateZ = 19.18 * mm;       // total plate height (the lid is separate, see below) (=well depth+floor thickness=17.78+1.4)
    G4double fFloor = 1.4 * mm;         // thickness of the rectangular bottom sheet (= well floor)
    G4double fWellDiameter = 10.0 * mm; // INNER diameter of each well tube
    G4double fWallThickness = 1.0 * mm; // wall of each well tube. Outer diameter = fWellDiameter +
                                         // 2*fWallThickness must not exceed fWellPitch (13 mm); 0.95 mm
                                         // leaves a 0.1 mm air gap so neighbouring tubes do not touch
    G4double fWellPitch = 12.9 * mm;
    G4double fFillHeight = 5.0 * mm;    // liquid depth (includes the cell layer)

    // Scoring region. The bottom `fCellLayerThickness` of the liquid is a separate volume
    // standing in for the adherent cell monolayer (water-equivalent).
    //   true  : dose = energy deposited in that layer only, divided by its mass
    //   false : dose = energy deposited in the whole liquid column, divided by its mass
    G4bool fScoreCellLayerOnly = true;
    G4double fCellLayerThickness = 50.0 * um;

    // Voxel scoring in the cell layer of the chosen wells (writes voxel_dose.csv). Off by default.
    // The voxels are fVoxelSize x fVoxelSize in x and y and fCellLayerThickness thick, laid out on
    // a grid centred on each well axis.
    G4bool fVoxelScoring = true;
    // Wells to divide into voxels: {"ALL"} for all 48, or a list such as {"C5"} or {"C4", "C5", "B5"}.
    // (C5 is the well right under the source.)
    std::vector<std::string> fVoxelWells = {"ALL"};
    G4double fVoxelSize = 0.5 * mm;

    // Lid (polystyrene): flat top held above the plate, plus a skirt around the outside.
    // Typical values - set them to match the datasheet of your plate.
    G4bool fHasLid = true;              // false = no lid at all
    G4double fLidThickness = 1.0 * mm;  // thickness of the flat top
    G4double fLidAirGap = 0.5 * mm;     // air between plate top and underside of the lid
    G4double fLidSkirtHeight = 4.0 * mm; // skirt length below the underside of the lid top (0 = no skirt)
    G4double fLidWall = 1.0 * mm;       // skirt wall thickness
    G4double fLidClearance = 0.5 * mm;  // air gap between plate side and skirt

    // Source: strontium disc, then a silver disc on each side, then a stainless-steel disc on each side
    G4double fSourceRadius = 4.1 * mm;       // strontium disc radius
    G4double fSourceThickness = 1.2 * um;    // strontium disc thickness
    G4double fSilverThickness = 50.0 * um;   // each of the two silver discs
    G4double fSteelThickness = 50.0 * um;    // each of the two outer stainless-steel discs  (placeholder value)

    // Spectrum of the electrons crossing a plane next to the steel filter
    G4bool fSpectrumTargetSide = true;       // true: the side of the source facing the target (what reaches it),
                                             // false: the far side
    G4double fSpectrumPlaneOffset = fFloor; // distance of the plane from the outer face of the steel on that
                                             // side; positive = further from the source, negative = inside the steel
    G4double fSpectrumRadius = 4.1 * mm;     // radius of the plane, centred on the source axis
    G4double fSpectrumEmax = 2.5 * MeV;      // Y-90 end point is 2.28 MeV; higher energies go in an overflow bin
    G4double fSpectrumBinWidth = 5.0 * keV;
    G4double fHousingRadius = 4.1 * mm;      // radius of the silver and steel discs
    G4double fGap = 30.0 * mm;               // source window to top of plate (NOT to the lid)
    G4double fSourceX = fWellPitch * 0.5;
    G4double fSourceY = fWellPitch * 0.5;
    // ----------------------------------------

    void SetupVoxels();   // fills the members below; called from the constructor

    G4LogicalVolume* fLiquidLV = nullptr;      // liquid above the cell layer
    G4LogicalVolume* fCellLayerLV = nullptr;   // bottom layer of the liquid
    G4ThreeVector fSourceCentre;
    G4double fSpectrumPlaneZ = 0.;          // set in Construct
    G4int fSpectrumDirection = 1;

    std::vector<G4int> fVoxelWellIndex;        // slot -> well index (row*kCols + col)
    std::vector<G4int> fVoxelSlot;             // well index -> slot, -1 if not voxelised
    G4int fNVox = 0;                           // voxels per side
    std::vector<G4double> fVoxelArea;          // area of each voxel inside a well (same for every well)
};

#endif
