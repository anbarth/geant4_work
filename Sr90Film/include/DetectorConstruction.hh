#ifndef DetectorConstruction_hh
#define DetectorConstruction_hh

#include "G4VUserDetectorConstruction.hh"
#include "G4ThreeVector.hh"
#include "G4SystemOfUnits.hh"
#include <cmath>

class G4LogicalVolume;
class G4VPhysicalVolume;

// Geometry: a three-layer radiochromic film in air, centred on the z axis, the active layer
// sandwiched between two base layers (top of the stack at z = 0):
//   top base layer      (z from -fTopBaseThickness to 0)
//   active layer        (below it)
//   bottom base layer   (below the active layer)
// A Sr-90/Y-90 applicator sits on top with its window facing down (-z), fGap above the top base
// layer (fGap = 0: window in direct contact with the film).
// Dose is scored in voxels of the ACTIVE layer: fVoxelSize x fVoxelSize in x and y, and the whole
// thickness of the layer in z. Voxel (ix, iy) has index iy*nx + ix, ix along +x, iy along +y.
class DetectorConstruction : public G4VUserDetectorConstruction
{
public:
    DetectorConstruction() = default;
    ~DetectorConstruction() override = default;

    G4VPhysicalVolume* Construct() override;

    G4bool IsActiveLayer(const G4LogicalVolume* lv) const { return lv == fActiveLV; }

    // ---- voxel grid over the active layer ----
    G4int GetNumVoxelsX() const { return static_cast<G4int>(std::lround(fFilmX / fVoxelSize)); }
    G4int GetNumVoxelsY() const { return static_cast<G4int>(std::lround(fFilmY / fVoxelSize)); }
    G4int GetNumVoxels() const { return GetNumVoxelsX() * GetNumVoxelsY(); }
    G4double GetVoxelSize() const { return fVoxelSize; }
    G4double GetVoxelArea() const { return fVoxelSize * fVoxelSize; }
    G4double GetActiveThickness() const { return fActiveThickness; }

    // Voxel containing the point (x, y) (global coordinates); -1 if outside the film
    G4int GetVoxelIndex(G4double x, G4double y) const;
    // Centre of a voxel relative to the centre of the film (what voxel_dose.csv contains)
    void GetVoxelCentre(G4int index, G4double& xRel, G4double& yRel) const;

    G4double GetVoxelMass() const;          // mass of one voxel of the active layer
    G4double GetActiveLayerMass() const;    // mass of the whole active layer

    // Active disc of the source (for primary generation)
    G4double GetSourceRadius() const { return fSourceRadius; }
    G4double GetSourceThickness() const { return fSourceThickness; }
    G4ThreeVector GetSourceCentre() const { return fSourceCentre; }

private:
    // ---- Edit these to change the setup ----
    // Film (both layers have the same footprint)
    G4double fFilmX = 50.0 * mm;
    G4double fFilmY = 50.0 * mm;
    G4double fTopBaseThickness = 0.125 * mm;  // base layer above the active layer (facing the source)
    G4double fActiveThickness = 0.025 * mm;   // middle layer
    G4double fBaseThickness = 0.125 * mm;     // base layer below the active layer
    G4double fFilmCentreX = 0.0 * mm;         // film centre in x, y (z: top of the film is at z = 0)
    G4double fFilmCentreY = 0.0 * mm;

    // Voxels in the active layer. fFilmX and fFilmY must be whole multiples of fVoxelSize.
    G4double fVoxelSize = 0.5 * mm;

    // Source applicator
    G4double fSourceRadius = 8.2 * mm;       // active disc radius
    G4double fSourceThickness = 1.2 * um;    // active disc thickness
    G4double fHousingRadius = 8.2 * mm;
    G4double fHousingHeight = 101.2 * um;
    G4double fWindow = 50.0 * um;            // thickness of housing between active disc and outside
    G4double fGap = 0.0 * mm;                // window to top of the film (0 = in contact)
    G4double fSourceX = 0.0 * mm;            // lateral offset of the source from the z axis
    G4double fSourceY = 0.0 * mm;
    // ----------------------------------------

    G4LogicalVolume* fActiveLV = nullptr;
    G4ThreeVector fSourceCentre;
};

#endif
