Simulation of the Sr90 source irradiating a 48-well plate. The simulation measures the spectrum in the plane of the well floor, as well as the dose deposited to the bottom 50um of each well.



The simulation was run on 11/13/2026 exactly as written here, except that the source position was moved to three different spots:



source position 1

G4double fSourceX = fWellPitch \* 0.5;

G4double fSourceY = fWellPitch \* 0.5;



source position 2

G4double fSourceX = fWellPitch \* 0;

G4double fSourceY = fWellPitch \* 0.5;



source position 3

G4double fSourceX = fWellPitch \* 0;

G4double fSourceY = fWellPitch \* 0;



