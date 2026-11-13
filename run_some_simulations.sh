#!/bin/bash
set -e

##### Film spectrum
echo "=========" >> timelog.txt
date >> timelog.txt

cd Sr90Film_spectrum
cd build
cmake ..
make -j32
./sim run.mac
cd ..
cd ..


##### Dose in open position
echo "=========" >> timelog.txt
date >> timelog.txt

cd Sr90Water
cd build
cmake ..
make -j32
./sim run.mac
cd ..
cd ..


##### Plate dosage and spectrum
echo "=========" >> timelog.txt
cd Sr90Plate_spectrum

# Source positions to scan (multiples of fWellPitch)
POSITIONSX=(0.5 0 0)
POSITIONSY=(0.5 0.5 0)

# Base name for the results folders
RUN_TAG="2026.11.13/moving_source"

for i in "${!POSITIONSX[@]}"; do
    posx=${POSITIONSX[$i]}
    posy=${POSITIONSY[$i]}
    outdir="results/${RUN_TAG}_pos${i}"

    # 1. Write down the time
    date >> timelog.txt

    # 2. Edit source position (line 120)
    sed -i '' "140c\\
    G4double fSourceX = fWellPitch * ${posx};
" include/DetectorConstruction.hh
    sed -i '' "141c\\
    G4double fSourceY = fWellPitch * ${posy};
" include/DetectorConstruction.hh

    # 3. Build and run the simulation
    (
        cd build
        cmake ..
        make -j32
        ./sim run.mac
    )



    # 4. Make a directory and move the results there
    mkdir -p "${outdir}"
    mv build/well_dose.csv "${outdir}/"
    mv build/voxel_dose.csv "${outdir}/"
 

done

cd ..
cd ..

