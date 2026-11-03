#!/bin/bash
set -e

# Source positions to scan (multiples of fWellPitch). Edit this list as needed.
POSITIONS=(0 0.25 0.5)

# Base name for the results folders
RUN_TAG="2026.11.02/moving_source"

echo "=========" >> timelog.txt

for pos in "${POSITIONS[@]}"; do
    echo "=== Source position: fWellPitch * ${pos} ==="
    outdir="results/${RUN_TAG}_pos${pos}"

    # 1. Write down the time
    date >> timelog.txt

    # 2. Edit source position (line 120)
    sed -i '' "120c\\
    G4double fSourceX = fWellPitch * ${pos};
" include/DetectorConstruction.hh

    # 3. Build and run the simulation
    (
        cd build
        cmake ..
        make -j32
        ./sim run.mac
    )

    # 4. Make heatmap and voxel profile
    python3 plot_voxel_heatmap.py
    python3 voxel_line_profile.py --xy 0 6.5 --angle 0

    # 4. Make a directory and move the results there
    mkdir -p "${outdir}"
    mv build/well_dose.csv "${outdir}/"
    mv build/voxel_dose.csv "${outdir}/"
    mv voxel_dose_heatmap.png "${outdir}/"
    mv line_x0_y6.5_0deg.png "${outdir}/"
    mv line_x0_y6.5_0deg.csv "${outdir}/"

    

done