#! /bin/bash

common="-no_partial -labels -byte -step 1 1 1 -start -20 -25 -30 -nelements 40 50 60"

# define boundary conditions
# voxels with non-zero value are set at fixed potential:
# voxels with zero value will be filled with field potential

# make a thin plate - will have high potential (=1.0) (label=2)
make_phantom -rectangle -center 0 0 -10 -width 20 20 1 -fill_value 2 -edge_value 2 \
    $common \
    pos_plate_1.mnc -clob

# make anothe thin plate - will have high potential (=1.0) (label=2)
make_phantom -rectangle -center 0 0 10 -width 20 20 1 -fill_value 2 -edge_value 2 \
    $common \
    pos_plate_2.mnc -clob


# make a small grounded stick - will have low potential (=0.0) (label=1)
make_phantom -rectangle -center 10 0 0 -width 2 20 2 -fill_value 1 -edge_value 1 \
    $common \
    neg_cube.mnc -clob

# combine into a single volume
mincmath -max -byte -labels \
    pos_plate_1.mnc pos_plate_2.mnc neg_cube.mnc \
    boundaries.mnc  -clob
