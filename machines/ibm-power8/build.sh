#!/bin/bash

module purge
module load tbb
module load compiler/gcc-5.3.0

set -o verbose
cp -r ../lammps-10Mar16 .

# apply patches
cp ../../src lammps-10Mar16/src/
cd lammps-10Mar16/src/

make yes-MANYBODY
make yes-USER-OMP
make yes-USER-INTEL
make no-KOKKOS

binaries=power8

make power8 -j10
mv lmp_power8 ../../
