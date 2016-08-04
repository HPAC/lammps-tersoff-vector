#!/bin/bash

#script build.script
set -o verbose

#[[ -e lammps.tar.gz ]] || wget http://lammps.sandia.gov/tars/lammps.tar.gz
#[[ -d lammps-14Mar16 ]] || tar xzf lammps.tar.gz
cp -r ../lammps-10Mar16 .
cd lammps-10Mar16/
patch -p1 -N < ../../ccflags.patch
cd src
cd STUBS
make
cd ..
make yes-USER-OMP
make yes-USER-INTEL
make no-KOKKOS
make no-GPU

binaries="serial"

function build_binaries() {
  set -o verbose
  for f in $binaries; do
    touch pair_tersoff_intel.h
    [[ -e lmp_${f}_$2 ]] && continue
    make $f "$1" -j4
    mv lmp_$f lmp_${f}_$2 || exit
  done
}

build_binaries "CCFLAGS=-DLMP_INTEL_NOOP -DLAMMPS_MEMALIGN=64"  "default_vector" 
build_binaries "CCFLAGS=-DLMP_INTEL_TERSOFF_PACK_I=true -DLAMMPS_MEMALIGN=64"  "packi_vector" 
build_binaries "CCFLAGS=-DLMP_INTEL_TERSOFF_PACK_I=false -DLAMMPS_MEMALIGN=64"  "nopacki_vector" 

build_binaries "CCFLAGS=-LMP_INTEL_VECTOR_MIC=NONE -DLMP_INTEL_VECTOR_HOST=NONE -DLAMMPS_MEMALIGN=64"  "default_scalar" 
build_binaries "CCFLAGS=-DLMP_INTEL_TERSOFF_PACK_I=true -LMP_INTEL_VECTOR_MIC=NONE -DLMP_INTEL_VECTOR_HOST=NONE -DLAMMPS_MEMALIGN=64"  "packi_scalar" 
build_binaries "CCFLAGS=-DLMP_INTEL_TERSOFF_PACK_I=false -LMP_INTEL_VECTOR_MIC=NONE -DLMP_INTEL_VECTOR_HOST=NONE -DLAMMPS_MEMALIGN=64"  "nopacki_scalar" 
