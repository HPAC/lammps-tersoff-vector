#!/bin/bash

. /etc/profile
. /etc/profile.d/modules.sh
#script build.script
set -o verbose

#module load gcc/4.9
module unload intel
module load intel/15.0
module load tbb

[[ -d lammps-10Mar16 ]] || cp -r ../lammps-10Mar16 .
cd lammps-10Mar16/
patch -p1 -N < ../../ccflags.patch
cd src
make yes-USER-OMP
make yes-USER-INTEL

binaries="intel_cpu intel_phi"

function build_binaries() {
  set -o verbose
  for f in $binaries; do
    touch pair_tersoff_intel.h
    [[ -e lmp_${f}_$2 ]] && continue
    make $f "$1" "LIB=$MKL_LIB -L$TBB_LIBDIR" -j
    mv lmp_$f lmp_${f}_$2 || exit
  done
}

build_binaries "CCFLAGS=-DLMP_INTEL_NOOP $MKL_INC"  "default_vector" 
build_binaries "CCFLAGS=-DLMP_INTEL_TERSOFF_PACK_I=true $MKL_INC"  "packi_vector" 
build_binaries "CCFLAGS=-DLMP_INTEL_TERSOFF_PACK_I=false $MKL_INC"  "nopacki_vector" 

build_binaries "CCFLAGS=-LMP_INTEL_VECTOR_MIC=NONE -DLMP_INTEL_VECTOR_HOST=NONE $MKL_INC"  "default_scalar" 
build_binaries "CCFLAGS=-DLMP_INTEL_TERSOFF_PACK_I=true -LMP_INTEL_VECTOR_MIC=NONE -DLMP_INTEL_VECTOR_HOST=NONE $MKL_INC"  "packi_scalar" 
build_binaries "CCFLAGS=-DLMP_INTEL_TERSOFF_PACK_I=false -LMP_INTEL_VECTOR_MIC=NONE -DLMP_INTEL_VECTOR_HOST=NONE $MKL_INC"  "nopacki_scalar" 
