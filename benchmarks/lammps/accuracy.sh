#!/usr/bin/zsh

module switch intel intel/16.0
module switch openmpi intelmpi
module unload python

mpirun -np 16 -host cluster-phi ../../machines/rwth-sb_phi/lammps-10Mar16/src/lmp_intel_phi_default_vector -in in.tersoff-acc -pk intel 1 mode $1 -sf intel -v p vanilla | tee acc-$1
