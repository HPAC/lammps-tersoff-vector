#!/bin/bash

module purge 
module load compiler/gcc-5.3.0
module load hpcx

export OMP_NUM_THREADS=4
# basic run
mpirun -np 20 ./lmp_power8 -in in.tersoff-acc -pk intel 1 mode single -sf intel -v p vanilla | tee acc-single > accuracy_mpi_results
#./lmp_power8 -in in.tersoff -pk omp 0 -pk intel 1 mode single -sf intel -v p vanilla
