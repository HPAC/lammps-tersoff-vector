#!/bin/bash

module purge 
module load compiler/gcc-5.3.0
module load hpcx

# basic run
./lmp_power8 -in in.tersoff -pk omp 0 -pk intel 1 mode single -sf intel -v p vanilla
