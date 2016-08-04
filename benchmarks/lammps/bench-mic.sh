#!/bin/bash
h=`hostname -s`
d=../../machines/$1/lammps-10Mar16/src

# <arch$1> <binary-type$2> <total_cores$3> <factor_1$4> <factor_2$5>

export OMP_NUM_THREADS=1
for i in `seq 0 5`; do
  e=lmp_intel_$2_default_vector
  echo $i
  for p in single double mixed; do
    $d/$e -in in.tersoff -v p vanilla -sf intel -pk intel 0 mode $p > results/$h-tersoff-$e-$p-cpu_intel-$i
  done
  $d/$e -in in.tersoff -v p vanilla > results/$h-tersoff-$e-double-cpu_normal-$i
done
# multi core run
for s in tersoff tersoff_bench; do
  for i in `seq 0 5`; do
    e=lmp_intel_$2_default_vector
    echo $i
    for p in single double mixed; do
      mpirun -np $3 -host $h $d/$e -in in.$s -v p vanilla -sf intel -pk intel 0 mode $p > results/$h-$s-$e-$p-cpu_intel_all_mpi-$i 2>&1
      OMP_NUM_THREADS=$5 mpirun -np $4 -host $h $d/$e -in in.$s -v p vanilla -sf intel -pk intel 0 mode $p > results/$h-$s-$e-$p-cpu_intel_some_mpi-$i 2>&1
      OMP_NUM_THREADS=$4 mpirun -np $5 -host $h $d/$e -in in.$s -v p vanilla -sf intel -pk intel 0 mode $p > results/$h-$s-$e-$p-cpu_intel_some_omp-$i 2>&1
      OMP_NUM_THREADS=$3 d/$e -in in.$s -v p vanilla -sf intel -pk intel 0 mode $p > results/$h-$s-$e-$p-cpu_intel_all_omp-$i 2>&1
    done
    mpirun -np $3 -host $h $d/$e -in in.$s -v p vanilla > results/$h-$s-$e-double-cpu_normal_all_mpi-$i 2>&1 
    OMP_NUM_THREADS=$5 mpirun -np $4 -host $h $d/$e -in in.$s -v p vanilla -sf omp -pk omp 0 > results/$h-$s-$e-double-cpu_normal_some_mpi-$i 2>&1
    OMP_NUM_THREADS=$4 mpirun -np $5 -host $h $d/$e -in in.$s -v p vanilla -sf omp -pk omp 0 > results/$h-$s-$e-double-cpu_normal_some_omp-$i 2>&1
    OMP_NUM_THREADS=$3 d/$e -in in.$s -v p vanilla -sf omp -pk omp 0 > results/$h-$s-$e-double-cpu_normal_all_omp-$i 2>&1
  done
done
