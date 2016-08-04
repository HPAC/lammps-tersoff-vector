#!/usr/bin/zsh
module load gcc/4.9
module switch openmpi intelmpi
module switch intel intel/16.0
module load cuda
h=`hostname -s`
d=../../machines/$1/lammps-10Mar16/src

for s in tersoff tersoff_gpu; do
  # GPU
  for p in single double mixed; do
    e=lmp_mpi_gpu_$p
    echo $e
    for i in `seq 0 5`; do
      $d/$e -in in.$s -v p vanilla -sf gpu > results/$h-$s-$e-$i
    done
  done
  for v in vect novect; do
    e=lmp_kokkos_cuda_$v
    echo $e
    for i in `seq 0 5`; do
      $d/$e -in in.$s -v p kokkos -sf kk -k on t 0 g 1 > results/$h-$s-$e-$i
    done
  done
done
