#!/usr/bin/zsh
module load gcc/4.9
module switch openmpi intelmpi
module switch intel intel/16.0
module load cuda
h=`hostname -s`
d=../../machines/$1/lammps-10Mar16/src

export OMP_NUM_THREADS=1
for s in tersoff_bench; do
  for p in single double mixed; do
    e=lmp_intel_phi_default_vector
    echo $s $p
    for i in `seq 0 5`; do
      echo $i
      mpirun -np $2 -host $h $d/$e -in in.$s -v p vanilla -sf intel -pk intel 1 balance 1 mode $p > results/$h-$s-$e-$p-only-$i
      mpirun -np $2 -host $h $d/$e -in in.$s -v p vanilla -sf intel -pk intel 1 balance 0 mode $p > results/$h-$s-$e-$p-host-$i
      mpirun -np $2 -host $h $d/$e -in in.$s -v p vanilla -sf intel -pk intel 1 balance -1 mode $p > results/$h-$s-$e-$p-auto-$i
    done
  done
done
