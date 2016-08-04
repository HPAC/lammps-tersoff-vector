import os
template = """#!/bin/bash
#@ wall_clock_limit = 00:30:00
#@ job_name = {node}-{name}
#@ job_type = parallel
#@ class = phi
#@ node = {node}
#@ tasks_per_node = 16
#@ node_usage = not_shared
#@ initialdir = {dir}/../../../benchmarks/lammps
#@ output = {dir}/{name}-{node}-$(jobid).out
#@ error  = {dir}/{name}-{node}-$(jobid).err
#@ queue

. /etc/profile
. /etc/profile.d/modules.sh

module unload intel
module load intel/15.0

module list

cat $LOADL_HOSTFILE

export I_MPI_DAPL_PROVIDER_LIST=ofa-v2-mlx4_0-1u
mpiexec ../../machines/lrz-ib_phi/lammps-10Mar16/src/lmp_intel_phi_default_vector -in in.tersoff_2e6 {pkg}
"""
pks=[("default","")]
dir=os.path.dirname(os.path.abspath(__file__))+'/run'
for prec in ["single", "double", "mixed"]:
    pks.append(("intel-%s" % prec, "-pk intel 0 mode %s -sf intel" % prec))
    for nphi in [1, 2]:
        for balance in [("auto",-1), ("spec",0.55), ("full",1)]:
            pks.append(("phi-%d-%s-%s" % (nphi,prec,balance[0]), "-pk intel %d mode %s balance %f -sf intel" % (nphi, prec, balance[1])))
for nodes in [1,2,4,8,12,16]:#,20,24,28,32]:
    for pkn, pk in pks:
        with open('%s/scale-%s-%d.job' % (dir,pkn, nodes), 'w') as f:
            f.write(template.format(dir=dir,node=nodes, pkg=pk, name="scale-"+pkn))

