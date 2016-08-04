/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator

   Original Version:
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   See the README file in the top-level LAMMPS directory.

   -----------------------------------------------------------------------

   USER-CUDA Package and associated modifications:
   https://sourceforge.net/projects/lammpscuda/

   Christian Trott, christian.trott@tu-ilmenau.de
   Lars Winterfeld, lars.winterfeld@tu-ilmenau.de
   Theoretical Physics II, University of Technology Ilmenau, Germany

   See the README file in the USER-CUDA directory.

   This software is distributed under the GNU General Public License.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

/* ----------------------------------------------------------------------
   Contributing author: Paul Crozier (SNL)
------------------------------------------------------------------------- */

#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include "pair_lj_sdk_coul_long_cuda.h"
#include "pair_lj_sdk_coul_long_cuda_cu.h"
#include "cuda_data.h"
#include "atom.h"
#include "comm.h"
#include "force.h"
#include "kspace.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "neigh_request.h"
#include "cuda_neigh_list.h"
#include "update.h"
#include "integrate.h"
#include "respa.h"
#include "memory.h"
#include "error.h"
#include "user_cuda.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

PairLJSDKCoulLongCuda::PairLJSDKCoulLongCuda(LAMMPS *lmp) : PairLJSDKCoulLong(lmp)
{
  cuda = lmp->cuda;
   if(cuda == NULL)
        error->all(FLERR,"You cannot use a /cuda class, without activating 'cuda' acceleration. Provide '-c on' as command-line argument to LAMMPS..");

        allocated2 = false;
        lj_type_double = NULL;
        cuda->shared_data.pair.cudable_force = 1;
        cuda->setSystemParams();
}

/* ----------------------------------------------------------------------
   remember pointer to arrays in cuda shared data
------------------------------------------------------------------------- */

void PairLJSDKCoulLongCuda::allocate()
{
        if(! allocated) PairLJSDKCoulLong::allocate();
        int n = atom->ntypes;
        if(! allocated2)
        {
                allocated2 = true;


                  memory->create(lj_type_double,n+1,n+1,"pairlj:ljtypedouble");

                cuda->shared_data.pair.cut     = cut_lj;
                cuda->shared_data.pair.cut_coul= NULL;
                cuda->shared_data.pair.coeff1  = lj1;
                cuda->shared_data.pair.coeff2  = lj2;
                cuda->shared_data.pair.coeff3  = lj3;
                cuda->shared_data.pair.coeff4  = lj4;
                cuda->shared_data.pair.coeff5  = lj_type_double;
                cuda->shared_data.pair.offset  = offset;
                cuda->shared_data.pair.special_lj  = force->special_lj;
                cuda->shared_data.pair.special_coul  = force->special_coul;

        }
          for (int i = 1; i <= n; i++) {
      for (int j = i; j <= n; j++) {
        lj_type_double[i][j] = lj_type[i][j];
        lj_type_double[j][i] = lj_type[i][j];
      }
    }
}

/* ---------------------------------------------------------------------- */

void PairLJSDKCoulLongCuda::compute(int eflag, int vflag)
{
        if (eflag || vflag) ev_setup(eflag,vflag);
        if(eflag) cuda->cu_eng_vdwl->upload();
        if(eflag) cuda->cu_eng_coul->upload();
        if(vflag) cuda->cu_virial->upload();

        Cuda_PairLJSDKCoulLongCuda(& cuda->shared_data, & cuda_neigh_list->sneighlist, eflag, vflag, eflag_atom, vflag_atom);

    if(not cuda->shared_data.pair.collect_forces_later)
    {
          if(eflag) cuda->cu_eng_vdwl->download();
          if(eflag) cuda->cu_eng_coul->download();
          if(vflag) cuda->cu_virial->download();
    }

}

/* ---------------------------------------------------------------------- */

void PairLJSDKCoulLongCuda::settings(int narg, char **arg)
{
        PairLJSDKCoulLong::settings(narg, arg);
        cuda->shared_data.pair.cut_global = (F_CFLOAT) cut_lj_global;
        cuda->shared_data.pair.cut_coul_global = (F_CFLOAT) cut_coul;
}

/* ---------------------------------------------------------------------- */

void PairLJSDKCoulLongCuda::coeff(int narg, char **arg)
{
        PairLJSDKCoulLong::coeff(narg, arg);
        allocate();
}

void PairLJSDKCoulLongCuda::init_style()
{
        MYDBG(printf("# CUDA PairLJSDKCoulLongCuda::init_style start\n"); )
  // request regular or rRESPA neighbor lists

  int irequest;

        irequest = neighbor->request(this,instance_me);
    neighbor->requests[irequest]->full = 1;
    neighbor->requests[irequest]->half = 0;
    neighbor->requests[irequest]->cudable = 1;

  g_ewald = force->kspace->g_ewald;
  cuda->shared_data.pair.g_ewald=g_ewald;
  cuda->shared_data.pppm.qqrd2e=force->qqrd2e;
  if (force->newton) error->warning(FLERR,"Pair style uses does not use \"newton\" setting. You might test if \"newton off\" makes the simulation run faster.");
  MYDBG(printf("# CUDA PairLJSDKCoulLongCuda::init_style end\n"); )
}

void PairLJSDKCoulLongCuda::init_list(int id, NeighList *ptr)
{
        MYDBG(printf("# CUDA PairLJSDKCoulLongCuda::init_list\n");)
        PairLJSDKCoulLong::init_list(id, ptr);
        #ifndef CUDA_USE_BINNING
        // right now we can only handle verlet (id 0), not respa
        if(id == 0) cuda_neigh_list = cuda->registerNeighborList(ptr);
        // see Neighbor::init() for details on lammps lists' logic
        #endif
        MYDBG(printf("# CUDA PairLJSDKCoulLongCuda::init_list end\n");)
}

void PairLJSDKCoulLongCuda::ev_setup(int eflag, int vflag)
{
        int maxeatomold=maxeatom;
        PairLJSDKCoulLong::ev_setup(eflag,vflag);

  if (eflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_eatom; cuda->cu_eatom = new cCudaData<double, ENERGY_CFLOAT, x > ((double*)eatom, & cuda->shared_data.atom.eatom , atom->nmax  );}

  if (vflag_atom && atom->nmax > maxeatomold)
        {delete cuda->cu_vatom; cuda->cu_vatom = new cCudaData<double, ENERGY_CFLOAT, yx > ((double*)vatom, & cuda->shared_data.atom.vatom , atom->nmax, 6  );}

}
