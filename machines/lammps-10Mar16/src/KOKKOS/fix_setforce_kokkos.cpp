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

#include <string.h>
#include <stdlib.h>
#include "fix_setforce_kokkos.h"
#include "atom_kokkos.h"
#include "update.h"
#include "modify.h"
#include "domain.h"
#include "region.h"
#include "respa.h"
#include "input.h"
#include "variable.h"
#include "memory.h"
#include "error.h"
#include "force.h"
#include "atom_masks.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{NONE,CONSTANT,EQUAL,ATOM};

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSetForceKokkos<DeviceType>::FixSetForceKokkos(LAMMPS *lmp, int narg, char **arg) :
  FixSetForce(lmp, narg, arg)
{
  atomKK = (AtomKokkos *) atom;
  execution_space = ExecutionSpaceFromDevice<DeviceType>::space;
  datamask_read = EMPTY_MASK;
  datamask_modify = EMPTY_MASK;

  memory->destroy(sforce);
  memory->create_kokkos(k_sforce,sforce,maxatom,3,"setforce:sforce");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
FixSetForceKokkos<DeviceType>::~FixSetForceKokkos()
{
  if (copymode) return;

  memory->destroy_kokkos(k_sforce,sforce);
  sforce = NULL;
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSetForceKokkos<DeviceType>::init()
{
  FixSetForce::init();

  if (strstr(update->integrate_style,"respa"))
    error->all(FLERR,"Cannot (yet) use respa with Kokkos");
}

/* ---------------------------------------------------------------------- */

template<class DeviceType>
void FixSetForceKokkos<DeviceType>::post_force(int vflag)
{
  atomKK->sync(execution_space, X_MASK | F_MASK | MASK_MASK);

  x = atomKK->k_x.view<DeviceType>();
  f = atomKK->k_f.view<DeviceType>();
  mask = atomKK->k_mask.view<DeviceType>();

  int nlocal = atom->nlocal;

  // update region if necessary

  region = NULL;
  if (iregion >= 0) {
    region = domain->regions[iregion];
    region->prematch();
    d_match = DAT::t_int_1d("setforce:d_match",nlocal);
    region->match_all_kokkos(groupbit,d_match);
  }

  // reallocate sforce array if necessary

  if (varflag == ATOM && nlocal > maxatom) {
    maxatom = atom->nmax;
    memory->destroy_kokkos(k_sforce,sforce);
    memory->create_kokkos(k_sforce,sforce,maxatom,3,"setforce:sforce");
  }

  foriginal[0] = foriginal[1] = foriginal[2] = 0.0;
  double_3 foriginal_kk;
  force_flag = 0;

  if (varflag == CONSTANT) {
    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixSetForceConstant>(0,nlocal),*this,foriginal_kk);
    DeviceType::fence();
    copymode = 0;

  // variable force, wrap with clear/add

  } else {

    atomKK->sync(Host,ALL_MASK); // this can be removed when variable class is ported to Kokkos

    modify->clearstep_compute();

    if (xstyle == EQUAL) xvalue = input->variable->compute_equal(xvar);
    else if (xstyle == ATOM)
      input->variable->compute_atom(xvar,igroup,&sforce[0][0],3,0);
    if (ystyle == EQUAL) yvalue = input->variable->compute_equal(yvar);
    else if (ystyle == ATOM)
      input->variable->compute_atom(yvar,igroup,&sforce[0][1],3,0);
    if (zstyle == EQUAL) zvalue = input->variable->compute_equal(zvar);
    else if (zstyle == ATOM)
      input->variable->compute_atom(zvar,igroup,&sforce[0][2],3,0);

    modify->addstep_compute(update->ntimestep + 1);

    if (varflag == ATOM) {  // this can be removed when variable class is ported to Kokkos
      k_sforce.modify<LMPHostType>();
      k_sforce.sync<DeviceType>();
    }

    copymode = 1;
    Kokkos::parallel_reduce(Kokkos::RangePolicy<DeviceType, TagFixSetForceNonConstant>(0,nlocal),*this,foriginal_kk);
    DeviceType::fence();
    copymode = 0;
  }

  atomKK->modified(execution_space, F_MASK);

  foriginal[0] = foriginal_kk.d0;
  foriginal[1] = foriginal_kk.d1;
  foriginal[2] = foriginal_kk.d2;
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixSetForceKokkos<DeviceType>::operator()(TagFixSetForceConstant, const int &i, double_3& foriginal_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    foriginal_kk.d0 += f(i,0);
    foriginal_kk.d1 += f(i,1);
    foriginal_kk.d2 += f(i,2);
    if (xstyle) f(i,0) = xvalue;
    if (ystyle) f(i,1) = yvalue;
    if (zstyle) f(i,2) = zvalue;
  }
}

template<class DeviceType>
KOKKOS_INLINE_FUNCTION
void FixSetForceKokkos<DeviceType>::operator()(TagFixSetForceNonConstant, const int &i, double_3& foriginal_kk) const {
  if (mask[i] & groupbit) {
    if (region && !d_match[i]) return;
    foriginal_kk.d0 += f(i,0);
    foriginal_kk.d1 += f(i,1);
    foriginal_kk.d2 += f(i,2);
    if (xstyle == ATOM) f(i,0) = d_sforce(i,0);
    else if (xstyle) f(i,0) = xvalue;
    if (ystyle == ATOM) f(i,1) = d_sforce(i,1);
    else if (ystyle) f(i,1) = yvalue;
    if (zstyle == ATOM) f(i,2) = d_sforce(i,2);
    else if (zstyle) f(i,2) = zvalue;
  }
}

template class FixSetForceKokkos<LMPDeviceType>;
#ifdef KOKKOS_HAVE_CUDA
template class FixSetForceKokkos<LMPHostType>;
#endif