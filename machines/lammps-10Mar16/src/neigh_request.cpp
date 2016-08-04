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

#include "neigh_request.h"
#include "atom.h"
#include "memory.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

NeighRequest::NeighRequest(LAMMPS *lmp) : Pointers(lmp)
{
  // default ID = 0

  id = 0;
  unprocessed = 1;

  // default is pair request

  pair = 1;
  fix = compute = command = 0;

  // default is half neighbor list

  half = 1;
  full = 0;
  full_cluster = 0;
  gran = granhistory = 0;
  respainner = respamiddle = respaouter = 0;
  half_from_full = 0;

  // default is every reneighboring
  // default is use newton_pair setting in force
  // default is encode special bond flags
  // default is no auxiliary floating point values
  // default is no neighbors of ghosts
  // default is no CUDA neighbor list build
  // default is no multi-threaded neighbor list build
  // default is no Kokkos neighbor list build
  // default is no Shardlow Splitting Algorithm (SSA) neighbor list build

  occasional = 0;
  newton = 0;
  special = 1;
  dnum = 0;
  ghost = 0;
  cudable = 0;
  omp = 0;
  intel = 0;
  kokkos_host = kokkos_device = 0;
  ssa = 0;

  // default is no copy or skip

  copy = 0;
  skip = 0;
  iskip = NULL;
  ijskip = NULL;
  otherlist = -1;
}

/* ---------------------------------------------------------------------- */

NeighRequest::~NeighRequest()
{
  delete [] iskip;
  memory->destroy(ijskip);
}

/* ----------------------------------------------------------------------
   archive request params that Neighbor may change after call to identical()
------------------------------------------------------------------------- */

void NeighRequest::archive()
{
  half_original = half;
  half_from_full_original = half_from_full;
  copy_original = copy;
  otherlist_original = otherlist;
}

/* ----------------------------------------------------------------------
   compare this request to other request
   identical means all params set by requestor are the same
   compare to original values in other if Neighbor may have changed them
   return 1 if identical, 0 if not
------------------------------------------------------------------------- */

int NeighRequest::identical(NeighRequest *other)
{
  int same = 1;

  // set same = 0 if old list was never processed
  // use of requestor_instance and instance counter
  //   prevents an old fix from being unfix/refix in same memory location
  //   and appearing to be old, when it is really new
  //   only needed for classes with persistent neigh lists: Fix, Compute, Pair

  if (other->unprocessed) same = 0;

  if (requestor != other->requestor) same = 0;
  if (requestor_instance != other->requestor_instance) same = 0;
  if (id != other->id) same = 0;

  if (pair != other->pair) same = 0;
  if (fix != other->fix) same = 0;
  if (compute != other->compute) same = 0;
  if (command != other->command) same = 0;

  if (half != other->half_original) same = 0;
  if (full != other->full) same = 0;
  if (gran != other->gran) same = 0;
  if (granhistory != other->granhistory) same = 0;
  if (respainner != other->respainner) same = 0;
  if (respamiddle != other->respamiddle) same = 0;
  if (respaouter != other->respaouter) same = 0;
  if (half_from_full != other->half_from_full_original) same = 0;

  if (newton != other->newton) same = 0;
  if (occasional != other->occasional) same = 0;
  if (special != other->special) same = 0;
  if (dnum != other->dnum) same = 0;
  if (ghost != other->ghost) same = 0;
  if (cudable != other->cudable) same = 0;
  if (omp != other->omp) same = 0;
  if (intel != other->intel) same = 0;
  if (ssa != other->ssa) same = 0;

  if (copy != other->copy_original) same = 0;
  if (same_skip(other) == 0) same = 0;
  if (otherlist != other->otherlist_original) same = 0;

  return same;
}

/* ----------------------------------------------------------------------
   compare kind of this request to other request
   return 1 if same, 0 if different
------------------------------------------------------------------------- */

int NeighRequest::same_kind(NeighRequest *other)
{
  int same = 1;

  if (half != other->half) same = 0;
  if (full != other->full) same = 0;
  if (gran != other->gran) same = 0;
  if (granhistory != other->granhistory) same = 0;
  if (respainner != other->respainner) same = 0;
  if (respamiddle != other->respamiddle) same = 0;
  if (respaouter != other->respaouter) same = 0;
  if (half_from_full != other->half_from_full) same = 0;
  if (newton != other->newton) same = 0;
  if (ghost != other->ghost) same = 0;
  if (cudable != other->cudable) same = 0;
  if (omp != other->omp) same = 0;
  if (intel != other->intel) same = 0;
  if (ssa != other->ssa) same = 0;

  return same;
}

/* ----------------------------------------------------------------------
   compare skip attributes of this request to other request
   return 1 if same, 0 if different
------------------------------------------------------------------------- */

int NeighRequest::same_skip(NeighRequest *other)
{
  int i,j;

  int same = 1;

  if (skip != other->skip) same = 0;
  if (skip && other->skip) {
    int ntypes = atom->ntypes;
    for (i = 1; i <= ntypes; i++)
      if (iskip[i] != other->iskip[i]) same = 0;
    for (i = 1; i <= ntypes; i++)
      for (j = 1; j <= ntypes; j++)
        if (ijskip[i][j] != other->ijskip[i][j]) same = 0;
  }

  return same;
}

/* ----------------------------------------------------------------------
   set kind and other values of this request to that of other request
------------------------------------------------------------------------- */

void NeighRequest::copy_request(NeighRequest *other)
{
  half = 0;

  if (other->half) half = 1;
  if (other->full) full = 1;
  if (other->gran) gran = 1;
  if (other->granhistory) granhistory = 1;
  if (other->respainner) respainner = 1;
  if (other->respamiddle) respamiddle = 1;
  if (other->respaouter) respaouter = 1;
  if (other->half_from_full) half_from_full = 1;

  newton = other->newton;
  dnum = other->dnum;
  ghost = other->ghost;
  cudable = other->cudable;
  omp = other->omp;
  intel = other->intel;
  ssa = other->ssa;
}
