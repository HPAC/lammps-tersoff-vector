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

#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <ctype.h>
#include "kokkos.h"
#include "lammps.h"
#include "force.h"
#include "neighbor_kokkos.h"
#include "neigh_list_kokkos.h"
#include "error.h"

using namespace LAMMPS_NS;

/* ---------------------------------------------------------------------- */

KokkosLMP::KokkosLMP(LAMMPS *lmp, int narg, char **arg) : Pointers(lmp)
{
  kokkos_exists = 1;
  lmp->kokkos = this;

  int me = 0;
  MPI_Comm_rank(world,&me);
  if (me == 0) error->message(FLERR,"KOKKOS mode is enabled");

  // process any command-line args that invoke Kokkos settings

  ngpu = 0;
  int device = 0;
  num_threads = 1;
  numa = 1;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"d") == 0 || strcmp(arg[iarg],"device") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Invalid Kokkos command-line args");
      device = atoi(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"g") == 0 ||
               strcmp(arg[iarg],"gpus") == 0) {
#ifndef KOKKOS_HAVE_CUDA
      error->all(FLERR,"GPUs are requested but Kokkos has not been compiled for CUDA");
#endif
      if (iarg+2 > narg) error->all(FLERR,"Invalid Kokkos command-line args");
      ngpu = atoi(arg[iarg+1]);

      int skip_gpu = 9999;
      if (iarg+2 < narg && isdigit(arg[iarg+2][0])) {
        skip_gpu = atoi(arg[iarg+2]);
        iarg++;
      }
      iarg += 2;

      char *str;
      if ((str = getenv("SLURM_LOCALID"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpu;
        if (device >= skip_gpu) device++;
      }
      if ((str = getenv("MV2_COMM_WORLD_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpu;
        if (device >= skip_gpu) device++;
      }
      if ((str = getenv("OMPI_COMM_WORLD_LOCAL_RANK"))) {
        int local_rank = atoi(str);
        device = local_rank % ngpu;
        if (device >= skip_gpu) device++;
      }

    } else if (strcmp(arg[iarg],"t") == 0 ||
               strcmp(arg[iarg],"threads") == 0) {
      num_threads = atoi(arg[iarg+1]);
      iarg += 2;

    } else if (strcmp(arg[iarg],"n") == 0 ||
               strcmp(arg[iarg],"numa") == 0) {
      numa = atoi(arg[iarg+1]);
      iarg += 2;

    } else error->all(FLERR,"Invalid Kokkos command-line args");
  }

  // initialize Kokkos

  if (me == 0) {
    if (screen) fprintf(screen,"  using %d GPU(s)\n",ngpu);
    if (logfile) fprintf(logfile,"  using %d GPU(s)\n",ngpu);
  }

#ifdef KOKKOS_HAVE_CUDA
  if (ngpu <= 0)
    error->all(FLERR,"Kokkos has been compiled for CUDA but no GPUs are requested");

  Kokkos::HostSpace::execution_space::initialize(num_threads,numa);
  Kokkos::Cuda::SelectDevice select_device(device);
  Kokkos::Cuda::initialize(select_device);
#else
  LMPHostType::initialize(num_threads,numa);
#endif

  // default settings for package kokkos command

  neighflag = FULL;
  exchange_comm_classic = 0;
  forward_comm_classic = 0;
  exchange_comm_on_host = 1;
  forward_comm_on_host = 1;
}

/* ---------------------------------------------------------------------- */

KokkosLMP::~KokkosLMP()
{
  // finalize Kokkos

#ifdef KOKKOS_HAVE_CUDA
  Kokkos::Cuda::finalize();
  Kokkos::HostSpace::execution_space::finalize();
#else
  LMPHostType::finalize();
#endif
}

/* ----------------------------------------------------------------------
   invoked by package kokkos command
------------------------------------------------------------------------- */

void KokkosLMP::accelerator(int narg, char **arg)
{
  // defaults

  neighflag = FULL;
  int newtonflag = 0;
  double binsize = 0.0;
  exchange_comm_classic = forward_comm_classic = 0;
  exchange_comm_on_host = forward_comm_on_host = 0;

  int iarg = 0;
  while (iarg < narg) {
    if (strcmp(arg[iarg],"neigh") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"full") == 0) neighflag = FULL;
      else if (strcmp(arg[iarg+1],"half/thread") == 0) neighflag = HALFTHREAD;
      else if (strcmp(arg[iarg+1],"half") == 0) neighflag = HALF;
      else if (strcmp(arg[iarg+1],"n2") == 0) neighflag = N2;
      else if (strcmp(arg[iarg+1],"full/cluster") == 0) neighflag = FULLCLUSTER;
      else error->all(FLERR,"Illegal package kokkos command");
      if (neighflag == HALF && (num_threads > 1 || ngpu > 0))
        error->all(FLERR,"Must use Kokkos half/thread or full neighbor list with threads or GPUs");
      iarg += 2;
    } else if (strcmp(arg[iarg],"binsize") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      binsize = force->numeric(FLERR,arg[iarg+1]);
      iarg += 2;
    } else if (strcmp(arg[iarg],"newton") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"off") == 0) newtonflag = 0;
      else if (strcmp(arg[iarg+1],"on") == 0) newtonflag = 1;
      else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"no") == 0) {
        exchange_comm_classic = forward_comm_classic = 1;
      } else if (strcmp(arg[iarg+1],"host") == 0) {
        exchange_comm_classic = forward_comm_classic = 0;
        exchange_comm_on_host = forward_comm_on_host = 1;
      } else if (strcmp(arg[iarg+1],"device") == 0) {
        exchange_comm_classic = forward_comm_classic = 0;
        exchange_comm_on_host = forward_comm_on_host = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm/exchange") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"no") == 0) exchange_comm_classic = 1;
      else if (strcmp(arg[iarg+1],"host") == 0) {
        exchange_comm_classic = 0;
        exchange_comm_on_host = 1;
      } else if (strcmp(arg[iarg+1],"device") == 0) {
        exchange_comm_classic = 0;
        exchange_comm_on_host = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else if (strcmp(arg[iarg],"comm/forward") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal package kokkos command");
      if (strcmp(arg[iarg+1],"no") == 0) forward_comm_classic = 1;
      else if (strcmp(arg[iarg+1],"host") == 0) {
        forward_comm_classic = 0;
        forward_comm_on_host = 1;
      } else if (strcmp(arg[iarg+1],"device") == 0) {
        forward_comm_classic = 0;
        forward_comm_on_host = 0;
      } else error->all(FLERR,"Illegal package kokkos command");
      iarg += 2;
    } else error->all(FLERR,"Illegal package kokkos command");
  }

  // set newton flags
  // set neighbor binsize, same as neigh_modify command

  force->newton = force->newton_pair = force->newton_bond = newtonflag;

  neighbor->binsize_user = binsize;
  if (binsize <= 0.0) neighbor->binsizeflag = 0;
  else neighbor->binsizeflag = 1;
}

/* ----------------------------------------------------------------------
   called by Finish
------------------------------------------------------------------------- */

int KokkosLMP::neigh_list_kokkos(int m)
{
  NeighborKokkos *nk = (NeighborKokkos *) neighbor;
  if (nk->lists_host[m] && nk->lists_host[m]->d_numneigh.dimension_0())
    return 1;
  if (nk->lists_device[m] && nk->lists_device[m]->d_numneigh.dimension_0())
    return 1;
  return 0;
}

/* ----------------------------------------------------------------------
   called by Finish
------------------------------------------------------------------------- */

int KokkosLMP::neigh_count(int m)
{
  int inum;
  int nneigh = 0;

  ArrayTypes<LMPHostType>::t_int_1d h_ilist;
  ArrayTypes<LMPHostType>::t_int_1d h_numneigh;

  NeighborKokkos *nk = (NeighborKokkos *) neighbor;
  if (nk->lists_host[m]) {
    inum = nk->lists_host[m]->inum;
#ifndef KOKKOS_USE_CUDA_UVM
    h_ilist = Kokkos::create_mirror_view(nk->lists_host[m]->d_ilist);
    h_numneigh = Kokkos::create_mirror_view(nk->lists_host[m]->d_numneigh);
#else
    h_ilist = nk->lists_host[m]->d_ilist;
    h_numneigh = nk->lists_host[m]->d_numneigh;
#endif
    Kokkos::deep_copy(h_ilist,nk->lists_host[m]->d_ilist);
    Kokkos::deep_copy(h_numneigh,nk->lists_host[m]->d_numneigh);
  } else if (nk->lists_device[m]) {
    inum = nk->lists_device[m]->inum;
#ifndef KOKKOS_USE_CUDA_UVM
    h_ilist = Kokkos::create_mirror_view(nk->lists_device[m]->d_ilist);
    h_numneigh = Kokkos::create_mirror_view(nk->lists_device[m]->d_numneigh);
#else
    h_ilist = nk->lists_device[m]->d_ilist;
    h_numneigh = nk->lists_device[m]->d_numneigh;
#endif
    Kokkos::deep_copy(h_ilist,nk->lists_device[m]->d_ilist);
    Kokkos::deep_copy(h_numneigh,nk->lists_device[m]->d_numneigh);
  }

  for (int i = 0; i < inum; i++) nneigh += h_numneigh[h_ilist[i]];

  return nneigh;
}
