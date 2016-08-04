/* -*- c++ -*- ----------------------------------------------------------
   LAMMPS - Large-scale Atomic/Molecular Massively Parallel Simulator
   http://lammps.sandia.gov, Sandia National Laboratories
   Steve Plimpton, sjplimp@sandia.gov

   Copyright (2003) Sandia Corporation.  Under the terms of Contract
   DE-AC04-94AL85000 with Sandia Corporation, the U.S. Government retains
   certain rights in this software.  This software is distributed under
   the GNU General Public License.

   See the README file in the top-level LAMMPS directory.
------------------------------------------------------------------------- */

#ifndef LMP_NEIGHBOR_OMP_H
#define LMP_NEIGHBOR_OMP_H

#if defined(_OPENMP)
#include <omp.h>
#endif

#include "modify.h"
#include "timer.h"
#include "fix_omp.h"
#include "thr_data.h"

namespace LAMMPS_NS {

// these macros hide some ugly and redundant OpenMP related stuff
#if defined(_OPENMP)

// make sure we have at least one page for each thread
#define NEIGH_OMP_INIT                             \
  const int nthreads = comm->nthreads;             \
  const int ifix = modify->find_fix("package_omp")

// get thread id and then assign each thread a fixed chunk of atoms
#define NEIGH_OMP_SETUP(num)                    \
  {                                             \
    const int tid = omp_get_thread_num();       \
    const int idelta = 1 + num/nthreads;        \
    const int ifrom = tid*idelta;               \
    const int ito   = ((ifrom + idelta) > num)  \
      ? num : (ifrom+idelta);                   \
    FixOMP *fix = static_cast<FixOMP *>(modify->fix[ifix]); \
    ThrData *thr = fix->get_thr(tid);           \
    thr->timer(Timer::START);

#define NEIGH_OMP_CLOSE                         \
      thr->timer(Timer::NEIGH);                 \
    }

#else /* !defined(_OPENMP) */

#define NEIGH_OMP_INIT

#define NEIGH_OMP_SETUP(num)                    \
  const int tid = 0;                            \
  const int ifrom = 0;                          \
  const int ito = num

#define NEIGH_OMP_CLOSE

#endif

}

#endif
