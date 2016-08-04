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
#include "fix_balance.h"
#include "balance.h"
#include "update.h"
#include "atom.h"
#include "comm.h"
#include "domain.h"
#include "neighbor.h"
#include "irregular.h"
#include "force.h"
#include "kspace.h"
#include "rcb.h"
#include "error.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{SHIFT,BISECTION};
enum{LAYOUT_UNIFORM,LAYOUT_NONUNIFORM,LAYOUT_TILED};    // several files

/* ---------------------------------------------------------------------- */

FixBalance::FixBalance(LAMMPS *lmp, int narg, char **arg) :
  Fix(lmp, narg, arg)
{
  if (narg < 6) error->all(FLERR,"Illegal fix balance command");

  box_change_domain = 1;
  scalar_flag = 1;
  extscalar = 0;
  vector_flag = 1;
  size_vector = 3;
  extvector = 0;
  global_freq = 1;

  // parse arguments

  int dimension = domain->dimension;

  nevery = force->inumeric(FLERR,arg[3]);
  if (nevery < 0) error->all(FLERR,"Illegal fix balance command");
  thresh = force->numeric(FLERR,arg[4]);

  if (strcmp(arg[5],"shift") == 0) lbstyle = SHIFT;
  else if (strcmp(arg[5],"rcb") == 0) lbstyle = BISECTION;
  else error->all(FLERR,"Illegal fix balance command");

  int iarg = 5;
  if (lbstyle == SHIFT) {
    if (iarg+4 > narg) error->all(FLERR,"Illegal fix balance command");
    if (strlen(arg[iarg+1]) > 3) error->all(FLERR,"Illegal fix balance command");
    strcpy(bstr,arg[iarg+1]);
    nitermax = force->inumeric(FLERR,arg[iarg+2]);
    if (nitermax <= 0) error->all(FLERR,"Illegal fix balance command");
    stopthresh = force->numeric(FLERR,arg[iarg+3]);
    if (stopthresh < 1.0) error->all(FLERR,"Illegal fix balance command");
    iarg += 4;
  } else if (lbstyle == BISECTION) {
    iarg++;
  }

  // optional args

  outflag = 0;
  int outarg = 0;
  fp = NULL;

  while (iarg < narg) {
    if (strcmp(arg[iarg],"out") == 0) {
      if (iarg+2 > narg) error->all(FLERR,"Illegal fix balance command");
      outflag = 1;
      outarg = iarg+1;
      iarg += 2;
    } else error->all(FLERR,"Illegal fix balance command");
  }

  // error check

  if (lbstyle == SHIFT) {
    int blen = strlen(bstr);
    for (int i = 0; i < blen; i++) {
      if (bstr[i] != 'x' && bstr[i] != 'y' && bstr[i] != 'z')
        error->all(FLERR,"Fix balance shift string is invalid");
      if (bstr[i] == 'z' && dimension == 2)
        error->all(FLERR,"Fix balance shift string is invalid");
      for (int j = i+1; j < blen; j++)
        if (bstr[i] == bstr[j])
          error->all(FLERR,"Fix balance shift string is invalid");
    }
  }

  if (lbstyle == BISECTION && comm->style == 0)
    error->all(FLERR,"Fix balance rcb cannot be used with comm_style brick");

  // create instance of Balance class
  // if SHIFT, initialize it with params

  balance = new Balance(lmp);
  if (lbstyle == SHIFT) balance->shift_setup(bstr,nitermax,thresh);

  // create instance of Irregular class

  irregular = new Irregular(lmp);

  // output file

  if (outflag && comm->me == 0) {
    fp = fopen(arg[outarg],"w");
    if (fp == NULL) error->one(FLERR,"Cannot open fix balance output file");
  }

  // only force reneighboring if nevery > 0

  if (nevery) force_reneighbor = 1;

  // compute initial outputs

  imbfinal = imbprev = balance->imbalance_nlocal(maxperproc);
  itercount = 0;
  pending = 0;
}

/* ---------------------------------------------------------------------- */

FixBalance::~FixBalance()
{
  if (fp) fclose(fp);
  delete balance;
  delete irregular;
}

/* ---------------------------------------------------------------------- */

int FixBalance::setmask()
{
  int mask = 0;
  mask |= PRE_EXCHANGE;
  mask |= PRE_NEIGHBOR;
  return mask;
}

/* ---------------------------------------------------------------------- */

void FixBalance::init()
{
  if (force->kspace) kspace_flag = 1;
  else kspace_flag = 0;
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup(int vflag)
{
  // compute final imbalance factor if setup_pre_exchange() invoked balancer
  // this is called at end of run setup, before output

  pre_neighbor();
}

/* ---------------------------------------------------------------------- */

void FixBalance::setup_pre_exchange()
{
  // insure atoms are in current box & update box via shrink-wrap
  // has to be be done before rebalance() invokes Irregular::migrate_atoms()
  //   since it requires atoms be inside simulation box
  //   even though pbc() will be done again in Verlet::run()
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // perform a rebalance if threshhold exceeded

  imbnow = balance->imbalance_nlocal(maxperproc);
  if (imbnow > thresh) rebalance();

  // next_reneighbor = next time to force reneighboring

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::pre_exchange()
{
  // return if not a rebalance timestep

  if (nevery && update->ntimestep < next_reneighbor) return;

  // insure atoms are in current box & update box via shrink-wrap
  // no exchange() since doesn't matter if atoms are assigned to correct procs

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  domain->pbc();
  domain->reset_box();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // return if imbalance < threshhold

  imbnow = balance->imbalance_nlocal(maxperproc);
  if (imbnow <= thresh) {
    if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
    return;
  }

  rebalance();

  // next timestep to rebalance

  if (nevery) next_reneighbor = (update->ntimestep/nevery)*nevery + nevery;
}

/* ----------------------------------------------------------------------
   compute final imbalance factor based on nlocal after comm->exchange()
   only do this if rebalancing just occured
------------------------------------------------------------------------- */

void FixBalance::pre_neighbor()
{
  if (!pending) return;
  imbfinal = balance->imbalance_nlocal(maxperproc);
  pending = 0;
}

/* ----------------------------------------------------------------------
   perform dynamic load balancing
------------------------------------------------------------------------- */

void FixBalance::rebalance()
{
  imbprev = imbnow;

  // invoke balancer and reset comm->uniform flag

  int *sendproc;
  if (lbstyle == SHIFT) {
    itercount = balance->shift();
    comm->layout = LAYOUT_NONUNIFORM;
  } else if (lbstyle == BISECTION) {
    sendproc = balance->bisection();
    comm->layout = LAYOUT_TILED;
  }

  // output of new decomposition

  if (outflag) balance->dumpout(update->ntimestep,fp);

  // reset proc sub-domains
  // check and warn if any proc's subbox is smaller than neigh skin
  //   since may lead to lost atoms in exchange()

  if (domain->triclinic) domain->set_lamda_box();
  domain->set_local_box();
  domain->subbox_too_small_check(neighbor->skin);

  // move atoms to new processors via irregular()
  // only needed if migrate_check() says an atom moves to far
  // else allow caller's comm->exchange() to do it

  if (domain->triclinic) domain->x2lamda(atom->nlocal);
  if (lbstyle == BISECTION) irregular->migrate_atoms(0,1,sendproc);
  else if (irregular->migrate_check()) irregular->migrate_atoms();
  if (domain->triclinic) domain->lamda2x(atom->nlocal);

  // invoke KSpace setup_grid() to adjust to new proc sub-domains

  if (kspace_flag) force->kspace->setup_grid();

  // pending triggers pre_neighbor() to compute final imbalance factor
  // can only be done after atoms migrate in caller's comm->exchange()

  pending = 1;
}

/* ----------------------------------------------------------------------
   return imbalance factor after last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_scalar()
{
  return imbfinal;
}

/* ----------------------------------------------------------------------
   return stats for last rebalance
------------------------------------------------------------------------- */

double FixBalance::compute_vector(int i)
{
  if (i == 0) return (double) maxperproc;
  if (i == 1) return (double) itercount;
  return imbprev;
}

/* ----------------------------------------------------------------------
   return # of bytes of allocated memory
------------------------------------------------------------------------- */

double FixBalance::memory_usage()
{
  double bytes = irregular->memory_usage();
  if (balance->rcb) bytes += balance->rcb->memory_usage();
  return bytes;
}
