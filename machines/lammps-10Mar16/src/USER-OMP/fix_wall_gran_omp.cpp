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
   Contributing author: Axel Kohlmeyer (Temple U)
------------------------------------------------------------------------- */

#include <math.h>
#include "fix_wall_gran_omp.h"
#include "atom.h"
#include "update.h"

using namespace LAMMPS_NS;
using namespace FixConst;

enum{XPLANE=0,YPLANE=1,ZPLANE=2,ZCYLINDER};    // XYZ PLANE need to be 0,1,2
enum{HOOKE,HOOKE_HISTORY,HERTZ_HISTORY};

#define BIG 1.0e20

/* ---------------------------------------------------------------------- */

FixWallGranOMP::FixWallGranOMP(LAMMPS *lmp, int narg, char **arg) :
  FixWallGran(lmp, narg, arg) { }

/* ---------------------------------------------------------------------- */

void FixWallGranOMP::post_force(int vflag)
{
  double vwall[3];

  // set position of wall to initial settings and velocity to 0.0
  // if wiggle or shear, set wall position and velocity accordingly

  double wlo = lo;
  double whi = hi;
  vwall[0] = vwall[1] = vwall[2] = 0.0;
  if (wiggle) {
    double arg = omega * (update->ntimestep - time_origin) * dt;
    if (wallstyle == axis) {
      wlo = lo + amplitude - amplitude*cos(arg);
      whi = hi + amplitude - amplitude*cos(arg);
    }
    vwall[axis] = amplitude*omega*sin(arg);
  } else if (wshear) vwall[axis] = vshear;

  // loop over all my atoms
  // rsq = distance from wall
  // dx,dy,dz = signed distance from wall
  // for rotating cylinder, reset vwall based on particle position
  // skip atom if not close enough to wall
  //   if wall was set to NULL, it's skipped since lo/hi are infinity
  // compute force and torque on atom if close enough to wall
  //   via wall potential matched to pair potential
  // set shear if pair potential stores history

  double * const * const x = atom->x;
  double * const * const v = atom->v;
  double * const * const f = atom->f;
  double * const * const omega = atom->omega;
  double * const * const torque = atom->torque;
  double * const radius = atom->radius;
  double * const rmass = atom->rmass;
  const int * const mask = atom->mask;
  const int nlocal = atom->nlocal;

  shearupdate = (update->setupflag) ? 0 : 1;

  int i;
#if defined(_OPENMP)
#pragma omp parallel for private(i) default(none) firstprivate(vwall,wlo,whi)
#endif
  for (i = 0; i < nlocal; i++) {

    if (mask[i] & groupbit) {
      double dx,dy,dz,del1,del2,delxy,delr,rsq;

      dx = dy = dz = 0.0;

      if (wallstyle == XPLANE) {
        del1 = x[i][0] - wlo;
        del2 = whi - x[i][0];
        if (del1 < del2) dx = del1;
        else dx = -del2;
      } else if (wallstyle == YPLANE) {
        del1 = x[i][1] - wlo;
        del2 = whi - x[i][1];
        if (del1 < del2) dy = del1;
        else dy = -del2;
      } else if (wallstyle == ZPLANE) {
        del1 = x[i][2] - wlo;
        del2 = whi - x[i][2];
        if (del1 < del2) dz = del1;
        else dz = -del2;
      } else if (wallstyle == ZCYLINDER) {
        delxy = sqrt(x[i][0]*x[i][0] + x[i][1]*x[i][1]);
        delr = cylradius - delxy;
        if (delr > radius[i]) dz = cylradius;
        else {
          dx = -delr/delxy * x[i][0];
          dy = -delr/delxy * x[i][1];
          if (wshear && axis != 2) {
            vwall[0] = vshear * x[i][1]/delxy;
            vwall[1] = -vshear * x[i][0]/delxy;
            vwall[2] = 0.0;
          }
        }
      }

      rsq = dx*dx + dy*dy + dz*dz;

      if (rsq > radius[i]*radius[i]) {
        if (pairstyle != HOOKE) {
          shear[i][0] = 0.0;
          shear[i][1] = 0.0;
          shear[i][2] = 0.0;
        }
      } else {
        if (pairstyle == HOOKE)
          hooke(rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],
                radius[i],rmass[i]);
        else if (pairstyle == HOOKE_HISTORY)
          hooke_history(rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],
                        radius[i],rmass[i],shear[i]);
        else if (pairstyle == HERTZ_HISTORY)
          hertz_history(rsq,dx,dy,dz,vwall,v[i],f[i],omega[i],torque[i],
                        radius[i],rmass[i],shear[i]);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */

void FixWallGranOMP::post_force_respa(int vflag, int ilevel, int iloop)
{
  if (ilevel == nlevels_respa-1) post_force(vflag);
}
