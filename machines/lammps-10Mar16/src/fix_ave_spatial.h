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

#ifdef FIX_CLASS

FixStyle(ave/spatial,FixAveSpatial)

#else

#ifndef LMP_FIX_AVE_SPATIAL_H
#define LMP_FIX_AVE_SPATIAL_H

#include <stdio.h>
#include "fix.h"

namespace LAMMPS_NS {

class FixAveSpatial : public Fix {
 public:
  FixAveSpatial(class LAMMPS *, int, char **);
  ~FixAveSpatial();
  int setmask();
  void init();
  void setup(int);
  void end_of_step();
  double compute_array(int,int);
  double memory_usage();

 private:
  int me,nvalues;
  int nrepeat,nfreq,irepeat;
  bigint nvalid,nvalid_last;
  int ndim,normflag,regionflag,overwrite,discard;
  char *tstring,*sstring,*idregion;
  int *which,*argindex,*value2index;
  char **ids;
  FILE *fp;
  class Region *region;

  int ave,nwindow,scaleflag;
  int norm,iwindow,window_limit;
  double xscale,yscale,zscale;
  double bin_volume;

  long filepos;
  int dim[3],originflag[3],nlayers[3];
  double origin[3],delta[3];
  double offset[3],invdelta[3];

  int maxvar;
  double *varatom;

  int maxatom;
  int *bin;

  int nbins,maxbin;
  double **coord;
  double *count_one,*count_many,*count_sum;
  double **values_one,**values_many,**values_sum;
  double *count_total,**count_list;
  double **values_total,***values_list;

  int minflag[3],maxflag[3];
  double minvalue[3],maxvalue[3];

  void setup_bins();
  void atom2bin1d();
  void atom2bin2d();
  void atom2bin3d();
  bigint nextvalid();
};

}

#endif
#endif

/* ERROR/WARNING messages:

W: The fix ave/spatial command has been replaced by the more flexible fix ave/chunk and compute chunk/atom commands -- fix ave/spatial will be removed in the summer of 2015

Self-explanatory.

E: Illegal ... command

Self-explanatory.  Check the input script syntax and compare to the
documentation for the command.  You can use -echo screen as a
command-line option when running LAMMPS to see the offending line.

E: Cannot use fix ave/spatial z for 2 dimensional model

Self-explanatory.

E: Same dimension twice in fix ave/spatial

Self-explanatory.

E: No input values for fix ave/spatial

Self-explanatory.

E: Region ID for fix ave/spatial does not exist

Self-explanatory.

E: Cannot open fix ave/spatial file %s

The specified file cannot be opened.  Check that the path and name are
correct.

E: Compute ID for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial compute does not calculate per-atom values

A compute used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial compute does not calculate a per-atom vector

A compute used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial compute does not calculate a per-atom array

Self-explanatory.

E: Fix ave/spatial compute vector is accessed out-of-range

The index for the vector is out of bounds.

E: Fix ID for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial fix does not calculate per-atom values

A fix used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial fix does not calculate a per-atom vector

A fix used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial fix does not calculate a per-atom array

Self-explanatory.

E: Fix ave/spatial fix vector is accessed out-of-range

The index for the vector is out of bounds.

E: Variable name for fix ave/spatial does not exist

Self-explanatory.

E: Fix ave/spatial variable is not atom-style variable

A variable used by fix ave/spatial must generate per-atom values.

E: Fix ave/spatial for triclinic boxes requires units reduced

Self-explanatory.

E: Fix ave/spatial settings invalid with changing box size

If the box size changes, only the units reduced option can be
used.

E: Fix for fix ave/spatial not computed at compatible time

Fixes generate their values on specific timesteps.  Fix ave/spatial is
requesting a value on a non-allowed timestep.

E: Invalid timestep reset for fix ave/spatial

Resetting the timestep has invalidated the sequence of timesteps this
fix needs to process.

E: Invalid bin bounds in fix ave/spatial

The lo/hi values are inconsistent.

*/
