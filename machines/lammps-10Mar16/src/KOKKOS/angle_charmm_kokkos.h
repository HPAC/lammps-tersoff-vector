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

#ifdef ANGLE_CLASS

AngleStyle(charmm/kk,AngleCharmmKokkos<LMPDeviceType>)
AngleStyle(charmm/kk/device,AngleCharmmKokkos<LMPDeviceType>)
AngleStyle(charmm/kk/host,AngleCharmmKokkos<LMPHostType>)

#else

#ifndef LMP_ANGLE_CHARMM_KOKKOS_H
#define LMP_ANGLE_CHARMM_KOKKOS_H

#include "angle_charmm.h"
#include "kokkos_type.h"

namespace LAMMPS_NS {

template<int NEWTON_BOND, int EVFLAG>
struct TagAngleCharmmCompute{};

template<class DeviceType>
class AngleCharmmKokkos : public AngleCharmm {
 public:
  typedef DeviceType device_type;
  typedef EV_FLOAT value_type;

  AngleCharmmKokkos(class LAMMPS *);
  virtual ~AngleCharmmKokkos();
  virtual void compute(int, int);
  virtual void coeff(int, char **);

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleCharmmCompute<NEWTON_BOND,EVFLAG>, const int&, EV_FLOAT&) const;

  template<int NEWTON_BOND, int EVFLAG>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagAngleCharmmCompute<NEWTON_BOND,EVFLAG>, const int&) const;

  //template<int NEWTON_BOND>
  KOKKOS_INLINE_FUNCTION
  void ev_tally(EV_FLOAT &ev, const int i, const int j, const int k,
                     F_FLOAT &eangle, F_FLOAT *f1, F_FLOAT *f3,
                     const F_FLOAT &delx1, const F_FLOAT &dely1, const F_FLOAT &delz1,
                     const F_FLOAT &delx2, const F_FLOAT &dely2, const F_FLOAT &delz2) const;

 protected:

  class NeighborKokkos *neighborKK;

  typename ArrayTypes<DeviceType>::t_x_array_randomread x;
  typename ArrayTypes<DeviceType>::t_f_array f;
  typename ArrayTypes<DeviceType>::t_int_2d anglelist;

  DAT::tdual_efloat_1d k_eatom;
  DAT::tdual_virial_array k_vatom;
  DAT::t_efloat_1d d_eatom;
  DAT::t_virial_array d_vatom;

  int nlocal,newton_bond;
  int eflag,vflag;

  DAT::tdual_ffloat_1d k_k;
  DAT::tdual_ffloat_1d k_theta0;
  DAT::tdual_ffloat_1d k_k_ub;
  DAT::tdual_ffloat_1d k_r_ub;

  DAT::t_ffloat_1d d_k;
  DAT::t_ffloat_1d d_theta0;
  DAT::t_ffloat_1d d_k_ub;
  DAT::t_ffloat_1d d_r_ub;

  virtual void allocate();
};

}

#endif
#endif

/* ERROR/WARNING messages:

*/
