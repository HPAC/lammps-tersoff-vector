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

#ifndef LMP_NEIGHBOR_KOKKOS_H
#define LMP_NEIGHBOR_KOKKOS_H

#include "neighbor.h"
#include "neigh_list_kokkos.h"
#include "neigh_bond_kokkos.h"
#include "kokkos_type.h"
#include <math.h>

namespace LAMMPS_NS {

template<class Device>
class NeighborKokkosExecute
{
  typedef ArrayTypes<Device> AT;

 public:
  NeighListKokkos<Device> neigh_list;
  const typename AT::t_xfloat_2d_randomread cutneighsq;
  const typename AT::t_int_1d bincount;
  const typename AT::t_int_1d_const c_bincount;
  typename AT::t_int_2d bins;
  typename AT::t_int_2d_const c_bins;
  const typename AT::t_x_array_randomread x;
  const typename AT::t_int_1d_const type,mask,molecule;

  const typename AT::t_tagint_1d_const tag;
  const typename AT::t_tagint_2d_const special;
  const typename AT::t_int_2d_const nspecial;
  const int molecular;
  int moltemplate;

  int special_flag[4];

  const int nbinx,nbiny,nbinz;
  const int mbinx,mbiny,mbinz;
  const int mbinxlo,mbinylo,mbinzlo;
  const X_FLOAT bininvx,bininvy,bininvz;
  X_FLOAT bboxhi[3],bboxlo[3];

  const int nlocal;

  const int exclude;

  const int nex_type;
  const int maxex_type;
  const typename AT::t_int_1d_const ex1_type,ex2_type;
  const typename AT::t_int_2d_const ex_type;

  const int nex_group;
  const int maxex_group;
  const typename AT::t_int_1d_const ex1_group,ex2_group;
  const typename AT::t_int_1d_const ex1_bit,ex2_bit;

  const int nex_mol;
  const int maxex_mol;
  const typename AT::t_int_1d_const ex_mol_group;
  const typename AT::t_int_1d_const ex_mol_bit;

  typename AT::t_int_scalar resize;
  typename AT::t_int_scalar new_maxneighs;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_resize;
  typename ArrayTypes<LMPHostType>::t_int_scalar h_new_maxneighs;

  const int xperiodic, yperiodic, zperiodic;
  const int xprd_half, yprd_half, zprd_half;

  NeighborKokkosExecute(
                        const NeighListKokkos<Device> &_neigh_list,
                        const typename AT::t_xfloat_2d_randomread &_cutneighsq,
                        const typename AT::t_int_1d &_bincount,
                        const typename AT::t_int_2d &_bins,
                        const int _nlocal,
                        const typename AT::t_x_array_randomread &_x,
                        const typename AT::t_int_1d_const &_type,
                        const typename AT::t_int_1d_const &_mask,
                        const typename AT::t_int_1d_const &_molecule,
                        const typename AT::t_tagint_1d_const &_tag,
                        const typename AT::t_tagint_2d_const &_special,
                        const typename AT::t_int_2d_const &_nspecial,
                        const int &_molecular,
                        const int & _nbinx,const int & _nbiny,const int & _nbinz,
                        const int & _mbinx,const int & _mbiny,const int & _mbinz,
                        const int & _mbinxlo,const int & _mbinylo,const int & _mbinzlo,
                        const X_FLOAT &_bininvx,const X_FLOAT &_bininvy,const X_FLOAT &_bininvz,
                        const int & _exclude,const int & _nex_type,const int & _maxex_type,
                        const typename AT::t_int_1d_const & _ex1_type,
                        const typename AT::t_int_1d_const & _ex2_type,
                        const typename AT::t_int_2d_const & _ex_type,
                        const int & _nex_group,const int & _maxex_group,
                        const typename AT::t_int_1d_const & _ex1_group,
                        const typename AT::t_int_1d_const & _ex2_group,
                        const typename AT::t_int_1d_const & _ex1_bit,
                        const typename AT::t_int_1d_const & _ex2_bit,
                        const int & _nex_mol,const int & _maxex_mol,
                        const typename AT::t_int_1d_const & _ex_mol_group,
                        const typename AT::t_int_1d_const & _ex_mol_bit,
                        const X_FLOAT *_bboxhi, const X_FLOAT* _bboxlo,
                        const int & _xperiodic, const int & _yperiodic, const int & _zperiodic,
                        const int & _xprd_half, const int & _yprd_half, const int & _zprd_half):
    neigh_list(_neigh_list), cutneighsq(_cutneighsq),
    bincount(_bincount),c_bincount(_bincount),bins(_bins),c_bins(_bins),
    nlocal(_nlocal),
    x(_x),type(_type),mask(_mask),molecule(_molecule),
    tag(_tag),special(_special),nspecial(_nspecial),molecular(_molecular),
    nbinx(_nbinx),nbiny(_nbiny),nbinz(_nbinz),
    mbinx(_mbinx),mbiny(_mbiny),mbinz(_mbinz),
    mbinxlo(_mbinxlo),mbinylo(_mbinylo),mbinzlo(_mbinzlo),
    bininvx(_bininvx),bininvy(_bininvy),bininvz(_bininvz),
    exclude(_exclude),nex_type(_nex_type),maxex_type(_maxex_type),
    ex1_type(_ex1_type),ex2_type(_ex2_type),ex_type(_ex_type),
    nex_group(_nex_group),maxex_group(_maxex_group),
    ex1_group(_ex1_group),ex2_group(_ex2_group),
    ex1_bit(_ex1_bit),ex2_bit(_ex2_bit),nex_mol(_nex_mol),maxex_mol(_maxex_mol),
    ex_mol_group(_ex_mol_group),ex_mol_bit(_ex_mol_bit),
    xperiodic(_xperiodic),yperiodic(_yperiodic),zperiodic(_zperiodic),
    xprd_half(_xprd_half),yprd_half(_yprd_half),zprd_half(_zprd_half){

    if (molecular == 2) moltemplate = 1;
    else moltemplate = 0;

    bboxlo[0] = _bboxlo[0]; bboxlo[1] = _bboxlo[1]; bboxlo[2] = _bboxlo[2];
    bboxhi[0] = _bboxhi[0]; bboxhi[1] = _bboxhi[1]; bboxhi[2] = _bboxhi[2];

    resize = typename AT::t_int_scalar("NeighborKokkosFunctor::resize");
#ifndef KOKKOS_USE_CUDA_UVM
    h_resize = Kokkos::create_mirror_view(resize);
#else
    h_resize = resize;
#endif
    h_resize() = 1;
    new_maxneighs = typename AT::
      t_int_scalar("NeighborKokkosFunctor::new_maxneighs");
#ifndef KOKKOS_USE_CUDA_UVM
    h_new_maxneighs = Kokkos::create_mirror_view(new_maxneighs);
#else
    h_new_maxneighs = new_maxneighs;
#endif
    h_new_maxneighs() = neigh_list.maxneighs;
  };

  ~NeighborKokkosExecute() {neigh_list.clean_copy();};

  template<int HalfNeigh, int GhostNewton>
  KOKKOS_FUNCTION
  void build_Item(const int &i) const;

  KOKKOS_FUNCTION
  void build_Item_Full_Ghost(const int &i) const;

  template<int ClusterSize>
  KOKKOS_FUNCTION
  void build_cluster_Item(const int &i) const;

#ifdef KOKKOS_HAVE_CUDA
  template<int HalfNeigh, int GhostNewton>
  __device__ inline
  void build_ItemCuda(typename Kokkos::TeamPolicy<Device>::member_type dev) const;
#endif

  KOKKOS_INLINE_FUNCTION
  void binatomsItem(const int &i) const;

  KOKKOS_INLINE_FUNCTION
  int coord2bin(const X_FLOAT & x,const X_FLOAT & y,const X_FLOAT & z) const
  {
    int ix,iy,iz;

    if (x >= bboxhi[0])
      ix = static_cast<int> ((x-bboxhi[0])*bininvx) + nbinx;
    else if (x >= bboxlo[0]) {
      ix = static_cast<int> ((x-bboxlo[0])*bininvx);
      ix = MIN(ix,nbinx-1);
    } else
      ix = static_cast<int> ((x-bboxlo[0])*bininvx) - 1;

    if (y >= bboxhi[1])
      iy = static_cast<int> ((y-bboxhi[1])*bininvy) + nbiny;
    else if (y >= bboxlo[1]) {
      iy = static_cast<int> ((y-bboxlo[1])*bininvy);
      iy = MIN(iy,nbiny-1);
    } else
      iy = static_cast<int> ((y-bboxlo[1])*bininvy) - 1;

    if (z >= bboxhi[2])
      iz = static_cast<int> ((z-bboxhi[2])*bininvz) + nbinz;
    else if (z >= bboxlo[2]) {
      iz = static_cast<int> ((z-bboxlo[2])*bininvz);
      iz = MIN(iz,nbinz-1);
    } else
      iz = static_cast<int> ((z-bboxlo[2])*bininvz) - 1;

    return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
  }

  KOKKOS_INLINE_FUNCTION
  int coord2bin(const X_FLOAT & x,const X_FLOAT & y,const X_FLOAT & z, int* i) const
  {
    int ix,iy,iz;

    if (x >= bboxhi[0])
      ix = static_cast<int> ((x-bboxhi[0])*bininvx) + nbinx;
    else if (x >= bboxlo[0]) {
      ix = static_cast<int> ((x-bboxlo[0])*bininvx);
      ix = MIN(ix,nbinx-1);
    } else
      ix = static_cast<int> ((x-bboxlo[0])*bininvx) - 1;

    if (y >= bboxhi[1])
      iy = static_cast<int> ((y-bboxhi[1])*bininvy) + nbiny;
    else if (y >= bboxlo[1]) {
      iy = static_cast<int> ((y-bboxlo[1])*bininvy);
      iy = MIN(iy,nbiny-1);
    } else
      iy = static_cast<int> ((y-bboxlo[1])*bininvy) - 1;

    if (z >= bboxhi[2])
      iz = static_cast<int> ((z-bboxhi[2])*bininvz) + nbinz;
    else if (z >= bboxlo[2]) {
      iz = static_cast<int> ((z-bboxlo[2])*bininvz);
      iz = MIN(iz,nbinz-1);
    } else
      iz = static_cast<int> ((z-bboxlo[2])*bininvz) - 1;

    i[0] = ix - mbinxlo;
    i[1] = iy - mbinylo;
    i[2] = iz - mbinzlo;

    return (iz-mbinzlo)*mbiny*mbinx + (iy-mbinylo)*mbinx + (ix-mbinxlo);
  }

  KOKKOS_INLINE_FUNCTION
  int exclusion(const int &i,const int &j, const int &itype,const int &jtype) const;

  KOKKOS_INLINE_FUNCTION
  int find_special(const int &i, const int &j) const;

  KOKKOS_INLINE_FUNCTION
  int minimum_image_check(double dx, double dy, double dz) const {
    if (xperiodic && fabs(dx) > xprd_half) return 1;
    if (yperiodic && fabs(dy) > yprd_half) return 1;
    if (zperiodic && fabs(dz) > zprd_half) return 1;
    return 0;
  }

};

template<class Device>
struct NeighborKokkosBinAtomsFunctor {
  typedef Device device_type;

  const NeighborKokkosExecute<Device> c;

  NeighborKokkosBinAtomsFunctor(const NeighborKokkosExecute<Device> &_c):
    c(_c) {};
  ~NeighborKokkosBinAtomsFunctor() {}
  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.binatomsItem(i);
  }
};

template<class Device,int HALF_NEIGH,int GHOST_NEWTON>
struct NeighborKokkosBuildFunctor {
  typedef Device device_type;

  const NeighborKokkosExecute<Device> c;
  const size_t sharedsize;

  NeighborKokkosBuildFunctor(const NeighborKokkosExecute<Device> &_c,
                             const size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_Item<HALF_NEIGH,GHOST_NEWTON>(i);
  }
#ifdef KOKKOS_HAVE_CUDA
  KOKKOS_INLINE_FUNCTION
  void operator() (typename Kokkos::TeamPolicy<Device>::member_type dev) const {
    c.template build_ItemCuda<HALF_NEIGH,GHOST_NEWTON>(dev);
  }
  size_t shmem_size(const int team_size) const { (void) team_size; return sharedsize; }
#endif
};

template<class Device>
struct NeighborKokkosBuildFunctorFullGhost {
  typedef Device device_type;

  const NeighborKokkosExecute<Device> c;
  const size_t sharedsize;

  NeighborKokkosBuildFunctorFullGhost(const NeighborKokkosExecute<Device> &_c,
                             const size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.build_Item_Full_Ghost(i);
  }
};

template<class Device,int ClusterSize>
struct NeighborClusterKokkosBuildFunctor {
  typedef Device device_type;

  const NeighborKokkosExecute<Device> c;
  const size_t sharedsize;

  NeighborClusterKokkosBuildFunctor(const NeighborKokkosExecute<Device> &_c,
                             const size_t _sharedsize):c(_c),
                             sharedsize(_sharedsize) {};

  KOKKOS_INLINE_FUNCTION
  void operator() (const int & i) const {
    c.template build_cluster_Item<ClusterSize>(i);
  }
};

template<class DeviceType>
struct TagNeighborCheckDistance{};

template<class DeviceType>
struct TagNeighborXhold{};

class NeighborKokkos : public Neighbor {
 public:
  typedef int value_type;



  int nlist_host;                       // pairwise neighbor lists on Host
  NeighListKokkos<LMPHostType> **lists_host;
  int nlist_device;                     // pairwise neighbor lists on Device
  NeighListKokkos<LMPDeviceType> **lists_device;

  NeighBondKokkos<LMPHostType> neighbond_host;
  NeighBondKokkos<LMPDeviceType> neighbond_device;

  DAT::tdual_int_2d k_bondlist;
  DAT::tdual_int_2d k_anglelist;
  DAT::tdual_int_2d k_dihedrallist;
  DAT::tdual_int_2d k_improperlist;

  NeighborKokkos(class LAMMPS *);
  ~NeighborKokkos();
  void init();

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighborCheckDistance<DeviceType>, const int&, int&) const;

  template<class DeviceType>
  KOKKOS_INLINE_FUNCTION
  void operator()(TagNeighborXhold<DeviceType>, const int&) const;

 private:
  int atoms_per_bin;
  DAT::tdual_xfloat_2d k_cutneighsq;
  DAT::tdual_int_1d k_bincount;
  DAT::tdual_int_2d k_bins;

  DAT::tdual_int_1d k_ex1_type,k_ex2_type;
  DAT::tdual_int_2d k_ex_type;
  DAT::tdual_int_1d k_ex1_group,k_ex2_group;
  DAT::tdual_int_1d k_ex1_bit,k_ex2_bit;
  DAT::tdual_int_1d k_ex_mol_group;
  DAT::tdual_int_1d k_ex_mol_bit;

  DAT::tdual_x_array x;
  DAT::tdual_x_array xhold;

  X_FLOAT deltasq;
  int device_flag;

  void init_cutneighsq_kokkos(int);
  int init_lists_kokkos();
  void init_list_flags1_kokkos(int);
  void init_list_flags2_kokkos(int);
  void init_list_grow_kokkos(int);
  void init_ex_type_kokkos(int);
  void init_ex_bit_kokkos();
  void init_ex_mol_bit_kokkos();
  void choose_build(int, NeighRequest *);
  virtual int check_distance();
  template<class DeviceType> int check_distance_kokkos();
  virtual void build(int);
  template<class DeviceType> void build_kokkos(int);
  void setup_bins_kokkos(int);
  void modify_ex_type_grow_kokkos();
  void modify_ex_group_grow_kokkos();
  void modify_mol_group_grow_kokkos();
  void init_topology_kokkos();
  void build_topology_kokkos();

  typedef void (NeighborKokkos::*PairPtrHost)
    (class NeighListKokkos<LMPHostType> *);
  PairPtrHost *pair_build_host;
  typedef void (NeighborKokkos::*PairPtrDevice)
    (class NeighListKokkos<LMPDeviceType> *);
  PairPtrDevice *pair_build_device;

  template<class DeviceType,int HALF_NEIGH, int GHOST>
  void full_bin_kokkos(NeighListKokkos<DeviceType> *list);
  template<class DeviceType>
  void full_bin_cluster_kokkos(NeighListKokkos<DeviceType> *list);

  typedef void (NeighborKokkos::*StencilPtrHost)
    (class NeighListKokkos<LMPHostType> *, int, int, int);
  StencilPtrHost *stencil_create_host;
  typedef void (NeighborKokkos::*StencilPtrDevice)
    (class NeighListKokkos<LMPDeviceType> *, int, int, int);
  StencilPtrDevice *stencil_create_device;
};

}

#endif

/* ERROR/WARNING messages:

E: Cannot (yet) request ghost atoms with Kokkos half neighbor list

This feature is not yet supported.

E: Too many local+ghost atoms for neighbor list

The number of nlocal + nghost atoms on a processor
is limited by the size of a 32-bit integer with 2 bits
removed for masking 1-2, 1-3, 1-4 neighbors.

*/
