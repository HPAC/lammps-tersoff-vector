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

#include <stdio.h>


#include "pair_tersoff_cuda_cu.h"
__device__ __constant__ Param_Float params[MANYBODY_NPAIR* MANYBODY_NPAIR* MANYBODY_NPAIR];
__device__ __constant__ F_CFLOAT* _glob_zeta_ij; //zeta_ij
__device__ __constant__ F_CFLOAT4* _glob_r_ij; //r_ij (x,y,z,r^2) for pairs within force cutoff
__device__ __constant__ bool _zbl; //is tersoff zbl?


#include "pair_tersoff_cuda_kernel_nc.cu"

#include <time.h>


void Cuda_PairTersoffCuda_Init(cuda_shared_data* sdata, Param_Float* params_host, void* map_host, void* elem2param_host, int nelements_h, bool zbl)
{
  unsigned cuda_ntypes = sdata->atom.ntypes + 1;
  X_CFLOAT box_size[3] = {
    sdata->domain.subhi[0] - sdata->domain.sublo[0],
    sdata->domain.subhi[1] - sdata->domain.sublo[1],
    sdata->domain.subhi[2] - sdata->domain.sublo[2]
  };

  cudaMemcpyToSymbol(MY_AP(box_size)     , box_size                      , sizeof(X_CFLOAT) * 3);
  cudaMemcpyToSymbol(MY_AP(cuda_ntypes)  , &cuda_ntypes                   , sizeof(unsigned));
  cudaMemcpyToSymbol(MY_AP(virial)       , &sdata->pair.virial.dev_data   , sizeof(ENERGY_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(eng_vdwl)     , &sdata->pair.eng_vdwl.dev_data , sizeof(ENERGY_CFLOAT*));
  cudaMemcpyToSymbol(MY_AP(periodicity)  , sdata->domain.periodicity     , sizeof(int) * 3);
  cudaMemcpyToSymbol(MY_AP(collect_forces_later), &sdata->pair.collect_forces_later  , sizeof(int));
  cudaMemcpyToSymbol(params, params_host  , sizeof(Param_Float)*nelements_h * nelements_h * nelements_h);
  cudaMemcpyToSymbol(elem2param, elem2param_host  , sizeof(int)*nelements_h * nelements_h * nelements_h);
  cudaMemcpyToSymbol(map, map_host  , sizeof(int)*cuda_ntypes);
  cudaMemcpyToSymbol(nelements, &nelements_h, sizeof(int));
  cudaMemcpyToSymbol(_zbl, &zbl, sizeof(bool));

}

void Cuda_PairTersoffCuda(cuda_shared_data* sdata, cuda_shared_neighlist* sneighlist, int eflag, int vflag, int eflag_atom, int vflag_atom)
{
  static F_CFLOAT* glob_zeta_ij = NULL;
  static int glob_zeta_ij_size = 0;
  static F_CFLOAT4* glob_r_ij = NULL;
  static int* glob_numneigh_red = NULL;
  static int* glob_neighbors_red = NULL;
  static int* glob_neightype_red = NULL;

  if(glob_zeta_ij_size < sdata->atom.nall * sneighlist->maxneighbors * sizeof(F_CFLOAT)) {
    glob_zeta_ij_size = sdata->atom.nall * sneighlist->maxneighbors * sizeof(F_CFLOAT);
    cudaFree(glob_zeta_ij);
    cudaFree(glob_r_ij);
    cudaFree(glob_numneigh_red);
    cudaFree(glob_neighbors_red);
    cudaFree(glob_neightype_red);
    cudaMalloc(&glob_zeta_ij, glob_zeta_ij_size);
    cudaMalloc(&glob_r_ij, glob_zeta_ij_size * 4);
    cudaMalloc(&glob_numneigh_red, sdata->atom.nall * sizeof(int));
    cudaMalloc(&glob_neighbors_red, sdata->atom.nall * sneighlist->maxneighbors * sizeof(int));
    cudaMalloc(&glob_neightype_red, sdata->atom.nall * sneighlist->maxneighbors * sizeof(int));
    cudaMemcpyToSymbol(_glob_numneigh_red, &glob_numneigh_red  , sizeof(int*));
    cudaMemcpyToSymbol(_glob_neighbors_red, &glob_neighbors_red  , sizeof(int*));
    cudaMemcpyToSymbol(_glob_neightype_red, &glob_neightype_red  , sizeof(int*));
    cudaMemcpyToSymbol(_glob_r_ij, &glob_r_ij  , sizeof(F_CFLOAT4*));
    cudaMemcpyToSymbol(_glob_zeta_ij, &glob_zeta_ij  , sizeof(F_CFLOAT*));
  }

  dim3 grid, threads;
  int sharedperproc;

  Cuda_Pair_PreKernel_AllStyles(sdata, sneighlist, eflag, vflag, grid, threads, sharedperproc, false, 64);
  cudaStream_t* streams = (cudaStream_t*) CudaWrapper_returnStreams();



  dim3 grid2;

  if(sdata->atom.nall <= 256 * 64000) {
    grid2.x = (sdata->atom.nall + 255) / 256;
    grid2.y = 1;
  } else {
    grid2.x = (sdata->atom.nall + 256 * 128 - 1) / (256 * 128);
    grid2.y = 128;
  }

  grid2.z = 1;
  dim3 threads2;
  threads2.x = 256;
  threads2.y = 1;
  threads2.z = 1;

  my_times time1, time2;

  //pre-calculate all neighbordistances and zeta_ij
  my_gettime(CLOCK_REALTIME, &time1);
  Pair_Tersoff_Kernel_TpA_RIJ <<< grid2, threads2, 0, streams[1]>>>
  ();
  cudaThreadSynchronize();
  Pair_Tersoff_Kernel_TpA_ZetaIJ <<< grid2, threads2, 0, streams[1]>>>
  ();
  cudaThreadSynchronize();
  my_gettime(CLOCK_REALTIME, &time2);
  sdata->cuda_timings.test1 +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;
  my_gettime(CLOCK_REALTIME, &time1);

  //actual force calculation
  unsigned int sharedsize = (sharedperproc * sizeof(ENERGY_CFLOAT) + 4 * sizeof(F_CFLOAT)) * threads.x; //extra 4 floats per thread used to reduce register pressure

  if(eflag) {
    if(vflag)
      Pair_Tersoff_Kernel_TpA<1, 1> <<< grid, threads, sharedsize, streams[1]>>>
      (eflag_atom, vflag_atom);
    else
      Pair_Tersoff_Kernel_TpA<1, 0> <<< grid, threads, sharedsize, streams[1]>>>
      (eflag_atom, vflag_atom);
  } else {
    if(vflag)
      Pair_Tersoff_Kernel_TpA<0, 1> <<< grid, threads, sharedsize, streams[1]>>>
      (eflag_atom, vflag_atom);
    else
      Pair_Tersoff_Kernel_TpA<0, 0> <<< grid, threads, sharedsize, streams[1]>>>
      (eflag_atom, vflag_atom);
  }
  cudaThreadSynchronize();
  my_gettime(CLOCK_REALTIME, &time2);
  sdata->cuda_timings.test2 +=
    time2.tv_sec - time1.tv_sec + 1.0 * (time2.tv_nsec - time1.tv_nsec) / 1000000000;

  Cuda_Pair_PostKernel_AllStyles(sdata, grid, sharedperproc, eflag, vflag);
}

