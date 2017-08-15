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
   Contributing author: W. Michael Brown (Intel)
------------------------------------------------------------------------- */

#include "pair_sw_intel.h"

#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(push,target(mic))
#endif
#include <cmath>
#ifdef _LMP_INTEL_OFFLOAD
#pragma offload_attribute(pop)
#endif

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "atom.h"
#include "neighbor.h"
#include "neigh_request.h"
#include "force.h"
#include "comm.h"
#include "memory.h"
#include "neighbor.h"
#include "neigh_list.h"
#include "memory.h"
#include "error.h"
#include "modify.h"
#include "suffix.h"

#ifdef LMP_USE_AVXCD
#define OUTER_CHUNK 1
#include "intel_simd.h"
using namespace ip_simd;
#endif

using namespace LAMMPS_NS;

#define FC_PACKED0_T typename ForceConst<flt_t>::fc_packed0
#define FC_PACKED1_T typename ForceConst<flt_t>::fc_packed1
#define FC_PACKED1p2_T typename ForceConst<flt_t>::fc_packed1p2
#define FC_PACKED2_T typename ForceConst<flt_t>::fc_packed2
#define FC_PACKED3_T typename ForceConst<flt_t>::fc_packed3

#define MAXLINE 1024
#define DELTA 4

/* ---------------------------------------------------------------------- */

PairSWIntel::PairSWIntel(LAMMPS *lmp) : PairSW(lmp)
{
  suffix_flag |= Suffix::INTEL;
}

/* ---------------------------------------------------------------------- */

PairSWIntel::~PairSWIntel()
{
}

/* ---------------------------------------------------------------------- */

void PairSWIntel::compute(int eflag, int vflag)
{
  if (fix->precision() == FixIntel::PREC_MODE_MIXED)
    compute<float,double>(eflag, vflag, fix->get_mixed_buffers(),
			  force_const_single);
  else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE)
    compute<double,double>(eflag, vflag, fix->get_double_buffers(),
                           force_const_double);
  else
    compute<float,float>(eflag, vflag, fix->get_single_buffers(),
                         force_const_single);

  fix->balance_stamp();
  vflag_fdotr = 0;
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void PairSWIntel::compute(int eflag, int vflag,
                          IntelBuffers<flt_t,acc_t> *buffers,
                          const ForceConst<flt_t> &fc)
{
  if (eflag || vflag) {
    ev_setup(eflag, vflag);
  } else evflag = vflag_fdotr = 0;

  const int inum = list->inum;
  const int nthreads = comm->nthreads;
  const int host_start = fix->host_start_pair();
  const int offload_end = fix->offload_end_pair();
  const int ago = neighbor->ago;

  if (ago != 0 && fix->separate_buffers() == 0) {
    fix->start_watch(TIME_PACK);

    IntelBuffers<flt_t,acc_t> *buffers_ = buffers;
    IntelBuffers<flt_t,acc_t> *buffers = buffers_;
    #if defined(_OPENMP)
    #pragma omp parallel default(none) shared(eflag,vflag,buffers,fc)
    #endif
    {
      int ifrom, ito, tid;
      IP_PRE_omp_range_id_align(ifrom, ito, tid, atom->nlocal + atom->nghost,
                                nthreads, sizeof(ATOM_T));
      buffers->thr_pack(ifrom, ito, ago);
    }

    fix->stop_watch(TIME_PACK);
  }

  if (_onetype) {
    if (_spq) {
      if (evflag || vflag_fdotr) {
	int ovflag = 0;
	if (vflag_fdotr) ovflag = 2;
	else if (vflag) ovflag = 1;
	if (eflag) {
	  eval<1,1,1,1>(1, ovflag, buffers, fc, 0, offload_end, _offload_pad);
	  eval<1,1,1,1>(0, ovflag, buffers, fc, host_start, inum, _host_pad);
	} else {
	  eval<1,1,1,0>(1, ovflag, buffers, fc, 0, offload_end, _offload_pad);
	  eval<1,1,1,0>(0, ovflag, buffers, fc, host_start, inum, _host_pad);
	}
      } else {
	eval<1,1,0,0>(1, 0, buffers, fc, 0, offload_end, _offload_pad);
	eval<1,1,0,0>(0, 0, buffers, fc, host_start, inum, _host_pad);
      }
    } else {
      if (evflag || vflag_fdotr) {
	int ovflag = 0;
	if (vflag_fdotr) ovflag = 2;
	else if (vflag) ovflag = 1;
	if (eflag) {
	  eval<0,1,1,1>(1, ovflag, buffers, fc, 0, offload_end, _offload_pad);
	  eval<0,1,1,1>(0, ovflag, buffers, fc, host_start, inum, _host_pad);
	} else {
	  eval<0,1,1,0>(1, ovflag, buffers, fc, 0, offload_end, _offload_pad);
	  eval<0,1,1,0>(0, ovflag, buffers, fc, host_start, inum, _host_pad);
	}
      } else {
	eval<0,1,0,0>(1, 0, buffers, fc, 0, offload_end, _offload_pad);
	eval<0,1,0,0>(0, 0, buffers, fc, host_start, inum, _host_pad);
      }
    }
  } else {
    if (_spq) {
      if (evflag || vflag_fdotr) {
	int ovflag = 0;
	if (vflag_fdotr) ovflag = 2;
	else if (vflag) ovflag = 1;
	if (eflag) {
	  eval<1,0,1,1>(1, ovflag, buffers, fc, 0, offload_end, _offload_pad);
	  eval<1,0,1,1>(0, ovflag, buffers, fc, host_start, inum, _host_pad);
	} else {
	  eval<1,0,1,0>(1, ovflag, buffers, fc, 0, offload_end, _offload_pad);
	  eval<1,0,1,0>(0, ovflag, buffers, fc, host_start, inum, _host_pad);
	}
      } else {
	eval<1,0,0,0>(1, 0, buffers, fc, 0, offload_end, _offload_pad);
	eval<1,0,0,0>(0, 0, buffers, fc, host_start, inum, _host_pad);
      }
    } else {
      if (evflag || vflag_fdotr) {
	int ovflag = 0;
	if (vflag_fdotr) ovflag = 2;
	else if (vflag) ovflag = 1;
	if (eflag) {
	  eval<0,0,1,1>(1, ovflag, buffers, fc, 0, offload_end, _offload_pad);
	  eval<0,0,1,1>(0, ovflag, buffers, fc, host_start, inum, _host_pad);
	} else {
	  eval<0,0,1,0>(1, ovflag, buffers, fc, 0, offload_end, _offload_pad);
	  eval<0,0,1,0>(0, ovflag, buffers, fc, host_start, inum, _host_pad);
	}
      } else {
	eval<0,0,0,0>(1, 0, buffers, fc, 0, offload_end, _offload_pad);
	eval<0,0,0,0>(0, 0, buffers, fc, host_start, inum, _host_pad);
      }
    }
  }
}

/* ---------------------------------------------------------------------- */
#ifndef LMP_USE_AVXCD

template <int SPQ,int ONETYPE,int EVFLAG,int EFLAG,class flt_t,class acc_t>
void PairSWIntel::eval(const int offload, const int vflag,
                       IntelBuffers<flt_t,acc_t> *buffers,
                       const ForceConst<flt_t> &fc, const int astart,
		       const int aend, const int pad_width)
{
  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  IP_PRE_pack_separate_buffers(fix, buffers, ago, offload, nlocal, nall);

  ATOM_T * _noalias const x = buffers->get_x(offload);
  const int * _noalias const numneighhalf = buffers->get_atombin();
  const int * _noalias const numneigh = list->numneigh;
  const int * _noalias const cnumneigh = buffers->cnumneigh(list);
  const int * _noalias const firstneigh = buffers->firstneigh(list);

  const FC_PACKED0_T * _noalias const p2 = fc.p2[0];
  const FC_PACKED1_T * _noalias const p2f = fc.p2f[0];
  const FC_PACKED1p2_T * _noalias const p2f2 = fc.p2f2[0];
  const FC_PACKED2_T * _noalias const p2e = fc.p2e[0];
  const FC_PACKED3_T * _noalias const p3 = fc.p3[0][0];

  flt_t * _noalias const ccachex = buffers->get_ccachex();
  flt_t * _noalias const ccachey = buffers->get_ccachey();
  flt_t * _noalias const ccachez = buffers->get_ccachez();
  flt_t * _noalias const ccachew = buffers->get_ccachew();
  int * _noalias const ccachei = buffers->get_ccachei();
  int * _noalias const ccachej = buffers->get_ccachej();
  const int ccache_stride = _ccache_stride;

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  // Determine how much data to transfer
  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, /* NEWTON_PAIR*/ 1, EVFLAG, EFLAG, vflag,
                       buffers, offload, fix, separate_flag,
                       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;

  #ifdef _LMP_INTEL_OFFLOAD
  double *timer_compute = fix->off_watch_pair();
  int *overflow = fix->get_off_overflow_flag();

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if(offload) \
    in(p2,p2f,p2f2,p2e,p3:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(numneighhalf:length(0) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(ccachex,ccachey,ccachez,ccachew:length(0) alloc_if(0) free_if(0)) \
    in(ccachei,ccachej:length(0) alloc_if(0) free_if(0)) \
    in(ccache_stride,nthreads,inum,nall,ntypes,vflag,eatom,offload) \
    in(astart,nlocal,f_stride,minlocal,separate_flag,pad_width) \
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(f_start)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(1, separate_flag, nlocal, nall,
                              f_stride, x, 0);

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EVFLAG) {
      oevdwl = (acc_t)0;
      if (vflag) ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;
    }

    #if defined(_OPENMP)
    #pragma omp parallel default(none) \
      shared(f_start,f_stride,nlocal,nall,minlocal) \
      reduction(+:oevdwl,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iito, tid;
      IP_PRE_omp_range_id(iifrom, iito, tid, inum, nthreads);
      iifrom += astart;
      iito += astart;

      FORCE_T * _noalias const f = f_start - minlocal + (tid * f_stride);
      memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      const int toffs = tid * ccache_stride;
      flt_t * _noalias const tdelx = ccachex + toffs;
      flt_t * _noalias const tdely = ccachey + toffs;
      flt_t * _noalias const tdelz = ccachez + toffs;
      flt_t * _noalias const trsq = ccachew + toffs;
      int * _noalias const tj = ccachei + toffs;
      int * _noalias const tjtype = ccachej + toffs;

      // loop over full neighbor list of my atoms
      flt_t cutsq, cut, powerp, powerq, sigma, c1, c2, c3, c4, c5, c6;
      flt_t sigma_gamma, costheta, lambda_epsilon, lambda_epsilon2;
      if (ONETYPE) {
        cutsq = p2[3].cutsq;
	cut = p2f[3].cut;
	sigma = p2f[3].sigma;
	c1 = p2f2[3].c1;
	c2 = p2f2[3].c2;
	c3 = p2f2[3].c3;
	c4 = p2f2[3].c4;
	sigma_gamma = p2[3].sigma_gamma;
	costheta = p3[7].costheta;
	lambda_epsilon = p3[7].lambda_epsilon;
	lambda_epsilon2 = p3[7].lambda_epsilon2;
	if (SPQ == 0) {
          powerp = p2f[3].powerp;
	  powerq = p2f[3].powerq;
        }
	if (EFLAG) {
          c5 = p2e[3].c5;
	  c6 = p2e[3].c6;
        }
      }

      for (int i = iifrom; i < iito; ++i) {
        int itype, itype_offset;
        const flt_t xtmp = x[i].x;
        const flt_t ytmp = x[i].y;
        const flt_t ztmp = x[i].z;

	if (!ONETYPE) {
          itype = x[i].w;
	  itype_offset = itype * ntypes;
        } 

        const int * _noalias const jlist = firstneigh + cnumneigh[i];
        const int jnum = numneigh[i];
	const int jnumhalf = numneighhalf[i];

        acc_t fxtmp, fytmp, fztmp, fwtmp;
        acc_t sevdwl, sv0, sv1, sv2, sv3, sv4, sv5;
        fxtmp = fytmp = fztmp = (acc_t)0.0;
        if (EVFLAG) {
          if (EFLAG) fwtmp = sevdwl = (acc_t)0;
          if (vflag==1) sv0 = sv1 = sv2 = sv3 = sv4 = sv5 = (acc_t)0;
        }

	int ejnum = 0, ejnumhalf = 0;
        for (int jj = 0; jj < jnum; jj++) {
          int j = jlist[jj];
          j &= NEIGHMASK;
          const flt_t delx = x[j].x - xtmp;
          const flt_t dely = x[j].y - ytmp;
          const flt_t delz = x[j].z - ztmp;
          int jtype, ijtype;
          if (!ONETYPE) {
	    jtype = x[j].w;
	    ijtype = itype_offset + jtype;
	    cutsq = p2[ijtype].cutsq;
	  } 
          const flt_t rsq1 = delx * delx + dely * dely + delz * delz;
          if (rsq1 < cutsq) {
	    tdelx[ejnum] = delx;
	    tdely[ejnum] = dely;
	    tdelz[ejnum] = delz;
	    trsq[ejnum] = rsq1;
	    tj[ejnum] = j;
	    if (!ONETYPE) tjtype[ejnum] = jtype;
	    ejnum++;
	    if (jj < jnumhalf) ejnumhalf++;
	  }
	}
	int ejnum_pad = ejnum;
	
	while ( (ejnum_pad % pad_width) != 0) {
	  tdelx[ejnum_pad] = (flt_t)0.0;
	  tdely[ejnum_pad] = (flt_t)0.0;
	  tdelz[ejnum_pad] = (flt_t)0.0;
	  trsq[ejnum_pad] = p2[3].cutsq + (flt_t)1.0;
	  tj[ejnum_pad] = nall;
	  if (!ONETYPE) tjtype[ejnum_pad] = 0;
	  ejnum_pad++;
	}
	
        #if defined(LMP_SIMD_COMPILER)
	#pragma vector aligned
		#pragma simd reduction(+:fxtmp, fytmp, fztmp, fwtmp, sevdwl, \
					 sv0, sv1, sv2, sv3, sv4, sv5)
	#endif
        for (int jj = 0; jj < ejnum_pad; jj++) {
          acc_t fjxtmp, fjytmp, fjztmp, fjtmp;
          fjxtmp = fjytmp = fjztmp = (acc_t)0.0;
          if (EFLAG) fjtmp = (acc_t)0.0;
	  int ijtype;

          const flt_t delx = tdelx[jj];
          const flt_t dely = tdely[jj];
          const flt_t delz = tdelz[jj];
	  if (!ONETYPE) ijtype = tjtype[jj] + itype_offset;
          const flt_t rsq1 = trsq[jj];

          const flt_t r1 = sqrt(rsq1);
          const flt_t rinvsq1 = (flt_t)1.0 / rsq1;
	  if (!ONETYPE) cut = p2f[ijtype].cut;
          const flt_t rainv1 = (flt_t)1.0 / (r1 - cut);
	  
	  // two-body interactions, skip half of them
	  flt_t rp, rq;
	  if (SPQ == 1) {
	    rp = r1 * r1;
	    rp *= rp;
	    rp = (flt_t)1.0 / rp;
	    rq = (flt_t)1.0;
	  } else {
            if (!ONETYPE) {
              powerp = p2f[ijtype].powerp;
	      powerq = p2f[ijtype].powerq;
            }
	    rp = std::pow(r1, powerp);
	    rq = std::pow(r1, powerq);
	  }

	  if (!ONETYPE) {
            sigma = p2f[ijtype].sigma;
	    c1 = p2f2[ijtype].c1;
	    c2 = p2f2[ijtype].c2;
	    c3 = p2f2[ijtype].c3; 
	    c4 = p2f2[ijtype].c4;
          }

	  const flt_t rainvsq = rainv1 * rainv1 * r1;
	  flt_t expsrainv = exp(sigma * rainv1);
	  if (jj >= ejnumhalf) expsrainv = (flt_t)0.0;
	  const flt_t fpair = (c1 * rp - c2 * rq + (c3 * rp - c4 * rq) * 
			       rainvsq) * expsrainv * rinvsq1;

	  fxtmp -= delx * fpair;
	  fytmp -= dely * fpair;
	  fztmp -= delz * fpair;
	  fjxtmp += delx * fpair;
	  fjytmp += dely * fpair;
	  fjztmp += delz * fpair;

	  if (EVFLAG) {
	    if (EFLAG) {
	      flt_t evdwl;
	      if (!ONETYPE) {
		c5 = p2e[ijtype].c5;
		c6 = p2e[ijtype].c6;
	      }
	      evdwl = (c5 * rp - c6 * rq) * expsrainv;
	      sevdwl += evdwl;
	      if (eatom) {
		fwtmp += (acc_t)0.5 * evdwl;
		fjtmp += (acc_t)0.5 * evdwl;
	      }
	    }
	    IP_PRE_ev_tally_nbor(vflag, (flt_t)1.0, fpair,
				 -delx, -dely, -delz);
	  }

	  /*---------------------------------------------*/

	  int ijkoff;
	  if (!ONETYPE) {
            sigma_gamma = p2[ijtype].sigma_gamma;
	    ijkoff = ijtype * ntypes;
          }

          flt_t gsrainv1 = sigma_gamma * rainv1;
          flt_t gsrainvsq1 = gsrainv1 * rainv1 / r1;
          flt_t expgsrainv1 = exp(gsrainv1);

          for (int kk = 0; kk < ejnum; kk++) {
            int iktype, ijktype;
	    if (!ONETYPE) {
              iktype = tjtype[kk];
	      ijktype = ijkoff + iktype;
	      iktype += itype_offset;
	      cut = p2[iktype].cut;
	      sigma_gamma = p2[iktype].sigma_gamma;
	      costheta = p3[ijktype].costheta;
	      lambda_epsilon = p3[ijktype].lambda_epsilon;
	      lambda_epsilon2 = p3[ijktype].lambda_epsilon2;
            }

            flt_t delr2[3];
            delr2[0] = tdelx[kk];
            delr2[1] = tdely[kk];
            delr2[2] = tdelz[kk];
            const flt_t rsq2 = trsq[kk];

	    const flt_t r2 = sqrt(rsq2);
	    const flt_t rinvsq2 = (flt_t)1.0 / rsq2;
	    const flt_t rainv2 = (flt_t)1.0 / (r2 - cut);
	    const flt_t gsrainv2 = sigma_gamma * rainv2;
	    const flt_t gsrainvsq2 = gsrainv2 * rainv2 / r2;
	    const flt_t expgsrainv2 = exp(gsrainv2);

	    const flt_t rinv12 = (flt_t)1.0 / (r1 * r2);
	    const flt_t cs = (delx * delr2[0] + dely * delr2[1] +
                              delz * delr2[2]) * rinv12;
	    const flt_t delcs = cs - costheta;
	    const flt_t delcssq = delcs*delcs;

	    flt_t kfactor;
	    if (jj == kk || jj >= ejnum) kfactor = (flt_t)0.0;
	    else kfactor = (flt_t)1.0;

	    const flt_t facexp = expgsrainv1*expgsrainv2*kfactor;
	    const flt_t facrad = lambda_epsilon * facexp * delcssq;
	    const flt_t frad1 = facrad*gsrainvsq1;
	    const flt_t frad2 = facrad*gsrainvsq2;
	    const flt_t facang = lambda_epsilon2 * facexp * delcs;
	    const flt_t facang12 = rinv12*facang;
	    const flt_t csfacang = cs*facang;
	    const flt_t csfac1 = rinvsq1*csfacang;

	    const flt_t fjx = delx*(frad1+csfac1)-delr2[0]*facang12;
	    const flt_t fjy = dely*(frad1+csfac1)-delr2[1]*facang12;
	    const flt_t fjz = delz*(frad1+csfac1)-delr2[2]*facang12;

	    fxtmp -= fjx;
	    fytmp -= fjy;
	    fztmp -= fjz;
	    fjxtmp += fjx;
	    fjytmp += fjy;
	    fjztmp += fjz;

	    if (EVFLAG) {
	      if (EFLAG) {
	        const flt_t evdwl = facrad * (flt_t)0.5;
		sevdwl += evdwl;
		if (eatom) {
		  fwtmp += (acc_t)0.33333333 * evdwl;
		  fjtmp += (acc_t)0.33333333 * facrad;
		}
	      }
	      IP_PRE_ev_tally_nbor3v(vflag, fjx, fjy, fjz,
				     delx, dely, delz);
	    }
	  } // for kk
	  const int j = tj[jj];
          f[j].x += fjxtmp;
          f[j].y += fjytmp;
          f[j].z += fjztmp;
          if (EFLAG) 
	    if (eatom) f[j].w += fjtmp;
        } // for jj

        f[i].x += fxtmp;
        f[i].y += fytmp;
        f[i].z += fztmp;
        IP_PRE_ev_tally_atom(EVFLAG, EFLAG, vflag, f, fwtmp);
      } // for ii

      #ifndef _LMP_INTEL_OFFLOAD
      if (vflag == 2)
      #endif
      {
        #if defined(_OPENMP)
        #pragma omp barrier
        #endif
        IP_PRE_fdotr_acc_force(1, EVFLAG,  EFLAG, vflag, eatom, nall,
			       nlocal, minlocal, nthreads, f_start, f_stride,
			       x, offload);
      }
    } // end omp
    if (EVFLAG) {
      if (EFLAG) {
        ev_global[0] = oevdwl;
        ev_global[1] = (acc_t)0.0;
      }
      if (vflag) {
        ev_global[2] = ov0;
        ev_global[3] = ov1;
        ev_global[4] = ov2;
        ev_global[5] = ov3;
        ev_global[6] = ov4;
        ev_global[7] = ov5;
      }
    }
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end offload
  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EVFLAG)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, vflag);
  else
    fix->add_result_array(f_start, 0, offload);
}

#else 

/* ----------------------------------------------------------------------

Vector intrinsics are temporarily being used for the Stillinger-Weber
potential to allow for advanced features in the AVX512 instruction set to
be exploited on early hardware. We hope to see compiler improvements for
AVX512 that will eliminate this requirement, so it is not recommended to
develop code based on the intrinsics implementation. Please e-mail the
authors for more details.

------------------------------------------------------------------------- */

template <int SPQ,int ONETYPE,int EVFLAG,int EFLAG,class flt_t,class acc_t>
void PairSWIntel::eval(const int offload, const int vflag,
                       IntelBuffers<flt_t,acc_t> *buffers,
                       const ForceConst<flt_t> &fc, const int astart, 
		       const int aend, const int pad_width)
{
  typedef typename SIMD_type<flt_t>::SIMD_vec SIMD_flt_t;
  typedef typename SIMD_type<acc_t>::SIMD_vec SIMD_acc_t;
  const int swidth = SIMD_type<flt_t>::width();

  const int inum = aend - astart;
  if (inum == 0) return;
  int nlocal, nall, minlocal;
  fix->get_buffern(offload, nlocal, nall, minlocal);

  const int ago = neighbor->ago;
  IP_PRE_pack_separate_buffers(fix, buffers, ago, offload, nlocal, nall);

  ATOM_T * _noalias const x = buffers->get_x(offload);
  const int * _noalias const numneighhalf = buffers->get_atombin();
  const int * _noalias const numneigh = list->numneigh;
  const int * _noalias const cnumneigh = buffers->cnumneigh(list);
  const int * _noalias const firstneigh = buffers->firstneigh(list);

  const FC_PACKED0_T * _noalias const p2 = fc.p2[0];
  const FC_PACKED1_T * _noalias const p2f = fc.p2f[0];
  const FC_PACKED1p2_T * _noalias const p2f2 = fc.p2f2[0];
  const FC_PACKED2_T * _noalias const p2e = fc.p2e[0];
  const FC_PACKED3_T * _noalias const p3 = fc.p3[0][0];

  flt_t * _noalias const ccachex = buffers->get_ccachex();
  flt_t * _noalias const ccachey = buffers->get_ccachey();
  flt_t * _noalias const ccachez = buffers->get_ccachez();
  flt_t * _noalias const ccachew = buffers->get_ccachew();
  int * _noalias const ccachei = buffers->get_ccachei();
  int * _noalias const ccachej = buffers->get_ccachej();
  acc_t * _noalias const ccachef = buffers->get_ccachef();
  const int ccache_stride = _ccache_stride;
  const int ccache_stride3 = _ccache_stride3;

  const int ntypes = atom->ntypes + 1;
  const int eatom = this->eflag_atom;

  // Determine how much data to transfer
  int x_size, q_size, f_stride, ev_size, separate_flag;
  IP_PRE_get_transfern(ago, /* NEWTON_PAIR*/ 1, EVFLAG, EFLAG, vflag,
                       buffers, offload, fix, separate_flag,
                       x_size, q_size, ev_size, f_stride);

  int tc;
  FORCE_T * _noalias f_start;
  acc_t * _noalias ev_global;
  IP_PRE_get_buffers(offload, buffers, fix, tc, f_start, ev_global);
  const int nthreads = tc;

  #ifdef _LMP_INTEL_OFFLOAD
  double *timer_compute = fix->off_watch_pair();
  int *overflow = fix->get_off_overflow_flag();

  if (offload) fix->start_watch(TIME_OFFLOAD_LATENCY);
  #pragma offload target(mic:_cop) if(offload) \
    in(p2,p2f,p2f2,p2e,p3:length(0) alloc_if(0) free_if(0)) \
    in(firstneigh:length(0) alloc_if(0) free_if(0)) \
    in(cnumneigh:length(0) alloc_if(0) free_if(0)) \
    in(numneigh:length(0) alloc_if(0) free_if(0)) \
    in(x:length(x_size) alloc_if(0) free_if(0)) \
    in(numneighhalf:length(0) alloc_if(0) free_if(0)) \
    in(overflow:length(0) alloc_if(0) free_if(0)) \
    in(ccachex,ccachey,ccachez,ccachew:length(0) alloc_if(0) free_if(0)) \
    in(ccachei,ccachej,ccachef:length(0) alloc_if(0) free_if(0)) \
    in(ccache_stride,nthreads,inum,nall,ntypes,vflag,eatom,offload) \
    in(astart,nlocal,f_stride,minlocal,separate_flag,pad_width) \
    in(ccache_stride3)						\
    out(f_start:length(f_stride) alloc_if(0) free_if(0)) \
    out(ev_global:length(ev_size) alloc_if(0) free_if(0)) \
    out(timer_compute:length(1) alloc_if(0) free_if(0)) \
    signal(f_start)
  #endif
  {
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime();
    #endif

    IP_PRE_repack_for_offload(1, separate_flag, nlocal, nall,
                              f_stride, x, 0);

    acc_t oevdwl, ov0, ov1, ov2, ov3, ov4, ov5;
    if (EVFLAG) {
      oevdwl = (acc_t)0;
      if (vflag) ov0 = ov1 = ov2 = ov3 = ov4 = ov5 = (acc_t)0;
    }

    #if defined(_OPENMP)
    #pragma omp parallel default(none) \
      shared(f_start,f_stride,nlocal,nall,minlocal) \
      reduction(+:oevdwl,ov0,ov1,ov2,ov3,ov4,ov5)
    #endif
    {
      int iifrom, iito, tid;
      IP_PRE_omp_range_id_vec(iifrom, iito, tid, inum, nthreads, swidth);
      iifrom += astart;
      iito += astart;

      FORCE_T * _noalias const f = f_start - minlocal + (tid * f_stride);
      memset(f + minlocal, 0, f_stride * sizeof(FORCE_T));

      const int toffs = tid * ccache_stride;
      flt_t * _noalias const tdelx = ccachex + toffs;
      flt_t * _noalias const tdely = ccachey + toffs;
      flt_t * _noalias const tdelz = ccachez + toffs;
      flt_t * _noalias const trsq = ccachew + toffs;
      int * _noalias const tj = ccachei + toffs;
      int * _noalias const tjtype = ccachej + toffs;
      acc_t * _noalias const tf = ccachef + tid * ccache_stride3;

      // loop over full neighbor list of my atoms

      SIMD_flt_t cutsq, cut, powerp, powerq, sigma, c1, c2, c3,c4, c5, c6;
      SIMD_flt_t sigma_gamma, costheta, lambda_epsilon, lambda_epsilon2;
      if (ONETYPE) {
	cutsq = SIMD_set(p2[3].cutsq);
	cut = SIMD_set(p2f[3].cut);
	sigma = SIMD_set(p2f[3].sigma);
	c1 = SIMD_set(p2f2[3].c1);
	c2 = SIMD_set(p2f2[3].c2);
	c3 = SIMD_set(p2f2[3].c3);
	c4 = SIMD_set(p2f2[3].c4);
	sigma_gamma = SIMD_set(p2[3].sigma_gamma);
	costheta = SIMD_set(p3[7].costheta);
	lambda_epsilon = SIMD_set(p3[7].lambda_epsilon);
	lambda_epsilon2 = SIMD_set(p3[7].lambda_epsilon2);
	if (SPQ == 0) {
	  powerp = SIMD_set(p2f[3].powerp);
	  powerq = SIMD_set(p2f[3].powerq);
	}
	if (EFLAG) {
          c5 = SIMD_set(p2e[3].c5);
          c6 = SIMD_set(p2e[3].c6);
        }
      }

      SIMD_int ilist = SIMD_set(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15);
      const SIMD_int goffset = SIMD_set(0,16,32,48,64,80,96,112,128,
					144,160,176,192,208,224,240);
      ilist = ilist + iifrom;
      for (int i = iifrom; i < iito; i += swidth) {
	SIMD_mask imask = ilist < iito;
	SIMD_flt_t xtmp, ytmp, ztmp;
	SIMD_int itype, itype_offset;

	if (ONETYPE)
	  SIMD_atom_gather(imask, &(x[i].x), goffset, xtmp, ytmp, ztmp);
	else {
	  SIMD_atom_gather(imask, &(x[i].x), goffset, xtmp, ytmp, ztmp, itype);
	  itype_offset = itype * ntypes;
	}

	#ifdef OUTER_CHUNK
	const int* ng = firstneigh + cnumneigh[i] - swidth;
	#else
        SIMD_int ng = SIMD_load(cnumneigh + i);
	ng = ng - 1;
	#endif
	const SIMD_int jnum = SIMD_loadz(imask, numneigh + i);
	const SIMD_int jnumhalf = SIMD_loadz(imask, numneighhalf + i);
	const int jnum_max = SIMD_max(jnum);

	SIMD_acc_t fxtmp = SIMD_set((acc_t)0);
	SIMD_acc_t fytmp = SIMD_set((acc_t)0);
	SIMD_acc_t fztmp = SIMD_set((acc_t)0);
	SIMD_acc_t fwtmp, fxtmp2, fytmp2, fztmp2, fwtmp2;
	if (is_same<flt_t,acc_t>::value == 0) {
	  fxtmp2 = SIMD_set((acc_t)0);
	  fytmp2 = SIMD_set((acc_t)0);
	  fztmp2 = SIMD_set((acc_t)0);
          if (EFLAG) fwtmp2 = SIMD_set((acc_t)0);
	}

        SIMD_acc_t sevdwl, sv0, sv1, sv2, sv3, sv4, sv5;
        if (EVFLAG) {
          if (EFLAG) {
            fwtmp = SIMD_set((acc_t)0);
	    sevdwl = SIMD_set((acc_t)0);
          }
          if (vflag==1) {
            sv0 = SIMD_set((acc_t)0);
	    sv1 = SIMD_set((acc_t)0);
	    sv2 = SIMD_set((acc_t)0);
	    sv3 = SIMD_set((acc_t)0);
	    sv4 = SIMD_set((acc_t)0);
	    sv5 = SIMD_set((acc_t)0);
          }
        }

	SIMD_int ejnum = SIMD_set(0);
	SIMD_int ejnumhalf = SIMD_set(0);
	SIMD_int coffset = SIMD_set(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 
				    11, 12, 13, 14, 15);
        for (int jj = 0; jj < jnum_max; jj++) {
          SIMD_mask jmask = jj < jnum;

	  #ifdef OUTER_CHUNK
	  ng += swidth;
	  SIMD_int j = SIMD_load(ng);
	  #else
	  ng = ng + 1;
	  SIMD_int j = SIMD_gather(jmask, firstneigh, ng);
	  #endif
          j = j & SIMD_set(NEIGHMASK);
	  const SIMD_int joffset = j << 4;

	  SIMD_flt_t delx, dely, delz;
	  SIMD_int jtype, ijtype;
	  if (ONETYPE)
	    SIMD_atom_gather(jmask, &(x[0].x), joffset, delx, dely, delz);
	  else {
	    SIMD_atom_gather(jmask, &(x[0].x), joffset, delx, dely, delz, 
			     jtype);
	    ijtype = (jtype + itype_offset) << 2;
	    cutsq = SIMD_gather(jmask, &(p2[0].cutsq), ijtype);
	  }

	  delx = delx - xtmp;
	  dely = dely - ytmp;
	  delz = delz - ztmp;
          SIMD_flt_t rsq1 = delx * delx;
	  rsq1 = SIMD_fma(dely, dely, rsq1);
	  rsq1 = SIMD_fma(delz, delz, rsq1);

	  const SIMD_mask rmask = SIMD_lt(jmask, rsq1, cutsq);
	  SIMD_scatter(rmask, tdelx, coffset, delx);
	  SIMD_scatter(rmask, tdely, coffset, dely);
	  SIMD_scatter(rmask, tdelz, coffset, delz);
	  SIMD_scatter(rmask, trsq, coffset, rsq1);
	  SIMD_scatter(rmask, tj, coffset, j);
	  if (!ONETYPE) SIMD_scatter(rmask, tjtype, coffset, jtype);
	  ejnum = SIMD_add(rmask, ejnum, 1);
	  coffset = SIMD_add(rmask, coffset, swidth);
	  const SIMD_mask hmask = SIMD_lt(rmask, SIMD_set(jj), jnumhalf);
	  ejnumhalf = SIMD_add(hmask, ejnumhalf, 1);
	}

	const int ejnum_max = SIMD_max(ejnum);
	const int ejnumhalf_max = SIMD_max(ejnumhalf);
	memset(tf, 0, ejnum_max * sizeof(acc_t) * swidth * 3);
        for (int jj = 0; jj < ejnum_max; jj++) {
          SIMD_int ijtype;
	  const int coffset = jj * swidth;
	  if (!ONETYPE) {
	    ijtype = SIMD_load(tjtype + coffset);
	    ijtype = (ijtype + itype_offset) << 2;
	    cut = SIMD_gather(&(p2f[0].cut), ijtype);
	  }

	  SIMD_acc_t fjxtmp = SIMD_set((acc_t)0);
	  SIMD_acc_t fjytmp = SIMD_set((acc_t)0);
	  SIMD_acc_t fjztmp = SIMD_set((acc_t)0);
	  SIMD_acc_t fjtmp, fjxtmp2, fjytmp2, fjztmp2, fjtmp2;
          if (EFLAG) fjtmp = SIMD_set((acc_t)0.0);

	  if (is_same<flt_t,acc_t>::value == 0) {
	    fjxtmp2 = SIMD_set((acc_t)0);
	    fjytmp2 = SIMD_set((acc_t)0);
	    fjztmp2 = SIMD_set((acc_t)0);
	    if (EFLAG) fjtmp2 = SIMD_set((acc_t)0.0);
	  }

          const SIMD_flt_t delx = SIMD_load(tdelx + coffset);
          const SIMD_flt_t dely = SIMD_load(tdely + coffset);
          const SIMD_flt_t delz = SIMD_load(tdelz + coffset);
          const SIMD_flt_t rsq1 = SIMD_load(trsq + coffset);

          const SIMD_flt_t r1 = SIMD_sqrt(rsq1); 
          const SIMD_flt_t rinvsq1 = SIMD_rcp(rsq1);
          const SIMD_flt_t rainv1 = SIMD_rcp(r1 - cut);
	  
	  // two-body interactions, skip half of them
	  if (jj < ejnumhalf_max) {
            SIMD_flt_t rp, rq;
	    if (SPQ == 1) {
              rp = r1 * r1;
	      rp = rp * rp;
	      rp = SIMD_rcp(rp);
	      rq = SIMD_set((flt_t)1.0);
            } else {
	      if (!ONETYPE) {
		powerp = SIMD_gather(&(p2f[0].powerp), ijtype);
		powerq = SIMD_gather(&(p2f[0].powerq), ijtype);
	      }
	      rp = SIMD_pow(r1, powerp);
	      rq = SIMD_pow(r1, powerq);
	    }

	    if (!ONETYPE) {
	      sigma = SIMD_gather(&(p2f[0].sigma), ijtype);
	      c1 = SIMD_gather(&(p2f2[0].c1), ijtype);
	      c2 = SIMD_gather(&(p2f2[0].c2), ijtype);
	      c3 = SIMD_gather(&(p2f2[0].c3), ijtype);
	      c4 = SIMD_gather(&(p2f2[0].c4), ijtype);
	    }

	    const SIMD_flt_t rainvsq = rainv1 * rainv1 * r1;
	    const SIMD_flt_t expsrainv = SIMD_exp(sigma * rainv1);
	    const SIMD_flt_t fpair = (c1 * rp - c2 * rq + (c3 * rp - c4 * rq) *
				      rainvsq) * expsrainv * rinvsq1;

	    const SIMD_flt_t fjx = delx * fpair;
	    const SIMD_flt_t fjy = dely * fpair;
	    const SIMD_flt_t fjz = delz * fpair;

	    const SIMD_mask hmask = jj < ejnumhalf;
	    SIMD_accumulate3(hmask, fjx, fjy, fjz, fxtmp, fytmp, fztmp,
			     fjxtmp, fjytmp, fjztmp, fxtmp2, fytmp2,
			     fztmp2, fjxtmp2, fjytmp2, fjztmp2); 
          
	    if (EVFLAG) {
	      if (EFLAG) {
	        if (!ONETYPE) {
	          c5 = SIMD_gather(&(p2e[0].c5), ijtype);
	          c6 = SIMD_gather(&(p2e[0].c6), ijtype);
                }            
	        SIMD_flt_t evdwl;
		evdwl = (c5 * rp - c6 * rq) * expsrainv;
		SIMD_acc_energy3(hmask, evdwl, eatom, sevdwl, fwtmp, fjtmp,
                                 fwtmp2, fjtmp2);
	      }
	      SIMD_ev_tally_nbor(hmask, vflag, (flt_t)1.0, fpair, delx, dely, 
                                 delz, sv0, sv1, sv2, sv3, sv4, sv5);
	    }
          }

	  /*---------------------------------------------*/
	  SIMD_int ijkoff;
	  if (!ONETYPE) {
	    sigma_gamma = SIMD_gather(&(p2[0].sigma_gamma), ijtype);
	    ijkoff = ijtype * ntypes;
	  }
          const SIMD_flt_t gsrainv1 = sigma_gamma * rainv1;
          const SIMD_flt_t gsrainvsq1 = gsrainv1 * rainv1 / r1;
          const SIMD_flt_t expgsrainv1 = SIMD_exp(gsrainv1);

	  const SIMD_mask jmask = jj < ejnum;
          for (int kk = jj+1; kk < ejnum_max; kk++) {
	    SIMD_int iktype, ijktype;
	    const int kcoffset = kk * swidth;
	    if (!ONETYPE) {
	      iktype = SIMD_load(tjtype + kcoffset);
	      ijktype = ijkoff + (iktype << 2);
	      iktype = (iktype + itype_offset) << 2;
	      cut = SIMD_gather(&(p2[0].cut), iktype);
	      sigma_gamma = SIMD_gather(&(p2[0].sigma_gamma), iktype);
	      costheta = SIMD_gather(&(p3[0].costheta), ijktype);
	      lambda_epsilon = SIMD_gather(&(p3[0].lambda_epsilon), ijktype);
	      lambda_epsilon2 = SIMD_gather(&(p3[0].lambda_epsilon2), ijktype);
	    }
	    const SIMD_flt_t delr2x = SIMD_load(tdelx + kcoffset);
	    const SIMD_flt_t delr2y = SIMD_load(tdely + kcoffset);
	    const SIMD_flt_t delr2z = SIMD_load(tdelz + kcoffset);
	    const SIMD_flt_t rsq2 = SIMD_load(trsq + kcoffset);

	    const SIMD_flt_t r2 = SIMD_sqrt(rsq2);
	    const SIMD_flt_t rinvsq2 = SIMD_rcp(rsq2);
	    const SIMD_flt_t rainv2 = SIMD_rcp(r2 - cut);
	    const SIMD_flt_t gsrainv2 = sigma_gamma * rainv2;
	    const SIMD_flt_t gsrainvsq2 = gsrainv2 * rainv2 / r2;
	    const SIMD_flt_t expgsrainv2 = SIMD_exp(gsrainv2);
	    const SIMD_flt_t rinv12 = SIMD_rcp(r1 * r2);
	    const SIMD_flt_t cs = (delx * delr2x + dely * delr2y + 
                              delz * delr2z) * rinv12;
	    const SIMD_flt_t delcs = cs - costheta;
	    const SIMD_flt_t delcssq = delcs*delcs;

	    const SIMD_flt_t facexp = expgsrainv1*expgsrainv2;
	    const SIMD_flt_t facrad = lambda_epsilon * facexp * delcssq;
	    const SIMD_flt_t frad1 = facrad * gsrainvsq1;
	    const SIMD_flt_t frad2 = facrad * gsrainvsq2;
	    const SIMD_flt_t facang = lambda_epsilon2 * facexp * delcs;
	    const SIMD_flt_t facang12 = rinv12 * facang;
	    const SIMD_flt_t csfacang = cs * facang;

	    const SIMD_flt_t csfac1 = rinvsq1 * csfacang;
	    const SIMD_flt_t fjx = delx * (frad1 + csfac1)-delr2x*facang12;
	    const SIMD_flt_t fjy = dely * (frad1 + csfac1)-delr2y*facang12;
	    const SIMD_flt_t fjz = delz * (frad1 + csfac1)-delr2z*facang12;

	    const SIMD_flt_t csfac2 = rinvsq2 * csfacang;
	    SIMD_flt_t fkx = delx * facang12 - delr2x * (frad2 + csfac2);
	    SIMD_flt_t fky = dely * facang12 - delr2y * (frad2 + csfac2);
	    SIMD_flt_t fkz = delz * facang12 - delr2z * (frad2 + csfac2);

	    const SIMD_mask kmask = SIMD_lt(jmask, kk, ejnum);

	    SIMD_acc_cache3(kmask, fjx, fjy, fjz, fkx, fky, fkz, fxtmp, fytmp,
			    fztmp, fjxtmp, fjytmp, fjztmp, fxtmp2, fytmp2,
			    fztmp2, fjxtmp2, fjytmp2, fjztmp2, 
			    tf + kcoffset * 3, swidth); 

	    if (EVFLAG) {
	      if (EFLAG) {
		SIMD_int k;
		if (eatom) {
		  k = SIMD_load(tj + kcoffset);
		  k = k << 4;
		}
		SIMD_acc_three(kmask, facrad, eatom, sevdwl, fwtmp, fjtmp,
			       fwtmp2, fjtmp2, k, &(f[0].x));
	      }
	      SIMD_ev_tally_nbor3v(kmask, vflag, fjx, fjy, fjz, fkx, fky, fkz,
				   delx, dely, delz, delr2x, delr2y, delr2z,
				   sv0, sv1, sv2, sv3, sv4, sv5);
	    }
	    
	  } // for kk
	  if (is_same<flt_t,acc_t>::value == 1)
	    SIMD_cache3(tf + coffset * 3, swidth, fjxtmp, fjytmp, fjztmp);
	  else
	    SIMD_cache3(tf + coffset * 3, swidth, fjxtmp, fjytmp, fjztmp, 
	                fjxtmp2, fjytmp2, fjztmp2);

	  if (EFLAG) {
	    if (eatom) { 
	      SIMD_int j = SIMD_load(tj + coffset);
	      j = j << 4;
	      SIMD_jeng_update(jmask, &(f[0].x) + 3, j, fjtmp);
	      if (is_same<flt_t,acc_t>::value == 0)
		SIMD_jeng_update_hi(jmask, &(f[0].x) + 3, j, fjtmp2);
	    }
	  }
        } // for jj first loop

        for (int jj = 0; jj < ejnum_max; jj++) {
	  const int coffset = jj * swidth;
	  const SIMD_mask jmask = jj < ejnum;
          const SIMD_int j = SIMD_load(tj + coffset);
	  const SIMD_int joffset = j << 4;

	  SIMD_acc_t fjxtmp, fjytmp, fjztmp, fjxtmp2, fjytmp2, fjztmp2;
	  int foffset = swidth;
	  if (is_same<flt_t,acc_t>::value == 0) foffset = foffset >> 1;
	  acc_t *p = tf + coffset * 3;
	  fjxtmp = SIMD_load(p);
	  if (is_same<flt_t,acc_t>::value == 0) {
	    p = p + foffset;
	    fjxtmp2 = SIMD_load(p);
	  }
	  p = p + foffset;
	  fjytmp = SIMD_load(p);
	  if (is_same<flt_t,acc_t>::value == 0) {
	    p = p + foffset;
	    fjytmp2 = SIMD_load(p);
	  }
	  p = p + foffset;
	  fjztmp = SIMD_load(p);
	  if (is_same<flt_t,acc_t>::value == 0) {
	    p = p + foffset;
	    fjztmp2 = SIMD_load(p);
	  }
	  
	  SIMD_conflict_pi_reduce3(jmask, joffset, fjxtmp, fjytmp, fjztmp);
	  SIMD_jforce_update(jmask, &(f[0].x), joffset, fjxtmp, fjytmp, 
			     fjztmp);
          if (is_same<flt_t,acc_t>::value == 0) {
	    SIMD_int joffset2 = _mm512_shuffle_i32x4(joffset, joffset, 238);
	    SIMD_mask jmask2 = jmask >> 8;
	    SIMD_conflict_pi_reduce3(jmask2, joffset2, fjxtmp2, fjytmp2, 
				     fjztmp2);
	    SIMD_jforce_update(jmask2, &(f[0].x), joffset2, fjxtmp2, fjytmp2, 
			       fjztmp2);
	  }
	} // for jj second loop

	SIMD_iforce_update(imask, &(f[i].x), goffset, fxtmp, fytmp, fztmp,
			   EVFLAG, eatom, fwtmp);
	if (is_same<flt_t,acc_t>::value == 0) {
	  imask = imask >> 8;
	  SIMD_iforce_update(imask, &(f[i+8].x), goffset, fxtmp2, fytmp2, 
			     fztmp2, EVFLAG, eatom, fwtmp2);
	}
	if (EVFLAG) {
	  if (EFLAG) oevdwl += SIMD_sum(sevdwl);
	  if (vflag == 1) {
	    ov0 += SIMD_sum(sv0);
	    ov1 += SIMD_sum(sv1);
	    ov2 += SIMD_sum(sv2);
	    ov3 += SIMD_sum(sv3);
	    ov4 += SIMD_sum(sv4);
	    ov5 += SIMD_sum(sv5);
	  }
	}
	ilist = ilist + swidth;
      } // for ii

      #ifndef _LMP_INTEL_OFFLOAD
      if (vflag == 2)
      #endif
      {
        #if defined(_OPENMP)
        #pragma omp barrier
        #endif
        IP_PRE_fdotr_acc_force(1, EVFLAG,  EFLAG, vflag, eatom, nall, nlocal, 
			       minlocal, nthreads, f_start, f_stride, x, 
			       offload);
      }
    } // end omp
  
    if (EVFLAG) {
      if (EFLAG) {
        ev_global[0] = oevdwl;
        ev_global[1] = (acc_t)0.0;
      }
      if (vflag) {
        ev_global[2] = ov0;
        ev_global[3] = ov1;
        ev_global[4] = ov2;
        ev_global[5] = ov3;
        ev_global[6] = ov4;
        ev_global[7] = ov5;
      }
    }
    #if defined(__MIC__) && defined(_LMP_INTEL_OFFLOAD)
    *timer_compute = MIC_Wtime() - *timer_compute;
    #endif
  } // end offload
  if (offload)
    fix->stop_watch(TIME_OFFLOAD_LATENCY);
  else
    fix->stop_watch(TIME_HOST_PAIR);

  if (EVFLAG)
    fix->add_result_array(f_start, ev_global, offload, eatom, 0, vflag);
  else
    fix->add_result_array(f_start, 0, offload);
}

#endif

/* ---------------------------------------------------------------------- */

void PairSWIntel::allocate()
{
  PairSW::allocate();
}

/* ----------------------------------------------------------------------
   init specific to this pair style
------------------------------------------------------------------------- */

void PairSWIntel::init_style()
{
  PairSW::init_style();
  neighbor->requests[neighbor->nrequest-1]->intel = 1;
  map[0] = map[1];

  int ifix = modify->find_fix("package_intel");
  if (ifix < 0)
    error->all(FLERR,
               "The 'package intel' command is required for /intel styles");
  fix = static_cast<FixIntel *>(modify->fix[ifix]);

  fix->pair_init_check(true);

  #ifdef _LMP_INTEL_OFFLOAD
  _cop = fix->coprocessor_number();
  #endif

  if (fix->precision() == FixIntel::PREC_MODE_MIXED) {
    pack_force_const(force_const_single, fix->get_mixed_buffers());
    fix->get_mixed_buffers()->need_tag(1);
  } else if (fix->precision() == FixIntel::PREC_MODE_DOUBLE) {
    pack_force_const(force_const_double, fix->get_double_buffers());
    fix->get_double_buffers()->need_tag(1);
  } else {
    pack_force_const(force_const_single, fix->get_single_buffers());
    fix->get_single_buffers()->need_tag(1);
  }

  #ifdef _LMP_INTEL_OFFLOAD
  if (fix->offload_noghost())
    error->all(FLERR,"The 'ghost no' option cannot be used with sw/intel.");
  #endif

  #if defined(__INTEL_COMPILER)
  if (__INTEL_COMPILER_BUILD_DATE < 20141023)
    error->all(FLERR, "Intel compiler versions before "
	       "15 Update 1 not supported for sw/intel");
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t, class acc_t>
void PairSWIntel::pack_force_const(ForceConst<flt_t> &fc,
                                   IntelBuffers<flt_t,acc_t> *buffers)
{
  #ifdef LMP_USE_AVXCD
  fix->nbor_pack_width(SIMD_type<flt_t>::width());
  #endif

  int off_ccache = 0;
  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop >= 0) off_ccache = 1;
  #endif

  #ifdef LMP_USE_AVXCD
  const int swidth = SIMD_type<flt_t>::width();
  #else
  const int swidth = 1;
  #endif

  buffers->grow_ccache(off_ccache, comm->nthreads, swidth);
  _ccache_stride = buffers->ccache_stride();
  #ifdef LMP_USE_AVXCD
  _ccache_stride3 = buffers->ccache_stride3();
  #endif

  int tp1 = atom->ntypes + 1;
  fc.set_ntypes(tp1,memory,_cop);
  buffers->set_ntypes(tp1);
  flt_t **cutneighsq = buffers->get_cutneighsq();

  // Repeat cutsq calculation because done after call to init_style
  double cut, cutneigh;
  for (int i = 1; i <= atom->ntypes; i++) {
    for (int j = i; j <= atom->ntypes; j++) {
      if (setflag[i][j] != 0 || (setflag[i][i] != 0 && setflag[j][j] != 0)) {
        cut = init_one(i,j);
        cutneigh = cut + neighbor->skin;
        cutsq[i][j] = cutsq[j][i] = cut*cut;
        cutneighsq[i][j] = cutneighsq[j][i] = cutneigh * cutneigh;
      }
    }
  }
  
  _onetype = 0;
  if (atom->ntypes == 1) _onetype = 1;

  _spq = 1;
  for (int ii = 0; ii < tp1; ii++) {
    int i = map[ii];
    for (int jj = 0; jj < tp1; jj++) {
      int j = map[jj];
      if (i < 0 || j < 0 || ii == 0 || jj == 0) {
	fc.p2[ii][jj].cutsq = 0;
	fc.p2[ii][jj].cut = 0;
	fc.p2[ii][jj].sigma_gamma = 0;
	fc.p2f[ii][jj].cut = 0;
	fc.p2f[ii][jj].powerp = 0;
	fc.p2f[ii][jj].powerq = 0;
	fc.p2f[ii][jj].sigma = 0;
	fc.p2f2[ii][jj].c1 = 0;
	fc.p2f2[ii][jj].c2 = 0;
	fc.p2f2[ii][jj].c3 = 0;
	fc.p2f2[ii][jj].c4 = 0;
	fc.p2e[ii][jj].c5 = 0;
	fc.p2e[ii][jj].c6 = 0;
      } else {
	int ijparam = elem2param[i][j][j];
	fc.p2[ii][jj].cutsq = params[ijparam].cutsq;
	fc.p2[ii][jj].cut = params[ijparam].cut;
	fc.p2[ii][jj].sigma_gamma = params[ijparam].sigma_gamma;
	fc.p2f[ii][jj].cut = params[ijparam].cut;
	fc.p2f[ii][jj].powerp = -params[ijparam].powerp;
	fc.p2f[ii][jj].powerq = -params[ijparam].powerq;
	fc.p2f[ii][jj].sigma = params[ijparam].sigma;
	fc.p2f2[ii][jj].c1 = params[ijparam].c1;
	fc.p2f2[ii][jj].c2 = params[ijparam].c2;
	fc.p2f2[ii][jj].c3 = params[ijparam].c3;
	fc.p2f2[ii][jj].c4 = params[ijparam].c4;
	fc.p2e[ii][jj].c5 = params[ijparam].c5;
	fc.p2e[ii][jj].c6 = params[ijparam].c6;

	double cutcut = params[ijparam].cut * params[ijparam].cut;
	if (params[ijparam].cutsq >= cutcut)
	  fc.p2[ii][jj].cutsq *= 0.98;

	if (params[ijparam].powerp != 4.0 || params[ijparam].powerq != 0.0)
	  _spq = 0;
      }

      for (int kk = 0; kk < tp1; kk++) {
        int k = map[kk];
	if (i < 0 || j < 0 || k < 0  || ii == 0 || jj == 0 || kk == 0) {
	  fc.p3[ii][jj][kk].costheta = 0;
	  fc.p3[ii][jj][kk].lambda_epsilon = 0;
	  fc.p3[ii][jj][kk].lambda_epsilon2 = 0;
	} else {
	  int ijkparam = elem2param[i][j][k];
	  fc.p3[ii][jj][kk].costheta = params[ijkparam].costheta;
	  fc.p3[ii][jj][kk].lambda_epsilon = params[ijkparam].lambda_epsilon;
	  fc.p3[ii][jj][kk].lambda_epsilon2 = params[ijkparam].lambda_epsilon2;
	}
      }
    }
  }

  _host_pad = 1;
  _offload_pad = 1;

  if (INTEL_NBOR_PAD > 1)
    _host_pad = INTEL_NBOR_PAD * sizeof(float) / sizeof(flt_t);

  #ifdef _LMP_INTEL_OFFLOAD
  if (_cop < 0) return;
  FC_PACKED0_T *op2 = fc.p2[0];
  FC_PACKED1_T *op2f = fc.p2f[0];
  FC_PACKED1p2_T *op2f2 = fc.p2f2[0];
  FC_PACKED2_T *op2e = fc.p2e[0];
  FC_PACKED3_T *op3 = fc.p3[0][0];
  flt_t * ocutneighsq = cutneighsq[0];
  int tp1sq = tp1 * tp1;
  int tp1cu = tp1sq * tp1;
  if (op2 != NULL && op2f != NULL && op2f2 != NULL && op2e != NULL && 
      op3 != NULL && ocutneighsq != NULL) {
    #pragma offload_transfer target(mic:_cop) \
      in(op2,op2f,op2f2,op2e: length(tp1sq) alloc_if(0) free_if(0))	\
      in(op3: length(tp1cu) alloc_if(0) free_if(0)) \
      in(ocutneighsq: length(tp1sq))
  }
  #endif
}

/* ---------------------------------------------------------------------- */

template <class flt_t>
void PairSWIntel::ForceConst<flt_t>::set_ntypes(const int ntypes,
                                                Memory *memory,
                                                const int cop) {
  if (ntypes != _ntypes) {
    if (_ntypes > 0) {
      fc_packed0 *op2 = p2[0];
      fc_packed1 *op2f = p2f[0];
      fc_packed1p2 *op2f2 = p2f2[0];
      fc_packed2 *op2e = p2e[0];
      fc_packed3 *op3 = p3[0][0];

      #ifdef _LMP_INTEL_OFFLOAD
      if (op2 != NULL && op2f != NULL && op2f2 != NULL && op2e != NULL && 
	  op3 != NULL && _cop >= 0) {
        #pragma offload_transfer target(mic:_cop) \
          nocopy(op2, op2f, op2f2, op2e, op3: alloc_if(0) free_if(1))
      }
      #endif

      _memory->destroy(op2);
      _memory->destroy(op2f);
      _memory->destroy(op2f2);
      _memory->destroy(op2e);
      _memory->destroy(op3);
    }
    if (ntypes > 0) {
      _cop = cop;
      memory->create(p2,ntypes,ntypes,"fc.p2");
      memory->create(p2f,ntypes,ntypes,"fc.p2f");
      memory->create(p2f2,ntypes,ntypes,"fc.p2f2");
      memory->create(p2e,ntypes,ntypes,"fc.p2e");
      memory->create(p3,ntypes,ntypes,ntypes,"fc.p3");

      #ifdef _LMP_INTEL_OFFLOAD
      fc_packed0 *op2 = p2[0];
      fc_packed1 *op2f = p2f[0];
      fc_packed1p2 *op2f2 = p2f2[0];
      fc_packed2 *op2e = p2e[0];
      fc_packed3 *op3 = p3[0][0];
      int tp1sq = ntypes * ntypes;
      int tp1cu = tp1sq * ntypes;
      if (op2 != NULL && op2f != NULL && op2f2 != NULL && op2e != NULL && 
	  op3 != NULL && cop >= 0) {
        #pragma offload_transfer target(mic:cop) \
          nocopy(op2,op2f,op2f2,op2e: length(tp1sq) alloc_if(1) free_if(0)) \
          nocopy(op3: length(tp1cu) alloc_if(1) free_if(0))
      }
      #endif
    }
  }
  _ntypes = ntypes;
  _memory = memory;
}
