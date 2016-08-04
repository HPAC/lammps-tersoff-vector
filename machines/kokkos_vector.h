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
   Contributing author: Markus Höhnerbach (RWTH)
------------------------------------------------------------------------- */

// This file provides an intrinsics abstraction that allows access to the
// underlying SIMD registers.
// It allows for algorithms that are templated over the used vector length
// The final interface is provided by vector_routines, which provides the
// support for different precision modes.
// The vector_ops interface provides routines specific to one floating point
// data type, and is specialized for various architectures.
// The routines work best with AVX-512 and AVX2, as both support the gather
// instructions.
// For both AVX and SSE we miss some optimization opportunities in the gather
// implementations.

// Vector classes provided with the intel compiler
namespace lmp_intel {

#ifdef __CUDA_ARCH__
#define KOKKOS_VECTOR_INTRINSIC static __device__ __host__
#else
#define KOKKOS_VECTOR_INTRINSIC static
#endif
// Self explanatory mostly, KNC=IMCI and AVX-512, NONE=Scalar, AN=Array Not.
enum CalculationMode { KNC, AVX, AVX2, SSE, NONE, AN, CUDA };
#ifdef __MIC__
  #ifdef LMP_INTEL_VECTOR_MIC
  KOKKOS_VECTOR_INTRINSIC const CalculationMode mode = LMP_INTEL_VECTOR_MIC;
  #else
  KOKKOS_VECTOR_INTRINSIC const CalculationMode mode = KNC;
  #endif
#else
  static const CalculationMode mode = NONE;
  #ifdef LMP_INTEL_VECTOR_HOST
  KOKKOS_VECTOR_INTRINSIC const CalculationMode mode = LMP_INTEL_VECTOR_HOST;
  #else
  //  #ifdef __AVX512F__
  //  KOKKOS_VECTOR_INTRINSIC const CalculationMode mode = KNC;
  //  #else
  //    #ifdef __AVX2__
  //    KOKKOS_VECTOR_INTRINSIC const CalculationMode mode = AVX2;
  //    #else
  //      #ifdef __AVX__
  //      KOKKOS_VECTOR_INTRINSIC const CalculationMode mode = AVX;
  //      #else
  //      KOKKOS_VECTOR_INTRINSIC const CalculationMode mode = SSE;
  //      #endif
  //    #endif
  //  #endif
  //#endif
  #endif
#endif

// This is used in the selection logic
template<CalculationMode mode>
struct vector_traits { 
    static const bool support_integer_and_gather_ops = true; 
};

template<>
struct vector_traits<AVX> { 
    static const bool support_integer_and_gather_ops = false; 
};

// This is the base template for all the different architectures
// It will get specialized
template<class flt_t, CalculationMode mode>
struct vector_ops {};

// Scalar implementation
template<class flt_t>
struct vector_ops<flt_t, NONE> {
    static const int VL = 1;
    static const int ALIGN = 4;
    typedef flt_t fscal;
    typedef flt_t fvec;
    typedef int ivec;
    typedef bool bvec;
    typedef flt_t farr[1];
    typedef int iarr[1];
    KOKKOS_VECTOR_INTRINSIC fvec recip(const fvec &a) {
      return ((flt_t) 1.) / a;
    }
    template<int scale>
    KOKKOS_VECTOR_INTRINSIC void gather_prefetch_t0(const ivec &idx, const bvec &mask, const void *base) {
      // nop
    }
    template<int scale>
    KOKKOS_VECTOR_INTRINSIC fvec gather(const fvec &from, const bvec &mask, const ivec &idx, const void *base) {
      return mask ? *reinterpret_cast<const flt_t*>(reinterpret_cast<const char*>(base) + scale * idx) : from;
    }
    template<class T>
    KOKKOS_VECTOR_INTRINSIC void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      *x = gather<1>(*x, mask, idxs, &base->x);
      *y = gather<1>(*y, mask, idxs, &base->y);
      *z = gather<1>(*z, mask, idxs, &base->z);
      *w = int_gather<1>(*w, mask, idxs, &base->w);
    }
    KOKKOS_VECTOR_INTRINSIC void gather_8(const ivec &idxs, const bvec &mask, const void *base, 
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      fvec a = zero(), b = zero(), c = zero(), d = zero();
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 4 * sizeof(fscal), r4, r5, r6, r7);
    }
    KOKKOS_VECTOR_INTRINSIC void gather_4(const ivec &idxs, const bvec &mask, const void *base, 
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char*>(base) +  0 * sizeof(fscal));
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char*>(base) +  1 * sizeof(fscal));
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char*>(base) +  2 * sizeof(fscal));
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char*>(base) +  3 * sizeof(fscal));
    }
    KOKKOS_VECTOR_INTRINSIC fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return mask ? b : a;
    }
    KOKKOS_VECTOR_INTRINSIC ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return mask ? b : a;
    }
    KOKKOS_VECTOR_INTRINSIC fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    KOKKOS_VECTOR_INTRINSIC fvec zero() {
      return 0.;
    }
    KOKKOS_VECTOR_INTRINSIC bvec cmpeq(const fvec &a, const fvec &b) {
      return a == b;
    }
    KOKKOS_VECTOR_INTRINSIC bvec cmpnle(const fvec &a, const fvec &b) {
      return !(a <= b);
    }
    KOKKOS_VECTOR_INTRINSIC bvec cmple(const fvec &a, const fvec &b) {
      return a <= b;
    }
    KOKKOS_VECTOR_INTRINSIC bvec cmplt(const fvec &a, const fvec &b) {
      return a < b;
    }
    KOKKOS_VECTOR_INTRINSIC bvec int_cmpneq(const ivec &a, const ivec &b) {
      return a != b;
    }
    KOKKOS_VECTOR_INTRINSIC bvec int_cmplt(const ivec &a, const ivec &b) {
      return a < b;
    }
    KOKKOS_VECTOR_INTRINSIC fvec invsqrt(const fvec &a) {
      return 1. / sqrt(a);
    }
    KOKKOS_VECTOR_INTRINSIC fvec sincos(fvec *c, const fvec &a) {
      *c = cos(a);
      return sin(a);
    }
    KOKKOS_VECTOR_INTRINSIC fscal reduce_add(const fvec &a) {
      return a;
    }
    KOKKOS_VECTOR_INTRINSIC ivec int_mullo(const ivec &a, const ivec &b) {
      return a * b;
    }
    KOKKOS_VECTOR_INTRINSIC ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return mask ? a + b : src;
    }
    template<int scale>
    KOKKOS_VECTOR_INTRINSIC ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      return mask ? *reinterpret_cast<const int*>(reinterpret_cast<const char*>(base) + scale * idx) : from;
    }
    KOKKOS_VECTOR_INTRINSIC fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return mask ? a + b : src;
    }
    KOKKOS_VECTOR_INTRINSIC void store(void *at, const fvec &a) {
      *reinterpret_cast<flt_t*>(at) = a;
    }
    KOKKOS_VECTOR_INTRINSIC void int_store(int *at, const ivec &a) {
      *reinterpret_cast<int*>(at) = a;
    }
    KOKKOS_VECTOR_INTRINSIC void mask_store(int *at, const bvec &a) {
      *at = a;
    }
    KOKKOS_VECTOR_INTRINSIC fvec min(const fvec &a, const fvec &b) {
      return a < b ? a : b;
    }
    KOKKOS_VECTOR_INTRINSIC bool mask_test_at(const bvec &mask, int at) {
      return mask;
    }
    KOKKOS_VECTOR_INTRINSIC bool mask_testz(const bvec &mask) {
#ifdef __CUDA_ARCH__
      return __all(! mask);
#else
      return ! mask;
#endif
    }
    KOKKOS_VECTOR_INTRINSIC bvec mask_enable_lower(int n) {
      return n > 0 ? true : false;
    }
    KOKKOS_VECTOR_INTRINSIC ivec int_load_vl(const int *a) {
      return *a;
    }
    KOKKOS_VECTOR_INTRINSIC void int_clear_arr(int *a) {
      *a = 0;
    }
    KOKKOS_VECTOR_INTRINSIC bvec full_mask() {
      return true;
    }
    KOKKOS_VECTOR_INTRINSIC void int_print(const ivec &a) {
    }
};

// Mixins to implement mixed precision and single/single and double/double
// This one is for single/single and double/double
template<class BASE_flt_t, CalculationMode BASE_mic>
struct AccumulatorSameMixin {
  typedef vector_ops<BASE_flt_t, BASE_mic> BASE;
  typedef typename BASE::fvec avec;
  typedef typename BASE::farr aarr;

  KOKKOS_VECTOR_INTRINSIC avec acc_mask_add(const avec &src, const typename BASE::bvec &m, const avec &a, const typename BASE::fvec &b) {
    return BASE::mask_add(src, m, a, b);
  }

  KOKKOS_VECTOR_INTRINSIC typename BASE::fscal acc_reduce_add(const avec &a) {
    return BASE::reduce_add(a);
  }

  KOKKOS_VECTOR_INTRINSIC avec acc_zero() {
    return BASE::zero();
  }

  KOKKOS_VECTOR_INTRINSIC void acc_store(aarr mem, const avec &a) {
    BASE::store(mem, a);
  }

};

// Mixed precision for cases where double vectors contain fewer elements
template<class BASE_flt_t, class HIGH_flt_t, CalculationMode mic>
struct AccumulatorTwiceMixin {
  typedef vector_ops<BASE_flt_t, mic> BASE;
  typedef vector_ops<HIGH_flt_t, mic> HIGH;

  struct avec_t {
    typename HIGH::fvec lo, hi;
    avec_t(const typename HIGH::fvec &alo, const typename HIGH::fvec &ahi) : lo(alo), hi(ahi) {}
    avec_t(const typename BASE::fvec &a) {
      lo = BASE::cvtup_lo(a);
      hi = BASE::cvtup_hi(a);
    }
    friend avec_t operator +(const avec_t &a, const avec_t &b) {
      return avec_t(a.lo + b.lo, a.hi + b.hi);
    }
    friend avec_t operator -(const avec_t &a, const avec_t &b) {
      return avec_t(a.lo - b.lo, a.hi - b.hi);
    }
    friend avec_t operator *(const avec_t &a, const avec_t &b) {
      return avec_t(a.lo * b.lo, a.hi * b.hi);
    }
    operator typename BASE::fvec() const {
      return BASE::cvtdown(lo, hi);
    }
  };

  typedef avec_t avec;
  typedef typename HIGH::fscal aarr[BASE::VL] __attribute__((aligned(BASE::ALIGN)));
  
  KOKKOS_VECTOR_INTRINSIC avec acc_mask_add(const avec &src, const typename BASE::bvec &m, const avec &a, const typename BASE::fvec &b) {
    typename HIGH::fvec blo = BASE::cvtup_lo(b);
    typename HIGH::fvec bhi = BASE::cvtup_hi(b);
    typename HIGH::bvec mlo, mhi;
    BASE::mask_cvtup(m, &mlo, &mhi);
    return avec(HIGH::mask_add(src.lo, mlo, a.lo, blo), HIGH::mask_add(src.hi, mhi, a.hi, bhi));
  }
  
  KOKKOS_VECTOR_INTRINSIC typename HIGH::fscal acc_reduce_add(const avec &a) {
    return HIGH::reduce_add(a.lo + a.hi);
  }

  KOKKOS_VECTOR_INTRINSIC avec acc_zero() {
    return avec(HIGH::zero(), HIGH::zero());
  }

  KOKKOS_VECTOR_INTRINSIC void acc_store(aarr mem, const avec &a) {
    HIGH::store(mem, a.lo);
    HIGH::store(mem + BASE::VL / 2, a.hi);
  }

};

// For cases where vector_ops<float,x>::VL == vector_ops<double,x>::VL
// i.e. scalar & AN
template<class BASE_flt_t, class HIGH_flt_t, CalculationMode mic>
struct AccumulatorTwiceMixinNone {
  typedef vector_ops<BASE_flt_t, mic> BASE;
  typedef vector_ops<HIGH_flt_t, mic> HIGH;
 
  typedef typename HIGH::fvec avec;
  typedef typename HIGH::fscal aarr[BASE::VL];
  
  KOKKOS_VECTOR_INTRINSIC avec acc_mask_add(const avec &src, const typename BASE::bvec &m, const avec &a, const typename BASE::fvec &b) {
     return HIGH::mask_add(src, m, a, static_cast<typename HIGH::fvec>(b));
  }  
  KOKKOS_VECTOR_INTRINSIC typename HIGH::fscal acc_reduce_add(const avec &a) {
    return HIGH::reduce_add(a);
  }

  KOKKOS_VECTOR_INTRINSIC avec acc_zero() {
    return HIGH::zero();
  }

  KOKKOS_VECTOR_INTRINSIC void acc_store(aarr mem, const avec &a) {
    HIGH::store(mem, a);
  }

};

// This is the interfact that the user will see in the end.
template<class flt_t, class acc_t, CalculationMode mic>
struct vector_routines {};

template<CalculationMode mic>
struct vector_routines<double,double,mic> : public vector_ops<double, mic>, public AccumulatorSameMixin<double, mic> {};

template<CalculationMode mic>
struct vector_routines<float,float,mic> : public vector_ops<float, mic>, public AccumulatorSameMixin<float, mic> {};

template<CalculationMode mic>
struct vector_routines<float,double,mic> : public vector_ops<float, mic>, public AccumulatorTwiceMixin<float,double, mic> {};

// Specialize for AN and scalar
template<>
struct vector_routines<float,double,NONE> : public vector_ops<float, NONE>, public AccumulatorTwiceMixinNone<float,double, NONE> {};

} // namespace lmp_intel
