#ifdef __ALTIVEC__
#include<altivec.h>
#include "vector_math_altivec.h"

namespace altivec {
#pragma pack(push, 16)
  struct fvec {
    vector float data;  
    fvec() {}
    explicit fvec(const float f) { data = vec_splats(f); }
    explicit fvec(const float f[4]) { data = vec_ld(0, f); }
    //fvec(const fvec &a) { data = a.data; }
    fvec(const vector float &a) { data = a; }
    operator vector float() const { return data; }
    //fvec & operator = (const fvec &a) { data = a.data; return *this; }
    const fvec operator +(const fvec &b) const {
      return vec_add(data, b.data);
    }
    const fvec operator -(const fvec &b) const {
      return vec_sub(data, b.data);
    }
    const fvec operator *(const fvec &b) const {
      return vec_mul(data, b.data);
    }
    const fvec operator /(const fvec &b) const {
      return vec_div(data, b.data);
    }
    fvec& operator +=(const fvec &a) {
      data = vec_add(data, a.data); return *this;
    }
    fvec& operator -=(const fvec &a) {
      data = vec_sub(data, a.data); return *this;
    }
    fvec& operator *=(const fvec &a) {
      data = vec_mul(data, a.data); return *this;
    }
    fvec& operator /=(const fvec &a) {
      data = data / a.data; return *this;
    }
    fvec operator - () const {
        return - data;
    }
  };
  fvec invsqrt(const fvec &a) {
    //vector float out;
    //out = vec_insert(1.0f/sqrtf(vec_extract(a.data, 0)), out, 0);
    //out = vec_insert(1.0f/sqrtf(vec_extract(a.data, 1)), out, 1);
    //out = vec_insert(1.0f/sqrtf(vec_extract(a.data, 2)), out, 2);
    //out = vec_insert(1.0f/sqrtf(vec_extract(a.data, 3)), out, 3);
    //return out;
    return vec_rsqrt(a.data);
  }
  fvec sqrt(const fvec &a) {
    return vec_sqrt(a.data);
  }
  fvec exp(const fvec &a) {
    return altivec::exp_ps(a.data);
  }
  fvec sin(const fvec &a) {
    return altivec::sin_ps(a.data);
  }
  fvec pow(const fvec &a, const fvec &b) {
    fvec c = altivec::log_ps(a.data);
    fvec d = c * b;
    return altivec::exp_ps(d);
  }
  struct ivec {
    vector int data;
    ivec () {}
    explicit ivec(const int f) { data = vec_splats(f); }
    explicit ivec(const int f[4]) { data = vec_ld(0, f); }
    ivec(const vector int &d) { data = d; }
    operator vector int() const {
      return data;
    }
    const ivec operator &(const ivec &b) const {
      return vec_and(data, b.data);
    }
    const ivec operator |(const ivec &b) const {
      return vec_or(data, b.data);
    }
    const ivec operator +(const ivec &b) const {
      return vec_add(data, b.data);
    }
  };
  struct bvec {
    vector bool int data;
    bvec() {}
    explicit bvec(int i) {
      vector int a;
      if (i == 0) data = vec_cmpgt(a, a);
      else data = vec_cmpeq(a, a);
    }
    explicit bvec(const vector bool int &d) { data = d; }
    friend bvec operator &(const bvec &a, const bvec &b) {
      return bvec(vec_and(a.data, b.data));
    }
    friend bvec operator |(const bvec &a, const bvec &b) {
      return bvec(vec_or(a.data, b.data));
    }
    friend bvec operator ~(const bvec &a) {
      return bvec(vec_andc(bvec(1).data, a.data));
    }
    bvec& operator &=(const bvec &a) {
      data = vec_and(data, a.data);
      return *this;
    }
  };
#pragma pack(pop)
}
template<>
struct vector_ops<float, ALTIVEC> {
    static const int VL = 4;
    static const int ALIGN = 16;
    typedef float fscal;
    typedef altivec::fvec fvec;
    typedef altivec::ivec ivec;
    typedef altivec::bvec bvec;
    typedef float farr[4] __attribute__((aligned(16)));
    typedef int iarr[4] __attribute__((aligned(16)));
    static fvec recip(const fvec &a) {
      fvec b = vec_re(a.data);
      b = b + b - a * b * b;
      return  b + b - a * b * b;
      //return vec_recipdiv(fvec(1).data, a.data);
    }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, const bvec &mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, const bvec &mask, const ivec &idx, const void *base) {
      vector float result = from.data;
      #define ALTIVEC_GATHER_UNROLL(i)                                               \
        if (vec_extract(mask.data, i) != 0) {                                        \
          const char * raw_base = reinterpret_cast<const char*>(base);               \
          int idx_at_i = vec_extract(idx.data, i);                                   \
          size_t offset = scale * idx_at_i;                                          \
          const float * dest_ptr = reinterpret_cast<const float*>(&raw_base[offset]); \
          result = vec_insert(*dest_ptr, result, i);                                 \
        }
      ALTIVEC_GATHER_UNROLL(0) 
      ALTIVEC_GATHER_UNROLL(1) 
      ALTIVEC_GATHER_UNROLL(2) 
      ALTIVEC_GATHER_UNROLL(3) 
      return result;
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      vector float vm_0, vm_1, vm_2, vm_3;
      vector float vr_0, vr_1, vr_2, vr_3;
      #define ALTIVEC_GATHER_X_UNROLL(i) \
        if (vec_extract(mask.data, i) != 0)\
          vm_##i = vec_ld(0, reinterpret_cast<const float*>(&reinterpret_cast<const char*>(base)[1 * (idxs.data[i])]));
      // http://www.freevec.org/function/matrix_4x4_transpose_floats
      ALTIVEC_GATHER_X_UNROLL(0)
      ALTIVEC_GATHER_X_UNROLL(1)
      ALTIVEC_GATHER_X_UNROLL(2)
      ALTIVEC_GATHER_X_UNROLL(3)
      vr_0 = vec_mergeh(vm_0, vm_2);
      vr_1 = vec_mergel(vm_0, vm_2);
      vr_2 = vec_mergeh(vm_1, vm_3);
      vr_3 = vec_mergel(vm_1, vm_3);
      // Get the resulting vectors
      vm_0 = vec_mergeh(vr_0, vr_2);
      vm_1 = vec_mergel(vr_0, vr_2);
      vm_2 = vec_mergeh(vr_1, vr_3);
      vm_3 = vec_mergel(vr_1, vr_3);
      *x = blend(mask, *x, vm_0);
      *y = blend(mask, *y, vm_1);
      *z = blend(mask, *z, vm_2);
      *w = int_blend(mask, *w, (vector int)vm_3);
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base, 
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 4 * sizeof(fscal), r4, r5, r6, r7);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base, 
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      vector float vm_0, vm_1, vm_2, vm_3;
      vector float vr_0, vr_1, vr_2, vr_3;
      #define ALTIVEC_GATHER_4_UNROLL(i) \
        if (vec_extract(mask.data, i) != 0)\
          vm_##i = vec_ld(0, reinterpret_cast<const float*>(reinterpret_cast<const char*>(base) + 4 * (vec_extract(idxs.data, i))));
      ALTIVEC_GATHER_4_UNROLL(0)
      ALTIVEC_GATHER_4_UNROLL(1)
      ALTIVEC_GATHER_4_UNROLL(2)
      ALTIVEC_GATHER_4_UNROLL(3)
      vr_0 = vec_mergeh(vm_0, vm_2);
      vr_1 = vec_mergel(vm_0, vm_2);
      vr_2 = vec_mergeh(vm_1, vm_3);
      vr_3 = vec_mergel(vm_1, vm_3);
      // Get the resulting vectors
      vm_0 = vec_mergeh(vr_0, vr_2);
      vm_1 = vec_mergel(vr_0, vr_2);
      vm_2 = vec_mergeh(vr_1, vr_3);
      vm_3 = vec_mergel(vr_1, vr_3);
      *r0 = blend(mask, *r0, fvec(vm_0));
      *r1 = blend(mask, *r1, fvec(vm_1));
      *r2 = blend(mask, *r2, fvec(vm_2));
      *r3 = blend(mask, *r3, fvec(vm_3));
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return vec_or(vec_and((vector float) mask.data, b.data), vec_andc(a.data, (vector float) mask.data));
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return vec_or(vec_and(mask.data, b.data), vec_andc(a.data, mask.data));
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      vector float a;
      return vec_xor(a, a);
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return bvec(vec_cmpeq(a.data, b.data));
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return bvec(vec_cmpgt(a.data, b.data));
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return bvec(vec_cmple(a.data, b.data));
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return bvec(vec_cmplt(a.data, b.data));
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      return bvec(vec_or(vec_cmpgt(a.data, b.data), vec_cmplt(a.data, b.data)));
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return bvec(vec_cmplt(a.data, b.data));
    }
    static fvec invsqrt(const fvec &a) {
      return altivec::invsqrt(a);
    }
    static fvec sincos(fvec *c, const fvec &a) {
      fvec ret;
      altivec::sincos_ps(a.data, &ret.data, &c->data);
      return ret;
    }
    static fscal reduce_add(const fvec &a) {
      fscal sum = 0;
      sum += vec_extract(a.data, 0);
      sum += vec_extract(a.data, 1);
      sum += vec_extract(a.data, 2);
      sum += vec_extract(a.data, 3);
      return sum;
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
#     if __GNUC__ < 7
      vector int result;
      result = vec_insert(vec_extract(a.data, 0) * vec_extract(b.data, 0), result, 0);
      result = vec_insert(vec_extract(a.data, 1) * vec_extract(b.data, 1), result, 1);
      result = vec_insert(vec_extract(a.data, 2) * vec_extract(b.data, 2), result, 2);
      result = vec_insert(vec_extract(a.data, 3) * vec_extract(b.data, 3), result, 3);
      return result;
#     else
      return vec_mul(a.data, b.data);
#     endif
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return int_blend(mask, src, a + b);
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      vector int result = from.data;
      #define ALTIVEC_INT_GATHER_UNROLL(i)                                           \
        if (vec_extract(mask.data, i) != 0) {                                        \
          const char * raw_base = reinterpret_cast<const char*>(base);               \
          int idx_at_i = vec_extract(idx.data, i);                                   \
          size_t offset = scale * idx_at_i;                                          \
          const int * dest_ptr = reinterpret_cast<const int*>(&raw_base[offset]);     \
          result = vec_insert(*dest_ptr, result, i);                                 \
        }

      ALTIVEC_INT_GATHER_UNROLL(0) 
      ALTIVEC_INT_GATHER_UNROLL(1) 
      ALTIVEC_INT_GATHER_UNROLL(2) 
      ALTIVEC_INT_GATHER_UNROLL(3) 
      return result;
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return blend(mask, src, a + b);
    }
    static void store(void *at, const fvec &a) {
      vec_st(a.data, 0, reinterpret_cast<float*>(at));
    }
    static void int_store(int *at, const ivec &a) {
      vec_st(a.data, 0, at);
    }
    static fvec min(const fvec &a, const fvec &b) {
      return vec_min(a.data, b.data);
    }
    static bool mask_test_at(const bvec &mask, int at) {
      switch (at) {
      case 0: return vec_extract(mask.data, 0) != 0;
      case 1: return vec_extract(mask.data, 1) != 0;
      case 2: return vec_extract(mask.data, 2) != 0;
      case 3: return vec_extract(mask.data, 3) != 0;
      }
      return 0;
    }
    static bool mask_testz(const bvec &mask) {
      return vec_all_eq(mask.data, vec_splats(0));
//      if (vec_extract(mask.data, 0)) return false;
//      if (vec_extract(mask.data, 1)) return false;
//      if (vec_extract(mask.data, 2)) return false;
//      if (vec_extract(mask.data, 3)) return false;
//      return true;
    }
    static bvec mask_enable_lower(int n) {
      vector int a[5] __attribute__((aligned(16))) = {
        { 0, 0, 0, 0},
        {-1, 0, 0, 0},
        {-1,-1, 0, 0},
        {-1,-1,-1, 0},
        {-1,-1,-1,-1},
      };
      return bvec((vector bool int) a[n > 5 ? 5 : n]);
      //if (n > 0) a = vec_insert(-1, a, 0);
      //if (n > 1) a = vec_insert(-1, a, 1);
      //if (n > 2) a = vec_insert(-1, a, 2);
      //if (n > 3) a = vec_insert(-1, a, 3);
      //return bvec((vector bool int) a);
    }
    static ivec int_load_vl(const int *a) {
      return vec_ld(0, a);
    }
    static void int_clear_arr(int *a) {
      vec_st(vec_splats(0), 0, a);
    }
    static bvec full_mask() {
      return bvec(1);
    }
};

template<>
struct vector_ops<double, ALTIVEC> {
    static const int VL = 1;
    static const int ALIGN = 4;
    typedef double flt_t;
    typedef double fscal;
    typedef double fvec;
    typedef int ivec;
    typedef bool bvec;
    typedef double farr[1];
    typedef int iarr[1];
    static fvec recip(const fvec &a) {
      return ((flt_t) 1.) / a;
    }
    template<int scale>
    static void gather_prefetch_t0(const ivec &idx, const bvec &mask, const void *base) {
      // nop
    }
    template<int scale>
    static fvec gather(const fvec &from, const bvec &mask, const ivec &idx, const void *base) {
      return mask ? *reinterpret_cast<const flt_t*>(reinterpret_cast<const char*>(base) + scale * idx) : from;
    }
    template<class T>
    static void gather_x(const ivec &idxs, const bvec &mask, const T *base, fvec *x, fvec *y, fvec *z, ivec *w) {
      *x = gather<1>(*x, mask, idxs, &base->x);
      *y = gather<1>(*y, mask, idxs, &base->y);
      *z = gather<1>(*z, mask, idxs, &base->z);
      *w = int_gather<1>(*w, mask, idxs, &base->w);
    }
    static void gather_8(const ivec &idxs, const bvec &mask, const void *base, 
        fvec *r0, fvec *r1, fvec *r2, fvec *r3, fvec *r4, fvec *r5, fvec *r6, fvec *r7) {
      gather_4(idxs, mask, base, r0, r1, r2, r3);
      gather_4(idxs, mask, reinterpret_cast<const char*>(base) + 4 * sizeof(fscal), r4, r5, r6, r7);
    }
    static void gather_4(const ivec &idxs, const bvec &mask, const void *base, 
        fvec *r0, fvec *r1, fvec *r2, fvec *r3) {
      *r0 = gather<4>(*r0, mask, idxs, reinterpret_cast<const char*>(base) +  0 * sizeof(fscal));
      *r1 = gather<4>(*r1, mask, idxs, reinterpret_cast<const char*>(base) +  1 * sizeof(fscal));
      *r2 = gather<4>(*r2, mask, idxs, reinterpret_cast<const char*>(base) +  2 * sizeof(fscal));
      *r3 = gather<4>(*r3, mask, idxs, reinterpret_cast<const char*>(base) +  3 * sizeof(fscal));
    }
    static fvec blend(const bvec &mask, const fvec &a, const fvec &b) {
      return mask ? b : a;
    }
    static ivec int_blend(const bvec &mask, const ivec &a, const ivec &b) {
      return mask ? b : a;
    }
    static fvec fmadd(const fvec &a, const fvec &b, const fvec &c) {
      return a*b + c;
    }
    static fvec zero() {
      return 0.;
    }
    static bvec cmpeq(const fvec &a, const fvec &b) {
      return a == b;
    }
    static bvec cmpnle(const fvec &a, const fvec &b) {
      return !(a <= b);
    }
    static bvec cmple(const fvec &a, const fvec &b) {
      return a <= b;
    }
    static bvec cmplt(const fvec &a, const fvec &b) {
      return a < b;
    }
    static bvec int_cmpneq(const ivec &a, const ivec &b) {
      return a != b;
    }
    static bvec int_cmplt(const ivec &a, const ivec &b) {
      return a < b;
    }
    static fvec invsqrt(const fvec &a) {
      return 1. / sqrt(a);
    }
    static fvec sincos(fvec *c, const fvec &a) {
      *c = cos(a);
      return sin(a);
    }
    static fscal reduce_add(const fvec &a) {
      return a;
    }
    static ivec int_mullo(const ivec &a, const ivec &b) {
      return a * b;
    }
    static ivec int_mask_add(const ivec &src, const bvec &mask, const ivec &a, const ivec &b) {
      return mask ? a + b : src;
    }
    template<int scale>
    static ivec int_gather(const ivec &from, bvec mask, const ivec &idx, const void *base) {
      return mask ? *reinterpret_cast<const int*>(reinterpret_cast<const char*>(base) + scale * idx) : from;
    }
    static fvec mask_add(const fvec &src, const bvec &mask, const fvec &a, const fvec &b) {
      return mask ? a + b : src;
    }
    static void store(void *at, const fvec &a) {
      *reinterpret_cast<flt_t*>(at) = a;
    }
    static void int_store(int *at, const ivec &a) {
      *reinterpret_cast<int*>(at) = a;
    }
    static void mask_store(int *at, const bvec &a) {
      *at = a;
    }
    static fvec min(const fvec &a, const fvec &b) {
      return a < b ? a : b;
    }
    static bool mask_test_at(const bvec &mask, int at) {
      return mask;
    }
    static bool mask_testz(const bvec &mask) {
      return ! mask;
    }
    static bvec mask_enable_lower(int n) {
      return n > 0 ? true : false;
    }
    static ivec int_load_vl(const int *a) {
      return *a;
    }
    static void int_clear_arr(int *a) {
      *a = 0;
    }
    static bvec full_mask() {
      return true;
    }
    static void int_print(const ivec &a) {
    }
};

#endif
