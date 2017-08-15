/* NEON implementation of sin, cos, exp and log

   Inspired by Intel Approximate Math library, and based on the
   corresponding algorithms of the cephes math library
*/

/* Copyright (C) 2011  Julien Pommier

  This software is provided 'as-is', without any express or implied
  warranty.  In no event will the authors be held liable for any damages
  arising from the use of this software.

  Permission is granted to anyone to use this software for any purpose,
  including commercial applications, and to alter it and redistribute it
  freely, subject to the following restrictions:

  1. The origin of this software must not be misrepresented; you must not
     claim that you wrote the original software. If you use this software
     in a product, an acknowledgment in the product documentation would be
     appreciated but is not required.
  2. Altered source versions must be plainly marked as such, and must not be
     misrepresented as being the original software.
  3. This notice may not be removed or altered from any source distribution.

  (this is the zlib license)
*/

// Modified to be usable as a C++ header
#ifdef __ALTIVEC__

#include <altivec.h>

namespace altivec {

typedef vector float v4sf;  // vector of 4 float
typedef vector unsigned int v4su;  // vector of 4 uint32
typedef vector signed int v4si;  // vector of 4 uint32

#define c_inv_mant_mask ~0x7f800000u
#define c_cephes_SQRTHF 0.707106781186547524
#define c_cephes_log_p0 7.0376836292E-2
#define c_cephes_log_p1 - 1.1514610310E-1
#define c_cephes_log_p2 1.1676998740E-1
#define c_cephes_log_p3 - 1.2420140846E-1
#define c_cephes_log_p4 + 1.4249322787E-1
#define c_cephes_log_p5 - 1.6668057665E-1
#define c_cephes_log_p6 + 2.0000714765E-1
#define c_cephes_log_p7 - 2.4999993993E-1
#define c_cephes_log_p8 + 3.3333331174E-1
#define c_cephes_log_q1 -2.12194440e-4
#define c_cephes_log_q2 0.693359375

/* natural logarithm computed for 4 simultaneous float 
   return NaN for x <= 0
*/
inline v4sf log_ps(v4sf x) {
  //v4sf out;
  //out = vec_insert(log(vec_extract(x, 0)), out, 0);
  //out = vec_insert(log(vec_extract(x, 1)), out, 1);
  //out = vec_insert(log(vec_extract(x, 2)), out, 2);
  //out = vec_insert(log(vec_extract(x, 3)), out, 3);
  //return out;

  v4sf one = vec_splats(1.0f);

  x = vec_max(x, vec_splats(0.0f)); /* force flush to zero on denormal values */
  v4su invalid_mask = (v4su) vec_cmplt(x, vec_splats(0.0f));

  v4si ux = (v4si)(x);
  //
  v4si emm0 = vec_sr(ux, vec_splats(23u));

  /* keep only the fractional part */
  ux = vec_and(ux, vec_splats((signed int) c_inv_mant_mask));
  ux = vec_or(ux, (v4si)(vec_splats(0.5f)));
  x = (v4sf)(ux);

  emm0 = vec_sub(emm0, vec_splats(0x7f));
  v4sf e = vec_ctf(emm0, 0);

  e = vec_add(e, one);

  /* part2: 
     if( x < SQRTHF ) {
       e -= 1;
       x = x + x - 1.0;
     } else { x = x - 1.0; }
  */
  v4su mask = (v4su) vec_cmplt(x, vec_splats((float) c_cephes_SQRTHF));
  v4sf tmp = (v4sf)(vec_and((v4su)(x), mask));
  x = vec_sub(x, one);
  e = vec_sub(e, (v4sf)(vec_and((v4su)(one), mask)));
  x = vec_add(x, tmp);

  v4sf z = vec_mul(x,x);

  v4sf y = vec_splats((float) c_cephes_log_p0);
  y = vec_mul(y, x);
  y = vec_add(y, vec_splats((float) c_cephes_log_p1));
  y = vec_mul(y, x);
  y = vec_add(y, vec_splats((float) c_cephes_log_p2));
  y = vec_mul(y, x);
  y = vec_add(y, vec_splats((float) c_cephes_log_p3));
  y = vec_mul(y, x);
  y = vec_add(y, vec_splats((float) c_cephes_log_p4));
  y = vec_mul(y, x);
  y = vec_add(y, vec_splats((float) c_cephes_log_p5));
  y = vec_mul(y, x);
  y = vec_add(y, vec_splats((float) c_cephes_log_p6));
  y = vec_mul(y, x);
  y = vec_add(y, vec_splats((float) c_cephes_log_p7));
  y = vec_mul(y, x);
  y = vec_add(y, vec_splats((float) c_cephes_log_p8));
  y = vec_mul(y, x);

  y = vec_mul(y, z);
  

  tmp = vec_mul(e, vec_splats((float) c_cephes_log_q1));
  y = vec_add(y, tmp);


  tmp = vec_mul(z, vec_splats(0.5f));
  y = vec_sub(y, tmp);

  tmp = vec_mul(e, vec_splats((float) c_cephes_log_q2));
  x = vec_add(x, y);
  x = vec_add(x, tmp);
  x = (v4sf)(vec_or((v4su)(x), invalid_mask)); // negative arg will be NAN
  return x;
}

#define c_exp_hi 88.3762626647949f
#define c_exp_lo -88.3762626647949f

#define c_cephes_LOG2EF 1.44269504088896341
#define c_cephes_exp_C1 0.693359375
#define c_cephes_exp_C2 -2.12194440e-4

#define c_cephes_exp_p0 1.9875691500E-4
#define c_cephes_exp_p1 1.3981999507E-3
#define c_cephes_exp_p2 8.3334519073E-3
#define c_cephes_exp_p3 4.1665795894E-2
#define c_cephes_exp_p4 1.6666665459E-1
#define c_cephes_exp_p5 5.0000001201E-1

/* exp() computed for 4 float at once */
inline v4sf exp_ps(v4sf x) {
  //v4sf out;
  //out = vec_insert(exp(vec_extract(x, 0)), out, 0);
  //out = vec_insert(exp(vec_extract(x, 1)), out, 1);
  //out = vec_insert(exp(vec_extract(x, 2)), out, 2);
  //out = vec_insert(exp(vec_extract(x, 3)), out, 3);
  //return out;

  v4sf tmp, fx;

  v4sf one = vec_splats(1.0f);
  x = vec_min(x, vec_splats(c_exp_hi));
  x = vec_max(x, vec_splats(c_exp_lo));

  /* express exp(x) as exp(g + n*log(2)) */
  //fx = vmlaq_f32(vec_splats(0.5f), x, vec_splats((float) c_cephes_LOG2EF));
  fx = vec_madd(x, vec_splats((float) c_cephes_LOG2EF), vec_splats(0.5f));

  /* perform a floorf */
  tmp = vec_ctf(vec_cts(fx, 0), 0);

  /* if greater, substract 1 */
  v4su mask = (v4su) vec_cmpgt(tmp, fx);    
  mask = vec_and(mask, (v4su)(one));


  fx = vec_sub(tmp, (v4sf)(mask));

  tmp = vec_mul(fx, vec_splats((float) c_cephes_exp_C1));
  v4sf z = vec_mul(fx, vec_splats((float) c_cephes_exp_C2));
  x = vec_sub(x, tmp);
  x = vec_sub(x, z);

  v4sf y =  vec_splats((float) c_cephes_exp_p0);
  v4sf c1 = vec_splats((float) c_cephes_exp_p1); 
  v4sf c2 = vec_splats((float) c_cephes_exp_p2); 
  v4sf c3 = vec_splats((float) c_cephes_exp_p3); 
  v4sf c4 = vec_splats((float) c_cephes_exp_p4); 
  v4sf c5 = vec_splats((float) c_cephes_exp_p5);

  y = vec_mul(y, x);
  z = vec_mul(x,x);
  y = vec_add(y, c1);
  y = vec_mul(y, x);
  y = vec_add(y, c2);
  y = vec_mul(y, x);
  y = vec_add(y, c3);
  y = vec_mul(y, x);
  y = vec_add(y, c4);
  y = vec_mul(y, x);
  y = vec_add(y, c5);
  
  y = vec_mul(y, z);
  y = vec_add(y, x);
  y = vec_add(y, one);

  /* build 2^n */
  v4si mm;
  mm = vec_cts(fx, 0);
  mm = vec_add(mm, vec_splats(0x7f));
  mm = vec_sl(mm, vec_splats(23u));
  v4sf pow2n = (v4sf)(mm);

  y = vec_mul(y, pow2n);
  return y;
}

#define c_minus_cephes_DP1 -0.78515625
#define c_minus_cephes_DP2 -2.4187564849853515625e-4
#define c_minus_cephes_DP3 -3.77489497744594108e-8
#define c_sincof_p0 -1.9515295891E-4
#define c_sincof_p1  8.3321608736E-3
#define c_sincof_p2 -1.6666654611E-1
#define c_coscof_p0  2.443315711809948E-005
#define c_coscof_p1 -1.388731625493765E-003
#define c_coscof_p2  4.166664568298827E-002
#define c_cephes_FOPI 1.27323954473516 // 4 / M_PI

/* evaluation of 4 sines & cosines at once.

   The code is the exact rewriting of the cephes sinf function.
   Precision is excellent as long as x < 8192 (I did not bother to
   take into account the special handling they have for greater values
   -- it does not return garbage for arguments over 8192, though, but
   the extra precision is missing).

   Note that it is such that sinf((float)M_PI) = 8.74e-8, which is the
   surprising but correct result.

   Note also that when you compute sin(x), cos(x) is available at
   almost no extra price so both sin_ps and cos_ps make use of
   sincos_ps..
  */
inline void sincos_ps(v4sf x, v4sf *ysin, v4sf *ycos) { // any x
  //v4sf out;
  //out = vec_insert(sin(vec_extract(x, 0)), out, 0);
  //out = vec_insert(sin(vec_extract(x, 1)), out, 1);
  //out = vec_insert(sin(vec_extract(x, 2)), out, 2);
  //out = vec_insert(sin(vec_extract(x, 3)), out, 3);
  //*ysin = out;
  //out = vec_insert(cos(vec_extract(x, 0)), out, 0);
  //out = vec_insert(cos(vec_extract(x, 1)), out, 1);
  //out = vec_insert(cos(vec_extract(x, 2)), out, 2);
  //out = vec_insert(cos(vec_extract(x, 3)), out, 3);
  //*ycos = out;
  //return;
  v4sf xmm1, xmm2, xmm3, y;

  v4su emm2;
  
  v4su sign_mask_sin, sign_mask_cos;
  sign_mask_sin = (v4su) vec_cmplt(x, vec_splats(0.0f));
  x = vec_abs(x);

  /* scale by 4/Pi */
  y = vec_mul(x, vec_splats((float) c_cephes_FOPI));

  /* store the integer part of y in mm0 */
  emm2 = vec_ctu(y, 0);
  /* j=(j+1) & (~1) (see the cephes sources) */
  emm2 = vec_add(emm2, vec_splats(1u));
  emm2 = vec_and(emm2, vec_splats(~1u));
  y = vec_ctf(emm2, 0);

  /* get the polynom selection mask 
     there is one polynom for 0 <= x <= Pi/4
     and another one for Pi/4<x<=Pi/2

     Both branches will be computed.
  */
  v4su suzero = vec_splats(0u);
  // workaround tst and vec_cmpne
  #define vec_tstq_u32(a, b) vec_nor(suzero, (v4su) vec_cmpeq(suzero, vec_and(a, b)))
  v4su poly_mask = vec_tstq_u32(emm2, vec_splats(2u));
  
  /* The magic pass: "Extended precision modular arithmetic" 
     x = ((x - y * DP1) - y * DP2) - y * DP3; */
  xmm1 = vec_mul(y, vec_splats((float) c_minus_cephes_DP1));
  xmm2 = vec_mul(y, vec_splats((float) c_minus_cephes_DP2));
  xmm3 = vec_mul(y, vec_splats((float) c_minus_cephes_DP3));
  x = vec_add(x, xmm1);
  x = vec_add(x, xmm2);
  x = vec_add(x, xmm3);

  sign_mask_sin = vec_xor(sign_mask_sin, vec_tstq_u32(emm2, vec_splats(4u)));
  sign_mask_cos = vec_tstq_u32(vec_sub(emm2, vec_splats(2u)), vec_splats(4u));

  /* Evaluate the first polynom  (0 <= x <= Pi/4) in y1, 
     and the second polynom      (Pi/4 <= x <= 0) in y2 */
  v4sf z = vec_mul(x,x);
  v4sf y1, y2;

  y1 = vec_mul(z, vec_splats((float) c_coscof_p0));
  y2 = vec_mul(z, vec_splats((float) c_sincof_p0));
  y1 = vec_add(y1, vec_splats((float) c_coscof_p1));
  y2 = vec_add(y2, vec_splats((float) c_sincof_p1));
  y1 = vec_mul(y1, z);
  y2 = vec_mul(y2, z);
  y1 = vec_add(y1, vec_splats((float) c_coscof_p2));
  y2 = vec_add(y2, vec_splats((float) c_sincof_p2));
  y1 = vec_mul(y1, z);
  y2 = vec_mul(y2, z);
  y1 = vec_mul(y1, z);
  y2 = vec_mul(y2, x);
  y1 = vec_sub(y1, vec_mul(z, vec_splats(0.5f)));
  y2 = vec_add(y2, x);
  y1 = vec_add(y1, vec_splats(1.0f));

  /* select the correct result from the two polynoms */  
  // gcc workaround for vec_neg
  #define vec_negq_f32(a) vec_sub(vec_splats(0.0f), a)
  v4sf ys = vec_sel(y2, y1, poly_mask);
  v4sf yc = vec_sel(y1, y2, poly_mask);
  *ysin = vec_sel(ys, vec_negq_f32(ys), sign_mask_sin);
  *ycos = vec_sel(vec_negq_f32(yc), yc, sign_mask_cos);
//  v4sf ys = vbslq_f32(poly_mask, y1, y2);
//  v4sf yc = vbslq_f32(poly_mask, y2, y1);
//  *ysin = vbslq_f32(sign_mask_sin, vec_neg(ys), ys);
//  *ycos = vbslq_f32(sign_mask_cos, yc, vec_neg(yc));
}

inline v4sf sin_ps(v4sf x) {
  v4sf ysin, ycos; 
  sincos_ps(x, &ysin, &ycos); 
  return ysin;
}

inline v4sf cos_ps(v4sf x) {
  v4sf ysin, ycos; 
  sincos_ps(x, &ysin, &ycos); 
  return ycos;
}

}

#endif
