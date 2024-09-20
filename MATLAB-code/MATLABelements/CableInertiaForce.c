/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableInertiaForce.c
 *
 * Code generation for function 'CableInertiaForce'
 *
 */

/* Include files */
#include "CableInertiaForce.h"
#include "CableForce.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_internal_types.h"
#include "CableForceRotBCinCoord_types.h"
#include "mtimes.h"
#include "norm.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void L_Wp(const captured_var *tau0_x_tau,
                 const b_captured_var *tau_dot_tau0, const c_captured_var *I3,
                 const b_captured_var *nxip, const c_captured_var *G1,
                 const c_captured_var *PP, const c_captured_var *W,
                 const captured_var *tau, const c_captured_var *R0_,
                 const c_captured_var *Theta, const double upsilon[3],
                 double LL[9]);

static void Q1T_Wpp(const captured_var *tau0_x_tau,
                    const b_captured_var *tau_dot_tau0,
                    const c_captured_var *I3, const b_captured_var *nxip,
                    const c_captured_var *G1, const c_captured_var *PP,
                    const c_captured_var *W, const captured_var *tau,
                    const c_captured_var *R0_, const c_captured_var *Theta,
                    const captured_var *tau0, const double upsilon[3],
                    const double pie[3], double QQ[9]);

static void Q1T_Wthetp(const captured_var *tau0_x_tau,
                       const b_captured_var *tau_dot_tau0,
                       const c_captured_var *I3, const b_captured_var *nxip,
                       const c_captured_var *G1, const c_captured_var *PP,
                       const c_captured_var *W, const captured_var *tau,
                       const c_captured_var *R0_, const c_captured_var *Theta,
                       const captured_var *e1, const double upsilon[3],
                       const double pie[3], double QQ[3]);

static void Q1_Wpp(const captured_var *tau0_x_tau,
                   const b_captured_var *tau_dot_tau0, const c_captured_var *I3,
                   const b_captured_var *nxip, const c_captured_var *G1,
                   const c_captured_var *PP, const c_captured_var *W,
                   const captured_var *tau, const c_captured_var *R0_,
                   const c_captured_var *Theta, const captured_var *tau0,
                   const double upsilon[3], const double pie[3], double QQ[9]);

static void Q1_Wpthet(const captured_var *e1, const c_captured_var *W,
                      const captured_var *tau0_x_tau,
                      const b_captured_var *tau_dot_tau0,
                      const c_captured_var *I3, const b_captured_var *nxip,
                      const c_captured_var *G1, const c_captured_var *PP,
                      const captured_var *tau, const c_captured_var *R0_,
                      const c_captured_var *Theta, const double upsilon[3],
                      const double pie[3], double QQ[3]);

static void Q1_Wthetp(const captured_var *tau0_x_tau,
                      const b_captured_var *tau_dot_tau0,
                      const c_captured_var *I3, const b_captured_var *nxip,
                      const c_captured_var *G1, const c_captured_var *PP,
                      const c_captured_var *W, const captured_var *tau,
                      const c_captured_var *R0_, const c_captured_var *Theta,
                      const captured_var *e1, const double upsilon[3],
                      double pie, double QQ[9]);

static void Q2_Wp(const c_captured_var *PP, const c_captured_var *G1,
                  const b_captured_var *nxip, const captured_var *tau0_x_tau,
                  const b_captured_var *tau_dot_tau0, const captured_var *tau,
                  const c_captured_var *W, const c_captured_var *R0_,
                  const c_captured_var *Theta, const double pie[3],
                  double QQ[9]);

static void binary_expand_op_16(emxArray_real_T *in1,
                                const emxArray_int32_T *in2,
                                const emxArray_int32_T *in3,
                                const emxArray_uint32_T *in4,
                                const emxArray_real_T *in5);

static void binary_expand_op_17(emxArray_real_T *in1,
                                const emxArray_int32_T *in2, int in3, int in4,
                                const emxArray_uint32_T *in5, int in6, int in7,
                                const emxArray_real_T *in8);

static void binary_expand_op_22(emxArray_real_T *in1, int in2, int in3,
                                const emxArray_int32_T *in4, int in5, int in6,
                                const emxArray_uint32_T *in7,
                                const emxArray_real_T *in8);

/* Function Definitions */
static void L_Wp(const captured_var *tau0_x_tau,
                 const b_captured_var *tau_dot_tau0, const c_captured_var *I3,
                 const b_captured_var *nxip, const c_captured_var *G1,
                 const c_captured_var *PP, const c_captured_var *W,
                 const captured_var *tau, const c_captured_var *R0_,
                 const c_captured_var *Theta, const double upsilon[3],
                 double LL[9])
{
  double b_G1[9];
  double b_y[9];
  double c_G1[9];
  double y[9];
  double b_W[3];
  double a_tmp;
  double b_tau0_x_tau;
  double d;
  double d1;
  double d2;
  int aoffset;
  int coffset;
  int i;
  int j;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  FUNCTIONS FOR DERIVATIVE COMPUTATIONS */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  b_tau0_x_tau = 0.0;
  a_tmp = 1.0 / nxip->contents;
  for (j = 0; j < 3; j++) {
    double d3;
    coffset = j * 3;
    b_tau0_x_tau += tau0_x_tau->contents[j] * upsilon[j];
    d = 0.0;
    d1 = Theta->contents[j];
    d2 = Theta->contents[j + 3];
    d3 = Theta->contents[j + 6];
    for (i = 0; i < 3; i++) {
      aoffset = i * 3;
      y[coffset + i] =
          (R0_->contents[aoffset] * d1 + R0_->contents[aoffset + 1] * d2) +
          R0_->contents[aoffset + 2] * d3;
      d += W->contents[j + 3 * i] * upsilon[i];
    }
    b_W[j] = -d;
  }
  b_tau0_x_tau /= tau_dot_tau0->contents + 1.0;
  b_y[0] = 0.0;
  b_y[3] = -upsilon[2];
  b_y[6] = upsilon[1];
  b_y[1] = upsilon[2];
  b_y[4] = 0.0;
  b_y[7] = -upsilon[0];
  b_y[2] = -upsilon[1];
  b_y[5] = upsilon[0];
  b_y[8] = 0.0;
  for (coffset = 0; coffset < 9; coffset++) {
    b_y[coffset] += b_tau0_x_tau * I3->contents[coffset];
  }
  for (coffset = 0; coffset < 3; coffset++) {
    d = G1->contents[coffset];
    d1 = G1->contents[coffset + 3];
    d2 = G1->contents[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_G1[coffset + 3 * aoffset] =
          (d * b_y[3 * aoffset] + d1 * b_y[3 * aoffset + 1]) +
          d2 * b_y[3 * aoffset + 2];
    }
    d = b_G1[coffset];
    d1 = b_G1[coffset + 3];
    d2 = b_G1[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      c_G1[coffset + 3 * aoffset] =
          (d * PP->contents[3 * aoffset] + d1 * PP->contents[3 * aoffset + 1]) +
          d2 * PP->contents[3 * aoffset + 2];
    }
  }
  for (coffset = 0; coffset < 9; coffset++) {
    c_G1[coffset] *= a_tmp;
  }
  for (coffset = 0; coffset < 3; coffset++) {
    d = y[coffset];
    d1 = y[coffset + 3];
    d2 = y[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      LL[aoffset + 3 * coffset] = b_W[aoffset] * tau->contents[coffset];
      b_y[coffset + 3 * aoffset] =
          (d * c_G1[3 * aoffset] + d1 * c_G1[3 * aoffset + 1]) +
          d2 * c_G1[3 * aoffset + 2];
    }
  }
  for (coffset = 0; coffset < 9; coffset++) {
    LL[coffset] = a_tmp * (LL[coffset] + b_y[coffset]);
  }
}

static void Q1T_Wpp(const captured_var *tau0_x_tau,
                    const b_captured_var *tau_dot_tau0,
                    const c_captured_var *I3, const b_captured_var *nxip,
                    const c_captured_var *G1, const c_captured_var *PP,
                    const c_captured_var *W, const captured_var *tau,
                    const c_captured_var *R0_, const c_captured_var *Theta,
                    const captured_var *tau0, const double upsilon[3],
                    const double pie[3], double QQ[9])
{
  double b_G1[9];
  double b_t_hat___tmp[9];
  double b_y[9];
  double e_a[9];
  double f_a[9];
  double mm[9];
  double t_hat___tmp[9];
  double y[9];
  double G1Tpie[3];
  double b_W[3];
  double c_y[3];
  double d_a[3];
  double a;
  double a_tmp;
  double a_tmp_tmp;
  double b_a;
  double b_a_tmp;
  double b_tau;
  double b_tau0_x_tau;
  double c_a;
  double d;
  double d1;
  double d2;
  double d3;
  int aoffset;
  int coffset;
  int i;
  int j;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  FUNCTIONS FOR DERIVATIVE COMPUTATIONS */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  b_tau0_x_tau = 0.0;
  t_hat___tmp[0] = 0.0;
  t_hat___tmp[3] = -upsilon[2];
  t_hat___tmp[6] = upsilon[1];
  t_hat___tmp[1] = upsilon[2];
  t_hat___tmp[4] = 0.0;
  t_hat___tmp[7] = -upsilon[0];
  t_hat___tmp[2] = -upsilon[1];
  t_hat___tmp[5] = upsilon[0];
  t_hat___tmp[8] = 0.0;
  a_tmp = 1.0 / nxip->contents;
  for (j = 0; j < 3; j++) {
    coffset = j * 3;
    b_tau0_x_tau += tau0_x_tau->contents[j] * upsilon[j];
    d = 0.0;
    d1 = Theta->contents[j];
    d2 = Theta->contents[j + 3];
    d3 = Theta->contents[j + 6];
    for (i = 0; i < 3; i++) {
      aoffset = i * 3;
      y[coffset + i] =
          (R0_->contents[aoffset] * d1 + R0_->contents[aoffset + 1] * d2) +
          R0_->contents[aoffset + 2] * d3;
      d += W->contents[j + 3 * i] * upsilon[i];
    }
    b_W[j] = -d;
  }
  a = b_tau0_x_tau / (tau_dot_tau0->contents + 1.0);
  for (j = 0; j < 9; j++) {
    b_t_hat___tmp[j] = t_hat___tmp[j] + a * I3->contents[j];
  }
  for (j = 0; j < 3; j++) {
    d = G1->contents[j];
    d1 = G1->contents[j + 3];
    d2 = G1->contents[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_G1[j + 3 * aoffset] = (d * b_t_hat___tmp[3 * aoffset] +
                               d1 * b_t_hat___tmp[3 * aoffset + 1]) +
                              d2 * b_t_hat___tmp[3 * aoffset + 2];
    }
  }
  for (j = 0; j < 3; j++) {
    d = b_G1[j];
    d1 = b_G1[j + 3];
    d2 = b_G1[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_t_hat___tmp[j + 3 * aoffset] =
          (d * PP->contents[3 * aoffset] + d1 * PP->contents[3 * aoffset + 1]) +
          d2 * PP->contents[3 * aoffset + 2];
    }
  }
  for (j = 0; j < 9; j++) {
    b_t_hat___tmp[j] *= a_tmp;
  }
  for (j = 0; j < 3; j++) {
    d = y[j];
    d1 = y[j + 3];
    d2 = y[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_G1[aoffset + 3 * j] = b_W[aoffset] * tau->contents[j];
      b_y[j + 3 * aoffset] = (d * b_t_hat___tmp[3 * aoffset] +
                              d1 * b_t_hat___tmp[3 * aoffset + 1]) +
                             d2 * b_t_hat___tmp[3 * aoffset + 2];
    }
  }
  d = tau->contents[0];
  d1 = tau->contents[1];
  d2 = tau->contents[2];
  for (j = 0; j < 3; j++) {
    d3 = pie[j];
    b_t_hat___tmp[3 * j] = d * d3;
    b_t_hat___tmp[3 * j + 1] = d1 * d3;
    b_t_hat___tmp[3 * j + 2] = d2 * d3;
  }
  for (j = 0; j < 9; j++) {
    b_G1[j] = a_tmp * (b_G1[j] + b_y[j]);
  }
  a = 0.0;
  b_tau0_x_tau = 0.0;
  for (j = 0; j < 3; j++) {
    d = b_t_hat___tmp[j];
    d1 = b_t_hat___tmp[j + 3];
    d2 = b_t_hat___tmp[j + 6];
    d3 = 0.0;
    for (aoffset = 0; aoffset < 3; aoffset++) {
      mm[j + 3 * aoffset] =
          (d * b_G1[3 * aoffset] + d1 * b_G1[3 * aoffset + 1]) +
          d2 * b_G1[3 * aoffset + 2];
      d3 += pie[aoffset] * W->contents[aoffset + 3 * j];
    }
    d = upsilon[j];
    a += d3 * d;
    d1 = Theta->contents[j];
    d2 = Theta->contents[j + 3];
    d3 = Theta->contents[j + 6];
    b_a = 0.0;
    for (aoffset = 0; aoffset < 3; aoffset++) {
      c_a = (d1 * R0_->contents[3 * aoffset] +
             d2 * R0_->contents[3 * aoffset + 1]) +
            d3 * R0_->contents[3 * aoffset + 2];
      b_t_hat___tmp[j + 3 * aoffset] = c_a;
      b_a += c_a * pie[aoffset];
    }
    b_W[j] = b_a;
    b_tau0_x_tau += tau0_x_tau->contents[j] * d;
  }
  a /= nxip->contents;
  b_a = b_tau0_x_tau / (tau_dot_tau0->contents + 1.0);
  b_tau0_x_tau = 0.0;
  d = b_W[0];
  d1 = b_W[1];
  d2 = b_W[2];
  for (j = 0; j < 3; j++) {
    G1Tpie[j] = (G1->contents[3 * j] * d + G1->contents[3 * j + 1] * d1) +
                G1->contents[3 * j + 2] * d2;
    b_tau0_x_tau += tau0_x_tau->contents[j] * upsilon[j];
  }
  c_a = b_tau0_x_tau / (tau_dot_tau0->contents + 1.0);
  for (j = 0; j < 3; j++) {
    y[3 * j] = t_hat___tmp[j] + b_a * I3->contents[j];
    y[3 * j + 1] = t_hat___tmp[j + 3] + b_a * I3->contents[j + 3];
    y[3 * j + 2] = t_hat___tmp[j + 6] + b_a * I3->contents[j + 6];
  }
  a_tmp_tmp = -1.0 / nxip->contents;
  b_a_tmp = a_tmp_tmp / (tau_dot_tau0->contents + 1.0);
  b_tau0_x_tau = 0.0;
  b_tau = 0.0;
  d = G1Tpie[0];
  d1 = G1Tpie[1];
  d2 = G1Tpie[2];
  for (j = 0; j < 3; j++) {
    d3 = (y[j] * d + y[j + 3] * d1) + y[j + 6] * d2;
    c_y[j] = d3;
    b_tau0_x_tau += tau0_x_tau->contents[j] * upsilon[j];
    b_tau += tau->contents[j] * d3;
  }
  b_a = b_tau0_x_tau / (tau_dot_tau0->contents + 1.0);
  for (j = 0; j < 9; j++) {
    b_t_hat___tmp[j] = t_hat___tmp[j] + c_a * I3->contents[j];
  }
  for (j = 0; j < 3; j++) {
    d = G1->contents[j];
    d1 = G1->contents[j + 3];
    d2 = G1->contents[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_G1[j + 3 * aoffset] = (d * b_t_hat___tmp[3 * aoffset] +
                               d1 * b_t_hat___tmp[3 * aoffset + 1]) +
                              d2 * b_t_hat___tmp[3 * aoffset + 2];
    }
  }
  for (j = 0; j < 3; j++) {
    d = PP->contents[3 * j];
    d1 = PP->contents[3 * j + 1];
    d2 = PP->contents[3 * j + 2];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_t_hat___tmp[j + 3 * aoffset] =
          (b_G1[aoffset] * d + b_G1[aoffset + 3] * d1) + b_G1[aoffset + 6] * d2;
    }
  }
  for (j = 0; j < 9; j++) {
    b_t_hat___tmp[j] = -(a_tmp * b_t_hat___tmp[j]);
  }
  d = b_W[0];
  d1 = b_W[1];
  d2 = b_W[2];
  for (j = 0; j < 3; j++) {
    d_a[j] = (b_t_hat___tmp[j] * d + b_t_hat___tmp[j + 3] * d1) +
             b_t_hat___tmp[j + 6] * d2;
  }
  d = tau->contents[0];
  d1 = tau->contents[1];
  d2 = tau->contents[2];
  for (j = 0; j < 3; j++) {
    d3 = c_y[j];
    b_t_hat___tmp[3 * j] = d * d3 + b_tau * I3->contents[3 * j];
    coffset = 3 * j + 1;
    b_t_hat___tmp[coffset] = d1 * d3 + b_tau * I3->contents[coffset];
    coffset = 3 * j + 2;
    b_t_hat___tmp[coffset] = d2 * d3 + b_tau * I3->contents[coffset];
  }
  for (j = 0; j < 9; j++) {
    b_t_hat___tmp[j] *= a_tmp_tmp;
  }
  for (j = 0; j < 3; j++) {
    d = b_t_hat___tmp[j];
    d1 = b_t_hat___tmp[j + 3];
    d2 = b_t_hat___tmp[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      coffset = aoffset + 3 * j;
      e_a[coffset] = d_a[aoffset] * tau->contents[j];
      f_a[j + 3 * aoffset] =
          (d * PP->contents[3 * aoffset] + d1 * PP->contents[3 * aoffset + 1]) +
          d2 * PP->contents[3 * aoffset + 2];
      b_G1[coffset] = b_a_tmp * G1Tpie[aoffset] * tau0->contents[j];
    }
  }
  for (j = 0; j < 9; j++) {
    t_hat___tmp[j] += b_a * I3->contents[j];
  }
  for (j = 0; j < 3; j++) {
    d = b_G1[j];
    d1 = b_G1[j + 3];
    d2 = b_G1[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_t_hat___tmp[j + 3 * aoffset] =
          (d * t_hat___tmp[3 * aoffset] + d1 * t_hat___tmp[3 * aoffset + 1]) +
          d2 * t_hat___tmp[3 * aoffset + 2];
    }
  }
  d = tau0->contents[0];
  d1 = tau0->contents[1];
  d2 = tau0->contents[2];
  for (j = 0; j < 3; j++) {
    d3 = b_W[j];
    b_G1[3 * j] = b_a_tmp * d * d3;
    b_G1[3 * j + 1] = b_a_tmp * d1 * d3;
    b_G1[3 * j + 2] = b_a_tmp * d2 * d3;
  }
  for (j = 0; j < 3; j++) {
    d = b_G1[j];
    d1 = b_G1[j + 3];
    d2 = b_G1[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      t_hat___tmp[j + 3 * aoffset] =
          (d * G1->contents[3 * aoffset] + d1 * G1->contents[3 * aoffset + 1]) +
          d2 * G1->contents[3 * aoffset + 2];
    }
    d = t_hat___tmp[j];
    d1 = t_hat___tmp[j + 3];
    d2 = t_hat___tmp[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_G1[j + 3 * aoffset] =
          (d * PP->contents[3 * aoffset] + d1 * PP->contents[3 * aoffset + 1]) +
          d2 * PP->contents[3 * aoffset + 2];
    }
    d = b_t_hat___tmp[j];
    d1 = b_t_hat___tmp[j + 3];
    d2 = b_t_hat___tmp[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      t_hat___tmp[j + 3 * aoffset] =
          (d * PP->contents[3 * aoffset] + d1 * PP->contents[3 * aoffset + 1]) +
          d2 * PP->contents[3 * aoffset + 2];
    }
  }
  for (j = 0; j < 3; j++) {
    d = y[j];
    d1 = y[j + 3];
    d2 = y[j + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_y[j + 3 * aoffset] =
          (d * b_G1[3 * aoffset] + d1 * b_G1[3 * aoffset + 1]) +
          d2 * b_G1[3 * aoffset + 2];
    }
  }
  for (j = 0; j < 9; j++) {
    t_hat___tmp[j] += b_y[j];
  }
  for (j = 0; j < 3; j++) {
    for (aoffset = 0; aoffset < 3; aoffset++) {
      coffset = j + 3 * aoffset;
      QQ[coffset] =
          ((-(mm[coffset] + mm[aoffset + 3 * j]) - a * PP->contents[coffset]) +
           ((e_a[coffset] + f_a[coffset]) +
            ((PP->contents[j] * t_hat___tmp[3 * aoffset] +
              PP->contents[j + 3] * t_hat___tmp[3 * aoffset + 1]) +
             PP->contents[j + 6] * t_hat___tmp[3 * aoffset + 2])) /
               nxip->contents) /
          nxip->contents;
    }
  }
}

static void Q1T_Wthetp(const captured_var *tau0_x_tau,
                       const b_captured_var *tau_dot_tau0,
                       const c_captured_var *I3, const b_captured_var *nxip,
                       const c_captured_var *G1, const c_captured_var *PP,
                       const c_captured_var *W, const captured_var *tau,
                       const c_captured_var *R0_, const c_captured_var *Theta,
                       const captured_var *e1, const double upsilon[3],
                       const double pie[3], double QQ[3])
{
  double b_G1[9];
  double c_G1[9];
  double c_W[9];
  double y[9];
  double b_W[3];
  double a_tmp;
  double b_tau0_x_tau;
  double d;
  double e1_idx_1;
  double e1_idx_2;
  int aoffset;
  int coffset;
  int i;
  int j;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  FUNCTIONS FOR DERIVATIVE COMPUTATIONS */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  b_tau0_x_tau = 0.0;
  a_tmp = 1.0 / nxip->contents;
  for (j = 0; j < 3; j++) {
    double d1;
    coffset = j * 3;
    b_tau0_x_tau += tau0_x_tau->contents[j] * upsilon[j];
    e1_idx_1 = 0.0;
    e1_idx_2 = Theta->contents[j];
    d = Theta->contents[j + 3];
    d1 = Theta->contents[j + 6];
    for (i = 0; i < 3; i++) {
      aoffset = i * 3;
      y[coffset + i] =
          (R0_->contents[aoffset] * e1_idx_2 + R0_->contents[aoffset + 1] * d) +
          R0_->contents[aoffset + 2] * d1;
      e1_idx_1 += W->contents[j + 3 * i] * upsilon[i];
    }
    b_W[j] = -e1_idx_1;
  }
  b_tau0_x_tau /= tau_dot_tau0->contents + 1.0;
  c_W[0] = 0.0;
  c_W[3] = -upsilon[2];
  c_W[6] = upsilon[1];
  c_W[1] = upsilon[2];
  c_W[4] = 0.0;
  c_W[7] = -upsilon[0];
  c_W[2] = -upsilon[1];
  c_W[5] = upsilon[0];
  c_W[8] = 0.0;
  for (coffset = 0; coffset < 9; coffset++) {
    c_W[coffset] += b_tau0_x_tau * I3->contents[coffset];
  }
  for (coffset = 0; coffset < 3; coffset++) {
    e1_idx_1 = G1->contents[coffset];
    e1_idx_2 = G1->contents[coffset + 3];
    d = G1->contents[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_G1[coffset + 3 * aoffset] =
          (e1_idx_1 * c_W[3 * aoffset] + e1_idx_2 * c_W[3 * aoffset + 1]) +
          d * c_W[3 * aoffset + 2];
    }
    e1_idx_1 = b_G1[coffset];
    e1_idx_2 = b_G1[coffset + 3];
    d = b_G1[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      c_G1[coffset + 3 * aoffset] = (e1_idx_1 * PP->contents[3 * aoffset] +
                                     e1_idx_2 * PP->contents[3 * aoffset + 1]) +
                                    d * PP->contents[3 * aoffset + 2];
    }
  }
  for (coffset = 0; coffset < 9; coffset++) {
    c_G1[coffset] *= a_tmp;
  }
  for (coffset = 0; coffset < 3; coffset++) {
    e1_idx_1 = y[coffset];
    e1_idx_2 = y[coffset + 3];
    d = y[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      c_W[aoffset + 3 * coffset] = b_W[aoffset] * tau->contents[coffset];
      b_G1[coffset + 3 * aoffset] =
          (e1_idx_1 * c_G1[3 * aoffset] + e1_idx_2 * c_G1[3 * aoffset + 1]) +
          d * c_G1[3 * aoffset + 2];
    }
  }
  b_tau0_x_tau = e1->contents[1] * pie[2] - pie[1] * e1->contents[2];
  e1_idx_1 = pie[0] * e1->contents[2] - e1->contents[0] * pie[2];
  e1_idx_2 = e1->contents[0] * pie[1] - pie[0] * e1->contents[1];
  for (coffset = 0; coffset < 9; coffset++) {
    c_W[coffset] = a_tmp * (c_W[coffset] + b_G1[coffset]);
  }
  for (coffset = 0; coffset < 3; coffset++) {
    QQ[coffset] =
        (b_tau0_x_tau * c_W[3 * coffset] + e1_idx_1 * c_W[3 * coffset + 1]) +
        e1_idx_2 * c_W[3 * coffset + 2];
  }
}

static void Q1_Wpp(const captured_var *tau0_x_tau,
                   const b_captured_var *tau_dot_tau0, const c_captured_var *I3,
                   const b_captured_var *nxip, const c_captured_var *G1,
                   const c_captured_var *PP, const c_captured_var *W,
                   const captured_var *tau, const c_captured_var *R0_,
                   const c_captured_var *Theta, const captured_var *tau0,
                   const double upsilon[3], const double pie[3], double QQ[9])
{
  double G2[9];
  double b_G1[9];
  double b_W[9];
  double b_a[9];
  double b_t_hat___tmp[9];
  double b_y[9];
  double c_G1[9];
  double c_a[9];
  double d_a[9];
  double t_hat___tmp[9];
  double y[9];
  double Ppie[3];
  double b_G2[3];
  double a;
  double a_tmp;
  double b_a_tmp;
  double b_tau;
  double b_tau0_x_tau;
  double c_a_tmp;
  double c_tau;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double d_a_tmp;
  double s_tmp;
  int aoffset_tmp;
  int b_i;
  int coffset_tmp;
  int i;
  int j;
  b_tau = (tau->contents[0] * pie[0] + tau->contents[1] * pie[1]) +
          tau->contents[2] * pie[2];
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  FUNCTIONS FOR DERIVATIVE COMPUTATIONS */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  b_tau0_x_tau = 0.0;
  t_hat___tmp[0] = 0.0;
  t_hat___tmp[3] = -upsilon[2];
  t_hat___tmp[6] = upsilon[1];
  t_hat___tmp[1] = upsilon[2];
  t_hat___tmp[4] = 0.0;
  t_hat___tmp[7] = -upsilon[0];
  t_hat___tmp[2] = -upsilon[1];
  t_hat___tmp[5] = upsilon[0];
  t_hat___tmp[8] = 0.0;
  a_tmp = 1.0 / nxip->contents;
  for (j = 0; j < 3; j++) {
    coffset_tmp = j * 3;
    b_tau0_x_tau += tau0_x_tau->contents[j] * upsilon[j];
    d = Theta->contents[j];
    d1 = Theta->contents[j + 3];
    d2 = Theta->contents[j + 6];
    for (i = 0; i < 3; i++) {
      aoffset_tmp = i * 3;
      s_tmp = (R0_->contents[aoffset_tmp] * d +
               R0_->contents[aoffset_tmp + 1] * d1) +
              R0_->contents[aoffset_tmp + 2] * d2;
      aoffset_tmp = coffset_tmp + i;
      y[aoffset_tmp] = s_tmp;
      b_y[aoffset_tmp] = s_tmp;
    }
  }
  b_a_tmp = b_tau0_x_tau / (tau_dot_tau0->contents + 1.0);
  for (i = 0; i < 9; i++) {
    G2[i] = t_hat___tmp[i] + b_a_tmp * I3->contents[i];
  }
  b_tau0_x_tau = 0.0;
  c_a_tmp = -1.0 / nxip->contents;
  a = c_a_tmp / (tau_dot_tau0->contents + 1.0);
  c_tau = 0.0;
  for (i = 0; i < 3; i++) {
    Ppie[i] = (PP->contents[i] * pie[0] + PP->contents[i + 3] * pie[1]) +
              PP->contents[i + 6] * pie[2];
    b_tau0_x_tau += tau0_x_tau->contents[i] * upsilon[i];
    c_tau += tau->contents[i] * pie[i];
  }
  d_a_tmp = b_tau0_x_tau / (tau_dot_tau0->contents + 1.0);
  s_tmp = 0.0;
  d = Ppie[0];
  d1 = Ppie[1];
  d2 = Ppie[2];
  d3 = upsilon[0];
  d4 = upsilon[1];
  d5 = upsilon[2];
  for (i = 0; i < 3; i++) {
    s_tmp +=
        tau0->contents[i] * ((G2[i] * d + G2[i + 3] * d1) + G2[i + 6] * d2);
    b_G2[i] = -((W->contents[i] * d3 + W->contents[i + 3] * d4) +
                W->contents[i + 6] * d5);
  }
  b_tau0_x_tau = -s_tmp / nxip->contents / (tau_dot_tau0->contents + 1.0);
  for (i = 0; i < 9; i++) {
    b_t_hat___tmp[i] = t_hat___tmp[i] + b_a_tmp * I3->contents[i];
  }
  for (i = 0; i < 3; i++) {
    d = G1->contents[i];
    d1 = G1->contents[i + 3];
    d2 = G1->contents[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      b_G1[i + 3 * b_i] =
          (d * b_t_hat___tmp[3 * b_i] + d1 * b_t_hat___tmp[3 * b_i + 1]) +
          d2 * b_t_hat___tmp[3 * b_i + 2];
    }
    d = b_G1[i];
    d1 = b_G1[i + 3];
    d2 = b_G1[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      c_G1[i + 3 * b_i] =
          (d * PP->contents[3 * b_i] + d1 * PP->contents[3 * b_i + 1]) +
          d2 * PP->contents[3 * b_i + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    c_G1[i] *= a_tmp;
  }
  for (i = 0; i < 3; i++) {
    d = y[i];
    d1 = y[i + 3];
    d2 = y[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      b_W[b_i + 3 * i] = b_G2[b_i] * tau->contents[i];
      b_t_hat___tmp[i + 3 * b_i] =
          (d * c_G1[3 * b_i] + d1 * c_G1[3 * b_i + 1]) + d2 * c_G1[3 * b_i + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    b_W[i] = -(a_tmp * (b_W[i] + b_t_hat___tmp[i]));
  }
  d = pie[0];
  d1 = pie[1];
  d2 = pie[2];
  d3 = upsilon[0];
  d4 = upsilon[1];
  d5 = upsilon[2];
  for (i = 0; i < 3; i++) {
    s_tmp = tau->contents[i];
    b_t_hat___tmp[3 * i] = d * s_tmp + b_tau * I3->contents[3 * i];
    aoffset_tmp = 3 * i + 1;
    b_t_hat___tmp[aoffset_tmp] = d1 * s_tmp + b_tau * I3->contents[aoffset_tmp];
    aoffset_tmp = 3 * i + 2;
    b_t_hat___tmp[aoffset_tmp] = d2 * s_tmp + b_tau * I3->contents[aoffset_tmp];
    b_G2[i] = a_tmp * ((W->contents[i] * d3 + W->contents[i + 3] * d4) +
                       W->contents[i + 6] * d5);
  }
  for (i = 0; i < 3; i++) {
    d = b_W[i];
    d1 = b_W[i + 3];
    d2 = b_W[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      c_G1[b_i + 3 * i] = b_G2[b_i] * pie[i];
      b_a[i + 3 * b_i] =
          (d * b_t_hat___tmp[3 * b_i] + d1 * b_t_hat___tmp[3 * b_i + 1]) +
          d2 * b_t_hat___tmp[3 * b_i + 2];
    }
  }
  for (i = 0; i < 3; i++) {
    d = c_G1[i];
    d1 = c_G1[i + 3];
    d2 = c_G1[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      c_a[i + 3 * b_i] =
          (d * PP->contents[3 * b_i] + d1 * PP->contents[3 * b_i + 1]) +
          d2 * PP->contents[3 * b_i + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    b_t_hat___tmp[i] = t_hat___tmp[i] + d_a_tmp * I3->contents[i];
  }
  for (i = 0; i < 3; i++) {
    d = G1->contents[i];
    d1 = G1->contents[i + 3];
    d2 = G1->contents[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      b_G1[i + 3 * b_i] =
          (d * b_t_hat___tmp[3 * b_i] + d1 * b_t_hat___tmp[3 * b_i + 1]) +
          d2 * b_t_hat___tmp[3 * b_i + 2];
    }
    d = b_G1[i];
    d1 = b_G1[i + 3];
    d2 = b_G1[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      c_G1[i + 3 * b_i] =
          (d * PP->contents[3 * b_i] + d1 * PP->contents[3 * b_i + 1]) +
          d2 * PP->contents[3 * b_i + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    c_G1[i] = -(a_tmp * c_G1[i]);
  }
  d = pie[0];
  d1 = pie[1];
  d2 = pie[2];
  for (i = 0; i < 3; i++) {
    b_G2[i] = (c_G1[i] * d + c_G1[i + 3] * d1) + c_G1[i + 6] * d2;
  }
  for (i = 0; i < 3; i++) {
    d = G1->contents[i];
    d1 = G1->contents[i + 3];
    d2 = G1->contents[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      aoffset_tmp = b_i + 3 * i;
      c_G1[aoffset_tmp] = b_G2[b_i] * tau->contents[i];
      y[i + 3 * b_i] = (b_tau0_x_tau * d * PP->contents[3 * b_i] +
                        b_tau0_x_tau * d1 * PP->contents[3 * b_i + 1]) +
                       b_tau0_x_tau * d2 * PP->contents[3 * b_i + 2];
      b_W[aoffset_tmp] = a * Ppie[b_i] * tau0->contents[i];
    }
  }
  for (i = 0; i < 9; i++) {
    t_hat___tmp[i] += d_a_tmp * I3->contents[i];
  }
  for (i = 0; i < 3; i++) {
    d = b_W[i];
    d1 = b_W[i + 3];
    d2 = b_W[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      d_a[i + 3 * b_i] =
          (d * t_hat___tmp[3 * b_i] + d1 * t_hat___tmp[3 * b_i + 1]) +
          d2 * t_hat___tmp[3 * b_i + 2];
    }
    d = d_a[i];
    d1 = d_a[i + 3];
    d2 = d_a[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      b_W[i + 3 * b_i] =
          (d * PP->contents[3 * b_i] + d1 * PP->contents[3 * b_i + 1]) +
          d2 * PP->contents[3 * b_i + 2];
    }
  }
  d = tau->contents[0];
  d1 = tau->contents[1];
  d2 = tau->contents[2];
  for (i = 0; i < 3; i++) {
    d3 = pie[i];
    b_t_hat___tmp[3 * i] = d * d3 + c_tau * I3->contents[3 * i];
    aoffset_tmp = 3 * i + 1;
    b_t_hat___tmp[aoffset_tmp] = d1 * d3 + c_tau * I3->contents[aoffset_tmp];
    aoffset_tmp = 3 * i + 2;
    b_t_hat___tmp[aoffset_tmp] = d2 * d3 + c_tau * I3->contents[aoffset_tmp];
  }
  for (i = 0; i < 9; i++) {
    b_t_hat___tmp[i] *= c_a_tmp;
  }
  for (i = 0; i < 3; i++) {
    d = b_t_hat___tmp[i];
    d1 = b_t_hat___tmp[i + 3];
    d2 = b_t_hat___tmp[i + 6];
    d3 = G1->contents[i];
    d4 = G1->contents[i + 3];
    d5 = G1->contents[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      aoffset_tmp = 3 * b_i + 1;
      coffset_tmp = 3 * b_i + 2;
      j = i + 3 * b_i;
      d_a[j] = (d * PP->contents[3 * b_i] + d1 * PP->contents[aoffset_tmp]) +
               d2 * PP->contents[coffset_tmp];
      b_G1[j] =
          (d3 * G2[3 * b_i] + d4 * G2[aoffset_tmp]) + d5 * G2[coffset_tmp];
    }
    for (b_i = 0; b_i < 3; b_i++) {
      aoffset_tmp = i + 3 * b_i;
      b_t_hat___tmp[aoffset_tmp] =
          (c_G1[aoffset_tmp] + y[aoffset_tmp]) +
          ((d3 * b_W[3 * b_i] + d4 * b_W[3 * b_i + 1]) + d5 * b_W[3 * b_i + 2]);
    }
  }
  for (i = 0; i < 3; i++) {
    d = b_G1[i];
    d1 = b_G1[i + 3];
    d2 = b_G1[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      c_G1[i + 3 * b_i] =
          (d * d_a[3 * b_i] + d1 * d_a[3 * b_i + 1]) + d2 * d_a[3 * b_i + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    b_t_hat___tmp[i] = (b_t_hat___tmp[i] + c_G1[i]) / nxip->contents;
  }
  for (i = 0; i < 3; i++) {
    d = b_y[i];
    d1 = b_y[i + 3];
    d2 = b_y[i + 6];
    for (b_i = 0; b_i < 3; b_i++) {
      aoffset_tmp = i + 3 * b_i;
      QQ[aoffset_tmp] =
          ((b_a[aoffset_tmp] - c_a[aoffset_tmp]) +
           ((d * b_t_hat___tmp[3 * b_i] + d1 * b_t_hat___tmp[3 * b_i + 1]) +
            d2 * b_t_hat___tmp[3 * b_i + 2])) /
          nxip->contents;
    }
  }
}

static void Q1_Wpthet(const captured_var *e1, const c_captured_var *W,
                      const captured_var *tau0_x_tau,
                      const b_captured_var *tau_dot_tau0,
                      const c_captured_var *I3, const b_captured_var *nxip,
                      const c_captured_var *G1, const c_captured_var *PP,
                      const captured_var *tau, const c_captured_var *R0_,
                      const c_captured_var *Theta, const double upsilon[3],
                      const double pie[3], double QQ[3])
{
  double b_G1[9];
  double b_y[9];
  double c_G1[9];
  double c_y[9];
  double d_y[3];
  double y[3];
  double a;
  double b_tau;
  double b_tau0_x_tau;
  double d;
  double d1;
  double d2;
  int aoffset;
  int coffset;
  int i;
  int j;
  b_tau = (tau->contents[0] * pie[0] + tau->contents[1] * pie[1]) +
          tau->contents[2] * pie[2];
  b_tau0_x_tau = 0.0;
  for (j = 0; j < 3; j++) {
    coffset = j * 3;
    d = 0.0;
    d1 = Theta->contents[j];
    d2 = Theta->contents[j + 3];
    a = Theta->contents[j + 6];
    for (i = 0; i < 3; i++) {
      d += W->contents[j + 3 * i] * upsilon[i];
      aoffset = i * 3;
      b_y[coffset + i] =
          (R0_->contents[aoffset] * d1 + R0_->contents[aoffset + 1] * d2) +
          R0_->contents[aoffset + 2] * a;
    }
    y[j] = d;
    b_tau0_x_tau += tau0_x_tau->contents[j] * upsilon[j];
  }
  b_tau0_x_tau /= tau_dot_tau0->contents + 1.0;
  a = 1.0 / nxip->contents;
  c_y[0] = 0.0;
  c_y[3] = -upsilon[2];
  c_y[6] = upsilon[1];
  c_y[1] = upsilon[2];
  c_y[4] = 0.0;
  c_y[7] = -upsilon[0];
  c_y[2] = -upsilon[1];
  c_y[5] = upsilon[0];
  c_y[8] = 0.0;
  for (coffset = 0; coffset < 9; coffset++) {
    c_y[coffset] += b_tau0_x_tau * I3->contents[coffset];
  }
  for (coffset = 0; coffset < 3; coffset++) {
    d = G1->contents[coffset];
    d1 = G1->contents[coffset + 3];
    d2 = G1->contents[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_G1[coffset + 3 * aoffset] =
          (d * c_y[3 * aoffset] + d1 * c_y[3 * aoffset + 1]) +
          d2 * c_y[3 * aoffset + 2];
    }
    d = b_G1[coffset];
    d1 = b_G1[coffset + 3];
    d2 = b_G1[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      c_G1[coffset + 3 * aoffset] =
          (d * PP->contents[3 * aoffset] + d1 * PP->contents[3 * aoffset + 1]) +
          d2 * PP->contents[3 * aoffset + 2];
    }
  }
  for (coffset = 0; coffset < 9; coffset++) {
    c_G1[coffset] *= a;
  }
  for (coffset = 0; coffset < 3; coffset++) {
    d = 0.0;
    d1 = b_y[coffset];
    d2 = b_y[coffset + 3];
    a = b_y[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      d += ((d1 * c_G1[3 * aoffset] + d2 * c_G1[3 * aoffset + 1]) +
            a * c_G1[3 * aoffset + 2]) *
           pie[aoffset];
    }
    d_y[coffset] = d;
  }
  b_tau0_x_tau = -1.0 / nxip->contents;
  QQ[0] = b_tau0_x_tau *
          (b_tau * -(e1->contents[1] * y[2] - y[1] * e1->contents[2]) +
           (e1->contents[1] * d_y[2] - d_y[1] * e1->contents[2]));
  QQ[1] = b_tau0_x_tau *
          (b_tau * -(y[0] * e1->contents[2] - e1->contents[0] * y[2]) +
           (d_y[0] * e1->contents[2] - e1->contents[0] * d_y[2]));
  QQ[2] = b_tau0_x_tau *
          (b_tau * -(e1->contents[0] * y[1] - y[0] * e1->contents[1]) +
           (e1->contents[0] * d_y[1] - d_y[0] * e1->contents[1]));
}

static void Q1_Wthetp(const captured_var *tau0_x_tau,
                      const b_captured_var *tau_dot_tau0,
                      const c_captured_var *I3, const b_captured_var *nxip,
                      const c_captured_var *G1, const c_captured_var *PP,
                      const c_captured_var *W, const captured_var *tau,
                      const c_captured_var *R0_, const c_captured_var *Theta,
                      const captured_var *e1, const double upsilon[3],
                      double pie, double QQ[9])
{
  double b_G1[9];
  double c_G1[9];
  double c_W[9];
  double y[9];
  double b_W[3];
  double a_tmp;
  double b_tau0_x_tau;
  double d;
  double d1;
  double d2;
  int aoffset;
  int coffset;
  int i;
  int j;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  FUNCTIONS FOR DERIVATIVE COMPUTATIONS */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  b_tau0_x_tau = 0.0;
  a_tmp = 1.0 / nxip->contents;
  for (j = 0; j < 3; j++) {
    double d3;
    coffset = j * 3;
    b_tau0_x_tau += tau0_x_tau->contents[j] * upsilon[j];
    d = 0.0;
    d1 = Theta->contents[j];
    d2 = Theta->contents[j + 3];
    d3 = Theta->contents[j + 6];
    for (i = 0; i < 3; i++) {
      aoffset = i * 3;
      y[coffset + i] =
          (R0_->contents[aoffset] * d1 + R0_->contents[aoffset + 1] * d2) +
          R0_->contents[aoffset + 2] * d3;
      d += W->contents[j + 3 * i] * upsilon[i];
    }
    b_W[j] = -d;
  }
  b_tau0_x_tau /= tau_dot_tau0->contents + 1.0;
  b_G1[0] = 0.0;
  b_G1[3] = -upsilon[2];
  b_G1[6] = upsilon[1];
  b_G1[1] = upsilon[2];
  b_G1[4] = 0.0;
  b_G1[7] = -upsilon[0];
  b_G1[2] = -upsilon[1];
  b_G1[5] = upsilon[0];
  b_G1[8] = 0.0;
  for (coffset = 0; coffset < 9; coffset++) {
    b_G1[coffset] += b_tau0_x_tau * I3->contents[coffset];
  }
  for (coffset = 0; coffset < 3; coffset++) {
    d = G1->contents[coffset];
    d1 = G1->contents[coffset + 3];
    d2 = G1->contents[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      c_G1[coffset + 3 * aoffset] =
          (d * b_G1[3 * aoffset] + d1 * b_G1[3 * aoffset + 1]) +
          d2 * b_G1[3 * aoffset + 2];
    }
  }
  for (coffset = 0; coffset < 3; coffset++) {
    d = c_G1[coffset];
    d1 = c_G1[coffset + 3];
    d2 = c_G1[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      b_G1[coffset + 3 * aoffset] =
          (d * PP->contents[3 * aoffset] + d1 * PP->contents[3 * aoffset + 1]) +
          d2 * PP->contents[3 * aoffset + 2];
    }
  }
  for (coffset = 0; coffset < 9; coffset++) {
    b_G1[coffset] *= a_tmp;
  }
  for (coffset = 0; coffset < 3; coffset++) {
    d = y[coffset];
    d1 = y[coffset + 3];
    d2 = y[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      c_W[aoffset + 3 * coffset] = b_W[aoffset] * tau->contents[coffset];
      c_G1[coffset + 3 * aoffset] =
          (d * b_G1[3 * aoffset] + d1 * b_G1[3 * aoffset + 1]) +
          d2 * b_G1[3 * aoffset + 2];
    }
  }
  b_G1[0] = -0.0;
  b_G1[3] = e1->contents[2];
  b_G1[6] = -e1->contents[1];
  b_G1[1] = -e1->contents[2];
  b_G1[4] = -0.0;
  b_G1[7] = e1->contents[0];
  b_G1[2] = e1->contents[1];
  b_G1[5] = -e1->contents[0];
  b_G1[8] = -0.0;
  for (coffset = 0; coffset < 9; coffset++) {
    c_W[coffset] = a_tmp * (c_W[coffset] + c_G1[coffset]);
  }
  for (coffset = 0; coffset < 3; coffset++) {
    d = b_G1[coffset];
    d1 = b_G1[coffset + 3];
    d2 = b_G1[coffset + 6];
    for (aoffset = 0; aoffset < 3; aoffset++) {
      QQ[coffset + 3 * aoffset] =
          (d * c_W[3 * aoffset] + d1 * c_W[3 * aoffset + 1]) +
          d2 * c_W[3 * aoffset + 2];
    }
  }
  for (coffset = 0; coffset < 9; coffset++) {
    QQ[coffset] *= pie;
  }
}

static void Q2_Wp(const c_captured_var *PP, const c_captured_var *G1,
                  const b_captured_var *nxip, const captured_var *tau0_x_tau,
                  const b_captured_var *tau_dot_tau0, const captured_var *tau,
                  const c_captured_var *W, const c_captured_var *R0_,
                  const c_captured_var *Theta, const double pie[3],
                  double QQ[9])
{
  double b_G1[9];
  double b_Ppie[9];
  double dv[9];
  double y[9];
  double Ppie[3];
  double b_tau;
  double b_tau_dot_tau0;
  double d;
  double d1;
  double d2;
  double d3;
  int aoffset;
  int coffset;
  int i;
  int j;
  b_tau = (tau->contents[0] * pie[0] + tau->contents[1] * pie[1]) +
          tau->contents[2] * pie[2];
  for (j = 0; j < 3; j++) {
    coffset = j * 3;
    d = Theta->contents[j];
    d1 = Theta->contents[j + 3];
    d2 = Theta->contents[j + 6];
    d3 = 0.0;
    for (i = 0; i < 3; i++) {
      aoffset = i * 3;
      y[coffset + i] =
          (R0_->contents[aoffset] * d + R0_->contents[aoffset + 1] * d1) +
          R0_->contents[aoffset + 2] * d2;
      d3 += PP->contents[j + 3 * i] * pie[i];
    }
    Ppie[j] = d3;
  }
  b_tau_dot_tau0 = tau_dot_tau0->contents + 1.0;
  dv[0] = -0.0;
  dv[3] = Ppie[2];
  dv[6] = -Ppie[1];
  dv[1] = -Ppie[2];
  dv[4] = -0.0;
  dv[7] = Ppie[0];
  dv[2] = Ppie[1];
  dv[5] = -Ppie[0];
  dv[8] = -0.0;
  d = Ppie[0];
  d1 = Ppie[1];
  d2 = Ppie[2];
  for (aoffset = 0; aoffset < 3; aoffset++) {
    d3 = tau0_x_tau->contents[aoffset];
    b_Ppie[3 * aoffset] = d * d3 / b_tau_dot_tau0;
    b_Ppie[3 * aoffset + 1] = d1 * d3 / b_tau_dot_tau0;
    b_Ppie[3 * aoffset + 2] = d2 * d3 / b_tau_dot_tau0;
  }
  for (aoffset = 0; aoffset < 9; aoffset++) {
    b_G1[aoffset] = G1->contents[aoffset] / nxip->contents;
    dv[aoffset] += b_Ppie[aoffset];
  }
  for (aoffset = 0; aoffset < 3; aoffset++) {
    d = b_G1[aoffset];
    d1 = b_G1[aoffset + 3];
    d2 = b_G1[aoffset + 6];
    for (j = 0; j < 3; j++) {
      b_Ppie[aoffset + 3 * j] =
          (d * dv[3 * j] + d1 * dv[3 * j + 1]) + d2 * dv[3 * j + 2];
    }
  }
  for (aoffset = 0; aoffset < 3; aoffset++) {
    d = y[aoffset];
    d1 = y[aoffset + 3];
    d2 = y[aoffset + 6];
    for (j = 0; j < 3; j++) {
      coffset = aoffset + 3 * j;
      QQ[coffset] = (-b_tau * W->contents[coffset] +
                     ((d * b_Ppie[3 * j] + d1 * b_Ppie[3 * j + 1]) +
                      d2 * b_Ppie[3 * j + 2])) /
                    nxip->contents;
    }
  }
}

static void binary_expand_op_16(emxArray_real_T *in1,
                                const emxArray_int32_T *in2,
                                const emxArray_int32_T *in3,
                                const emxArray_uint32_T *in4,
                                const emxArray_real_T *in5)
{
  emxArray_real_T *b_in1;
  const double *in5_data;
  double *b_in1_data;
  double *in1_data;
  const int *in2_data;
  const int *in3_data;
  const unsigned int *in4_data;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0_tmp;
  int stride_1_0;
  int stride_1_1;
  in5_data = in5->data;
  in4_data = in4->data;
  in3_data = in3->data;
  in2_data = in2->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  if (in5->size[0] == 1) {
    loop_ub = in4->size[0];
  } else {
    loop_ub = in5->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  if (in5->size[1] == 1) {
    b_loop_ub = in4->size[0];
  } else {
    b_loop_ub = in5->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0_tmp = (in4->size[0] != 1);
  stride_1_0 = (in5->size[0] != 1);
  stride_1_1 = (in5->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[((int)in4_data[i1 * stride_0_0_tmp] +
                    in1->size[0] * ((int)in4_data[aux_0_1] - 1)) -
                   1] +
          in5_data[i1 * stride_1_0 + in5->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_0_tmp;
  }
  loop_ub = in2->size[0];
  b_loop_ub = in3->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      in1_data[in2_data[i1] + in1->size[0] * in3_data[i]] =
          b_in1_data[i1 + loop_ub * i];
    }
  }
  emxFree_real_T(&b_in1);
}

static void binary_expand_op_17(emxArray_real_T *in1,
                                const emxArray_int32_T *in2, int in3, int in4,
                                const emxArray_uint32_T *in5, int in6, int in7,
                                const emxArray_real_T *in8)
{
  emxArray_real_T *b_in1;
  const double *in8_data;
  double *b_in1_data;
  double *in1_data;
  const int *in2_data;
  const unsigned int *in5_data;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  in8_data = in8->data;
  in5_data = in5->data;
  in2_data = in2->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  if (in8->size[0] == 1) {
    loop_ub = in5->size[0];
  } else {
    loop_ub = in8->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  i1 = (in7 - in6) + 1;
  if (in8->size[1] == 1) {
    b_loop_ub = i1;
  } else {
    b_loop_ub = in8->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in5->size[0] != 1);
  stride_0_1 = (i1 != 1);
  stride_1_0 = (in8->size[0] != 1);
  stride_1_1 = (in8->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[((int)in5_data[i1 * stride_0_0] +
                    in1->size[0] * (in6 + aux_0_1)) -
                   1] +
          in8_data[i1 * stride_1_0 + in8->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  loop_ub = in2->size[0];
  b_loop_ub = in4 - in3;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      in1_data[in2_data[i1] + in1->size[0] * (in3 + i)] =
          b_in1_data[i1 + loop_ub * i];
    }
  }
  emxFree_real_T(&b_in1);
}

static void binary_expand_op_22(emxArray_real_T *in1, int in2, int in3,
                                const emxArray_int32_T *in4, int in5, int in6,
                                const emxArray_uint32_T *in7,
                                const emxArray_real_T *in8)
{
  emxArray_real_T *b_in1;
  const double *in8_data;
  double *b_in1_data;
  double *in1_data;
  const int *in4_data;
  const unsigned int *in7_data;
  int aux_0_1;
  int aux_1_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  in8_data = in8->data;
  in7_data = in7->data;
  in4_data = in4->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  i = (in6 - in5) + 1;
  if (in8->size[0] == 1) {
    loop_ub = i;
  } else {
    loop_ub = in8->size[0];
  }
  i1 = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  if (in8->size[1] == 1) {
    b_loop_ub = in7->size[0];
  } else {
    b_loop_ub = in8->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(b_in1, i1);
  b_in1_data = b_in1->data;
  stride_0_0 = (i != 1);
  stride_0_1 = (in7->size[0] != 1);
  stride_1_0 = (in8->size[0] != 1);
  stride_1_1 = (in8->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[(in5 + i1 * stride_0_0) +
                   in1->size[0] * ((int)in7_data[aux_0_1] - 1)] +
          in8_data[i1 * stride_1_0 + in8->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  loop_ub = in3 - in2;
  b_loop_ub = in4->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      in1_data[(in2 + i1) + in1->size[0] * in4_data[i]] =
          b_in1_data[i1 + loop_ub * i];
    }
  }
  emxFree_real_T(&b_in1);
}

void CableInertiaForce(
    const emxArray_real_T *P, const emxArray_real_T *P0,
    const emxArray_real_T *Pdot, const emxArray_real_T *Pddot,
    const emxArray_real_T *varTheta, const emxArray_real_T *varThetadot,
    const emxArray_real_T *varThetaddot, const emxArray_real_T *R0, double rho,
    const double II[9], const emxArray_real_T *wg, const emxArray_real_T *nel,
    const emxArray_real_T *colmat, const emxArray_real_T *colmat_brev, double d,
    double dbrev, double alph0, emxArray_real_T *F, emxArray_real_T *mu,
    emxArray_real_T *F_ij, emxArray_real_T *F_ib, emxArray_real_T *F_ijd,
    emxArray_real_T *F_ibd, emxArray_real_T *F_ijdd, emxArray_real_T *F_ibdd,
    emxArray_real_T *mu_aj, emxArray_real_T *mu_ab, emxArray_real_T *mu_ajd,
    emxArray_real_T *mu_abd, emxArray_real_T *mu_ajdd, emxArray_real_T *mu_abdd)
{
  b_captured_var nxip;
  b_captured_var tau_dot_tau0;
  c_captured_var G1;
  c_captured_var I3;
  c_captured_var PP;
  c_captured_var R0_;
  c_captured_var Theta;
  c_captured_var b_W;
  captured_var b_tau;
  captured_var e1;
  captured_var tau0;
  captured_var tau0_x_tau;
  emxArray_int32_T *r;
  emxArray_int32_T *r1;
  emxArray_int32_T *r2;
  emxArray_real_T *B;
  emxArray_real_T *Bbrev;
  emxArray_real_T *Bp;
  emxArray_real_T *F_;
  emxArray_real_T *F_ib_;
  emxArray_real_T *F_ibd_;
  emxArray_real_T *F_ibdd_;
  emxArray_real_T *F_ij_;
  emxArray_real_T *F_ijd_;
  emxArray_real_T *F_ijdd_;
  emxArray_real_T *b_F;
  emxArray_real_T *b_F_ij;
  emxArray_real_T *b_mu;
  emxArray_real_T *mu_ab_;
  emxArray_real_T *mu_abd_;
  emxArray_real_T *mu_abdd_;
  emxArray_real_T *mu_aj_;
  emxArray_real_T *mu_ajd_;
  emxArray_real_T *mu_ajdd_;
  emxArray_real_T *y;
  emxArray_uint32_T *r3;
  double dv1[9];
  double dv2[9];
  double htau0[9];
  double L_Wthet_xipd[3];
  double dv[3];
  const double *R0_data;
  const double *colmat_brev_data;
  const double *colmat_data;
  const double *nel_data;
  const double *varTheta_data;
  const double *varThetaddot_data;
  const double *varThetadot_data;
  const double *wg_data;
  double M2_contents;
  double a;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double *B_data;
  double *Bbrev_data;
  double *Bp_data;
  double *F__data;
  double *F_data;
  double *F_ib__data;
  double *F_ib_data;
  double *F_ibd__data;
  double *F_ibd_data;
  double *F_ibdd__data;
  double *F_ibdd_data;
  double *F_ij__data;
  double *F_ij_data;
  double *F_ijd__data;
  double *F_ijd_data;
  double *F_ijdd__data;
  double *F_ijdd_data;
  double *b_F_data;
  double *b_mu_data;
  double *mu_ab__data;
  double *mu_ab_data;
  double *mu_abd__data;
  double *mu_abd_data;
  double *mu_abdd__data;
  double *mu_abdd_data;
  double *mu_aj__data;
  double *mu_aj_data;
  double *mu_ajd__data;
  double *mu_ajd_data;
  double *mu_ajdd__data;
  double *mu_ajdd_data;
  double *mu_data;
  double *y_data;
  int b_i;
  int b_loop_ub;
  int b_loop_ub_tmp;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int loop_ub;
  int loop_ub_tmp;
  int n;
  int *r4;
  unsigned int *r5;
  int *r6;
  colmat_brev_data = colmat_brev->data;
  colmat_data = colmat->data;
  nel_data = nel->data;
  wg_data = wg->data;
  R0_data = R0->data;
  varThetaddot_data = varThetaddot->data;
  varThetadot_data = varThetadot->data;
  varTheta_data = varTheta->data;
  /*  INPUTS */
  /*  P = spline position DOF (represented as 3*N matrix), \mathscr{P} */
  /*  P0 = reference configuration */
  /*  Pdot = \dot{P} */
  /*  Pddot = \ddot{P} */
  /*  varTheta = spline twist DOF (\vartheta) */
  /*  varThetadot = \dot{vartheta} */
  /*  varThetaddot = \ddot{vartheta} */
  /*  R0 = cell array of orientations in the reference configuration at the */
  /*       quadrature points (+ the two end points) */
  /*  rho = mass per unit length */
  /*  II = body frame mass moment of inertia (3x3 matrix) per unit length */
  /*  sg = quadrature points */
  /*  wg = quadrature weights */
  /*  nel = nel(n) = index of element to which quadrature pt sg(n) belongs to */
  /*  colmat = B-spline basis for position and its derivatives (B) */
  /*  colmat_brev = B-spline basis for twist and its derivative (Brev) */
  /*  d = degree of B-spline basis for position (B) */
  /*  dbrev = degree of B-spline basis for twist (Bbrev) */
  /*  alph0 = mass-proportional damping coefficient (JB Model) */
  /*  OUTPUTS */
  /*  F = Vector of inertia forces length 3*N */
  /*  mu = Vector of inertia moments length Nbrev */
  /*  F_ij, F_ib, F_ijd, F_ibd, F_ijdd, F_ibdd, ... */
  /*  mu_aj, mu_ab, mu_ajd, mu_abd, mu_ajdd, mu_abdd = derivatives */
  /*  number of basis functions (or control points) */
  /*  number of quadrature (or collocation) points  */
  e1.contents[0] = 1.0;
  e1.contents[1] = 0.0;
  e1.contents[2] = 0.0;
  for (i = 0; i < 9; i++) {
    I3.contents[i] = iv[i];
  }
  i = F->size[0];
  F->size[0] = 3 * P->size[1];
  emxEnsureCapacity_real_T(F, i);
  F_data = F->data;
  loop_ub = 3 * P->size[1];
  for (i = 0; i < loop_ub; i++) {
    F_data[i] = 0.0;
  }
  /*  Generalized centrifugal+coriolis force */
  emxInit_real_T(&F_, 1);
  loop_ub = (int)(3.0 * (d + 1.0));
  i = F_->size[0];
  F_->size[0] = loop_ub;
  emxEnsureCapacity_real_T(F_, i);
  F__data = F_->data;
  for (i = 0; i < loop_ub; i++) {
    F__data[i] = 0.0;
  }
  /*  temp storage when computing F */
  i = mu->size[0];
  mu->size[0] = varTheta->size[0];
  emxEnsureCapacity_real_T(mu, i);
  mu_data = mu->data;
  b_loop_ub = varTheta->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    mu_data[i] = 0.0;
  }
  /*  Generalized centrifugal+coriolis moment */
  i = F_ij->size[0] * F_ij->size[1];
  F_ij->size[0] = 3 * colmat->size[1];
  F_ij->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(F_ij, i);
  F_ij_data = F_ij->data;
  b_loop_ub = 3 * colmat->size[1] * (3 * colmat->size[1]);
  for (i = 0; i < b_loop_ub; i++) {
    F_ij_data[i] = 0.0;
  }
  emxInit_real_T(&F_ij_, 2);
  i = F_ij_->size[0] * F_ij_->size[1];
  F_ij_->size[0] = loop_ub;
  F_ij_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(F_ij_, i);
  F_ij__data = F_ij_->data;
  b_i = loop_ub * loop_ub;
  for (i = 0; i < b_i; i++) {
    F_ij__data[i] = 0.0;
  }
  /*  temp storage */
  i = F_ib->size[0] * F_ib->size[1];
  F_ib->size[0] = 3 * colmat->size[1];
  F_ib->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(F_ib, i);
  F_ib_data = F_ib->data;
  loop_ub_tmp = 3 * colmat->size[1] * colmat_brev->size[1];
  for (i = 0; i < loop_ub_tmp; i++) {
    F_ib_data[i] = 0.0;
  }
  emxInit_real_T(&F_ib_, 2);
  i = F_ib_->size[0] * F_ib_->size[1];
  F_ib_->size[0] = loop_ub;
  i1 = (int)(dbrev + 1.0);
  F_ib_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(F_ib_, i);
  F_ib__data = F_ib_->data;
  b_loop_ub_tmp = loop_ub * (int)(dbrev + 1.0);
  for (i = 0; i < b_loop_ub_tmp; i++) {
    F_ib__data[i] = 0.0;
  }
  /*  temp storage */
  i = F_ijd->size[0] * F_ijd->size[1];
  F_ijd->size[0] = 3 * colmat->size[1];
  F_ijd->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(F_ijd, i);
  F_ijd_data = F_ijd->data;
  for (i = 0; i < b_loop_ub; i++) {
    F_ijd_data[i] = 0.0;
  }
  emxInit_real_T(&F_ijd_, 2);
  i = F_ijd_->size[0] * F_ijd_->size[1];
  F_ijd_->size[0] = loop_ub;
  F_ijd_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(F_ijd_, i);
  F_ijd__data = F_ijd_->data;
  for (i = 0; i < b_i; i++) {
    F_ijd__data[i] = 0.0;
  }
  /*  temp storage */
  i = F_ibd->size[0] * F_ibd->size[1];
  F_ibd->size[0] = 3 * colmat->size[1];
  F_ibd->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(F_ibd, i);
  F_ibd_data = F_ibd->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    F_ibd_data[i] = 0.0;
  }
  emxInit_real_T(&F_ibd_, 2);
  i = F_ibd_->size[0] * F_ibd_->size[1];
  F_ibd_->size[0] = loop_ub;
  F_ibd_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(F_ibd_, i);
  F_ibd__data = F_ibd_->data;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    F_ibd__data[i] = 0.0;
  }
  /*  temp storage */
  i = F_ijdd->size[0] * F_ijdd->size[1];
  F_ijdd->size[0] = 3 * colmat->size[1];
  F_ijdd->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(F_ijdd, i);
  F_ijdd_data = F_ijdd->data;
  for (i = 0; i < b_loop_ub; i++) {
    F_ijdd_data[i] = 0.0;
  }
  emxInit_real_T(&F_ijdd_, 2);
  i = F_ijdd_->size[0] * F_ijdd_->size[1];
  F_ijdd_->size[0] = loop_ub;
  F_ijdd_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(F_ijdd_, i);
  F_ijdd__data = F_ijdd_->data;
  for (i = 0; i < b_i; i++) {
    F_ijdd__data[i] = 0.0;
  }
  /*  temp storage */
  i = F_ibdd->size[0] * F_ibdd->size[1];
  F_ibdd->size[0] = 3 * colmat->size[1];
  F_ibdd->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(F_ibdd, i);
  F_ibdd_data = F_ibdd->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    F_ibdd_data[i] = 0.0;
  }
  emxInit_real_T(&F_ibdd_, 2);
  i = F_ibdd_->size[0] * F_ibdd_->size[1];
  F_ibdd_->size[0] = loop_ub;
  F_ibdd_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(F_ibdd_, i);
  F_ibdd__data = F_ibdd_->data;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    F_ibdd__data[i] = 0.0;
  }
  /*  temp storage */
  i = mu_aj->size[0] * mu_aj->size[1];
  mu_aj->size[0] = colmat_brev->size[1];
  mu_aj->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(mu_aj, i);
  mu_aj_data = mu_aj->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    mu_aj_data[i] = 0.0;
  }
  emxInit_real_T(&mu_aj_, 2);
  i = mu_aj_->size[0] * mu_aj_->size[1];
  mu_aj_->size[0] = (int)(dbrev + 1.0);
  mu_aj_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(mu_aj_, i);
  mu_aj__data = mu_aj_->data;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    mu_aj__data[i] = 0.0;
  }
  /*  temp storage */
  i = mu_ab->size[0] * mu_ab->size[1];
  mu_ab->size[0] = colmat_brev->size[1];
  mu_ab->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(mu_ab, i);
  mu_ab_data = mu_ab->data;
  b_loop_ub = colmat_brev->size[1] * colmat_brev->size[1];
  for (i = 0; i < b_loop_ub; i++) {
    mu_ab_data[i] = 0.0;
  }
  emxInit_real_T(&mu_ab_, 2);
  i = mu_ab_->size[0] * mu_ab_->size[1];
  mu_ab_->size[0] = (int)(dbrev + 1.0);
  mu_ab_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(mu_ab_, i);
  mu_ab__data = mu_ab_->data;
  b_i = (int)(dbrev + 1.0) * (int)(dbrev + 1.0);
  for (i = 0; i < b_i; i++) {
    mu_ab__data[i] = 0.0;
  }
  /*  temp storage */
  i = mu_ajd->size[0] * mu_ajd->size[1];
  mu_ajd->size[0] = colmat_brev->size[1];
  mu_ajd->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(mu_ajd, i);
  mu_ajd_data = mu_ajd->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    mu_ajd_data[i] = 0.0;
  }
  emxInit_real_T(&mu_ajd_, 2);
  i = mu_ajd_->size[0] * mu_ajd_->size[1];
  mu_ajd_->size[0] = (int)(dbrev + 1.0);
  mu_ajd_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(mu_ajd_, i);
  mu_ajd__data = mu_ajd_->data;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    mu_ajd__data[i] = 0.0;
  }
  /*  temp storage */
  i = mu_abd->size[0] * mu_abd->size[1];
  mu_abd->size[0] = colmat_brev->size[1];
  mu_abd->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(mu_abd, i);
  mu_abd_data = mu_abd->data;
  for (i = 0; i < b_loop_ub; i++) {
    mu_abd_data[i] = 0.0;
  }
  emxInit_real_T(&mu_abd_, 2);
  i = mu_abd_->size[0] * mu_abd_->size[1];
  mu_abd_->size[0] = (int)(dbrev + 1.0);
  mu_abd_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(mu_abd_, i);
  mu_abd__data = mu_abd_->data;
  for (i = 0; i < b_i; i++) {
    mu_abd__data[i] = 0.0;
  }
  /*  temp storage */
  i = mu_ajdd->size[0] * mu_ajdd->size[1];
  mu_ajdd->size[0] = colmat_brev->size[1];
  mu_ajdd->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(mu_ajdd, i);
  mu_ajdd_data = mu_ajdd->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    mu_ajdd_data[i] = 0.0;
  }
  emxInit_real_T(&mu_ajdd_, 2);
  i = mu_ajdd_->size[0] * mu_ajdd_->size[1];
  mu_ajdd_->size[0] = (int)(dbrev + 1.0);
  mu_ajdd_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(mu_ajdd_, i);
  mu_ajdd__data = mu_ajdd_->data;
  for (i = 0; i < b_loop_ub_tmp; i++) {
    mu_ajdd__data[i] = 0.0;
  }
  /*  temp storage */
  i = mu_abdd->size[0] * mu_abdd->size[1];
  mu_abdd->size[0] = colmat_brev->size[1];
  mu_abdd->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(mu_abdd, i);
  mu_abdd_data = mu_abdd->data;
  for (i = 0; i < b_loop_ub; i++) {
    mu_abdd_data[i] = 0.0;
  }
  emxInit_real_T(&mu_abdd_, 2);
  i = mu_abdd_->size[0] * mu_abdd_->size[1];
  mu_abdd_->size[0] = (int)(dbrev + 1.0);
  mu_abdd_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(mu_abdd_, i);
  mu_abdd__data = mu_abdd_->data;
  for (i = 0; i < b_i; i++) {
    mu_abdd__data[i] = 0.0;
  }
  /*  temp storage */
  i = wg->size[0];
  emxInit_real_T(&B, 1);
  emxInit_real_T(&Bp, 1);
  emxInit_real_T(&Bbrev, 1);
  emxInit_int32_T(&r, 2);
  emxInit_int32_T(&r1, 1);
  emxInit_int32_T(&r2, 1);
  emxInit_real_T(&y, 2);
  emxInit_uint32_T(&r3);
  emxInit_real_T(&b_F, 2);
  emxInit_real_T(&b_mu, 2);
  emxInit_real_T(&b_F_ij, 2);
  if (i - 1 >= 0) {
    L_Wthet_xipd[0] = -0.0;
    a = alph0 * rho;
    i3 = (int)(d + 1.0);
    M2_contents = 0.0;
    d1 = II[0];
    d2 = II[1];
    d3 = II[2];
    d4 = II[0];
    d5 = II[1];
    d6 = II[2];
    d7 = II[0];
    d8 = II[3];
    d9 = II[6];
  }
  for (n = 0; n < i; n++) {
    double LL[9];
    double L_WTp_l[9];
    double L_Wp_xipd[9];
    double L_lp[9];
    double L_lpd[9];
    double M1_contents[9];
    double M3_contents[9];
    double Q2_WTp_xipd[9];
    double Q2_WTthet_thetd[9];
    double Q2_Wp_xipd[9];
    double W[9];
    double Wd_xipd_p[9];
    double b[9];
    double b_I3[9];
    double b_M1_contents[9];
    double b_a[9];
    double b_htau0[9];
    double b_y_tmp[9];
    double c_y[9];
    double e_a[9];
    double tau[9];
    double xi0p[9];
    double F1[3];
    double IWd_xipd[3];
    double Ialph[3];
    double L_WTthet_l[3];
    double L_lthet[3];
    double Q1[3];
    double Wd_xipd_thet[3];
    double b_e1[3];
    double b_xi0p[3];
    double b_y[3];
    double c_W[3];
    double c_e1[3];
    double d_W[3];
    double d_y[3];
    double f_a[3];
    double l[3];
    double xipd[3];
    double xipdd[3];
    double M2_contents_tmp;
    double a_tmp;
    double b_M2_contents;
    double b_Q1;
    double b_d;
    double b_varThetaddot;
    double b_varThetadot;
    double c_a;
    double c_tau;
    double c_xi0p;
    double c_y_tmp;
    double d10;
    double d11;
    double d12;
    double d13;
    double d14;
    double d15;
    double d16;
    double d17;
    double d18;
    double d19;
    double d20;
    double d21;
    double d_a;
    double d_tau;
    double d_xi0p;
    double nel_;
    double nxi0p;
    double y_tmp;
    nel_ = nel_data[n];
    b_d = 3.0 * (((double)n + 1.0) - 1.0);
    loop_ub = colmat->size[1];
    i2 = B->size[0];
    B->size[0] = colmat->size[1];
    emxEnsureCapacity_real_T(B, i2);
    B_data = B->data;
    i2 = Bp->size[0];
    Bp->size[0] = colmat->size[1];
    emxEnsureCapacity_real_T(Bp, i2);
    Bp_data = Bp->data;
    for (i2 = 0; i2 < loop_ub; i2++) {
      B_data[i2] = colmat_data[((int)(b_d + 1.0) + colmat->size[0] * i2) - 1];
      Bp_data[i2] = colmat_data[((int)(b_d + 2.0) + colmat->size[0] * i2) - 1];
    }
    i2 = (int)((unsigned int)n << 1);
    loop_ub = colmat_brev->size[1];
    i4 = Bbrev->size[0];
    Bbrev->size[0] = colmat_brev->size[1];
    emxEnsureCapacity_real_T(Bbrev, i4);
    Bbrev_data = Bbrev->data;
    for (i4 = 0; i4 < loop_ub; i4++) {
      Bbrev_data[i4] = colmat_brev_data[i2 + colmat_brev->size[0] * i4];
    }
    mtimes(P0, Bp, b_xi0p);
    nxi0p = b_norm(b_xi0p);
    tau0.contents[0] = b_xi0p[0] / nxi0p;
    tau0.contents[1] = b_xi0p[1] / nxi0p;
    tau0.contents[2] = b_xi0p[2] / nxi0p;
    mtimes(P, Bp, b_xi0p);
    mtimes(Pdot, Bp, xipd);
    mtimes(Pddot, Bp, xipdd);
    nxip.contents = b_norm(b_xi0p);
    b_d = b_xi0p[0] / nxip.contents;
    dv[0] = b_d;
    b_tau.contents[0] = b_d;
    d10 = b_xi0p[1] / nxip.contents;
    dv[1] = d10;
    b_tau.contents[1] = d10;
    d11 = b_xi0p[2] / nxip.contents;
    dv[2] = d11;
    b_tau.contents[2] = d11;
    b_Q1 = 0.0;
    loop_ub = varTheta->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      b_Q1 += varTheta_data[i2] * Bbrev_data[i2];
    }
    b_varThetadot = 0.0;
    loop_ub = varThetadot->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      b_varThetadot += varThetadot_data[i2] * Bbrev_data[i2];
    }
    b_varThetaddot = 0.0;
    loop_ub = varThetaddot->size[0];
    for (i2 = 0; i2 < loop_ub; i2++) {
      b_varThetaddot += varThetaddot_data[i2] * Bbrev_data[i2];
    }
    c_a = sin(b_Q1);
    htau0[0] = 0.0;
    htau0[3] = -tau0.contents[2];
    htau0[6] = tau0.contents[1];
    htau0[1] = tau0.contents[2];
    htau0[4] = 0.0;
    htau0[7] = -tau0.contents[0];
    htau0[2] = -tau0.contents[1];
    htau0[5] = tau0.contents[0];
    htau0[8] = 0.0;
    d12 = 3.0 * ((double)n + 1.0) + 1.0;
    for (i2 = 0; i2 < 3; i2++) {
      i4 = (int)(d12 + (double)i2) - 1;
      dv2[3 * i2] = R0_data[3 * i4];
      dv2[3 * i2 + 1] = R0_data[3 * i4 + 1];
      dv2[3 * i2 + 2] = R0_data[3 * i4 + 2];
    }
    memcpy(&R0_.contents[0], &dv2[0], 9U * sizeof(double));
    d_a = 1.0 - cos(b_Q1);
    for (i2 = 0; i2 < 3; i2++) {
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] =
            (htau0[i2] * htau0[3 * i4] + htau0[i2 + 3] * htau0[3 * i4 + 1]) +
            htau0[i2 + 6] * htau0[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      Theta.contents[i2] =
          (I3.contents[i2] + c_a * htau0[i2]) + d_a * b_htau0[i2];
    }
    c_tau = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d12 = dv[i2];
      PP.contents[3 * i2] = I3.contents[3 * i2] - b_d * d12;
      i4 = 3 * i2 + 1;
      PP.contents[i4] = I3.contents[i4] - d10 * d12;
      i4 = 3 * i2 + 2;
      PP.contents[i4] = I3.contents[i4] - d11 * d12;
      c_tau += b_tau.contents[i2] * tau0.contents[i2];
    }
    tau_dot_tau0.contents = c_tau;
    tau0_x_tau.contents[0] = tau0.contents[1] * d11 - d10 * tau0.contents[2];
    tau0_x_tau.contents[1] = b_d * tau0.contents[2] - tau0.contents[0] * d11;
    tau0_x_tau.contents[2] = tau0.contents[0] * d10 - b_d * tau0.contents[1];
    y_tmp = 1.0 / nxip.contents;
    for (i2 = 0; i2 < 3; i2++) {
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] =
            (htau0[i4] * htau0[3 * i2] + htau0[i4 + 3] * htau0[3 * i2 + 1]) +
            htau0[i4 + 6] * htau0[3 * i2 + 2];
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      b_I3[3 * i2] =
          (I3.contents[i2] + c_a * htau0[i2]) + d_a * b_htau0[3 * i2];
      b_loop_ub = 3 * i2 + 1;
      b_I3[b_loop_ub] = (I3.contents[i2 + 3] + c_a * htau0[i2 + 3]) +
                        d_a * b_htau0[b_loop_ub];
      b_loop_ub = 3 * i2 + 2;
      b_I3[b_loop_ub] = (I3.contents[i2 + 6] + c_a * htau0[i2 + 6]) +
                        d_a * b_htau0[b_loop_ub];
      c_W[i2] = 2.0 * tau0.contents[i2] + dv[i2];
    }
    dv1[0] = 0.0;
    dv1[3] = -d11;
    dv1[6] = d10;
    dv1[1] = d11;
    dv1[4] = 0.0;
    dv1[7] = -b_d;
    dv1[2] = -d10;
    dv1[5] = b_d;
    dv1[8] = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d12 = dv2[3 * i2];
      d13 = dv2[3 * i2 + 1];
      d14 = dv2[3 * i2 + 2];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i4 + 3 * i2] =
            c_W[i4] * tau0_x_tau.contents[i2] / (c_tau + 1.0);
        c_y[i2 + 3 * i4] =
            (y_tmp * d12 * b_I3[3 * i4] + y_tmp * d13 * b_I3[3 * i4 + 1]) +
            y_tmp * d14 * b_I3[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      dv1[i2] -= b_htau0[i2];
    }
    for (i2 = 0; i2 < 3; i2++) {
      d12 = 0.0;
      d13 = c_y[i2];
      d14 = c_y[i2 + 3];
      d15 = c_y[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        d16 =
            (d13 * dv1[3 * i4] + d14 * dv1[3 * i4 + 1]) + d15 * dv1[3 * i4 + 2];
        b_W.contents[i2 + 3 * i4] = d16;
        d12 += d16 * xipd[i4];
      }
      d_W[i2] = d12 + e1.contents[i2] * b_varThetadot;
    }
    /*  angular momentum */
    d12 = d_W[0];
    d13 = d_W[1];
    d14 = d_W[2];
    for (i2 = 0; i2 < 3; i2++) {
      l[i2] = (II[i2] * d12 + II[i2 + 3] * d13) + II[i2 + 6] * d14;
      c_W[i2] = 2.0 * tau0.contents[i2] + dv[i2];
    }
    d12 = c_W[0];
    d13 = c_W[1];
    d14 = c_W[2];
    for (i2 = 0; i2 < 3; i2++) {
      d15 = tau0.contents[i2];
      G1.contents[3 * i2] = d12 * d15 / (c_tau + 1.0) - I3.contents[3 * i2];
      i4 = 3 * i2 + 1;
      G1.contents[i4] = d13 * d15 / (c_tau + 1.0) - I3.contents[i4];
      i4 = 3 * i2 + 2;
      G1.contents[i4] = d14 * d15 / (c_tau + 1.0) - I3.contents[i4];
    }
    for (i2 = 0; i2 < 9; i2++) {
      d12 = htau0[i2] / (c_tau + 1.0);
      dv1[i2] = d12;
      M1_contents[i2] = d12;
    }
    d12 = xipdd[0];
    d13 = xipdd[1];
    d14 = xipdd[2];
    for (b_i = 0; b_i < 3; b_i++) {
      d15 = tau0.contents[b_i];
      c_e1[b_i] = 2.0 * d15 + b_tau.contents[b_i];
      b_Q1 = b_d * d15 / (c_tau + 1.0) - I3.contents[3 * b_i];
      LL[3 * b_i] = b_Q1;
      M3_contents[3 * b_i] = b_Q1;
      i2 = 3 * b_i + 1;
      b_Q1 = d10 * d15 / (c_tau + 1.0) - I3.contents[i2];
      LL[i2] = b_Q1;
      M3_contents[i2] = b_Q1;
      i2 = 3 * b_i + 2;
      b_Q1 = d11 * d15 / (c_tau + 1.0) - I3.contents[i2];
      LL[i2] = b_Q1;
      M3_contents[i2] = b_Q1;
      d_W[b_i] = ((b_W.contents[b_i] * d12 + b_W.contents[b_i + 3] * d13) +
                  b_W.contents[b_i + 6] * d14) +
                 e1.contents[b_i] * b_varThetaddot;
    }
    b_M2_contents = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d12 = 0.0;
      d13 = 0.0;
      d14 = Theta.contents[i2];
      d15 = Theta.contents[i2 + 3];
      d16 = Theta.contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        d12 += II[i2 + 3 * i4] * d_W[i4];
        d13 += ((d14 * dv2[3 * i4] + d15 * dv2[3 * i4 + 1]) +
                d16 * dv2[3 * i4 + 2]) *
               l[i4];
      }
      b_y[i2] = d13;
      Ialph[i2] = d12;
      b_M2_contents += c_e1[i2] * d13;
    }
    d12 = l[0];
    d13 = l[1];
    d14 = l[2];
    for (i2 = 0; i2 < 3; i2++) {
      d15 = b_y[i2];
      e_a[3 * i2] = b_M2_contents * LL[3 * i2] - b_d * d15;
      i4 = 3 * i2 + 1;
      e_a[i4] = b_M2_contents * LL[i4] - d10 * d15;
      i5 = 3 * i2 + 2;
      e_a[i5] = b_M2_contents * LL[i5] - d11 * d15;
      d_W[i2] = -((b_W.contents[3 * i2] * d12 + b_W.contents[i4] * d13) +
                  b_W.contents[i5] * d14);
    }
    dv2[0] = 0.0;
    dv2[3] = -b_y[2];
    dv2[6] = b_y[1];
    dv2[1] = b_y[2];
    dv2[4] = 0.0;
    dv2[7] = -b_y[0];
    dv2[2] = -b_y[1];
    dv2[5] = b_y[0];
    dv2[8] = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d12 = M1_contents[i2];
      d13 = M1_contents[i2 + 3];
      d14 = M1_contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] =
            (d12 * e_a[3 * i4] + d13 * e_a[3 * i4 + 1]) + d14 * e_a[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      dv2[i2] += b_htau0[i2];
    }
    for (i2 = 0; i2 < 3; i2++) {
      d12 = b_tau.contents[i2];
      b_I3[3 * i2] = I3.contents[3 * i2] - b_d * d12;
      b_loop_ub = 3 * i2 + 1;
      b_I3[b_loop_ub] = I3.contents[b_loop_ub] - d10 * d12;
      b_loop_ub = 3 * i2 + 2;
      b_I3[b_loop_ub] = I3.contents[b_loop_ub] - d11 * d12;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d12 = dv2[i2];
      d13 = dv2[i2 + 3];
      d14 = dv2[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] = (d12 * b_I3[3 * i4] + d13 * b_I3[3 * i4 + 1]) +
                               d14 * b_I3[3 * i4 + 2];
      }
    }
    d12 = d_W[0];
    d13 = d_W[1];
    d14 = d_W[2];
    for (i2 = 0; i2 < 3; i2++) {
      d15 = dv[i2];
      L_WTp_l[3 * i2] = d12 * d15 + y_tmp * b_htau0[3 * i2];
      b_loop_ub = 3 * i2 + 1;
      L_WTp_l[b_loop_ub] = d13 * d15 + y_tmp * b_htau0[b_loop_ub];
      b_loop_ub = 3 * i2 + 2;
      L_WTp_l[b_loop_ub] = d14 * d15 + y_tmp * b_htau0[b_loop_ub];
    }
    for (i2 = 0; i2 < 9; i2++) {
      L_WTp_l[i2] *= y_tmp;
    }
    d12 = 0.0 - l[2];
    d13 = l[1];
    for (i2 = 0; i2 < 3; i2++) {
      L_WTthet_l[i2] =
          b_W.contents[3 * i2 + 1] * d12 + b_W.contents[3 * i2 + 2] * d13;
    }
    L_Wp(&tau0_x_tau, &tau_dot_tau0, &I3, &nxip, &G1, &PP, &b_W, &b_tau, &R0_,
         &Theta, xipd, L_Wp_xipd);
    d12 = xipd[0];
    d13 = xipd[1];
    d14 = xipd[2];
    for (i2 = 0; i2 < 3; i2++) {
      b_y[i2] = (b_W.contents[i2] * d12 + b_W.contents[i2 + 3] * d13) +
                b_W.contents[i2 + 6] * d14;
    }
    L_Wthet_xipd[1] = -(0.0 - b_y[2]);
    L_Wthet_xipd[2] = -b_y[1];
    /*  Wdot*xipd */
    d12 = xipd[0];
    d13 = xipd[1];
    d14 = xipd[2];
    for (i2 = 0; i2 < 3; i2++) {
      b_y[i2] = ((L_Wp_xipd[i2] * d12 + L_Wp_xipd[i2 + 3] * d13) +
                 L_Wp_xipd[i2 + 6] * d14) +
                L_Wthet_xipd[i2] * b_varThetadot;
    }
    mtimes(Pddot, B, b_xi0p);
    d12 = b_y[0];
    d13 = b_y[1];
    d14 = b_y[2];
    for (i2 = 0; i2 < 3; i2++) {
      b_xi0p[i2] *= rho;
      b_y_tmp[3 * i2] = L_Wp_xipd[i2];
      b_y_tmp[3 * i2 + 1] = L_Wp_xipd[i2 + 3];
      b_y_tmp[3 * i2 + 2] = L_Wp_xipd[i2 + 6];
      IWd_xipd[i2] = (II[i2] * d12 + II[i2 + 3] * d13) + II[i2 + 6] * d14;
    }
    mtimes(Pdot, B, c_W);
    c_y_tmp = Ialph[0] + IWd_xipd[0];
    d12 = Ialph[1] + IWd_xipd[1];
    d13 = Ialph[2] + IWd_xipd[2];
    for (i2 = 0; i2 < 3; i2++) {
      d_W[i2] =
          (b_W.contents[3 * i2] * c_y_tmp + b_W.contents[3 * i2 + 1] * d12) +
          b_W.contents[3 * i2 + 2] * d13;
    }
    d12 = xipd[0];
    d13 = xipd[1];
    d14 = xipd[2];
    d15 = l[0];
    d16 = l[1];
    d17 = l[2];
    for (i2 = 0; i2 < 3; i2++) {
      F1[i2] = (((d_W[i2] + ((L_WTp_l[i2] * d12 + L_WTp_l[i2 + 3] * d13) +
                             L_WTp_l[i2 + 6] * d14)) +
                 L_WTthet_l[i2] * b_varThetadot) -
                ((b_y_tmp[i2] * d15 + b_y_tmp[i2 + 3] * d16) +
                 b_y_tmp[i2 + 6] * d17)) +
               a * c_W[i2];
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i3; loop_ub_tmp++) {
      b_i = (int)(nel_ + (double)loop_ub_tmp) - 1;
      c_a = B_data[b_i];
      d_a = Bp_data[b_i];
      b_Q1 = (double)loop_ub_tmp * 3.0;
      F__data[(int)(b_Q1 + 1.0) - 1] =
          (c_a * b_xi0p[0] + d_a * F1[0]) * nxi0p * wg_data[n];
      F__data[(int)(b_Q1 + 2.0) - 1] =
          (c_a * b_xi0p[1] + d_a * F1[1]) * nxi0p * wg_data[n];
      F__data[(int)(b_Q1 + 3.0) - 1] =
          (c_a * b_xi0p[2] + d_a * F1[2]) * nxi0p * wg_data[n];
    }
    d12 = nel_data[n];
    d13 = (d12 - 1.0) * 3.0 + 1.0;
    d12 = (d12 + d) * 3.0;
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
    }
    if (i4 - i2 == F_->size[0]) {
      b_loop_ub = i6 - i5;
      i4 = b_F->size[0] * b_F->size[1];
      b_F->size[0] = 1;
      b_F->size[1] = b_loop_ub;
      emxEnsureCapacity_real_T(b_F, i4);
      b_F_data = b_F->data;
      for (i4 = 0; i4 < b_loop_ub; i4++) {
        b_F_data[i4] = F_data[i2 + i4] + F__data[i4];
      }
      loop_ub = b_F->size[1];
      for (i2 = 0; i2 < loop_ub; i2++) {
        F_data[i5 + i2] = b_F_data[i2];
      }
    } else {
      binary_expand_op_14(F, i5, i6, i2, i4 - 1, F_);
      F_data = F->data;
    }
    b_varThetaddot =
        c_y_tmp - (L_Wthet_xipd[1] * l[1] + L_Wthet_xipd[2] * l[2]);
    if (dbrev < 0.0) {
      b_F->size[1] = 0;
    } else {
      i2 = b_F->size[0] * b_F->size[1];
      b_F->size[0] = 1;
      b_F->size[1] = (int)dbrev + 1;
      emxEnsureCapacity_real_T(b_F, i2);
      b_F_data = b_F->data;
      loop_ub = (int)dbrev;
      for (i2 = 0; i2 <= loop_ub; i2++) {
        b_F_data[i2] = i2;
      }
    }
    i2 = b_F->size[0] * b_F->size[1];
    b_F->size[0] = 1;
    emxEnsureCapacity_real_T(b_F, i2);
    b_F_data = b_F->data;
    b_Q1 = nel_data[n];
    loop_ub = b_F->size[1] - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
      b_F_data[i2] += b_Q1;
    }
    if (dbrev < 0.0) {
      y->size[1] = 0;
    } else {
      i2 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)dbrev + 1;
      emxEnsureCapacity_real_T(y, i2);
      y_data = y->data;
      loop_ub = (int)dbrev;
      for (i2 = 0; i2 <= loop_ub; i2++) {
        y_data[i2] = i2;
      }
    }
    i2 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_real_T(y, i2);
    y_data = y->data;
    b_Q1 = nel_data[n];
    loop_ub = y->size[1] - 1;
    for (i2 = 0; i2 <= loop_ub; i2++) {
      y_data[i2] += b_Q1;
    }
    i2 = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = y->size[1];
    emxEnsureCapacity_int32_T(r, i2);
    r4 = r->data;
    loop_ub = y->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      r4[i2] = (int)y_data[i2];
    }
    i2 = b_mu->size[0] * b_mu->size[1];
    b_mu->size[0] = 1;
    b_mu->size[1] = r->size[1];
    emxEnsureCapacity_real_T(b_mu, i2);
    b_mu_data = b_mu->data;
    loop_ub = r->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      d14 = b_F_data[i2];
      b_mu_data[i2] = mu_data[(int)d14 - 1] + Bbrev_data[(int)d14 - 1] *
                                                  b_varThetaddot * nxi0p *
                                                  wg_data[n];
    }
    loop_ub = b_mu->size[1];
    for (i2 = 0; i2 < loop_ub; i2++) {
      mu_data[r4[i2] - 1] = b_mu_data[i2];
    }
    d_tau = 0.0;
    b_Q1 = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = 0.0;
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = i2 + 3 * i4;
        d14 += II[b_loop_ub] * L_Wthet_xipd[i4];
        d15 = II[i2];
        d16 = d15 * L_Wp_xipd[3 * i4];
        d17 = d15 * b_W.contents[3 * i4];
        d15 = II[i2 + 3];
        i5 = 3 * i4 + 1;
        d16 += d15 * L_Wp_xipd[i5];
        d17 += d15 * b_W.contents[i5];
        d15 = II[i2 + 6];
        i5 = 3 * i4 + 2;
        d16 += d15 * L_Wp_xipd[i5];
        d17 += d15 * b_W.contents[i5];
        L_lpd[b_loop_ub] = d17;
        L_lp[b_loop_ub] = d16;
      }
      L_lthet[i2] = d14;
      d_tau += b_tau.contents[i2] * xipd[i2];
      d14 = (PP.contents[i2] * xipd[0] + PP.contents[i2 + 3] * xipd[1]) +
            PP.contents[i2 + 6] * xipd[2];
      b_xi0p[i2] = d14;
      b_Q1 += tau0.contents[i2] * d14;
    }
    c_a = b_Q1 / (c_tau + 1.0);
    for (i2 = 0; i2 < 3; i2++) {
      f_a[i2] = c_a * b_tau.contents[i2] - b_xi0p[i2];
      d14 = xipd[i2];
      tau[3 * i2] = b_d * d14;
      tau[3 * i2 + 1] = d10 * d14;
      tau[3 * i2 + 2] = d11 * d14;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = tau[i2];
      d15 = tau[i2 + 3];
      d16 = tau[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        e_a[i4 + 3 * i2] = f_a[i4] * c_e1[i2];
        b_htau0[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
            d16 * PP.contents[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      e_a[i2] -= b_htau0[i2];
    }
    dv2[0] = -0.0;
    dv2[3] = b_xi0p[2];
    dv2[6] = -b_xi0p[1];
    dv2[1] = -b_xi0p[2];
    dv2[4] = -0.0;
    dv2[7] = b_xi0p[0];
    dv2[2] = b_xi0p[1];
    dv2[5] = -b_xi0p[0];
    dv2[8] = -0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = dv1[i2];
      d15 = dv1[i2 + 3];
      d16 = dv1[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] =
            (d14 * e_a[3 * i4] + d15 * e_a[3 * i4 + 1]) + d16 * e_a[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      dv2[i2] = (dv2[i2] + b_htau0[i2]) / nxip.contents;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = dv2[i2];
      d15 = dv2[i2 + 3];
      d16 = dv2[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] =
            (d14 * Theta.contents[3 * i4] + d15 * Theta.contents[3 * i4 + 1]) +
            d16 * Theta.contents[3 * i4 + 2];
      }
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        Q2_WTp_xipd[i2 + 3 * i4] =
            -d_tau * b_W.contents[i4 + 3 * i2] +
            ((d14 * R0_.contents[3 * i4] + d15 * R0_.contents[3 * i4 + 1]) +
             d16 * R0_.contents[3 * i4 + 2]);
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      Q2_WTp_xipd[i2] *= y_tmp;
    }
    dv2[0] = 0.0;
    dv2[3] = -0.0;
    dv2[6] = 0.0;
    dv2[1] = 0.0;
    dv2[4] = 0.0;
    dv2[7] = -1.0;
    dv2[2] = -0.0;
    dv2[5] = 1.0;
    dv2[8] = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_W.contents[3 * i2];
      d15 = b_W.contents[3 * i2 + 1];
      d16 = b_W.contents[3 * i2 + 2];
      for (i4 = 0; i4 < 3; i4++) {
        Q2_WTthet_thetd[i2 + 3 * i4] =
            (d14 * dv2[3 * i4] + d15 * dv2[3 * i4 + 1]) + d16 * dv2[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      Q2_WTthet_thetd[i2] *= b_varThetadot;
    }
    Q1_Wpp(&tau0_x_tau, &tau_dot_tau0, &I3, &nxip, &G1, &PP, &b_W, &b_tau, &R0_,
           &Theta, &tau0, xipd, xipd, Wd_xipd_p);
    Q1_Wthetp(&tau0_x_tau, &tau_dot_tau0, &I3, &nxip, &G1, &PP, &b_W, &b_tau,
              &R0_, &Theta, &e1, xipd, b_varThetadot, dv2);
    for (i2 = 0; i2 < 9; i2++) {
      Wd_xipd_p[i2] += dv2[i2];
    }
    htau0[0] = 0.0;
    htau0[3] = -0.0;
    htau0[6] = 0.0;
    htau0[1] = 0.0;
    htau0[4] = 0.0;
    htau0[7] = -1.0;
    htau0[2] = -0.0;
    htau0[5] = 1.0;
    htau0[8] = 0.0;
    d14 = xipd[0];
    d15 = xipd[1];
    d16 = xipd[2];
    for (i2 = 0; i2 < 3; i2++) {
      b_y[i2] = (b_W.contents[i2] * d14 + b_W.contents[i2 + 3] * d15) +
                b_W.contents[i2 + 6] * d16;
    }
    F1[1] = 0.0 - b_y[2];
    F1[2] = b_y[1];
    Q1_Wpthet(&e1, &b_W, &tau0_x_tau, &tau_dot_tau0, &I3, &nxip, &G1, &PP,
              &b_tau, &R0_, &Theta, xipd, xipd, c_W);
    for (i2 = 0; i2 < 9; i2++) {
      htau0[i2] = -htau0[i2];
    }
    d14 = -F1[1];
    d15 = -F1[2];
    for (i2 = 0; i2 < 3; i2++) {
      Wd_xipd_thet[i2] =
          c_W[i2] + (htau0[i2 + 3] * d14 + htau0[i2 + 6] * d15) * b_varThetadot;
    }
    Q2_Wp(&PP, &G1, &nxip, &tau0_x_tau, &tau_dot_tau0, &b_tau, &b_W, &R0_,
          &Theta, xipd, Q2_Wp_xipd);
    /*  QQ = -mycross(e1,W*pie); */
    /* %%% F_ij %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    memcpy(&dv2[0], &b_W.contents[0], 9U * sizeof(double));
    L_Wp(&tau0_x_tau, &tau_dot_tau0, &I3, &nxip, &G1, &PP, &b_W, &b_tau, &R0_,
         &Theta, xipdd, b);
    b_M2_contents = 0.0;
    for (b_i = 0; b_i < 3; b_i++) {
      c_W[b_i] = b_tau.contents[b_i];
      d14 = 0.0;
      d15 = 0.0;
      d16 = Theta.contents[b_i];
      d17 = Theta.contents[b_i + 3];
      d18 = Theta.contents[b_i + 6];
      for (i2 = 0; i2 < 3; i2++) {
        d19 = Ialph[i2];
        d14 += b_W.contents[i2 + 3 * b_i] * d19;
        d15 += ((d16 * R0_.contents[3 * i2] + d17 * R0_.contents[3 * i2 + 1]) +
                d18 * R0_.contents[3 * i2 + 2]) *
               d19;
      }
      b_y[b_i] = d15;
      Q1[b_i] = -d14;
      b_M2_contents += c_e1[b_i] * d15;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_y[i2];
      e_a[3 * i2] = b_M2_contents * M3_contents[3 * i2] - b_d * d14;
      b_i = 3 * i2 + 1;
      e_a[b_i] = b_M2_contents * M3_contents[b_i] - d10 * d14;
      b_i = 3 * i2 + 2;
      e_a[b_i] = b_M2_contents * M3_contents[b_i] - d11 * d14;
    }
    b_htau0[0] = 0.0;
    b_htau0[3] = -b_y[2];
    b_htau0[6] = b_y[1];
    b_htau0[1] = b_y[2];
    b_htau0[4] = 0.0;
    b_htau0[7] = -b_y[0];
    b_htau0[2] = -b_y[1];
    b_htau0[5] = b_y[0];
    b_htau0[8] = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = dv1[i2];
      d15 = dv1[i2 + 3];
      d16 = dv1[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * e_a[3 * i4] + d15 * e_a[3 * i4 + 1]) + d16 * e_a[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      b_htau0[i2] += b_I3[i2];
    }
    d_a = 0.0;
    b_M2_contents = 0.0;
    M2_contents_tmp = 0.0;
    htau0[0] = 0.0;
    htau0[4] = 0.0;
    htau0[8] = 0.0;
    d_tau = 0.0;
    a_tmp = -1.0 / nxip.contents;
    c_xi0p = 0.0;
    d_xi0p = 0.0;
    b_varThetaddot = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      c_a = xipd[i2] * b_tau.contents[i2];
      d_a += c_a;
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      d17 = Theta.contents[i2];
      d18 = Theta.contents[i2 + 3];
      d19 = Theta.contents[i2 + 6];
      d20 = 0.0;
      d21 = 0.0;
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = 3 * i4 + 1;
        b_i = 3 * i4 + 2;
        c_y[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[b_loop_ub]) +
            d16 * PP.contents[b_i];
        c_y_tmp =
            ((d17 * R0_.contents[3 * i4] + d18 * R0_.contents[b_loop_ub]) +
             d19 * R0_.contents[b_i]) *
            l[i4];
        d20 += c_y_tmp;
        d21 += c_y_tmp;
      }
      d_y[i2] = d21;
      b_y[i2] = d20;
      d14 = c_e1[i2];
      b_M2_contents += d14 * d20;
      b_Q1 = d14 * d21;
      M2_contents_tmp += b_Q1;
      M2_contents = M2_contents_tmp;
      d14 = (PP.contents[i2] * xipd[0] + PP.contents[i2 + 3] * xipd[1]) +
            PP.contents[i2 + 6] * xipd[2];
      b_xi0p[i2] = d14;
      d_tau += c_a;
      c_xi0p += d14 * d21;
      d_xi0p += d14 * tau0.contents[i2];
      b_varThetaddot += b_Q1;
    }
    htau0[3] = -d_y[2];
    htau0[6] = d_y[1];
    htau0[1] = d_y[2];
    htau0[7] = -d_y[0];
    htau0[2] = -d_y[1];
    htau0[5] = d_y[0];
    c_a = b_varThetaddot / (c_tau + 1.0);
    d14 = l[0];
    d15 = l[1];
    d16 = l[2];
    for (i2 = 0; i2 < 3; i2++) {
      d17 = b_y[i2];
      e_a[3 * i2] = b_M2_contents * M3_contents[3 * i2] - b_d * d17;
      i4 = 3 * i2 + 1;
      e_a[i4] = b_M2_contents * M3_contents[i4] - d10 * d17;
      i5 = 3 * i2 + 2;
      e_a[i5] = b_M2_contents * M3_contents[i5] - d11 * d17;
      d_W[i2] = -((b_W.contents[3 * i2] * d14 + b_W.contents[i4] * d15) +
                  b_W.contents[i5] * d16);
    }
    b_htau0[0] = 0.0;
    b_htau0[3] = -b_y[2];
    b_htau0[6] = b_y[1];
    b_htau0[1] = b_y[2];
    b_htau0[4] = 0.0;
    b_htau0[7] = -b_y[0];
    b_htau0[2] = -b_y[1];
    b_htau0[5] = b_y[0];
    b_htau0[8] = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = dv1[i2];
      d15 = dv1[i2 + 3];
      d16 = dv1[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * e_a[3 * i4] + d15 * e_a[3 * i4 + 1]) + d16 * e_a[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      b_htau0[i2] += b_I3[i2];
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
            d16 * PP.contents[3 * i4 + 2];
      }
    }
    d14 = d_W[0];
    d15 = d_W[1];
    d16 = d_W[2];
    for (i2 = 0; i2 < 3; i2++) {
      d17 = b_tau.contents[i2];
      W[3 * i2] = d14 * d17 + y_tmp * b_I3[3 * i2];
      b_loop_ub = 3 * i2 + 1;
      W[b_loop_ub] = d15 * d17 + y_tmp * b_I3[b_loop_ub];
      b_loop_ub = 3 * i2 + 2;
      W[b_loop_ub] = d16 * d17 + y_tmp * b_I3[b_loop_ub];
    }
    for (i2 = 0; i2 < 9; i2++) {
      W[i2] = -(y_tmp * W[i2]);
    }
    d14 = xipd[0];
    d15 = xipd[1];
    d16 = xipd[2];
    d17 = l[0];
    d18 = l[1];
    d19 = l[2];
    for (i2 = 0; i2 < 3; i2++) {
      d20 = b_tau.contents[i2];
      b_I3[3 * i2] = d14 * d20 + d_a * I3.contents[3 * i2];
      b_loop_ub = 3 * i2 + 1;
      b_I3[b_loop_ub] = d15 * d20 + d_a * I3.contents[b_loop_ub];
      b_i = 3 * i2 + 2;
      b_I3[b_i] = d16 * d20 + d_a * I3.contents[b_i];
      d_W[i2] = (b_W.contents[3 * i2] * d17 + b_W.contents[b_loop_ub] * d18) +
                b_W.contents[b_i] * d19;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = W[i2];
      d15 = W[i2 + 3];
      d16 = W[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i4 + 3 * i2] = d_W[i4] * xipd[i2];
        e_a[i2 + 3 * i4] = (d14 * b_I3[3 * i4] + d15 * b_I3[3 * i4 + 1]) +
                           d16 * b_I3[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        W[i2 + 3 * i4] =
            ((d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
             d16 * PP.contents[3 * i4 + 2]) /
            nxip.contents;
        b_i = i4 + 3 * i2;
        b_a[b_i] = M2_contents * LL[b_i] - b_tau.contents[i4] * d_y[i2];
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = M1_contents[i2];
      d15 = M1_contents[i2 + 3];
      d16 = M1_contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = i2 + 3 * i4;
        b_htau0[b_loop_ub] =
            htau0[b_loop_ub] + ((d14 * b_a[3 * i4] + d15 * b_a[3 * i4 + 1]) +
                                d16 * b_a[3 * i4 + 2]);
      }
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
            d16 * PP.contents[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      b_I3[i2] = -(y_tmp * b_I3[i2]);
    }
    d14 = xipd[0];
    d15 = xipd[1];
    d16 = xipd[2];
    for (i2 = 0; i2 < 3; i2++) {
      f_a[i2] = (b_I3[i2] * d14 + b_I3[i2 + 3] * d15) + b_I3[i2 + 6] * d16;
    }
    for (i2 = 0; i2 < 9; i2++) {
      b_htau0[i2] = M1_contents[i2] / nxip.contents;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        i5 = i4 + 3 * i2;
        i6 = (int)I3.contents[i5];
        b_M1_contents[i2 + 3 * i4] =
            (d14 * M3_contents[3 * i4] + d15 * M3_contents[3 * i4 + 1]) +
            d16 * M3_contents[3 * i4 + 2];
        d17 = b_xi0p[i4];
        b_I3[i5] = (d17 * d_y[i2] + c_xi0p * (double)i6) -
                   c_a * (d17 * tau0.contents[i2] + d_xi0p * (double)i6);
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_M1_contents[i2];
      d15 = b_M1_contents[i2 + 3];
      d16 = b_M1_contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] = (d14 * b_I3[3 * i4] + d15 * b_I3[3 * i4 + 1]) +
                               d16 * b_I3[3 * i4 + 2];
        b_a[i4 + 3 * i2] = f_a[i4] * b_tau.contents[i2];
      }
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_M1_contents[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
            d16 * PP.contents[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = d_y[i2];
      b_I3[3 * i2] = M2_contents_tmp * LL[3 * i2] - b_d * d14;
      d15 = xipd[i2];
      tau[3 * i2] = b_d * d15 + d_tau * I3.contents[3 * i2];
      b_loop_ub = 3 * i2 + 1;
      b_I3[b_loop_ub] = M2_contents_tmp * LL[b_loop_ub] - d10 * d14;
      tau[b_loop_ub] = d10 * d15 + d_tau * I3.contents[b_loop_ub];
      b_loop_ub = 3 * i2 + 2;
      b_I3[b_loop_ub] = M2_contents_tmp * LL[b_loop_ub] - d11 * d14;
      tau[b_loop_ub] = d11 * d15 + d_tau * I3.contents[b_loop_ub];
    }
    for (i2 = 0; i2 < 9; i2++) {
      tau[i2] *= a_tmp;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = M1_contents[i2];
      d15 = M1_contents[i2 + 3];
      d16 = M1_contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = i2 + 3 * i4;
        b_htau0[b_loop_ub] =
            htau0[b_loop_ub] + ((d14 * b_I3[3 * i4] + d15 * b_I3[3 * i4 + 1]) +
                                d16 * b_I3[3 * i4 + 2]);
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = tau[i2];
      d15 = tau[i2 + 3];
      d16 = tau[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
            d16 * PP.contents[3 * i4 + 2];
      }
    }
    F1[0] = 0.0;
    F1[1] = 0.0 - l[2];
    F1[2] = l[1];
    for (i2 = 0; i2 < 3; i2++) {
      d14 = 0.0;
      d15 = b_htau0[i2];
      d16 = b_htau0[i2 + 3];
      d17 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = i2 + 3 * i4;
        xi0p[b_loop_ub] = ((e_a[b_loop_ub] - W[b_loop_ub]) +
                           ((b_a[b_loop_ub] + b_M1_contents[b_loop_ub]) +
                            ((d15 * b_I3[3 * i4] + d16 * b_I3[3 * i4 + 1]) +
                             d17 * b_I3[3 * i4 + 2])) /
                               nxip.contents) /
                          nxip.contents;
        d14 += b_W.contents[i4 + 3 * i2] * F1[i4];
      }
      d_W[i2] = -d14;
    }
    b_M2_contents = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = 0.0;
      d15 = Theta.contents[i2];
      d16 = Theta.contents[i2 + 3];
      d17 = Theta.contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        htau0[i4 + 3 * i2] = d_W[i4] * b_tau.contents[i2];
        d14 += ((d15 * R0_.contents[3 * i4] + d16 * R0_.contents[3 * i4 + 1]) +
                d17 * R0_.contents[3 * i4 + 2]) *
               F1[i4];
      }
      b_y[i2] = d14;
      b_M2_contents += c_e1[i2] * d14;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_y[i2];
      e_a[3 * i2] = b_M2_contents * M3_contents[3 * i2] - b_d * d14;
      b_i = 3 * i2 + 1;
      e_a[b_i] = b_M2_contents * M3_contents[b_i] - d10 * d14;
      b_i = 3 * i2 + 2;
      e_a[b_i] = b_M2_contents * M3_contents[b_i] - d11 * d14;
    }
    b_htau0[0] = 0.0;
    b_htau0[3] = -b_y[2];
    b_htau0[6] = b_y[1];
    b_htau0[1] = b_y[2];
    b_htau0[4] = 0.0;
    b_htau0[7] = -b_y[0];
    b_htau0[2] = -b_y[1];
    b_htau0[5] = b_y[0];
    b_htau0[8] = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = dv1[i2];
      d15 = dv1[i2 + 3];
      d16 = dv1[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * e_a[3 * i4] + d15 * e_a[3 * i4 + 1]) + d16 * e_a[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      b_htau0[i2] += b_I3[i2];
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        LL[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
            d16 * PP.contents[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      LL[i2] *= y_tmp;
    }
    b_M2_contents = 0.0;
    for (b_i = 0; b_i < 3; b_i++) {
      F1[b_i] = b_tau.contents[b_i];
      d14 = 0.0;
      d15 = 0.0;
      d16 = Theta.contents[b_i];
      d17 = Theta.contents[b_i + 3];
      d18 = Theta.contents[b_i + 6];
      for (i2 = 0; i2 < 3; i2++) {
        d19 = IWd_xipd[i2];
        d14 += b_W.contents[i2 + 3 * b_i] * d19;
        d15 += ((d16 * R0_.contents[3 * i2] + d17 * R0_.contents[3 * i2 + 1]) +
                d18 * R0_.contents[3 * i2 + 2]) *
               d19;
      }
      b_y[b_i] = d15;
      b_xi0p[b_i] = -d14;
      b_M2_contents += c_e1[b_i] * d15;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_y[i2];
      e_a[3 * i2] = b_M2_contents * M3_contents[3 * i2] - b_d * d14;
      b_i = 3 * i2 + 1;
      e_a[b_i] = b_M2_contents * M3_contents[b_i] - d10 * d14;
      b_i = 3 * i2 + 2;
      e_a[b_i] = b_M2_contents * M3_contents[b_i] - d11 * d14;
    }
    b_htau0[0] = 0.0;
    b_htau0[3] = -b_y[2];
    b_htau0[6] = b_y[1];
    b_htau0[1] = b_y[2];
    b_htau0[4] = 0.0;
    b_htau0[7] = -b_y[0];
    b_htau0[2] = -b_y[1];
    b_htau0[5] = b_y[0];
    b_htau0[8] = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = dv1[i2];
      d15 = dv1[i2 + 3];
      d16 = dv1[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * e_a[3 * i4] + d15 * e_a[3 * i4 + 1]) + d16 * e_a[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      b_htau0[i2] += b_I3[i2];
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
            d16 * PP.contents[3 * i4 + 2];
      }
    }
    Q1T_Wpp(&tau0_x_tau, &tau_dot_tau0, &I3, &nxip, &G1, &PP, &b_W, &b_tau,
            &R0_, &Theta, &tau0, xipd, l, dv1);
    for (i2 = 0; i2 < 3; i2++) {
      d14 = dv2[3 * i2];
      d15 = dv2[3 * i2 + 1];
      d16 = dv2[3 * i2 + 2];
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = i4 + 3 * i2;
        tau[b_loop_ub] = Q1[i4] * c_W[i2] + y_tmp * c_y[b_loop_ub];
        b_htau0[i2 + 3 * i4] =
            (d14 * II[3 * i4] + d15 * II[3 * i4 + 1]) + d16 * II[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = Q2_WTp_xipd[i2];
      d15 = Q2_WTp_xipd[i2 + 3];
      d16 = Q2_WTp_xipd[i2 + 6];
      d17 = b_htau0[i2];
      d18 = b_htau0[i2 + 3];
      d19 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        i5 = 3 * i4 + 1;
        i6 = 3 * i4 + 2;
        b_loop_ub = i2 + 3 * i4;
        b_M1_contents[b_loop_ub] =
            (d14 * L_lp[3 * i4] + d15 * L_lp[i5]) + d16 * L_lp[i6];
        e_a[b_loop_ub] = (y_tmp * tau[b_loop_ub] +
                          ((d17 * b[3 * i4] + d18 * b[i5]) + d19 * b[i6])) +
                         xi0p[b_loop_ub];
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_W.contents[3 * i2];
      d15 = b_W.contents[3 * i2 + 1];
      d16 = b_W.contents[3 * i2 + 2];
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = i4 + 3 * i2;
        xi0p[b_loop_ub] = b_xi0p[i4] * F1[i2] + y_tmp * b_I3[b_loop_ub];
        W[i2 + 3 * i4] =
            (d14 * II[3 * i4] + d15 * II[3 * i4 + 1]) + d16 * II[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = Q2_WTthet_thetd[i2];
      d15 = Q2_WTthet_thetd[i2 + 3];
      d16 = Q2_WTthet_thetd[i2 + 6];
      d17 = b_y_tmp[i2];
      d18 = b_y_tmp[i2 + 3];
      d19 = b_y_tmp[i2 + 6];
      d20 = W[i2];
      d21 = W[i2 + 3];
      b_varThetaddot = W[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        d_a = L_lp[3 * i4];
        c_a = d14 * d_a;
        b_M2_contents = d17 * d_a;
        i5 = 3 * i4 + 1;
        d_a = L_lp[i5];
        c_a += d15 * d_a;
        b_M2_contents += d18 * d_a;
        i6 = 3 * i4 + 2;
        d_a = L_lp[i6];
        c_a += d16 * d_a;
        b_M2_contents += d19 * d_a;
        b_loop_ub = i2 + 3 * i4;
        tau[b_loop_ub] =
            ((((((e_a[b_loop_ub] + b_M1_contents[b_loop_ub]) +
                 y_tmp * (htau0[b_loop_ub] + LL[b_loop_ub]) * b_varThetadot) +
                c_a) +
               y_tmp * xi0p[b_loop_ub]) +
              ((d20 * Wd_xipd_p[3 * i4] + d21 * Wd_xipd_p[i5]) +
               b_varThetaddot * Wd_xipd_p[i6])) -
             dv1[b_loop_ub]) -
            b_M2_contents;
      }
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i3; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i3; b_loop_ub_tmp++) {
        c_a = Bp_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
              Bp_data[(int)(nel_ + (double)loop_ub_tmp) - 1];
        b_Q1 = (double)b_loop_ub_tmp * 3.0;
        b_varThetaddot = (double)loop_ub_tmp * 3.0;
        for (i2 = 0; i2 < 3; i2++) {
          b_loop_ub = (int)(b_varThetaddot + ((double)i2 + 1.0)) - 1;
          F_ij__data[((int)(b_Q1 + 1.0) + F_ij_->size[0] * b_loop_ub) - 1] =
              c_a * tau[3 * i2] * nxi0p * wg_data[n];
          F_ij__data[((int)(b_Q1 + 2.0) + F_ij_->size[0] * b_loop_ub) - 1] =
              c_a * tau[3 * i2 + 1] * nxi0p * wg_data[n];
          F_ij__data[((int)(b_Q1 + 3.0) + F_ij_->size[0] * b_loop_ub) - 1] =
              c_a * tau[3 * i2 + 2] * nxi0p * wg_data[n];
        }
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
      i7 = 0;
      b_i = 0;
      b_loop_ub_tmp = 0;
      loop_ub_tmp = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
      i7 = (int)d13 - 1;
      b_i = (int)d12;
      b_loop_ub_tmp = (int)d13 - 1;
      loop_ub_tmp = (int)d12;
    }
    loop_ub = i4 - i2;
    b_loop_ub = i6 - i5;
    if ((loop_ub == F_ij_->size[0]) && (b_loop_ub == F_ij_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = loop_ub;
      b_F_ij->size[1] = b_loop_ub;
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      for (i4 = 0; i4 < b_loop_ub; i4++) {
        for (i6 = 0; i6 < loop_ub; i6++) {
          b_F_data[i6 + b_F_ij->size[0] * i4] =
              F_ij_data[(i2 + i6) + F_ij->size[0] * (i5 + i4)] +
              F_ij__data[i6 + F_ij_->size[0] * i4];
        }
      }
      b_i -= i7;
      b_loop_ub = loop_ub_tmp - b_loop_ub_tmp;
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          F_ij_data[(i7 + i4) + F_ij->size[0] * (b_loop_ub_tmp + i2)] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_9(F_ij, i7, b_i, b_loop_ub_tmp, loop_ub_tmp, i2, i4 - 1,
                         i5, i6 - 1, F_ij_);
      F_ij_data = F_ij->data;
    }
    /* %%% F_ib %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    d_tau = 0.0;
    d14 = xipdd[0];
    d15 = xipdd[1];
    d16 = xipdd[2];
    for (i2 = 0; i2 < 3; i2++) {
      b_y[i2] = (b_W.contents[i2] * d14 + b_W.contents[i2 + 3] * d15) +
                b_W.contents[i2 + 6] * d16;
      d_tau += b_tau.contents[i2] * xipd[i2];
    }
    b_Q1 = 0.0;
    d_a = 0.0;
    Q1[0] = 0.0;
    Q1[1] = 0.0 - Ialph[2];
    Q1[2] = Ialph[1];
    Ialph[0] = -0.0;
    Ialph[1] = -(0.0 - b_y[2]);
    Ialph[2] = -b_y[1];
    for (i2 = 0; i2 < 3; i2++) {
      d_a += tau0_x_tau.contents[i2] * xipd[i2];
      d14 = 0.0;
      d15 = 0.0;
      d16 = 0.0;
      d17 = 0.0;
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = i2 + 3 * i4;
        d18 = xipd[i4];
        d14 += PP.contents[b_loop_ub] * d18;
        d15 += b_W.contents[b_loop_ub] * d18;
        d16 += b_W.contents[i4 + 3 * i2] * Q1[i4];
        d17 += ((b_W.contents[3 * i2] * II[3 * i4] +
                 b_W.contents[3 * i2 + 1] * II[3 * i4 + 1]) +
                b_W.contents[3 * i2 + 2] * II[3 * i4 + 2]) *
               Ialph[i4];
      }
      c_W[i2] = d17;
      d_W[i2] = d16;
      d_y[i2] = d15;
      b_xi0p[i2] = d14;
      b_Q1 += tau0.contents[i2] * d14;
    }
    c_a = b_Q1 / (c_tau + 1.0);
    d_a /= c_tau + 1.0;
    d14 = 0.0 - l[2];
    d15 = l[1];
    for (i2 = 0; i2 < 3; i2++) {
      f_a[i2] = c_a * b_tau.contents[i2] - b_xi0p[i2];
      b_e1[i2] = 2.0 * tau0.contents[i2] + dv[i2];
      d16 = xipd[i2];
      tau[3 * i2] = b_d * d16;
      i4 = 3 * i2 + 1;
      tau[i4] = d10 * d16;
      i5 = 3 * i2 + 2;
      tau[i5] = d11 * d16;
      Ialph[i2] = b_W.contents[i4] * d14 + b_W.contents[i5] * d15;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = tau[i2];
      d15 = tau[i2 + 3];
      d16 = tau[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        e_a[i4 + 3 * i2] = f_a[i4] * b_e1[i2];
        b_htau0[i2 + 3 * i4] =
            (d14 * PP.contents[3 * i4] + d15 * PP.contents[3 * i4 + 1]) +
            d16 * PP.contents[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      e_a[i2] -= b_htau0[i2];
    }
    dv2[0] = -0.0;
    dv2[3] = b_xi0p[2];
    dv2[6] = -b_xi0p[1];
    dv2[1] = -b_xi0p[2];
    dv2[4] = -0.0;
    dv2[7] = b_xi0p[0];
    dv2[2] = b_xi0p[1];
    dv2[5] = -b_xi0p[0];
    dv2[8] = -0.0;
    for (i2 = 0; i2 < 3; i2++) {
      d14 = M1_contents[i2];
      d15 = M1_contents[i2 + 3];
      d16 = M1_contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] =
            (d14 * e_a[3 * i4] + d15 * e_a[3 * i4 + 1]) + d16 * e_a[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      dv2[i2] = (dv2[i2] + b_htau0[i2]) / nxip.contents;
    }
    Q1[0] = 0.0;
    Q1[1] = 0.0 - l[2];
    Q1[2] = l[1];
    for (i2 = 0; i2 < 3; i2++) {
      d14 = dv2[i2];
      d15 = dv2[i2 + 3];
      d16 = dv2[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        dv1[i2 + 3 * i4] =
            (d14 * Theta.contents[3 * i4] + d15 * Theta.contents[3 * i4 + 1]) +
            d16 * Theta.contents[3 * i4 + 2];
      }
      d14 = 0.0;
      d15 = dv1[i2];
      d16 = dv1[i2 + 3];
      d17 = dv1[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        d18 = (d15 * R0_.contents[3 * i4] + d16 * R0_.contents[3 * i4 + 1]) +
              d17 * R0_.contents[3 * i4 + 2];
        dv2[i2 + 3 * i4] = d18;
        d14 += d18 * Q1[i4];
      }
      f_a[i2] = -d_tau * Ialph[i2] + d14;
    }
    Q1[1] = 0.0 - l[1];
    Q1[2] = 0.0 - l[2];
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_W.contents[3 * i2];
      d15 = b_W.contents[3 * i2 + 1];
      d16 = b_W.contents[3 * i2 + 2];
      b_y[i2] =
          (((d_W[i2] + c_W[i2]) + y_tmp * f_a[i2]) +
           ((Q2_WTp_xipd[i2] * L_lthet[0] + Q2_WTp_xipd[i2 + 3] * L_lthet[1]) +
            Q2_WTp_xipd[i2 + 6] * L_lthet[2])) +
          (d15 * Q1[1] + d16 * Q1[2]) * b_varThetadot;
      d17 = 0.0;
      for (i4 = 0; i4 < 3; i4++) {
        b_loop_ub = i2 + 3 * i4;
        d17 += Q2_WTthet_thetd[b_loop_ub] * L_lthet[i4];
        W[b_loop_ub] =
            (d14 * II[3 * i4] + d15 * II[3 * i4 + 1]) + d16 * II[3 * i4 + 2];
      }
      Ialph[i2] = d17;
    }
    d14 = Wd_xipd_thet[0];
    d15 = Wd_xipd_thet[1];
    d16 = Wd_xipd_thet[2];
    d17 = 0.0 - IWd_xipd[2];
    d18 = IWd_xipd[1];
    for (i2 = 0; i2 < 3; i2++) {
      d19 = l[i2];
      tau[3 * i2] = b_d * d19;
      i4 = 3 * i2 + 1;
      tau[i4] = d10 * d19;
      i5 = 3 * i2 + 2;
      tau[i5] = d11 * d19;
      c_W[i2] = (W[i2] * d14 + W[i2 + 3] * d15) + W[i2 + 6] * d16;
      d_W[i2] = (b_y[i2] + Ialph[i2]) +
                (b_W.contents[i4] * d17 + b_W.contents[i5] * d18);
    }
    Q1[0] = -0.0;
    Q1[1] = -(0.0 - d_y[2]);
    Q1[2] = -d_y[1];
    dv2[0] = 0.0;
    dv2[3] = -xipd[2];
    dv2[6] = xipd[1];
    dv2[1] = xipd[2];
    dv2[4] = 0.0;
    dv2[7] = -xipd[0];
    dv2[2] = -xipd[1];
    dv2[5] = xipd[0];
    dv2[8] = 0.0;
    for (i2 = 0; i2 < 9; i2++) {
      dv2[i2] += d_a * I3.contents[i2];
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = G1.contents[i2];
      d15 = G1.contents[i2 + 3];
      d16 = G1.contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        b_I3[i2 + 3 * i4] =
            (d14 * dv2[3 * i4] + d15 * dv2[3 * i4 + 1]) + d16 * dv2[3 * i4 + 2];
      }
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = PP.contents[3 * i2];
      d15 = PP.contents[3 * i2 + 1];
      d16 = PP.contents[3 * i2 + 2];
      for (i4 = 0; i4 < 3; i4++) {
        b_htau0[i2 + 3 * i4] =
            (b_I3[i4] * d14 + b_I3[i4 + 3] * d15) + b_I3[i4 + 6] * d16;
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      b_htau0[i2] *= y_tmp;
    }
    Ialph[0] = 0.0;
    Ialph[1] = 0.0 - l[2];
    Ialph[2] = l[1];
    for (i2 = 0; i2 < 3; i2++) {
      d14 = b_htau0[i2];
      d15 = b_htau0[i2 + 3];
      d16 = b_htau0[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        e_a[i2 + 3 * i4] =
            (d14 * Theta.contents[3 * i4] + d15 * Theta.contents[3 * i4 + 1]) +
            d16 * Theta.contents[3 * i4 + 2];
      }
      d14 = 0.0;
      d15 = 0.0;
      d16 = e_a[i2];
      d17 = e_a[i2 + 3];
      d18 = e_a[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        d14 += tau[i2 + 3 * i4] * Q1[i4];
        d15 += ((d16 * R0_.contents[3 * i4] + d17 * R0_.contents[3 * i4 + 1]) +
                d18 * R0_.contents[3 * i4 + 2]) *
               Ialph[i4];
      }
      f_a[i2] = d15;
      b_y[i2] = d14;
    }
    d14 = L_lthet[0];
    d15 = L_lthet[1];
    d16 = L_lthet[2];
    for (i2 = 0; i2 < 3; i2++) {
      Q1[i2] =
          ((d_W[i2] + c_W[i2]) - a_tmp * (b_y[i2] - f_a[i2])) -
          ((b_y_tmp[i2] * d14 + b_y_tmp[i2 + 3] * d15) + b_y_tmp[i2 + 6] * d16);
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i1; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i3; b_loop_ub_tmp++) {
        c_a = Bp_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
              Bbrev_data[(int)(nel_ + (double)loop_ub_tmp) - 1];
        b_Q1 = (double)b_loop_ub_tmp * 3.0;
        F_ib__data[((int)(b_Q1 + 1.0) + F_ib_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[0] * nxi0p * wg_data[n];
        F_ib__data[((int)(b_Q1 + 2.0) + F_ib_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[1] * nxi0p * wg_data[n];
        F_ib__data[((int)(b_Q1 + 3.0) + F_ib_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[2] * nxi0p * wg_data[n];
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
    }
    i5 = r3->size[0];
    r3->size[0] = y->size[1];
    emxEnsureCapacity_uint32_T(r3, i5);
    r5 = r3->data;
    loop_ub = y->size[1];
    for (i5 = 0; i5 < loop_ub; i5++) {
      r5[i5] = (unsigned int)y_data[i5];
    }
    if (d13 > d12) {
      i5 = 0;
      i6 = 0;
    } else {
      i5 = (int)d13 - 1;
      i6 = (int)d12;
    }
    i7 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i7);
    r4 = r1->data;
    loop_ub = r3->size[0];
    for (i7 = 0; i7 < loop_ub; i7++) {
      r4[i7] = (int)r5[i7] - 1;
    }
    loop_ub = i4 - i2;
    if ((loop_ub == F_ib_->size[0]) && (r3->size[0] == F_ib_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = loop_ub;
      b_F_ij->size[1] = r3->size[0];
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      b_loop_ub = r3->size[0];
      for (i4 = 0; i4 < b_loop_ub; i4++) {
        for (i7 = 0; i7 < loop_ub; i7++) {
          b_F_data[i7 + b_F_ij->size[0] * i4] =
              F_ib_data[(i2 + i7) + F_ib->size[0] * ((int)r5[i4] - 1)] +
              F_ib__data[i7 + F_ib_->size[0] * i4];
        }
      }
      b_i = i6 - i5;
      loop_ub = r1->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          F_ib_data[(i5 + i4) + F_ib->size[0] * r4[i2]] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_22(F_ib, i5, i6, r1, i2, i4 - 1, r3, F_ib_);
      F_ib_data = F_ib->data;
    }
    /* %%% F_ijd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    dv2[0] = -0.0;
    dv2[3] = 0.0;
    dv2[6] = -0.0;
    dv2[1] = -0.0;
    dv2[4] = -0.0;
    dv2[7] = 1.0;
    dv2[2] = 0.0;
    dv2[5] = -1.0;
    dv2[8] = -0.0;
    for (i2 = 0; i2 < 3; i2++) {
      i4 = (int)dv2[i2];
      i5 = (int)dv2[i2 + 3];
      i6 = (int)dv2[i2 + 6];
      for (i7 = 0; i7 < 3; i7++) {
        LL[i2 + 3 * i7] = ((double)i4 * b_W.contents[3 * i7] +
                           (double)i5 * b_W.contents[3 * i7 + 1]) +
                          (double)i6 * b_W.contents[3 * i7 + 2];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      LL[i2] = (Q2_Wp_xipd[i2] + L_Wp_xipd[i2]) + LL[i2] * b_varThetadot;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = 0.0;
      d15 = Theta.contents[i2];
      d16 = Theta.contents[i2 + 3];
      d17 = Theta.contents[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        d14 += ((d15 * R0_.contents[3 * i4] + d16 * R0_.contents[3 * i4 + 1]) +
                d17 * R0_.contents[3 * i4 + 2]) *
               l[i4];
      }
      b_y[i2] = d14;
    }
    for (i2 = 0; i2 < 3; i2++) {
      d14 = 0.0;
      d15 = Q2_WTp_xipd[i2];
      d16 = Q2_WTp_xipd[i2 + 3];
      d17 = Q2_WTp_xipd[i2 + 6];
      d18 = Q2_WTthet_thetd[i2];
      d19 = Q2_WTthet_thetd[i2 + 3];
      d20 = Q2_WTthet_thetd[i2 + 6];
      d21 = b_W.contents[3 * i2];
      b_varThetaddot = b_W.contents[3 * i2 + 1];
      d_a = b_W.contents[3 * i2 + 2];
      for (i4 = 0; i4 < 3; i4++) {
        d14 += G1.contents[i4 + 3 * i2] * b_y[i4];
        c_a = L_lpd[3 * i4];
        b_M2_contents = d15 * c_a;
        b_Q1 = d18 * c_a;
        i5 = 3 * i4 + 1;
        c_a = L_lpd[i5];
        b_M2_contents += d16 * c_a;
        b_Q1 += d19 * c_a;
        i6 = 3 * i4 + 2;
        c_a = L_lpd[i6];
        b_M2_contents += d17 * c_a;
        b_Q1 += d20 * c_a;
        b_loop_ub = i2 + 3 * i4;
        W[b_loop_ub] =
            (d21 * II[3 * i4] + b_varThetaddot * II[i5]) + d_a * II[i6];
        b_I3[b_loop_ub] = b_Q1;
        b_M1_contents[b_loop_ub] = b_M2_contents + L_WTp_l[b_loop_ub];
      }
      b_xi0p[i2] = d14;
    }
    dv2[0] = 0.0;
    dv2[3] = -b_xi0p[2];
    dv2[6] = b_xi0p[1];
    dv2[1] = b_xi0p[2];
    dv2[4] = 0.0;
    dv2[7] = -b_xi0p[0];
    dv2[2] = -b_xi0p[1];
    dv2[5] = b_xi0p[0];
    dv2[8] = 0.0;
    d14 = b_xi0p[0];
    d15 = b_xi0p[1];
    d16 = b_xi0p[2];
    for (i2 = 0; i2 < 3; i2++) {
      d17 = l[i2];
      tau[3 * i2] = -b_d * d17;
      d18 = tau0_x_tau.contents[i2];
      xi0p[3 * i2] = d14 * d18 / (c_tau + 1.0);
      b_loop_ub = 3 * i2 + 1;
      tau[b_loop_ub] = -d10 * d17;
      xi0p[b_loop_ub] = d15 * d18 / (c_tau + 1.0);
      b_loop_ub = 3 * i2 + 2;
      tau[b_loop_ub] = -d11 * d17;
      xi0p[b_loop_ub] = d16 * d18 / (c_tau + 1.0);
    }
    for (i2 = 0; i2 < 9; i2++) {
      dv2[i2] += xi0p[i2];
    }
    for (i2 = 0; i2 < 3; i2++) {
      b_d = PP.contents[i2];
      d10 = PP.contents[i2 + 3];
      d11 = PP.contents[i2 + 6];
      d14 = tau[i2];
      d15 = tau[i2 + 3];
      d16 = tau[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        i5 = 3 * i4 + 1;
        i6 = 3 * i4 + 2;
        b_i = i2 + 3 * i4;
        e_a[b_i] = (y_tmp * b_d * dv2[3 * i4] + y_tmp * d10 * dv2[i5]) +
                   y_tmp * d11 * dv2[i6];
        b_htau0[b_i] = (d14 * b_W.contents[3 * i4] + d15 * b_W.contents[i5]) +
                       d16 * b_W.contents[i6];
      }
      b_d = W[i2];
      d10 = W[i2 + 3];
      d11 = W[i2 + 6];
      d14 = b_y_tmp[i2];
      d15 = b_y_tmp[i2 + 3];
      d16 = b_y_tmp[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        i5 = 3 * i4 + 1;
        i6 = 3 * i4 + 2;
        b_loop_ub = i2 + 3 * i4;
        tau[b_loop_ub] = ((b_M1_contents[b_loop_ub] + b_I3[b_loop_ub]) +
                          ((b_d * LL[3 * i4] + d10 * LL[i5]) + d11 * LL[i6])) -
                         y_tmp * (b_htau0[b_loop_ub] + e_a[b_loop_ub]);
        b_I3[b_loop_ub] =
            (d14 * L_lpd[3 * i4] + d15 * L_lpd[i5]) + d16 * L_lpd[i6];
      }
    }
    for (i2 = 0; i2 < 9; i2++) {
      tau[i2] -= b_I3[i2];
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i3; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i3; b_loop_ub_tmp++) {
        b_i = (int)(nel_ + (double)b_loop_ub_tmp) - 1;
        b_loop_ub = (int)(nel_ + (double)loop_ub_tmp) - 1;
        c_a = B_data[b_i] * B_data[b_loop_ub] * alph0 * rho;
        b_Q1 = Bp_data[b_i] * Bp_data[b_loop_ub];
        for (i2 = 0; i2 < 9; i2++) {
          htau0[i2] = c_a * I3.contents[i2] + b_Q1;
        }
        for (i2 = 0; i2 < 3; i2++) {
          b_d = htau0[i2];
          d10 = htau0[i2 + 3];
          d11 = htau0[i2 + 6];
          for (i4 = 0; i4 < 3; i4++) {
            b_htau0[i2 + 3 * i4] = (b_d * tau[3 * i4] + d10 * tau[3 * i4 + 1]) +
                                   d11 * tau[3 * i4 + 2];
          }
        }
        b_Q1 = (double)b_loop_ub_tmp * 3.0;
        b_varThetaddot = (double)loop_ub_tmp * 3.0;
        for (i2 = 0; i2 < 3; i2++) {
          b_loop_ub = (int)(b_varThetaddot + ((double)i2 + 1.0)) - 1;
          F_ijd__data[((int)(b_Q1 + 1.0) + F_ijd_->size[0] * b_loop_ub) - 1] =
              b_htau0[3 * i2] * nxi0p * wg_data[n];
          F_ijd__data[((int)(b_Q1 + 2.0) + F_ijd_->size[0] * b_loop_ub) - 1] =
              b_htau0[3 * i2 + 1] * nxi0p * wg_data[n];
          F_ijd__data[((int)(b_Q1 + 3.0) + F_ijd_->size[0] * b_loop_ub) - 1] =
              b_htau0[3 * i2 + 2] * nxi0p * wg_data[n];
        }
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
      i7 = 0;
      b_i = 0;
      b_loop_ub_tmp = 0;
      loop_ub_tmp = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
      i7 = (int)d13 - 1;
      b_i = (int)d12;
      b_loop_ub_tmp = (int)d13 - 1;
      loop_ub_tmp = (int)d12;
    }
    loop_ub = i4 - i2;
    b_loop_ub = i6 - i5;
    if ((loop_ub == F_ijd_->size[0]) && (b_loop_ub == F_ijd_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = loop_ub;
      b_F_ij->size[1] = b_loop_ub;
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      for (i4 = 0; i4 < b_loop_ub; i4++) {
        for (i6 = 0; i6 < loop_ub; i6++) {
          b_F_data[i6 + b_F_ij->size[0] * i4] =
              F_ijd_data[(i2 + i6) + F_ijd->size[0] * (i5 + i4)] +
              F_ijd__data[i6 + F_ijd_->size[0] * i4];
        }
      }
      b_i -= i7;
      b_loop_ub = loop_ub_tmp - b_loop_ub_tmp;
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          F_ijd_data[(i7 + i4) + F_ijd->size[0] * (b_loop_ub_tmp + i2)] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_9(F_ijd, i7, b_i, b_loop_ub_tmp, loop_ub_tmp, i2, i4 - 1,
                         i5, i6 - 1, F_ijd_);
      F_ijd_data = F_ijd->data;
    }
    /* %%% F_ibd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    b_d = L_lthet[0];
    d10 = L_lthet[1];
    d11 = L_lthet[2];
    for (i2 = 0; i2 < 3; i2++) {
      Q1[i2] =
          (((((Q2_WTp_xipd[i2] * d1 + Q2_WTp_xipd[i2 + 3] * d2) +
              Q2_WTp_xipd[i2 + 6] * d3) +
             ((Q2_WTthet_thetd[i2] * d1 + Q2_WTthet_thetd[i2 + 3] * d2) +
              Q2_WTthet_thetd[i2 + 6] * d3)) +
            L_WTthet_l[i2]) +
           ((b_W.contents[3 * i2] * b_d + b_W.contents[3 * i2 + 1] * d10) +
            b_W.contents[3 * i2 + 2] * d11)) -
          ((b_y_tmp[i2] * d1 + b_y_tmp[i2 + 3] * d2) + b_y_tmp[i2 + 6] * d3);
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i1; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i3; b_loop_ub_tmp++) {
        c_a = Bp_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
              Bbrev_data[(int)(nel_ + (double)loop_ub_tmp) - 1];
        b_Q1 = (double)b_loop_ub_tmp * 3.0;
        F_ibd__data[((int)(b_Q1 + 1.0) + F_ibd_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[0] * nxi0p * wg_data[n];
        F_ibd__data[((int)(b_Q1 + 2.0) + F_ibd_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[1] * nxi0p * wg_data[n];
        F_ibd__data[((int)(b_Q1 + 3.0) + F_ibd_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[2] * nxi0p * wg_data[n];
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
    }
    i7 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i7);
    r4 = r1->data;
    loop_ub = r3->size[0];
    for (i7 = 0; i7 < loop_ub; i7++) {
      r4[i7] = (int)r5[i7] - 1;
    }
    loop_ub = i4 - i2;
    if ((loop_ub == F_ibd_->size[0]) && (r3->size[0] == F_ibd_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = loop_ub;
      b_F_ij->size[1] = r3->size[0];
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      b_loop_ub = r3->size[0];
      for (i4 = 0; i4 < b_loop_ub; i4++) {
        for (i7 = 0; i7 < loop_ub; i7++) {
          b_F_data[i7 + b_F_ij->size[0] * i4] =
              F_ibd_data[(i2 + i7) + F_ibd->size[0] * ((int)r5[i4] - 1)] +
              F_ibd__data[i7 + F_ibd_->size[0] * i4];
        }
      }
      b_i = i6 - i5;
      loop_ub = r1->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          F_ibd_data[(i5 + i4) + F_ibd->size[0] * r4[i2]] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_22(F_ibd, i5, i6, r1, i2, i4 - 1, r3, F_ibd_);
      F_ibd_data = F_ibd->data;
    }
    /* %%% F_ijdd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    for (i2 = 0; i2 < 3; i2++) {
      b_d = b_W.contents[3 * i2];
      d10 = b_W.contents[3 * i2 + 1];
      d11 = b_W.contents[3 * i2 + 2];
      for (i4 = 0; i4 < 3; i4++) {
        W[i2 + 3 * i4] =
            (b_d * II[3 * i4] + d10 * II[3 * i4 + 1]) + d11 * II[3 * i4 + 2];
      }
      b_d = W[i2];
      d10 = W[i2 + 3];
      d11 = W[i2 + 6];
      for (i4 = 0; i4 < 3; i4++) {
        tau[i2 + 3 * i4] =
            (b_d * b_W.contents[3 * i4] + d10 * b_W.contents[3 * i4 + 1]) +
            d11 * b_W.contents[3 * i4 + 2];
      }
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i3; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i3; b_loop_ub_tmp++) {
        b_i = (int)(nel_ + (double)b_loop_ub_tmp) - 1;
        b_loop_ub = (int)(nel_ + (double)loop_ub_tmp) - 1;
        c_a = B_data[b_i] * B_data[b_loop_ub] * rho;
        d_a = Bp_data[b_i] * Bp_data[b_loop_ub];
        b_Q1 = (double)b_loop_ub_tmp * 3.0;
        b_varThetaddot = (double)loop_ub_tmp * 3.0;
        for (i2 = 0; i2 < 3; i2++) {
          b_loop_ub = (int)(b_varThetaddot + ((double)i2 + 1.0)) - 1;
          F_ijdd__data[((int)(b_Q1 + 1.0) + F_ijdd_->size[0] * b_loop_ub) - 1] =
              (c_a * I3.contents[3 * i2] + d_a * tau[3 * i2]) * nxi0p *
              wg_data[n];
          b_i = 3 * i2 + 1;
          F_ijdd__data[((int)(b_Q1 + 2.0) + F_ijdd_->size[0] * b_loop_ub) - 1] =
              (c_a * I3.contents[b_i] + d_a * tau[b_i]) * nxi0p * wg_data[n];
          b_i = 3 * i2 + 2;
          F_ijdd__data[((int)(b_Q1 + 3.0) + F_ijdd_->size[0] * b_loop_ub) - 1] =
              (c_a * I3.contents[b_i] + d_a * tau[b_i]) * nxi0p * wg_data[n];
        }
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
      i7 = 0;
      b_i = 0;
      b_loop_ub_tmp = 0;
      loop_ub_tmp = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
      i7 = (int)d13 - 1;
      b_i = (int)d12;
      b_loop_ub_tmp = (int)d13 - 1;
      loop_ub_tmp = (int)d12;
    }
    loop_ub = i4 - i2;
    b_loop_ub = i6 - i5;
    if ((loop_ub == F_ijdd_->size[0]) && (b_loop_ub == F_ijdd_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = loop_ub;
      b_F_ij->size[1] = b_loop_ub;
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      for (i4 = 0; i4 < b_loop_ub; i4++) {
        for (i6 = 0; i6 < loop_ub; i6++) {
          b_F_data[i6 + b_F_ij->size[0] * i4] =
              F_ijdd_data[(i2 + i6) + F_ijdd->size[0] * (i5 + i4)] +
              F_ijdd__data[i6 + F_ijdd_->size[0] * i4];
        }
      }
      b_i -= i7;
      b_loop_ub = loop_ub_tmp - b_loop_ub_tmp;
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          F_ijdd_data[(i7 + i4) + F_ijdd->size[0] * (b_loop_ub_tmp + i2)] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_9(F_ijdd, i7, b_i, b_loop_ub_tmp, loop_ub_tmp, i2,
                         i4 - 1, i5, i6 - 1, F_ijdd_);
      F_ijdd_data = F_ijdd->data;
    }
    /* %%% F_ibdd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    for (i2 = 0; i2 < 3; i2++) {
      Q1[i2] = (b_W.contents[3 * i2] * d4 + b_W.contents[3 * i2 + 1] * d5) +
               b_W.contents[3 * i2 + 2] * d6;
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i1; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i3; b_loop_ub_tmp++) {
        c_a = Bp_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
              Bbrev_data[(int)(nel_ + (double)loop_ub_tmp) - 1];
        b_Q1 = (double)b_loop_ub_tmp * 3.0;
        F_ibdd__data[((int)(b_Q1 + 1.0) + F_ibdd_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[0] * nxi0p * wg_data[n];
        F_ibdd__data[((int)(b_Q1 + 2.0) + F_ibdd_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[1] * nxi0p * wg_data[n];
        F_ibdd__data[((int)(b_Q1 + 3.0) + F_ibdd_->size[0] * loop_ub_tmp) - 1] =
            c_a * Q1[2] * nxi0p * wg_data[n];
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
    }
    i7 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i7);
    r4 = r1->data;
    loop_ub = r3->size[0];
    for (i7 = 0; i7 < loop_ub; i7++) {
      r4[i7] = (int)r5[i7] - 1;
    }
    loop_ub = i4 - i2;
    if ((loop_ub == F_ibdd_->size[0]) && (r3->size[0] == F_ibdd_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = loop_ub;
      b_F_ij->size[1] = r3->size[0];
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      b_loop_ub = r3->size[0];
      for (i4 = 0; i4 < b_loop_ub; i4++) {
        for (i7 = 0; i7 < loop_ub; i7++) {
          b_F_data[i7 + b_F_ij->size[0] * i4] =
              F_ibdd_data[(i2 + i7) + F_ibdd->size[0] * ((int)r5[i4] - 1)] +
              F_ibdd__data[i7 + F_ibdd_->size[0] * i4];
        }
      }
      b_i = i6 - i5;
      loop_ub = r1->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          F_ibdd_data[(i5 + i4) + F_ibdd->size[0] * r4[i2]] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_22(F_ibdd, i5, i6, r1, i2, i4 - 1, r3, F_ibdd_);
      F_ibdd_data = F_ibdd->data;
    }
    /* %%% mu_aj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    L_Wp(&tau0_x_tau, &tau_dot_tau0, &I3, &nxip, &G1, &PP, &b_W, &b_tau, &R0_,
         &Theta, xipdd, b);
    Q1T_Wthetp(&tau0_x_tau, &tau_dot_tau0, &I3, &nxip, &G1, &PP, &b_W, &b_tau,
               &R0_, &Theta, &e1, xipd, l, F1);
    for (i2 = 0; i2 < 3; i2++) {
      b_e1[i2] = II[3 * i2];
    }
    b_d = L_Wthet_xipd[1];
    d10 = L_Wthet_xipd[2];
    d11 = b_e1[0];
    d14 = b_e1[1];
    d15 = b_e1[2];
    for (i2 = 0; i2 < 3; i2++) {
      i4 = 3 * i2 + 1;
      i5 = 3 * i2 + 2;
      c_W[i2] = b_d * L_lp[i4] + d10 * L_lp[i5];
      c_e1[i2] = ((d11 * Wd_xipd_p[3 * i2] + d14 * Wd_xipd_p[i4]) +
                  d15 * Wd_xipd_p[i5]) -
                 F1[i2];
    }
    for (i2 = 0; i2 < 3; i2++) {
      b_e1[i2] = II[3 * i2];
    }
    b_d = b_e1[0];
    d10 = b_e1[1];
    d11 = b_e1[2];
    for (i2 = 0; i2 < 3; i2++) {
      Q1[i2] = (c_e1[i2] - c_W[i2]) +
               ((b_d * b[3 * i2] + d10 * b[3 * i2 + 1]) + d11 * b[3 * i2 + 2]);
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i3; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i1; b_loop_ub_tmp++) {
        c_a = Bbrev_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
              Bp_data[(int)(nel_ + (double)loop_ub_tmp) - 1];
        b_varThetaddot = (double)loop_ub_tmp * 3.0;
        mu_aj__data[b_loop_ub_tmp +
                    mu_aj_->size[0] * ((int)(b_varThetaddot + 1.0) - 1)] =
            c_a * Q1[0] * nxi0p * wg_data[n];
        mu_aj__data[b_loop_ub_tmp +
                    mu_aj_->size[0] * ((int)(b_varThetaddot + 2.0) - 1)] =
            c_a * Q1[1] * nxi0p * wg_data[n];
        mu_aj__data[b_loop_ub_tmp +
                    mu_aj_->size[0] * ((int)(b_varThetaddot + 3.0) - 1)] =
            c_a * Q1[2] * nxi0p * wg_data[n];
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
    }
    i7 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i7);
    r4 = r1->data;
    loop_ub = r3->size[0];
    for (i7 = 0; i7 < loop_ub; i7++) {
      r4[i7] = (int)r5[i7] - 1;
    }
    loop_ub = i4 - i2;
    if ((r3->size[0] == mu_aj_->size[0]) && (loop_ub == mu_aj_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = r3->size[0];
      b_F_ij->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      for (i4 = 0; i4 < loop_ub; i4++) {
        b_loop_ub = r3->size[0];
        for (i7 = 0; i7 < b_loop_ub; i7++) {
          b_F_data[i7 + b_F_ij->size[0] * i4] =
              mu_aj_data[((int)r5[i7] + mu_aj->size[0] * (i2 + i4)) - 1] +
              mu_aj__data[i7 + mu_aj_->size[0] * i4];
        }
      }
      b_i = r1->size[0];
      b_loop_ub = i6 - i5;
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          mu_aj_data[r4[i4] + mu_aj->size[0] * (i5 + i2)] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_17(mu_aj, r1, i5, i6, r3, i2, i4 - 1, mu_aj_);
      mu_aj_data = mu_aj->data;
    }
    /* %%% mu_ab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    d_a = 0.0;
    b_d = xipdd[0];
    d10 = xipdd[1];
    d11 = xipdd[2];
    d14 = xipd[0];
    d15 = xipd[1];
    d16 = xipd[2];
    for (i2 = 0; i2 < 3; i2++) {
      d17 = b_W.contents[i2];
      d18 = d17 * b_d;
      d19 = d17 * d14;
      d17 = b_W.contents[i2 + 3];
      d18 += d17 * d10;
      d19 += d17 * d15;
      d17 = b_W.contents[i2 + 6];
      d18 += d17 * d11;
      d19 += d17 * d16;
      d_y[i2] = d19;
      b_y[i2] = d18;
      d_a += II[3 * i2] * Wd_xipd_thet[i2];
    }
    b_e1[0] = 0.0;
    b_e1[1] = 0.0 - l[2];
    b_e1[2] = l[1];
    Q1[0] = -0.0;
    Q1[1] = -(0.0 - d_y[2]);
    Q1[2] = -d_y[1];
    b_Q1 = 0.0;
    b_varThetaddot = 0.0;
    for (i2 = 0; i2 < 3; i2++) {
      b_Q1 += b_e1[i2] * Q1[i2];
      b_varThetaddot += L_Wthet_xipd[i2] * L_lthet[i2];
      b_e1[i2] = II[3 * i2];
    }
    b_Q1 = ((d_a - b_Q1) - b_varThetaddot) +
           (b_e1[1] * -(0.0 - b_y[2]) + b_e1[2] * -b_y[1]);
    for (loop_ub_tmp = 0; loop_ub_tmp < i1; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i1; b_loop_ub_tmp++) {
        mu_ab__data[b_loop_ub_tmp + mu_ab_->size[0] * loop_ub_tmp] =
            Bbrev_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
            Bbrev_data[(int)(nel_ + (double)loop_ub_tmp) - 1] * b_Q1 * nxi0p *
            wg_data[n];
      }
    }
    i2 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i2);
    r4 = r1->data;
    loop_ub = r3->size[0];
    i2 = r2->size[0];
    r2->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r2, i2);
    r6 = r2->data;
    for (i2 = 0; i2 < loop_ub; i2++) {
      i4 = (int)r5[i2] - 1;
      r4[i2] = i4;
      r6[i2] = i4;
    }
    if ((r3->size[0] == mu_ab_->size[0]) && (r3->size[0] == mu_ab_->size[1])) {
      i2 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = r3->size[0];
      b_F_ij->size[1] = r3->size[0];
      emxEnsureCapacity_real_T(b_F_ij, i2);
      b_F_data = b_F_ij->data;
      loop_ub = r3->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_loop_ub = r3->size[0];
        for (i4 = 0; i4 < b_loop_ub; i4++) {
          b_F_data[i4 + b_F_ij->size[0] * i2] =
              mu_ab_data[((int)r5[i4] + mu_ab->size[0] * ((int)r5[i2] - 1)) -
                         1] +
              mu_ab__data[i4 + mu_ab_->size[0] * i2];
        }
      }
      b_i = r1->size[0];
      loop_ub = r2->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          mu_ab_data[r4[i4] + mu_ab->size[0] * r6[i2]] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_16(mu_ab, r1, r2, r3, mu_ab_);
      mu_ab_data = mu_ab->data;
    }
    /* %%% mu_ajd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    for (i2 = 0; i2 < 3; i2++) {
      b_e1[i2] = II[3 * i2];
    }
    b_d = b_e1[0];
    d10 = b_e1[1];
    d11 = b_e1[2];
    for (i2 = 0; i2 < 3; i2++) {
      c_W[i2] =
          (b_d * LL[3 * i2] + d10 * LL[3 * i2 + 1]) + d11 * LL[3 * i2 + 2];
    }
    b_d = 0.0 - l[2];
    d10 = l[1];
    for (i2 = 0; i2 < 3; i2++) {
      b_e1[i2] =
          b_d * b_W.contents[3 * i2 + 1] + d10 * b_W.contents[3 * i2 + 2];
    }
    b_d = L_Wthet_xipd[1];
    d10 = L_Wthet_xipd[2];
    for (i2 = 0; i2 < 3; i2++) {
      Q1[i2] = (c_W[i2] - b_e1[i2]) -
               (b_d * L_lpd[3 * i2 + 1] + d10 * L_lpd[3 * i2 + 2]);
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i3; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i1; b_loop_ub_tmp++) {
        c_a = Bbrev_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
              Bp_data[(int)(nel_ + (double)loop_ub_tmp) - 1];
        b_varThetaddot = (double)loop_ub_tmp * 3.0;
        mu_ajd__data[b_loop_ub_tmp +
                     mu_ajd_->size[0] * ((int)(b_varThetaddot + 1.0) - 1)] =
            c_a * Q1[0] * nxi0p * wg_data[n];
        mu_ajd__data[b_loop_ub_tmp +
                     mu_ajd_->size[0] * ((int)(b_varThetaddot + 2.0) - 1)] =
            c_a * Q1[1] * nxi0p * wg_data[n];
        mu_ajd__data[b_loop_ub_tmp +
                     mu_ajd_->size[0] * ((int)(b_varThetaddot + 3.0) - 1)] =
            c_a * Q1[2] * nxi0p * wg_data[n];
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
    }
    i7 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i7);
    r4 = r1->data;
    loop_ub = r3->size[0];
    for (i7 = 0; i7 < loop_ub; i7++) {
      r4[i7] = (int)r5[i7] - 1;
    }
    loop_ub = i4 - i2;
    if ((r3->size[0] == mu_ajd_->size[0]) && (loop_ub == mu_ajd_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = r3->size[0];
      b_F_ij->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      for (i4 = 0; i4 < loop_ub; i4++) {
        b_loop_ub = r3->size[0];
        for (i7 = 0; i7 < b_loop_ub; i7++) {
          b_F_data[i7 + b_F_ij->size[0] * i4] =
              mu_ajd_data[((int)r5[i7] + mu_ajd->size[0] * (i2 + i4)) - 1] +
              mu_ajd__data[i7 + mu_ajd_->size[0] * i4];
        }
      }
      b_i = r1->size[0];
      b_loop_ub = i6 - i5;
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          mu_ajd_data[r4[i4] + mu_ajd->size[0] * (i5 + i2)] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_17(mu_ajd, r1, i5, i6, r3, i2, i4 - 1, mu_ajd_);
      mu_ajd_data = mu_ajd->data;
    }
    /* %%% mu_abd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    b_Q1 = L_lthet[0] - (L_Wthet_xipd[1] * II[1] + L_Wthet_xipd[2] * II[2]);
    for (loop_ub_tmp = 0; loop_ub_tmp < i1; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i1; b_loop_ub_tmp++) {
        mu_abd__data[b_loop_ub_tmp + mu_abd_->size[0] * loop_ub_tmp] =
            Bbrev_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
            Bbrev_data[(int)(nel_ + (double)loop_ub_tmp) - 1] * b_Q1 * nxi0p *
            wg_data[n];
      }
    }
    i2 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i2);
    r4 = r1->data;
    loop_ub = r3->size[0];
    i2 = r2->size[0];
    r2->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r2, i2);
    r6 = r2->data;
    for (i2 = 0; i2 < loop_ub; i2++) {
      i4 = (int)r5[i2] - 1;
      r4[i2] = i4;
      r6[i2] = i4;
    }
    if ((r3->size[0] == mu_abd_->size[0]) &&
        (r3->size[0] == mu_abd_->size[1])) {
      i2 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = r3->size[0];
      b_F_ij->size[1] = r3->size[0];
      emxEnsureCapacity_real_T(b_F_ij, i2);
      b_F_data = b_F_ij->data;
      loop_ub = r3->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_loop_ub = r3->size[0];
        for (i4 = 0; i4 < b_loop_ub; i4++) {
          b_F_data[i4 + b_F_ij->size[0] * i2] =
              mu_abd_data[((int)r5[i4] + mu_abd->size[0] * ((int)r5[i2] - 1)) -
                          1] +
              mu_abd__data[i4 + mu_abd_->size[0] * i2];
        }
      }
      b_i = r1->size[0];
      loop_ub = r2->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          mu_abd_data[r4[i4] + mu_abd->size[0] * r6[i2]] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_16(mu_abd, r1, r2, r3, mu_abd_);
      mu_abd_data = mu_abd->data;
    }
    /* %%% mu_ajdd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    for (i2 = 0; i2 < 3; i2++) {
      Q1[i2] = (d7 * b_W.contents[3 * i2] + d8 * b_W.contents[3 * i2 + 1]) +
               d9 * b_W.contents[3 * i2 + 2];
    }
    for (loop_ub_tmp = 0; loop_ub_tmp < i3; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i1; b_loop_ub_tmp++) {
        c_a = Bbrev_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
              Bp_data[(int)(nel_ + (double)loop_ub_tmp) - 1];
        b_varThetaddot = (double)loop_ub_tmp * 3.0;
        mu_ajdd__data[b_loop_ub_tmp +
                      mu_ajdd_->size[0] * ((int)(b_varThetaddot + 1.0) - 1)] =
            c_a * Q1[0] * nxi0p * wg_data[n];
        mu_ajdd__data[b_loop_ub_tmp +
                      mu_ajdd_->size[0] * ((int)(b_varThetaddot + 2.0) - 1)] =
            c_a * Q1[1] * nxi0p * wg_data[n];
        mu_ajdd__data[b_loop_ub_tmp +
                      mu_ajdd_->size[0] * ((int)(b_varThetaddot + 3.0) - 1)] =
            c_a * Q1[2] * nxi0p * wg_data[n];
      }
    }
    if (d13 > d12) {
      i2 = 0;
      i4 = 0;
      i5 = 0;
      i6 = 0;
    } else {
      i2 = (int)d13 - 1;
      i4 = (int)d12;
      i5 = (int)d13 - 1;
      i6 = (int)d12;
    }
    i7 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i7);
    r4 = r1->data;
    loop_ub = r3->size[0];
    for (i7 = 0; i7 < loop_ub; i7++) {
      r4[i7] = (int)r5[i7] - 1;
    }
    loop_ub = i4 - i2;
    if ((r3->size[0] == mu_ajdd_->size[0]) && (loop_ub == mu_ajdd_->size[1])) {
      i4 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = r3->size[0];
      b_F_ij->size[1] = loop_ub;
      emxEnsureCapacity_real_T(b_F_ij, i4);
      b_F_data = b_F_ij->data;
      for (i4 = 0; i4 < loop_ub; i4++) {
        b_loop_ub = r3->size[0];
        for (i7 = 0; i7 < b_loop_ub; i7++) {
          b_F_data[i7 + b_F_ij->size[0] * i4] =
              mu_ajdd_data[((int)r5[i7] + mu_ajdd->size[0] * (i2 + i4)) - 1] +
              mu_ajdd__data[i7 + mu_ajdd_->size[0] * i4];
        }
      }
      b_i = r1->size[0];
      b_loop_ub = i6 - i5;
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          mu_ajdd_data[r4[i4] + mu_ajdd->size[0] * (i5 + i2)] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_17(mu_ajdd, r1, i5, i6, r3, i2, i4 - 1, mu_ajdd_);
      mu_ajdd_data = mu_ajdd->data;
    }
    /* %%% mu_abdd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
    for (loop_ub_tmp = 0; loop_ub_tmp < i1; loop_ub_tmp++) {
      for (b_loop_ub_tmp = 0; b_loop_ub_tmp < i1; b_loop_ub_tmp++) {
        mu_abdd__data[b_loop_ub_tmp + mu_abdd_->size[0] * loop_ub_tmp] =
            Bbrev_data[(int)(nel_ + (double)b_loop_ub_tmp) - 1] *
            Bbrev_data[(int)(nel_ + (double)loop_ub_tmp) - 1] * II[0] * nxi0p *
            wg_data[n];
      }
    }
    i2 = r1->size[0];
    r1->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r1, i2);
    r4 = r1->data;
    loop_ub = r3->size[0];
    i2 = r2->size[0];
    r2->size[0] = r3->size[0];
    emxEnsureCapacity_int32_T(r2, i2);
    r6 = r2->data;
    for (i2 = 0; i2 < loop_ub; i2++) {
      i4 = (int)r5[i2] - 1;
      r4[i2] = i4;
      r6[i2] = i4;
    }
    if ((r3->size[0] == mu_abdd_->size[0]) &&
        (r3->size[0] == mu_abdd_->size[1])) {
      i2 = b_F_ij->size[0] * b_F_ij->size[1];
      b_F_ij->size[0] = r3->size[0];
      b_F_ij->size[1] = r3->size[0];
      emxEnsureCapacity_real_T(b_F_ij, i2);
      b_F_data = b_F_ij->data;
      loop_ub = r3->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        b_loop_ub = r3->size[0];
        for (i4 = 0; i4 < b_loop_ub; i4++) {
          b_F_data[i4 + b_F_ij->size[0] * i2] =
              mu_abdd_data[((int)r5[i4] +
                            mu_abdd->size[0] * ((int)r5[i2] - 1)) -
                           1] +
              mu_abdd__data[i4 + mu_abdd_->size[0] * i2];
        }
      }
      b_i = r1->size[0];
      loop_ub = r2->size[0];
      for (i2 = 0; i2 < loop_ub; i2++) {
        for (i4 = 0; i4 < b_i; i4++) {
          mu_abdd_data[r4[i4] + mu_abdd->size[0] * r6[i2]] =
              b_F_data[i4 + b_i * i2];
        }
      }
    } else {
      binary_expand_op_16(mu_abdd, r1, r2, r3, mu_abdd_);
      mu_abdd_data = mu_abdd->data;
    }
  }
  emxFree_real_T(&b_F_ij);
  emxFree_real_T(&b_mu);
  emxFree_real_T(&b_F);
  emxFree_uint32_T(&r3);
  emxFree_real_T(&y);
  emxFree_int32_T(&r2);
  emxFree_int32_T(&r1);
  emxFree_int32_T(&r);
  emxFree_real_T(&Bbrev);
  emxFree_real_T(&Bp);
  emxFree_real_T(&B);
  emxFree_real_T(&mu_abdd_);
  emxFree_real_T(&mu_ajdd_);
  emxFree_real_T(&mu_abd_);
  emxFree_real_T(&mu_ajd_);
  emxFree_real_T(&mu_ab_);
  emxFree_real_T(&mu_aj_);
  emxFree_real_T(&F_ibdd_);
  emxFree_real_T(&F_ijdd_);
  emxFree_real_T(&F_ibd_);
  emxFree_real_T(&F_ijd_);
  emxFree_real_T(&F_ib_);
  emxFree_real_T(&F_ij_);
  emxFree_real_T(&F_);
}

/* End of code generation (CableInertiaForce.c) */
