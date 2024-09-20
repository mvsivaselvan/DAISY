/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableForce.c
 *
 * Code generation for function 'CableForce'
 *
 */

/* Include files */
#include "CableForce.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_internal_types.h"
#include "CableForceRotBCinCoord_types.h"
#include "diag.h"
#include "mldivide.h"
#include "mtimes.h"
#include "norm.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void Q_rr_I(const c_captured_var *L_KI, const captured_var *K,
                   const b_captured_var *nxip, const c_captured_var *I3,
                   const captured_var *tau, const c_captured_var *PP,
                   const captured_var *tau0, const c_captured_var *L_tauI,
                   const captured_var *K0, const c_captured_var *L_GI,
                   const captured_var *G, const double eta[3], double Q[9]);

static void Q_rr_III(const b_captured_var *nxip, const c_captured_var *L_KII,
                     const captured_var *tau, const captured_var *tau0,
                     const c_captured_var *L_tauI, const c_captured_var *L_GII,
                     const double eta[3], double Q[9]);

static void Qtil_rr_I(const b_captured_var *nxip, const captured_var *K,
                      const c_captured_var *L_KI, const captured_var *tau,
                      const c_captured_var *I3, const c_captured_var *PP,
                      const captured_var *tau0, const c_captured_var *L_tauI,
                      const captured_var *K0, const captured_var *G,
                      const c_captured_var *L_GI, const double nu[3],
                      double Q[9]);

static void Qtil_rr_III(const b_captured_var *nxip, const captured_var *tau,
                        const c_captured_var *L_KII,
                        const c_captured_var *L_tauI, const captured_var *tau0,
                        const c_captured_var *L_GII, const double nu[3],
                        double Q[9]);

static void Qtil_rr_IV(const b_captured_var *nxip, const c_captured_var *L_KII,
                       const captured_var *tau, const captured_var *tau0,
                       const c_captured_var *L_tauI,
                       const c_captured_var *L_GII, const double nu[3],
                       double Q[9]);

static void binary_expand_op_15(emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3);

static void binary_expand_op_4(emxArray_real_T *in1,
                               const emxArray_int32_T *in2,
                               const emxArray_real_T *in3,
                               const emxArray_real_T *in4, double in5,
                               double in6, const emxArray_real_T *in7, int in8);

static void binary_expand_op_5(emxArray_real_T *in1,
                               const emxArray_int32_T *in2,
                               const emxArray_int32_T *in3,
                               const emxArray_real_T *in4,
                               const emxArray_real_T *in5);

static void binary_expand_op_7(emxArray_real_T *in1,
                               const emxArray_int32_T *in2, int in3, int in4,
                               const emxArray_real_T *in5, int in6, int in7,
                               const emxArray_real_T *in8);

static void binary_expand_op_8(emxArray_real_T *in1, int in2, int in3,
                               const emxArray_int32_T *in4, int in5, int in6,
                               const emxArray_real_T *in7,
                               const emxArray_real_T *in8);

static void plus(emxArray_real_T *in1, const emxArray_real_T *in2);

/* Function Definitions */
static void Q_rr_I(const c_captured_var *L_KI, const captured_var *K,
                   const b_captured_var *nxip, const c_captured_var *I3,
                   const captured_var *tau, const c_captured_var *PP,
                   const captured_var *tau0, const c_captured_var *L_tauI,
                   const captured_var *K0, const c_captured_var *L_GI,
                   const captured_var *G, const double eta[3], double Q[9])
{
  double b_LGIeta[9];
  double b_LKItau0[9];
  double b_LtauIeta[9];
  double b_LtauItau0[9];
  double c_LKIeta[9];
  double c_LtauItau0[9];
  double d_LKIeta[9];
  double d_Peta[9];
  double e_Peta[9];
  double e_a[9];
  double e_tau[9];
  double f_Peta[9];
  double f_tau[9];
  double g_tau[9];
  double LGIeta[3];
  double LKIeta[3];
  double LKItau0[3];
  double LtauIeta[3];
  double LtauItau0[3];
  double Peta[3];
  double b_LKIeta[3];
  double b_Peta[3];
  double c_Peta[3];
  double a;
  double a_tmp;
  double b_a;
  double b_eta;
  double b_tau;
  double b_tau0;
  double b_tau0_tmp;
  double b_tmp;
  double c_a;
  double c_eta;
  double c_tau;
  double c_tau0;
  double c_tmp;
  double d;
  double d_a;
  double d_eta;
  double d_tau;
  double e_eta;
  double tau0_tmp;
  int i;
  /*  functions related to rr */
  tau0_tmp = 0.0;
  b_tau0 = 0.0;
  b_eta = 0.0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  UTILITY FUNCTIONS FOR SECOND DERIVATIVE COMPUTATIONS */
  /*  It is assumed in these functions that quantities xip, xipp,  */
  /*  their norms (nxip, nxipp), tau0, tau, PP, L_tauI, */
  /*  K0, K, L_KI, L_KII, G, L_GI, L_GII have been computed and are available */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  function related to tau */
  b_tau = 0.0;
  a_tmp = nxip->contents;
  c_tmp = a_tmp * a_tmp;
  /*  functions related to K */
  c_eta = 0.0;
  b_tmp = 2.0 / nxip->contents;
  /*  functions related to G */
  d_eta = 0.0;
  /*  functions related to K */
  b_tau0_tmp = 0.0;
  c_tau0 = 0.0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  UTILITY FUNCTIONS FOR SECOND DERIVATIVE COMPUTATIONS */
  /*  It is assumed in these functions that quantities xip, xipp,  */
  /*  their norms (nxip, nxipp), tau0, tau, PP, L_tauI, */
  /*  K0, K, L_KI, L_KII, G, L_GI, L_GII have been computed and are available */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  function related to tau */
  c_tau = 0.0;
  e_eta = 0.0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  UTILITY FUNCTIONS FOR SECOND DERIVATIVE COMPUTATIONS */
  /*  It is assumed in these functions that quantities xip, xipp,  */
  /*  their norms (nxip, nxipp), tau0, tau, PP, L_tauI, */
  /*  K0, K, L_KI, L_KII, G, L_GI, L_GII have been computed and are available */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  function related to tau */
  d_tau = 0.0;
  for (i = 0; i < 3; i++) {
    double d1;
    double d2;
    double d3;
    double d4;
    d = tau->contents[i];
    tau0_tmp += tau0->contents[i] * d;
    b_tau0 = tau0_tmp;
    d_a = L_tauI->contents[i];
    c_a = L_tauI->contents[i + 3];
    b_a = L_tauI->contents[i + 6];
    LtauItau0[i] = (d_a * tau0->contents[0] + c_a * tau0->contents[1]) +
                   b_a * tau0->contents[2];
    LGIeta[i] = (L_GI->contents[i] * eta[0] + L_GI->contents[i + 3] * eta[1]) +
                L_GI->contents[i + 6] * eta[2];
    b_eta += eta[i] * G->contents[i];
    a_tmp = PP->contents[i];
    a = PP->contents[i + 3];
    d1 = PP->contents[i + 6];
    Peta[i] = (a_tmp * tau0->contents[0] + a * tau0->contents[1]) +
              d1 * tau0->contents[2];
    d2 = tau0->contents[i];
    b_tau += d * d2;
    d3 = K->contents[i];
    c_eta += eta[i] * d3;
    d_eta += eta[i] * (2.0 * d2 + d);
    b_tau0_tmp += d2 * d3;
    c_tau0 = b_tau0_tmp;
    c_tau += d * eta[i];
    e_eta += eta[i] * d2;
    d2 = L_KI->contents[i];
    d3 = d2 * eta[0];
    d4 = d2 * tau0->contents[0];
    d2 = L_KI->contents[i + 3];
    d3 += d2 * eta[1];
    d4 += d2 * tau0->contents[1];
    d2 = L_KI->contents[i + 6];
    d3 += d2 * eta[2];
    d2 = d4 + d2 * tau0->contents[2];
    c_Peta[i] =
        (a_tmp * K0->contents[0] + a * K0->contents[1]) + d1 * K0->contents[2];
    b_Peta[i] = (a_tmp * eta[0] + a * eta[1]) + d1 * eta[2];
    b_LKIeta[i] = d2;
    LtauIeta[i] = (d_a * eta[0] + c_a * eta[1]) + b_a * eta[2];
    LKItau0[i] = d2;
    LKIeta[i] = d3;
    d_tau += d * K0->contents[i];
  }
  a = 2.0 * b_eta / pow(b_tau0 + 1.0, 3.0);
  a_tmp = (b_tau0 + 1.0) * (b_tau0 + 1.0);
  b_a = b_eta / a_tmp;
  c_a = c_eta / nxip->contents;
  d_a = b_tau0_tmp / nxip->contents;
  for (i = 0; i < 3; i++) {
    int i1;
    d = I3->contents[3 * i];
    Q[3 * i] = c_a * d + tau->contents[0] * LKIeta[i];
    c_LKIeta[3 * i] = LKIeta[0] * tau->contents[i];
    e_a[3 * i] = d_a * d + tau->contents[0] * b_LKIeta[i];
    d_LKIeta[3 * i] = b_LKIeta[0] * tau->contents[i];
    b_LKItau0[3 * i] = LKItau0[0] * LtauIeta[i];
    b_LtauIeta[3 * i] = LtauIeta[0] * LKItau0[i];
    e_tau[3 * i] = tau->contents[0] * b_Peta[i];
    d_Peta[3 * i] = b_Peta[0] * tau->contents[i];
    f_tau[3 * i] = tau->contents[0] * c_Peta[i];
    e_Peta[3 * i] = c_Peta[0] * tau->contents[i];
    b_LtauItau0[3 * i] = LtauItau0[0] * LtauItau0[i];
    g_tau[3 * i] = tau->contents[0] * Peta[i];
    f_Peta[3 * i] = Peta[0] * tau->contents[i];
    c_LtauItau0[3 * i] = LtauItau0[0] * LGIeta[i];
    b_LGIeta[3 * i] = LGIeta[0] * LtauItau0[i];
    i1 = 3 * i + 1;
    d = I3->contents[i1];
    Q[i1] = c_a * d + tau->contents[1] * LKIeta[i];
    c_LKIeta[i1] = LKIeta[1] * tau->contents[i];
    e_a[i1] = d_a * d + tau->contents[1] * b_LKIeta[i];
    d_LKIeta[i1] = b_LKIeta[1] * tau->contents[i];
    b_LKItau0[i1] = LKItau0[1] * LtauIeta[i];
    b_LtauIeta[i1] = LtauIeta[1] * LKItau0[i];
    e_tau[i1] = tau->contents[1] * b_Peta[i];
    d_Peta[i1] = b_Peta[1] * tau->contents[i];
    f_tau[i1] = tau->contents[1] * c_Peta[i];
    e_Peta[i1] = c_Peta[1] * tau->contents[i];
    b_LtauItau0[i1] = LtauItau0[1] * LtauItau0[i];
    g_tau[i1] = tau->contents[1] * Peta[i];
    f_Peta[i1] = Peta[1] * tau->contents[i];
    c_LtauItau0[i1] = LtauItau0[1] * LGIeta[i];
    b_LGIeta[i1] = LGIeta[1] * LtauItau0[i];
    i1 = 3 * i + 2;
    d = I3->contents[i1];
    Q[i1] = c_a * d + tau->contents[2] * LKIeta[i];
    c_LKIeta[i1] = LKIeta[2] * tau->contents[i];
    e_a[i1] = d_a * d + tau->contents[2] * b_LKIeta[i];
    d_LKIeta[i1] = b_LKIeta[2] * tau->contents[i];
    b_LKItau0[i1] = LKItau0[2] * LtauIeta[i];
    b_LtauIeta[i1] = LtauIeta[2] * LKItau0[i];
    e_tau[i1] = tau->contents[2] * b_Peta[i];
    d_Peta[i1] = b_Peta[2] * tau->contents[i];
    f_tau[i1] = tau->contents[2] * c_Peta[i];
    e_Peta[i1] = c_Peta[2] * tau->contents[i];
    b_LtauItau0[i1] = LtauItau0[2] * LtauItau0[i];
    g_tau[i1] = tau->contents[2] * Peta[i];
    f_Peta[i1] = Peta[2] * tau->contents[i];
    c_LtauItau0[i1] = LtauItau0[2] * LGIeta[i];
    b_LGIeta[i1] = LGIeta[2] * LtauItau0[i];
  }
  for (i = 0; i < 9; i++) {
    d = PP->contents[i];
    Q[i] = (((-(Q[i] + c_LKIeta[i]) * b_tmp -
              (((d_eta * (-(e_a[i] + d_LKIeta[i]) * b_tmp) +
                 (b_LKItau0[i] + b_LtauIeta[i])) +
                c_tau0 * (-((e_tau[i] + d_Peta[i]) + c_tau * d) / c_tmp)) -
               e_eta * (-((f_tau[i] + e_Peta[i]) + d_tau * d) / c_tmp)) /
                  (tau0_tmp + 1.0)) -
             a * b_LtauItau0[i]) +
            b_a * (-((g_tau[i] + f_Peta[i]) + b_tau * d) / c_tmp)) +
           (c_LtauItau0[i] + b_LGIeta[i]) / a_tmp;
  }
}

static void Q_rr_III(const b_captured_var *nxip, const c_captured_var *L_KII,
                     const captured_var *tau, const captured_var *tau0,
                     const c_captured_var *L_tauI, const c_captured_var *L_GII,
                     const double eta[3], double Q[9])
{
  double c_L_KII[9];
  double dv[9];
  double dv1[9];
  double e_L_KII[9];
  double b_L_KII[3];
  double b_L_tauI[3];
  double d_L_KII[3];
  double a_tmp;
  double b_a_tmp;
  double b_eta;
  double b_tau0;
  double c;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  int i;
  b_tau0 = 0.0;
  a_tmp = -(2.0 / nxip->contents);
  b_a_tmp = nxip->contents;
  b_a_tmp *= b_a_tmp;
  b_eta = 0.0;
  for (i = 0; i < 3; i++) {
    d = tau0->contents[i];
    d1 = tau->contents[i];
    b_tau0 += d * d1;
    b_eta += eta[i] * (2.0 * d + d1);
    b_L_KII[i] =
        (L_KII->contents[i] * eta[0] + L_KII->contents[i + 3] * eta[1]) +
        L_KII->contents[i + 6] * eta[2];
  }
  c = (b_tau0 + 1.0) * (b_tau0 + 1.0);
  dv[0] = 0.0 / b_a_tmp;
  dv[3] = -eta[2] / b_a_tmp;
  dv[6] = eta[1] / b_a_tmp;
  dv[1] = eta[2] / b_a_tmp;
  dv[4] = 0.0 / b_a_tmp;
  dv[7] = -eta[0] / b_a_tmp;
  dv[2] = -eta[1] / b_a_tmp;
  dv[5] = eta[0] / b_a_tmp;
  dv[8] = 0.0 / b_a_tmp;
  d = b_L_KII[0];
  d1 = b_L_KII[1];
  d2 = b_L_KII[2];
  d3 = eta[0];
  d4 = eta[1];
  d5 = eta[2];
  for (i = 0; i < 3; i++) {
    d6 = tau->contents[i];
    c_L_KII[3 * i] = d * d6;
    c_L_KII[3 * i + 1] = d1 * d6;
    c_L_KII[3 * i + 2] = d2 * d6;
    b_L_tauI[i] = (L_tauI->contents[i] * d3 + L_tauI->contents[i + 3] * d4) +
                  L_tauI->contents[i + 6] * d5;
  }
  d = tau0->contents[0];
  d1 = tau0->contents[1];
  d2 = tau0->contents[2];
  for (i = 0; i < 3; i++) {
    d3 = (L_KII->contents[i] * d + L_KII->contents[i + 3] * d1) +
         L_KII->contents[i + 6] * d2;
    d_L_KII[i] = d3;
    b_L_KII[i] = d3;
  }
  d = d_L_KII[0];
  d1 = d_L_KII[1];
  d2 = d_L_KII[2];
  for (i = 0; i < 3; i++) {
    d3 = tau->contents[i];
    e_L_KII[3 * i] = d * d3;
    e_L_KII[3 * i + 1] = d1 * d3;
    e_L_KII[3 * i + 2] = d2 * d3;
  }
  dv1[0] = 0.0 / b_a_tmp;
  dv1[3] = -tau0->contents[2] / b_a_tmp;
  dv1[6] = tau0->contents[1] / b_a_tmp;
  dv1[1] = tau0->contents[2] / b_a_tmp;
  dv1[4] = 0.0 / b_a_tmp;
  dv1[7] = -tau0->contents[0] / b_a_tmp;
  dv1[2] = -tau0->contents[1] / b_a_tmp;
  dv1[5] = tau0->contents[0] / b_a_tmp;
  dv1[8] = 0.0 / b_a_tmp;
  d = b_L_KII[0];
  d1 = b_L_KII[1];
  d2 = b_L_KII[2];
  d3 = eta[0];
  d4 = eta[1];
  d5 = eta[2];
  d6 = tau0->contents[0];
  b_a_tmp = tau0->contents[1];
  d7 = tau0->contents[2];
  for (i = 0; i < 3; i++) {
    double d8;
    int Q_tmp;
    d8 = b_L_tauI[i];
    Q[3 * i] = (a_tmp * c_L_KII[3 * i] + dv[3 * i]) -
               (d * d8 + b_eta * (a_tmp * e_L_KII[3 * i] + dv1[3 * i])) /
                   (b_tau0 + 1.0);
    Q_tmp = 3 * i + 1;
    Q[Q_tmp] = (a_tmp * c_L_KII[Q_tmp] + dv[Q_tmp]) -
               (d1 * d8 + b_eta * (a_tmp * e_L_KII[Q_tmp] + dv1[Q_tmp])) /
                   (b_tau0 + 1.0);
    Q_tmp = 3 * i + 2;
    Q[Q_tmp] = (a_tmp * c_L_KII[Q_tmp] + dv[Q_tmp]) -
               (d2 * d8 + b_eta * (a_tmp * e_L_KII[Q_tmp] + dv1[Q_tmp])) /
                   (b_tau0 + 1.0);
    b_L_KII[i] = (L_GII->contents[i] * d3 + L_GII->contents[i + 3] * d4) +
                 L_GII->contents[i + 6] * d5;
    d_L_KII[i] =
        (L_tauI->contents[i] * d6 + L_tauI->contents[i + 3] * b_a_tmp) +
        L_tauI->contents[i + 6] * d7;
  }
  d = b_L_KII[0];
  d1 = b_L_KII[1];
  d2 = b_L_KII[2];
  for (i = 0; i < 3; i++) {
    d3 = d_L_KII[i];
    c_L_KII[3 * i] = d * d3 / c;
    c_L_KII[3 * i + 1] = d1 * d3 / c;
    c_L_KII[3 * i + 2] = d2 * d3 / c;
  }
  for (i = 0; i < 9; i++) {
    Q[i] += c_L_KII[i];
  }
}

static void Qtil_rr_I(const b_captured_var *nxip, const captured_var *K,
                      const c_captured_var *L_KI, const captured_var *tau,
                      const c_captured_var *I3, const c_captured_var *PP,
                      const captured_var *tau0, const c_captured_var *L_tauI,
                      const captured_var *K0, const captured_var *G,
                      const c_captured_var *L_GI, const double nu[3],
                      double Q[9])
{
  double b_G[9];
  double b_K[9];
  double b_L_KI[9];
  double b_nu[9];
  double b_y_tmp[9];
  double c_a[9];
  double d_Peta[9];
  double d_y_tmp[9];
  double dv[9];
  double g_tau[9];
  double y_tmp[9];
  double LKIeta[3];
  double LKItau0[3];
  double LtauItau0[3];
  double Peta[3];
  double b_Peta[3];
  double c_Peta[3];
  double c_y_tmp[3];
  double f_tau[3];
  double a;
  double a_tmp;
  double b_LKItau0;
  double b_LtauItau0;
  double b_a;
  double b_tau;
  double b_tau0;
  double c;
  double c_tau;
  double c_tau0;
  double c_tmp;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double d_tau;
  double e_tau;
  int i;
  int i1;
  int y_tmp_tmp;
  b_tau0 = (tau0->contents[0] * tau->contents[0] +
            tau0->contents[1] * tau->contents[1]) +
           tau0->contents[2] * tau->contents[2];
  b_LtauItau0 = 0.0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  UTILITY FUNCTIONS FOR SECOND DERIVATIVE COMPUTATIONS */
  /*  It is assumed in these functions that quantities xip, xipp,  */
  /*  their norms (nxip, nxipp), tau0, tau, PP, L_tauI, */
  /*  K0, K, L_KI, L_KII, G, L_GI, L_GII have been computed and are available */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  function related to tau */
  b_tau = 0.0;
  a_tmp = nxip->contents;
  c_tmp = a_tmp * a_tmp;
  c = (b_tau0 + 1.0) * (b_tau0 + 1.0);
  c_tau = 0.0;
  a_tmp = 2.0 / nxip->contents;
  /*  functions related to K */
  c_tau0 = 0.0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  UTILITY FUNCTIONS FOR SECOND DERIVATIVE COMPUTATIONS */
  /*  It is assumed in these functions that quantities xip, xipp,  */
  /*  their norms (nxip, nxipp), tau0, tau, PP, L_tauI, */
  /*  K0, K, L_KI, L_KII, G, L_GI, L_GII have been computed and are available */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  function related to tau */
  d_tau = 0.0;
  for (i = 0; i < 3; i++) {
    d = tau0->contents[0];
    d1 = L_tauI->contents[i] * d;
    d2 = PP->contents[i];
    d3 = d2 * d;
    d = tau0->contents[1];
    d1 += L_tauI->contents[i + 3] * d;
    d4 = PP->contents[i + 3];
    d3 += d4 * d;
    d = tau0->contents[2];
    d1 += L_tauI->contents[i + 6] * d;
    d5 = PP->contents[i + 6];
    d3 += d5 * d;
    Peta[i] = d3;
    LtauItau0[i] = d1;
    d = nu[i];
    b_LtauItau0 += d1 * d;
    d1 = tau->contents[i];
    b_tau += d1 * tau0->contents[i];
    y_tmp[3 * i] = L_GI->contents[i];
    y_tmp[3 * i + 1] = L_GI->contents[i + 3];
    y_tmp[3 * i + 2] = L_GI->contents[i + 6];
    c_tau += d1 * d;
    d = (L_KI->contents[i] * tau0->contents[0] +
         L_KI->contents[i + 3] * tau0->contents[1]) +
        L_KI->contents[i + 6] * tau0->contents[2];
    LKIeta[i] = d;
    LKItau0[i] = d;
    c_tau0 += tau0->contents[i] * K->contents[i];
    b_Peta[i] =
        (d2 * K0->contents[0] + d4 * K0->contents[1]) + d5 * K0->contents[2];
    d_tau += d1 * K0->contents[i];
  }
  a = -(2.0 / (b_tau0 + 1.0)) * b_LtauItau0;
  b_a = c_tau0 / nxip->contents;
  b_LKItau0 = 0.0;
  c_tau0 = 0.0;
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  UTILITY FUNCTIONS FOR SECOND DERIVATIVE COMPUTATIONS */
  /*  It is assumed in these functions that quantities xip, xipp,  */
  /*  their norms (nxip, nxipp), tau0, tau, PP, L_tauI, */
  /*  K0, K, L_KI, L_KII, G, L_GI, L_GII have been computed and are available */
  /* %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
   */
  /*  function related to tau */
  e_tau = 0.0;
  for (i = 0; i < 3; i++) {
    d = nu[i];
    b_LKItau0 += LKItau0[i] * d;
    c_tau0 += tau0->contents[i] * K->contents[i];
    d1 = tau->contents[i];
    e_tau += d1 * d;
    b_y_tmp[3 * i] = L_tauI->contents[i];
    b_nu[3 * i] = nu[0] * d1 + c_tau * I3->contents[3 * i];
    b_K[3 * i] = K->contents[0] * nu[i] / nxip->contents;
    y_tmp_tmp = 3 * i + 1;
    b_y_tmp[y_tmp_tmp] = L_tauI->contents[i + 3];
    b_nu[y_tmp_tmp] = nu[1] * d1 + c_tau * I3->contents[y_tmp_tmp];
    b_K[y_tmp_tmp] = K->contents[1] * nu[i] / nxip->contents;
    y_tmp_tmp = 3 * i + 2;
    b_y_tmp[y_tmp_tmp] = L_tauI->contents[i + 6];
    b_nu[y_tmp_tmp] = nu[2] * d1 + c_tau * I3->contents[y_tmp_tmp];
    b_K[y_tmp_tmp] = K->contents[2] * nu[i] / nxip->contents;
    c_Peta[i] = (PP->contents[i] * nu[0] + PP->contents[i + 3] * nu[1]) +
                PP->contents[i + 6] * nu[2];
  }
  for (i = 0; i < 3; i++) {
    d = L_KI->contents[3 * i];
    d1 = L_KI->contents[3 * i + 1];
    d2 = L_KI->contents[3 * i + 2];
    for (i1 = 0; i1 < 3; i1++) {
      b_L_KI[i + 3 * i1] =
          (d * b_nu[3 * i1] + d1 * b_nu[3 * i1 + 1]) + d2 * b_nu[3 * i1 + 2];
      y_tmp_tmp = i1 + 3 * i;
      c_a[y_tmp_tmp] =
          b_a * I3->contents[y_tmp_tmp] + tau->contents[i1] * LKIeta[i];
    }
  }
  d = LKIeta[0];
  d1 = LKIeta[1];
  d2 = LKIeta[2];
  for (i = 0; i < 3; i++) {
    d3 = tau->contents[i];
    b_nu[3 * i] = d * d3;
    b_nu[3 * i + 1] = d1 * d3;
    b_nu[3 * i + 2] = d2 * d3;
  }
  for (i = 0; i < 9; i++) {
    c_a[i] = -(c_a[i] + b_nu[i]) * a_tmp;
  }
  d = nu[0];
  d1 = nu[1];
  d2 = nu[2];
  for (i = 0; i < 3; i++) {
    f_tau[i] = 2.0 * tau0->contents[i] + tau->contents[i];
    c_y_tmp[i] = (b_y_tmp[i] * d + b_y_tmp[i + 3] * d1) + b_y_tmp[i + 6] * d2;
    LKIeta[i] = (c_a[i] * d + c_a[i + 3] * d1) + c_a[i + 6] * d2;
  }
  d = f_tau[0];
  d1 = f_tau[1];
  d2 = f_tau[2];
  d3 = c_y_tmp[0];
  d4 = c_y_tmp[1];
  d5 = c_y_tmp[2];
  for (i = 0; i < 3; i++) {
    b_a = LKIeta[i];
    dv[3 * i] = d * b_a;
    c_tau = LKItau0[i];
    d_y_tmp[3 * i] = d3 * c_tau;
    g_tau[3 * i] = tau->contents[0] * b_Peta[i];
    d_Peta[3 * i] = b_Peta[0] * tau->contents[i];
    i1 = 3 * i + 1;
    dv[i1] = d1 * b_a;
    d_y_tmp[i1] = d4 * c_tau;
    g_tau[i1] = tau->contents[1] * b_Peta[i];
    d_Peta[i1] = b_Peta[1] * tau->contents[i];
    i1 = 3 * i + 2;
    dv[i1] = d2 * b_a;
    d_y_tmp[i1] = d5 * c_tau;
    g_tau[i1] = tau->contents[2] * b_Peta[i];
    d_Peta[i1] = b_Peta[2] * tau->contents[i];
  }
  for (i = 0; i < 9; i++) {
    g_tau[i] = -((g_tau[i] + d_Peta[i]) + d_tau * PP->contents[i]) / c_tmp;
  }
  d = nu[0];
  d1 = nu[1];
  d2 = nu[2];
  for (i = 0; i < 3; i++) {
    f_tau[i] = (g_tau[i] * d + g_tau[i + 3] * d1) + g_tau[i + 6] * d2;
  }
  d = G->contents[0];
  d1 = G->contents[1];
  d2 = G->contents[2];
  for (i = 0; i < 3; i++) {
    g_tau[3 * i] = tau->contents[0] * c_Peta[i];
    d_Peta[3 * i] = c_Peta[0] * tau->contents[i];
    d3 = LtauItau0[i];
    b_G[3 * i] = d * d3;
    b_nu[3 * i] = tau->contents[0] * Peta[i];
    c_a[3 * i] = Peta[0] * tau->contents[i];
    y_tmp_tmp = 3 * i + 1;
    g_tau[y_tmp_tmp] = tau->contents[1] * c_Peta[i];
    d_Peta[y_tmp_tmp] = c_Peta[1] * tau->contents[i];
    b_G[y_tmp_tmp] = d1 * d3;
    b_nu[y_tmp_tmp] = tau->contents[1] * Peta[i];
    c_a[y_tmp_tmp] = Peta[1] * tau->contents[i];
    y_tmp_tmp = 3 * i + 2;
    g_tau[y_tmp_tmp] = tau->contents[2] * c_Peta[i];
    d_Peta[y_tmp_tmp] = c_Peta[2] * tau->contents[i];
    b_G[y_tmp_tmp] = d2 * d3;
    b_nu[y_tmp_tmp] = tau->contents[2] * Peta[i];
    c_a[y_tmp_tmp] = Peta[2] * tau->contents[i];
  }
  for (i = 0; i < 9; i++) {
    b_nu[i] = -((b_nu[i] + c_a[i]) + b_tau * PP->contents[i]) / c_tmp;
  }
  d = nu[0];
  d1 = nu[1];
  d2 = nu[2];
  d3 = G->contents[0];
  d4 = G->contents[1];
  d5 = G->contents[2];
  for (i = 0; i < 3; i++) {
    c_y_tmp[i] = (y_tmp[i] * d + y_tmp[i + 3] * d1) + y_tmp[i + 6] * d2;
    b_a = (b_nu[i] * d + b_nu[i + 3] * d1) + b_nu[i + 6] * d2;
    c_a[3 * i] = a * b_G[3 * i] + d3 * b_a;
    y_tmp_tmp = 3 * i + 1;
    c_a[y_tmp_tmp] = a * b_G[y_tmp_tmp] + d4 * b_a;
    y_tmp_tmp = 3 * i + 2;
    c_a[y_tmp_tmp] = a * b_G[y_tmp_tmp] + d5 * b_a;
  }
  d = tau0->contents[0];
  d1 = tau0->contents[1];
  d2 = tau0->contents[2];
  d3 = c_y_tmp[0];
  d4 = c_y_tmp[1];
  d5 = c_y_tmp[2];
  for (i = 0; i < 3; i++) {
    b_a = f_tau[i];
    c_tau = LtauItau0[i];
    Q[3 * i] =
        (-a_tmp * (b_K[3 * i] + b_L_KI[3 * i]) -
         ((((dv[3 * i] + d_y_tmp[3 * i]) + b_LKItau0 * b_y_tmp[3 * i]) -
           d * b_a) +
          c_tau0 *
              (-((g_tau[3 * i] + d_Peta[3 * i]) + e_tau * PP->contents[3 * i]) /
               c_tmp)) /
             (b_tau0 + 1.0)) +
        ((c_a[3 * i] + d3 * c_tau) + b_LtauItau0 * y_tmp[3 * i]) / c;
    y_tmp_tmp = 3 * i + 1;
    Q[y_tmp_tmp] =
        (-a_tmp * (b_K[y_tmp_tmp] + b_L_KI[y_tmp_tmp]) -
         ((((dv[y_tmp_tmp] + d_y_tmp[y_tmp_tmp]) +
            b_LKItau0 * b_y_tmp[y_tmp_tmp]) -
           d1 * b_a) +
          c_tau0 * (-((g_tau[y_tmp_tmp] + d_Peta[y_tmp_tmp]) +
                      e_tau * PP->contents[y_tmp_tmp]) /
                    c_tmp)) /
             (b_tau0 + 1.0)) +
        ((c_a[y_tmp_tmp] + d4 * c_tau) + b_LtauItau0 * y_tmp[y_tmp_tmp]) / c;
    y_tmp_tmp = 3 * i + 2;
    Q[y_tmp_tmp] =
        (-a_tmp * (b_K[y_tmp_tmp] + b_L_KI[y_tmp_tmp]) -
         ((((dv[y_tmp_tmp] + d_y_tmp[y_tmp_tmp]) +
            b_LKItau0 * b_y_tmp[y_tmp_tmp]) -
           d2 * b_a) +
          c_tau0 * (-((g_tau[y_tmp_tmp] + d_Peta[y_tmp_tmp]) +
                      e_tau * PP->contents[y_tmp_tmp]) /
                    c_tmp)) /
             (b_tau0 + 1.0)) +
        ((c_a[y_tmp_tmp] + d5 * c_tau) + b_LtauItau0 * y_tmp[y_tmp_tmp]) / c;
  }
}

static void Qtil_rr_III(const b_captured_var *nxip, const captured_var *tau,
                        const c_captured_var *L_KII,
                        const c_captured_var *L_tauI, const captured_var *tau0,
                        const c_captured_var *L_GII, const double nu[3],
                        double Q[9])
{
  double d_L_KII[9];
  double dv[9];
  double dv1[9];
  double b_L_KII[3];
  double b_L_tauI[3];
  double c_L_KII[3];
  double dv2[3];
  double a_tmp;
  double b_tau;
  double b_y;
  double c_L_tauI;
  double d;
  double d1;
  double d2;
  double d3;
  double tau0_tmp;
  double y;
  double y_tmp;
  int i;
  tau0_tmp = (tau0->contents[0] * tau->contents[0] +
              tau0->contents[1] * tau->contents[1]) +
             tau0->contents[2] * tau->contents[2];
  c_L_tauI = 0.0;
  b_tau = 0.0;
  y_tmp = -(2.0 / nxip->contents);
  a_tmp = nxip->contents;
  a_tmp *= a_tmp;
  dv[0] = 0.0 / a_tmp;
  dv[3] = -nu[2] / a_tmp;
  dv[6] = nu[1] / a_tmp;
  dv[1] = nu[2] / a_tmp;
  dv[4] = 0.0 / a_tmp;
  dv[7] = -nu[0] / a_tmp;
  dv[2] = -nu[1] / a_tmp;
  dv[5] = nu[0] / a_tmp;
  dv[8] = 0.0 / a_tmp;
  d = tau0->contents[0];
  d1 = tau0->contents[1];
  d2 = tau0->contents[2];
  for (i = 0; i < 3; i++) {
    d3 = nu[i];
    c_L_tauI += ((L_tauI->contents[i] * d + L_tauI->contents[i + 3] * d1) +
                 L_tauI->contents[i + 6] * d2) *
                d3;
    b_tau += tau->contents[i] * d3;
    d3 = (L_KII->contents[i] * d + L_KII->contents[i + 3] * d1) +
         L_KII->contents[i + 6] * d2;
    c_L_KII[i] = d3;
    b_L_tauI[i] = (L_tauI->contents[3 * i] * nu[0] +
                   L_tauI->contents[3 * i + 1] * nu[1]) +
                  L_tauI->contents[3 * i + 2] * nu[2];
    b_L_KII[i] = d3;
  }
  y = c_L_tauI / ((tau0_tmp + 1.0) * (tau0_tmp + 1.0));
  b_y = y_tmp * b_tau;
  d = c_L_KII[0];
  d1 = c_L_KII[1];
  d2 = c_L_KII[2];
  for (i = 0; i < 3; i++) {
    d3 = tau->contents[i];
    d_L_KII[3 * i] = d * d3;
    d_L_KII[3 * i + 1] = d1 * d3;
    d_L_KII[3 * i + 2] = d2 * d3;
  }
  dv1[0] = 0.0 / a_tmp;
  dv1[3] = -tau0->contents[2] / a_tmp;
  dv1[6] = tau0->contents[1] / a_tmp;
  dv1[1] = tau0->contents[2] / a_tmp;
  dv1[4] = 0.0 / a_tmp;
  dv1[7] = -tau0->contents[0] / a_tmp;
  dv1[2] = -tau0->contents[1] / a_tmp;
  dv1[5] = tau0->contents[0] / a_tmp;
  dv1[8] = 0.0 / a_tmp;
  for (i = 0; i < 9; i++) {
    d_L_KII[i] = y_tmp * d_L_KII[i] + dv1[i];
  }
  d = nu[0];
  d1 = nu[1];
  d2 = nu[2];
  for (i = 0; i < 3; i++) {
    c_L_KII[i] = (d_L_KII[i] * d + d_L_KII[i + 3] * d1) + d_L_KII[i + 6] * d2;
    dv2[i] = 2.0 * tau0->contents[i] + tau->contents[i];
  }
  d = b_L_tauI[0];
  d1 = b_L_tauI[1];
  d2 = b_L_tauI[2];
  d3 = dv2[0];
  a_tmp = dv2[1];
  c_L_tauI = dv2[2];
  for (i = 0; i < 3; i++) {
    int Q_tmp;
    b_tau = b_L_KII[i];
    y_tmp = c_L_KII[i];
    Q[3 * i] = ((b_y * L_KII->contents[i] + dv[3 * i]) -
                (d * b_tau + d3 * y_tmp) / (tau0_tmp + 1.0)) +
               y * L_GII->contents[i];
    Q_tmp = 3 * i + 1;
    Q[Q_tmp] = ((b_y * L_KII->contents[i + 3] + dv[Q_tmp]) -
                (d1 * b_tau + a_tmp * y_tmp) / (tau0_tmp + 1.0)) +
               y * L_GII->contents[i + 3];
    Q_tmp = 3 * i + 2;
    Q[Q_tmp] = ((b_y * L_KII->contents[i + 6] + dv[Q_tmp]) -
                (d2 * b_tau + c_L_tauI * y_tmp) / (tau0_tmp + 1.0)) +
               y * L_GII->contents[i + 6];
  }
}

static void Qtil_rr_IV(const b_captured_var *nxip, const c_captured_var *L_KII,
                       const captured_var *tau, const captured_var *tau0,
                       const c_captured_var *L_tauI,
                       const c_captured_var *L_GII, const double nu[3],
                       double Q[9])
{
  double a[9];
  double d_L_KII[9];
  double dv[9];
  double dv1[9];
  double b_L_GII[3];
  double b_L_KII[3];
  double b_L_tauI[3];
  double dv2[3];
  double a_tmp;
  double b_a_tmp;
  double b_tau0;
  double c;
  double c_L_KII;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  int i;
  b_tau0 = (tau0->contents[0] * tau->contents[0] +
            tau0->contents[1] * tau->contents[1]) +
           tau0->contents[2] * tau->contents[2];
  c = (b_tau0 + 1.0) * (b_tau0 + 1.0);
  a_tmp = -(2.0 / nxip->contents);
  b_a_tmp = nxip->contents;
  b_a_tmp *= b_a_tmp;
  c_L_KII = 0.0;
  d = tau0->contents[0];
  d1 = tau0->contents[1];
  d2 = tau0->contents[2];
  for (i = 0; i < 3; i++) {
    c_L_KII += ((L_KII->contents[i] * d + L_KII->contents[i + 3] * d1) +
                L_KII->contents[i + 6] * d2) *
               nu[i];
    b_L_KII[i] =
        a_tmp *
        ((L_KII->contents[3 * i] * nu[0] + L_KII->contents[3 * i + 1] * nu[1]) +
         L_KII->contents[3 * i + 2] * nu[2]);
  }
  d = b_L_KII[0];
  d1 = b_L_KII[1];
  d2 = b_L_KII[2];
  for (i = 0; i < 3; i++) {
    d3 = tau->contents[i];
    a[3 * i] = d * d3;
    a[3 * i + 1] = d1 * d3;
    a[3 * i + 2] = d2 * d3;
  }
  dv[0] = 0.0 / b_a_tmp;
  dv[3] = -nu[2] / b_a_tmp;
  dv[6] = nu[1] / b_a_tmp;
  dv[1] = nu[2] / b_a_tmp;
  dv[4] = 0.0 / b_a_tmp;
  dv[7] = -nu[0] / b_a_tmp;
  dv[2] = -nu[1] / b_a_tmp;
  dv[5] = nu[0] / b_a_tmp;
  dv[8] = 0.0 / b_a_tmp;
  d = tau0->contents[0];
  d1 = tau0->contents[1];
  d2 = tau0->contents[2];
  d3 = tau->contents[0];
  d4 = tau->contents[1];
  d5 = tau->contents[2];
  for (i = 0; i < 3; i++) {
    d6 = (L_KII->contents[i] * d + L_KII->contents[i + 3] * d1) +
         L_KII->contents[i + 6] * d2;
    d_L_KII[3 * i] = d6 * d3;
    d_L_KII[3 * i + 1] = d6 * d4;
    d_L_KII[3 * i + 2] = d6 * d5;
  }
  dv1[0] = 0.0 / b_a_tmp;
  dv1[1] = -tau0->contents[2] / b_a_tmp;
  dv1[2] = tau0->contents[1] / b_a_tmp;
  dv1[3] = tau0->contents[2] / b_a_tmp;
  dv1[4] = 0.0 / b_a_tmp;
  dv1[5] = -tau0->contents[0] / b_a_tmp;
  dv1[6] = -tau0->contents[1] / b_a_tmp;
  dv1[7] = tau0->contents[0] / b_a_tmp;
  dv1[8] = 0.0 / b_a_tmp;
  for (i = 0; i < 9; i++) {
    d_L_KII[i] = a_tmp * d_L_KII[i] + dv1[i];
  }
  d = nu[0];
  d1 = nu[1];
  d2 = nu[2];
  for (i = 0; i < 3; i++) {
    dv2[i] = 2.0 * tau0->contents[i] + tau->contents[i];
    b_L_GII[i] =
        (L_GII->contents[3 * i] * d + L_GII->contents[3 * i + 1] * d1) +
        L_GII->contents[3 * i + 2] * d2;
    b_L_tauI[i] = (L_tauI->contents[i] * tau0->contents[0] +
                   L_tauI->contents[i + 3] * tau0->contents[1]) +
                  L_tauI->contents[i + 6] * tau0->contents[2];
    b_L_KII[i] = (d_L_KII[i] * d + d_L_KII[i + 3] * d1) + d_L_KII[i + 6] * d2;
  }
  d = dv2[0];
  d1 = dv2[1];
  d2 = dv2[2];
  d3 = b_L_GII[0];
  d4 = b_L_GII[1];
  d5 = b_L_GII[2];
  for (i = 0; i < 3; i++) {
    int Q_tmp;
    d6 = b_L_KII[i];
    Q[3 * i] = (a[3 * i] - dv[3 * i]) -
               (c_L_KII * L_tauI->contents[i] + d * d6) / (b_tau0 + 1.0);
    b_a_tmp = b_L_tauI[i];
    d_L_KII[3 * i] = d3 * b_a_tmp / c;
    Q_tmp = 3 * i + 1;
    Q[Q_tmp] = (a[Q_tmp] - dv[Q_tmp]) -
               (c_L_KII * L_tauI->contents[i + 3] + d1 * d6) / (b_tau0 + 1.0);
    d_L_KII[Q_tmp] = d4 * b_a_tmp / c;
    Q_tmp = 3 * i + 2;
    Q[Q_tmp] = (a[Q_tmp] - dv[Q_tmp]) -
               (c_L_KII * L_tauI->contents[i + 6] + d2 * d6) / (b_tau0 + 1.0);
    d_L_KII[Q_tmp] = d5 * b_a_tmp / c;
  }
  for (i = 0; i < 9; i++) {
    Q[i] += d_L_KII[i];
  }
}

static void binary_expand_op_15(emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3)
{
  emxArray_real_T *b_in1;
  const double *in2_data;
  const double *in3_data;
  double *b_in1_data;
  double *in1_data;
  int aux_0_1;
  int aux_1_1;
  int aux_2_1;
  int b_loop_ub;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  int stride_2_0;
  int stride_2_1;
  in3_data = in3->data;
  in2_data = in2->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  if (in3->size[0] == 1) {
    if (in2->size[0] == 1) {
      loop_ub = in1->size[0];
    } else {
      loop_ub = in2->size[0];
    }
  } else {
    loop_ub = in3->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  if (in3->size[1] == 1) {
    if (in2->size[1] == 1) {
      b_loop_ub = in1->size[1];
    } else {
      b_loop_ub = in2->size[1];
    }
  } else {
    b_loop_ub = in3->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_0_1 = (in1->size[1] != 1);
  stride_1_0 = (in2->size[0] != 1);
  stride_1_1 = (in2->size[1] != 1);
  stride_2_0 = (in3->size[0] != 1);
  stride_2_1 = (in3->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  aux_2_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          (in1_data[i1 * stride_0_0 + in1->size[0] * aux_0_1] +
           in2_data[i1 * stride_1_0 + in2->size[0] * aux_1_1]) +
          in3_data[i1 * stride_2_0 + in3->size[0] * aux_2_1];
    }
    aux_2_1 += stride_2_1;
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = b_in1->size[0];
  in1->size[1] = b_in1->size[1];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  loop_ub = b_in1->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = b_in1->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real_T(&b_in1);
}

static void binary_expand_op_4(emxArray_real_T *in1,
                               const emxArray_int32_T *in2,
                               const emxArray_real_T *in3,
                               const emxArray_real_T *in4, double in5,
                               double in6, const emxArray_real_T *in7, int in8)
{
  const double *in3_data;
  const double *in4_data;
  const double *in7_data;
  double b_in7;
  double *in1_data;
  const int *in2_data;
  int i;
  int loop_ub;
  int stride_0_1;
  int stride_1_1;
  in7_data = in7->data;
  in4_data = in4->data;
  in3_data = in3->data;
  in2_data = in2->data;
  in1_data = in1->data;
  b_in7 = in7_data[in8];
  stride_0_1 = (in3->size[1] != 1);
  stride_1_1 = (in4->size[1] != 1);
  loop_ub = in2->size[1];
  for (i = 0; i < loop_ub; i++) {
    in1_data[in2_data[i] - 1] =
        in3_data[i * stride_0_1] + in4_data[i * stride_1_1] * in5 * in6 * b_in7;
  }
}

static void binary_expand_op_5(emxArray_real_T *in1,
                               const emxArray_int32_T *in2,
                               const emxArray_int32_T *in3,
                               const emxArray_real_T *in4,
                               const emxArray_real_T *in5)
{
  emxArray_real_T *b_in1;
  const double *in4_data;
  const double *in5_data;
  double *b_in1_data;
  double *in1_data;
  const int *in2_data;
  const int *in3_data;
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

static void binary_expand_op_7(emxArray_real_T *in1,
                               const emxArray_int32_T *in2, int in3, int in4,
                               const emxArray_real_T *in5, int in6, int in7,
                               const emxArray_real_T *in8)
{
  emxArray_real_T *b_in1;
  const double *in5_data;
  const double *in8_data;
  double *b_in1_data;
  double *in1_data;
  const int *in2_data;
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

static void binary_expand_op_8(emxArray_real_T *in1, int in2, int in3,
                               const emxArray_int32_T *in4, int in5, int in6,
                               const emxArray_real_T *in7,
                               const emxArray_real_T *in8)
{
  emxArray_real_T *b_in1;
  const double *in7_data;
  const double *in8_data;
  double *b_in1_data;
  double *in1_data;
  const int *in4_data;
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

static void plus(emxArray_real_T *in1, const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const double *in2_data;
  double *b_in1_data;
  double *in1_data;
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
  in2_data = in2->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  if (in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = in2->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  if (in2->size[1] == 1) {
    b_loop_ub = in1->size[1];
  } else {
    b_loop_ub = in2->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_0_1 = (in1->size[1] != 1);
  stride_1_0 = (in2->size[0] != 1);
  stride_1_1 = (in2->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * aux_0_1] +
          in2_data[i1 * stride_1_0 + in2->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = b_in1->size[0];
  in1->size[1] = b_in1->size[1];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  loop_ub = b_in1->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = b_in1->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real_T(&b_in1);
}

void CableForce(const emxArray_real_T *P, const emxArray_real_T *P0,
                const emxArray_real_T *Pdot, const emxArray_real_T *varTheta,
                const emxArray_real_T *varThetadot, const emxArray_real_T *R0,
                double rho, double EA, double EI, double GJ, double betAX,
                double betBEND, double betTOR, const emxArray_real_T *wg,
                const emxArray_real_T *nel, const emxArray_real_T *colmat,
                const emxArray_real_T *colmat_brev,
                const emxArray_real_T *colmat_bar, double d, double dbrev,
                double dbar, const emxArray_real_T *Mbar, const double u[3],
                const emxArray_real_T *Kbar11, const emxArray_real_T *Dbar11,
                emxArray_real_T *F, emxArray_real_T *mu, emxArray_real_T *F_ij,
                emxArray_real_T *F_ib, emxArray_real_T *F_ijd,
                emxArray_real_T *F_ibd, emxArray_real_T *mu_aj,
                emxArray_real_T *mu_ab, emxArray_real_T *mu_ajd,
                emxArray_real_T *mu_abd)
{
  b_captured_var nxip;
  c_captured_var I3;
  c_captured_var L_GI;
  c_captured_var L_GII;
  c_captured_var L_KI;
  c_captured_var L_KII;
  c_captured_var L_tauI;
  c_captured_var PP;
  captured_var G;
  captured_var K;
  captured_var K0;
  captured_var tau;
  captured_var tau0;
  emxArray_int32_T *r;
  emxArray_int32_T *r3;
  emxArray_int32_T *r4;
  emxArray_real_T *B;
  emxArray_real_T *Bbar;
  emxArray_real_T *Bbrev;
  emxArray_real_T *Bbrevp;
  emxArray_real_T *Bp;
  emxArray_real_T *Bpp;
  emxArray_real_T *F_;
  emxArray_real_T *F_ib_;
  emxArray_real_T *F_ibd_;
  emxArray_real_T *F_ij_;
  emxArray_real_T *F_ijd_;
  emxArray_real_T *bNbar;
  emxArray_real_T *b_bepsbar;
  emxArray_real_T *b_bepscomma;
  emxArray_real_T *b_y;
  emxArray_real_T *bepsbar;
  emxArray_real_T *bepscomma;
  emxArray_real_T *bepscomma_;
  emxArray_real_T *bepsdbar;
  emxArray_real_T *bepsdcomma;
  emxArray_real_T *bepsdcomma_;
  emxArray_real_T *cNbar;
  emxArray_real_T *mu_ab_;
  emxArray_real_T *mu_abd_;
  emxArray_real_T *mu_aj_;
  emxArray_real_T *r2;
  emxArray_real_T *y;
  double Dm22[9];
  double Em22[9];
  double R0_[9];
  double XX_contents[3];
  double mm[3];
  double rrd[3];
  double xi0p[3];
  double xip[3];
  const double *R0_data;
  const double *colmat_bar_data;
  const double *colmat_brev_data;
  const double *colmat_data;
  const double *nel_data;
  const double *varTheta_data;
  const double *varThetadot_data;
  const double *wg_data;
  double Dm11;
  double Nbar;
  double b_d;
  double ct;
  double st;
  double *B_data;
  double *Bbar_data;
  double *Bbrev_data;
  double *Bbrevp_data;
  double *Bp_data;
  double *Bpp_data;
  double *F__data;
  double *F_data;
  double *F_ib__data;
  double *F_ib_data;
  double *F_ibd__data;
  double *F_ibd_data;
  double *F_ij__data;
  double *F_ij_data;
  double *F_ijd__data;
  double *F_ijd_data;
  double *bNbar_data;
  double *b_bepsbar_data;
  double *b_y_data;
  double *bepsbar_data;
  double *bepscomma__data;
  double *bepscomma_data;
  double *bepsdbar_data;
  double *bepsdcomma__data;
  double *bepsdcomma_data;
  double *cNbar_data;
  double *mu_ab__data;
  double *mu_ab_data;
  double *mu_abd__data;
  double *mu_abd_data;
  double *mu_aj__data;
  double *mu_aj_data;
  double *mu_data;
  double *y_data;
  int N;
  int Nbrev;
  int b_i;
  int b_loop_ub;
  int b_loop_ub_tmp;
  int c_loop_ub;
  int c_loop_ub_tmp;
  int d_loop_ub;
  int d_loop_ub_tmp;
  int ddc;
  int e_loop_ub;
  int e_loop_ub_tmp;
  int f_loop_ub;
  int f_loop_ub_tmp;
  int g_loop_ub;
  int g_loop_ub_tmp;
  int h_loop_ub;
  int h_loop_ub_tmp;
  int i;
  int i1;
  int i2;
  int i3;
  int i4;
  int i5;
  int i6;
  int i7;
  int i8;
  int i_loop_ub;
  int j_loop_ub;
  int k_loop_ub;
  int l_loop_ub;
  int loop_ub;
  int loop_ub_tmp;
  int n;
  int pass;
  int *r1;
  int *r5;
  colmat_bar_data = colmat_bar->data;
  colmat_brev_data = colmat_brev->data;
  colmat_data = colmat->data;
  nel_data = nel->data;
  wg_data = wg->data;
  R0_data = R0->data;
  varThetadot_data = varThetadot->data;
  varTheta_data = varTheta->data;
  /*  INPUTS */
  /*  P = spline position DOF (represented as 3*N matrix), \mathscr{P} */
  /*  P0 = reference configuration */
  /*  Pd = \dot{P} */
  /*  varTheta = spline twist DOF (\vartheta) */
  /*  varThetdot = \dot{vartheta} */
  /*  R0 = cell array of orientations in the reference configuration at the */
  /*       quadrature points (+ the two end points) */
  /*  rho = mass per unit length */
  /*  EA = section axial rigidity */
  /*  EI = section bending rigidity */
  /*  GJ = section torsional rigidity */
  /*  betAX = damping coeff associated with rate of axial deformation */
  /*  betBEND = damping coeff associated with rate of bending deformation */
  /*  betTOR = damping coeff associated with rate of torsional deformation */
  /*  sg = quadrature points */
  /*  wg = quadrature weights */
  /*  nel = nel(n) = index of element to which quadrature pt sg(n) belongs to */
  /*  colmat = B-spline basis for position and its derivatives (B) */
  /*  colmat_brev = B-spline basis for twist and its derivative (Brev) */
  /*  colmat_bar = B-spline basis for strain projection (Bbar) */
  /*  d = degree of B-spline basis for position (B) */
  /*  dbrev = degree of B-spline basis for twist (Bbrev) */
  /*  dbar = degree of B-spline basis for strain projection (Bbar) */
  /*  Mbar = "mass matrix" for strain projection */
  /*  u = force per unit mass (acceleration due to gravity and ground)  */
  /*  Kbar11, Dbar11 = as defined in equation (33) [EQ NUM needs to be updated
   */
  /*                   to be consistent with document] */
  /*  OUTPUTS */
  /*  F = Force vector of length 3*N */
  /*  mu = Moment vector of length Nbrev */
  /*  F_ij, F_ib, F_ijd, F_ibd, mu_aj, mu_ab, mu_ajd, mu_abd = stiffness terms
   */
  /*  number of basis functions (or control points) */
  N = colmat->size[1];
  Nbrev = colmat_brev->size[1];
  Nbar = colmat_bar->size[1];
  /*  Constitutive matrices */
  mm[0] = GJ;
  mm[1] = EI;
  mm[2] = EI;
  diag(mm, Em22);
  Dm11 = betAX * EA;
  mm[0] = betTOR * GJ;
  ct = betBEND * EI;
  mm[1] = ct;
  mm[2] = ct;
  diag(mm, Dm22);
  /*  number of quadrature (or collocation) points  */
  for (i = 0; i < 9; i++) {
    I3.contents[i] = iv[i];
  }
  /*  Compute coefficients for the projected strain */
  emxInit_real_T(&bepsbar, 1);
  i = bepsbar->size[0];
  bepsbar->size[0] = colmat_bar->size[1];
  emxEnsureCapacity_real_T(bepsbar, i);
  bepsbar_data = bepsbar->data;
  loop_ub = colmat_bar->size[1];
  emxInit_real_T(&bepsdbar, 1);
  i = bepsdbar->size[0];
  bepsdbar->size[0] = colmat_bar->size[1];
  emxEnsureCapacity_real_T(bepsdbar, i);
  bepsdbar_data = bepsdbar->data;
  for (i = 0; i < loop_ub; i++) {
    bepsbar_data[i] = 0.0;
    bepsdbar_data[i] = 0.0;
  }
  i = wg->size[0];
  emxInit_real_T(&y, 2);
  y_data = y->data;
  if (wg->size[0] - 1 >= 0) {
    b_loop_ub = colmat->size[1];
    c_loop_ub = colmat_bar->size[1];
  }
  emxInit_real_T(&Bp, 1);
  emxInit_real_T(&Bbar, 1);
  emxInit_int32_T(&r, 2);
  emxInit_real_T(&b_y, 2);
  b_y_data = b_y->data;
  emxInit_real_T(&b_bepsbar, 2);
  b_bepsbar_data = b_bepsbar->data;
  for (n = 0; n < i; n++) {
    i1 = (int)(3.0 * (((double)n + 1.0) - 1.0) + 2.0);
    i2 = Bp->size[0];
    Bp->size[0] = colmat->size[1];
    emxEnsureCapacity_real_T(Bp, i2);
    Bp_data = Bp->data;
    for (i2 = 0; i2 < b_loop_ub; i2++) {
      Bp_data[i2] = colmat_data[(i1 + colmat->size[0] * i2) - 1];
    }
    i1 = Bbar->size[0];
    Bbar->size[0] = colmat_bar->size[1];
    emxEnsureCapacity_real_T(Bbar, i1);
    Bbar_data = Bbar->data;
    for (i1 = 0; i1 < c_loop_ub; i1++) {
      Bbar_data[i1] = colmat_bar_data[n + colmat_bar->size[0] * i1];
    }
    mtimes(P, Bp, xip);
    nxip.contents = b_norm(xip);
    if (dbar < 0.0) {
      b_y->size[1] = 0;
    } else {
      i1 = b_y->size[0] * b_y->size[1];
      b_y->size[0] = 1;
      b_y->size[1] = (int)dbar + 1;
      emxEnsureCapacity_real_T(b_y, i1);
      b_y_data = b_y->data;
      loop_ub = (int)dbar;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        b_y_data[i1] = i1;
      }
    }
    i1 = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    emxEnsureCapacity_real_T(b_y, i1);
    b_y_data = b_y->data;
    ct = nel_data[n];
    loop_ub = b_y->size[1] - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      b_y_data[i1] += ct;
    }
    mtimes(P0, Bp, rrd);
    st = (nxip.contents - b_norm(rrd)) * wg_data[n];
    if (dbar < 0.0) {
      y->size[1] = 0;
    } else {
      i1 = y->size[0] * y->size[1];
      y->size[0] = 1;
      y->size[1] = (int)dbar + 1;
      emxEnsureCapacity_real_T(y, i1);
      y_data = y->data;
      loop_ub = (int)dbar;
      for (i1 = 0; i1 <= loop_ub; i1++) {
        y_data[i1] = i1;
      }
    }
    i1 = y->size[0] * y->size[1];
    y->size[0] = 1;
    emxEnsureCapacity_real_T(y, i1);
    y_data = y->data;
    ct = nel_data[n];
    loop_ub = y->size[1] - 1;
    for (i1 = 0; i1 <= loop_ub; i1++) {
      y_data[i1] += ct;
    }
    i1 = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = y->size[1];
    emxEnsureCapacity_int32_T(r, i1);
    r1 = r->data;
    loop_ub = y->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r1[i1] = (int)y_data[i1];
    }
    i1 = b_bepsbar->size[0] * b_bepsbar->size[1];
    b_bepsbar->size[0] = 1;
    b_bepsbar->size[1] = r->size[1];
    emxEnsureCapacity_real_T(b_bepsbar, i1);
    b_bepsbar_data = b_bepsbar->data;
    loop_ub = r->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_d = b_y_data[i1];
      b_bepsbar_data[i1] =
          bepsbar_data[(int)b_d - 1] + Bbar_data[(int)b_d - 1] * st;
    }
    loop_ub = b_bepsbar->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      bepsbar_data[r1[i1] - 1] = b_bepsbar_data[i1];
    }
    mtimes(Pdot, Bp, xi0p);
    st =
        ((xip[0] / nxip.contents * xi0p[0] + xip[1] / nxip.contents * xi0p[1]) +
         xip[2] / nxip.contents * xi0p[2]) *
        wg_data[n];
    i1 = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = y->size[1];
    emxEnsureCapacity_int32_T(r, i1);
    r1 = r->data;
    loop_ub = y->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      r1[i1] = (int)y_data[i1];
    }
    i1 = b_bepsbar->size[0] * b_bepsbar->size[1];
    b_bepsbar->size[0] = 1;
    b_bepsbar->size[1] = r->size[1];
    emxEnsureCapacity_real_T(b_bepsbar, i1);
    b_bepsbar_data = b_bepsbar->data;
    loop_ub = r->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_d = b_y_data[i1];
      b_bepsbar_data[i1] =
          bepsdbar_data[(int)b_d - 1] + Bbar_data[(int)b_d - 1] * st;
    }
    loop_ub = b_bepsbar->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      bepsdbar_data[r1[i1] - 1] = b_bepsbar_data[i1];
    }
  }
  b_mldivide(Mbar, bepsbar);
  bepsbar_data = bepsbar->data;
  b_mldivide(Mbar, bepsdbar);
  bepsdbar_data = bepsdbar->data;
  /*  In the first pass, projection Nbar is computed, */
  /*  and in the second pass, force vector and stiffness are  */
  /*  assembled. The code is organized in this way because  */
  /*  there is significant repetition of code between the two */
  /*  passes */
  /* ---------------------------------------------------------------------- */
  /*  The following two allocations are to define size for code generation */
  emxInit_real_T(&bNbar, 1);
  i = bNbar->size[0];
  bNbar->size[0] = colmat_bar->size[1];
  emxEnsureCapacity_real_T(bNbar, i);
  bNbar_data = bNbar->data;
  loop_ub = colmat_bar->size[1];
  /*  allocate, computed in pass 1 */
  emxInit_real_T(&cNbar, 1);
  i = cNbar->size[0];
  cNbar->size[0] = colmat_bar->size[1];
  emxEnsureCapacity_real_T(cNbar, i);
  cNbar_data = cNbar->data;
  for (i = 0; i < loop_ub; i++) {
    bNbar_data[i] = 0.0;
    cNbar_data[i] = 0.0;
  }
  /*  allocate, computed in pass 1 */
  emxInit_real_T(&F_, 1);
  loop_ub = (int)(3.0 * (d + 1.0));
  i = F_->size[0];
  F_->size[0] = loop_ub;
  emxEnsureCapacity_real_T(F_, i);
  F__data = F_->data;
  for (i = 0; i < loop_ub; i++) {
    F__data[i] = 0.0;
  }
  /*  allocate, computed in pass 2 */
  i = F->size[0];
  F->size[0] = 3 * P->size[1];
  emxEnsureCapacity_real_T(F, i);
  F_data = F->data;
  loop_ub_tmp = 3 * P->size[1];
  for (i = 0; i < loop_ub_tmp; i++) {
    F_data[i] = 0.0;
  }
  /*  Generalized force */
  i = mu->size[0];
  mu->size[0] = varTheta->size[0];
  emxEnsureCapacity_real_T(mu, i);
  mu_data = mu->data;
  b_loop_ub = varTheta->size[0];
  for (i = 0; i < b_loop_ub; i++) {
    mu_data[i] = 0.0;
  }
  /*  Generalized moment */
  emxInit_real_T(&bepscomma, 2);
  i = bepscomma->size[0] * bepscomma->size[1];
  bepscomma->size[0] = colmat_bar->size[1];
  bepscomma->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(bepscomma, i);
  bepscomma_data = bepscomma->data;
  b_loop_ub_tmp = colmat_bar->size[1] * (3 * colmat->size[1]);
  for (i = 0; i < b_loop_ub_tmp; i++) {
    bepscomma_data[i] = 0.0;
  }
  emxInit_real_T(&bepscomma_, 2);
  i = (int)(dbar + 1.0);
  i1 = bepscomma_->size[0] * bepscomma_->size[1];
  bepscomma_->size[0] = (int)(dbar + 1.0);
  bepscomma_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(bepscomma_, i1);
  bepscomma__data = bepscomma_->data;
  c_loop_ub_tmp = (int)(dbar + 1.0) * loop_ub;
  for (i1 = 0; i1 < c_loop_ub_tmp; i1++) {
    bepscomma__data[i1] = 0.0;
  }
  /*  temp storage */
  emxInit_real_T(&bepsdcomma, 2);
  i1 = bepsdcomma->size[0] * bepsdcomma->size[1];
  bepsdcomma->size[0] = colmat_bar->size[1];
  bepsdcomma->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(bepsdcomma, i1);
  bepsdcomma_data = bepsdcomma->data;
  for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
    bepsdcomma_data[i1] = 0.0;
  }
  emxInit_real_T(&bepsdcomma_, 2);
  i1 = bepsdcomma_->size[0] * bepsdcomma_->size[1];
  bepsdcomma_->size[0] = (int)(dbar + 1.0);
  bepsdcomma_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(bepsdcomma_, i1);
  bepsdcomma__data = bepsdcomma_->data;
  for (i1 = 0; i1 < c_loop_ub_tmp; i1++) {
    bepsdcomma__data[i1] = 0.0;
  }
  /*  temp storage */
  i1 = F_ij->size[0] * F_ij->size[1];
  F_ij->size[0] = 3 * colmat->size[1];
  F_ij->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(F_ij, i1);
  F_ij_data = F_ij->data;
  b_loop_ub_tmp = 3 * colmat->size[1] * (3 * colmat->size[1]);
  for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
    F_ij_data[i1] = 0.0;
  }
  emxInit_real_T(&F_ij_, 2);
  i1 = F_ij_->size[0] * F_ij_->size[1];
  F_ij_->size[0] = loop_ub;
  F_ij_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(F_ij_, i1);
  F_ij__data = F_ij_->data;
  d_loop_ub_tmp = loop_ub * loop_ub;
  for (i1 = 0; i1 < d_loop_ub_tmp; i1++) {
    F_ij__data[i1] = 0.0;
  }
  /*  temp storage */
  i1 = F_ib->size[0] * F_ib->size[1];
  F_ib->size[0] = 3 * colmat->size[1];
  F_ib->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(F_ib, i1);
  F_ib_data = F_ib->data;
  e_loop_ub_tmp = 3 * colmat->size[1] * colmat_brev->size[1];
  for (i1 = 0; i1 < e_loop_ub_tmp; i1++) {
    F_ib_data[i1] = 0.0;
  }
  emxInit_real_T(&F_ib_, 2);
  i1 = F_ib_->size[0] * F_ib_->size[1];
  F_ib_->size[0] = loop_ub;
  i2 = (int)(dbrev + 1.0);
  F_ib_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(F_ib_, i1);
  F_ib__data = F_ib_->data;
  f_loop_ub_tmp = loop_ub * (int)(dbrev + 1.0);
  for (i1 = 0; i1 < f_loop_ub_tmp; i1++) {
    F_ib__data[i1] = 0.0;
  }
  /*  temp storage */
  i1 = F_ijd->size[0] * F_ijd->size[1];
  F_ijd->size[0] = 3 * colmat->size[1];
  F_ijd->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(F_ijd, i1);
  F_ijd_data = F_ijd->data;
  for (i1 = 0; i1 < b_loop_ub_tmp; i1++) {
    F_ijd_data[i1] = 0.0;
  }
  emxInit_real_T(&F_ijd_, 2);
  i1 = F_ijd_->size[0] * F_ijd_->size[1];
  F_ijd_->size[0] = loop_ub;
  F_ijd_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(F_ijd_, i1);
  F_ijd__data = F_ijd_->data;
  for (i1 = 0; i1 < d_loop_ub_tmp; i1++) {
    F_ijd__data[i1] = 0.0;
  }
  /*  temp storage */
  i1 = F_ibd->size[0] * F_ibd->size[1];
  F_ibd->size[0] = 3 * colmat->size[1];
  F_ibd->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(F_ibd, i1);
  F_ibd_data = F_ibd->data;
  for (i1 = 0; i1 < e_loop_ub_tmp; i1++) {
    F_ibd_data[i1] = 0.0;
  }
  emxInit_real_T(&F_ibd_, 2);
  i1 = F_ibd_->size[0] * F_ibd_->size[1];
  F_ibd_->size[0] = loop_ub;
  F_ibd_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(F_ibd_, i1);
  F_ibd__data = F_ibd_->data;
  for (i1 = 0; i1 < f_loop_ub_tmp; i1++) {
    F_ibd__data[i1] = 0.0;
  }
  /*  temp storage */
  i1 = mu_aj->size[0] * mu_aj->size[1];
  mu_aj->size[0] = colmat_brev->size[1];
  mu_aj->size[1] = 3 * colmat->size[1];
  emxEnsureCapacity_real_T(mu_aj, i1);
  mu_aj_data = mu_aj->data;
  for (i1 = 0; i1 < e_loop_ub_tmp; i1++) {
    mu_aj_data[i1] = 0.0;
  }
  emxInit_real_T(&mu_aj_, 2);
  i1 = mu_aj_->size[0] * mu_aj_->size[1];
  mu_aj_->size[0] = (int)(dbrev + 1.0);
  mu_aj_->size[1] = loop_ub;
  emxEnsureCapacity_real_T(mu_aj_, i1);
  mu_aj__data = mu_aj_->data;
  for (i1 = 0; i1 < f_loop_ub_tmp; i1++) {
    mu_aj__data[i1] = 0.0;
  }
  /*  temp storage */
  i1 = mu_ab->size[0] * mu_ab->size[1];
  mu_ab->size[0] = colmat_brev->size[1];
  mu_ab->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(mu_ab, i1);
  mu_ab_data = mu_ab->data;
  g_loop_ub_tmp = colmat_brev->size[1] * colmat_brev->size[1];
  for (i1 = 0; i1 < g_loop_ub_tmp; i1++) {
    mu_ab_data[i1] = 0.0;
  }
  emxInit_real_T(&mu_ab_, 2);
  i1 = mu_ab_->size[0] * mu_ab_->size[1];
  mu_ab_->size[0] = (int)(dbrev + 1.0);
  mu_ab_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(mu_ab_, i1);
  mu_ab__data = mu_ab_->data;
  h_loop_ub_tmp = (int)(dbrev + 1.0) * (int)(dbrev + 1.0);
  for (i1 = 0; i1 < h_loop_ub_tmp; i1++) {
    mu_ab__data[i1] = 0.0;
  }
  /*  temp storage */
  i1 = mu_abd->size[0] * mu_abd->size[1];
  mu_abd->size[0] = colmat_brev->size[1];
  mu_abd->size[1] = colmat_brev->size[1];
  emxEnsureCapacity_real_T(mu_abd, i1);
  mu_abd_data = mu_abd->data;
  for (i1 = 0; i1 < g_loop_ub_tmp; i1++) {
    mu_abd_data[i1] = 0.0;
  }
  emxInit_real_T(&mu_abd_, 2);
  i1 = mu_abd_->size[0] * mu_abd_->size[1];
  mu_abd_->size[0] = (int)(dbrev + 1.0);
  mu_abd_->size[1] = (int)(dbrev + 1.0);
  emxEnsureCapacity_real_T(mu_abd_, i1);
  mu_abd__data = mu_abd_->data;
  for (i1 = 0; i1 < h_loop_ub_tmp; i1++) {
    mu_abd__data[i1] = 0.0;
  }
  /*  temp storage */
  /* ---------------------------------------------------------------------- */
  i1 = wg->size[0];
  if (wg->size[0] - 1 >= 0) {
    d_loop_ub = colmat->size[1];
    e_loop_ub = colmat_brev->size[1];
    f_loop_ub = colmat_bar->size[1];
    L_KII.contents[0] = -0.0;
    L_KII.contents[4] = -0.0;
    L_KII.contents[8] = -0.0;
    g_loop_ub = varTheta->size[0];
    h_loop_ub = varThetadot->size[0];
    i_loop_ub = varTheta->size[0];
    j_loop_ub = varThetadot->size[0];
    k_loop_ub = bepsbar->size[0];
    l_loop_ub = bepsdbar->size[0];
  }
  emxInit_real_T(&B, 1);
  emxInit_real_T(&Bpp, 1);
  emxInit_real_T(&Bbrev, 1);
  emxInit_real_T(&Bbrevp, 1);
  emxInit_real_T(&r2, 2);
  emxInit_int32_T(&r3, 1);
  emxInit_int32_T(&r4, 1);
  emxInit_real_T(&b_bepscomma, 2);
  for (pass = 0; pass < 2; pass++) {
    if (pass + 1 == 1) {
      c_loop_ub = (int)Nbar;
      i3 = bNbar->size[0];
      bNbar->size[0] = (int)Nbar;
      emxEnsureCapacity_real_T(bNbar, i3);
      bNbar_data = bNbar->data;
      for (i3 = 0; i3 < c_loop_ub; i3++) {
        bNbar_data[i3] = 0.0;
      }
    } else {
      /*  pass == 2 */
      i3 = F->size[0];
      F->size[0] = 3 * P->size[1];
      emxEnsureCapacity_real_T(F, i3);
      F_data = F->data;
      for (i3 = 0; i3 < loop_ub_tmp; i3++) {
        F_data[i3] = 0.0;
      }
      /*  Generalized force */
      i3 = F_->size[0];
      F_->size[0] = loop_ub;
      emxEnsureCapacity_real_T(F_, i3);
      F__data = F_->data;
      for (i3 = 0; i3 < loop_ub; i3++) {
        F__data[i3] = 0.0;
      }
      /*  temp storage when computing F */
      b_loop_ub = varTheta->size[0];
      i3 = mu->size[0];
      mu->size[0] = varTheta->size[0];
      emxEnsureCapacity_real_T(mu, i3);
      mu_data = mu->data;
      for (i3 = 0; i3 < b_loop_ub; i3++) {
        mu_data[i3] = 0.0;
      }
      /*  Generalized moment */
      i3 = bepscomma->size[0] * bepscomma->size[1];
      bepscomma->size[0] = (int)Nbar;
      bepscomma->size[1] = 3 * N;
      emxEnsureCapacity_real_T(bepscomma, i3);
      bepscomma_data = bepscomma->data;
      c_loop_ub = (int)Nbar * (3 * N);
      for (i3 = 0; i3 < c_loop_ub; i3++) {
        bepscomma_data[i3] = 0.0;
      }
      i3 = bepscomma_->size[0] * bepscomma_->size[1];
      bepscomma_->size[0] = (int)(dbar + 1.0);
      bepscomma_->size[1] = loop_ub;
      emxEnsureCapacity_real_T(bepscomma_, i3);
      bepscomma__data = bepscomma_->data;
      for (i3 = 0; i3 < c_loop_ub_tmp; i3++) {
        bepscomma__data[i3] = 0.0;
      }
      /*  temp storage */
      i3 = bepsdcomma->size[0] * bepsdcomma->size[1];
      bepsdcomma->size[0] = (int)Nbar;
      bepsdcomma->size[1] = 3 * N;
      emxEnsureCapacity_real_T(bepsdcomma, i3);
      bepsdcomma_data = bepsdcomma->data;
      for (i3 = 0; i3 < c_loop_ub; i3++) {
        bepsdcomma_data[i3] = 0.0;
      }
      i3 = bepsdcomma_->size[0] * bepsdcomma_->size[1];
      bepsdcomma_->size[0] = (int)(dbar + 1.0);
      bepsdcomma_->size[1] = loop_ub;
      emxEnsureCapacity_real_T(bepsdcomma_, i3);
      bepsdcomma__data = bepsdcomma_->data;
      for (i3 = 0; i3 < c_loop_ub_tmp; i3++) {
        bepsdcomma__data[i3] = 0.0;
      }
      /*  temp storage */
      i3 = F_ij->size[0] * F_ij->size[1];
      F_ij->size[0] = 3 * N;
      F_ij->size[1] = 3 * N;
      emxEnsureCapacity_real_T(F_ij, i3);
      F_ij_data = F_ij->data;
      for (i3 = 0; i3 < b_loop_ub_tmp; i3++) {
        F_ij_data[i3] = 0.0;
      }
      i3 = F_ij_->size[0] * F_ij_->size[1];
      F_ij_->size[0] = loop_ub;
      F_ij_->size[1] = loop_ub;
      emxEnsureCapacity_real_T(F_ij_, i3);
      F_ij__data = F_ij_->data;
      for (i3 = 0; i3 < d_loop_ub_tmp; i3++) {
        F_ij__data[i3] = 0.0;
      }
      /*  temp storage */
      i3 = F_ib->size[0] * F_ib->size[1];
      F_ib->size[0] = 3 * N;
      F_ib->size[1] = Nbrev;
      emxEnsureCapacity_real_T(F_ib, i3);
      F_ib_data = F_ib->data;
      for (i3 = 0; i3 < e_loop_ub_tmp; i3++) {
        F_ib_data[i3] = 0.0;
      }
      i3 = F_ib_->size[0] * F_ib_->size[1];
      F_ib_->size[0] = loop_ub;
      F_ib_->size[1] = (int)(dbrev + 1.0);
      emxEnsureCapacity_real_T(F_ib_, i3);
      F_ib__data = F_ib_->data;
      for (i3 = 0; i3 < f_loop_ub_tmp; i3++) {
        F_ib__data[i3] = 0.0;
      }
      /*  temp storage */
      i3 = F_ijd->size[0] * F_ijd->size[1];
      F_ijd->size[0] = 3 * N;
      F_ijd->size[1] = 3 * N;
      emxEnsureCapacity_real_T(F_ijd, i3);
      F_ijd_data = F_ijd->data;
      for (i3 = 0; i3 < b_loop_ub_tmp; i3++) {
        F_ijd_data[i3] = 0.0;
      }
      i3 = F_ijd_->size[0] * F_ijd_->size[1];
      F_ijd_->size[0] = loop_ub;
      F_ijd_->size[1] = loop_ub;
      emxEnsureCapacity_real_T(F_ijd_, i3);
      F_ijd__data = F_ijd_->data;
      for (i3 = 0; i3 < d_loop_ub_tmp; i3++) {
        F_ijd__data[i3] = 0.0;
      }
      /*  temp storage */
      i3 = F_ibd->size[0] * F_ibd->size[1];
      F_ibd->size[0] = 3 * N;
      F_ibd->size[1] = Nbrev;
      emxEnsureCapacity_real_T(F_ibd, i3);
      F_ibd_data = F_ibd->data;
      for (i3 = 0; i3 < e_loop_ub_tmp; i3++) {
        F_ibd_data[i3] = 0.0;
      }
      i3 = F_ibd_->size[0] * F_ibd_->size[1];
      F_ibd_->size[0] = loop_ub;
      F_ibd_->size[1] = (int)(dbrev + 1.0);
      emxEnsureCapacity_real_T(F_ibd_, i3);
      F_ibd__data = F_ibd_->data;
      for (i3 = 0; i3 < f_loop_ub_tmp; i3++) {
        F_ibd__data[i3] = 0.0;
      }
      /*  temp storage */
      i3 = mu_aj->size[0] * mu_aj->size[1];
      mu_aj->size[0] = Nbrev;
      mu_aj->size[1] = 3 * N;
      emxEnsureCapacity_real_T(mu_aj, i3);
      mu_aj_data = mu_aj->data;
      for (i3 = 0; i3 < e_loop_ub_tmp; i3++) {
        mu_aj_data[i3] = 0.0;
      }
      i3 = mu_aj_->size[0] * mu_aj_->size[1];
      mu_aj_->size[0] = (int)(dbrev + 1.0);
      mu_aj_->size[1] = loop_ub;
      emxEnsureCapacity_real_T(mu_aj_, i3);
      mu_aj__data = mu_aj_->data;
      for (i3 = 0; i3 < f_loop_ub_tmp; i3++) {
        mu_aj__data[i3] = 0.0;
      }
      /*  temp storage */
      i3 = mu_ab->size[0] * mu_ab->size[1];
      mu_ab->size[0] = Nbrev;
      mu_ab->size[1] = Nbrev;
      emxEnsureCapacity_real_T(mu_ab, i3);
      mu_ab_data = mu_ab->data;
      for (i3 = 0; i3 < g_loop_ub_tmp; i3++) {
        mu_ab_data[i3] = 0.0;
      }
      i3 = mu_ab_->size[0] * mu_ab_->size[1];
      mu_ab_->size[0] = (int)(dbrev + 1.0);
      mu_ab_->size[1] = (int)(dbrev + 1.0);
      emxEnsureCapacity_real_T(mu_ab_, i3);
      mu_ab__data = mu_ab_->data;
      for (i3 = 0; i3 < h_loop_ub_tmp; i3++) {
        mu_ab__data[i3] = 0.0;
      }
      /*  temp storage */
      i3 = mu_abd->size[0] * mu_abd->size[1];
      mu_abd->size[0] = Nbrev;
      mu_abd->size[1] = Nbrev;
      emxEnsureCapacity_real_T(mu_abd, i3);
      mu_abd_data = mu_abd->data;
      for (i3 = 0; i3 < g_loop_ub_tmp; i3++) {
        mu_abd_data[i3] = 0.0;
      }
      i3 = mu_abd_->size[0] * mu_abd_->size[1];
      mu_abd_->size[0] = (int)(dbrev + 1.0);
      mu_abd_->size[1] = (int)(dbrev + 1.0);
      emxEnsureCapacity_real_T(mu_abd_, i3);
      mu_abd__data = mu_abd_->data;
      for (i3 = 0; i3 < h_loop_ub_tmp; i3++) {
        mu_abd__data[i3] = 0.0;
      }
      /*  temp storage */
    }
    for (n = 0; n < i1; n++) {
      double D22UU[9];
      double L_rrI[9];
      double L_rrII[9];
      double Q2[9];
      double Theta[9];
      double ThetaR0[9];
      double b_y_tmp[9];
      double htau0[9];
      double y_tmp[9];
      double D22WW[3];
      double Q20[3];
      double b_L_chiO[3];
      double rr[3];
      double xi0pp[3];
      double xipd[3];
      double xippd[3];
      double D22WW_tmp;
      double L_chiO;
      double L_chiO_tmp;
      double a;
      double b_D22WW_tmp;
      double b_L_chiO_tmp;
      double b_tau0;
      double b_varTheta;
      double b_varThetadot;
      double c_D22WW_tmp;
      double c_L_chiO_tmp;
      double c_tmp;
      double d1;
      double d2;
      double d3;
      double d4;
      double d5;
      double d6;
      double nel_;
      double nxi0p;
      double tau0_tmp;
      nel_ = nel_data[n];
      b_d = 3.0 * (((double)n + 1.0) - 1.0);
      i3 = B->size[0];
      B->size[0] = colmat->size[1];
      emxEnsureCapacity_real_T(B, i3);
      B_data = B->data;
      i3 = Bp->size[0];
      Bp->size[0] = colmat->size[1];
      emxEnsureCapacity_real_T(Bp, i3);
      Bp_data = Bp->data;
      i3 = Bpp->size[0];
      Bpp->size[0] = colmat->size[1];
      emxEnsureCapacity_real_T(Bpp, i3);
      Bpp_data = Bpp->data;
      for (i3 = 0; i3 < d_loop_ub; i3++) {
        B_data[i3] = colmat_data[((int)(b_d + 1.0) + colmat->size[0] * i3) - 1];
        Bp_data[i3] =
            colmat_data[((int)(b_d + 2.0) + colmat->size[0] * i3) - 1];
        Bpp_data[i3] =
            colmat_data[((int)(b_d + 3.0) + colmat->size[0] * i3) - 1];
      }
      i3 = (int)((unsigned int)n << 1);
      i4 = Bbrev->size[0];
      Bbrev->size[0] = colmat_brev->size[1];
      emxEnsureCapacity_real_T(Bbrev, i4);
      Bbrev_data = Bbrev->data;
      i4 = Bbrevp->size[0];
      Bbrevp->size[0] = colmat_brev->size[1];
      emxEnsureCapacity_real_T(Bbrevp, i4);
      Bbrevp_data = Bbrevp->data;
      for (i4 = 0; i4 < e_loop_ub; i4++) {
        Bbrev_data[i4] = colmat_brev_data[i3 + colmat_brev->size[0] * i4];
        Bbrevp_data[i4] =
            colmat_brev_data[(i3 + colmat_brev->size[0] * i4) + 1];
      }
      i3 = Bbar->size[0];
      Bbar->size[0] = colmat_bar->size[1];
      emxEnsureCapacity_real_T(Bbar, i3);
      Bbar_data = Bbar->data;
      for (i3 = 0; i3 < f_loop_ub; i3++) {
        Bbar_data[i3] = colmat_bar_data[n + colmat_bar->size[0] * i3];
      }
      mtimes(P0, Bp, xi0p);
      mtimes(P0, Bpp, xi0pp);
      nxi0p = b_norm(xi0p);
      tau0.contents[0] = xi0p[0] / nxi0p;
      tau0.contents[1] = xi0p[1] / nxi0p;
      tau0.contents[2] = xi0p[2] / nxi0p;
      ct = nxi0p * nxi0p;
      K0.contents[0] = (xi0p[1] * xi0pp[2] - xi0pp[1] * xi0p[2]) / ct;
      K0.contents[1] = (xi0pp[0] * xi0p[2] - xi0p[0] * xi0pp[2]) / ct;
      K0.contents[2] = (xi0p[0] * xi0pp[1] - xi0pp[0] * xi0p[1]) / ct;
      mtimes(P, Bp, xip);
      mtimes(Pdot, Bp, xipd);
      mtimes(P, Bpp, xi0p);
      mtimes(Pdot, Bpp, xippd);
      nxip.contents = b_norm(xip);
      c_tmp = nxip.contents * nxip.contents;
      K.contents[0] = (xip[1] * xi0p[2] - xi0p[1] * xip[2]) / c_tmp;
      K.contents[1] = (xi0p[0] * xip[2] - xip[0] * xi0p[2]) / c_tmp;
      K.contents[2] = (xip[0] * xi0p[1] - xi0p[0] * xip[1]) / c_tmp;
      b_d = xip[0] / nxip.contents;
      tau.contents[0] = b_d;
      ct = K0.contents[0] * b_d;
      b_d = xip[1] / nxip.contents;
      tau.contents[1] = b_d;
      ct += K0.contents[1] * b_d;
      b_d = xip[2] / nxip.contents;
      tau.contents[2] = b_d;
      b_tau0 = (tau0.contents[0] * K.contents[0] +
                tau0.contents[1] * K.contents[1]) +
               tau0.contents[2] * K.contents[2];
      ct += K0.contents[2] * b_d;
      G.contents[0] = b_tau0 * (2.0 * tau0.contents[0] + tau.contents[0]) -
                      ct * tau0.contents[0];
      G.contents[1] = b_tau0 * (2.0 * tau0.contents[1] + tau.contents[1]) -
                      ct * tau0.contents[1];
      G.contents[2] =
          b_tau0 * (2.0 * tau0.contents[2] + b_d) - ct * tau0.contents[2];
      tau0_tmp = (tau0.contents[0] * tau.contents[0] +
                  tau0.contents[1] * tau.contents[1]) +
                 tau0.contents[2] * b_d;
      for (b_i = 0; b_i < 3; b_i++) {
        rr[b_i] = (K.contents[b_i] - K0.contents[b_i]) -
                  G.contents[b_i] / (tau0_tmp + 1.0);
        PP.contents[3 * b_i] =
            I3.contents[3 * b_i] - tau.contents[0] * tau.contents[b_i];
        i3 = 3 * b_i + 1;
        PP.contents[i3] = I3.contents[i3] - tau.contents[1] * tau.contents[b_i];
        i3 = 3 * b_i + 2;
        PP.contents[i3] = I3.contents[i3] - b_d * tau.contents[b_i];
      }
      for (i3 = 0; i3 < 9; i3++) {
        L_tauI.contents[i3] = PP.contents[i3] / nxip.contents;
      }
      a = -2.0 / nxip.contents;
      d1 = tau.contents[0];
      d2 = tau.contents[1];
      for (i3 = 0; i3 < 3; i3++) {
        xi0p[i3] /= c_tmp;
        d3 = K.contents[i3];
        Q2[3 * i3] = a * d1 * d3;
        Q2[3 * i3 + 1] = a * d2 * d3;
        Q2[3 * i3 + 2] = a * b_d * d3;
      }
      D22UU[0] = 0.0;
      D22UU[3] = -xi0p[2];
      D22UU[6] = xi0p[1];
      D22UU[1] = xi0p[2];
      D22UU[4] = 0.0;
      D22UU[7] = -xi0p[0];
      D22UU[2] = -xi0p[1];
      D22UU[5] = xi0p[0];
      D22UU[8] = 0.0;
      for (i3 = 0; i3 < 9; i3++) {
        L_KI.contents[i3] = Q2[i3] + D22UU[i3];
      }
      b_tau0 = 0.0;
      for (i3 = 0; i3 < 3; i3++) {
        b_tau0 += tau0.contents[i3] * K.contents[i3];
        Q20[i3] = (L_KI.contents[i3] * tau0.contents[0] +
                   L_KI.contents[i3 + 3] * tau0.contents[1]) +
                  L_KI.contents[i3 + 6] * tau0.contents[2];
        mm[i3] = 2.0 * tau0.contents[i3] + tau.contents[i3];
      }
      d1 = Q20[0];
      d2 = Q20[1];
      d3 = Q20[2];
      L_chiO = K0.contents[0];
      d4 = K0.contents[1];
      d5 = K0.contents[2];
      for (i3 = 0; i3 < 3; i3++) {
        d6 = mm[i3];
        Q2[3 * i3] = d1 * d6 + b_tau0 * L_tauI.contents[3 * i3];
        c_loop_ub = 3 * i3 + 1;
        Q2[c_loop_ub] = d2 * d6 + b_tau0 * L_tauI.contents[c_loop_ub];
        c_loop_ub = 3 * i3 + 2;
        Q2[c_loop_ub] = d3 * d6 + b_tau0 * L_tauI.contents[c_loop_ub];
        D22WW[i3] =
            (L_tauI.contents[i3] * L_chiO + L_tauI.contents[i3 + 3] * d4) +
            L_tauI.contents[i3 + 6] * d5;
      }
      d1 = D22WW[0];
      d2 = D22WW[1];
      d3 = D22WW[2];
      for (i3 = 0; i3 < 3; i3++) {
        L_chiO = tau0.contents[i3];
        D22UU[3 * i3] = d1 * L_chiO;
        D22UU[3 * i3 + 1] = d2 * L_chiO;
        D22UU[3 * i3 + 2] = d3 * L_chiO;
      }
      for (i3 = 0; i3 < 9; i3++) {
        L_GI.contents[i3] = Q2[i3] - D22UU[i3];
      }
      ct = (tau0_tmp + 1.0) * (tau0_tmp + 1.0);
      d1 = tau0.contents[0];
      d2 = tau0.contents[1];
      d3 = tau0.contents[2];
      for (i3 = 0; i3 < 3; i3++) {
        D22WW[i3] = ((L_tauI.contents[i3] * d1 + L_tauI.contents[i3 + 3] * d2) +
                     L_tauI.contents[i3 + 6] * d3) /
                    ct;
      }
      d1 = D22WW[0];
      d2 = D22WW[1];
      d3 = D22WW[2];
      for (i3 = 0; i3 < 3; i3++) {
        L_chiO = G.contents[i3];
        L_rrI[3 * i3] = (L_KI.contents[3 * i3] + d1 * L_chiO) -
                        L_GI.contents[3 * i3] / (tau0_tmp + 1.0);
        c_loop_ub = 3 * i3 + 1;
        L_rrI[c_loop_ub] = (L_KI.contents[c_loop_ub] + d2 * L_chiO) -
                           L_GI.contents[c_loop_ub] / (tau0_tmp + 1.0);
        c_loop_ub = 3 * i3 + 2;
        L_rrI[c_loop_ub] = (L_KI.contents[c_loop_ub] + d3 * L_chiO) -
                           L_GI.contents[c_loop_ub] / (tau0_tmp + 1.0);
        xip[i3] /= c_tmp;
      }
      L_KII.contents[3] = xip[2];
      L_KII.contents[6] = -xip[1];
      L_KII.contents[1] = -xip[2];
      L_KII.contents[7] = xip[0];
      L_KII.contents[2] = xip[1];
      L_KII.contents[5] = -xip[0];
      for (i3 = 0; i3 < 3; i3++) {
        Q20[i3] = (L_KII.contents[i3] * tau0.contents[0] +
                   L_KII.contents[i3 + 3] * tau0.contents[1]) +
                  L_KII.contents[i3 + 6] * tau0.contents[2];
        mm[i3] = 2.0 * tau0.contents[i3] + tau.contents[i3];
      }
      b_tau0 = 0.0;
      d1 = Q20[0];
      d2 = Q20[1];
      d3 = Q20[2];
      for (i3 = 0; i3 < 3; i3++) {
        L_chiO = mm[i3];
        L_GII.contents[3 * i3] = d1 * L_chiO;
        L_GII.contents[3 * i3 + 1] = d2 * L_chiO;
        L_GII.contents[3 * i3 + 2] = d3 * L_chiO;
        b_tau0 += tau0.contents[i3] * tau.contents[i3];
      }
      for (i3 = 0; i3 < 9; i3++) {
        L_rrII[i3] = L_KII.contents[i3] - L_GII.contents[i3] / (b_tau0 + 1.0);
      }
      ct = 0.0;
      for (i3 = 0; i3 < g_loop_ub; i3++) {
        ct += varTheta_data[i3] * Bbrev_data[i3];
      }
      b_varThetadot = 0.0;
      for (i3 = 0; i3 < h_loop_ub; i3++) {
        b_varThetadot += varThetadot_data[i3] * Bbrev_data[i3];
      }
      b_varTheta = 0.0;
      for (i3 = 0; i3 < i_loop_ub; i3++) {
        b_varTheta += varTheta_data[i3] * Bbrevp_data[i3];
      }
      c_tmp = 0.0;
      for (i3 = 0; i3 < j_loop_ub; i3++) {
        c_tmp += varThetadot_data[i3] * Bbrevp_data[i3];
      }
      st = sin(ct);
      ct = cos(ct);
      htau0[0] = 0.0;
      htau0[3] = -tau0.contents[2];
      htau0[6] = tau0.contents[1];
      htau0[1] = tau0.contents[2];
      htau0[4] = 0.0;
      htau0[7] = -tau0.contents[0];
      htau0[2] = -tau0.contents[1];
      htau0[5] = tau0.contents[0];
      htau0[8] = 0.0;
      b_tau0 = (tau0.contents[0] * xi0pp[0] + tau0.contents[1] * xi0pp[1]) +
               tau0.contents[2] * xi0pp[2];
      d1 = 3.0 * ((double)n + 1.0) + 1.0;
      for (i3 = 0; i3 < 3; i3++) {
        xi0pp[i3] = (xi0pp[i3] - b_tau0 * tau0.contents[i3]) / nxi0p;
        for (i4 = 0; i4 < 3; i4++) {
          R0_[i4 + 3 * i3] = R0_data[i4 + 3 * ((int)(d1 + (double)i3) - 1)];
          Theta[i3 + 3 * i4] =
              (htau0[i3] * htau0[3 * i4] + htau0[i3 + 3] * htau0[3 * i4 + 1]) +
              htau0[i3 + 6] * htau0[3 * i4 + 2];
        }
      }
      for (i3 = 0; i3 < 9; i3++) {
        Theta[i3] = (I3.contents[i3] + st * htau0[i3]) + (1.0 - ct) * Theta[i3];
      }
      L_chiO_tmp = rr[1] * tau0.contents[2] - tau0.contents[1] * rr[2];
      b_L_chiO_tmp = tau0.contents[0] * rr[2] - rr[0] * tau0.contents[2];
      c_L_chiO_tmp = rr[0] * tau0.contents[1] - tau0.contents[0] * rr[1];
      b_L_chiO[0] = L_chiO_tmp + xi0pp[0];
      b_L_chiO[1] = b_L_chiO_tmp + xi0pp[1];
      b_L_chiO[2] = c_L_chiO_tmp + xi0pp[2];
      /*  calculate epsilon, epsilond, kappa, kappad */
      a = 0.0;
      for (i3 = 0; i3 < k_loop_ub; i3++) {
        a += bepsbar_data[i3] * Bbar_data[i3];
      }
      D22WW_tmp = tau0.contents[1] * xi0pp[2] - xi0pp[1] * tau0.contents[2];
      D22WW[0] = (1.0 - ct) * D22WW_tmp;
      b_D22WW_tmp = xi0pp[0] * tau0.contents[2] - tau0.contents[0] * xi0pp[2];
      D22WW[1] = (1.0 - ct) * b_D22WW_tmp;
      c_D22WW_tmp = tau0.contents[0] * xi0pp[1] - xi0pp[0] * tau0.contents[1];
      D22WW[2] = (1.0 - ct) * c_D22WW_tmp;
      d1 = rr[0];
      d2 = rr[1];
      d3 = rr[2];
      for (i3 = 0; i3 < 3; i3++) {
        Q20[i3] =
            ((Theta[3 * i3] * d1 + Theta[3 * i3 + 1] * d2) +
             Theta[3 * i3 + 2] * d3) +
            ((b_varTheta * tau0.contents[i3] + st * xi0pp[i3]) - D22WW[i3]);
      }
      d1 = Q20[0];
      d2 = Q20[1];
      d3 = Q20[2];
      for (i3 = 0; i3 < 3; i3++) {
        xi0p[i3] =
            ((R0_[3 * i3] * d1 + R0_[3 * i3 + 1] * d2) + R0_[3 * i3 + 2] * d3) /
            nxi0p;
      }
      tau0_tmp = 0.0;
      for (i3 = 0; i3 < l_loop_ub; i3++) {
        tau0_tmp += bepsdbar_data[i3] * Bbar_data[i3];
      }
      for (i3 = 0; i3 < 3; i3++) {
        y_tmp[3 * i3] = L_rrI[i3];
        b_y_tmp[3 * i3] = L_rrII[i3];
        c_loop_ub = 3 * i3 + 1;
        y_tmp[c_loop_ub] = L_rrI[i3 + 3];
        b_y_tmp[c_loop_ub] = L_rrII[i3 + 3];
        c_loop_ub = 3 * i3 + 2;
        y_tmp[c_loop_ub] = L_rrI[i3 + 6];
        b_y_tmp[c_loop_ub] = L_rrII[i3 + 6];
      }
      b_mtimes(R0_, Theta, ThetaR0);
      d1 = xipd[0];
      d2 = xipd[1];
      d3 = xipd[2];
      L_chiO = xippd[0];
      d4 = xippd[1];
      d5 = xippd[2];
      for (i3 = 0; i3 < 3; i3++) {
        d6 = ((y_tmp[i3] * d1 + y_tmp[i3 + 3] * d2) + y_tmp[i3 + 6] * d3) +
             ((b_y_tmp[i3] * L_chiO + b_y_tmp[i3 + 3] * d4) +
              b_y_tmp[i3 + 6] * d5);
        rrd[i3] = d6;
        Q20[i3] =
            (d6 + b_L_chiO[i3] * b_varThetadot) + tau0.contents[i3] * c_tmp;
      }
      d1 = Q20[0];
      d2 = Q20[1];
      d3 = Q20[2];
      for (i3 = 0; i3 < 3; i3++) {
        xip[i3] =
            ((ThetaR0[i3] * d1 + ThetaR0[i3 + 3] * d2) + ThetaR0[i3 + 6] * d3) /
            nxi0p;
      }
      if (pass + 1 == 1) {
        st = EA * a + Dm11 * tau0_tmp;
        if (dbar < 0.0) {
          b_y->size[0] = 1;
          b_y->size[1] = 0;
          b_bepsbar->size[0] = 1;
          b_bepsbar->size[1] = 0;
        } else {
          i3 = b_y->size[0] * b_y->size[1];
          b_y->size[0] = 1;
          b_y->size[1] = (int)dbar + 1;
          emxEnsureCapacity_real_T(b_y, i3);
          b_y_data = b_y->data;
          c_loop_ub = (int)dbar;
          for (i3 = 0; i3 <= c_loop_ub; i3++) {
            b_y_data[i3] = i3;
          }
          i3 = b_bepsbar->size[0] * b_bepsbar->size[1];
          b_bepsbar->size[0] = 1;
          b_bepsbar->size[1] = (int)dbar + 1;
          emxEnsureCapacity_real_T(b_bepsbar, i3);
          b_bepsbar_data = b_bepsbar->data;
          for (i3 = 0; i3 <= c_loop_ub; i3++) {
            b_bepsbar_data[i3] = i3;
          }
        }
        i3 = r2->size[0] * r2->size[1];
        r2->size[0] = 1;
        r2->size[1] = b_bepsbar->size[1];
        emxEnsureCapacity_real_T(r2, i3);
        B_data = r2->data;
        ct = nel_data[n];
        b_loop_ub = b_bepsbar->size[1];
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          B_data[i3] = Bbar_data[(int)(ct + b_bepsbar_data[i3]) - 1];
        }
        i3 = b_bepsbar->size[0] * b_bepsbar->size[1];
        b_bepsbar->size[0] = 1;
        b_bepsbar->size[1] = b_y->size[1];
        emxEnsureCapacity_real_T(b_bepsbar, i3);
        b_bepsbar_data = b_bepsbar->data;
        ct = nel_data[n];
        b_loop_ub = b_y->size[1];
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          b_bepsbar_data[i3] = bNbar_data[(int)(ct + b_y_data[i3]) - 1];
        }
        if (dbar < 0.0) {
          y->size[0] = 1;
          y->size[1] = 0;
        } else {
          i3 = y->size[0] * y->size[1];
          y->size[0] = 1;
          y->size[1] = (int)dbar + 1;
          emxEnsureCapacity_real_T(y, i3);
          y_data = y->data;
          b_loop_ub = (int)dbar;
          for (i3 = 0; i3 <= b_loop_ub; i3++) {
            y_data[i3] = i3;
          }
        }
        i3 = r->size[0] * r->size[1];
        r->size[0] = 1;
        r->size[1] = y->size[1];
        emxEnsureCapacity_int32_T(r, i3);
        r1 = r->data;
        ct = nel_data[n];
        b_loop_ub = y->size[1];
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          r1[i3] = (int)(ct + y_data[i3]);
        }
        if (b_bepsbar->size[1] == r2->size[1]) {
          b_loop_ub = r->size[1];
          for (i3 = 0; i3 < b_loop_ub; i3++) {
            bNbar_data[r1[i3] - 1] =
                b_bepsbar_data[i3] + B_data[i3] * st * nxi0p * wg_data[n];
          }
        } else {
          binary_expand_op_4(bNbar, r, b_bepsbar, r2, st, nxi0p, wg, n);
          bNbar_data = bNbar->data;
        }
      } else {
        double D22[9];
        double D22VV[9];
        double E22[9];
        double Q1[9];
        double VV[9];
        double F3[3];
        double b_a;
        double c_a;
        int i9;
        /*  pass == 2 */
        Nbar = 0.0;
        b_loop_ub = cNbar->size[0];
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          Nbar += cNbar_data[i3] * Bbar_data[i3];
        }
        for (i3 = 0; i3 < 3; i3++) {
          d1 = 0.0;
          d2 = Theta[i3];
          d3 = Theta[i3 + 3];
          L_chiO = Theta[i3 + 6];
          d4 = 0.0;
          for (i4 = 0; i4 < 3; i4++) {
            c_loop_ub = i3 + 3 * i4;
            ThetaR0[c_loop_ub] = (d2 * R0_[3 * i4] + d3 * R0_[3 * i4 + 1]) +
                                 L_chiO * R0_[3 * i4 + 2];
            d1 += Em22[c_loop_ub] * xi0p[i4];
            d4 += Dm22[c_loop_ub] * xip[i4];
          }
          D22WW[i3] = d1 + d4;
        }
        a = -(rho * nxi0p);
        d1 = D22WW[0];
        d2 = D22WW[1];
        d3 = D22WW[2];
        for (b_i = 0; b_i < 3; b_i++) {
          mm[b_i] = (ThetaR0[b_i] * d1 + ThetaR0[b_i + 3] * d2) +
                    ThetaR0[b_i + 6] * d3;
          xi0p[b_i] = a * u[b_i];
        }
        L_chiO = 0.0;
        b_tau0 = 0.0;
        for (i3 = 0; i3 < 3; i3++) {
          F3[i3] = (L_rrII[i3] * mm[0] + L_rrII[i3 + 3] * mm[1]) +
                   L_rrII[i3 + 6] * mm[2];
          xip[i3] = tau.contents[i3] * Nbar +
                    ((L_rrI[i3] * mm[0] + L_rrI[i3 + 3] * mm[1]) +
                     L_rrI[i3 + 6] * mm[2]);
          L_chiO += b_L_chiO[i3] * mm[i3];
          b_tau0 += tau0.contents[i3] * mm[i3];
        }
        i3 = (int)(d + 1.0);
        for (c_loop_ub = 0; c_loop_ub < i3; c_loop_ub++) {
          b_i = (int)(nel_ + (double)c_loop_ub) - 1;
          a = B_data[b_i];
          b_a = Bp_data[b_i];
          c_a = Bpp_data[b_i];
          ct = (double)c_loop_ub * 3.0;
          F__data[(int)(ct + 1.0) - 1] =
              ((a * xi0p[0] + b_a * xip[0]) + c_a * F3[0]) * wg_data[n];
          F__data[(int)(ct + 2.0) - 1] =
              ((a * xi0p[1] + b_a * xip[1]) + c_a * F3[1]) * wg_data[n];
          F__data[(int)(ct + 3.0) - 1] =
              ((a * xi0p[2] + b_a * xip[2]) + c_a * F3[2]) * wg_data[n];
        }
        d1 = nel_data[n];
        d2 = (d1 - 1.0) * 3.0 + 1.0;
        d1 = (d1 + d) * 3.0;
        if (d2 > d1) {
          i4 = 0;
          i5 = 0;
          i6 = 0;
          i7 = 0;
        } else {
          i4 = (int)d2 - 1;
          i5 = (int)d1;
          i6 = (int)d2 - 1;
          i7 = (int)d1;
        }
        if (i5 - i4 == F_->size[0]) {
          c_loop_ub = i7 - i6;
          i5 = b_bepsbar->size[0] * b_bepsbar->size[1];
          b_bepsbar->size[0] = 1;
          b_bepsbar->size[1] = c_loop_ub;
          emxEnsureCapacity_real_T(b_bepsbar, i5);
          b_bepsbar_data = b_bepsbar->data;
          for (i5 = 0; i5 < c_loop_ub; i5++) {
            b_bepsbar_data[i5] = F_data[i4 + i5] + F__data[i5];
          }
          b_loop_ub = b_bepsbar->size[1];
          for (i4 = 0; i4 < b_loop_ub; i4++) {
            F_data[i6 + i4] = b_bepsbar_data[i4];
          }
        } else {
          binary_expand_op_14(F, i6, i7, i4, i5 - 1, F_);
          F_data = F->data;
        }
        if (dbrev < 0.0) {
          b_y->size[1] = 0;
        } else {
          i4 = b_y->size[0] * b_y->size[1];
          b_y->size[0] = 1;
          b_y->size[1] = (int)dbrev + 1;
          emxEnsureCapacity_real_T(b_y, i4);
          b_y_data = b_y->data;
          b_loop_ub = (int)dbrev;
          for (i4 = 0; i4 <= b_loop_ub; i4++) {
            b_y_data[i4] = i4;
          }
        }
        i4 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        emxEnsureCapacity_real_T(b_y, i4);
        b_y_data = b_y->data;
        ct = nel_data[n];
        b_loop_ub = b_y->size[1] - 1;
        for (i4 = 0; i4 <= b_loop_ub; i4++) {
          b_y_data[i4] += ct;
        }
        if (dbrev < 0.0) {
          y->size[1] = 0;
        } else {
          i4 = y->size[0] * y->size[1];
          y->size[0] = 1;
          y->size[1] = (int)dbrev + 1;
          emxEnsureCapacity_real_T(y, i4);
          y_data = y->data;
          b_loop_ub = (int)dbrev;
          for (i4 = 0; i4 <= b_loop_ub; i4++) {
            y_data[i4] = i4;
          }
        }
        i4 = y->size[0] * y->size[1];
        y->size[0] = 1;
        emxEnsureCapacity_real_T(y, i4);
        y_data = y->data;
        ct = nel_data[n];
        b_loop_ub = y->size[1] - 1;
        for (i4 = 0; i4 <= b_loop_ub; i4++) {
          y_data[i4] += ct;
        }
        i4 = r->size[0] * r->size[1];
        r->size[0] = 1;
        r->size[1] = y->size[1];
        emxEnsureCapacity_int32_T(r, i4);
        r1 = r->data;
        b_loop_ub = y->size[1];
        for (i4 = 0; i4 < b_loop_ub; i4++) {
          r1[i4] = (int)y_data[i4];
        }
        i4 = b_bepsbar->size[0] * b_bepsbar->size[1];
        b_bepsbar->size[0] = 1;
        b_bepsbar->size[1] = r->size[1];
        emxEnsureCapacity_real_T(b_bepsbar, i4);
        b_bepsbar_data = b_bepsbar->data;
        b_loop_ub = r->size[1];
        for (i4 = 0; i4 < b_loop_ub; i4++) {
          d3 = b_y_data[i4];
          b_bepsbar_data[i4] =
              mu_data[(int)d3 - 1] + (Bbrev_data[(int)d3 - 1] * L_chiO +
                                      Bbrevp_data[(int)d3 - 1] * b_tau0) *
                                         wg_data[n];
        }
        b_loop_ub = b_bepsbar->size[1];
        for (i4 = 0; i4 < b_loop_ub; i4++) {
          mu_data[r1[i4] - 1] = b_bepsbar_data[i4];
        }
        for (i4 = 0; i4 < 9; i4++) {
          htau0[i4] *= b_varThetadot;
        }
        XX_contents[0] = (b_L_chiO_tmp * tau0.contents[2] -
                          tau0.contents[1] * c_L_chiO_tmp) -
                         D22WW_tmp;
        XX_contents[1] =
            (tau0.contents[0] * c_L_chiO_tmp - L_chiO_tmp * tau0.contents[2]) -
            b_D22WW_tmp;
        XX_contents[2] =
            (L_chiO_tmp * tau0.contents[1] - tau0.contents[0] * b_L_chiO_tmp) -
            c_D22WW_tmp;
        for (i4 = 0; i4 < 3; i4++) {
          d3 = 0.0;
          for (i5 = 0; i5 < 3; i5++) {
            i6 = i4 + 3 * i5;
            d3 += PP.contents[i6] * xipd[i5];
            R0_[i5 + 3 * i4] = ThetaR0[i6];
            Q2[i6] = (ThetaR0[i4] * Em22[3 * i5] +
                      ThetaR0[i4 + 3] * Em22[3 * i5 + 1]) +
                     ThetaR0[i4 + 6] * Em22[3 * i5 + 2];
          }
          xip[i4] = d3 / nxip.contents;
        }
        for (i4 = 0; i4 < 3; i4++) {
          d3 = Q2[i4];
          L_chiO = Q2[i4 + 3];
          d4 = Q2[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            E22[i4 + 3 * i5] = ((d3 * R0_[3 * i5] + L_chiO * R0_[3 * i5 + 1]) +
                                d4 * R0_[3 * i5 + 2]) /
                               nxi0p;
          }
          d3 = ThetaR0[i4];
          L_chiO = ThetaR0[i4 + 3];
          d4 = ThetaR0[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            Q2[i4 + 3 * i5] = (d3 * Dm22[3 * i5] + L_chiO * Dm22[3 * i5 + 1]) +
                              d4 * Dm22[3 * i5 + 2];
          }
          d3 = Q2[i4];
          L_chiO = Q2[i4 + 3];
          d4 = Q2[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            D22[i4 + 3 * i5] = ((d3 * R0_[3 * i5] + L_chiO * R0_[3 * i5 + 1]) +
                                d4 * R0_[3 * i5 + 2]) /
                               nxi0p;
          }
        }
        Qtil_rr_I(&nxip, &K, &L_KI, &tau, &I3, &PP, &tau0, &L_tauI, &K0, &G,
                  &L_GI, xipd, ThetaR0);
        Qtil_rr_IV(&nxip, &L_KII, &tau, &tau0, &L_tauI, &L_GII, xippd, R0_);
        for (i4 = 0; i4 < 3; i4++) {
          d3 = htau0[i4];
          L_chiO = htau0[i4 + 3];
          d4 = htau0[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            c_loop_ub = i4 + 3 * i5;
            Theta[c_loop_ub] =
                (ThetaR0[c_loop_ub] + R0_[c_loop_ub]) -
                ((d3 * y_tmp[3 * i5] + L_chiO * y_tmp[3 * i5 + 1]) +
                 d4 * y_tmp[3 * i5 + 2]);
          }
        }
        Qtil_rr_III(&nxip, &tau, &L_KII, &L_tauI, &tau0, &L_GII, xipd, ThetaR0);
        F3[0] = tau0.contents[1] * mm[2] - mm[1] * tau0.contents[2];
        F3[1] = mm[0] * tau0.contents[2] - tau0.contents[0] * mm[2];
        F3[2] = tau0.contents[0] * mm[1] - mm[0] * tau0.contents[1];
        for (i4 = 0; i4 < 3; i4++) {
          d3 = htau0[i4];
          L_chiO = htau0[i4 + 3];
          d4 = htau0[i4 + 6];
          d5 = 0.0;
          for (i5 = 0; i5 < 3; i5++) {
            c_loop_ub = i4 + 3 * i5;
            VV[c_loop_ub] =
                ThetaR0[c_loop_ub] -
                ((d3 * b_y_tmp[3 * i5] + L_chiO * b_y_tmp[3 * i5 + 1]) +
                 d4 * b_y_tmp[3 * i5 + 2]);
            d5 += E22[c_loop_ub] * b_L_chiO[i5];
          }
          F3[i4] += d5;
        }
        /* %%% the b terms %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        for (b_loop_ub = 0; b_loop_ub < i; b_loop_ub++) {
          for (ddc = 0; ddc < i3; ddc++) {
            b_varTheta = (double)ddc * 3.0;
            ct = Bbar_data[(int)(nel_ + (double)b_loop_ub) - 1] *
                 Bp_data[(int)(nel_ + (double)ddc) - 1] * wg_data[n];
            bepscomma__data[b_loop_ub + bepscomma_->size[0] *
                                            ((int)(b_varTheta + 1.0) - 1)] =
                ct * tau.contents[0];
            bepsdcomma__data[b_loop_ub + bepsdcomma_->size[0] *
                                             ((int)(b_varTheta + 1.0) - 1)] =
                ct * xip[0];
            bepscomma__data[b_loop_ub + bepscomma_->size[0] *
                                            ((int)(b_varTheta + 2.0) - 1)] =
                ct * tau.contents[1];
            bepsdcomma__data[b_loop_ub + bepsdcomma_->size[0] *
                                             ((int)(b_varTheta + 2.0) - 1)] =
                ct * xip[1];
            bepscomma__data[b_loop_ub + bepscomma_->size[0] *
                                            ((int)(b_varTheta + 3.0) - 1)] =
                ct * b_d;
            bepsdcomma__data[b_loop_ub + bepsdcomma_->size[0] *
                                             ((int)(b_varTheta + 3.0) - 1)] =
                ct * xip[2];
          }
        }
        if (d2 > d1) {
          i4 = 0;
          i5 = 0;
        } else {
          i4 = (int)d2 - 1;
          i5 = (int)d1;
        }
        if (dbar < 0.0) {
          b_y->size[1] = 0;
        } else {
          i6 = b_y->size[0] * b_y->size[1];
          b_y->size[0] = 1;
          b_y->size[1] = (int)dbar + 1;
          emxEnsureCapacity_real_T(b_y, i6);
          b_y_data = b_y->data;
          b_loop_ub = (int)dbar;
          for (i6 = 0; i6 <= b_loop_ub; i6++) {
            b_y_data[i6] = i6;
          }
        }
        i6 = b_y->size[0] * b_y->size[1];
        b_y->size[0] = 1;
        emxEnsureCapacity_real_T(b_y, i6);
        b_y_data = b_y->data;
        ct = nel_data[n];
        b_loop_ub = b_y->size[1] - 1;
        for (i6 = 0; i6 <= b_loop_ub; i6++) {
          b_y_data[i6] += ct;
        }
        i6 = Bbar->size[0];
        Bbar->size[0] = b_y->size[1];
        emxEnsureCapacity_real_T(Bbar, i6);
        Bbar_data = Bbar->data;
        b_loop_ub = b_y->size[1];
        for (i6 = 0; i6 < b_loop_ub; i6++) {
          Bbar_data[i6] = b_y_data[i6];
        }
        if (d2 > d1) {
          i6 = 0;
          i7 = 0;
        } else {
          i6 = (int)d2 - 1;
          i7 = (int)d1;
        }
        i8 = r3->size[0];
        r3->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r3, i8);
        r1 = r3->data;
        b_loop_ub = Bbar->size[0];
        for (i8 = 0; i8 < b_loop_ub; i8++) {
          r1[i8] = (int)Bbar_data[i8] - 1;
        }
        b_loop_ub = i5 - i4;
        if ((Bbar->size[0] == bepscomma_->size[0]) &&
            (b_loop_ub == bepscomma_->size[1])) {
          i5 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = Bbar->size[0];
          b_bepscomma->size[1] = b_loop_ub;
          emxEnsureCapacity_real_T(b_bepscomma, i5);
          B_data = b_bepscomma->data;
          for (i5 = 0; i5 < b_loop_ub; i5++) {
            c_loop_ub = Bbar->size[0];
            for (i8 = 0; i8 < c_loop_ub; i8++) {
              B_data[i8 + b_bepscomma->size[0] * i5] =
                  bepscomma_data[((int)Bbar_data[i8] +
                                  bepscomma->size[0] * (i4 + i5)) -
                                 1] +
                  bepscomma__data[i8 + bepscomma_->size[0] * i5];
            }
          }
          b_i = r3->size[0];
          c_loop_ub = i7 - i6;
          for (i4 = 0; i4 < c_loop_ub; i4++) {
            for (i5 = 0; i5 < b_i; i5++) {
              bepscomma_data[r1[i5] + bepscomma->size[0] * (i6 + i4)] =
                  B_data[i5 + b_i * i4];
            }
          }
        } else {
          binary_expand_op_7(bepscomma, r3, i6, i7, Bbar, i4, i5 - 1,
                             bepscomma_);
          bepscomma_data = bepscomma->data;
        }
        if (d2 > d1) {
          i4 = 0;
          i5 = 0;
          i6 = 0;
          i7 = 0;
        } else {
          i4 = (int)d2 - 1;
          i5 = (int)d1;
          i6 = (int)d2 - 1;
          i7 = (int)d1;
        }
        i8 = r3->size[0];
        r3->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r3, i8);
        r1 = r3->data;
        b_loop_ub = Bbar->size[0];
        for (i8 = 0; i8 < b_loop_ub; i8++) {
          r1[i8] = (int)Bbar_data[i8] - 1;
        }
        b_loop_ub = i5 - i4;
        if ((Bbar->size[0] == bepsdcomma_->size[0]) &&
            (b_loop_ub == bepsdcomma_->size[1])) {
          i5 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = Bbar->size[0];
          b_bepscomma->size[1] = b_loop_ub;
          emxEnsureCapacity_real_T(b_bepscomma, i5);
          B_data = b_bepscomma->data;
          for (i5 = 0; i5 < b_loop_ub; i5++) {
            c_loop_ub = Bbar->size[0];
            for (i8 = 0; i8 < c_loop_ub; i8++) {
              B_data[i8 + b_bepscomma->size[0] * i5] =
                  bepsdcomma_data[((int)Bbar_data[i8] +
                                   bepsdcomma->size[0] * (i4 + i5)) -
                                  1] +
                  bepsdcomma__data[i8 + bepsdcomma_->size[0] * i5];
            }
          }
          b_i = r3->size[0];
          c_loop_ub = i7 - i6;
          for (i4 = 0; i4 < c_loop_ub; i4++) {
            for (i5 = 0; i5 < b_i; i5++) {
              bepsdcomma_data[r1[i5] + bepsdcomma->size[0] * (i6 + i4)] =
                  B_data[i5 + b_i * i4];
            }
          }
        } else {
          binary_expand_op_7(bepsdcomma, r3, i6, i7, Bbar, i4, i5 - 1,
                             bepsdcomma_);
          bepsdcomma_data = bepsdcomma->data;
        }
        /* %%% F_ij %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        a = Nbar / nxip.contents;
        memcpy(&R0_[0], &PP.contents[0], 9U * sizeof(double));
        Q_rr_I(&L_KI, &K, &nxip, &I3, &tau, &PP, &tau0, &L_tauI, &K0, &L_GI, &G,
               mm, ThetaR0);
        for (i4 = 0; i4 < 3; i4++) {
          b_d = L_rrI[i4];
          d3 = L_rrI[i4 + 3];
          L_chiO = L_rrI[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            Q2[i4 + 3 * i5] = (b_d * E22[3 * i5] + d3 * E22[3 * i5 + 1]) +
                              L_chiO * E22[3 * i5 + 2];
          }
          b_d = Q2[i4];
          d3 = Q2[i4 + 3];
          L_chiO = Q2[i4 + 6];
          d4 = L_rrII[i4];
          d5 = L_rrII[i4 + 3];
          d6 = L_rrII[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            i6 = 3 * i5 + 1;
            i7 = 3 * i5 + 2;
            c_loop_ub = i4 + 3 * i5;
            Q1[c_loop_ub] =
                (a * R0_[c_loop_ub] + ThetaR0[c_loop_ub]) +
                ((b_d * y_tmp[3 * i5] + d3 * y_tmp[i6]) + L_chiO * y_tmp[i7]);
            R0_[c_loop_ub] = (d4 * E22[3 * i5] + d5 * E22[i6]) + d6 * E22[i7];
          }
        }
        Q_rr_III(&nxip, &L_KII, &tau, &tau0, &L_tauI, &L_GII, mm, ThetaR0);
        for (i4 = 0; i4 < 3; i4++) {
          b_d = R0_[i4];
          d3 = R0_[i4 + 3];
          L_chiO = R0_[i4 + 6];
          d4 = D22[i4];
          d5 = D22[i4 + 3];
          d6 = D22[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            i6 = 3 * i5 + 1;
            i7 = 3 * i5 + 2;
            c_loop_ub = i4 + 3 * i5;
            D22VV[c_loop_ub] = (d4 * VV[3 * i5] + d5 * VV[i6]) + d6 * VV[i7];
            D22UU[c_loop_ub] =
                (d4 * Theta[3 * i5] + d5 * Theta[i6]) + d6 * Theta[i7];
            htau0[c_loop_ub] = (b_d * b_y_tmp[3 * i5] + d3 * b_y_tmp[i6]) +
                               L_chiO * b_y_tmp[i7];
            Q2[c_loop_ub] =
                ThetaR0[c_loop_ub] +
                ((b_d * y_tmp[3 * i5] + d3 * y_tmp[i6]) + L_chiO * y_tmp[i7]);
          }
        }
        for (i4 = 0; i4 < 3; i4++) {
          b_d = L_rrI[i4];
          d3 = L_rrI[i4 + 3];
          L_chiO = L_rrI[i4 + 6];
          d4 = L_rrII[i4];
          d5 = L_rrII[i4 + 3];
          d6 = L_rrII[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            tau0_tmp = D22UU[3 * i5];
            a = b_d * tau0_tmp;
            c_a = D22VV[3 * i5];
            c_tmp = b_d * c_a;
            ct = d4 * tau0_tmp;
            b_a = d4 * c_a;
            i6 = 3 * i5 + 1;
            tau0_tmp = D22UU[i6];
            a += d3 * tau0_tmp;
            c_a = D22VV[i6];
            c_tmp += d3 * c_a;
            ct += d5 * tau0_tmp;
            b_a += d5 * c_a;
            i6 = 3 * i5 + 2;
            tau0_tmp = D22UU[i6];
            a += L_chiO * tau0_tmp;
            c_a = D22VV[i6];
            c_tmp += L_chiO * c_a;
            ct += d6 * tau0_tmp;
            b_a += d6 * c_a;
            c_loop_ub = i4 + 3 * i5;
            VV[c_loop_ub] = b_a;
            ThetaR0[c_loop_ub] = ct;
            Theta[c_loop_ub] = c_tmp;
            R0_[c_loop_ub] = a;
          }
        }
        for (ddc = 0; ddc < i3; ddc++) {
          for (b_loop_ub = 0; b_loop_ub < i3; b_loop_ub++) {
            c_loop_ub = (int)(nel_ + (double)b_loop_ub) - 1;
            ct = Bp_data[c_loop_ub];
            b_i = (int)(nel_ + (double)ddc) - 1;
            tau0_tmp = Bp_data[b_i];
            a = ct * tau0_tmp;
            c_tmp = Bpp_data[b_i];
            b_a = ct * c_tmp;
            ct = Bpp_data[c_loop_ub];
            c_a = ct * tau0_tmp;
            ct *= c_tmp;
            tau0_tmp = (double)b_loop_ub * 3.0;
            b_varTheta = (double)ddc * 3.0;
            for (i4 = 0; i4 < 3; i4++) {
              c_loop_ub = (int)(b_varTheta + ((double)i4 + 1.0)) - 1;
              F_ij__data[((int)(tau0_tmp + 1.0) + F_ij_->size[0] * c_loop_ub) -
                         1] = (((a * (Q1[3 * i4] + R0_[3 * i4]) +
                                 b_a * (Q2[i4] + Theta[3 * i4])) +
                                c_a * (Q2[3 * i4] + ThetaR0[3 * i4])) +
                               ct * (htau0[3 * i4] + VV[3 * i4])) *
                              wg_data[n];
              b_i = 3 * i4 + 1;
              F_ij__data[((int)(tau0_tmp + 2.0) + F_ij_->size[0] * c_loop_ub) -
                         1] = (((a * (Q1[b_i] + R0_[b_i]) +
                                 b_a * (Q2[i4 + 3] + Theta[b_i])) +
                                c_a * (Q2[b_i] + ThetaR0[b_i])) +
                               ct * (htau0[b_i] + VV[b_i])) *
                              wg_data[n];
              b_i = 3 * i4 + 2;
              F_ij__data[((int)(tau0_tmp + 3.0) + F_ij_->size[0] * c_loop_ub) -
                         1] = (((a * (Q1[b_i] + R0_[b_i]) +
                                 b_a * (Q2[i4 + 6] + Theta[b_i])) +
                                c_a * (Q2[b_i] + ThetaR0[b_i])) +
                               ct * (htau0[b_i] + VV[b_i])) *
                              wg_data[n];
            }
          }
        }
        if (d2 > d1) {
          i4 = 0;
          i5 = 0;
          i6 = 0;
          i7 = 0;
          i8 = 0;
          b_i = 0;
          i9 = 0;
          ddc = 0;
        } else {
          i4 = (int)d2 - 1;
          i5 = (int)d1;
          i6 = (int)d2 - 1;
          i7 = (int)d1;
          i8 = (int)d2 - 1;
          b_i = (int)d1;
          i9 = (int)d2 - 1;
          ddc = (int)d1;
        }
        b_loop_ub = i5 - i4;
        c_loop_ub = i7 - i6;
        if ((b_loop_ub == F_ij_->size[0]) && (c_loop_ub == F_ij_->size[1])) {
          i5 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = b_loop_ub;
          b_bepscomma->size[1] = c_loop_ub;
          emxEnsureCapacity_real_T(b_bepscomma, i5);
          B_data = b_bepscomma->data;
          for (i5 = 0; i5 < c_loop_ub; i5++) {
            for (i7 = 0; i7 < b_loop_ub; i7++) {
              B_data[i7 + b_bepscomma->size[0] * i5] =
                  F_ij_data[(i4 + i7) + F_ij->size[0] * (i6 + i5)] +
                  F_ij__data[i7 + F_ij_->size[0] * i5];
            }
          }
          b_i -= i8;
          c_loop_ub = ddc - i9;
          for (i4 = 0; i4 < c_loop_ub; i4++) {
            for (i5 = 0; i5 < b_i; i5++) {
              F_ij_data[(i8 + i5) + F_ij->size[0] * (i9 + i4)] =
                  B_data[i5 + b_i * i4];
            }
          }
        } else {
          binary_expand_op_9(F_ij, i8, b_i, i9, ddc, i4, i5 - 1, i6, i7 - 1,
                             F_ij_);
          F_ij_data = F_ij->data;
        }
        /* %%% F_ib %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        Q20[0] = rrd[1] * tau0.contents[2] - tau0.contents[1] * rrd[2];
        Q20[1] = tau0.contents[0] * rrd[2] - rrd[0] * tau0.contents[2];
        Q20[2] = rrd[0] * tau0.contents[1] - tau0.contents[0] * rrd[1];
        b_d = tau0.contents[0];
        d3 = tau0.contents[1];
        L_chiO = tau0.contents[2];
        for (i4 = 0; i4 < 3; i4++) {
          rrd[i4] = (E22[i4] * b_d + E22[i4 + 3] * d3) + E22[i4 + 6] * L_chiO;
          Q20[i4] += b_varThetadot * XX_contents[i4];
        }
        b_d = F3[0];
        d3 = F3[1];
        L_chiO = F3[2];
        d4 = rrd[0];
        d5 = rrd[1];
        d6 = rrd[2];
        tau0_tmp = Q20[0];
        a = Q20[1];
        c_a = Q20[2];
        for (i4 = 0; i4 < 3; i4++) {
          c_tmp = L_rrI[i4];
          ct = c_tmp * b_d;
          b_a = c_tmp * d4;
          c_tmp = L_rrII[i4];
          b_varTheta = c_tmp * b_d;
          st = c_tmp * d4;
          c_tmp = L_rrI[i4 + 3];
          ct += c_tmp * d3;
          b_a += c_tmp * d5;
          c_tmp = L_rrII[i4 + 3];
          b_varTheta += c_tmp * d3;
          st += c_tmp * d5;
          c_tmp = L_rrI[i4 + 6];
          ct += c_tmp * L_chiO;
          b_a += c_tmp * d6;
          c_tmp = L_rrII[i4 + 6];
          b_varTheta += c_tmp * L_chiO;
          st += c_tmp * d6;
          xipd[i4] = st;
          xi0pp[i4] = b_varTheta;
          rr[i4] = b_a;
          xippd[i4] = ct;
          D22WW[i4] =
              (D22[i4] * tau0_tmp + D22[i4 + 3] * a) + D22[i4 + 6] * c_a;
        }
        b_d = D22WW[0];
        d3 = D22WW[1];
        L_chiO = D22WW[2];
        for (i4 = 0; i4 < 3; i4++) {
          xip[i4] = (L_rrII[i4] * b_d + L_rrII[i4 + 3] * d3) +
                    L_rrII[i4 + 6] * L_chiO;
          xi0p[i4] =
              (L_rrI[i4] * b_d + L_rrI[i4 + 3] * d3) + L_rrI[i4 + 6] * L_chiO;
        }
        for (ddc = 0; ddc < i2; ddc++) {
          for (b_loop_ub = 0; b_loop_ub < i3; b_loop_ub++) {
            c_loop_ub = (int)(nel_ + (double)b_loop_ub) - 1;
            ct = Bp_data[c_loop_ub];
            b_i = (int)(nel_ + (double)ddc) - 1;
            tau0_tmp = Bbrev_data[b_i];
            a = ct * tau0_tmp;
            c_tmp = Bbrevp_data[b_i];
            b_a = ct * c_tmp;
            ct = Bpp_data[c_loop_ub];
            c_a = ct * tau0_tmp;
            ct *= c_tmp;
            tau0_tmp = (double)b_loop_ub * 3.0;
            F_ib__data[((int)(tau0_tmp + 1.0) + F_ib_->size[0] * ddc) - 1] =
                (((a * (xippd[0] + xi0p[0]) + b_a * rr[0]) +
                  c_a * (xi0pp[0] + xip[0])) +
                 ct * xipd[0]) *
                wg_data[n];
            F_ib__data[((int)(tau0_tmp + 2.0) + F_ib_->size[0] * ddc) - 1] =
                (((a * (xippd[1] + xi0p[1]) + b_a * rr[1]) +
                  c_a * (xi0pp[1] + xip[1])) +
                 ct * xipd[1]) *
                wg_data[n];
            F_ib__data[((int)(tau0_tmp + 3.0) + F_ib_->size[0] * ddc) - 1] =
                (((a * (xippd[2] + xi0p[2]) + b_a * rr[2]) +
                  c_a * (xi0pp[2] + xip[2])) +
                 ct * xipd[2]) *
                wg_data[n];
          }
        }
        if (d2 > d1) {
          i4 = 0;
          i5 = 0;
        } else {
          i4 = (int)d2 - 1;
          i5 = (int)d1;
        }
        i6 = Bbar->size[0];
        Bbar->size[0] = y->size[1];
        emxEnsureCapacity_real_T(Bbar, i6);
        Bbar_data = Bbar->data;
        b_loop_ub = y->size[1];
        for (i6 = 0; i6 < b_loop_ub; i6++) {
          Bbar_data[i6] = y_data[i6];
        }
        if (d2 > d1) {
          i6 = 0;
          i7 = 0;
        } else {
          i6 = (int)d2 - 1;
          i7 = (int)d1;
        }
        i8 = r3->size[0];
        r3->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r3, i8);
        r1 = r3->data;
        b_loop_ub = Bbar->size[0];
        for (i8 = 0; i8 < b_loop_ub; i8++) {
          r1[i8] = (int)Bbar_data[i8] - 1;
        }
        b_loop_ub = i5 - i4;
        if ((b_loop_ub == F_ib_->size[0]) &&
            (Bbar->size[0] == F_ib_->size[1])) {
          i5 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = b_loop_ub;
          b_bepscomma->size[1] = Bbar->size[0];
          emxEnsureCapacity_real_T(b_bepscomma, i5);
          B_data = b_bepscomma->data;
          c_loop_ub = Bbar->size[0];
          for (i5 = 0; i5 < c_loop_ub; i5++) {
            for (i8 = 0; i8 < b_loop_ub; i8++) {
              B_data[i8 + b_bepscomma->size[0] * i5] =
                  F_ib_data[(i4 + i8) +
                            F_ib->size[0] * ((int)Bbar_data[i5] - 1)] +
                  F_ib__data[i8 + F_ib_->size[0] * i5];
            }
          }
          b_i = i7 - i6;
          b_loop_ub = r3->size[0];
          for (i4 = 0; i4 < b_loop_ub; i4++) {
            for (i5 = 0; i5 < b_i; i5++) {
              F_ib_data[(i6 + i5) + F_ib->size[0] * r1[i4]] =
                  B_data[i5 + b_i * i4];
            }
          }
        } else {
          binary_expand_op_8(F_ib, i6, i7, r3, i4, i5 - 1, Bbar, F_ib_);
          F_ib_data = F_ib->data;
        }
        /* %%% F_ijd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        for (i4 = 0; i4 < 3; i4++) {
          b_d = L_rrI[i4];
          d3 = L_rrI[i4 + 3];
          L_chiO = L_rrI[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            VV[i4 + 3 * i5] = (b_d * D22[3 * i5] + d3 * D22[3 * i5 + 1]) +
                              L_chiO * D22[3 * i5 + 2];
          }
          b_d = VV[i4];
          d3 = VV[i4 + 3];
          L_chiO = VV[i4 + 6];
          d4 = L_rrII[i4];
          d5 = L_rrII[i4 + 3];
          d6 = L_rrII[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            i6 = 3 * i5 + 1;
            i7 = 3 * i5 + 2;
            c_loop_ub = i4 + 3 * i5;
            ThetaR0[c_loop_ub] =
                (d4 * D22[3 * i5] + d5 * D22[i6]) + d6 * D22[i7];
            Theta[c_loop_ub] = (b_d * b_y_tmp[3 * i5] + d3 * b_y_tmp[i6]) +
                               L_chiO * b_y_tmp[i7];
            htau0[c_loop_ub] =
                (b_d * y_tmp[3 * i5] + d3 * y_tmp[i6]) + L_chiO * y_tmp[i7];
          }
          b_d = ThetaR0[i4];
          d3 = ThetaR0[i4 + 3];
          L_chiO = ThetaR0[i4 + 6];
          for (i5 = 0; i5 < 3; i5++) {
            R0_[i4 + 3 * i5] =
                (b_d * b_y_tmp[3 * i5] + d3 * b_y_tmp[3 * i5 + 1]) +
                L_chiO * b_y_tmp[3 * i5 + 2];
          }
        }
        for (ddc = 0; ddc < i3; ddc++) {
          for (b_loop_ub = 0; b_loop_ub < i3; b_loop_ub++) {
            c_loop_ub = (int)(nel_ + (double)b_loop_ub) - 1;
            ct = Bp_data[c_loop_ub];
            b_i = (int)(nel_ + (double)ddc) - 1;
            tau0_tmp = Bp_data[b_i];
            a = ct * tau0_tmp;
            c_tmp = Bpp_data[b_i];
            b_a = ct * c_tmp;
            ct = Bpp_data[c_loop_ub];
            st = ct * tau0_tmp;
            c_a = ct * c_tmp;
            tau0_tmp = (double)b_loop_ub * 3.0;
            b_varTheta = (double)ddc * 3.0;
            for (i4 = 0; i4 < 3; i4++) {
              c_loop_ub = (int)(b_varTheta + ((double)i4 + 1.0)) - 1;
              F_ijd__data[((int)(tau0_tmp + 1.0) +
                           F_ijd_->size[0] * c_loop_ub) -
                          1] = (((a * htau0[3 * i4] + b_a * Theta[3 * i4]) +
                                 st * Theta[i4]) +
                                c_a * R0_[3 * i4]) *
                               wg_data[n];
              b_i = 3 * i4 + 1;
              F_ijd__data[((int)(tau0_tmp + 2.0) +
                           F_ijd_->size[0] * c_loop_ub) -
                          1] =
                  (((a * htau0[b_i] + b_a * Theta[b_i]) + st * Theta[i4 + 3]) +
                   c_a * R0_[b_i]) *
                  wg_data[n];
              b_i = 3 * i4 + 2;
              F_ijd__data[((int)(tau0_tmp + 3.0) +
                           F_ijd_->size[0] * c_loop_ub) -
                          1] =
                  (((a * htau0[b_i] + b_a * Theta[b_i]) + st * Theta[i4 + 6]) +
                   c_a * R0_[b_i]) *
                  wg_data[n];
            }
          }
        }
        if (d2 > d1) {
          i4 = 0;
          i5 = 0;
          i6 = 0;
          i7 = 0;
          i8 = 0;
          b_i = 0;
          i9 = 0;
          ddc = 0;
        } else {
          i4 = (int)d2 - 1;
          i5 = (int)d1;
          i6 = (int)d2 - 1;
          i7 = (int)d1;
          i8 = (int)d2 - 1;
          b_i = (int)d1;
          i9 = (int)d2 - 1;
          ddc = (int)d1;
        }
        b_loop_ub = i5 - i4;
        c_loop_ub = i7 - i6;
        if ((b_loop_ub == F_ijd_->size[0]) && (c_loop_ub == F_ijd_->size[1])) {
          i5 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = b_loop_ub;
          b_bepscomma->size[1] = c_loop_ub;
          emxEnsureCapacity_real_T(b_bepscomma, i5);
          B_data = b_bepscomma->data;
          for (i5 = 0; i5 < c_loop_ub; i5++) {
            for (i7 = 0; i7 < b_loop_ub; i7++) {
              B_data[i7 + b_bepscomma->size[0] * i5] =
                  F_ijd_data[(i4 + i7) + F_ijd->size[0] * (i6 + i5)] +
                  F_ijd__data[i7 + F_ijd_->size[0] * i5];
            }
          }
          b_i -= i8;
          c_loop_ub = ddc - i9;
          for (i4 = 0; i4 < c_loop_ub; i4++) {
            for (i5 = 0; i5 < b_i; i5++) {
              F_ijd_data[(i8 + i5) + F_ijd->size[0] * (i9 + i4)] =
                  B_data[i5 + b_i * i4];
            }
          }
        } else {
          binary_expand_op_9(F_ijd, i8, b_i, i9, ddc, i4, i5 - 1, i6, i7 - 1,
                             F_ijd_);
          F_ijd_data = F_ijd->data;
        }
        /* %%% F_ibd and mu_ajd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        b_d = b_L_chiO[0];
        d3 = b_L_chiO[1];
        L_chiO = b_L_chiO[2];
        d4 = tau0.contents[0];
        d5 = tau0.contents[1];
        d6 = tau0.contents[2];
        for (i4 = 0; i4 < 3; i4++) {
          tau0_tmp = VV[i4];
          a = tau0_tmp * b_d;
          c_a = tau0_tmp * d4;
          tau0_tmp = ThetaR0[i4];
          c_tmp = tau0_tmp * b_d;
          ct = tau0_tmp * d4;
          tau0_tmp = VV[i4 + 3];
          a += tau0_tmp * d3;
          c_a += tau0_tmp * d5;
          tau0_tmp = ThetaR0[i4 + 3];
          c_tmp += tau0_tmp * d3;
          ct += tau0_tmp * d5;
          tau0_tmp = VV[i4 + 6];
          a += tau0_tmp * L_chiO;
          c_a += tau0_tmp * d6;
          tau0_tmp = ThetaR0[i4 + 6];
          c_tmp += tau0_tmp * L_chiO;
          ct += tau0_tmp * d6;
          Q20[i4] = ct;
          F3[i4] = c_tmp;
          xip[i4] = c_a;
          xi0p[i4] = a;
        }
        for (ddc = 0; ddc < i2; ddc++) {
          for (b_loop_ub = 0; b_loop_ub < i3; b_loop_ub++) {
            c_loop_ub = (int)(nel_ + (double)b_loop_ub) - 1;
            ct = Bp_data[c_loop_ub];
            b_i = (int)(nel_ + (double)ddc) - 1;
            tau0_tmp = Bbrev_data[b_i];
            a = ct * tau0_tmp;
            c_tmp = Bbrevp_data[b_i];
            b_a = ct * c_tmp;
            ct = Bpp_data[c_loop_ub];
            c_a = ct * tau0_tmp;
            ct *= c_tmp;
            tau0_tmp = (double)b_loop_ub * 3.0;
            F_ibd__data[((int)(tau0_tmp + 1.0) + F_ibd_->size[0] * ddc) - 1] =
                (((a * xi0p[0] + b_a * xip[0]) + c_a * F3[0]) + ct * Q20[0]) *
                wg_data[n];
            F_ibd__data[((int)(tau0_tmp + 2.0) + F_ibd_->size[0] * ddc) - 1] =
                (((a * xi0p[1] + b_a * xip[1]) + c_a * F3[1]) + ct * Q20[1]) *
                wg_data[n];
            F_ibd__data[((int)(tau0_tmp + 3.0) + F_ibd_->size[0] * ddc) - 1] =
                (((a * xi0p[2] + b_a * xip[2]) + c_a * F3[2]) + ct * Q20[2]) *
                wg_data[n];
          }
        }
        if (d2 > d1) {
          i4 = 0;
          i5 = 0;
          i6 = 0;
          i7 = 0;
        } else {
          i4 = (int)d2 - 1;
          i5 = (int)d1;
          i6 = (int)d2 - 1;
          i7 = (int)d1;
        }
        i8 = r3->size[0];
        r3->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r3, i8);
        r1 = r3->data;
        b_loop_ub = Bbar->size[0];
        for (i8 = 0; i8 < b_loop_ub; i8++) {
          r1[i8] = (int)Bbar_data[i8] - 1;
        }
        b_loop_ub = i5 - i4;
        if ((b_loop_ub == F_ibd_->size[0]) &&
            (Bbar->size[0] == F_ibd_->size[1])) {
          i5 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = b_loop_ub;
          b_bepscomma->size[1] = Bbar->size[0];
          emxEnsureCapacity_real_T(b_bepscomma, i5);
          B_data = b_bepscomma->data;
          c_loop_ub = Bbar->size[0];
          for (i5 = 0; i5 < c_loop_ub; i5++) {
            for (i8 = 0; i8 < b_loop_ub; i8++) {
              B_data[i8 + b_bepscomma->size[0] * i5] =
                  F_ibd_data[(i4 + i8) +
                             F_ibd->size[0] * ((int)Bbar_data[i5] - 1)] +
                  F_ibd__data[i8 + F_ibd_->size[0] * i5];
            }
          }
          b_i = i7 - i6;
          b_loop_ub = r3->size[0];
          for (i4 = 0; i4 < b_loop_ub; i4++) {
            for (i5 = 0; i5 < b_i; i5++) {
              F_ibd_data[(i6 + i5) + F_ibd->size[0] * r1[i4]] =
                  B_data[i5 + b_i * i4];
            }
          }
        } else {
          binary_expand_op_8(F_ibd, i6, i7, r3, i4, i5 - 1, Bbar, F_ibd_);
          F_ibd_data = F_ibd->data;
        }
        /* %%% mu_aj %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        b_d = b_L_chiO[0];
        d3 = b_L_chiO[1];
        L_chiO = b_L_chiO[2];
        d4 = tau0.contents[0];
        d5 = tau0.contents[1];
        d6 = tau0.contents[2];
        for (i4 = 0; i4 < 3; i4++) {
          tau0_tmp = D22UU[3 * i4];
          i5 = 3 * i4 + 1;
          a = D22UU[i5];
          i6 = 3 * i4 + 2;
          c_a = D22UU[i6];
          xi0p[i4] = (b_d * tau0_tmp + d3 * a) + L_chiO * c_a;
          xip[i4] = (d4 * tau0_tmp + d5 * a) + d6 * c_a;
          tau0_tmp = D22VV[3 * i4];
          a = D22VV[i5];
          c_a = D22VV[i6];
          F3[i4] = (b_d * tau0_tmp + d3 * a) + L_chiO * c_a;
          Q20[i4] = (d4 * tau0_tmp + d5 * a) + d6 * c_a;
        }
        for (ddc = 0; ddc < i3; ddc++) {
          for (b_loop_ub = 0; b_loop_ub < i2; b_loop_ub++) {
            b_i = (int)(nel_ + (double)b_loop_ub) - 1;
            c_loop_ub = (int)(nel_ + (double)ddc) - 1;
            ct = Bp_data[c_loop_ub];
            tau0_tmp = Bbrev_data[b_i];
            a = tau0_tmp * ct;
            c_tmp = Bbrevp_data[b_i];
            b_a = c_tmp * ct;
            ct = Bpp_data[c_loop_ub];
            c_a = tau0_tmp * ct;
            ct *= c_tmp;
            b_varTheta = (double)ddc * 3.0;
            mu_aj__data[b_loop_ub +
                        mu_aj_->size[0] * ((int)(b_varTheta + 1.0) - 1)] =
                (((a * (xippd[0] + xi0p[0]) + b_a * (rr[0] + xip[0])) +
                  c_a * (xi0pp[0] + F3[0])) +
                 ct * (xipd[0] + Q20[0])) *
                wg_data[n];
            mu_aj__data[b_loop_ub +
                        mu_aj_->size[0] * ((int)(b_varTheta + 2.0) - 1)] =
                (((a * (xippd[1] + xi0p[1]) + b_a * (rr[1] + xip[1])) +
                  c_a * (xi0pp[1] + F3[1])) +
                 ct * (xipd[1] + Q20[1])) *
                wg_data[n];
            mu_aj__data[b_loop_ub +
                        mu_aj_->size[0] * ((int)(b_varTheta + 3.0) - 1)] =
                (((a * (xippd[2] + xi0p[2]) + b_a * (rr[2] + xip[2])) +
                  c_a * (xi0pp[2] + F3[2])) +
                 ct * (xipd[2] + Q20[2])) *
                wg_data[n];
          }
        }
        if (d2 > d1) {
          i3 = 0;
          i4 = 0;
          i5 = 0;
          i6 = 0;
        } else {
          i3 = (int)d2 - 1;
          i4 = (int)d1;
          i5 = (int)d2 - 1;
          i6 = (int)d1;
        }
        i7 = r3->size[0];
        r3->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r3, i7);
        r1 = r3->data;
        b_loop_ub = Bbar->size[0];
        for (i7 = 0; i7 < b_loop_ub; i7++) {
          r1[i7] = (int)Bbar_data[i7] - 1;
        }
        b_loop_ub = i4 - i3;
        if ((Bbar->size[0] == mu_aj_->size[0]) &&
            (b_loop_ub == mu_aj_->size[1])) {
          i4 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = Bbar->size[0];
          b_bepscomma->size[1] = b_loop_ub;
          emxEnsureCapacity_real_T(b_bepscomma, i4);
          B_data = b_bepscomma->data;
          for (i4 = 0; i4 < b_loop_ub; i4++) {
            c_loop_ub = Bbar->size[0];
            for (i7 = 0; i7 < c_loop_ub; i7++) {
              B_data[i7 + b_bepscomma->size[0] * i4] =
                  mu_aj_data[((int)Bbar_data[i7] + mu_aj->size[0] * (i3 + i4)) -
                             1] +
                  mu_aj__data[i7 + mu_aj_->size[0] * i4];
            }
          }
          b_i = r3->size[0];
          c_loop_ub = i6 - i5;
          for (i3 = 0; i3 < c_loop_ub; i3++) {
            for (i4 = 0; i4 < b_i; i4++) {
              mu_aj_data[r1[i4] + mu_aj->size[0] * (i5 + i3)] =
                  B_data[i4 + b_i * i3];
            }
          }
        } else {
          binary_expand_op_7(mu_aj, r3, i5, i6, Bbar, i3, i4 - 1, mu_aj_);
          mu_aj_data = mu_aj->data;
        }
        /* %%% mu_ab %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        /*  functions related to chi */
        b_a = 0.0;
        ct = 0.0;
        L_chiO = 0.0;
        b_tau0 = 0.0;
        a = 0.0;
        c_a = 0.0;
        for (i3 = 0; i3 < 3; i3++) {
          b_a += ((b_L_chiO[0] * E22[3 * i3] + b_L_chiO[1] * E22[3 * i3 + 1]) +
                  b_L_chiO[2] * E22[3 * i3 + 2]) *
                 b_L_chiO[i3];
          ct += mm[i3] * XX_contents[i3];
          b_d = b_L_chiO[i3];
          d1 = rrd[i3];
          L_chiO += b_d * d1;
          d2 = tau0.contents[i3];
          b_tau0 += d2 * d1;
          d1 = D22WW[i3];
          a += b_d * d1;
          c_a += d2 * d1;
        }
        ct += b_a;
        for (ddc = 0; ddc < i2; ddc++) {
          for (b_loop_ub = 0; b_loop_ub < i2; b_loop_ub++) {
            c_loop_ub = (int)(nel_ + (double)b_loop_ub) - 1;
            c_tmp = Bbrev_data[c_loop_ub];
            b_i = (int)(nel_ + (double)ddc) - 1;
            b_varTheta = Bbrev_data[b_i];
            tau0_tmp = Bbrevp_data[c_loop_ub];
            st = Bbrevp_data[b_i];
            mu_ab__data[b_loop_ub + mu_ab_->size[0] * ddc] =
                (((c_tmp * b_varTheta * (ct + a) + c_tmp * st * L_chiO) +
                  tau0_tmp * b_varTheta * (L_chiO + c_a)) +
                 tau0_tmp * st * b_tau0) *
                wg_data[n];
          }
        }
        i3 = r3->size[0];
        r3->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r3, i3);
        r1 = r3->data;
        b_loop_ub = Bbar->size[0];
        i3 = r4->size[0];
        r4->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r4, i3);
        r5 = r4->data;
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          i4 = (int)Bbar_data[i3] - 1;
          r1[i3] = i4;
          r5[i3] = i4;
        }
        if ((Bbar->size[0] == mu_ab_->size[0]) &&
            (Bbar->size[0] == mu_ab_->size[1])) {
          i3 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = Bbar->size[0];
          b_bepscomma->size[1] = Bbar->size[0];
          emxEnsureCapacity_real_T(b_bepscomma, i3);
          B_data = b_bepscomma->data;
          b_loop_ub = Bbar->size[0];
          for (i3 = 0; i3 < b_loop_ub; i3++) {
            c_loop_ub = Bbar->size[0];
            for (i4 = 0; i4 < c_loop_ub; i4++) {
              B_data[i4 + b_bepscomma->size[0] * i3] =
                  mu_ab_data[((int)Bbar_data[i4] +
                              mu_ab->size[0] * ((int)Bbar_data[i3] - 1)) -
                             1] +
                  mu_ab__data[i4 + mu_ab_->size[0] * i3];
            }
          }
          b_i = r3->size[0];
          b_loop_ub = r4->size[0];
          for (i3 = 0; i3 < b_loop_ub; i3++) {
            for (i4 = 0; i4 < b_i; i4++) {
              mu_ab_data[r1[i4] + mu_ab->size[0] * r5[i3]] =
                  B_data[i4 + b_i * i3];
            }
          }
        } else {
          binary_expand_op_5(mu_ab, r3, r4, Bbar, mu_ab_);
          mu_ab_data = mu_ab->data;
        }
        /* %%% mu_ajd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        /*  computed at the end as simply F_ibd' */
        /* %%% mu_abd %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% */
        st = 0.0;
        a = 0.0;
        b_a = 0.0;
        for (i3 = 0; i3 < 3; i3++) {
          b_d = D22[3 * i3];
          d1 = D22[3 * i3 + 1];
          d2 = D22[3 * i3 + 2];
          d3 = (b_L_chiO[0] * b_d + b_L_chiO[1] * d1) + b_L_chiO[2] * d2;
          st += d3 * b_L_chiO[i3];
          a += d3 * tau0.contents[i3];
          b_a += ((tau0.contents[0] * b_d + tau0.contents[1] * d1) +
                  tau0.contents[2] * d2) *
                 tau0.contents[i3];
        }
        for (ddc = 0; ddc < i2; ddc++) {
          for (b_loop_ub = 0; b_loop_ub < i2; b_loop_ub++) {
            c_loop_ub = (int)(nel_ + (double)b_loop_ub) - 1;
            ct = Bbrev_data[c_loop_ub];
            b_i = (int)(nel_ + (double)ddc) - 1;
            tau0_tmp = Bbrev_data[b_i];
            c_tmp = Bbrevp_data[c_loop_ub];
            b_varTheta = Bbrevp_data[b_i];
            mu_abd__data[b_loop_ub + mu_abd_->size[0] * ddc] =
                (((ct * tau0_tmp * st + ct * b_varTheta * a) +
                  c_tmp * tau0_tmp * a) +
                 c_tmp * b_varTheta * b_a) *
                wg_data[n];
          }
        }
        i3 = r3->size[0];
        r3->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r3, i3);
        r1 = r3->data;
        b_loop_ub = Bbar->size[0];
        i3 = r4->size[0];
        r4->size[0] = Bbar->size[0];
        emxEnsureCapacity_int32_T(r4, i3);
        r5 = r4->data;
        for (i3 = 0; i3 < b_loop_ub; i3++) {
          i4 = (int)Bbar_data[i3] - 1;
          r1[i3] = i4;
          r5[i3] = i4;
        }
        if ((Bbar->size[0] == mu_abd_->size[0]) &&
            (Bbar->size[0] == mu_abd_->size[1])) {
          i3 = b_bepscomma->size[0] * b_bepscomma->size[1];
          b_bepscomma->size[0] = Bbar->size[0];
          b_bepscomma->size[1] = Bbar->size[0];
          emxEnsureCapacity_real_T(b_bepscomma, i3);
          B_data = b_bepscomma->data;
          b_loop_ub = Bbar->size[0];
          for (i3 = 0; i3 < b_loop_ub; i3++) {
            c_loop_ub = Bbar->size[0];
            for (i4 = 0; i4 < c_loop_ub; i4++) {
              B_data[i4 + b_bepscomma->size[0] * i3] =
                  mu_abd_data[((int)Bbar_data[i4] +
                               mu_abd->size[0] * ((int)Bbar_data[i3] - 1)) -
                              1] +
                  mu_abd__data[i4 + mu_abd_->size[0] * i3];
            }
          }
          b_i = r3->size[0];
          b_loop_ub = r4->size[0];
          for (i3 = 0; i3 < b_loop_ub; i3++) {
            for (i4 = 0; i4 < b_i; i4++) {
              mu_abd_data[r1[i4] + mu_abd->size[0] * r5[i3]] =
                  B_data[i4 + b_i * i3];
            }
          }
        } else {
          binary_expand_op_5(mu_abd, r3, r4, Bbar, mu_abd_);
          mu_abd_data = mu_abd->data;
        }
      }
    }
    if (pass + 1 == 1) {
      i3 = cNbar->size[0];
      cNbar->size[0] = bNbar->size[0];
      emxEnsureCapacity_real_T(cNbar, i3);
      cNbar_data = cNbar->data;
      b_loop_ub = bNbar->size[0];
      for (i3 = 0; i3 < b_loop_ub; i3++) {
        cNbar_data[i3] = bNbar_data[i3];
      }
      b_mldivide(Mbar, cNbar);
      cNbar_data = cNbar->data;
    }
  }
  emxFree_real_T(&b_bepscomma);
  emxFree_real_T(&b_bepsbar);
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
  emxFree_int32_T(&r4);
  emxFree_int32_T(&r3);
  emxFree_real_T(&r2);
  emxFree_int32_T(&r);
  emxFree_real_T(&Bbrevp);
  emxFree_real_T(&Bbrev);
  emxFree_real_T(&Bpp);
  emxFree_real_T(&B);
  emxFree_real_T(&mu_abd_);
  emxFree_real_T(&mu_ab_);
  emxFree_real_T(&mu_aj_);
  emxFree_real_T(&F_ibd_);
  emxFree_real_T(&F_ijd_);
  emxFree_real_T(&F_ib_);
  emxFree_real_T(&F_);
  emxFree_real_T(&cNbar);
  emxFree_real_T(&bNbar);
  emxFree_real_T(&Bbar);
  emxFree_real_T(&Bp);
  emxFree_real_T(&bepsdbar);
  emxFree_real_T(&bepsbar);
  c_mtimes(bepscomma, Dbar11, bepscomma_);
  c_mtimes(bepscomma, Kbar11, bepsdcomma_);
  d_mtimes(bepsdcomma_, bepscomma, F_ij_);
  F_ij__data = F_ij_->data;
  d_mtimes(bepscomma_, bepsdcomma, bepsdcomma_);
  bepsdcomma__data = bepsdcomma_->data;
  emxFree_real_T(&bepsdcomma);
  if (F_ij->size[0] == 1) {
    i = F_ij_->size[0];
  } else {
    i = F_ij->size[0];
  }
  if (F_ij->size[1] == 1) {
    i1 = F_ij_->size[1];
  } else {
    i1 = F_ij->size[1];
  }
  if ((F_ij->size[0] == F_ij_->size[0]) && (F_ij->size[1] == F_ij_->size[1]) &&
      (i == bepsdcomma_->size[0]) && (i1 == bepsdcomma_->size[1])) {
    loop_ub = F_ij->size[0] * F_ij->size[1];
    for (i = 0; i < loop_ub; i++) {
      F_ij_data[i] = (F_ij_data[i] + F_ij__data[i]) + bepsdcomma__data[i];
    }
  } else {
    binary_expand_op_15(F_ij, F_ij_, bepsdcomma_);
  }
  emxFree_real_T(&F_ij_);
  d_mtimes(bepscomma_, bepscomma, bepsdcomma_);
  bepsdcomma__data = bepsdcomma_->data;
  emxFree_real_T(&bepscomma_);
  emxFree_real_T(&bepscomma);
  if ((F_ijd->size[0] == bepsdcomma_->size[0]) &&
      (F_ijd->size[1] == bepsdcomma_->size[1])) {
    loop_ub = F_ijd->size[0] * F_ijd->size[1];
    for (i = 0; i < loop_ub; i++) {
      F_ijd_data[i] += bepsdcomma__data[i];
    }
  } else {
    plus(F_ijd, bepsdcomma_);
  }
  emxFree_real_T(&bepsdcomma_);
  i = mu_ajd->size[0] * mu_ajd->size[1];
  mu_ajd->size[0] = F_ibd->size[1];
  mu_ajd->size[1] = F_ibd->size[0];
  emxEnsureCapacity_real_T(mu_ajd, i);
  B_data = mu_ajd->data;
  loop_ub = F_ibd->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = F_ibd->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      B_data[i1 + mu_ajd->size[0] * i] = F_ibd_data[i + F_ibd->size[0] * i1];
    }
  }
}

void binary_expand_op_14(emxArray_real_T *in1, int in2, int in3, int in4,
                         int in5, const emxArray_real_T *in6)
{
  emxArray_real_T *b_in1;
  const double *in6_data;
  double *b_in1_data;
  double *in1_data;
  int i;
  int stride_0_1;
  int stride_1_1;
  int unnamed_idx_1;
  in6_data = in6->data;
  in1_data = in1->data;
  unnamed_idx_1 = in3 - in2;
  emxInit_real_T(&b_in1, 2);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = 1;
  b_in1->size[1] = unnamed_idx_1;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_1 = ((in5 - in4) + 1 != 1);
  stride_1_1 = (in6->size[0] != 1);
  for (i = 0; i < unnamed_idx_1; i++) {
    b_in1_data[i] = in1_data[in4 + i * stride_0_1] + in6_data[i * stride_1_1];
  }
  unnamed_idx_1 = b_in1->size[1];
  for (i = 0; i < unnamed_idx_1; i++) {
    in1_data[in2 + i] = b_in1_data[i];
  }
  emxFree_real_T(&b_in1);
}

void binary_expand_op_9(emxArray_real_T *in1, int in2, int in3, int in4,
                        int in5, int in6, int in7, int in8, int in9,
                        const emxArray_real_T *in10)
{
  emxArray_real_T *b_in1;
  const double *in10_data;
  double *b_in1_data;
  double *in1_data;
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
  in10_data = in10->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 2);
  i = (in7 - in6) + 1;
  if (in10->size[0] == 1) {
    loop_ub = i;
  } else {
    loop_ub = in10->size[0];
  }
  i1 = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  stride_0_1 = (in9 - in8) + 1;
  if (in10->size[1] == 1) {
    b_loop_ub = stride_0_1;
  } else {
    b_loop_ub = in10->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(b_in1, i1);
  b_in1_data = b_in1->data;
  stride_0_0 = (i != 1);
  stride_0_1 = (stride_0_1 != 1);
  stride_1_0 = (in10->size[0] != 1);
  stride_1_1 = (in10->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[(in6 + i1 * stride_0_0) + in1->size[0] * (in8 + aux_0_1)] +
          in10_data[i1 * stride_1_0 + in10->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  stride_0_1 = in3 - in2;
  stride_1_0 = in5 - in4;
  for (i = 0; i < stride_1_0; i++) {
    for (i1 = 0; i1 < stride_0_1; i1++) {
      in1_data[(in2 + i1) + in1->size[0] * (in4 + i)] =
          b_in1_data[i1 + stride_0_1 * i];
    }
  }
  emxFree_real_T(&b_in1);
}

/* End of code generation (CableForce.c) */
