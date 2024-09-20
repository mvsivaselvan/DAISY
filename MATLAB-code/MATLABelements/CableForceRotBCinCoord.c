/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableForceRotBCinCoord.c
 *
 * Code generation for function 'CableForceRotBCinCoord'
 *
 */

/* Include files */
#include "CableForceRotBCinCoord.h"
#include "CableBCTransinCoord.h"
#include "CableForce.h"
#include "CableForceInputDerivative.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_initialize.h"
#include "CableForceRotBCinCoord_types.h"
#include "CableInertiaForce.h"
#include "blkdiag.h"
#include "mtimes.h"

/* Function Declarations */
static void binary_expand_op(emxArray_real_T *in1, const emxArray_real_T *in2,
                             const emxArray_real_T *in3,
                             const emxArray_real_T *in4);

static void binary_expand_op_1(emxArray_real_T *in1, const unsigned int in2[2],
                               const emxArray_real_T *in3, const int in4[2],
                               int in5, const emxArray_real_T *in6,
                               const int in7[2]);

static void binary_expand_op_2(emxArray_real_T *in1, const emxArray_real_T *in2,
                               const int in3[2], int in4,
                               const emxArray_real_T *in5, const int in6[2]);

/* Function Definitions */
static void binary_expand_op(emxArray_real_T *in1, const emxArray_real_T *in2,
                             const emxArray_real_T *in3,
                             const emxArray_real_T *in4)
{
  const double *in2_data;
  const double *in3_data;
  const double *in4_data;
  double *in1_data;
  int i;
  int i1;
  int loop_ub;
  int stride_0_0;
  in4_data = in4->data;
  in3_data = in3->data;
  in2_data = in2->data;
  in1_data = in1->data;
  stride_0_0 = (in2->size[0] != 1);
  loop_ub = in1->size[0];
  for (i = 0; i < 14; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      int i2;
      i2 = i1 * stride_0_0;
      in1_data[i1 + in1->size[0] * i] =
          (in2_data[i2 + in2->size[0] * i] + in3_data[i2 + in3->size[0] * i]) +
          in4_data[i1 + in4->size[0] * i];
    }
  }
}

static void binary_expand_op_1(emxArray_real_T *in1, const unsigned int in2[2],
                               const emxArray_real_T *in3, const int in4[2],
                               int in5, const emxArray_real_T *in6,
                               const int in7[2])
{
  emxArray_real_T *r;
  const double *in3_data;
  const double *in6_data;
  double *in1_data;
  double *r1;
  int aux_1_1;
  int i;
  int i1;
  int in4_idx_0;
  int in7_idx_0;
  int stride_1_0;
  int stride_1_1;
  in6_data = in6->data;
  in3_data = in3->data;
  in4_idx_0 = in4[0];
  in7_idx_0 = in7[0];
  i = in1->size[0] * in1->size[1];
  in1->size[0] = in4_idx_0 + in7_idx_0;
  in1->size[1] = in5;
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  for (i = 0; i < in5; i++) {
    for (i1 = 0; i1 < in4_idx_0; i1++) {
      in1_data[i1 + in1->size[0] * i] = in3_data[i1 + in4_idx_0 * i];
    }
    for (i1 = 0; i1 < in7_idx_0; i1++) {
      in1_data[(i1 + in4_idx_0) + in1->size[0] * i] =
          in6_data[i1 + in7_idx_0 * i];
    }
  }
  emxInit_real_T(&r, 2);
  if (in1->size[0] == 1) {
    in4_idx_0 = (int)in2[0];
  } else {
    in4_idx_0 = in1->size[0];
  }
  i = r->size[0] * r->size[1];
  r->size[0] = in4_idx_0;
  if (in1->size[1] == 1) {
    in7_idx_0 = (int)in2[1];
  } else {
    in7_idx_0 = in1->size[1];
  }
  r->size[1] = in7_idx_0;
  emxEnsureCapacity_real_T(r, i);
  r1 = r->data;
  stride_1_0 = (in1->size[0] != 1);
  stride_1_1 = (in1->size[1] != 1);
  aux_1_1 = 0;
  for (i = 0; i < in7_idx_0; i++) {
    for (i1 = 0; i1 < in4_idx_0; i1++) {
      r1[i1 + r->size[0] * i] =
          in1_data[i1 * stride_1_0 + in1->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
  }
  i = in1->size[0] * in1->size[1];
  in1->size[0] = r->size[0];
  in1->size[1] = r->size[1];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  in4_idx_0 = r->size[1];
  for (i = 0; i < in4_idx_0; i++) {
    in7_idx_0 = r->size[0];
    for (i1 = 0; i1 < in7_idx_0; i1++) {
      in1_data[i1 + in1->size[0] * i] = r1[i1 + r->size[0] * i];
    }
  }
  emxFree_real_T(&r);
}

static void binary_expand_op_2(emxArray_real_T *in1, const emxArray_real_T *in2,
                               const int in3[2], int in4,
                               const emxArray_real_T *in5, const int in6[2])
{
  emxArray_real_T *b_in1;
  emxArray_real_T *b_in2;
  const double *in2_data;
  const double *in5_data;
  double *b_in1_data;
  double *b_in2_data;
  double *in1_data;
  int aux_0_1;
  int aux_1_1;
  int i;
  int i1;
  int in3_idx_0;
  int in6_idx_0;
  int stride_0_0;
  int stride_0_1;
  int stride_1_0;
  int stride_1_1;
  in5_data = in5->data;
  in2_data = in2->data;
  in1_data = in1->data;
  in3_idx_0 = in3[0];
  in6_idx_0 = in6[0];
  emxInit_real_T(&b_in2, 2);
  i = b_in2->size[0] * b_in2->size[1];
  b_in2->size[0] = in3_idx_0 + in6_idx_0;
  b_in2->size[1] = in4;
  emxEnsureCapacity_real_T(b_in2, i);
  b_in2_data = b_in2->data;
  for (i = 0; i < in4; i++) {
    for (i1 = 0; i1 < in3_idx_0; i1++) {
      b_in2_data[i1 + b_in2->size[0] * i] = in2_data[i1 + in3_idx_0 * i];
    }
    for (i1 = 0; i1 < in6_idx_0; i1++) {
      b_in2_data[(i1 + in3_idx_0) + b_in2->size[0] * i] =
          in5_data[i1 + in6_idx_0 * i];
    }
  }
  emxInit_real_T(&b_in1, 2);
  if (b_in2->size[0] == 1) {
    in3_idx_0 = in1->size[0];
  } else {
    in3_idx_0 = b_in2->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = in3_idx_0;
  if (b_in2->size[1] == 1) {
    in6_idx_0 = in1->size[1];
  } else {
    in6_idx_0 = b_in2->size[1];
  }
  b_in1->size[1] = in6_idx_0;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_0_1 = (in1->size[1] != 1);
  stride_1_0 = (b_in2->size[0] != 1);
  stride_1_1 = (b_in2->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < in6_idx_0; i++) {
    for (i1 = 0; i1 < in3_idx_0; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * aux_0_1] +
          b_in2_data[i1 * stride_1_0 + b_in2->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_1;
    aux_0_1 += stride_0_1;
  }
  emxFree_real_T(&b_in2);
  i = in1->size[0] * in1->size[1];
  in1->size[0] = b_in1->size[0];
  in1->size[1] = b_in1->size[1];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  in3_idx_0 = b_in1->size[1];
  for (i = 0; i < in3_idx_0; i++) {
    in6_idx_0 = b_in1->size[0];
    for (i1 = 0; i1 < in6_idx_0; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real_T(&b_in1);
}

void CableForceRotBCinCoord(
    const double d1[3], const double phi1[3], double gamm1, const double d2[3],
    const double phi2[3], double gamm2, const emxArray_real_T *Pmid,
    const emxArray_real_T *varThetamid, const emxArray_real_T *P0,
    const double d1dot[3], const double phi1dot[3], double gamm1dot,
    const double d2dot[3], const double phi2dot[3], double gamm2dot,
    const emxArray_real_T *Pmiddot, const emxArray_real_T *varThetamiddot,
    const double d1ddot[3], const double phi1ddot[3], double gamm1ddot,
    const double d2ddot[3], const double phi2ddot[3], double gamm2ddot,
    const emxArray_real_T *Pmidddot, const emxArray_real_T *varThetamidddot,
    const double x01[3], const double RJ1[9], const double RE1[9],
    const double r1[3], const double x02[3], const double RJ2[9],
    const double RE2[9], const double r2[3], const emxArray_real_T *R0,
    const double II[9], double rho, double EA, double EI, double GJ,
    double betAX, double betBEND, double betTOR, const emxArray_real_T *sg,
    const emxArray_real_T *wg, const emxArray_real_T *nel,
    const emxArray_real_T *colmat, const emxArray_real_T *colmat_brev,
    const emxArray_real_T *colmat_bar, double d, double dbrev, double dbar,
    const emxArray_real_T *Mbar, const double u[3],
    const emxArray_real_T *Kbar11, const emxArray_real_T *Dbar11,
    double dynamic, double alph0, emxArray_real_T *Fb, emxArray_real_T *Kb,
    emxArray_real_T *Cb, emxArray_real_T *Mb, emxArray_real_T *Bb)
{
  emxArray_int32_T *c_y;
  emxArray_real_T *Bb2__;
  emxArray_real_T *F;
  emxArray_real_T *F_ib;
  emxArray_real_T *F_ibd;
  emxArray_real_T *F_ij;
  emxArray_real_T *F_ijd;
  emxArray_real_T *Fi;
  emxArray_real_T *P;
  emxArray_real_T *Pdot;
  emxArray_real_T *b_Kb;
  emxArray_real_T *b_q1dot;
  emxArray_real_T *b_r1;
  emxArray_real_T *b_r2;
  emxArray_real_T *b_y;
  emxArray_real_T *c_q1dot;
  emxArray_real_T *dFdu;
  emxArray_real_T *d_Kb;
  emxArray_real_T *d_y;
  emxArray_real_T *e_y;
  emxArray_real_T *mu;
  emxArray_real_T *mu_ab;
  emxArray_real_T *mu_abd;
  emxArray_real_T *mu_aj;
  emxArray_real_T *mu_ajd;
  emxArray_real_T *mui;
  emxArray_real_T *mui_abd;
  emxArray_real_T *mui_abdd;
  emxArray_real_T *mui_ajd;
  emxArray_real_T *mui_ajdd;
  emxArray_real_T *r;
  emxArray_real_T *r3;
  emxArray_real_T *result;
  emxArray_real_T *varTheta;
  emxArray_real_T *varThetadot;
  emxArray_real_T *y;
  double J[196];
  double Qtilde[196];
  double c_Kb[196];
  double dv[196];
  double C1[49];
  double C2[49];
  double Fb1__tmp[49];
  double Fb2__tmp[49];
  double J1[49];
  double J2[49];
  double Q1[49];
  double Q2[49];
  double Qtilde1[49];
  double Qtilde1a[49];
  double Qtilde2[49];
  double Qtilde2a[49];
  double a__1[49];
  double b_Fb1__tmp[21];
  double b_dFdu[21];
  double R01[9];
  double R0end[9];
  double q1dot[7];
  double q2dot[7];
  double qbar1[7];
  double qbar1ddot[7];
  double qbar1dot[7];
  double qbar2[7];
  double qbar2ddot[7];
  double qbar2dot[7];
  const double *Pmid_data;
  const double *Pmidddot_data;
  const double *Pmiddot_data;
  const double *R0_data;
  const double *varThetamid_data;
  const double *varThetamidddot_data;
  const double *varThetamiddot_data;
  double a_tmp;
  double b_d;
  double b_d1;
  double b_sg;
  double *Cb_data;
  double *F_data;
  double *F_ib_data;
  double *F_ibd_data;
  double *F_ij_data;
  double *F_ijd_data;
  double *Fi_data;
  double *Kb_data;
  double *Mb_data;
  double *P_data;
  double *mu_ab_data;
  double *mu_abd_data;
  double *mu_aj_data;
  double *mu_data;
  double *mui_ajd_data;
  double *mui_data;
  double *result_data;
  double *varTheta_data;
  double *varThetadot_data;
  int b_input_sizes[2];
  int c_input_sizes[2];
  int input_sizes[2];
  int sizes[2];
  unsigned int uv[2];
  int b_i;
  int i;
  int i1;
  int i2;
  int input_sizes_idx_0;
  int input_sizes_idx_1;
  int loop_ub;
  int loop_ub_tmp;
  int *y_data;
  boolean_T b;
  boolean_T empty_non_axis_sizes;
  if (!isInitialized_CableForceRotBCinCoord) {
    CableForceRotBCinCoord_initialize();
  }
  R0_data = R0->data;
  varThetamidddot_data = varThetamidddot->data;
  Pmidddot_data = Pmidddot->data;
  varThetamiddot_data = varThetamiddot->data;
  Pmiddot_data = Pmiddot->data;
  varThetamid_data = varThetamid->data;
  Pmid_data = Pmid->data;
  /*  INPUTS */
  /*  d1 = displacement of joint 1  */
  /*  phi1 = exponential coordinates of rotation of joint 1 */
  /*  gamm1 = distance between control points p1 and p2 */
  /*  d2, phi2, gamm2 - same at end 2 (s=1) */
  /*  Pmid = positions of control points 3:(N-2); 3x(N-4) matrix */
  /*  varThetamid = twist DOF 2:(Nbrev-1) */
  /*  P0 = positions of control points in reference configuration */
  /*  d1dot = velocity of end 1 */
  /*  phi1dot = rate of change of phi1 */
  /*  gamm1dot = rate of change of gamma1 */
  /*  d2dot, phi2dot, gamm2dot - same at end 2 (s=1) */
  /*  Pmiddot = rate of change of control points 3:(N-2)  */
  /*  varThetamiddot = rate of change of twist DOF 2:(Nbrev-1) */
  /*  d1ddot ... varThetamidddot = corresponding acceleration terms */
  /*  x01 = reference position of joint 1; (3x1) vector */
  /*  RJ1 = rotation of joint 1 coordinate frame with respect to global */
  /*  RE1 = rotation of cable end 1 with respect to joint 1 coordinate frame */
  /*  r1 = position of cable end 1 relative to joint 1 in joint 1 coordinate  */
  /*        system (end offset, 3x1 vector) */
  /*  x02, RJ2, RE2, r2 - same at end 2 (s=1)  */
  /*  R0 = cell array of orientations in the reference configuration; */
  /*       R0{1} = orientation at end1 (s=0); */
  /*       R0{end} = orientation at end2 (s=1); */
  /*       R0(2:end-1) = orientations at the quadrature points (this is what is
   */
  /*                     passed to CableForce and CableInertiaForce) */
  /*  II = body frame mass moment of inertia (3x3 matrix) per unit length */
  /*  Inputs rho ... Dbar11 are same as in the CableForce function */
  /*  dynamic = binary argument to define static(0)/dynamic(1) state  */
  /*  alph0 = mass-proportional damping coefficient (JB Model) */
  /*  OUTPUTS */
  /*  Fb = force vector (including INERTIA FORCES) corresponding to the DOF  */
  /*      (d1, phi1, gamm1, d2, phi2, gamm2, Pmid(:), varThetamid) */
  /*      in that order; length of the vector is 3N+Nbrev */
  /*  Kb = stiffness matrix Fb_{i,j} with i,j spanning the above DOF */
  /*  Cb = damping matrix Fb_{i,jdot} */
  /*  Mb = inertia matrix Fb_{i,jddot} */
  /*  Bb = derivate of force w.r.t. input */
  b_sg = ((double)sg->size[0] + 1.0) * 3.0;
  /*  Get the control point DOF */
  qbar1[6] = gamm1;
  qbar2[6] = gamm2;
  qbar1dot[6] = gamm1dot;
  for (i = 0; i < 3; i++) {
    R01[3 * i] = R0_data[3 * i];
    input_sizes_idx_1 = (int)(b_sg + ((double)i + 1.0)) - 1;
    R0end[3 * i] = R0_data[3 * input_sizes_idx_1];
    input_sizes_idx_0 = 3 * i + 1;
    R01[input_sizes_idx_0] = R0_data[input_sizes_idx_0];
    R0end[input_sizes_idx_0] = R0_data[3 * input_sizes_idx_1 + 1];
    input_sizes_idx_0 = 3 * i + 2;
    R01[input_sizes_idx_0] = R0_data[input_sizes_idx_0];
    R0end[input_sizes_idx_0] = R0_data[3 * input_sizes_idx_1 + 2];
    qbar1[i] = d1[i];
    qbar1[i + 3] = phi1[i];
    qbar2[i] = d2[i];
    qbar2[i + 3] = phi2[i];
    qbar1dot[i] = d1dot[i];
    qbar1dot[i + 3] = phi1dot[i];
    qbar2dot[i] = d2dot[i];
    qbar2dot[i + 3] = phi2dot[i];
  }
  qbar2dot[6] = gamm2dot;
  CableBCTransinCoord(qbar1, x01, RJ1, RE1, r1, R01, qbar1dot, qbar1dot, q1dot,
                      J1, a__1, Qtilde1);
  CableBCTransinCoord(qbar2, x02, RJ2, RE2, r2, R0end, qbar2dot, qbar2dot,
                      q2dot, J2, a__1, Qtilde2);
  emxInit_real_T(&P, 2);
  b_i = P->size[0] * P->size[1];
  P->size[0] = 3;
  P->size[1] = Pmid->size[1] + 4;
  emxEnsureCapacity_real_T(P, b_i);
  P_data = P->data;
  P_data[0] = q1dot[0];
  P_data[3] = q1dot[3];
  P_data[1] = q1dot[1];
  P_data[4] = q1dot[4];
  P_data[2] = q1dot[2];
  P_data[5] = q1dot[5];
  loop_ub = Pmid->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    i1 = 3 * (b_i + 2);
    P_data[i1] = Pmid_data[3 * b_i];
    P_data[i1 + 1] = Pmid_data[3 * b_i + 1];
    P_data[i1 + 2] = Pmid_data[3 * b_i + 2];
  }
  P_data[3 * (Pmid->size[1] + 2)] = q2dot[3];
  P_data[3 * (Pmid->size[1] + 3)] = q2dot[0];
  P_data[3 * (Pmid->size[1] + 2) + 1] = q2dot[4];
  P_data[3 * (Pmid->size[1] + 3) + 1] = q2dot[1];
  P_data[3 * (Pmid->size[1] + 2) + 2] = q2dot[5];
  P_data[3 * (Pmid->size[1] + 3) + 2] = q2dot[2];
  emxInit_real_T(&varTheta, 1);
  b_i = varTheta->size[0];
  varTheta->size[0] = varThetamid->size[0] + 2;
  emxEnsureCapacity_real_T(varTheta, b_i);
  varTheta_data = varTheta->data;
  varTheta_data[0] = q1dot[6];
  loop_ub = varThetamid->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    varTheta_data[b_i + 1] = varThetamid_data[b_i];
  }
  varTheta_data[varThetamid->size[0] + 1] = q2dot[6];
  /*  Get the corresponding velocities  */
  for (b_i = 0; b_i < 7; b_i++) {
    b_d = 0.0;
    b_d1 = 0.0;
    for (i1 = 0; i1 < 7; i1++) {
      i = b_i + 7 * i1;
      b_d += J1[i] * qbar1dot[i1];
      b_d1 += J2[i] * qbar2dot[i1];
    }
    q2dot[b_i] = b_d1;
    q1dot[b_i] = b_d;
  }
  emxInit_real_T(&Pdot, 2);
  b_i = Pdot->size[0] * Pdot->size[1];
  Pdot->size[0] = 3;
  Pdot->size[1] = Pmiddot->size[1] + 4;
  emxEnsureCapacity_real_T(Pdot, b_i);
  varTheta_data = Pdot->data;
  varTheta_data[0] = q1dot[0];
  varTheta_data[3] = q1dot[3];
  varTheta_data[1] = q1dot[1];
  varTheta_data[4] = q1dot[4];
  varTheta_data[2] = q1dot[2];
  varTheta_data[5] = q1dot[5];
  loop_ub = Pmiddot->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    i1 = 3 * (b_i + 2);
    varTheta_data[i1] = Pmiddot_data[3 * b_i];
    varTheta_data[i1 + 1] = Pmiddot_data[3 * b_i + 1];
    varTheta_data[i1 + 2] = Pmiddot_data[3 * b_i + 2];
  }
  varTheta_data[3 * (Pmiddot->size[1] + 2)] = q2dot[3];
  varTheta_data[3 * (Pmiddot->size[1] + 3)] = q2dot[0];
  varTheta_data[3 * (Pmiddot->size[1] + 2) + 1] = q2dot[4];
  varTheta_data[3 * (Pmiddot->size[1] + 3) + 1] = q2dot[1];
  varTheta_data[3 * (Pmiddot->size[1] + 2) + 2] = q2dot[5];
  varTheta_data[3 * (Pmiddot->size[1] + 3) + 2] = q2dot[2];
  emxInit_real_T(&varThetadot, 1);
  b_i = varThetadot->size[0];
  varThetadot->size[0] = varThetamiddot->size[0] + 2;
  emxEnsureCapacity_real_T(varThetadot, b_i);
  varThetadot_data = varThetadot->data;
  varThetadot_data[0] = q1dot[6];
  loop_ub = varThetamiddot->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    varThetadot_data[b_i + 1] = varThetamiddot_data[b_i];
  }
  varThetadot_data[varThetamiddot->size[0] + 1] = q2dot[6];
  /*  Get the corresponding accelerations (dynamic case) */
  qbar1ddot[6] = gamm1ddot;
  qbar1ddot[0] = d1ddot[0];
  qbar1ddot[3] = phi1ddot[0];
  qbar2ddot[0] = d2ddot[0];
  qbar2ddot[3] = phi2ddot[0];
  qbar1ddot[1] = d1ddot[1];
  qbar1ddot[4] = phi1ddot[1];
  qbar2ddot[1] = d2ddot[1];
  qbar2ddot[4] = phi2ddot[1];
  qbar1ddot[2] = d1ddot[2];
  qbar1ddot[5] = phi1ddot[2];
  qbar2ddot[2] = d2ddot[2];
  qbar2ddot[5] = phi2ddot[2];
  qbar2ddot[6] = gamm2ddot;
  for (b_i = 0; b_i < 7; b_i++) {
    b_d = 0.0;
    b_d1 = 0.0;
    b_sg = 0.0;
    a_tmp = 0.0;
    for (i1 = 0; i1 < 7; i1++) {
      i = b_i + 7 * i1;
      b_d += J1[i] * qbar1ddot[i1];
      b_d1 += Qtilde1[i] * qbar1dot[i1];
      b_sg += J2[i] * qbar2ddot[i1];
      a_tmp += Qtilde2[i] * qbar2dot[i1];
    }
    q1dot[b_i] = b_d + b_d1;
    q2dot[b_i] = b_sg + a_tmp;
  }
  /*  Get forces and stiffness in control point coordinates */
  /*  also get stiffness, damping and inertia */
  emxInit_real_T(&F, 1);
  emxInit_real_T(&mu, 1);
  emxInit_real_T(&F_ij, 2);
  emxInit_real_T(&F_ib, 2);
  emxInit_real_T(&F_ijd, 2);
  emxInit_real_T(&F_ibd, 2);
  emxInit_real_T(&mu_aj, 2);
  emxInit_real_T(&mu_ab, 2);
  emxInit_real_T(&mu_ajd, 2);
  emxInit_real_T(&mu_abd, 2);
  CableForce(P, P0, Pdot, varTheta, varThetadot, R0, rho, EA, EI, GJ, betAX,
             betBEND, betTOR, wg, nel, colmat, colmat_brev, colmat_bar, d,
             dbrev, dbar, Mbar, u, Kbar11, Dbar11, F, mu, F_ij, F_ib, F_ijd,
             F_ibd, mu_aj, mu_ab, mu_ajd, mu_abd);
  mu_abd_data = mu_abd->data;
  Mb_data = mu_ajd->data;
  mu_ab_data = mu_ab->data;
  mu_aj_data = mu_aj->data;
  F_ibd_data = F_ibd->data;
  F_ijd_data = F_ijd->data;
  F_ib_data = F_ib->data;
  F_ij_data = F_ij->data;
  mu_data = mu->data;
  F_data = F->data;
  emxInit_real_T(&Fi, 1);
  b_i = Fi->size[0];
  Fi->size[0] = F->size[0];
  emxEnsureCapacity_real_T(Fi, b_i);
  Fi_data = Fi->data;
  loop_ub = F->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    Fi_data[b_i] = 0.0;
  }
  emxInit_real_T(&mui, 1);
  b_i = mui->size[0];
  mui->size[0] = mu->size[0];
  emxEnsureCapacity_real_T(mui, b_i);
  mui_data = mui->data;
  loop_ub = mu->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    mui_data[b_i] = 0.0;
  }
  b = ((F_ij->size[0] != 0) && (F_ij->size[1] != 0));
  if (b) {
    i = F_ij->size[0];
  } else if ((F_ib->size[0] != 0) && (F_ib->size[1] != 0)) {
    i = F_ib->size[0];
  } else {
    i = F_ij->size[0];
    if (F_ib->size[0] > F_ij->size[0]) {
      i = F_ib->size[0];
    }
  }
  empty_non_axis_sizes = (i == 0);
  if (empty_non_axis_sizes || b) {
    input_sizes_idx_1 = F_ij->size[1];
  } else {
    input_sizes_idx_1 = 0;
  }
  if (empty_non_axis_sizes || ((F_ib->size[0] != 0) && (F_ib->size[1] != 0))) {
    sizes[1] = F_ib->size[1];
  } else {
    sizes[1] = 0;
  }
  emxInit_real_T(&result, 2);
  b_i = result->size[0] * result->size[1];
  result->size[0] = i;
  result->size[1] = input_sizes_idx_1 + sizes[1];
  emxEnsureCapacity_real_T(result, b_i);
  result_data = result->data;
  for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
    for (i1 = 0; i1 < i; i1++) {
      result_data[i1 + result->size[0] * b_i] = F_ij_data[i1 + i * b_i];
    }
  }
  loop_ub = sizes[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < i; i1++) {
      result_data[i1 + result->size[0] * (b_i + input_sizes_idx_1)] =
          F_ib_data[i1 + i * b_i];
    }
  }
  b = ((mu_aj->size[0] != 0) && (mu_aj->size[1] != 0));
  if (b) {
    i = mu_aj->size[0];
  } else if ((mu_ab->size[0] != 0) && (mu_ab->size[1] != 0)) {
    i = mu_ab->size[0];
  } else {
    i = mu_aj->size[0];
    if (mu_ab->size[0] > mu_aj->size[0]) {
      i = mu_ab->size[0];
    }
  }
  empty_non_axis_sizes = (i == 0);
  if (empty_non_axis_sizes || b) {
    input_sizes_idx_1 = mu_aj->size[1];
  } else {
    input_sizes_idx_1 = 0;
  }
  if (empty_non_axis_sizes ||
      ((mu_ab->size[0] != 0) && (mu_ab->size[1] != 0))) {
    sizes[1] = mu_ab->size[1];
  } else {
    sizes[1] = 0;
  }
  b_i = F_ij->size[0] * F_ij->size[1];
  F_ij->size[0] = i;
  F_ij->size[1] = input_sizes_idx_1 + sizes[1];
  emxEnsureCapacity_real_T(F_ij, b_i);
  F_ij_data = F_ij->data;
  for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
    for (i1 = 0; i1 < i; i1++) {
      F_ij_data[i1 + F_ij->size[0] * b_i] = mu_aj_data[i1 + i * b_i];
    }
  }
  loop_ub = sizes[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < i; i1++) {
      F_ij_data[i1 + F_ij->size[0] * (b_i + input_sizes_idx_1)] =
          mu_ab_data[i1 + i * b_i];
    }
  }
  b = ((result->size[0] != 0) && (result->size[1] != 0));
  if (b) {
    i = result->size[1];
  } else if ((F_ij->size[0] != 0) && (F_ij->size[1] != 0)) {
    i = F_ij->size[1];
  } else {
    i = result->size[1];
    if (F_ij->size[1] > result->size[1]) {
      i = F_ij->size[1];
    }
  }
  empty_non_axis_sizes = (i == 0);
  if (empty_non_axis_sizes || b) {
    input_sizes_idx_0 = result->size[0];
  } else {
    input_sizes_idx_0 = 0;
  }
  if (empty_non_axis_sizes || ((F_ij->size[0] != 0) && (F_ij->size[1] != 0))) {
    sizes[0] = F_ij->size[0];
  } else {
    sizes[0] = 0;
  }
  input_sizes_idx_1 = input_sizes_idx_0;
  input_sizes_idx_0 = sizes[0];
  b_i = Kb->size[0] * Kb->size[1];
  Kb->size[0] = input_sizes_idx_1 + sizes[0];
  Kb->size[1] = i;
  emxEnsureCapacity_real_T(Kb, b_i);
  Kb_data = Kb->data;
  for (b_i = 0; b_i < i; b_i++) {
    for (i1 = 0; i1 < input_sizes_idx_1; i1++) {
      Kb_data[i1 + Kb->size[0] * b_i] =
          result_data[i1 + input_sizes_idx_1 * b_i];
    }
    for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
      Kb_data[(i1 + input_sizes_idx_1) + Kb->size[0] * b_i] =
          F_ij_data[i1 + input_sizes_idx_0 * b_i];
    }
  }
  b = ((F_ijd->size[0] != 0) && (F_ijd->size[1] != 0));
  if (b) {
    i = F_ijd->size[0];
  } else if ((F_ibd->size[0] != 0) && (F_ibd->size[1] != 0)) {
    i = F_ibd->size[0];
  } else {
    i = F_ijd->size[0];
    if (F_ibd->size[0] > F_ijd->size[0]) {
      i = F_ibd->size[0];
    }
  }
  empty_non_axis_sizes = (i == 0);
  if (empty_non_axis_sizes || b) {
    input_sizes_idx_1 = F_ijd->size[1];
  } else {
    input_sizes_idx_1 = 0;
  }
  if (empty_non_axis_sizes ||
      ((F_ibd->size[0] != 0) && (F_ibd->size[1] != 0))) {
    sizes[1] = F_ibd->size[1];
  } else {
    sizes[1] = 0;
  }
  b_i = result->size[0] * result->size[1];
  result->size[0] = i;
  result->size[1] = input_sizes_idx_1 + sizes[1];
  emxEnsureCapacity_real_T(result, b_i);
  result_data = result->data;
  for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
    for (i1 = 0; i1 < i; i1++) {
      result_data[i1 + result->size[0] * b_i] = F_ijd_data[i1 + i * b_i];
    }
  }
  loop_ub = sizes[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < i; i1++) {
      result_data[i1 + result->size[0] * (b_i + input_sizes_idx_1)] =
          F_ibd_data[i1 + i * b_i];
    }
  }
  b = ((mu_ajd->size[0] != 0) && (mu_ajd->size[1] != 0));
  if (b) {
    i = mu_ajd->size[0];
  } else if ((mu_abd->size[0] != 0) && (mu_abd->size[1] != 0)) {
    i = mu_abd->size[0];
  } else {
    i = mu_ajd->size[0];
    if (mu_abd->size[0] > mu_ajd->size[0]) {
      i = mu_abd->size[0];
    }
  }
  empty_non_axis_sizes = (i == 0);
  if (empty_non_axis_sizes || b) {
    input_sizes_idx_1 = mu_ajd->size[1];
  } else {
    input_sizes_idx_1 = 0;
  }
  if (empty_non_axis_sizes ||
      ((mu_abd->size[0] != 0) && (mu_abd->size[1] != 0))) {
    sizes[1] = mu_abd->size[1];
  } else {
    sizes[1] = 0;
  }
  b_i = F_ij->size[0] * F_ij->size[1];
  F_ij->size[0] = i;
  F_ij->size[1] = input_sizes_idx_1 + sizes[1];
  emxEnsureCapacity_real_T(F_ij, b_i);
  F_ij_data = F_ij->data;
  for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
    for (i1 = 0; i1 < i; i1++) {
      F_ij_data[i1 + F_ij->size[0] * b_i] = Mb_data[i1 + i * b_i];
    }
  }
  loop_ub = sizes[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < i; i1++) {
      F_ij_data[i1 + F_ij->size[0] * (b_i + input_sizes_idx_1)] =
          mu_abd_data[i1 + i * b_i];
    }
  }
  b = ((result->size[0] != 0) && (result->size[1] != 0));
  if (b) {
    i = result->size[1];
  } else if ((F_ij->size[0] != 0) && (F_ij->size[1] != 0)) {
    i = F_ij->size[1];
  } else {
    i = result->size[1];
    if (F_ij->size[1] > result->size[1]) {
      i = F_ij->size[1];
    }
  }
  empty_non_axis_sizes = (i == 0);
  if (empty_non_axis_sizes || b) {
    input_sizes_idx_1 = result->size[0];
  } else {
    input_sizes_idx_1 = 0;
  }
  if (empty_non_axis_sizes || ((F_ij->size[0] != 0) && (F_ij->size[1] != 0))) {
    sizes[0] = F_ij->size[0];
  } else {
    sizes[0] = 0;
  }
  input_sizes_idx_0 = sizes[0];
  b_i = Cb->size[0] * Cb->size[1];
  Cb->size[0] = input_sizes_idx_1 + sizes[0];
  Cb->size[1] = i;
  emxEnsureCapacity_real_T(Cb, b_i);
  Cb_data = Cb->data;
  for (b_i = 0; b_i < i; b_i++) {
    for (i1 = 0; i1 < input_sizes_idx_1; i1++) {
      Cb_data[i1 + Cb->size[0] * b_i] =
          result_data[i1 + input_sizes_idx_1 * b_i];
    }
    for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
      Cb_data[(i1 + input_sizes_idx_1) + Cb->size[0] * b_i] =
          F_ij_data[i1 + input_sizes_idx_0 * b_i];
    }
  }
  uv[0] = (unsigned int)Kb->size[0];
  uv[1] = (unsigned int)Kb->size[1];
  b_i = Mb->size[0] * Mb->size[1];
  Mb->size[0] = Kb->size[0];
  Mb->size[1] = Kb->size[1];
  emxEnsureCapacity_real_T(Mb, b_i);
  Mb_data = Mb->data;
  loop_ub_tmp = Kb->size[0] * Kb->size[1];
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    Mb_data[b_i] = 0.0;
  }
  if (dynamic != 0.0) {
    emxInit_real_T(&b_q1dot, 2);
    b_i = b_q1dot->size[0] * b_q1dot->size[1];
    b_q1dot->size[0] = 3;
    b_q1dot->size[1] = Pmidddot->size[1] + 4;
    emxEnsureCapacity_real_T(b_q1dot, b_i);
    varTheta_data = b_q1dot->data;
    varTheta_data[0] = q1dot[0];
    varTheta_data[3] = q1dot[3];
    varTheta_data[1] = q1dot[1];
    varTheta_data[4] = q1dot[4];
    varTheta_data[2] = q1dot[2];
    varTheta_data[5] = q1dot[5];
    loop_ub = Pmidddot->size[1];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      i1 = 3 * (b_i + 2);
      varTheta_data[i1] = Pmidddot_data[3 * b_i];
      varTheta_data[i1 + 1] = Pmidddot_data[3 * b_i + 1];
      varTheta_data[i1 + 2] = Pmidddot_data[3 * b_i + 2];
    }
    varTheta_data[3 * (Pmidddot->size[1] + 2)] = q2dot[3];
    varTheta_data[3 * (Pmidddot->size[1] + 3)] = q2dot[0];
    varTheta_data[3 * (Pmidddot->size[1] + 2) + 1] = q2dot[4];
    varTheta_data[3 * (Pmidddot->size[1] + 3) + 1] = q2dot[1];
    varTheta_data[3 * (Pmidddot->size[1] + 2) + 2] = q2dot[5];
    varTheta_data[3 * (Pmidddot->size[1] + 3) + 2] = q2dot[2];
    emxInit_real_T(&c_q1dot, 1);
    b_i = c_q1dot->size[0];
    c_q1dot->size[0] = varThetamidddot->size[0] + 2;
    emxEnsureCapacity_real_T(c_q1dot, b_i);
    varTheta_data = c_q1dot->data;
    varTheta_data[0] = q1dot[6];
    loop_ub = varThetamidddot->size[0];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      varTheta_data[b_i + 1] = varThetamidddot_data[b_i];
    }
    varTheta_data[varThetamidddot->size[0] + 1] = q2dot[6];
    emxInit_real_T(&mui_ajd, 2);
    emxInit_real_T(&mui_abd, 2);
    emxInit_real_T(&mui_ajdd, 2);
    emxInit_real_T(&mui_abdd, 2);
    CableInertiaForce(P, P0, Pdot, b_q1dot, varTheta, varThetadot, c_q1dot, R0,
                      rho, II, wg, nel, colmat, colmat_brev, d, dbrev, alph0,
                      Fi, mui, F_ij, F_ib, mu_aj, mu_ab, F_ijd, F_ibd, mu_ajd,
                      mu_abd, mui_ajd, mui_abd, mui_ajdd, mui_abdd);
    varTheta_data = mui_abdd->data;
    P_data = mui_ajdd->data;
    varThetadot_data = mui_abd->data;
    mui_ajd_data = mui_ajd->data;
    mu_abd_data = mu_abd->data;
    Mb_data = mu_ajd->data;
    F_ibd_data = F_ibd->data;
    F_ijd_data = F_ijd->data;
    mu_ab_data = mu_ab->data;
    mu_aj_data = mu_aj->data;
    F_ib_data = F_ib->data;
    F_ij_data = F_ij->data;
    mui_data = mui->data;
    Fi_data = Fi->data;
    emxFree_real_T(&c_q1dot);
    emxFree_real_T(&b_q1dot);
    b = ((F_ij->size[0] != 0) && (F_ij->size[1] != 0));
    if (b) {
      i = F_ij->size[0];
    } else if ((F_ib->size[0] != 0) && (F_ib->size[1] != 0)) {
      i = F_ib->size[0];
    } else {
      i = F_ij->size[0];
      if (F_ib->size[0] > F_ij->size[0]) {
        i = F_ib->size[0];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      input_sizes_idx_1 = F_ij->size[1];
    } else {
      input_sizes_idx_1 = 0;
    }
    if (empty_non_axis_sizes ||
        ((F_ib->size[0] != 0) && (F_ib->size[1] != 0))) {
      sizes[1] = F_ib->size[1];
    } else {
      sizes[1] = 0;
    }
    b_i = result->size[0] * result->size[1];
    result->size[0] = i;
    result->size[1] = input_sizes_idx_1 + sizes[1];
    emxEnsureCapacity_real_T(result, b_i);
    result_data = result->data;
    for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        result_data[i1 + result->size[0] * b_i] = F_ij_data[i1 + i * b_i];
      }
    }
    loop_ub = sizes[1];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        result_data[i1 + result->size[0] * (b_i + input_sizes_idx_1)] =
            F_ib_data[i1 + i * b_i];
      }
    }
    b = ((mu_ajd->size[0] != 0) && (mu_ajd->size[1] != 0));
    if (b) {
      i = mu_ajd->size[0];
    } else if ((mu_abd->size[0] != 0) && (mu_abd->size[1] != 0)) {
      i = mu_abd->size[0];
    } else {
      i = mu_ajd->size[0];
      if (mu_abd->size[0] > mu_ajd->size[0]) {
        i = mu_abd->size[0];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      input_sizes_idx_1 = mu_ajd->size[1];
    } else {
      input_sizes_idx_1 = 0;
    }
    if (empty_non_axis_sizes ||
        ((mu_abd->size[0] != 0) && (mu_abd->size[1] != 0))) {
      sizes[1] = mu_abd->size[1];
    } else {
      sizes[1] = 0;
    }
    b_i = F_ij->size[0] * F_ij->size[1];
    F_ij->size[0] = i;
    F_ij->size[1] = input_sizes_idx_1 + sizes[1];
    emxEnsureCapacity_real_T(F_ij, b_i);
    F_ij_data = F_ij->data;
    for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        F_ij_data[i1 + F_ij->size[0] * b_i] = Mb_data[i1 + i * b_i];
      }
    }
    loop_ub = sizes[1];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        F_ij_data[i1 + F_ij->size[0] * (b_i + input_sizes_idx_1)] =
            mu_abd_data[i1 + i * b_i];
      }
    }
    b = ((result->size[0] != 0) && (result->size[1] != 0));
    if (b) {
      i = result->size[1];
    } else if ((F_ij->size[0] != 0) && (F_ij->size[1] != 0)) {
      i = F_ij->size[1];
    } else {
      i = result->size[1];
      if (F_ij->size[1] > result->size[1]) {
        i = F_ij->size[1];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      input_sizes[0] = result->size[0];
    } else {
      input_sizes[0] = 0;
    }
    if (empty_non_axis_sizes ||
        ((F_ij->size[0] != 0) && (F_ij->size[1] != 0))) {
      sizes[0] = F_ij->size[0];
    } else {
      sizes[0] = 0;
    }
    b_i = input_sizes[0] + sizes[0];
    if ((Kb->size[0] == b_i) && (Kb->size[1] == i)) {
      input_sizes_idx_1 = input_sizes[0];
      input_sizes_idx_0 = sizes[0];
      i1 = F_ib->size[0] * F_ib->size[1];
      F_ib->size[0] = b_i;
      F_ib->size[1] = i;
      emxEnsureCapacity_real_T(F_ib, i1);
      F_ib_data = F_ib->data;
      for (b_i = 0; b_i < i; b_i++) {
        for (i1 = 0; i1 < input_sizes_idx_1; i1++) {
          F_ib_data[i1 + F_ib->size[0] * b_i] =
              result_data[i1 + input_sizes_idx_1 * b_i];
        }
        for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
          F_ib_data[(i1 + input_sizes_idx_1) + F_ib->size[0] * b_i] =
              F_ij_data[i1 + input_sizes_idx_0 * b_i];
        }
      }
      loop_ub = Kb->size[0] * Kb->size[1];
      for (b_i = 0; b_i < loop_ub; b_i++) {
        Kb_data[b_i] += F_ib_data[b_i];
      }
    } else {
      binary_expand_op_2(Kb, result, input_sizes, i, F_ij, sizes);
      Kb_data = Kb->data;
    }
    b = ((mu_aj->size[0] != 0) && (mu_aj->size[1] != 0));
    if (b) {
      i = mu_aj->size[0];
    } else if ((mu_ab->size[0] != 0) && (mu_ab->size[1] != 0)) {
      i = mu_ab->size[0];
    } else {
      i = mu_aj->size[0];
      if (mu_ab->size[0] > mu_aj->size[0]) {
        i = mu_ab->size[0];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      input_sizes_idx_1 = mu_aj->size[1];
    } else {
      input_sizes_idx_1 = 0;
    }
    if (empty_non_axis_sizes ||
        ((mu_ab->size[0] != 0) && (mu_ab->size[1] != 0))) {
      sizes[1] = mu_ab->size[1];
    } else {
      sizes[1] = 0;
    }
    b_i = result->size[0] * result->size[1];
    result->size[0] = i;
    result->size[1] = input_sizes_idx_1 + sizes[1];
    emxEnsureCapacity_real_T(result, b_i);
    result_data = result->data;
    for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        result_data[i1 + result->size[0] * b_i] = mu_aj_data[i1 + i * b_i];
      }
    }
    loop_ub = sizes[1];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        result_data[i1 + result->size[0] * (b_i + input_sizes_idx_1)] =
            mu_ab_data[i1 + i * b_i];
      }
    }
    b = ((mui_ajd->size[0] != 0) && (mui_ajd->size[1] != 0));
    if (b) {
      i = mui_ajd->size[0];
    } else if ((mui_abd->size[0] != 0) && (mui_abd->size[1] != 0)) {
      i = mui_abd->size[0];
    } else {
      i = mui_ajd->size[0];
      if (mui_abd->size[0] > mui_ajd->size[0]) {
        i = mui_abd->size[0];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      input_sizes_idx_1 = mui_ajd->size[1];
    } else {
      input_sizes_idx_1 = 0;
    }
    if (empty_non_axis_sizes ||
        ((mui_abd->size[0] != 0) && (mui_abd->size[1] != 0))) {
      sizes[1] = mui_abd->size[1];
    } else {
      sizes[1] = 0;
    }
    b_i = F_ij->size[0] * F_ij->size[1];
    F_ij->size[0] = i;
    F_ij->size[1] = input_sizes_idx_1 + sizes[1];
    emxEnsureCapacity_real_T(F_ij, b_i);
    F_ij_data = F_ij->data;
    for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        F_ij_data[i1 + F_ij->size[0] * b_i] = mui_ajd_data[i1 + i * b_i];
      }
    }
    emxFree_real_T(&mui_ajd);
    loop_ub = sizes[1];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        F_ij_data[i1 + F_ij->size[0] * (b_i + input_sizes_idx_1)] =
            varThetadot_data[i1 + i * b_i];
      }
    }
    emxFree_real_T(&mui_abd);
    b = ((result->size[0] != 0) && (result->size[1] != 0));
    if (b) {
      i = result->size[1];
    } else if ((F_ij->size[0] != 0) && (F_ij->size[1] != 0)) {
      i = F_ij->size[1];
    } else {
      i = result->size[1];
      if (F_ij->size[1] > result->size[1]) {
        i = F_ij->size[1];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      b_input_sizes[0] = result->size[0];
    } else {
      b_input_sizes[0] = 0;
    }
    if (empty_non_axis_sizes ||
        ((F_ij->size[0] != 0) && (F_ij->size[1] != 0))) {
      sizes[0] = F_ij->size[0];
    } else {
      sizes[0] = 0;
    }
    b_i = b_input_sizes[0] + sizes[0];
    if ((Cb->size[0] == b_i) && (Cb->size[1] == i)) {
      input_sizes_idx_1 = b_input_sizes[0];
      input_sizes_idx_0 = sizes[0];
      i1 = F_ib->size[0] * F_ib->size[1];
      F_ib->size[0] = b_i;
      F_ib->size[1] = i;
      emxEnsureCapacity_real_T(F_ib, i1);
      F_ib_data = F_ib->data;
      for (b_i = 0; b_i < i; b_i++) {
        for (i1 = 0; i1 < input_sizes_idx_1; i1++) {
          F_ib_data[i1 + F_ib->size[0] * b_i] =
              result_data[i1 + input_sizes_idx_1 * b_i];
        }
        for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
          F_ib_data[(i1 + input_sizes_idx_1) + F_ib->size[0] * b_i] =
              F_ij_data[i1 + input_sizes_idx_0 * b_i];
        }
      }
      loop_ub = Cb->size[0] * Cb->size[1];
      for (b_i = 0; b_i < loop_ub; b_i++) {
        Cb_data[b_i] += F_ib_data[b_i];
      }
    } else {
      binary_expand_op_2(Cb, result, b_input_sizes, i, F_ij, sizes);
      Cb_data = Cb->data;
    }
    b = ((F_ijd->size[0] != 0) && (F_ijd->size[1] != 0));
    if (b) {
      i = F_ijd->size[0];
    } else if ((F_ibd->size[0] != 0) && (F_ibd->size[1] != 0)) {
      i = F_ibd->size[0];
    } else {
      i = F_ijd->size[0];
      if (F_ibd->size[0] > F_ijd->size[0]) {
        i = F_ibd->size[0];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      input_sizes_idx_1 = F_ijd->size[1];
    } else {
      input_sizes_idx_1 = 0;
    }
    if (empty_non_axis_sizes ||
        ((F_ibd->size[0] != 0) && (F_ibd->size[1] != 0))) {
      sizes[1] = F_ibd->size[1];
    } else {
      sizes[1] = 0;
    }
    b_i = result->size[0] * result->size[1];
    result->size[0] = i;
    result->size[1] = input_sizes_idx_1 + sizes[1];
    emxEnsureCapacity_real_T(result, b_i);
    result_data = result->data;
    for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        result_data[i1 + result->size[0] * b_i] = F_ijd_data[i1 + i * b_i];
      }
    }
    loop_ub = sizes[1];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        result_data[i1 + result->size[0] * (b_i + input_sizes_idx_1)] =
            F_ibd_data[i1 + i * b_i];
      }
    }
    b = ((mui_ajdd->size[0] != 0) && (mui_ajdd->size[1] != 0));
    if (b) {
      i = mui_ajdd->size[0];
    } else if ((mui_abdd->size[0] != 0) && (mui_abdd->size[1] != 0)) {
      i = mui_abdd->size[0];
    } else {
      i = mui_ajdd->size[0];
      if (mui_abdd->size[0] > mui_ajdd->size[0]) {
        i = mui_abdd->size[0];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      input_sizes_idx_1 = mui_ajdd->size[1];
    } else {
      input_sizes_idx_1 = 0;
    }
    if (empty_non_axis_sizes ||
        ((mui_abdd->size[0] != 0) && (mui_abdd->size[1] != 0))) {
      sizes[1] = mui_abdd->size[1];
    } else {
      sizes[1] = 0;
    }
    b_i = F_ij->size[0] * F_ij->size[1];
    F_ij->size[0] = i;
    F_ij->size[1] = input_sizes_idx_1 + sizes[1];
    emxEnsureCapacity_real_T(F_ij, b_i);
    F_ij_data = F_ij->data;
    for (b_i = 0; b_i < input_sizes_idx_1; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        F_ij_data[i1 + F_ij->size[0] * b_i] = P_data[i1 + i * b_i];
      }
    }
    emxFree_real_T(&mui_ajdd);
    loop_ub = sizes[1];
    for (b_i = 0; b_i < loop_ub; b_i++) {
      for (i1 = 0; i1 < i; i1++) {
        F_ij_data[i1 + F_ij->size[0] * (b_i + input_sizes_idx_1)] =
            varTheta_data[i1 + i * b_i];
      }
    }
    emxFree_real_T(&mui_abdd);
    b = ((result->size[0] != 0) && (result->size[1] != 0));
    if (b) {
      i = result->size[1];
    } else if ((F_ij->size[0] != 0) && (F_ij->size[1] != 0)) {
      i = F_ij->size[1];
    } else {
      i = result->size[1];
      if (F_ij->size[1] > result->size[1]) {
        i = F_ij->size[1];
      }
    }
    empty_non_axis_sizes = (i == 0);
    if (empty_non_axis_sizes || b) {
      c_input_sizes[0] = result->size[0];
    } else {
      c_input_sizes[0] = 0;
    }
    if (empty_non_axis_sizes ||
        ((F_ij->size[0] != 0) && (F_ij->size[1] != 0))) {
      sizes[0] = F_ij->size[0];
    } else {
      sizes[0] = 0;
    }
    b_i = c_input_sizes[0] + sizes[0];
    if (((int)uv[0] == b_i) && ((int)uv[1] == i)) {
      input_sizes_idx_1 = c_input_sizes[0];
      input_sizes_idx_0 = sizes[0];
      i1 = F_ib->size[0] * F_ib->size[1];
      F_ib->size[0] = b_i;
      F_ib->size[1] = i;
      emxEnsureCapacity_real_T(F_ib, i1);
      F_ib_data = F_ib->data;
      for (b_i = 0; b_i < i; b_i++) {
        for (i1 = 0; i1 < input_sizes_idx_1; i1++) {
          F_ib_data[i1 + F_ib->size[0] * b_i] =
              result_data[i1 + input_sizes_idx_1 * b_i];
        }
        for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
          F_ib_data[(i1 + input_sizes_idx_1) + F_ib->size[0] * b_i] =
              F_ij_data[i1 + input_sizes_idx_0 * b_i];
        }
      }
      b_i = Mb->size[0] * Mb->size[1];
      Mb->size[0] = (int)uv[0];
      Mb->size[1] = (int)uv[1];
      emxEnsureCapacity_real_T(Mb, b_i);
      Mb_data = Mb->data;
      for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
        Mb_data[b_i] = F_ib_data[b_i];
      }
    } else {
      binary_expand_op_1(Mb, uv, result, c_input_sizes, i, F_ij, sizes);
      Mb_data = Mb->data;
    }
  }
  emxFree_real_T(&result);
  emxFree_real_T(&mu_abd);
  emxFree_real_T(&mu_ajd);
  emxFree_real_T(&mu_ab);
  emxFree_real_T(&mu_aj);
  emxFree_real_T(&F_ibd);
  emxFree_real_T(&F_ijd);
  emxFree_real_T(&F_ij);
  emxFree_real_T(&Pdot);
  loop_ub = F->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    F_data[b_i] += Fi_data[b_i];
  }
  emxFree_real_T(&Fi);
  loop_ub = mu->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    mu_data[b_i] += mui_data[b_i];
  }
  emxFree_real_T(&mui);
  for (b_i = 0; b_i < 6; b_i++) {
    q2dot[b_i] = F_data[b_i];
  }
  q2dot[6] = mu_data[0];
  for (b_i = 0; b_i < 7; b_i++) {
    for (i1 = 0; i1 < 7; i1++) {
      Fb1__tmp[i1 + 7 * b_i] = J1[b_i + 7 * i1];
    }
  }
  emxInit_real_T(&y, 2);
  mui_ajd_data = y->data;
  a_tmp = 3.0 * (double)P->size[1];
  emxFree_real_T(&P);
  if (a_tmp < a_tmp - 2.0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    b_i = y->size[0] * y->size[1];
    y->size[0] = 1;
    loop_ub = (int)(a_tmp - (a_tmp - 2.0));
    y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(y, b_i);
    mui_ajd_data = y->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      mui_ajd_data[b_i] = (a_tmp - 2.0) + (double)b_i;
    }
  }
  emxInit_real_T(&b_y, 2);
  mu_abd_data = b_y->data;
  if (a_tmp - 3.0 < a_tmp - 5.0) {
    b_y->size[0] = 1;
    b_y->size[1] = 0;
  } else {
    b_i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    loop_ub = (int)((a_tmp - 3.0) - (a_tmp - 5.0));
    b_y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(b_y, b_i);
    mu_abd_data = b_y->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      mu_abd_data[b_i] = (a_tmp - 5.0) + (double)b_i;
    }
  }
  emxInit_int32_T(&c_y, 1);
  b_i = c_y->size[0];
  c_y->size[0] = y->size[1] + b_y->size[1];
  emxEnsureCapacity_int32_T(c_y, b_i);
  y_data = c_y->data;
  loop_ub = y->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    y_data[b_i] = (int)mui_ajd_data[b_i] - 1;
  }
  loop_ub = b_y->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    y_data[b_i + y->size[1]] = (int)mu_abd_data[b_i] - 1;
  }
  b_i = varThetadot->size[0];
  varThetadot->size[0] = c_y->size[0] + 1;
  emxEnsureCapacity_real_T(varThetadot, b_i);
  varThetadot_data = varThetadot->data;
  loop_ub = c_y->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    varThetadot_data[b_i] = F_data[y_data[b_i]];
  }
  varThetadot_data[c_y->size[0]] = mu_data[varTheta->size[0] - 1];
  for (b_i = 0; b_i < 7; b_i++) {
    for (i1 = 0; i1 < 7; i1++) {
      Fb2__tmp[i1 + 7 * b_i] = J2[b_i + 7 * i1];
    }
  }
  if (a_tmp - 6.0 < 7.0) {
    b_i = 0;
    i1 = 0;
  } else {
    b_i = 6;
    i1 = (int)(a_tmp - 6.0);
  }
  if (varTheta->size[0] - 1 < 2) {
    i = 0;
    i2 = 1;
  } else {
    i = 1;
    i2 = varTheta->size[0];
  }
  loop_ub = i1 - b_i;
  input_sizes_idx_1 = Fb->size[0];
  Fb->size[0] = ((loop_ub + i2) - i) + 13;
  emxEnsureCapacity_real_T(Fb, input_sizes_idx_1);
  varTheta_data = Fb->data;
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < 7; input_sizes_idx_1++) {
    b_d = 0.0;
    b_d1 = 0.0;
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < 7; input_sizes_idx_0++) {
      loop_ub_tmp = input_sizes_idx_1 + 7 * input_sizes_idx_0;
      b_d += Fb1__tmp[loop_ub_tmp] * q2dot[input_sizes_idx_0];
      b_d1 += Fb2__tmp[loop_ub_tmp] * varThetadot_data[input_sizes_idx_0];
    }
    varTheta_data[input_sizes_idx_1] = b_d;
    varTheta_data[input_sizes_idx_1 + 7] = b_d1;
  }
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < loop_ub;
       input_sizes_idx_1++) {
    varTheta_data[input_sizes_idx_1 + 14] = F_data[b_i + input_sizes_idx_1];
  }
  emxFree_real_T(&F);
  loop_ub = i2 - i;
  for (i2 = 0; i2 <= loop_ub - 2; i2++) {
    varTheta_data[((i2 + i1) - b_i) + 14] = mu_data[i + i2];
  }
  double a__10[49];
  double a__9[49];
  emxFree_real_T(&mu);
  /*  Rearrangement of DOFs in both Force and Stiffness (CableForceRotBC.m) */
  /*  are not applied here. If needed, Force vector has to be permuted as */
  /*  Fb = Fb([1:6 8:13 7 14:end]); */
  /*  stiffness, damping and inertia needed */
  /*  Calculate auxiliary terms */
  CableBCTransinCoord(qbar1, x01, RJ1, RE1, r1, R01, q2dot, qbar1ddot, q1dot,
                      a__1, Q1, Qtilde1a);
  b_CableBCTransinCoord(qbar2, x02, RJ2, RE2, r2, R0end, varThetadot, qbar2ddot,
                        q1dot, a__1, Q2, Qtilde2a);
  c_CableBCTransinCoord(qbar1, x01, RJ1, RE1, r1, R01, q2dot, qbar1dot,
                        qbar1dot, q1dot, a__1, a__9, a__10, C1);
  d_CableBCTransinCoord(qbar2, x02, RJ2, RE2, r2, R0end, varThetadot, qbar2dot,
                        qbar2dot, q1dot, a__1, a__9, a__10, C2);
  blkdiag(J1, J2, J);
  /*  J = q^i_{ib}  */
  /*  Q_{ib,jb} = Fb_i q^i_{ib,jb} */
  blkdiag(Qtilde1, Qtilde2, Qtilde);
  /*  Qtilde^{i}_{ib} = */
  /*  q^i_{ib,jb} qbardot^{jb} */
  /*  Qtildea^{i}_{ib}  */
  /*  = q^i_{ib,jb} qbarddot^{jb} */
  /*  C^{i}_{ib} = q^i_{ib,jb,kb}qdot^{jb}qdot^{kb} */
  /*  Permute to apply transformation */
  if (a_tmp < a_tmp - 2.0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    b_i = y->size[0] * y->size[1];
    y->size[0] = 1;
    loop_ub = (int)(a_tmp - (a_tmp - 2.0));
    y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(y, b_i);
    mui_ajd_data = y->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      mui_ajd_data[b_i] = (a_tmp - 2.0) + (double)b_i;
    }
  }
  if (a_tmp - 3.0 < a_tmp - 5.0) {
    b_y->size[0] = 1;
    b_y->size[1] = 0;
  } else {
    b_i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    loop_ub = (int)((a_tmp - 3.0) - (a_tmp - 5.0));
    b_y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(b_y, b_i);
    mu_abd_data = b_y->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      mu_abd_data[b_i] = (a_tmp - 5.0) + (double)b_i;
    }
  }
  emxInit_real_T(&d_y, 2);
  varTheta_data = d_y->data;
  if (a_tmp - 6.0 < 7.0) {
    d_y->size[0] = 1;
    d_y->size[1] = 0;
  } else {
    b_i = d_y->size[0] * d_y->size[1];
    d_y->size[0] = 1;
    d_y->size[1] = (int)((a_tmp - 6.0) - 7.0) + 1;
    emxEnsureCapacity_real_T(d_y, b_i);
    varTheta_data = d_y->data;
    loop_ub = (int)((a_tmp - 6.0) - 7.0);
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      varTheta_data[b_i] = (double)b_i + 7.0;
    }
  }
  emxInit_real_T(&e_y, 2);
  P_data = e_y->data;
  b_sg = a_tmp + (double)varTheta->size[0];
  if (b_sg - 1.0 < a_tmp + 2.0) {
    e_y->size[0] = 1;
    e_y->size[1] = 0;
  } else {
    b_i = e_y->size[0] * e_y->size[1];
    e_y->size[0] = 1;
    loop_ub = (int)((b_sg - 1.0) - (a_tmp + 2.0));
    e_y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(e_y, b_i);
    P_data = e_y->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      P_data[b_i] = (a_tmp + 2.0) + (double)b_i;
    }
  }
  /*  Permute Stiffness, Damping, and Mass matrices  */
  b_i = y->size[1] + b_y->size[1];
  i1 = varThetadot->size[0];
  varThetadot->size[0] = ((b_i + d_y->size[1]) + e_y->size[1]) + 8;
  emxEnsureCapacity_real_T(varThetadot, i1);
  varThetadot_data = varThetadot->data;
  for (i1 = 0; i1 < 6; i1++) {
    varThetadot_data[i1] = (double)i1 + 1.0;
  }
  varThetadot_data[6] = a_tmp + 1.0;
  loop_ub = y->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    varThetadot_data[i1 + 7] = mui_ajd_data[i1];
  }
  loop_ub = b_y->size[1];
  for (i1 = 0; i1 < loop_ub; i1++) {
    varThetadot_data[(i1 + y->size[1]) + 7] = mu_abd_data[i1];
  }
  varThetadot_data[b_i + 7] = b_sg;
  loop_ub = d_y->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    varThetadot_data[((b_i + y->size[1]) + b_y->size[1]) + 8] =
        varTheta_data[b_i];
  }
  loop_ub = e_y->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    varThetadot_data[(((b_i + y->size[1]) + b_y->size[1]) + d_y->size[1]) + 8] =
        P_data[b_i];
  }
  emxFree_real_T(&e_y);
  emxFree_real_T(&d_y);
  b_i = F_ib->size[0] * F_ib->size[1];
  F_ib->size[0] = varThetadot->size[0];
  F_ib->size[1] = varThetadot->size[0];
  emxEnsureCapacity_real_T(F_ib, b_i);
  F_ib_data = F_ib->data;
  loop_ub = varThetadot->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    input_sizes_idx_0 = varThetadot->size[0];
    for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
      F_ib_data[i1 + F_ib->size[0] * b_i] =
          Kb_data[((int)varThetadot_data[i1] +
                   Kb->size[0] * ((int)varThetadot_data[b_i] - 1)) -
                  1];
    }
  }
  b_i = Kb->size[0] * Kb->size[1];
  Kb->size[0] = F_ib->size[0];
  Kb->size[1] = F_ib->size[1];
  emxEnsureCapacity_real_T(Kb, b_i);
  Kb_data = Kb->data;
  loop_ub_tmp = F_ib->size[0] * F_ib->size[1];
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    Kb_data[b_i] = F_ib_data[b_i];
  }
  b_i = F_ib->size[0] * F_ib->size[1];
  F_ib->size[0] = varThetadot->size[0];
  F_ib->size[1] = varThetadot->size[0];
  emxEnsureCapacity_real_T(F_ib, b_i);
  F_ib_data = F_ib->data;
  loop_ub = varThetadot->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    input_sizes_idx_0 = varThetadot->size[0];
    for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
      F_ib_data[i1 + F_ib->size[0] * b_i] =
          Cb_data[((int)varThetadot_data[i1] +
                   Cb->size[0] * ((int)varThetadot_data[b_i] - 1)) -
                  1];
    }
  }
  b_i = Cb->size[0] * Cb->size[1];
  Cb->size[0] = F_ib->size[0];
  Cb->size[1] = F_ib->size[1];
  emxEnsureCapacity_real_T(Cb, b_i);
  Cb_data = Cb->data;
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    Cb_data[b_i] = F_ib_data[b_i];
  }
  b_i = F_ib->size[0] * F_ib->size[1];
  F_ib->size[0] = varThetadot->size[0];
  F_ib->size[1] = varThetadot->size[0];
  emxEnsureCapacity_real_T(F_ib, b_i);
  F_ib_data = F_ib->data;
  loop_ub = varThetadot->size[0];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    input_sizes_idx_0 = varThetadot->size[0];
    for (i1 = 0; i1 < input_sizes_idx_0; i1++) {
      F_ib_data[i1 + F_ib->size[0] * b_i] =
          Mb_data[((int)varThetadot_data[i1] +
                   Mb->size[0] * ((int)varThetadot_data[b_i] - 1)) -
                  1];
    }
  }
  emxFree_real_T(&varThetadot);
  b_i = Mb->size[0] * Mb->size[1];
  Mb->size[0] = F_ib->size[0];
  Mb->size[1] = F_ib->size[1];
  emxEnsureCapacity_real_T(Mb, b_i);
  Mb_data = Mb->data;
  for (b_i = 0; b_i < loop_ub_tmp; b_i++) {
    Mb_data[b_i] = F_ib_data[b_i];
  }
  emxFree_real_T(&F_ib);
  /*  Stiffness, Damping and Mass matrices in terms of new DOF */
  emxInit_real_T(&b_Kb, 2);
  b_i = b_Kb->size[0] * b_Kb->size[1];
  b_Kb->size[0] = Kb->size[0];
  b_Kb->size[1] = 14;
  emxEnsureCapacity_real_T(b_Kb, b_i);
  varTheta_data = b_Kb->data;
  loop_ub = Kb->size[0];
  for (b_i = 0; b_i < 14; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      varTheta_data[i1 + b_Kb->size[0] * b_i] = Kb_data[i1 + Kb->size[0] * b_i];
    }
  }
  emxInit_real_T(&r, 2);
  e_mtimes(b_Kb, J, r);
  P_data = r->data;
  b_i = b_Kb->size[0] * b_Kb->size[1];
  b_Kb->size[0] = Cb->size[0];
  b_Kb->size[1] = 14;
  emxEnsureCapacity_real_T(b_Kb, b_i);
  varTheta_data = b_Kb->data;
  loop_ub = Cb->size[0];
  for (b_i = 0; b_i < 14; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      varTheta_data[i1 + b_Kb->size[0] * b_i] = Cb_data[i1 + Cb->size[0] * b_i];
    }
  }
  emxInit_real_T(&b_r1, 2);
  e_mtimes(b_Kb, Qtilde, b_r1);
  varThetadot_data = b_r1->data;
  blkdiag(Qtilde1a, Qtilde2a, dv);
  blkdiag(C1, C2, c_Kb);
  b_i = b_Kb->size[0] * b_Kb->size[1];
  b_Kb->size[0] = Mb->size[0];
  b_Kb->size[1] = 14;
  emxEnsureCapacity_real_T(b_Kb, b_i);
  varTheta_data = b_Kb->data;
  loop_ub = Mb->size[0];
  for (b_i = 0; b_i < 14; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      varTheta_data[i1 + b_Kb->size[0] * b_i] = Mb_data[i1 + Mb->size[0] * b_i];
    }
  }
  for (b_i = 0; b_i < 196; b_i++) {
    dv[b_i] += c_Kb[b_i];
  }
  emxInit_real_T(&b_r2, 2);
  e_mtimes(b_Kb, dv, b_r2);
  varTheta_data = b_r2->data;
  if (r->size[0] == b_r2->size[0]) {
    loop_ub = Kb->size[0];
    for (b_i = 0; b_i < 14; b_i++) {
      for (i1 = 0; i1 < loop_ub; i1++) {
        Kb_data[i1 + Kb->size[0] * b_i] =
            (P_data[i1 + r->size[0] * b_i] +
             varThetadot_data[i1 + b_r1->size[0] * b_i]) +
            varTheta_data[i1 + b_r2->size[0] * b_i];
      }
    }
  } else {
    binary_expand_op(Kb, r, b_r1, b_r2);
    Kb_data = Kb->data;
  }
  emxInit_real_T(&d_Kb, 2);
  b_i = d_Kb->size[0] * d_Kb->size[1];
  d_Kb->size[0] = 14;
  d_Kb->size[1] = Kb->size[1];
  emxEnsureCapacity_real_T(d_Kb, b_i);
  varTheta_data = d_Kb->data;
  loop_ub = Kb->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < 14; i1++) {
      varTheta_data[i1 + 14 * b_i] = Kb_data[i1 + Kb->size[0] * b_i];
    }
  }
  emxInit_real_T(&r3, 2);
  f_mtimes(J, d_Kb, r3);
  P_data = r3->data;
  loop_ub = Kb->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < 14; i1++) {
      Kb_data[i1 + Kb->size[0] * b_i] = P_data[i1 + 14 * b_i];
    }
  }
  blkdiag(Q1, Q2, dv);
  /*  add geometric stiffness */
  b_i = b_Kb->size[0] * b_Kb->size[1];
  b_Kb->size[0] = Cb->size[0];
  b_Kb->size[1] = 14;
  emxEnsureCapacity_real_T(b_Kb, b_i);
  varTheta_data = b_Kb->data;
  loop_ub = Cb->size[0];
  for (b_i = 0; b_i < 14; b_i++) {
    for (i1 = 0; i1 < 14; i1++) {
      Kb_data[i1 + Kb->size[0] * b_i] += dv[i1 + 14 * b_i];
    }
    for (i1 = 0; i1 < loop_ub; i1++) {
      varTheta_data[i1 + b_Kb->size[0] * b_i] = Cb_data[i1 + Cb->size[0] * b_i];
    }
  }
  e_mtimes(b_Kb, J, r);
  P_data = r->data;
  b_i = b_r2->size[0] * b_r2->size[1];
  b_r2->size[0] = Mb->size[0];
  b_r2->size[1] = 14;
  emxEnsureCapacity_real_T(b_r2, b_i);
  varTheta_data = b_r2->data;
  loop_ub = Mb->size[0];
  for (b_i = 0; b_i < 14; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      varTheta_data[i1 + b_r2->size[0] * b_i] =
          2.0 * Mb_data[i1 + Mb->size[0] * b_i];
    }
  }
  e_mtimes(b_r2, Qtilde, b_r1);
  varThetadot_data = b_r1->data;
  emxFree_real_T(&b_r2);
  loop_ub = Cb->size[0];
  for (b_i = 0; b_i < 14; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      Cb_data[i1 + Cb->size[0] * b_i] =
          P_data[i1 + r->size[0] * b_i] +
          varThetadot_data[i1 + b_r1->size[0] * b_i];
    }
  }
  emxFree_real_T(&b_r1);
  b_i = d_Kb->size[0] * d_Kb->size[1];
  d_Kb->size[0] = 14;
  d_Kb->size[1] = Cb->size[1];
  emxEnsureCapacity_real_T(d_Kb, b_i);
  varTheta_data = d_Kb->data;
  loop_ub = Cb->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < 14; i1++) {
      varTheta_data[i1 + 14 * b_i] = Cb_data[i1 + Cb->size[0] * b_i];
    }
  }
  f_mtimes(J, d_Kb, r3);
  P_data = r3->data;
  loop_ub = Cb->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < 14; i1++) {
      Cb_data[i1 + Cb->size[0] * b_i] = P_data[i1 + 14 * b_i];
    }
  }
  b_i = b_Kb->size[0] * b_Kb->size[1];
  b_Kb->size[0] = Mb->size[0];
  b_Kb->size[1] = 14;
  emxEnsureCapacity_real_T(b_Kb, b_i);
  varTheta_data = b_Kb->data;
  loop_ub = Mb->size[0];
  for (b_i = 0; b_i < 14; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      varTheta_data[i1 + b_Kb->size[0] * b_i] = Mb_data[i1 + Mb->size[0] * b_i];
    }
  }
  e_mtimes(b_Kb, J, r);
  P_data = r->data;
  emxFree_real_T(&b_Kb);
  loop_ub = Mb->size[0];
  for (b_i = 0; b_i < 14; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      Mb_data[i1 + Mb->size[0] * b_i] = P_data[i1 + r->size[0] * b_i];
    }
  }
  emxFree_real_T(&r);
  b_i = d_Kb->size[0] * d_Kb->size[1];
  d_Kb->size[0] = 14;
  d_Kb->size[1] = Mb->size[1];
  emxEnsureCapacity_real_T(d_Kb, b_i);
  varTheta_data = d_Kb->data;
  loop_ub = Mb->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < 14; i1++) {
      varTheta_data[i1 + 14 * b_i] = Mb_data[i1 + Mb->size[0] * b_i];
    }
  }
  f_mtimes(J, d_Kb, r3);
  P_data = r3->data;
  emxFree_real_T(&d_Kb);
  loop_ub = Mb->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    for (i1 = 0; i1 < 14; i1++) {
      Mb_data[i1 + Mb->size[0] * b_i] = P_data[i1 + 14 * b_i];
    }
  }
  emxFree_real_T(&r3);
  /*  derivative w.r.t. input needed */
  /*  compute Binp = dF/du */
  emxInit_real_T(&dFdu, 2);
  CableForceInputDerivative(P0, rho, wg, nel, colmat, d, dFdu);
  P_data = dFdu->data;
  /*  dmu/du = 0 */
  /*  Now transform each column of this matrix in the same way as F and mu */
  if (a_tmp < a_tmp - 2.0) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    b_i = y->size[0] * y->size[1];
    y->size[0] = 1;
    loop_ub = (int)(a_tmp - (a_tmp - 2.0));
    y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(y, b_i);
    mui_ajd_data = y->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      mui_ajd_data[b_i] = (a_tmp - 2.0) + (double)b_i;
    }
  }
  if (a_tmp - 3.0 < a_tmp - 5.0) {
    b_y->size[0] = 1;
    b_y->size[1] = 0;
  } else {
    b_i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = 1;
    loop_ub = (int)((a_tmp - 3.0) - (a_tmp - 5.0));
    b_y->size[1] = loop_ub + 1;
    emxEnsureCapacity_real_T(b_y, b_i);
    mu_abd_data = b_y->data;
    for (b_i = 0; b_i <= loop_ub; b_i++) {
      mu_abd_data[b_i] = (a_tmp - 5.0) + (double)b_i;
    }
  }
  b_i = c_y->size[0];
  c_y->size[0] = y->size[1] + b_y->size[1];
  emxEnsureCapacity_int32_T(c_y, b_i);
  y_data = c_y->data;
  loop_ub = y->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    y_data[b_i] = (int)mui_ajd_data[b_i] - 1;
  }
  loop_ub = b_y->size[1];
  for (b_i = 0; b_i < loop_ub; b_i++) {
    y_data[b_i + y->size[1]] = (int)mu_abd_data[b_i] - 1;
  }
  emxFree_real_T(&b_y);
  emxFree_real_T(&y);
  emxInit_real_T(&Bb2__, 2);
  b_i = Bb2__->size[0] * Bb2__->size[1];
  Bb2__->size[0] = c_y->size[0] + 1;
  Bb2__->size[1] = 3;
  emxEnsureCapacity_real_T(Bb2__, b_i);
  varTheta_data = Bb2__->data;
  loop_ub = c_y->size[0];
  for (b_i = 0; b_i < 3; b_i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      varTheta_data[i1 + Bb2__->size[0] * b_i] =
          P_data[y_data[i1] + dFdu->size[0] * b_i];
    }
    varTheta_data[c_y->size[0] + Bb2__->size[0] * b_i] = 0.0;
  }
  emxFree_int32_T(&c_y);
  if (a_tmp - 6.0 < 7.0) {
    b_i = 0;
    i1 = 0;
  } else {
    b_i = 6;
    i1 = (int)(a_tmp - 6.0);
  }
  if (varTheta->size[0] - 1 < 2) {
    i = -15;
    i2 = -16;
  } else {
    i = -14;
    i2 = varTheta->size[0] - 17;
  }
  emxFree_real_T(&varTheta);
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < 3; input_sizes_idx_1++) {
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < 6; input_sizes_idx_0++) {
      b_dFdu[input_sizes_idx_0 + 7 * input_sizes_idx_1] =
          P_data[input_sizes_idx_0 + dFdu->size[0] * input_sizes_idx_1];
    }
    b_dFdu[7 * input_sizes_idx_1 + 6] = 0.0;
  }
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < 7; input_sizes_idx_1++) {
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < 3; input_sizes_idx_0++) {
      b_d = 0.0;
      for (loop_ub_tmp = 0; loop_ub_tmp < 7; loop_ub_tmp++) {
        b_d += Fb1__tmp[input_sizes_idx_1 + 7 * loop_ub_tmp] *
               b_dFdu[loop_ub_tmp + 7 * input_sizes_idx_0];
      }
      b_Fb1__tmp[input_sizes_idx_1 + 7 * input_sizes_idx_0] = b_d;
    }
  }
  for (input_sizes_idx_1 = 0; input_sizes_idx_1 < 7; input_sizes_idx_1++) {
    for (input_sizes_idx_0 = 0; input_sizes_idx_0 < 3; input_sizes_idx_0++) {
      b_d = 0.0;
      for (loop_ub_tmp = 0; loop_ub_tmp < 7; loop_ub_tmp++) {
        b_d += Fb2__tmp[input_sizes_idx_1 + 7 * loop_ub_tmp] *
               varTheta_data[loop_ub_tmp + 7 * input_sizes_idx_0];
      }
      b_dFdu[input_sizes_idx_1 + 7 * input_sizes_idx_0] = b_d;
    }
  }
  emxFree_real_T(&Bb2__);
  loop_ub = i1 - b_i;
  input_sizes_idx_1 = Bb->size[0] * Bb->size[1];
  Bb->size[0] = ((loop_ub + i2) - i) + 15;
  Bb->size[1] = 3;
  emxEnsureCapacity_real_T(Bb, input_sizes_idx_1);
  varTheta_data = Bb->data;
  input_sizes_idx_0 = i2 - i;
  for (i = 0; i < 3; i++) {
    for (i2 = 0; i2 < 7; i2++) {
      input_sizes_idx_1 = i2 + 7 * i;
      varTheta_data[i2 + Bb->size[0] * i] = b_Fb1__tmp[input_sizes_idx_1];
      varTheta_data[(i2 + Bb->size[0] * i) + 7] = b_dFdu[input_sizes_idx_1];
    }
    for (i2 = 0; i2 < loop_ub; i2++) {
      varTheta_data[(i2 + Bb->size[0] * i) + 14] =
          P_data[(b_i + i2) + dFdu->size[0] * i];
    }
    for (i2 = 0; i2 <= input_sizes_idx_0; i2++) {
      varTheta_data[(((i2 + i1) - b_i) + Bb->size[0] * i) + 14] = 0.0;
    }
  }
  emxFree_real_T(&dFdu);
}

/* End of code generation (CableForceRotBCinCoord.c) */
