/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getCableDeformedShape.c
 *
 * Code generation for function 'getCableDeformedShape'
 *
 */

/* Include files */
#include "getCableDeformedShape.h"
#include "CableBCTransinCoord.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_initialize.h"
#include "CableForceRotBCinCoord_types.h"

/* Function Definitions */
void getCableDeformedShape(const double d1[3], const double phi1[3],
                           double gamm1, const double d2[3],
                           const double phi2[3], double gamm2,
                           const emxArray_real_T *Pmid, const double x01[3],
                           const double RJ1[9], const double RE1[9],
                           const double r1[3], const double Rb10[9],
                           const double x02[3], const double RJ2[9],
                           const double RE2[9], const double r2[3],
                           const double Rb20[9], emxArray_real_T *P)
{
  double b_d1[7];
  double q1[7];
  double q2[7];
  const double *Pmid_data;
  double *P_data;
  int i;
  int loop_ub;
  if (!isInitialized_CableForceRotBCinCoord) {
    CableForceRotBCinCoord_initialize();
  }
  Pmid_data = Pmid->data;
  /*  INPUTS  */
  /*  d1 = displacement of joint1  */
  /*  phi1 = exponential coordinates of rotation of joint */
  /*  gamm1 = distance between first and second control points */
  /*  d2, phi2, gamm2 = same of joint 2 (end) */
  /*  Pmid = displaced control points 3:N-2 (3*(N-4) matrix) */
  /*  x01 = reference position of joint 1 (3x1) vector */
  /*  RJ1 = rotation of joint 1 coordinate frame with respect to global */
  /*  RE1 = rotation of cable end with respect to joint 1 coordinate frame */
  /*  r1 = position of cable end relative to joint 1 in joint coordinate system
   */
  /*      (end offset, 3x1 vector) */
  /*  x02, RJ2, RE2, r2 = same for joint 2 */
  /*  OUTPUT */
  /*  P = displaced control points */
  b_d1[0] = d1[0];
  b_d1[3] = phi1[0];
  b_d1[1] = d1[1];
  b_d1[4] = phi1[1];
  b_d1[2] = d1[2];
  b_d1[5] = phi1[2];
  b_d1[6] = gamm1;
  e_CableBCTransinCoord(b_d1, x01, RJ1, RE1, r1, Rb10, q1);
  b_d1[0] = d2[0];
  b_d1[3] = phi2[0];
  b_d1[1] = d2[1];
  b_d1[4] = phi2[1];
  b_d1[2] = d2[2];
  b_d1[5] = phi2[2];
  b_d1[6] = gamm2;
  e_CableBCTransinCoord(b_d1, x02, RJ2, RE2, r2, Rb20, q2);
  i = P->size[0] * P->size[1];
  P->size[0] = 3;
  P->size[1] = Pmid->size[1] + 4;
  emxEnsureCapacity_real_T(P, i);
  P_data = P->data;
  P_data[0] = q1[0];
  P_data[3] = q1[3];
  P_data[1] = q1[1];
  P_data[4] = q1[4];
  P_data[2] = q1[2];
  P_data[5] = q1[5];
  loop_ub = Pmid->size[1];
  for (i = 0; i < loop_ub; i++) {
    int i1;
    i1 = 3 * (i + 2);
    P_data[i1] = Pmid_data[3 * i];
    P_data[i1 + 1] = Pmid_data[3 * i + 1];
    P_data[i1 + 2] = Pmid_data[3 * i + 2];
  }
  P_data[3 * (Pmid->size[1] + 2)] = q2[3];
  P_data[3 * (Pmid->size[1] + 3)] = q2[0];
  P_data[3 * (Pmid->size[1] + 2) + 1] = q2[4];
  P_data[3 * (Pmid->size[1] + 3) + 1] = q2[1];
  P_data[3 * (Pmid->size[1] + 2) + 2] = q2[5];
  P_data[3 * (Pmid->size[1] + 3) + 2] = q2[2];
}

/* End of code generation (getCableDeformedShape.c) */
