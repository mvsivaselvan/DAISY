/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableForceInputDerivative.c
 *
 * Code generation for function 'CableForceInputDerivative'
 *
 */

/* Include files */
#include "CableForceInputDerivative.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_types.h"
#include <math.h>

/* Function Declarations */
static void binary_expand_op_29(emxArray_real_T *in1, int in2, int in3, int in4,
                                int in5, const emxArray_real_T *in6);

/* Function Definitions */
static void binary_expand_op_29(emxArray_real_T *in1, int in2, int in3, int in4,
                                int in5, const emxArray_real_T *in6)
{
  emxArray_real_T *b_in1;
  const double *in6_data;
  double *b_in1_data;
  double *in1_data;
  int i;
  int i1;
  int in3_idx_0;
  int stride_0_0;
  int stride_1_0;
  in6_data = in6->data;
  in1_data = in1->data;
  in3_idx_0 = in3 - in2;
  emxInit_real_T(&b_in1, 2);
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = in3_idx_0;
  b_in1->size[1] = 3;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = ((in5 - in4) + 1 != 1);
  stride_1_0 = (in6->size[0] != 1);
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < in3_idx_0; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[(in4 + i1 * stride_0_0) + in1->size[0] * i] +
          in6_data[i1 * stride_1_0 + in6->size[0] * i];
    }
  }
  in3_idx_0 = b_in1->size[0];
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < in3_idx_0; i1++) {
      in1_data[(in2 + i1) + in1->size[0] * i] =
          b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real_T(&b_in1);
}

void CableForceInputDerivative(const emxArray_real_T *P0, double rho,
                               const emxArray_real_T *wg,
                               const emxArray_real_T *nel,
                               const emxArray_real_T *colmat, double d,
                               emxArray_real_T *dFdu)
{
  emxArray_real_T *B;
  emxArray_real_T *b_dFdu;
  emxArray_real_T *dFdu_;
  emxArray_real_T *r;
  const double *P0_data;
  const double *colmat_data;
  const double *nel_data;
  const double *wg_data;
  double *B_data;
  double *dFdu__data;
  double *dFdu_data;
  double *r1;
  int aoffset;
  int i;
  int i1;
  int i2;
  int k;
  int loop_ub;
  int n;
  colmat_data = colmat->data;
  nel_data = nel->data;
  wg_data = wg->data;
  P0_data = P0->data;
  /*  Computes the derivative of cable force, F, w.r.t. input, u, i.e. dF/du */
  /*  Since dmu/du = 0, it is not computed. */
  /*  INPUTS */
  /*  P0 = reference configuration */
  /*  rho = mass per unit length */
  /*  wg = quadrature weights */
  /*  nel = nel(n) = index of element to which quadrature pt sg(n) belongs to */
  /*  colmat = B-spline basis for position and its derivatives (B) */
  /*  d = degree of B-spline basis for position (B) */
  /*  OUTPUT */
  /*  dFdu = dF/du (a 3Nx3 matrix, since there are 3N force components and */
  /*                3 input components) */
  /*  number of basis functions (or control points) */
  /*  number of quadrature (or collocation) points  */
  i = dFdu->size[0] * dFdu->size[1];
  dFdu->size[0] = 3 * colmat->size[1];
  dFdu->size[1] = 3;
  emxEnsureCapacity_real_T(dFdu, i);
  dFdu_data = dFdu->data;
  loop_ub = 3 * colmat->size[1] * 3;
  for (i = 0; i < loop_ub; i++) {
    dFdu_data[i] = 0.0;
  }
  /*  3DOF per control point, 3 inputs */
  emxInit_real_T(&dFdu_, 2);
  i = (int)(3.0 * (d + 1.0));
  i1 = dFdu_->size[0] * dFdu_->size[1];
  dFdu_->size[0] = i;
  dFdu_->size[1] = 3;
  emxEnsureCapacity_real_T(dFdu_, i1);
  dFdu__data = dFdu_->data;
  loop_ub = i * 3;
  for (i = 0; i < loop_ub; i++) {
    dFdu__data[i] = 0.0;
  }
  /*  temp storage when computing dFdu */
  i = wg->size[0];
  emxInit_real_T(&B, 1);
  emxInit_real_T(&r, 2);
  emxInit_real_T(&b_dFdu, 2);
  for (n = 0; n < i; n++) {
    double absxk;
    double nxi0p;
    double scale;
    double xi0p_idx_0;
    double xi0p_idx_1;
    double xi0p_idx_2;
    xi0p_idx_0 = 3.0 * (((double)n + 1.0) - 1.0);
    loop_ub = colmat->size[1];
    i1 = B->size[0];
    B->size[0] = colmat->size[1];
    emxEnsureCapacity_real_T(B, i1);
    B_data = B->data;
    i1 = r->size[0] * r->size[1];
    r->size[0] = 1;
    r->size[1] = colmat->size[1];
    emxEnsureCapacity_real_T(r, i1);
    r1 = r->data;
    for (i1 = 0; i1 < loop_ub; i1++) {
      B_data[i1] =
          colmat_data[((int)(xi0p_idx_0 + 1.0) + colmat->size[0] * i1) - 1];
      r1[i1] =
          colmat_data[((int)(xi0p_idx_0 + 2.0) + colmat->size[0] * i1) - 1];
    }
    loop_ub = P0->size[1];
    xi0p_idx_0 = 0.0;
    xi0p_idx_1 = 0.0;
    xi0p_idx_2 = 0.0;
    for (k = 0; k < loop_ub; k++) {
      aoffset = k * 3;
      xi0p_idx_0 += P0_data[aoffset] * r1[k];
      xi0p_idx_1 += P0_data[aoffset + 1] * r1[k];
      xi0p_idx_2 += P0_data[aoffset + 2] * r1[k];
    }
    scale = 3.3121686421112381E-170;
    absxk = fabs(xi0p_idx_0);
    if (absxk > 3.3121686421112381E-170) {
      nxi0p = 1.0;
      scale = absxk;
    } else {
      xi0p_idx_0 = absxk / 3.3121686421112381E-170;
      nxi0p = xi0p_idx_0 * xi0p_idx_0;
    }
    absxk = fabs(xi0p_idx_1);
    if (absxk > scale) {
      xi0p_idx_0 = scale / absxk;
      nxi0p = nxi0p * xi0p_idx_0 * xi0p_idx_0 + 1.0;
      scale = absxk;
    } else {
      xi0p_idx_0 = absxk / scale;
      nxi0p += xi0p_idx_0 * xi0p_idx_0;
    }
    absxk = fabs(xi0p_idx_2);
    if (absxk > scale) {
      xi0p_idx_0 = scale / absxk;
      nxi0p = nxi0p * xi0p_idx_0 * xi0p_idx_0 + 1.0;
      scale = absxk;
    } else {
      xi0p_idx_0 = absxk / scale;
      nxi0p += xi0p_idx_0 * xi0p_idx_0;
    }
    nxi0p = scale * sqrt(nxi0p);
    i1 = (int)(d + 1.0);
    for (aoffset = 0; aoffset < i1; aoffset++) {
      xi0p_idx_0 = -(B_data[(int)(nel_data[n] + (double)aoffset) - 1] * rho *
                     nxi0p * wg_data[n]);
      absxk = (double)aoffset * 3.0;
      for (i2 = 0; i2 < 3; i2++) {
        dFdu__data[((int)(absxk + 1.0) + dFdu_->size[0] * i2) - 1] =
            xi0p_idx_0 * (double)iv[3 * i2];
        dFdu__data[((int)(absxk + 2.0) + dFdu_->size[0] * i2) - 1] =
            xi0p_idx_0 * (double)iv[3 * i2 + 1];
        dFdu__data[((int)(absxk + 3.0) + dFdu_->size[0] * i2) - 1] =
            xi0p_idx_0 * (double)iv[3 * i2 + 2];
      }
    }
    xi0p_idx_0 = nel_data[n];
    absxk = (xi0p_idx_0 - 1.0) * 3.0 + 1.0;
    xi0p_idx_0 = (xi0p_idx_0 + d) * 3.0;
    if (absxk > xi0p_idx_0) {
      i1 = 0;
      i2 = 0;
      k = 0;
      loop_ub = 0;
    } else {
      i1 = (int)absxk - 1;
      i2 = (int)xi0p_idx_0;
      k = (int)absxk - 1;
      loop_ub = (int)xi0p_idx_0;
    }
    if (i2 - i1 == dFdu_->size[0]) {
      aoffset = loop_ub - k;
      i2 = b_dFdu->size[0] * b_dFdu->size[1];
      b_dFdu->size[0] = aoffset;
      b_dFdu->size[1] = 3;
      emxEnsureCapacity_real_T(b_dFdu, i2);
      B_data = b_dFdu->data;
      for (i2 = 0; i2 < 3; i2++) {
        for (loop_ub = 0; loop_ub < aoffset; loop_ub++) {
          B_data[loop_ub + b_dFdu->size[0] * i2] =
              dFdu_data[(i1 + loop_ub) + dFdu->size[0] * i2] +
              dFdu__data[loop_ub + dFdu_->size[0] * i2];
        }
      }
      loop_ub = b_dFdu->size[0];
      for (i1 = 0; i1 < 3; i1++) {
        for (i2 = 0; i2 < loop_ub; i2++) {
          dFdu_data[(k + i2) + dFdu->size[0] * i1] =
              B_data[i2 + b_dFdu->size[0] * i1];
        }
      }
    } else {
      binary_expand_op_29(dFdu, k, loop_ub, i1, i2 - 1, dFdu_);
      dFdu_data = dFdu->data;
    }
  }
  emxFree_real_T(&b_dFdu);
  emxFree_real_T(&r);
  emxFree_real_T(&B);
  emxFree_real_T(&dFdu_);
}

/* End of code generation (CableForceInputDerivative.c) */
