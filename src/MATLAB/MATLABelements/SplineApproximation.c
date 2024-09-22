/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SplineApproximation.c
 *
 * Code generation for function 'SplineApproximation'
 *
 */

/* Include files */
#include "SplineApproximation.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_initialize.h"
#include "CableForceRotBCinCoord_types.h"
#include "mldivide.h"
#include <math.h>

/* Function Declarations */
static void binary_expand_op_33(emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3, int in4,
                                const emxArray_real_T *in5);

static void binary_expand_op_34(emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3,
                                const emxArray_real_T *in4, int in5,
                                const emxArray_real_T *in6);

static void minus(emxArray_real_T *in1, const emxArray_real_T *in2);

/* Function Definitions */
static void binary_expand_op_33(emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3, int in4,
                                const emxArray_real_T *in5)
{
  emxArray_real_T *b_in1;
  const double *in2_data;
  const double *in3_data;
  const double *in5_data;
  double b_in3;
  double b_in5;
  double *b_in1_data;
  double *in1_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  in5_data = in5->data;
  in3_data = in3->data;
  in2_data = in2->data;
  in1_data = in1->data;
  b_in3 = in3_data[in4 + 1];
  b_in5 = in5_data[in4];
  emxInit_real_T(&b_in1, 1);
  if (in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = in2->size[0];
  }
  i = b_in1->size[0];
  b_in1->size[0] = loop_ub;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_1_0 = (in2->size[0] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] =
        in1_data[i * stride_0_0] + in2_data[i * stride_1_0] * b_in3 * b_in5;
  }
  i = in1->size[0];
  in1->size[0] = b_in1->size[0];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  loop_ub = b_in1->size[0];
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(&b_in1);
}

static void binary_expand_op_34(emxArray_real_T *in1,
                                const emxArray_real_T *in2,
                                const emxArray_real_T *in3,
                                const emxArray_real_T *in4, int in5,
                                const emxArray_real_T *in6)
{
  emxArray_real_T *b_in1;
  emxArray_real_T *b_in2;
  const double *in2_data;
  const double *in3_data;
  const double *in4_data;
  const double *in6_data;
  double b_in4;
  double b_in6;
  double *b_in1_data;
  double *b_in2_data;
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
  in6_data = in6->data;
  in4_data = in4->data;
  in3_data = in3->data;
  in2_data = in2->data;
  in1_data = in1->data;
  b_in4 = in4_data[in5 + 1];
  emxInit_real_T(&b_in2, 2);
  i = b_in2->size[0] * b_in2->size[1];
  b_in2->size[0] = in2->size[0];
  b_in2->size[1] = in3->size[1];
  emxEnsureCapacity_real_T(b_in2, i);
  b_in2_data = b_in2->data;
  loop_ub = in3->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = in2->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in2_data[i1 + b_in2->size[0] * i] = in2_data[i1] * in3_data[i];
    }
  }
  b_in6 = in6_data[in5];
  emxInit_real_T(&b_in1, 2);
  if (b_in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = b_in2->size[0];
  }
  i = b_in1->size[0] * b_in1->size[1];
  b_in1->size[0] = loop_ub;
  if (b_in2->size[1] == 1) {
    b_loop_ub = in1->size[1];
  } else {
    b_loop_ub = b_in2->size[1];
  }
  b_in1->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_0_1 = (in1->size[1] != 1);
  stride_1_0 = (b_in2->size[0] != 1);
  stride_1_1 = (b_in2->size[1] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * aux_0_1] +
          b_in2_data[i1 * stride_1_0 + b_in2->size[0] * aux_1_1] * b_in4 *
              b_in6;
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
  loop_ub = b_in1->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = b_in1->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      in1_data[i1 + in1->size[0] * i] = b_in1_data[i1 + b_in1->size[0] * i];
    }
  }
  emxFree_real_T(&b_in1);
}

static void minus(emxArray_real_T *in1, const emxArray_real_T *in2)
{
  emxArray_real_T *b_in1;
  const double *in2_data;
  double *b_in1_data;
  double *in1_data;
  int i;
  int loop_ub;
  int stride_0_0;
  int stride_1_0;
  in2_data = in2->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in1, 1);
  if (in2->size[0] == 1) {
    loop_ub = in1->size[0];
  } else {
    loop_ub = in2->size[0];
  }
  i = b_in1->size[0];
  b_in1->size[0] = loop_ub;
  emxEnsureCapacity_real_T(b_in1, i);
  b_in1_data = b_in1->data;
  stride_0_0 = (in1->size[0] != 1);
  stride_1_0 = (in2->size[0] != 1);
  for (i = 0; i < loop_ub; i++) {
    b_in1_data[i] = in1_data[i * stride_0_0] - in2_data[i * stride_1_0];
  }
  i = in1->size[0];
  in1->size[0] = b_in1->size[0];
  emxEnsureCapacity_real_T(in1, i);
  in1_data = in1->data;
  loop_ub = b_in1->size[0];
  for (i = 0; i < loop_ub; i++) {
    in1_data[i] = b_in1_data[i];
  }
  emxFree_real_T(&b_in1);
}

void SplineApproximation(const emxArray_real_T *gamm_,
                         const emxArray_real_T *J_, double N,
                         const emxArray_real_T *xg, const emxArray_real_T *wg,
                         const emxArray_real_T *colmat, emxArray_real_T *p,
                         double *err)
{
  emxArray_int8_T *C;
  emxArray_int8_T *b_varargin_2;
  emxArray_int8_T *varargin_2;
  emxArray_real_T *A;
  emxArray_real_T *B;
  emxArray_real_T *M;
  emxArray_real_T *a;
  emxArray_real_T *b;
  emxArray_real_T *b_x;
  emxArray_real_T *x;
  const double *J__data;
  const double *colmat_data;
  const double *gamm__data;
  const double *wg_data;
  double scale;
  double *B_data;
  double *M_data;
  double *b_data;
  double *b_x_data;
  double *x_data;
  int b_i1;
  int i;
  int i1;
  int kidx;
  int loop_ub;
  int ma;
  int na;
  int sizes_idx_0;
  int sizes_idx_1;
  signed char *C_data;
  signed char *varargin_2_data;
  boolean_T empty_non_axis_sizes;
  boolean_T sizes_idx_1_tmp;
  if (!isInitialized_CableForceRotBCinCoord) {
    CableForceRotBCinCoord_initialize();
  }
  colmat_data = colmat->data;
  wg_data = wg->data;
  J__data = J_->data;
  gamm__data = gamm_->data;
  /*  Approximate a smooth curve by a splne curve */
  /*  Inputs */
  /*  gamm = function handle for a parametric curve */
  /*  N = number of control points in approximating spline curve */
  /*  xg = quadrature points to sample desired curve */
  /*  wg = quadrature weights */
  /*  colmat = B-Spline basis function evaluated at quadrature points */
  /*           The matrix also contains first and second derivatives */
  /*  verbose - 0: no plots; 1: plot gamm and spline spproximation */
  /*  Outputs */
  /*  p = control points (3xN) array */
  /*  knots = uniform knots */
  emxInit_real_T(&M, 2);
  i = M->size[0] * M->size[1];
  M->size[0] = (int)N;
  M->size[1] = (int)N;
  emxEnsureCapacity_real_T(M, i);
  M_data = M->data;
  loop_ub = (int)N * (int)N;
  for (i = 0; i < loop_ub; i++) {
    M_data[i] = 0.0;
  }
  emxInit_real_T(&b, 1);
  loop_ub = (int)(3.0 * N);
  i = b->size[0];
  b->size[0] = loop_ub;
  emxEnsureCapacity_real_T(b, i);
  b_data = b->data;
  for (i = 0; i < loop_ub; i++) {
    b_data[i] = 0.0;
  }
  i = xg->size[0];
  emxInit_real_T(&B, 2);
  emxInit_real_T(&x, 1);
  emxInit_real_T(&a, 1);
  emxInit_real_T(&b_x, 2);
  for (sizes_idx_1 = 0; sizes_idx_1 < i; sizes_idx_1++) {
    i1 = (int)(3.0 * (((double)sizes_idx_1 + 1.0) - 1.0) + 1.0);
    na = B->size[0] * B->size[1];
    B->size[0] = 1;
    kidx = colmat->size[1];
    B->size[1] = colmat->size[1];
    emxEnsureCapacity_real_T(B, na);
    B_data = B->data;
    for (na = 0; na < kidx; na++) {
      B_data[na] = colmat_data[(i1 + colmat->size[0] * na) - 1];
    }
    i1 = x->size[0];
    x->size[0] = B->size[1];
    emxEnsureCapacity_real_T(x, i1);
    b_x_data = x->data;
    kidx = B->size[1];
    for (i1 = 0; i1 < kidx; i1++) {
      b_x_data[i1] = B_data[i1];
    }
    if ((M->size[0] == x->size[0]) && (M->size[1] == B->size[1])) {
      scale = J__data[sizes_idx_1 + 1];
      i1 = b_x->size[0] * b_x->size[1];
      b_x->size[0] = x->size[0];
      b_x->size[1] = B->size[1];
      emxEnsureCapacity_real_T(b_x, i1);
      x_data = b_x->data;
      kidx = B->size[1];
      for (i1 = 0; i1 < kidx; i1++) {
        sizes_idx_0 = x->size[0];
        for (na = 0; na < sizes_idx_0; na++) {
          x_data[na + b_x->size[0] * i1] = b_x_data[na] * B_data[i1];
        }
      }
      kidx = M->size[0] * M->size[1];
      for (i1 = 0; i1 < kidx; i1++) {
        M_data[i1] += x_data[i1] * scale * wg_data[sizes_idx_1];
      }
    } else {
      binary_expand_op_34(M, x, B, J_, sizes_idx_1, wg);
      M_data = M->data;
    }
    /*  index n+1 for gamm_ and J because the  */
    /*  first element is for s=0 (start point) */
    ma = x->size[0];
    i1 = a->size[0];
    a->size[0] = x->size[0] * 3;
    emxEnsureCapacity_real_T(a, i1);
    B_data = a->data;
    kidx = -1;
    for (b_i1 = 0; b_i1 < ma; b_i1++) {
      i1 = 3 * (sizes_idx_1 + 1);
      B_data[kidx + 1] = b_x_data[b_i1] * gamm__data[i1];
      B_data[kidx + 2] = b_x_data[b_i1] * gamm__data[i1 + 1];
      B_data[kidx + 3] = b_x_data[b_i1] * gamm__data[i1 + 2];
      kidx += 3;
    }
    kidx = b->size[0];
    if (b->size[0] == a->size[0]) {
      scale = J__data[sizes_idx_1 + 1];
      for (i1 = 0; i1 < kidx; i1++) {
        b_data[i1] += B_data[i1] * scale * wg_data[sizes_idx_1];
      }
    } else {
      binary_expand_op_33(b, a, J_, sizes_idx_1, wg);
      b_data = b->data;
    }
  }
  emxFree_real_T(&B);
  i = b_x->size[0] * b_x->size[1];
  b_x->size[0] = M->size[0];
  b_x->size[1] = M->size[1];
  emxEnsureCapacity_real_T(b_x, i);
  x_data = b_x->data;
  kidx = M->size[0] * M->size[1];
  for (i = 0; i < kidx; i++) {
    x_data[i] = M_data[i];
  }
  ma = M->size[0];
  na = M->size[1];
  sizes_idx_0 = M->size[0] * 3;
  sizes_idx_1 = M->size[1] * 3;
  i = M->size[0] * M->size[1];
  M->size[0] = sizes_idx_0;
  M->size[1] = sizes_idx_1;
  emxEnsureCapacity_real_T(M, i);
  M_data = M->data;
  kidx = -1;
  for (sizes_idx_0 = 0; sizes_idx_0 < na; sizes_idx_0++) {
    for (sizes_idx_1 = 0; sizes_idx_1 < 3; sizes_idx_1++) {
      for (b_i1 = 0; b_i1 < ma; b_i1++) {
        M_data[kidx + 1] = x_data[b_i1 + b_x->size[0] * sizes_idx_0] *
                           (double)iv[3 * sizes_idx_1];
        M_data[kidx + 2] = x_data[b_i1 + b_x->size[0] * sizes_idx_0] *
                           (double)iv[3 * sizes_idx_1 + 1];
        M_data[kidx + 3] = x_data[b_i1 + b_x->size[0] * sizes_idx_0] *
                           (double)iv[3 * sizes_idx_1 + 2];
        kidx += 3;
      }
    }
  }
  emxInit_int8_T(&C);
  i = C->size[0] * C->size[1];
  C->size[0] = 6;
  C->size[1] = loop_ub;
  emxEnsureCapacity_int8_T(C, i);
  C_data = C->data;
  loop_ub *= 6;
  for (i = 0; i < loop_ub; i++) {
    C_data[i] = 0;
  }
  scale = 3.0 * N - 3.0;
  for (i = 0; i < 3; i++) {
    signed char i2;
    i2 = iv[3 * i];
    C_data[6 * i] = i2;
    kidx = (int)(scale + ((double)i + 1.0)) - 1;
    C_data[6 * kidx + 3] = i2;
    i2 = iv[3 * i + 1];
    C_data[6 * i + 1] = i2;
    C_data[6 * kidx + 4] = i2;
    i2 = iv[3 * i + 2];
    C_data[6 * i + 2] = i2;
    C_data[6 * kidx + 5] = i2;
  }
  emxInit_int8_T(&varargin_2);
  i = varargin_2->size[0] * varargin_2->size[1];
  varargin_2->size[0] = C->size[1];
  varargin_2->size[1] = 6;
  emxEnsureCapacity_int8_T(varargin_2, i);
  varargin_2_data = varargin_2->data;
  loop_ub = C->size[1];
  for (i = 0; i < 6; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      varargin_2_data[i1 + varargin_2->size[0] * i] = C_data[i + 6 * i1];
    }
  }
  sizes_idx_1_tmp = ((M->size[0] != 0) && (M->size[1] != 0));
  if (sizes_idx_1_tmp) {
    kidx = M->size[0];
  } else if (varargin_2->size[0] != 0) {
    kidx = varargin_2->size[0];
  } else {
    kidx = M->size[0];
  }
  empty_non_axis_sizes = (kidx == 0);
  if (empty_non_axis_sizes || sizes_idx_1_tmp) {
    sizes_idx_0 = M->size[1];
  } else {
    sizes_idx_0 = 0;
  }
  if (empty_non_axis_sizes || (varargin_2->size[0] != 0)) {
    sizes_idx_1 = 6;
  } else {
    sizes_idx_1 = 0;
  }
  i = b_x->size[0] * b_x->size[1];
  b_x->size[0] = kidx;
  b_x->size[1] = sizes_idx_0 + sizes_idx_1;
  emxEnsureCapacity_real_T(b_x, i);
  x_data = b_x->data;
  for (i = 0; i < sizes_idx_0; i++) {
    for (i1 = 0; i1 < kidx; i1++) {
      x_data[i1 + b_x->size[0] * i] = M_data[i1 + kidx * i];
    }
  }
  for (i = 0; i < sizes_idx_1; i++) {
    for (i1 = 0; i1 < kidx; i1++) {
      x_data[i1 + b_x->size[0] * (i + sizes_idx_0)] =
          varargin_2_data[i1 + kidx * i];
    }
  }
  emxFree_int8_T(&varargin_2);
  emxInit_int8_T(&b_varargin_2);
  i = b_varargin_2->size[0] * b_varargin_2->size[1];
  b_varargin_2->size[0] = 6;
  b_varargin_2->size[1] = C->size[1] + 6;
  emxEnsureCapacity_int8_T(b_varargin_2, i);
  varargin_2_data = b_varargin_2->data;
  loop_ub = C->size[1];
  for (i = 0; i < loop_ub; i++) {
    for (i1 = 0; i1 < 6; i1++) {
      na = i1 + 6 * i;
      varargin_2_data[na] = C_data[na];
    }
  }
  for (i = 0; i < 6; i++) {
    for (i1 = 0; i1 < 6; i1++) {
      varargin_2_data[i1 + 6 * (i + C->size[1])] = 0;
    }
  }
  emxFree_int8_T(&C);
  sizes_idx_1_tmp = ((b_x->size[0] != 0) && (b_x->size[1] != 0));
  if (sizes_idx_1_tmp) {
    sizes_idx_1 = b_x->size[1];
  } else {
    sizes_idx_1 = b_varargin_2->size[1];
  }
  i = x->size[0];
  x->size[0] = b->size[0] + 6;
  emxEnsureCapacity_real_T(x, i);
  b_x_data = x->data;
  loop_ub = b->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_x_data[i] = b_data[i];
  }
  b_x_data[b->size[0]] = gamm__data[0];
  b_x_data[b->size[0] + 3] = gamm__data[3 * (gamm_->size[1] - 1)];
  b_x_data[b->size[0] + 1] = gamm__data[1];
  b_x_data[b->size[0] + 4] = gamm__data[3 * (gamm_->size[1] - 1) + 1];
  b_x_data[b->size[0] + 2] = gamm__data[2];
  b_x_data[b->size[0] + 5] = gamm__data[3 * (gamm_->size[1] - 1) + 2];
  if (sizes_idx_1_tmp) {
    kidx = b_x->size[0];
  } else {
    kidx = 0;
  }
  emxInit_real_T(&A, 2);
  i = A->size[0] * A->size[1];
  A->size[0] = kidx + 6;
  A->size[1] = sizes_idx_1;
  emxEnsureCapacity_real_T(A, i);
  B_data = A->data;
  for (i = 0; i < sizes_idx_1; i++) {
    for (i1 = 0; i1 < kidx; i1++) {
      B_data[i1 + A->size[0] * i] = x_data[i1 + kidx * i];
    }
    for (i1 = 0; i1 < 6; i1++) {
      B_data[(i1 + kidx) + A->size[0] * i] = varargin_2_data[i1 + 6 * i];
    }
  }
  emxFree_real_T(&b_x);
  emxFree_int8_T(&b_varargin_2);
  b_mldivide(A, x);
  b_x_data = x->data;
  emxFree_real_T(&A);
  sizes_idx_0 = M->size[0] - 1;
  na = M->size[1];
  i = a->size[0];
  a->size[0] = M->size[0];
  emxEnsureCapacity_real_T(a, i);
  B_data = a->data;
  for (sizes_idx_1 = 0; sizes_idx_1 <= sizes_idx_0; sizes_idx_1++) {
    B_data[sizes_idx_1] = 0.0;
  }
  for (ma = 0; ma < na; ma++) {
    kidx = ma * M->size[0];
    for (sizes_idx_1 = 0; sizes_idx_1 <= sizes_idx_0; sizes_idx_1++) {
      B_data[sizes_idx_1] += M_data[kidx + sizes_idx_1] * b_x_data[ma];
    }
  }
  emxFree_real_T(&M);
  if (a->size[0] == b->size[0]) {
    loop_ub = a->size[0];
    for (i = 0; i < loop_ub; i++) {
      B_data[i] -= b_data[i];
    }
  } else {
    minus(a, b);
    B_data = a->data;
  }
  emxFree_real_T(&b);
  if (a->size[0] == 0) {
    *err = 0.0;
  } else {
    *err = 0.0;
    if (a->size[0] == 1) {
      *err = fabs(B_data[0]);
    } else {
      scale = 3.3121686421112381E-170;
      kidx = a->size[0];
      for (ma = 0; ma < kidx; ma++) {
        double absxk;
        absxk = fabs(B_data[ma]);
        if (absxk > scale) {
          double t;
          t = scale / absxk;
          *err = *err * t * t + 1.0;
          scale = absxk;
        } else {
          double t;
          t = absxk / scale;
          *err += t * t;
        }
      }
      *err = scale * sqrt(*err);
    }
  }
  emxFree_real_T(&a);
  i = p->size[0] * p->size[1];
  p->size[0] = 3;
  p->size[1] = (int)N;
  emxEnsureCapacity_real_T(p, i);
  B_data = p->data;
  loop_ub = 3 * (int)N;
  for (i = 0; i < loop_ub; i++) {
    B_data[i] = b_x_data[i];
  }
  emxFree_real_T(&x);
}

/* End of code generation (SplineApproximation.c) */
