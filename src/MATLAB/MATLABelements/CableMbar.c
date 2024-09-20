/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableMbar.c
 *
 * Code generation for function 'CableMbar'
 *
 */

/* Include files */
#include "CableMbar.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_initialize.h"
#include "CableForceRotBCinCoord_types.h"
#include "mldivide.h"
#include "mrdivide_helper.h"
#include <math.h>

/* Function Declarations */
static void binary_expand_op_30(emxArray_real_T *in1,
                                const emxArray_real_T *in2, double in3);

/* Function Definitions */
static void binary_expand_op_30(emxArray_real_T *in1,
                                const emxArray_real_T *in2, double in3)
{
  emxArray_real_T *b_in1;
  emxArray_real_T *b_in2;
  const double *in2_data;
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
  int stride_1_0_tmp;
  in2_data = in2->data;
  in1_data = in1->data;
  emxInit_real_T(&b_in2, 2);
  i = b_in2->size[0] * b_in2->size[1];
  b_in2->size[0] = in2->size[0];
  b_in2->size[1] = in2->size[0];
  emxEnsureCapacity_real_T(b_in2, i);
  b_in2_data = b_in2->data;
  loop_ub = in2->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = in2->size[0];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_in2_data[i1 + b_in2->size[0] * i] = in2_data[i1] * in3 * in2_data[i];
    }
  }
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
  stride_1_0_tmp = (b_in2->size[0] != 1);
  aux_0_1 = 0;
  aux_1_1 = 0;
  for (i = 0; i < b_loop_ub; i++) {
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_in1_data[i1 + b_in1->size[0] * i] =
          in1_data[i1 * stride_0_0 + in1->size[0] * aux_0_1] +
          b_in2_data[i1 * stride_1_0_tmp + b_in2->size[0] * aux_1_1];
    }
    aux_1_1 += stride_1_0_tmp;
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

void CableMbar(const emxArray_real_T *P0, double EA, double betAX,
               const emxArray_real_T *colmat, const emxArray_real_T *colmat_bar,
               const emxArray_real_T *wg, emxArray_real_T *Mbar,
               emxArray_real_T *Kbar11, emxArray_real_T *Dbar11)
{
  emxArray_real_T *Bbar;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  const double *P0_data;
  const double *colmat_bar_data;
  const double *colmat_data;
  const double *wg_data;
  double Dm11;
  double *Bbar_data;
  double *Mbar_data;
  int i;
  int i1;
  int k;
  int loop_ub_tmp;
  int n;
  if (!isInitialized_CableForceRotBCinCoord) {
    CableForceRotBCinCoord_initialize();
  }
  wg_data = wg->data;
  colmat_bar_data = colmat_bar->data;
  colmat_data = colmat->data;
  P0_data = P0->data;
  /*  INPUTS */
  /*  P0 = reference configuration */
  /*  EA = section axial rigidity */
  /*  betAX = damping coeff associated with rate of axial deformation */
  /*  colmat = B-spline basis for position and its derivatives (B) */
  /*  colmat_bar = B-spline basis for strain projection (Bbar) */
  /*  wg = quadrature weights */
  /*  OUTPUTS */
  /*  Mbar = "mass matrix" for strain projection */
  /*  Kbar11, Dbar11 = as defined in equation (33) [EQ NUM needs to be updated
   */
  /*                   to be consistent with document] */
  Dm11 = betAX * EA;
  /*  number of strain projection basis functions */
  /*  number of quadrature (or collocation) points  */
  i = Mbar->size[0] * Mbar->size[1];
  Mbar->size[0] = colmat_bar->size[1];
  Mbar->size[1] = colmat_bar->size[1];
  emxEnsureCapacity_real_T(Mbar, i);
  Mbar_data = Mbar->data;
  loop_ub_tmp = colmat_bar->size[1] * colmat_bar->size[1];
  for (i = 0; i < loop_ub_tmp; i++) {
    Mbar_data[i] = 0.0;
  }
  i = Kbar11->size[0] * Kbar11->size[1];
  Kbar11->size[0] = colmat_bar->size[1];
  Kbar11->size[1] = colmat_bar->size[1];
  emxEnsureCapacity_real_T(Kbar11, i);
  Mbar_data = Kbar11->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    Mbar_data[i] = 0.0;
  }
  i = Dbar11->size[0] * Dbar11->size[1];
  Dbar11->size[0] = colmat_bar->size[1];
  Dbar11->size[1] = colmat_bar->size[1];
  emxEnsureCapacity_real_T(Dbar11, i);
  Mbar_data = Dbar11->data;
  for (i = 0; i < loop_ub_tmp; i++) {
    Mbar_data[i] = 0.0;
  }
  i = wg->size[0];
  emxInit_real_T(&Bbar, 1);
  emxInit_real_T(&r, 2);
  for (n = 0; n < i; n++) {
    double nxi0p;
    double scale;
    double t;
    double xi0p_idx_0;
    double xi0p_idx_1;
    double xi0p_idx_2;
    int aoffset;
    k = (int)(3.0 * (((double)n + 1.0) - 1.0) + 2.0);
    i1 = r->size[0] * r->size[1];
    r->size[0] = 1;
    loop_ub_tmp = colmat->size[1];
    r->size[1] = colmat->size[1];
    emxEnsureCapacity_real_T(r, i1);
    Mbar_data = r->data;
    for (i1 = 0; i1 < loop_ub_tmp; i1++) {
      Mbar_data[i1] = colmat_data[(k + colmat->size[0] * i1) - 1];
    }
    loop_ub_tmp = colmat_bar->size[1];
    k = Bbar->size[0];
    Bbar->size[0] = colmat_bar->size[1];
    emxEnsureCapacity_real_T(Bbar, k);
    Bbar_data = Bbar->data;
    for (k = 0; k < loop_ub_tmp; k++) {
      Bbar_data[k] = colmat_bar_data[n + colmat_bar->size[0] * k];
    }
    loop_ub_tmp = P0->size[1];
    xi0p_idx_0 = 0.0;
    xi0p_idx_1 = 0.0;
    xi0p_idx_2 = 0.0;
    for (k = 0; k < loop_ub_tmp; k++) {
      aoffset = k * 3;
      xi0p_idx_0 += P0_data[aoffset] * Mbar_data[k];
      xi0p_idx_1 += P0_data[aoffset + 1] * Mbar_data[k];
      xi0p_idx_2 += P0_data[aoffset + 2] * Mbar_data[k];
    }
    scale = 3.3121686421112381E-170;
    xi0p_idx_0 = fabs(xi0p_idx_0);
    if (xi0p_idx_0 > 3.3121686421112381E-170) {
      nxi0p = 1.0;
      scale = xi0p_idx_0;
    } else {
      t = xi0p_idx_0 / 3.3121686421112381E-170;
      nxi0p = t * t;
    }
    xi0p_idx_0 = fabs(xi0p_idx_1);
    if (xi0p_idx_0 > scale) {
      t = scale / xi0p_idx_0;
      nxi0p = nxi0p * t * t + 1.0;
      scale = xi0p_idx_0;
    } else {
      t = xi0p_idx_0 / scale;
      nxi0p += t * t;
    }
    xi0p_idx_0 = fabs(xi0p_idx_2);
    if (xi0p_idx_0 > scale) {
      t = scale / xi0p_idx_0;
      nxi0p = nxi0p * t * t + 1.0;
      scale = xi0p_idx_0;
    } else {
      t = xi0p_idx_0 / scale;
      nxi0p += t * t;
    }
    nxi0p = scale * sqrt(nxi0p);
    xi0p_idx_0 = nxi0p * wg_data[n];
    if ((Mbar->size[0] == Bbar->size[0]) && (Bbar->size[0] == Mbar->size[1])) {
      k = Mbar->size[0] * Mbar->size[1];
      Mbar->size[0] = Bbar->size[0];
      Mbar->size[1] = Bbar->size[0];
      emxEnsureCapacity_real_T(Mbar, k);
      Mbar_data = Mbar->data;
      loop_ub_tmp = Bbar->size[0];
      for (k = 0; k < loop_ub_tmp; k++) {
        aoffset = Bbar->size[0];
        for (i1 = 0; i1 < aoffset; i1++) {
          Mbar_data[i1 + Mbar->size[0] * k] +=
              Bbar_data[i1] * xi0p_idx_0 * Bbar_data[k];
        }
      }
    } else {
      binary_expand_op_30(Mbar, Bbar, xi0p_idx_0);
    }
    xi0p_idx_0 = EA * nxi0p * wg_data[n];
    if ((Kbar11->size[0] == Bbar->size[0]) &&
        (Bbar->size[0] == Kbar11->size[1])) {
      k = Kbar11->size[0] * Kbar11->size[1];
      Kbar11->size[0] = Bbar->size[0];
      Kbar11->size[1] = Bbar->size[0];
      emxEnsureCapacity_real_T(Kbar11, k);
      Mbar_data = Kbar11->data;
      loop_ub_tmp = Bbar->size[0];
      for (k = 0; k < loop_ub_tmp; k++) {
        aoffset = Bbar->size[0];
        for (i1 = 0; i1 < aoffset; i1++) {
          Mbar_data[i1 + Kbar11->size[0] * k] +=
              Bbar_data[i1] * xi0p_idx_0 * Bbar_data[k];
        }
      }
    } else {
      binary_expand_op_30(Kbar11, Bbar, xi0p_idx_0);
    }
    xi0p_idx_0 = Dm11 * nxi0p * wg_data[n];
    if ((Dbar11->size[0] == Bbar->size[0]) &&
        (Bbar->size[0] == Dbar11->size[1])) {
      k = Dbar11->size[0] * Dbar11->size[1];
      Dbar11->size[0] = Bbar->size[0];
      Dbar11->size[1] = Bbar->size[0];
      emxEnsureCapacity_real_T(Dbar11, k);
      Mbar_data = Dbar11->data;
      loop_ub_tmp = Bbar->size[0];
      for (k = 0; k < loop_ub_tmp; k++) {
        aoffset = Bbar->size[0];
        for (i1 = 0; i1 < aoffset; i1++) {
          Mbar_data[i1 + Dbar11->size[0] * k] +=
              Bbar_data[i1] * xi0p_idx_0 * Bbar_data[k];
        }
      }
    } else {
      binary_expand_op_30(Dbar11, Bbar, xi0p_idx_0);
    }
  }
  emxFree_real_T(&r);
  emxFree_real_T(&Bbar);
  emxInit_real_T(&r1, 2);
  mrdiv(Kbar11, Mbar, r1);
  mldivide(Mbar, r1, Kbar11);
  mrdiv(Dbar11, Mbar, r1);
  mldivide(Mbar, r1, Dbar11);
  emxFree_real_T(&r1);
}

/* End of code generation (CableMbar.c) */
