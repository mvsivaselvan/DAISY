/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mldivide.c
 *
 * Code generation for function 'mldivide'
 *
 */

/* Include files */
#include "mldivide.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_types.h"
#include "qrsolve.h"
#include "xgeqp3.h"
#include "xgetrf.h"

/* Function Definitions */
void b_mldivide(const emxArray_real_T *A, emxArray_real_T *B)
{
  emxArray_int32_T *jpvt;
  emxArray_real_T *b_A;
  emxArray_real_T *b_B;
  emxArray_real_T *tau;
  const double *A_data;
  double *B_data;
  double *b_A_data;
  double *b_B_data;
  double *tau_data;
  int LDA;
  int b_i;
  int i;
  int rankA;
  int *jpvt_data;
  B_data = B->data;
  A_data = A->data;
  emxInit_real_T(&b_A, 2);
  emxInit_real_T(&tau, 1);
  emxInit_int32_T(&jpvt, 2);
  emxInit_real_T(&b_B, 1);
  if ((A->size[0] == 0) || (A->size[1] == 0) || (B->size[0] == 0)) {
    i = B->size[0];
    B->size[0] = A->size[1];
    emxEnsureCapacity_real_T(B, i);
    B_data = B->data;
    LDA = A->size[1];
    for (i = 0; i < LDA; i++) {
      B_data[i] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    double wj;
    int m;
    int mn;
    LDA = A->size[0];
    mn = A->size[1];
    if (LDA <= mn) {
      mn = LDA;
    }
    LDA = B->size[0];
    if (LDA <= mn) {
      mn = LDA;
    }
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    b_A_data = b_A->data;
    LDA = A->size[0] * A->size[1];
    for (i = 0; i < LDA; i++) {
      b_A_data[i] = A_data[i];
    }
    xgetrf(mn, mn, b_A, A->size[0], jpvt);
    jpvt_data = jpvt->data;
    b_A_data = b_A->data;
    LDA = b_A->size[0];
    for (b_i = 0; b_i <= mn - 2; b_i++) {
      i = jpvt_data[b_i];
      if (i != b_i + 1) {
        wj = B_data[b_i];
        B_data[b_i] = B_data[i - 1];
        B_data[i - 1] = wj;
      }
    }
    for (rankA = 0; rankA < mn; rankA++) {
      m = LDA * rankA;
      if (B_data[rankA] != 0.0) {
        i = rankA + 2;
        for (b_i = i; b_i <= mn; b_i++) {
          B_data[b_i - 1] -= B_data[rankA] * b_A_data[(b_i + m) - 1];
        }
      }
    }
    for (rankA = mn; rankA >= 1; rankA--) {
      m = LDA * (rankA - 1);
      wj = B_data[rankA - 1];
      if (wj != 0.0) {
        wj /= b_A_data[(rankA + m) - 1];
        B_data[rankA - 1] = wj;
        for (b_i = 0; b_i <= rankA - 2; b_i++) {
          B_data[b_i] -= B_data[rankA - 1] * b_A_data[b_i + m];
        }
      }
    }
  } else {
    int m;
    int mn;
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    b_A_data = b_A->data;
    LDA = A->size[0] * A->size[1];
    for (i = 0; i < LDA; i++) {
      b_A_data[i] = A_data[i];
    }
    xgeqp3(b_A, tau, jpvt);
    jpvt_data = jpvt->data;
    tau_data = tau->data;
    b_A_data = b_A->data;
    rankA = rankFromQR(b_A);
    i = b_B->size[0];
    b_B->size[0] = B->size[0];
    emxEnsureCapacity_real_T(b_B, i);
    b_B_data = b_B->data;
    LDA = B->size[0];
    for (i = 0; i < LDA; i++) {
      b_B_data[i] = B_data[i];
    }
    i = B->size[0];
    B->size[0] = b_A->size[1];
    emxEnsureCapacity_real_T(B, i);
    B_data = B->data;
    LDA = b_A->size[1];
    for (i = 0; i < LDA; i++) {
      B_data[i] = 0.0;
    }
    m = b_A->size[0];
    LDA = b_A->size[0];
    mn = b_A->size[1];
    if (LDA <= mn) {
      mn = LDA;
    }
    for (LDA = 0; LDA < mn; LDA++) {
      if (tau_data[LDA] != 0.0) {
        double wj;
        wj = b_B_data[LDA];
        i = LDA + 2;
        for (b_i = i; b_i <= m; b_i++) {
          wj += b_A_data[(b_i + b_A->size[0] * LDA) - 1] * b_B_data[b_i - 1];
        }
        wj *= tau_data[LDA];
        if (wj != 0.0) {
          b_B_data[LDA] -= wj;
          for (b_i = i; b_i <= m; b_i++) {
            b_B_data[b_i - 1] -= b_A_data[(b_i + b_A->size[0] * LDA) - 1] * wj;
          }
        }
      }
    }
    for (b_i = 0; b_i < rankA; b_i++) {
      B_data[jpvt_data[b_i] - 1] = b_B_data[b_i];
    }
    for (LDA = rankA; LDA >= 1; LDA--) {
      i = jpvt_data[LDA - 1];
      B_data[i - 1] /= b_A_data[(LDA + b_A->size[0] * (LDA - 1)) - 1];
      for (b_i = 0; b_i <= LDA - 2; b_i++) {
        B_data[jpvt_data[b_i] - 1] -=
            B_data[i - 1] * b_A_data[b_i + b_A->size[0] * (LDA - 1)];
      }
    }
  }
  emxFree_real_T(&b_B);
  emxFree_int32_T(&jpvt);
  emxFree_real_T(&tau);
  emxFree_real_T(&b_A);
}

void mldivide(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *Y)
{
  emxArray_int32_T *ipiv;
  emxArray_real_T *b_A;
  const double *A_data;
  const double *B_data;
  double *Y_data;
  double *b_A_data;
  int b_i;
  int i;
  int j;
  int k;
  int *ipiv_data;
  B_data = B->data;
  A_data = A->data;
  emxInit_real_T(&b_A, 2);
  emxInit_int32_T(&ipiv, 2);
  if ((A->size[0] == 0) || (A->size[1] == 0) ||
      ((B->size[0] == 0) || (B->size[1] == 0))) {
    int LDA;
    i = Y->size[0] * Y->size[1];
    Y->size[0] = A->size[1];
    Y->size[1] = B->size[1];
    emxEnsureCapacity_real_T(Y, i);
    Y_data = Y->data;
    LDA = A->size[1] * B->size[1];
    for (i = 0; i < LDA; i++) {
      Y_data[i] = 0.0;
    }
  } else if (A->size[0] == A->size[1]) {
    double temp;
    int LDA;
    int LDB;
    int i1;
    int jBcol;
    int kAcol;
    int n;
    int nrhs;
    i = Y->size[0] * Y->size[1];
    Y->size[0] = B->size[0];
    Y->size[1] = B->size[1];
    emxEnsureCapacity_real_T(Y, i);
    Y_data = Y->data;
    LDA = B->size[0] * B->size[1];
    for (i = 0; i < LDA; i++) {
      Y_data[i] = B_data[i];
    }
    LDA = A->size[0];
    n = A->size[1];
    if (LDA <= n) {
      n = LDA;
    }
    LDA = B->size[0];
    if (LDA <= n) {
      n = LDA;
    }
    nrhs = B->size[1] - 1;
    i = b_A->size[0] * b_A->size[1];
    b_A->size[0] = A->size[0];
    b_A->size[1] = A->size[1];
    emxEnsureCapacity_real_T(b_A, i);
    b_A_data = b_A->data;
    LDA = A->size[0] * A->size[1];
    for (i = 0; i < LDA; i++) {
      b_A_data[i] = A_data[i];
    }
    xgetrf(n, n, b_A, A->size[0], ipiv);
    ipiv_data = ipiv->data;
    b_A_data = b_A->data;
    LDA = b_A->size[0];
    LDB = B->size[0];
    for (b_i = 0; b_i <= n - 2; b_i++) {
      i = ipiv_data[b_i];
      if (i != b_i + 1) {
        for (j = 0; j <= nrhs; j++) {
          temp = Y_data[b_i + Y->size[0] * j];
          Y_data[b_i + Y->size[0] * j] = Y_data[(i + Y->size[0] * j) - 1];
          Y_data[(i + Y->size[0] * j) - 1] = temp;
        }
      }
    }
    for (j = 0; j <= nrhs; j++) {
      jBcol = LDB * j;
      for (k = 0; k < n; k++) {
        kAcol = LDA * k;
        i = k + jBcol;
        if (Y_data[i] != 0.0) {
          i1 = k + 2;
          for (b_i = i1; b_i <= n; b_i++) {
            int i2;
            i2 = (b_i + jBcol) - 1;
            Y_data[i2] -= Y_data[i] * b_A_data[(b_i + kAcol) - 1];
          }
        }
      }
    }
    for (j = 0; j <= nrhs; j++) {
      jBcol = LDB * j - 1;
      for (k = n; k >= 1; k--) {
        kAcol = LDA * (k - 1) - 1;
        i = k + jBcol;
        temp = Y_data[i];
        if (temp != 0.0) {
          Y_data[i] = temp / b_A_data[k + kAcol];
          for (b_i = 0; b_i <= k - 2; b_i++) {
            i1 = (b_i + jBcol) + 1;
            Y_data[i1] -= Y_data[i] * b_A_data[(b_i + kAcol) + 1];
          }
        }
      }
    }
  } else {
    qrsolve(A, B, Y);
  }
  emxFree_int32_T(&ipiv);
  emxFree_real_T(&b_A);
}

/* End of code generation (mldivide.c) */
