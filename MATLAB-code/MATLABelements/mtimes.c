/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mtimes.c
 *
 * Code generation for function 'mtimes'
 *
 */

/* Include files */
#include "mtimes.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_types.h"

/* Function Definitions */
void b_mtimes(const double A[9], const double B[9], double C[9])
{
  int i;
  int j;
  for (j = 0; j < 3; j++) {
    double d;
    double d1;
    double d2;
    int coffset;
    coffset = j * 3;
    d = B[j];
    d1 = B[j + 3];
    d2 = B[j + 6];
    for (i = 0; i < 3; i++) {
      int aoffset;
      aoffset = i * 3;
      C[coffset + i] =
          (A[aoffset] * d + A[aoffset + 1] * d1) + A[aoffset + 2] * d2;
    }
  }
}

void c_mtimes(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *C)
{
  const double *A_data;
  const double *B_data;
  double *C_data;
  int b_i;
  int i;
  int inner;
  int j;
  int k;
  int mc;
  int nc;
  B_data = B->data;
  A_data = A->data;
  mc = A->size[1];
  inner = A->size[0];
  nc = B->size[1];
  i = C->size[0] * C->size[1];
  C->size[0] = A->size[1];
  C->size[1] = B->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  for (j = 0; j < nc; j++) {
    int boffset;
    int coffset;
    coffset = j * mc;
    boffset = j * B->size[0];
    for (b_i = 0; b_i < mc; b_i++) {
      C_data[coffset + b_i] = 0.0;
    }
    for (k = 0; k < inner; k++) {
      double bkj;
      bkj = B_data[boffset + k];
      for (b_i = 0; b_i < mc; b_i++) {
        i = coffset + b_i;
        C_data[i] += A_data[b_i * A->size[0] + k] * bkj;
      }
    }
  }
}

void d_mtimes(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *C)
{
  const double *A_data;
  const double *B_data;
  double *C_data;
  int b_i;
  int i;
  int inner;
  int j;
  int k;
  int mc;
  int nc;
  B_data = B->data;
  A_data = A->data;
  mc = A->size[0];
  inner = A->size[1];
  nc = B->size[1];
  i = C->size[0] * C->size[1];
  C->size[0] = A->size[0];
  C->size[1] = B->size[1];
  emxEnsureCapacity_real_T(C, i);
  C_data = C->data;
  for (j = 0; j < nc; j++) {
    int boffset;
    int coffset;
    coffset = j * mc;
    boffset = j * B->size[0];
    for (b_i = 0; b_i < mc; b_i++) {
      C_data[coffset + b_i] = 0.0;
    }
    for (k = 0; k < inner; k++) {
      double bkj;
      int aoffset;
      aoffset = k * A->size[0];
      bkj = B_data[boffset + k];
      for (b_i = 0; b_i < mc; b_i++) {
        i = coffset + b_i;
        C_data[i] += A_data[aoffset + b_i] * bkj;
      }
    }
  }
}

void e_mtimes(const emxArray_real_T *A, const double B[196], emxArray_real_T *C)
{
  const double *A_data;
  double *C_data;
  int coffset;
  int i;
  int j;
  int k;
  int m;
  A_data = A->data;
  m = A->size[0];
  coffset = C->size[0] * C->size[1];
  C->size[0] = A->size[0];
  C->size[1] = 14;
  emxEnsureCapacity_real_T(C, coffset);
  C_data = C->data;
  for (j = 0; j < 14; j++) {
    int boffset;
    coffset = j * m;
    boffset = j * 14;
    for (i = 0; i < m; i++) {
      double s;
      s = 0.0;
      for (k = 0; k < 14; k++) {
        s += A_data[k * A->size[0] + i] * B[boffset + k];
      }
      C_data[coffset + i] = s;
    }
  }
}

void f_mtimes(const double A[196], const emxArray_real_T *B, emxArray_real_T *C)
{
  const double *B_data;
  double *C_data;
  int coffset_tmp;
  int i;
  int j;
  int k;
  int n;
  B_data = B->data;
  n = B->size[1];
  coffset_tmp = C->size[0] * C->size[1];
  C->size[0] = 14;
  C->size[1] = B->size[1];
  emxEnsureCapacity_real_T(C, coffset_tmp);
  C_data = C->data;
  for (j = 0; j < n; j++) {
    coffset_tmp = j * 14;
    for (i = 0; i < 14; i++) {
      double s;
      int aoffset;
      aoffset = i * 14;
      s = 0.0;
      for (k = 0; k < 14; k++) {
        s += A[aoffset + k] * B_data[coffset_tmp + k];
      }
      C_data[coffset_tmp + i] = s;
    }
  }
}

void mtimes(const emxArray_real_T *A, const emxArray_real_T *B, double C[3])
{
  const double *A_data;
  const double *B_data;
  int inner;
  int k;
  B_data = B->data;
  A_data = A->data;
  inner = A->size[1];
  C[0] = 0.0;
  C[1] = 0.0;
  C[2] = 0.0;
  for (k = 0; k < inner; k++) {
    int aoffset;
    aoffset = k * 3;
    C[0] += A_data[aoffset] * B_data[k];
    C[1] += A_data[aoffset + 1] * B_data[k];
    C[2] += A_data[aoffset + 2] * B_data[k];
  }
}

/* End of code generation (mtimes.c) */
