/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * blkdiag.c
 *
 * Code generation for function 'blkdiag'
 *
 */

/* Include files */
#include "blkdiag.h"
#include <string.h>

/* Function Definitions */
void blkdiag(const double varargin_1[49], const double varargin_2[49],
             double y[196])
{
  int i;
  int i1;
  memset(&y[0], 0, 196U * sizeof(double));
  for (i = 0; i < 7; i++) {
    for (i1 = 0; i1 < 7; i1++) {
      int y_tmp;
      y_tmp = i1 + 7 * i;
      y[i1 + 14 * i] = varargin_1[y_tmp];
      y[(i1 + 14 * (i + 7)) + 7] = varargin_2[y_tmp];
    }
  }
}

/* End of code generation (blkdiag.c) */
