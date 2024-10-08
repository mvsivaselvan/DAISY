/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * mldivide.h
 *
 * Code generation for function 'mldivide'
 *
 */

#ifndef MLDIVIDE_H
#define MLDIVIDE_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void b_mldivide(const emxArray_real_T *A, emxArray_real_T *B);

void mldivide(const emxArray_real_T *A, const emxArray_real_T *B,
              emxArray_real_T *Y);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (mldivide.h) */
