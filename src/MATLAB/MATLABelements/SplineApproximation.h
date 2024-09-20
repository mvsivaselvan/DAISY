/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * SplineApproximation.h
 *
 * Code generation for function 'SplineApproximation'
 *
 */

#ifndef SPLINEAPPROXIMATION_H
#define SPLINEAPPROXIMATION_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void SplineApproximation(const emxArray_real_T *gamm_,
                                const emxArray_real_T *J_, double N,
                                const emxArray_real_T *xg,
                                const emxArray_real_T *wg,
                                const emxArray_real_T *colmat,
                                emxArray_real_T *p, double *err);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (SplineApproximation.h) */
