/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableForceInputDerivative.h
 *
 * Code generation for function 'CableForceInputDerivative'
 *
 */

#ifndef CABLEFORCEINPUTDERIVATIVE_H
#define CABLEFORCEINPUTDERIVATIVE_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void CableForceInputDerivative(const emxArray_real_T *P0, double rho,
                               const emxArray_real_T *wg,
                               const emxArray_real_T *nel,
                               const emxArray_real_T *colmat, double d,
                               emxArray_real_T *dFdu);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (CableForceInputDerivative.h) */
