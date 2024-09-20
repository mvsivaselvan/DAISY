/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getBishopFrame.h
 *
 * Code generation for function 'getBishopFrame'
 *
 */

#ifndef GETBISHOPFRAME_H
#define GETBISHOPFRAME_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void getBishopFrame(const emxArray_real_T *P,
                           const emxArray_real_T *knots, double d,
                           const emxArray_real_T *xg, emxArray_real_T *R);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (getBishopFrame.h) */
