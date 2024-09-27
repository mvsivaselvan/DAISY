/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getCableDeformedShape.h
 *
 * Code generation for function 'getCableDeformedShape'
 *
 */

#ifndef GETCABLEDEFORMEDSHAPE_H
#define GETCABLEDEFORMEDSHAPE_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void getCableDeformedShape(
    const double d1[3], const double phi1[3], double gamm1, const double d2[3],
    const double phi2[3], double gamm2, const emxArray_real_T *Pmid,
    const double x01[3], const double RJ1[9], const double RE1[9],
    const double r1[3], const double Rb10[9], const double x02[3],
    const double RJ2[9], const double RE2[9], const double r2[3],
    const double Rb20[9], emxArray_real_T *P);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (getCableDeformedShape.h) */
