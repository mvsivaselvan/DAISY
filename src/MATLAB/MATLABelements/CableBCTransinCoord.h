/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableBCTransinCoord.h
 *
 * Code generation for function 'CableBCTransinCoord'
 *
 */

#ifndef CABLEBCTRANSINCOORD_H
#define CABLEBCTRANSINCOORD_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void CableBCTransinCoord(const double qbar[7], const double x0[3],
                         const double RJ[9], const double RE[9],
                         const double r[3], const double R0[9],
                         const double eta[7], const double rho1[7], double q[7],
                         double J[49], double Q[49], double Qtilde[49]);

void b_CableBCTransinCoord(const double qbar[7], const double x0[3],
                           const double RJ[9], const double RE[9],
                           const double r[3], const double R0[9],
                           const emxArray_real_T *eta, const double rho1[7],
                           double q[7], double J[49], double Q[49],
                           double Qtilde[49]);

void c_CableBCTransinCoord(const double qbar[7], const double x0[3],
                           const double RJ[9], const double RE[9],
                           const double r[3], const double R0[9],
                           const double eta[7], const double rho1[7],
                           const double rho2[7], double q[7], double J[49],
                           double Q[49], double Qtilde[49], double C[49]);

void d_CableBCTransinCoord(const double qbar[7], const double x0[3],
                           const double RJ[9], const double RE[9],
                           const double r[3], const double R0[9],
                           const emxArray_real_T *eta, const double rho1[7],
                           const double rho2[7], double q[7], double J[49],
                           double Q[49], double Qtilde[49], double C[49]);

void e_CableBCTransinCoord(const double qbar[7], const double x0[3],
                           const double RJ[9], const double RE[9],
                           const double r[3], const double R0[9], double q[7]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (CableBCTransinCoord.h) */
