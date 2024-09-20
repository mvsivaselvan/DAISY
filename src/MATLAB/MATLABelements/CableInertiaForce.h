/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableInertiaForce.h
 *
 * Code generation for function 'CableInertiaForce'
 *
 */

#ifndef CABLEINERTIAFORCE_H
#define CABLEINERTIAFORCE_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void CableInertiaForce(
    const emxArray_real_T *P, const emxArray_real_T *P0,
    const emxArray_real_T *Pdot, const emxArray_real_T *Pddot,
    const emxArray_real_T *varTheta, const emxArray_real_T *varThetadot,
    const emxArray_real_T *varThetaddot, const emxArray_real_T *R0, double rho,
    const double II[9], const emxArray_real_T *wg, const emxArray_real_T *nel,
    const emxArray_real_T *colmat, const emxArray_real_T *colmat_brev, double d,
    double dbrev, double alph0, emxArray_real_T *F, emxArray_real_T *mu,
    emxArray_real_T *F_ij, emxArray_real_T *F_ib, emxArray_real_T *F_ijd,
    emxArray_real_T *F_ibd, emxArray_real_T *F_ijdd, emxArray_real_T *F_ibdd,
    emxArray_real_T *mu_aj, emxArray_real_T *mu_ab, emxArray_real_T *mu_ajd,
    emxArray_real_T *mu_abd, emxArray_real_T *mu_ajdd,
    emxArray_real_T *mu_abdd);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (CableInertiaForce.h) */
