/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableForceRotBCinCoord.h
 *
 * Code generation for function 'CableForceRotBCinCoord'
 *
 */

#ifndef CABLEFORCEROTBCINCOORD_H
#define CABLEFORCEROTBCINCOORD_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void CableForceRotBCinCoord(
    const double d1[3], const double phi1[3], double gamm1, const double d2[3],
    const double phi2[3], double gamm2, const emxArray_real_T *Pmid,
    const emxArray_real_T *varThetamid, const emxArray_real_T *P0,
    const double d1dot[3], const double phi1dot[3], double gamm1dot,
    const double d2dot[3], const double phi2dot[3], double gamm2dot,
    const emxArray_real_T *Pmiddot, const emxArray_real_T *varThetamiddot,
    const double d1ddot[3], const double phi1ddot[3], double gamm1ddot,
    const double d2ddot[3], const double phi2ddot[3], double gamm2ddot,
    const emxArray_real_T *Pmidddot, const emxArray_real_T *varThetamidddot,
    const double x01[3], const double RJ1[9], const double RE1[9],
    const double r1[3], const double x02[3], const double RJ2[9],
    const double RE2[9], const double r2[3], const emxArray_real_T *R0,
    const double II[9], double rho, double EA, double EI, double GJ,
    double betAX, double betBEND, double betTOR, const emxArray_real_T *sg,
    const emxArray_real_T *wg, const emxArray_real_T *nel,
    const emxArray_real_T *colmat, const emxArray_real_T *colmat_brev,
    const emxArray_real_T *colmat_bar, double d, double dbrev, double dbar,
    const emxArray_real_T *Mbar, const double u[3],
    const emxArray_real_T *Kbar11, const emxArray_real_T *Dbar11,
    double dynamic, double alph0, emxArray_real_T *Fb, emxArray_real_T *Kb,
    emxArray_real_T *Cb, emxArray_real_T *Mb, emxArray_real_T *Bb);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (CableForceRotBCinCoord.h) */
