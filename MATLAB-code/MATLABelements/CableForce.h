/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableForce.h
 *
 * Code generation for function 'CableForce'
 *
 */

#ifndef CABLEFORCE_H
#define CABLEFORCE_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
void CableForce(const emxArray_real_T *P, const emxArray_real_T *P0,
                const emxArray_real_T *Pdot, const emxArray_real_T *varTheta,
                const emxArray_real_T *varThetadot, const emxArray_real_T *R0,
                double rho, double EA, double EI, double GJ, double betAX,
                double betBEND, double betTOR, const emxArray_real_T *wg,
                const emxArray_real_T *nel, const emxArray_real_T *colmat,
                const emxArray_real_T *colmat_brev,
                const emxArray_real_T *colmat_bar, double d, double dbrev,
                double dbar, const emxArray_real_T *Mbar, const double u[3],
                const emxArray_real_T *Kbar11, const emxArray_real_T *Dbar11,
                emxArray_real_T *F, emxArray_real_T *mu, emxArray_real_T *F_ij,
                emxArray_real_T *F_ib, emxArray_real_T *F_ijd,
                emxArray_real_T *F_ibd, emxArray_real_T *mu_aj,
                emxArray_real_T *mu_ab, emxArray_real_T *mu_ajd,
                emxArray_real_T *mu_abd);

void binary_expand_op_14(emxArray_real_T *in1, int in2, int in3, int in4,
                         int in5, const emxArray_real_T *in6);

void binary_expand_op_9(emxArray_real_T *in1, int in2, int in3, int in4,
                        int in5, int in6, int in7, int in8, int in9,
                        const emxArray_real_T *in10);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (CableForce.h) */
