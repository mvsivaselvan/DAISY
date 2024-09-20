/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RigidBodyForce.h
 *
 * Code generation for function 'RigidBodyForce'
 *
 */

#ifndef RIGIDBODYFORCE_H
#define RIGIDBODYFORCE_H

/* Include files */
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void RigidBodyForce(const double d[3], const double phi[3],
                           const double dd[3], const double phid[3],
                           const double ddd[3], const double phidd[3], double m,
                           const double II[9], const double KT[9],
                           const double CT[9], const double KR[9],
                           const double CR[9], const double rr[3],
                           const double RJ[9], const double u[3], double F[6],
                           double K[36], double C[36], double M[36],
                           double B[18]);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (RigidBodyForce.h) */
