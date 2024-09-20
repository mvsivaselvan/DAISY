/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableMbar.h
 *
 * Code generation for function 'CableMbar'
 *
 */

#ifndef CABLEMBAR_H
#define CABLEMBAR_H

/* Include files */
#include "CableForceRotBCinCoord_types.h"
#include "rtwtypes.h"
#include <stddef.h>
#include <stdlib.h>

#ifdef __cplusplus
extern "C" {
#endif

/* Function Declarations */
extern void CableMbar(const emxArray_real_T *P0, double EA, double betAX,
                      const emxArray_real_T *colmat,
                      const emxArray_real_T *colmat_bar,
                      const emxArray_real_T *wg, emxArray_real_T *Mbar,
                      emxArray_real_T *Kbar11, emxArray_real_T *Dbar11);

#ifdef __cplusplus
}
#endif

#endif
/* End of code generation (CableMbar.h) */
