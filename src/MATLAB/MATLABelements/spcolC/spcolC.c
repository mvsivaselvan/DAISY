/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * spcolC.c
 *
 * Code generation for function 'spcolC'
 *
 */

/* Include files */
#include "../spcolC.h"
#include "../CableForceRotBCinCoord_emxutil.h"
#include "../CableForceRotBCinCoord_types.h"

extern void mySspcol(double *knots, int nknots, int k, double *colpts,
                     int ncolpts, int nderiv, double *colmat, int ncolmat);

/* Function Definitions */
void spcolC(const emxArray_real_T *knots, double k, double colpt,
            emxArray_real_T *colmat)
{
  int nknots = knots->size[0];
  int nbasis = nknots - (int)k;
  colmat->data = (double *)malloc(3 * nbasis * sizeof(double));
  colmat->size[0] = 3;
  colmat->size[1] = nbasis;
  mySpcol(knots->data, nknots, (int)k, &colpt, 1, 3, colmat->data, 3*nbasis);
}

/* End of code generation (spcolC.c) */
