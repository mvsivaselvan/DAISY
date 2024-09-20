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
#include "spcolC.h"
#include "CableForceRotBCinCoord_types.h"
#include "rand.h"

/* Function Definitions */
void spcolC(const emxArray_real_T *knots, double k, emxArray_real_T *colmat)
{
  /*  This is a dummy function, since MATLAB coder cannot generate for spcol */
  /*  This function is called instead by the generator code.  */
  /*  It is forced to never inline, so that the this function is explicitly */
  /*  called in the nested function Bishop */
  /*  Code generated for this function will be deleted. The function spcolC in
   */
  /*  a Spline library will actually compute the necessary */
  /*  This is the C prototype: */
  /*  spcolC(double* knots, int nknots, int k, double* colpts, int ncolpts,  */
  /*         int nderiv, double* colmat, int ncolmat); */
  b_rand((double)knots->size[1] - k, colmat);
}

/* End of code generation (spcolC.c) */
