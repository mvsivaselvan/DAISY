/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * getBishopFrame.c
 *
 * Code generation for function 'getBishopFrame'
 *
 */

/* Include files */
#include "getBishopFrame.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_emxutil.h"
#include "CableForceRotBCinCoord_initialize.h"
#include "CableForceRotBCinCoord_internal_types.h"
#include "CableForceRotBCinCoord_types.h"
#include "spcolC.h"
#include <math.h>
#include <string.h>

/* Function Declarations */
static void Bishop(const d_captured_var *knots, const b_captured_var *d,
                   const d_captured_var *P, const c_captured_var *I3,
                   const double q[4], double qdot[4]);

/* Function Definitions */
static void Bishop(const d_captured_var *knots, const b_captured_var *d,
                   const d_captured_var *P, const c_captured_var *I3,
                   const double q[4], double qdot[4])
{
  emxArray_real_T *colmat;
  double dv[16];
  double b_I3[9];
  double taup[3];
  double xip[3];
  double nxip;
  double xipp_idx_0;
  double xipp_idx_1;
  double xipp_idx_2;
  double *colmat_data;
  int aoffset;
  int inner;
  int k;
  int xip_tmp;
  /*  generating code */
  emxInit_real_T(&colmat, 2);
  spcolC(knots->contents, d->contents + 1.0, colmat);
  colmat_data = colmat->data;
  /*  B = colmat(1,:)'; % not used */
  inner = P->contents->size[1];
  xip[0] = 0.0;
  xip[1] = 0.0;
  xip[2] = 0.0;
  for (k = 0; k < inner; k++) {
    aoffset = k * 3;
    xip_tmp = 3 * k + 1;
    xip[0] += P->contents->data[aoffset] * colmat_data[xip_tmp];
    xip[1] += P->contents->data[aoffset + 1] * colmat_data[xip_tmp];
    xip[2] += P->contents->data[aoffset + 2] * colmat_data[xip_tmp];
  }
  inner = P->contents->size[1];
  xipp_idx_0 = 0.0;
  xipp_idx_1 = 0.0;
  xipp_idx_2 = 0.0;
  for (k = 0; k < inner; k++) {
    aoffset = k * 3;
    xip_tmp = 3 * k + 2;
    xipp_idx_0 += P->contents->data[aoffset] * colmat_data[xip_tmp];
    xipp_idx_1 += P->contents->data[aoffset + 1] * colmat_data[xip_tmp];
    xipp_idx_2 += P->contents->data[aoffset + 2] * colmat_data[xip_tmp];
  }
  emxFree_real_T(&colmat);
  nxip = sqrt((xip[0] * xip[0] + xip[1] * xip[1]) + xip[2] * xip[2]);
  /*  if (nxip == 0) % This situation will never heppen, it is to avoid */
  /*                 % nan and inf checks in code generation */
  /*      qdot = [0; 0; 0; 0]; */
  /*      return */
  /*  end */
  xip[0] /= nxip;
  xip[1] /= nxip;
  xip[2] /= nxip;
  for (inner = 0; inner < 3; inner++) {
    b_I3[3 * inner] = I3->contents[3 * inner] - xip[0] * xip[inner];
    xip_tmp = 3 * inner + 1;
    b_I3[xip_tmp] = I3->contents[xip_tmp] - xip[1] * xip[inner];
    xip_tmp = 3 * inner + 2;
    b_I3[xip_tmp] = I3->contents[xip_tmp] - xip[2] * xip[inner];
  }
  for (inner = 0; inner < 3; inner++) {
    taup[inner] = ((b_I3[inner] * xipp_idx_0 + b_I3[inner + 3] * xipp_idx_1) +
                   b_I3[inner + 6] * xipp_idx_2) /
                  nxip;
  }
  double b_d;
  double d1;
  double d2;
  xipp_idx_0 = xip[1] * taup[2] - taup[1] * xip[2];
  xipp_idx_1 = taup[0] * xip[2] - xip[0] * taup[2];
  xipp_idx_2 = xip[0] * taup[1] - taup[0] * xip[1];
  /*  feedback gain to prevent drift of q from unit norm */
  nxip = ((q[0] * q[0] + q[1] * q[1]) + q[2] * q[2]) + q[3] * q[3];
  nxip = 0.05 / nxip * (1.0 - nxip);
  dv[0] = 0.0;
  b_d = 0.5 * -xipp_idx_0;
  dv[4] = b_d;
  dv[1] = 0.5 * xipp_idx_0;
  d1 = 0.5 * -xipp_idx_1;
  dv[8] = d1;
  dv[2] = 0.5 * xipp_idx_1;
  d2 = 0.5 * -xipp_idx_2;
  dv[12] = d2;
  dv[3] = 0.5 * xipp_idx_2;
  dv[5] = 0.0;
  dv[9] = d2;
  dv[13] = 0.5 * xipp_idx_1;
  dv[6] = 0.5 * xipp_idx_2;
  dv[10] = 0.0;
  dv[14] = b_d;
  dv[7] = d1;
  dv[11] = 0.5 * xipp_idx_0;
  dv[15] = 0.0;
  for (inner = 0; inner < 4; inner++) {
    qdot[inner] =
        (((dv[inner] * q[0] + dv[inner + 4] * q[1]) + dv[inner + 8] * q[2]) +
         dv[inner + 12] * q[3]) +
        nxip * q[inner];
  }
}

void getBishopFrame(const emxArray_real_T *P, const emxArray_real_T *knots,
                    double d, const emxArray_real_T *xg, emxArray_real_T *R)
{
  static const double x[21] = {0.2,
                               0.075,
                               0.225,
                               0.97777777777777775,
                               -3.7333333333333334,
                               3.5555555555555554,
                               2.9525986892242035,
                               -11.595793324188385,
                               9.8228928516994358,
                               -0.29080932784636487,
                               2.8462752525252526,
                               -10.757575757575758,
                               8.9064227177434727,
                               0.27840909090909088,
                               -0.2735313036020583,
                               0.091145833333333329,
                               0.0,
                               0.44923629829290207,
                               0.65104166666666663,
                               -0.322376179245283,
                               0.13095238095238096};
  static const double b[7] = {0.0012326388888888888,
                              0.0,
                              -0.0042527702905061394,
                              0.036979166666666667,
                              -0.05086379716981132,
                              0.0419047619047619,
                              -0.025};
  static const double b_b[7] = {-2.859375,
                                0.0,
                                4.0431266846361185,
                                -3.90625,
                                2.7939268867924527,
                                -1.5714285714285714,
                                1.5};
  static const double c_b[7] = {3.0833333333333335,
                                0.0,
                                -6.2893081761006293,
                                10.416666666666666,
                                -6.8773584905660377,
                                3.6666666666666665,
                                -4.0};
  static const double d_b[7] = {-1.1328125,
                                0.0,
                                2.6954177897574123,
                                -5.859375,
                                3.7610554245283021,
                                -1.9642857142857142,
                                2.5};
  static const signed char b_iv[9] = {1, 0, 0, 0, 1, 0, 0, 0, 0};
  b_captured_var b_d;
  c_captured_var I3;
  d_captured_var b_P;
  d_captured_var b_knots;
  emxArray_real_T *Q;
  emxArray_real_T *yout;
  double f[28];
  double f0[4];
  double q0[4];
  const double *P_data;
  const double *knots_data;
  const double *xg_data;
  double absh;
  double c_d;
  double d1;
  double d2;
  double d3;
  double h;
  double hmax;
  double hmin;
  double normtau0;
  double t;
  double tau0_idx_0;
  double tau0_idx_1;
  double tau0_idx_2;
  double tfinal_tmp;
  double tnew;
  double *Q_data;
  double *yout_data;
  int Bcolidx;
  int exponent;
  int i;
  int i1;
  int ia;
  int iac;
  int j;
  int next;
  int nnxt;
  int nout;
  int tdir;
  boolean_T Done;
  boolean_T MinStepExit;
  if (!isInitialized_CableForceRotBCinCoord) {
    CableForceRotBCinCoord_initialize();
  }
  xg_data = xg->data;
  knots_data = knots->data;
  P_data = P->data;
  emxInitStruct_captured_var(&b_P);
  i = b_P.contents->size[0] * b_P.contents->size[1];
  b_P.contents->size[0] = 3;
  b_P.contents->size[1] = P->size[1];
  emxEnsureCapacity_real_T(b_P.contents, i);
  nnxt = 3 * P->size[1];
  for (i = 0; i < nnxt; i++) {
    b_P.contents->data[i] = P_data[i];
  }
  emxInitStruct_captured_var(&b_knots);
  i = b_knots.contents->size[0] * b_knots.contents->size[1];
  b_knots.contents->size[0] = 1;
  b_knots.contents->size[1] = knots->size[1];
  emxEnsureCapacity_real_T(b_knots.contents, i);
  nnxt = knots->size[1];
  for (i = 0; i < nnxt; i++) {
    b_knots.contents->data[i] = knots_data[i];
  }
  b_d.contents = d;
  /*  Computes the (twist-free) Bishop framing of a given curve */
  /*  INPUTS */
  /*  P = control points of curve arranged as 3XN matrix */
  /*  knots = knot vector for B-spline basis functions */
  /*  d = degree of B-spline basis functions used */
  /*  xg = collocation points (in arclength parametrization) where Bishop */
  /*       frames are to be computed */
  /*  q0 = initial condition for integration (quaternion representation  */
  /*       of R at s=0) */
  /*  OUTPUT */
  /*  R = cell array with length(xg) rotation matrices defining the Bishop */
  /*      frame */
  c_d = P_data[3] - P_data[0];
  tau0_idx_0 = c_d;
  normtau0 = c_d * c_d;
  c_d = P_data[4] - P_data[1];
  tau0_idx_1 = c_d;
  normtau0 += c_d * c_d;
  c_d = P_data[5] - P_data[2];
  normtau0 += c_d * c_d;
  normtau0 = sqrt(normtau0);
  tau0_idx_0 /= normtau0;
  tau0_idx_1 /= normtau0;
  tau0_idx_2 = c_d / normtau0;
  normtau0 = acos(tau0_idx_0);
  /*  tau0(1) = dot product of tau0 with [1;0;0] */
  if (fabs(normtau0) < 1.0E-15) {
    /*  tau0 is oriented along [1;0;0] */
    q0[0] = 1.0;
    q0[1] = 0.0;
    q0[2] = 0.0;
    q0[3] = 0.0;
  } else {
    tau0_idx_0 = sqrt(tau0_idx_1 * tau0_idx_1 + tau0_idx_2 * tau0_idx_2);
    /*  i.e. unit vector along cross([1;0;0],tau0) */
    hmax = sin(normtau0 / 2.0);
    q0[0] = cos(normtau0 / 2.0);
    q0[1] = hmax * (0.0 / tau0_idx_0);
    q0[2] = hmax * (-tau0_idx_2 / tau0_idx_0);
    q0[3] = hmax * (tau0_idx_1 / tau0_idx_0);
  }
  for (i = 0; i < 9; i++) {
    I3.contents[i] = b_iv[i];
  }
  tfinal_tmp = xg_data[xg->size[0] - 1];
  Bishop(&b_knots, &b_d, &b_P, &I3, q0, f0);
  emxInit_real_T(&yout, 2);
  i = yout->size[0] * yout->size[1];
  yout->size[0] = 4;
  yout->size[1] = xg->size[0];
  emxEnsureCapacity_real_T(yout, i);
  yout_data = yout->data;
  nnxt = xg->size[0] << 2;
  for (i = 0; i < nnxt; i++) {
    yout_data[i] = 0.0;
  }
  nout = 1;
  yout_data[0] = q0[0];
  yout_data[1] = q0[1];
  yout_data[2] = q0[2];
  yout_data[3] = q0[3];
  tau0_idx_1 = tfinal_tmp - xg_data[0];
  normtau0 = fabs(tau0_idx_1);
  tau0_idx_0 = fabs(xg_data[0]);
  hmax = fmin(normtau0,
              fmax(0.1 * normtau0, 3.5527136788005009E-15 *
                                       fmax(tau0_idx_0, fabs(tfinal_tmp))));
  if (tau0_idx_0 < 4.4501477170144028E-308) {
    tau0_idx_2 = 4.94065645841247E-324;
  } else {
    frexp(tau0_idx_0, &Bcolidx);
    tau0_idx_2 = ldexp(1.0, Bcolidx - 53);
  }
  absh = fmin(hmax, fabs(xg_data[1] - xg_data[0]));
  normtau0 = 0.0;
  tau0_idx_0 = fabs(f0[0] / fmax(fabs(q0[0]), 0.004));
  if (tau0_idx_0 > 0.0) {
    normtau0 = tau0_idx_0;
  }
  tau0_idx_0 = fabs(f0[1] / 0.004);
  if (tau0_idx_0 > normtau0) {
    normtau0 = tau0_idx_0;
  }
  tau0_idx_0 = fabs(f0[2] / fmax(fabs(q0[2]), 0.004));
  if (tau0_idx_0 > normtau0) {
    normtau0 = tau0_idx_0;
  }
  tau0_idx_0 = fabs(f0[3] / fmax(fabs(q0[3]), 0.004));
  if (tau0_idx_0 > normtau0) {
    normtau0 = tau0_idx_0;
  }
  normtau0 /= 0.0015229231509727023;
  if (absh * normtau0 > 1.0) {
    absh = 1.0 / normtau0;
  }
  absh = fmax(absh, 16.0 * tau0_idx_2);
  t = xg_data[0];
  memset(&f[0], 0, 28U * sizeof(double));
  f[0] = f0[0];
  f[1] = f0[1];
  f[2] = f0[2];
  f[3] = f0[3];
  if (tau0_idx_1 < 0.0) {
    tdir = -1;
  } else {
    tdir = (tau0_idx_1 > 0.0);
  }
  next = 0;
  MinStepExit = false;
  Done = false;
  int exitg1;
  do {
    double ynew[4];
    boolean_T NoFailedAttempts;
    exitg1 = 0;
    tau0_idx_0 = fabs(t);
    if (tau0_idx_0 < 4.4501477170144028E-308) {
      tau0_idx_2 = 4.94065645841247E-324;
    } else {
      frexp(tau0_idx_0, &exponent);
      tau0_idx_2 = ldexp(1.0, exponent - 53);
    }
    hmin = 16.0 * tau0_idx_2;
    absh = fmin(hmax, fmax(hmin, absh));
    h = (double)tdir * absh;
    c_d = tfinal_tmp - t;
    d1 = fabs(c_d);
    if (1.1 * absh >= d1) {
      h = c_d;
      absh = d1;
      Done = true;
    }
    NoFailedAttempts = true;
    int exitg2;
    do {
      exitg2 = 0;
      Bcolidx = 0;
      for (j = 0; j < 5; j++) {
        Bcolidx += j;
        f0[0] = q0[0];
        f0[1] = q0[1];
        f0[2] = q0[2];
        f0[3] = q0[3];
        if (h != 0.0) {
          i = (j << 2) + 1;
          for (iac = 1; iac <= i; iac += 4) {
            normtau0 = h * x[Bcolidx + ((iac - 1) >> 2)];
            i1 = iac + 3;
            for (ia = iac; ia <= i1; ia++) {
              nnxt = ia - iac;
              f0[nnxt] += f[ia - 1] * normtau0;
            }
          }
        }
        Bishop(&b_knots, &b_d, &b_P, &I3, f0,
               *(double(*)[4]) & f[(j + 1) << 2]);
      }
      tnew = t + h;
      ynew[0] = q0[0];
      ynew[1] = q0[1];
      ynew[2] = q0[2];
      ynew[3] = q0[3];
      if (h != 0.0) {
        for (iac = 0; iac <= 20; iac += 4) {
          normtau0 = h * x[(Bcolidx + (iac >> 2)) + 5];
          i = iac + 4;
          for (ia = iac + 1; ia <= i; ia++) {
            nnxt = (ia - iac) - 1;
            ynew[nnxt] += f[ia - 1] * normtau0;
          }
        }
      }
      Bishop(&b_knots, &b_d, &b_P, &I3, ynew, *(double(*)[4]) & f[24]);
      for (i = 0; i < 4; i++) {
        c_d = 0.0;
        for (i1 = 0; i1 < 7; i1++) {
          c_d += f[i + (i1 << 2)] * b[i1];
        }
        f0[i] = c_d;
      }
      if (Done) {
        tnew = tfinal_tmp;
      }
      h = tnew - t;
      tau0_idx_1 = 0.0;
      normtau0 = fabs(f0[0]);
      tau0_idx_0 = fabs(q0[0]);
      tau0_idx_2 = fabs(ynew[0]);
      if (tau0_idx_0 > tau0_idx_2) {
        if (tau0_idx_0 > 0.004) {
          normtau0 /= tau0_idx_0;
        } else {
          normtau0 /= 0.004;
        }
      } else if (tau0_idx_2 > 0.004) {
        normtau0 /= tau0_idx_2;
      } else {
        normtau0 /= 0.004;
      }
      if (normtau0 > 0.0) {
        tau0_idx_1 = normtau0;
      }
      normtau0 = fabs(f0[1]);
      tau0_idx_0 = fabs(q0[1]);
      tau0_idx_2 = fabs(ynew[1]);
      if (tau0_idx_0 > tau0_idx_2) {
        if (tau0_idx_0 > 0.004) {
          normtau0 /= tau0_idx_0;
        } else {
          normtau0 /= 0.004;
        }
      } else if (tau0_idx_2 > 0.004) {
        normtau0 /= tau0_idx_2;
      } else {
        normtau0 /= 0.004;
      }
      if (normtau0 > tau0_idx_1) {
        tau0_idx_1 = normtau0;
      }
      normtau0 = fabs(f0[2]);
      tau0_idx_0 = fabs(q0[2]);
      tau0_idx_2 = fabs(ynew[2]);
      if (tau0_idx_0 > tau0_idx_2) {
        if (tau0_idx_0 > 0.004) {
          normtau0 /= tau0_idx_0;
        } else {
          normtau0 /= 0.004;
        }
      } else if (tau0_idx_2 > 0.004) {
        normtau0 /= tau0_idx_2;
      } else {
        normtau0 /= 0.004;
      }
      if (normtau0 > tau0_idx_1) {
        tau0_idx_1 = normtau0;
      }
      normtau0 = fabs(f0[3]);
      tau0_idx_0 = fabs(q0[3]);
      tau0_idx_2 = fabs(ynew[3]);
      if (tau0_idx_0 > tau0_idx_2) {
        if (tau0_idx_0 > 0.004) {
          normtau0 /= tau0_idx_0;
        } else {
          normtau0 /= 0.004;
        }
      } else if (tau0_idx_2 > 0.004) {
        normtau0 /= tau0_idx_2;
      } else {
        normtau0 /= 0.004;
      }
      if (normtau0 > tau0_idx_1) {
        tau0_idx_1 = normtau0;
      }
      tau0_idx_0 = absh * tau0_idx_1;
      if (tau0_idx_0 > 2.5E-14) {
        if (absh <= hmin) {
          MinStepExit = true;
          exitg2 = 1;
        } else {
          if (NoFailedAttempts) {
            NoFailedAttempts = false;
            absh = fmax(hmin,
                        absh * fmax(0.1, 0.8 * pow(2.5E-14 / tau0_idx_0, 0.2)));
          } else {
            absh = fmax(hmin, 0.5 * absh);
          }
          h = (double)tdir * absh;
          Done = false;
        }
      } else {
        exitg2 = 1;
      }
    } while (exitg2 == 0);
    if (MinStepExit) {
      exitg1 = 1;
    } else {
      nnxt = next;
      while ((nnxt + 2 <= xg->size[0]) &&
             ((double)tdir * (tnew - xg_data[nnxt + 1]) >= 0.0)) {
        nnxt++;
      }
      Bcolidx = nnxt - next;
      if (Bcolidx > 0) {
        for (j = next + 2; j <= nnxt; j++) {
          normtau0 = (xg_data[j - 1] - t) / h;
          for (iac = 0; iac < 4; iac++) {
            c_d = 0.0;
            d1 = 0.0;
            d2 = 0.0;
            for (i = 0; i < 7; i++) {
              d3 = f[iac + (i << 2)];
              c_d += d3 * (h * b_b[i]);
              d1 += d3 * (h * c_b[i]);
              d2 += d3 * (h * d_b[i]);
            }
            yout_data[iac + 4 * (j - 1)] =
                (((d2 * normtau0 + d1) * normtau0 + c_d) * normtau0 +
                 f[iac] * h) *
                    normtau0 +
                q0[iac];
          }
        }
        if (xg_data[nnxt] == tnew) {
          yout_data[4 * nnxt] = ynew[0];
          yout_data[4 * nnxt + 1] = ynew[1];
          yout_data[4 * nnxt + 2] = ynew[2];
          yout_data[4 * nnxt + 3] = ynew[3];
        } else {
          normtau0 = (xg_data[nnxt] - t) / h;
          for (j = 0; j < 4; j++) {
            c_d = 0.0;
            d1 = 0.0;
            d2 = 0.0;
            for (i = 0; i < 7; i++) {
              d3 = f[j + (i << 2)];
              c_d += d3 * (h * b_b[i]);
              d1 += d3 * (h * c_b[i]);
              d2 += d3 * (h * d_b[i]);
            }
            yout_data[j + 4 * nnxt] =
                (((d2 * normtau0 + d1) * normtau0 + c_d) * normtau0 +
                 f[j] * h) *
                    normtau0 +
                q0[j];
          }
        }
        nout += Bcolidx;
        next = nnxt;
      }
      if (Done) {
        exitg1 = 1;
      } else {
        if (NoFailedAttempts) {
          normtau0 = 1.25 * pow(tau0_idx_0 / 2.5E-14, 0.2);
          if (normtau0 > 0.2) {
            absh /= normtau0;
          } else {
            absh *= 5.0;
          }
        }
        t = tnew;
        q0[0] = ynew[0];
        f[0] = f[24];
        q0[1] = ynew[1];
        f[1] = f[25];
        q0[2] = ynew[2];
        f[2] = f[26];
        q0[3] = ynew[3];
        f[3] = f[27];
      }
    }
  } while (exitg1 == 0);
  emxFreeStruct_captured_var(&b_knots);
  emxFreeStruct_captured_var(&b_P);
  emxInit_real_T(&Q, 2);
  i = Q->size[0] * Q->size[1];
  Q->size[0] = nout;
  Q->size[1] = 4;
  emxEnsureCapacity_real_T(Q, i);
  Q_data = Q->data;
  for (i = 0; i < 4; i++) {
    for (i1 = 0; i1 < nout; i1++) {
      Q_data[i1 + Q->size[0] * i] = yout_data[i + 4 * i1];
    }
  }
  emxFree_real_T(&yout);
  i = R->size[0] * R->size[1];
  R->size[0] = 3;
  R->size[1] = 3 * xg->size[0];
  emxEnsureCapacity_real_T(R, i);
  yout_data = R->data;
  nnxt = 3 * (3 * xg->size[0]);
  for (i = 0; i < nnxt; i++) {
    yout_data[i] = 0.0;
  }
  i = xg->size[0];
  for (tdir = 0; tdir < i; tdir++) {
    c_d = Q_data[tdir];
    hmin = c_d * c_d;
    c_d = Q_data[tdir + Q->size[0]];
    h = c_d * c_d;
    normtau0 = hmin + h;
    c_d = Q_data[tdir + Q->size[0] * 2];
    tnew = c_d * c_d;
    c_d = Q_data[tdir + Q->size[0] * 3];
    t = c_d * c_d;
    hmax = sqrt((normtau0 + tnew) + t);
    hmax = 1.0 / (hmax * hmax);
    c_d = 3.0 * (((double)tdir + 1.0) - 1.0) + 1.0;
    d1 = Q_data[tdir];
    d2 = Q_data[tdir + Q->size[0]];
    d3 = Q_data[tdir + Q->size[0] * 2];
    absh = Q_data[tdir + Q->size[0] * 3];
    i1 = 3 * ((int)c_d - 1);
    yout_data[i1] = hmax * ((normtau0 - tnew) - t);
    normtau0 = d2 * d3;
    tau0_idx_0 = d1 * absh;
    nnxt = 3 * ((int)(c_d + 1.0) - 1);
    yout_data[nnxt] = hmax * (2.0 * (normtau0 - tau0_idx_0));
    tau0_idx_2 = d2 * absh;
    tau0_idx_1 = d1 * d3;
    Bcolidx = 3 * ((int)(c_d + 2.0) - 1);
    yout_data[Bcolidx] = hmax * (2.0 * (tau0_idx_2 + tau0_idx_1));
    yout_data[i1 + 1] = hmax * (2.0 * (normtau0 + tau0_idx_0));
    yout_data[nnxt + 1] = hmax * (((hmin + tnew) - h) - t);
    normtau0 = d3 * absh;
    tau0_idx_0 = d1 * d2;
    yout_data[Bcolidx + 1] = hmax * (2.0 * (normtau0 - tau0_idx_0));
    yout_data[i1 + 2] = hmax * (2.0 * (tau0_idx_2 - tau0_idx_1));
    yout_data[nnxt + 2] = hmax * (2.0 * (normtau0 + tau0_idx_0));
    yout_data[Bcolidx + 2] = hmax * (((hmin + t) - h) - tnew);
  }
  emxFree_real_T(&Q);
}

/* End of code generation (getBishopFrame.c) */
