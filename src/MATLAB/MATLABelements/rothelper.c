/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * rothelper.c
 *
 * Code generation for function 'rothelper'
 *
 */

/* Include files */
#include "rothelper.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void b_rothelper(const double phi[3], double a[11])
{
  double absxk;
  double nphi;
  double s;
  double scale;
  double t;
  /*  Computes helper functions associated with rotation matrix */
  /*  phi = exponential coordinates of rotation */
  /*  if flag == 0, compute a(1:2) */
  /*  if flag == 1, compute a(1:4,9) */
  /*  if flag == 2, compute a(1:4,9:10) */
  /*  if flag == 3, compute a(1:11) */
  memset(&a[0], 0, 11U * sizeof(double));
  scale = 3.3121686421112381E-170;
  absxk = fabs(phi[0]);
  if (absxk > 3.3121686421112381E-170) {
    nphi = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    nphi = t * t;
  }
  absxk = fabs(phi[1]);
  if (absxk > scale) {
    t = scale / absxk;
    nphi = nphi * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nphi += t * t;
  }
  absxk = fabs(phi[2]);
  if (absxk > scale) {
    t = scale / absxk;
    nphi = nphi * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nphi += t * t;
  }
  nphi = scale * sqrt(nphi);
  t = cos(nphi);
  s = sin(nphi);
  if (nphi > 1.0E-5) {
    double nphi2_tmp;
    nphi2_tmp = nphi * nphi;
    a[0] = s / nphi;
    a[1] = (1.0 - t) / nphi2_tmp;
    scale = nphi * nphi2_tmp;
    absxk = nphi * scale;
    a[2] = (nphi * t - s) / scale;
    a[3] = (nphi * s - 2.0 * (1.0 - t)) / absxk;
    a[8] = (nphi - s) / scale;
    scale = nphi * absxk;
    a[9] = (3.0 * s - nphi * (t + 2.0)) / scale;
    absxk = nphi * scale;
    a[4] = -(3.0 * nphi * t + (nphi2_tmp - 3.0) * s) / scale;
    a[5] = (((nphi2_tmp - 8.0) * t + 8.0) - 5.0 * nphi * s) / absxk;
    a[10] = ((8.0 * nphi + 7.0 * nphi * t) + (nphi2_tmp - 15.0) * s) /
            (nphi * absxk);
  } else {
    a[0] = nphi * (nphi * (nphi * (nphi * 0.0083333333333333332) -
                           0.16666666666666666)) +
           1.0;
    a[1] = nphi * (nphi * (nphi * (nphi * 0.0013888888888888889) -
                           0.041666666666666664)) +
           0.5;
    a[2] = nphi * (nphi * (nphi * (nphi * -0.0011904761904761906) +
                           0.033333333333333333)) -
           0.33333333333333331;
    a[3] = nphi * (nphi * (nphi * (nphi * -0.00014880952380952382) +
                           0.0055555555555555558)) -
           0.083333333333333329;
    a[8] = nphi * (nphi * (nphi * (nphi * 0.00019841269841269841) -
                           0.0083333333333333332)) +
           0.16666666666666666;
    a[9] = nphi * (nphi * (nphi * (nphi * -1.6534391534391536E-5) +
                           0.00079365079365079365)) -
           0.016666666666666666;
    a[4] = nphi * (nphi * (nphi * (nphi * 0.00013227513227513228) -
                           0.0047619047619047623)) +
           0.066666666666666666;
    a[5] = nphi * (nphi * (nphi * (nphi * 1.3227513227513228E-5) -
                           0.00059523809523809529)) +
           0.011111111111111112;
    a[10] = nphi * (nphi * (nphi * (nphi * 1.2025012025012026E-6) -
                            6.6137566137566142E-5)) +
            0.0015873015873015873;
  }
}

void rothelper(const double phi[3], double a[11])
{
  double absxk;
  double nphi;
  double s;
  double scale;
  double t;
  /*  Computes helper functions associated with rotation matrix */
  /*  phi = exponential coordinates of rotation */
  /*  if flag == 0, compute a(1:2) */
  /*  if flag == 1, compute a(1:4,9) */
  /*  if flag == 2, compute a(1:4,9:10) */
  /*  if flag == 3, compute a(1:11) */
  memset(&a[0], 0, 11U * sizeof(double));
  scale = 3.3121686421112381E-170;
  absxk = fabs(phi[0]);
  if (absxk > 3.3121686421112381E-170) {
    nphi = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    nphi = t * t;
  }
  absxk = fabs(phi[1]);
  if (absxk > scale) {
    t = scale / absxk;
    nphi = nphi * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nphi += t * t;
  }
  absxk = fabs(phi[2]);
  if (absxk > scale) {
    t = scale / absxk;
    nphi = nphi * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nphi += t * t;
  }
  nphi = scale * sqrt(nphi);
  t = cos(nphi);
  s = sin(nphi);
  if (nphi > 1.0E-5) {
    scale = nphi * nphi;
    a[0] = s / nphi;
    a[1] = (1.0 - t) / scale;
    scale *= nphi;
    absxk = nphi * scale;
    a[2] = (nphi * t - s) / scale;
    a[3] = (nphi * s - 2.0 * (1.0 - t)) / absxk;
    a[8] = (nphi - s) / scale;
    a[9] = (3.0 * s - nphi * (t + 2.0)) / (nphi * absxk);
  } else {
    a[0] = nphi * (nphi * (nphi * (nphi * 0.0083333333333333332) -
                           0.16666666666666666)) +
           1.0;
    a[1] = nphi * (nphi * (nphi * (nphi * 0.0013888888888888889) -
                           0.041666666666666664)) +
           0.5;
    a[2] = nphi * (nphi * (nphi * (nphi * -0.0011904761904761906) +
                           0.033333333333333333)) -
           0.33333333333333331;
    a[3] = nphi * (nphi * (nphi * (nphi * -0.00014880952380952382) +
                           0.0055555555555555558)) -
           0.083333333333333329;
    a[8] = nphi * (nphi * (nphi * (nphi * 0.00019841269841269841) -
                           0.0083333333333333332)) +
           0.16666666666666666;
    a[9] = nphi * (nphi * (nphi * (nphi * -1.6534391534391536E-5) +
                           0.00079365079365079365)) -
           0.016666666666666666;
  }
}

/* End of code generation (rothelper.c) */
