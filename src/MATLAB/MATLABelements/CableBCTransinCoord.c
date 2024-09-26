/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * CableBCTransinCoord.c
 *
 * Code generation for function 'CableBCTransinCoord'
 *
 */

/* Include files */
#include "CableBCTransinCoord.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_types.h"
#include "rothelper.h"
#include <math.h>
#include <string.h>

/* Function Definitions */
void CableBCTransinCoord(const double qbar[7], const double x0[3],
                         const double RJ[9], const double RE[9],
                         const double r[3], const double R0[9],
                         const double eta[7], const double rho1[7], double q[7],
                         double J[49], double Q[49], double Qtilde[49])
{
  double a[11];
  double D2rotr[9];
  double Drotr[9];
  double Dtau[9];
  double Dws1[9];
  double RJdexpT_[9];
  double Rb[9];
  double b_a[9];
  double dv[9];
  double e_a[9];
  double hphi[9];
  double hphi2[9];
  double tau0_dot_tau[9];
  double tauhat[9];
  double Thet_tilde_e2[3];
  double p2_tmp[3];
  double rotr[3];
  double tau0_cross_tau[3];
  double v[3];
  double absxk;
  double c_a;
  double c_tmp;
  double d;
  double d1;
  double d10;
  double d11;
  double d12;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double d_a;
  double nu_;
  double scale;
  double t;
  int Qtilde_tmp;
  int b_Qtilde_tmp;
  int b_i;
  int c_Qtilde_tmp;
  int i;
  /*  Computes the transformation from (d, phi, gamma) -> (p1, p2, vartheta) for
   */
  /*  one end of the cable; NOTE: the main difference from CableBCtrans.m is */
  /*  that this function is based on exponential coordinates for the rotation,
   */
  /*  rather than the rotation matrix; also, a joint coordinate system and end
   */
  /*  offset are included, making it a more general. */
  /*  Inputs */
  /*  qbar = [d; phi; gamm] */
  /*    where d = displacement of joint,  */
  /*          phi = exponential coordinates of rotation of joint */
  /*          gamm = distance between first and second control points */
  /*  x0 = reference position of joint (3x1) vector */
  /*  RJ = rotation of joint coordinate frame with respect to global */
  /*  RE = rotation of cable end with respect to joint coordinate frame */
  /*  r = position of cable end relative to joint in joint coordinate system */
  /*      (end offset, 3x1 vector) */
  /*  R0 = orientation of cable end in reference configuration */
  /*  eta = vector to contract with for second derivative, Q */
  /*  rho1 = vector to contract with for second derivative, Qtilde */
  /*  rho2 = vector to contract with for third derivative, C */
  /*  Outputs: */
  /*  q = [p1; p2; vartheta1] */
  /*    where p1, p2 = positions of first and second control points */
  /*          vartheta1 = (twist) rotation of first control point */
  /*  J = first derivative of transformation */
  /*  Q = second derivative contracted over upper index,  */
  /*      Q_{ib,jb} = eta_i q^i_{ib,jb} */
  /*  Qtilde = second derivative contracted over one of the lower indices */
  /*      Qtilde^{i}_{ib} = rho1^{jb} q^i_{ib,jb} */
  /*  C = third derivative contracted over two lower indices */
  /*      C^{i}_{ib} = q^i_{ib,jb,kb}rho1^{jb}rho2^{kb} */
  rothelper(&qbar[3], a);
  hphi[0] = 0.0;
  hphi[3] = -qbar[5];
  hphi[6] = qbar[4];
  hphi[1] = qbar[5];
  hphi[4] = 0.0;
  hphi[7] = -qbar[3];
  hphi[2] = -qbar[4];
  hphi[5] = qbar[3];
  hphi[8] = 0.0;
  for (i = 0; i < 3; i++) {
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      hphi2[i + 3 * Qtilde_tmp] = (hphi[i] * hphi[3 * Qtilde_tmp] +
                                   hphi[i + 3] * hphi[3 * Qtilde_tmp + 1]) +
                                  hphi[i + 6] * hphi[3 * Qtilde_tmp + 2];
    }
  }
  d = a[0];
  d1 = a[1];
  for (i = 0; i < 9; i++) {
    b_a[i] = (d * hphi[i] + (double)iv[i]) + d1 * hphi2[i];
  }
  absxk = 0.0;
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      Drotr[i + 3 * Qtilde_tmp] =
          (d * b_a[3 * Qtilde_tmp] + d1 * b_a[3 * Qtilde_tmp + 1]) +
          d2 * b_a[3 * Qtilde_tmp + 2];
    }
    d = Drotr[i];
    d1 = Drotr[i + 3];
    d2 = Drotr[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      Rb[i + 3 * Qtilde_tmp] =
          (d * RE[3 * Qtilde_tmp] + d1 * RE[3 * Qtilde_tmp + 1]) +
          d2 * RE[3 * Qtilde_tmp + 2];
    }
    absxk += R0[i] * Rb[i];
  }
  tau0_cross_tau[0] = R0[1] * Rb[2] - Rb[1] * R0[2];
  tau0_cross_tau[1] = Rb[0] * R0[2] - R0[0] * Rb[2];
  tau0_cross_tau[2] = R0[0] * Rb[1] - Rb[0] * R0[1];
  dv[0] = 0.0;
  dv[1] = -tau0_cross_tau[2];
  dv[2] = tau0_cross_tau[1];
  dv[3] = tau0_cross_tau[2];
  dv[4] = 0.0;
  dv[5] = -tau0_cross_tau[0];
  dv[6] = -tau0_cross_tau[1];
  dv[7] = tau0_cross_tau[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d = tau0_cross_tau[i] / (absxk + 1.0);
    tau0_dot_tau[3 * i] =
        (absxk * (double)iv[i] + dv[3 * i]) + d * tau0_cross_tau[0];
    b_i = 3 * i + 1;
    tau0_dot_tau[b_i] =
        (absxk * (double)iv[i + 3] + dv[b_i]) + d * tau0_cross_tau[1];
    b_i = 3 * i + 2;
    tau0_dot_tau[b_i] =
        (absxk * (double)iv[i + 6] + dv[b_i]) + d * tau0_cross_tau[2];
  }
  d = Rb[3];
  d1 = Rb[4];
  d2 = Rb[5];
  for (i = 0; i < 3; i++) {
    tau0_cross_tau[i] = (tau0_dot_tau[i] * d + tau0_dot_tau[i + 3] * d1) +
                        tau0_dot_tau[i + 6] * d2;
  }
  /*  first derivative needed */
  nu_ = 0.0;
  scale = 3.3121686421112381E-170;
  d = r[0];
  d1 = r[1];
  d2 = r[2];
  d3 = tau0_cross_tau[0];
  d4 = tau0_cross_tau[1];
  d5 = tau0_cross_tau[2];
  d6 = qbar[0];
  d7 = qbar[1];
  d8 = qbar[2];
  d9 = qbar[6];
  for (b_i = 0; b_i < 3; b_i++) {
    d10 = (Drotr[b_i] * d + Drotr[b_i + 3] * d1) + Drotr[b_i + 6] * d2;
    rotr[b_i] = d10;
    Thet_tilde_e2[b_i] =
        (R0[3 * b_i] * d3 + R0[3 * b_i + 1] * d4) + R0[3 * b_i + 2] * d5;
    d10 += x0[b_i] + ((RJ[b_i] * d6 + RJ[b_i + 3] * d7) + RJ[b_i + 6] * d8);
    d11 = Rb[b_i];
    d12 = d9 * d11;
    p2_tmp[b_i] = d12;
    q[b_i] = d10;
    q[b_i + 3] = d10 + d12;
    d10 = d11 + R0[b_i];
    tau0_cross_tau[b_i] = d10;
    absxk = fabs(d10);
    if (absxk > scale) {
      t = scale / absxk;
      nu_ = nu_ * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      nu_ += t * t;
    }
  }
  q[6] = atan2(Thet_tilde_e2[2], Thet_tilde_e2[1]);
  nu_ = scale * sqrt(nu_);
  c_tmp = nu_ * nu_;
  d = a[1];
  d1 = a[8];
  for (b_i = 0; b_i < 3; b_i++) {
    v[b_i] = 2.0 * tau0_cross_tau[b_i] / c_tmp;
    tau0_dot_tau[3 * b_i] = ((double)iv[b_i] - d * hphi[b_i]) + d1 * hphi2[b_i];
    tau0_dot_tau[3 * b_i + 1] =
        ((double)iv[b_i + 3] - d * hphi[b_i + 3]) + d1 * hphi2[b_i + 3];
    tau0_dot_tau[3 * b_i + 2] =
        ((double)iv[b_i + 6] - d * hphi[b_i + 6]) + d1 * hphi2[b_i + 6];
  }
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      RJdexpT_[i + 3 * Qtilde_tmp] = (d * tau0_dot_tau[3 * Qtilde_tmp] +
                                      d1 * tau0_dot_tau[3 * Qtilde_tmp + 1]) +
                                     d2 * tau0_dot_tau[3 * Qtilde_tmp + 2];
    }
  }
  hphi[0] = 0.0;
  hphi[3] = -rotr[2];
  hphi[6] = rotr[1];
  hphi[1] = rotr[2];
  hphi[4] = 0.0;
  hphi[7] = -rotr[0];
  hphi[2] = -rotr[1];
  hphi[5] = rotr[0];
  hphi[8] = 0.0;
  tauhat[0] = 0.0;
  tauhat[3] = -Rb[2];
  tauhat[6] = Rb[1];
  tauhat[1] = Rb[2];
  tauhat[4] = 0.0;
  tauhat[7] = -Rb[0];
  tauhat[2] = -Rb[1];
  tauhat[5] = Rb[0];
  tauhat[8] = 0.0;
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] = -hphi[i];
  }
  for (i = 0; i < 3; i++) {
    d = tau0_dot_tau[i];
    d1 = tau0_dot_tau[i + 3];
    d2 = tau0_dot_tau[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      Drotr[i + 3 * Qtilde_tmp] =
          (d * RJdexpT_[3 * Qtilde_tmp] + d1 * RJdexpT_[3 * Qtilde_tmp + 1]) +
          d2 * RJdexpT_[3 * Qtilde_tmp + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] = -tauhat[i];
  }
  for (i = 0; i < 3; i++) {
    d = tau0_dot_tau[i];
    d1 = tau0_dot_tau[i + 3];
    d2 = tau0_dot_tau[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      Dtau[i + 3 * Qtilde_tmp] =
          (d * RJdexpT_[3 * Qtilde_tmp] + d1 * RJdexpT_[3 * Qtilde_tmp + 1]) +
          d2 * RJdexpT_[3 * Qtilde_tmp + 2];
    }
  }
  memset(&J[0], 0, 49U * sizeof(double));
  /*  second derivative needed */
  /*  ws1 is short form for w(phi;dphi1) */
  hphi2[0] = 0.0;
  hphi2[4] = 0.0;
  hphi2[8] = 0.0;
  /*  Dws1 is short form for Dw(phi;dphi1) */
  scale = 0.0;
  d = v[0];
  d1 = v[1];
  d2 = v[2];
  d3 = a[2];
  d4 = a[3] * (rho1[4] * qbar[5] - qbar[4] * rho1[5]);
  d5 = a[3] * (qbar[3] * rho1[5] - rho1[3] * qbar[5]);
  d6 = a[3] * (rho1[3] * qbar[4] - qbar[3] * rho1[4]);
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    absxk = qbar[i + 3];
    scale += absxk * rho1[i + 3];
    d7 = RJ[3 * i];
    J[7 * i] = d7;
    J[7 * i + 3] = d7;
    d7 = Drotr[3 * i];
    J[b_i] = d7;
    J[b_i + 3] = d7 + qbar[6] * Dtau[3 * i];
    Drotr[3 * i] = qbar[3] * absxk;
    b_a[3 * i] = d3 * rho1[3] * absxk;
    Dws1[3 * i] = d4 * absxk;
    Qtilde_tmp = 3 * i + 1;
    d7 = RJ[Qtilde_tmp];
    J[7 * i + 1] = d7;
    J[7 * i + 4] = d7;
    d7 = Drotr[Qtilde_tmp];
    J[b_i + 1] = d7;
    J[b_i + 4] = d7 + qbar[6] * Dtau[Qtilde_tmp];
    Drotr[Qtilde_tmp] = qbar[4] * absxk;
    b_a[Qtilde_tmp] = d3 * rho1[4] * absxk;
    Dws1[Qtilde_tmp] = d5 * absxk;
    b_Qtilde_tmp = 3 * i + 2;
    d7 = RJ[b_Qtilde_tmp];
    J[7 * i + 2] = d7;
    J[7 * i + 5] = d7;
    d7 = Drotr[b_Qtilde_tmp];
    J[b_i + 2] = d7;
    J[b_i + 5] = d7 + qbar[6] * Dtau[b_Qtilde_tmp];
    J[b_i + 6] = (d * RJdexpT_[3 * i] + d1 * RJdexpT_[Qtilde_tmp]) +
                 d2 * RJdexpT_[b_Qtilde_tmp];
    Drotr[b_Qtilde_tmp] = qbar[5] * absxk;
    b_a[b_Qtilde_tmp] = d3 * rho1[5] * absxk;
    Dws1[b_Qtilde_tmp] = d6 * absxk;
    Thet_tilde_e2[i] = (RJdexpT_[i] * rho1[3] + RJdexpT_[i + 3] * rho1[4]) +
                       RJdexpT_[i + 6] * rho1[5];
    J[i + 45] = Rb[i];
  }
  hphi2[3] = -Thet_tilde_e2[2];
  hphi2[6] = Thet_tilde_e2[1];
  hphi2[1] = Thet_tilde_e2[2];
  hphi2[7] = -Thet_tilde_e2[0];
  hphi2[2] = -Thet_tilde_e2[1];
  hphi2[5] = Thet_tilde_e2[0];
  c_a = a[9] * scale;
  d_a = a[8] * scale;
  dv[0] = 0.0;
  dv[3] = a[1] * -rho1[5];
  dv[6] = a[1] * rho1[4];
  dv[1] = a[1] * rho1[5];
  dv[4] = 0.0;
  dv[7] = a[1] * -rho1[3];
  dv[2] = a[1] * -rho1[4];
  dv[5] = a[1] * rho1[3];
  dv[8] = 0.0;
  d = qbar[3];
  d1 = qbar[4];
  d2 = qbar[5];
  d3 = a[8];
  for (i = 0; i < 3; i++) {
    absxk = rho1[i + 3];
    e_a[3 * i] =
        ((((b_a[3 * i] - Dws1[3 * i]) - dv[3 * i]) + c_a * Drotr[3 * i]) +
         d * d3 * absxk) +
        d_a * (double)iv[3 * i];
    b_i = 3 * i + 1;
    e_a[b_i] = ((((b_a[b_i] - Dws1[b_i]) - dv[b_i]) + c_a * Drotr[b_i]) +
                d1 * d3 * absxk) +
               d_a * (double)iv[b_i];
    b_i = 3 * i + 2;
    e_a[b_i] = ((((b_a[b_i] - Dws1[b_i]) - dv[b_i]) + c_a * Drotr[b_i]) +
                d2 * d3 * absxk) +
               d_a * (double)iv[b_i];
  }
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      Dws1[i + 3 * Qtilde_tmp] =
          (d * e_a[3 * Qtilde_tmp] + d1 * e_a[3 * Qtilde_tmp + 1]) +
          d2 * e_a[3 * Qtilde_tmp + 2];
    }
  }
  /*  D2rotr is short form for D^2rotr(phi)(dphi1,dot) */
  for (i = 0; i < 3; i++) {
    d = hphi2[i];
    d1 = hphi2[i + 3];
    d2 = hphi2[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      b_Qtilde_tmp = 3 * Qtilde_tmp + 1;
      c_Qtilde_tmp = 3 * Qtilde_tmp + 2;
      b_i = i + 3 * Qtilde_tmp;
      D2rotr[b_i] =
          (hphi[i] * Dws1[3 * Qtilde_tmp] + hphi[i + 3] * Dws1[b_Qtilde_tmp]) +
          hphi[i + 6] * Dws1[c_Qtilde_tmp];
      e_a[b_i] = (d * hphi[3 * Qtilde_tmp] + d1 * hphi[b_Qtilde_tmp]) +
                 d2 * hphi[c_Qtilde_tmp];
    }
    d = e_a[i];
    d1 = e_a[i + 3];
    d2 = e_a[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      tau0_dot_tau[i + 3 * Qtilde_tmp] =
          (d * RJdexpT_[3 * Qtilde_tmp] + d1 * RJdexpT_[3 * Qtilde_tmp + 1]) +
          d2 * RJdexpT_[3 * Qtilde_tmp + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    D2rotr[i] = -(D2rotr[i] + tau0_dot_tau[i]);
  }
  /*  D2tau1 is short form for D^2tau(phi)(dphi1,dot) */
  tau0_cross_tau[0] /= nu_;
  tau0_cross_tau[1] /= nu_;
  tau0_cross_tau[2] /= nu_;
  /*  unit vector in the direction of u */
  c_a = -(2.0 / c_tmp);
  for (i = 0; i < 3; i++) {
    tau0_dot_tau[3 * i] =
        (double)iv[3 * i] - tau0_cross_tau[0] * tau0_cross_tau[i];
    b_i = 3 * i + 1;
    tau0_dot_tau[b_i] = (double)iv[b_i] - tau0_cross_tau[1] * tau0_cross_tau[i];
    b_i = 3 * i + 2;
    tau0_dot_tau[b_i] = (double)iv[b_i] - tau0_cross_tau[2] * tau0_cross_tau[i];
  }
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] = c_a * (2.0 * tau0_dot_tau[i] - (double)iv[i]);
  }
  for (i = 0; i < 3; i++) {
    d = tau0_dot_tau[i];
    d1 = tau0_dot_tau[i + 3];
    d2 = tau0_dot_tau[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      b_a[i + 3 * Qtilde_tmp] =
          (d * tauhat[3 * Qtilde_tmp] + d1 * tauhat[3 * Qtilde_tmp + 1]) +
          d2 * tauhat[3 * Qtilde_tmp + 2];
    }
    d = b_a[i];
    d1 = b_a[i + 3];
    d2 = b_a[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      hphi[i + 3 * Qtilde_tmp] =
          (d * RJdexpT_[3 * Qtilde_tmp] + d1 * RJdexpT_[3 * Qtilde_tmp + 1]) +
          d2 * RJdexpT_[3 * Qtilde_tmp + 2];
    }
  }
  memset(&Qtilde[0], 0, 49U * sizeof(double));
  for (i = 0; i < 3; i++) {
    d = hphi2[i];
    d1 = hphi2[i + 3];
    d2 = hphi2[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      Qtilde[Qtilde_tmp + 7 * (i + 3)] = D2rotr[Qtilde_tmp + 3 * i];
      b_Qtilde_tmp = 3 * Qtilde_tmp + 1;
      c_Qtilde_tmp = 3 * Qtilde_tmp + 2;
      b_i = i + 3 * Qtilde_tmp;
      tau0_dot_tau[b_i] = (tauhat[i] * Dws1[3 * Qtilde_tmp] +
                           tauhat[i + 3] * Dws1[b_Qtilde_tmp]) +
                          tauhat[i + 6] * Dws1[c_Qtilde_tmp];
      e_a[b_i] = (d * tauhat[3 * Qtilde_tmp] + d1 * tauhat[b_Qtilde_tmp]) +
                 d2 * tauhat[c_Qtilde_tmp];
    }
    d = e_a[i];
    d1 = e_a[i + 3];
    d2 = e_a[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      hphi2[i + 3 * Qtilde_tmp] =
          (d * RJdexpT_[3 * Qtilde_tmp] + d1 * RJdexpT_[3 * Qtilde_tmp + 1]) +
          d2 * RJdexpT_[3 * Qtilde_tmp + 2];
    }
  }
  d = qbar[6];
  d1 = rho1[6];
  d2 = rho1[3];
  d3 = rho1[4];
  d4 = rho1[5];
  d5 = Thet_tilde_e2[0];
  d6 = Thet_tilde_e2[1];
  d7 = Thet_tilde_e2[2];
  d8 = v[0];
  d9 = v[1];
  d10 = v[2];
  for (b_i = 0; b_i < 3; b_i++) {
    c_Qtilde_tmp = 7 * (b_i + 3);
    Qtilde[c_Qtilde_tmp + 3] =
        (D2rotr[3 * b_i] + d * -(tau0_dot_tau[3 * b_i] + hphi2[3 * b_i])) +
        d1 * Dtau[3 * b_i];
    b_Qtilde_tmp = 3 * b_i + 1;
    Qtilde[c_Qtilde_tmp + 4] =
        (D2rotr[b_Qtilde_tmp] +
         d * -(tau0_dot_tau[b_Qtilde_tmp] + hphi2[b_Qtilde_tmp])) +
        d1 * Dtau[b_Qtilde_tmp];
    Qtilde_tmp = 3 * b_i + 2;
    Qtilde[c_Qtilde_tmp + 5] =
        (D2rotr[Qtilde_tmp] +
         d * -(tau0_dot_tau[Qtilde_tmp] + hphi2[Qtilde_tmp])) +
        d1 * Dtau[Qtilde_tmp];
    Qtilde[c_Qtilde_tmp + 6] = ((d5 * hphi[3 * b_i] + d6 * hphi[b_Qtilde_tmp]) +
                                d7 * hphi[Qtilde_tmp]) +
                               ((d8 * Dws1[3 * b_i] + d9 * Dws1[b_Qtilde_tmp]) +
                                d10 * Dws1[Qtilde_tmp]);
    Qtilde[b_i + 45] =
        (Dtau[b_i] * d2 + Dtau[b_i + 3] * d3) + Dtau[b_i + 6] * d4;
    tau0_cross_tau[b_i] = eta[b_i] + eta[b_i + 3];
  }
  absxk = ((rotr[1] * tau0_cross_tau[2] - tau0_cross_tau[1] * rotr[2]) +
           qbar[6] * (Rb[1] * eta[5] - Rb[2] * eta[4])) +
          v[0] * eta[6];
  scale = ((tau0_cross_tau[0] * rotr[2] - rotr[0] * tau0_cross_tau[2]) +
           qbar[6] * (Rb[2] * eta[3] - Rb[0] * eta[5])) +
          v[1] * eta[6];
  t = ((rotr[0] * tau0_cross_tau[1] - tau0_cross_tau[0] * rotr[1]) +
       qbar[6] * (Rb[0] * eta[4] - Rb[1] * eta[3])) +
      v[2] * eta[6];
  c_tmp = 0.0;
  for (i = 0; i < 3; i++) {
    d = (RJ[3 * i] * absxk + RJ[3 * i + 1] * scale) + RJ[3 * i + 2] * t;
    Thet_tilde_e2[i] = d;
    c_tmp += d * qbar[i + 3];
  }
  c_a = a[9] * c_tmp;
  d_a = a[8] * c_tmp;
  memset(&Q[0], 0, 49U * sizeof(double));
  absxk = ((rotr[0] * tau0_cross_tau[0] + rotr[1] * tau0_cross_tau[1]) +
           rotr[2] * tau0_cross_tau[2]) +
          qbar[6] * ((Rb[0] * eta[3] + Rb[1] * eta[4]) + Rb[2] * eta[5]);
  dv[0] = 0.0;
  dv[3] = a[1] * -Thet_tilde_e2[2];
  dv[6] = a[1] * Thet_tilde_e2[1];
  dv[1] = a[1] * Thet_tilde_e2[2];
  dv[4] = 0.0;
  dv[7] = a[1] * -Thet_tilde_e2[0];
  dv[2] = a[1] * -Thet_tilde_e2[1];
  dv[5] = Thet_tilde_e2[0] * a[1];
  dv[8] = 0.0;
  d = Thet_tilde_e2[0];
  d1 = Thet_tilde_e2[1];
  d2 = Thet_tilde_e2[2];
  d3 = a[2];
  d4 = a[3] * (Thet_tilde_e2[2] * qbar[4] - Thet_tilde_e2[1] * qbar[5]);
  d5 = a[3] * (Thet_tilde_e2[0] * qbar[5] - Thet_tilde_e2[2] * qbar[3]);
  d6 = a[3] * (Thet_tilde_e2[1] * qbar[3] - Thet_tilde_e2[0] * qbar[4]);
  d7 = rotr[0];
  d8 = rotr[1];
  d9 = rotr[2];
  d10 = p2_tmp[0];
  d11 = p2_tmp[1];
  d12 = p2_tmp[2];
  for (i = 0; i < 3; i++) {
    t = qbar[i + 3];
    b_a[3 * i] = d * d3 * t;
    Dws1[3 * i] = d4 * t;
    c_tmp = tau0_cross_tau[i];
    e_a[3 * i] = d7 * c_tmp;
    scale = eta[i + 3];
    tau0_dot_tau[3 * i] = d10 * scale;
    b_i = 3 * i + 1;
    b_a[b_i] = d1 * d3 * t;
    Dws1[b_i] = d5 * t;
    e_a[b_i] = d8 * c_tmp;
    tau0_dot_tau[b_i] = d11 * scale;
    b_i = 3 * i + 2;
    b_a[b_i] = d3 * d2 * t;
    Dws1[b_i] = d6 * t;
    e_a[b_i] = d9 * c_tmp;
    tau0_dot_tau[b_i] = d12 * scale;
  }
  for (i = 0; i < 9; i++) {
    e_a[i] = (e_a[i] + tau0_dot_tau[i]) - absxk * (double)iv[i];
  }
  for (i = 0; i < 3; i++) {
    d = e_a[i];
    d1 = e_a[i + 3];
    d2 = e_a[i + 6];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      b_i = i + 3 * Qtilde_tmp;
      D2rotr[b_i] =
          ((d * RJdexpT_[3 * Qtilde_tmp] + d1 * RJdexpT_[3 * Qtilde_tmp + 1]) +
           d2 * RJdexpT_[3 * Qtilde_tmp + 2]) +
          eta[6] * hphi[b_i];
    }
  }
  for (i = 0; i < 3; i++) {
    d = RJdexpT_[3 * i];
    d1 = RJdexpT_[3 * i + 1];
    d2 = RJdexpT_[3 * i + 2];
    for (Qtilde_tmp = 0; Qtilde_tmp < 3; Qtilde_tmp++) {
      b_i = Qtilde_tmp + 3 * i;
      e_a[b_i] = ((((b_a[b_i] - Dws1[b_i]) + dv[b_i]) + c_a * Drotr[b_i]) +
                  d_a * (double)iv[b_i]) +
                 a[8] * qbar[Qtilde_tmp + 3] * Thet_tilde_e2[i];
      tau0_dot_tau[i + 3 * Qtilde_tmp] =
          (d * D2rotr[3 * Qtilde_tmp] + d1 * D2rotr[3 * Qtilde_tmp + 1]) +
          d2 * D2rotr[3 * Qtilde_tmp + 2];
    }
  }
  d = eta[3];
  d1 = eta[4];
  d2 = eta[5];
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    Q[b_i + 3] = e_a[3 * i] + tau0_dot_tau[3 * i];
    c_Qtilde_tmp = 3 * i + 1;
    Q[b_i + 4] = e_a[c_Qtilde_tmp] + tau0_dot_tau[c_Qtilde_tmp];
    b_Qtilde_tmp = 3 * i + 2;
    Q[b_i + 5] = e_a[b_Qtilde_tmp] + tau0_dot_tau[b_Qtilde_tmp];
    d3 = (d * Dtau[3 * i] + d1 * Dtau[c_Qtilde_tmp]) + d2 * Dtau[b_Qtilde_tmp];
    Q[b_i + 6] = d3;
    Q[i + 45] = d3;
  }
}

void b_CableBCTransinCoord(const double qbar[7], const double x0[3],
                           const double RJ[9], const double RE[9],
                           const double r[3], const double R0[9],
                           const emxArray_real_T *eta, const double rho1[7],
                           double q[7], double J[49], double Q[49],
                           double Qtilde[49])
{
  double a[11];
  double D2rotr[9];
  double Drotr[9];
  double Dtau[9];
  double Dws1[9];
  double RJdexpT_[9];
  double Rb[9];
  double b_a[9];
  double dv[9];
  double e_a[9];
  double hphi[9];
  double hphi2[9];
  double tau0_dot_tau[9];
  double tauhat[9];
  double Thet_tilde_e2[3];
  double p2_tmp[3];
  double rotr[3];
  double tau0_cross_tau[3];
  double v[3];
  const double *eta_data;
  double absxk;
  double c_a;
  double c_tmp;
  double d;
  double d1;
  double d10;
  double d11;
  double d12;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double d_a;
  double nu_;
  double scale;
  double t;
  int Qtilde_tmp;
  int b_Qtilde_tmp;
  int b_i;
  int i;
  int i1;
  eta_data = eta->data;
  /*  Computes the transformation from (d, phi, gamma) -> (p1, p2, vartheta) for
   */
  /*  one end of the cable; NOTE: the main difference from CableBCtrans.m is */
  /*  that this function is based on exponential coordinates for the rotation,
   */
  /*  rather than the rotation matrix; also, a joint coordinate system and end
   */
  /*  offset are included, making it a more general. */
  /*  Inputs */
  /*  qbar = [d; phi; gamm] */
  /*    where d = displacement of joint,  */
  /*          phi = exponential coordinates of rotation of joint */
  /*          gamm = distance between first and second control points */
  /*  x0 = reference position of joint (3x1) vector */
  /*  RJ = rotation of joint coordinate frame with respect to global */
  /*  RE = rotation of cable end with respect to joint coordinate frame */
  /*  r = position of cable end relative to joint in joint coordinate system */
  /*      (end offset, 3x1 vector) */
  /*  R0 = orientation of cable end in reference configuration */
  /*  eta = vector to contract with for second derivative, Q */
  /*  rho1 = vector to contract with for second derivative, Qtilde */
  /*  rho2 = vector to contract with for third derivative, C */
  /*  Outputs: */
  /*  q = [p1; p2; vartheta1] */
  /*    where p1, p2 = positions of first and second control points */
  /*          vartheta1 = (twist) rotation of first control point */
  /*  J = first derivative of transformation */
  /*  Q = second derivative contracted over upper index,  */
  /*      Q_{ib,jb} = eta_i q^i_{ib,jb} */
  /*  Qtilde = second derivative contracted over one of the lower indices */
  /*      Qtilde^{i}_{ib} = rho1^{jb} q^i_{ib,jb} */
  /*  C = third derivative contracted over two lower indices */
  /*      C^{i}_{ib} = q^i_{ib,jb,kb}rho1^{jb}rho2^{kb} */
  rothelper(&qbar[3], a);
  hphi[0] = 0.0;
  hphi[3] = -qbar[5];
  hphi[6] = qbar[4];
  hphi[1] = qbar[5];
  hphi[4] = 0.0;
  hphi[7] = -qbar[3];
  hphi[2] = -qbar[4];
  hphi[5] = qbar[3];
  hphi[8] = 0.0;
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      hphi2[i + 3 * i1] =
          (hphi[i] * hphi[3 * i1] + hphi[i + 3] * hphi[3 * i1 + 1]) +
          hphi[i + 6] * hphi[3 * i1 + 2];
    }
  }
  d = a[0];
  d1 = a[1];
  for (i = 0; i < 9; i++) {
    b_a[i] = (d * hphi[i] + (double)iv[i]) + d1 * hphi2[i];
  }
  absxk = 0.0;
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Drotr[i + 3 * i1] =
          (d * b_a[3 * i1] + d1 * b_a[3 * i1 + 1]) + d2 * b_a[3 * i1 + 2];
    }
    d = Drotr[i];
    d1 = Drotr[i + 3];
    d2 = Drotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Rb[i + 3 * i1] =
          (d * RE[3 * i1] + d1 * RE[3 * i1 + 1]) + d2 * RE[3 * i1 + 2];
    }
    absxk += R0[i] * Rb[i];
  }
  tau0_cross_tau[0] = R0[1] * Rb[2] - Rb[1] * R0[2];
  tau0_cross_tau[1] = Rb[0] * R0[2] - R0[0] * Rb[2];
  tau0_cross_tau[2] = R0[0] * Rb[1] - Rb[0] * R0[1];
  dv[0] = 0.0;
  dv[1] = -tau0_cross_tau[2];
  dv[2] = tau0_cross_tau[1];
  dv[3] = tau0_cross_tau[2];
  dv[4] = 0.0;
  dv[5] = -tau0_cross_tau[0];
  dv[6] = -tau0_cross_tau[1];
  dv[7] = tau0_cross_tau[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d = tau0_cross_tau[i] / (absxk + 1.0);
    tau0_dot_tau[3 * i] =
        (absxk * (double)iv[i] + dv[3 * i]) + d * tau0_cross_tau[0];
    b_i = 3 * i + 1;
    tau0_dot_tau[b_i] =
        (absxk * (double)iv[i + 3] + dv[b_i]) + d * tau0_cross_tau[1];
    b_i = 3 * i + 2;
    tau0_dot_tau[b_i] =
        (absxk * (double)iv[i + 6] + dv[b_i]) + d * tau0_cross_tau[2];
  }
  d = Rb[3];
  d1 = Rb[4];
  d2 = Rb[5];
  for (i = 0; i < 3; i++) {
    tau0_cross_tau[i] = (tau0_dot_tau[i] * d + tau0_dot_tau[i + 3] * d1) +
                        tau0_dot_tau[i + 6] * d2;
  }
  /*  first derivative needed */
  nu_ = 0.0;
  scale = 3.3121686421112381E-170;
  d = r[0];
  d1 = r[1];
  d2 = r[2];
  d3 = tau0_cross_tau[0];
  d4 = tau0_cross_tau[1];
  d5 = tau0_cross_tau[2];
  d6 = qbar[0];
  d7 = qbar[1];
  d8 = qbar[2];
  d9 = qbar[6];
  for (b_i = 0; b_i < 3; b_i++) {
    d10 = (Drotr[b_i] * d + Drotr[b_i + 3] * d1) + Drotr[b_i + 6] * d2;
    rotr[b_i] = d10;
    Thet_tilde_e2[b_i] =
        (R0[3 * b_i] * d3 + R0[3 * b_i + 1] * d4) + R0[3 * b_i + 2] * d5;
    d10 += x0[b_i] + ((RJ[b_i] * d6 + RJ[b_i + 3] * d7) + RJ[b_i + 6] * d8);
    d11 = Rb[b_i];
    d12 = d9 * d11;
    p2_tmp[b_i] = d12;
    q[b_i] = d10;
    q[b_i + 3] = d10 + d12;
    d10 = d11 + R0[b_i];
    tau0_cross_tau[b_i] = d10;
    absxk = fabs(d10);
    if (absxk > scale) {
      t = scale / absxk;
      nu_ = nu_ * t * t + 1.0;
      scale = absxk;
    } else {
      t = absxk / scale;
      nu_ += t * t;
    }
  }
  q[6] = atan2(Thet_tilde_e2[2], Thet_tilde_e2[1]);
  nu_ = scale * sqrt(nu_);
  c_tmp = nu_ * nu_;
  d = a[1];
  d1 = a[8];
  for (b_i = 0; b_i < 3; b_i++) {
    v[b_i] = 2.0 * tau0_cross_tau[b_i] / c_tmp;
    tau0_dot_tau[3 * b_i] = ((double)iv[b_i] - d * hphi[b_i]) + d1 * hphi2[b_i];
    tau0_dot_tau[3 * b_i + 1] =
        ((double)iv[b_i + 3] - d * hphi[b_i + 3]) + d1 * hphi2[b_i + 3];
    tau0_dot_tau[3 * b_i + 2] =
        ((double)iv[b_i + 6] - d * hphi[b_i + 6]) + d1 * hphi2[b_i + 6];
  }
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      RJdexpT_[i + 3 * i1] =
          (d * tau0_dot_tau[3 * i1] + d1 * tau0_dot_tau[3 * i1 + 1]) +
          d2 * tau0_dot_tau[3 * i1 + 2];
    }
  }
  hphi[0] = 0.0;
  hphi[3] = -rotr[2];
  hphi[6] = rotr[1];
  hphi[1] = rotr[2];
  hphi[4] = 0.0;
  hphi[7] = -rotr[0];
  hphi[2] = -rotr[1];
  hphi[5] = rotr[0];
  hphi[8] = 0.0;
  tauhat[0] = 0.0;
  tauhat[3] = -Rb[2];
  tauhat[6] = Rb[1];
  tauhat[1] = Rb[2];
  tauhat[4] = 0.0;
  tauhat[7] = -Rb[0];
  tauhat[2] = -Rb[1];
  tauhat[5] = Rb[0];
  tauhat[8] = 0.0;
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] = -hphi[i];
  }
  for (i = 0; i < 3; i++) {
    d = tau0_dot_tau[i];
    d1 = tau0_dot_tau[i + 3];
    d2 = tau0_dot_tau[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Drotr[i + 3 * i1] = (d * RJdexpT_[3 * i1] + d1 * RJdexpT_[3 * i1 + 1]) +
                          d2 * RJdexpT_[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] = -tauhat[i];
  }
  for (i = 0; i < 3; i++) {
    d = tau0_dot_tau[i];
    d1 = tau0_dot_tau[i + 3];
    d2 = tau0_dot_tau[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Dtau[i + 3 * i1] = (d * RJdexpT_[3 * i1] + d1 * RJdexpT_[3 * i1 + 1]) +
                         d2 * RJdexpT_[3 * i1 + 2];
    }
  }
  memset(&J[0], 0, 49U * sizeof(double));
  /*  second derivative needed */
  /*  ws1 is short form for w(phi;dphi1) */
  hphi2[0] = 0.0;
  hphi2[4] = 0.0;
  hphi2[8] = 0.0;
  /*  Dws1 is short form for Dw(phi;dphi1) */
  scale = 0.0;
  d = v[0];
  d1 = v[1];
  d2 = v[2];
  d3 = a[2];
  d4 = a[3] * (rho1[4] * qbar[5] - qbar[4] * rho1[5]);
  d5 = a[3] * (qbar[3] * rho1[5] - rho1[3] * qbar[5]);
  d6 = a[3] * (rho1[3] * qbar[4] - qbar[3] * rho1[4]);
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    absxk = qbar[i + 3];
    scale += absxk * rho1[i + 3];
    d7 = RJ[3 * i];
    J[7 * i] = d7;
    J[7 * i + 3] = d7;
    d7 = Drotr[3 * i];
    J[b_i] = d7;
    J[b_i + 3] = d7 + qbar[6] * Dtau[3 * i];
    Drotr[3 * i] = qbar[3] * absxk;
    b_a[3 * i] = d3 * rho1[3] * absxk;
    Dws1[3 * i] = d4 * absxk;
    i1 = 3 * i + 1;
    d7 = RJ[i1];
    J[7 * i + 1] = d7;
    J[7 * i + 4] = d7;
    d7 = Drotr[i1];
    J[b_i + 1] = d7;
    J[b_i + 4] = d7 + qbar[6] * Dtau[i1];
    Drotr[i1] = qbar[4] * absxk;
    b_a[i1] = d3 * rho1[4] * absxk;
    Dws1[i1] = d5 * absxk;
    Qtilde_tmp = 3 * i + 2;
    d7 = RJ[Qtilde_tmp];
    J[7 * i + 2] = d7;
    J[7 * i + 5] = d7;
    d7 = Drotr[Qtilde_tmp];
    J[b_i + 2] = d7;
    J[b_i + 5] = d7 + qbar[6] * Dtau[Qtilde_tmp];
    J[b_i + 6] =
        (d * RJdexpT_[3 * i] + d1 * RJdexpT_[i1]) + d2 * RJdexpT_[Qtilde_tmp];
    Drotr[Qtilde_tmp] = qbar[5] * absxk;
    b_a[Qtilde_tmp] = d3 * rho1[5] * absxk;
    Dws1[Qtilde_tmp] = d6 * absxk;
    Thet_tilde_e2[i] = (RJdexpT_[i] * rho1[3] + RJdexpT_[i + 3] * rho1[4]) +
                       RJdexpT_[i + 6] * rho1[5];
    J[i + 45] = Rb[i];
  }
  hphi2[3] = -Thet_tilde_e2[2];
  hphi2[6] = Thet_tilde_e2[1];
  hphi2[1] = Thet_tilde_e2[2];
  hphi2[7] = -Thet_tilde_e2[0];
  hphi2[2] = -Thet_tilde_e2[1];
  hphi2[5] = Thet_tilde_e2[0];
  c_a = a[9] * scale;
  d_a = a[8] * scale;
  dv[0] = 0.0;
  dv[3] = a[1] * -rho1[5];
  dv[6] = a[1] * rho1[4];
  dv[1] = a[1] * rho1[5];
  dv[4] = 0.0;
  dv[7] = a[1] * -rho1[3];
  dv[2] = a[1] * -rho1[4];
  dv[5] = a[1] * rho1[3];
  dv[8] = 0.0;
  d = qbar[3];
  d1 = qbar[4];
  d2 = qbar[5];
  d3 = a[8];
  for (i = 0; i < 3; i++) {
    absxk = rho1[i + 3];
    e_a[3 * i] =
        ((((b_a[3 * i] - Dws1[3 * i]) - dv[3 * i]) + c_a * Drotr[3 * i]) +
         d * d3 * absxk) +
        d_a * (double)iv[3 * i];
    b_i = 3 * i + 1;
    e_a[b_i] = ((((b_a[b_i] - Dws1[b_i]) - dv[b_i]) + c_a * Drotr[b_i]) +
                d1 * d3 * absxk) +
               d_a * (double)iv[b_i];
    b_i = 3 * i + 2;
    e_a[b_i] = ((((b_a[b_i] - Dws1[b_i]) - dv[b_i]) + c_a * Drotr[b_i]) +
                d2 * d3 * absxk) +
               d_a * (double)iv[b_i];
  }
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Dws1[i + 3 * i1] =
          (d * e_a[3 * i1] + d1 * e_a[3 * i1 + 1]) + d2 * e_a[3 * i1 + 2];
    }
  }
  /*  D2rotr is short form for D^2rotr(phi)(dphi1,dot) */
  for (i = 0; i < 3; i++) {
    d = hphi2[i];
    d1 = hphi2[i + 3];
    d2 = hphi2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Qtilde_tmp = 3 * i1 + 1;
      b_Qtilde_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      D2rotr[b_i] = (hphi[i] * Dws1[3 * i1] + hphi[i + 3] * Dws1[Qtilde_tmp]) +
                    hphi[i + 6] * Dws1[b_Qtilde_tmp];
      e_a[b_i] =
          (d * hphi[3 * i1] + d1 * hphi[Qtilde_tmp]) + d2 * hphi[b_Qtilde_tmp];
    }
    d = e_a[i];
    d1 = e_a[i + 3];
    d2 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      tau0_dot_tau[i + 3 * i1] =
          (d * RJdexpT_[3 * i1] + d1 * RJdexpT_[3 * i1 + 1]) +
          d2 * RJdexpT_[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    D2rotr[i] = -(D2rotr[i] + tau0_dot_tau[i]);
  }
  /*  D2tau1 is short form for D^2tau(phi)(dphi1,dot) */
  tau0_cross_tau[0] /= nu_;
  tau0_cross_tau[1] /= nu_;
  tau0_cross_tau[2] /= nu_;
  /*  unit vector in the direction of u */
  c_a = -(2.0 / c_tmp);
  for (i = 0; i < 3; i++) {
    tau0_dot_tau[3 * i] =
        (double)iv[3 * i] - tau0_cross_tau[0] * tau0_cross_tau[i];
    b_i = 3 * i + 1;
    tau0_dot_tau[b_i] = (double)iv[b_i] - tau0_cross_tau[1] * tau0_cross_tau[i];
    b_i = 3 * i + 2;
    tau0_dot_tau[b_i] = (double)iv[b_i] - tau0_cross_tau[2] * tau0_cross_tau[i];
  }
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] = c_a * (2.0 * tau0_dot_tau[i] - (double)iv[i]);
  }
  for (i = 0; i < 3; i++) {
    d = tau0_dot_tau[i];
    d1 = tau0_dot_tau[i + 3];
    d2 = tau0_dot_tau[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_a[i + 3 * i1] = (d * tauhat[3 * i1] + d1 * tauhat[3 * i1 + 1]) +
                        d2 * tauhat[3 * i1 + 2];
    }
    d = b_a[i];
    d1 = b_a[i + 3];
    d2 = b_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      hphi[i + 3 * i1] = (d * RJdexpT_[3 * i1] + d1 * RJdexpT_[3 * i1 + 1]) +
                         d2 * RJdexpT_[3 * i1 + 2];
    }
  }
  memset(&Qtilde[0], 0, 49U * sizeof(double));
  for (i = 0; i < 3; i++) {
    d = hphi2[i];
    d1 = hphi2[i + 3];
    d2 = hphi2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Qtilde[i1 + 7 * (i + 3)] = D2rotr[i1 + 3 * i];
      Qtilde_tmp = 3 * i1 + 1;
      b_Qtilde_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      tau0_dot_tau[b_i] =
          (tauhat[i] * Dws1[3 * i1] + tauhat[i + 3] * Dws1[Qtilde_tmp]) +
          tauhat[i + 6] * Dws1[b_Qtilde_tmp];
      e_a[b_i] = (d * tauhat[3 * i1] + d1 * tauhat[Qtilde_tmp]) +
                 d2 * tauhat[b_Qtilde_tmp];
    }
    d = e_a[i];
    d1 = e_a[i + 3];
    d2 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      hphi2[i + 3 * i1] = (d * RJdexpT_[3 * i1] + d1 * RJdexpT_[3 * i1 + 1]) +
                          d2 * RJdexpT_[3 * i1 + 2];
    }
  }
  d = qbar[6];
  d1 = rho1[6];
  d2 = rho1[3];
  d3 = rho1[4];
  d4 = rho1[5];
  d5 = Thet_tilde_e2[0];
  d6 = Thet_tilde_e2[1];
  d7 = Thet_tilde_e2[2];
  d8 = v[0];
  d9 = v[1];
  d10 = v[2];
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    Qtilde[b_i + 3] =
        (D2rotr[3 * i] + d * -(tau0_dot_tau[3 * i] + hphi2[3 * i])) +
        d1 * Dtau[3 * i];
    b_Qtilde_tmp = 3 * i + 1;
    Qtilde[b_i + 4] = (D2rotr[b_Qtilde_tmp] + d * -(tau0_dot_tau[b_Qtilde_tmp] +
                                                    hphi2[b_Qtilde_tmp])) +
                      d1 * Dtau[b_Qtilde_tmp];
    Qtilde_tmp = 3 * i + 2;
    Qtilde[b_i + 5] = (D2rotr[Qtilde_tmp] +
                       d * -(tau0_dot_tau[Qtilde_tmp] + hphi2[Qtilde_tmp])) +
                      d1 * Dtau[Qtilde_tmp];
    Qtilde[b_i + 6] =
        ((d5 * hphi[3 * i] + d6 * hphi[b_Qtilde_tmp]) + d7 * hphi[Qtilde_tmp]) +
        ((d8 * Dws1[3 * i] + d9 * Dws1[b_Qtilde_tmp]) + d10 * Dws1[Qtilde_tmp]);
    Qtilde[i + 45] = (Dtau[i] * d2 + Dtau[i + 3] * d3) + Dtau[i + 6] * d4;
    tau0_cross_tau[i] = eta_data[i + 3];
  }
  Thet_tilde_e2[0] = eta_data[0] + tau0_cross_tau[0];
  Thet_tilde_e2[1] = eta_data[1] + tau0_cross_tau[1];
  Thet_tilde_e2[2] = eta_data[2] + tau0_cross_tau[2];
  absxk = ((rotr[1] * Thet_tilde_e2[2] - Thet_tilde_e2[1] * rotr[2]) +
           qbar[6] * (Rb[1] * eta_data[5] - Rb[2] * eta_data[4])) +
          v[0] * eta_data[6];
  scale = ((Thet_tilde_e2[0] * rotr[2] - rotr[0] * Thet_tilde_e2[2]) +
           qbar[6] * (Rb[2] * eta_data[3] - Rb[0] * eta_data[5])) +
          v[1] * eta_data[6];
  t = ((rotr[0] * Thet_tilde_e2[1] - Thet_tilde_e2[0] * rotr[1]) +
       qbar[6] * (Rb[0] * eta_data[4] - Rb[1] * eta_data[3])) +
      v[2] * eta_data[6];
  c_tmp = 0.0;
  for (i = 0; i < 3; i++) {
    d = (RJ[3 * i] * absxk + RJ[3 * i + 1] * scale) + RJ[3 * i + 2] * t;
    v[i] = d;
    c_tmp += d * qbar[i + 3];
  }
  c_a = a[9] * c_tmp;
  d_a = a[8] * c_tmp;
  memset(&Q[0], 0, 49U * sizeof(double));
  absxk = ((rotr[0] * Thet_tilde_e2[0] + rotr[1] * Thet_tilde_e2[1]) +
           rotr[2] * Thet_tilde_e2[2]) +
          qbar[6] * ((Rb[0] * tau0_cross_tau[0] + Rb[1] * tau0_cross_tau[1]) +
                     Rb[2] * tau0_cross_tau[2]);
  dv[0] = 0.0;
  dv[3] = a[1] * -v[2];
  dv[6] = a[1] * v[1];
  dv[1] = a[1] * v[2];
  dv[4] = 0.0;
  dv[7] = a[1] * -v[0];
  dv[2] = a[1] * -v[1];
  dv[5] = v[0] * a[1];
  dv[8] = 0.0;
  d = v[0];
  d1 = v[1];
  d2 = v[2];
  d3 = a[2];
  d4 = a[3] * (v[2] * qbar[4] - v[1] * qbar[5]);
  d5 = a[3] * (v[0] * qbar[5] - v[2] * qbar[3]);
  d6 = a[3] * (v[1] * qbar[3] - v[0] * qbar[4]);
  d7 = rotr[0];
  d8 = rotr[1];
  d9 = rotr[2];
  d10 = p2_tmp[0];
  d11 = p2_tmp[1];
  d12 = p2_tmp[2];
  for (i = 0; i < 3; i++) {
    t = qbar[i + 3];
    b_a[3 * i] = d * d3 * t;
    Dws1[3 * i] = d4 * t;
    c_tmp = Thet_tilde_e2[i];
    e_a[3 * i] = d7 * c_tmp;
    scale = tau0_cross_tau[i];
    tau0_dot_tau[3 * i] = d10 * scale;
    b_i = 3 * i + 1;
    b_a[b_i] = d1 * d3 * t;
    Dws1[b_i] = d5 * t;
    e_a[b_i] = d8 * c_tmp;
    tau0_dot_tau[b_i] = d11 * scale;
    b_i = 3 * i + 2;
    b_a[b_i] = d3 * d2 * t;
    Dws1[b_i] = d6 * t;
    e_a[b_i] = d9 * c_tmp;
    tau0_dot_tau[b_i] = d12 * scale;
  }
  for (i = 0; i < 9; i++) {
    e_a[i] = (e_a[i] + tau0_dot_tau[i]) - absxk * (double)iv[i];
  }
  for (i = 0; i < 3; i++) {
    d = e_a[i];
    d1 = e_a[i + 3];
    d2 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      D2rotr[b_i] = ((d * RJdexpT_[3 * i1] + d1 * RJdexpT_[3 * i1 + 1]) +
                     d2 * RJdexpT_[3 * i1 + 2]) +
                    eta_data[6] * hphi[b_i];
    }
  }
  for (i = 0; i < 3; i++) {
    d = RJdexpT_[3 * i];
    d1 = RJdexpT_[3 * i + 1];
    d2 = RJdexpT_[3 * i + 2];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i1 + 3 * i;
      e_a[b_i] = ((((b_a[b_i] - Dws1[b_i]) + dv[b_i]) + c_a * Drotr[b_i]) +
                  d_a * (double)iv[b_i]) +
                 a[8] * qbar[i1 + 3] * v[i];
      tau0_dot_tau[i + 3 * i1] =
          (d * D2rotr[3 * i1] + d1 * D2rotr[3 * i1 + 1]) +
          d2 * D2rotr[3 * i1 + 2];
    }
  }
  d = tau0_cross_tau[0];
  d1 = tau0_cross_tau[1];
  d2 = tau0_cross_tau[2];
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    Q[b_i + 3] = e_a[3 * i] + tau0_dot_tau[3 * i];
    b_Qtilde_tmp = 3 * i + 1;
    Q[b_i + 4] = e_a[b_Qtilde_tmp] + tau0_dot_tau[b_Qtilde_tmp];
    Qtilde_tmp = 3 * i + 2;
    Q[b_i + 5] = e_a[Qtilde_tmp] + tau0_dot_tau[Qtilde_tmp];
    d3 = (d * Dtau[3 * i] + d1 * Dtau[b_Qtilde_tmp]) + d2 * Dtau[Qtilde_tmp];
    Q[b_i + 6] = d3;
    Q[i + 45] = d3;
  }
}

void c_CableBCTransinCoord(const double qbar[7], const double x0[3],
                           const double RJ[9], const double RE[9],
                           const double r[3], const double R0[9],
                           const double eta[7], const double rho1[7],
                           const double rho2[7], double q[7], double J[49],
                           double Q[49], double Qtilde[49], double C[49])
{
  double a[11];
  double D2rotr_tmp[9];
  double D2tau1[9];
  double D2tau1_tmp[9];
  double D2ws1[9];
  double Drotr[9];
  double Dtau[9];
  double Dv[9];
  double Dv_tmp[9];
  double Dws1[9];
  double Dws2[9];
  double PP[9];
  double RJdexpT_[9];
  double Rb[9];
  double a_tmp[9];
  double b_a[9];
  double b_rotr[9];
  double c_a[9];
  double dv[9];
  double e_a[9];
  double hphi[9];
  double hphi2[9];
  double rotrhat[9];
  double tau0_dot_tau[9];
  double tauhat[9];
  double ws1hat[9];
  double ws2hat[9];
  double Thet_tilde_e2[3];
  double rotr[3];
  double tau0_cross_tau[3];
  double u[3];
  double v[3];
  double ws1[3];
  double ws2[3];
  double absxk;
  double b_a_tmp;
  double b_qbar;
  double c_a_tmp;
  double d;
  double d1;
  double d10;
  double d11;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double d_a;
  double d_a_tmp;
  double dphi1xphi_idx_0;
  double dphi1xphi_idx_1;
  double dphi1xphi_idx_2;
  double f_a;
  double g_a;
  double h_a;
  double nu_;
  double scale;
  double t;
  double tau0_dot_tau_tmp;
  int J_tmp;
  int Q_tmp;
  int b_i;
  int i;
  int i1;
  /*  Computes the transformation from (d, phi, gamma) -> (p1, p2, vartheta) for
   */
  /*  one end of the cable; NOTE: the main difference from CableBCtrans.m is */
  /*  that this function is based on exponential coordinates for the rotation,
   */
  /*  rather than the rotation matrix; also, a joint coordinate system and end
   */
  /*  offset are included, making it a more general. */
  /*  Inputs */
  /*  qbar = [d; phi; gamm] */
  /*    where d = displacement of joint,  */
  /*          phi = exponential coordinates of rotation of joint */
  /*          gamm = distance between first and second control points */
  /*  x0 = reference position of joint (3x1) vector */
  /*  RJ = rotation of joint coordinate frame with respect to global */
  /*  RE = rotation of cable end with respect to joint coordinate frame */
  /*  r = position of cable end relative to joint in joint coordinate system */
  /*      (end offset, 3x1 vector) */
  /*  R0 = orientation of cable end in reference configuration */
  /*  eta = vector to contract with for second derivative, Q */
  /*  rho1 = vector to contract with for second derivative, Qtilde */
  /*  rho2 = vector to contract with for third derivative, C */
  /*  Outputs: */
  /*  q = [p1; p2; vartheta1] */
  /*    where p1, p2 = positions of first and second control points */
  /*          vartheta1 = (twist) rotation of first control point */
  /*  J = first derivative of transformation */
  /*  Q = second derivative contracted over upper index,  */
  /*      Q_{ib,jb} = eta_i q^i_{ib,jb} */
  /*  Qtilde = second derivative contracted over one of the lower indices */
  /*      Qtilde^{i}_{ib} = rho1^{jb} q^i_{ib,jb} */
  /*  C = third derivative contracted over two lower indices */
  /*      C^{i}_{ib} = q^i_{ib,jb,kb}rho1^{jb}rho2^{kb} */
  b_rothelper(&qbar[3], a);
  hphi[0] = 0.0;
  hphi[3] = -qbar[5];
  hphi[6] = qbar[4];
  hphi[1] = qbar[5];
  hphi[4] = 0.0;
  hphi[7] = -qbar[3];
  hphi[2] = -qbar[4];
  hphi[5] = qbar[3];
  hphi[8] = 0.0;
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      hphi2[i + 3 * i1] =
          (hphi[i] * hphi[3 * i1] + hphi[i + 3] * hphi[3 * i1 + 1]) +
          hphi[i + 6] * hphi[3 * i1 + 2];
    }
  }
  d = a[0];
  d1 = a[1];
  for (i = 0; i < 9; i++) {
    b_a[i] = (d * hphi[i] + (double)iv[i]) + d1 * hphi2[i];
  }
  scale = 0.0;
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Drotr[i + 3 * i1] =
          (d * b_a[3 * i1] + d1 * b_a[3 * i1 + 1]) + d2 * b_a[3 * i1 + 2];
    }
    d = Drotr[i];
    d1 = Drotr[i + 3];
    d2 = Drotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Rb[i + 3 * i1] =
          (d * RE[3 * i1] + d1 * RE[3 * i1 + 1]) + d2 * RE[3 * i1 + 2];
    }
    scale += R0[i] * Rb[i];
  }
  tau0_cross_tau[0] = R0[1] * Rb[2] - Rb[1] * R0[2];
  tau0_cross_tau[1] = Rb[0] * R0[2] - R0[0] * Rb[2];
  tau0_cross_tau[2] = R0[0] * Rb[1] - Rb[0] * R0[1];
  dv[0] = 0.0;
  dv[1] = -tau0_cross_tau[2];
  dv[2] = tau0_cross_tau[1];
  dv[3] = tau0_cross_tau[2];
  dv[4] = 0.0;
  dv[5] = -tau0_cross_tau[0];
  dv[6] = -tau0_cross_tau[1];
  dv[7] = tau0_cross_tau[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d = tau0_cross_tau[i] / (scale + 1.0);
    tau0_dot_tau[3 * i] =
        (scale * (double)iv[i] + dv[3 * i]) + d * tau0_cross_tau[0];
    b_i = 3 * i + 1;
    tau0_dot_tau[b_i] =
        (scale * (double)iv[i + 3] + dv[b_i]) + d * tau0_cross_tau[1];
    b_i = 3 * i + 2;
    tau0_dot_tau[b_i] =
        (scale * (double)iv[i + 6] + dv[b_i]) + d * tau0_cross_tau[2];
  }
  d = Rb[3];
  d1 = Rb[4];
  d2 = Rb[5];
  for (i = 0; i < 3; i++) {
    ws2[i] = (tau0_dot_tau[i] * d + tau0_dot_tau[i + 3] * d1) +
             tau0_dot_tau[i + 6] * d2;
  }
  d = r[0];
  d1 = r[1];
  d2 = r[2];
  d3 = ws2[0];
  d4 = ws2[1];
  d5 = ws2[2];
  d6 = qbar[0];
  d7 = qbar[1];
  d8 = qbar[2];
  d9 = qbar[6];
  for (i = 0; i < 3; i++) {
    d10 = (Drotr[i] * d + Drotr[i + 3] * d1) + Drotr[i + 6] * d2;
    rotr[i] = d10;
    Thet_tilde_e2[i] =
        (R0[3 * i] * d3 + R0[3 * i + 1] * d4) + R0[3 * i + 2] * d5;
    d10 += x0[i] + ((RJ[i] * d6 + RJ[i + 3] * d7) + RJ[i + 6] * d8);
    d11 = d9 * Rb[i];
    ws2[i] = d11;
    q[i] = d10;
    q[i + 3] = d10 + d11;
  }
  q[6] = atan2(Thet_tilde_e2[2], Thet_tilde_e2[1]);
  /*  first derivative needed */
  scale = 3.3121686421112381E-170;
  d = Rb[0] + R0[0];
  u[0] = d;
  absxk = fabs(d);
  if (absxk > 3.3121686421112381E-170) {
    nu_ = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    nu_ = t * t;
  }
  d = Rb[1] + R0[1];
  u[1] = d;
  absxk = fabs(d);
  if (absxk > scale) {
    t = scale / absxk;
    nu_ = nu_ * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nu_ += t * t;
  }
  d = Rb[2] + R0[2];
  u[2] = d;
  absxk = fabs(d);
  if (absxk > scale) {
    t = scale / absxk;
    nu_ = nu_ * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nu_ += t * t;
  }
  nu_ = scale * sqrt(nu_);
  absxk = nu_ * nu_;
  d1 = a[1];
  d2 = a[8];
  for (b_i = 0; b_i < 3; b_i++) {
    v[b_i] = 2.0 * u[b_i] / absxk;
    tau0_dot_tau[3 * b_i] =
        ((double)iv[b_i] - d1 * hphi[b_i]) + d2 * hphi2[b_i];
    tau0_dot_tau[3 * b_i + 1] =
        ((double)iv[b_i + 3] - d1 * hphi[b_i + 3]) + d2 * hphi2[b_i + 3];
    tau0_dot_tau[3 * b_i + 2] =
        ((double)iv[b_i + 6] - d1 * hphi[b_i + 6]) + d2 * hphi2[b_i + 6];
  }
  for (i = 0; i < 3; i++) {
    d1 = RJ[i];
    d2 = RJ[i + 3];
    d3 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      RJdexpT_[i + 3 * i1] =
          (d1 * tau0_dot_tau[3 * i1] + d2 * tau0_dot_tau[3 * i1 + 1]) +
          d3 * tau0_dot_tau[3 * i1 + 2];
    }
  }
  rotrhat[0] = 0.0;
  rotrhat[3] = -rotr[2];
  rotrhat[6] = rotr[1];
  rotrhat[1] = rotr[2];
  rotrhat[4] = 0.0;
  rotrhat[7] = -rotr[0];
  rotrhat[2] = -rotr[1];
  rotrhat[5] = rotr[0];
  rotrhat[8] = 0.0;
  tauhat[0] = 0.0;
  tauhat[3] = -Rb[2];
  tauhat[6] = Rb[1];
  tauhat[1] = Rb[2];
  tauhat[4] = 0.0;
  tauhat[7] = -Rb[0];
  tauhat[2] = -Rb[1];
  tauhat[5] = Rb[0];
  tauhat[8] = 0.0;
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] = -rotrhat[i];
  }
  for (i = 0; i < 3; i++) {
    d1 = tau0_dot_tau[i];
    d2 = tau0_dot_tau[i + 3];
    d3 = tau0_dot_tau[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Drotr[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                          d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    a_tmp[i] = -tauhat[i];
  }
  for (i = 0; i < 3; i++) {
    d1 = a_tmp[i];
    d2 = a_tmp[i + 3];
    d3 = a_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Dtau[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                         d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  memset(&J[0], 0, 49U * sizeof(double));
  /*  second derivative needed */
  /*  ws1 is short form for w(phi;dphi1) */
  ws1hat[0] = 0.0;
  ws1hat[4] = 0.0;
  ws1hat[8] = 0.0;
  /*  Dws1 is short form for Dw(phi;dphi1) */
  b_qbar = 0.0;
  dphi1xphi_idx_0 = rho1[4] * qbar[5] - qbar[4] * rho1[5];
  dphi1xphi_idx_1 = qbar[3] * rho1[5] - rho1[3] * qbar[5];
  dphi1xphi_idx_2 = rho1[3] * qbar[4] - qbar[3] * rho1[4];
  hphi2[0] = 0.0;
  hphi2[3] = -rho1[5];
  hphi2[6] = rho1[4];
  hphi2[1] = rho1[5];
  hphi2[4] = 0.0;
  hphi2[7] = -rho1[3];
  hphi2[2] = -rho1[4];
  hphi2[5] = rho1[3];
  hphi2[8] = 0.0;
  d1 = v[0];
  d2 = v[1];
  d3 = v[2];
  d4 = a[2];
  d5 = a[3];
  for (i = 0; i < 3; i++) {
    J_tmp = 7 * (i + 3);
    scale = qbar[i + 3];
    b_qbar += scale * rho1[i + 3];
    d6 = RJ[3 * i];
    J[7 * i] = d6;
    J[7 * i + 3] = d6;
    d6 = Drotr[3 * i];
    J[J_tmp] = d6;
    J[J_tmp + 3] = d6 + qbar[6] * Dtau[3 * i];
    Drotr[3 * i] = qbar[3] * scale;
    b_a[3 * i] = d4 * rho1[3] * scale;
    c_a[3 * i] = d5 * dphi1xphi_idx_0 * scale;
    i1 = 3 * i + 1;
    d6 = RJ[i1];
    J[7 * i + 1] = d6;
    J[7 * i + 4] = d6;
    d6 = Drotr[i1];
    J[J_tmp + 1] = d6;
    J[J_tmp + 4] = d6 + qbar[6] * Dtau[i1];
    Drotr[i1] = qbar[4] * scale;
    b_a[i1] = d4 * rho1[4] * scale;
    c_a[i1] = d5 * dphi1xphi_idx_1 * scale;
    Q_tmp = 3 * i + 2;
    d6 = RJ[Q_tmp];
    J[7 * i + 2] = d6;
    J[7 * i + 5] = d6;
    d6 = Drotr[Q_tmp];
    J[J_tmp + 2] = d6;
    J[J_tmp + 5] = d6 + qbar[6] * Dtau[Q_tmp];
    J[J_tmp + 6] =
        (d1 * RJdexpT_[3 * i] + d2 * RJdexpT_[i1]) + d3 * RJdexpT_[Q_tmp];
    Drotr[Q_tmp] = qbar[5] * scale;
    b_a[Q_tmp] = d4 * rho1[5] * scale;
    c_a[Q_tmp] = d5 * dphi1xphi_idx_2 * scale;
    ws1[i] = (RJdexpT_[i] * rho1[3] + RJdexpT_[i + 3] * rho1[4]) +
             RJdexpT_[i + 6] * rho1[5];
    J[i + 45] = Rb[i];
  }
  ws1hat[3] = -ws1[2];
  ws1hat[6] = ws1[1];
  ws1hat[1] = ws1[2];
  ws1hat[7] = -ws1[0];
  ws1hat[2] = -ws1[1];
  ws1hat[5] = ws1[0];
  b_a_tmp = a[9] * b_qbar;
  d_a = a[8] * b_qbar;
  d1 = a[1];
  d2 = a[8];
  d3 = qbar[3];
  d4 = qbar[4];
  d5 = qbar[5];
  for (i = 0; i < 3; i++) {
    c_a_tmp = rho1[i + 3];
    e_a[3 * i] = ((((b_a[3 * i] - c_a[3 * i]) - d1 * hphi2[3 * i]) +
                   b_a_tmp * Drotr[3 * i]) +
                  d3 * d2 * c_a_tmp) +
                 d_a * (double)iv[3 * i];
    b_i = 3 * i + 1;
    e_a[b_i] =
        ((((b_a[b_i] - c_a[b_i]) - d1 * hphi2[b_i]) + b_a_tmp * Drotr[b_i]) +
         d4 * d2 * c_a_tmp) +
        d_a * (double)iv[b_i];
    b_i = 3 * i + 2;
    e_a[b_i] =
        ((((b_a[b_i] - c_a[b_i]) - d1 * hphi2[b_i]) + b_a_tmp * Drotr[b_i]) +
         d5 * d2 * c_a_tmp) +
        d_a * (double)iv[b_i];
  }
  /*  D2rotr is short form for D^2rotr(phi)(dphi1,dot) */
  for (i = 0; i < 3; i++) {
    d1 = ws1hat[i];
    d2 = ws1hat[i + 3];
    d3 = ws1hat[i + 6];
    d4 = RJ[i];
    d5 = RJ[i + 3];
    d6 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      D2rotr_tmp[b_i] =
          (d1 * rotrhat[3 * i1] + d2 * rotrhat[Q_tmp]) + d3 * rotrhat[J_tmp];
      Dws1[b_i] = (d4 * e_a[3 * i1] + d5 * e_a[Q_tmp]) + d6 * e_a[J_tmp];
    }
  }
  for (i = 0; i < 3; i++) {
    d1 = D2rotr_tmp[i];
    d2 = D2rotr_tmp[i + 3];
    d3 = D2rotr_tmp[i + 6];
    d4 = rotrhat[i];
    d5 = rotrhat[i + 3];
    d6 = rotrhat[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      tau0_dot_tau[b_i] =
          (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[Q_tmp]) + d3 * RJdexpT_[J_tmp];
      hphi[b_i] = (d4 * Dws1[3 * i1] + d5 * Dws1[Q_tmp]) + d6 * Dws1[J_tmp];
    }
  }
  for (i = 0; i < 9; i++) {
    hphi[i] = -(hphi[i] + tau0_dot_tau[i]);
  }
  /*  D2tau1 is short form for D^2tau(phi)(dphi1,dot) */
  for (i = 0; i < 3; i++) {
    d1 = ws1hat[i];
    d2 = ws1hat[i + 3];
    d3 = ws1hat[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      D2tau1[b_i] = (tauhat[i] * Dws1[3 * i1] + tauhat[i + 3] * Dws1[Q_tmp]) +
                    tauhat[i + 6] * Dws1[J_tmp];
      D2tau1_tmp[b_i] =
          (d1 * tauhat[3 * i1] + d2 * tauhat[Q_tmp]) + d3 * tauhat[J_tmp];
    }
    d1 = D2tau1_tmp[i];
    d2 = D2tau1_tmp[i + 3];
    d3 = D2tau1_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      tau0_dot_tau[i + 3 * i1] =
          (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
          d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    D2tau1[i] = -(D2tau1[i] + tau0_dot_tau[i]);
  }
  tau0_cross_tau[0] = u[0] / nu_;
  tau0_cross_tau[1] = u[1] / nu_;
  tau0_cross_tau[2] = d / nu_;
  /*  unit vector in the direction of u */
  for (i = 0; i < 3; i++) {
    PP[3 * i] = (double)iv[3 * i] - tau0_cross_tau[0] * tau0_cross_tau[i];
    b_i = 3 * i + 1;
    PP[b_i] = (double)iv[b_i] - tau0_cross_tau[1] * tau0_cross_tau[i];
    b_i = 3 * i + 2;
    PP[b_i] = (double)iv[b_i] - tau0_cross_tau[2] * tau0_cross_tau[i];
  }
  for (i = 0; i < 9; i++) {
    Dv_tmp[i] = 2.0 * PP[i] - (double)iv[i];
  }
  c_a_tmp = 2.0 / absxk;
  for (i = 0; i < 3; i++) {
    d1 = Dv_tmp[i];
    d2 = Dv_tmp[i + 3];
    d3 = Dv_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_a[i + 3 * i1] = (-c_a_tmp * d1 * tauhat[3 * i1] +
                         -c_a_tmp * d2 * tauhat[3 * i1 + 1]) +
                        -c_a_tmp * d3 * tauhat[3 * i1 + 2];
    }
    d1 = b_a[i];
    d2 = b_a[i + 3];
    d3 = b_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Dv[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                       d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  memset(&Qtilde[0], 0, 49U * sizeof(double));
  d1 = qbar[6];
  d2 = rho1[6];
  d3 = rho1[3];
  d4 = rho1[4];
  d5 = rho1[5];
  d6 = ws1[0];
  d7 = ws1[1];
  d8 = ws1[2];
  d9 = v[0];
  d10 = v[1];
  d11 = v[2];
  for (b_i = 0; b_i < 3; b_i++) {
    scale = hphi[3 * b_i];
    J_tmp = 7 * (b_i + 3);
    Qtilde[J_tmp] = scale;
    Qtilde[J_tmp + 3] = (scale + d1 * D2tau1[3 * b_i]) + d2 * Dtau[3 * b_i];
    i = 3 * b_i + 1;
    scale = hphi[i];
    Qtilde[J_tmp + 1] = scale;
    Qtilde[J_tmp + 4] = (scale + d1 * D2tau1[i]) + d2 * Dtau[i];
    i1 = 3 * b_i + 2;
    scale = hphi[i1];
    Qtilde[J_tmp + 2] = scale;
    Qtilde[J_tmp + 5] = (scale + d1 * D2tau1[i1]) + d2 * Dtau[i1];
    Qtilde[J_tmp + 6] = ((d6 * Dv[3 * b_i] + d7 * Dv[i]) + d8 * Dv[i1]) +
                        ((d9 * Dws1[3 * b_i] + d10 * Dws1[i]) + d11 * Dws1[i1]);
    Qtilde[b_i + 45] =
        (Dtau[b_i] * d3 + Dtau[b_i + 3] * d4) + Dtau[b_i + 6] * d5;
    tau0_cross_tau[b_i] = eta[b_i] + eta[b_i + 3];
  }
  scale = 0.0;
  d1 = ((rotr[1] * tau0_cross_tau[2] - tau0_cross_tau[1] * rotr[2]) +
        qbar[6] * (Rb[1] * eta[5] - Rb[2] * eta[4])) +
       v[0] * eta[6];
  d2 = ((tau0_cross_tau[0] * rotr[2] - rotr[0] * tau0_cross_tau[2]) +
        qbar[6] * (Rb[2] * eta[3] - Rb[0] * eta[5])) +
       v[1] * eta[6];
  d3 = ((rotr[0] * tau0_cross_tau[1] - tau0_cross_tau[0] * rotr[1]) +
        qbar[6] * (Rb[0] * eta[4] - Rb[1] * eta[3])) +
       v[2] * eta[6];
  for (i = 0; i < 3; i++) {
    d4 = (RJ[3 * i] * d1 + RJ[3 * i + 1] * d2) + RJ[3 * i + 2] * d3;
    Thet_tilde_e2[i] = d4;
    scale += d4 * qbar[i + 3];
  }
  d_a = a[9] * scale;
  f_a = a[8] * scale;
  memset(&Q[0], 0, 49U * sizeof(double));
  g_a = ((rotr[0] * tau0_cross_tau[0] + rotr[1] * tau0_cross_tau[1]) +
         rotr[2] * tau0_cross_tau[2]) +
        qbar[6] * ((Rb[0] * eta[3] + Rb[1] * eta[4]) + Rb[2] * eta[5]);
  dv[0] = 0.0;
  dv[3] = a[1] * -Thet_tilde_e2[2];
  dv[6] = a[1] * Thet_tilde_e2[1];
  dv[1] = a[1] * Thet_tilde_e2[2];
  dv[4] = 0.0;
  dv[7] = a[1] * -Thet_tilde_e2[0];
  dv[2] = a[1] * -Thet_tilde_e2[1];
  dv[5] = Thet_tilde_e2[0] * a[1];
  dv[8] = 0.0;
  d1 = Thet_tilde_e2[0];
  d2 = Thet_tilde_e2[1];
  d3 = Thet_tilde_e2[2];
  d4 = a[2];
  d5 = a[3] * (Thet_tilde_e2[2] * qbar[4] - Thet_tilde_e2[1] * qbar[5]);
  d6 = a[3] * (Thet_tilde_e2[0] * qbar[5] - Thet_tilde_e2[2] * qbar[3]);
  d7 = a[3] * (Thet_tilde_e2[1] * qbar[3] - Thet_tilde_e2[0] * qbar[4]);
  d8 = rotr[0];
  d9 = rotr[1];
  d10 = rotr[2];
  d11 = ws2[0];
  scale = ws2[1];
  absxk = ws2[2];
  for (i = 0; i < 3; i++) {
    t = qbar[i + 3];
    b_a[3 * i] = d1 * d4 * t;
    c_a[3 * i] = d5 * t;
    h_a = tau0_cross_tau[i];
    e_a[3 * i] = d8 * h_a;
    tau0_dot_tau_tmp = eta[i + 3];
    tau0_dot_tau[3 * i] = d11 * tau0_dot_tau_tmp;
    b_i = 3 * i + 1;
    b_a[b_i] = d2 * d4 * t;
    c_a[b_i] = d6 * t;
    e_a[b_i] = d9 * h_a;
    tau0_dot_tau[b_i] = scale * tau0_dot_tau_tmp;
    b_i = 3 * i + 2;
    b_a[b_i] = d4 * d3 * t;
    c_a[b_i] = d7 * t;
    e_a[b_i] = d10 * h_a;
    tau0_dot_tau[b_i] = absxk * tau0_dot_tau_tmp;
  }
  for (i = 0; i < 9; i++) {
    e_a[i] = (e_a[i] + tau0_dot_tau[i]) - g_a * (double)iv[i];
  }
  for (i = 0; i < 3; i++) {
    d1 = e_a[i];
    d2 = e_a[i + 3];
    d3 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      b_rotr[b_i] = ((d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                     d3 * RJdexpT_[3 * i1 + 2]) +
                    eta[6] * Dv[b_i];
    }
  }
  for (i = 0; i < 3; i++) {
    d1 = RJdexpT_[3 * i];
    d2 = RJdexpT_[3 * i + 1];
    d3 = RJdexpT_[3 * i + 2];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i1 + 3 * i;
      e_a[b_i] = ((((b_a[b_i] - c_a[b_i]) + dv[b_i]) + d_a * Drotr[b_i]) +
                  f_a * (double)iv[b_i]) +
                 a[8] * qbar[i1 + 3] * Thet_tilde_e2[i];
      tau0_dot_tau[i + 3 * i1] =
          (d1 * b_rotr[3 * i1] + d2 * b_rotr[3 * i1 + 1]) +
          d3 * b_rotr[3 * i1 + 2];
    }
  }
  d1 = eta[3];
  d2 = eta[4];
  d3 = eta[5];
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    Q[b_i + 3] = e_a[3 * i] + tau0_dot_tau[3 * i];
    J_tmp = 3 * i + 1;
    Q[b_i + 4] = e_a[J_tmp] + tau0_dot_tau[J_tmp];
    Q_tmp = 3 * i + 2;
    Q[b_i + 5] = e_a[Q_tmp] + tau0_dot_tau[Q_tmp];
    d4 = (d1 * Dtau[3 * i] + d2 * Dtau[J_tmp]) + d3 * Dtau[Q_tmp];
    Q[b_i + 6] = d4;
    Q[i + 45] = d4;
  }
  /*  third derivtive needed */
  /*  ws2 is short form for ws(phi;dphi2) */
  ws2hat[0] = 0.0;
  ws2hat[4] = 0.0;
  ws2hat[8] = 0.0;
  /*  Dws2 is short form for Dws(phi;dphi2) */
  scale = 0.0;
  d1 = a[2];
  d2 = a[3] * (rho2[4] * qbar[5] - qbar[4] * rho2[5]);
  d3 = a[3] * (qbar[3] * rho2[5] - rho2[3] * qbar[5]);
  d4 = a[3] * (rho2[3] * qbar[4] - qbar[3] * rho2[4]);
  for (i = 0; i < 3; i++) {
    d5 = qbar[i + 3];
    scale += d5 * rho2[i + 3];
    b_a[3 * i] = d1 * rho2[3] * d5;
    c_a[3 * i] = d2 * d5;
    b_i = 3 * i + 1;
    b_a[b_i] = d1 * rho2[4] * d5;
    c_a[b_i] = d3 * d5;
    b_i = 3 * i + 2;
    b_a[b_i] = d1 * rho2[5] * d5;
    c_a[b_i] = d4 * d5;
    ws2[i] = (RJdexpT_[i] * rho2[3] + RJdexpT_[i + 3] * rho2[4]) +
             RJdexpT_[i + 6] * rho2[5];
  }
  ws2hat[3] = -ws2[2];
  ws2hat[6] = ws2[1];
  ws2hat[1] = ws2[2];
  ws2hat[7] = -ws2[0];
  ws2hat[2] = -ws2[1];
  ws2hat[5] = ws2[0];
  d_a_tmp = a[9] * scale;
  d_a = a[8] * scale;
  dv[0] = 0.0;
  dv[3] = a[1] * -rho2[5];
  dv[6] = a[1] * rho2[4];
  dv[1] = a[1] * rho2[5];
  dv[4] = 0.0;
  dv[7] = a[1] * -rho2[3];
  dv[2] = a[1] * -rho2[4];
  dv[5] = a[1] * rho2[3];
  dv[8] = 0.0;
  d1 = qbar[3];
  d2 = qbar[4];
  d3 = qbar[5];
  d4 = a[8];
  for (i = 0; i < 3; i++) {
    t = rho2[i + 3];
    e_a[3 * i] =
        ((((b_a[3 * i] - c_a[3 * i]) - dv[3 * i]) + d_a_tmp * Drotr[3 * i]) +
         d1 * d4 * t) +
        d_a * (double)iv[3 * i];
    b_i = 3 * i + 1;
    e_a[b_i] = ((((b_a[b_i] - c_a[b_i]) - dv[b_i]) + d_a_tmp * Drotr[b_i]) +
                d2 * d4 * t) +
               d_a * (double)iv[b_i];
    b_i = 3 * i + 2;
    e_a[b_i] = ((((b_a[b_i] - c_a[b_i]) - dv[b_i]) + d_a_tmp * Drotr[b_i]) +
                d3 * d4 * t) +
               d_a * (double)iv[b_i];
  }
  /*  D2ws1 is short form for D^2 ws(phi;dphi1)(dphi2,dot) */
  absxk = 0.0;
  d_a = a[4] * scale;
  f_a = a[5] * scale;
  g_a = a[3] * scale;
  tau0_dot_tau_tmp = a[10] * b_qbar * scale;
  h_a = b_a_tmp * scale;
  for (i = 0; i < 3; i++) {
    scale = rho2[i + 3];
    absxk += rho1[i + 3] * scale;
    d1 = RJ[i];
    d2 = RJ[i + 3];
    d3 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Dws2[i + 3 * i1] =
          (d1 * e_a[3 * i1] + d2 * e_a[3 * i1 + 1]) + d3 * e_a[3 * i1 + 2];
      t = rho1[i1 + 3];
      b_i = i1 + 3 * i;
      b_a[b_i] = d_a * t * qbar[i + 3];
      c_a[b_i] = a[2] * t * scale;
    }
  }
  d_a = a[9] * absxk;
  scale = a[8] * absxk;
  d1 = a[3];
  d2 = a[8];
  d3 = a[3] * (rho1[4] * rho2[5] - rho2[4] * rho1[5]);
  d4 = a[3] * (rho2[3] * rho1[5] - rho1[3] * rho2[5]);
  d5 = a[3] * (rho1[3] * rho2[4] - rho2[3] * rho1[4]);
  for (i = 0; i < 3; i++) {
    signed char i2;
    d6 = rho1[i + 3];
    d7 = Drotr[3 * i];
    d8 = qbar[i + 3];
    d9 = rho2[i + 3];
    d10 = (((((b_a[3 * i] + c_a[3 * i]) - f_a * dphi1xphi_idx_0 * d8) -
             d1 * dphi1xphi_idx_0 * d9) -
            g_a * hphi2[3 * i]) -
           d3 * d8) +
          tau0_dot_tau_tmp * d7;
    i2 = iv[3 * i];
    d11 = b_a_tmp * rho2[3] * d8;
    b_a[3 * i] = d11;
    c_a[3 * i] = ((((((d10 + d_a_tmp * qbar[3] * d6) + b_a_tmp * qbar[3] * d9) +
                     h_a * (double)i2) +
                    d_a * d7) +
                   scale * (double)i2) +
                  d11) +
                 rho2[3] * d2 * d6;
    i1 = 3 * i + 1;
    d7 = Drotr[i1];
    d10 = (((((b_a[i1] + c_a[i1]) - f_a * dphi1xphi_idx_1 * d8) -
             d1 * dphi1xphi_idx_1 * d9) -
            g_a * hphi2[i1]) -
           d4 * d8) +
          tau0_dot_tau_tmp * d7;
    i2 = iv[i1];
    d11 = b_a_tmp * rho2[4] * d8;
    b_a[i1] = d11;
    c_a[i1] = ((((((d10 + d_a_tmp * qbar[4] * d6) + b_a_tmp * qbar[4] * d9) +
                  h_a * (double)i2) +
                 d_a * d7) +
                scale * (double)i2) +
               d11) +
              rho2[4] * d2 * d6;
    i1 = 3 * i + 2;
    d7 = Drotr[i1];
    d10 = (((((b_a[i1] + c_a[i1]) - f_a * dphi1xphi_idx_2 * d8) -
             d1 * dphi1xphi_idx_2 * d9) -
            g_a * hphi2[i1]) -
           d5 * d8) +
          tau0_dot_tau_tmp * d7;
    i2 = iv[i1];
    d11 = b_a_tmp * rho2[5] * d8;
    b_a[i1] = d11;
    c_a[i1] = ((((((d10 + d_a_tmp * qbar[5] * d6) + b_a_tmp * qbar[5] * d9) +
                  h_a * (double)i2) +
                 d_a * d7) +
                scale * (double)i2) +
               d11) +
              rho2[5] * d2 * d6;
  }
  /*  D2v is short form for D^2 v(phi)(dphi2,dot) */
  absxk = 0.0;
  d_a = -(4.0 / pow(nu_, 4.0));
  /*  D3rotr is short form for D^3 rotr(phi)(dphi1,dphi2,dot) */
  for (i = 0; i < 3; i++) {
    d1 = RJ[i];
    d2 = RJ[i + 3];
    d3 = RJ[i + 6];
    d4 = 0.0;
    d5 = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      D2ws1[b_i] =
          (d1 * c_a[3 * i1] + d2 * c_a[3 * i1 + 1]) + d3 * c_a[3 * i1 + 2];
      scale = rho2[i1 + 3];
      d4 += Dtau[b_i] * scale;
      d5 += Dws1[b_i] * scale;
    }
    Thet_tilde_e2[i] = d5;
    tau0_cross_tau[i] = d4;
    absxk += u[i] * d4;
  }
  hphi2[0] = 0.0;
  hphi2[3] = -Thet_tilde_e2[2];
  hphi2[6] = Thet_tilde_e2[1];
  hphi2[1] = Thet_tilde_e2[2];
  hphi2[4] = 0.0;
  hphi2[7] = -Thet_tilde_e2[0];
  hphi2[2] = -Thet_tilde_e2[1];
  hphi2[5] = Thet_tilde_e2[0];
  hphi2[8] = 0.0;
  dphi1xphi_idx_0 = ws2[1] * rotr[2] - rotr[1] * ws2[2];
  dphi1xphi_idx_1 = rotr[0] * ws2[2] - ws2[0] * rotr[2];
  dphi1xphi_idx_2 = ws2[0] * rotr[1] - rotr[0] * ws2[1];
  dv[0] = 0.0;
  dv[3] = -dphi1xphi_idx_2;
  dv[6] = dphi1xphi_idx_1;
  dv[1] = dphi1xphi_idx_2;
  dv[4] = 0.0;
  dv[7] = -dphi1xphi_idx_0;
  dv[2] = -dphi1xphi_idx_1;
  dv[5] = dphi1xphi_idx_0;
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d1 = hphi2[i];
    d2 = hphi2[i + 3];
    d3 = hphi2[i + 6];
    d4 = ws1hat[i];
    d5 = ws1hat[i + 3];
    d6 = ws1hat[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      tau0_dot_tau[b_i] =
          (rotrhat[i] * D2ws1[3 * i1] + rotrhat[i + 3] * D2ws1[Q_tmp]) +
          rotrhat[i + 6] * D2ws1[J_tmp];
      e_a[b_i] =
          (d1 * rotrhat[3 * i1] + d2 * rotrhat[Q_tmp]) + d3 * rotrhat[J_tmp];
      Drotr[b_i] =
          (d4 * ws2hat[3 * i1] + d5 * ws2hat[Q_tmp]) + d6 * ws2hat[J_tmp];
    }
    d1 = e_a[i];
    d2 = e_a[i + 3];
    d3 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_rotr[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                           d3 * RJdexpT_[3 * i1 + 2];
    }
    d1 = dv[i];
    d2 = dv[i + 3];
    d3 = dv[i + 6];
    d4 = Drotr[i];
    d5 = Drotr[i + 3];
    d6 = Drotr[i + 6];
    d7 = D2rotr_tmp[i];
    d8 = D2rotr_tmp[i + 3];
    d9 = D2rotr_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      e_a[b_i] = (tau0_dot_tau[b_i] + b_rotr[b_i]) +
                 ((d1 * Dws1[3 * i1] + d2 * Dws1[Q_tmp]) + d3 * Dws1[J_tmp]);
      b_rotr[b_i] =
          (d4 * rotrhat[3 * i1] + d5 * rotrhat[Q_tmp]) + d6 * rotrhat[J_tmp];
      tau0_dot_tau[b_i] =
          (d7 * Dws2[3 * i1] + d8 * Dws2[Q_tmp]) + d9 * Dws2[J_tmp];
    }
    d1 = b_rotr[i];
    d2 = b_rotr[i + 3];
    d3 = b_rotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      hphi[b_i] = -((e_a[b_i] + tau0_dot_tau[b_i]) +
                    ((d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                     d3 * RJdexpT_[3 * i1 + 2]));
    }
  }
  /*  D3tau is short form for D^3 tau(phi)(dphi1,dphi2,dot) */
  dphi1xphi_idx_0 = ws2[1] * Rb[2] - Rb[1] * ws2[2];
  dphi1xphi_idx_1 = Rb[0] * ws2[2] - ws2[0] * Rb[2];
  dphi1xphi_idx_2 = ws2[0] * Rb[1] - Rb[0] * ws2[1];
  /*  D3varthet1 is short form for D^3 varthet1(phi)(dphi1,dphi2,dot) */
  /*  D2tau2 is short form for D^2tau(phi)(dphi2,dot) */
  memset(&C[0], 0, 49U * sizeof(double));
  dv[0] = 0.0;
  dv[3] = -dphi1xphi_idx_2;
  dv[6] = dphi1xphi_idx_1;
  dv[1] = dphi1xphi_idx_2;
  dv[4] = 0.0;
  dv[7] = -dphi1xphi_idx_0;
  dv[2] = -dphi1xphi_idx_1;
  dv[5] = dphi1xphi_idx_0;
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d1 = hphi2[i];
    d2 = hphi2[i + 3];
    d3 = hphi2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      C[i1 + 7 * (i + 3)] = hphi[i1 + 3 * i];
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_a[b_i] = (tauhat[i] * D2ws1[3 * i1] + tauhat[i + 3] * D2ws1[Q_tmp]) +
                 tauhat[i + 6] * D2ws1[J_tmp];
      e_a[b_i] =
          (d1 * tauhat[3 * i1] + d2 * tauhat[Q_tmp]) + d3 * tauhat[J_tmp];
    }
    d1 = e_a[i];
    d2 = e_a[i + 3];
    d3 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      hphi2[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                          d3 * RJdexpT_[3 * i1 + 2];
    }
    d1 = Drotr[i];
    d2 = Drotr[i + 3];
    d3 = Drotr[i + 6];
    d4 = D2tau1_tmp[i];
    d5 = D2tau1_tmp[i + 3];
    d6 = D2tau1_tmp[i + 6];
    d7 = dv[i];
    d8 = dv[i + 3];
    d9 = dv[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      b_rotr[b_i] =
          (d1 * tauhat[3 * i1] + d2 * tauhat[Q_tmp]) + d3 * tauhat[J_tmp];
      tau0_dot_tau[b_i] =
          (d4 * Dws2[3 * i1] + d5 * Dws2[Q_tmp]) + d6 * Dws2[J_tmp];
      e_a[b_i] = (c_a[b_i] + hphi2[b_i]) +
                 ((d7 * Dws1[3 * i1] + d8 * Dws1[Q_tmp]) + d9 * Dws1[J_tmp]);
    }
    d1 = b_rotr[i];
    d2 = b_rotr[i + 3];
    d3 = b_rotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      c_a[b_i] = -((e_a[b_i] + tau0_dot_tau[b_i]) +
                   ((d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                    d3 * RJdexpT_[3 * i1 + 2]));
    }
    d1 = ws2hat[i];
    d2 = ws2hat[i + 3];
    d3 = ws2hat[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      e_a[b_i] = (tauhat[i] * Dws2[3 * i1] + tauhat[i + 3] * Dws2[Q_tmp]) +
                 tauhat[i + 6] * Dws2[J_tmp];
      b_rotr[b_i] =
          (d1 * tauhat[3 * i1] + d2 * tauhat[Q_tmp]) + d3 * tauhat[J_tmp];
    }
    d1 = b_rotr[i];
    d2 = b_rotr[i + 3];
    d3 = b_rotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      tau0_dot_tau[i + 3 * i1] =
          (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
          d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  d1 = qbar[6];
  d2 = rho1[6];
  d3 = rho2[6];
  d4 = rho2[3];
  d5 = rho2[4];
  d6 = rho2[5];
  d7 = u[0];
  d8 = u[1];
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    C[b_i + 3] = ((hphi[3 * i] + d1 * c_a[3 * i]) +
                  d2 * -(e_a[3 * i] + tau0_dot_tau[3 * i])) +
                 d3 * D2tau1[3 * i];
    d9 = tau0_cross_tau[i];
    b_a[3 * i] = absxk * (double)iv[3 * i] + d7 * d9;
    J_tmp = 3 * i + 1;
    C[b_i + 4] = ((hphi[J_tmp] + d1 * c_a[J_tmp]) +
                  d2 * -(e_a[J_tmp] + tau0_dot_tau[J_tmp])) +
                 d3 * D2tau1[J_tmp];
    b_a[J_tmp] = absxk * (double)iv[J_tmp] + d8 * d9;
    J_tmp = 3 * i + 2;
    C[b_i + 5] = ((hphi[J_tmp] + d1 * c_a[J_tmp]) +
                  d2 * -(e_a[J_tmp] + tau0_dot_tau[J_tmp])) +
                 d3 * D2tau1[J_tmp];
    b_a[J_tmp] = absxk * (double)iv[J_tmp] + d * d9;
    rotr[i] =
        (Dv_tmp[i] * tau0_cross_tau[0] + Dv_tmp[i + 3] * tau0_cross_tau[1]) +
        Dv_tmp[i + 6] * tau0_cross_tau[2];
    ws2[i] = (Dv[i] * d4 + Dv[i + 3] * d5) + Dv[i + 6] * d6;
  }
  for (i = 0; i < 3; i++) {
    d = b_a[i];
    d1 = b_a[i + 3];
    d2 = b_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      e_a[i1 + 3 * i] = rotr[i1] * u[i];
      c_a[i + 3 * i1] =
          (d * PP[3 * i1] + d1 * PP[3 * i1 + 1]) + d2 * PP[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    e_a[i] = d_a * (e_a[i] + c_a[i]);
  }
  for (i = 0; i < 3; i++) {
    d = ws2hat[i];
    d1 = ws2hat[i + 3];
    d2 = ws2hat[i + 6];
    d3 = a_tmp[i];
    d4 = a_tmp[i + 3];
    d5 = a_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      b_rotr[b_i] = (d * Dtau[3 * i1] + d1 * Dtau[Q_tmp]) + d2 * Dtau[J_tmp];
      tau0_dot_tau[b_i] =
          (d3 * Dws2[3 * i1] + d4 * Dws2[Q_tmp]) + d5 * Dws2[J_tmp];
    }
  }
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] += b_rotr[i];
  }
  for (i = 0; i < 3; i++) {
    d = Dv_tmp[i];
    d1 = Dv_tmp[i + 3];
    d2 = Dv_tmp[i + 6];
    d3 = e_a[i];
    d4 = e_a[i + 3];
    d5 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      J_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_a[b_i] = (c_a_tmp * d * tau0_dot_tau[3 * i1] +
                  c_a_tmp * d1 * tau0_dot_tau[Q_tmp]) +
                 c_a_tmp * d2 * tau0_dot_tau[J_tmp];
      b_a[b_i] = (d3 * Dtau[3 * i1] + d4 * Dtau[Q_tmp]) + d5 * Dtau[J_tmp];
    }
  }
  for (i = 0; i < 9; i++) {
    b_a[i] += c_a[i];
  }
  d = ws2[0];
  d1 = ws2[1];
  d2 = ws2[2];
  d3 = ws1[0];
  d4 = ws1[1];
  d5 = ws1[2];
  d6 = Thet_tilde_e2[0];
  d7 = Thet_tilde_e2[1];
  d8 = Thet_tilde_e2[2];
  d9 = v[0];
  d10 = v[1];
  d11 = v[2];
  scale = rho2[3];
  absxk = rho2[4];
  t = rho2[5];
  for (i = 0; i < 3; i++) {
    i1 = 3 * i + 1;
    Q_tmp = 3 * i + 2;
    C[7 * (i + 3) + 6] =
        ((((d * Dws1[3 * i] + d1 * Dws1[i1]) + d2 * Dws1[Q_tmp]) +
          ((d3 * b_a[3 * i] + d4 * b_a[i1]) + d5 * b_a[Q_tmp])) +
         ((d6 * Dv[3 * i] + d7 * Dv[i1]) + d8 * Dv[Q_tmp])) +
        ((d9 * D2ws1[3 * i] + d10 * D2ws1[i1]) + d11 * D2ws1[Q_tmp]);
    C[i + 45] = (D2tau1[i] * scale + D2tau1[i + 3] * absxk) + D2tau1[i + 6] * t;
  }
}

void d_CableBCTransinCoord(const double qbar[7], const double x0[3],
                           const double RJ[9], const double RE[9],
                           const double r[3], const double R0[9],
                           const emxArray_real_T *eta, const double rho1[7],
                           const double rho2[7], double q[7], double J[49],
                           double Q[49], double Qtilde[49], double C[49])
{
  double a[11];
  double D2rotr_tmp[9];
  double D2tau1[9];
  double D2tau1_tmp[9];
  double D2ws1[9];
  double Drotr[9];
  double Dtau[9];
  double Dv[9];
  double Dv_tmp[9];
  double Dws1[9];
  double Dws2[9];
  double PP[9];
  double RJdexpT_[9];
  double Rb[9];
  double a_tmp[9];
  double b_a[9];
  double c_a[9];
  double c_rotr[9];
  double dv[9];
  double e_a[9];
  double hphi[9];
  double hphi2[9];
  double rotrhat[9];
  double tau0_dot_tau[9];
  double tauhat[9];
  double ws1hat[9];
  double ws2hat[9];
  double Thet_tilde_e2[3];
  double b_rotr[3];
  double m[3];
  double rotr[3];
  double tau0_cross_tau[3];
  double u[3];
  double v[3];
  double ws1[3];
  double ws2[3];
  const double *eta_data;
  double absxk;
  double b_a_tmp;
  double b_qbar;
  double c_a_tmp;
  double d;
  double d1;
  double d10;
  double d11;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d9;
  double d_a;
  double d_a_tmp;
  double dphi1xphi_idx_0;
  double dphi1xphi_idx_1;
  double dphi1xphi_idx_2;
  double f_a;
  double g_a;
  double h_a;
  double i_a;
  double nu_;
  double scale;
  double t;
  int Q_tmp;
  int b_Q_tmp;
  int b_i;
  int i;
  int i1;
  eta_data = eta->data;
  /*  Computes the transformation from (d, phi, gamma) -> (p1, p2, vartheta) for
   */
  /*  one end of the cable; NOTE: the main difference from CableBCtrans.m is */
  /*  that this function is based on exponential coordinates for the rotation,
   */
  /*  rather than the rotation matrix; also, a joint coordinate system and end
   */
  /*  offset are included, making it a more general. */
  /*  Inputs */
  /*  qbar = [d; phi; gamm] */
  /*    where d = displacement of joint,  */
  /*          phi = exponential coordinates of rotation of joint */
  /*          gamm = distance between first and second control points */
  /*  x0 = reference position of joint (3x1) vector */
  /*  RJ = rotation of joint coordinate frame with respect to global */
  /*  RE = rotation of cable end with respect to joint coordinate frame */
  /*  r = position of cable end relative to joint in joint coordinate system */
  /*      (end offset, 3x1 vector) */
  /*  R0 = orientation of cable end in reference configuration */
  /*  eta = vector to contract with for second derivative, Q */
  /*  rho1 = vector to contract with for second derivative, Qtilde */
  /*  rho2 = vector to contract with for third derivative, C */
  /*  Outputs: */
  /*  q = [p1; p2; vartheta1] */
  /*    where p1, p2 = positions of first and second control points */
  /*          vartheta1 = (twist) rotation of first control point */
  /*  J = first derivative of transformation */
  /*  Q = second derivative contracted over upper index,  */
  /*      Q_{ib,jb} = eta_i q^i_{ib,jb} */
  /*  Qtilde = second derivative contracted over one of the lower indices */
  /*      Qtilde^{i}_{ib} = rho1^{jb} q^i_{ib,jb} */
  /*  C = third derivative contracted over two lower indices */
  /*      C^{i}_{ib} = q^i_{ib,jb,kb}rho1^{jb}rho2^{kb} */
  b_rothelper(&qbar[3], a);
  hphi[0] = 0.0;
  hphi[3] = -qbar[5];
  hphi[6] = qbar[4];
  hphi[1] = qbar[5];
  hphi[4] = 0.0;
  hphi[7] = -qbar[3];
  hphi[2] = -qbar[4];
  hphi[5] = qbar[3];
  hphi[8] = 0.0;
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      hphi2[i + 3 * i1] =
          (hphi[i] * hphi[3 * i1] + hphi[i + 3] * hphi[3 * i1 + 1]) +
          hphi[i + 6] * hphi[3 * i1 + 2];
    }
  }
  d = a[0];
  d1 = a[1];
  for (i = 0; i < 9; i++) {
    b_a[i] = (d * hphi[i] + (double)iv[i]) + d1 * hphi2[i];
  }
  scale = 0.0;
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Drotr[i + 3 * i1] =
          (d * b_a[3 * i1] + d1 * b_a[3 * i1 + 1]) + d2 * b_a[3 * i1 + 2];
    }
    d = Drotr[i];
    d1 = Drotr[i + 3];
    d2 = Drotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Rb[i + 3 * i1] =
          (d * RE[3 * i1] + d1 * RE[3 * i1 + 1]) + d2 * RE[3 * i1 + 2];
    }
    scale += R0[i] * Rb[i];
  }
  tau0_cross_tau[0] = R0[1] * Rb[2] - Rb[1] * R0[2];
  tau0_cross_tau[1] = Rb[0] * R0[2] - R0[0] * Rb[2];
  tau0_cross_tau[2] = R0[0] * Rb[1] - Rb[0] * R0[1];
  dv[0] = 0.0;
  dv[1] = -tau0_cross_tau[2];
  dv[2] = tau0_cross_tau[1];
  dv[3] = tau0_cross_tau[2];
  dv[4] = 0.0;
  dv[5] = -tau0_cross_tau[0];
  dv[6] = -tau0_cross_tau[1];
  dv[7] = tau0_cross_tau[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d = tau0_cross_tau[i] / (scale + 1.0);
    tau0_dot_tau[3 * i] =
        (scale * (double)iv[i] + dv[3 * i]) + d * tau0_cross_tau[0];
    b_i = 3 * i + 1;
    tau0_dot_tau[b_i] =
        (scale * (double)iv[i + 3] + dv[b_i]) + d * tau0_cross_tau[1];
    b_i = 3 * i + 2;
    tau0_dot_tau[b_i] =
        (scale * (double)iv[i + 6] + dv[b_i]) + d * tau0_cross_tau[2];
  }
  d = Rb[3];
  d1 = Rb[4];
  d2 = Rb[5];
  for (i = 0; i < 3; i++) {
    m[i] = (tau0_dot_tau[i] * d + tau0_dot_tau[i + 3] * d1) +
           tau0_dot_tau[i + 6] * d2;
  }
  d = r[0];
  d1 = r[1];
  d2 = r[2];
  d3 = m[0];
  d4 = m[1];
  d5 = m[2];
  d6 = qbar[0];
  d7 = qbar[1];
  d8 = qbar[2];
  d9 = qbar[6];
  for (i = 0; i < 3; i++) {
    d10 = (Drotr[i] * d + Drotr[i + 3] * d1) + Drotr[i + 6] * d2;
    rotr[i] = d10;
    Thet_tilde_e2[i] =
        (R0[3 * i] * d3 + R0[3 * i + 1] * d4) + R0[3 * i + 2] * d5;
    d10 += x0[i] + ((RJ[i] * d6 + RJ[i + 3] * d7) + RJ[i + 6] * d8);
    d11 = d9 * Rb[i];
    ws2[i] = d11;
    q[i] = d10;
    q[i + 3] = d10 + d11;
  }
  q[6] = atan2(Thet_tilde_e2[2], Thet_tilde_e2[1]);
  /*  first derivative needed */
  scale = 3.3121686421112381E-170;
  d = Rb[0] + R0[0];
  u[0] = d;
  absxk = fabs(d);
  if (absxk > 3.3121686421112381E-170) {
    nu_ = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    nu_ = t * t;
  }
  d = Rb[1] + R0[1];
  u[1] = d;
  absxk = fabs(d);
  if (absxk > scale) {
    t = scale / absxk;
    nu_ = nu_ * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nu_ += t * t;
  }
  d = Rb[2] + R0[2];
  u[2] = d;
  absxk = fabs(d);
  if (absxk > scale) {
    t = scale / absxk;
    nu_ = nu_ * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nu_ += t * t;
  }
  nu_ = scale * sqrt(nu_);
  absxk = nu_ * nu_;
  d1 = a[1];
  d2 = a[8];
  for (b_i = 0; b_i < 3; b_i++) {
    v[b_i] = 2.0 * u[b_i] / absxk;
    tau0_dot_tau[3 * b_i] =
        ((double)iv[b_i] - d1 * hphi[b_i]) + d2 * hphi2[b_i];
    tau0_dot_tau[3 * b_i + 1] =
        ((double)iv[b_i + 3] - d1 * hphi[b_i + 3]) + d2 * hphi2[b_i + 3];
    tau0_dot_tau[3 * b_i + 2] =
        ((double)iv[b_i + 6] - d1 * hphi[b_i + 6]) + d2 * hphi2[b_i + 6];
  }
  for (i = 0; i < 3; i++) {
    d1 = RJ[i];
    d2 = RJ[i + 3];
    d3 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      RJdexpT_[i + 3 * i1] =
          (d1 * tau0_dot_tau[3 * i1] + d2 * tau0_dot_tau[3 * i1 + 1]) +
          d3 * tau0_dot_tau[3 * i1 + 2];
    }
  }
  rotrhat[0] = 0.0;
  rotrhat[3] = -rotr[2];
  rotrhat[6] = rotr[1];
  rotrhat[1] = rotr[2];
  rotrhat[4] = 0.0;
  rotrhat[7] = -rotr[0];
  rotrhat[2] = -rotr[1];
  rotrhat[5] = rotr[0];
  rotrhat[8] = 0.0;
  tauhat[0] = 0.0;
  tauhat[3] = -Rb[2];
  tauhat[6] = Rb[1];
  tauhat[1] = Rb[2];
  tauhat[4] = 0.0;
  tauhat[7] = -Rb[0];
  tauhat[2] = -Rb[1];
  tauhat[5] = Rb[0];
  tauhat[8] = 0.0;
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] = -rotrhat[i];
  }
  for (i = 0; i < 3; i++) {
    d1 = tau0_dot_tau[i];
    d2 = tau0_dot_tau[i + 3];
    d3 = tau0_dot_tau[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Drotr[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                          d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    a_tmp[i] = -tauhat[i];
  }
  for (i = 0; i < 3; i++) {
    d1 = a_tmp[i];
    d2 = a_tmp[i + 3];
    d3 = a_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Dtau[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                         d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  memset(&J[0], 0, 49U * sizeof(double));
  /*  second derivative needed */
  /*  ws1 is short form for w(phi;dphi1) */
  ws1hat[0] = 0.0;
  ws1hat[4] = 0.0;
  ws1hat[8] = 0.0;
  /*  Dws1 is short form for Dw(phi;dphi1) */
  b_qbar = 0.0;
  dphi1xphi_idx_0 = rho1[4] * qbar[5] - qbar[4] * rho1[5];
  dphi1xphi_idx_1 = qbar[3] * rho1[5] - rho1[3] * qbar[5];
  dphi1xphi_idx_2 = rho1[3] * qbar[4] - qbar[3] * rho1[4];
  hphi2[0] = 0.0;
  hphi2[3] = -rho1[5];
  hphi2[6] = rho1[4];
  hphi2[1] = rho1[5];
  hphi2[4] = 0.0;
  hphi2[7] = -rho1[3];
  hphi2[2] = -rho1[4];
  hphi2[5] = rho1[3];
  hphi2[8] = 0.0;
  d1 = v[0];
  d2 = v[1];
  d3 = v[2];
  d4 = a[2];
  d5 = a[3];
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    scale = qbar[i + 3];
    b_qbar += scale * rho1[i + 3];
    d6 = RJ[3 * i];
    J[7 * i] = d6;
    J[7 * i + 3] = d6;
    d6 = Drotr[3 * i];
    J[b_i] = d6;
    J[b_i + 3] = d6 + qbar[6] * Dtau[3 * i];
    Drotr[3 * i] = qbar[3] * scale;
    b_a[3 * i] = d4 * rho1[3] * scale;
    c_a[3 * i] = d5 * dphi1xphi_idx_0 * scale;
    i1 = 3 * i + 1;
    d6 = RJ[i1];
    J[7 * i + 1] = d6;
    J[7 * i + 4] = d6;
    d6 = Drotr[i1];
    J[b_i + 1] = d6;
    J[b_i + 4] = d6 + qbar[6] * Dtau[i1];
    Drotr[i1] = qbar[4] * scale;
    b_a[i1] = d4 * rho1[4] * scale;
    c_a[i1] = d5 * dphi1xphi_idx_1 * scale;
    Q_tmp = 3 * i + 2;
    d6 = RJ[Q_tmp];
    J[7 * i + 2] = d6;
    J[7 * i + 5] = d6;
    d6 = Drotr[Q_tmp];
    J[b_i + 2] = d6;
    J[b_i + 5] = d6 + qbar[6] * Dtau[Q_tmp];
    J[b_i + 6] =
        (d1 * RJdexpT_[3 * i] + d2 * RJdexpT_[i1]) + d3 * RJdexpT_[Q_tmp];
    Drotr[Q_tmp] = qbar[5] * scale;
    b_a[Q_tmp] = d4 * rho1[5] * scale;
    c_a[Q_tmp] = d5 * dphi1xphi_idx_2 * scale;
    ws1[i] = (RJdexpT_[i] * rho1[3] + RJdexpT_[i + 3] * rho1[4]) +
             RJdexpT_[i + 6] * rho1[5];
    J[i + 45] = Rb[i];
  }
  ws1hat[3] = -ws1[2];
  ws1hat[6] = ws1[1];
  ws1hat[1] = ws1[2];
  ws1hat[7] = -ws1[0];
  ws1hat[2] = -ws1[1];
  ws1hat[5] = ws1[0];
  b_a_tmp = a[9] * b_qbar;
  d_a = a[8] * b_qbar;
  d1 = a[1];
  d2 = a[8];
  d3 = qbar[3];
  d4 = qbar[4];
  d5 = qbar[5];
  for (i = 0; i < 3; i++) {
    c_a_tmp = rho1[i + 3];
    e_a[3 * i] = ((((b_a[3 * i] - c_a[3 * i]) - d1 * hphi2[3 * i]) +
                   b_a_tmp * Drotr[3 * i]) +
                  d3 * d2 * c_a_tmp) +
                 d_a * (double)iv[3 * i];
    b_i = 3 * i + 1;
    e_a[b_i] =
        ((((b_a[b_i] - c_a[b_i]) - d1 * hphi2[b_i]) + b_a_tmp * Drotr[b_i]) +
         d4 * d2 * c_a_tmp) +
        d_a * (double)iv[b_i];
    b_i = 3 * i + 2;
    e_a[b_i] =
        ((((b_a[b_i] - c_a[b_i]) - d1 * hphi2[b_i]) + b_a_tmp * Drotr[b_i]) +
         d5 * d2 * c_a_tmp) +
        d_a * (double)iv[b_i];
  }
  /*  D2rotr is short form for D^2rotr(phi)(dphi1,dot) */
  for (i = 0; i < 3; i++) {
    d1 = ws1hat[i];
    d2 = ws1hat[i + 3];
    d3 = ws1hat[i + 6];
    d4 = RJ[i];
    d5 = RJ[i + 3];
    d6 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      D2rotr_tmp[b_i] =
          (d1 * rotrhat[3 * i1] + d2 * rotrhat[Q_tmp]) + d3 * rotrhat[b_Q_tmp];
      Dws1[b_i] = (d4 * e_a[3 * i1] + d5 * e_a[Q_tmp]) + d6 * e_a[b_Q_tmp];
    }
  }
  for (i = 0; i < 3; i++) {
    d1 = D2rotr_tmp[i];
    d2 = D2rotr_tmp[i + 3];
    d3 = D2rotr_tmp[i + 6];
    d4 = rotrhat[i];
    d5 = rotrhat[i + 3];
    d6 = rotrhat[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      tau0_dot_tau[b_i] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[Q_tmp]) +
                          d3 * RJdexpT_[b_Q_tmp];
      hphi[b_i] = (d4 * Dws1[3 * i1] + d5 * Dws1[Q_tmp]) + d6 * Dws1[b_Q_tmp];
    }
  }
  for (i = 0; i < 9; i++) {
    hphi[i] = -(hphi[i] + tau0_dot_tau[i]);
  }
  /*  D2tau1 is short form for D^2tau(phi)(dphi1,dot) */
  for (i = 0; i < 3; i++) {
    d1 = ws1hat[i];
    d2 = ws1hat[i + 3];
    d3 = ws1hat[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      D2tau1[b_i] = (tauhat[i] * Dws1[3 * i1] + tauhat[i + 3] * Dws1[Q_tmp]) +
                    tauhat[i + 6] * Dws1[b_Q_tmp];
      D2tau1_tmp[b_i] =
          (d1 * tauhat[3 * i1] + d2 * tauhat[Q_tmp]) + d3 * tauhat[b_Q_tmp];
    }
    d1 = D2tau1_tmp[i];
    d2 = D2tau1_tmp[i + 3];
    d3 = D2tau1_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      tau0_dot_tau[i + 3 * i1] =
          (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
          d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    D2tau1[i] = -(D2tau1[i] + tau0_dot_tau[i]);
  }
  tau0_cross_tau[0] = u[0] / nu_;
  tau0_cross_tau[1] = u[1] / nu_;
  tau0_cross_tau[2] = d / nu_;
  /*  unit vector in the direction of u */
  for (i = 0; i < 3; i++) {
    PP[3 * i] = (double)iv[3 * i] - tau0_cross_tau[0] * tau0_cross_tau[i];
    b_i = 3 * i + 1;
    PP[b_i] = (double)iv[b_i] - tau0_cross_tau[1] * tau0_cross_tau[i];
    b_i = 3 * i + 2;
    PP[b_i] = (double)iv[b_i] - tau0_cross_tau[2] * tau0_cross_tau[i];
  }
  for (i = 0; i < 9; i++) {
    Dv_tmp[i] = 2.0 * PP[i] - (double)iv[i];
  }
  c_a_tmp = 2.0 / absxk;
  for (i = 0; i < 3; i++) {
    d1 = Dv_tmp[i];
    d2 = Dv_tmp[i + 3];
    d3 = Dv_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_a[i + 3 * i1] = (-c_a_tmp * d1 * tauhat[3 * i1] +
                         -c_a_tmp * d2 * tauhat[3 * i1 + 1]) +
                        -c_a_tmp * d3 * tauhat[3 * i1 + 2];
    }
    d1 = b_a[i];
    d2 = b_a[i + 3];
    d3 = b_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Dv[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                       d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  memset(&Qtilde[0], 0, 49U * sizeof(double));
  d1 = qbar[6];
  d2 = rho1[6];
  d3 = rho1[3];
  d4 = rho1[4];
  d5 = rho1[5];
  d6 = ws1[0];
  d7 = ws1[1];
  d8 = ws1[2];
  d9 = v[0];
  d10 = v[1];
  d11 = v[2];
  for (i = 0; i < 3; i++) {
    scale = hphi[3 * i];
    b_i = 7 * (i + 3);
    Qtilde[b_i] = scale;
    Qtilde[b_i + 3] = (scale + d1 * D2tau1[3 * i]) + d2 * Dtau[3 * i];
    i1 = 3 * i + 1;
    scale = hphi[i1];
    Qtilde[b_i + 1] = scale;
    Qtilde[b_i + 4] = (scale + d1 * D2tau1[i1]) + d2 * Dtau[i1];
    Q_tmp = 3 * i + 2;
    scale = hphi[Q_tmp];
    Qtilde[b_i + 2] = scale;
    Qtilde[b_i + 5] = (scale + d1 * D2tau1[Q_tmp]) + d2 * Dtau[Q_tmp];
    Qtilde[b_i + 6] = ((d6 * Dv[3 * i] + d7 * Dv[i1]) + d8 * Dv[Q_tmp]) +
                      ((d9 * Dws1[3 * i] + d10 * Dws1[i1]) + d11 * Dws1[Q_tmp]);
    Qtilde[i + 45] = (Dtau[i] * d3 + Dtau[i + 3] * d4) + Dtau[i + 6] * d5;
    absxk = eta_data[i + 3];
    tau0_cross_tau[i] = absxk;
    Thet_tilde_e2[i] = eta_data[i] + absxk;
  }
  b_rotr[0] = ((rotr[1] * Thet_tilde_e2[2] - Thet_tilde_e2[1] * rotr[2]) +
               qbar[6] * (Rb[1] * eta_data[5] - Rb[2] * eta_data[4])) +
              v[0] * eta_data[6];
  b_rotr[1] = ((Thet_tilde_e2[0] * rotr[2] - rotr[0] * Thet_tilde_e2[2]) +
               qbar[6] * (Rb[2] * eta_data[3] - Rb[0] * eta_data[5])) +
              v[1] * eta_data[6];
  b_rotr[2] = ((rotr[0] * Thet_tilde_e2[1] - Thet_tilde_e2[0] * rotr[1]) +
               qbar[6] * (Rb[0] * eta_data[4] - Rb[1] * eta_data[3])) +
              v[2] * eta_data[6];
  scale = 0.0;
  d1 = b_rotr[0];
  d2 = b_rotr[1];
  d3 = b_rotr[2];
  for (i = 0; i < 3; i++) {
    d4 = (RJ[3 * i] * d1 + RJ[3 * i + 1] * d2) + RJ[3 * i + 2] * d3;
    m[i] = d4;
    scale += d4 * qbar[i + 3];
  }
  d_a = a[9] * scale;
  f_a = a[8] * scale;
  memset(&Q[0], 0, 49U * sizeof(double));
  g_a = ((rotr[0] * Thet_tilde_e2[0] + rotr[1] * Thet_tilde_e2[1]) +
         rotr[2] * Thet_tilde_e2[2]) +
        qbar[6] * ((Rb[0] * tau0_cross_tau[0] + Rb[1] * tau0_cross_tau[1]) +
                   Rb[2] * tau0_cross_tau[2]);
  dv[0] = 0.0;
  dv[3] = a[1] * -m[2];
  dv[6] = a[1] * m[1];
  dv[1] = a[1] * m[2];
  dv[4] = 0.0;
  dv[7] = a[1] * -m[0];
  dv[2] = a[1] * -m[1];
  dv[5] = m[0] * a[1];
  dv[8] = 0.0;
  d1 = m[0];
  d2 = m[1];
  d3 = m[2];
  d4 = a[2];
  d5 = a[3] * (m[2] * qbar[4] - m[1] * qbar[5]);
  d6 = a[3] * (m[0] * qbar[5] - m[2] * qbar[3]);
  d7 = a[3] * (m[1] * qbar[3] - m[0] * qbar[4]);
  d8 = rotr[0];
  d9 = rotr[1];
  d10 = rotr[2];
  d11 = ws2[0];
  scale = ws2[1];
  absxk = ws2[2];
  for (i = 0; i < 3; i++) {
    t = qbar[i + 3];
    b_a[3 * i] = d1 * d4 * t;
    c_a[3 * i] = d5 * t;
    h_a = Thet_tilde_e2[i];
    e_a[3 * i] = d8 * h_a;
    i_a = tau0_cross_tau[i];
    tau0_dot_tau[3 * i] = d11 * i_a;
    b_i = 3 * i + 1;
    b_a[b_i] = d2 * d4 * t;
    c_a[b_i] = d6 * t;
    e_a[b_i] = d9 * h_a;
    tau0_dot_tau[b_i] = scale * i_a;
    b_i = 3 * i + 2;
    b_a[b_i] = d4 * d3 * t;
    c_a[b_i] = d7 * t;
    e_a[b_i] = d10 * h_a;
    tau0_dot_tau[b_i] = absxk * i_a;
  }
  for (i = 0; i < 9; i++) {
    e_a[i] = (e_a[i] + tau0_dot_tau[i]) - g_a * (double)iv[i];
  }
  for (i = 0; i < 3; i++) {
    d1 = e_a[i];
    d2 = e_a[i + 3];
    d3 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      c_rotr[b_i] = ((d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                     d3 * RJdexpT_[3 * i1 + 2]) +
                    eta_data[6] * Dv[b_i];
    }
  }
  for (i = 0; i < 3; i++) {
    d1 = RJdexpT_[3 * i];
    d2 = RJdexpT_[3 * i + 1];
    d3 = RJdexpT_[3 * i + 2];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i1 + 3 * i;
      e_a[b_i] = ((((b_a[b_i] - c_a[b_i]) + dv[b_i]) + d_a * Drotr[b_i]) +
                  f_a * (double)iv[b_i]) +
                 a[8] * qbar[i1 + 3] * m[i];
      tau0_dot_tau[i + 3 * i1] =
          (d1 * c_rotr[3 * i1] + d2 * c_rotr[3 * i1 + 1]) +
          d3 * c_rotr[3 * i1 + 2];
    }
  }
  d1 = tau0_cross_tau[0];
  d2 = tau0_cross_tau[1];
  d3 = tau0_cross_tau[2];
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    Q[b_i + 3] = e_a[3 * i] + tau0_dot_tau[3 * i];
    b_Q_tmp = 3 * i + 1;
    Q[b_i + 4] = e_a[b_Q_tmp] + tau0_dot_tau[b_Q_tmp];
    Q_tmp = 3 * i + 2;
    Q[b_i + 5] = e_a[Q_tmp] + tau0_dot_tau[Q_tmp];
    d4 = (d1 * Dtau[3 * i] + d2 * Dtau[b_Q_tmp]) + d3 * Dtau[Q_tmp];
    Q[b_i + 6] = d4;
    Q[i + 45] = d4;
  }
  /*  third derivtive needed */
  /*  ws2 is short form for ws(phi;dphi2) */
  ws2hat[0] = 0.0;
  ws2hat[4] = 0.0;
  ws2hat[8] = 0.0;
  /*  Dws2 is short form for Dws(phi;dphi2) */
  scale = 0.0;
  d1 = a[2];
  d2 = a[3] * (rho2[4] * qbar[5] - qbar[4] * rho2[5]);
  d3 = a[3] * (qbar[3] * rho2[5] - rho2[3] * qbar[5]);
  d4 = a[3] * (rho2[3] * qbar[4] - qbar[3] * rho2[4]);
  for (i = 0; i < 3; i++) {
    d5 = qbar[i + 3];
    scale += d5 * rho2[i + 3];
    b_a[3 * i] = d1 * rho2[3] * d5;
    c_a[3 * i] = d2 * d5;
    b_i = 3 * i + 1;
    b_a[b_i] = d1 * rho2[4] * d5;
    c_a[b_i] = d3 * d5;
    b_i = 3 * i + 2;
    b_a[b_i] = d1 * rho2[5] * d5;
    c_a[b_i] = d4 * d5;
    ws2[i] = (RJdexpT_[i] * rho2[3] + RJdexpT_[i + 3] * rho2[4]) +
             RJdexpT_[i + 6] * rho2[5];
  }
  ws2hat[3] = -ws2[2];
  ws2hat[6] = ws2[1];
  ws2hat[1] = ws2[2];
  ws2hat[7] = -ws2[0];
  ws2hat[2] = -ws2[1];
  ws2hat[5] = ws2[0];
  d_a_tmp = a[9] * scale;
  d_a = a[8] * scale;
  dv[0] = 0.0;
  dv[3] = a[1] * -rho2[5];
  dv[6] = a[1] * rho2[4];
  dv[1] = a[1] * rho2[5];
  dv[4] = 0.0;
  dv[7] = a[1] * -rho2[3];
  dv[2] = a[1] * -rho2[4];
  dv[5] = a[1] * rho2[3];
  dv[8] = 0.0;
  d1 = qbar[3];
  d2 = qbar[4];
  d3 = qbar[5];
  d4 = a[8];
  for (i = 0; i < 3; i++) {
    t = rho2[i + 3];
    e_a[3 * i] =
        ((((b_a[3 * i] - c_a[3 * i]) - dv[3 * i]) + d_a_tmp * Drotr[3 * i]) +
         d1 * d4 * t) +
        d_a * (double)iv[3 * i];
    b_i = 3 * i + 1;
    e_a[b_i] = ((((b_a[b_i] - c_a[b_i]) - dv[b_i]) + d_a_tmp * Drotr[b_i]) +
                d2 * d4 * t) +
               d_a * (double)iv[b_i];
    b_i = 3 * i + 2;
    e_a[b_i] = ((((b_a[b_i] - c_a[b_i]) - dv[b_i]) + d_a_tmp * Drotr[b_i]) +
                d3 * d4 * t) +
               d_a * (double)iv[b_i];
  }
  /*  D2ws1 is short form for D^2 ws(phi;dphi1)(dphi2,dot) */
  absxk = 0.0;
  d_a = a[4] * scale;
  f_a = a[5] * scale;
  g_a = a[3] * scale;
  i_a = a[10] * b_qbar * scale;
  h_a = b_a_tmp * scale;
  for (i = 0; i < 3; i++) {
    scale = rho2[i + 3];
    absxk += rho1[i + 3] * scale;
    d1 = RJ[i];
    d2 = RJ[i + 3];
    d3 = RJ[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Dws2[i + 3 * i1] =
          (d1 * e_a[3 * i1] + d2 * e_a[3 * i1 + 1]) + d3 * e_a[3 * i1 + 2];
      t = rho1[i1 + 3];
      b_i = i1 + 3 * i;
      b_a[b_i] = d_a * t * qbar[i + 3];
      c_a[b_i] = a[2] * t * scale;
    }
  }
  d_a = a[9] * absxk;
  scale = a[8] * absxk;
  d1 = a[3];
  d2 = a[8];
  d3 = a[3] * (rho1[4] * rho2[5] - rho2[4] * rho1[5]);
  d4 = a[3] * (rho2[3] * rho1[5] - rho1[3] * rho2[5]);
  d5 = a[3] * (rho1[3] * rho2[4] - rho2[3] * rho1[4]);
  for (i = 0; i < 3; i++) {
    signed char i2;
    d6 = rho1[i + 3];
    d7 = Drotr[3 * i];
    d8 = qbar[i + 3];
    d9 = rho2[i + 3];
    d10 = (((((b_a[3 * i] + c_a[3 * i]) - f_a * dphi1xphi_idx_0 * d8) -
             d1 * dphi1xphi_idx_0 * d9) -
            g_a * hphi2[3 * i]) -
           d3 * d8) +
          i_a * d7;
    i2 = iv[3 * i];
    d11 = b_a_tmp * rho2[3] * d8;
    b_a[3 * i] = d11;
    c_a[3 * i] = ((((((d10 + d_a_tmp * qbar[3] * d6) + b_a_tmp * qbar[3] * d9) +
                     h_a * (double)i2) +
                    d_a * d7) +
                   scale * (double)i2) +
                  d11) +
                 rho2[3] * d2 * d6;
    i1 = 3 * i + 1;
    d7 = Drotr[i1];
    d10 = (((((b_a[i1] + c_a[i1]) - f_a * dphi1xphi_idx_1 * d8) -
             d1 * dphi1xphi_idx_1 * d9) -
            g_a * hphi2[i1]) -
           d4 * d8) +
          i_a * d7;
    i2 = iv[i1];
    d11 = b_a_tmp * rho2[4] * d8;
    b_a[i1] = d11;
    c_a[i1] = ((((((d10 + d_a_tmp * qbar[4] * d6) + b_a_tmp * qbar[4] * d9) +
                  h_a * (double)i2) +
                 d_a * d7) +
                scale * (double)i2) +
               d11) +
              rho2[4] * d2 * d6;
    i1 = 3 * i + 2;
    d7 = Drotr[i1];
    d10 = (((((b_a[i1] + c_a[i1]) - f_a * dphi1xphi_idx_2 * d8) -
             d1 * dphi1xphi_idx_2 * d9) -
            g_a * hphi2[i1]) -
           d5 * d8) +
          i_a * d7;
    i2 = iv[i1];
    d11 = b_a_tmp * rho2[5] * d8;
    b_a[i1] = d11;
    c_a[i1] = ((((((d10 + d_a_tmp * qbar[5] * d6) + b_a_tmp * qbar[5] * d9) +
                  h_a * (double)i2) +
                 d_a * d7) +
                scale * (double)i2) +
               d11) +
              rho2[5] * d2 * d6;
  }
  /*  D2v is short form for D^2 v(phi)(dphi2,dot) */
  scale = 0.0;
  d_a = -(4.0 / pow(nu_, 4.0));
  /*  D3rotr is short form for D^3 rotr(phi)(dphi1,dphi2,dot) */
  for (i = 0; i < 3; i++) {
    d1 = RJ[i];
    d2 = RJ[i + 3];
    d3 = RJ[i + 6];
    d4 = 0.0;
    d5 = 0.0;
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      D2ws1[b_i] =
          (d1 * c_a[3 * i1] + d2 * c_a[3 * i1 + 1]) + d3 * c_a[3 * i1 + 2];
      absxk = rho2[i1 + 3];
      d4 += Dtau[b_i] * absxk;
      d5 += Dws1[b_i] * absxk;
    }
    Thet_tilde_e2[i] = d5;
    tau0_cross_tau[i] = d4;
    scale += u[i] * d4;
  }
  hphi2[0] = 0.0;
  hphi2[3] = -Thet_tilde_e2[2];
  hphi2[6] = Thet_tilde_e2[1];
  hphi2[1] = Thet_tilde_e2[2];
  hphi2[4] = 0.0;
  hphi2[7] = -Thet_tilde_e2[0];
  hphi2[2] = -Thet_tilde_e2[1];
  hphi2[5] = Thet_tilde_e2[0];
  hphi2[8] = 0.0;
  m[0] = ws2[1] * rotr[2] - rotr[1] * ws2[2];
  m[1] = rotr[0] * ws2[2] - ws2[0] * rotr[2];
  m[2] = ws2[0] * rotr[1] - rotr[0] * ws2[1];
  dv[0] = 0.0;
  dv[3] = -m[2];
  dv[6] = m[1];
  dv[1] = m[2];
  dv[4] = 0.0;
  dv[7] = -m[0];
  dv[2] = -m[1];
  dv[5] = m[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d1 = hphi2[i];
    d2 = hphi2[i + 3];
    d3 = hphi2[i + 6];
    d4 = ws1hat[i];
    d5 = ws1hat[i + 3];
    d6 = ws1hat[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      tau0_dot_tau[b_i] =
          (rotrhat[i] * D2ws1[3 * i1] + rotrhat[i + 3] * D2ws1[Q_tmp]) +
          rotrhat[i + 6] * D2ws1[b_Q_tmp];
      e_a[b_i] =
          (d1 * rotrhat[3 * i1] + d2 * rotrhat[Q_tmp]) + d3 * rotrhat[b_Q_tmp];
      Drotr[b_i] =
          (d4 * ws2hat[3 * i1] + d5 * ws2hat[Q_tmp]) + d6 * ws2hat[b_Q_tmp];
    }
    d1 = e_a[i];
    d2 = e_a[i + 3];
    d3 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      c_rotr[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                           d3 * RJdexpT_[3 * i1 + 2];
    }
    d1 = dv[i];
    d2 = dv[i + 3];
    d3 = dv[i + 6];
    d4 = Drotr[i];
    d5 = Drotr[i + 3];
    d6 = Drotr[i + 6];
    d7 = D2rotr_tmp[i];
    d8 = D2rotr_tmp[i + 3];
    d9 = D2rotr_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      e_a[b_i] = (tau0_dot_tau[b_i] + c_rotr[b_i]) +
                 ((d1 * Dws1[3 * i1] + d2 * Dws1[Q_tmp]) + d3 * Dws1[b_Q_tmp]);
      c_rotr[b_i] =
          (d4 * rotrhat[3 * i1] + d5 * rotrhat[Q_tmp]) + d6 * rotrhat[b_Q_tmp];
      tau0_dot_tau[b_i] =
          (d7 * Dws2[3 * i1] + d8 * Dws2[Q_tmp]) + d9 * Dws2[b_Q_tmp];
    }
    d1 = c_rotr[i];
    d2 = c_rotr[i + 3];
    d3 = c_rotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      hphi[b_i] = -((e_a[b_i] + tau0_dot_tau[b_i]) +
                    ((d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                     d3 * RJdexpT_[3 * i1 + 2]));
    }
  }
  /*  D3tau is short form for D^3 tau(phi)(dphi1,dphi2,dot) */
  m[0] = ws2[1] * Rb[2] - Rb[1] * ws2[2];
  m[1] = Rb[0] * ws2[2] - ws2[0] * Rb[2];
  m[2] = ws2[0] * Rb[1] - Rb[0] * ws2[1];
  /*  D3varthet1 is short form for D^3 varthet1(phi)(dphi1,dphi2,dot) */
  /*  D2tau2 is short form for D^2tau(phi)(dphi2,dot) */
  memset(&C[0], 0, 49U * sizeof(double));
  dv[0] = 0.0;
  dv[3] = -m[2];
  dv[6] = m[1];
  dv[1] = m[2];
  dv[4] = 0.0;
  dv[7] = -m[0];
  dv[2] = -m[1];
  dv[5] = m[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d1 = hphi2[i];
    d2 = hphi2[i + 3];
    d3 = hphi2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      C[i1 + 7 * (i + 3)] = hphi[i1 + 3 * i];
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_a[b_i] = (tauhat[i] * D2ws1[3 * i1] + tauhat[i + 3] * D2ws1[Q_tmp]) +
                 tauhat[i + 6] * D2ws1[b_Q_tmp];
      e_a[b_i] =
          (d1 * tauhat[3 * i1] + d2 * tauhat[Q_tmp]) + d3 * tauhat[b_Q_tmp];
    }
    d1 = e_a[i];
    d2 = e_a[i + 3];
    d3 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      hphi2[i + 3 * i1] = (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                          d3 * RJdexpT_[3 * i1 + 2];
    }
    d1 = Drotr[i];
    d2 = Drotr[i + 3];
    d3 = Drotr[i + 6];
    d4 = D2tau1_tmp[i];
    d5 = D2tau1_tmp[i + 3];
    d6 = D2tau1_tmp[i + 6];
    d7 = dv[i];
    d8 = dv[i + 3];
    d9 = dv[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_rotr[b_i] =
          (d1 * tauhat[3 * i1] + d2 * tauhat[Q_tmp]) + d3 * tauhat[b_Q_tmp];
      tau0_dot_tau[b_i] =
          (d4 * Dws2[3 * i1] + d5 * Dws2[Q_tmp]) + d6 * Dws2[b_Q_tmp];
      e_a[b_i] = (c_a[b_i] + hphi2[b_i]) +
                 ((d7 * Dws1[3 * i1] + d8 * Dws1[Q_tmp]) + d9 * Dws1[b_Q_tmp]);
    }
    d1 = c_rotr[i];
    d2 = c_rotr[i + 3];
    d3 = c_rotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      c_a[b_i] = -((e_a[b_i] + tau0_dot_tau[b_i]) +
                   ((d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
                    d3 * RJdexpT_[3 * i1 + 2]));
    }
    d1 = ws2hat[i];
    d2 = ws2hat[i + 3];
    d3 = ws2hat[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      e_a[b_i] = (tauhat[i] * Dws2[3 * i1] + tauhat[i + 3] * Dws2[Q_tmp]) +
                 tauhat[i + 6] * Dws2[b_Q_tmp];
      c_rotr[b_i] =
          (d1 * tauhat[3 * i1] + d2 * tauhat[Q_tmp]) + d3 * tauhat[b_Q_tmp];
    }
    d1 = c_rotr[i];
    d2 = c_rotr[i + 3];
    d3 = c_rotr[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      tau0_dot_tau[i + 3 * i1] =
          (d1 * RJdexpT_[3 * i1] + d2 * RJdexpT_[3 * i1 + 1]) +
          d3 * RJdexpT_[3 * i1 + 2];
    }
  }
  d1 = qbar[6];
  d2 = rho1[6];
  d3 = rho2[6];
  d4 = rho2[3];
  d5 = rho2[4];
  d6 = rho2[5];
  d7 = u[0];
  d8 = u[1];
  for (i = 0; i < 3; i++) {
    b_i = 7 * (i + 3);
    C[b_i + 3] = ((hphi[3 * i] + d1 * c_a[3 * i]) +
                  d2 * -(e_a[3 * i] + tau0_dot_tau[3 * i])) +
                 d3 * D2tau1[3 * i];
    d9 = tau0_cross_tau[i];
    b_a[3 * i] = scale * (double)iv[3 * i] + d7 * d9;
    b_Q_tmp = 3 * i + 1;
    C[b_i + 4] = ((hphi[b_Q_tmp] + d1 * c_a[b_Q_tmp]) +
                  d2 * -(e_a[b_Q_tmp] + tau0_dot_tau[b_Q_tmp])) +
                 d3 * D2tau1[b_Q_tmp];
    b_a[b_Q_tmp] = scale * (double)iv[b_Q_tmp] + d8 * d9;
    b_Q_tmp = 3 * i + 2;
    C[b_i + 5] = ((hphi[b_Q_tmp] + d1 * c_a[b_Q_tmp]) +
                  d2 * -(e_a[b_Q_tmp] + tau0_dot_tau[b_Q_tmp])) +
                 d3 * D2tau1[b_Q_tmp];
    b_a[b_Q_tmp] = scale * (double)iv[b_Q_tmp] + d * d9;
    b_rotr[i] =
        (Dv_tmp[i] * tau0_cross_tau[0] + Dv_tmp[i + 3] * tau0_cross_tau[1]) +
        Dv_tmp[i + 6] * tau0_cross_tau[2];
    m[i] = (Dv[i] * d4 + Dv[i + 3] * d5) + Dv[i + 6] * d6;
  }
  for (i = 0; i < 3; i++) {
    d = b_a[i];
    d1 = b_a[i + 3];
    d2 = b_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      e_a[i1 + 3 * i] = b_rotr[i1] * u[i];
      c_a[i + 3 * i1] =
          (d * PP[3 * i1] + d1 * PP[3 * i1 + 1]) + d2 * PP[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    e_a[i] = d_a * (e_a[i] + c_a[i]);
  }
  for (i = 0; i < 3; i++) {
    d = ws2hat[i];
    d1 = ws2hat[i + 3];
    d2 = ws2hat[i + 6];
    d3 = a_tmp[i];
    d4 = a_tmp[i + 3];
    d5 = a_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_rotr[b_i] = (d * Dtau[3 * i1] + d1 * Dtau[Q_tmp]) + d2 * Dtau[b_Q_tmp];
      tau0_dot_tau[b_i] =
          (d3 * Dws2[3 * i1] + d4 * Dws2[Q_tmp]) + d5 * Dws2[b_Q_tmp];
    }
  }
  for (i = 0; i < 9; i++) {
    tau0_dot_tau[i] += c_rotr[i];
  }
  for (i = 0; i < 3; i++) {
    d = Dv_tmp[i];
    d1 = Dv_tmp[i + 3];
    d2 = Dv_tmp[i + 6];
    d3 = e_a[i];
    d4 = e_a[i + 3];
    d5 = e_a[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      Q_tmp = 3 * i1 + 1;
      b_Q_tmp = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_a[b_i] = (c_a_tmp * d * tau0_dot_tau[3 * i1] +
                  c_a_tmp * d1 * tau0_dot_tau[Q_tmp]) +
                 c_a_tmp * d2 * tau0_dot_tau[b_Q_tmp];
      b_a[b_i] = (d3 * Dtau[3 * i1] + d4 * Dtau[Q_tmp]) + d5 * Dtau[b_Q_tmp];
    }
  }
  for (i = 0; i < 9; i++) {
    b_a[i] += c_a[i];
  }
  d = m[0];
  d1 = m[1];
  d2 = m[2];
  d3 = ws1[0];
  d4 = ws1[1];
  d5 = ws1[2];
  d6 = Thet_tilde_e2[0];
  d7 = Thet_tilde_e2[1];
  d8 = Thet_tilde_e2[2];
  d9 = v[0];
  d10 = v[1];
  d11 = v[2];
  scale = rho2[3];
  absxk = rho2[4];
  t = rho2[5];
  for (i = 0; i < 3; i++) {
    i1 = 3 * i + 1;
    Q_tmp = 3 * i + 2;
    C[7 * (i + 3) + 6] =
        ((((d * Dws1[3 * i] + d1 * Dws1[i1]) + d2 * Dws1[Q_tmp]) +
          ((d3 * b_a[3 * i] + d4 * b_a[i1]) + d5 * b_a[Q_tmp])) +
         ((d6 * Dv[3 * i] + d7 * Dv[i1]) + d8 * Dv[Q_tmp])) +
        ((d9 * D2ws1[3 * i] + d10 * D2ws1[i1]) + d11 * D2ws1[Q_tmp]);
    C[i + 45] = (D2tau1[i] * scale + D2tau1[i + 3] * absxk) + D2tau1[i + 6] * t;
  }
}

void e_CableBCTransinCoord(const double qbar[7], const double x0[3],
                           const double RJ[9], const double RE[9],
                           const double r[3], const double R0[9], double q[7])
{
  double b_hphi[9];
  double dv[9];
  double hphi[9];
  double tau0_dot_tau[9];
  double tau0_cross_tau[3];
  double absxk;
  double d;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double nphi;
  double scale;
  double t;
  int i;
  int tau0_dot_tau_tmp;
  /*  Computes the transformation from (d, phi, gamma) -> (p1, p2, vartheta) for
   */
  /*  one end of the cable; NOTE: the main difference from CableBCtrans.m is */
  /*  that this function is based on exponential coordinates for the rotation,
   */
  /*  rather than the rotation matrix; also, a joint coordinate system and end
   */
  /*  offset are included, making it a more general. */
  /*  Inputs */
  /*  qbar = [d; phi; gamm] */
  /*    where d = displacement of joint,  */
  /*          phi = exponential coordinates of rotation of joint */
  /*          gamm = distance between first and second control points */
  /*  x0 = reference position of joint (3x1) vector */
  /*  RJ = rotation of joint coordinate frame with respect to global */
  /*  RE = rotation of cable end with respect to joint coordinate frame */
  /*  r = position of cable end relative to joint in joint coordinate system */
  /*      (end offset, 3x1 vector) */
  /*  R0 = orientation of cable end in reference configuration */
  /*  eta = vector to contract with for second derivative, Q */
  /*  rho1 = vector to contract with for second derivative, Qtilde */
  /*  rho2 = vector to contract with for third derivative, C */
  /*  Outputs: */
  /*  q = [p1; p2; vartheta1] */
  /*    where p1, p2 = positions of first and second control points */
  /*          vartheta1 = (twist) rotation of first control point */
  /*  J = first derivative of transformation */
  /*  Q = second derivative contracted over upper index,  */
  /*      Q_{ib,jb} = eta_i q^i_{ib,jb} */
  /*  Qtilde = second derivative contracted over one of the lower indices */
  /*      Qtilde^{i}_{ib} = rho1^{jb} q^i_{ib,jb} */
  /*  C = third derivative contracted over two lower indices */
  /*      C^{i}_{ib} = q^i_{ib,jb,kb}rho1^{jb}rho2^{kb} */
  /*  Computes helper functions associated with rotation matrix */
  /*  phi = exponential coordinates of rotation */
  /*  if flag == 0, compute a(1:2) */
  /*  if flag == 1, compute a(1:4,9) */
  /*  if flag == 2, compute a(1:4,9:10) */
  /*  if flag == 3, compute a(1:11) */
  scale = 3.3121686421112381E-170;
  absxk = fabs(qbar[3]);
  if (absxk > 3.3121686421112381E-170) {
    nphi = 1.0;
    scale = absxk;
  } else {
    t = absxk / 3.3121686421112381E-170;
    nphi = t * t;
  }
  absxk = fabs(qbar[4]);
  if (absxk > scale) {
    t = scale / absxk;
    nphi = nphi * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nphi += t * t;
  }
  absxk = fabs(qbar[5]);
  if (absxk > scale) {
    t = scale / absxk;
    nphi = nphi * t * t + 1.0;
    scale = absxk;
  } else {
    t = absxk / scale;
    nphi += t * t;
  }
  nphi = scale * sqrt(nphi);
  if (nphi > 1.0E-5) {
    absxk = sin(nphi) / nphi;
    scale = (1.0 - cos(nphi)) / (nphi * nphi);
  } else {
    absxk = nphi * (nphi * (nphi * (nphi * 0.0083333333333333332) -
                            0.16666666666666666)) +
            1.0;
    scale = nphi * (nphi * (nphi * (nphi * 0.0013888888888888889) -
                            0.041666666666666664)) +
            0.5;
  }
  hphi[0] = 0.0;
  hphi[3] = -qbar[5];
  hphi[6] = qbar[4];
  hphi[1] = qbar[5];
  hphi[4] = 0.0;
  hphi[7] = -qbar[3];
  hphi[2] = -qbar[4];
  hphi[5] = qbar[3];
  hphi[8] = 0.0;
  for (i = 0; i < 3; i++) {
    for (tau0_dot_tau_tmp = 0; tau0_dot_tau_tmp < 3; tau0_dot_tau_tmp++) {
      b_hphi[i + 3 * tau0_dot_tau_tmp] =
          (hphi[i] * hphi[3 * tau0_dot_tau_tmp] +
           hphi[i + 3] * hphi[3 * tau0_dot_tau_tmp + 1]) +
          hphi[i + 6] * hphi[3 * tau0_dot_tau_tmp + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    hphi[i] = (absxk * hphi[i] + (double)iv[i]) + scale * b_hphi[i];
  }
  for (i = 0; i < 3; i++) {
    d = RJ[i];
    d1 = RJ[i + 3];
    d2 = RJ[i + 6];
    for (tau0_dot_tau_tmp = 0; tau0_dot_tau_tmp < 3; tau0_dot_tau_tmp++) {
      b_hphi[i + 3 * tau0_dot_tau_tmp] = (d * hphi[3 * tau0_dot_tau_tmp] +
                                          d1 * hphi[3 * tau0_dot_tau_tmp + 1]) +
                                         d2 * hphi[3 * tau0_dot_tau_tmp + 2];
    }
  }
  scale = 0.0;
  for (i = 0; i < 3; i++) {
    d = b_hphi[i];
    d1 = b_hphi[i + 3];
    d2 = b_hphi[i + 6];
    for (tau0_dot_tau_tmp = 0; tau0_dot_tau_tmp < 3; tau0_dot_tau_tmp++) {
      hphi[i + 3 * tau0_dot_tau_tmp] =
          (d * RE[3 * tau0_dot_tau_tmp] + d1 * RE[3 * tau0_dot_tau_tmp + 1]) +
          d2 * RE[3 * tau0_dot_tau_tmp + 2];
    }
    scale += R0[i] * hphi[i];
  }
  tau0_cross_tau[0] = R0[1] * hphi[2] - hphi[1] * R0[2];
  tau0_cross_tau[1] = hphi[0] * R0[2] - R0[0] * hphi[2];
  tau0_cross_tau[2] = R0[0] * hphi[1] - hphi[0] * R0[1];
  dv[0] = 0.0;
  dv[1] = -tau0_cross_tau[2];
  dv[2] = tau0_cross_tau[1];
  dv[3] = tau0_cross_tau[2];
  dv[4] = 0.0;
  dv[5] = -tau0_cross_tau[0];
  dv[6] = -tau0_cross_tau[1];
  dv[7] = tau0_cross_tau[0];
  dv[8] = 0.0;
  for (i = 0; i < 3; i++) {
    d = tau0_cross_tau[i] / (scale + 1.0);
    tau0_dot_tau[3 * i] =
        (scale * (double)iv[i] + dv[3 * i]) + d * tau0_cross_tau[0];
    tau0_dot_tau_tmp = 3 * i + 1;
    tau0_dot_tau[tau0_dot_tau_tmp] =
        (scale * (double)iv[i + 3] + dv[tau0_dot_tau_tmp]) +
        d * tau0_cross_tau[1];
    tau0_dot_tau_tmp = 3 * i + 2;
    tau0_dot_tau[tau0_dot_tau_tmp] =
        (scale * (double)iv[i + 6] + dv[tau0_dot_tau_tmp]) +
        d * tau0_cross_tau[2];
  }
  d = hphi[3];
  d1 = hphi[4];
  d2 = hphi[5];
  for (i = 0; i < 3; i++) {
    tau0_cross_tau[i] = (tau0_dot_tau[i] * d + tau0_dot_tau[i + 3] * d1) +
                        tau0_dot_tau[i + 6] * d2;
  }
  d = tau0_cross_tau[0];
  d1 = tau0_cross_tau[1];
  d2 = tau0_cross_tau[2];
  scale = qbar[0];
  absxk = qbar[1];
  t = qbar[2];
  nphi = qbar[6];
  d3 = r[0];
  d4 = r[1];
  d5 = r[2];
  for (i = 0; i < 3; i++) {
    double d6;
    tau0_cross_tau[i] =
        (R0[3 * i] * d + R0[3 * i + 1] * d1) + R0[3 * i + 2] * d2;
    d6 = (x0[i] + ((RJ[i] * scale + RJ[i + 3] * absxk) + RJ[i + 6] * t)) +
         ((b_hphi[i] * d3 + b_hphi[i + 3] * d4) + b_hphi[i + 6] * d5);
    q[i] = d6;
    q[i + 3] = d6 + nphi * hphi[i];
  }
  q[6] = atan2(tau0_cross_tau[2], tau0_cross_tau[1]);
}

/* End of code generation (CableBCTransinCoord.c) */
