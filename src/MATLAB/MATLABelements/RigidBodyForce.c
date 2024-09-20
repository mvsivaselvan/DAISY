/*
 * Academic License - for use in teaching, academic research, and meeting
 * course requirements at degree granting institutions only.  Not for
 * government, commercial, or other organizational use.
 *
 * RigidBodyForce.c
 *
 * Code generation for function 'RigidBodyForce'
 *
 */

/* Include files */
#include "RigidBodyForce.h"
#include "CableForceRotBCinCoord_data.h"
#include "CableForceRotBCinCoord_initialize.h"
#include "rothelper.h"

/* Function Definitions */
void RigidBodyForce(const double d[3], const double phi[3], const double dd[3],
                    const double phid[3], const double ddd[3],
                    const double phidd[3], double m, const double II[9],
                    const double KT[9], const double CT[9], const double KR[9],
                    const double CR[9], const double rr[3], const double RJ[9],
                    const double u[3], double F[6], double K[36], double C[36],
                    double M[36], double B[18])
{
  double a_contents[11];
  double dv[11];
  double D1w_phid[9];
  double F_tmp[9];
  double R[9];
  double b_F_tmp[9];
  double b_a_contents[9];
  double b_tmp[9];
  double c_a_contents[9];
  double c_phi[9];
  double dv1[9];
  double dv2[9];
  double dv3[9];
  double hphi[9];
  double hphi2[9];
  double rhat[9];
  double tmp1[9];
  double tmp2[9];
  double tmp3[9];
  double tmp4[9];
  double what[9];
  double wtBody_tmp[9];
  double y_tmp[9];
  double RTddd[3];
  double b_w[3];
  double mubar[3];
  double mubar_tmp[3];
  double w[3];
  double wd[3];
  double wmubar[3];
  double wtBody[3];
  double wtJoint[3];
  double a;
  double a_tmp;
  double abar_idx_0;
  double abar_idx_1;
  double abar_idx_2;
  double b_a;
  double b_d;
  double b_phi;
  double b_phid;
  double b_w_tmp;
  double c_a;
  double c_w_tmp;
  double d1;
  double d2;
  double d3;
  double d4;
  double d5;
  double d6;
  double d7;
  double d8;
  double d_phi;
  double phi_tmp;
  double w_tmp;
  double wxr_idx_0;
  double wxr_idx_1;
  double wxr_idx_2;
  int a_contents_tmp_tmp;
  int b_i;
  int i;
  int i1;
  int i5;
  signed char i2;
  if (!isInitialized_CableForceRotBCinCoord) {
    CableForceRotBCinCoord_initialize();
  }
  /*  Inputs */
  /*  d = displacement of joint in joint coordinate system */
  /*  phi = exponential coordinates of rotation in joint coordinate system */
  /*  dd = \dot{d} */
  /*  phid = \dot{phi} */
  /*  ddd = \ddot(d} */
  /*  phidd = \ddot(phi} */
  /*  m = mass */
  /*  II = moment of inertia about joint */
  /*  KT = translational stiffness matrix (3x3) */
  /*  CT = translational damping matrix (3x3) */
  /*  KR = rotational stiffness matrix (3x3) */
  /*  CR = rotational damping matrix (3x3) */
  /*  rr = position vector of center of mass relative to joint in joint */
  /*       coordinate system */
  /*  RJ = oritentation of joint coordinate system relative to global */
  /*  u = force per unit mass (acceleration due to gravity and ground)  */
  /*  Outputs */
  /*  F = joint force (6x1 - force and moment) */
  /*  K = stiffness matrix (derivative of force w.r.t. (d,phi)) */
  /*  C = derivative of force w.r.t. (ddot,phidot) */
  /*  M = mass matrix (derivative of force w.r.t. (dddot,phiddot)) */
  /*  B = derivative of force w.r.t. input */
  b_rothelper(phi, dv);
  b_rothelper(phi, a_contents);
  hphi[0] = 0.0;
  hphi[3] = -phi[2];
  hphi[6] = phi[1];
  hphi[1] = phi[2];
  hphi[4] = 0.0;
  hphi[7] = -phi[0];
  hphi[2] = -phi[1];
  hphi[5] = phi[0];
  hphi[8] = 0.0;
  for (i = 0; i < 3; i++) {
    for (i1 = 0; i1 < 3; i1++) {
      hphi2[i + 3 * i1] =
          (hphi[i] * hphi[3 * i1] + hphi[i + 3] * hphi[3 * i1 + 1]) +
          hphi[i + 6] * hphi[3 * i1 + 2];
    }
  }
  b_d = a_contents[0];
  d1 = a_contents[1];
  d2 = a_contents[8];
  for (i = 0; i < 9; i++) {
    b_i = iv[i];
    d3 = hphi[i];
    d4 = hphi2[i];
    R[i] = ((double)b_i + b_d * d3) + d1 * d4;
    d3 = ((double)b_i - d1 * d3) + d2 * d4;
    hphi[i] = d3;
  }
  for (i = 0; i < 3; i++) {
    y_tmp[3 * i] = RJ[i];
    y_tmp[3 * i + 1] = RJ[i + 3];
    y_tmp[3 * i + 2] = RJ[i + 6];
  }
  /*  Weight in the joint frame */
  b_d = u[0];
  d1 = u[1];
  d2 = u[2];
  for (i = 0; i < 3; i++) {
    wtBody_tmp[3 * i] = R[i];
    wtBody_tmp[3 * i + 1] = R[i + 3];
    wtBody_tmp[3 * i + 2] = R[i + 6];
    wtJoint[i] = m * ((y_tmp[i] * b_d + y_tmp[i + 3] * d1) + y_tmp[i + 6] * d2);
  }
  /*  weight in the body frame */
  b_phi = 0.0;
  b_tmp[0] = 0.0;
  b_tmp[3] = -phid[2];
  b_tmp[6] = phid[1];
  b_tmp[1] = phid[2];
  b_tmp[4] = 0.0;
  b_tmp[7] = -phid[0];
  b_tmp[2] = -phid[1];
  b_tmp[5] = phid[0];
  b_tmp[8] = 0.0;
  w_tmp = phid[1] * phi[2] - phi[1] * phid[2];
  b_w_tmp = phi[0] * phid[2] - phid[0] * phi[2];
  c_w_tmp = phid[0] * phi[1] - phi[0] * phid[1];
  b_d = dv[2];
  d1 = dv[3] * w_tmp;
  d2 = dv[3] * b_w_tmp;
  d3 = dv[3] * c_w_tmp;
  d4 = wtJoint[0];
  d5 = wtJoint[1];
  d6 = wtJoint[2];
  for (i = 0; i < 3; i++) {
    d7 = phi[i];
    b_phi += d7 * phid[i];
    dv1[3 * i] = phid[0] * b_d * d7;
    dv2[3 * i] = d1 * d7;
    c_phi[3 * i] = phi[0] * phi[i];
    i1 = 3 * i + 1;
    d7 = phi[i];
    dv1[i1] = phid[1] * b_d * d7;
    dv2[i1] = d2 * d7;
    c_phi[i1] = phi[1] * phi[i];
    i1 = 3 * i + 2;
    d7 = phi[i];
    dv1[i1] = b_d * phid[2] * d7;
    dv2[i1] = d3 * d7;
    c_phi[i1] = phi[2] * phi[i];
    w[i] = (hphi[i] * phid[0] + hphi[i + 3] * phid[1]) + hphi[i + 6] * phid[2];
    wtBody[i] =
        (wtBody_tmp[i] * d4 + wtBody_tmp[i + 3] * d5) + wtBody_tmp[i + 6] * d6;
  }
  a = dv[9] * b_phi;
  b_a = dv[8] * b_phi;
  b_d = a_contents[1];
  d1 = phi[0];
  d2 = phi[1];
  d3 = phi[2];
  d4 = dv[8];
  d5 = phidd[0];
  d6 = phidd[1];
  d7 = phidd[2];
  for (i = 0; i < 3; i++) {
    d8 = phid[i];
    D1w_phid[3 * i] =
        ((((dv1[3 * i] + dv2[3 * i]) + b_d * b_tmp[3 * i]) + a * c_phi[3 * i]) +
         d1 * d4 * d8) +
        b_a * (double)iv[3 * i];
    b_i = 3 * i + 1;
    D1w_phid[b_i] =
        ((((dv1[b_i] + dv2[b_i]) + b_d * b_tmp[b_i]) + a * c_phi[b_i]) +
         d2 * d4 * d8) +
        b_a * (double)iv[b_i];
    b_i = 3 * i + 2;
    D1w_phid[b_i] =
        ((((dv1[b_i] + dv2[b_i]) + b_d * b_tmp[b_i]) + a * c_phi[b_i]) +
         d3 * d4 * d8) +
        b_a * (double)iv[b_i];
    wd[i] = (hphi[i] * d5 + hphi[i + 3] * d6) + hphi[i + 6] * d7;
  }
  b_d = phid[0];
  d1 = phid[1];
  d2 = phid[2];
  for (i = 0; i < 3; i++) {
    wd[i] += (D1w_phid[i] * b_d + D1w_phid[i + 3] * d1) + D1w_phid[i + 6] * d2;
  }
  wxr_idx_0 = w[1] * rr[2] - rr[1] * w[2];
  wxr_idx_1 = rr[0] * w[2] - w[0] * rr[2];
  wxr_idx_2 = w[0] * rr[1] - rr[0] * w[1];
  abar_idx_0 =
      (wd[1] * rr[2] - rr[1] * wd[2]) + (w[1] * wxr_idx_2 - wxr_idx_1 * w[2]);
  abar_idx_1 =
      (rr[0] * wd[2] - wd[0] * rr[2]) + (wxr_idx_0 * w[2] - w[0] * wxr_idx_2);
  abar_idx_2 =
      (wd[0] * rr[1] - rr[0] * wd[1]) + (w[0] * wxr_idx_1 - wxr_idx_0 * w[1]);
  b_d = ddd[0];
  d1 = ddd[1];
  d2 = ddd[2];
  d3 = w[0];
  d4 = w[1];
  d5 = w[2];
  for (i = 0; i < 3; i++) {
    RTddd[i] =
        (wtBody_tmp[i] * b_d + wtBody_tmp[i + 3] * d1) + wtBody_tmp[i + 6] * d2;
    mubar_tmp[i] = (II[i] * d3 + II[i + 3] * d4) + II[i + 6] * d5;
  }
  mubar[0] = m * (rr[1] * RTddd[2] - RTddd[1] * rr[2]);
  mubar[1] = m * (RTddd[0] * rr[2] - rr[0] * RTddd[2]);
  mubar[2] = m * (rr[0] * RTddd[1] - RTddd[0] * rr[1]);
  b_w[0] = w[1] * mubar_tmp[2] - mubar_tmp[1] * w[2];
  b_w[1] = mubar_tmp[0] * w[2] - w[0] * mubar_tmp[2];
  b_w[2] = w[0] * mubar_tmp[1] - mubar_tmp[0] * w[1];
  wmubar[0] = rr[1] * wtBody[2] - wtBody[1] * rr[2];
  wmubar[1] = wtBody[0] * rr[2] - rr[0] * wtBody[2];
  wmubar[2] = rr[0] * wtBody[1] - wtBody[0] * rr[1];
  b_d = wd[0];
  d1 = wd[1];
  d2 = wd[2];
  for (i = 0; i < 3; i++) {
    mubar[i] = ((mubar[i] + ((II[i] * b_d + II[i + 3] * d1) + II[i + 6] * d2)) +
                b_w[i]) -
               wmubar[i];
  }
  b_d = mubar[0];
  d1 = mubar[1];
  d2 = mubar[2];
  for (i = 0; i < 3; i++) {
    wmubar[i] = (hphi[i] * b_d + hphi[i + 3] * d1) + hphi[i + 6] * d2;
  }
  for (b_i = 0; b_i < 6; b_i++) {
    F[b_i] = 0.0;
  }
  for (i = 0; i < 9; i++) {
    F_tmp[i] = m * R[i];
  }
  /*  need derivatives K, C and M */
  rhat[0] = 0.0;
  rhat[3] = -rr[2];
  rhat[6] = rr[1];
  rhat[1] = rr[2];
  rhat[4] = 0.0;
  rhat[7] = -rr[0];
  rhat[2] = -rr[1];
  rhat[5] = rr[0];
  rhat[8] = 0.0;
  what[0] = 0.0;
  what[3] = -w[2];
  what[6] = w[1];
  what[1] = w[2];
  what[4] = 0.0;
  what[7] = -w[0];
  what[2] = -w[1];
  what[5] = w[0];
  what[8] = 0.0;
  /*  ocmputes matrix D_1^2 w(phi,v)(dphi1,.) */
  phi_tmp = 0.0;
  b_phid = 0.0;
  wd[0] = w_tmp;
  wd[1] = b_w_tmp;
  wd[2] = c_w_tmp;
  b_phi = 0.0;
  d_phi = 0.0;
  b_d = d[0];
  d1 = d[1];
  d2 = d[2];
  d3 = dd[0];
  d4 = dd[1];
  d5 = dd[2];
  d6 = wmubar[0];
  d7 = wmubar[1];
  d8 = wmubar[2];
  for (i = 0; i < 3; i++) {
    F[i] =
        (((m * ddd[i] + ((F_tmp[i] * abar_idx_0 + F_tmp[i + 3] * abar_idx_1) +
                         F_tmp[i + 6] * abar_idx_2)) +
          ((KT[i] * b_d + KT[i + 3] * d1) + KT[i + 6] * d2)) +
         ((CT[i] * d3 + CT[i + 3] * d4) + CT[i + 6] * d5)) -
        wtJoint[i];
    phi_tmp += phi[i] * phid[i];
    hphi2[3 * i] = phi[0] * phi[i];
    hphi2[3 * i + 1] = phi[1] * phi[i];
    hphi2[3 * i + 2] = phi[2] * phi[i];
    F[i + 3] = (((R[i] * d6 + R[i + 3] * d7) + R[i + 6] * d8) +
                ((KR[i] * phi[0] + KR[i + 3] * phi[1]) + KR[i + 6] * phi[2])) +
               ((CR[i] * phid[0] + CR[i + 3] * phid[1]) + CR[i + 6] * phid[2]);
    c_a = phid[i];
    b_phid += c_a * c_a;
    b_phi = phi_tmp;
    d_phi += phi[i] * phidd[i];
  }
  a_tmp = dv[9] * phi_tmp;
  a = dv[10] * phi_tmp;
  b_a = dv[9] * b_phid;
  phi_tmp = dv[8] * b_phid;
  b_phid = dv[9] * d_phi;
  c_a = dv[8] * d_phi;
  b_d = a_contents[2];
  d1 = a_contents[3];
  d2 = a_contents[4];
  d3 = a_contents[5];
  d4 = a_contents[9];
  for (i = 0; i < 3; i++) {
    b_w[i] = (b_d * phid[i] + d1 * wd[i]) + a_tmp * phi[i];
    c_a_contents[3 * i] =
        ((((phid[0] * d2 * phi[i] + w_tmp * d3 * phi[i]) + d1 * b_tmp[3 * i]) +
          a * hphi2[3 * i]) +
         phi[0] * d4 * phid[i]) +
        a_tmp * (double)iv[3 * i];
    b_i = 3 * i + 1;
    c_a_contents[b_i] =
        ((((phid[1] * d2 * phi[i] + b_w_tmp * d3 * phi[i]) + d1 * b_tmp[b_i]) +
          a * hphi2[b_i]) +
         phi[1] * d4 * phid[i]) +
        a_tmp * (double)iv[b_i];
    b_i = 3 * i + 2;
    c_a_contents[b_i] =
        ((((phid[2] * d2 * phi[i] + c_w_tmp * d3 * phi[i]) + d1 * b_tmp[b_i]) +
          a * hphi2[b_i]) +
         phi[2] * d4 * phid[i]) +
        a_tmp * (double)iv[b_i];
  }
  b_a_contents[0] = 0.0;
  b_a_contents[3] = dv[1] * -phidd[2];
  b_a_contents[6] = dv[1] * phidd[1];
  b_a_contents[1] = dv[1] * phidd[2];
  b_a_contents[4] = 0.0;
  b_a_contents[7] = dv[1] * -phidd[0];
  b_a_contents[2] = dv[1] * -phidd[1];
  b_a_contents[5] = phidd[0] * dv[1];
  b_a_contents[8] = 0.0;
  b_d = b_w[0];
  d1 = b_w[1];
  d2 = b_w[2];
  d3 = dv[8];
  d4 = dv[2];
  d5 = dv[3] * (phidd[1] * phi[2] - phi[1] * phidd[2]);
  d6 = dv[3] * (phi[0] * phidd[2] - phidd[0] * phi[2]);
  d7 = dv[3] * (phidd[0] * phi[1] - phi[0] * phidd[1]);
  for (i = 0; i < 3; i++) {
    signed char i3;
    signed char i4;
    d8 = phi[i];
    i2 = iv[3 * i];
    b_i = 3 * i + 1;
    i3 = iv[b_i];
    a_contents_tmp_tmp = 3 * i + 2;
    i4 = iv[a_contents_tmp_tmp];
    tmp1[3 * i] =
        (((((b_d * phid[i] + b_phi * c_a_contents[3 * i]) +
            b_a * hphi2[3 * i]) +
           phi_tmp * (double)i2) +
          a_tmp * phid[0] * d8) +
         phid[0] * d3 * phid[i]) +
        (((((phidd[0] * d4 * phi[i] + d5 * phi[i]) + b_a_contents[3 * i]) +
           b_phid * (phi[0] * phi[i])) +
          phi[0] * d3 * phidd[i]) +
         c_a * (double)i2);
    tmp1[b_i] =
        (((((d1 * phid[i] + b_phi * c_a_contents[b_i]) + b_a * hphi2[b_i]) +
           phi_tmp * (double)i3) +
          a_tmp * phid[1] * d8) +
         phid[1] * d3 * phid[i]) +
        (((((phidd[1] * d4 * phi[i] + d6 * phi[i]) + b_a_contents[b_i]) +
           b_phid * (phi[1] * phi[i])) +
          phi[1] * d3 * phidd[i]) +
         c_a * (double)i3);
    tmp1[a_contents_tmp_tmp] =
        (((((d2 * phid[i] + b_phi * c_a_contents[a_contents_tmp_tmp]) +
            b_a * hphi2[a_contents_tmp_tmp]) +
           phi_tmp * (double)i4) +
          a_tmp * phid[2] * d8) +
         phid[2] * d3 * phid[i]) +
        (((((d4 * phidd[2] * phi[i] + d7 * phi[i]) +
            b_a_contents[a_contents_tmp_tmp]) +
           b_phid * (phi[2] * phi[i])) +
          phi[2] * d3 * phidd[i]) +
         c_a * (double)i4);
  }
  tmp2[0] = 0.0;
  tmp2[3] = -wxr_idx_2;
  tmp2[6] = wxr_idx_1;
  tmp2[1] = wxr_idx_2;
  tmp2[4] = 0.0;
  tmp2[7] = -wxr_idx_0;
  tmp2[2] = -wxr_idx_1;
  tmp2[5] = wxr_idx_0;
  tmp2[8] = 0.0;
  for (i = 0; i < 3; i++) {
    b_d = what[i];
    d1 = what[i + 3];
    d2 = what[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      c_a_contents[i + 3 * i1] =
          (b_d * rhat[3 * i1] + d1 * rhat[3 * i1 + 1]) + d2 * rhat[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    tmp2[i] += c_a_contents[i];
  }
  /*  computes D_2D_1 w(phi,~)(dphi1,.) */
  b_phi = 0.0;
  tmp3[0] = 0.0;
  tmp3[3] = a_contents[3] * -phi[2];
  tmp3[6] = phi[1] * a_contents[3];
  tmp3[1] = phi[2] * a_contents[3];
  tmp3[4] = 0.0;
  tmp3[7] = a_contents[3] * -phi[0];
  tmp3[2] = a_contents[3] * -phi[1];
  tmp3[5] = phi[0] * a_contents[3];
  tmp3[8] = 0.0;
  for (i = 0; i < 3; i++) {
    b_d = phid[i];
    b_phi += phi[i] * b_d;
    c_phi[3 * i] = phi[0] * phi[i];
    c_a_contents[3 * i] = phi[0] * b_d;
    b_a_contents[3 * i] = phid[0] * phi[i];
    b_i = 3 * i + 1;
    c_phi[b_i] = phi[1] * phi[i];
    c_a_contents[b_i] = phi[1] * phid[i];
    b_a_contents[b_i] = phid[1] * phi[i];
    b_i = 3 * i + 2;
    c_phi[b_i] = phi[2] * phi[i];
    c_a_contents[b_i] = phi[2] * phid[i];
    b_a_contents[b_i] = phid[2] * phi[i];
  }
  b_d = a_contents[2];
  d1 = a_contents[9];
  d2 = dv[1];
  d3 = dv[8];
  for (i = 0; i < 9; i++) {
    tmp3[i] = D1w_phid[i] +
              ((b_phi * ((b_d * (double)iv[i] - tmp3[i]) + d1 * c_phi[i]) -
                d2 * b_tmp[i]) +
               d3 * (c_a_contents[i] + b_a_contents[i]));
  }
  tmp4[0] = -0.0;
  tmp4[3] = mubar_tmp[2];
  tmp4[6] = -mubar_tmp[1];
  tmp4[1] = -mubar_tmp[2];
  tmp4[4] = -0.0;
  tmp4[7] = mubar_tmp[0];
  tmp4[2] = mubar_tmp[1];
  tmp4[5] = -mubar_tmp[0];
  tmp4[8] = -0.0;
  for (i = 0; i < 3; i++) {
    b_d = what[i];
    d1 = what[i + 3];
    d2 = what[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      c_a_contents[i + 3 * i1] =
          (b_d * II[3 * i1] + d1 * II[3 * i1 + 1]) + d2 * II[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    tmp4[i] += c_a_contents[i];
    hphi2[i] = -rhat[i];
  }
  c_a = (rr[0] * RTddd[0] + rr[1] * RTddd[1]) + rr[2] * RTddd[2];
  phi_tmp = (rr[0] * wtBody[0] + rr[1] * wtBody[1]) + rr[2] * wtBody[2];
  b_phi = 0.0;
  dv1[0] = -0.0;
  dv1[3] = abar_idx_2;
  dv1[6] = -abar_idx_1;
  dv1[1] = -abar_idx_2;
  dv1[4] = -0.0;
  dv1[7] = abar_idx_0;
  dv1[2] = abar_idx_1;
  dv1[5] = -abar_idx_0;
  dv1[8] = -0.0;
  for (i = 0; i < 3; i++) {
    b_tmp[3 * i] = hphi[i];
    b_tmp[3 * i + 1] = hphi[i + 3];
    b_tmp[3 * i + 2] = hphi[i + 6];
    b_phi += phi[i] * mubar[i];
    b_d = dv1[i];
    d1 = dv1[i + 3];
    d2 = dv1[i + 6];
    d3 = hphi2[i];
    d4 = hphi2[i + 3];
    d5 = hphi2[i + 6];
    d6 = tmp2[i];
    d7 = tmp2[i + 3];
    d8 = tmp2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      a_contents_tmp_tmp = 3 * i1 + 1;
      i5 = 3 * i1 + 2;
      dv2[i + 3 * i1] =
          ((b_d * hphi[3 * i1] + d1 * hphi[a_contents_tmp_tmp]) +
           d2 * hphi[i5]) +
          (((d3 * tmp1[3 * i1] + d4 * tmp1[a_contents_tmp_tmp]) +
            d5 * tmp1[i5]) -
           ((d6 * D1w_phid[3 * i1] + d7 * D1w_phid[a_contents_tmp_tmp]) +
            d8 * D1w_phid[i5]));
    }
  }
  a = dv[9] * b_phi;
  b_a = dv[8] * b_phi;
  for (i = 0; i < 3; i++) {
    b_d = F_tmp[i];
    d1 = F_tmp[i + 3];
    d2 = F_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_F_tmp[i + 3 * i1] =
          (b_d * dv2[3 * i1] + d1 * dv2[3 * i1 + 1]) + d2 * dv2[3 * i1 + 2];
    }
  }
  dv1[0] = -0.0;
  dv1[3] = wmubar[2];
  dv1[6] = -wmubar[1];
  dv1[1] = -wmubar[2];
  dv1[4] = -0.0;
  dv1[7] = wmubar[0];
  dv1[2] = wmubar[1];
  dv1[5] = -wmubar[0];
  dv1[8] = -0.0;
  b_w[0] = dv[3] * (mubar[1] * phi[2] - phi[1] * mubar[2]);
  b_w[1] = dv[3] * (phi[0] * mubar[2] - mubar[0] * phi[2]);
  b_w[2] = dv[3] * (mubar[0] * phi[1] - phi[0] * mubar[1]);
  b_a_contents[0] = 0.0;
  b_a_contents[3] = a_contents[1] * -mubar[2];
  b_a_contents[6] = a_contents[1] * mubar[1];
  b_a_contents[1] = a_contents[1] * mubar[2];
  b_a_contents[4] = 0.0;
  b_a_contents[7] = a_contents[1] * -mubar[0];
  b_a_contents[2] = a_contents[1] * -mubar[1];
  b_a_contents[5] = mubar[0] * a_contents[1];
  b_a_contents[8] = 0.0;
  for (i = 0; i < 3; i++) {
    b_d = dv1[i];
    d1 = dv1[i + 3];
    d2 = dv1[i + 6];
    d3 = phi[i];
    for (i1 = 0; i1 < 3; i1++) {
      a_contents_tmp_tmp = i1 + 3 * i;
      dv2[a_contents_tmp_tmp] = dv[2] * mubar[i1] * d3;
      c_a_contents[a_contents_tmp_tmp] = b_w[i1] * d3;
      c_phi[a_contents_tmp_tmp] = phi[i1] * phi[i];
      dv3[i + 3 * i1] =
          (b_d * hphi[3 * i1] + d1 * hphi[3 * i1 + 1]) + d2 * hphi[3 * i1 + 2];
    }
  }
  b_d = phi[0];
  d1 = phi[1];
  d2 = phi[2];
  d3 = dv[8];
  d4 = RTddd[0];
  d5 = RTddd[1];
  d6 = RTddd[2];
  for (i = 0; i < 3; i++) {
    i2 = iv[3 * i];
    d7 = mubar[i];
    dv1[3 * i] = ((((dv2[3 * i] + c_a_contents[3 * i]) + b_a_contents[3 * i]) +
                   a * c_phi[3 * i]) +
                  b_d * d3 * d7) +
                 b_a * (double)i2;
    d8 = rr[i];
    c_a_contents[3 * i] = d4 * d8 - c_a * (double)i2;
    i1 = 3 * i + 1;
    i2 = iv[i1];
    dv1[i1] =
        ((((dv2[i1] + c_a_contents[i1]) + b_a_contents[i1]) + a * c_phi[i1]) +
         d1 * d3 * d7) +
        b_a * (double)i2;
    c_a_contents[i1] = d5 * d8 - c_a * (double)i2;
    i1 = 3 * i + 2;
    i2 = iv[i1];
    dv1[i1] =
        ((((dv2[i1] + c_a_contents[i1]) + b_a_contents[i1]) + a * c_phi[i1]) +
         d2 * d3 * d7) +
        b_a * (double)i2;
    c_a_contents[i1] = d6 * d8 - c_a * (double)i2;
  }
  for (i = 0; i < 9; i++) {
    c_a_contents[i] *= m;
  }
  for (i = 0; i < 3; i++) {
    b_d = c_a_contents[i];
    d1 = c_a_contents[i + 3];
    d2 = c_a_contents[i + 6];
    d3 = II[i];
    d4 = II[i + 3];
    d5 = II[i + 6];
    d6 = tmp4[i];
    d7 = tmp4[i + 3];
    d8 = tmp4[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i1 + 3 * i;
      b_a_contents[b_i] = phi_tmp * (double)iv[b_i] - wtBody[i1] * rr[i];
      a_contents_tmp_tmp = 3 * i1 + 1;
      i5 = 3 * i1 + 2;
      what[i + 3 * i1] =
          (((b_d * hphi[3 * i1] + d1 * hphi[a_contents_tmp_tmp]) +
            d2 * hphi[i5]) +
           ((d3 * tmp1[3 * i1] + d4 * tmp1[a_contents_tmp_tmp]) +
            d5 * tmp1[i5])) +
          ((d6 * D1w_phid[3 * i1] + d7 * D1w_phid[a_contents_tmp_tmp]) +
           d8 * D1w_phid[i5]);
    }
  }
  for (i = 0; i < 3; i++) {
    b_d = b_a_contents[i];
    d1 = b_a_contents[i + 3];
    d2 = b_a_contents[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      c_a_contents[i + 3 * i1] =
          (b_d * hphi[3 * i1] + d1 * hphi[3 * i1 + 1]) + d2 * hphi[3 * i1 + 2];
    }
  }
  for (i = 0; i < 9; i++) {
    what[i] += c_a_contents[i];
  }
  for (i = 0; i < 3; i++) {
    b_d = hphi[i];
    d1 = hphi[i + 3];
    d2 = hphi[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      a_contents_tmp_tmp = i + 3 * i1;
      dv2[a_contents_tmp_tmp] =
          (dv3[a_contents_tmp_tmp] + dv1[a_contents_tmp_tmp]) +
          ((b_d * what[3 * i1] + d1 * what[3 * i1 + 1]) +
           d2 * what[3 * i1 + 2]);
    }
  }
  for (i = 0; i < 3; i++) {
    b_d = R[i];
    d1 = R[i + 3];
    d2 = R[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      c_a_contents[b_i] =
          ((b_d * dv2[3 * i1] + d1 * dv2[3 * i1 + 1]) + d2 * dv2[3 * i1 + 2]) +
          KR[b_i];
      b_i = i1 + 3 * i;
      a_contents_tmp_tmp = i1 + 6 * i;
      K[a_contents_tmp_tmp] = KT[b_i];
      K[i1 + 6 * (i + 3)] = b_F_tmp[b_i];
      K[a_contents_tmp_tmp + 3] = 0.0;
    }
  }
  for (i = 0; i < 3; i++) {
    b_i = 6 * (i + 3);
    K[b_i + 3] = c_a_contents[3 * i];
    K[b_i + 4] = c_a_contents[3 * i + 1];
    K[b_i + 5] = c_a_contents[3 * i + 2];
  }
  for (i = 0; i < 3; i++) {
    b_d = tmp2[i];
    d1 = tmp2[i + 3];
    d2 = tmp2[i + 6];
    d3 = hphi2[i];
    d4 = hphi2[i + 3];
    d5 = hphi2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      a_contents_tmp_tmp = 3 * i1 + 1;
      i5 = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_a_contents[b_i] =
          (b_d * hphi[3 * i1] + d1 * hphi[a_contents_tmp_tmp]) + d2 * hphi[i5];
      b_a_contents[b_i] =
          (d3 * tmp3[3 * i1] + d4 * tmp3[a_contents_tmp_tmp]) + d5 * tmp3[i5];
    }
  }
  for (i = 0; i < 9; i++) {
    b_a_contents[i] -= c_a_contents[i];
  }
  for (i = 0; i < 3; i++) {
    b_d = tmp4[i];
    d1 = tmp4[i + 3];
    d2 = tmp4[i + 6];
    d3 = II[i];
    d4 = II[i + 3];
    d5 = II[i + 6];
    d6 = F_tmp[i];
    d7 = F_tmp[i + 3];
    d8 = F_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      a_contents_tmp_tmp = 3 * i1 + 1;
      i5 = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_a_contents[b_i] =
          (b_d * hphi[3 * i1] + d1 * hphi[a_contents_tmp_tmp]) + d2 * hphi[i5];
      c_phi[b_i] =
          (d3 * tmp3[3 * i1] + d4 * tmp3[a_contents_tmp_tmp]) + d5 * tmp3[i5];
      b_F_tmp[b_i] =
          (d6 * b_a_contents[3 * i1] + d7 * b_a_contents[a_contents_tmp_tmp]) +
          d8 * b_a_contents[i5];
    }
  }
  for (i = 0; i < 9; i++) {
    c_phi[i] += c_a_contents[i];
  }
  for (i = 0; i < 3; i++) {
    b_d = b_tmp[i];
    d1 = b_tmp[i + 3];
    d2 = b_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_i = i + 3 * i1;
      c_a_contents[b_i] = ((b_d * c_phi[3 * i1] + d1 * c_phi[3 * i1 + 1]) +
                           d2 * c_phi[3 * i1 + 2]) +
                          CR[b_i];
      b_i = i1 + 3 * i;
      a_contents_tmp_tmp = i1 + 6 * i;
      C[a_contents_tmp_tmp] = CT[b_i];
      C[i1 + 6 * (i + 3)] = b_F_tmp[b_i];
      C[a_contents_tmp_tmp + 3] = 0.0;
    }
  }
  for (i = 0; i < 3; i++) {
    b_i = 6 * (i + 3);
    C[b_i + 3] = c_a_contents[3 * i];
    C[b_i + 4] = c_a_contents[3 * i + 1];
    C[b_i + 5] = c_a_contents[3 * i + 2];
  }
  for (i = 0; i < 3; i++) {
    b_d = hphi2[i];
    d1 = hphi2[i + 3];
    d2 = hphi2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      b_a_contents[i + 3 * i1] =
          (b_d * hphi[3 * i1] + d1 * hphi[3 * i1 + 1]) + d2 * hphi[3 * i1 + 2];
    }
  }
  for (i = 0; i < 3; i++) {
    b_d = II[i];
    d1 = II[i + 3];
    d2 = II[i + 6];
    d3 = b_tmp[i];
    d4 = b_tmp[i + 3];
    d5 = b_tmp[i + 6];
    d6 = F_tmp[i];
    d7 = F_tmp[i + 3];
    d8 = F_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      a_contents_tmp_tmp = 3 * i1 + 1;
      i5 = 3 * i1 + 2;
      b_i = i + 3 * i1;
      c_phi[b_i] =
          (b_d * hphi[3 * i1] + d1 * hphi[a_contents_tmp_tmp]) + d2 * hphi[i5];
      hphi2[b_i] = (m * d3 * rhat[3 * i1] + m * d4 * rhat[a_contents_tmp_tmp]) +
                   m * d5 * rhat[i5];
      b_F_tmp[b_i] =
          (d6 * b_a_contents[3 * i1] + d7 * b_a_contents[a_contents_tmp_tmp]) +
          d8 * b_a_contents[i5];
    }
    b_d = hphi2[i];
    d1 = hphi2[i + 3];
    d2 = hphi2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      what[i + 3 * i1] =
          (b_d * wtBody_tmp[3 * i1] + d1 * wtBody_tmp[3 * i1 + 1]) +
          d2 * wtBody_tmp[3 * i1 + 2];
    }
  }
  for (i = 0; i < 3; i++) {
    b_d = b_tmp[i];
    d1 = b_tmp[i + 3];
    d2 = b_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      c_a_contents[i + 3 * i1] =
          (b_d * c_phi[3 * i1] + d1 * c_phi[3 * i1 + 1]) +
          d2 * c_phi[3 * i1 + 2];
      b_i = i1 + 3 * i;
      a_contents_tmp_tmp = i1 + 6 * i;
      M[a_contents_tmp_tmp] = m * (double)iv[b_i];
      M[i1 + 6 * (i + 3)] = b_F_tmp[b_i];
      M[a_contents_tmp_tmp + 3] = what[b_i];
    }
  }
  /*  need B, derivative of force w.r.t input */
  for (i = 0; i < 3; i++) {
    b_d = b_tmp[i];
    d1 = b_tmp[i + 3];
    d2 = b_tmp[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      M[(i1 + 6 * (i + 3)) + 3] = c_a_contents[i1 + 3 * i];
      hphi2[i + 3 * i1] =
          (-m * b_d * rhat[3 * i1] + -m * d1 * rhat[3 * i1 + 1]) +
          -m * d2 * rhat[3 * i1 + 2];
    }
    b_d = hphi2[i];
    d1 = hphi2[i + 3];
    d2 = hphi2[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      what[i + 3 * i1] =
          (b_d * wtBody_tmp[3 * i1] + d1 * wtBody_tmp[3 * i1 + 1]) +
          d2 * wtBody_tmp[3 * i1 + 2];
    }
    b_d = what[i];
    d1 = what[i + 3];
    d2 = what[i + 6];
    for (i1 = 0; i1 < 3; i1++) {
      hphi2[i + 3 * i1] = (b_d * y_tmp[3 * i1] + d1 * y_tmp[3 * i1 + 1]) +
                          d2 * y_tmp[3 * i1 + 2];
      B[i1 + 6 * i] = -m * y_tmp[i1 + 3 * i];
    }
  }
  for (i = 0; i < 3; i++) {
    B[6 * i + 3] = hphi2[3 * i];
    B[6 * i + 4] = hphi2[3 * i + 1];
    B[6 * i + 5] = hphi2[3 * i + 2];
  }
}

/* End of code generation (RigidBodyForce.c) */
