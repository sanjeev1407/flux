// This file is part of the Finite-element soLver for Unsteady electromagnetiX
// (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of
// Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=====================================================================================================

#include "flux/flux2D_solver.hpp"

using namespace std;
using namespace mfem;

double JouleHeatingCoefficient::Eval(ElementTransformation& T,
                                     const IntegrationPoint& ip) {
  double E_re, E_im;
  double thisSigma_coefficient;
  E_re = E.real().GetValue(T, ip);
  E_im = E.imag().GetValue(T, ip);
  thisSigma_coefficient = sigma_gf_coef.Eval(T, ip);
  return 0.5 * thisSigma_coefficient * (E_re * E_re + E_im * E_im);
}

double TotalJouleHeatingCoefficient::Eval(ElementTransformation& T,
                                          const IntegrationPoint& ip) {
  double E_re, E_im;
  double thisSigma_coefficient;
  double twopir;
  E_re = E.real().GetValue(T, ip);
  E_im = E.imag().GetValue(T, ip);
  thisSigma_coefficient = sigma_gf_coef.Eval(T, ip);
  twopir = two_pi_r_coef.Eval(T, ip);

  return twopir * 0.5 * thisSigma_coefficient * (E_re * E_re + E_im * E_im);
}

double rEreCoefficient::Eval(ElementTransformation& T,
                             const IntegrationPoint& ip) {
  double E_re;
  double r;
  E_re = E.real().GetValue(T, ip);
  r = r_coef.Eval(T, ip);
  return r * E_re;
}

double rEimCoefficient::Eval(ElementTransformation& T,
                             const IntegrationPoint& ip) {
  double E_im;
  double r;
  E_im = E.imag().GetValue(T, ip);
  r = r_coef.Eval(T, ip);
  return r * E_im;
}

void BreCoefficient::Eval(Vector& V, ElementTransformation& T,
                          const IntegrationPoint& ip) {
  Vector grad_E_im, grad_rE_im;
  double inv_r;
  grad_E_im_gf.GetVectorValue(T, ip, grad_E_im);
  grad_rEim_gf.GetVectorValue(T, ip, grad_rE_im);
  inv_r = r_inv_coef.Eval(T, ip);
  V[0] = -inv_r * grad_rE_im[1];
  V[1] = grad_E_im[0];

  if (inv_r > 1.0e6) {
    V[0] = 0;
  }

  if (inv_r > 1.0e6) {
    V[1] = 0;
  }
}

void BimCoefficient::Eval(Vector& V, ElementTransformation& T,
                          const IntegrationPoint& ip) {
  Vector grad_E_re, grad_rE_re;
  double inv_r;
  grad_E_re_gf.GetVectorValue(T, ip, grad_E_re);
  grad_rEre_gf.GetVectorValue(T, ip, grad_rE_re);
  inv_r = r_inv_coef.Eval(T, ip);
  V[0] = inv_r * grad_rE_re[1];
  V[1] = -grad_E_re[0];

  if (inv_r > 1.0e6) {
    V[0] = 0;
  }

  if (inv_r > 1.0e6) {
    V[1] = 0;
  }
}

double AxialLorentzForceCoefficient::Eval(ElementTransformation& T,
                                          const IntegrationPoint& ip) {
  double sig;
  double Ereal, Eimag;
  Vector Breal, Bimag;
  Ereal = E.real().GetValue(T, ip);
  Eimag = E.imag().GetValue(T, ip);
  B_re.GetVectorValue(T, ip, Breal);
  B_im.GetVectorValue(T, ip, Bimag);
  sig = sigma_gf_coef.Eval(T, ip);

  // Axial force
  return -0.5 * sig * (Eimag * Bimag[1] + Ereal * Breal[1]);
}

double RadialLorentzForceCoefficient::Eval(ElementTransformation& T,
                                           const IntegrationPoint& ip) {
  double sig;
  double Ereal, Eimag;
  Vector Breal, Bimag;
  Ereal = E.real().GetValue(T, ip);
  Eimag = E.imag().GetValue(T, ip);
  B_re.GetVectorValue(T, ip, Breal);
  B_im.GetVectorValue(T, ip, Bimag);
  sig = sigma_gf_coef.Eval(T, ip);

  // Radial force
  return 0.5 * sig * (Ereal * Breal[0] + Eimag * Bimag[0]);
}

