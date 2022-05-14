// This file is part of the Finite-element soLver for Unsteady electromagnetiX
// (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of
// Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=================================================================================

#ifndef FLUX_FLUX_2D_SOLVER_HPP_
#define FLUX_FLUX_2D_SOLVER_HPP_

#include <fstream>
#include <iostream>
#include <memory>
#include "mfem.hpp"

using namespace mfem;

// Some global variable for convenience
const double SOLVER_TOL = 1.0e-6;
const int SOLVER_MAX_IT = 1000;
const int SOLVER_PRINT_LEVEL = 0;

// These are defined in emsolve.cpp
// If these are defined in flux2D.cpp,
// they should be moved into flux2D.hpp
double current_source(const Vector& x);
double voltage_bc(const Vector& x);
double unit_coeff(const Vector& x);
double sigmaFunc(const Vector& x);
double w_mu_r_sigmaFunc(const Vector& x);
double sigma_conductor(const Vector& x);
double sigma_non_conductor(const Vector& x);
double rFunc(const Vector& x);
double rinvFunc(const Vector& x);
double two_pi_rFunc(const Vector& x);


// A Coefficient is an object with a function Eval that returns a double. The
// JouleHeatingCoefficient object will contain a reference to the electric field
// grid function, and the conductivity sigma, and returns sigma E dot E at a
// point.
class JouleHeatingCoefficient : public Coefficient {
 private:
  GridFunctionCoefficient& sigma_gf_coef;
  ParComplexGridFunction& E;

 public:
  JouleHeatingCoefficient(GridFunctionCoefficient& sigma_gf_coef_,
                          ParComplexGridFunction& E_)
      : sigma_gf_coef(sigma_gf_coef_), E(E_) {}
  virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip);
  virtual ~JouleHeatingCoefficient() {}
};

class TotalJouleHeatingCoefficient : public Coefficient {
 private:
  FunctionCoefficient& two_pi_r_coef;
  GridFunctionCoefficient& sigma_gf_coef;
  ParComplexGridFunction& E;

 public:
  TotalJouleHeatingCoefficient(FunctionCoefficient& two_pi_r_coef_,
                               GridFunctionCoefficient& sigma_gf_coef_,
                               ParComplexGridFunction& E_)
      : two_pi_r_coef(two_pi_r_coef_), sigma_gf_coef(sigma_gf_coef_), E(E_) {}
  virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip);
  virtual ~TotalJouleHeatingCoefficient() {}
};

class rEreCoefficient : public Coefficient {
 private:
  FunctionCoefficient& r_coef;
  ParComplexGridFunction& E;

 public:
  rEreCoefficient(FunctionCoefficient& r_coef_, ParComplexGridFunction& E_)
      : r_coef(r_coef_), E(E_) {}
  virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip);
  virtual ~rEreCoefficient() {}
};

class rEimCoefficient : public Coefficient {
 private:
  FunctionCoefficient& r_coef;
  ParComplexGridFunction& E;

 public:
  rEimCoefficient(FunctionCoefficient& r_coef_, ParComplexGridFunction& E_)
      : r_coef(r_coef_), E(E_) {}
  virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip);
  virtual ~rEimCoefficient() {}
};

class BreCoefficient : public VectorCoefficient {
 private:
  FunctionCoefficient& r_inv_coef;
  ParGridFunction& grad_E_im_gf;
  ParGridFunction& grad_rEim_gf;

 public:
  BreCoefficient(int dim, FunctionCoefficient& r_inv_coef_,
                 ParGridFunction& grad_E_im_gf_, ParGridFunction& grad_rEim_gf_)
      : VectorCoefficient(dim),
        r_inv_coef(r_inv_coef_),
        grad_E_im_gf(grad_E_im_gf_),
        grad_rEim_gf(grad_rEim_gf_) {}
  virtual void Eval(Vector& V, ElementTransformation& T,
                    const IntegrationPoint& ip);
  virtual ~BreCoefficient() {}
};

class BimCoefficient : public VectorCoefficient {
 private:
  FunctionCoefficient& r_inv_coef;
  ParGridFunction& grad_E_re_gf;
  ParGridFunction& grad_rEre_gf;

 public:
  BimCoefficient(int dim, FunctionCoefficient& r_inv_coef_,
                 ParGridFunction& grad_E_re_gf_, ParGridFunction& grad_rEre_gf_)
      : VectorCoefficient(dim),
        r_inv_coef(r_inv_coef_),
        grad_E_re_gf(grad_E_re_gf_),
        grad_rEre_gf(grad_rEre_gf_) {}
  virtual void Eval(Vector& V, ElementTransformation& T,
                    const IntegrationPoint& ip);
  virtual ~BimCoefficient() {}
};

class AxialLorentzForceCoefficient : public Coefficient {
 private:
  GridFunctionCoefficient& sigma_gf_coef;
  ParComplexGridFunction& E;
  ParGridFunction& B_re;
  ParGridFunction& B_im;

 public:
  AxialLorentzForceCoefficient(GridFunctionCoefficient& sigma_gf_coef_,
                               ParComplexGridFunction& E_,
                               ParGridFunction& B_re_, ParGridFunction& B_im_)
      : sigma_gf_coef(sigma_gf_coef_), E(E_), B_re(B_re_), B_im(B_im_) {}
  virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip);
  virtual ~AxialLorentzForceCoefficient() {}
};

class RadialLorentzForceCoefficient : public Coefficient {
 private:
  GridFunctionCoefficient& sigma_gf_coef;
  ParComplexGridFunction& E;
  ParGridFunction& B_re;
  ParGridFunction& B_im;

 public:
  RadialLorentzForceCoefficient(GridFunctionCoefficient& sigma_gf_coef_,
                                ParComplexGridFunction& E_,
                                ParGridFunction& B_re_, ParGridFunction& B_im_)
      : sigma_gf_coef(sigma_gf_coef_), E(E_), B_re(B_re_), B_im(B_im_) {}
  virtual double Eval(ElementTransformation& T, const IntegrationPoint& ip);
  virtual ~RadialLorentzForceCoefficient() {}
};


#endif  // FLUX_FLUX_2D_SOLVER_HPP_
