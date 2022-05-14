// This file is part of the Finite-element soLver for Unsteady electromagnetiX
// (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of
// Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//============================================================================

#ifndef FLUX_MFEM_COEFFICIENTS_HPP_
#define FLUX_MFEM_COEFFICIENTS_HPP_

#include <functional>
#include <vector>

#include <mfem/mfem.hpp>

#include "flux/debug_assert.hpp"

class LambdaCoefficient : public mfem::Coefficient {
 public:
  LambdaCoefficient(void) = default;

  LambdaCoefficient(std::function<double(mfem::ElementTransformation& a_T,
                                         const mfem::IntegrationPoint& a_ip)>
                        a_function)
      : action_m(a_function) {}

  void SetLambda(std::function<double(mfem::ElementTransformation& a_T,
                                      const mfem::IntegrationPoint& a_ip)>
                     a_function) {
    action_m = a_function;
  }

  virtual double Eval(mfem::ElementTransformation& a_T,
                      const mfem::IntegrationPoint& a_ip) override final {
    return action_m(a_T, a_ip);
  }

 private:
  std::function<double(mfem::ElementTransformation& a_T,
                       const mfem::IntegrationPoint& a_ip)>
      action_m;
};

class LambdaVectorCoefficient : public mfem::VectorCoefficient {
 public:
  LambdaVectorCoefficient(void) = delete;

  LambdaVectorCoefficient(
      const int a_size,
      std::function<void(mfem::Vector& a_vec, mfem::ElementTransformation& a_T,
                         const mfem::IntegrationPoint& a_ip)>
          a_function)
      : mfem::VectorCoefficient(a_size), action_m(a_function) {}

  void SetLambda(
      const int a_size,
      std::function<void(mfem::Vector& a_vec, mfem::ElementTransformation& a_T,
                         const mfem::IntegrationPoint& a_ip)>
          a_function) {
    vdim = a_size;
    action_m = a_function;
  }

  virtual void Eval(mfem::Vector& a_vec, mfem::ElementTransformation& a_T,
                    const mfem::IntegrationPoint& a_ip) override final {
    a_vec.SetSize(vdim);
    return action_m(a_vec, a_T, a_ip);
  }

 private:
  std::function<void(mfem::Vector& a_vec, mfem::ElementTransformation& a_T,
                     const mfem::IntegrationPoint& a_ip)>
      action_m;
};

class LambdaMatrixCoefficient : public mfem::MatrixCoefficient {
 public:
  LambdaMatrixCoefficient(void) = delete;

  LambdaMatrixCoefficient(
      const int a_dim, std::function<void(mfem::DenseMatrix& a_mat,
                                          mfem::ElementTransformation& a_T,
                                          const mfem::IntegrationPoint& a_ip)>
                           a_function)
      : mfem::MatrixCoefficient(a_dim, false), action_m(a_function) {}

  void SetLambda(const int a_dim,
                 std::function<void(mfem::DenseMatrix& a_mat,
                                    mfem::ElementTransformation& a_T,
                                    const mfem::IntegrationPoint& a_ip)>
                     a_function) {
    height = a_dim;
    width = a_dim;
    action_m = a_function;
  }

  virtual void Eval(mfem::DenseMatrix& a_mat, mfem::ElementTransformation& a_T,
                    const mfem::IntegrationPoint& a_ip) override final {
    a_mat.SetSize(height);
    return action_m(a_mat, a_T, a_ip);
  }

 private:
  std::function<void(mfem::DenseMatrix& a_mat, mfem::ElementTransformation& a_T,
                     const mfem::IntegrationPoint& a_ip)>
      action_m;
};

class MatrixGridFunctionCoefficient : public mfem::MatrixCoefficient {
 public:
  MatrixGridFunctionCoefficient(void) = delete;

  MatrixGridFunctionCoefficient(const int a_dim,
                                const mfem::ParGridFunction* a_gf)
      : mfem::MatrixCoefficient(a_dim, false),
        grid_function_m(a_gf),
        tmp_vector_m(a_dim * a_dim) {}

  virtual void Eval(mfem::DenseMatrix& a_mat, mfem::ElementTransformation& a_T,
                    const mfem::IntegrationPoint& a_ip) override final {
    grid_function_m->GetVectorValue(a_T, a_ip, tmp_vector_m);
    a_mat.Reset(tmp_vector_m.GetData(), height, width);
  }

 private:
  const mfem::ParGridFunction* grid_function_m;
  mfem::Vector tmp_vector_m;
};

class MatrixQuadratureFunctionCoefficient : public mfem::MatrixCoefficient {
 public:
  MatrixQuadratureFunctionCoefficient(void) = delete;

  MatrixQuadratureFunctionCoefficient(const int a_dim,
                                      const mfem::QuadratureFunction* a_qf)
      : mfem::MatrixCoefficient(a_dim, false),
        quadrature_function_m(a_qf),
        tmp_vector_m(a_dim * a_dim) {}

  virtual void Eval(mfem::DenseMatrix& a_mat, mfem::ElementTransformation& a_T,
                    const mfem::IntegrationPoint& a_ip) override final {
    quadrature_function_m->GetElementValues(a_T.ElementNo, a_ip.index,
                                            tmp_vector_m);
    a_mat.Reset(tmp_vector_m.GetData(), height, width);
  }

 private:
  const mfem::QuadratureFunction* quadrature_function_m;
  mfem::Vector tmp_vector_m;
};

#endif  // FLUX_MFEM_COEFFICIENTS_HPP_
