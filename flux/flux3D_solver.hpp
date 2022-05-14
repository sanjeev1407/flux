// This file is part of the Finite-element soLver for Unsteady electromagnetiX (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=====================================================================================================

#ifdef MFEM_USE_MPI

#include <memory>
#include <iostream>
#include <fstream>

namespace mfem
{

  namespace electromagnetics
  {

    // Some global variable for convenience
    const double       SOLVER_TOL = 1.0e-9;
    const int       SOLVER_MAX_IT = 1000;
    // Initialized in electric.cpp and used in electric_solver.cpp:
    const int SOLVER_PRINT_LEVEL = 0;

    // These are defined in emsolve.cpp
    void current_source(const Vector &x, Vector &Js);
    double voltage_bc(const Vector &x);
    double unit_coeff(const Vector &x);
    double sigma_func(const Vector &x);
    double sigma_conductor(const Vector &x);
    double sigma_non_conductor(const Vector &x);
    double sigma_profile(const Vector &x);
    double Pe_func(const Vector &x);
    double ne_func(const Vector &x);
    
    
    // A Coefficient is an object with a function Eval that returns a double.  A
    // MeshDependentCoefficient returns a different value depending upon the given
    // mesh attribute, i.e. a "material property".
    // Somewhat inefficiently, this is achieved using a GridFunction.
    class MeshDependentCoefficient: public Coefficient
    {
    private:
      std::map<int, std::function<double(const Vector &)>> *materialMap;
      double scaleFactor;
    public:
      MeshDependentCoefficient(std::map<int, std::function<double(const Vector &)>> &inputMap,
			       double scale = 1.0);
      MeshDependentCoefficient(MeshDependentCoefficient &cloneMe);
      virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
      void SetScaleFactor(const double &scale) { scaleFactor = scale; }
      virtual ~MeshDependentCoefficient()
      {
	if (materialMap != NULL) { delete materialMap; }
      }
    };

    // A Coefficient is an object with a function Eval that returns a double. The
    // JouleHeatingCoefficient object will contain a reference to the electric field
    // grid function, and the conductivity sigma, and returns sigma E dot E at a
    // point.
    class JouleHeatingCoefficient: public Coefficient
    {
    private:
      ParComplexGridFunction &E;
      MeshDependentCoefficient sigma_coefficient;
    public:
      JouleHeatingCoefficient(MeshDependentCoefficient &sigma_coefficient_,
			      ParComplexGridFunction &E_)
	: E(E_), sigma_coefficient(sigma_coefficient_) {}
      virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
      virtual ~JouleHeatingCoefficient() {}
    };

    class AmbipolarCoefficient: public VectorCoefficient
    {
    private:
      ParGridFunction &grad_Pe;
      ParGridFunction &ne;
    public:
      AmbipolarCoefficient(int dim, ParGridFunction &grad_Pe_,
			   ParGridFunction &ne_)
	: VectorCoefficient(dim), grad_Pe(grad_Pe_), ne(ne_) {}
      virtual void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip);
      virtual ~AmbipolarCoefficient() {}
      };

    class sigmaEreCoefficient: public VectorCoefficient
    {
    private:
      MeshDependentCoefficient sigma_coefficient;
      ParComplexGridFunction &E;
      
    public:
      sigmaEreCoefficient(int dim, MeshDependentCoefficient sigma_coefficient_, ParComplexGridFunction &E_)
	: VectorCoefficient(dim), sigma_coefficient(sigma_coefficient_), E(E_) {}
      virtual void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip);
      virtual ~sigmaEreCoefficient() {}
      };

    class sigmaEimCoefficient: public VectorCoefficient
    {
    private:
      MeshDependentCoefficient sigma_coefficient;
      ParComplexGridFunction &E;
      
    public:
      sigmaEimCoefficient(int dim, MeshDependentCoefficient sigma_coefficient_, ParComplexGridFunction &E_)
	: VectorCoefficient(dim), sigma_coefficient(sigma_coefficient_), E(E_) {}
      virtual void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip);
      virtual ~sigmaEimCoefficient() {}
      };

    class LorentzForceCoefficient: public VectorCoefficient
    {
    private:
      MeshDependentCoefficient sigma_coefficient;
      ParComplexGridFunction &E;
      ParGridFunction &B_re;
      ParGridFunction &B_im;
      
    public:
      LorentzForceCoefficient(int dim, MeshDependentCoefficient sigma_coefficient_, ParComplexGridFunction &E_, ParGridFunction &B_re_, ParGridFunction &B_im_)
	: VectorCoefficient(dim), sigma_coefficient(sigma_coefficient_), E(E_), B_re(B_re_), B_im(B_im_) {}
      virtual void Eval(Vector &V, ElementTransformation &T, const IntegrationPoint &ip);
      virtual ~LorentzForceCoefficient() {}
      };
    
    
  } // namespace electromagnetics

} // namespace mfem

#endif // MFEM_USE_MPI

