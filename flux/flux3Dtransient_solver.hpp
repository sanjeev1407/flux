// This file is part of the Finite-element soLver for Unsteady electromagnetiX (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=====================================================================================================

#ifndef MFEM_ELECTRIC_SOLVER
#define MFEM_ELECTRIC_SOLVER

#include "../common/pfem_extras.hpp"

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
extern int SOLVER_PRINT_LEVEL;
extern int        STATIC_COND;

// These are defined in electric.cpp
void e_tan_zero_bc(const Vector &x, Vector &E);
void e_tan_cosine_bc(const Vector &x, double t, Vector &E);
void current_source(const Vector &x, double t, Vector &Js);
double cosine_p_bc(const Vector &x, double t);
double zero_p_bc(const Vector &x, double t);
void e_exact(const Vector &x, double t, Vector &E);
void b_exact(const Vector &x, double t, Vector &B);

// A Coefficient is an object with a function Eval that returns a double.  A
// MeshDependentCoefficient returns a different value depending upon the given
// mesh attribute, i.e. a "material property".
// Somewhat inefficiently, this is achieved using a GridFunction.
class MeshDependentCoefficient: public Coefficient
{
private:
   std::map<int, double> *materialMap;
   double scaleFactor;
public:
   MeshDependentCoefficient(const std::map<int, double> &inputMap,
                            double scale = 1.0);
   MeshDependentCoefficient(const MeshDependentCoefficient &cloneMe);
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   void SetScaleFactor(const double &scale) { scaleFactor = scale; }
   virtual ~MeshDependentCoefficient()
   {
      if (materialMap != NULL) { delete materialMap; }
   }
};

// This Coefficient is a product of a GridFunction and MeshDependentCoefficient
// for example if T (temperature) is a GridFunction and c (heat capacity) is a
// MeshDependentCoefficient, this function can compute c*T.
class ScaledGFCoefficient: public GridFunctionCoefficient
{
private:
   MeshDependentCoefficient mdc;
public:
   ScaledGFCoefficient(GridFunction *gf, MeshDependentCoefficient &input_mdc);
   virtual double Eval(ElementTransformation &T, const IntegrationPoint &ip);
   void SetMDC(const MeshDependentCoefficient &input_mdc) { mdc = input_mdc; }
   virtual ~ScaledGFCoefficient() {}
};

/**
   After spatial discretization, the magnetic diffusion equation can be written
   as a system of ODEs:

      dE/dt       = - (M1(sigma) + dt S1(1/mu))^{-1}*(S1(1/mu)*E + M1(sigma)*grad dP/dt )
      dB/dt       = - Curl(E)
      
   where

     E is the 1-form electric field,
     B is the 2-form magnetic flux,
     P is the 0-form electrostatic potential
   
   Class ElectricDiffusionEOperator represents the right-hand side of the above
   system of ODEs.
*/
class ElectricDiffusionEOperator : public TimeDependentOperator
{
protected:
   // These ParFiniteElementSpace objects provide degree-of-freedom mappings.
   // To create these you must provide the mesh and the definition of the FE
   // space. These objects are used to create hypre vectors to store the DOFs,
   // they are used to create grid functions to perform FEM interpolation, and
   // they are used by bilinear forms.
   ParFiniteElementSpace &HCurlFESpace;
   ParFiniteElementSpace &HDivFESpace;
   ParFiniteElementSpace &HGradFESpace;
   
   // ParBilinearForms are used to create sparse matrices representing discrete
   // linear operators.
   ParBilinearForm *a0,*a1, *m1, *s1;
   ParDiscreteLinearOperator  *curl, *grad;
   ParMixedBilinearForm  *weakCurl;
   ParLinearForm *b;

  // Hypre matrices and vectors for 1-form systems A1 X1 = B1
   HypreParMatrix *A0, *A1, *M1;
   Vector *X0, *X1, *B0, *B1;

   // temporary work vectors
   ParGridFunction *v0, *v1;

   // HypreSolver is derived from Solver, which is derived from Operator. So a
   // HypreSolver object has a Mult() operator, which is actually the solver
   // operation y = A^-1 x i.e. multiplication by A^-1.
   // HyprePCG is a wrapper for hypre's preconditioned conjugate gradient.
   mutable HypreSolver * amg_a0;
   mutable HyprePCG    * pcg_a0;
   mutable HypreSolver * ams_a1;
   mutable HyprePCG    * pcg_a1;
   mutable HypreSolver * dsp_m1;
   mutable HyprePCG    * pcg_m1;
   
   mutable Array<int> ess_bdr_zero;
   mutable Array<int> ess_bdr_zero_vdofs;
   mutable Array<int> ess_bdr_cosine;
   mutable Array<int> ess_bdr_cosine_vdofs;
   mutable Array<int> phi_ess_bdr_zero;
   mutable Array<int> phi_ess_bdr_zero_vdofs;
   mutable Array<int> phi_ess_bdr_cosine;
   mutable Array<int> phi_ess_bdr_cosine_vdofs;
  
   
   MeshDependentCoefficient *sigma;
   double mu, dt_A1;

   // The method builA2 creates the ParBilinearForm a2, the HypreParMatrix A2,
   // and the solver and preconditioner pcg_a2 and amg_a2. The other build
   // functions do similar things.
   void buildA0(MeshDependentCoefficient &sigma);
   void buildA1(double muInv, MeshDependentCoefficient &sigma, double dt);
   void buildM1(MeshDependentCoefficient &sigma);
   void buildS1(double muInv);
   void buildCurl(double muInv);
   void buildGrad();
   
public:
   ElectricDiffusionEOperator(int len,
                              ParFiniteElementSpace &HCurlFES,
                              ParFiniteElementSpace &HDivFES,
			      ParFiniteElementSpace &HGradFES,
			      Array<int> &ess_bdr_zero,
			      Array<int> &ess_bdr_cosine,
			      Array<int> &phi_ess_bdr_zero,
			      Array<int> &phi_ess_bdr_cosine,
                              double mu,
                              std::map<int, double> sigmaAttMap
                              );

   // Initialize the fields. This is where restart would go to.
   void Init(Vector &vx);

   // class TimeDependentOperator is derived from Operator, and class Operator
   // has the virtual function Mult(x,y) which computes y = A x for some matrix
   // A, or more generally a nonlinear operator y = A(x).
   virtual void Mult(const Vector &vx, Vector &dvx_dt) const;

   // Solve the Backward-Euler equation: k = f(x + dt*k, t), for the unknown
   // slope k. This is the only requirement for high-order SDIRK implicit
   // integration. This is a virtual function of class TimeDependentOperator.
   virtual void ImplicitSolve(const double dt, const Vector &x, Vector &k);

   
   void SetTime(const double _t);

   // Write all the hypre matrices and vectors to disk.
   void Debug(const char *basefilename, double time);

   virtual ~ElectricDiffusionEOperator();
};


} // namespace electromagnetics

} // namespace mfem

#endif // MFEM_USE_MPI

#endif // MFEM_JOULE_SOLVER
