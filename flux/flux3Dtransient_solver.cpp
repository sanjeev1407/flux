// This file is part of the Finite-element soLver for Unsteady electromagnetiX (FLUX)
//
// Copyright (C) 2021 Sanjeev Kumar <sanjeev4@illinois.edu> (University of Illinois at Urbana-Champaign)
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/
//=====================================================================================================

#include "electric_solver.hpp"

#ifdef MFEM_USE_MPI

using namespace std;

namespace mfem
{

using namespace common;

namespace electromagnetics
{

ElectricDiffusionEOperator::ElectricDiffusionEOperator(
   int stateVectorLen,
   ParFiniteElementSpace &HCurlFES,
   ParFiniteElementSpace &HDivFES,
   ParFiniteElementSpace &HGradFES,
   Array<int> &ess_bdr_zero_arg,
   Array<int> &ess_bdr_cosine_arg,
   Array<int> &phi_ess_bdr_zero_arg,
   Array<int> &phi_ess_bdr_cosine_arg,
   double mu_coef,
   std::map<int, double> sigmaAttMap)
   
   : TimeDependentOperator(stateVectorLen, 0.0),
     HCurlFESpace(HCurlFES), HDivFESpace(HDivFES), HGradFESpace(HGradFES),
     a0(NULL), a1(NULL), m1(NULL), s1(NULL), curl(NULL), 
     weakCurl(NULL), grad(NULL), A0(NULL),
     A1(NULL), M1(NULL), X0(NULL), X1(NULL), B0(NULL), B1(NULL), 
     v0(NULL), v1(NULL), ams_a1(NULL), amg_a0(NULL), pcg_a0(NULL),
     pcg_a1(NULL), dsp_m1(NULL), pcg_m1(NULL), 
     mu(mu_coef), dt_A1(-1.0)
   {
   ess_bdr_zero.SetSize(ess_bdr_zero_arg.Size());
   for (int i=0; i<ess_bdr_zero_arg.Size(); i++)
   {
      ess_bdr_zero[i] = ess_bdr_zero_arg[i];
   }

   ess_bdr_cosine.SetSize(ess_bdr_cosine_arg.Size());
   for (int i=0; i<ess_bdr_cosine_arg.Size(); i++)
   {
      ess_bdr_cosine[i] = ess_bdr_cosine_arg[i];
   }

   phi_ess_bdr_zero.SetSize(phi_ess_bdr_zero_arg.Size());
   for (int i=0; i<phi_ess_bdr_zero_arg.Size(); i++)
   {
      phi_ess_bdr_zero[i] = phi_ess_bdr_zero_arg[i];
   }

   phi_ess_bdr_cosine.SetSize(phi_ess_bdr_cosine_arg.Size());
   for (int i=0; i<phi_ess_bdr_cosine_arg.Size(); i++)
   {
      phi_ess_bdr_cosine[i] = phi_ess_bdr_cosine_arg[i];
   }

   
   sigma     = new MeshDependentCoefficient(sigmaAttMap);

   this->buildA0(*sigma);
   this->buildM1(*sigma);
   this->buildS1(1.0/mu);
   this->buildCurl(1.0/mu);
   this->buildGrad();

   v0 = new ParGridFunction(&HGradFESpace);
   v1 = new ParGridFunction(&HCurlFESpace);
   A0 = new HypreParMatrix;
   A1 = new HypreParMatrix;
   X0 = new Vector;
   X1 = new Vector;
   B0 = new Vector;
   B1 = new Vector;

}

void ElectricDiffusionEOperator::Init(Vector &X)
{
   Vector zero_vec(3); zero_vec = 0.0;
   VectorConstantCoefficient Zero_vec(zero_vec);
   ConstantCoefficient Zero(0.0);

   // The big BlockVector stores the fields as follows:
   //    E field
   //    B field
   //    P field

   int Vsize_nd = HCurlFESpace.GetVSize();
   int Vsize_rt = HDivFESpace.GetVSize();
   int Vsize_h1 = HGradFESpace.GetVSize();

   Array<int> true_offset(4);
   true_offset[0] = 0;
   true_offset[1] = true_offset[0] + Vsize_nd;
   true_offset[2] = true_offset[1] + Vsize_rt;
   true_offset[3] = true_offset[2] + Vsize_h1;

   Vector* xptr = (Vector*) &X;
   ParGridFunction E, B, P;
   E.MakeRef(&HCurlFESpace,*xptr,true_offset[0]);
   B.MakeRef(&HDivFESpace, *xptr,true_offset[1]);
   P.MakeRef(&HGradFESpace, *xptr,true_offset[2]);
  
   E.ProjectCoefficient(Zero_vec);
   B.ProjectCoefficient(Zero_vec);
   P.ProjectCoefficient(Zero);
 
}

/*
This is an experimental Mult() method for explicit integration.
Not recommended for actual use.
S0 P  = 0
M1 dEdt = -S1*E + M1*grad P => M1 * E = weakCurl^T B + M1*grad P
   dBdt = -Curl E

Boundary conditions are applied to E and P.  No boundary conditions are applied to B.
*/
void ElectricDiffusionEOperator::Mult(const Vector &X, Vector &dX_dt) const
{
   dX_dt = 0.0;

   // The big BlockVector stores the fields as follows:
   //    E field
   //    B field
   //    P field
-p precice_config.xml -cr 0.08 -cl 0.5 -pdt
   int Vsize_nd = HCurlFESpace.GetVSize();
   int Vsize_rt = HDivFESpace.GetVSize();
   int Vsize_h1 = HGradFESpace.GetVSize();

   Array<int> true_offset(4);
   true_offset[0] = 0;
   true_offset[1] = true_offset[0] + Vsize_nd;
   true_offset[2] = true_offset[1] + Vsize_rt;
   true_offset[3] = true_offset[2] + Vsize_h1;

   Vector* xptr = (Vector*) &X;
   ParGridFunction E, B, P;
   E.MakeRef(&HCurlFESpace,*xptr,true_offset[0]);
   B.MakeRef(&HDivFESpace, *xptr,true_offset[1]);
   P.MakeRef(&HGradFESpace, *xptr,true_offset[2]);

   ParGridFunction dE, dB, dP;
   dE.MakeRef(&HCurlFESpace,dX_dt,true_offset[0]);
   dB.MakeRef(&HDivFESpace, dX_dt,true_offset[1]);
   dP.MakeRef(&HGradFESpace, dX_dt,true_offset[2]);

   // db = - Curl E
   curl->Mult(E, dB);
   dB *= -1.0;

   // form the Laplacian and solve it
   ParGridFunction P_gf(&HGradFESpace);

   // cosine_p_bc is given function defining electrostatic potential on surface
   FunctionCoefficient cosine_voltage(cosine_p_bc);
   cosine_voltage.SetTime(this->GetTime());
   P_gf = 0.0;
   P_gf.ProjectBdrCoefficient(cosine_voltage,phi_ess_bdr_cosine);

   FunctionCoefficient zero_voltage(zero_p_bc);
   P_gf.ProjectBdrCoefficient(zero_voltage,phi_ess_bdr_zero);


   Array<int> phi_ess_tdof_list;
   Array<int> phi_ess_bdr(4);

   for (int i=0; i<4; i++)
     {
       phi_ess_bdr[i] = phi_ess_bdr_zero[i] + phi_ess_bdr_cosine[i];
     }

   HGradFESpace.GetEssentialTrueDofs(phi_ess_bdr, phi_ess_tdof_list);

   *v0 = 0.0;
   a0->FormLinearSystem(phi_ess_tdof_list,P_gf,*v0,*A0,*X0,*B0);

   if (amg_a0 == NULL) { amg_a0 = new HypreBoomerAMG(*A0); }
   if (pcg_a0 == NULL)
     {
       pcg_a0 = new HyprePCG(*A0);
       pcg_a0->SetTol(SOLVER_TOL);
       pcg_a0->SetMaxIter(SOLVER_MAX_IT);
       pcg_a0->SetPrintLevel(SOLVER_PRINT_LEVEL);
       pcg_a0->SetPreconditioner(*amg_a0);
     }

   // pcg "Mult" operation is a solve
   // X0 = A0^-1 * B0
   pcg_a0->Mult(*B0, *X0);

   a0->RecoverFEMSolution(*X0,*v0,P);
   dP = 0.0;

   // Solving for Electric field
   // RHS v1 = <1/mu v, curl u> B + M1 grad P
   weakCurl->MultTranspose(B, *v1);

   grad->Mult(P,E);
   m1->AddMult(E,*v1,1.0);

   // OK now v1 is the right hand side, just need to add essential BC's

   ParGridFunction E_gf(&HCurlFESpace);

   //n x E = 0 Dirichlet BC 
   VectorFunctionCoefficient Ezero(3, e_tan_zero_bc);
   E_gf = 0.0;
   E_gf.ProjectBdrCoefficientTangent(Ezero,ess_bdr_zero);

   //n x E = cos(w*t) dirichlet BC 
   VectorFunctionCoefficient E_cosine(3, e_tan_cosine_bc);
   E_cosine.SetTime(this->GetTime());
   E_gf.ProjectBdrCoefficientTangent(E_cosine,ess_bdr_cosine);

   // apply essential BC's
   // the new system to solve is M1 X1 = B1
   Array<int> ess_tdof_list;
   Array<int> ess_bdr(4);

   for (int i=0; i<4; i++)
     {
       ess_bdr[i] = ess_bdr_zero[i] + ess_bdr_cosine[i];
     }

   HCurlFESpace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);   

   m1->FormLinearSystem(ess_tdof_list,E_gf,*v1,*A1,*X1,*B1);

   if (dsp_m1 == NULL) { dsp_m1 = new HypreDiagScale(*A1); }
   if (pcg_m1 == NULL)
   {
      pcg_m1 = new HyprePCG(*A1);
      pcg_m1->SetTol(SOLVER_TOL);
      pcg_m1->SetMaxIter(SOLVER_MAX_IT);
      pcg_m1->SetPrintLevel(SOLVER_PRINT_LEVEL);
      pcg_m1->SetPreconditioner(*dsp_m1);
   }
   // pcg "Mult" operation is a solve
   // X1 = M1^-1 * B1 = M1^-1 [(curl x 1/mu B) + M1 * grad P]
   pcg_m1->Mult(*B1, *X1);

   // "undo" the static condensation and fill in grid function dE
   m1->RecoverFEMSolution(*X1,*v1,E);
   dE = 0.0;

   // the total field is E_tot = E_ind - Grad Phi
   // so we need to subtract out Grad Phi
   // E = E - grad (P)
   grad->AddMult(P,E,-1.0);
}

/*
This is the main computational code that computes dX/dt implicitly
where X is the state vector containing E, B

        S0 P = 0
(M1+dt S1)E = WeakCurl^T B + M1 * grad P
          dBdt = -Curl E_{n+1}

*/
void ElectricDiffusionEOperator::ImplicitSolve(const double dt,
                                               const Vector &X, Vector &dX_dt)
{
   if ( A1 == NULL || fabs(dt-dt_A1) > 1.0e-12*dt )
   {
      this->buildA1(1.0/mu, *sigma, dt);
   }

   dX_dt = 0.0;

   // The big BlockVector stores the fields as follows:
   //    E field
   //    B field
   //    P field

   int Vsize_nd = HCurlFESpace.GetVSize();
   int Vsize_rt = HDivFESpace.GetVSize();
   int Vsize_h1 = HGradFESpace.GetVSize();

   Array<int> true_offset(4);
   true_offset[0] = 0;
   true_offset[1] = true_offset[0] + Vsize_nd;
   true_offset[2] = true_offset[1] + Vsize_rt;
   true_offset[3] = true_offset[2] + Vsize_h1;

   Vector* xptr  = (Vector*) &X;
   ParGridFunction E, B, P;
   E.MakeRef(&HCurlFESpace,*xptr,true_offset[0]);
   B.MakeRef(&HDivFESpace, *xptr,true_offset[1]);
   P.MakeRef(&HGradFESpace, *xptr,true_offset[2]);

   ParGridFunction dE, dB, dP;
   dE.MakeRef(&HCurlFESpace,dX_dt,true_offset[0]);
   dB.MakeRef(&HDivFESpace, dX_dt,true_offset[1]);
   dP.MakeRef(&HGradFESpace, dX_dt,true_offset[2]);

   // Solving for electrostatuc potential
   // form the Laplacian and solve it
   ParGridFunction P_gf(&HGradFESpace);

   // cosine_p_bc is given function defining electrostatic potential on surface
   FunctionCoefficient cosine_voltage(cosine_p_bc);
   cosine_voltage.SetTime(this->GetTime());
   P_gf = 0.0;
   P_gf.ProjectBdrCoefficient(cosine_voltage,phi_ess_bdr_cosine);

   FunctionCoefficient zero_voltage(zero_p_bc);
   P_gf.ProjectBdrCoefficient(zero_voltage,phi_ess_bdr_zero);

   Array<int> phi_ess_tdof_list;
   Array<int> phi_ess_bdr(4);

   for (int i=0; i<4; i++)
     {
       phi_ess_bdr[i] = phi_ess_bdr_zero[i] + phi_ess_bdr_cosine[i];
     }

   HGradFESpace.GetEssentialTrueDofs(phi_ess_bdr, phi_ess_tdof_list);

   *v0 = 0.0;
   a0->FormLinearSystem(phi_ess_tdof_list,P_gf,*v0,*A0,*X0,*B0);

   if (amg_a0 == NULL) { amg_a0 = new HypreBoomerAMG(*A0); }
   if (pcg_a0 == NULL)
     {
       pcg_a0 = new HyprePCG(*A0);
       pcg_a0->SetTol(SOLVER_TOL);
       pcg_a0->SetMaxIter(SOLVER_MAX_IT);
       pcg_a0->SetPrintLevel(SOLVER_PRINT_LEVEL);
       pcg_a0->SetPreconditioner(*amg_a0);
     }

   // pcg "Mult" operation is a solve
   // X0 = A0^-1 * B0
   pcg_a0->Mult(*B0, *X0);

   a0->RecoverFEMSolution(*X0,*v0,P);
   dP = 0.0;


   // Solving for Electric field
   // RHS : v1 = <1/mu v, curl u> B + M1 * grad P
   weakCurl->MultTranspose(B, *v1);

   /*VectorFunctionCoefficient source(3, current_source);
   source.SetTime(this->GetTime());

   ParLinearForm *b = new LinearForm(&HCurlFESpace);
   b->AddDomainIntegrator(new VectorFEDomainLFIntegrator(source));   
   b->Assemble();*/ 

   grad->Mult(P,E);
   m1->AddMult(E,*v1,1.0);

   //Subtract b from v1
   //v1->Add(-1.0, *b);

   ParGridFunction E_gf(&HCurlFESpace);

   //n x E = 0 Dirichlet BC 
   VectorFunctionCoefficient Ezero(3, e_tan_zero_bc);
   E_gf = 0.0;
   E_gf.ProjectBdrCoefficientTangent(Ezero,ess_bdr_zero);

   //n x E = cos(w*t) dirichlet BC 
   VectorFunctionCoefficient E_cosine(3, e_tan_cosine_bc);
   E_cosine.SetTime(this->GetTime());
   E_gf.ProjectBdrCoefficientTangent(E_cosine,ess_bdr_cosine);

   // apply essential BC's
   // the new system to solve is A1 X1 = B1
   Array<int> ess_tdof_list;
   Array<int> ess_bdr(4);

   for (int i=0; i<4; i++)
    {
      ess_bdr[i] = ess_bdr_zero[i] + ess_bdr_cosine[i];
    }

   HCurlFESpace.GetEssentialTrueDofs(ess_bdr, ess_tdof_list);


   a1->FormLinearSystem(ess_tdof_list,E_gf,*v1,*A1,*X1,*B1);

   // We only need to create the solver and preconditioner once
   if ( ams_a1 == NULL )
   {
     ParFiniteElementSpace *prec_fespace = &HCurlFESpace;
	
      ams_a1 = new HypreAMS(*A1, prec_fespace);
   }
   if ( pcg_a1 == NULL )
   {
      pcg_a1 = new HyprePCG(*A1);
      pcg_a1->SetTol(SOLVER_TOL);
      pcg_a1->SetMaxIter(SOLVER_MAX_IT);
      pcg_a1->SetPrintLevel(SOLVER_PRINT_LEVEL);
      pcg_a1->SetPreconditioner(*ams_a1);
   }
   // solve the system
   // dE = (A1)^-1 [-S1 E + M1 * grad dP] => E = A1^(-1)[weakCurl^T B + M1*grad P]
   pcg_a1->Mult(*B1, *X1);

   // E is a grid function
   a1->RecoverFEMSolution(*X1,*v1,E);
   dE = 0.0;

   // the total field is E_tot = E_ind - Grad Phi
   // so we need to subtract out Grad Phi
   // E = E - grad (P)
   // note grad maps GF to GF
   grad->AddMult(P,E,-1.0);

    // Compute dB/dt = -Curl(E_{n+1})
    curl->Mult(E, dB);
    dB *= -1.0;
  
   
 }
  
   void ElectricDiffusionEOperator::buildA0(MeshDependentCoefficient &Sigma)
   {
     if ( a0 != NULL ) { delete a0; }

     // First create and assemble the bilinear form.  For now we assume the mesh
     // isn't moving, the materials are time independent, and dt is constant. So
     // we only need to do this once.

     // ConstantCoefficient Sigma(sigma);
     a0 = new ParBilinearForm(&HGradFESpace);
     a0->AddDomainIntegrator(new DiffusionIntegrator(Sigma));
     a0->Assemble();

     // Don't finalize or parallel assemble this is done in FormLinearSystem.
   }
  

 void ElectricDiffusionEOperator::buildA1(double muInv,
                                          MeshDependentCoefficient &Sigma,
                                          double dt)
 {
    if ( a1 != NULL ) { delete a1; }

    // First create and assemble the bilinear form.  For now we assume the mesh
    // isn't moving, the materials are time independent, and dt is constant. So
    // we only need to do this once.

    ConstantCoefficient dtMuInv(dt*muInv);
    a1 = new ParBilinearForm(&HCurlFESpace);
    a1->AddDomainIntegrator(new VectorFEMassIntegrator(Sigma));
    a1->AddDomainIntegrator(new CurlCurlIntegrator(dtMuInv));
    a1->Assemble();

    // Don't finalize or parallel assemble this is done in FormLinearSystem.

    dt_A1 = dt;
 }


 void ElectricDiffusionEOperator::buildM1(MeshDependentCoefficient &Sigma)
 {
    if ( m1 != NULL ) { delete m1; }

    m1 = new ParBilinearForm(&HCurlFESpace);
    m1->AddDomainIntegrator(new VectorFEMassIntegrator(Sigma));
    m1->Assemble();

    // Don't finalize or parallel assemble this is done in FormLinearSystem.
 }


 void ElectricDiffusionEOperator::buildS1(double muInv)
 {
    if ( s1 != NULL ) { delete s1; }

    ConstantCoefficient MuInv(muInv);
    s1 = new ParBilinearForm(&HCurlFESpace);
    s1->AddDomainIntegrator(new CurlCurlIntegrator(MuInv));
    s1->Assemble();
 }


 void ElectricDiffusionEOperator::buildCurl(double muInv)
 {
    if ( curl != NULL ) { delete curl; }
   if ( weakCurl != NULL ) { delete weakCurl; }

   curl = new ParDiscreteLinearOperator(&HCurlFESpace, &HDivFESpace);
   curl->AddDomainInterpolator(new CurlInterpolator);
   curl->Assemble();

   ConstantCoefficient MuInv(muInv);
   weakCurl = new ParMixedBilinearForm(&HCurlFESpace, &HDivFESpace);
   weakCurl->AddDomainIntegrator(new VectorFECurlIntegrator(MuInv));
   weakCurl->Assemble();

   // no ParallelAssemble since this will be applied to GridFunctions
}

  void ElectricDiffusionEOperator::buildGrad()
  {
    if ( grad != NULL ) { delete grad; }

    grad = new ParDiscreteLinearOperator(&HGradFESpace, &HCurlFESpace);
    grad->AddDomainInterpolator(new GradientInterpolator());
    grad->Assemble();

    // no ParallelAssemble since this will be applied to GridFunctions
  }


void ElectricDiffusionEOperator::SetTime(const double _t)
{ t = _t; }

ElectricDiffusionEOperator::~ElectricDiffusionEOperator()
{
   if ( ams_a1 != NULL ) { delete ams_a1; }
   if ( pcg_a1 != NULL ) { delete pcg_a1; }

   if ( dsp_m1 != NULL ) { delete dsp_m1; }
   if ( pcg_m1 != NULL ) { delete pcg_m1; }


   if ( curl != NULL ) { delete curl; }
   if ( weakCurl != NULL ) { delete weakCurl; }
   if ( grad != NULL ) { delete grad; }

   if ( a0 != NULL ) { delete a0; }
   if ( a1 != NULL ) { delete a1; }
   if ( m1 != NULL ) { delete m1; }
   if ( s1 != NULL ) { delete s1; }

   if ( A0 != NULL ) { delete A0; }
   if ( X0 != NULL ) { delete X0; }
   if ( B0 != NULL ) { delete B0; }

   if ( A1 != NULL ) { delete A1; }
   if ( X1 != NULL ) { delete X1; }
   if ( B1 != NULL ) { delete B1; }


   if ( v1 != NULL ) { delete v1; }

   if (sigma     != NULL) { delete sigma; }

   delete amg_a0;
   delete pcg_a0;
   delete M1;
   delete v0;
}

void ElectricDiffusionEOperator::Debug(const char *base, double)
{
   {
      hypre_ParCSRMatrixPrint(*A1,"A1_");
      HypreParVector tempB1(A1->GetComm(),A1->N(),B1->GetData(),A1->ColPart());
      tempB1.Print("B1_");
       HypreParVector tempX1(A1->GetComm(),A1->N(),X1->GetData(),A1->ColPart());
       tempX1.Print("X1_");
    }


 }


 MeshDependentCoefficient::MeshDependentCoefficient(
    const std::map<int, double> &inputMap, double scale)
    : Coefficient()
 {
    // make a copy of the magic attribute-value map for later use
    materialMap = new std::map<int, double>(inputMap);
    scaleFactor = scale;
 }

 MeshDependentCoefficient::MeshDependentCoefficient(
    const MeshDependentCoefficient &cloneMe)
    : Coefficient()
 {
    // make a copy of the magic attribute-value map for later use
    materialMap = new std::map<int, double>(*(cloneMe.materialMap));
    scaleFactor = cloneMe.scaleFactor;
 }

 double MeshDependentCoefficient::Eval(ElementTransformation &T,
                                       const IntegrationPoint &ip)
{
   // given the attribute, extract the coefficient value from the map
   std::map<int, double>::iterator it;
   int thisAtt = T.Attribute;
   double value;
   it = materialMap->find(thisAtt);
   if (it != materialMap->end())
   {
      value = it->second;
   }
   else
   {
      value = 0.0; // avoid compile warning
      std::cerr << "MeshDependentCoefficient attribute " << thisAtt
                << " not found" << std::endl;
      mfem_error();
   }

   return value*scaleFactor;
}

ScaledGFCoefficient::ScaledGFCoefficient(GridFunction *gf,
                                         MeshDependentCoefficient &input_mdc)
   : GridFunctionCoefficient(gf), mdc(input_mdc) {}

double ScaledGFCoefficient::Eval(ElementTransformation &T,
                                 const IntegrationPoint &ip)
{
   return mdc.Eval(T,ip) * GridFunctionCoefficient::Eval(T,ip);
}

} // namespace electromagnetics

} // namespace mfem

#endif // MFEM_USE_MPI
